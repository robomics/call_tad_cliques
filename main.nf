#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    cis_coolers = Channel.fromPath(params.mcools)
                         .map { tuple(it.getBaseName(), it, params.cis_bin_size) }
    trans_coolers = Channel.fromPath(params.mcools)
                           .map { tuple(it.getBaseName(), it, params.trans_bin_size) }

    if (params.tads) {
        if (file(params.tads, glob: true) instanceof java.util.LinkedList) {
            // TADs are expected to be in the same order as the cooler files
            tads = Channel.fromPath(params.tads)
	                      .merge(cis_coolers)
                          .map { tuple(it[1], it[0]) }
        } else {
            // A single TAD file was provided: use this annotation for all .cool files
            tads = cis_coolers.map { tuple(it[0], file(params.tads)) }
        }
    } else {
        hicexplorer_find_tads(cis_coolers)
        tads = hicexplorer_find_tads.out.bed
    }

    extract_chrom_sizes_from_cooler(trans_coolers.first())
    generate_bed_mask(file(params.assembly_gaps), file(params.cytoband))

    chrom_sizes = extract_chrom_sizes_from_cooler.out.chrom_sizes
    mask = generate_bed_mask.out.bed

    tads_to_domains(tads, chrom_sizes)
    bedtools_bed_setdiff(tads_to_domains.out.bed, mask)
    domains = bedtools_bed_setdiff.out.bed

    map_intrachrom_interactions(
        cis_coolers.join(domains,
                         failOnMismatch: true,
                         failOnDuplicate: true))

    collect_interchrom_interactions(trans_coolers, mask)

    select_significant_intrachrom_interactions(map_intrachrom_interactions.out.bedpe,
                                               params.cis_bin_size,
                                               params.cis_log_ratio,
                                               params.cis_fdr)

    select_significant_interchrom_interactions(collect_interchrom_interactions.out.bedpe,
                                               params.trans_log_ratio,
                                               params.trans_fdr)

    map_interchrom_interactions_to_domains(select_significant_interchrom_interactions
              .out
              .bedpe
              .join(domains,
                    failOnMismatch: true,
                    failOnDuplicate: true)
        )

    merge_interactions(select_significant_intrachrom_interactions
              .out
              .bedpe
              .join(map_interchrom_interactions_to_domains.out.bedpe,
                    failOnMismatch: true,
                    failOnDuplicate: true)
        )

    call_cliques(domains
              .join(merge_interactions.out.bedpe,
                    failOnMismatch: true,
                    failOnDuplicate: true),
                 params.clique_size_thresh)
}

process extract_chrom_sizes_from_cooler {
    label 'very_short'

    input:
        tuple val(sample), path(cooler), val(resolution)

    output:
        path "*.chrom.sizes", emit: chrom_sizes

    shell:
        outname="${sample}.chrom.sizes"
        '''
        set -o pipefail

        if [ !{resolution} -eq 0 ]; then
            cooler='!{cooler}'
        else
            cooler='!{cooler}::/resolutions/!{resolution}'
        fi

        cooler dump --table chroms "$cooler" > '!{outname}'
        '''
}

process generate_bed_mask {
    label 'very_short'

    input:
        path gaps
        path cytoband

    output:
        path "mask.bed.zst", emit: bed

    shell:
        '''
        set -o pipefail

        # Concatenate and sort intervals from gaps and cytoband files,
        # then merge overlapping intervals

        mkfifo gaps.bed.tmp cytoband.bed.tmp
        trap 'rm -f mkfifo gaps.bed.tmp cytoband.bed.tmp' EXIT

        # Decompress BED files and select BED3 columns
        zcat '!{gaps}' | cut -f 2-4 >> gaps.bed.tmp &
        zcat '!{cytoband}' | cut -f 1-3 | grep 'acen$' >> cytoband.bed.tmp &

        # Concatenate, sort and merge intervals
        cat *.bed.tmp |
        sort -k1,1V -k2,2n --parallel=!{task.cpus} |
        bedtools merge -i stdin |
        zstd -T!{task.cpus}                  \
             -!{params.zstd_compression_lvl} \
             -o mask.bed.zst
        '''
}

process hicexplorer_find_tads {
    label 'process_medium'

    input:
        tuple val(sample), path(cooler), val(resolution)

    output:
        tuple val(sample), path("*_tads.bed.zst"), emit: bed

    shell:
        outprefix="${cooler.baseName}"
        '''
        if [ '!{resolution}' -ne 0 ]; then
            cooler='!{cooler}::/resolutions/!{resolution}'
        else
            cooler='!{cooler}'
        fi

        NUMEXPR_MAX_THREADS='!{task.cpus}'
        export NUMEXPR_MAX_THREADS

        hicFindTADs -p '!{task.cpus}'          \
                    --matrix "$cooler"         \
                    --outPrefix '!{outprefix}' \
                    --correctForMultipleTesting fdr

        zstd -T!{task.cpus}                  \
             -!{params.zstd_compression_lvl} \
             '!{outprefix}_domains.bed'

        mv '!{outprefix}_domains.bed.zst' '!{outprefix}_tads.bed.zst'
        '''
}

process tads_to_domains {
    label 'very_short'

    input:
        tuple val(sample), path(tads)
        path chrom_sizes

    output:
        tuple val(sample), path("*_domains.bed.zst"), emit: bed

    shell:
        outname="${sample}_domains.bed.zst"
        '''
        set -o pipefail

        convert_tads_to_domains.py  \
            '!{chrom_sizes}'        \
            '!{tads}'               |
            zstd -T!{task.cpus}                  \
                 -!{params.zstd_compression_lvl} \
                 -o '!{outname}'
        '''
}

process map_intrachrom_interactions {
    label 'very_short'

    input:
        tuple val(sample), path(cooler), val(resolution), path(domains)

    output:
        tuple val(sample), path("*.bedpe.zst"), emit: bedpe

    shell:
        outname="${sample}_domains_with_cis_interactions.bedpe.zst"
        '''
        set -o pipefail

        if [ !{resolution} -eq 0 ]; then
            cooler='!{cooler}'
        else
            cooler='!{cooler}::/resolutions/!{resolution}'
        fi

        map_intrachrom_interactions_to_domains.py  \
            "$cooler"                              \
            '!{domains}'                           |
            zstd -T!{task.cpus}                    \
                 -!{params.zstd_compression_lvl}   \
                 -o '!{outname}'
        '''
}

process bedtools_bedpe_setdiff {
    label 'very_short'

    input:
        tuple val(sample), path(bedpe)
        path mask

    output:
        tuple val(sample), path("*_filtered.bedpe.zst"), emit: bedpe

    shell:
        outname=bed.toString().replace(".bedpe.zst", "_filtered.bedpe.zst")
        '''
        set -o pipefail

        bedtools pairtobed -a <(zstdcat '!{bedpe}') \
                           -b <(zstdcat '!{mask}')  \
                           -type neither            |
            zstd -T!{task.cpus}                     \
                 -!{params.zstd_compression_lvl}    \
                 -o '!{outname}'
        '''
}

process bedtools_bed_setdiff {
    label 'very_short'

    input:
        tuple val(sample), path(bed)
        path mask

    output:
        tuple val(sample), path("*_filtered.bed.zst"), emit: bed

    shell:
        outname=bed.toString().replace(".bed.zst", "_filtered.bed.zst")
        '''
        set -o pipefail

        bedtools subtract -A                      \
                          -a <(zstdcat '!{bed}')  \
                          -b <(zstdcat '!{mask}') |
            zstd -T!{task.cpus}                   \
                 -!{params.zstd_compression_lvl}  \
                 -o '!{outname}'
        '''
}

process collect_interchrom_interactions {
    label 'long'

    input:
        tuple val(sample), path(cooler), val(resolution)
        path mask

    output:
        tuple val(sample), path("*.bedpe.zst"), emit: bedpe

    shell:
        outname="${cooler.baseName}_trans_interactions.bedpe.zst"
        '''
        set -o pipefail

        if [ !{resolution} -eq 0 ]; then
            cooler='!{cooler}'
        else
            cooler='!{cooler}::/resolutions/!{resolution}'
        fi

        cooler dump --join "$cooler" |
            awk '$1!=$4{print}'      |
            bedtools pairtobed -a stdin                 \
                               -b <(zstdcat '!{mask}')  \
                               -type neither            |
            zstd -T!{task.cpus}                  \
                 -!{params.zstd_compression_lvl} \
                 -o '!{outname}'
        '''
}

process select_significant_intrachrom_interactions {
    label 'short'

    input:
        tuple val(sample), path(bedpe)
        val resolution
        val log_ratio
        val fdr

    output:
        tuple val(sample), path("*_cis_significant.bedpe.zst"), emit: bedpe

    shell:
        outname=bedpe.toString().replace(".bedpe.zst", "_cis_significant.bedpe.zst")
        '''
        set -o pipefail

        run_nchg.py                        \
            --fdr='!{fdr}'                 \
            --log-ratio='!{log_ratio}'     \
            --resolution='!{resolution}'   \
            --drop-not-significant         \
            <(zstdcat '!{bedpe}')          \
            intra                          |
        cut -f 1-6,8 |
        zstd -T!{task.cpus}                  \
             -!{params.zstd_compression_lvl} \
             -o '!{outname}'
        '''
}

process select_significant_interchrom_interactions {
    label 'short'

    input:
        tuple val(sample), path(bedpe)
        val log_ratio
        val fdr

    output:
        tuple val(sample), path("*_trans_significant.bedpe.zst"), emit: bedpe

    shell:
        outname=bedpe.toString().replace(".bedpe.zst", "_trans_significant.bedpe.zst")
        '''
        set -o pipefail

        run_nchg.py                        \
            --fdr='!{fdr}'                 \
            --log-ratio='!{log_ratio}'     \
            --drop-not-significant         \
            <(zstdcat '!{bedpe}')          \
            inter                          |
        cut -f 1-6,8 |
        zstd -T!{task.cpus}                  \
             -!{params.zstd_compression_lvl} \
             -o '!{outname}'
        '''
}

process map_interchrom_interactions_to_domains {
    label 'short'

    input:
        tuple val(sample), path(interchrom_interactions), path(domains)

    output:
        tuple val(sample), path("*.bedpe.zst"), emit: bedpe

    shell:
        outname=domains.toString().replace(".bed.zst", "_overlapping_domains.bedpe.zst")
        '''
        set -o pipefail

        mkfifo {left,right}.fifo

        zstdcat '!{interchrom_interactions}' |
            awk '{printf("%s\\t%s\\t%s\\n", $1,($2+$3)/2,1+($2+$3)/2)}' |
            bedtools intersect -wao -a stdin -b <(zstdcat '!{domains}') |
            cut -f 4-6 >> left.fifo &

        zstdcat '!{interchrom_interactions}' |
            awk '{printf("%s\\t%s\\t%s\\t%s\\n", $4,($5+$6)/2,1+($5+$6)/2,$7)}' |
            bedtools intersect -wao -a stdin -b <(zstdcat '!{domains}')         |
            awk '{printf("%s\\t%s\\t%s\\t%s\\n", $5,$6,$7,$4)}' >> right.fifo &

        # We're discarding rows without match because the domains passed as input to this step
        # do not cover the entire genome
        paste left.fifo right.fifo                  |
            awk '$1!="." && $4!="."'                |
            sort -u -k1,1V -k2,2n -k4,4V -k5,5n     |
        zstd -T!{task.cpus}                         \
             -!{params.zstd_compression_lvl}        \
             -o '!{outname}'
        '''
}

process merge_interactions {
    label 'very_short'

    input:
        tuple val(sample), path(cis), path(trans)

    output:
        tuple val(sample), path("*_interactions.bedpe.zst"), emit: bedpe

    shell:
        outname="${sample}_domains_with_significant_interactions.bedpe.zst"
        '''
        set -o pipefail

        zstdcat *.bedpe.zst                  |
        awk '$1!="." && $4!="."'             |
        sort -u -k1,1V -k2,2n -k4,4V -k5,5n  |
        zstd -T!{task.cpus}                  \
             -!{params.zstd_compression_lvl} \
             -o '!{outname}'
        '''
}

process call_cliques {
    publishDir params.outdir, mode: 'copy'

    label 'very_short'

    input:
        tuple val(sample), path(domains), path(significant_interactions)
        val clique_size

    output:
        tuple val(sample), path("*_clique_interactions.bedpe"), emit: clique_interactions
        tuple val(sample), path("*_clique_sizes.bedGraph"), emit: clique_sizes
        tuple val(sample), path("*_clique_stats.tsv"), emit: clique_stats
        tuple val(sample), path("*_tad_interactions.bedpe"), emit: tad_interactions

    shell:
        outprefix="${sample}"
        '''
        call_cliques.py                   \
            '!{domains}'                  \
            '!{significant_interactions}' \
            '!{outprefix}'
        '''
}