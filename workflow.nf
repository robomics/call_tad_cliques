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

    extract_chrom_sizes_from_cooler(trans_coolers)
    generate_bed_mask(file(params.assembly_gaps), file(params.cytoband))

    chrom_sizes = extract_chrom_sizes_from_cooler.out.chrom_sizes
    mask = generate_bed_mask.out.bed.collect()

    tads_to_domains(chrom_sizes, file(params.tads))
    domains = bedtools_bed_setdiff(tads_to_domains.out.bed, mask)

    map_intrachrom_interactions(cis_coolers, domains)
    collect_interchrom_interactions(trans_coolers, mask)

    select_significant_intrachrom_interactions(map_intrachrom_interactions.out.bedpe,
                                               params.cis_bin_size,
                                               params.cis_log_ratio,
                                               params.cis_fdr)

    select_significant_interchrom_interactions(collect_interchrom_interactions.out.bedpe,
                                               params.trans_log_ratio,
                                               params.trans_fdr)

    map_interchrom_interactions_to_domains(domains,
                                           select_significant_interchrom_interactions.out.bedpe)

    // TODO need to pair samples before processing
    // merge_interactions()
}

process extract_chrom_sizes_from_cooler {
    label 'very_short'

    input:
        tuple val(sample), path(cooler), val(resolution)

    output:
        val sample, emit: sample
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
    publishDir params.output_dir, mode: 'copy'
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

process tads_to_domains {
    publishDir params.output_dir, mode: 'copy'
    label 'very_short'

    input:
        path chrom_sizes
        path tads

    output:
        path "domains.bed.zst", emit: bed

    shell:
        '''
        set -o pipefail

        '!{params.script_dir}/convert_tads_to_domains.py' \
            '!{chrom_sizes}' \
            '!{tads}'        |
            zstd -T!{task.cpus}                  \
                 -!{params.zstd_compression_lvl} \
                 -o domains.bed.zst
        '''
}

process map_intrachrom_interactions {
    publishDir params.output_dir, mode: 'copy'
    label 'very_short'

    input:
        tuple val(sample), path(cooler), val(resolution)
        path domains

    output:
        path "*.bedpe.zst", emit: bedpe

    shell:
        outname="${sample}_domains_with_cis_interactions.bedpe.zst"
        '''
        set -o pipefail

        if [ !{resolution} -eq 0 ]; then
            cooler='!{cooler}'
        else
            cooler='!{cooler}::/resolutions/!{resolution}'
        fi

        '!{params.script_dir}/map_intrachrom_interactions_to_domains.py' \
            "$cooler"         \
            '!{domains}'      |
            zstd -T!{task.cpus}                  \
                 -!{params.zstd_compression_lvl} \
                 -o '!{outname}'
        '''
}

process bedtools_bedpe_setdiff {
    publishDir params.output_dir, mode: 'copy'
    label 'very_short'

    input:
        path bedpe
        path mask

    output:
        path "*_filtered.bedpe.zst", emit: bedpe

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
    publishDir params.output_dir, mode: 'copy'
    label 'very_short'

    input:
        path bed
        path mask

    output:
        path "*_filtered.bed.zst", emit: bed

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
    publishDir params.output_dir, mode: 'copy'

    input:
        tuple val(sample), path(cooler), val(resolution)
        path mask

    output:
        path "*.bedpe.zst", emit: bedpe

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
    publishDir params.output_dir, mode: 'copy'
    label 'very_short'

    input:
        path bedpe
        val resolution
        val log_ratio
        val fdr

    output:
        path "*_cis_significant.bedpe.zst", emit: bedpe

    shell:
        outname=bedpe.toString().replace(".bedpe.zst", "_cis_significant.bedpe.zst")
        '''
        set -o pipefail

        '!{params.script_dir}/run_nchg.py' \
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
    publishDir params.output_dir, mode: 'copy'

    input:
        path bedpe
        val log_ratio
        val fdr

    output:
        path "*_trans_significant.bedpe.zst", emit: bedpe

    shell:
        outname=bedpe.toString().replace(".bedpe.zst", "_trans_significant.bedpe.zst")
        '''
        set -o pipefail

        '!{params.script_dir}/run_nchg.py' \
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
    publishDir params.output_dir, mode: 'copy'

    input:
        path domains
        path interchrom_interactions

    output:
        path "*.bedpe.zst", emit: bedpe

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
    publishDir params.output_dir, mode: 'copy'

    input:
        path cis
        path trans

    output:
        path "*_interactions.bedpe.zst", emit: bedpe

    shell:
        outname=cis.toString().replace(".bed.zst", "_interactions.bedpe.zst")
        '''
        zstdcat *.bedpe.zst                  |
        awk '$1!="." && $4!="."'             |
        sort -u -k1,1V -k2,2n -k4,4V -k5,5n  |
        zstd -T!{task.cpus}                  \
             -!{params.zstd_compression_lvl} \
             -o '!{outname}'
        '''
}