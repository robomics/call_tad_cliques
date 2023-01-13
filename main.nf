#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

def strip_resolution_from_cooler_uri(uri) {
    file(uri.replaceFirst(/::\/resolutions\/\d+$/, ""), checkIfExists: true)
}

// Workaround for optional input files: https://github.com/nextflow-io/nextflow/issues/1694
def make_optional_input(path) {
    if (path?.trim()) {
        return [file(path)]
    }
    return []
}

def parse_sample_sheet_row(row) {
    cooler_cis_fname = strip_resolution_from_cooler_uri(row.cooler_cis)
    cooler_trans_fname = strip_resolution_from_cooler_uri(row.cooler_trans)
    if (cooler_cis_fname == cooler_trans_fname) {
        // Make file optional to workaround input file name collisions
        coolers = [cooler_cis_fname]
    } else {
        coolers = [cooler_cis_fname, cooler_trans_fname]
    }

    tads = make_optional_input(row.tads)

    tuple(row.sample,
          coolers,
          row.cooler_cis,
          row.cooler_trans,
          tads)
}

def row_to_tuple(row) {
    tuple(row.id,
          row.cooler_cis,
          row.cooler_trans,
          row.tads,
          row.cis_resolution,
          row.trans_resolution)
}


workflow {
    if (params.sample_sheet) {
        check_sample_sheet(file(params.sample_sheet, checkIfExists: true))
    } else {
        generate_sample_sheet(params.sample,
                              params.cooler_cis,
                              params.cooler_trans,
                              params.tads ? params.tads : "")
        check_sample_sheet(generate_sample_sheet.out.tsv)
    }
    sample_sheet = check_sample_sheet.out.tsv

    sample_sheet.splitCsv(sep: "\t", header: true)
                .map {
                        it = parse_sample_sheet_row(it)
                        it[1] + it[4]  // Concatenate path to coolers and tads (when available)
                    }
                .set { files_from_sample_sheet }

    process_sample_sheet(sample_sheet, files_from_sample_sheet.collect())

    process_sample_sheet.out.tsv
           .splitCsv(sep: "\t", header: true)
           .map { row -> tuple(row.id, row.cooler_cis, row.cis_resolution) }
           .set { cis_coolers }

    process_sample_sheet.out.tsv
           .splitCsv(sep: "\t", header: true)
           .map { row -> tuple(row.id, row.cooler_trans, row.trans_resolution) }
           .set { trans_coolers }

    extract_chrom_sizes_from_cooler(cis_coolers.first())
    generate_bed_mask(make_optional_input(params.assembly_gaps),
                      make_optional_input(params.cytoband))

    chrom_sizes = extract_chrom_sizes_from_cooler.out.chrom_sizes
    mask = generate_bed_mask.out.bed

    process_tads(
        process_sample_sheet.out.tsv
            .splitCsv(sep: "\t", header: true)
            .map { row ->
                   tuple(row.id, row.cooler_cis, row.cis_resolution,
                         make_optional_input(row.tads))
                 }
    )


    fill_gaps_between_tads(process_tads.out.bed.map { tuple(it[0], it[1]) },
                           chrom_sizes)

    bedtools_bed_setdiff(fill_gaps_between_tads.out.bed,
                         mask)

    process_sample_sheet.out.tsv
            .splitCsv(sep: "\t", header: true)
            .map { row_to_tuple(it) }
            .join(bedtools_bed_setdiff.out.bed)
            .map { // replace old (empty or unfiltered) tads with new ones
                   it.set(3, it[6])
                   it[0..5]
                 }
            .set { sample_table }

    // Sample table schema:
    // 0. id
    // 1. cis_cooler
    // 2. trans_cooler
    // 3. tads
    // 4. cis_resolution
    // 5. trans_resolution

    sample_table.map { it -> tuple(it[0], it[3]) }
                .set { tads }

    // Process intra-chromosomal interactions
    map_intrachrom_interactions(
        // id, cis_cooler, tads, cis_resolution
        sample_table.map { tuple(it[0], it[1], it[3], it[4]) }
    )

    select_significant_intrachrom_interactions(map_intrachrom_interactions.out.bedpe,
                                               params.cis_log_ratio,
                                               params.cis_fdr)

    // Process inter-chromosomal interactions
    collect_interchrom_interactions(trans_coolers,
                                    mask)

    select_significant_interchrom_interactions(collect_interchrom_interactions.out.bedpe,
                                               params.trans_log_ratio,
                                               params.trans_fdr)

    map_interchrom_interactions_to_tads(
        select_significant_interchrom_interactions.out.bedpe
            .join(tads)
            // id, significant_interactions, tads
            .map { tuple(it[0], it[1], it[3]) }
    )

    // Merge inter and intra-chromosomal interactions
    merge_interactions(
        select_significant_intrachrom_interactions.out.bedpe
            .join(map_interchrom_interactions_to_tads.out.bedpe)
            // id, cis_interactions, trans_interactions
            .map { tuple(it[0], it[1], it[3]) }
    )

    call_cliques(tads.join(merge_interactions.out.bedpe),
                 params.clique_size_thresh)
}

process generate_sample_sheet {
    label 'very_short'

    cpus 1

    input:
        val sample
        val cooler_cis
        val cooler_trans
        val tads

    output:
        path "sample_sheet.tsv", emit: tsv

    shell:
        '''
        printf 'sample\\tcooler_cis\\tcooler_trans\\ttads\\n' > sample_sheet.tsv
        printf '%s\\t%s\\t%s\\t%s\\n' '!{sample}' \
                                      '!{cooler_cis}' '!{cooler_trans}' \
                                      '!{tads}' >> sample_sheet.tsv
        '''
}

process check_sample_sheet {
    label 'very_short'

    cpus 1

    input:
        path sample_sheet

    output:
        path "${sample_sheet}", includeInputs: true, emit: tsv

    shell:
        '''
        parse_samplesheet.py --detached '!{sample_sheet}' > /dev/null
        '''
}

process process_sample_sheet {
    label 'very_short'

    cpus 1

    input:
        path sample_sheet
        path files

    output:
        path "*.ok", emit: tsv

    shell:
        '''
        parse_samplesheet.py '!{sample_sheet}' > '!{sample_sheet}.ok'
        '''
}

process extract_chrom_sizes_from_cooler {
    // publishDir "${params.outdir}/dbg", mode: 'copy'
    label 'very_short'

    cpus 1

    input:
        tuple val(id),
              path(cooler),
              val(resolution)

    output:
        path "*.chrom.sizes", emit: chrom_sizes

    shell:
        outname="${id}.chrom.sizes"
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
    // publishDir "${params.outdir}/dbg", mode: 'copy'
    label 'very_short'

    cpus 1

    input:
        path gaps
        path cytoband

    output:
        path "mask.bed.zst", emit: bed

    shell:
        '''
        set -o pipefail

        # Implement logic to support optional gaps/cytoband input files
        touch empty.tmp

        if [[ '!{gaps}' == "" ]]; then
            ln -s empty.tmp gaps
        else
            ln -s '!{gaps}' gaps
        fi

        if [[ '!{cytoband}' == "" ]]; then
            ln -s empty.tmp cytoband
        else
            ln -s '!{cytoband}' cytoband
        fi

        # Concatenate and sort intervals from gaps and cytoband files,
        # then merge overlapping intervals

        mkfifo gaps.bed.tmp cytoband.bed.tmp
        trap 'rm -f mkfifo gaps.bed.tmp cytoband.bed.tmp empty.tmp' EXIT

        # Decompress BED files and select BED3 columns
        zcat -f gaps | cut -f 2-4 >> gaps.bed.tmp &
        zcat -f cytoband | cut -f 1-3 | grep 'acen$' >> cytoband.bed.tmp &

        # Concatenate, sort and merge intervals
        cat *.bed.tmp |
        sort -k1,1V -k2,2n --parallel=!{task.cpus} |
        bedtools merge -i stdin |
        zstd -T!{task.cpus}                  \
             -!{params.zstd_compression_lvl} \
             -o mask.bed.zst
        '''
}

process bedtools_bed_setdiff {
    // publishDir "${params.outdir}/dbg", mode: 'copy'
    label 'very_short'

    cpus 1

    input:
        tuple val(id),
              path(bed)
        path mask

    output:
        tuple val(id),
              path("*_filtered.bed.zst"),
        emit: bed

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

process process_tads {
    // publishDir "${params.outdir}/dbg", mode: 'copy'
    label 'process_medium'
    
    cpus 1

    input:
        tuple val(id),
              path(cooler),
              val(resolution),
              path(tads)

    output:
        tuple val(id),
              path("*_tads.bed.zst"),
              val(resolution),
        emit: bed

    shell:
        outprefix="${cooler.baseName}"
        '''
        NUMEXPR_MAX_THREADS='!{task.cpus}'
        export NUMEXPR_MAX_THREADS

        if [[ '!{tads}' == "" ]]; then
            if [ '!{resolution}' -ne 0 ]; then
                cooler='!{cooler}::/resolutions/!{resolution}'
            else
                cooler='!{cooler}'
            fi

            hicFindTADs -p '!{task.cpus}'          \
                        --matrix "$cooler"         \
                        --outPrefix '!{outprefix}' \
                        --correctForMultipleTesting fdr
        else
            zcat -f '!{tads}' > '!{outprefix}_domains.bed'
        fi

        zstd -T!{task.cpus}                  \
             -!{params.zstd_compression_lvl} \
             '!{outprefix}_domains.bed'

        mv '!{outprefix}_domains.bed.zst' '!{outprefix}_tads.bed.zst'
        '''
}

process fill_gaps_between_tads {
    // publishDir "${params.outdir}/dbg", mode: 'copy'
    label 'very_short'

    cpus 1

    input:
        tuple val(id),
              path(tads)
        path chrom_sizes

    output:
        tuple val(id),
              path("*_tads.bed.zst"),
        emit: bed

    shell:
        outname="${id}_tads.bed.zst"
        '''
        set -o pipefail

        fill_gaps_between_tads.py  \
            '!{chrom_sizes}'        \
            '!{tads}'               |
            zstd -T!{task.cpus}                  \
                 -!{params.zstd_compression_lvl} \
                 -o '!{outname}'
        '''
}

process map_intrachrom_interactions {
    // publishDir "${params.outdir}/dbg", mode: 'copy'
    label 'very_short'

    cpus 1

    input:
        tuple val(id),
              path(cooler),
              path(tads),
              val(resolution)

    output:
        tuple val(id),
              path("*.bedpe.zst"),
              val(resolution),
        emit: bedpe

    shell:
        outname="${id}_tads_with_cis_interactions.bedpe.zst"
        '''
        set -o pipefail

        if [ !{resolution} -eq 0 ]; then
            cooler='!{cooler}'
        else
            cooler='!{cooler}::/resolutions/!{resolution}'
        fi

        map_intrachrom_interactions_to_tads.py     \
            "$cooler"                              \
            '!{tads}'                              |
            zstd -T!{task.cpus}                    \
                 -!{params.zstd_compression_lvl}   \
                 -o '!{outname}'
        '''
}

process select_significant_intrachrom_interactions {
    // publishDir "${params.outdir}/dbg", mode: 'copy'
    label 'short'

    cpus 1

    input:
        tuple val(id),
              path(bedpe),
              val(resolution)
        val log_ratio
        val fdr

    output:
        tuple val(id),
              path("*_significant_cis_interactions.bedpe.zst"),
              val(resolution),
        emit: bedpe

    shell:
        outname=bedpe.toString()
                     .replace("_cis_interactions.bedpe.zst",
                              "_significant_cis_interactions.bedpe.zst")
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

process collect_interchrom_interactions {
    // publishDir "${params.outdir}/dbg", mode: 'copy'
    label 'long'

    cpus 1

    input:
        tuple val(id),
              path(cooler),
              val(resolution)
        path mask

    output:
        tuple val(id),
              path("*.bedpe.zst"),
              val(resolution),
        emit: bedpe

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

process select_significant_interchrom_interactions {
    // publishDir "${params.outdir}/dbg", mode: 'copy'
    label 'short'

    cpus 1

    input:
        tuple val(id),
              path(bedpe),
              val(resolution)
        val log_ratio
        val fdr

    output:
        tuple val(id),
              path("*_significant_trans_significant.bedpe.zst"),
              val(resolution),
        emit: bedpe

    shell:
        outname=bedpe.toString()
                     .replace("_trans_interactions.bedpe.zst",
                              "_significant_trans_significant.bedpe.zst")
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

process map_interchrom_interactions_to_tads {
    // publishDir "${params.outdir}/dbg", mode: 'copy'
    label 'short'

    cpus 1

    input:
        tuple val(id),
              path(interchrom_interactions),
              path(tads)

    output:
        tuple val(id),
              path("*.bedpe.zst"),
        emit: bedpe

    shell:
        outname=interchrom_interactions.toString()
                                       .replace(".bedpe.zst",
                                                "_overlapping_tads.bedpe.zst")
        '''
        set -o pipefail

        mkfifo {left,right}.fifo

        zstdcat '!{interchrom_interactions}' |
            awk '{printf("%s\\t%s\\t%s\\n", $1,($2+$3)/2,1+($2+$3)/2)}' |
            bedtools intersect -wao -a stdin -b <(zstdcat '!{tads}')    |
            cut -f 4-6 >> left.fifo &

        zstdcat '!{interchrom_interactions}' |
            awk '{printf("%s\\t%s\\t%s\\t%s\\n", $4,($5+$6)/2,1+($5+$6)/2,$7)}' |
            bedtools intersect -wao -a stdin -b <(zstdcat '!{tads}')            |
            awk '{printf("%s\\t%s\\t%s\\t%s\\n", $5,$6,$7,$4)}' >> right.fifo &

        # We're discarding rows without match because the tads passed as input to this step
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
    // publishDir "${params.outdir}/dbg", mode: 'copy'
    label 'very_short'

    cpus 1

    input:
        tuple val(id),
              path(cis_bedpe),
              path(trans_bedpe)

    output:
        tuple val(id),
              path("*_interactions.bedpe.zst"),
        emit: bedpe

    shell:
        outname="${id}_tads_with_significant_interactions.bedpe.zst"
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
    label 'short'
    
    cpus 1

    input:
        tuple val(id),
              path(tads),
              path(significant_interactions)
        val min_clique_size

    output:
        tuple val(id), path("*_clique_interactions.bedpe"), emit: clique_interactions
        tuple val(id), path("*_clique_sizes.bedGraph"), emit: clique_sizes
        tuple val(id), path("*_clique_stats.tsv"), emit: clique_stats
        tuple val(id), path("*_tad_interactions.bedpe"), emit: tad_interactions

    shell:
        outprefix="${id}"
        '''
        call_cliques.py                   \
            '!{tads}'                     \
            '!{significant_interactions}' \
            '!{outprefix}'                \
            --clique-size-threshold=!{min_clique_size}
        '''
}
