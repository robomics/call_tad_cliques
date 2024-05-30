// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT


workflow NCHG {

    take:
        sample_sheet
        domains
        mad_max
        bad_bin_fraction

        fdr
        log_ratio

        zstd_compression_lvl

    main:

        sample_sheet
            .splitCsv(sep: "\t", header: true)
            .map { row -> tuple(row.sample,
                                file(row.hic_file, checkIfExists: true),
                                row.resolution)
            }
            .join(domains,
                  failOnDuplicate: true,
                  failOnMismatch: true)
            .set { nchg_tasks }

        COMPUTE(
            nchg_tasks,
            mad_max,
            bad_bin_fraction
        )

       MERGE(
            COMPUTE.out.parquet
       )

       FILTER(
            MERGE.out.parquet,
            fdr,
            log_ratio
       )

       VIEW(
           FILTER.out.parquet,
           zstd_compression_lvl
       )

    emit:
        tsv = VIEW.out.tsv

}


process COMPUTE {
    label 'process_long'
    label 'process_very_high'
    tag "$sample"

    cpus 16

    input:
        tuple val(sample),
              path(hic),
              val(resolution),
              path(domains)

        val mad_max
        val bad_bin_fraction

    output:
        tuple val(sample),
              path("*.chrom.sizes"),
              path("*.parquet"),
        emit: parquet

    shell:
        '''
        set -o pipefail

        zstdcat '!{domains}' > domains.bed

        NCHG compute \\
            '!{hic}' \\
            '!{sample}' \\
            --resolution='!{resolution}' \\
            --domains=domains.bed \\
            --mad-max='!{mad_max}' \\
            --bad-bin-fraction='!{bad_bin_fraction}' \\
            --threads='!{task.cpus}'
        '''
}

process MERGE {
    tag "$sample"

    cpus 2

    input:
        tuple val(sample),
              path(chrom_sizes),
              path(parquets)

    output:
        tuple val(sample),
              path(outname),
        emit: parquet

    shell:
        outname="${sample}.parquet"
        '''
        NCHG merge '!{sample}' '!{outname}' \\
            --threads='!{task.cpus}'
        '''
}

process FILTER {
    tag "$sample"

    cpus 2

    input:
        tuple val(sample),
              path(parquet)

        val fdr
        val log_ratio

    output:
        tuple val(sample),
              path("*.parquet"),
        emit: parquet

    shell:
        '''
        NCHG filter \\
            '!{parquet}' \\
            '!{sample}.filtered.parquet' \\
            --fdr='!{fdr}' \\
            --log-ratio='!{log_ratio}' \\
            --threads='!{task.cpus}'
        '''
}

process VIEW {
    publishDir "${params.publish_dir}/nchg/",
        enabled: !!params.publish_dir,
        mode: params.publish_dir_mode

    label 'process_medium'
    tag "$sample"

    input:
        tuple val(sample),
              path(parquet)

        val zstd_compression_lvl

    output:
        tuple val(sample),
              path("*.tsv.zst"),
        emit: tsv

    shell:
        '''
        set -o pipefail

        NCHG view \\
            '!{parquet}' |
            zstd '-!{zstd_compression_lvl}' \\
                 -T!{task.cpus} \\
            > '!{sample}.filtered.tsv.zst'
        '''
}
