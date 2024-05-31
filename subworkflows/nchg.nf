// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// Workaround for optional input files: https://github.com/nextflow-io/nextflow/issues/1694
def make_optional_input(path) {
    if (path?.trim()) {
        return [file(path)]
    }
    return []
}

workflow NCHG {

    take:
        sample_sheet
        domains
        mad_max
        bad_bin_fraction

        fdr
        log_ratio

        cytoband
        gaps
        custom_mask

    main:

        GENERATE_MASK(
            make_optional_input(cytoband),
            make_optional_input(gaps),
            make_optional_input(custom_mask),
            params.zstd_compression_lvl
        )

        MASK_DOMAINS(
            domains,
            GENERATE_MASK.out.bed,
            params.zstd_compression_lvl
        )

        MASK_DOMAINS.out.set { filtered_domains }

        sample_sheet
            .splitCsv(sep: "\t", header: true)
            .map { row -> tuple(row.sample,
                                file(row.hic_file, checkIfExists: true),
                                row.resolution)
            }
            .join(filtered_domains,
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
           FILTER.out.parquet
       )

    emit:
        tsv = VIEW.out.tsv

}


process GENERATE_MASK {
    label 'duration_very_short'

    cpus 1

    input:
        path cytoband
        path gaps
        path other
        val zstd_compression_lvl

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

        if [[ '!{other}' == "" ]]; then
            ln -s empty.tmp other
        else
            ln -s '!{other}' other
        fi

        # Concatenate and sort intervals from gaps and cytoband files,
        # then merge overlapping intervals

        mkfifo gaps.bed.tmp cytoband.bed.tmp other.bed.tmp
        trap 'rm -f mkfifo gaps.bed.tmp cytoband.bed.tmp other.bed.tmp empty.tmp' EXIT

        # Decompress BED files and select BED3 columns
        zcat -f gaps | cut -f 2-4 >> gaps.bed.tmp &
        zcat -f cytoband | grep 'acen$' | cut -f 1-3 >> cytoband.bed.tmp &
        zcat -f other | cut -f 1-3 >> other.bed.tmp &

        # Concatenate, sort and merge intervals
        cat *.bed.tmp |
        sort -k1,1V -k2,2n --parallel=!{task.cpus} |
        bedtools merge -i stdin |
        zstd -T!{task.cpus} \\
             -!{zstd_compression_lvl} \\
             -o mask.bed.zst
        '''
}

process MASK_DOMAINS {
    label 'duration_very_short'
    tag "$sample"

    cpus 1

    input:
        tuple val(sample),
              path(domains)

        path mask
        val zstd_compression_lvl

    output:
        tuple val(sample),
              path(outname)

    shell:
        outname=domains.toString().replace(".bed.gz", "_filtered.bed.zst")
        '''
        set -o pipefail

        bedtools subtract -A \\
                          -a <(zcat -f '!{domains}') \\
                          -b <(zstdcat -f '!{mask}') |
            zstd -T!{task.cpus} \\
                 -!{zstd_compression_lvl} \\
                 -o '!{outname}'
        '''
}


process COMPUTE {
    label 'process_long'
    label 'process_high'
    tag "$sample"


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

        zstdcat -f '!{domains}' > domains.bed

        NCHG compute \\
            '!{hic}' \\
            '!{sample}' \\
            --resolution='!{resolution}' \\
            --domains=domains.bed \\
            --mad-max='!{mad_max}' \\
            --bad-bin-fraction='!{bad_bin_fraction}' \\
            --threads='!{task.cpus}'

        ls -lah
        exit 1
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

    output:
        tuple val(sample),
              path("*.tsv.gz"),
        emit: tsv

    shell:
        '''
        set -o pipefail

        NCHG view '!{parquet}' |
            pigz -9 -p !{task.cpus} > '!{sample}.filtered.tsv.gz'
        '''
}
