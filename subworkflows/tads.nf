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


workflow TADS {

    take:
        sample_sheet
        hic_norm
        cool_norm
        zstd_compression_lvl

    main:

        sample_sheet
            .splitCsv(sep: "\t", header: true)
            .map { row ->
                   tuple(row.sample,
                         file(row.hic_file, checkIfExists: true),
                         row.resolution,
                         make_optional_input(row.tads))
            }
            .branch { call: !it[3]
                      copy: !!it[3]
            }
            .set { tasks }

        SELECT_NORMALIZATION_METHOD(
            tasks.call.map { tuple(it[0], it[1]) },
            hic_norm,
            cool_norm
        )

        tasks.call
            .join(SELECT_NORMALIZATION_METHOD.out.norm)
            .set { hicexplorer_tasks }



        APPLY_NORMALIZATION(
            hicexplorer_tasks
        )

        HICEXPLORER_FIND_TADS(
            APPLY_NORMALIZATION.out,
            zstd_compression_lvl
        )

        COPY(
            tasks.copy,
            zstd_compression_lvl
        )

    emit:
        tsv = HICEXPLORER_FIND_TADS.out.bed.mix(COPY.out.bed)

}


process SELECT_NORMALIZATION_METHOD {
    label 'process_very_short'
    tag "$id"

    input:
        tuple val(id),
              path(file)

        val hic_norm
        val cool_norm

    output:
        tuple val(id),
              stdout,
        emit: norm

    shell:
        '''
        #!/usr/bin/env python3

        import hictkpy

        if hictkpy.is_hic("!{file}"):
            print("!{hic_norm}", end="")
        else:
            print("!{cool_norm}", end="")
        '''
}


process APPLY_NORMALIZATION {
    tag "$id ($resolution; $normalization)"

    input:
        tuple val(id),
              path(file),
              val(resolution),
              path(tads),
              val(normalization)

    output:
        tuple val(id),
              path(outname)

    shell:
        outname="${id}.cool"
        '''
        apply_normalization.py '!{file}' \\
                               '!{outname}' \\
                               --norm-name='!{normalization}' \\
                               --resolution='!{resolution}'
        '''
}

process HICEXPLORER_FIND_TADS {
    publishDir "${params.publish_dir}/tads",
        enabled: !!params.publish_dir,
        mode: params.publish_dir_mode

    label 'process_medium'
    tag "$id"

    input:
        tuple val(id),
              path(cooler)

        val zstd_compression_lvl

    output:
        tuple val(id),
              path("*_tads.bed.zst"),
        emit: bed

    shell:
        outprefix="${cooler.baseName}"
        '''
        NUMEXPR_MAX_THREADS='!{task.cpus}'
        export NUMEXPR_MAX_THREADS

        hicFindTADs -p '!{task.cpus}' \\
                    --matrix '!{cooler}' \\
                    --outPrefix '!{outprefix}' \\
                    --correctForMultipleTesting fdr

        zstd -T!{task.cpus} \\
             -'!{zstd_compression_lvl}' \\
             '!{outprefix}_domains.bed'

        mv '!{outprefix}_domains.bed.zst' '!{outprefix}_tads.bed.zst'
        '''
}

process COPY {
    publishDir "${params.publish_dir}/tads",
        enabled: !!params.publish_dir,
        mode: params.publish_dir_mode

    label 'process_very_short'
    tag "$id"

    input:
        tuple val(id),
              path(file),
              val(resolution),
              path(tads)

        val zstd_compression_lvl

    output:
        tuple val(id),
              path("*_tads.bed.zst"),
        emit: bed

    shell:
        outprefix="${file.baseName}"
        '''
        zcat -f '!{tads}' |
        zstd -T!{task.cpus} \\
             -'!{zstd_compression_lvl}' \\
             > '!{outprefix}_tads.bed.zst'
        '''
}
