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
            APPLY_NORMALIZATION.out
        )

        COPY(
            tasks.copy
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

    output:
        tuple val(id),
              path("*_tads.bed.gz"),
        emit: bed

    shell:
        outprefix="${cooler.baseName}"
        '''
        set -o pipefail

        NUMEXPR_MAX_THREADS='!{task.cpus}'
        export NUMEXPR_MAX_THREADS

        hicFindTADs -p '!{task.cpus}' \\
                    --matrix '!{cooler}' \\
                    --outPrefix '!{outprefix}' \\
                    --correctForMultipleTesting fdr

        pigz -9 -p !{task.cpus} < '!{outprefix}_domains.bed' > '!{outprefix}_tads.bed.gz'
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

    output:
        tuple val(id),
              path("*_tads.bed.gz"),
        emit: bed

    shell:
        outprefix="${file.baseName}"
        '''
        set -o pipefail

        zcat -f '!{tads}' |
        pigz -9 -p !{task.cpus} > '!{outprefix}_tads.bed.gz'
        '''
}
