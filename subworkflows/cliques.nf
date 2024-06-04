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

workflow CLIQUES {

    take:
        significant_interactions
        tads
        clique_size_thresh


    main:

        interaction_types = ["cis-only", "trans-only", "all"]

        significant_interactions
            .join(tads)
            .combine(interaction_types)
            .set { tasks }

        CALL(
            tasks,
            clique_size_thresh
        )

        CALL.out.cliques
            .map { tuple(it[1], it[3]) }
            .groupTuple()
            .set { cliques }

        PLOT_MAXIMAL_CLIQUE_SIZE(
            cliques
        )

}


process CALL {
    publishDir "${params.publish_dir}/cliques",
        enabled: !!params.publish_dir,
        mode: params.publish_dir_mode
    tag "$id"

    cpus 1

    input:
        tuple val(id),
              path(significant_interactions),
              path(tads),
              val(interaction_type)

        val min_clique_size

    output:
        tuple val(id),
              val(interaction_type),
              path("*.bed.gz"),
              path("*.tsv.gz"),
        emit: cliques

    shell:
        suffix=interaction_type.replaceAll('-', '_').replaceAll('_only', '')
        outprefix="${id}_${suffix}"
        '''
        call_cliques.py                                \\
            '!{tads}'                                  \\
            '!{significant_interactions}'              \\
            '!{outprefix}'                             \\
            --interaction-type='!{interaction_type}'   \\
            --clique-size-threshold=!{min_clique_size}

        gzip -9 *.{bed,tsv}
        '''
}

process PLOT_MAXIMAL_CLIQUE_SIZE {
    publishDir "${params.publish_dir}/plots/cliques",
        enabled: !!params.publish_dir,
        mode: params.publish_dir_mode

    label 'duration_very_short'
    tag "$interaction_type"

    cpus 1

    input:
        tuple val(interaction_type),
              path(cliques)

    output:
        tuple val(interaction_type),
              path("*.png"), emit: png

        tuple val(interaction_type),
              path("*.svg"), emit: svg

    shell:
        '''
        plot_maximal_clique_sizes.py \\
            *_cliques.tsv.gz         \\
            -o '!{interaction_type}_maximal_clique_size_abs'

        plot_maximal_clique_sizes.py \\
            *_cliques.tsv.gz         \\
            --stat='density'         \\
            -o '!{interaction_type}_maximal_clique_size_rel'
        '''
}
