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
        mask
        clique_size_thresh


    main:

        interaction_types = []
        if (params.call_cis_cliques) {
            interaction_types.add("cis")
        }
        if (params.call_trans_cliques) {
            interaction_types.add("trans")
        }

        significant_interactions
            .join(tads)
            .combine(interaction_types)
            .set { tasks }

        CALL(
            tasks,
            clique_size_thresh
        )

        CALL.out.cliques.join(mask)
            .set { mask_tasks }

        MASK(
            mask_tasks
        )

        MASK.out.cliques
            .map { tuple(it[1], it[3]) }
            .groupTuple()
            .set { cliques }

        PLOT_MAXIMAL_CLIQUE_SIZE_DISTRIBUTION_BY_TAD(
            cliques
        )

        PLOT_CLIQUE_SIZE_DISTRIBUTION(
            cliques
        )

}


process CALL {
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
        outprefix="${id}_${interaction_type}.unfiltered"
        '''
        call_cliques.py \\
            '!{tads}' \\
            '!{significant_interactions}' \\
            '!{outprefix}' \\
            --interaction-type='!{interaction_type}-only' \\
            --clique-size-threshold=!{min_clique_size}

        gzip -9 *.{bed,tsv}
        '''
}

process MASK {
    publishDir "${params.publish_dir}/cliques",
        enabled: !!params.publish_dir,
        mode: params.publish_dir_mode
    tag "$id"

    cpus 1

    input:
        tuple val(id),
              val(interaction_type),
              path(domains),
              path(cliques),
              path(mask)

    output:
        tuple val(id),
              val(interaction_type),
              path("*.bed.gz"),
              path("*.tsv.gz"),
        emit: cliques

    shell:
        opts=[]
        if (!mask.toString().isEmpty()) {
            opts.push("--mask='${mask}'")
        }

        opts=opts.join(" ")
        outprefix="${id}_${interaction_type}"
        '''
        set -o pipefail

        mask_cliques.py '!{cliques}' '!{domains}' !{opts} |
            gzip -9 > '!{outprefix}_cliques.tsv.gz'

        cp '!{domains}' '!{outprefix}_domains.bed.gz'
        '''
}

process PLOT_MAXIMAL_CLIQUE_SIZE_DISTRIBUTION_BY_TAD {
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
        plot_tad_maximal_clique_size_distribution.py \\
            *_cliques.tsv.gz \\
            -o '!{interaction_type}_tad_max_clique_size_distribution_abs'

        plot_tad_maximal_clique_size_distribution.py \\
            *_cliques.tsv.gz \\
            --stat='density' \\
            -o '!{interaction_type}_tad_max_clique_size_distribution_rel'
        '''
}

process PLOT_CLIQUE_SIZE_DISTRIBUTION {
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
            *_cliques.tsv.gz \\
            -o '!{interaction_type}_clique_size_distribution_abs'

        plot_maximal_clique_sizes.py \\
            *_cliques.tsv.gz \\
            --stat='density' \\
            -o '!{interaction_type}_clique_size_distribution_rel'
        '''
}
