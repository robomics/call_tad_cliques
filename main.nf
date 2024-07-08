#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

params.publish_dir = params.outdir

include { SAMPLESHEET } from './subworkflows/samplesheet.nf'
include { TADS } from './subworkflows/tads.nf'
include { NCHG } from './subworkflows/nchg.nf'
include { CLIQUES } from './subworkflows/cliques.nf'


// Workaround for optional input files: https://github.com/nextflow-io/nextflow/issues/1694
def make_optional_input(path) {
    if (path?.trim()) {
        return [file(path)]
    }
    return []
}


workflow {

    log.info("-- PARAMETERS")
    log.info("")
    if (params.sample_sheet) {
        log.info("-- sample_sheet: ${params.sample_sheet}")
    } else {
        log.info("-- sample: ${params.sample}")
        log.info("-- hic_file: ${params.hic_file}")
        log.info("-- resolution: ${params.resolution}")
        log.info("-- tads: ${params.tads}")
        log.info("-- mask: ${params.mask}")
    }
    log.info("-- outdir: ${params.outdir}")
    log.info("-- publish_dir_mode: ${params.publish_dir_mode}")
    log.info("-- cytoband: ${params.cytoband}")
    log.info("-- assembly_gaps: ${params.assembly_gaps}")

    log.info("-- hicexplorer_hic_norm: ${params.hicexplorer_hic_norm}")
    log.info("-- hicexplorer_cool_norm: ${params.hicexplorer_cool_norm}")

    log.info("-- nchg_mad_max: ${params.nchg_mad_max}")
    log.info("-- nchg_bad_bin_fraction: ${params.nchg_bad_bin_fraction}")

    log.info("-- nchg_fdr_cis: ${params.nchg_fdr_cis}")
    log.info("-- nchg_log_ratio_cis: ${params.nchg_log_ratio_cis}")
    log.info("-- nchg_fdr_trans: ${params.nchg_fdr_trans}")
    log.info("-- nchg_log_ratio_trans: ${params.nchg_log_ratio_trans}")

    log.info("-- clique_size_thresh: ${params.clique_size_thresh}")
    log.info("-- call_cis_cliques: ${params.call_cis_cliques}")
    log.info("-- call_trans_cliques: ${params.call_trans_cliques}")

    log.info("-- plot_format: ${params.plot_format}")
    log.info("-- hic_tgt_resolution_plots: ${params.hic_tgt_resolution_plots}")
    log.info("-- plot_sig_interactions_cmap_lb: ${params.plot_sig_interactions_cmap_lb}")
    log.info("-- plot_sig_interactions_cmap_ub: ${params.plot_sig_interactions_cmap_ub}")
    log.info("-- skip_expected_plots: ${params.skip_expected_plots}")
    log.info("-- skip_sign_interaction_plots: ${params.skip_sign_interaction_plots}")

    log.info("-- zstd_compression_lvl: ${params.zstd_compression_lvl}")
    log.info("")

    SAMPLESHEET(
        params.sample_sheet,
        params.sample,
        params.hic_file,
        params.resolution,
        params.tads,
        params.mask
    )

    SAMPLESHEET.out.tsv.set { sample_sheet }

    TADS(
        sample_sheet,
        params.hicexplorer_hic_norm,
        params.hicexplorer_cool_norm
    )

    NCHG(
        sample_sheet,
        TADS.out.tsv,
        params.nchg_mad_max,
        params.nchg_bad_bin_fraction,
        params.cytoband,
        params.assembly_gaps
    )


    sample_sheet
        .splitCsv(sep: "\t", header: true)
        .map { row -> tuple(row.sample, make_optional_input(row.mask))
        }
        .set { masks }

    CLIQUES(
        NCHG.out.tsv,
        TADS.out.tsv,
        masks,
        params.clique_size_thresh
    )
}
