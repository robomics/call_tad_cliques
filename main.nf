#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2


def row_to_tuple(row) {
    tuple(row.sample,
          file(row.cooler_cis, checkIfExists: true),
          file(row.cooler_trans, checkIfExists: true),
          row.tads,
          row.cis_resolution,
          row.trans_resolution)
}

params.publish_dir = params.outdir

include { SAMPLESHEET } from './subworkflows/samplesheet.nf'
include { TADS } from './subworkflows/tads.nf'
include { NCHG } from './subworkflows/nchg.nf'
include { CLIQUES } from './subworkflows/cliques.nf'


workflow {

    SAMPLESHEET(
        params.sample_sheet,
        params.sample,
        params.hic_file,
        params.resolution,
        params.tads
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
        params.assembly_gaps,
        params.custom_mask,
    )

    CLIQUES(
        NCHG.out.tsv,
        TADS.out.tsv,
        params.clique_size_thresh
    )
}
