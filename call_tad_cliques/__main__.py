#!/usr/bin/env python

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from __future__ import absolute_import

import argparse
import concurrent.futures
import logging
import pathlib
import shutil
import sys

import bioframe as bf
import cooler
import multiprocess as mp
import pandas as pd

from . import io, cis, trans, preprocessing, call_cliques


def find_executable(name):
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(f"Unable to find command \"{name}\"")
    return path


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument("low-resolution-cooler",
                     type=str,
                     help="Low resolution (e.g. 1Mbp) contact matrix in Cooler format.")
    cli.add_argument("high-resolution-cooler",
                     type=str,
                     help="High resolution (e.g. 50kbp) contact matrix in Cooler format.")
    cli.add_argument("tads",
                     type=str,
                     help="Path to a BED file with the list of TADs to be processed.")
    cli.add_argument("--masked-regions",
                     nargs="+",
                     type=str,
                     help="Path to one or more BED files with the list of regions to be masked out.")
    cli.add_argument("--nchg-bin",
                     type=str,
                     default="./NCHG",  # TODO changeme
                     help="Path to NCHG executable.")
    return cli


def main():
    assert len(sys.argv) == 2
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s")
    logging.getLogger().setLevel(logging.INFO)
    nchg_bin = find_executable("./NCHG")

    data_dir = sys.argv[1]
    path_to_cooler = pathlib.Path(f"{data_dir}/HiC_001_10A.mcool")
    path_to_tads = pathlib.Path(f"{data_dir}/HiC_001_10A_domains.bed")
    path_to_centromeres = pathlib.Path(f"{data_dir}/cytoBand.txt.gz")
    path_to_unmappable_regions = pathlib.Path(f"{data_dir}/unmappable_regions.json")

    chrom_sizes = io.import_chrom_sizes(f"{path_to_cooler}::/resolutions/50000")

    logging.info("Importing centromeric and masked regions...")
    masked_regions = bf.merge(pd.concat([io.import_centromeric_regions(path_to_centromeres),
                                         io.import_assembly_gaps(path_to_unmappable_regions, format="json")]))
    logging.info(f"Imported {len(masked_regions)} regions")

    tads = io.import_tads(path_to_tads)
    segments = preprocessing.fill_gaps_between_tads(tads, chrom_sizes)
    segments = preprocessing.mask_cis_regions(segments, masked_regions)

    with mp.Pool() as mp_pool, concurrent.futures.ThreadPoolExecutor() as tpool:  # noqa
        cis_df = tpool.submit(cis.process_cis_interactions, path_to_cooler, 50_000, segments, 2.0, 0.01, nchg_bin, mp_pool)
        trans_df = tpool.submit(trans.process_trans_interactions_dbg, path_to_cooler, 1_000_000, segments, masked_regions, 2.0, 0.01,
                                nchg_bin, mp_pool)

        cis_df = cis_df.result()
        trans_df = trans_df.result()

    cis_df = cis_df[cis_df["significant"]]

    cis_df.to_csv("/tmp/cis_df.tsv", sep="\t", index=False)
    trans_df.to_csv("/tmp/trans_df.tsv", sep="\t", index=False)

    columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
    interactions = pd.concat([cis_df[columns], trans_df[columns]]).sort_values(by=columns).reset_index()
    interactions.to_csv("/tmp/interactions.tsv", sep="\t", index=False)

    for chrom in cooler.Cooler(f"{path_to_cooler}::/resolutions/50000").chromnames:
        clique_stats, clique_interactions, tad_interactions, clique_sizes = call_cliques.call_cliques(segments, interactions, chrom, 5)
        clique_stats.to_csv("/tmp/clique_stats.tsv", sep="\t", index=False)
        clique_interactions.to_csv("/tmp/clique_interactions.tsv", sep="\t", index=False)
        tad_interactions.to_csv("/tmp/tad_interactions.tsv", sep="\t", index=False)
        clique_sizes.to_csv("/tmp/clique_sizes.tsv", sep="\t", index=False)
        break

    # cis_df.to_csv("/tmp/cis.tsv", sep="\t", header=False, index=False)
    # trans_df.to_csv("/tmp/trans.tsv", sep="\t", header=False, index=False)

    return 0


if __name__ == "__main__":
    main()
