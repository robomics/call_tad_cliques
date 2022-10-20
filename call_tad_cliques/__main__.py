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

import bioframe as bf
import cooler
import multiprocess as mp
import pandas as pd
import numpy as np
from . import io, cis, trans, preprocessing, call_cliques


def find_executable(name):
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(f'Unable to find command "{name}"')
    return path


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "low-resolution-cooler", type=str, help="Low resolution (e.g. 1Mbp) contact matrix in Cooler format."
    )
    cli.add_argument(
        "high-resolution-cooler", type=str, help="High resolution (e.g. 50kbp) contact matrix in Cooler format."
    )
    cli.add_argument(
        "domains", type=str, help="Path to a BED file with the list of domains (usually TADs) to be processed."
    )
    cli.add_argument("--cytoband-file", type=str, help="Path to a cytoBand file. Used to mask centromeric regions.")
    cli.add_argument(
        "--unmappable-regions",
        type=str,
        help="Path to a JSON file with a list of assembly gaps/unmappable regions.\n"
        "Example: https://api.genome.ucsc.edu/getData/track?genome=hg38;track=gap;jsonOutputArrays=1",
    )
    cli.add_argument("--nchg-bin", type=str, default="./NCHG", help="Path to NCHG executable.")  # TODO changeme
    cli.add_argument("--odd-log-ratio-cis", type=float, default=2.0)
    cli.add_argument("--odd-log-ratio-trans", type=float, default=2.0)
    cli.add_argument("--fdr-cis", type=float, default=0.01)
    cli.add_argument("--fdr-trans", type=float, default=0.01)
    cli.add_argument("--outprefix", type=str, required=True, help="TODO")
    return cli


def setup_logger(level=logging.INFO):
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s")
    logging.getLogger().setLevel(level)


def validate_coolers(path_to_lowres_cooler, path_to_highres_cooler):
    c1 = cooler.Cooler(path_to_lowres_cooler)
    c2 = cooler.Cooler(path_to_highres_cooler)

    if c1.binsize < c2.binsize:
        raise RuntimeError("TODO: low res < high res")

    if not c1.chromsizes.equals(c2.chromsizes):
        raise RuntimeError("TODO: different assembly")

    return c1, c2


def annotate_masked_regions(path_to_cytoband, path_to_assembly_gaps):
    masked_regions = []
    if path_to_cytoband is not None:
        logging.info(f"Importing centromeric regions from {path_to_cytoband}...")
        masked_regions.append(io.import_centromeric_regions(path_to_cytoband))
        logging.info(f"DONE! Imported {len(masked_regions[-1])} regions!")

    if path_to_assembly_gaps is not None:
        logging.info(f"Importing assembly gaps from {path_to_assembly_gaps}...")

        format = "json" if path_to_assembly_gaps.endswith(".json") else "bed"
        masked_regions.append(io.import_assembly_gaps(path_to_assembly_gaps, format=format))

        logging.info(f"DONE! Imported {len(masked_regions[-1])} regions!")

    if len(masked_regions) == 0:
        return pd.DataFrame({"chrom": [], "start": [], "end": []})

    logging.info(f"Generating list of regions to be masked...")
    masked_regions = bf.merge(pd.concat(masked_regions))

    masked_region_size_mbp = np.sum(np.abs(masked_regions["end"] - masked_regions["start"])) / 1.0e6
    logging.info(
        f"DONE! Generated {len(masked_regions)} intervals for masking ({masked_region_size_mbp:.2f} Mbp in total)"
    )

    return masked_regions


def import_domains(chrom_sizes: pd.DataFrame, path_to_tads, masked_regions: pd.DataFrame):
    tads = io.import_tads(path_to_tads)
    domains = preprocessing.fill_gaps_between_tads(tads, chrom_sizes)
    domains = preprocessing.mask_cis_regions(domains, masked_regions)

    return domains


def import_interactions(path_to_low_res_cooler, path_to_high_res_cooler, domains, masked_regions, nchg_binary):
    with mp.Pool() as mp_pool, concurrent.futures.ThreadPoolExecutor() as tpool:  # noqa
        cis_df = tpool.submit(
            cis.process_cis_interactions,
            path_to_high_res_cooler,
            domains,
            odd_log_ratio_thresh=2.0,
            fdr_thresh=0.01,
            nchg_bin=nchg_binary,
            exec_pool=mp_pool,
        )

        trans_df = tpool.submit(
            trans.process_trans_interactions,
            path_to_low_res_cooler,
            domains,
            masked_regions,
            odd_log_ratio_thresh=2.0,
            fdr_thresh=0.01,
            nchg_bin=nchg_binary,
            exec_pool=mp_pool,
        )

        return cis_df.result(), trans_df.result()


def main():
    setup_logger()

    args = vars(make_cli().parse_args())

    nchg_bin = find_executable(args["nchg_bin"])
    low_res_cooler, high_res_cooler = validate_coolers(args["low-resolution-cooler"], args["high-resolution-cooler"])

    chrom_sizes = io.import_chrom_sizes(low_res_cooler)

    masked_regions = annotate_masked_regions(args.get("cytoband_file"), args.get("unmappable_regions"))
    domains = import_domains(chrom_sizes, args["domains"], masked_regions)

    assert len(domains) > 0

    # path_to_cooler = pathlib.Path(f"{data_dir}/HiC_001_10A.mcool")
    # path_to_tads = pathlib.Path(f"{data_dir}/HiC_001_10A_domains.bed")
    # path_to_centromeres = pathlib.Path(f"{data_dir}/cytoBand.txt.gz")
    # path_to_unmappable_regions = pathlib.Path(f"{data_dir}/unmappable_regions.json")

    cis_df, trans_df = import_interactions(low_res_cooler, high_res_cooler, domains, masked_regions, nchg_bin)

    cis_df = cis_df[cis_df["significant"]]

    cis_df.to_csv("/tmp/cis_df.tsv", sep="\t", index=False)
    trans_df.to_csv("/tmp/trans_df.tsv", sep="\t", index=False)

    columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
    interactions = pd.concat([cis_df[columns], trans_df[columns]]).sort_values(by=columns).reset_index()
    interactions.to_csv("/tmp/interactions.tsv", sep="\t", index=False)

    outprefix = pathlib.Path(args["outprefix"])
    outprefix.parent.mkdir(exist_ok=True)
    for chrom in chrom_sizes["chrom"]:
        clique_stats, clique_interactions, tad_interactions, clique_sizes = call_cliques.call_cliques(
            domains, interactions, chrom, 5
        )
        clique_stats.to_csv(f"{outprefix}clique_stats.tsv", sep="\t", index=False)
        clique_interactions.to_csv(f"{outprefix}clique_interactions.tsv", sep="\t", index=False)
        tad_interactions.to_csv(f"{outprefix}tad_interactions.tsv", sep="\t", index=False)
        clique_sizes.to_csv(f"{outprefix}clique_sizes.tsv", sep="\t", index=False)

    # outprefix = pathlib.Path(args["outprefix"])
    # outprefix.parent.mkdir(exist_ok=True)
    # cis_df.to_csv(f"{outprefix}cis.tsv", sep="\t", header=False, index=False)
    # trans_df.to_csv(f"{outprefix}trans.tsv", sep="\t", header=False, index=False)

    return 0


if __name__ == "__main__":
    main()
