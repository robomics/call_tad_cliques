#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import argparse

import pandas as pd
import pathlib
import warnings
import bioframe as bf
import logging
import itertools
import numpy as np


def make_cli() -> argparse.ArgumentParser:
    def existing_file(arg):
        if (file := pathlib.Path(arg)).exists():
            return file

        raise FileNotFoundError(arg)

    cli = argparse.ArgumentParser(
        description="Map interchromosomal (trans) interactions between a list of domains in BED format"
    )

    cli.add_argument("domains", type=existing_file, help="Path to a BED file with the domain list.")
    cli.add_argument("bedpe", type=existing_file, help="Path to a BEDPE file with the list of significant trans interactions.")

    return cli


def import_domains(path_to_domains: pathlib.Path, schema: str = "bed3") -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # Ignore warnings about data loss when parsing BED4+
        return bf.read_table(str(path_to_domains), schema=schema)


def import_interactions(path_to_interactions: pathlib.Path, schema: str = "bedpe") -> pd.DataFrame:
    return bf.read_table(str(path_to_interactions), schema=schema)


def map_intrachrom_interactions_to_domains_serial(cf: cooler.Cooler, domains: pd.DataFrame) -> pd.DataFrame:  # noqa
    logging.info(f"Mapping pairwise cis-interactions between {len(domains)} regions...")

    nnz_regions = 0
    sel = cf.matrix(balance=False, as_pixels=False)
    bin_size = cf.binsize
    for chrom, df in domains.groupby("chrom"):
        logging.info(f"Processing {chrom}...")
        pixels = np.triu(sel.fetch(chrom))
        qstarts = df["start"].to_numpy()
        qends = df["end"].to_numpy()

        idx0 = (qstarts / bin_size).astype(int)
        idx1 = (qends / bin_size).astype(int)

        num_queries = len(idx0)
        for i in range(num_queries):
            i0, i1 = idx0[i], idx1[i]
            for j in range(i + 1, num_queries):
                i2, i3 = idx0[j], idx1[j]
                if (tot := np.sum(pixels[i0:i1, i2:i3])) != 0:
                    print(f"{chrom}\t{qstarts[i]}\t{qends[i]}\t{chrom}\t{qstarts[j]}\t{qends[j]}\t{tot}")
                    nnz_regions += 1

    logging.info(f"Found {nnz_regions} regions with one or more interactions")


def setup_logger(level=logging.INFO):
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s")
    logging.getLogger().setLevel(level)


if __name__ == "__main__":
    setup_logger()

    args = vars(make_cli().parse_args())

    cooler_file = cooler.Cooler(args["cooler"])

    if (bs := cooler_file.binsize) < 20000:
        logging.warning(f"Mapping interactions using a bin size of {bs} bp, processing of large chromosomes may require a lot of memory.")

    chrom_sizes = import_chrom_sizes(cooler_file)
    if len(chrom_sizes) == 0:
        raise RuntimeError(f"Unable to import any chromosome from Cooler at URI {cooler_file.uri}")

    domains = import_domains(args["domains"])
    if len(domains) == 0:
        raise RuntimeError("Unable to import any domain from BED file " + str(args["domains"]))

    map_intrachrom_interactions_to_domains_serial(cooler_file, domains)