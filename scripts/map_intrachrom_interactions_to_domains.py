#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import argparse
import itertools
import logging
import pathlib
import warnings

import bioframe as bf
import cooler
import numpy as np
import pandas as pd


def make_cli() -> argparse.ArgumentParser:
    def existing_file(arg):
        if (file := pathlib.Path(arg)).exists():
            return file

        raise FileNotFoundError(arg)

    def existing_cooler(arg):
        return cooler.Cooler(arg).uri

    cli = argparse.ArgumentParser(
        description="Map intrachromosomal (cis) interactions between a list of domains in BED format"
    )

    cli.add_argument(
        "cooler",
        type=existing_cooler,
        help="Path to a Cooler file (URI syntax supported).",
    )
    cli.add_argument("domains", type=existing_file, help="Path to a BED file with a list of TADs.")
    return cli


def import_chrom_sizes(cf: cooler.Cooler) -> pd.DataFrame:
    return bf.from_dict(cf.chromsizes)


def import_domains(path_to_domains: pathlib.Path, schema: str = "bed3") -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # Ignore warnings about data loss when parsing BED4+
        return bf.read_table(str(path_to_domains), schema=schema)


def map_intrachrom_interactions_to_domains_serial_slow(
    cf: cooler.Cooler, domains: pd.DataFrame
) -> pd.DataFrame:  # noqa
    logging.info("Mapping pairwise cis-interactions between %d regions...", len(domains))

    def interval_to_str(interval):  # noqa
        chrom, start, end = interval
        return f"{chrom}:{start}-{end}"

    num_regions = 0
    nnz_regions = 0
    sel = cf.matrix(balance=False, as_pixels=True, join=True)
    for chrom, df in domains.groupby("chrom"):
        queries = tuple(df.itertuples(index=False))
        for q1, q2 in itertools.product(queries, repeat=2):
            if q1 == q2:
                continue

            tot = sel.fetch(interval_to_str(q1), interval_to_str(q2))["count"].sum()
            num_regions += 1
            if tot != 0:
                print(f"{chrom}\t{q1.start}\t{q1.end}\t{chrom}\t{q2.start}\t{q2.end}\t{tot}")
                nnz_regions += 1

    logging.info(
        "Mapping produced %d/%d regions with one or more interactions",
        nnz_regions,
        num_regions,
    )


def map_intrachrom_interactions_to_domains_serial(cf: cooler.Cooler, domains: pd.DataFrame) -> pd.DataFrame:  # noqa
    logging.info("Mapping pairwise cis-interactions between %d regions...", len(domains))

    nnz_regions = 0
    sel = cf.matrix(balance=False, as_pixels=False)
    bin_size = cf.binsize
    for chrom, df in domains.groupby("chrom"):
        logging.info("Processing %s...", chrom)
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

    logging.info("Found %d regions with one or more interactions", nnz_regions)


def setup_logger(level=logging.INFO):
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s")
    logging.getLogger().setLevel(level)


def main():
    args = vars(make_cli().parse_args())

    cooler_file = cooler.Cooler(args["cooler"])

    if (bs := cooler_file.binsize) < 20_000:
        logging.warning(
            "Mapping interactions using a bin size of %d bp. "
            "Processing of large chromosomes may require a lot of memory.",
            bs,
        )

    chrom_sizes = import_chrom_sizes(cooler_file)
    if len(chrom_sizes) == 0:
        raise RuntimeError(f"Unable to import any chromosome from Cooler at URI {cooler_file.uri}")

    domains = import_domains(args["domains"])
    if len(domains) == 0:
        raise RuntimeError("Unable to import any domain from BED file " + str(args["domains"]))

    map_intrachrom_interactions_to_domains_serial(cooler_file, domains)


if __name__ == "__main__":
    setup_logger()
    main()
