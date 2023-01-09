#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import sys
import warnings

import bioframe as bf
import numpy as np
import pandas as pd


def make_cli() -> argparse.ArgumentParser:
    def existing_file(arg):
        if (file := pathlib.Path(arg)).exists():
            return file

        raise FileNotFoundError(arg)

    cli = argparse.ArgumentParser(
        description="Fill gaps between TADs and output domains without gaps in BED format to stdout"
    )

    cli.add_argument("chrom-sizes", type=existing_file, help="Path to a chrom.sizes file.")
    cli.add_argument("tads", type=existing_file, help="Path to a BED file with a list of TADs.")

    return cli


def import_tads(path_to_tads: pathlib.Path) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # Ignore warnings about data loss when parsing BED4+
        return bf.read_table(str(path_to_tads), schema="bed3")


def import_chrom_sizes(path_to_chrom_sizes: pathlib.Path) -> pd.DataFrame:
    df = pd.read_table(str(path_to_chrom_sizes), names=["chrom", "end"], usecols=list(range(2)))
    df["start"] = 0

    return df[bf.SCHEMAS["bed3"]]


def fill_gaps_between_tads(chrom_sizes: pd.DataFrame, tads: pd.DataFrame) -> pd.DataFrame:
    """
    Generate a dataframe of contiguous domains where gaps between TAD domains are filled with non-TAD domains
    """
    df = pd.concat([tads, bf.complement(tads, chrom_sizes)])

    # Map chromosome names to their rank in the chrom_sizes df
    mappings = {row.chrom: i for i, row in enumerate(chrom_sizes.itertuples(index=False))}  # noqa

    # Assign a chromosome rank to each entry in the segment df
    df["chrom_rank"] = np.array([mappings.get(chrom) for chrom in df["chrom"]], dtype=float)
    df[~np.isfinite(df["chrom_rank"])] = len(chrom_sizes) + 1

    # Sort domains such that chromosomes appear in the same order as in the chrom_sizes dataframe
    return (
        df.sort_values(by=["chrom_rank", "start", "end"])
        .drop(columns=["chrom_rank", "view_region"])
        .reset_index(drop=True)
    )


if __name__ == "__main__":
    args = vars(make_cli().parse_args())

    domains = fill_gaps_between_tads(import_chrom_sizes(args["chrom-sizes"]), import_tads(args["tads"]))
    domains.to_csv(sys.stdout, sep="\t", header=False, index=False)
