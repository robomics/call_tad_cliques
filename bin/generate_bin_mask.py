#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import sys

import bioframe as bf
import pandas as pd


def make_cli() -> argparse.ArgumentParser:
    def existing_file(arg):
        if (file := pathlib.Path(arg)).exists():
            return file

        raise FileNotFoundError(arg)

    cli = argparse.ArgumentParser(description="Generate a unified list of bins to be masked")

    cli.add_argument(
        "masks",
        nargs="+",
        type=existing_file,
        help="Path to one or more BED files with the regions to be masked.",
    )
    cli.add_argument(
        "--cytoband",
        type=existing_file,
        help="Path to a cytoband file.\n" "Regions marked ad 'acen' will be used for masking.",
    )

    return cli


if __name__ == "__main__":
    args = vars(make_cli().parse_args())

    dfs = [pd.read_table(p, usecols=[0, 1, 2], names=["chrom", "start", "end"]) for p in args["masks"]]

    if args["cytoband"] is not None:
        df = pd.read_table(args["cytoband"], usecols=[0, 1, 2, 3, 4], names=["chrom", "start", "end", "arm", "type"])
        dfs.append(df.loc[df["type"] == "acen", ["chrom", "start", "end"]])

    df = bf.merge(pd.concat(dfs))[["chrom", "start", "end"]]
    df.to_csv(sys.stdout, sep="\t", header=False, index=False)
