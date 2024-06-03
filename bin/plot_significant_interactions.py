#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT


import argparse
import pathlib
from typing import Tuple

import hictkpy
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
from matplotlib.colors import LogNorm


def make_cli():
    cli = argparse.ArgumentParser("Plot the significant interactions for a given chromosome pair.")

    cli.add_argument(
        "hic-file",
        type=pathlib.Path,
        help="Path to a Hi-C matrix in .hic, .mcool, or .cool format",
    )

    cli.add_argument(
        "parquet",
        type=pathlib.Path,
        help="Path to the .parquet file produced by NCHG filter.",
    )

    cli.add_argument("chrom1", type=str, help="Name of the first chromosome.")
    cli.add_argument("chrom2", type=str, help="Name of the second chromosome.")

    cli.add_argument("output-path", type=pathlib.Path, help="Path where to store the output plot.")

    cli.add_argument(
        "--resolution",
        default=0,
        type=int,
        help="Hi-C matrix resolution.\n" "Required when input matrix is in .hic or .mcool format.",
    )

    cli.add_argument(
        "--pvalue",
        type=float,
        default=1.0,
        help="Pvalue threshold used to classify significant interactions.",
    )
    cli.add_argument(
        "--log-ratio",
        type=float,
        default=-np.inf,
        help="Log-ratio threshold used to classify significant interactions.",
    )

    cli.add_argument(
        "--min-value",
        type=float,
        default=None,
        help="Lower bound for the color scale of the log-ratio heatmap.",
    )
    cli.add_argument(
        "--max-value", type=float, default=None, help="Upper bound for the color scale of the log-ratio heatmap."
    )

    return cli


def preprocess_data(df: pd.DataFrame) -> pd.DataFrame:
    df["chrom1"] = df["chrom1"].astype("string")
    df["chrom2"] = df["chrom2"].astype("string")

    chroms = pd.Series(df["chrom1"].unique().tolist() + df["chrom2"].unique().tolist()).unique()

    df.sort_values(["chrom1", "start1", "chrom2", "start2"], inplace=True)
    df["chrom1"] = pd.Categorical(df["chrom1"], categories=chroms)
    df["chrom2"] = pd.Categorical(df["chrom2"], categories=chroms)

    df.reset_index(drop=True, inplace=True)

    return df


def import_data(path: str, chrom1: str, chrom2: str) -> pd.DataFrame:
    df = preprocess_data(pd.read_parquet(path))
    return df[(df["chrom1"] == chrom1) & (df["chrom2"] == chrom2)]


def fetch_hic_matrix(path: pathlib.Path, resolution: int, chrom1: str, chrom2: str) -> Tuple[int, npt.NDArray]:
    path = str(path)
    if hictkpy.is_hic(path) or hictkpy.is_mcool_file(path):
        if resolution is None:
            raise RuntimeError("--resolution is required when Hi-C matrix is in .hic or .mcool format.")
    f = hictkpy.File(path, resolution)
    f.fetch(chrom1, chrom2).to_df()
    return f.resolution(), f.fetch(chrom1, chrom2).to_numpy()


def df_to_matrix(df: pd.DataFrame, bin_size: int, shape: Tuple[int, int]) -> npt.NDArray:
    m = np.full(shape, np.nan)

    cols = ["start1", "end1", "start2", "end2", "log_ratio"]
    for _, (start1, end1, start2, end2, log_ratio) in df[cols].iterrows():
        i0 = int(start1 // bin_size)
        i1 = int(end1 // bin_size)
        j0 = int(start2 // bin_size)
        j1 = int(end2 // bin_size)

        if np.isfinite(log_ratio):
            m[i0:i1, j0:j1] = np.maximum(np.nan_to_num(m[i0:i1, j0:j1], nan=-np.inf), log_ratio)

    return m


def main():
    args = vars(make_cli().parse_args())

    df = import_data(args["parquet"], args["chrom1"], args["chrom2"])

    if "pvalue_corrected" in df:
        pval_col = "pvalue_corrected"
    else:
        pval_col = "pvalue"

    df = df[(df[pval_col] < args["pvalue"]) & (df["log_ratio"] >= args["log_ratio"]) & (~df["log_ratio"].isna())]

    resolution, obs_matrix = fetch_hic_matrix(args["hic-file"], args["resolution"], args["chrom1"], args["chrom2"])

    sig_matrix = df_to_matrix(df, resolution, obs_matrix.shape)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(2 * 6.4, 6.4))

    img0 = ax0.imshow(obs_matrix, norm=LogNorm(), interpolation="nearest")
    img1 = ax1.imshow(sig_matrix, vmin=args["min_value"], vmax=args["max_value"], interpolation="nearest", cmap="Reds")

    ax0.set(xlabel=args["chrom2"], ylabel=args["chrom1"], title="Observed matrix")
    ax1.set(
        xlabel=args["chrom2"],
        ylabel=args["chrom1"],
        title="Significant interactions (Log-ratio)",
    )

    plt.colorbar(img0, ax=ax0)
    plt.colorbar(img1, ax=ax1)

    plt.tight_layout()
    fig.savefig(args["output-path"], dpi=300)


if __name__ == "__main__":
    main()
