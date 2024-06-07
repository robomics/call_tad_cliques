#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT


import argparse
import pathlib

import h5py
import matplotlib.pyplot as plt
import numpy as np


def make_cli():
    cli = argparse.ArgumentParser("Plot the expected value profiles computed by NCHG expected.")

    cli.add_argument(
        "path",
        type=pathlib.Path,
        help="Path to a .h5 file with the expected value profiles computed by NCHG expected.",
    )

    cli.add_argument("output-path", type=pathlib.Path, help="Path where to store the output plot.")

    grp = cli.add_mutually_exclusive_group()
    grp.add_argument(
        "--plot-cis",
        default=False,
        action="store_true",
        help="Plot profile(s) for cis matrices.",
    )

    grp.add_argument(
        "--plot-trans",
        default=False,
        action="store_true",
        help="Plot profile(s) for trans matrices.",
    )

    cli.add_argument("--xscale-log", action="store_true")
    cli.add_argument("--yscale-log", action="store_true")

    return cli


def plot_cis(path: pathlib.Path, x_scale_log: bool, y_scale_log: bool) -> plt.Figure:
    with h5py.File(path) as f:
        profile = f["profile/values"][:]

        file_name = f.attrs["source-file"][:]
        resolution = f.attrs["resolution"]

    fig, ax = plt.subplots(1, 1)
    ax.plot([x * resolution for x in range(len(profile))], profile)

    ax.set(
        title=f"{file_name} ({resolution} bp): expected profile (cis)",
        xlabel="Genomic distance (bp)",
        ylabel="Interaction frequency",
    )

    if x_scale_log:
        ax.set(xscale="log")
    if y_scale_log:
        ax.set(yscale="log")

    return fig


def plot_trans(path: pathlib.Path) -> plt.Figure:
    with h5py.File(path) as f:
        chrom1 = [s.decode("utf-8") for s in f["avg-values/chrom1"]]
        chrom2 = [s.decode("utf-8") for s in f["avg-values/chrom2"]]
        values = f["avg-values/value"][:]

        file_name = f.attrs["source-file"][:]
        resolution = f.attrs["resolution"]

    chroms = {c: i for i, c in enumerate(dict.fromkeys(chrom1))}
    m = np.full([len(chroms), len(chroms)], np.nan)

    for c1, c2, x in zip(chrom1, chrom2, values):
        i1 = chroms.get(c1)
        i2 = chroms.get(c2)

        if i1 is None or i2 is None:
            continue

        m[i1, i2] = x
        m[i2, i1] = x

    fig, ax = plt.subplots(1, 1)

    lb = np.nanmin(m[m != 0])
    img = ax.imshow(m, vmin=lb)
    fig.colorbar(img, ax=ax)

    ax.set(
        title=f"{file_name} ({resolution} bp): expected values (trans)",
        xlabel="Chrom2",
        ylabel="Chrom1",
    )

    ax.set_xticks(list(range(len(chroms))), list(chroms.keys()), rotation=40)
    ax.set_yticks(list(range(len(chroms))), list(chroms.keys()))

    return fig


def main():
    args = vars(make_cli().parse_args())

    if args["plot_cis"]:
        fig = plot_cis(args["path"], args["xscale_log"], args["yscale_log"])
    else:
        fig = plot_trans(args["path"])

    plt.tight_layout()
    fig.savefig(args["output-path"], dpi=300)


if __name__ == "__main__":
    main()
