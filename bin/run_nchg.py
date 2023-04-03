#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import argparse
import logging
import pathlib
import shutil
import subprocess as sp
import sys
import tempfile
import warnings
from typing import Union

import numpy as np
import pandas as pd
from statsmodels.stats import multitest


def make_cli() -> argparse.ArgumentParser:
    def file(arg):
        if (file_ := pathlib.Path(arg)).exists():
            return file_

        raise FileNotFoundError(arg)

    def probability(arg):
        if 0 <= (n := float(arg)) <= 1:
            return n

        raise ValueError("Not a valid probability")

    def nonnegative_float(arg):
        if (n := float(arg)) >= 0:
            return n

        raise ValueError("Not a non-negative float")

    def executable_file(name):
        if (path := shutil.which(name)) is not None:
            return pathlib.Path(path)

        raise RuntimeError(f'Unable to find command "{name}"')

    cli = argparse.ArgumentParser(description="Run NCHG and identify significant interactions")

    cli.add_argument(
        "bedpe",
        type=file,
        help="Path to a BEDPE with the interactions to be processed.",
    )
    cli.add_argument("method", choices={"intra", "inter"}, help="Processing method.")
    cli.add_argument(
        "--resolution",
        type=int,
        required="intra" in sys.argv,
        help="Resolution in bp to be passed to NCHG. Required when method=intra.",
    )
    cli.add_argument(
        "--fdr",
        type=probability,
        default=0.01,
        help="FDR threshold used to identify significant interactions.",
    )
    cli.add_argument(
        "--log-ratio",
        type=nonnegative_float,
        default=2.0,
        help="Log-ratio cutoff used to identify significant interactions.",
    )
    cli.add_argument(
        "--nchg-bin",
        type=executable_file,
        default="NCHG",
        help="Path to NCHG executable.",
    )
    cli.add_argument(
        "--drop-not-significant",
        action="store_true",
        default=False,
        help="Only return significant interactions.",
    )
    cli.add_argument("--write-header", action="store_true", default=False, help="Write BEDPE header.")

    return cli


def ensure_file_is_plaintext(path: Union[str, pathlib.Path]) -> bool:
    try:
        with open(path, "r", encoding="utf-8") as f:
            for i, _ in enumerate(f):
                if i >= 128:
                    break
    except UnicodeDecodeError:
        return False

    return True


def run_nchg(bedpe: pathlib.Path, nchg_bin: pathlib.Path, method: str, resolution: int) -> pd.DataFrame:
    assert method in {"intra", "inter"}
    if method == "intra":
        cmd = [str(nchg_bin), "--mindelta", str(resolution), "--printcounts"]
    else:
        cmd = [str(nchg_bin), "--useInter", "--printcounts"]

    # NCHG is picky about file types and returns incomplete outputs if e.g. the input file is a FIFO
    # To be on the safe side, if bedpe is not a regular file, we make a temporary copy
    with tempfile.NamedTemporaryFile(mode="w+b") as tmpbedpe:
        if bedpe.is_file():
            cmd.append(str(bedpe))
        else:
            with open(bedpe, "rb") as f:
                shutil.copyfileobj(f, tmpbedpe.file)
            cmd.append(tmpbedpe.name)

        # NCHG returns an empty output when the input file is not in plaintext
        if not ensure_file_is_plaintext(cmd[-1]):
            raise RuntimeError(f"File {bedpe} is not a valid UTF-8 text file")

        logging.debug("Spawning subprocess for %s...", cmd)
        with sp.Popen(cmd, stdin=sp.DEVNULL, stdout=sp.PIPE, stderr=None, encoding="utf-8") as proc:
            columns = [
                "chrom1",
                "start1",
                "end1",
                "chrom2",
                "start2",
                "end2",
                "pvalue",
                "observed_count",
                "expected_count",
            ]

            logging.debug("Capturing stdout for %s...", cmd)
            # Asyncronously read NCHG output into a dataframe
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")  # Ignore warnings about data loss header and data shape mismatch
                df = pd.read_table(proc.stdout, names=columns, index_col=False, delim_whitespace=True)
            logging.debug(
                "DONE capturing stdout for %s.\nSuccessfully parsed %d records",
                cmd,
                len(df),
            )

            proc.wait()
            if proc.returncode != 0:
                raise RuntimeError(f"NCHG terminated with exit code {proc.returncode}")

    logging.debug("Command %s terminated without errors", cmd)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df["log_ratio"] = np.nan_to_num(np.log2(df["observed_count"]) - np.log2(df["expected_count"]))

    return df.sort_values(by=["chrom1", "start1", "chrom2", "start2"])


def correct_pval(df: pd.DataFrame, correction_method: str, log_ratio_thresh: float, fdr_thresh: float) -> pd.DataFrame:
    assert len(df) != 0
    columns = list(df.columns.copy())
    # Correct pvalues
    df["pvalue_corrected"] = multitest.multipletests(df["pvalue"], method=correction_method)[1]

    # Mark significant interactions
    if log_ratio_thresh > 0 or fdr_thresh < 1:
        df = df[(df["log_ratio"] >= log_ratio_thresh) & (df["pvalue_corrected"] <= fdr_thresh)]

    # Reorder columns and return df
    columns.insert(columns.index("pvalue") + 1, "pvalue_corrected")

    return df[columns]


def setup_logger(level=logging.INFO):
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s")
    logging.getLogger().setLevel(level)


def main():
    args = vars(make_cli().parse_args())

    if not args["drop_not_significant"]:
        args["log_ratio"] = 0.0
        args["fdr_thresh"] = 1.0

    df = run_nchg(
        args["bedpe"],
        nchg_bin=args["nchg_bin"],
        method=args["method"],
        resolution=args["resolution"],
    )
    df = correct_pval(
        df,
        correction_method="fdr_bh",
        log_ratio_thresh=args["log_ratio"],
        fdr_thresh=args["fdr"],
    )

    df = df.sort_values(by=["chrom1", "start1", "chrom2", "start2"])

    if args["write_header"]:
        df.columns = [f"#{df.columns[0]}"] + df.columns[1:].tolist()

    df.to_csv(sys.stdout, sep="\t", index=False, header=args["write_header"])


if __name__ == "__main__":
    setup_logger()
    main()
