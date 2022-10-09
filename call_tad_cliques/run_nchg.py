# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import logging
import pathlib
import subprocess as sp
from typing import Union

import pandas as pd
from statsmodels.stats import multitest


def run_nchg_and_generate_df(bedpe: Union[pathlib.Path, str],
                             nchg_bin: Union[pathlib.Path, str],
                             method: str,
                             resolution: int) -> pd.DataFrame:
    assert method in {"intra", "inter"}
    assert pathlib.Path(bedpe).is_file()

    if method == "intra":
        cmd = [str(nchg_bin),
               "--mindelta", str(resolution),
               "--printcounts",
               str(bedpe)]
    else:
        cmd = [str(nchg_bin),
               "--useInter",
               "--printcounts",
               str(bedpe)]

    logging.debug(f"Running {cmd} as subprocess...")

    proc = sp.Popen(cmd,
                    stdin=sp.DEVNULL,
                    stdout=sp.PIPE,
                    stderr=None,
                    encoding="utf-8")

    columns = ["chrom1", "start1", "end1",
               "chrom2", "start2", "end2",
               "pvalue",
               "observed_count", "expected_count",
               "odds_ratio", "omega"]

    logging.debug(f"Capturing stdout for {cmd}...")
    # Asyncronously read NCHG output into a dataframe
    df = pd.read_table(proc.stdout,
                       names=columns,
                       index_col=False,
                       delim_whitespace=True)
    logging.debug(f"DONE capturing stdout for {cmd}."
                  f"Successfully parsed {len(df)} records")

    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"NCHG terminated with exit code {proc.returncode}")

    logging.debug(f"Command {cmd} terminated without errors")

    return df


def correct_nchg_pvalue_and_mark_significant_interactions(df: pd.DataFrame,
                                                          correction_method: str,
                                                          log_ratio_thresh: float,
                                                          fdr_thresh: float) -> pd.DataFrame:
    assert len(df) != 0
    columns = list(df.columns.copy())
    # Correct pvalues
    df["pvalue_corrected"] = multitest.multipletests(df["pvalue"], method=correction_method)[1]

    # Mark significant interactions
    df["significant"] = (df["odds_ratio"] >= log_ratio_thresh) & (df["pvalue_corrected"] <= fdr_thresh)

    # Reorder columns and return df
    columns.insert(columns.index("pvalue") + 1, "pvalue_corrected")
    columns.insert(columns.index("pvalue_corrected") + 1, "significant")

    return df[columns]
