# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import bioframe as bf
import numpy as np
import pandas as pd
import pathlib
from typing import Union
import cooler
import logging
import tempfile
import itertools

from .run_nchg import run_nchg_and_generate_df, correct_nchg_pvalue_and_mark_significant_interactions


def __map_interchrom_interactions_to_domains(trans_interactions: pd.DataFrame, domains: pd.DataFrame) -> pd.DataFrame:
    df1 = trans_interactions.copy()
    df1["start1"] = np.array((df1["start1"] + df1["end1"]) / 2, dtype=int)
    df1["end1"] = df1["start1"] + 1

    df1["start2"] = np.array((df1["start2"] + df1["end2"]) / 2, dtype=int)
    df1["end2"] = df1["start2"] + 1

    df2 = bf.overlap(domains, df1, how="left", cols2=["chrom1", "start1", "end1"], suffixes=["1", "__"])  # noqa
    df3 = bf.overlap(domains, df1, how="left", cols2=["chrom2", "start2", "end2"], suffixes=["2", "__"])  # noqa

    df4 = df2[["chrom1", "start1", "end1"]].copy()
    df4[["chrom2", "start2", "end2"]] = df3[["chrom2", "start2", "end2"]]

    return (
        df4.sort_values(by=["chrom1", "start1", "chrom2", "start2"]).dropna().drop_duplicates().reset_index(drop=True)
    )


def __identify_significant_trans_interactions(
    cooler_uri: str,
    masked_regions: pd.DataFrame,
    log_ratio_thresh: float,
    fdr_thresh: float,
    multiple_test_correction_method: str,
    nchg_bin: Union[pathlib.Path, str],
) -> pd.DataFrame:
    cooler_file = cooler.Cooler(cooler_uri)
    pixel_selector = cooler_file.matrix(balance=False, as_pixels=True, join=True)

    # We are using a temp file because NCHG returns incomplete results when the input BEDPE
    # is not a regular file
    logging.info("Identifying significant trans-regions...")
    with tempfile.NamedTemporaryFile(mode="w+t") as tmpbedpe:
        logging.debug(f"Writing nnz pixels for trans-regions to file {tmpbedpe.name}...")
        num_nnz_pixels = 0
        for chrom1, chrom2 in itertools.product(cooler_file.chromnames, repeat=2):
            if chrom1 == chrom2:
                continue

            pixels = pixel_selector.fetch(chrom1, chrom2)
            pixels = bf.setdiff(pixels, masked_regions, cols1=["chrom1", "start1", "end1"])  # noqa
            pixels = bf.setdiff(pixels, masked_regions, cols1=["chrom2", "start2", "end2"])  # noqa
            num_nnz_pixels += len(pixels)

            pixels.to_csv(tmpbedpe.name, mode="a", sep="\t", header=False, index=False)

        tmpbedpe.flush()
        logging.debug(f"Written {num_nnz_pixels} nnz for trans-regions to file {tmpbedpe.name}")

        df = run_nchg_and_generate_df(tmpbedpe.name, nchg_bin, method="inter", resolution=cooler_file.binsize)

    assert len(df) > 0
    df = correct_nchg_pvalue_and_mark_significant_interactions(
        df, multiple_test_correction_method, log_ratio_thresh, fdr_thresh
    )

    num_significant_interactions = df["significant"].sum()
    logging.info(f"DONE! Identified {num_significant_interactions}/{len(df)} significant trans-regions")
    return df


def process_trans_interactions(
    cf: cooler.Cooler,
    domains: pd.DataFrame,
    masked_regions: Union[pd.DataFrame, None],
    odd_log_ratio_thresh: float,
    fdr_thresh: float,
    nchg_bin: Union[str, pathlib.Path],
    exec_pool,
) -> pd.DataFrame:
    trans_df = exec_pool.apply(
        __identify_significant_trans_interactions,
        args=[cf.uri, masked_regions, odd_log_ratio_thresh, fdr_thresh, "fdr_bh", nchg_bin],
    )
    trans_df.to_csv("/tmp/trans_001.tsv", sep="\t", header=False, index=False)

    trans_df = exec_pool.apply(
        __map_interchrom_interactions_to_domains, args=[trans_df[trans_df["significant"]], domains]
    )
    trans_df.to_csv("/tmp/trans_002.tsv", sep="\t", header=False, index=False)
    return trans_df


def process_trans_interactions_dbg(
    cf: cooler.Cooler,
    domains: pd.DataFrame,
    masked_regions: Union[pd.DataFrame, None],
    odd_log_ratio_thresh,
    fdr_thresh,
    nchg_bin,
    pool,
):
    nchg_out = pd.read_table("/tmp/nchg_output.tsv", index_col=None)
    trans_df = correct_nchg_pvalue_and_mark_significant_interactions(
        nchg_out, "fdr_bh", odd_log_ratio_thresh, fdr_thresh
    )

    trans_df = __map_interchrom_interactions_to_domains(trans_df[trans_df["significant"]], domains)
    trans_df.to_csv("/tmp/trans_002.tsv", sep="\t", header=False, index=False)
    return trans_df
