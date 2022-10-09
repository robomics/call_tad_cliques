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

def __map_interchrom_interactions_to_segments(trans_interactions: pd.DataFrame, segments: pd.DataFrame) -> pd.DataFrame:
    df1 = trans_interactions.copy()
    df1["start1"] = np.array((df1["start1"] + df1["end1"]) / 2, dtype=int)
    df1["end1"] = df1["start1"] + 1

    df1["start2"] = np.array((df1["start2"] + df1["end2"]) / 2, dtype=int)
    df1["end2"] = df1["start2"] + 1

    df2 = bf.overlap(df1, segments, how="left", cols1=["chrom1", "start1", "end1"], suffixes=["__", "1"])  # noqa
    df3 = bf.overlap(df1, segments, how="left", cols1=["chrom2", "start2", "end2"], suffixes=["__", "2"])  # noqa

    df4 = df2[["chrom1", "start1", "end1"]].copy()
    df4[["chrom2", "start2", "end2"]] = df3[["chrom2", "start2", "end2"]]
    df4["pvalue_corrected"] = df3["pvalue_corrected__"].tolist()

    return df4.sort_values(by=["chrom1", "start1", "chrom2", "start2"]).dropna().reset_index(drop=True)


def __identify_significant_trans_interactions(path_to_cooler: Union[pathlib.Path, str],
                                            masked_regions: pd.DataFrame,
                                            log_ratio_thresh: float,
                                            fdr_thresh: float,
                                            multiple_test_correction_method: str,
                                            nchg_bin: Union[pathlib.Path, str]) -> pd.DataFrame:
    cooler_file = cooler.Cooler(str(path_to_cooler))
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
    df = correct_nchg_pvalue_and_mark_significant_interactions(df, multiple_test_correction_method, log_ratio_thresh,
                                                          fdr_thresh)

    num_significant_interactions = df["significant"].sum()
    logging.info(f"DONE! Identified {num_significant_interactions}/{len(df)} significant trans-regions")
    return df


def process_trans_interactions(path_to_cooler, resolution, segments, masked_regions, odd_log_ratio_thresh, fdr_thresh, nchg_bin, mp_pool):
    trans_df = mp_pool.apply_async(__identify_significant_trans_interactions,
                                   args=[f"{path_to_cooler}::/resolutions/{resolution}",
                                         masked_regions,
                                         odd_log_ratio_thresh,
                                         fdr_thresh,
                                         "fdr_bh",
                                         nchg_bin])
    trans_df = trans_df.get()
    trans_df.to_csv("/tmp/trans_001.tsv", sep="\t", header=False, index=False)

    trans_df = mp_pool.apply_async(__map_interchrom_interactions_to_segments, trans_df[trans_df["significant"]], segments, masked_regions)
    trans_df = trans_df.get()
    trans_df.to_csv("/tmp/trans_002.tsv", sep="\t", header=False, index=False)
    return trans_df


def process_trans_interactions_dbg(path_to_cooler, resolution, segments, masked_regions, odd_log_ratio_thresh, fdr_thresh, nchg_bin, mp_pool):
    nchg_out = pd.read_table("/tmp/nchg_output.tsv", index_col=None)
    trans_df = correct_nchg_pvalue_and_mark_significant_interactions(nchg_out, "fdr_bh", odd_log_ratio_thresh, fdr_thresh)

    trans_df = __map_interchrom_interactions_to_segments(trans_df[trans_df["significant"]], segments)
    trans_df.to_csv("/tmp/trans_002.tsv", sep="\t", header=False, index=False)
    return trans_df
