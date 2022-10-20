# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import logging

import bioframe as bf
import numpy as np
import pandas as pd


def fill_gaps_between_tads(tad_df: pd.DataFrame, chrom_sizes_df: pd.DataFrame) -> pd.DataFrame:
    """
    Generate a dataframe of contiguous domains where gaps between TAD domains are filled with non-TAD domains
    """
    logging.info("Filling gaps between TADs...")
    df1 = tad_df.copy()
    df1["name"] = "tad"

    df2 = bf.complement(tad_df, chrom_sizes_df)
    df2["name"] = "non-tad"

    domains = pd.concat([df1, df2])

    # Map chromosome names to their rank in the chrom_sizes_df
    mappings = {row.chrom: i for i, row in enumerate(chrom_sizes_df.itertuples(index=False))}  # noqa

    # Assign a chromosome rank to each entry in the segment df
    domains["chrom_rank"] = np.array([mappings.get(chrom) for chrom in domains["chrom"]], dtype=float)
    domains[~np.isfinite(domains["chrom_rank"])] = len(chrom_sizes_df) + 1

    # Sort domains such that chromosomes appear in the same order as in the chrom_sizes_df
    domains = (
        domains.sort_values(by=["chrom_rank", "start", "end"])
        .drop(columns=["chrom_rank", "view_region"])
        .reset_index(drop=True)
    )

    logging.info(f"Filled {len(domains) - len(tad_df)} gaps")
    return domains


def mask_cis_regions(regions: pd.DataFrame, mask: pd.DataFrame):
    logging.info(f"Masking {len(regions)} cis-regions...")
    df = bf.setdiff(regions, mask)
    logging.info(f"Masked {len(regions) - len(df)} regions")
    return df
