# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import itertools
import logging
import pathlib
import tempfile
import time
from typing import Union

import cooler
import multiprocess as mp
import pandas as pd

from .run_nchg import run_nchg_and_generate_df, correct_nchg_pvalue_and_mark_significant_interactions


def __map_intrachrom_interactions_to_domains_worker(args: list) -> list:
    def segment_to_str(segment):  # noqa
        return f"{segment.chrom}:{segment.start}-{segment.end}"

    path_to_cooler, i, domains = args  # Unpack params
    logging.debug(
        f"Mapping pairwise cis-interactions for genomic regions between {segment_to_str(domains[i])} and {segment_to_str(domains[-1])}"
    )

    query = domains[i]

    cooler_file = cooler.Cooler(str(path_to_cooler))
    chrom_size = cooler_file.chromsizes[query.chrom]

    # Fetch and index pixel corresponding to rows overlapping the query segment
    # Pixels are indexed by start2 and end2 for performance reasons
    pixels = (
        cooler_file.matrix(balance=False, as_pixels=True, join=True)
        .fetch(f"{query.chrom}:{query.start}-{query.end}", f"{query.chrom}:{query.start}-{chrom_size}")
        .set_index(["start2", "end2"], drop=False)
        .sort_index()
        .drop(columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2"])
    )

    results = []  # Loop over domains located downstream of the query segment
    for segment in domains[i + 1 :]:
        # Slice pixels df and aggregate contacts
        count = pixels.loc[
            (pixels.index.get_level_values("start2") >= segment.start)
            & (pixels.index.get_level_values("end2") < segment.end),
            "count",
        ].sum()
        if count != 0:
            results.append([query.chrom, query.start, query.end, query.chrom, segment.start, segment.end, count])
    logging.debug(
        f"DONE mapping pairwise cis-interactions for genomic regions between {segment_to_str(domains[i])} and {segment_to_str(domains[-1])} "
        f"({len(results)} nnz out of {len(domains[i + 1:])} regions)"
    )
    return results


def __map_intrachrom_interactions_to_domains(
    cooler_uri: str, domains: pd.DataFrame, thread_pool: mp.Pool
) -> pd.DataFrame:  # noqa
    """
    Generate a BEDPE dataframe with the number of cis-interactions between all domains
    """
    t0 = time.time()
    logging.info(f"Mapping pairwise cis-interactions between {len(domains)} regions...")
    async_results = []
    for chrom, df in domains.groupby("chrom"):
        # Asyncronously loop over domains overlapping a given chromosome
        queries = tuple(df.itertuples(index=False))
        async_results.append(
            thread_pool.map_async(
                __map_intrachrom_interactions_to_domains_worker,
                zip(itertools.repeat(cooler_uri), range(len(queries)), itertools.repeat(queries)),
            )
        )

    records = []
    for result in async_results:
        for chunk in result.get():
            records.extend((res for res in chunk if len(res) != 0))

    df = pd.DataFrame(records, columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2", "score"])
    t1 = time.time()
    nnz_regions = len(df)
    tot_regions = (len(domains) ** 2) - len(domains)
    elapsed_time = time.strftime("%Hh%Mm%Ss", time.gmtime(t1 - t0))
    logging.info(
        f"Mapping of pairwise cis-interactions took {elapsed_time}. "
        f"Found {nnz_regions}/{tot_regions} regions with one or more interactions"
    )
    return df


def __identify_significant_cis_interactions(
    domains: pd.DataFrame,
    resolution: int,
    log_ratio_thresh: float,
    fdr_thresh: float,
    multiple_test_correction_method: str,
    nchg_bin: Union[pathlib.Path, str],
) -> pd.DataFrame:
    logging.info("Identifying significant cis-regions...")
    # We are using a temp file because NCHG returns incomplete results when the input BEDPE
    # is not a regular file
    with tempfile.NamedTemporaryFile(mode="w+t") as tmpbedpe:
        # Write domains to tmp file
        domains.to_csv(tmpbedpe, sep="\t", header=False, index=False)
        tmpbedpe.flush()

        df = run_nchg_and_generate_df(tmpbedpe.name, nchg_bin, method="intra", resolution=resolution)

    assert len(df) > 0
    df = correct_nchg_pvalue_and_mark_significant_interactions(
        df, multiple_test_correction_method, log_ratio_thresh, fdr_thresh
    )

    num_significant_interactions = df["significant"].sum()
    logging.info(f"DONE! Identified {num_significant_interactions}/{len(df)} significant cis-regions")
    return df


def process_cis_interactions(
    cf: cooler.Cooler,
    domains: pd.DataFrame,
    odd_log_ratio_thresh: float,
    fdr_thresh: float,
    nchg_bin: Union[pathlib.Path, str],
    exec_pool,
):
    cis_df = __map_intrachrom_interactions_to_domains(cf.uri, domains, exec_pool)

    return exec_pool.apply(
        __identify_significant_cis_interactions,
        args=[cis_df, cf.binsize, odd_log_ratio_thresh, fdr_thresh, "fdr_bh", nchg_bin],
    )
