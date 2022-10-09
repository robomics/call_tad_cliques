# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import json
import logging
import pathlib
import warnings
from typing import Union

import bioframe as bf
import cooler
import pandas as pd


# https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
def import_centromeric_regions(path: Union[pathlib.Path, str]) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df = bf.read_table(path, schema="bed4")
        return df[df["name"] == "acen"].drop(columns="name")


# https://api.genome.ucsc.edu/getData/track?genome=hg38;track=gap;jsonOutputArrays=1
def __array_json_to_bed(path: Union[pathlib.Path, str]) -> list:
    # Parse JSON array(s) into a flat list of BED3 records
    records = []
    with open(path, "r") as f:
        for array in json.load(f)["gap"].values():
            records.extend([[str(r[1]), int(r[2]), int(r[3])] for r in array])
    return records


def import_assembly_gaps(path: Union[pathlib.Path, str], format: str = "json") -> pd.DataFrame:  # noqa
    supported_formats = {"json", "bed"}

    if format not in supported_formats:
        raise RuntimeError(
            f"Unknown format \"{format}\". Supported formats: " + ", ".join([f for f in supported_formats]))

    if format == "json":
        return bf.from_list(__array_json_to_bed(path), cols=["chrom", "start", "end"])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return bf.read_table(path, schema="bed3")


def import_chrom_sizes(path_to_cooler: Union[pathlib.Path, str]) -> pd.DataFrame:
    logging.info(f"Reading chrom sizes from {path_to_cooler}...")
    chroms = bf.from_dict(cooler.Cooler(str(path_to_cooler)).chromsizes)
    logging.info(f"Imported {len(chroms)} chromosomes")
    return chroms


def import_tads(path: Union[pathlib.Path, str]) -> pd.DataFrame:
    logging.info(f"Reading TADs from {path}...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tads = bf.read_table(str(path), schema="bed3")

    logging.info(f"Imported {len(tads)}")
    return tads
