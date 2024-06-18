#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import sys
from collections import namedtuple
from typing import Tuple, Union

import hictkpy
import pandas as pd

FileURI = namedtuple("FileURI", ["parent", "name", "suffix"])


def make_cli() -> argparse.ArgumentParser:
    def existing_file(arg):
        if (file := pathlib.Path(arg)).exists():
            return file

        raise FileNotFoundError(arg)

    cli = argparse.ArgumentParser()

    cli.add_argument("tsv", type=existing_file, help="Path to a sample sheet.")
    cli.add_argument(
        "--detached",
        action="store_true",
        default=False,
        help="Only validate the sample sheet's syntax (i.e. skip steps that require accessing files listed in the sample sheet)",
    )

    return cli


def parse_file_uri(uri: str) -> FileURI:
    if "::" in uri:
        path, _, suffix = uri.partition("::")
    else:
        path = uri
        suffix = None

    return FileURI(pathlib.Path(path).parent, pathlib.Path(path).name, suffix)


def file_uri_basename(uri: str) -> str:
    _, name, suffix = parse_file_uri(uri)

    if suffix is None:
        return name

    return f"{name}::{suffix}"


def check_sample_ids(df: pd.DataFrame):
    # Look for empty/missing sample IDs
    invalid_ids_mask = (df["sample"].str.strip().str.len() == 0) | (df["sample"].isnull())
    if invalid_ids_mask.sum() != 0:
        affected_rows = "\n- ".join((str(i + 1) for i, invalid in enumerate(invalid_ids_mask) if invalid))
        raise RuntimeError(f"Found {invalid_ids_mask.sum()} missing sample ID(s).\nAffected row(s):\n- {affected_rows}")

    # Look for duplicate sample IDs
    df1 = df.groupby("sample").count()
    duplicates = df1.loc[df1[EXPECTED_COLUMNS[1]] != 1, EXPECTED_COLUMNS[1]]
    if len(duplicates) != 0:
        affected_samples = "\n- ".join(
            (f"{sample}: {num_duplicates} entries" for sample, num_duplicates in duplicates.iteritems())
        )
        raise RuntimeError(f"Found {len(duplicates)} rows with duplicate sample IDs:\n - {affected_samples}")


def check_column_is_integral(col: pd.Series, col_name: Union[str, None] = None):
    if not pd.api.types.is_integer_dtype(col):
        if col_name is None:
            col_name = col.name
        raise RuntimeError(f'column "{col_name}" contains invalid values (expected int, found {col.dtype})')


def check_is_valid_hic_file(path: pathlib.Path, resolution: int, strip_parent_dirs: bool = True):
    if strip_parent_dirs:
        path = file_uri_basename(str(path))
    try:
        f = hictkpy.File(str(path), resolution)
        if f.resolution() != resolution:
            raise RuntimeError(f'expected file "{path}" to have resolution {resolution}, but found {f.resolution()}')
    except OSError as e:
        raise RuntimeError(f'unable to open file "{path}": {e}')
    except Exception as e:
        raise RuntimeError(f'failed to validate file "{path}" at resolution {resolution}: {e}')


def check_is_valid_bed3(path: pathlib.Path, strip_parent_dirs: bool = True):
    if strip_parent_dirs:
        path = file_uri_basename(str(path))
    try:
        df = pd.read_table(path, names=["chrom", "start", "end"], usecols=[0, 1, 2])
        if len(df) == 0:
            return
        if df.isnull().values.any():
            raise RuntimeError("found one or more null/nan value(s) in one of the first three columns")

        check_column_is_integral(df["start"])
        check_column_is_integral(df["end"])

        if (df["start"] >= df["end"]).any():
            raise RuntimeError("found one or more invalid intervals: start position >= end position")
    except RuntimeError as e:
        raise RuntimeError(f'validation failed for file "{path}": {e}')


def main():
    df = pd.read_table(sample_sheet)

    if tuple(df.columns) != EXPECTED_COLUMNS:
        expected = "\n- ".join(EXPECTED_COLUMNS)
        found = "\n- ".join(df.columns)
        raise RuntimeError(f"\nexpected columns:\n- {expected}\nfound:\n- {found}")

    check_sample_ids(df)
    check_column_is_integral(df["resolution"])

    if not detached:
        for _, row in df.iterrows():
            check_is_valid_hic_file(row["hic_file"], row["resolution"])

        for path in df["tads"].fillna(""):
            if path != "":
                check_is_valid_bed3(path)

        for path in df["mask"].fillna(""):
            if path != "":
                check_is_valid_bed3(path)

        paths = [uri.partition("::")[0] for uri in df["hic_file"]]

        df["hic_file"] = paths

    df.to_csv(sys.stdout, sep="\t", index=False, header=True)


if __name__ == "__main__":
    EXPECTED_COLUMNS = tuple(["sample", "hic_file", "resolution", "tads", "mask"])
    args = vars(make_cli().parse_args())
    sample_sheet = args["tsv"]
    detached = args["detached"]

    main()

    print(f'Validation of samplesheet "{sample_sheet}" was successful!', file=sys.stderr)
