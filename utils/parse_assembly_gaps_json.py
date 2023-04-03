#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import json
import sys

import pandas as pd

# This script reads assembly gaps in JSON format from stdin and outputs the gaps in BED format to stdout.
# Gaps in JSON format can be downloaded from UCSC.
# e.g. https://api.genome.ucsc.edu/getData/track?genome=hg38;track=gap;jsonOutputArrays=1;maxItemsOutput=-1

if __name__ == "__main__":
    records = []
    for array in json.load(sys.stdin)["gap"].values():
        records.extend([[str(r[1]), int(r[2]), int(r[3])] for r in array])

    pd.DataFrame(records, columns=["chrom", "start", "end"]).sort_values(
        by=["chrom", "start", "end"]
    ).to_csv(sys.stdout, sep="\t", header=False, index=False)
