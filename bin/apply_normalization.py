#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import logging
import os
import pathlib

import hictkpy
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "hic-file",
        type=pathlib.Path,
        help="Path to a .cool or .hic file (URI syntax supported).",
    )
    cli.add_argument(
        "output-cooler",
        type=pathlib.Path,
        help="Name of the output .cool file.",
    )
    cli.add_argument(
        "--norm-name",
        type=str,
        default="weight",
        help="Name of the normalization to be applied.",
    )
    cli.add_argument(
        "--resolution",
        type=int,
        help="File resolution. Required when input file is in .hic or .mcool format.",
    )

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing files (if any).",
    )
    return cli


def main():
    args = vars(make_cli().parse_args())
    input_file = str(args["hic-file"])
    resolution = args["resolution"]
    output_file = str(args["output-cooler"])

    if (hictkpy.is_hic(input_file) or hictkpy.is_mcool_file(input_file)) and resolution is None:
        raise RuntimeError("--resolution is a mandatory argument when input file is in .mcool or .hic format.")

    fin = hictkpy.File(input_file, resolution)
    sel = fin.fetch(normalization=args["norm_name"], join=False)

    if args["force"]:
        os.remove(output_file)

    writer = hictkpy.cooler.FileWriter(output_file, fin.chromosomes(), fin.resolution())

    chunk_size = int(5e6)

    bin1_ids = []
    bin2_ids = []
    counts = []
    i = 0
    for pixel in sel:
        bin1_ids.append(pixel.bin1_id)
        bin2_ids.append(pixel.bin2_id)
        counts.append(pixel.count)

        if len(bin1_ids) == chunk_size:
            logging.info(f"Writing chunk #{i}")
            i += 1
            writer.add_pixels(pd.DataFrame({"bin1_id": bin1_ids, "bin2_id": bin2_ids, "count": counts}))

            bin1_ids = []
            bin2_ids = []
            counts = []

    if len(bin1_ids) != 0:
        logging.info(f"Writing chunk #{i}")
        writer.add_pixels(pd.DataFrame({"bin1_id": bin1_ids, "bin2_id": bin2_ids, "count": counts}))

    logging.info(f"Finalizing cooler file...")
    writer.finalize()


if __name__ == "__main__":
    main()
