#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import argparse
import itertools
import logging
import pathlib
import re
import warnings
from typing import Tuple

import bioframe as bf
import networkx
import pandas as pd


class Bead3D:
    def __init__(self, chrom, start, end, radius, id):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.pos = abs(self.end + self.start) / 2
        self.radius = float(radius)
        self.id = str(id)

    def __eq__(self, other):
        return self.chrom == other.chrom and self.start == other.start and self.end == other.end

    def __lt__(self, other):
        if self.chrom == other.chrom:
            return self.pos < other.pos

        return self.chrom < other.chrom

    def __gt__(self, other):
        if self.chrom == other.chrom:
            return self.pos > other.pos

        return self.chrom > other.chrom

    def __repr__(self):
        return self.id

    def __str__(self):
        return self.id

    def __hash__(self):
        return hash(f"{self.chrom}:{self.start}-{self.end}")

    def __sub__(self, other):
        return abs(self.pos - other.pos)

    def get_BED(self, strip_copy=False):  # noqa
        return f"{self.get_chrom(strip_copy)}\t{self.start}\t{self.end}"

    def get_coords(self, strip_copy=False):
        return self.get_chrom(strip_copy), self.start, self.end

    def get_chrom(self, strip_copy=False):
        if strip_copy:
            return self.chrom.split("_")[0]
        return self.chrom

    def get_chrom_copy(self):
        if "_" in self.chrom:
            return self.chrom.split("_")[1]
        return "A"


class Interaction3D:
    def __init__(self, bead1, bead2):
        self.bead1, self.bead2 = sorted([bead1, bead2])

    def __iter__(self):
        return iter((self.bead1, self.bead2))

    def __repr__(self):
        return str(self)

    def __str__(self):
        return str(self.bead1) + " " + str(self.bead2)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)


def build_tad_graph(beads: list, interactions: set, chrom: str) -> networkx.Graph:
    graph = networkx.Graph()
    for bead in beads:
        if bead.chrom == chrom:
            graph.add_node(bead)

    for bead1, bead2 in interactions:
        if bead1.chrom == bead2.chrom == chrom:
            graph.add_edge(bead1, bead2)

    return graph


def compute_clique_stats(cliques, clique_size_thresh: int) -> pd.DataFrame:
    records = []
    for i, clique in enumerate(cliques):
        if len(clique) < clique_size_thresh:
            continue

        clique_str = ";".join((str(x) for x in sorted(clique)))
        records.append([clique[0].get_chrom(True), i, len(clique), f"CLICSTAT # {clique_str}"])

    return pd.DataFrame(records, columns=["chrom", "bead1_id", "bead2_id", "comment"]).sort_values(
        by=["chrom", "bead1_id", "bead2_id"]
    )


def map_clique_interactions(cliques, clique_size_thresh: int) -> Tuple[set, pd.DataFrame]:
    clique_interactions = set()
    for clique in cliques:
        if len(clique) < clique_size_thresh:
            continue

        for node1, node2 in itertools.product(clique, repeat=2):
            clique_interactions.add(Interaction3D(node1, node2))
    records = [list(n1.get_coords()) + list(n2.get_coords()) for n1, n2 in clique_interactions]

    columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
    return clique_interactions, pd.DataFrame(records, columns=columns).sort_values(
        by=["chrom1", "start1", "chrom2", "start2"]
    )


def map_tad_interactions(interactions: set, clique_interactions: set, chrom: str) -> pd.DataFrame:
    records = []
    for bead1, bead2 in interactions:
        if bead1.chrom == chrom and bead2.chrom == chrom:
            chrom1, start1, end1 = bead1.get_coords()
            chrom2, start2, end2 = bead2.get_coords()
            bead_belongs_to_clique = Interaction3D(bead1, bead2) in clique_interactions

            records.append([chrom1, start1, end1, chrom2, start2, end2, bead_belongs_to_clique, "INTERACT"])

    columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "bead_part_of_clique", "comment"]
    return pd.DataFrame(records, columns=columns).sort_values(by=["chrom1", "start1", "chrom2", "start2"])


def compute_clique_sizes(tad_graph, clique_size_thresh: int) -> pd.DataFrame:
    cliquenum = networkx.node_clique_number(tad_graph)
    records = [list(tad.get_coords(True)) + [cliquenum[tad]] for tad in cliquenum]

    columns = ["chrom", "start", "end", "size"]
    df = pd.DataFrame(records, columns=columns).sort_values(by=["chrom", "start"])

    return df[df["size"] >= clique_size_thresh]


def preprocess_data(domains: pd.DataFrame, interactions: pd.DataFrame) -> Tuple[set, list]:
    segment_dict = {}
    for interaction in interactions.itertuples(index=False):
        id1 = f"{interaction.chrom1}:{interaction.start1}-{interaction.end1}"
        id2 = f"{interaction.chrom2}:{interaction.start2}-{interaction.end2}"

        segment_dict.setdefault(id1, list()).append(id2)
        segment_dict.setdefault(id2, list()).append(id1)

    combined_dict = {}

    for segment in domains.itertuples(index=False):
        segment_id = f"{segment.chrom}:{segment.start}-{segment.end}"
        combined_dict[segment_id] = segment_dict[segment_id] if segment_id in segment_dict else "."

    records = []
    pattern = re.compile(r"[:-]")
    for segment_id, edges in combined_dict.items():
        chrom, start, end = pattern.split(segment_id)
        records.append([str(chrom), int(start), int(end), str(segment_id), 1, edges])

    columns = ["chrom", "start", "end", "id", "radius", "edges"]
    data = pd.DataFrame(records, columns=columns).sort_values(by=["chrom", "start"])

    res = {}
    idpairs = []
    id2bead = {}
    beads = []

    for row in data.itertuples(index=False):
        bead = Bead3D(row.chrom, row.start, row.end, row.radius, row.id)

        id2bead[row.id] = bead
        beads.append(bead)

        if row.edges != ".":
            idpairs.extend([tuple([row.id, id2]) for id2 in row.edges])

        res.setdefault(row.chrom, []).append(bead)

    for beads in res.values():
        beads.sort()

    interacting_beads = {Interaction3D(id2bead[id1], id2bead[id2]) for id1, id2 in idpairs}

    return interacting_beads, beads


def import_domains(path_to_bed: pathlib.Path, schema: str = "bed3") -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # Ignore warnings about data loss when parsing BED4+
        return bf.read_table(str(path_to_bed), schema=schema)


def import_interactions(path_to_bedpe: pathlib.Path, schema: str = "bedpe") -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # Ignore warnings about data loss due to column truncation
        return bf.read_table(str(path_to_bedpe), schema=schema)


def setup_logger(level=logging.INFO):
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s")
    logging.getLogger().setLevel(level)


def make_cli() -> argparse.ArgumentParser:
    def existing_file(arg):
        if (file := pathlib.Path(arg)).exists():
            return file

        raise FileNotFoundError(arg)

    def nonnegative_int(arg):
        if (n := int(arg)) >= 0:
            return n

        raise ValueError("Not a non-negative integer")

    cli = argparse.ArgumentParser(
        description="Map interchromosomal (trans) interactions between a list of domains in BED format"
    )

    cli.add_argument("domains", type=existing_file, help="Path to a BED file with a list of TADs.")
    cli.add_argument(
        "interactions",
        type=existing_file,
        help="Path to a BEDPE with the set of significant cis and trans interactions.",
    )
    cli.add_argument("output_prefix", type=pathlib.Path, help="Path to output prefix (including parent folder(s)).")
    cli.add_argument("--force", action="store_true", default=False, help="Force overwrite existing files.")
    cli.add_argument(
        "--clique-size-threshold",
        type=nonnegative_int,
        help="Minimum clique size. Cliques smaller than this threshold will be dropped.",
        default=5,
    )

    return cli


def main():
    args = vars(make_cli().parse_args())

    domains = import_domains(args["domains"])
    interactions = import_interactions(args["interactions"])
    clique_size_thresh = args["clique_size_threshold"]
    out_prefix = args["output_prefix"]

    if (parent := out_prefix.parent) != "" and parent != ".":
        parent.mkdir(exist_ok=True)

    interactions, beads = preprocess_data(domains, interactions)

    for suffix in [
        "_clique_stats.tsv",
        "_clique_sizes.bedGraph",
        "_clique_interactions.bedpe",
        "_tad_interactions.bedpe",
    ]:
        if (file := pathlib.Path(f"{out_prefix}{suffix}")).exists():
            if args["force"]:
                file.unlink()
            else:
                raise RuntimeError(f"Refusing to overwrite existing file {file}")

    print_header = True
    for chrom in domains["chrom"].unique():
        logging.info("Processing %s...", chrom)
        tad_graph = build_tad_graph(beads, interactions, chrom)
        cliques = tuple(networkx.find_cliques(tad_graph))

        clique_stats_df = compute_clique_stats(cliques, clique_size_thresh)
        clique_interactions, clique_interactions_df = map_clique_interactions(cliques, clique_size_thresh)
        tad_interactions_df = map_tad_interactions(interactions, clique_interactions, chrom)
        clique_sizes_df = compute_clique_sizes(tad_graph, clique_size_thresh)

        clique_stats_df.to_csv(f"{out_prefix}_clique_stats.tsv", index=False, header=print_header, sep="\t", mode="a")
        clique_sizes_df.to_csv(f"{out_prefix}_clique_sizes.bedGraph", index=False, header=False, sep="\t", mode="a")
        clique_interactions_df.to_csv(
            f"{out_prefix}_clique_interactions.bedpe", index=False, header=False, sep="\t", mode="a"
        )
        tad_interactions_df.to_csv(
            f"{out_prefix}_tad_interactions.bedpe", index=False, header=print_header, sep="\t", mode="a"
        )

        print_header = False


if __name__ == "__main__":
    setup_logger()
    main()
