#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import re
import sys
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
            return self.chrom.partition("_")[0]
        return self.chrom

    def get_chrom_copy(self):
        if "_" in self.chrom:
            return self.chrom.partition("_")[2]
        return "A"


class Interaction3D:
    def __init__(self, bead1, bead2):
        self.bead1, self.bead2 = sorted([bead1, bead2])

    def __iter__(self):
        return iter((self.bead1, self.bead2))

    def __repr__(self):
        return str(self)

    def __str__(self):
        return f"{self.bead1} {self.bead2}"

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)


def build_tad_graph(beads: list, interactions: set) -> networkx.Graph:
    graph = networkx.Graph()
    for bead in beads:
        graph.add_node(bead)

    for bead1, bead2 in interactions:
        graph.add_edge(bead1, bead2)

    return graph


def collect_clique_interactions(cliques, clique_size_thresh: int) -> pd.DataFrame:
    records = []
    clique_id = 0
    for clique in cliques:
        clique_size = len(clique)
        if clique_size < clique_size_thresh:
            continue

        clique_name = f"CLIQUE_#{clique_id}"
        for node in clique:
            records.append(list(node.get_coords()) + [clique_name, clique_size])
        clique_id += 1

    columns = bf.SCHEMAS["bed"][:5]
    return pd.DataFrame(records, columns=columns).sort_values(by=["chrom", "start", "name"])


def preprocess_data(domains: pd.DataFrame, interactions: pd.DataFrame) -> Tuple[set, list]:
    segment_dict = {}
    for interaction in interactions.itertuples(index=False):
        id1 = f"{interaction.chrom1}:{interaction.start1}-{interaction.end1}"
        id2 = f"{interaction.chrom2}:{interaction.start2}-{interaction.end2}"

        segment_dict.setdefault(id1, []).append(id2)
        segment_dict.setdefault(id2, []).append(id1)

    combined_dict = {}

    for segment in domains.itertuples(index=False):
        segment_id = f"{segment.chrom}:{segment.start}-{segment.end}"
        combined_dict[segment_id] = segment_dict.get(segment_id, ".")

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


def import_interactions(path_to_bedpe: pathlib.Path, interaction_type: str, schema: str = "bedpe") -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # Ignore warnings about data loss due to column truncation
        df = bf.read_table(str(path_to_bedpe), schema=schema)

        if interaction_type == "cis-only":
            return df[df["chrom1"] == df["chrom2"]]

        if interaction_type == "trans-only":
            return df[df["chrom1"] != df["chrom2"]]

        assert interaction_type == "all"
        return df


def generate_list_of_domains_with_id(clique_interactions: pd.DataFrame) -> pd.DataFrame:
    domains = (
        clique_interactions[["chrom", "start", "end"]]
        .drop_duplicates()
        .sort_values(by=["chrom", "start"])
        .reset_index(drop=True)
    )
    domains.index.name = "id"

    return domains


def generate_list_of_cliques(clique_interactions: pd.DataFrame, tads: pd.DataFrame) -> pd.DataFrame:
    cliques = clique_interactions.merge(tads.reset_index(), on=["chrom", "start", "end"], how="left")
    cliques = cliques[["id", "name"]].groupby(by=["name"]).aggregate(lambda ids: ",".join(str(i) for i in sorted(ids)))

    cliques.rename(columns={"id": "tad_ids"}, inplace=True)
    if len(cliques) == 0:
        cliques["size"] = pd.Series(dtype=int)
        return cliques

    cliques["size"] = cliques["tad_ids"].str.count(",") + 1

    return cliques


def make_cli() -> argparse.ArgumentParser:
    def existing_file(arg):
        if (file := pathlib.Path(arg)).exists():
            return file

        raise FileNotFoundError(arg)

    def nonnegative_int(arg):
        if (n := int(arg)) >= 0:
            return n

        raise ValueError("Not a non-negative integer")

    cli = argparse.ArgumentParser(description="Call TAD cliques given a list of TADs and significant interactions")

    cli.add_argument("domains", type=existing_file, help="Path to a BED file with a list of TADs.")
    cli.add_argument(
        "interactions",
        type=existing_file,
        help="Path to a BEDPE with the set of significant cis and trans interactions.",
    )
    cli.add_argument("output-prefix", type=pathlib.Path, help="Output path prefix.")
    cli.add_argument(
        "--interaction-type",
        type=str,
        default="all",
        choices={"cis-only", "trans-only", "all"},
        help="Type of interactions to consider when calling cliques.",
    )
    cli.add_argument(
        "--clique-size-threshold",
        type=nonnegative_int,
        help="Minimum clique size. Cliques smaller than this threshold will be dropped.",
        default=5,
    )
    cli.add_argument("--force", action="store_true", default=False, help="Overwrite existing files (if any).")

    return cli


def handle_path_collisions(*paths):
    collisions = [f for f in paths if f.exists()]
    if len(collisions) != 0:
        collisions = "\n - ".join((str(p) for p in collisions))
        raise RuntimeError(
            f"Refusing to overwrite file(s):\n - {collisions}\nPass --force to overwrite existing file(s)."
        )


def main():
    args = vars(make_cli().parse_args())

    domain_file = pathlib.Path(str(args["output-prefix"]) + "_domains.bed")
    clique_file = pathlib.Path(str(args["output-prefix"]) + "_cliques.tsv")

    if not args["force"]:
        handle_path_collisions(domain_file, clique_file)

    domains = import_domains(args["domains"])
    interactions = import_interactions(args["interactions"], args["interaction_type"])
    clique_size_thresh = args["clique_size_threshold"]

    interactions, beads = preprocess_data(domains, interactions)

    tad_graph = build_tad_graph(beads, interactions)
    cliques = tuple(networkx.find_cliques(tad_graph))

    clique_interactions = collect_clique_interactions(cliques, clique_size_thresh)

    domains = generate_list_of_domains_with_id(clique_interactions)
    cliques = generate_list_of_cliques(clique_interactions, domains)

    domains = domains.reset_index()[["chrom", "start", "end", "id"]]
    domains.to_csv(domain_file, sep="\t", index=False, header=False)
    cliques.to_csv(clique_file, sep="\t", index=True, header=True)


if __name__ == "__main__":
    main()
