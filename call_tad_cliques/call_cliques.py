# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT
import itertools
import pathlib
import re
from collections import namedtuple
from typing import Union, Tuple

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

    def getBED(self, stripCopy=False):
        return f"{self.getChrom(stripCopy)}\t{self.start}\t{self.end}"

    def getCoords(self, stripCopy=False):
        return self.getChrom(stripCopy), self.start, self.end

    def getChrom(self, stripCopy=False):
        if stripCopy:
            return self.chrom.split("_")[0]
        return self.chrom

    def getChromCopy(self):
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


def __parse_gtrack_header(line: str,
                          required_fields=tuple(["seqid", "start", "end",
                                                 "id", "radius", "edges"])):
    header = tuple(line.replace("###", "").split())
    for field in required_fields:
        if field not in header:
            raise RuntimeError(f"header is missing mandatory field {field}.")

    return header


__GTrackRecord = namedtuple("__GTrackRecord", ["chrom", "start", "end", "id", "radius"])


def __parse_gtrack_record(line: str, header: tuple) -> Tuple[__GTrackRecord,
                                                             Bead3D,
                                                             Union[tuple, None]]:
    record = dict(zip(header, line.split()))

    edges = record["edges"]
    if edges != ".":
        edges = tuple(edges.split(";"))
    else:
        edges = None

    record = __GTrackRecord(record["seqid"],
                            record["start"],
                            record["end"],
                            record["id"],
                            record["radius"])

    bead = Bead3D(record.chrom,
                  record.start,
                  record.end,
                  record.radius,
                  record.id)

    return record, bead, edges


def __read_gtrack_file(path_to_gtrack_file: Union[pathlib.Path, str]) -> Tuple[set, list]:
    """ Reads data as GTrack file
        OBS: Does no actual sanity check on the file format
    """
    res = {}
    idpairs = []
    id2bead = {}
    beads = []
    header_ids = None
    with open(path_to_gtrack_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                if line.startswith("###"):  # File header
                    assert header_ids is None
                    header_ids = __parse_gtrack_header(line)
                continue

            if not header_ids:
                raise RuntimeError(f"GTrack file {path_to_gtrack_file} is missing the file header")

            record, bead, edges = __parse_gtrack_record(line, header_ids)

            id2bead[record.id] = bead
            beads.append(bead)

            if edges is not None:
                idpairs.extend([tuple([record.id, id2]) for id2 in edges])

            if record.chrom not in res:
                res[record.chrom] = [bead]
            else:
                res[record.chrom].append(bead)

    [beads.sort() for beads in res.values()]
    interacting_beads = {Interaction3D(id2bead[id1], id2bead[id2]) for id1, id2 in idpairs}

    return interacting_beads, beads


def __build_tad_graph(beads: list, interactions: set, chrom: str) -> networkx.Graph:
    graph = networkx.Graph()
    [graph.add_node(bead) for bead in beads if bead.chrom == chrom]
    [graph.add_edge(bead1, bead2) for bead1, bead2 in interactions if bead1.chrom == bead2.chrom == chrom]

    return graph


def __compute_clique_stats(cliques) -> pd.DataFrame:
    records = []
    for i, clique in enumerate(cliques):
        clique_str = ";".join((str(x) for x in sorted(clique)))
        records.append([clique[0].getChrom(True), i, len(clique), f"CLICSTAT # {clique_str}"])

    return pd.DataFrame(records, columns=["chrom", "bead1_id", "bead2_id", "comment"]) \
        .sort_values(by=["chrom", "bead1_id", "bead2_id"])


def __map_clique_interactions(cliques, clique_size_thresh: int) -> Tuple[set, pd.DataFrame]:
    clique_interactions = set()
    for i, clique in enumerate(cliques):
        if len(clique) < clique_size_thresh:
            continue
        for node1, node2 in itertools.product(clique, repeat=2):
            clique_interactions.add(Interaction3D(node1, node2))
    records = [list(n1.getCoords()) + list(n2.getCoords()) for n1, n2 in clique_interactions]

    columns = ["chrom1", "start1", "end1",
               "chrom2", "start2", "end2"]
    return clique_interactions, \
           pd.DataFrame(records, columns=columns).sort_values(by=["chrom1", "start1", "chrom2", "start2"])


def __map_tad_interactions(interactions: set, clique_interactions: set, chrom: str) -> pd.DataFrame:
    records = []
    for bead1, bead2 in interactions:
        if bead1.chrom == chrom and bead2.chrom == chrom:
            bead_belongs_to_clique = Interaction3D(bead1, bead2) in clique_interactions
            records.append([bead1.id, bead2.id, bead_belongs_to_clique, "INTERACT"])

    columns = ["bead1_id", "bead2_id", "bead_part_of_clique", "comment"]
    return pd.DataFrame(records, columns=columns).sort_values(by=["bead1_id", "bead2_id"])


def __compute_clique_sizes(tad_graph) -> pd.DataFrame:
    cliquenum = networkx.node_clique_number(tad_graph)
    records = [list(tad.getCoords(True)) + [cliquenum[tad]] for tad in cliquenum]

    columns = ["chrom", "start", "end", "size"]
    return pd.DataFrame(records, columns=columns).sort_values(by=["chrom", "start"])


def __preprocess_data(segments: pd.DataFrame, interactions: pd.DataFrame) -> Tuple[set, list]:
    segment_dict = {}
    for interaction in interactions.itertuples(index=False):
        id1 = f"{interaction.chrom1}:{interaction.start1}-{interaction.end1}"
        id2 = f"{interaction.chrom2}:{interaction.start2}-{interaction.end2}"

        segment_dict.setdefault(id1, list()).append(id2)
        segment_dict.setdefault(id2, list()).append(id1)

    combined_dict = {}

    for segment in segments.itertuples(index=False):
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
        bead = Bead3D(row.chrom,
                      row.start,
                      row.end,
                      row.radius,
                      row.id)

        id2bead[row.id] = bead
        beads.append(bead)

        if row.edges != ".":
            idpairs.extend([tuple([row.id, id2]) for id2 in row.edges])

        res.setdefault(row.chrom, []).append(bead)

    [beads.sort() for beads in res.values()]
    interacting_beads = {Interaction3D(id2bead[id1], id2bead[id2]) for id1, id2 in idpairs}

    return interacting_beads, beads


def call_cliques(segments: pd.DataFrame, interactions: pd.DataFrame, chrom: str, clique_size_thresh: int) -> Tuple[pd.DataFrame,
                                                                                                                   pd.DataFrame,
                                                                                                                   pd.DataFrame,
                                                                                                                   pd.DataFrame]:
    interactions, beads = __preprocess_data(segments, interactions)
    # interactions, beads = __read_gtrack_file(path_to_gtrack)

    tad_graph = __build_tad_graph(beads, interactions, chrom)
    cliques = tuple(networkx.find_cliques(tad_graph))

    clique_stats_df = __compute_clique_stats(cliques)
    clique_interactions, clique_interactions_df = __map_clique_interactions(cliques, clique_size_thresh)
    tad_interactions_df = __map_tad_interactions(interactions, clique_interactions, chrom)
    clique_sizes_df = __compute_clique_sizes(tad_graph)

    return clique_stats_df, clique_interactions_df, tad_interactions_df, clique_sizes_df
