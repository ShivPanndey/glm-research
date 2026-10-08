"""Parse functional labels and attach them to genome coordinates.

Inputs we expect from the lab's pipeline:
  * eggNOG-mapper  *.emapper.annotations  (21 TSV columns, KEGG_ko is column 12)
  * Prodigal       *.faa  (header carries "# start # end # strand # ID=...")
                   or *.gff (CDS rows with coordinates)

The emapper file names each ORF (e.g. "contig7_12" = 12th gene on contig7) but
has no coordinates. Prodigal output supplies the coordinates, so joining the two
gives ORF -> (contig, start, end, KOs), which is what window-level retrieval needs.
"""
from __future__ import annotations

import gzip
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

KEGG_KO_COL = 11  # 0-based index of column 12


def _open(path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def parse_emapper(path: str | Path) -> dict[str, set[str]]:
    """ORF id -> set of KOs (e.g. {'K00001'}). ORFs with KEGG_ko '-' are kept
    with an empty set so we can count unannotated genes."""
    orf_kos: dict[str, set[str]] = {}
    ko_col = KEGG_KO_COL
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("#query"):
                    cols = line[1:].rstrip("\n").split("\t")
                    if "KEGG_ko" in cols:
                        ko_col = cols.index("KEGG_ko")
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) <= ko_col:
                continue
            raw = f[ko_col].strip()
            kos = set()
            if raw and raw != "-":
                kos = {k.replace("ko:", "").strip() for k in raw.split(",") if k.strip()}
            orf_kos[f[0]] = kos
    return orf_kos


@dataclass
class ORF:
    orf_id: str
    contig: str
    start: int   # 0-based inclusive
    end: int     # 0-based exclusive
    strand: int


def parse_prodigal_faa(path: str | Path) -> list[ORF]:
    """Coordinates from Prodigal .faa headers:
    >contig7_12 # 10432 # 11301 # -1 # ID=7_12;partial=00;..."""
    orfs = []
    with _open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = [p.strip() for p in line[1:].split("#")]
            orf_id = parts[0].split()[0]
            start, end, strand = int(parts[1]), int(parts[2]), int(parts[3])
            contig = orf_id.rsplit("_", 1)[0]
            orfs.append(ORF(orf_id, contig, start - 1, end, strand))
    return orfs


def parse_prodigal_gff(path: str | Path) -> list[ORF]:
    """Coordinates from Prodigal GFF (ID=<contig_idx>_<gene_idx>; ORF name is
    rebuilt as <seqid>_<gene_idx> to match emapper query names)."""
    orfs = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            attrs = dict(a.split("=", 1) for a in f[8].strip(";").split(";") if "=" in a)
            gene_idx = attrs.get("ID", "").rsplit("_", 1)[-1]
            orfs.append(ORF(f"{f[0]}_{gene_idx}", f[0], int(f[3]) - 1, int(f[4]),
                            1 if f[6] == "+" else -1))
    return orfs


def genome_ko_set(orf_kos: dict[str, set[str]]) -> set[str]:
    out: set[str] = set()
    for kos in orf_kos.values():
        out |= kos
    return out


def window_ko_sets(contigs: np.ndarray, starts: np.ndarray, ends: np.ndarray,
                   orfs: list[ORF], orf_kos: dict[str, set[str]],
                   min_overlap: float = 0.5) -> list[set[str]]:
    """KOs assigned to each window: an ORF contributes its KOs to a window if at
    least `min_overlap` of the ORF's length lies inside the window."""
    by_contig: dict[str, list[ORF]] = defaultdict(list)
    for o in orfs:
        if orf_kos.get(o.orf_id):
            by_contig[o.contig].append(o)
    index = {}
    for c, lst in by_contig.items():
        lst.sort(key=lambda o: o.start)
        index[c] = (np.array([o.start for o in lst]), np.array([o.end for o in lst]), lst)
    out = []
    for c, s, e in zip(contigs, starts, ends):
        kos: set[str] = set()
        if c in index:
            os_, oe, lst = index[c]
            hit = np.nonzero((os_ < e) & (oe > s))[0]
            for i in hit:
                ov = min(oe[i], e) - max(os_[i], s)
                if ov / (oe[i] - os_[i]) >= min_overlap:
                    kos |= orf_kos[lst[i].orf_id]
        out.append(kos)
    return out


class KOVocab:
    """Maps KO strings to column indices of a multi-label matrix."""

    def __init__(self, ko_sets):
        vocab = sorted(set().union(*ko_sets)) if ko_sets else []
        self.kos = vocab
        self.idx = {k: i for i, k in enumerate(vocab)}

    def __len__(self):
        return len(self.kos)

    def to_matrix(self, ko_sets) -> np.ndarray:
        Y = np.zeros((len(ko_sets), len(self.kos)), dtype=np.float32)
        for r, s in enumerate(ko_sets):
            for k in s:
                j = self.idx.get(k)
                if j is not None:
                    Y[r, j] = 1.0
        return Y
