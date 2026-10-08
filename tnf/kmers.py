"""Tetranucleotide frequency (TNF) featurization.

A DNA strand is read 5'->3', but a genome is double stranded: the k-mer AAGT on
one strand is ACTT on the other. Counting each 4-mer together with its reverse
complement ("canonical" k-mers) makes the vector strand-independent. Of the
256 possible 4-mers, 16 are their own reverse complement (palindromes), so the
canonical space has (256 - 16) / 2 + 16 = 136 dimensions.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Iterator

import numpy as np

K = 4
_BASES = "ACGT"
_ENC = np.full(256, -1, dtype=np.int8)
for _i, _b in enumerate(_BASES):
    _ENC[ord(_b)] = _i
    _ENC[ord(_b.lower())] = _i


def _revcomp_code(code: int) -> int:
    """Reverse complement of a 4-mer stored as a base-4 integer (A0 C1 G2 T3)."""
    digits = [(code >> (2 * (K - 1 - j))) & 3 for j in range(K)]
    rc = [3 - d for d in reversed(digits)]
    out = 0
    for d in rc:
        out = out * 4 + d
    return out


def _build_canonical_table():
    canon_of = np.empty(4 ** K, dtype=np.int16)
    names: list[str] = []
    seen: dict[int, int] = {}
    for code in range(4 ** K):
        rep = min(code, _revcomp_code(code))
        if rep not in seen:
            seen[rep] = len(names)
            names.append("".join(_BASES[(rep >> (2 * (K - 1 - j))) & 3] for j in range(K)))
        canon_of[code] = seen[rep]
    return canon_of, names


CANON_INDEX, CANON_KMERS = _build_canonical_table()
N_FEATURES = len(CANON_KMERS)  # 136
assert N_FEATURES == 136


# --------------------------------------------------------------------------- io
def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield (record_id, sequence) from a plain or gzipped FASTA file."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    name, chunks = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


# ----------------------------------------------------------------- counting
def kmer_index_track(seq: str) -> np.ndarray:
    """Canonical 4-mer index (0..135) starting at every position; -1 where the
    4-mer contains a non-ACGT base (N, IUPAC codes) or runs off the end."""
    s = _ENC[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
    n = len(s) - K + 1
    if n <= 0:
        return np.empty(0, dtype=np.int16)
    valid = s >= 0
    ok = valid[:n].copy()
    for j in range(1, K):
        ok &= valid[j : j + n]
    sc = np.where(valid, s, 0).astype(np.int32)
    code = np.zeros(n, dtype=np.int32)
    for j in range(K):
        code = code * 4 + sc[j : j + n]
    track = CANON_INDEX[code].astype(np.int16)
    track[~ok] = -1
    return track


def count_kmers(seq: str) -> np.ndarray:
    """136-dim canonical 4-mer counts for a whole sequence."""
    t = kmer_index_track(seq)
    t = t[t >= 0]
    return np.bincount(t, minlength=N_FEATURES).astype(np.float64)


def window_counts(seq: str, window: int, stride: int):
    """Counts for sliding windows along one contig.

    Uses a block trick for speed: k-mers are binned into stride-sized blocks by
    start position, and each window is the sum of window/stride consecutive
    blocks. (The last 3 k-mers of a window may extend up to 3 bp past its end,
    which is negligible.) Returns (starts, counts[n_windows, 136]).
    """
    if window % stride != 0:
        raise ValueError("window must be a multiple of stride")
    L = len(seq)
    if L < window:
        return np.empty(0, dtype=np.int64), np.empty((0, N_FEATURES))
    t = kmer_index_track(seq)
    pos = np.nonzero(t >= 0)[0]
    n_blocks = L // stride
    blk = pos // stride
    keep = blk < n_blocks
    flat = blk[keep] * N_FEATURES + t[pos[keep]]
    blocks = np.bincount(flat, minlength=n_blocks * N_FEATURES).reshape(n_blocks, N_FEATURES)
    m = window // stride
    csum = np.vstack([np.zeros((1, N_FEATURES)), np.cumsum(blocks, axis=0)])
    n_win = n_blocks - m + 1
    counts = csum[m : m + n_win] - csum[:n_win]
    starts = np.arange(n_win, dtype=np.int64) * stride
    return starts, counts


# ------------------------------------------------------------ normalization
def clr(counts: np.ndarray, pseudocount: float = 0.5) -> np.ndarray:
    """Centered log-ratio transform. TNF is compositional (only proportions
    carry information), so we compare log-ratios rather than raw frequencies."""
    x = np.log(np.asarray(counts, dtype=np.float64) + pseudocount)
    return x - x.mean(axis=-1, keepdims=True)


# -------------------------------------------------------------- genome level
@dataclass
class WindowTable:
    genome: str
    contigs: np.ndarray   # contig id per window
    starts: np.ndarray    # 0-based start
    ends: np.ndarray      # exclusive end
    counts: np.ndarray    # [n, 136]


def featurize_genome(path: str | Path, window: int | None = None, stride: int | None = None,
                     genome_id: str | None = None):
    """Return the whole-genome 136-dim count vector and, if window is given,
    a WindowTable of per-window counts (windows never span contig boundaries)."""
    path = Path(path)
    gid = genome_id or path.name.split(".")[0]
    total = np.zeros(N_FEATURES)
    cids, sts, cnts = [], [], []
    for cid, seq in read_fasta(path):
        total += count_kmers(seq)
        if window:
            st, c = window_counts(seq, window, stride or window)
            cids.append(np.full(len(st), cid, dtype=object))
            sts.append(st)
            cnts.append(c)
    if not window:
        return total, None
    starts = np.concatenate(sts) if sts else np.empty(0, dtype=np.int64)
    table = WindowTable(
        genome=gid,
        contigs=np.concatenate(cids) if cids else np.empty(0, dtype=object),
        starts=starts,
        ends=starts + window,
        counts=np.vstack(cnts) if cnts else np.empty((0, N_FEATURES)),
    )
    return total, table


__all__ = [
    "CANON_KMERS", "N_FEATURES", "read_fasta", "kmer_index_track", "count_kmers",
    "window_counts", "clr", "featurize_genome", "WindowTable",
]
