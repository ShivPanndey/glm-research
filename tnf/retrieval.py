"""Vector database over TNF vectors, plus KO label transfer by nearest neighbors.

Retrieval idea: a query genome (or a window of it) is embedded as a TNF vector;
we look up its k most similar labeled references and predict each KO with a
similarity-weighted vote over the references that carry it.
"""
from __future__ import annotations

import faiss
import numpy as np

from .kmers import clr


class TNFIndex:
    """CLR -> z-score (fit on references) -> L2 normalize -> exact cosine search.

    IndexFlatIP is exact. For millions of windows on Sapelo2, swap in
    faiss.IndexHNSWFlat or IndexIVFFlat; the interface stays the same.
    """

    def __init__(self, pseudocount: float = 0.5):
        self.pseudocount = pseudocount
        self.mu = self.sd = None
        self.index = None
        self.meta = None

    def _embed(self, counts: np.ndarray) -> np.ndarray:
        x = clr(counts, self.pseudocount)
        x = (x - self.mu) / self.sd
        x = np.ascontiguousarray(x, dtype=np.float32)
        faiss.normalize_L2(x)
        return x

    def build(self, counts: np.ndarray, meta=None):
        x = clr(counts, self.pseudocount)
        self.mu = x.mean(axis=0)
        self.sd = x.std(axis=0) + 1e-6
        emb = self._embed(counts)
        self.index = faiss.IndexFlatIP(emb.shape[1])
        self.index.add(emb)
        self.meta = meta
        return self

    def search(self, counts: np.ndarray, k: int = 10):
        """Returns (similarities [n,k], reference row ids [n,k])."""
        return self.index.search(self._embed(counts), k)


def vote(sims: np.ndarray, ids: np.ndarray, Y_ref: np.ndarray,
         temperature: float = 0.05) -> np.ndarray:
    """Similarity-weighted KO vote. Weights are softmax(sim / T) over the k
    neighbors, so P(KO) = sum of weights of neighbors carrying that KO.
    Lower T trusts the single nearest neighbor more."""
    w = np.exp((sims - sims.max(axis=1, keepdims=True)) / temperature)
    w /= w.sum(axis=1, keepdims=True)
    return np.einsum("nk,nkj->nj", w, Y_ref[ids])


def transfer_genome_level(index: TNFIndex, Y_ref: np.ndarray, query_counts: np.ndarray,
                          k: int = 5, temperature: float = 0.05):
    sims, ids = index.search(query_counts, k)
    return vote(sims, ids, Y_ref, temperature), sims, ids


def transfer_window_level(index: TNFIndex, Y_ref_windows: np.ndarray,
                          query_window_counts: np.ndarray, k: int = 10,
                          temperature: float = 0.05, pool: str = "max"):
    """Traverse a query genome window by window. Each window gets a KO
    probability vector from its neighbors; pooling across windows gives the
    genome-level prediction ("max": a KO is present if any window supports it;
    "noisy_or": 1 - prod(1 - p) across windows)."""
    sims, ids = index.search(query_window_counts, k)
    P_win = vote(sims, ids, Y_ref_windows, temperature)
    if pool == "max":
        P = P_win.max(axis=0)
    elif pool == "noisy_or":
        P = 1.0 - np.exp(np.log1p(-np.clip(P_win, 0, 1 - 1e-6)).sum(axis=0))
    else:
        raise ValueError(pool)
    return P, P_win, sims, ids
