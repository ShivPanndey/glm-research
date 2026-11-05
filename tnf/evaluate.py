"""Multi-label metrics and leakage-safe splits."""
from __future__ import annotations

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform


def multilabel_metrics(Y_true: np.ndarray, P: np.ndarray, threshold: float = 0.5) -> dict:
    """Micro/macro precision, recall, F1. Macro averages only over KOs present
    in the ground truth (absent KOs would otherwise inflate or deflate it)."""
    Yp = (P >= threshold).astype(np.float32)
    tp = (Yp * Y_true).sum(0)
    fp = (Yp * (1 - Y_true)).sum(0)
    fn = ((1 - Yp) * Y_true).sum(0)
    micro_p = tp.sum() / max(tp.sum() + fp.sum(), 1e-9)
    micro_r = tp.sum() / max(tp.sum() + fn.sum(), 1e-9)
    micro_f1 = 2 * micro_p * micro_r / max(micro_p + micro_r, 1e-9)
    present = Y_true.sum(0) > 0
    p_j = tp / np.maximum(tp + fp, 1e-9)
    r_j = tp / np.maximum(tp + fn, 1e-9)
    f_j = 2 * p_j * r_j / np.maximum(p_j + r_j, 1e-9)
    return {
        "micro_precision": float(micro_p), "micro_recall": float(micro_r),
        "micro_f1": float(micro_f1),
        "macro_f1": float(f_j[present].mean()) if present.any() else float("nan"),
        "threshold": threshold,
    }


def best_threshold(Y_true, P, grid=np.linspace(0.05, 0.95, 19)) -> dict:
    """Pick the threshold on VALIDATION data only, then apply it to test."""
    scores = [multilabel_metrics(Y_true, P, t) for t in grid]
    return max(scores, key=lambda d: d["micro_f1"])


def cluster_genomes(dist: np.ndarray, threshold: float, method: str = "average") -> np.ndarray:
    """Average-linkage clusters from a genome-genome distance matrix.

    On the real data, use 1 - ANI/100 (skani or fastANI) with threshold 0.05
    (95% ANI ~ species boundary). TNF distance works as a quick fallback.
    """
    d = squareform(dist, checks=False)
    return fcluster(linkage(d, method=method), t=threshold, criterion="distance")


def grouped_split(groups: np.ndarray, frac=(0.7, 0.15, 0.15), seed: int = 0):
    """Assign whole clusters to train/val/test so close relatives never
    straddle the split. Clusters are placed largest first (ties shuffled) into
    whichever split is furthest below its target share. Returns row indices."""
    rng = np.random.default_rng(seed)
    uniq, sizes = np.unique(groups, return_counts=True)
    order = sorted(rng.permutation(len(uniq)), key=lambda i: -sizes[i])
    target = np.array(frac) * len(groups)
    filled = np.zeros(3)
    split_of = {}
    for i in order:
        s = int(np.argmax((target - filled) / target))
        split_of[uniq[i]] = s
        filled[s] += sizes[i]
    lab = np.array([split_of[g] for g in groups])
    return [np.nonzero(lab == i)[0] for i in range(3)]
