#!/usr/bin/env python
"""TNF retrieval baseline on the lab dataset.

Step 1 (parallel on Sapelo2, one SLURM array task per chunk of genomes):
  python scripts/run_baseline.py featurize --genomes DIR --out cache/ [--window 5000 --stride 2500]

Step 2 (single job):
  python scripts/run_baseline.py evaluate --cache cache/ --annotations DIR [--faa DIR]
         [--clusters clusters.tsv] --out results/

Genome <-> annotation files are matched on a normalized id: the file name up to
the first '.', with a trailing '_genome', '_proteins' or '_genes' removed.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tnf import (KOVocab, TNFIndex, best_threshold, clr, cluster_genomes, featurize_genome,  # noqa: E402
                 genome_ko_set, grouped_split, multilabel_metrics, parse_emapper,
                 parse_prodigal_faa, transfer_genome_level, transfer_window_level,
                 window_ko_sets)

FASTA_RE = re.compile(r"\.(fa|fna|fasta)(\.gz)?$")


def norm_id(p: Path) -> str:
    return re.sub(r"_(genome|proteins|genes)$", "", p.name.split(".")[0])


# ------------------------------------------------------------------ featurize
def cmd_featurize(a):
    files = sorted(p for p in Path(a.genomes).iterdir() if FASTA_RE.search(p.name))
    task = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
    n_tasks = int(os.environ.get("SLURM_ARRAY_TASK_COUNT", 1))
    mine = files[task::n_tasks]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for f in mine:
        gid = norm_id(f)
        dest = out / f"{gid}.npz"
        if dest.exists() and not a.overwrite:
            continue
        total, tab = featurize_genome(f, a.window, a.stride, genome_id=gid)
        if tab is None:
            np.savez_compressed(dest, total=total)
        else:
            np.savez_compressed(dest, total=total, contigs=tab.contigs.astype(str),
                                starts=tab.starts, ends=tab.ends,
                                counts=tab.counts.astype(np.uint16))
    print(f"task {task}/{n_tasks}: featurized {len(mine)} genomes -> {out}")


# ------------------------------------------------------------------- evaluate
def load_labels(ann_dir: Path):
    out = {}
    for p in sorted(ann_dir.glob("*.emapper.annotations*")):
        out[norm_id(p)] = parse_emapper(p)
    return out


def frequency_prior(Y_train: np.ndarray, n_query: int) -> np.ndarray:
    """Baseline: P(KO) = fraction of training genomes that carry it. Predicts
    the same 'typical genome' for every query. Retrieval must beat this."""
    return np.tile(Y_train.mean(0), (n_query, 1))


def per_genome_f1(Y, P, t):
    Yp = P >= t
    tp = (Yp & (Y > 0)).sum(1)
    prec = tp / np.maximum(Yp.sum(1), 1)
    rec = tp / np.maximum((Y > 0).sum(1), 1)
    return 2 * prec * rec / np.maximum(prec + rec, 1e-9)


def cmd_evaluate(a):
    cache = {p.stem: p for p in Path(a.cache).glob("*.npz")}
    labels = load_labels(Path(a.annotations))
    ids = sorted(set(cache) & set(labels))
    missing = sorted(set(cache) ^ set(labels))
    print(f"{len(ids)} genomes with both TNF and labels; {len(missing)} unmatched")
    if len(ids) < 10:
        sys.exit("need at least ~10 matched genomes for a train/val/test split")

    G = np.vstack([np.load(cache[g])["total"] for g in ids])
    ko_sets = [genome_ko_set(labels[g]) for g in ids]

    # clusters: ANI-based file if given, else TNF-distance fallback
    if a.clusters:
        cmap = dict(l.split()[:2] for l in open(a.clusters) if l.strip())
        groups = np.array([cmap.get(g, g) for g in ids])
    else:
        Z = clr(G)
        Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-6)
        Z /= np.linalg.norm(Z, axis=1, keepdims=True)
        groups = cluster_genomes(np.clip(1 - Z @ Z.T, 0, None), a.tnf_cluster_threshold)
    tr, va, te = grouped_split(groups, seed=a.seed)
    vocab = KOVocab([ko_sets[i] for i in tr])  # vocabulary from train only
    Y = vocab.to_matrix(ko_sets)
    print(f"split: train={len(tr)} val={len(va)} test={len(te)} genomes, "
          f"{len(np.unique(groups))} clusters, {len(vocab)} KOs")

    index = TNFIndex().build(G[tr])
    report = {"n_genomes": len(ids), "n_kos": len(vocab),
              "split_sizes": [len(tr), len(va), len(te)], "methods": {}}

    # ---- genome-level retrieval: tune k and temperature on validation
    best = None
    for k in (1, 3, 5, 10, 20):
        for T in (0.01, 0.05, 0.2):
            kk = min(k, len(tr))
            P, _, _ = transfer_genome_level(index, Y[tr], G[va], kk, T)
            m = best_threshold(Y[va], P)
            if best is None or m["micro_f1"] > best[0]["micro_f1"]:
                best = (m, kk, T)
    m_val, k, T = best
    P_te, sims, _ = transfer_genome_level(index, Y[tr], G[te], k, T)
    report["methods"]["tnf_genome_knn"] = {
        "k": k, "temperature": T, "val": m_val,
        "test": multilabel_metrics(Y[te], P_te, m_val["threshold"])}

    # ---- frequency prior baseline (threshold tuned on val)
    m_val_p = best_threshold(Y[va], frequency_prior(Y[tr], len(va)))
    P_prior = frequency_prior(Y[tr], len(te))
    report["methods"]["frequency_prior"] = {
        "val": m_val_p, "test": multilabel_metrics(Y[te], P_prior, m_val_p["threshold"])}

    # ---- OOD view: accuracy vs similarity to nearest training genome
    f1_g = per_genome_f1(Y[te], P_te, m_val["threshold"])
    report["ood_test_genomes"] = [{"genome": ids[i], "nearest_train_cosine": float(s), "f1": float(f)}
                                  for i, s, f in zip(te, sims[:, 0], f1_g)]

    # ---- window-level traversal (needs Prodigal coordinates)
    if a.faa:
        faa = {norm_id(p): p for p in Path(a.faa).iterdir() if ".faa" in p.name}
        win_ids = [i for i in range(len(ids)) if ids[i] in faa and "counts" in np.load(cache[ids[i]])]
        if win_ids:
            W = {}
            for i in win_ids:
                z = np.load(cache[ids[i]], allow_pickle=False)
                orfs = parse_prodigal_faa(faa[ids[i]])
                wk = window_ko_sets(z["contigs"], z["starts"], z["ends"], orfs, labels[ids[i]])
                W[i] = (z["counts"].astype(np.float64), vocab.to_matrix(wk))
            trw = [i for i in tr if i in W]
            Cw = np.vstack([W[i][0] for i in trw])
            Yw = np.vstack([W[i][1] for i in trw])
            widx = TNFIndex().build(Cw)
            preds = {}
            for split, members in (("val", va), ("test", te)):
                rows = [i for i in members if i in W]
                preds[split] = (rows, np.vstack([
                    transfer_window_level(widx, Yw, W[i][0], k=a.window_k, pool="max")[0]
                    for i in rows]))
            rv, Pv = preds["val"]
            rt, Pt = preds["test"]
            m_v = best_threshold(Y[rv], Pv)
            report["methods"]["tnf_window_traversal"] = {
                "k": a.window_k, "pool": "max", "val": m_v,
                "test": multilabel_metrics(Y[rt], Pt, m_v["threshold"])}

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "baseline_report.json").write_text(json.dumps(report, indent=2))
    for name, r in report["methods"].items():
        t = r["test"]
        print(f"{name:24s} test micro-F1={t['micro_f1']:.3f}  macro-F1={t['macro_f1']:.3f}  "
              f"P={t['micro_precision']:.3f} R={t['micro_recall']:.3f}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("featurize")
    f.add_argument("--genomes", required=True)
    f.add_argument("--out", required=True)
    f.add_argument("--window", type=int, default=None)
    f.add_argument("--stride", type=int, default=None)
    f.add_argument("--overwrite", action="store_true")
    e = sub.add_parser("evaluate")
    e.add_argument("--cache", required=True)
    e.add_argument("--annotations", required=True)
    e.add_argument("--faa", default=None)
    e.add_argument("--clusters", default=None, help="TSV: genome_id<TAB>cluster_id (e.g. from skani at 95% ANI)")
    e.add_argument("--tnf-cluster-threshold", type=float, default=0.1)
    e.add_argument("--window-k", type=int, default=10)
    e.add_argument("--seed", type=int, default=0)
    e.add_argument("--out", required=True)
    a = ap.parse_args()
    {"featurize": cmd_featurize, "evaluate": cmd_evaluate}[a.cmd](a)


if __name__ == "__main__":
    main()
