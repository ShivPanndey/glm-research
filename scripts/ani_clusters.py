#!/usr/bin/env python
"""Turn skani all-vs-all output into 95% ANI clusters for leakage-safe splits.

ANI (average nucleotide identity) is how alike two genomes are along the parts they share.
At about 95% they are the same species: near-duplicate data points that must stay on one side
of a train/test split. skani only reports pairs above roughly 80% ANI; missing pairs count as 0.

  skani triangle -E -l genomes.txt -o skani_edges.tsv
  python scripts/ani_clusters.py --skani skani_edges.tsv --genomes genomes/ --out ani/

Writes clusters_ani95.tsv (id<TAB>cluster, the format `evaluate --clusters` reads) and
ani_summary.json (cluster sizes and each genome's best ANI to any other genome; if most are
below 80%, the set is close to one genome per family).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_baseline import FASTA_RE, norm_id  # noqa: E402
from tnf import cluster_genomes  # noqa: E402

CUT = 0.05  # distance 1 - ANI/100, i.e. 95% ANI
NEEDED = ("Ref_file", "Query_file", "ANI")


def read_skani_sparse(path: Path) -> list[tuple[str, str, float]]:
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for col in NEEDED:
            if col not in (reader.fieldnames or []):
                raise ValueError(f"skani output {path} has no '{col}' column; got {reader.fieldnames}. "
                                 "Was it made with `skani triangle -E`?")
        return [(norm_id(Path(r["Ref_file"])), norm_id(Path(r["Query_file"])), float(r["ANI"]))
                for r in reader]


def ani_distance(ids: list[str], pairs) -> np.ndarray:
    pos = {g: i for i, g in enumerate(ids)}
    D = np.ones((len(ids), len(ids)))
    np.fill_diagonal(D, 0.0)
    for a, b, ani in pairs:
        if a in pos and b in pos and a != b:
            d = 1 - ani / 100
            D[pos[a], pos[b]] = D[pos[b], pos[a]] = min(d, D[pos[a], pos[b]])
    return D


def nearest_ani(ids: list[str], pairs) -> dict[str, float]:
    best = dict.fromkeys(ids, 0.0)
    for a, b, ani in pairs:
        for g in (a, b):
            if g in best and a != b:
                best[g] = max(best[g], ani)
    return best


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skani", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    ids = sorted(norm_id(p) for p in Path(a.genomes).iterdir() if FASTA_RE.search(p.name))
    pairs = read_skani_sparse(Path(a.skani))
    labels = cluster_genomes(ani_distance(ids, pairs), CUT) if len(ids) > 1 else np.ones(len(ids), int)
    near = nearest_ani(ids, pairs)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "clusters_ani95.tsv").write_text("".join(f"{g}\t{c}\n" for g, c in zip(ids, labels)))
    sizes = Counter(Counter(labels.tolist()).values())
    vals = np.array(list(near.values()))
    summary = {
        "n_genomes": len(ids),
        "n_clusters": len(set(labels.tolist())),
        "cluster_size_counts": {str(k): v for k, v in sorted(sizes.items())},
        "nearest_ani_quantiles": {str(q): float(np.quantile(vals, q)) for q in (0, .1, .25, .5, .75, .9, 1)},
        "n_nearest_below_80": int((vals < 80).sum()),
        "n_pairs_reported": len(pairs),
    }
    (out / "ani_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"ani: {summary['n_genomes']} genomes, {summary['n_clusters']} clusters at 95% ANI, "
          f"{summary['n_nearest_below_80']} genomes with no relative above 80% ANI")
    return 0


if __name__ == "__main__":
    sys.exit(main())
