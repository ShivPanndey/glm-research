"""Does window-level TNF retrieval land on the homologous locus in a relative?

Query every 5 kb window of E. coli K-12 against all windows of E. coli EC590
(and, as a control, against H. influenzae). If retrieval tracks sequence
homology, best-hit positions form a synteny line in a dot plot.
"""
import sys
from pathlib import Path
import json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from tnf import TNFIndex, featurize_genome

def gp(s): return ROOT / "data/genomes" / f"{s}.fna.gz"
W, STRIDE = 5000, 2500
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

_, q = featurize_genome(gp("ecoli_K12"), W, STRIDE)
out = {}
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.4))
for ax, target, label in [(axes[0], "ecoli_EC590", "E. coli EC590 (same species)"),
                          (axes[1], "haemophilus_influenzae", "H. influenzae (different genus)")]:
    _, r = featurize_genome(gp(target), W, STRIDE)
    idx = TNFIndex().build(r.counts)
    sims, ids = idx.search(q.counts, 1)
    hit_pos = r.starts[ids[:, 0]]
    # synteny score: fraction of query windows whose best hit lies within 50 kb of
    # the hit predicted by a single best rotation/reflection fit
    L = r.starts.max() + W
    # fit the best offset separately for each orientation, then count a window as
    # syntenic if it sits within 50 kb of either line (handles inversions)
    syn = np.zeros(len(hit_pos), dtype=bool)
    for flip in (1, -1):
        d = (flip * hit_pos - q.starts) % L
        hist = np.bincount((d // 10000).astype(int))
        off = hist.argmax() * 10000 + 5000
        resid = np.minimum((d - off) % L, (off - d) % L)
        syn |= resid <= 50000
    best = float(syn.mean())
    out[target] = {"median_best_hit_cosine": float(np.median(sims[:, 0])), "syntenic_fraction": best}
    ax.scatter(q.starts / 1e6, hit_pos / 1e6, s=2, alpha=0.5, color="tab:blue", rasterized=True)
    ax.set_xlabel("E. coli K-12 query window position (Mb)")
    ax.set_ylabel(f"Best-hit position in {label.split(' (')[0]} (Mb)")
    ax.set_title(f"vs {label}\n{100*best:.0f}% of windows hit the syntenic locus", fontsize=10)
fig.suptitle("5 kb TNF retrieval: nearest window across genomes", y=1.0)
fig.tight_layout()
fig.savefig(ROOT / "results/fig5_locus_recovery_dotplot.png", dpi=200, bbox_inches="tight")
print(json.dumps(out, indent=1))
(ROOT / "results/locus_recovery.json").write_text(json.dumps(out, indent=1))
