"""Preliminary experiment on six public bacterial genomes.

Questions:
  Q1. Do whole-genome TNF vectors separate organisms? (sanity check)
  Q2. How large must a window be before its TNF vector reliably retrieves
      windows from the right organism? This sets the floor for window-level
      KO transfer: if a window can't find its relatives, it can't inherit
      their functions.
  Q3. Can TNF tell apart two strains of the same species (E. coli K-12 vs EC590)?

Leakage control: each contig is cut in half. Windows from the first half form
the reference database; windows from the second half are queries. Query and
reference windows never overlap.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.metrics import balanced_accuracy_score, confusion_matrix

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tnf import TNFIndex, clr, featurize_genome, read_fasta  # noqa: E402

GENOMES = {  # file stem: (display name, species, phylum)
    "ecoli_K12":              ("E. coli K-12",       "E. coli",        "Pseudomonadota"),
    "ecoli_EC590":            ("E. coli EC590",      "E. coli",        "Pseudomonadota"),
    "haemophilus_influenzae": ("H. influenzae",      "H. influenzae",  "Pseudomonadota"),
    "portiera_aleyrodidarum": ("Ca. Portiera",       "Ca. Portiera",   "Pseudomonadota"),
    "bacteroides_fragilis":   ("B. fragilis",        "B. fragilis",    "Bacteroidota"),
    "NCTC11397":              ("C. diphtheriae",     "C. diphtheriae", "Actinomycetota"),
}
WINDOWS = [500, 1000, 2000, 5000, 10000, 20000]
K = 10
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})


def gpath(stem):
    return ROOT / "data" / "genomes" / f"{stem}.fna.gz"


def split_halves(table, contig_len):
    half = np.array([contig_len[c] / 2 for c in table.contigs])
    ref = table.ends <= half
    qry = table.starts >= half
    return ref, qry


def knn_vote(sims, ids, ref_labels, n_classes, T=0.05):
    w = np.exp((sims - sims.max(1, keepdims=True)) / T)
    scores = np.zeros((len(sims), n_classes))
    for j in range(ids.shape[1]):
        np.add.at(scores, (np.arange(len(sims)), ref_labels[ids[:, j]]), w[:, j])
    return scores.argmax(1)


def main():
    stems = list(GENOMES)
    contig_len = {}
    for s in stems:
        for cid, seq in read_fasta(gpath(s)):
            contig_len[cid] = len(seq)

    # ---------------- Q1: whole-genome similarity -----------------------
    G = np.vstack([featurize_genome(gpath(s))[0] for s in stems])
    Z = clr(G)
    Zc = Z - Z.mean(0)
    Zn = Zc / np.linalg.norm(Zc, axis=1, keepdims=True)
    S = Zn @ Zn.T
    names = [GENOMES[s][0] for s in stems]
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    im = ax.imshow(S, cmap="viridis", vmin=-1, vmax=1)
    ax.set_xticks(range(len(stems)), names, rotation=40, ha="right")
    ax.set_yticks(range(len(stems)), names)
    for i in range(len(stems)):
        for j in range(len(stems)):
            ax.text(j, i, f"{S[i, j]:.2f}", ha="center", va="center", fontsize=8,
                    color="black" if S[i, j] > 0.3 else "white")
    fig.colorbar(im, ax=ax, label="cosine similarity (centered CLR TNF)")
    ax.set_title("Whole-genome TNF similarity")
    fig.tight_layout()
    fig.savefig(OUT / "fig1_genome_similarity.png", dpi=200)
    plt.close(fig)

    # ---------------- Q2/Q3: window retrieval vs window size -------------
    species_list = sorted({GENOMES[s][1] for s in stems})
    sp_idx = {sp: i for i, sp in enumerate(species_list)}
    results = []
    pca_payload = None
    for w in WINDOWS:
        stride = w // 2
        counts, genome_lab, ref_mask, qry_mask = [], [], [], []
        for gi, s in enumerate(stems):
            _, tab = featurize_genome(gpath(s), window=w, stride=stride)
            r, q = split_halves(tab, contig_len)
            counts.append(tab.counts)
            genome_lab.append(np.full(len(tab.starts), gi))
            ref_mask.append(r)
            qry_mask.append(q)
        C = np.vstack(counts)
        g = np.concatenate(genome_lab)
        r = np.concatenate(ref_mask)
        q = np.concatenate(qry_mask)
        sp = np.array([sp_idx[GENOMES[stems[i]][1]] for i in g])

        # Species test uses ONE genome per species (EC590 excluded), so a correct
        # hit can only come from shared composition, never from a homologous
        # locus in a sister strain (see locus_recovery.py for that effect).
        one = g != stems.index("ecoli_EC590")
        rs, qs = r & one, q & one
        idx = TNFIndex().build(C[rs])
        sims, ids = idx.search(C[qs], K)
        pred_sp = knn_vote(sims, ids, sp[rs], len(species_list))
        bal_acc = balanced_accuracy_score(sp[qs], pred_sp)
        per_species = {species_list[c]: float((pred_sp[sp[qs] == c] == c).mean())
                       for c in range(len(species_list))}

        # strain test: E. coli query windows, E. coli reference only, 2 classes
        ec = [stems.index("ecoli_K12"), stems.index("ecoli_EC590")]
        r_ec = r & np.isin(g, ec)
        q_ec = q & np.isin(g, ec)
        idx_ec = TNFIndex().build(C[r_ec])
        s2, i2 = idx_ec.search(C[q_ec], K)
        strain_ref = (g[r_ec] == ec[1]).astype(int)
        pred_strain = knn_vote(s2, i2, strain_ref, 2)
        strain_acc = balanced_accuracy_score((g[q_ec] == ec[1]).astype(int), pred_strain)

        row = {"window_bp": w, "stride_bp": stride, "n_reference": int(rs.sum()),
               "n_query": int(qs.sum()), "species_balanced_accuracy": float(bal_acc),
               "per_species_recall": per_species, "ecoli_strain_balanced_accuracy": float(strain_acc),
               "chance_species": 1 / len(species_list)}
        results.append(row)
        print(f"w={w:>6} bp  species bal.acc={bal_acc:.3f}  strain bal.acc={strain_acc:.3f}  "
              f"(ref={r.sum()}, query={q.sum()})")
        if w == 5000:
            cm = confusion_matrix(sp[qs], pred_sp, normalize="true")
            pca_payload = (C, g)
            cm_payload = cm

    (OUT / "window_taxonomy_results.json").write_text(json.dumps(results, indent=2))

    # Fig 2: accuracy vs window size
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    ws = [r_["window_bp"] / 1000 for r_ in results]
    ax.plot(ws, [r_["species_balanced_accuracy"] for r_ in results], "o-", lw=2,
            label="Species (5 classes)")
    for x, r_ in zip(ws, results):
        ax.annotate(f"{r_['species_balanced_accuracy']:.2f}", (x, r_["species_balanced_accuracy"]),
                    textcoords="offset points", xytext=(0, -14), ha="center", fontsize=8)
    ax.axhline(0.2, color="gray", lw=1, ls=":")
    ax.text(ws[-1], 0.215, "chance (0.20)", ha="right", fontsize=8, color="gray")
    ax.set_xscale("log")
    ax.set_xticks(ws, [f"{x:g}" for x in ws])
    ax.set_xlabel("Window size (kb)")
    ax.set_ylabel("Balanced accuracy (held-out windows)")
    ax.set_ylim(0, 1.02)
    ax.set_title(f"Held-out windows: species identified by TNF retrieval ({K}-NN)")
    ax.legend(loc="lower right", frameon=False)
    ax.text(ws[0], 0.05, "one genome per species; reference = first half of each contig,\nqueries = second half",
            fontsize=7.5, color="dimgray")
    fig.tight_layout()
    fig.savefig(OUT / "fig2_accuracy_vs_window.png", dpi=200)
    plt.close(fig)

    # Fig 3: PCA of 5 kb windows
    C, g = pca_payload
    Z = clr(C)
    Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-6)
    pc = PCA(n_components=2, random_state=0).fit(Z)
    X = pc.transform(Z)
    fig, ax = plt.subplots(figsize=(5.8, 4.6))
    rng = np.random.default_rng(0)
    order = rng.permutation(len(X))
    colors = plt.cm.tab10(np.arange(len(stems)))
    for gi, s in enumerate(stems):
        m = g[order] == gi
        ax.scatter(X[order][m, 0], X[order][m, 1], s=3, alpha=0.35, color=colors[gi],
                   label=GENOMES[s][0], rasterized=True)
    ax.set_xlabel(f"PC1 ({100 * pc.explained_variance_ratio_[0]:.0f}% var.)")
    ax.set_ylabel(f"PC2 ({100 * pc.explained_variance_ratio_[1]:.0f}% var.)")
    ax.set_title("5 kb windows embedded by TNF (PCA)")
    ax.legend(markerscale=4, frameon=False, fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(OUT / "fig3_pca_5kb_windows.png", dpi=200)
    plt.close(fig)

    # Fig 4: confusion at 5 kb
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    im = ax.imshow(cm_payload, cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(len(species_list)), species_list, rotation=40, ha="right")
    ax.set_yticks(range(len(species_list)), species_list)
    for i in range(len(species_list)):
        for j in range(len(species_list)):
            ax.text(j, i, f"{cm_payload[i, j]:.2f}", ha="center", va="center", fontsize=8,
                    color="white" if cm_payload[i, j] > 0.5 else "black")
    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    ax.set_title("Species confusion, 5 kb windows")
    fig.colorbar(im, ax=ax, label="row-normalized")
    fig.tight_layout()
    fig.savefig(OUT / "fig4_confusion_5kb.png", dpi=200)
    plt.close(fig)
    print(json.dumps({"genome_similarity": dict(zip(names, S.round(3).tolist()))}, indent=1))


if __name__ == "__main__":
    main()
