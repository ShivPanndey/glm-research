# TNF retrieval baseline for KEGG_ko prediction

Predicts KEGG_ko functional labels for a prokaryotic genome by retrieving
compositionally similar, already-annotated genomes (or genome windows) from a
vector database of tetranucleotide frequency (TNF) vectors. This is the
retrieval baseline and retrieval component of the GLM project: it bypasses gene
calling and homology search at prediction time.

## How it works

1. **Featurize.** Count canonical 4-mers (136 dimensions; each 4-mer merged with
   its reverse complement) over the whole genome and over sliding windows that
   never cross contig boundaries.
2. **Embed.** Centered log-ratio transform (TNF is compositional), z-score using
   reference statistics, L2 normalize.
3. **Index.** Exact cosine search with FAISS (`IndexFlatIP`; swap in HNSW or IVF
   for millions of windows).
4. **Transfer.** Each KO gets a similarity-weighted vote from the k nearest
   references. Window mode votes per window, then pools across the genome.
5. **Evaluate.** Genome-clustered train/val/test splits; k, temperature and
   decision threshold tuned on validation only; compared with a KO-frequency
   prior that predicts the same "typical genome" for every query.

## Layout

```
tnf/kmers.py        canonical 4-mer counting, windows, CLR
tnf/labels.py       eggNOG-mapper KEGG_ko parser, Prodigal .faa/.gff coordinates, window -> KO
tnf/retrieval.py    FAISS index, genome-level and window-level label transfer
tnf/evaluate.py     multi-label metrics, genome clustering, leakage-safe splits
scripts/run_baseline.py        featurize (SLURM-array aware) and evaluate
scripts/sapelo2/               Sapelo2 workflow: config, setup, preflight, jobs, RUNBOOK.md
experiments/window_taxonomy.py preliminary experiment, figures 1 to 4
experiments/locus_recovery.py  preliminary experiment, figure 5
tests/                         unit tests and smoke-test data generator
```

## Running on the lab dataset

On UGA's Sapelo2 cluster, follow `scripts/sapelo2/RUNBOOK.md`. In short:

```bash
bash scripts/sapelo2/setup_env.sh                       # once: venv, requirements, skani
source scripts/sapelo2/config.sh && python scripts/sapelo2/preflight.py
bash scripts/sapelo2/submit_all.sh --pilot              # 20 genomes; then without --pilot
bash scripts/sapelo2/collect_results.sh /scratch/$USER/glm/pilot
```

`submit_all.sh` chains gene re-calling, the qend check, featurization, skani ANI clusters,
and `evaluate --seeds 0,1,2`. Run `make check` before committing.

Prefer ANI-based clusters (for example skani at 95% ANI) over the TNF fallback,
because TNF clustering would make the split partly depend on the feature being
tested.

## Preliminary results (6 public genomes, before lab data)

Genomes: *E. coli* K-12 and EC590, *H. influenzae*, *Ca.* Portiera
aleyrodidarum, *B. fragilis*, *C. diphtheriae* (3 phyla, GC 26 to 54%).
Reproduce with `python experiments/window_taxonomy.py` and
`python experiments/locus_recovery.py`.

| Window | Species balanced accuracy (5 classes, chance 0.20) |
|---|---|
| 0.5 kb | 0.89 |
| 1 kb | 0.95 |
| 2 kb | 0.98 |
| 5 kb | 0.99 |
| 10 kb | 1.00 |

Across strains, 93% of *E. coli* K-12 5 kb windows retrieve the homologous
locus in EC590 (median cosine 0.88). Across genera (*H. influenzae*), 12% do,
which is chance level for that test, and hits collapse onto a few hub windows.

Caveat: these five species differ strongly in GC content, so species
identification here is an easy case. The lab's environmental MAGs will contain
many related lineages and will be harder.
