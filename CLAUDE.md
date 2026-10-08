# CLAUDE.md: GLM Research, TNF retrieval baseline

Context file for Claude Code. Read this fully before doing anything in this repo.

## Who I am and how to work with me

- Shiv (Shivank Pandey), undergraduate AI/ML student at UGA, first author on this project.
- I am an AI engineer with no prior genomics or microbiology background. When biology
  matters, explain it from first principles in ML terms (features, labels, leakage,
  distribution shift) before using jargon.
- House style for any prose you write for me (README, reports, comments, commit
  messages): **no em dashes**. Use commas, colons, semicolons, parentheses, or new sentences.
- Never report a number you did not actually compute in this repo. Mark public-data
  results as sanity checks, not project results.
- I want project progress to be based on the **lab's own data**, not public genomes.

## Project

Predict microbial genome function directly from genome sequence with an interpretable
genome language model (GLM), with calibrated uncertainty. Co-authors and advisors:
Dr. Srimanti Duttagupta and Dr. Avishek Dutta (UGA Geology; Dutta also Savannah River
Ecology Lab). Timeline: Sept 2026 to Jan 6, 2027. Emphasis on carbon, nitrogen, sulfur,
and metal-cycling functions.

### Scope (decided Sept 2026, do not drift)

- **"Function" means only the KEGG_ko field.** The model predicts KEGG_ko labels from
  genome `.fa` files.
- **Standard pipeline to beat:** genome -> Prodigal gene calls -> amino acid sequences ->
  eggNOG-mapper -> KEGG_ko.
- **In scope (Approach 1 only):** whole genome sequence -> model -> KEGG_ko, with no gene
  calling or homology search at prediction time. Components: a TNF retrieval baseline
  (this repo), a windowed per-region model, a genome-level multi-label model, and a
  hybrid (retrieved KOs fed to the GLM as prior or context).
- **Out of scope:** Approach 2 (genome -> 16S rRNA -> function). Do not bring up 16S,
  Barrnap, PICRUSt2, or PAPRICA.

### Core idea of this repo

Embed genomes (and genome windows) as tetranucleotide frequency (TNF) vectors, store them
in a vector database, retrieve the nearest annotated references for a query, and
transfer their KOs by similarity-weighted vote.

## Data

- **Location:** UGA SharePoint site `DUTTA-LAB-Ashutosh-project`, folder `Data_for_Shivank`,
  single file `family_1000_dataset.tar.gz` (1.35 GB compressed, uploaded 2026-09-22).
  **Contents not yet inspected.** Plan: move it to Sapelo2 scratch.
- **Known conventions from earlier samples:**
  - Genomes are multi-contig `.fa` files. Example: `2013515001_1_genome.fa`, a Yellowstone
    (YNP8) Crenarchaeota MAG, 81 contigs, 1.81 Mb, GC 32%.
  - Labels are eggNOG-mapper `.emapper.annotations` TSVs, 21 columns, **KEGG_ko in column
    12**, values like `ko:K00001,ko:K00002` or `-`. Produced with emapper 2.1.9 (diamond
    mode) from Prodigal `.faa` on a SLURM cluster (original paths under `/scratch/ak03777/`,
    directory `all_dereplicated_genomes`).
  - The annotation files name each ORF (for example `contig7_12`, the 12th gene on
    contig7) but carry **no genome coordinates**. Prodigal `.faa` headers or `.gff` files
    supply coordinates. Without them, window-level labels cannot be built.
- **Open data request (blocking window-level work):** Prodigal `.faa` or `.gff` files, if
  not already inside the tarball.

### First thing to do with the tarball

```bash
# Inventory and a small matched sample, without full extraction:
python scripts/peek_dataset.py family_1000_dataset.tar.gz --n 20
# -> dataset_summary.txt and dataset_sample.zip next to the archive
```

Then check `tnf_retrieval` file matching: `scripts/run_baseline.py::norm_id` matches
genome and annotation files on the name up to the first `.`, minus a trailing
`_genome`, `_proteins`, or `_genes`. **Adapt `norm_id` to the real naming scheme** and
confirm the "unmatched" count printed by `evaluate` is near zero.

## Compute

UGA GACRC **Sapelo2** (SLURM). Partition `batch`. Use `/scratch/$USER/` for data. A
template array job is in `scripts/featurize_array.sbatch` (module name and env path are
placeholders: check `module spider Python`).

## Repo layout

```
tnf/kmers.py          canonical 4-mer counting (136 dims), sliding windows, CLR, FASTA reader
tnf/labels.py         emapper KEGG_ko parser, Prodigal .faa/.gff coordinate parsers,
                      ORF -> window KO mapping, KOVocab (multi-label matrix)
tnf/retrieval.py      TNFIndex (CLR -> z-score -> L2 -> FAISS IndexFlatIP),
                      softmax similarity vote, genome-level and window-level transfer
tnf/evaluate.py       micro/macro P/R/F1, threshold tuning, genome clustering, grouped splits
scripts/run_baseline.py         `featurize` (SLURM-array aware) and `evaluate` subcommands
scripts/featurize_array.sbatch  Sapelo2 array template
scripts/peek_dataset.py         inventory + sample extraction from the lab tarball
experiments/window_taxonomy.py  public-genome sanity check, figs 1 to 4
experiments/locus_recovery.py   public-genome synteny test, fig 5
tests/test_tnf.py               unit tests (run: python -m pytest -q tests)
tests/make_smoke_data.py        builds a fake dataset in lab formats (synthetic KOs) for smoke tests
data/genomes/                   6 public genomes used only for sanity checks
results/                        figures and JSON from the experiments
```

## Design decisions (and why)

- **Canonical k-mers, k=4, 136 dims.** Each 4-mer merged with its reverse complement so
  the vector does not depend on which DNA strand was sequenced. 16 of 256 4-mers are
  palindromes, hence (256 - 16) / 2 + 16 = 136.
- **CLR transform** (pseudocount 0.5). TNF is compositional (only proportions matter).
- **z-score on reference statistics, then L2 norm, then inner product** = cosine search.
  `IndexFlatIP` is exact; switch to HNSW or IVF when the window index reaches millions.
- **Window counting uses a block trick:** k-mers binned into stride-sized blocks, windows
  summed from blocks. The last 3 k-mers can overhang the window end by up to 3 bp
  (tested and accepted). Requires `window % stride == 0`. Windows never cross contigs.
- **Window counts cached as uint16.** Safe for windows up to about 65 kb; widen the dtype
  if you go larger.
- **KO vote:** weights = softmax(similarity / T) over k neighbors; P(KO) = sum of weights of
  neighbors carrying it. Window mode pools per-window probabilities by max (or noisy-OR).
- **Leakage control:** whole clusters go to train, val, or test (`grouped_split`, largest
  clusters placed first into the split furthest below target). **Prefer ANI clusters**
  (for example skani at 95% ANI, passed via `--clusters`). The TNF-distance fallback is
  circular because the split then depends on the feature being evaluated.
- **Baselines:** the KO frequency prior (same "typical genome" for every query) must
  always be reported next to retrieval. k, T, and the decision threshold are tuned on
  validation only; test is touched once.
- **KO vocabulary is built from training genomes only.**

## Results so far (public genomes, sanity checks only)

Six public genomes: E. coli K-12 and EC590, H. influenzae, Ca. Portiera aleyrodidarum,
B. fragilis, C. diphtheriae (3 phyla, GC 26 to 54%). Each contig cut in half: first half
is the reference, second half the queries.

- **Species from held-out windows** (one genome per species, 10-NN, chance 0.20): balanced
  accuracy 0.89 at 0.5 kb, 0.95 at 1 kb, 0.98 at 2 kb, 0.99 at 5 kb, 1.00 at 10 to 20 kb.
- **Locus recovery:** 93% of E. coli K-12 5 kb windows retrieve the homologous locus in
  EC590 (median cosine 0.88); the dot plot shows the EC590 assembly is reverse
  complemented and rotated, plus a real inversion. Against H. influenzae only 12% (chance
  level, median cosine 0.64), and hits collapse onto a few **hub windows**.
- **Interpretation:** TNF retrieval should transfer function when a close relative is in
  the reference set and degrade toward a genome-average guess with phylogenetic
  distance. Key evaluation: accuracy versus similarity to the nearest training genome.
- **Caveat:** these species differ strongly in GC, so this is an easy case. No KO
  accuracy yet, because that needs the lab annotations.
- Gotcha found: a naive strain-level test scored below chance (0.11) because query windows
  matched the homologous locus in the sister strain. Any window-level split must account
  for homology across related genomes, not just overlap within one genome.

## Next steps (in order)

1. Inspect the tarball (`peek_dataset.py`), record its real structure in this file, adapt
   `norm_id`.
2. Run `run_baseline.py evaluate` on a sample of lab genomes: first real KEGG_ko numbers,
   TNF genome-kNN versus frequency prior.
3. Compute ANI clusters (skani or fastANI) for proper splits; rerun.
4. If `.faa`/`.gff` are available: window-level traversal; sweep window size (1, 2, 5,
   10 kb) and k; isotonic calibration of window probabilities.
5. OOD plot: per-genome F1 versus nearest-training-genome cosine (data already written to
   `baseline_report.json` under `ood_test_genomes`).
6. Hubness mitigation for window retrieval (CSLS or similar); compare.
7. Break results out by carbon, nitrogen, sulfur, and metal-cycling KO subsets.
8. Hybrid hook: retrieved KOs as prior or context for the GLM.

Original 15-week plan (from Sept 24): weeks 1 to 2 data and labels; 3 to 4 featurizer
(done); 5 to 6 genome-level baseline (code done, needs lab data); 7 to 9 window traversal
and calibration; week 10 Thanksgiving buffer; 11 to 12 functional-category and OOD
evaluation (light, finals); 13 to 15 hybrid prototype and results summary for the advisors.

## Commands

```bash
pip install -r requirements.txt
python -m pytest -q tests
python experiments/window_taxonomy.py && python experiments/locus_recovery.py

# smoke test of the full pipeline (synthetic labels, scores meaningless)
python tests/make_smoke_data.py /tmp/smoke
python scripts/run_baseline.py featurize --genomes /tmp/smoke/genomes --out /tmp/smoke/cache --window 5000 --stride 2500
python scripts/run_baseline.py evaluate --cache /tmp/smoke/cache --annotations /tmp/smoke/ann --faa /tmp/smoke/faa --out /tmp/smoke/results

# real data on Sapelo2
sbatch scripts/featurize_array.sbatch
python scripts/run_baseline.py evaluate --cache /scratch/$USER/glm/tnf_cache_w5k \
  --annotations $DATA/annotations --faa $DATA/prodigal --clusters clusters_ani95.tsv --out results/
```
