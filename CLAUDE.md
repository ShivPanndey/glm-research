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
- **I am a visual, diagram-first person.** When planning, explaining, or asking me to approve
  a design, lead with diagrams (flow charts, maps of where code and data live, before/after
  pictures), then short text. Keep a running visual plan page updated as we discuss.

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
  single file `family_1000_dataset.tar.gz` (1,351,413,479 bytes, uploaded 2026-09-22). Some
  browsers (Safari) auto-decompress it to a plain `.tar`; that is the same data.
  Plan: move it to Sapelo2 scratch.

### Dataset structure (inspected 2026-10-08 with `peek_dataset.py`)

Re-inspect without full extraction (.tar or .tar.gz): `python scripts/peek_dataset.py
family_1000_dataset.tar --n 20` writes `dataset_summary.txt` and `dataset_sample.zip` next
to the archive (both gitignored; lab data stays off GitHub).

```
family_1000_dataset_final_v2/            4,013 files, 4.60 GB uncompressed
  genomes/<id>.fa                        1,000 genomes, flat folder, 422 KB to 12.5 MB (median 2.55 MB)
  eggnog/<id>/<id>.emapper.annotations   labels (KEGG_ko is column 12)
  eggnog/<id>/<id>.emapper.seed_orthologs, .emapper.hits   intermediate emapper outputs
  logs/                                  12 assembly logs (see below)
```

- All 1,000 genome ids have a matching annotation; 0 unmatched.
- Id families seen: `GCMeta_NNNNNNNN`, `spire_mag_NNNNNNNN`, `specI_v4_NNNNN`, IMG-style
  `3300017444_117`, and NCBI assemblies such as `GCF_003590795.1_ASM359079v1_genomic` (dots
  inside the id, so `norm_id` strips known suffixes instead of splitting at the first `.`).
- Contig headers vary by source: `GCMeta_01254050_00001`, `k141_126446`,
  `3300017444.a:Ga0185300_10004342`, `NODE_26_length_263618_cov_24.617230`.
- Annotation query ids are Prodigal style, `<contig>_<n>` = n-th gene on that contig.
  Only genes with an emapper hit get a row. Produced by emapper 2.1.9, diamond mode, from
  `chunk_*/<id>.faa` on Sapelo2 (`/scratch/ak03777/objective2/analysis/eggnog/...`).
- **No Prodigal `.faa` or `.gff` in the tarball.** Workaround: `scripts/call_genes.py`
  re-calls genes with pyrodigal in single-genome mode (Prodigal's default). On the 20-genome
  sample this reproduced 100% of annotated gene ids in every genome, with no annotated index
  beyond the predicted gene count; metagenome mode did not (95 to 99.8%). Id agreement does
  not yet prove identical coordinates: still to check that `.seed_orthologs` alignment ends
  (qend) fit inside the re-called protein lengths.
- **Logs to read before evaluating:** `selected_1000_names.txt`, `retained_995_names.txt`,
  `replacement_5_names.txt`, `final_genome_names.txt`, `copy_failures.tsv`, plus
  metadata/eggnog/fasta id lists. Suggests 5 of the original 1,000 picks were replaced. A
  metadata table (taxonomy?) was used but is not in the tarball.
- **Open question for the advisors:** does "family_1000" mean one genome per taxonomic
  family? If so, no test genome has a close relative in the reference set, which is the
  hardest case for retrieval (see results below), and the evaluation should say so.

### Matching files

`scripts/run_baseline.py::norm_id` matches files on the name minus its known suffix
(`.fa`, `.faa`, `.emapper.annotations`, ...) and a trailing `_genome`, `_proteins`, or
`_genes`. `--annotations` and `--faa` are searched recursively, so point `--annotations`
at `eggnog/`. Checked on the 20-genome sample: "0 unmatched".

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
scripts/call_genes.py           re-creates Prodigal .faa (pyrodigal, single mode), checks ids vs emapper
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

1. Done: tarball inspected, structure recorded above, `norm_id` adapted, `call_genes.py`
   added. Pipeline runs end to end on the 20-genome sample (too few genomes for numbers).
   Note: the sample is the 20 smallest genomes (0.42 to 0.64 MB), not representative.
2. Read the logs; confirm with the advisors what "family_1000" means. Verify re-called
   coordinates against `.seed_orthologs` qend.
3. On Sapelo2: `call_genes.py` and `featurize` as array jobs over all 1,000 genomes, then
   `evaluate`: first real KEGG_ko numbers, TNF genome-kNN versus frequency prior.
4. Compute ANI clusters (skani or fastANI) for proper splits; rerun.
5. Window-level traversal; sweep window size (1, 2, 5,
   10 kb) and k; isotonic calibration of window probabilities.
6. OOD plot: per-genome F1 versus nearest-training-genome cosine (data already written to
   `baseline_report.json` under `ood_test_genomes`).
7. Hubness mitigation for window retrieval (CSLS or similar); compare.
8. Break results out by carbon, nitrogen, sulfur, and metal-cycling KO subsets.
9. Hybrid hook: retrieved KOs as prior or context for the GLM.

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

# real data on Sapelo2 ($DATA = family_1000_dataset_final_v2)
python scripts/call_genes.py --genomes $DATA/genomes --annotations $DATA/eggnog --out $DATA/prodigal
sbatch scripts/featurize_array.sbatch
python scripts/run_baseline.py evaluate --cache /scratch/$USER/glm/tnf_cache_w5k \
  --annotations $DATA/eggnog --faa $DATA/prodigal --clusters clusters_ani95.tsv --out results/
```
