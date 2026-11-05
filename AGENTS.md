# AGENTS.md

Map of this repo for coding agents. The full context file is `CLAUDE.md`; read it first.

## What this repo is

TNF retrieval baseline for predicting KEGG_ko labels from whole genome sequence. Genomes become
136-dim tetranucleotide vectors, a FAISS index finds the nearest annotated genomes, and their KOs
are transferred by a similarity-weighted vote.

## Where things are

```
tnf/                  library: kmers, labels, retrieval, evaluate
scripts/              CLIs: run_baseline (featurize, evaluate), call_genes, check_coords,
                      ani_clusters, peek_dataset, check_style, jobutil
scripts/sapelo2/      cluster workflow: config, setup, preflight, sbatch jobs, RUNBOOK.md
experiments/          public-genome sanity checks (figures 1 to 5)
tests/                pytest suite, smoke data builder, smoke runner
docs/decisions.md     design decisions and why
docs/results-log.md   every reported number, with commit and command
docs/harness-log.md   what tripped us up and what changed in the harness
docs/superpowers/     specs and implementation plans
```

## Rules

1. Run `make check` before every commit and push. It runs ruff, pytest, the smoke pipeline, the
   em-dash check, and `bash -n` on shell scripts.
2. No em dashes in any text: use commas, colons, semicolons, parentheses, or new sentences.
   In code that must produce one, write `chr(0x2014)`.
3. Lab data never goes to git: no `.fa`, `.faa`, `.emapper.*`, `.tar*`, `.npz` outside
   `data/genomes/` and `tests/`. The commit hook enforces this.
4. Never report a number that was not computed in this repo. Log every number in
   `docs/results-log.md`, copied from a committed JSON.
5. Leakage rules: KO vocabulary from training genomes only; k, temperature, and threshold tuned
   on validation only; test touched once per seed; always report the frequency prior beside
   retrieval; prefer ANI clusters for splits.
6. Explain biology in ML terms for the project owner, and lead with diagrams when planning.
