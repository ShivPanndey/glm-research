# Design decisions (and why)

Moved from CLAUDE.md on 2026-11-05. Add new decisions here with a date.


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
