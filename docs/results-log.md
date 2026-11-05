# Results log

Append-only. Every number reported anywhere (README, CLAUDE.md, reports, emails to advisors)
must appear here first, copied from a committed JSON file, never retyped from memory.

Entry format:

```
## YYYY-MM-DD: <title>
- commit: <git commit that produced the JSON>
- command: <exact command>
- source: <path to committed JSON>
- split: <how train and test were separated>
- kind: project result | sanity check (public genomes)
<metrics table>
```

## 2026-10-08: Window taxonomy and locus recovery on six public genomes

- commit: 9f01bff
- command: `python experiments/window_taxonomy.py && python experiments/locus_recovery.py`
- source: `results/window_taxonomy_results.json`, `results/locus_recovery.json`
- split: each contig cut in half; first half reference, second half queries
- kind: sanity check (public genomes), not a project result

Species from held-out windows (10-NN, chance 0.20):

| window | queries | species balanced accuracy | E. coli strain balanced accuracy |
|---|---|---|---|
| 0.5 kb | 29,202 | 0.894 | 0.111 |
| 1 kb | 14,595 | 0.949 | 0.112 |
| 2 kb | 7,291 | 0.976 | 0.105 |
| 5 kb | 2,910 | 0.990 | 0.104 |
| 10 kb | 1,449 | 0.997 | 0.106 |
| 20 kb | 719 | 1.000 | 0.102 |

Locus recovery of E. coli K-12 5 kb windows:

| target genome | median best-hit cosine | syntenic fraction |
|---|---|---|
| E. coli EC590 | 0.880 | 0.932 |
| H. influenzae | 0.639 | 0.125 |
