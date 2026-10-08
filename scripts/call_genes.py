#!/usr/bin/env python
"""Re-create the Prodigal .faa files the lab's eggNOG-mapper run started from.

The lab tarball ships genomes and emapper outputs but no Prodigal files, and the
annotations name genes (contig_N = N-th gene on a contig) without coordinates.
Prodigal in single-genome mode (its default) reproduces those ids, so its .faa
headers give each annotated gene a position, which window-level labels need.

  python scripts/call_genes.py --genomes DIR --out faa/ [--annotations DIR]

SLURM-array aware like `run_baseline.py featurize`. With --annotations, every
genome is checked: each annotated gene id must exist in the new .faa, and a
genome that fails is reported (its .faa is still written, for inspection).
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tnf import read_fasta  # noqa: E402
from run_baseline import FASTA_RE, norm_id  # noqa: E402


def annotated_ids(path: Path) -> list[str]:
    with open(path) as fh:
        return [line.split("\t", 1)[0] for line in fh if not line.startswith("#")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--annotations", default=None, help="checked recursively for <id>.emapper.annotations")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()

    files = sorted(p for p in Path(a.genomes).iterdir() if FASTA_RE.search(p.name))
    task = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
    n_tasks = int(os.environ.get("SLURM_ARRAY_TASK_COUNT", 1))
    mine = files[task::n_tasks]
    ann = {norm_id(p): p for p in Path(a.annotations).rglob("*.emapper.annotations")} if a.annotations else {}
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    failed = []
    for f in mine:
        gid = norm_id(f)
        dest = out / f"{gid}.faa"
        if dest.exists() and not a.overwrite:
            continue
        contigs = [(name.split()[0], seq) for name, seq in read_fasta(f)]
        gf = pyrodigal.GeneFinder(meta=False)
        gf.train(*[s.encode() for _, s in contigs])
        ids = set()
        with open(dest, "w") as fh:
            for name, seq in contigs:
                genes = gf.find_genes(seq.encode())
                genes.write_translations(fh, sequence_id=name)
                ids.update(f"{name}_{i}" for i in range(1, len(genes) + 1))
        if gid in ann:
            want = annotated_ids(ann[gid])
            missing = sum(q not in ids for q in want)
            if missing:
                failed.append(gid)
                print(f"{gid}: {missing}/{len(want)} annotated gene ids not reproduced")
    print(f"task {task}/{n_tasks}: called genes for {len(mine)} genomes -> {out}"
          + (f"; {len(failed)} failed the id check" if ann else ""))
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
