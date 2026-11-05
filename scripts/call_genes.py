#!/usr/bin/env python
"""Re-create the Prodigal .faa files the lab's eggNOG-mapper run started from.

The lab tarball ships genomes and emapper outputs but no Prodigal files, and the
annotations name genes (contig_N = N-th gene on a contig) without coordinates.
Prodigal in single-genome mode (its default) reproduces those ids, so its .faa
headers give each annotated gene a position, which window-level labels need.

  python scripts/call_genes.py --genomes DIR --out faa/ [--annotations DIR]

SLURM-array aware like `run_baseline.py featurize`. With --annotations, every
genome is checked: each annotated gene id must exist in the new .faa, and a
genome that fails is listed in failures.task<ID>.tsv and its proteins are kept as <id>.faa.failed
for inspection, so reruns retry it and later steps never mistake it for a good .faa.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tnf import read_fasta  # noqa: E402
from jobutil import FailureLog, my_shard, task_count, task_id  # noqa: E402
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
    mine = my_shard(files)
    ann = {norm_id(p): p for p in Path(a.annotations).rglob("*.emapper.annotations")} if a.annotations else {}
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    with FailureLog(out) as log:
        for f in mine:
            gid = norm_id(f)
            dest = out / f"{gid}.faa"
            if dest.exists() and not a.overwrite:
                continue
            try:
                ids = write_genes(f, dest)
                if gid in ann:
                    want = annotated_ids(ann[gid])
                    missing = sum(q not in ids for q in want)
                    if missing:
                        # keep it for inspection, but not under the name later steps and
                        # reruns treat as done
                        dest.replace(dest.with_name(dest.name + ".failed"))
                        log.record(gid, f"ids: {missing}/{len(want)} not reproduced")
            except Exception as e:  # one bad genome must not stop the others
                log.record(gid, f"{type(e).__name__}: {e}")
    print(f"task {task_id()}/{task_count()}: called genes for {len(mine)} genomes -> {out}"
          + (f"; {log.count} failed, see {log.path.name}" if log.count else ""))
    if log.count:
        sys.exit(1)


def write_genes(fasta: Path, dest: Path) -> set[str]:
    """Train Prodigal on the whole genome (single mode), write the .faa, return gene ids."""
    contigs = [(name.split()[0], seq) for name, seq in read_fasta(fasta)]
    if not contigs:
        raise ValueError("no sequence")
    gf = pyrodigal.GeneFinder(meta=False)
    gf.train(*[s.encode() for _, s in contigs])
    ids = set()
    with open(dest, "w") as fh:
        for name, seq in contigs:
            genes = gf.find_genes(seq.encode())
            genes.write_translations(fh, sequence_id=name)
            ids.update(f"{name}_{i}" for i in range(1, len(genes) + 1))
    return ids

if __name__ == "__main__":
    main()
