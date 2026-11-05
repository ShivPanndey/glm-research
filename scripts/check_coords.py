#!/usr/bin/env python
"""Check that re-called proteins match the ones the lab's eggNOG-mapper run aligned.

The lab tarball has no Prodigal files, so `call_genes.py` re-creates them. Gene ids already
match; this checks lengths too. For every annotated gene, the alignment end on the query
(qend, in amino acids, column 6 of <id>.emapper.seed_orthologs) must fit inside the
re-called protein. A gene longer in the original run would show qend > our length.

  python scripts/check_coords.py --faa prodigal/ --annotations eggnog/ --out coords/

Writes coords_check.tsv (one row per genome) and coords_check.json (summary). Report only:
always exits 0, so it never blocks genome-level evaluation.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_baseline import norm_id  # noqa: E402

SEED_SUFFIX = re.compile(r"\.emapper\.seed_orthologs(\.gz)?$")
QEND_COL = 5  # emapper 2.1.9: qseqid sseqid evalue bitscore qstart qend ...


def protein_lengths(faa: Path) -> dict[str, int]:
    lengths, name, n = {}, None, 0
    for line in open(faa):
        if line.startswith(">"):
            if name is not None:
                lengths[name] = n
            name, n = line[1:].split()[0], 0
        else:
            n += len(line.strip().rstrip("*"))
    if name is not None:
        lengths[name] = n
    return lengths


def seed_qends(path: Path) -> dict[str, int]:
    qends: dict[str, int] = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        cols = line.rstrip("\n").split("\t")
        q, qend = cols[0], int(cols[QEND_COL])
        qends[q] = max(qend, qends.get(q, 0))
    return qends


def check_genome(lengths: dict[str, int], qends: dict[str, int]) -> dict:
    missing = [q for q in qends if q not in lengths]
    too_long = [q for q in qends if q in lengths and qends[q] > lengths[q]]
    return {"n_genes": len(qends), "n_ok": len(qends) - len(missing) - len(too_long),
            "n_missing": len(missing), "n_too_long": len(too_long),
            "examples": (too_long + missing)[:5]}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--faa", required=True)
    ap.add_argument("--annotations", required=True, help="searched recursively for *.emapper.seed_orthologs")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    faa = {norm_id(p): p for p in Path(a.faa).rglob("*.faa")}
    seeds = {SEED_SUFFIX.sub("", p.name): p for p in Path(a.annotations).rglob("*.emapper.seed_orthologs*")}
    rows, missing_faa = {}, []
    for gid in sorted(seeds):
        if gid not in faa:
            missing_faa.append(gid)
            continue
        rows[gid] = check_genome(protein_lengths(faa[gid]), seed_qends(seeds[gid]))

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "coords_check.tsv", "w") as fh:
        fh.write("genome\tn_genes\tn_ok\tn_missing\tn_too_long\n")
        for gid, r in rows.items():
            fh.write(f"{gid}\t{r['n_genes']}\t{r['n_ok']}\t{r['n_missing']}\t{r['n_too_long']}\n")
    n_genes = sum(r["n_genes"] for r in rows.values())
    summary = {
        "n_genomes": len(rows),
        "n_genomes_all_ok": sum(r["n_ok"] == r["n_genes"] for r in rows.values()),
        "gene_pass_rate": sum(r["n_ok"] for r in rows.values()) / n_genes if n_genes else None,
        "n_missing_faa": len(missing_faa),
        "missing_faa": missing_faa[:50],
        "per_genome_failures": [{"genome": g, **r} for g, r in rows.items() if r["n_ok"] < r["n_genes"]][:50],
    }
    (out / "coords_check.json").write_text(json.dumps(summary, indent=2))
    print(f"coords: {summary['n_genomes_all_ok']}/{summary['n_genomes']} genomes fully consistent, "
          f"gene pass rate {summary['gene_pass_rate']}, {len(missing_faa)} genomes without a .faa")
    return 0


if __name__ == "__main__":
    sys.exit(main())
