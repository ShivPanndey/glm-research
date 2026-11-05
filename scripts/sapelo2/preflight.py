#!/usr/bin/env python
"""Check the Sapelo2 setup before submitting any job. Runs in seconds on a login node.

  source scripts/sapelo2/config.sh && python scripts/sapelo2/preflight.py [--expect 1000]

Reads DATA, WORK_OUT, SCRATCH from the environment (set by config.sh). Prints one line per
check and exits 1 if any check fails. Paste the output to Claude after the first run.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_baseline import FASTA_RE, norm_id  # noqa: E402

Check = tuple[str, bool, str]


def writable(d: Path) -> tuple[bool, str]:
    try:
        d.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=d):
            pass
        return True, str(d)
    except OSError as e:
        return False, f"{d}: {e}"


def seed_format(eggnog: Path) -> tuple[bool, str]:
    for p in eggnog.rglob("*.emapper.seed_orthologs"):
        for line in open(p):
            if line.strip() and not line.startswith("#"):
                cols = line.rstrip("\n").split("\t")
                return len(cols) >= 6, f"{p.name}: {len(cols)} columns, qend (col 6) = {cols[5:6]}"
        return False, f"{p.name}: no data rows"
    return False, "no .emapper.seed_orthologs files found"


def run_checks(data: Path, work_out: Path, scratch: Path, expect: int = 1000) -> list[Check]:
    genomes = data / "genomes"
    eggnog = data / "eggnog"
    g_ids = {norm_id(p) for p in genomes.iterdir() if FASTA_RE.search(p.name)} if genomes.is_dir() else set()
    a_ids = {norm_id(p) for p in eggnog.rglob("*.emapper.annotations")} if eggnog.is_dir() else set()
    checks: list[Check] = [
        ("data folder", data.is_dir(), str(data)),
        ("genome count", len(g_ids) == expect, f"{len(g_ids)} (expected {expect})"),
        ("annotation count", len(a_ids) == expect, f"{len(a_ids)} (expected {expect})"),
        ("genomes without annotations", not (g_ids - a_ids), f"{len(g_ids - a_ids)}: {sorted(g_ids - a_ids)[:5]}"),
        ("annotations without genomes", not (a_ids - g_ids), f"{len(a_ids - g_ids)}: {sorted(a_ids - g_ids)[:5]}"),
        ("logs folder", (data / "logs").is_dir(), str(data / "logs")),
        ("seed_orthologs format", *seed_format(eggnog)),
        ("work_out writable", *writable(work_out)),
        ("scratch writable", *writable(scratch)),
        ("skani on PATH", shutil.which("skani") is not None, shutil.which("skani") or "not found (use --no-ani or rerun setup_env.sh)"),
    ]
    present = [d for d in ("/work", "/project") if Path(d).is_dir()]
    checks.append(("group storage roots (info)", True, ", ".join(present) or "neither /work nor /project found"))
    return checks


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--expect", type=int, default=1000, help="number of genomes (20 for a pilot check)")
    a = ap.parse_args()
    missing = [v for v in ("DATA", "WORK_OUT", "SCRATCH") if not os.environ.get(v)]
    if missing:
        sys.exit(f"{', '.join(missing)} not set: run `source scripts/sapelo2/config.sh` first")
    checks = run_checks(Path(os.environ["DATA"]), Path(os.environ["WORK_OUT"]),
                        Path(os.environ["SCRATCH"]), a.expect)
    for name, ok, detail in checks:
        print(f"{'PASS' if ok else 'FAIL'}  {name:30s} {detail}")
    failed = [c for c in checks if not c[1]]
    print(f"\n{len(checks) - len(failed)}/{len(checks)} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
