"""Smoke run of the whole pipeline on synthetic data, in a self-cleaning temp dir.

Used by `make smoke`. Scores from this data are meaningless; it only proves that
make_smoke_data -> featurize -> evaluate runs end to end.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EVAL_ARGS = sys.argv[1:]  # extra evaluate flags, e.g. --seeds 0


def run(*args: str, quiet: bool = True) -> str:
    out = subprocess.run([sys.executable, *args], cwd=ROOT, capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f"smoke step failed: {' '.join(args)}\n{out.stdout}\n{out.stderr}")
    return out.stdout


with tempfile.TemporaryDirectory() as d:
    run("tests/make_smoke_data.py", d)
    run("scripts/run_baseline.py", "featurize", "--genomes", f"{d}/genomes", "--out", f"{d}/cache",
        "--window", "5000", "--stride", "2500")
    print(run("scripts/run_baseline.py", "evaluate", "--cache", f"{d}/cache", "--annotations",
              f"{d}/ann", "--out", f"{d}/results", *EVAL_ARGS).strip().splitlines()[-1])
    print("smoke: ok")
