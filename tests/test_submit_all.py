"""submit_all.sh wiring, checked with a fake `sbatch` that records its arguments."""
import os
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SUBMIT = ROOT / "scripts/sapelo2/submit_all.sh"

FAKE_SBATCH = """#!/bin/bash
n=$(( $(cat "$FAKE_DIR/n" 2>/dev/null || echo 100) + 1 ))
echo $n > "$FAKE_DIR/n"
echo "$n $*" >> "$FAKE_DIR/calls.log"
echo $n
"""


@pytest.fixture
def env(tmp_path):
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "sbatch").write_text(FAKE_SBATCH)
    (bindir / "sbatch").chmod(0o755)
    data = tmp_path / "data"
    (data / "genomes").mkdir(parents=True)
    (data / "eggnog").mkdir()
    for i in range(25):
        (data / "genomes" / f"g{i:02d}.fa").write_text(">c\nACGT\n")
    cfg = tmp_path / "config.sh"
    cfg.write_text(f'export REPO="{ROOT}"\nexport DATA="{data}"\nexport WORK_OUT="{tmp_path}/prodigal"\n'
                   f'export SCRATCH="{tmp_path}/scratch"\nexport VENV="{tmp_path}/venv"\n'
                   'export PY_MODULE=Python\nexport MAIL=nobody@example.com\n')
    e = dict(os.environ, PATH=f"{bindir}:{os.environ['PATH']}", GLM_CONFIG=str(cfg), FAKE_DIR=str(tmp_path))
    return tmp_path, e


def submit(env, *args):
    tmp, e = env
    run = subprocess.run(["bash", str(SUBMIT), *args], env=e, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    calls = {}
    for line in (tmp / "calls.log").read_text().splitlines():
        jid, rest = line.split(" ", 1)
        name = rest.split("--job-name=")[1].split()[0]
        calls[name] = (jid, rest)
    return calls


def test_pilot_chain_with_ani(env):
    calls = submit(env, "--pilot")
    assert set(calls) == {"call_genes", "check_coords", "featurize", "skani", "evaluate"}
    fe, sk, cg = calls["featurize"][0], calls["skani"][0], calls["call_genes"][0]
    assert f"--dependency=afterok:{fe}:{sk}" in calls["evaluate"][1]
    assert f"--dependency=afterany:{cg}" in calls["check_coords"][1]
    assert "--array=0-3" in calls["featurize"][1]
    assert "CLUSTERS=" in calls["evaluate"][1] and "pilot/ani/clusters_ani95.tsv" in calls["evaluate"][1]
    links = list((env[0] / "scratch/pilot/genomes").iterdir())
    assert len(links) == 20 and all(p.is_symlink() for p in links)
    assert sorted(p.name for p in links)[-1] == "g19.fa"


def test_pilot_no_ani(env):
    calls = submit(env, "--pilot", "--no-ani")
    assert "skani" not in calls and len(calls) == 4
    assert f"--dependency=afterok:{calls['featurize'][0]} " in calls["evaluate"][1] + " "
    assert "CLUSTERS=" not in calls["evaluate"][1]


def test_full_run_uses_50_tasks(env):
    calls = submit(env)
    assert "--array=0-49" in calls["call_genes"][1]
    assert not (env[0] / "scratch/pilot").exists()
