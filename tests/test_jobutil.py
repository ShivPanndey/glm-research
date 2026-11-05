import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from jobutil import FailureLog, failures_path, my_shard  # noqa: E402

from tnf import read_fasta  # noqa: E402


def test_shards_cover_each_item_once(monkeypatch):
    items = list(range(50))
    seen = []
    monkeypatch.setenv("SLURM_ARRAY_TASK_COUNT", "7")
    for task in range(7):
        monkeypatch.setenv("SLURM_ARRAY_TASK_ID", str(task))
        seen += my_shard(items)
    assert sorted(seen) == items


def test_failure_log_fresh_and_cleanup(tmp_path, monkeypatch):
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "0")
    stale = failures_path(tmp_path)
    stale.write_text("old\tstale\n")
    with FailureLog(tmp_path) as log:
        pass
    assert log.count == 0 and not stale.exists()
    with FailureLog(tmp_path) as log:
        log.record("g1", "boom")
    assert log.count == 1
    assert stale.read_text() == "g1\tboom\n"


def test_featurize_skips_bad_genome(tmp_path):
    genomes, out = tmp_path / "genomes", tmp_path / "cache"
    genomes.mkdir()
    seq = next(read_fasta(ROOT / "data/genomes/ecoli_K12.fna.gz"))[1][:20000]
    (genomes / "good.fa").write_text(f">c1\n{seq}\n")
    (genomes / "bad.fa").write_text("not a fasta\n")
    cmd = [sys.executable, str(ROOT / "scripts/run_baseline.py"), "featurize",
           "--genomes", str(genomes), "--out", str(out)]
    env = {"PATH": "/usr/bin:/bin", "SLURM_ARRAY_TASK_ID": "0", "SLURM_ARRAY_TASK_COUNT": "1"}
    first = subprocess.run(cmd, capture_output=True, text=True, env=env)
    assert first.returncode == 1, first.stderr
    assert (out / "good.npz").exists()
    assert (out / "failures.task0.tsv").read_text().startswith("bad\t")
    mtime = (out / "good.npz").stat().st_mtime_ns
    subprocess.run(cmd, capture_output=True, text=True, env=env)
    assert (out / "good.npz").stat().st_mtime_ns == mtime


def test_call_genes_id_failure_survives_rerun(tmp_path):
    genomes, ann, out = tmp_path / "genomes", tmp_path / "ann", tmp_path / "faa"
    genomes.mkdir()
    ann.mkdir()
    seq = next(read_fasta(ROOT / "data/genomes/ecoli_K12.fna.gz"))[1][:60000]
    (genomes / "g.fa").write_text(f">c1\n{seq}\n")
    (ann / "g.emapper.annotations").write_text("## emapper\n#query\tx\nc1_99999\tx\n")
    cmd = [sys.executable, str(ROOT / "scripts/call_genes.py"), "--genomes", str(genomes),
           "--annotations", str(ann), "--out", str(out)]
    env = {"PATH": "/usr/bin:/bin"}
    for _ in range(2):
        run = subprocess.run(cmd, capture_output=True, text=True, env=env)
        assert run.returncode == 1, run.stdout + run.stderr
        assert (out / "failures.task0.tsv").read_text().startswith("g\tids: 1/1")
        assert not (out / "g.faa").exists() and (out / "g.faa.failed").exists()
