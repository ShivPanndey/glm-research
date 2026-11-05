"""evaluate CLI: several split seeds, a summary, and run metadata (synthetic labels)."""
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tnf import read_fasta  # noqa: E402

RUN = [sys.executable, str(ROOT / "scripts/run_baseline.py")]


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    """12 pseudo-genomes (2 slices of each public genome) with KOs tied to the source genome."""
    d = tmp_path_factory.mktemp("eval")
    (d / "genomes").mkdir()
    (d / "ann").mkdir()
    for s, src in enumerate(sorted((ROOT / "data/genomes").glob("*.fna.gz"))):
        seq = next(read_fasta(src))[1]
        for part in range(2):
            gid = f"{src.name.split('.')[0]}_s{part}"
            (d / "genomes" / f"{gid}.fa").write_text(f">{gid}_c\n{seq[part * 40000:(part + 1) * 40000]}\n")
            kos = [f"ko:K{s:03d}{j:02d}" for j in range(5)] + ["ko:K99999"]
            rows = "".join(f"{gid}_c_{i}\t" + "\t".join(["x"] * 10) + f"\t{ko}\n" for i, ko in enumerate(kos, 1))
            (d / "ann" / f"{gid}.emapper.annotations").write_text("## emapper\n#query\n" + rows)
    subprocess.run(RUN + ["featurize", "--genomes", str(d / "genomes"), "--out", str(d / "cache")],
                   check=True, capture_output=True)
    return d


def evaluate(d, out, *extra):
    subprocess.run(RUN + ["evaluate", "--cache", str(d / "cache"), "--annotations", str(d / "ann"),
                          "--out", str(out), *extra], check=True, capture_output=True, text=True)
    return json.loads((out / "baseline_report.json").read_text())


def test_evaluate_two_seeds(dataset, tmp_path):
    report = evaluate(dataset, tmp_path, "--seeds", "0,1")
    assert set(report["per_seed"]) == {"0", "1"}
    meta = report["meta"]
    assert meta["cluster_source"] == "tnf_fallback (circular)"
    assert meta["n_failures"] == 0 and meta["n_matched"] == 12 and meta["seeds"] == [0, 1]
    assert len(meta["git_commit"]) >= 7
    for m in ("tnf_genome_knn", "frequency_prior"):
        f1 = report["summary"][m]["micro_f1"]
        assert f1["min"] <= f1["mean"] <= f1["max"]


def test_evaluate_ani_clusters_label(dataset, tmp_path):
    clusters = tmp_path / "clusters.tsv"
    ids = sorted(p.stem for p in (dataset / "cache").glob("*.npz"))
    clusters.write_text("".join(f"{g}\t{i}\n" for i, g in enumerate(ids)))
    report = evaluate(dataset, tmp_path / "out", "--clusters", str(clusters))
    assert report["meta"]["cluster_source"] == "ani95"
    assert set(report["per_seed"]) == {"0"}


def test_clusters_file_must_match_ids(dataset, tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("nope1\t1\nnope2\t2\n")
    run = subprocess.run(RUN + ["evaluate", "--cache", str(dataset / "cache"), "--annotations",
                                str(dataset / "ann"), "--clusters", str(bad), "--out", str(tmp_path / "o")],
                         capture_output=True, text=True)
    assert run.returncode != 0 and "clusters file" in (run.stdout + run.stderr)


def test_meta_counts_ids_in_clusters_file(dataset, tmp_path):
    ids = sorted(p.stem for p in (dataset / "cache").glob("*.npz"))
    part = tmp_path / "part.tsv"
    part.write_text("".join(f"{g}\t{i}\n" for i, g in enumerate(ids[:6])))
    report = evaluate(dataset, tmp_path / "o", "--clusters", str(part))
    assert report["meta"]["n_in_clusters_file"] == 6
