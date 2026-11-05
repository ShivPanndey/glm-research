import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "sapelo2"))
sys.path.insert(0, str(ROOT / "scripts"))
from preflight import run_checks  # noqa: E402

SEED_LINE = "c1_1\ts\t1e-5\t50\t1\t3\t1\t3\t90\t100\t100\n"


def make_data(root: Path) -> Path:
    data = root / "family_1000_dataset_final_v2"
    (data / "logs").mkdir(parents=True)
    (data / "genomes").mkdir()
    for gid in ("a", "b", "GCF_1.1_x_genomic"):
        (data / "genomes" / f"{gid}.fa").write_text(">c\nACGT\n")
        (data / "eggnog" / gid).mkdir(parents=True)
        (data / "eggnog" / gid / f"{gid}.emapper.annotations").write_text("#q\n")
    (data / "eggnog" / "a" / "a.emapper.seed_orthologs").write_text("#qseqid\n" + SEED_LINE)
    return data


def results(tmp_path, data, expect=3):
    return {name: ok for name, ok, _ in run_checks(data, tmp_path / "work", tmp_path / "scratch", expect=expect)}


def test_all_pass_except_missing_skani(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))  # no skani here
    r = results(tmp_path, make_data(tmp_path))
    assert r.pop("skani on PATH") is False
    assert all(r.values()), r


def test_missing_annotation_fails_counts(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    data = make_data(tmp_path)
    (data / "eggnog" / "b" / "b.emapper.annotations").unlink()
    r = results(tmp_path, data)
    assert r["annotation count"] is False and r["genomes without annotations"] is False
    assert r["genome count"] is True
