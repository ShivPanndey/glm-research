import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from check_coords import check_genome, protein_lengths, seed_qends  # noqa: E402

FAA = (">c1_1 # 1 # 12 # 1 # ID=1_1\nMKV*\n"
       ">c1_2 # 20 # 37 # 1 # ID=1_2\nMKV\nLA\n")
SEEDS = ("#qseqid\tsseqid\tevalue\tbitscore\tqstart\tqend\tsstart\tsend\tpident\tqcov\tscov\n"
         "c1_1\ts\t1e-5\t50\t1\t3\t1\t3\t90\t100\t100\n"
         "c1_1\ts\t1e-5\t50\t1\t2\t1\t2\t90\t100\t100\n"
         "c1_2\ts\t1e-5\t50\t1\t6\t1\t6\t90\t100\t100\n"
         "c1_3\ts\t1e-5\t50\t1\t2\t1\t2\t90\t100\t100\n")


def write(tmp_path):
    (tmp_path / "g.faa").write_text(FAA)
    (tmp_path / "g.emapper.seed_orthologs").write_text(SEEDS)


def test_lengths_exclude_stop_and_join_lines(tmp_path):
    write(tmp_path)
    assert protein_lengths(tmp_path / "g.faa") == {"c1_1": 3, "c1_2": 5}


def test_qend_is_max_per_query(tmp_path):
    write(tmp_path)
    assert seed_qends(tmp_path / "g.emapper.seed_orthologs") == {"c1_1": 3, "c1_2": 6, "c1_3": 2}


def test_check_genome_counts(tmp_path):
    write(tmp_path)
    r = check_genome(protein_lengths(tmp_path / "g.faa"),
                     seed_qends(tmp_path / "g.emapper.seed_orthologs"))
    assert (r["n_genes"], r["n_ok"], r["n_missing"], r["n_too_long"]) == (3, 1, 1, 1)
    assert set(r["examples"]) == {"c1_2", "c1_3"}


def test_cli_reports_missing_faa(tmp_path):
    faa, ann, out = tmp_path / "faa", tmp_path / "eggnog", tmp_path / "out"
    faa.mkdir()
    for gid in ("g", "h"):
        (ann / gid).mkdir(parents=True)
        (ann / gid / f"{gid}.emapper.seed_orthologs").write_text(SEEDS)
    (faa / "g.faa").write_text(FAA)
    run = subprocess.run([sys.executable, str(ROOT / "scripts/check_coords.py"), "--faa", str(faa),
                          "--annotations", str(ann), "--out", str(out)], capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    summary = json.loads((out / "coords_check.json").read_text())
    assert summary["n_genomes"] == 1 and summary["n_missing_faa"] == 1
    assert summary["n_genomes_all_ok"] == 0
    assert (out / "coords_check.tsv").read_text().splitlines()[1].startswith("g\t3\t1\t1\t1")
