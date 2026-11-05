import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from ani_clusters import ani_distance, nearest_ani, read_skani_sparse  # noqa: E402

HEADER = "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
ROWS = "g/a.fa\tg/b.fa\t97.5\t90\t90\tx\ty\ng/a.fa\tg/c.fa\t85.0\t50\t50\tx\tz\n"
IDS = ["a", "b", "c", "d"]


def test_parse_by_header_and_nearest(tmp_path):
    p = tmp_path / "e.tsv"
    p.write_text(HEADER + ROWS)
    pairs = read_skani_sparse(p)
    assert pairs == [("a", "b", 97.5), ("a", "c", 85.0)]
    assert nearest_ani(IDS, pairs) == {"a": 97.5, "b": 97.5, "c": 85.0, "d": 0.0}


def test_distance_missing_pairs_are_far():
    D = ani_distance(IDS, [("a", "b", 97.5), ("a", "c", 85.0)])
    assert D[0, 1] == D[1, 0] == pytest.approx(0.025)
    assert D[0, 3] == 1.0 and D[2, 2] == 0.0


def test_missing_column_named(tmp_path):
    p = tmp_path / "e.tsv"
    p.write_text("Ref_file\tQuery_file\tidentity\n")
    with pytest.raises(ValueError, match="ANI"):
        read_skani_sparse(p)


def test_cli_clusters_and_summary(tmp_path):
    genomes = tmp_path / "g"
    genomes.mkdir()
    for gid in IDS:
        (genomes / f"{gid}.fa").write_text(">c\nACGT\n")
    edges = tmp_path / "e.tsv"
    edges.write_text(HEADER + ROWS)
    out = tmp_path / "out"
    run = subprocess.run([sys.executable, str(ROOT / "scripts/ani_clusters.py"), "--skani", str(edges),
                          "--genomes", str(genomes), "--out", str(out)], capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    cl = dict(line.split("\t") for line in (out / "clusters_ani95.tsv").read_text().splitlines())
    assert cl["a"] == cl["b"] and len({cl["a"], cl["c"], cl["d"]}) == 3
    s = json.loads((out / "ani_summary.json").read_text())
    assert s["n_genomes"] == 4 and s["n_clusters"] == 3 and s["n_nearest_below_80"] == 1
    assert s["cluster_size_counts"] == {"1": 2, "2": 1}


def test_chained_relatives_share_a_cluster(tmp_path):
    from ani_clusters import cluster_ani
    pairs = [("a", "b", 96.0), ("b", "c", 96.0)]  # a~c not reported by skani
    labels = cluster_ani(["a", "b", "c", "d"], pairs)
    assert labels[0] == labels[1] == labels[2] != labels[3]
