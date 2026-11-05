import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from check_style import find_em_dashes  # noqa: E402

EM = "—"


def test_flags_em_dash(tmp_path):
    p = tmp_path / "a.md"
    p.write_text(f"ok\na {EM} b\n", encoding="utf-8")
    assert find_em_dashes([p]) == [(str(p), 2)]


def test_ignores_binary(tmp_path):
    p = tmp_path / "b.bin"
    p.write_bytes(b"\xff\xfe\x14\x20")
    assert find_em_dashes([p]) == []


def test_clean_file(tmp_path):
    p = tmp_path / "c.md"
    p.write_text("a, b\n", encoding="utf-8")
    assert find_em_dashes([p]) == []
