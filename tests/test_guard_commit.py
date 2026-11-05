import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HOOK = ROOT / ".claude" / "hooks" / "guard_commit.py"
sys.path.insert(0, str(HOOK.parent))
from guard_commit import blocked_paths  # noqa: E402

EM = chr(0x2014)


def run_hook(command: str) -> int:
    payload = json.dumps({"tool_input": {"command": command}})
    return subprocess.run([sys.executable, str(HOOK)], input=payload, text=True, cwd=ROOT,
                          capture_output=True).returncode


def test_blocked_paths_lab_data_only():
    paths = ["genomes/x.fa", "data/genomes/y.fna.gz", "tests/fixtures/a.faa",
             "eggnog/x/x.emapper.annotations", "cache/x.npz", "README.md"]
    assert blocked_paths(paths) == ["genomes/x.fa", "eggnog/x/x.emapper.annotations", "cache/x.npz"]


def test_blocks_em_dash_in_commit_message():
    assert run_hook(f'git commit -m "a {EM} b"') == 2


def test_ignores_non_commit_commands():
    assert run_hook("ls") == 0
