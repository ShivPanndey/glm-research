#!/usr/bin/env python
"""PreToolUse hook for Bash: block a `git commit` that would break house rules.

Blocks (exit 2, reason on stderr) when the commit command contains an em dash, when a staged
path looks like lab data (genomes, proteins, emapper outputs, archives, TNF caches) outside
data/genomes/ and tests/, or when a staged text file contains an em dash.
"""
from __future__ import annotations

import fnmatch
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from check_style import EM_DASH, find_em_dashes  # noqa: E402

LAB_PATTERNS = ("*.fa", "*.fa.gz", "*.faa", "*.faa.gz", "*.emapper.*", "*.tar", "*.tar.gz", "*.npz")
ALLOWED_PREFIXES = ("data/genomes/", "tests/")


def blocked_paths(paths: list[str]) -> list[str]:
    return [p for p in paths
            if not p.startswith(ALLOWED_PREFIXES)
            and any(fnmatch.fnmatch(Path(p).name, pat) for pat in LAB_PATTERNS)]


def staged() -> list[str]:
    out = subprocess.run(["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR"],
                         cwd=ROOT, capture_output=True, text=True)
    return out.stdout.split()


def main() -> int:
    try:
        command = json.load(sys.stdin).get("tool_input", {}).get("command", "")
    except json.JSONDecodeError:
        return 0
    if "git commit" not in command:
        return 0
    problems = []
    if EM_DASH in command:
        problems.append("the commit command contains an em dash (house style)")
    paths = staged()
    bad = blocked_paths(paths)
    if bad:
        problems.append("lab data staged, keep it off GitHub: " + ", ".join(bad[:10]))
    hits = find_em_dashes([ROOT / p for p in paths])
    if hits:
        problems.append("em dash in staged files: " + ", ".join(f"{Path(p).relative_to(ROOT)}:{n}"
                                                             for p, n in hits[:10]))
    if problems:
        print("Commit blocked: " + "; ".join(problems), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
