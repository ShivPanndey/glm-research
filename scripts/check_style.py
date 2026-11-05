#!/usr/bin/env python
"""House-style gate: no em dashes in tracked text files or in a given string.

  python scripts/check_style.py              # every file `git ls-files` lists
  python scripts/check_style.py a.md b.py    # only these files
  python scripts/check_style.py --text "commit message"

Files that are not valid UTF-8 (gzipped genomes, images) are skipped.
Exits 1 and prints path:line for every hit.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

EM_DASH = "—"


def find_em_dashes(paths: list[Path]) -> list[tuple[str, int]]:
    hits = []
    for p in paths:
        try:
            text = Path(p).read_text(encoding="utf-8")
        except (UnicodeDecodeError, IsADirectoryError, FileNotFoundError):
            continue
        hits += [(str(p), i) for i, line in enumerate(text.splitlines(), 1) if EM_DASH in line]
    return hits


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="*")
    ap.add_argument("--text", default=None)
    a = ap.parse_args()
    if a.text is not None:
        if EM_DASH in a.text:
            print("em dash in text")
            return 1
        return 0
    paths = a.paths or subprocess.run(["git", "ls-files"], capture_output=True, text=True,
                                      check=True).stdout.split()
    hits = find_em_dashes([Path(p) for p in paths])
    for path, line in hits:
        print(f"{path}:{line}: em dash (house style: use commas, colons, or parentheses)")
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
