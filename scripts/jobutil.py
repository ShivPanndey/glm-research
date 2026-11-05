"""Helpers shared by the SLURM array scripts (featurize, call_genes).

Each array task handles every n-th genome and records per-genome failures in its own
`failures.task<ID>.tsv` (id<TAB>reason) inside the output folder, so concurrent tasks never
write the same file. The file is rewritten on every run and removed when nothing failed.
"""
from __future__ import annotations

import os
from pathlib import Path


def task_id() -> int:
    return int(os.environ.get("SLURM_ARRAY_TASK_ID") or 0)


def task_count() -> int:
    return int(os.environ.get("SLURM_ARRAY_TASK_COUNT") or 1)


def my_shard(items: list) -> list:
    return items[task_id()::task_count()]


def failures_path(out: Path) -> Path:
    return Path(out) / f"failures.task{task_id()}.tsv"


class FailureLog:
    """Context manager: `with FailureLog(out) as log: log.record(gid, reason)`."""

    def __init__(self, out: Path):
        self.path = failures_path(out)
        self.count = 0

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = open(self.path, "w")
        return self

    def record(self, gid: str, reason: str) -> None:
        self._fh.write(f"{gid}\t{' '.join(str(reason).split())}\n")
        self._fh.flush()
        self.count += 1

    def __exit__(self, *exc):
        self._fh.close()
        if self.count == 0:
            self.path.unlink()
        return False
