#!/usr/bin/env python3
"""Summarize family_1000_dataset.tar.gz (or the plain .tar a browser leaves) and pull a small, matched sample.

Usage (Python 3.8+, no extra packages):
    python peek_dataset.py family_1000_dataset.tar.gz            # default 20 genomes
    python peek_dataset.py family_1000_dataset.tar.gz --n 30

Writes, next to the archive:
    dataset_summary.txt   structure overview (paste or upload this)
    dataset_sample.zip    up to N genomes, each with its matching annotation and
                          Prodigal files when present, spread across folders
Nothing else is extracted. Two streaming passes over the archive (a few minutes).
"""
import argparse
import collections
import io
import re
import tarfile
import zipfile
from pathlib import Path

GENOME_RE = re.compile(r"\.(fa|fna|fasta)(\.gz)?$", re.I)
KINDS = [("annotation", re.compile(r"\.emapper\.annotations", re.I)),
         ("proteins", re.compile(r"\.faa(\.gz)?$", re.I)),
         ("gff", re.compile(r"\.gff3?(\.gz)?$", re.I)),
         ("genes_nt", re.compile(r"\.(ffn|fnn)(\.gz)?$", re.I)),
         ("genome", GENOME_RE)]
MAX_SAMPLE_BYTES = 28 * 1024 * 1024  # keep the zip uploadable


def kind_of(name):
    for k, rx in KINDS:
        if rx.search(name):
            return k
    return "other"


def norm_id(name):
    base = Path(name).name.split(".")[0]
    return re.sub(r"_(genome|proteins|genes|prodigal|annotations?)$", "", base)


def ext_of(name):
    n = Path(name).name
    m = re.search(r"(\.[A-Za-z0-9_]+)(\.gz)?$", n)
    return (m.group(0) if m else "(none)").lower()


def head_text(data, n_lines=6, width=160):
    try:
        if data[:2] == b"\x1f\x8b":
            import gzip
            data = gzip.decompress(data)
        lines = data.decode("utf-8", "replace").splitlines()[:n_lines]
        return "\n".join("    " + l[:width] for l in lines)
    except Exception as e:  # noqa: BLE001
        return f"    (could not read: {e})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("archive")
    ap.add_argument("--n", type=int, default=20)
    a = ap.parse_args()
    arc = Path(a.archive)
    out_dir = arc.parent

    # ---------------- pass 1: inventory
    files = []  # (name, size)
    with tarfile.open(arc, "r|*") as tf:
        for m in tf:
            if m.isfile():
                files.append((m.name, m.size))
    by_ext = collections.Counter(ext_of(n) for n, _ in files)
    by_kind = collections.Counter(kind_of(n) for n, _ in files)
    dirs = collections.Counter(str(Path(n).parent) for n, _ in files)
    depth1 = collections.Counter(Path(n).parts[0] if len(Path(n).parts) > 1 else "." for n, _ in files)
    depth2 = collections.Counter("/".join(Path(n).parts[:2]) for n, _ in files if len(Path(n).parts) > 2)

    groups = collections.defaultdict(dict)  # id -> kind -> (name, size)
    for n, s in files:
        k = kind_of(n)
        if k != "other":
            groups[norm_id(n)].setdefault(k, (n, s))
    complete = {g: v for g, v in groups.items() if "genome" in v and "annotation" in v}

    # spread the sample across parent folders (round robin), smallest genomes first
    by_parent = collections.defaultdict(list)
    for g, v in complete.items():
        by_parent[str(Path(v["genome"][0]).parent)].append(g)
    for p in by_parent:
        by_parent[p].sort(key=lambda g: complete[g]["genome"][1])
    chosen, budget = [], 0
    parents = sorted(by_parent)
    while len(chosen) < a.n and any(by_parent.values()):
        for p in parents:
            if by_parent[p] and len(chosen) < a.n:
                g = by_parent[p].pop(0)
                size = sum(s for _, s in complete[g].values())
                if budget + size <= MAX_SAMPLE_BYTES:
                    chosen.append(g)
                    budget += size
        if budget > MAX_SAMPLE_BYTES * 0.95:
            break
    want = {name for g in chosen for name, _ in complete[g].values()}
    preview_kind = {}  # one example file per kind for the summary
    for n, s in files:
        k = kind_of(n)
        if k not in preview_kind and s < 50_000_000:
            preview_kind[k] = n
    want_preview = set(preview_kind.values())

    # ---------------- pass 2: extract sample + previews
    previews = {}
    zpath = out_dir / "dataset_sample.zip"
    with tarfile.open(arc, "r|*") as tf, zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as zf:
        for m in tf:
            if not m.isfile() or (m.name not in want and m.name not in want_preview):
                continue
            data = tf.extractfile(m).read()
            if m.name in want:
                zf.writestr(m.name, data)
            if m.name in want_preview:
                previews[m.name] = head_text(data)

    # ---------------- summary
    total = sum(s for _, s in files)
    L = [f"Archive: {arc.name}",
         f"Files: {len(files):,}   Uncompressed total: {total / 1e9:.2f} GB",
         "", "Top-level folders (file counts):"]
    L += [f"  {d}: {c:,}" for d, c in depth1.most_common(30)]
    L += ["", "Second-level folders (first 40):"]
    L += [f"  {d}: {c:,}" for d, c in sorted(depth2.items())[:40]]
    L += [f"  ... {len(depth2):,} second-level folders total", "",
          f"Distinct parent folders: {len(dirs):,}", "", "File types (by extension):"]
    L += [f"  {e}: {c:,}" for e, c in by_ext.most_common(25)]
    L += ["", "Recognized kinds:"] + [f"  {k}: {c:,}" for k, c in by_kind.most_common()]
    n_faa = sum(1 for v in complete.values() if "proteins" in v)
    n_gff = sum(1 for v in complete.values() if "gff" in v)
    L += ["", f"Genome ids with genome + annotation: {len(complete):,}",
          f"  ...of which also have Prodigal .faa: {n_faa:,}   .gff: {n_gff:,}",
          f"Genome ids with only some files: {len(groups) - len(complete):,}", ""]
    gsizes = sorted(v["genome"][1] for v in complete.values())
    if gsizes:
        L.append(f"Genome file size (bytes): min {gsizes[0]:,}  median {gsizes[len(gsizes)//2]:,}  "
                 f"max {gsizes[-1]:,}")
    L += ["", "Example paths (first 25):"] + [f"  {n}" for n, _ in files[:25]]
    L += ["", "File previews (first lines):"]
    for k, n in preview_kind.items():
        L += [f"  [{k}] {n}", previews.get(n, "    (not read)"), ""]
    L += [f"Sample: {len(chosen)} genome ids, {budget / 1e6:.1f} MB uncompressed -> {zpath.name}"]
    L += [f"  {g}: {', '.join(sorted(complete[g]))}" for g in chosen]
    (out_dir / "dataset_summary.txt").write_text("\n".join(L))
    print("\n".join(L[:40]))
    print(f"\nWrote {out_dir / 'dataset_summary.txt'} and {zpath}")


if __name__ == "__main__":
    main()
