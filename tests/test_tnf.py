import random
from collections import Counter

import numpy as np

from tnf import (CANON_KMERS, N_FEATURES, KOVocab, count_kmers, genome_ko_set,
                 parse_emapper, parse_prodigal_faa, window_counts, window_ko_sets)

COMP = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(COMP)[::-1]

def brute_counts(seq):
    c = Counter()
    for i in range(len(seq) - 3):
        k = seq[i:i+4]
        if set(k) <= set("ACGT"):
            c[min(k, rc(k))] += 1
    return np.array([c[k] for k in CANON_KMERS], dtype=float)

def rand_seq(n, seed=0):
    r = random.Random(seed); return "".join(r.choice("ACGT") for _ in range(n))

def test_dimension():
    assert N_FEATURES == 136 and len(set(CANON_KMERS)) == 136
    assert sum(1 for k in CANON_KMERS if k == rc(k)) == 16

def test_matches_bruteforce_with_Ns():
    s = rand_seq(5000, 1)
    s = s[:1000] + "NNRY" + s[1000:]
    assert np.array_equal(count_kmers(s), brute_counts(s))

def test_strand_invariance():
    s = rand_seq(3000, 2)
    assert np.array_equal(count_kmers(s), count_kmers(rc(s)))

def test_window_counts_close_to_bruteforce():
    s = rand_seq(20000, 3)
    starts, C = window_counts(s, 2000, 500)
    for i in [0, 7, len(starts) - 2]:
        st = starts[i]
        exact = brute_counts(s[st:st + 2000 + 3])  # block trick includes 3 overhanging k-mers
        assert np.array_equal(C[i], exact)

def test_emapper_and_faa(tmp_path):
    em = tmp_path / "g.emapper.annotations"
    header = ["query","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl","COG_category",
              "Description","Preferred_name","GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module",
              "KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs"]
    rows = [["c1_1"] + ["x"]*10 + ["ko:K00001,ko:K00002"] + ["-"]*9,
            ["c1_2"] + ["x"]*10 + ["-"] + ["-"]*9,
            ["c2_1"] + ["x"]*10 + ["ko:K00003"] + ["-"]*9]
    em.write_text("## emapper-2.1.9\n#" + "\t".join(header) + "\n" +
                  "\n".join("\t".join(r) for r in rows) + "\n")
    kos = parse_emapper(em)
    assert kos["c1_1"] == {"K00001", "K00002"} and kos["c1_2"] == set()
    assert genome_ko_set(kos) == {"K00001", "K00002", "K00003"}
    faa = tmp_path / "g.faa"
    faa.write_text(">c1_1 # 101 # 400 # 1 # ID=1_1;partial=00\nMK*\n"
                   ">c1_2 # 1201 # 1500 # -1 # ID=1_2;partial=00\nMK*\n"
                   ">c2_1 # 51 # 350 # 1 # ID=2_1;partial=00\nMK*\n")
    orfs = parse_prodigal_faa(faa)
    assert (orfs[0].contig, orfs[0].start, orfs[0].end) == ("c1", 100, 400)
    w = window_ko_sets(np.array(["c1", "c1", "c2"], dtype=object),
                       np.array([0, 1000, 0]), np.array([1000, 2000, 1000]), orfs, kos)
    assert w == [{"K00001", "K00002"}, set(), {"K00003"}]
    Y = KOVocab(w).to_matrix(w)
    assert Y.shape == (3, 3) and Y.sum() == 3

def test_grouped_split_no_leakage_and_balanced():
    from tnf import grouped_split
    groups = np.repeat(np.arange(12), 8)  # 12 clusters of 8
    tr, va, te = grouped_split(groups, seed=1)
    assert set(groups[tr]).isdisjoint(groups[va]) and set(groups[tr]).isdisjoint(groups[te])
    assert set(groups[va]).isdisjoint(groups[te])
    assert len(te) >= 8 and len(va) >= 8 and len(tr) + len(va) + len(te) == 96
