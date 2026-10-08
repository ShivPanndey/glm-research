"""Build a small pseudo-dataset in the lab's file formats to smoke-test the runner.

Each real genome is cut into 8 pseudo-MAGs; genes are called with Prodigal
(pyrodigal) and written as real Prodigal .faa; KO labels are SYNTHETIC (a hash
of each protein's first 25 aa), so any scores from this data are meaningless.
It only checks that featurize -> split -> retrieve -> evaluate runs end to end.
"""
import gzip, hashlib, sys
from pathlib import Path
import pyrodigal
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from tnf import read_fasta

out = Path(sys.argv[1]); [ (out/d).mkdir(parents=True, exist_ok=True) for d in ("genomes","faa","ann") ]
HDR = "\t".join(["query","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl","COG_category",
  "Description","Preferred_name","GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction",
  "KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs"])
gf = pyrodigal.GeneFinder(meta=True)
for g in sorted((ROOT/"data/genomes").glob("*.fna.gz")):
    seq = "".join(s for _, s in read_fasta(g))
    n = 8; L = len(seq) // n
    for c in range(n):
        gid = f"{g.name.split('.')[0]}_p{c}"
        sub = seq[c*L:(c+1)*L]
        (out/"genomes"/f"{gid}_genome.fa").write_text(f">{gid}_ctg1\n{sub}\n")
        genes = gf.find_genes(sub)
        with open(out/"faa"/f"{gid}.faa", "w") as fh:
            genes.write_translations(fh, sequence_id=f"{gid}_ctg1")
        rows = []
        for i, gene in enumerate(genes, 1):
            prot = gene.translate()
            h = int(hashlib.md5(prot[:25].encode()).hexdigest(), 16)
            ko = f"ko:K{h % 3000:05d}" if h % 4 else "-"
            rows.append("\t".join([f"{gid}_ctg1_{i}"] + ["x"]*10 + [ko] + ["-"]*9))
        (out/"ann"/f"{gid}.emapper.annotations").write_text("## emapper-2.1.9\n#" + HDR + "\n" + "\n".join(rows) + "\n")
print("smoke data ->", out)
