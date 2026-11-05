#!/bin/bash
#SBATCH --job-name=skani
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
# Submitted by submit_all.sh, which exports the paths used below (see config.sh).
set -euo pipefail
module load "$PY_MODULE" 2>/dev/null || true
# shellcheck source=/dev/null
source "$VENV/bin/activate"
cd "$REPO"
module load skani 2>/dev/null || true
command -v skani >/dev/null || { echo "skani not found: run setup_env.sh, or submit with --no-ani" >&2; exit 1; }
find -L "$GENOMES" -maxdepth 1 \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) | sort > "$ANI_OUT/genomes.txt"
skani triangle -E -t "${SLURM_CPUS_PER_TASK:-1}" -l "$ANI_OUT/genomes.txt" -o "$ANI_OUT/skani_edges.tsv"
python scripts/ani_clusters.py --skani "$ANI_OUT/skani_edges.tsv" --genomes "$GENOMES" --out "$ANI_OUT"
