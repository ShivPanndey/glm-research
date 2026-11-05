#!/bin/bash
#SBATCH --job-name=call_genes
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --array=0-49
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
# Submitted by submit_all.sh, which exports the paths used below (see config.sh).
set -euo pipefail
module load "$PY_MODULE" 2>/dev/null || true
# shellcheck source=/dev/null
source "$VENV/bin/activate"
cd "$REPO"
python scripts/call_genes.py --genomes "$GENOMES" --annotations "$ANN" --out "$FAA"
