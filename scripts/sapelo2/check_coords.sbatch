#!/bin/bash
#SBATCH --job-name=check_coords
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
# Submitted by submit_all.sh, which exports the paths used below (see config.sh).
set -euo pipefail
module load "$PY_MODULE" 2>/dev/null || true
# shellcheck source=/dev/null
source "$VENV/bin/activate"
cd "$REPO"
python scripts/check_coords.py --faa "$FAA" --annotations "$ANN" --out "$COORDS"
