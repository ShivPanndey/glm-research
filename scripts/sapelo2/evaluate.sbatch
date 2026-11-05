#!/bin/bash
#SBATCH --job-name=evaluate
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
# genome-level only for the first score: no --faa (window traversal has its own spec)
args=(--cache "$CACHE" --annotations "$ANN" --seeds "${SEEDS:-0,1,2}" --out "$RESULTS")
if [ -n "${CLUSTERS:-}" ]; then args+=(--clusters "$CLUSTERS"); fi
python scripts/run_baseline.py evaluate "${args[@]}"
