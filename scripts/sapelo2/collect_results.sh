#!/bin/bash
# Copy the small outputs of a run into the repo so they can be committed and pushed:
#   bash scripts/sapelo2/collect_results.sh /scratch/$USER/glm/pilot     (or .../run)
# Copies JSON/TSV summaries and a SLURM job summary only. Refuses any file over 5 MB, so
# genomes, proteins, and caches can never be swept in by mistake.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${GLM_CONFIG:-$HERE/config.sh}"

RUN_DIR="${1:?usage: collect_results.sh RUN_DIR (for example $SCRATCH/pilot)}"
RUN_DIR="${RUN_DIR%/}"
cd "$REPO"
DEST="results/lab/$(date +%F)_$(git rev-parse --short HEAD)_$(basename "$RUN_DIR")"
mkdir -p "$DEST"

copy() {  # src dest-name
  [ -f "$1" ] || { echo "  missing: $1"; return 0; }
  if [ "$(stat -c %s "$1")" -gt 5242880 ]; then
    echo "refusing to copy $1: over 5 MB, not a summary file" >&2
    exit 1
  fi
  cp "$1" "$DEST/$2"
  echo "  copied: $2"
}

copy "$RUN_DIR/results/baseline_report.json" baseline_report.json
copy "$RUN_DIR/ani/ani_summary.json" ani_summary.json
copy "$RUN_DIR/ani/clusters_ani95.tsv" clusters_ani95.tsv
copy "$RUN_DIR/coords/coords_check.json" coords_check.json
copy "$RUN_DIR/coords/coords_check.tsv" coords_check.tsv
FAA_DIR="$RUN_DIR/prodigal"            # pilot runs keep proteins inside the run folder
[ -d "$FAA_DIR" ] || FAA_DIR="$WORK_OUT"  # full runs keep them in group storage
for dir in "$RUN_DIR/tnf_cache" "$FAA_DIR"; do
  for f in "$dir"/failures.task*.tsv; do
    [ -e "$f" ] && copy "$f" "$(basename "$dir")_$(basename "$f")"
  done
done
if [ -f "$RUN_DIR/job_ids.txt" ] && command -v sacct >/dev/null; then
  sacct -j "$(tr ' ' ',' < "$RUN_DIR/job_ids.txt")" --format=JobID,JobName,State,Elapsed,MaxRSS > "$DEST/sacct.txt"
  echo "  copied: sacct.txt"
fi
echo "Now: git add $DEST && git commit -m 'Add $(basename "$RUN_DIR") results' && git push"
