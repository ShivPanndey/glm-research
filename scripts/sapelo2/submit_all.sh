#!/bin/bash
# Submit the whole pipeline with one command. Each job waits for the ones it needs.
#
#   bash scripts/sapelo2/submit_all.sh            full run, 1,000 genomes
#   bash scripts/sapelo2/submit_all.sh --pilot    first 20 genomes by name, small arrays
#   add --no-ani to skip skani (the split then falls back to TNF distance, which is circular)
#
#   call_genes (array) --afterany--> check_coords            (report only)
#   featurize  (array) --+
#   skani      (1 job) --+--afterok--> evaluate (seeds 0,1,2)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${GLM_CONFIG:-$HERE/config.sh}"

PILOT=0
ANI=1
for arg in "$@"; do
  case "$arg" in
    --pilot) PILOT=1 ;;
    --no-ani) ANI=0 ;;
    *) echo "unknown option: $arg (use --pilot and/or --no-ani)" >&2; exit 2 ;;
  esac
done

GENOMES="$DATA/genomes"
ANN="$DATA/eggnog"
RUN="$SCRATCH/run"
FAA="$WORK_OUT"
ARRAY="0-49"
if [ "$PILOT" = 1 ]; then
  RUN="$SCRATCH/pilot"
  FAA="$RUN/prodigal"
  ARRAY="0-3"
  GENOMES="$RUN/genomes"
  mkdir -p "$GENOMES"
  find "$GENOMES" -maxdepth 1 -type l -delete
  find "$DATA/genomes" -maxdepth 1 \( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' \) -printf '%f\n' \
    | sort | head -n 20 | while read -r f; do ln -s "$DATA/genomes/$f" "$GENOMES/$f"; done
fi
CACHE="$RUN/tnf_cache"
ANI_OUT="$RUN/ani"
COORDS="$RUN/coords"
RESULTS="$RUN/results"
LOGS="$RUN/logs"
mkdir -p "$LOGS" "$CACHE" "$FAA" "$ANI_OUT" "$COORDS" "$RESULTS"
export GENOMES ANN FAA CACHE ANI_OUT COORDS RESULTS
unset CLUSTERS

sub() {  # sbatch with our defaults; prints the job id
  local id
  id=$(sbatch --parsable --mail-user="$MAIL" "$@")
  echo "${id%%;*}"
}

CG=$(sub --job-name=call_genes --array="$ARRAY" --output="$LOGS/%x_%A_%a.out" --export=ALL "$HERE/call_genes.sbatch")
CC=$(sub --job-name=check_coords --dependency="afterany:$CG" --output="$LOGS/%x_%j.out" --export=ALL "$HERE/check_coords.sbatch")
FE=$(sub --job-name=featurize --array="$ARRAY" --output="$LOGS/%x_%A_%a.out" --export=ALL "$HERE/featurize.sbatch")
if [ "$ANI" = 1 ]; then
  SK=$(sub --job-name=skani --output="$LOGS/%x_%j.out" --export=ALL "$HERE/skani.sbatch")
  EV=$(sub --job-name=evaluate --dependency="afterok:$FE:$SK" --output="$LOGS/%x_%j.out" \
           --export="ALL,CLUSTERS=$ANI_OUT/clusters_ani95.tsv" "$HERE/evaluate.sbatch")
  IDS="$CG $CC $FE $SK $EV"
else
  EV=$(sub --job-name=evaluate --dependency="afterok:$FE" --output="$LOGS/%x_%j.out" --export=ALL "$HERE/evaluate.sbatch")
  IDS="$CG $CC $FE $EV"
fi
echo "$IDS" > "$RUN/job_ids.txt"

echo "Submitted to $RUN"
echo "  call_genes $CG, check_coords $CC, featurize $FE${SK:+, skani $SK}, evaluate $EV"
echo "Watch:   squeue -u \$USER"
echo "Logs:    $LOGS"
echo "Collect: bash scripts/sapelo2/collect_results.sh $RUN"
