# Every Sapelo2 path lives here. On first login edit GROUP_DIR and PY_MODULE, nothing else.
# Sourced by setup_env.sh, submit_all.sh, and collect_results.sh. Set GLM_CONFIG to use
# another file instead (the tests do this).

# Your lab's group storage. Find it with: ls /work /project
export GROUP_DIR="${GROUP_DIR:-/work/CHANGE_ME}"
# Python module name. Find it with: module spider Python
export PY_MODULE="${PY_MODULE:-Python}"   # pick the full name module spider shows

export REPO="${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
export DATA="$GROUP_DIR/$USER/glm/family_1000_dataset_final_v2"   # extracted tarball (kept)
export WORK_OUT="$GROUP_DIR/$USER/glm/prodigal"                   # re-called proteins (kept)
export SCRATCH="/scratch/$USER/glm"                               # caches, logs (disposable)
export VENV="$HOME/envs/tnf"                                      # python venv
export MAIL="sp74395@uga.edu"                                     # SLURM end/fail emails
