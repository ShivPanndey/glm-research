#!/bin/bash
# One-time setup on Sapelo2 (run on a login node, from the repo root):
#   bash scripts/sapelo2/setup_env.sh
# Creates the Python venv, installs requirements, tries to make skani available, runs the tests.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${GLM_CONFIG:-$HERE/config.sh}"

module load "$PY_MODULE"
if [ ! -x "$VENV/bin/python" ]; then
  python -m venv "$VENV"
fi
# shellcheck source=/dev/null
source "$VENV/bin/activate"
pip install --upgrade pip -q
pip install -q -r "$REPO/requirements.txt"
echo "python: $(python --version), venv: $VENV"

# skani, first option that works: cluster module, conda (bioconda), prebuilt release binary.
if module load skani 2>/dev/null && command -v skani >/dev/null; then
  echo "skani: cluster module"
elif command -v skani >/dev/null; then
  echo "skani: already on PATH ($(command -v skani))"
elif command -v conda >/dev/null && conda create -y -q -p "$VENV/skani-conda" -c conda-forge -c bioconda skani >/dev/null; then
  ln -sf "$VENV/skani-conda/bin/skani" "$VENV/bin/skani"
  echo "skani: conda (bioconda)"
elif curl -fsSL -o "$VENV/bin/skani" "https://github.com/bluenote-1577/skani/releases/download/v0.2.2/skani" \
     && chmod +x "$VENV/bin/skani" && "$VENV/bin/skani" --version >/dev/null; then
  echo "skani: release binary v0.2.2"
else
  rm -f "$VENV/bin/skani"
  echo "skani: NOT installed. Try 'module spider skani' or 'module load Miniforge3' then rerun," >&2
  echo "       or submit with --no-ani (split then falls back to TNF distance, which is circular)." >&2
fi
command -v skani >/dev/null && skani --version

cd "$REPO" && python -m pytest -q tests
echo "setup done. Next: python scripts/sapelo2/preflight.py"
