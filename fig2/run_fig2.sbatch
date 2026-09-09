#!/bin/bash
#SBATCH -A bw25j003
#SBATCH -p cpu-single
#SBATCH -n 4
#SBATCH --cpus-per-task=1
#SBATCH --mem=48GB
#SBATCH -t 2:00:00
#SBATCH -J fig2_R
#SBATCH -o logs/fig2_rerun_%j.out
#SBATCH -e logs/fig2_rerun_%j.err

# Full rerun of Figure 2 (panels S1, B-F) to validate the path fixes and
# the Panel B prediction-score bug fix.
#
# Run from the repository root (github/):
#   sbatch fig2/run_fig2.sbatch
#
# fig2.R sets ROOT <- normalizePath(getwd()) and derives all paths from it,
# so the submit working directory must be the repository root.

set -euo pipefail

REPO_ROOT="$(pwd)"
if [[ ! -f "${REPO_ROOT}/fig2/fig2.R" ]]; then
    echo "ERROR: run this script from the repository root (github/)" >&2
    exit 1
fi

mkdir -p "${REPO_ROOT}/logs"

module load devel/miniforge
conda activate /home/fr/fr_fr/fr_se1100/.conda/envs/ai
conda info >/dev/null 2>&1 || { echo "ERROR: conda env 'ai' not activatable" >&2; exit 1; }

cd "${REPO_ROOT}"
echo "[$(date +%H:%M:%S)] Starting fig2.R (full rerun: S1,B,C,D,E,F)"
Rscript fig2/fig2.R
echo "[$(date +%H:%M:%S)] fig2.R finished with exit code $?"