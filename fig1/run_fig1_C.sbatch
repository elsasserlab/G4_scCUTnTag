#!/bin/bash
#SBATCH -A bw25j003
#SBATCH -p cpu-single
#SBATCH -n 4
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH -t 0:30:00
#SBATCH -J fig1_C
#SBATCH -o logs/fig1_C_%j.out
#SBATCH -e logs/fig1_C_%j.err

# Figure 1 Panel C: PQS overlap bar plot for MEF/ESC scG4 cluster peaks
# (with shared peaks bar) plus a regioneR permutation test (union of all
# peaks vs PQS, chr1, 1000 permutations) shown in the plot subtitle.
#
# Run from the repository root (github/):
#   sbatch fig1/run_fig1_C.sbatch
#
# fig1.R sets ROOT <- normalizePath(getwd()) and derives all paths from it,
# so the submit working directory must be the repository root.

set -euo pipefail

REPO_ROOT="$(pwd)"
if [[ ! -f "${REPO_ROOT}/fig1/fig1.R" ]]; then
    echo "ERROR: run this script from the repository root (github/)" >&2
    exit 1
fi

mkdir -p "${REPO_ROOT}/logs"

module load devel/miniforge
conda activate /home/fr/fr_fr/fr_se1100/.conda/envs/ai
conda info >/dev/null 2>&1 || { echo "ERROR: conda env 'ai' not activatable" >&2; exit 1; }

cd "${REPO_ROOT}"
echo "[$(date +%H:%M:%S)] Starting fig1.R panel C"
Rscript fig1/fig1.R --only=C
echo "[$(date +%H:%M:%S)] fig1.R finished with exit code $?"