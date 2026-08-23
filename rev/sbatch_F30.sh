#!/bin/bash
#SBATCH --job-name=F30_scG4_vs_scATAC
#SBATCH --partition=cpu-single
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=github/rev/logs/F30_scG4_vs_scATAC_%j.log
#SBATCH --error=github/rev/logs/F30_scG4_vs_scATAC_%j.err

# Run F30: scG4 vs scATAC scatter plot (Union + Intersection)
echo "Starting F30 analysis at $(date)"
Rscript github/rev/rev_F30_scG4_vs_scATAC_scatter.R
echo "Finished F30 analysis at $(date)"

# Run F30: scG4 vs scATAC scatter plot (Union + Intersection)
echo "Starting F30 analysis at $(date)"
Rscript github/rev/rev_F30_scG4_vs_scATAC_scatter.R
echo "Finished F30 analysis at $(date)"