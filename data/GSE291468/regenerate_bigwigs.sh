#!/bin/bash
# OPTIONAL regeneration script for the GSM8836086 bigwigs.
#
# The corrected, duplicate-free bigwigs already ship in this directory and
# replace the GEO originals (which used a wrong effective genome size, kept
# PCR duplicates and were built from cut sites instead of fragments). You
# normally do NOT need to run anything here.
#
# Primary use case: regenerating the smoothed variants of the Fig 3C track
# pair (.smooth150 = 150 bp moving average on 10 bp bins). Only those two
# jobs are active; all other tracks are commented out for reference.
#
set -euo pipefail

DATA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${DATA_DIR}"

submit() {
    local label="$1" dst="$2" smooth="${3:-}"
    echo "Submitting: ${dst##*/}  [${label}]"
    sbatch regenerate_bigwigs.sbatch "${label}" "${dst}" ${smooth}
}

# --- Active: smoothed Fig 3C pair -----------------------------------------
submit Predicted_Astrocytes \
    GSM8836086_Predicted_Astrocytes_RPGC.smooth150.bw 150

submit Predicted_non-Astrocytes \
    GSM8836086_Predicted_non-Astrocytes_RPGC.smooth150.bw 150

# --- Reference: unsmoothed tracks (files are supplied; uncomment to redo) --
# submit Predicted_Astrocytes \
#     GSM8836086_Predicted_Astrocytes_RPGC.bw
#
# submit Predicted_non-Astrocytes \
#     GSM8836086_Predicted_non-Astrocytes_RPGC.bw
#
# submit Unreliable_predictions \
#     GSM8836086_Unreliable_predictions_RPGC.bw
#
# submit Predicted_OPC \
#     GSM8836086_Predicted_Oligodendrocyte_progenitor_cells_RPGC.bw
#
# submit Predicted_COP-NFOL \
#     GSM8836086_Predicted_Committed_and_Newly_formed_oligodendr_RPGC.bw
#
# submit Predicted_Mature_oligodendrocytes \
#     GSM8836086_Predicted_Mature_oligodendrocytes_RPGC.bw
#
# submit Predicted_Olfactory_ensheathing_cells \
#     GSM8836086_Predicted_Olfactory_ensheathing_cells_RPGC.bw
#
# submit Predicted_Pericytes \
#     GSM8836086_Predicted_Pericytes_RPGC.bw
#
# submit Predicted_Vascular_endothelial_cells \
#     GSM8836086_Predicted_Vascular_endothelial_cells_RPGC.bw
#
# submit cluster_0 \
#     GSM8836086_scG4CnT_MouseBrain_GFPpos_cluster_0_RPGC.bw
#
# submit cluster_1 \
#     GSM8836086_scG4CnT_MouseBrain_GFPpos_cluster_1_RPGC.bw
#
# submit cluster_2 \
#     GSM8836086_scG4CnT_MouseBrain_GFPpos_cluster_2_RPGC.bw
#
# submit cluster_3 \
#     GSM8836086_scG4CnT_MouseBrain_GFPpos_cluster_3_RPGC.bw

echo ""
echo "All jobs submitted. Check progress with: squeue -u \$USER"
