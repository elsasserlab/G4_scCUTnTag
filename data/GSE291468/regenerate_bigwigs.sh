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
# Usage (from the repository root):
#   bash data/GSE291468/regenerate_bigwigs.sh
#
# Requires bedtools, samtools and deepTools (bamCoverage) on PATH, or an
# activated conda environment named by CONDA_ENV (default: scg4).
set -euo pipefail

DATA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${DATA_DIR}"

CONDA_ENV="${CONDA_ENV:-scg4}"
if command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "${CONDA_ENV}"
fi
for tool in bedtools samtools bamCoverage; do
  command -v "${tool}" >/dev/null 2>&1 || { echo "ERROR: ${tool} not on PATH (see README.md)" >&2; exit 1; }
done

FRAGMENTS="${DATA_DIR}/GSM8836086_GFP_sorted_mouse_brain/CellRanger/fragments.tsv.gz"
BC_DIR="${DATA_DIR}/scBridge_predictions/barcodes"
CLUSTER_DIR="${DATA_DIR}/Seurat_cluster_barcodes"
CHROM_SIZES="${DATA_DIR}/../genome/mm10.chrom.sizes.txt"
EFFECTIVE_GENOME_SIZE=2652783500
THREADS="${THREADS:-2}"

for f in "${FRAGMENTS}" "${CHROM_SIZES}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing file: $f" >&2; exit 1; }
done

WORKROOT="$(mktemp -d "${DATA_DIR}/.tmp_make_bigwig_XXXXXX")"
trap 'rm -rf "${WORKROOT}"' EXIT

run_track() {
  local LABEL="$1"
  local OUTPUT_BW="$2"
  local SMOOTH="${3:-}"
  local BC_FILE
  local WORKDIR

  case "${LABEL}" in
    Predicted_Astrocytes)
      BC_FILE="${BC_DIR}/barcodes_Astrocytes.tsv" ;;
    Unreliable_predictions)
      BC_FILE="${BC_DIR}/barcodes_Novel.tsv" ;;
    Predicted_non-Astrocytes)
      BC_FILE="${BC_DIR}/barcodes_non_AST.tsv" ;;
    Predicted_OPC)
      BC_FILE="${BC_DIR}/barcodes_OPC.tsv" ;;
    Predicted_COP-NFOL)
      BC_FILE="${BC_DIR}/barcodes_COP-NFOL.tsv" ;;
    Predicted_Mature_oligodendrocytes)
      BC_FILE="${BC_DIR}/barcodes_Oligodendrocytes.tsv" ;;
    Predicted_Olfactory_ensheathing_cells)
      BC_FILE="${BC_DIR}/barcodes_OEC.tsv" ;;
    Predicted_Pericytes)
      BC_FILE="${BC_DIR}/barcodes_Pericytes.tsv" ;;
    Predicted_Vascular_endothelial_cells)
      BC_FILE="${BC_DIR}/barcodes_VEC.tsv" ;;
    cluster_[0-3])
      BC_FILE="${CLUSTER_DIR}/barcodes_cluster_${LABEL##*_}.tsv" ;;
    *)
      echo "ERROR: unknown label '${LABEL}'." >&2
      return 1 ;;
  esac

  [[ -s "${BC_FILE}" ]] || { echo "ERROR: missing barcode file: ${BC_FILE}" >&2; return 1; }

  local BARCODES
  BARCODES="$(cut -f1 "${BC_FILE}")"
  local N_BC
  N_BC=$(printf '%s\n' "${BARCODES}" | grep -c . || true)

  WORKDIR="${WORKROOT}/${LABEL}"
  mkdir -p "${WORKDIR}" "$(dirname "${OUTPUT_BW}")"

  echo "[${LABEL}] $(date +%H:%M:%S) Filtering fragments (${N_BC} barcodes from ${BC_FILE##*/})"
  gzip -dc "${FRAGMENTS}" \
    | awk -F'\t' 'NR==FNR { b[$1]=1; next } ($4 in b) { print $1"\t"$2"\t"$3"\t"$4"\t0\t+\t"$2"\t"$3"\t255,0,0\t"$5 }' \
        <(printf '%s\n' "${BARCODES}") - \
    | bedtools bedtobam -i stdin -g "${CHROM_SIZES}" \
    > "${WORKDIR}/unsorted.bam"

  echo "[${LABEL}] $(date +%H:%M:%S) Sorting ($(samtools view -c "${WORKDIR}/unsorted.bam") reads)"
  samtools sort -@ "${THREADS}" -o "${WORKDIR}/coord.bam" "${WORKDIR}/unsorted.bam"

  echo "[${LABEL}] $(date +%H:%M:%S) Fixmate + duplicate removal"
  samtools collate -@ "${THREADS}" -O "${WORKDIR}/coord.bam" "${WORKDIR}/collate" \
    | samtools fixmate -@ "${THREADS}" -m - "${WORKDIR}/fixmate.bam"
  samtools sort -@ "${THREADS}" -o "${WORKDIR}/fixmate.coord.bam" "${WORKDIR}/fixmate.bam"
  samtools markdup -@ "${THREADS}" -r "${WORKDIR}/fixmate.coord.bam" "${WORKDIR}/dedup.bam"
  samtools index "${WORKDIR}/dedup.bam"

  local BEFORE AFTER DUP_PCT
  BEFORE=$(samtools view -c "${WORKDIR}/fixmate.coord.bam")
  AFTER=$(samtools view -c "${WORKDIR}/dedup.bam")
  DUP_PCT=$(awk "BEGIN {printf \"%.2f\", (1-${AFTER}/${BEFORE})*100}")
  echo "[${LABEL}] $(date +%H:%M:%S) Dedup done: ${BEFORE} -> ${AFTER} reads (${DUP_PCT}% removed)"

  echo "[${LABEL}] $(date +%H:%M:%S) bamCoverage -> ${OUTPUT_BW}"
  bamCoverage \
    --bam "${WORKDIR}/dedup.bam" \
    --outFileName "${OUTPUT_BW}" \
    --outFileFormat bigwig \
    --binSize 10 \
    ${SMOOTH:+--smoothLength ${SMOOTH}} \
    --normalizeUsing RPGC \
    --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE} \
    --numberOfProcessors "${THREADS}"

  echo "[${LABEL}] $(date +%H:%M:%S) Created ${OUTPUT_BW}"
}

# --- Active: smoothed Fig 3C pair -----------------------------------------
run_track Predicted_Astrocytes \
  GSM8836086_Predicted_Astrocytes_RPGC.smooth150.bw 150

run_track Predicted_non-Astrocytes \
  GSM8836086_Predicted_non-Astrocytes_RPGC.smooth150.bw 150

# --- Reference: unsmoothed tracks (files are supplied; uncomment to redo) --
# run_track Predicted_Astrocytes \
#   GSM8836086_Predicted_Astrocytes_RPGC.bw
#
# run_track Predicted_non-Astrocytes \
#   GSM8836086_Predicted_non-Astrocytes_RPGC.bw
#
# run_track Unreliable_predictions \
#   GSM8836086_Unreliable_predictions_RPGC.bw
