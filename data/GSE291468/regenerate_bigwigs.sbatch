#!/bin/bash
#SBATCH -A bw25j003
#SBATCH -p cpu-single
#SBATCH -n 2
#SBATCH --mem=8GB
#SBATCH -t 2:00:00
#SBATCH -J make_bigwig
#SBATCH -o ../../logs/make_bigwig_%j.out
#SBATCH -e ../../logs/make_bigwig_%j.err

# Regenerate one GSM8836086 bigwig from CellRanger fragments.
# Usage (run from this directory, data/GSE291468/):
#   sbatch regenerate_bigwigs.sbatch LABEL OUTPUT_BW [SMOOTH_BP]
# Barcode definitions (LABEL -> file) are hardcoded below; they are the
# original authors' cell-group lists used for the GEO-deposited files.
#
# This script is OPTIONAL: the corrected, duplicate-free bigwigs already ship
# with the repository. Its primary purpose is regenerating the .smooth150
# variants of the Fig 3C track pair (see regenerate_bigwigs.sh).
#
set -euo pipefail

module load bio/bedtools/2.30.0
module load bio/samtools/1.19.2
source /gpfs/bwfor/software/common/devel/miniforge/24.9.2-0/etc/profile.d/conda.sh
conda activate deeptools

LABEL="${1:?Usage: sbatch $0 LABEL OUTPUT_BW [SMOOTH_BP]}"
OUTPUT_BW="${2:?Usage: sbatch $0 LABEL OUTPUT_BW [SMOOTH_BP]}"
SMOOTH="${3:-}"

DATA_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
FRAGMENTS="${DATA_DIR}/GSM8836086_GFP_sorted_mouse_brain/CellRanger/fragments.tsv.gz"
BC_DIR="${DATA_DIR}/scBridge_predictions/barcodes"
CLUSTER_DIR="${DATA_DIR}/Seurat_cluster_barcodes"
CHROM_SIZES="${DATA_DIR}/../genome/mm10.chrom.sizes.txt"
EFFECTIVE_GENOME_SIZE=2652783500

for f in "${FRAGMENTS}" "${CHROM_SIZES}"; do
    [[ -s "$f" ]] || { echo "ERROR: missing file: $f" >&2; exit 1; }
done

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
        C="${LABEL##*_}"
        BC_FILE="${CLUSTER_DIR}/barcodes_cluster_${C}.tsv" ;;
    *)
        echo "ERROR: unknown label '${LABEL}'. Valid labels are listed in this script." >&2
        exit 1 ;;
esac

[[ -s "${BC_FILE}" ]] || { echo "ERROR: missing barcode file: ${BC_FILE}" >&2; exit 1; }
BARCODES="$(cut -f1 "${BC_FILE}")"
N_BC=$(printf '%s\n' "${BARCODES}" | grep -c . || true)

OUTPUT_BW_DIR="$(dirname "${OUTPUT_BW}")"
mkdir -p "${OUTPUT_BW_DIR}"
WORKDIR="$(mktemp -d "${OUTPUT_BW_DIR}/.tmp_make_bigwig_XXXXXX")"
trap 'rm -rf "${WORKDIR}"' EXIT

echo "[${LABEL}] $(date +%H:%M:%S) Filtering fragments (${N_BC} barcodes from ${BC_FILE##*/})"
zcat "${FRAGMENTS}" \
    | awk -F'\t' 'NR==FNR { b[$1]=1; next } ($4 in b) { print $1"\t"$2"\t"$3"\t"$4"\t0\t+\t"$2"\t"$3"\t255,0,0\t"$5 }' \
        <(printf '%s\n' "${BARCODES}") - \
    | bedtools bedtobam -i stdin -g "${CHROM_SIZES}" \
    > "${WORKDIR}/unsorted.bam"

echo "[${LABEL}] $(date +%H:%M:%S) Sorting ($(samtools view -c "${WORKDIR}/unsorted.bam") reads)"
samtools sort -@ 2 -o "${WORKDIR}/coord.bam" "${WORKDIR}/unsorted.bam"

echo "[${LABEL}] $(date +%H:%M:%S) Fixmate + duplicate removal"
samtools collate -@ 2 -O "${WORKDIR}/coord.bam" "${WORKDIR}/collate" \
    | samtools fixmate -@ 2 -m - "${WORKDIR}/fixmate.bam"
samtools sort -@ 2 -o "${WORKDIR}/fixmate.coord.bam" "${WORKDIR}/fixmate.bam"
samtools markdup -@ 2 -r "${WORKDIR}/fixmate.coord.bam" "${WORKDIR}/dedup.bam"
samtools index "${WORKDIR}/dedup.bam"

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
    --numberOfProcessors 2

echo "[${LABEL}] $(date +%H:%M:%S) Created ${OUTPUT_BW}"
