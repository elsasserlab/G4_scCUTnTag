#!/usr/bin/env bash
# Build a continuous fragment-coverage bigwig (RPGC) for GSE198467 from the
# raw CellRanger-ATAC fragments file.
#
# Modeled on data/GSE291468/regenerate_bigwigs.sh. Unlike the scG4 tracks we
# do NOT markdup: identical genomic coordinates across different cells are
# genuine multi-cell signal, not PCR duplicates.
#
# Usage (from the repository root):
#   bash data/GSE198467/make_atac_bigwig.sh
#
# Output:
#   data/GSE198467/GSE198467_ATAC_fragments.RPGC.bw  (10 bp bins, RPGC)
#
# Requires bedtools, samtools and deepTools (bamCoverage) on PATH, or an
# activated conda environment named by CONDA_ENV (default: scg4).
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-scg4}"
if command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "${CONDA_ENV}"
fi
for tool in bedtools samtools bamCoverage; do
  command -v "${tool}" >/dev/null 2>&1 || { echo "ERROR: ${tool} not on PATH (see README.md)" >&2; exit 1; }
done

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FRAG="${HERE}/GSE198467_ATAC_fragments.tsv.gz"
CHROM_SIZES="${HERE}/../genome/mm10.chrom.sizes.txt"
OUTBW="${HERE}/GSE198467_ATAC_fragments.RPGC.bw"
EFFECTIVE_GENOME_SIZE=2652783500
THREADS="${THREADS:-4}"

for f in "${FRAG}" "${CHROM_SIZES}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing file: $f" >&2; exit 1; }
done

WORKDIR="$(mktemp -d "${HERE}/.tmp_make_atac_bw_XXXXXX")"
trap 'rm -rf "${WORKDIR}"' EXIT

echo "[1/3] $(date +%H:%M:%S) Building BED6 (canonical chr) from fragments"
gzip -dc "${FRAG}" \
  | awk -F'\t' '$1 ~ /^chr([1-9]|1[0-9]|X|Y)$/ {print $1"\t"$2"\t"$3"\t"$4"\t0\t+"}' \
  > "${WORKDIR}/fragments.bed"

echo "[2/3] $(date +%H:%M:%S) bedtobam + sort ($(wc -l < "${WORKDIR}/fragments.bed") fragments)"
bedtools bedtobam -i "${WORKDIR}/fragments.bed" -g "${CHROM_SIZES}" \
  > "${WORKDIR}/unsorted.bam"
samtools sort -@ "${THREADS}" -o "${WORKDIR}/sorted.bam" "${WORKDIR}/unsorted.bam"
samtools index -@ "${THREADS}" "${WORKDIR}/sorted.bam"

echo "[3/3] $(date +%H:%M:%S) bamCoverage RPGC -> ${OUTBW}"
bamCoverage \
  --bam "${WORKDIR}/sorted.bam" \
  --outFileName "${OUTBW}" \
  --outFileFormat bigwig \
  --binSize 10 \
  --normalizeUsing RPGC \
  --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE} \
  --numberOfProcessors "${THREADS}"

echo "[done] $(date +%H:%M:%S) Created ${OUTBW}"
