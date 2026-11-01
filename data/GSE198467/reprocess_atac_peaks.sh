#!/usr/bin/env bash
# ===========================================================================
# Reproduce summit-containing consensus scATAC peaks for the GSE198467
# nano-CT ATAC modality from the raw fragment file.
#
# The author-provided Seurat object contains only broad peak windows
# (median ~3.7 kb, no summit offsets). MACS2 peak calling on the raw
# fragments recovers conventional narrowPeak intervals WITH summit offsets
# (column 10), so the scATAC reference can be compared with the scG4
# datasets using the same "original peak width" policy.
#
# Peak-calling parameters follow the CellRanger-ATAC / ENCODE convention for
# scATAC fragments (Tn5 shift/extsize) with --call-summits.
#
# Stringency is tunable via QVALUE (default 1e-6 yields ~93k peaks, in the
# realistic 50-100k range).
#
# Requires: macs2 (conda env 'scg4'), gzip. Run from anywhere:
#   QVALUE=1e-6 bash data/GSE198467/reprocess_atac_peaks.sh
#
# Outputs:
#   data/GSE198467/macs2_out/fragments_peaks.narrowPeak   (canonical chr only)
#   data/GSE198467/GSE198467_atac_consensus_peaks.narrowPeak  (final, indexable)
# ===========================================================================
set -euo pipefail

QVALUE="${QVALUE:-1e-6}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FRAG="${HERE}/GSE198467_ATAC_fragments.tsv.gz"
OUTDIR="${HERE}/macs2_out"
FINAL="${HERE}/GSE198467_atac_consensus_peaks.narrowPeak"
LOG="${OUTDIR}/macs2.log"

if [[ ! -s "${FRAG}" ]]; then
  echo "ERROR: fragments file not found: ${FRAG}" >&2
  echo "Download it first (see data/download_data.sh)." >&2
  exit 1
fi
if ! command -v macs2 >/dev/null 2>&1; then
  echo "ERROR: macs2 not on PATH (activate the scg4 environment first)." >&2
  exit 1
fi
gzip -t "${FRAG}"
mkdir -p "${OUTDIR}"

CANONICAL_RE='^chr([1-9]|1[0-9]|X|Y)$'

echo "[1/3] Building canonical-chromosome 3-column BED from fragments (large)..."
time gzip -dc "${FRAG}" \
  | awk -F'\t' '$1 ~ /^chr([1-9]|1[0-9]|X|Y)$/ {print $1"\t"$2"\t"$3}' \
  > "${OUTDIR}/fragments.bed"
echo "      BED lines: $(wc -l < "${OUTDIR}/fragments.bed")"

echo "[2/3] MACS2 peak calling (--nomodel --shift -75 --extsize 150 --call-summits, q=${QVALUE})..."
time macs2 callpeak -t "${OUTDIR}/fragments.bed" -f BED \
  --nomodel --shift -75 --extsize 150 --keep-dup all --call-summits \
  -g mm -q "${QVALUE}" --outdir "${OUTDIR}" -n fragments > "${LOG}" 2>&1
tail -5 "${LOG}"

echo "[3/3] Move canonical narrowPeak to final location..."
mv "${OUTDIR}/fragments_peaks.narrowPeak" "${FINAL}"

echo
echo "Summary (final file): ${FINAL}"
echo "  peaks:            $(wc -l < "${FINAL}")"
echo "  median width:     $(awk -F'\t' '{s+=($3-$2)} END {print int((s/NR)+0.5)}' "${FINAL}")"
echo "  summit col range: $(awk -F'\t' 'NR==1{min=$10;max=$10} {if($10<min)min=$10; if($10>max)max=$10} END {print min"-"max}' "${FINAL}")"

rm -f "${OUTDIR}/fragments.bed"
echo "Done."