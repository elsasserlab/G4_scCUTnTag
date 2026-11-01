#!/usr/bin/env bash
#
# download_data.sh - Download all data required to reproduce figures 1-4
#                    and revision supplementary figures.
#
# Data layout (After download):
#   data/
#   ├── genome/              Reference genome files (chrom.sizes, GTF, eG4, PQS)
#   ├── cCRE/                Candidate cis-regulatory elements (bundled)
#   ├── GSE291468/           Our scG4 CUT&Tag + bulk CUT&Tag (primary accession)
#   │   ├── GSM8836086_GFP_sorted_mouse_brain/CellRanger/   Cell Ranger outputs
#   │   ├── GSM8836087_unsorted_mouse_brain/CellRanger/     Cell Ranger outputs
#   │   ├── GSM8836088_mESC_MEF/CellRanger/                 Cell Ranger outputs
#   │   └── (bulk bigwigs, broadPeaks, Seurat objects, etc.)
#   ├── GSE163484/           Bartosovic et al. mouse brain scRNA-seq
#   ├── GSE198467/           Brain snATAC-seq reference + cluster-specific peaks
#   ├── GSE149080/           Published bulk mESC ATAC-seq
#   ├── GSE211123/           Published bulk 3T3 ATAC-seq
#   ├── GSE217860/           Published bulk 3T3 G4 CUT&Tag
#   └── GSE90894/            Chronis et al. mESC/MEF RNA-seq RPKM table
#
# Usage (from the repository root):
#   bash data/download_data.sh              # download everything
#   bash data/download_data.sh --fig1       # only Fig 1 data
#   bash data/download_data.sh --fig2       # only Fig 2 data
#   bash data/download_data.sh --fig4       # only Fig 4 data
#   bash data/download_data.sh --revisions  # only revision data
#
# Requirements: curl, gzip, tabix (for fragment indexes), git-lfs
#               (repository files tracked via Git LFS are fetched with
#               'git lfs pull', see Quick Start in README.md)
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${ROOT}")"
DATA="${ROOT}"

# Files that failed to download (the script continues and reports at the end)
FAILED=()

is_lfs_pointer() {
  # Detect an un-materialized Git LFS pointer file (checked out without LFS)
  [[ -f "$1" ]] && head -c 130 "$1" 2>/dev/null | grep -q "git-lfs.github.com/spec"
}

# ---------------------------------------------------------------------
# Preflight: warn about repository files that are still Git LFS pointers.
# Author-provided files (Seurat objects, bigwigs) are only shipped through
# LFS and cannot be re-downloaded from GEO - they MUST be materialized
# with 'git lfs pull' before any figure can be reproduced.
# ---------------------------------------------------------------------
if command -v git >/dev/null 2>&1 && git -C "${REPO_ROOT}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  if command -v git-lfs >/dev/null 2>&1; then
    POINTERS=$(git -C "${REPO_ROOT}" lfs ls-files 2>/dev/null | awk '$2 == "-" {print $3}' | tr '\n' ' ')
    if [[ -n "${POINTERS}" ]]; then
      echo "[warn] The following repository files are un-materialized Git LFS pointers:"
      echo "       ${POINTERS}"
      echo "       Run 'git lfs install && git lfs pull' in ${REPO_ROOT} first."
      echo "       Files that exist on GEO will be re-downloaded by this script;"
      echo "       author-provided files can ONLY be obtained via 'git lfs pull'."
      echo ""
    fi
  else
    echo "[warn] git-lfs is not installed. Several required data files (bigwigs,"
    echo "       author-provided Seurat objects) are tracked via Git LFS and will"
    echo "       appear as tiny pointer files. Install git-lfs and run"
    echo "       'git lfs install && git lfs pull' in ${REPO_ROOT}."
    echo ""
  fi
fi

download() {
  local url="$1" dest="$2"
  mkdir -p "$(dirname "${dest}")"
  if [[ "${FORCE_DOWNLOAD:-0}" != "1" && -s "${dest}" ]] && ! is_lfs_pointer "${dest}"; then
    echo "[skip] $(basename "${dest}") already exists"
    return 0
  fi
  if is_lfs_pointer "${dest}"; then
    echo "[warn] $(basename "${dest}") is an un-materialized Git LFS pointer;"
    echo "       prefer 'git lfs pull', re-downloading from GEO anyway"
  fi
  echo "[download] $(basename "${dest}")"
  if curl -fL --retry 5 --retry-delay 2 --continue-at - -o "${dest}.part" "${url}"; then
    mv "${dest}.part" "${dest}"
  else
    rm -f "${dest}.part"
    echo "[error] failed to download $(basename "${dest}") - continuing"
    FAILED+=("${dest}")
  fi
}

download_gzip_as_plain() {
  local url="$1" dest="$2"
  local compressed="${dest}.gz.part"
  mkdir -p "$(dirname "${dest}")"
  if [[ "${FORCE_DOWNLOAD:-0}" != "1" && -s "${dest}" ]] && ! is_lfs_pointer "${dest}"; then
    echo "[skip] $(basename "${dest}") already exists"
    return 0
  fi
  echo "[download] $(basename "${dest}")"
  if curl -fL --retry 5 --retry-delay 2 --continue-at - -o "${compressed}" "${url}"; then
    if gzip -dc "${compressed}" > "${dest}.part"; then
      rm -f "${compressed}"
      mv "${dest}.part" "${dest}"
    else
      rm -f "${compressed}" "${dest}.part"
      echo "[error] failed to decompress $(basename "${dest}") - continuing"
      FAILED+=("${dest}")
    fi
  else
    rm -f "${compressed}"
    echo "[error] failed to download $(basename "${dest}") - continuing"
    FAILED+=("${dest}")
  fi
}

download_and_unzip() {
  local url="$1" dest="$2"
  local compressed="${dest}.part"
  mkdir -p "$(dirname "${dest}")"
  if [[ "${FORCE_DOWNLOAD:-0}" != "1" && -s "${dest}" ]] && ! is_lfs_pointer "${dest}"; then
    echo "[skip] $(basename "${dest}") already exists"
    return 0
  fi
  echo "[download] $(basename "${dest}")"
  if ! curl -fL --retry 5 --retry-delay 2 --continue-at - -o "${compressed}" "${url}"; then
    rm -f "${compressed}"
    echo "[error] failed to download $(basename "${dest}") - continuing"
    FAILED+=("${dest}")
    return 0
  fi
  if [[ "${dest}" == *.gz ]]; then
    # The destination is itself a gzip file: store it as-is (do NOT gunzip,
    # otherwise the file would be corrupt/unreadable by tabix and R).
    mv "${compressed}" "${dest}"
  else
    # The destination is plain text: gunzip the downloaded file into place.
    mv "${compressed}" "${dest}.gz"
    gzip -d "${dest}.gz"
  fi
}

# Index fragment files with tabix (required by several figure scripts)
index_fragments() {
  local frag
  if command -v tabix >/dev/null 2>&1; then
    for frag in "$@"; do
      if [[ -s "${frag}" && ! -s "${frag}.tbi" ]]; then
        echo "[tabix] indexing ${frag}"
        tabix -f -p bed "${frag}" || echo "[warn] tabix failed for ${frag}"
      fi
    done
  else
    echo "[warn] tabix not found - fragment indexes (.tbi) not built."
    echo "       Install tabix and run: tabix -p bed <fragments.tsv.gz>"
    echo "       (the .tbi files shipped with the repository can also be used)"
  fi
}

# Parse filter
FIG="${1:-all}"
echo "============================================"
echo "  scG4 data download"
echo "  Destination: ${DATA}"
echo "  Filter: ${FIG}"
echo "============================================"

# ======================================================================
# GSE291468: Our scG4 CUT&Tag + bulk CUT&Tag (primary GEO accession)
# ======================================================================
if [[ "$FIG" == "all" || "$FIG" == "--fig1" || "$FIG" == "--fig2" || \
      "$FIG" == "--fig3" || "$FIG" == "--fig4" || "$FIG" == "--revisions" ]]; then
  echo ""
  echo "--- GSE291468: scG4 CUT&Tag + bulk CUT&Tag ---"

  # === GSM8836088: mESC-MEF mixture processed files (Figure 1) ===
  # These are the authoritative GEO outputs for Figure 1. Keep the local
  # filenames stable for the figure scripts while preserving the accession
  # and sample identity in this mapping.
  echo ""
  echo "[GSM8836088 mESC-MEF processed outputs for Figure 1]"
  GSM8836088_URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836088/suppl"
  download "${GSM8836088_URL}/GSM8836088_scG4CnT_mixture_mES_3T3_Seurat_object.Rds" \
           "${DATA}/GSE291468/GSM8836088_mESCMEF_Seurat_object.Rds"
  download "${GSM8836088_URL}/GSM8836088_scG4CnT_mixture_mES_3T3_cluster_0.bam_RPGC.bw" \
           "${DATA}/GSE291468/GSM8836088_cluster0_RPGC.bw"
  download "${GSM8836088_URL}/GSM8836088_scG4CnT_mixture_mES_3T3_cluster_1.bam_RPGC.bw" \
           "${DATA}/GSE291468/GSM8836088_cluster1_RPGC.bw"
  download_gzip_as_plain \
    "${GSM8836088_URL}/GSM8836088_scG4CnT_mixture_mES_3T3_cluster_0_peaks.narrowPeak.gz" \
    "${DATA}/GSE291468/GSM8836088_cluster_0_peaks.narrowPeak"
  download_gzip_as_plain \
    "${GSM8836088_URL}/GSM8836088_scG4CnT_mixture_mES_3T3_cluster_1_peaks.narrowPeak.gz" \
    "${DATA}/GSE291468/GSM8836088_cluster_1_peaks.narrowPeak"
  download "${GSM8836088_URL}/GSM8836088_scG4CnT_mixture_mES_3T3_fragments.tsv.gz" \
           "${DATA}/GSE291468/GSM8836088_mESC_MEF/CellRanger/fragments.tsv.gz"
  download "${GSM8836088_URL}/GSM8836088_scG4CnT_mixture_mES_3T3_raw_peak_matrix.h5" \
           "${DATA}/GSE291468/GSM8836088_mESC_MEF/CellRanger/raw_peak_matrix.h5"

  # === Cell Ranger outputs (fragments, peaks, barcodes) ===
  echo ""
  echo "[Cell Ranger outputs]"

  # GFP+ sorted mouse brain (GSM8836086)
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_scG4CnT_MouseBrain_GFPpos_fragments.tsv.gz" \
           "${DATA}/GSE291468/GSM8836086_GFP_sorted_mouse_brain/CellRanger/fragments.tsv.gz"
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_scG4CnT_MouseBrain_GFPpos_raw_peak_bc_matrix.h5" \
           "${DATA}/GSE291468/GSM8836086_GFP_sorted_mouse_brain/CellRanger/raw_peak_bc_matrix.h5"

  # Unsorted mouse brain (GSM8836087)
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836087/suppl/GSM8836087_scG4CnT_MouseBrain_unsorted_fragments.tsv.gz" \
           "${DATA}/GSE291468/GSM8836087_unsorted_mouse_brain/CellRanger/fragments.tsv.gz"
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836087/suppl/GSM8836087_scG4CnT_MouseBrain_unsorted_raw_peak_matrix.h5" \
           "${DATA}/GSE291468/GSM8836087_unsorted_mouse_brain/CellRanger/raw_peak_matrix.h5"

  # === Seurat objects (author-processed) ===
  echo ""
  echo "[Seurat objects]"
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_scG4CnT_MouseBrain_GFPpos_Seurat_object.Rds" \
           "${DATA}/GSE291468/GSM8836086_GFPpos_Seurat_object.Rds"
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836087/suppl/GSM8836087_scG4CnT_MouseBrain_unsorted_Seurat_object.Rds" \
           "${DATA}/GSE291468/GSM8836087_unsorted_Seurat_object.Rds"
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836088/suppl/GSM8836088_scG4CnT_mixture_mES_3T3_Seurat_object.Rds" \
           "${DATA}/GSE291468/GSM8836088_mESCMEF_Seurat_object.Rds"

  # === Cluster-specific narrowPeak files ===
  echo ""
  echo "[Cluster-specific peaks]"
  # Unsorted cluster 0 and 1 (used by Fig 4 panels B-D). These are GSM8836087
  # supplementary files from GSE291468, not files from GSE198467.
  download_and_unzip "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836087/suppl/GSM8836087_scG4CnT_MouseBrain_unsorted_cluster_0_peaks.narrowPeak.gz" \
           "${DATA}/GSE291468/GSM8836087_unsorted_cluster_0_peaks.narrowPeak"
  download_and_unzip "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836087/suppl/GSM8836087_scG4CnT_MouseBrain_unsorted_cluster_1_peaks.narrowPeak.gz" \
           "${DATA}/GSE291468/GSM8836087_unsorted_cluster_1_peaks.narrowPeak"

  # mESC/MEF cluster peaks (used by Fig 1)
  download_and_unzip "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836088/suppl/GSM8836088_scG4CnT_mixture_mES_3T3_cluster_0_peaks.narrowPeak.gz" \
           "${DATA}/GSE291468/GSM8836088_cluster_0_peaks.narrowPeak"
  download_and_unzip "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836088/suppl/GSM8836088_scG4CnT_mixture_mES_3T3_cluster_1_peaks.narrowPeak.gz" \
           "${DATA}/GSE291468/GSM8836088_cluster_1_peaks.narrowPeak"

  # GFP+ sorted cluster 0-3 narrowPeaks
  download_and_unzip "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_scG4CnT_MouseBrain_GFPpos_cluster_0_peaks.narrowPeak.gz" \
           "${DATA}/GSE291468/GSM8836086_GFPpos_cluster_0_peaks.narrowPeak"
  download_and_unzip "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_scG4CnT_MouseBrain_GFPpos_cluster_1_peaks.narrowPeak.gz" \
           "${DATA}/GSE291468/GSM8836086_GFPpos_cluster_1_peaks.narrowPeak"
  download_and_unzip "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_scG4CnT_MouseBrain_GFPpos_cluster_2_peaks.narrowPeak.gz" \
           "${DATA}/GSE291468/GSM8836086_GFPpos_cluster_2_peaks.narrowPeak"
  download_and_unzip "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_scG4CnT_MouseBrain_GFPpos_cluster_3_peaks.narrowPeak.gz" \
           "${DATA}/GSE291468/GSM8836086_GFPpos_cluster_3_peaks.narrowPeak"

  # === Cluster-specific RPGC bigWigs ===
  echo ""
  echo "[Cluster-specific bigWigs]"
  # GSM8836086 does not publish a whole-population RPGC bigWig in its
  # supplementary directory. Do not request the former 404 URL. The
  # cluster-specific tracks below are the published processed signal files.
  #
  # MASKED: the GSM8836086 single-cell bigWigs below are NOT downloaded from
  # GEO anymore. The published versions were built with a wrong effective
  # genome size, retained PCR duplicates and used cut-site coverage; we ship
  # regenerated, duplicate-free replacements (correct mm10 eGS, fragment
  # coverage) directly in data/GSE291468/ under the same filenames — see
  # data/GSE291468/JOB_STATUS.md and regenerate_bigwigs.sh there.
  # Unsorted clusters
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836087/suppl/GSM8836087_scG4CnT_MouseBrain_unsorted_cluster_0.bam_RPGC.bw" \
           "${DATA}/GSE291468/GSM8836087_cluster0_RPGC.bw"
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836087/suppl/GSM8836087_scG4CnT_MouseBrain_unsorted_cluster_1.bam_RPGC.bw" \
           "${DATA}/GSE291468/GSM8836087_cluster1_RPGC.bw"
  # GFP+ sorted clusters (masked — replaced by regenerated files)
  # download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_scG4CnT_MouseBrain_GFPpos_cluster_0_RPGC.bw" \
  #          "${DATA}/GSE291468/GSM8836086_scG4CnT_MouseBrain_GFPpos_cluster_0_RPGC.bw"
  # download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_scG4CnT_MouseBrain_GFPpos_cluster_1_RPGC.bw" \
  #          "${DATA}/GSE291468/GSM8836086_scG4CnT_MouseBrain_GFPpos_cluster_1_RPGC.bw"
  # download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_scG4CnT_MouseBrain_GFPpos_cluster_2_RPGC.bw" \
  #          "${DATA}/GSE291468/GSM8836086_scG4CnT_MouseBrain_GFPpos_cluster_2_RPGC.bw"
  # download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_scG4CnT_MouseBrain_GFPpos_cluster_3_RPGC.bw" \
  #          "${DATA}/GSE291468/GSM8836086_scG4CnT_MouseBrain_GFPpos_cluster_3_RPGC.bw"

  # === Cell-type predicted bigWigs (GFP+ sorted) ===
  # MASKED: all ten files in this block are replaced by our regenerated,
  # duplicate-free versions shipped in data/GSE291468/ (same filenames).
  #
  # download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_Predicted_Astrocytes_RPGC.bw" \
  #          "${DATA}/GSE291468/GSM8836086_Predicted_Astrocytes_RPGC.bw"
  # download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_Predicted_non-Astrocytes_RPGC.bw" \
  #          "${DATA}/GSE291468/GSM8836086_Predicted_non-Astrocytes_RPGC.bw"
  # download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_Predicted_Oligodendrocyte_progenitor_cells_RPGC.bw" \
  #          "${DATA}/GSE291468/GSM8836086_Predicted_Oligodendrocyte_progenitor_cells_RPGC.bw"
  # download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_Predicted_Committed_and_Newly_formed_oligodendr_RPGC.bw" \
  #          "${DATA}/GSE291468/GSM8836086_Predicted_Committed_and_Newly_formed_oligodendr_RPGC.bw"
  # download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_Predicted_Mature_oligodendrocytes_RPGC.bw" \
  #          "${DATA}/GSE291468/GSM8836086_Predicted_Mature_oligodendrocytes_RPGC.bw"
  # download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_Predicted_Olfactory_ensheathing_cells_RPGC.bw" \
  #          "${DATA}/GSE291468/GSM8836086_Predicted_Olfactory_ensheathing_cells_RPGC.bw"
  # download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_Predicted_Pericytes_RPGC.bw" \
  #          "${DATA}/GSE291468/GSM8836086_Predicted_Pericytes_RPGC.bw"
  # download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_Predicted_Vascular_endothelial_cells_RPGC.bw" \
  #          "${DATA}/GSE291468/GSM8836086_Predicted_Vascular_endothelial_cells_RPGC.bw"
  # download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_highly_predicted_Astrocytes_RPGC.bw" \
  #          "${DATA}/GSE291468/GSM8836086_highly_predicted_Astrocytes_RPGC.bw"
  # download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836086/suppl/GSM8836086_Unreliable_predictions_RPGC.bw" \
  #          "${DATA}/GSE291468/GSM8836086_Unreliable_predictions_RPGC.bw"

  # === Bulk G4 CUT&Tag ===
  echo ""
  echo "[Bulk G4 CUT&Tag]"
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836082/suppl/GSM8836082_bulkG4CnT_mESC_rep1.bw" \
           "${DATA}/GSE291468/GSM8836082_bulkG4CnT_mESC_rep1.bw"
  download_and_unzip "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836082/suppl/GSM8836082_bulkG4CnT_mESC_rep1.broadPeak.gz" \
           "${DATA}/GSE291468/GSM8836082_bulkG4CnT_mESC_rep1.broadPeak"
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836084/suppl/GSM8836084_bulkG4CnT_3T3_rep1.bw" \
           "${DATA}/GSE291468/GSM8836084_bulkG4CnT_3T3_rep1.bw"
  download_and_unzip "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836084/suppl/GSM8836084_bulkG4CnT_3T3_rep1.broadPeak.gz" \
           "${DATA}/GSE291468/GSM8836084_bulkG4CnT_3T3_rep1.broadPeak"

  # Rep 2 bigwigs
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836083/suppl/GSM8836083_bulkG4CnT_mESC_rep2.bw" \
           "${DATA}/GSE291468/GSM8836083_bulkG4CnT_mESC_rep2.bw"
  download_and_unzip "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836083/suppl/GSM8836083_bulkG4CnT_mESC_rep2.broadPeak.gz" \
           "${DATA}/GSE291468/GSM8836083_bulkG4CnT_mESC_rep2.broadPeak"
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836085/suppl/GSM8836085_bulkG4CnT_3T3_rep2.bw" \
           "${DATA}/GSE291468/GSM8836085_bulkG4CnT_3T3_rep2.bw"
  download_and_unzip "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8836nnn/GSM8836085/suppl/GSM8836085_bulkG4CnT_3T3_rep2.broadPeak.gz" \
           "${DATA}/GSE291468/GSM8836085_bulkG4CnT_3T3_rep2.broadPeak"

  # Tabix indexes for the fragment files (required by several scripts that
  # query fragments in specific regions; .tbi files are also shipped with
  # the repository)
  index_fragments "${DATA}"/GSE291468/GSM*/CellRanger/fragments.tsv.gz

  # Generate the PQS score-filtered subsets (PQS_scores.min40/min50.mm10.bed)
  # required by fig1/fig1.R and the revision scripts. Pure awk, no extra tools.
  if [[ ! -s "${DATA}/pqsfinder/PQS_scores.min40.mm10.bed" || \
        ! -s "${DATA}/pqsfinder/PQS_scores.min50.mm10.bed" ]]; then
    echo ""
    echo "[filter_pqs] generating PQS_scores.min40/min50.mm10.bed"
    bash "${DATA}/pqsfinder/filter_pqs.sh" || {
      echo "[warn] filter_pqs.sh failed - generate the subsets manually with:"
      echo "       bash data/pqsfinder/filter_pqs.sh"
    }
  fi
fi

# ======================================================================
# GSE163484: Bartosovic et al. mouse brain scRNA-seq (reference)
# ======================================================================
if [[ "$FIG" == "all" || "$FIG" == "--fig2" || "$FIG" == "--fig3" || "$FIG" == "--fig4" ]]; then
  echo ""
  echo "--- GSE163484: Bartosovic scRNA-seq ---"

  # The Seurat reference object `scRNA_Seq-mouse_brain.Rds` is NOT available
  # on GEO (GEO only provides raw count matrices). It is a processed Seurat
  # object provided directly by the study's authors, and it is SHIPPED with
  # this repository via Git LFS. If the file is missing or still an LFS
  # pointer, run 'git lfs pull' in the repository root.
  if [[ -f "${DATA}/GSE163484/scRNA_Seq-mouse_brain.Rds" ]] && ! is_lfs_pointer "${DATA}/GSE163484/scRNA_Seq-mouse_brain.Rds"; then
    echo "[ok] scRNA_Seq-mouse_brain.Rds present (shipped via Git LFS)"
  else
    echo "[warn] scRNA_Seq-mouse_brain.Rds missing or still a Git LFS pointer."
    echo "       Run 'git lfs pull' in the repository root to materialize it."
    echo "       (This author-provided file cannot be downloaded from GEO.)"
    FAILED+=("${DATA}/GSE163484/scRNA_Seq-mouse_brain.Rds")
  fi
fi

# ======================================================================
# Zeisel2018: Zeisel et al. mouse nervous system neuron scRNA-seq
# (Fig 4 neuron reference)
# ======================================================================
if [[ "$FIG" == "all" || "$FIG" == "--fig4" ]]; then
  echo ""
  echo "--- Zeisel2018: Zeisel neuron scRNA-seq ---"

  # The Seurat neuron object `scRNA_Seq-Zeisel_et_al-neuron.Rds` is NOT
  # available on GEO. It is a processed Seurat object built by
  # `data/Zeisel2018/processing_neuron_data.R` from the level-2 neuron loom
  # files (l2_neurons_*.agg.loom) of the Zeisel et al. (2018) mouse nervous
  # system atlas (Cell, DOI: 10.1016/j.cell.2018.06.021; http://mousebrain.org).
  # Those loom files are no longer publicly hosted (mousebrain.org is
  # offline), so the Rds was received directly from the study's authors.
  # It is SHIPPED with this repository via Git LFS. If the file is missing
  # or still an LFS pointer, run 'git lfs pull' in the repository root.
  if [[ -f "${DATA}/Zeisel2018/scRNA_Seq-Zeisel_et_al-neuron.Rds" ]] && ! is_lfs_pointer "${DATA}/Zeisel2018/scRNA_Seq-Zeisel_et_al-neuron.Rds"; then
    echo "[ok] scRNA_Seq-Zeisel_et_al-neuron.Rds present (shipped via Git LFS)"
  else
    echo "[warn] scRNA_Seq-Zeisel_et_al-neuron.Rds missing or still a Git LFS"
    echo "       pointer. Run 'git lfs pull' in the repository root to"
    echo "       materialize it (this author-provided file cannot be"
    echo "       downloaded from GEO)."
    FAILED+=("${DATA}/Zeisel2018/scRNA_Seq-Zeisel_et_al-neuron.Rds")
  fi
fi

# ======================================================================
# GSE198467: Brain snATAC-seq reference
# ======================================================================
if [[ "$FIG" == "all" || "$FIG" == "--fig4" || "$FIG" == "--revisions" ]]; then
  echo ""
  echo "--- GSE198467: Brain snATAC-seq ---"
  download_and_unzip "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198467/suppl/GSE198467_ATAC_Seurat_object_clustered_renamed.Rds.gz" \
           "${DATA}/GSE198467/GSE198467_ATAC_Seurat_object_clustered_renamed.Rds"

  # Raw ATAC fragments used to re-call consensus peaks with MACS2
  # (summit-containing narrowPeak), see data/GSE198467/reprocess_atac_peaks.sh
  download_and_unzip "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198467/suppl/GSE198467_ATAC_fragments.tsv.gz" \
           "${DATA}/GSE198467/GSE198467_ATAC_fragments.tsv.gz"
  index_fragments "${DATA}/GSE198467/GSE198467_ATAC_fragments.tsv.gz"

  # The summit-containing consensus narrowPeak and the fragment-coverage
  # RPGC bigWig are SHIPPED with the repository (the bigWig via Git LFS),
  # so they are ready after 'git lfs pull'. The two scripts below are the
  # original regeneration pipeline - they are only run here if the shipped
  # files are missing and the required tools are available.
  if [[ ! -s "${DATA}/GSE198467/GSE198467_atac_consensus_peaks.narrowPeak" ]]; then
    if command -v macs2 >/dev/null 2>&1; then
      echo ""
      echo "[reprocess_atac_peaks] re-calling consensus peaks with MACS2"
      bash "${DATA}/GSE198467/reprocess_atac_peaks.sh" || {
        echo "[warn] reprocess_atac_peaks.sh failed - rerun manually with:"
        echo "       bash data/GSE198467/reprocess_atac_peaks.sh"
      }
    else
      echo "[warn] GSE198467_atac_consensus_peaks.narrowPeak missing and macs2"
      echo "       not on PATH - regenerate with (after 'git lfs pull' it should"
      echo "       be present): bash data/GSE198467/reprocess_atac_peaks.sh"
    fi
  else
    echo "[ok] GSE198467_atac_consensus_peaks.narrowPeak present (shipped)"
  fi

  if [[ ! -s "${DATA}/GSE198467/GSE198467_ATAC_fragments.RPGC.bw" ]] || \
     is_lfs_pointer "${DATA}/GSE198467/GSE198467_ATAC_fragments.RPGC.bw"; then
    if command -v bamCoverage >/dev/null 2>&1; then
      echo ""
      echo "[make_atac_bigwig] building RPGC bigWig from fragments"
      bash "${DATA}/GSE198467/make_atac_bigwig.sh" || {
        echo "[warn] make_atac_bigwig.sh failed - rerun manually with:"
        echo "       bash data/GSE198467/make_atac_bigwig.sh"
      }
    else
      echo "[warn] GSE198467_ATAC_fragments.RPGC.bw missing (still a Git LFS"
      echo "       pointer?) - run 'git lfs pull' first, or regenerate with:"
      echo "       bash data/GSE198467/make_atac_bigwig.sh"
    fi
  else
    echo "[ok] GSE198467_ATAC_fragments.RPGC.bw present (shipped via Git LFS)"
  fi
fi

# ======================================================================
# Published bulk ATAC-seq references
# ======================================================================
if [[ "$FIG" == "all" || "$FIG" == "--revisions" ]]; then
  echo ""
  echo "--- Published bulk ATAC BigWigs (shipped via Git LFS) ---"

  # The bulk ATAC bigWigs and narrowPeak files are SHIPPED with this
  # repository via Git LFS and are NOT downloaded here, because the GEO
  # files do not match the data used for the published figures:
  #   - GSE149080 publishes GSM4661960 only as an mm9 bigWig
  #     (GSM4661960_ATAC-B2-H33WT.mm9.bw), while everything else in this
  #     repository is mm10. The shipped bigWig and narrowPeak are
  #     author-reprocessed mm10 versions (re-aligned and re-called with
  #     MACS2; GEO publishes no peak calls at all).
  #   - GSE211123 (GSM6451000) publishes an mm10 bigWig, but the shipped
  #     copy was re-processed/normalized by the same pipeline, so it also
  #     differs from the GEO file.
  # Downloading the GEO files instead would silently replace mm10 data
  # with wrong-genome (mm9) or differently processed signal.
  for f in "${DATA}/GSE149080/GSM4661960_ATAC_ESC_WT_batch2.rpgc.bw" \
           "${DATA}/GSE211123/GSM6451000_ATAC_3T3.rpgc.bw" \
           "${DATA}/GSE149080/GSM4661960_ESC_WT_batch2_peaks.narrowPeak" \
           "${DATA}/GSE211123/GSM6451000_3T3_control_ATAC_peaks.narrowPeak"; do
    if [[ -f "$f" ]] && ! is_lfs_pointer "$f"; then
      echo "[ok] $(basename "$f") present (shipped via Git LFS)"
    else
      echo "[warn] $(basename "$f") missing or still a Git LFS pointer."
      echo "       Run 'git lfs pull' in the repository root to materialize it."
      echo "       (These author-processed mm10 files cannot be downloaded from GEO.)"
      FAILED+=("$f")
    fi
  done

  echo ""
  echo "--- Published bulk 3T3 G4 CUT&Tag BigWigs ---"
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6729nnn/GSM6729104/suppl/GSM6729104_3T3-rep1_R1_val_1.bw" \
           "${DATA}/GSE217860/GSM6729104_3T3-rep1_R1_val_1.bw"
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6729nnn/GSM6729105/suppl/GSM6729105_3T3-rep2_R1_val_1.bw" \
           "${DATA}/GSE217860/GSM6729105_3T3-rep2_R1_val_1.bw"
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6729nnn/GSM6729106/suppl/GSM6729106_3T3-rep3_R1_val_1.bw" \
           "${DATA}/GSE217860/GSM6729106_3T3-rep3_R1_val_1.bw"

  echo ""
  echo "--- Chronis et al. mESC/MEF RNA-seq (GSE90894) ---"
  download "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE90nnn/GSE90894/suppl/GSE90894_RPKM_mRNAseq_table.xlsx" \
           "${DATA}/GSE90894/GSE90894_RPKM_mRNAseq_table.xlsx"
fi

# ======================================================================
# GENCODE mm10 annotation (Fig 1 feature distribution and revision F01)
# ======================================================================
if [[ "$FIG" == "all" || "$FIG" == "--fig1" || "$FIG" == "--revisions" ]]; then
  echo ""
  echo "--- GENCODE mm10 annotation ---"
  download "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz" \
           "${DATA}/genome/mm10_.annotation.gtf.gz"
fi

echo ""
echo "============================================"
if [[ "${#FAILED[@]}" -gt 0 ]]; then
  echo "  Download finished with ${#FAILED[@]} failed file(s):"
  for f in "${FAILED[@]}"; do
    echo "    - ${f}"
  done
  echo "  Fix the issues above and rerun this script (already-downloaded"
  echo "  files are skipped)."
  exit 1
fi
echo "  Download complete."
echo "  See README.md for preprocessing steps."
echo "============================================"
