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
# Usage:
#   cd data
#   bash download_data.sh              # download everything
#   bash download_data.sh --fig1       # only Fig 1 data
#   bash download_data.sh --fig2       # only Fig 2 data
#   bash download_data.sh --fig4       # only Fig 4 data
#   bash download_data.sh --revisions  # only revision data
#
# Requirements: curl, gzip
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${ROOT}"

download() {
  local url="$1" dest="$2"
  mkdir -p "$(dirname "${dest}")"
  if [[ "${FORCE_DOWNLOAD:-0}" != "1" && -s "${dest}" ]]; then
    echo "[skip] $(basename "${dest}") already exists"
    return
  fi
  echo "[download] $(basename "${dest}")"
  curl -fL --retry 5 --retry-delay 2 --continue-at - -o "${dest}.part" "${url}"
  mv "${dest}.part" "${dest}"
}

download_gzip_as_plain() {
  local url="$1" dest="$2"
  local compressed="${dest}.gz.part"
  mkdir -p "$(dirname "${dest}")"
  if [[ "${FORCE_DOWNLOAD:-0}" != "1" && -s "${dest}" ]]; then
    echo "[skip] $(basename "${dest}") already exists"
    return
  fi
  echo "[download] $(basename "${dest}")"
  curl -fL --retry 5 --retry-delay 2 --continue-at - -o "${compressed}" "${url}"
  gzip -dc "${compressed}" > "${dest}.part"
  rm -f "${compressed}"
  mv "${dest}.part" "${dest}"
}

download_and_unzip() {
  local url="$1" dest="$2"
  local compressed="${dest}.gz.part"
  mkdir -p "$(dirname "${dest}")"
  if [[ "${FORCE_DOWNLOAD:-0}" != "1" && -s "${dest}" ]]; then
    echo "[skip] $(basename "${dest}") already exists"
    return
  fi
  echo "[download] $(basename "${dest}")"
  curl -fL --retry 5 --retry-delay 2 --continue-at - -o "${compressed}" "${url}"
  mv "${compressed}" "${dest}.gz"
  gzip -d "${dest}.gz"
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
fi

# ======================================================================
# GSE163484: Bartosovic et al. mouse brain scRNA-seq (reference)
# ======================================================================
if [[ "$FIG" == "all" || "$FIG" == "--fig2" || "$FIG" == "--fig3" || "$FIG" == "--fig4" ]]; then
  echo ""
  echo "--- GSE163484: Bartosovic scRNA-seq ---"
  download "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163484/suppl/GSE163484_Bartosovic_et_al_counts.csv.gz" \
           "${DATA}/GSE163484/GSE163484_Bartosovic_et_al_counts.csv.gz"
  download "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163484/suppl/GSE163484_Bartosovic_et_al_annot.csv.gz" \
           "${DATA}/GSE163484/GSE163484_Bartosovic_et_al_annot.csv.gz"

  echo ""
  echo "[note] Zeisel et al. neuron data: download loom files from http://mousebrain.org/"
  echo "       Place at data/Zeisel_et_al/neuron_scRNA_Seq/l2_neurons_*.agg.loom"
fi

# ======================================================================
# GSE198467: Brain snATAC-seq reference
# ======================================================================
if [[ "$FIG" == "all" || "$FIG" == "--fig4" || "$FIG" == "--revisions" ]]; then
  echo ""
  echo "--- GSE198467: Brain snATAC-seq ---"
  download_and_unzip "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198467/suppl/GSE198467_ATAC_Seurat_object_clustered_renamed.Rds.gz" \
           "${DATA}/GSE198467/GSE198467_ATAC_Seurat_object_clustered_renamed.Rds"
fi

# ======================================================================
# Published bulk ATAC-seq references
# ======================================================================
if [[ "$FIG" == "all" || "$FIG" == "--revisions" ]]; then
  echo ""
  echo "--- Published bulk ATAC BigWigs ---"
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4661nnn/GSM4661960/suppl/GSM4661960_ATAC-B2-H33WT.mm9.bw" \
           "${DATA}/GSE149080/GSM4661960_ATAC_ESC_WT_batch2.rpgc.bw"
  download "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6451nnn/GSM6451000/suppl/GSM6451000_3T3-Control-ATAC.bw" \
           "${DATA}/GSE211123/GSM6451000_ATAC_3T3.rpgc.bw"

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

echo ""
echo "============================================"
echo "  Download complete."
echo "  See README.md for preprocessing steps."
echo "============================================"
