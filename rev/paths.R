# rev/paths.R — Shared path definitions for revision scripts
# Source this at the top of each revision script:
#   source("rev/paths.R")
#
# Run all scripts from the repo root directory.

ROOT <- normalizePath(getwd())
DATA <- file.path(ROOT, "data")
RESULTS <- file.path(ROOT, "results")

# Data directories
GENOME_DIR   <- file.path(DATA, "genome")
GSE291468    <- file.path(DATA, "GSE291468")
CELLRANGER   <- GSE291468
GSE163484    <- file.path(DATA, "GSE163484")
GSE198467    <- file.path(DATA, "GSE198467")
BULK_PEAKS   <- GSE291468
GSE149080    <- file.path(DATA, "GSE149080")
GSE211123    <- file.path(DATA, "GSE211123")
GSE217860    <- file.path(DATA, "GSE217860")
GSE90894     <- file.path(DATA, "GSE90894")
cCRE_DIR     <- file.path(DATA, "cCRE")

# Common reference files
CHROM_SIZES  <- file.path(GENOME_DIR, "mm10.chrom.sizes.txt")
GTF_PATH     <- file.path(GENOME_DIR, "mm10_.annotation.gtf.gz")
TPM_PATH     <- file.path(GENOME_DIR, "gene_level_TPM.tsv")
EG4_PATH     <- file.path(GENOME_DIR, "Mouse_eG4.txt")
REFGENE_TSS  <- file.path(GENOME_DIR, "refGene.tss.1kb.mm10.bed")
PQS_BED      <- file.path(DATA, "pqsfinder/PQS_scores.mm10.bed")
PQS_BW       <- file.path(DATA, "pqsfinder/PQS_scores.mm10.bw")

# Chronis et al. mESC/MEF RNA-seq RPKM table (GSE90894)
RNASEQ_XLSX  <- file.path(GSE90894, "GSE90894_RPKM_mRNAseq_table.xlsx")

# GEO supplementary signal tracks (always prefer GEO versions)
BW_CL0       <- file.path(GSE291468, "GSM8836087_cluster0_RPGC.bw")
BW_CL1       <- file.path(GSE291468, "GSM8836087_cluster1_RPGC.bw")
BW_AST       <- file.path(GSE291468, "GSM8836086_Predicted_Astrocytes_RPGC.bw")
BW_NONAST    <- file.path(GSE291468, "GSM8836086_Predicted_non-Astrocytes_RPGC.bw")
BW_GFPPOS    <- file.path(GSE291468, "GSM8836086_GFP_sorted_mousebrain.rpgc.bw")

# GEO mESC_MEF supplementary objects
SEURAT_MESC_MEF <- file.path(GSE291468, "GSM8836088_mESCMEF_Seurat_object.Rds")
# Cluster peaks: GEO files in data/GSE291468/ — use individual paths below
# or file.path(GSE291468, "GSM8836088_cluster_X_peaks.narrowPeak")
CLUSTER_PEAKS_0 <- file.path(GSE291468, "GSM8836088_cluster_0_peaks.narrowPeak")
CLUSTER_PEAKS_1 <- file.path(GSE291468, "GSM8836088_cluster_1_peaks.narrowPeak")

# Output directories
OUT_ROOT     <- file.path(ROOT, "rev", "outputs")
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

# Canonical chromosomes
CANONICAL <- paste0("chr", c(1:19, "X", "Y"))

# scBridge outputs (de-novo, not on GEO)
RESULTS_SCBRIDGE <- file.path(RESULTS, "scBridge/output")
