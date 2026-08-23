#!/usr/bin/env Rscript
# ===========================================================================
# Figure 1 - deeptools heatmaps
#
# Compute deepTools computeMatrix + plotHeatmap for 6 bigwigs
# over the 3 peak categories from Panel F (cl0-only, shared, cl1-only).
#
# Requires the deeptools conda env:
#   /home/fr/fr_fr/fr_se1100/.conda/envs/deeptools
#
# Run from repo root:
#   Rscript fig1/fig1_deeptools.R
# ===========================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

ROOT <- normalizePath(getwd())
DATA <- file.path(ROOT, "data")
OUT  <- file.path(ROOT, "fig1", "outputs", "deeptools")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

DEEPTOOLS <- "/home/fr/fr_fr/fr_se1100/.conda/envs/deeptools/bin"

# --- Step 1: Generate 3 peak categories (same as Panel F) ---
cat("== Generating peak categories ==\n")
cl0 <- rtracklayer::import(file.path(DATA, "GSE291468/GSM8836088_cluster_0_peaks.narrowPeak"))
cl1 <- rtracklayer::import(file.path(DATA, "GSE291468/GSM8836088_cluster_1_peaks.narrowPeak"))

hits <- findOverlaps(cl0, cl1)
shared0 <- unique(queryHits(hits))
shared1 <- unique(subjectHits(hits))

cl0_only <- cl0[-shared0]
cl1_only <- cl1[-shared1]
both <- GenomicRanges::reduce(sort(c(cl0[shared0], cl1[shared1])))

cat("  cl0-only:", length(cl0_only), "\n")
cat("  shared:  ", length(both), "\n")
cat("  cl1-only:", length(cl1_only), "\n")

to_bed <- function(gr) {
  data.frame(
    chrom  = as.character(seqnames(gr)),
    start  = start(gr) - 1,
    end    = end(gr),
    name   = ".",
    score  = 0,
    strand = "."
  )
}

bed_cl0   <- file.path(OUT, "cl0_only.bed")
bed_shared <- file.path(OUT, "shared.bed")
bed_cl1   <- file.path(OUT, "cl1_only.bed")
write.table(to_bed(cl0_only),   bed_cl0,    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(to_bed(both),       bed_shared, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(to_bed(cl1_only),   bed_cl1,    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("  Wrote BED files to", OUT, "\n")

# --- Step 2: Run computeMatrix + plotHeatmap per bigwig ---
cat("\n== Running deepTools ==\n")

bw_info <- list(
  list(label = "cl0",    bw = file.path(DATA, "GSE291468/GSM8836088_cluster0_RPGC.bw")),
  list(label = "cl1",    bw = file.path(DATA, "GSE291468/GSM8836088_cluster1_RPGC.bw")),
  list(label = "mESC",   bw = file.path(DATA, "GSE291468/GSM8836082_bulkG4CnT_mESC_rep1.bw")),
  list(label = "MEF",    bw = file.path(DATA, "GSE291468/GSM8836084_bulkG4CnT_3T3_rep1.bw")),
  list(label = "ATAC",   bw = file.path(DATA, "GSE149080/GSM4661960_ATAC_ESC_WT_batch2.rpgc.bw")),
  # TODO: replace with new PQS_scores.mm10.bw when available
  list(label = "PQS",    bw = file.path(DATA, "pqsfinder/PQS_scores.mm10.bw"))
)

for (info in bw_info) {
  if (!file.exists(info$bw)) {
    cat("  SKIP", info$label, ": file not found\n")
    next
  }
  cat("  ", info$label, "\n")

  mat_file <- file.path(OUT, paste0("matrix_", info$label, ".mat.gz"))
  pdf_file <- file.path(OUT, paste0("heatmap_", info$label, ".pdf"))

  # computeMatrix
  cmd_mat <- sprintf(
    paste(
      "%s/computeMatrix reference-point",
      "-S '%s'",
      "-R '%s' '%s' '%s'",
      "-b 3000 -a 3000",
      "--skipZeros --missingDataAsZero",
      "--samplesLabel '%s'",
      "-o '%s'"
    ),
    DEEPTOOLS, info$bw, bed_cl0, bed_shared, bed_cl1,
    info$label, mat_file
  )
  system(cmd_mat, intern = TRUE)

  # Compute 90th percentile zmax from the matrix for color scale capping
  zmax_cmd <- sprintf(
    "%s/Rscript -e 'm <- read.table(gzfile(\"%s\"), skip=1, comment.char=\"@\"); cat(quantile(as.matrix(m), 0.90, na.rm=TRUE), \"\\n\")'",
    DEEPTOOLS, mat_file
  )
  zmax_val <- as.numeric(system(zmax_cmd, intern = TRUE)[1])
  zmax_val <- round(zmax_val, 1)
  cat("    zmax (p90):", zmax_val, "\n")

  # plotHeatmap
  cmd_plot <- sprintf(
    paste(
      "%s/plotHeatmap -m '%s'",
      "-out '%s'",
      "--refPointLabel 'peak'",
      "--heatmapHeight 14",
      "--whatToShow 'heatmap and colorbar'",
      "--colorMap 'Blues'",
      "--zMax %s",
      "-z 'cl0' 'shared' 'cl1'",
      "--yAxisLabel '' --xAxisLabel ''",
      "--legendLocation 'none'"
    ),
    DEEPTOOLS, mat_file, pdf_file, zmax_val
  )
  system(cmd_plot, intern = TRUE)
  cat("    Saved", basename(pdf_file), "\n")
}

cat("\n== Done ==\n")
cat("Outputs in:", OUT, "\n")
