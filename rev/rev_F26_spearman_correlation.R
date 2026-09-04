#!/usr/bin/env Rscript
# =============================================================
# F26: Correlation between scG4 clusters and matched bulk G4
# CUT&Tag tracks at RNA-seq-matched differential G4 peaks.
#
# The original source script was lost, but its 992-region BED and
# correlation outputs survived. This script regenerates that peak
# set from F24 and writes both the legacy and explicit all-diff
# output names.
#
# Run from the repository root:
#   Rscript rev/rev_F26_spearman_correlation.R
# =============================================================

source("rev/paths.R")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(wigglescout)
  library(ComplexHeatmap)
  library(circlize)
  library(data.table)
})

OUT_DIR <- file.path(OUT_ROOT, "S1_marker_correlation")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

F24_CATEGORIES <- file.path(
  OUT_ROOT, "RNAseq_lfc_violin", "F24_peak_categories_rnaseq.csv"
)
MARKER_BED <- file.path(OUT_DIR, "F26_marker_regions.mm10.bed")
EXPECTED_N <- 992L

BIGWIGS <- c(
  "MEF bulk" = file.path(GSE291468, "GSM8836084_bulkG4CnT_3T3_rep1.bw"),
  "cluster0" = file.path(GSE291468, "GSM8836088_cluster0_RPGC.bw"),
  "cluster1" = file.path(GSE291468, "GSM8836088_cluster1_RPGC.bw"),
  "mESC bulk" = file.path(GSE291468, "GSM8836082_bulkG4CnT_mESC_rep1.bw")
)

missing_inputs <- c(F24_CATEGORIES, BIGWIGS)[
  !file.exists(c(F24_CATEGORIES, BIGWIGS))
]
if (length(missing_inputs)) {
  stop("Missing required input(s):\n  ", paste(missing_inputs, collapse = "\n  "))
}

cat("Loading F24 RNA-seq-matched differential peak categories...\n")
peak_table <- fread(F24_CATEGORIES)
required_columns <- c("seqnames", "start", "end", "diff")
if (!all(required_columns %in% names(peak_table))) {
  stop("F24 table is missing columns: ",
       paste(setdiff(required_columns, names(peak_table)), collapse = ", "))
}

all_diff <- peak_table[diff %chin% c("MEF", "ESC")]
setorder(all_diff, seqnames, start, end)
all_diff <- unique(all_diff, by = c("seqnames", "start", "end"))

if (nrow(all_diff) != EXPECTED_N) {
  stop("Expected ", EXPECTED_N, " differential peaks, found ", nrow(all_diff),
       ". Regenerate F24 before running F26 or inspect its thresholds.")
}

marker_gr <- GRanges(
  seqnames = all_diff$seqnames,
  ranges = IRanges(start = all_diff$start, end = all_diff$end),
  diff = all_diff$diff
)
marker_gr <- sort(marker_gr)

# Confirm that regeneration still recovers the surviving canonical BED.
if (file.exists(MARKER_BED)) {
  previous_gr <- sort(rtracklayer::import(MARKER_BED))
  if (!identical(
    as.character(marker_gr),
    as.character(previous_gr)
  )) {
    stop("Regenerated F26 regions do not match ", MARKER_BED)
  }
}
rtracklayer::export.bed(marker_gr, MARKER_BED)
cat(sprintf(
  "  Regenerated %d regions (%d MEF, %d ESC)\n",
  length(marker_gr), sum(marker_gr$diff == "MEF"), sum(marker_gr$diff == "ESC")
))

cat("Scoring four G4 CUT&Tag tracks with bw_loci...\n")
signal_gr <- bw_loci(
  BIGWIGS,
  loci = marker_gr,
  labels = names(BIGWIGS),
  default_na = 0
)
signal_table <- as.data.frame(mcols(signal_gr))[, make.names(names(BIGWIGS))]
names(signal_table) <- names(BIGWIGS)
signal_matrix <- as.matrix(signal_table)

write_correlation_csv <- function(cor_matrix, filename) {
  output <- as.data.table(cor_matrix, keep.rownames = "track")
  fwrite(output, file.path(OUT_DIR, filename))
}

draw_correlation_heatmap <- function(cor_matrix, method, filename) {
  label <- if (method == "spearman") "Spearman rho" else "Pearson r"
  color_function <- colorRamp2(
    c(-1, 0, 1),
    c("#2166ac", "white", "#b2182b")
  )
  heatmap <- Heatmap(
    cor_matrix,
    name = label,
    col = color_function,
    column_title = sprintf(
      "%s correlation (all differential G4 peaks, n = %d)",
      tools::toTitleCase(method), nrow(signal_matrix)
    ),
    rect_gp = grid::gpar(col = "black", lwd = 1),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.text(
        sprintf("%.2f", cor_matrix[i, j]),
        x, y,
        gp = grid::gpar(fontsize = 14)
      )
    },
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    heatmap_width = grid::unit(6, "cm"),
    heatmap_height = grid::unit(6, "cm"),
    row_names_gp = grid::gpar(fontsize = 12),
    column_names_gp = grid::gpar(fontsize = 12)
  )

  pdf(file.path(OUT_DIR, paste0(filename, ".pdf")), width = 7.5, height = 7.5)
  draw(heatmap)
  dev.off()

  png(
    file.path(OUT_DIR, paste0(filename, ".png")),
    width = 2250, height = 2250, res = 300
  )
  draw(heatmap)
  dev.off()
}

cat("Computing correlation matrices...\n")
correlations <- list(
  spearman = cor(signal_matrix, method = "spearman", use = "complete.obs"),
  pearson = cor(signal_matrix, method = "pearson", use = "complete.obs")
)

for (method in names(correlations)) {
  cor_matrix <- correlations[[method]]

  # Legacy filenames retained for the existing F26 panels.
  write_correlation_csv(cor_matrix, paste0("F26_", method, "_correlation.csv"))
  draw_correlation_heatmap(
    cor_matrix, method, paste0("panel_F26_", method, "_correlation")
  )

  # Explicit alternate filenames requested for the complete 992-peak set.
  write_correlation_csv(
    cor_matrix, paste0("F26_", method, "_correlation_all_diff.csv")
  )
  draw_correlation_heatmap(
    cor_matrix, method, paste0("panel_F26_", method, "_correlation_all_diff")
  )
}

cat("\nCorrelation matrices:\n")
print(lapply(correlations, round, digits = 3))
cat("\nWrote F26 outputs to ", OUT_DIR, "\n", sep = "")
