#!/usr/bin/env Rscript
# =============================================================
# F26: Correlation between scG4 clusters and matched bulk G4
# CUT&Tag tracks at differential G4 peaks.
#
# The original source script was lost, but its 992-region BED and
# correlation outputs survived. This script regenerates that legacy
# RNA-seq-matched set and an alternate set containing all differential
# peaks from F24 without the RNA-seq/promoter restriction.
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

F24_CATEGORIES_RNASEQ <- file.path(
  OUT_ROOT, "RNAseq_lfc_violin", "F24_peak_categories_rnaseq.csv"
)
F24_CATEGORIES_ALL <- file.path(
  OUT_ROOT, "RNAseq_lfc_violin", "F24_peak_categories.csv"
)
LEGACY_BED <- file.path(OUT_DIR, "F26_marker_regions.mm10.bed")
ALL_DIFF_BED <- file.path(OUT_DIR, "F26_marker_regions_all_diff.mm10.bed")
EXPECTED_LEGACY_N <- 992L
EXPECTED_ALL_DIFF_N <- 10713L

BIGWIGS <- c(
  "MEF bulk" = file.path(GSE291468, "GSM8836084_bulkG4CnT_3T3_rep1.bw"),
  "cluster0" = file.path(GSE291468, "GSM8836088_cluster0_RPGC.bw"),
  "cluster1" = file.path(GSE291468, "GSM8836088_cluster1_RPGC.bw"),
  "mESC bulk" = file.path(GSE291468, "GSM8836082_bulkG4CnT_mESC_rep1.bw")
)

required_inputs <- c(F24_CATEGORIES_RNASEQ, F24_CATEGORIES_ALL, BIGWIGS)
missing_inputs <- required_inputs[
  !file.exists(required_inputs)
]
if (length(missing_inputs)) {
  stop("Missing required input(s):\n  ", paste(missing_inputs, collapse = "\n  "))
}

load_diff_peaks <- function(path, expected_n, label) {
  peak_table <- fread(path)
  required_columns <- c("seqnames", "start", "end", "diff")
  if (!all(required_columns %in% names(peak_table))) {
    stop("F24 table is missing columns: ",
         paste(setdiff(required_columns, names(peak_table)), collapse = ", "))
  }

  peak_table <- peak_table[diff %chin% c("MEF", "ESC")]
  setorder(peak_table, seqnames, start, end)
  peak_table <- unique(peak_table, by = c("seqnames", "start", "end"))
  if (nrow(peak_table) != expected_n) {
    stop("Expected ", expected_n, " ", label, " peaks, found ",
         nrow(peak_table), ". Regenerate F24 or inspect its thresholds.")
  }

  sort(GRanges(
    seqnames = peak_table$seqnames,
    ranges = IRanges(start = peak_table$start, end = peak_table$end),
    diff = peak_table$diff
  ))
}

cat("Loading F24 differential peak categories...\n")
legacy_gr <- load_diff_peaks(
  F24_CATEGORIES_RNASEQ, EXPECTED_LEGACY_N, "RNA-seq-matched differential"
)
all_diff_gr <- load_diff_peaks(
  F24_CATEGORIES_ALL, EXPECTED_ALL_DIFF_N, "all differential"
)

# Confirm that regeneration still recovers the surviving canonical BED.
if (file.exists(LEGACY_BED)) {
  previous_gr <- sort(rtracklayer::import(LEGACY_BED))
  if (!identical(
    as.character(legacy_gr),
    as.character(previous_gr)
  )) {
    stop("Regenerated F26 regions do not match ", LEGACY_BED)
  }
}
rtracklayer::export.bed(legacy_gr, LEGACY_BED)
rtracklayer::export.bed(all_diff_gr, ALL_DIFF_BED)
cat(sprintf(
  "  Legacy: %d regions (%d MEF, %d ESC)\n",
  length(legacy_gr), sum(legacy_gr$diff == "MEF"), sum(legacy_gr$diff == "ESC")
))
cat(sprintf(
  "  All differential: %d regions (%d MEF, %d ESC)\n",
  length(all_diff_gr), sum(all_diff_gr$diff == "MEF"),
  sum(all_diff_gr$diff == "ESC")
))

score_tracks <- function(regions, label) {
  cat("Scoring ", label, " with bw_loci...\n", sep = "")
  signal_gr <- bw_loci(
    BIGWIGS,
    loci = regions,
    labels = names(BIGWIGS),
    default_na = 0
  )
  signal_table <- as.data.frame(mcols(signal_gr))[, make.names(names(BIGWIGS))]
  names(signal_table) <- names(BIGWIGS)
  as.matrix(signal_table)
}

write_correlation_csv <- function(cor_matrix, filename) {
  output <- as.data.table(cor_matrix, keep.rownames = "track")
  fwrite(output, file.path(OUT_DIR, filename))
}

draw_correlation_heatmap <- function(cor_matrix, method, filename, scope, n_peaks) {
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
      "%s correlation (%s, n = %s)",
      tools::toTitleCase(method), scope, format(n_peaks, big.mark = ",")
    ),
    rect_gp = grid::gpar(col = "black", lwd = 1),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.text(
        sprintf("%.2f", cor_matrix[i, j]),
        x, y,
        gp = grid::gpar(fontsize = 12)
      )
    },
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    heatmap_width = grid::unit(8, "cm"),
    heatmap_height = grid::unit(8, "cm"),
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

analyses <- list(
  legacy = list(
    regions = legacy_gr,
    suffix = "",
    scope = "RNA-seq-matched differential G4 peaks"
  ),
  all_diff = list(
    regions = all_diff_gr,
    suffix = "_all_diff",
    scope = "all differential G4 peaks"
  )
)

all_correlations <- list()
for (analysis_name in names(analyses)) {
  analysis <- analyses[[analysis_name]]
  signal_matrix <- score_tracks(analysis$regions, analysis$scope)
  correlations <- list(
    spearman = cor(signal_matrix, method = "spearman", use = "complete.obs"),
    pearson = cor(signal_matrix, method = "pearson", use = "complete.obs")
  )
  all_correlations[[analysis_name]] <- correlations

  for (method in names(correlations)) {
    cor_matrix <- correlations[[method]]
    suffix <- analysis$suffix
    write_correlation_csv(
      cor_matrix, paste0("F26_", method, "_correlation", suffix, ".csv")
    )
    draw_correlation_heatmap(
      cor_matrix,
      method,
      paste0("panel_F26_", method, "_correlation", suffix),
      analysis$scope,
      length(analysis$regions)
    )
  }
}

cat("\nCorrelation matrices:\n")
print(lapply(all_correlations, function(x) lapply(x, round, digits = 3)))
cat("\nWrote F26 outputs to ", OUT_DIR, "\n", sep = "")
