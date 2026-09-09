# =============================================================
# Fig 3C rebuttal panels:
#   (1) PQS-overlapping G4 peaks per gene, AST vs non-AST, at three
#       stringent thresholds: 30, 40, 50 RPGC (width >= 100 bp).
#   (2) AST / non-AST peak count ratio vs threshold (robustness check).
#
# To run in RStudio: open this file, set DATA_DIR if needed, then Source.
# =============================================================


source("rev/paths.R")
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(bedscout)
  library(ggplot2)
})

# ----- Paths -----
inp <- function(fn) {
  hits <- list.files(DATA, pattern = fn, recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("File not found in data/: ", fn)
  hits[1]
}
OUT_DIR  <- file.path(OUT_ROOT, "fig3C_peak_validation")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# BW_AST and BW_NONAST are defined in rev/paths.R
PQS_BED   <- file.path(DATA, "pqsfinder/PQS_scores.mm10.bed")
GTF_PATH  <- file.path(GENOME_DIR, "mm10_.annotation.gtf.gz")

# ----- Parameters -----
PANEL_THRESHOLDS <- c(30, 40, 50)              # for the headline PQS-overlap panel
SWEEP_THRESHOLDS <- c(5, 10, 15, 20, 30, 40, 50)  # for the ratio sweep panel
MIN_WIDTH        <- 100                        # bp
FIG3C_GENES      <- c("Tnik", "Pitpnc1", "Pbx1", "Nwd1")
CANONICAL        <- c(paste0("chr", 1:19), "chrX", "chrY")

# ----- Peak caller -----
call_peaks <- function(bw_rle, chrom, gs, ge, threshold, min_width = MIN_WIDTH) {
  if (!(chrom %in% names(bw_rle))) return(GRanges())
  cov <- bw_rle[[chrom]]
  if (ge > length(cov)) ge <- length(cov)
  v <- IRanges::slice(cov[gs:ge], lower = threshold, includeLower = TRUE)
  if (length(v) == 0) return(GRanges())
  rng <- ranges(v)
  rng <- rng[width(rng) >= min_width]
  if (length(rng) == 0) return(GRanges())
  rng <- IRanges::shift(rng, gs - 1L)
  GRanges(seqnames = chrom, ranges = rng)
}

# ----- Load data -----
cat("Loading BigWigs, PQS, gene coords...\n")
bw_ast    <- import(BW_AST,    as = "RleList")
bw_nonast <- import(BW_NONAST, as = "RleList")
pqs <- import(PQS_BED, format = "BED")
pqs <- pqs[as.character(seqnames(pqs)) %in% CANONICAL]
pqs <- pqs[!is.na(mcols(pqs)$score) & mcols(pqs)$score >= 50]

gtf <- rtracklayer::import(GTF_PATH)
genes <- gtf[gtf$type == "gene"]
gene_coord <- as.data.table(genes)[gene_name %in% FIG3C_GENES,
                                    .(gene_name, chrom = as.character(seqnames),
                                      start, end)]
gene_coord <- gene_coord[, .(chrom = data.table::first(chrom),
                              start = min(start), end = max(end)),
                         by = gene_name]
print(gene_coord)

# =============================================================
# Panel 1: PQS-overlapping peak counts per gene at 30, 40, 50 RPGC
# =============================================================
cat(sprintf("\nCalling peaks + PQS overlap at thresholds %s RPGC...\n",
            paste(PANEL_THRESHOLDS, collapse = ", ")))

panel_results <- data.table()
for (thr in PANEL_THRESHOLDS) {
  for (i in seq_len(nrow(gene_coord))) {
    g <- gene_coord[i]
    for (info in list(list(label = "AST",     bw = bw_ast),
                       list(label = "non-AST", bw = bw_nonast))) {
      peaks <- call_peaks(info$bw, g$chrom, g$start, g$end, threshold = thr)
      n_peaks <- length(peaks)
      n_PQS   <- if (n_peaks > 0) sum(overlapsAny(peaks, pqs, ignore.strand = TRUE)) else 0
      panel_results <- rbind(panel_results, data.table(
        threshold = thr,
        gene      = g$gene_name,
        cluster   = info$label,
        n_peaks   = n_peaks,
        n_PQS     = n_PQS,
        pct_PQS   = if (n_peaks > 0) round(100 * n_PQS / n_peaks, 1) else 0
      ))
    }
  }
}
panel_results$gene      <- factor(panel_results$gene,    levels = FIG3C_GENES)
panel_results$cluster   <- factor(panel_results$cluster, levels = c("AST", "non-AST"))
panel_results$threshold <- factor(panel_results$threshold,
                                   levels = PANEL_THRESHOLDS,
                                   labels = paste0(">= ", PANEL_THRESHOLDS, " RPGC"))
fwrite(panel_results, file.path(OUT_DIR, "fig3C_PQS_overlap_thresholds.csv"))
cat("\nPanel results:\n"); print(panel_results)

cat("\nDrawing PQS-overlap panel (faceted by threshold)...\n")
p1 <- ggplot(panel_results, aes(x = gene, y = n_PQS, fill = cluster)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d (%.0f%%)", n_PQS, pct_PQS)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 2.8, fontface = "bold") +
  facet_wrap(~ threshold, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c(AST = "#d62728", `non-AST` = "#7f7f7f"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = NULL, y = "Number of PQS-overlapping peaks",
       title = "PQS-overlapping G4 peaks per gene by signal threshold",
       subtitle = sprintf("AST vs non-AST, three stringency levels (min width = %d bp)",
                          MIN_WIDTH)) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "top")

print(p1)
ggsave(file.path(OUT_DIR, "fig3C_PQS_overlap_thresholds.pdf"), p1,
       width = 13, height = 5)

# =============================================================
# Panel 2: AST / non-AST peak count ratio vs threshold
# =============================================================
cat(sprintf("\nSweeping thresholds %s for the ratio plot...\n",
            paste(SWEEP_THRESHOLDS, collapse = ", ")))

ratio_results <- data.table()
for (thr in SWEEP_THRESHOLDS) {
  for (i in seq_len(nrow(gene_coord))) {
    g <- gene_coord[i]
    n_ast <- length(call_peaks(bw_ast,    g$chrom, g$start, g$end, threshold = thr))
    n_non <- length(call_peaks(bw_nonast, g$chrom, g$start, g$end, threshold = thr))
    ratio_results <- rbind(ratio_results, data.table(
      threshold = thr,
      gene      = g$gene_name,
      n_AST     = n_ast,
      n_nonAST  = n_non,
      ratio     = if (n_non > 0) n_ast / n_non else NA_real_
    ))
  }
}
ratio_results$gene <- factor(ratio_results$gene, levels = FIG3C_GENES)
fwrite(ratio_results, file.path(OUT_DIR, "fig3C_threshold_ratio.csv"))
cat("\nThreshold sweep:\n"); print(ratio_results)

# (Threshold-ratio line plot removed — not used in the slides. The
#  fig3C_threshold_ratio.csv supporting table above is still written.)

cat("\nDone. Outputs:\n",
    "  ", file.path(OUT_DIR, "fig3C_PQS_overlap_thresholds.csv"), "\n",
    "  ", file.path(OUT_DIR, "fig3C_PQS_overlap_thresholds.pdf"), "\n",
    "  ", file.path(OUT_DIR, "fig3C_PQS_overlap_thresholds.pdf"), "\n",
    "  ", file.path(OUT_DIR, "fig3C_threshold_ratio.csv"),         "\n", sep = "")
