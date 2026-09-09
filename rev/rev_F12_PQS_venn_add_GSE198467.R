#!/usr/bin/env Rscript
# ===========================================================================
# Add GSE198467 (GFP+ sorted scATAC) PQS overlap to existing analysis
# ===========================================================================

source("rev/paths.R")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(ggplot2)
})

# ----- Paths -----
OUT_DIR  <- file.path(OUT_ROOT, "PQS_venn_diagrams")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CANONICAL <- c(paste0("chr", 1:19), "chrX", "chrY")
TARGET_WIDTH <- 270

# ----- Load PQS -----
cat("Loading PQS BED...\n")
PQS_BED <- file.path(DATA, "pqsfinder/PQS_scores.mm10.bed")
pqs <- import(PQS_BED, format = "BED")
pqs <- pqs[as.character(seqnames(pqs)) %in% CANONICAL]
pqs <- pqs[!is.na(mcols(pqs)$score) & mcols(pqs)$score >= 20]
cat(sprintf("  PQS sites: %s\n", format(length(pqs), big.mark = ",")))

# ----- Overlap function -----
overlap_counts <- function(peaks_gr, pqs_gr, peaks_label, genome_bp = 2.65e9) {
  n_peaks  <- length(peaks_gr)
  n_pqs    <- length(pqs_gr)
  n_int    <- sum(IRanges::overlapsAny(peaks_gr, pqs_gr, ignore.strand = TRUE))
  n_pqs_in <- sum(IRanges::overlapsAny(pqs_gr, peaks_gr, ignore.strand = TRUE))
  
  pqs_reduced  <- GenomicRanges::reduce(pqs_gr, ignore.strand = TRUE)
  pqs_bp_total <- sum(as.numeric(width(pqs_reduced)))
  pqs_frac     <- pqs_bp_total / genome_bp
  med_peak_w   <- median(width(peaks_gr))
  pqs_frac_adj <- min(1, pqs_frac * (1 + med_peak_w / 30))
  exp_overlap  <- n_peaks * pqs_frac_adj
  fold_enrich  <- n_int / exp_overlap
  bin_p        <- binom.test(n_int, n_peaks, p = pqs_frac_adj, alternative = "greater")$p.value
  
  list(
    label             = peaks_label,
    n_peaks           = n_peaks,
    n_PQS_sites       = n_pqs,
    overlap           = n_int,
    pct_peaks_overlap = 100 * n_int / n_peaks,
    pct_PQS_overlap   = 100 * n_pqs_in / n_pqs,
    expected_overlap  = exp_overlap,
    fold_enrichment   = fold_enrich,
    binom_pvalue      = bin_p
  )
}

# ----- GSE198467: GFP+ sorted scATAC -----
cat("\nLoading GSE198467 GFP+ sorted scATAC peaks...\n")
GSE198467_RDS <- file.path(DATA, "GSE198467/GSE198467_ATAC_Seurat_object_clustered_renamed.Rds")
atac_gfp <- readRDS(GSE198467_RDS)

# Extract peak coordinates from Seurat object
atac_peak_names <- rownames(atac_gfp)
atac_peak_parts <- strsplit(atac_peak_names, "[-:]")
atac_peaks_gr <- GRanges(
  seqnames = sapply(atac_peak_parts, `[`, 1),
  ranges = IRanges(
    start = as.integer(sapply(atac_peak_parts, `[`, 2)),
    end = as.integer(sapply(atac_peak_parts, `[`, 3))
  )
)
atac_peaks_gr <- atac_peaks_gr[as.character(seqnames(atac_peaks_gr)) %in% CANONICAL]
cat(sprintf("  GSE198467 peaks (canonical): %s (median width %d bp)\n",
            format(length(atac_peaks_gr), big.mark = ","),
            median(width(atac_peaks_gr))))

# Resize to 270 bp for consistency
atac_resized <- resize(atac_peaks_gr, width = TARGET_WIDTH, fix = "center")
cat(sprintf("  Resized to %d bp: %s peaks\n",
            TARGET_WIDTH, format(length(atac_resized), big.mark = ",")))

oc_gse198467 <- overlap_counts(atac_resized, pqs, "GSE198467 GFP+ sorted scATAC")
cat(sprintf("  GSE198467 overlap: %.1f%% (%s / %s)  %.2fx  p=%g\n",
            oc_gse198467$pct_peaks_overlap,
            format(oc_gse198467$overlap, big.mark = ","),
            format(oc_gse198467$n_peaks, big.mark = ","),
            oc_gse198467$fold_enrichment,
            oc_gse198467$binom_pvalue))

# ----- Bar chart for GSE198467 -----
bar_gse198467 <- data.frame(
  cluster = "GFP+ sorted scATAC",
  pct     = oc_gse198467$pct_peaks_overlap,
  overlap = oc_gse198467$overlap,
  n_peaks = oc_gse198467$n_peaks,
  fold    = oc_gse198467$fold_enrichment,
  pval    = oc_gse198467$binom_pvalue
)
bar_gse198467$label <- sprintf("%.1f%%\n(%s / %s)",
                                bar_gse198467$pct,
                                format(bar_gse198467$overlap, big.mark = ","),
                                format(bar_gse198467$n_peaks, big.mark = ","))

p_gse198467 <- ggplot(bar_gse198467, aes(x = cluster, y = pct)) +
  geom_col(width = 0.6, fill = "#9ecae1", color = "grey30") +
  geom_text(aes(label = label), vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_y_continuous(limits = c(0, 105),
                     breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL,
       y = "% of scATAC peaks overlapping >=1 PQS site",
       title = "GSE198467 GFP+ sorted scATAC overlap with PQS",
       subtitle = sprintf("%.2fx enrichment, binomial p=%s",
                          oc_gse198467$fold_enrichment,
                          format.pval(oc_gse198467$binom_pvalue, digits = 2))) +
  theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 11))

ggsave(file.path(OUT_DIR, "bar_pct_peaks_overlap_PQS_scATAC_GSE198467_GFPsorted.pdf"), p_gse198467,
       width = 6, height = 6)
cat("  Saved bar_pct_peaks_overlap_PQS_scATAC_GSE198467_GFPsorted.pdf\n")

cat("\nDone.\n")
