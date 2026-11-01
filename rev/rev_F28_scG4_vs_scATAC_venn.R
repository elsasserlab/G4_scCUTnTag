#!/usr/bin/env Rscript
# ===========================================================================
# Figure F28: scG4 vs scATAC peak overlap
#
# Tests whether scG4 signal reflects chromatin accessibility by comparing
# GFP+ scG4 (GSM8836086) with GFP+ scATAC (GSE198467) from same cell population.
#
# scG4 peaks are the merged UNION of the per-cluster MACS2 narrowPeaks
# (GSM8836086_GFPpos_cluster_{0..3}_peaks.narrowPeak), reduced and
# restricted to canonical chromosomes (median ~0.5 kb), NOT the pre-merged
# Seurat rownames (~0.9 kb). Consistent with F12 and F30.
#
# scATAC peaks are the MACS2 re-called consensus narrowPeak
# (data/GSE198467/GSE198467_atac_consensus_peaks.narrowPeak), not the broad
# author-provided windows (median ~3.7 kb) stored in the Seurat rownames.
#
# Output (github/rev/outputs/scG4_vs_scATAC/):
#   - F28_peak_overlap_venn.pdf
#   - F28_peak_overlap_stats.csv
# ===========================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(ggplotify)
  library(eulerr)
  library(data.table)
})

# ----- Paths -----
ROOT <- normalizePath(getwd())
DATA <- file.path(ROOT, "data")
OUT_DIR <- file.path(ROOT, "rev/outputs/scG4_vs_scATAC")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# scG4 GFP+ data: merged union of the per-cluster MACS2 narrowPeaks
# (consistent with F12/F30; NOT the pre-merged Seurat rownames)
G4_PEAK_CLUST <- file.path(DATA, "GSE291468",
                           paste0("GSM8836086_GFPpos_cluster_", 0:3, "_peaks.narrowPeak"))
for (f in G4_PEAK_CLUST) stopifnot(file.exists(f))

# scATAC data: MACS2 re-called consensus narrowPeak
ATAC_NARROWPEAK <- file.path(DATA, "GSE198467/GSE198467_atac_consensus_peaks.narrowPeak")

cat("=== Loading data ===\n")

# Merge GFP+ cluster narrowPeaks into a union (canonical only)
cat("Merging GFP+ cluster narrowPeaks into a union...\n")
g4_clusters <- lapply(G4_PEAK_CLUST, function(f) import(f, format = "narrowPeak"))
for (i in seq_along(g4_clusters)) {
  cat(sprintf("  GFP+ cluster %d peaks: %s\n", i - 1,
              format(length(g4_clusters[[i]]), big.mark = ",")))
}
CANONICAL <- c(paste0("chr", 1:19), "chrX", "chrY")
g4_peaks_gr <- reduce(do.call(c, g4_clusters))
g4_peaks_gr <- g4_peaks_gr[as.character(seqnames(g4_peaks_gr)) %in% CANONICAL]
cat(sprintf("  scG4 peaks (canonical union): %s (median width %d bp)\n",
            format(length(g4_peaks_gr), big.mark = ","),
            median(width(g4_peaks_gr))))

# scATAC peak coordinates: use the MACS2 re-called consensus narrowPeak
# (the Seurat rownames are broad windows with no summit offsets; see
# data/GSE198467/reprocess_atac_peaks.sh).
cat("Loading re-called scATAC consensus narrowPeak...\n")
ATAC_NARROWPEAK <- file.path(DATA, "GSE198467/GSE198467_atac_consensus_peaks.narrowPeak")
stopifnot(file.exists(ATAC_NARROWPEAK))
atac_peaks_gr <- import(ATAC_NARROWPEAK, format = "narrowPeak")
# Filter to canonical chromosomes
CANONICAL <- c(paste0("chr", 1:19), "chrX", "chrY")
canonical_mask <- as.character(seqnames(atac_peaks_gr)) %in% CANONICAL
atac_peaks_gr <- atac_peaks_gr[canonical_mask]
cat(sprintf("  Total scATAC peaks (canonical): %s (median width %d bp)\n",
            format(length(atac_peaks_gr), big.mark = ","),
            median(width(atac_peaks_gr))))

# ===========================================================================
# Peak Overlap Analysis (using reduce for symmetric Venn diagram)
# ===========================================================================
cat("\n=== Peak Overlap Analysis ===\n")

# Combine all peaks and reduce to non-overlapping regions
cat("Merging overlapping peaks with reduce()...\n")
all_peaks <- c(g4_peaks_gr, atac_peaks_gr)
all_peaks_reduced <- reduce(all_peaks, ignore.strand = TRUE)
cat(sprintf("  Original peaks: %s (scG4) + %s (scATAC) = %s total\n",
             format(length(g4_peaks_gr), big.mark = ","),
             format(length(atac_peaks_gr), big.mark = ","),
             format(length(all_peaks), big.mark = ",")))
cat(sprintf("  After reduce(): %s unique regions\n",
             format(length(all_peaks_reduced), big.mark = ",")))

# Count signal in each reduced region
cat("Counting signal per reduced region...\n")
g4_hits <- countOverlaps(all_peaks_reduced, g4_peaks_gr, ignore.strand = TRUE) > 0
atac_hits <- countOverlaps(all_peaks_reduced, atac_peaks_gr, ignore.strand = TRUE) > 0

n_g4_only <- sum(g4_hits & !atac_hits)
n_atac_only <- sum(!g4_hits & atac_hits)
n_overlap <- sum(g4_hits & atac_hits)
n_g4_total <- n_g4_only + n_overlap
n_atac_total <- n_atac_only + n_overlap

cat(sprintf("  scG4 peaks: %s\n", format(n_g4_total, big.mark = ",")))
cat(sprintf("  scATAC peaks: %s\n", format(n_atac_total, big.mark = ",")))
cat(sprintf("  Overlap: %s (%.1f%% of scG4, %.1f%% of scATAC)\n",
             format(n_overlap, big.mark = ","),
             100 * n_overlap / n_g4_total,
             100 * n_overlap / n_atac_total))
cat(sprintf("  scG4-only: %s (%.1f%%)\n", format(n_g4_only, big.mark = ","), 100 * n_g4_only / n_g4_total))
cat(sprintf("  scATAC-only: %s (%.1f%%)\n", format(n_atac_only, big.mark = ","), 100 * n_atac_only / n_atac_total))

# Create Euler diagram using eulerr (proportional to actual counts)
# Show region-specific counts within each area
cat("\nCreating Euler diagram...\n")
fit <- euler(c("scG4" = n_g4_only, "scATAC" = n_atac_only, "scG4&scATAC" = n_overlap))

# plot.euler() returns a gTree grob (not a ggplot), so wrap it with
# ggplotify before adding ggplot annotations.
p_venn <- ggplotify::as.ggplot(
  plot(fit,
       fills = list(fill = c("#9ecae1", "#fc9272")),
       edges = list(col = "black", lwd = 1.5),
       labels = FALSE,  # Remove set labels (use totals outside instead)
       quantities = list(font = 2, cex = 2.0, digits = 0))  # Show counts in regions
)

# Add original peak counts as text outside the circles (panel coords are 0-1)
orig_g4 <- length(g4_peaks_gr)
orig_atac <- length(atac_peaks_gr)
p_venn <- p_venn +
  annotate("text", x = 0.03, y = 0.75, label = sprintf("scG4\n%s", format(orig_g4, big.mark = ",")),
           hjust = 0, vjust = 0.5, size = 4.5, fontface = "bold", color = "#9ecae1") +
  annotate("text", x = 0.97, y = 0.28, label = sprintf("scATAC\n%s", format(orig_atac, big.mark = ",")),
           hjust = 1, vjust = 0.5, size = 4.5, fontface = "bold", color = "#fc9272")

ggsave(file.path(OUT_DIR, "F28_peak_overlap_venn.pdf"), p_venn, width = 5.5, height = 4.5, dpi = 300)
cat("  Saved F28_peak_overlap_venn.pdf\n")

# Save overlap statistics (reduced regions for Venn diagram)
overlap_stats <- data.frame(
  category = c("scG4_total_reduced", "scATAC_total_reduced", "overlap", "scG4_only", "scATAC_only"),
  count = c(n_g4_total, n_atac_total, n_overlap, n_g4_only, n_atac_only),
  pct_of_g4 = c(100, NA, 100 * n_overlap / n_g4_total, 100 * n_g4_only / n_g4_total, NA),
  pct_of_atac = c(NA, 100, 100 * n_overlap / n_atac_total, NA, 100 * n_atac_only / n_atac_total)
)
fwrite(overlap_stats, file.path(OUT_DIR, "F28_peak_overlap_stats.csv"))

# Also save original peak counts for reference
original_stats <- data.frame(
  category = c("scG4_original", "scATAC_original"),
  count = c(length(g4_peaks_gr), length(atac_peaks_gr)),
  note = c("Original peak count (before reduce)", "Original peak count (canonical chromosomes)")
)
fwrite(original_stats, file.path(OUT_DIR, "F28_peak_overlap_stats_original.csv"), append = FALSE)
cat("  Saved F28_peak_overlap_stats.csv\n")

cat("\n=== Done ===\n")
cat("Output directory:", OUT_DIR, "\n")
