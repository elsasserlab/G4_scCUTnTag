#!/usr/bin/env Rscript
# ===========================================================================
# Figure F28: scG4 vs scATAC peak overlap
#
# Tests whether scG4 signal reflects chromatin accessibility by comparing
# GFP+ scG4 (GSM8836086) with GFP+ scATAC (GSE198467) from same cell population.
#
# Output (github/rev/outputs/scG4_vs_scATAC/):
#   - F28_peak_overlap_venn.pdf
#   - F28_peak_overlap_stats.csv
# ===========================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(ggplot2)
  library(eulerr)
  library(data.table)
})

# ----- Paths -----
ROOT <- normalizePath(getwd())
DATA <- file.path(ROOT, "data")
OUT_DIR <- file.path(ROOT, "rev/outputs/scG4_vs_scATAC")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# scG4 GFP+ data
G4_RDS <- file.path(DATA, "GSE291468/GSM8836086_GFPpos_Seurat_object.Rds")

# scATAC data
ATAC_RDS <- file.path(DATA, "GSE198467/GSE198467_ATAC_Seurat_object_clustered_renamed.Rds")

cat("=== Loading data ===\n")

# Load scG4 GFP+ object
cat("Loading scG4 GFP+ object...\n")
g4 <- readRDS(G4_RDS)
cat(sprintf("  Cells: %s, Peaks: %s, Clusters: %s\n",
            ncol(g4), nrow(g4), paste(levels(Idents(g4)), collapse=", ")))

# Extract peak coordinates from Seurat object rownames
cat("Extracting scG4 peak coordinates from Seurat object...\n")
g4_peak_names <- rownames(g4)
g4_peak_parts <- strsplit(g4_peak_names, "[-:]")
g4_peaks_gr <- GRanges(
  seqnames = sapply(g4_peak_parts, `[`, 1),
  ranges = IRanges(
    start = as.integer(sapply(g4_peak_parts, `[`, 2)),
    end = as.integer(sapply(g4_peak_parts, `[`, 3))
  )
)
cat(sprintf("  Total scG4 peaks: %s\n", format(length(g4_peaks_gr), big.mark = ",")))

# Load scATAC object
cat("Loading scATAC object...\n")
atac <- readRDS(ATAC_RDS)
cat(sprintf("  Cells: %s, Peaks: %s, Clusters: %s\n",
            ncol(atac), nrow(atac), paste(levels(Idents(atac)), collapse=", ")))

# Extract scATAC peak coordinates
cat("Extracting scATAC peak coordinates...\n")
atac_peak_names <- rownames(atac)
atac_peak_parts <- strsplit(atac_peak_names, "[-:]")
atac_peaks_gr <- GRanges(
  seqnames = sapply(atac_peak_parts, `[`, 1),
  ranges = IRanges(
    start = as.integer(sapply(atac_peak_parts, `[`, 2)),
    end = as.integer(sapply(atac_peak_parts, `[`, 3))
  )
)
# Filter to canonical chromosomes
CANONICAL <- c(paste0("chr", 1:19), "chrX", "chrY")
canonical_mask <- as.character(seqnames(atac_peaks_gr)) %in% CANONICAL
atac_peaks_gr <- atac_peaks_gr[canonical_mask]
cat(sprintf("  Total scATAC peaks (canonical): %s\n", format(length(atac_peaks_gr), big.mark = ",")))

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

p_venn <- plot(fit,
               fills = list(fill = c("#9ecae1", "#fc9272")),
               edges = list(col = "black", lwd = 1.5),
               labels = NULL,  # Remove set labels (use totals outside instead)
               quantities = list(font = 2, cex = 2.0, digits = 0))  # Show counts in regions

# Add original peak counts as text outside the circles
orig_g4 <- length(g4_peaks_gr)
orig_atac <- length(atac_peaks_gr)
p_venn <- p_venn +
  annotate("text", x = -1.8, y = 0.3, label = sprintf("scG4\n%s", format(orig_g4, big.mark = ",")),
           hjust = 1, vjust = 0.5, size = 4.5, fontface = "bold", color = "#9ecae1") +
  annotate("text", x = 1.8, y = -0.3, label = sprintf("scATAC\n%s", format(orig_atac, big.mark = ",")),
           hjust = 0, vjust = 0.5, size = 4.5, fontface = "bold", color = "#fc9272")

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
