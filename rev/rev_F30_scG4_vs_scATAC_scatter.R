#!/usr/bin/env Rscript
# ===========================================================================
# Figure F30: scG4 vs scATAC scatter plot (Union + Intersection)
#
# Compares GFP+ scG4 (GSM8836086) with GFP+ scATAC (GSE198467)
# 
# Approach:
#   1. Find G4 peaks that overlap with ATAC peaks
#   2. Merge overlapping peaks into unified intervals  
#   3. Count reads from BOTH assays in these merged intervals
#   4. Compare UNION (all peaks) vs INTERSECTION (overlapping peaks)
#
# Output (github/rev/outputs/scG4_vs_scATAC/):
#   - F30_scatter_union_intersection.pdf
#   - F30_intersection_data.csv
#   - F30_union_data.csv
# ===========================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(ggplot2)
  library(data.table)
  library(ggrastr)
})

set.seed(42)

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

# Load objects
cat("Loading scG4 GFP+ object...\n")
g4 <- readRDS(G4_RDS)
cat(sprintf("  Cells: %s, Peaks: %s\n", ncol(g4), nrow(g4)))

cat("Loading scATAC object...\n")
atac <- readRDS(ATAC_RDS)
cat(sprintf("  Cells: %s, Peaks: %s\n", ncol(atac), nrow(atac)))

# Extract count matrices (sum across all cells for total signal)
cat("\nExtracting count matrices...\n")
g4_counts <- GetAssayData(g4, assay = "peaks", layer = "counts")
atac_counts <- GetAssayData(atac, assay = "peaks", layer = "counts")

# Sum across all cells to get total signal per peak
g4_total <- Matrix::rowSums(g4_counts)
atac_total <- Matrix::rowSums(atac_counts)
cat(sprintf("  scG4 total counts: %s\n", format(sum(g4_total), big.mark = ",")))
cat(sprintf("  scATAC total counts: %s\n", format(sum(atac_total), big.mark = ",")))

# Extract peak coordinates as GRanges
cat("Extracting peak coordinates...\n")

g4_peak_names <- rownames(g4_counts)
g4_peak_parts <- strsplit(g4_peak_names, "[-:]")
g4_peaks_gr <- GRanges(
  seqnames = sapply(g4_peak_parts, `[`, 1),
  ranges = IRanges(
    start = as.integer(sapply(g4_peak_parts, `[`, 2)),
    end = as.integer(sapply(g4_peak_parts, `[`, 3))
  )
)

atac_peak_names <- rownames(atac_counts)
atac_peak_parts <- strsplit(atac_peak_names, "[-:]")
atac_peaks_gr <- GRanges(
  seqnames = sapply(atac_peak_parts, `[`, 1),
  ranges = IRanges(
    start = as.integer(sapply(atac_peak_parts, `[`, 2)),
    end = as.integer(sapply(atac_peak_parts, `[`, 3))
  )
)

# ===========================================================================
# Find overlapping peaks and create merged intervals
# ===========================================================================
cat("\n=== Finding overlapping peaks ===\n")

# Find all overlaps
ov <- findOverlaps(g4_peaks_gr, atac_peaks_gr, ignore.strand = TRUE)
ov_df <- as.data.frame(ov)
cat(sprintf("  Total overlap pairs: %s\n", format(nrow(ov_df), big.mark = ",")))

# Create merged intervals: union of each overlapping pair
cat("Creating merged intervals (INTERSECTION)...\n")

merged_intervals <- GRanges(
  seqnames = seqnames(g4_peaks_gr[ov_df$queryHits]),
  ranges = IRanges(
    start = pmin(start(g4_peaks_gr[ov_df$queryHits]), start(atac_peaks_gr[ov_df$subjectHits])),
    end = pmax(end(g4_peaks_gr[ov_df$queryHits]), end(atac_peaks_gr[ov_df$subjectHits]))
  )
)

# Add overlap pair info as metadata
mcols(merged_intervals) <- DataFrame(
  g4_peak_idx = ov_df$queryHits,
  atac_peak_idx = ov_df$subjectHits
)

# Reduce to merge overlapping merged intervals, keeping track of which original pairs contributed
cat("  Reducing merged intervals...\n")
merged_reduced <- reduce(merged_intervals, ignore.strand = TRUE, with.revmap = TRUE)
cat(sprintf("  Merged intervals (reduced, INTERSECTION): %s\n", format(length(merged_reduced), big.mark = ",")))

# ===========================================================================
# Create UNION peak set (all peaks from both datasets)
# ===========================================================================
cat("\nCreating UNION peak set...\n")

all_peaks <- c(g4_peaks_gr, atac_peaks_gr)
cat(sprintf("  Total peaks (concatenated): %s\n", format(length(all_peaks), big.mark = ",")))

union_peaks <- reduce(all_peaks, ignore.strand = TRUE, with.revmap = TRUE)
cat(sprintf("  Union regions (reduced): %s\n", format(length(union_peaks), big.mark = ",")))

# ===========================================================================
# Count signal in merged intervals (INTERSECTION)
# ===========================================================================
cat("\n=== Counting signal in INTERSECTION intervals ===\n")

# For each reduced interval, find which original overlap pairs it contains
# Then sum the counts from those peaks

g4_signal_intersect <- numeric(length(merged_reduced))
atac_signal_intersect <- numeric(length(merged_reduced))

for (i in seq_along(merged_reduced)) {
  # Get the indices of original overlap pairs that contributed to this merged interval
  pair_indices <- merged_reduced$revmap[[i]]
  
  if (length(pair_indices) > 0) {
    # Get the G4 and ATAC peak indices from these pairs
    g4_peaks_here <- unique(ov_df$queryHits[pair_indices])
    atac_peaks_here <- unique(ov_df$subjectHits[pair_indices])
    
    # Sum counts from all contributing peaks
    g4_signal_intersect[i] <- sum(g4_total[g4_peaks_here])
    atac_signal_intersect[i] <- sum(atac_total[atac_peaks_here])
  }
}

cat(sprintf("  INTERSECTION - G4 signal: min=%s, max=%s, mean=%.1f\n",
            min(g4_signal_intersect), format(max(g4_signal_intersect), big.mark = ","), mean(g4_signal_intersect)))
cat(sprintf("  INTERSECTION - ATAC signal: min=%s, max=%s, mean=%.1f\n",
            min(atac_signal_intersect), format(max(atac_signal_intersect), big.mark = ","), mean(atac_signal_intersect)))

# ===========================================================================
# Count signal in UNION intervals
# ===========================================================================
cat("\n=== Counting signal in UNION intervals ===\n")

# Map original peaks to union regions
g4_to_union <- findOverlaps(g4_peaks_gr, union_peaks, ignore.strand = TRUE)
atac_to_union <- findOverlaps(atac_peaks_gr, union_peaks, ignore.strand = TRUE)

g4_signal_union <- numeric(length(union_peaks))
atac_signal_union <- numeric(length(union_peaks))

# Convert overlaps to data frames (these are 1-based)
g4_ov_df <- as.data.frame(g4_to_union)
atac_ov_df <- as.data.frame(atac_to_union)

# Sum G4 signal per union region using aggregate
g4_signal_by_region <- aggregate(g4_total[g4_ov_df$queryHits] ~ g4_ov_df$subjectHits, FUN = sum)
colnames(g4_signal_by_region) <- c("region_idx", "signal")
g4_signal_union[g4_signal_by_region$region_idx] <- g4_signal_by_region$signal

# Sum ATAC signal per union region
atac_signal_by_region <- aggregate(atac_total[atac_ov_df$queryHits] ~ atac_ov_df$subjectHits, FUN = sum)
colnames(atac_signal_by_region) <- c("region_idx", "signal")
atac_signal_union[atac_signal_by_region$region_idx] <- atac_signal_by_region$signal

cat(sprintf("  UNION - G4 signal: min=%s, max=%s, mean=%.1f\n",
            min(g4_signal_union), format(max(g4_signal_union), big.mark = ","), mean(g4_signal_union)))
cat(sprintf("  UNION - ATAC signal: min=%s, max=%s, mean=%.1f\n",
            min(atac_signal_union), format(max(atac_signal_union), big.mark = ","), mean(atac_signal_union)))

# ===========================================================================
# Create scatter plot with UNION (gray) and INTERSECTION (dark blue)
# ===========================================================================
cat("\n=== Creating scatter plot ===\n")

# Create data frames for both sets
intersect_df <- data.frame(
  G4 = g4_signal_intersect,
  ATAC = atac_signal_intersect,
  category = "Intersection (overlapping peaks)"
)

union_df <- data.frame(
  G4 = g4_signal_union,
  ATAC = atac_signal_union,
  category = "Union (all peaks)"
)

# Filter to regions with signal in at least one assay
intersect_df <- intersect_df[intersect_df$G4 > 0 | intersect_df$ATAC > 0, ]
union_df <- union_df[union_df$G4 > 0 | union_df$ATAC > 0, ]

cat(sprintf("  INTERSECTION regions with signal: %s\n", format(nrow(intersect_df), big.mark = ",")))
cat(sprintf("  UNION regions with signal: %s\n", format(nrow(union_df), big.mark = ",")))

# Calculate correlations
# INTERSECTION: all have signal in both by definition
intersect_both <- intersect_df[intersect_df$G4 > 0 & intersect_df$ATAC > 0, ]
# UNION: calculate on ALL regions (including G4-only and ATAC-only)
union_all <- union_df[union_df$G4 > 0 | union_df$ATAC > 0, ]

cat(sprintf("  INTERSECTION with signal in both: %s\n", format(nrow(intersect_both), big.mark = ",")))
cat(sprintf("  UNION total regions: %s\n", format(nrow(union_all), big.mark = ",")))

pearson_intersect <- NA
spearman_intersect <- NA
pearson_union <- NA
spearman_union <- NA

if (nrow(intersect_both) > 10) {
  pearson_intersect <- cor(log10(intersect_both$G4 + 1), log10(intersect_both$ATAC + 1), method = "pearson")
  spearman_intersect <- cor(log10(intersect_both$G4 + 1), log10(intersect_both$ATAC + 1), method = "spearman")
  cat(sprintf("  INTERSECTION correlation: Pearson r=%.3f, Spearman ρ=%.3f\n", pearson_intersect, spearman_intersect))
}

if (nrow(union_all) > 10) {
  pearson_union <- cor(log10(union_all$G4 + 1), log10(union_all$ATAC + 1), method = "pearson")
  spearman_union <- cor(log10(union_all$G4 + 1), log10(union_all$ATAC + 1), method = "spearman")
  cat(sprintf("  UNION correlation (all regions): Pearson r=%.3f, Spearman ρ=%.3f\n", pearson_union, spearman_union))
}

# Generate plot with both overlaid
cat("Generating plot...\n")

# Use log scale
intersect_df$G4_log <- log10(intersect_df$G4 + 1)
intersect_df$ATAC_log <- log10(intersect_df$ATAC + 1)
union_df$G4_log <- log10(union_df$G4 + 1)
union_df$ATAC_log <- log10(union_df$ATAC + 1)

p_scatter <- ggplot() +
  # UNION (gray, low alpha) - plotted first so it's in the background
  ggrastr::rasterise(
    geom_point(data = union_df, aes(x = G4_log, y = ATAC_log),
               alpha = 0.2, size = 0.08, color = "gray60"),
    dpi = 300
  ) +
  # INTERSECTION (dark blue, high alpha) - plotted on top
  ggrastr::rasterise(
    geom_point(data = intersect_df, aes(x = G4_log, y = ATAC_log),
               alpha = 0.8, size = 0.12, color = "#1f5999"),
    dpi = 300
  ) +
  # Regression line for intersection
  geom_smooth(data = intersect_df, aes(x = G4_log, y = ATAC_log),
              method = "lm", se = FALSE, color = "#b2182b", linewidth = 0.35) +
  # Regression line for union
  geom_smooth(data = union_df, aes(x = G4_log, y = ATAC_log),
              method = "lm", se = FALSE, color = "gray40", linewidth = 0.2, linetype = "dashed") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text = element_text(size = 8.5),
        axis.title = element_text(size = 9.5)) +
  labs(x = "scG4 (log10(counts+1))",
       y = "scATAC (log10(counts+1))",
       title = "",
       subtitle = "") +
  annotate("text", x = 1, y = max(union_df$ATAC_log, na.rm = TRUE) * 0.95,
           label = sprintf("Intersection\nn = %s\nr = %.2f",
                          format(nrow(intersect_both), big.mark = "."),
                          pearson_intersect),
           hjust = 0, size = 3.2, fontface = "bold", color = "#1f5999") +
  annotate("text", x = 1, y = max(union_df$ATAC_log, na.rm = TRUE) * 0.80,
           label = sprintf("Union\nn = %s\nr = %.2f",
                          format(nrow(union_all), big.mark = "."),
                          pearson_union),
           hjust = 0, size = 2.8, fontface = "plain", color = "gray40")

ggsave(file.path(OUT_DIR, "F30_scatter_union_intersection.pdf"), p_scatter, width = 4.5, height = 3.8, dpi = 300)
cat("  Saved F30_scatter_union_intersection.pdf\n")

# Save data (both datasets)
fwrite(intersect_df, file.path(OUT_DIR, "F30_intersection_data.csv"))
fwrite(union_df, file.path(OUT_DIR, "F30_union_data.csv"))
cat("  Saved F30_intersection_data.csv and F30_union_data.csv\n")

cat("\n=== Done ===\n")
cat("Output directory:", OUT_DIR, "\n")
