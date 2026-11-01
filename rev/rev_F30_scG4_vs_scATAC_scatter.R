#!/usr/bin/env Rscript
# ===========================================================================
# Figure F30: scG4 vs scATAC scatter plot (Union + Intersection)
#
# Compares GFP+ scG4 (GSM8836086) with GFP+ scATAC (GSE198467)
# 
# Approach:
#   1. Find G4 peaks that overlap with ATAC peaks
#   2. Merge overlapping peaks into unified intervals  
#   3. Score BOTH assays with continuous coverage (bigWig + bw_loci) in these
#      intervals
#   4. Compare UNION (all peaks) vs INTERSECTION (overlapping peaks)
#
# Signal: continuous fragment coverage (RPGC-normalized, 10 bp bins) scored
# over the WHOLE interval with bw_loci, so single-modality regions can show
# genuine non-zero signal in the other assay (per-peak read counts forced 0
# there by construction). Both tracks are RPGC mean coverage:
#   - scG4:  GSM8836086_GFP_sorted_mousebrain.rpgc.bw (files in GSE291468/)
#   - scATAC: GSE198467_ATAC_fragments.RPGC.bw (data/GSE198467/,
#             built by data/GSE198467/make_atac_bigwig.sh from the raw fragments)
#
# scATAC peak coordinates come from the MACS2 re-called consensus narrowPeak
# + raw fragments file (data/GSE198467/), not the broad author-provided
# windows (median ~3.7 kb) stored in the Seurat rownames.
# scG4 peak coordinates are the merged UNION of the per-cluster MACS2
# narrowPeaks (GSM8836086_GFPpos_cluster_{0..3}_peaks.narrowPeak), reduced
# and restricted to canonical chromosomes (median width ~0.4 kb), NOT the
# pre-merged Seurat rownames (~0.9 kb) - consistent with F12/F28.
# See data/GSE198467/reprocess_atac_peaks.sh.
#
# Output (github/rev/outputs/scG4_vs_scATAC/):
#   - F30_scatter_union_intersection.pdf
#   - F30_intersection_data.csv
#   - F30_union_data.csv
# ===========================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(data.table)
  library(ggrastr)
  library(wigglescout)
})

set.seed(42)

# ----- Paths -----
ROOT <- normalizePath(getwd())
DATA <- file.path(ROOT, "data")
OUT_DIR <- file.path(ROOT, "rev/outputs/scG4_vs_scATAC")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# scG4 GFP+ data: merged union of the per-cluster MACS2 narrowPeaks
# (consistent with F12; NOT the pre-merged Seurat rownames)

G4_PEAK_CLUST <- file.path(DATA, "GSE291468",
                           paste0("GSM8836086_GFPpos_cluster_", 0:3, "_peaks.narrowPeak"))
G4_BW <- file.path(DATA, "GSE291468/GSM8836086_GFP_sorted_mousebrain.rpgc.bw")

# scATAC data
ATAC_NARROWPEAK <- file.path(DATA, "GSE198467/GSE198467_atac_consensus_peaks.narrowPeak")
ATAC_BW <- file.path(DATA, "GSE198467/GSE198467_ATAC_fragments.RPGC.bw")
CANONICAL <- c(paste0("chr", 1:19), "chrX", "chrY")

for (f in c(G4_PEAK_CLUST, G4_BW, ATAC_NARROWPEAK, ATAC_BW)) {
  stopifnot(file.exists(f))
}

cat("=== Loading data ===\n")

cat("Merging GFP+ cluster narrowPeaks into a union...\n")
g4_clusters <- lapply(G4_PEAK_CLUST, function(f) {
  import(f, format = "narrowPeak")
})
for (i in seq_along(g4_clusters)) {
  cat(sprintf("  GFP+ cluster %d peaks: %s\n", i - 1,
              format(length(g4_clusters[[i]]), big.mark = ",")))
}
g4_peaks_gr <- reduce(do.call(c, g4_clusters))
g4_peaks_gr <- g4_peaks_gr[as.character(seqnames(g4_peaks_gr)) %in% CANONICAL]
cat(sprintf("  scG4 peaks (canonical union): %s (median width %d bp)\n",
            format(length(g4_peaks_gr), big.mark = ","),
            median(width(g4_peaks_gr))))

# scATAC: MACS2 re-called consensus narrowPeak (canonical only)
cat("\nLoading re-called scATAC consensus narrowPeak...\n")
atac_peaks_gr <- import(ATAC_NARROWPEAK, format = "narrowPeak")
atac_peaks_gr <- atac_peaks_gr[as.character(seqnames(atac_peaks_gr)) %in% CANONICAL]
cat(sprintf("  scATAC peaks (canonical): %s (median width %d bp)\n",
            format(length(atac_peaks_gr), big.mark = ","),
            median(width(atac_peaks_gr))))

# ===========================================================================
# Find overlapping peaks and create merged intervals
# ===========================================================================
cat("\n=== Finding overlapping peaks ===\n")

ov <- findOverlaps(g4_peaks_gr, atac_peaks_gr, ignore.strand = TRUE)
ov_df <- as.data.frame(ov)
cat(sprintf("  Total overlap pairs: %s\n", format(nrow(ov_df), big.mark = ",")))

cat("Creating merged intervals (INTERSECTION)...\n")

# this is a reimplementation of reduce?
merged_intervals <- GRanges(
  seqnames = seqnames(g4_peaks_gr[ov_df$queryHits]),
  ranges = IRanges(
    start = pmin(start(g4_peaks_gr[ov_df$queryHits]), start(atac_peaks_gr[ov_df$subjectHits])),
    end = pmax(end(g4_peaks_gr[ov_df$queryHits]), end(atac_peaks_gr[ov_df$subjectHits]))
  )
)
mcols(merged_intervals) <- DataFrame(
  g4_peak_idx = ov_df$queryHits,
  atac_peak_idx = ov_df$subjectHits
)

cat("  Reducing merged intervals...\n")
# ======
# Hacking this here for the moment, to plot the scG4 peaks in blue, not the intersection
# merged_reduced <-  g4_peaks_gr
# ======
#merged_reduced <- reduce(merged_intervals, ignore.strand = TRUE, with.revmap = TRUE)

cat(sprintf("  Merged intervals (reduced, INTERSECTION): %s\n", format(length(merged_reduced), big.mark = ",")))

# ===========================================================================
# Create UNION peak set (all peaks from both datasets)
# ===========================================================================
cat("\nCreating UNION peak set...\n")

all_peaks <- c(g4_peaks_gr, atac_peaks_gr)
union_peaks <- reduce(all_peaks, ignore.strand = TRUE, with.revmap = TRUE)
cat(sprintf("  Union regions (reduced): %s\n", format(length(union_peaks), big.mark = ",")))

# ===========================================================================
# Score continuous coverage (bw_loci) in both interval sets
# ===========================================================================

score_loci <- function(regions, label) {
  cat("Scoring ", label, " with bw_loci...\n", sep = "")
  sig <- mcols(bw_loci(
    c(scG4 = G4_BW, scATAC = ATAC_BW),
    loci = regions,
    labels = c("scG4", "scATAC"),
    per_locus_stat = "mean",
    default_na = 0
  ))
  as.data.frame(sig)
}

cat("\n=== Scoring INTERSECTION intervals ===\n")
intersect_sig <- score_loci(merged_reduced, "INTERSECTION intervals")
cat(sprintf("  Mean scG4 coverage: %.3f ; mean scATAC coverage: %.3f\n",
            mean(intersect_sig$scG4), mean(intersect_sig$scATAC)))
cat(sprintf("  Zero-coverage counts: scG4=%s, scATAC=%s\n",
            sum(intersect_sig$scG4 == 0), sum(intersect_sig$scATAC == 0)))

cat("\n=== Scoring UNION intervals ===\n")
union_sig <- score_loci(union_peaks, "UNION intervals")
cat(sprintf("  Mean scG4 coverage: %.3f ; mean scATAC coverage: %.3f\n",
            mean(union_sig$scG4), mean(union_sig$scATAC)))
cat(sprintf("  Zero-coverage counts: scG4=%s, scATAC=%s\n",
            sum(union_sig$scG4 == 0), sum(union_sig$scATAC == 0)))

# ===========================================================================
# Create scatter plot with UNION (gray) and INTERSECTION (dark blue overlay)
# ===========================================================================
cat("\n=== Creating scatter plot ===\n")

intersect_df <- data.frame(
  G4 = intersect_sig$scG4,
  ATAC = intersect_sig$scATAC,
  category = "Intersection (overlapping peaks)"
)

union_df <- data.frame(
  G4 = union_sig$scG4,
  ATAC = union_sig$scATAC,
  category = "Union (all peaks)"
)

cat(sprintf("  N regions: INTERSECTION=%s, UNION=%s\n",
            format(nrow(intersect_df), big.mark = ","),
            format(nrow(union_df), big.mark = ",")))
cat(sprintf("  INTERSECTION with signal in both: %s\n",
            format(sum(intersect_df$G4 > 0 & intersect_df$ATAC > 0), big.mark = ",")))
cat(sprintf("  UNION with signal in both: %s\n",
            format(sum(union_df$G4 > 0 & union_df$ATAC > 0), big.mark = ",")))
cat(sprintf("  UNION with signal in only one: %s\n",
            format(sum((union_df$G4 > 0) != (union_df$ATAC > 0)), big.mark = ",")))

# log10(x + 1) transform
intersect_df$G4_log <- log10(intersect_df$G4 + 1)
intersect_df$ATAC_log <- log10(intersect_df$ATAC + 1)
union_df$G4_log <- log10(union_df$G4 + 1)
union_df$ATAC_log <- log10(union_df$ATAC + 1)

pearson_intersect <- cor(intersect_df$G4_log, intersect_df$ATAC_log, method = "pearson")
spearman_intersect <- cor(intersect_df$G4_log, intersect_df$ATAC_log, method = "spearman")
pearson_union <- cor(union_df$G4_log, union_df$ATAC_log, method = "pearson")
spearman_union <- cor(union_df$G4_log, union_df$ATAC_log, method = "spearman")

cat(sprintf("  INTERSECTION correlation: Pearson r=%.3f, Spearman ρ=%.3f\n",
            pearson_intersect, spearman_intersect))
cat(sprintf("  UNION correlation (all regions): Pearson r=%.3f, Spearman ρ=%.3f\n",
            pearson_union, spearman_union))

# Figure: build ggplot directly from the bw_loci data frames above (avoids a
# second bigwig read that plot_bw_loci_scatter would trigger internally).
cat("Generating plot...\n")

p_scatter <- ggplot() +
  ggrastr::rasterise(
    geom_point(data = union_df, aes(x = G4, y = ATAC),
               alpha = 0.2, size = 0.2, color = "gray60"),
    dpi = 300
  ) +
  ggrastr::rasterise(
    geom_point(data = intersect_df, aes(x = G4, y = ATAC),
               alpha = 0.5, size = 0.2, color = "#dc7252"),
    dpi = 300
  ) +
  scale_x_continuous(trans = "log1p") +
  scale_y_continuous(trans = "log1p") +
  theme_classic(base_size = 10) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text = element_text(size = 8.5),
        axis.title = element_text(size = 9.5)) +
  labs(x = "scG4 (coverage, log scale)",
       y = "scATAC (coverage, log scale)",
       title = "",
       subtitle = "") +
  annotate("text", x = 10^0.6 - 1, y = 10^(max(union_df$ATAC_log) * 0.95) - 1,
           label = sprintf("Intersection\nn = %s\nr = %.2f",
                          format(nrow(intersect_df), big.mark = "."),
                          pearson_intersect),
           hjust = 0, size = 3.2, fontface = "bold", color = "#dc7252") +
  annotate("text", x = 10^0.6 - 1, y = 10^(max(union_df$ATAC_log) * 0.80) - 1,
           label = sprintf("Union\nn = %s\nr = %.2f",
                          format(nrow(union_df), big.mark = "."),
                          pearson_union),
           hjust = 0, size = 2.8, fontface = "plain", color = "gray40")

ggsave(file.path(OUT_DIR, "F30_scatter_union_intersection.pdf"), p_scatter, width = 4.5, height = 3.8, dpi = 300)
cat("  Saved F30_scatter_union_intersection.pdf\n")

fwrite(intersect_df, file.path(OUT_DIR, "F30_intersection_data.csv"))
fwrite(union_df, file.path(OUT_DIR, "F30_union_data.csv"))
cat("  Saved F30_intersection_data.csv and F30_union_data.csv\n")

cat("\n=== Done ===\n")
cat("Output directory:", OUT_DIR, "\n")