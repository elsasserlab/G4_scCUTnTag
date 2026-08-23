# =============================================================
# F12 extension: PQS overlap analysis for GSE198467 scATAC-seq
# mouse brain clusters, with width-adjusted peaks.
#
# scATAC peaks are resized to 270 bp (matching GFP+ sorted
# brain scG4 median width) centered on their midpoint, to enable
# fair PQS overlap comparison. Original widths are ~879 bp median,
# which would inflate overlap to ~100%.
#
# Extracts per-cluster peak sets from the GSE198467 Seurat object
# using stringent peak calling (peaks present in >=5% of cells)
# to approximate scG4 peak counts (~10k-30k per cluster).
#
# Output (rev/outputs/PQS_venn_diagrams/):
#   - bar_pct_peaks_overlap_PQS_scATAC.pdf (width-adjusted)
#   - scATAC_PQS_overlap_counts.csv (width-adjusted)
# =============================================================

source("rev/paths.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(data.table)
  library(ggplot2)
})

# ----- Paths -----
OUT_DIR <- file.path(OUT_ROOT, "PQS_venn_diagrams")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GSE198467_OBJ <- file.path(GSE198467, "GSE198467_ATAC_Seurat_object_clustered_renamed.Rds")
PQS_BED <- file.path(DATA, "pqsfinder/PQS_scores.mm10.bed")

stopifnot(file.exists(GSE198467_OBJ), file.exists(PQS_BED))

CANONICAL <- c(paste0("chr", 1:19), "chrX", "chrY")
MM10_GENOME_BP <- 2.65e9

# ----- Load scATAC object -----
cat("Loading GSE198467 scATAC-seq object...\n")
obj <- readRDS(GSE198467_OBJ)
cat(sprintf("  Object: %s features x %s cells\n",
            format(nrow(obj), big.mark = ","),
            format(ncol(obj), big.mark = ",")))
cat(sprintf("  Clusters: %s\n", paste(levels(Idents(obj)), collapse = ", ")))

# ----- Extract peak ranges -----
cat("\nExtracting peak ranges...\n")
peak_names <- rownames(obj)
# Parse peak names (format: chr1-1000-2000 or chr1:1000-2000)
peak_parts <- strsplit(peak_names, "[-:]")
all_peak_ranges <- GRanges(
  seqnames = sapply(peak_parts, `[`, 1),
  ranges = IRanges(
    start = as.integer(sapply(peak_parts, `[`, 2)),
    end = as.integer(sapply(peak_parts, `[`, 3))
  )
)
# Filter to canonical chromosomes
canonical_mask <- as.character(seqnames(all_peak_ranges)) %in% CANONICAL
peak_ranges <- all_peak_ranges[canonical_mask]
# Also subset the peak matrix
peak_matrix <- GetAssayData(obj, assay = "peaks", layer = "counts")[canonical_mask, , drop = FALSE]
cat(sprintf("  Total peaks (canonical chr): %s\n", format(length(peak_ranges), big.mark = ",")))
cat(sprintf("  Original median width: %d bp\n", median(width(peak_ranges))))

# ----- Resize scATAC peaks to match scG4 width -----
# GFP+ sorted brain scG4 peaks have median width 263-285 bp
# Resize scATAC peaks to 270 bp centered on midpoint for fair PQS comparison
TARGET_WIDTH <- 270L
cat(sprintf("\nResizing scATAC peaks to %d bp (matching GFP+ scG4 median width)...\n", TARGET_WIDTH))
peak_ranges_resized <- resize(peak_ranges, width = TARGET_WIDTH, fix = "center")
cat(sprintf("  Resized median width: %d bp\n", median(width(peak_ranges_resized))))

# ----- Get per-cluster peak sets (using resized peaks) -----
# Use more stringent peak calling: peaks must be present in >=5% of cells
# to approximate scG4 peak counts (~10k-30k per cluster)
cat("\nExtracting per-cluster peak sets (stringent: >=5% cells)...\n")

cluster_peaks <- lapply(levels(Idents(obj)), function(cluster_id) {
  cells_in_cluster <- WhichCells(obj, idents = cluster_id)
  n_cells <- length(cells_in_cluster)
  min_cells <- max(1, floor(0.05 * n_cells))  # >=5% of cells
  cluster_counts <- peak_matrix[, cells_in_cluster, drop = FALSE]
  # Peaks with at least 1 count in >=5% of cells
  peak_present <- Matrix::rowSums(cluster_counts > 0) >= min_cells
  peak_ranges_resized[peak_present]  # Use resized peaks
})
names(cluster_peaks) <- paste0("Cluster_", levels(Idents(obj)))

cat("  Peak counts per cluster:\n")
for(cl in names(cluster_peaks)) {
  cat(sprintf("    %-15s %s\n", cl, format(length(cluster_peaks[[cl]]), big.mark = ",")))
}

# ----- Load PQS sites -----
cat("\nLoading PQS sites (min score 20)...\n")
pqs <- import(PQS_BED, format = "BED")
pqs <- pqs[as.character(seqnames(pqs)) %in% CANONICAL]
cat(sprintf("  PQS sites: %s\n", format(length(pqs), big.mark = ",")))

# ----- Overlap counting function -----
overlap_counts <- function(peaks_gr, pqs_gr, label, genome_bp = MM10_GENOME_BP) {
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
  
  # Binomial test for enrichment (one-sided, alternative="greater")
  bin_p        <- binom.test(n_int, n_peaks, p = pqs_frac_adj, alternative = "greater")$p.value
  
  # Significance stars
  sig_label <- if(bin_p < 0.001) "***" else if(bin_p < 0.01) "**" else if(bin_p < 0.05) "*" else ""
  
  list(
    peaks_only        = n_peaks - n_int,
    pqs_only          = n_pqs   - n_int,
    overlap           = n_int,
    n_peaks           = n_peaks,
    n_pqs             = n_pqs,
    pct_peaks_overlap = 100 * n_int / n_peaks,
    pct_pqs_overlap   = 100 * n_pqs_in / n_pqs,
    expected_overlap  = exp_overlap,
    fold_enrichment   = fold_enrich,
    binom_pvalue      = bin_p,
    sig_label         = sig_label,
    label             = label
  )
}

# ----- Calculate overlap for ALL peaks combined (not per-cluster) -----
cat("\nCalculating PQS overlap for all scATAC peaks combined...\n")

# Merge all cluster peaks (union across all cells)
all_cluster_peaks <- peak_ranges_resized[Matrix::rowSums(peak_matrix) >= 1]
cat(sprintf("  Total scATAC peaks (union across all cells): %s\n",
            format(length(all_cluster_peaks), big.mark = ",")))

oc_all <- overlap_counts(all_cluster_peaks, pqs, "scATAC (all peaks)")
cat(sprintf("  scATAC (all)    %6.1f%%  (%s / %s)  %.2fx  %s\n",
            oc_all$pct_peaks_overlap,
            format(oc_all$overlap, big.mark = ","),
            format(oc_all$n_peaks, big.mark = ","),
            oc_all$fold_enrichment,
            oc_all$sig_label))

# ----- Save overlap table -----
overlap_tbl <- data.frame(
  comparison = "scATAC brain (all peaks, width-adjusted to 270 bp)",
  n_peaks = oc_all$n_peaks,
  n_overlap_PQS = oc_all$overlap,
  pct_overlap_PQS = round(oc_all$pct_peaks_overlap, 1),
  pct_PQS_overlap_peaks = round(oc_all$pct_pqs_overlap, 2),
  expected_overlap = round(oc_all$expected_overlap),
  fold_enrichment = round(oc_all$fold_enrichment, 2),
  binomial_pvalue = oc_all$binom_pvalue,
  sig_label = oc_all$sig_label,
  row.names = NULL
)

write.csv(overlap_tbl, file.path(OUT_DIR, "scATAC_PQS_overlap_counts.csv"), row.names = FALSE)
cat("\n  Saved scATAC_PQS_overlap_counts.csv\n")

# ----- Generate bar plot (single bar) -----
cat("\nGenerating bar plot...\n")

bar_df <- data.frame(
  assay = "scATAC brain",
  pct = oc_all$pct_peaks_overlap,
  overlap = oc_all$overlap,
  n_peaks = oc_all$n_peaks,
  fold = oc_all$fold_enrichment,
  sig = oc_all$sig_label,
  stringsAsFactors = FALSE
)
bar_df$label <- sprintf("%.1f%%\n(%s / %s)\n%.2fx%s",
                        bar_df$pct,
                        format(bar_df$overlap, big.mark = ","),
                        format(bar_df$n_peaks, big.mark = ","),
                        bar_df$fold,
                        bar_df$sig)

p_scATAC <- ggplot(bar_df, aes(x = assay, y = pct)) +
  geom_col(width = 0.3, fill = "#b3de69", color = "grey30") +
  geom_text(aes(label = label), vjust = -0.3, size = 4, fontface = "bold") +
  scale_y_continuous(limits = c(0, 105),
                     breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL,
       y = "% of scATAC peaks overlapping >=1 PQS site",
       title = "Mouse brain scATAC-seq overlap with PQS (width-adjusted to 270 bp)",
       subtitle = "Fold-enrichment vs random expectation (binomial test, one-sided); ***p<0.001, **p<0.01, *p<0.05") +
  theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 11))

ggsave(file.path(OUT_DIR, "bar_pct_peaks_overlap_PQS_scATAC.pdf"), p_scATAC,
       width = 5, height = 6)
cat("  Saved bar_pct_peaks_overlap_PQS_scATAC.pdf\n")

cat("\n=== Done ===\n")
cat("Output directory:", OUT_DIR, "\n")