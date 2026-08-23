#!/usr/bin/env Rscript
# ===========================================================================
# Figure F29: scG4 vs scATAC cluster-cluster correlation (Pearson)
#
# Tests whether scG4 signal reflects chromatin accessibility by computing
# correlations between scG4 clusters and scATAC clusters at 1:1 matched peaks.
#
# Output (github/rev/outputs/scG4_vs_scATAC/):
#   - F29_cluster_correlation_pearson.pdf
#   - F29_cluster_correlation_pearson.csv
# ===========================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(ggplot2)
  library(reshape2)
  library(data.table)
})

# ----- Paths -----
ROOT <- normalizePath(getwd())
DATA <- file.path(ROOT, "github/data")
OUT_DIR <- file.path(ROOT, "github/rev/outputs/scG4_vs_scATAC")
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

# Load scATAC object
cat("Loading scATAC object...\n")
atac <- readRDS(ATAC_RDS)
cat(sprintf("  Cells: %s, Peaks: %s, Clusters: %s\n",
            ncol(atac), nrow(atac), paste(levels(Idents(atac)), collapse=", ")))

# Extract count matrices
cat("\nExtracting count matrices...\n")
g4_counts <- GetAssayData(g4, assay = "peaks", layer = "counts")
atac_counts <- GetAssayData(atac, assay = "peaks", layer = "counts")
cat(sprintf("  scG4: %s peaks x %s cells\n", nrow(g4_counts), ncol(g4_counts)))
cat(sprintf("  scATAC: %s peaks x %s cells\n", nrow(atac_counts), ncol(atac_counts)))

# Extract peak coordinates
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

# Find 1:1 matched peaks
cat("\nFinding 1:1 matched peaks...\n")
ov <- findOverlaps(g4_peaks_gr, atac_peaks_gr, ignore.strand = TRUE)
ov_df <- as.data.frame(ov)

g4_overlap_count <- table(ov_df$queryHits)
atac_overlap_count <- table(ov_df$subjectHits)

g4_unique <- as.integer(names(g4_overlap_count[g4_overlap_count == 1]))
atac_unique <- as.integer(names(atac_overlap_count[atac_overlap_count == 1]))

ov_filtered <- ov_df[ov_df$queryHits %in% g4_unique & ov_df$subjectHits %in% atac_unique, ]
g4_idx <- ov_filtered$queryHits
atac_idx <- ov_filtered$subjectHits

cat(sprintf("  1:1 matched peaks: %s\n", format(length(g4_idx), big.mark = ",")))

# Calculate pseudobulk for each cluster
cat("\nCalculating pseudobulk profiles per cluster...\n")

g4_clusters <- levels(Idents(g4))
g4_pseudobulk <- lapply(g4_clusters, function(cl) {
  cells <- WhichCells(g4, idents = cl)
  counts <- g4_counts[g4_idx, cells, drop = FALSE]
  lib_size <- colSums(counts)
  cpm <- sweep(counts, 2, lib_size, "/") * 1e6
  rowMeans(cpm)
})
names(g4_pseudobulk) <- paste0("G4_cl", g4_clusters)

atac_clusters <- levels(Idents(atac))
atac_pseudobulk <- lapply(atac_clusters, function(cl) {
  cells <- WhichCells(atac, idents = cl)
  counts <- atac_counts[atac_idx, cells, drop = FALSE]
  lib_size <- colSums(counts)
  cpm <- sweep(counts, 2, lib_size, "/") * 1e6
  rowMeans(cpm)
})
names(atac_pseudobulk) <- paste0("ATAC_cl", atac_clusters)

cat(sprintf("  scG4 clusters: %s\n", paste(names(g4_pseudobulk), collapse = ", ")))
cat(sprintf("  scATAC clusters: %s\n", paste(names(atac_pseudobulk), collapse = ", ")))

# Compute Pearson correlation matrix
cat("\nComputing Pearson correlation matrix...\n")
pearson_cor <- cor(
  do.call(cbind, c(g4_pseudobulk, atac_pseudobulk)),
  method = "pearson"
)

# Plot correlation matrix
p_pearson <- ggplot(reshape2::melt(pearson_cor), aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-1, 1),
                       name = "r") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8)) +
  labs(x = NULL, y = NULL, title = "Pearson", subtitle = NULL)

ggsave(file.path(OUT_DIR, "F29_cluster_correlation_pearson.pdf"), p_pearson, width = 5.5, height = 4.5, dpi = 300)
cat("  Saved F29_cluster_correlation_pearson.pdf\n")

# Save correlation matrix
fwrite(as.data.frame(pearson_cor), file.path(OUT_DIR, "F29_cluster_correlation_pearson.csv"))
cat("  Saved F29_cluster_correlation_pearson.csv\n")

cat("\n=== Done ===\n")
cat("Output directory:", OUT_DIR, "\n")