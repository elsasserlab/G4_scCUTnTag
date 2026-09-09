#!/usr/bin/env Rscript
# ===========================================================================
# QC Comparison: GSM5501245 (SINAG1) vs ESC-MEF (GSM8836088) from Fig 1
#
# Compares per-cell QC metrics between:
# - GSM5501245: U2OS+MCF7 mixed (from paper: Sci Rep 11, 23641, 2021)
# - GSM8836088: mESC+MEF mixture (from our Fig 1 Panel C)
#
# Output: rev/outputs/QC_GSM5501245/panel_QC_GSM5501245.pdf
# ===========================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(cowplot)
  library(data.table)
  library(dplyr)
})

# ----- Paths -----
ROOT <- normalizePath(getwd())
DATA <- file.path(ROOT, "data")
OUT_DIR <- file.path(ROOT, "rev/outputs/QC_GSM5501245")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ESC-MEF data paths
SEURAT_ESCMEF <- file.path(DATA, "GSE291468/GSM8836088_mESCMEF_Seurat_object.Rds")
FRAGMENTS_ESCMEF <- file.path(DATA, "GSE291468/GSM8836088_mESC_MEF/CellRanger/fragments.tsv.gz")

cat("=== Loading ESC-MEF (GSM8836088) data ===\n")
if (!file.exists(SEURAT_ESCMEF)) {
  stop("ESC-MEF Seurat object not found: ", SEURAT_ESCMEF)
}
seurat_escmef <- readRDS(SEURAT_ESCMEF)
cat(sprintf("  Loaded %s cells, %s features\n", ncol(seurat_escmef), nrow(seurat_escmef)))

# Compute QC metrics if not present
if (!"FRiP" %in% colnames(seurat_escmef@meta.data)) {
  cat("  Computing FRiP...\n")
  seurat_escmef$FRiP <- seurat_escmef$peak_region_fragments / seurat_escmef$passed_filters
}

if (!"TSS_enrichment" %in% colnames(seurat_escmef@meta.data)) {
  cat("  Computing TSS enrichment...\n")
  frag_list <- Fragments(seurat_escmef[["peaks"]])
  if (length(frag_list) > 0) {
    for (i in seq_along(frag_list)) frag_list[[i]]@path <- FRAGMENTS_ESCMEF
    seurat_escmef[["peaks"]]@fragments <- frag_list
    seurat_escmef <- TSSEnrichment(seurat_escmef, fast = TRUE)
    if ("TSS.enrichment" %in% colnames(seurat_escmef@meta.data)) {
      seurat_escmef$TSS_enrichment <- seurat_escmef$TSS.enrichment
    }
  }
}

if (!"nucleosome_signal" %in% colnames(seurat_escmef@meta.data)) {
  cat("  Computing nucleosome signal...\n")
  seurat_escmef <- NucleosomeSignal(seurat_escmef)
}

# Extract ESC-MEF metrics by cluster
cat("\n=== Extracting ESC-MEF metrics ===\n")
escmef_meta <- seurat_escmef@meta.data
escmef_meta$cluster <- seurat_escmef$seurat_clusters
escmef_meta$dataset <- "ESC-MEF\n(GSM8836088)"

escmef_summary <- escmef_meta %>%
  group_by(cluster, dataset) %>%
  summarise(
    n_cells = n(),
    median_nFeature = median(nFeature_peaks),
    median_nCount = median(nCount_peaks),
    median_FRiP = median(FRiP),
    median_TSS = median(TSS_enrichment, na.rm = TRUE),
    median_nuc = median(nucleosome_signal, na.rm = TRUE),
    .groups = "drop"
  )

cat("  ESC-MEF Cluster 0 (MEF):", escmef_summary$n_cells[1], "cells\n")
cat("  ESC-MEF Cluster 1 (mESC):", escmef_summary$n_cells[2], "cells\n")

# GSM5501245 metrics from paper (Sci Rep 11, 23641, 2021)
# SINAG1 = mixed U2OS+MCF7
cat("\n=== Creating GSM5501245 metrics from paper ===\n")
gsm5501245_data <- data.frame(
  cluster = c("MCF7", "U2OS"),
  dataset = "U2OS/MCF7\n(GSM5501245)",
  n_cells = c(593, 2071),
  median_nFeature = c(NA, NA),  # Not reported in paper
  median_nCount = c(NA, NA),    # Not reported in paper
  median_FRiP = c(NA, NA),      # Not reported in paper
  median_TSS = c(NA, NA),       # Not reported in paper
  median_nuc = c(NA, NA)        # Not reported in paper
)

# Paper reports median fragments/cell: 739 (MCF7), 939 (U2OS)
# We'll use these as proxy for nCount
gsm5501245_data$median_nCount <- c(739, 939)

cat("  GSM5501245 MCF7: 593 cells, 739 median fragments/cell\n")
cat("  GSM5501245 U2OS: 2,071 cells, 939 median fragments/cell\n")

# Combine datasets
combined <- rbind(
  escmef_summary %>% select(cluster, dataset, n_cells, median_nFeature, median_nCount, median_FRiP, median_TSS, median_nuc),
  gsm5501245_data
)

# Save summary table
fwrite(combined, file.path(OUT_DIR, "QC_comparison_summary.csv"))
cat("\n  Saved QC_comparison_summary.csv\n")

# ===========================================================================
# Create comparison plots
# ===========================================================================
cat("\n=== Creating comparison plots ===\n")

# Plot 1: Number of cells
p_cells <- ggplot(combined, aes(x = dataset, y = n_cells, fill = cluster)) +
  geom_col(position = "dodge", color = "grey30", width = 0.7) +
  geom_text(aes(label = format(n_cells, big.mark = ",")), 
            position = position_dodge(width = 0.7), vjust = -0.5, size = 4, fontface = "bold") +
  scale_y_log10() +
  scale_fill_manual(values = c("#9ecae1", "#fc9272", "#a1d99b", "#fdb462")) +
  labs(x = NULL, y = "Number of cells (log10)",
       title = "Cell Recovery",
       subtitle = "Cells passing QC filters") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        axis.text.x = element_text(size = 11, angle = 0),
        axis.text.y = element_text(size = 11),
        legend.position = "bottom",
        legend.title = element_blank())

# Plot 2: Median fragments per cell (nCount)
p_nCount <- ggplot(combined, aes(x = dataset, y = median_nCount, fill = cluster)) +
  geom_col(position = "dodge", color = "grey30", width = 0.7) +
  geom_text(aes(label = ifelse(is.na(median_nCount), "NA", format(round(median_nCount), big.mark = ","))), 
            position = position_dodge(width = 0.7), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("#9ecae1", "#fc9272", "#a1d99b", "#fdb462"), na.value = "grey80") +
  labs(x = NULL, y = "Median fragments/cell",
       title = "Library Complexity",
       subtitle = "Median high-quality fragments per cell") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        axis.text.x = element_text(size = 11, angle = 0),
        axis.text.y = element_text(size = 11),
        legend.position = "none")

# Plot 3: Total peaks (aggregate)
# ESC-MEF peaks from Fig 1
cl0_peaks <- 31009  # MEF
cl1_peaks <- 34133  # mESC
total_escmef <- cl0_peaks + cl1_peaks

# GSM5501245 peaks
gsm5501245_peaks <- 17956  # From peak.bed file

peak_data <- data.frame(
  dataset = c("ESC-MEF", "U2OS/MCF7"),
  total_peaks = c(total_escmef, gsm5501245_peaks),
  cluster0 = c(cl0_peaks, NA),
  cluster1 = c(cl1_peaks, NA),
  combined = c(NA, gsm5501245_peaks)
)

p_peaks <- ggplot(peak_data, aes(x = dataset, y = total_peaks, fill = dataset)) +
  geom_col(color = "grey30", width = 0.6) +
  geom_text(aes(label = format(total_peaks, big.mark = ",")), 
            vjust = -0.5, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("#fc9272", "#b3de69")) +
  labs(x = NULL, y = "Total peaks",
       title = "Peak Numbers",
       subtitle = "Aggregate peaks across all cells") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        axis.text.x = element_text(size = 11, angle = 0),
        axis.text.y = element_text(size = 11),
        legend.position = "none")

# Plot 4: TSS enrichment (ESC-MEF only, GSM5501245 not reported)
p_tss <- ggplot(escmef_summary, aes(x = dataset, y = median_TSS, fill = cluster)) +
  geom_col(position = "dodge", color = "grey30", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", median_TSS)), 
            position = position_dodge(width = 0.7), vjust = -0.5, size = 4, fontface = "bold") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 1.5, y = 1.1, label = "QC threshold", color = "red", size = 3, fontface = "italic") +
  scale_fill_manual(values = c("#9ecae1", "#fc9272")) +
  labs(x = NULL, y = "TSS enrichment",
       title = "TSS Enrichment",
       subtitle = "Median per cluster (GSM5501245: not reported)") +
  scale_y_continuous(limits = c(0, max(escmef_summary$median_TSS) * 1.3)) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        axis.text.x = element_text(size = 11, angle = 0),
        axis.text.y = element_text(size = 11),
        legend.position = "none")

# ===========================================================================
# Combine all panels
# ===========================================================================
cat("\n=== Combining panels ===\n")
p_combined <- plot_grid(
  p_cells, p_nCount,
  p_peaks, p_tss,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D"),
  label_size = 14,
  align = "hv"
)

ggsave(file.path(OUT_DIR, "panel_QC_GSM5501245.pdf"), p_combined,
       width = 14, height = 10, dpi = 300)
cat("  Saved panel_QC_GSM5501245.pdf\n")

cat("\n=== Summary ===\n")
cat("ESC-MEF (GSM8836088):\n")
cat(sprintf("  Cluster 0 (MEF):  %s cells, median %.0f fragments/cell, TSS=%.1f\n",
            escmef_summary$n_cells[1], escmef_summary$median_nCount[1], escmef_summary$median_TSS[1]))
cat(sprintf("  Cluster 1 (mESC): %s cells, median %.0f fragments/cell, TSS=%.1f\n",
            escmef_summary$n_cells[2], escmef_summary$median_nCount[2], escmef_summary$median_TSS[2]))
cat(sprintf("  Total peaks: %s (MEF) + %s (mESC) = %s\n",
            format(cl0_peaks, big.mark = ","), format(cl1_peaks, big.mark = ","), format(total_escmef, big.mark = ",")))
cat("\nU2OS/MCF7 (GSM5501245, from paper):\n")
cat(sprintf("  MCF7:  593 cells, 739 median fragments/cell\n"))
cat(sprintf("  U2OS:  2,071 cells, 939 median fragments/cell\n"))
cat(sprintf("  Total peaks: %s\n", format(gsm5501245_peaks, big.mark = ",")))
cat("\nDone. Outputs in:", OUT_DIR, "\n")
