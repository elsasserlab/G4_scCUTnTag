#!/usr/bin/env Rscript
# ===========================================================================
# QC Panel for GSM5501245 (SINAG1 scG4 CUT&Tag - U2OS/MCF7 mixed)
#
# Third-party data from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5501245
# 
# Analogous to panel_S1_qc_GFPpos.pdf but adapted for aggregate data
# Metrics from: Scientific Reports 11, 23641 (2021) - Supplementary Table S2
#
# Output: rev/outputs/QC_GSM5501245/panel_QC_GSM5501245.pdf
# ===========================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(ggplot2)
  library(cowplot)
  library(data.table)
})

# ----- Paths -----
ROOT <- normalizePath(getwd())
DATA <- file.path(ROOT, "data")
OUT_DIR <- file.path(ROOT, "rev/outputs/QC_GSM5501245")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PEAKS_FILE <- file.path(DATA, "GSM5501245/GSM5501245_SINAG1.cellranger.peak.bed")

cat("=== Loading GSM5501245 peaks ===\n")
peaks <- import(PEAKS_FILE, format = "BED")
cat(sprintf("  Total peaks: %s\n", format(length(peaks), big.mark = ",")))

# Filter to canonical chromosomes
CANONICAL <- c(paste0("chr", 1:22), "chrX", "chrY")
peaks <- peaks[as.character(seqnames(peaks)) %in% CANONICAL]
cat(sprintf("  Canonical chromosomes: %s\n", format(length(peaks), big.mark = ",")))

# ----- Published metrics from paper (Sci Rep 11, 23641, 2021) -----
cat("\n=== Using published QC metrics ===\n")
# From Supplementary Table S2:
# SINAG1 (mixed U2OS+MCF7): 593 MCF7 + 2071 U2OS cells
# Median fragments/cell: 739 (MCF7), 939 (U2OS)
# Range: 652-1226 high-quality fragments/cell

published_metrics <- data.frame(
  metric = c("Total cells", "MCF7 cells", "U2OS cells",
             "Median fragments/cell (MCF7)", "Median fragments/cell (U2OS)",
             "Fragment range"),
  value = c("2,664", "593", "2,071", "739", "939", "652-1,226")
)

# ----- QC Metric 1: Peak width distribution -----
cat("\n=== Computing peak width statistics ===\n")
peak_widths <- width(peaks)

p_width_vln <- ggplot(data.frame(width = peak_widths), aes(x = "", y = width)) +
  geom_violin(fill = "#2171b5", color = "grey30", scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, 
               color = "red", linewidth = 1) +
  scale_y_log10() +
  labs(x = "", y = "Peak width (log10 bp)",
       title = "Peak Width Distribution",
       subtitle = sprintf("Median: %d bp, IQR: %d-%d bp",
                         median(peak_widths),
                         quantile(peak_widths, 0.25),
                         quantile(peak_widths, 0.75))) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# ----- QC Metric 2: Peak distribution by chromosome -----
cat("\n=== Computing chromosome distribution ===\n")
chr_counts <- as.data.frame(table(seqnames(peaks)))
setnames(chr_counts, c("Chromosome", "Count"))
chr_counts$Chromosome <- as.character(chr_counts$Chromosome)
chr_counts <- chr_counts[order(match(chr_counts$Chromosome, CANONICAL)),]

p_chr <- ggplot(chr_counts, aes(x = Chromosome, y = Count)) +
  geom_col(fill = "#fc9272", color = "grey30") +
  labs(x = "Chromosome", y = "Peak count",
       title = "Peaks per Chromosome",
       subtitle = sprintf("Total: %s peaks", format(length(peaks), big.mark = ","))) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

# ----- QC Metric 3: Published per-cell metrics -----
cat("\n=== Creating per-cell metrics panel ===\n")
p_cells <- ggplot(published_metrics[1:3,], aes(x = metric, y = 1, fill = metric)) +
  geom_col(color = "grey30") +
  geom_text(aes(label = value), size = 5, fontface = "bold", vjust = 0.5) +
  scale_fill_manual(values = c("#a1d99b", "#9ecae1", "#fc9272")) +
  labs(x = NULL, y = NULL,
       title = "Cell Recovery (from paper)",
       subtitle = "GSM5501245 SINAG1 (mixed U2OS+MCF7)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

# ----- QC Metric 4: Published fragment metrics -----
frag_metrics <- published_metrics[4:6,]
p_fragments <- ggplot(frag_metrics, aes(x = metric, y = 1, fill = metric)) +
  geom_col(color = "grey30") +
  geom_text(aes(label = value), size = 4.5, fontface = "bold", vjust = 0.5) +
  scale_fill_manual(values = c("#b3de69", "#b3de69", "#fdb462")) +
  labs(x = NULL, y = NULL,
       title = "Fragments per Cell (from paper)",
       subtitle = "Median high-quality fragments") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

# ----- Combine all panels -----
cat("\n=== Combining panels ===\n")
p_combined <- plot_grid(
  p_width_vln, p_chr,
  p_cells, p_fragments,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D"),
  label_size = 14,
  rel_heights = c(1, 0.8)
)

ggsave(file.path(OUT_DIR, "panel_QC_GSM5501245.pdf"), p_combined,
       width = 14, height = 10, dpi = 300)
cat("  Saved panel_QC_GSM5501245.pdf\n")

# Save summary stats
summary_stats <- data.frame(
  sample = "GSM5501245",
  cell_type = "U2OS_MCF7_SINAG1 (mixed)",
  total_cells = 2664,
  mcf7_cells = 593,
  u2os_cells = 2071,
  median_fragments_mcf7 = 739,
  median_fragments_u2os = 939,
  fragment_range = "652-1226",
  total_peaks = length(peaks),
  median_peak_width = median(peak_widths),
  mean_peak_width = mean(peak_widths)
)
fwrite(summary_stats, file.path(OUT_DIR, "GSM5501245_summary_stats.csv"))
cat("  Saved GSM5501245_summary_stats.csv\n")

cat("\nDone. Outputs in:", OUT_DIR, "\n")
cat("\nNote: Per-cell QC metrics (nFeature, nCount, TSS enrichment, FRiP, nucleosome signal)\n")
cat("are from the published paper (Sci Rep 11, 23641, 2021 - Supplementary Table S2).\n")
cat("The fragments.tsv.gz file is not available from GEO for this sample.\n")
