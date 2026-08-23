# =============================================================
# F27: G4 peak overlap and signal at differentially expressed
# genes between mESC and MEF.
#
# Pipeline:
#   1. Build consensus G4 peaks from cluster0 + cluster1 narrowPeak
#      files (resize ±500 bp), score with RPGC bigWigs.
#   2. From RNA-seq (GSE90894), select DE genes by |lfc| > 4.
#   3. Intersect DE gene TSS (RefSeq ±1kb) with G4 peaks.
#   4. Stacked bar: fraction of DE genes with/without G4 overlap.
#   5. Violin of G4 lfc (log2 cl1/cl0) per DE category.
#   6. RPGC dodged violins (MEF/ESC, y 0–100) with mean+SD.
#
# Output (rev/outputs/G4_vs_DE_genes/):
#     panel_F27.pdf / .png              (stacked bar + lfc violin)
#     panel_F27_rpgc_violin.pdf / .png  (MEF/ESC RPGC dodged violins)
#     F27_DE_gene_overlap.csv
#     F27_paired_ttest_g4_signal.csv
# =============================================================


source("rev/paths.R")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(bedscout)
  library(wigglescout)
  library(data.table)
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
})

OUT_DIR <- file.path(OUT_ROOT, "G4_vs_DE_genes")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

LFC_CUTOFF <- 4

# ----- 1. Build consensus peaks (cluster0 + cluster1) -----
cat("Building consensus G4 peaks (cluster0 + cluster1)...\n")
MIX_PEAKS <- c(
  file.path(GSE291468, "GSM8836088_cluster_0_peaks.narrowPeak"),
  file.path(GSE291468, "GSM8836088_cluster_1_peaks.narrowPeak")
)
BW_CLUSTERS <- c(
  file.path(GSE291468, "GSM8836088_cluster0_RPGC.bw"),
  file.path(GSE291468, "GSM8836088_cluster1_RPGC.bw")
)
stopifnot(all(file.exists(MIX_PEAKS)), all(file.exists(BW_CLUSTERS)))

pkl_mix <- lapply(MIX_PEAKS, function(f) {
  peaks_universe(list(f), min_consensus = 1, resize_pre = 500, resize_post = 500)
})
names(pkl_mix) <- c("cluster_0", "cluster_1")
cat(sprintf("  scG4 peaks: cluster_0 %s, cluster_1 %s\n",
            format(length(pkl_mix[[1]]), big.mark = ","),
            format(length(pkl_mix[[2]]), big.mark = ",")))

pkl_all <- GRangesList(cl0 = pkl_mix[[1]], cl1 = pkl_mix[[2]])
for (n in names(pkl_all)) mcols(pkl_all[[n]])$source <- n
pkl_all <- unlist(pkl_all, use.names = TRUE)
pk_all <- GenomicRanges::reduce(pkl_all, with.revmap = TRUE)
mcols(pk_all)$source <- sapply(mcols(pk_all)$revmap, function(idx) {
  paste(sort(unique(mcols(pkl_all)$source[idx])), collapse = ",")
})
pk_all <- pk_all[grepl("chr", seqnames(pk_all)), ]
pk_all$revmap <- NULL
cat(sprintf("  Union of peaks (canonical chr): %s\n",
            format(length(pk_all), big.mark = ",")))

# ----- 2. Score peaks with RPGC bigWigs -----
cat("Scoring peaks with bw_loci...\n")
cov.pk_all <- bw_loci(BW_CLUSTERS, loci = pk_all, labels = c("cl0", "cl1"))
mcols(pk_all) <- cbind(mcols(pk_all), mcols(cov.pk_all))
pk_all$lfc  <- log2((pk_all$cl1 + 0.01) / (pk_all$cl0 + 0.01))
pk_all$mean <- (pk_all$cl1 + pk_all$cl0) / 2
cat(sprintf("  Scored %s peaks\n", format(length(pk_all), big.mark = ",")))

# ----- 3. RNA-seq: select DE genes by lfc cutoff -----
cat("Loading RNA-seq table and selecting DE genes...\n")
rnaseq <- read_xlsx(RNASEQ_XLSX, sheet = 1, skip = 4)
rnaseq$gn <- sub(";.+", "", rnaseq$name)
rnaseq <- rnaseq[!is.na(rnaseq$MEFs) & !is.na(rnaseq$ESCs), ]
rnaseq$lfc_rna <- log2((rnaseq$ESCs + 0.001) / (rnaseq$MEFs + 0.001))

top_up   <- rnaseq$gn[rnaseq$lfc_rna > LFC_CUTOFF]       # ESC-up
top_down <- rnaseq$gn[rnaseq$lfc_rna < -LFC_CUTOFF]       # MEF-up
n_up <- length(top_up)
n_down <- length(top_down)
cat(sprintf("  |lfc| > %d: %d ESC-up, %d MEF-up (%d total)\n",
            LFC_CUTOFF, n_up, n_down, n_up + n_down))

# ----- 4. Load RefSeq TSS and find overlaps -----
cat("Loading RefSeq TSS windows and finding peak overlaps...\n")
tss_gr <- rtracklayer::import(REFGENE_TSS)
tss_gr$gene_sym <- sub(";.*", "", tss_gr$name)
cat(sprintf("  RefSeq TSS windows: %s (%s gene symbols)\n",
            format(length(tss_gr), big.mark = ","),
            format(uniqueN(tss_gr$gene_sym), big.mark = ",")))

ov <- findOverlaps(pk_all, tss_gr, ignore.strand = TRUE)
peak_gene <- data.table(
  peak_idx  = queryHits(ov),
  lfc       = pk_all$lfc[queryHits(ov)],
  cl0       = pk_all$cl0[queryHits(ov)],
  cl1       = pk_all$cl1[queryHits(ov)],
  gene      = tss_gr$gene_sym[subjectHits(ov)]
)
cat(sprintf("  Overlaps: %s (mapping to %s genes)\n",
            format(nrow(ov), big.mark = ","),
            format(uniqueN(peak_gene$gene), big.mark = ",")))

# ----- 5. Overlap counts -----
cat("Computing overlap with DE gene sets...\n")
has_g4_up   <- sum(top_up %in% peak_gene$gene)
has_g4_down <- sum(top_down %in% peak_gene$gene)

overlap_dt <- data.table(
  de_set = rep(factor(c("ESC-up", "MEF-up"), levels = c("ESC-up", "MEF-up")), each = 2),
  status = rep(c("G4 overlapping", "Non-overlapping"), 2),
  n      = c(has_g4_up, n_up - has_g4_up, has_g4_down, n_down - has_g4_down),
  pct    = c(100 * has_g4_up / n_up, 100 * (n_up - has_g4_up) / n_up,
             100 * has_g4_down / n_down, 100 * (n_down - has_g4_down) / n_down)
)
fwrite(overlap_dt, file.path(OUT_DIR, "F27_DE_gene_overlap.csv"))

cat(sprintf("  ESC-up: %d / %d (%.1f%%) with G4 peak\n",
            has_g4_up, n_up, 100 * has_g4_up / n_up))
cat(sprintf("  MEF-up: %d / %d (%.1f%%) with G4 peak\n",
            has_g4_down, n_down, 100 * has_g4_down / n_down))

# ----- 6. G4 signal for overlapping genes (one peak/gene, max |lfc|) -----
peak_gene[, best_lfc := lfc[which.max(abs(lfc))], by = gene]
peak_gene[, best_cl0 := cl0[which.max(abs(lfc))], by = gene]
peak_gene[, best_cl1 := cl1[which.max(abs(lfc))], by = gene]
peak_gene_best <- unique(peak_gene[, .(gene, lfc = best_lfc, cl0 = best_cl0, cl1 = best_cl1)])

lfc_up <- peak_gene_best[gene %in% top_up, ]
lfc_up$set <- "ESC-up"
lfc_down <- peak_gene_best[gene %in% top_down, ]
lfc_down$set <- "MEF-up"
sig_all <- rbind(lfc_up, lfc_down)
sig_all$set <- factor(sig_all$set, levels = c("ESC-up", "MEF-up"))

cat(sprintf("  Genes with G4 signal (one peak/gene): ESC-up %d, MEF-up %d\n",
            nrow(lfc_up), nrow(lfc_down)))

# Paired t-tests (MEF cl0 vs ESC cl1 within each DE category)
ttests <- do.call(rbind, lapply(levels(sig_all$set), function(s) {
  sub <- sig_all[set == s]
  tt  <- t.test(sub$cl0, sub$cl1, paired = TRUE)
  p_lab <- ifelse(tt$p.value < 0.001, sprintf("p=%.1e", tt$p.value),
                  sprintf("p=%.3f", tt$p.value))
  data.frame(set = s, p_value = tt$p.value, p_label = p_lab,
             mean_MEF = mean(sub$cl0), mean_ESC = mean(sub$cl1),
             stringsAsFactors = FALSE)
}))
fwrite(ttests, file.path(OUT_DIR, "F27_paired_ttest_g4_signal.csv"))
cat("  Paired t-tests (MEF cl0 vs ESC cl1 G4 signal):\n")
print(ttests[, c("set", "p_label", "mean_MEF", "mean_ESC")])

# Reshape for dodged violins
sig_long <- tidyr::pivot_longer(
  sig_all[, .(gene, set, cl0, cl1)],
  cols = c(cl0, cl1),
  names_to = "signal", values_to = "g4_signal"
)
sig_long$signal <- factor(sig_long$signal,
                          levels = c("cl0", "cl1"),
                          labels = c("MEF", "ESC"))

# ----- 7. Theme -----
theme_fig1 <- theme_classic() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black")
  )

# ----- 8. Stacked bar -----
overlap_dt$status <- factor(overlap_dt$status,
                            levels = c("Non-overlapping", "G4 overlapping"))

p_bar <- ggplot(overlap_dt, aes(x = de_set, y = pct, fill = status)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
            position = position_stack(vjust = 0.5),
            size = 4, fontface = "bold") +
  scale_fill_manual(values = c("G4 overlapping" = "#d62728",
                                "Non-overlapping" = "#bdbdbd")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Fraction of DE genes (%)",
       title = "G4 peak overlap at DE genes",
       fill = NULL) +
  theme_fig1 +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))

# ----- 9. G4 lfc violin -----
p_lfc <- ggviolin(
  sig_all, "set", "lfc",
  fill = "set",
  add = "mean_sd",
  add.params = list(size = 0.6)
) +
  scale_fill_manual(values = c("ESC-up" = "#d62728", "MEF-up" = "#1f77b4")) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
  stat_compare_means(comparisons = list(c("ESC-up", "MEF-up")),
                     method = "wilcox.test", label = "p.signif",
                     size = 5, bracket.size = 0.3) +
  labs(x = NULL, y = "G4 signal log2(cluster1 / cluster0)",
       title = "G4 signal lfc at DE genes") +
  theme_fig1 +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

# ----- 10. RPGC dodged violins (y 0–100) -----
y_pos <- 100 * 1.02

p_rpgc <- ggviolin(
  sig_long, "set", "g4_signal",
  fill = "signal",
  add = "mean_sd",
  add.params = list(size = 0.6)
) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_fill_manual(values = c("MEF" = "#1f77b4", "ESC" = "#d62728")) +
  annotate("text", x = seq_along(levels(sig_all$set)), y = y_pos,
           label = ttests$p_label, size = 4.5, fontface = "italic") +
  labs(x = NULL, y = "G4 signal (RPGC)",
       title = "G4 RPGC at DE genes",
       fill = "Cell type") +
  theme_fig1 +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top")

# ----- 11. Save -----
p_main <- p_bar + p_lfc + plot_layout(widths = c(1, 1))

ggsave(file.path(OUT_DIR, "panel_F27.pdf"),
       p_main, width = 12, height = 6, device = "pdf")
ggsave(file.path(OUT_DIR, "panel_F27.png"),
       p_main, width = 12, height = 6, dpi = 300)

ggsave(file.path(OUT_DIR, "panel_F27_rpgc_violin.pdf"),
       p_rpgc, width = 7, height = 6, device = "pdf")
ggsave(file.path(OUT_DIR, "panel_F27_rpgc_violin.png"),
       p_rpgc, width = 7, height = 6, dpi = 300)

cat("\nWrote outputs to", OUT_DIR, "\n")