# =============================================================
# F24: RNA-seq expression (log2 ESC/MEF) of G4 peaks by
# differential category (violin) + split violin (RPKMs).
#
# Uses RefSeq TSS annotation (refGene.tss.1kb.mm10.bed) for
# gene annotation (consistent with original author's script).
#
# Pipeline:
#   1. Consensus G4 peaks: scG4 mixture clusters 0/1 only, resized
#      to +/- 500 bp around the summit.
#   2. Score with bw_loci on cluster 0/1 RPGC bigWigs; classify
#      peaks into low / unchanged / MEF / ESC.
#   3. Annotate to genes (refGene TSS +/- 1 kb).
#   4. Export per-peak table.
#   5. Violin plot of log2(ESC/MEF) RPKM per category.
#   6. Split violin of MEF and ESC RPKM per category (log10 axis)
#      with paired t-tests.
#
# Output (rev/outputs/RNAseq_lfc_violin/):
#     panel_F24_rnaseq_violin.pdf / .png
#     panel_F24_split_violin_rpkms.pdf / .png
#     F24_peak_categories.csv
#     F24_peak_categories_rnaseq.csv
#     F24_paired_ttest_rnaseq.csv
#     rnaseq_violin_stats.csv
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
})

# ----- Paths -----
OUT_DIR <- file.path(OUT_ROOT, "RNAseq_lfc_violin")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MIX_PEAKS <- c(
  file.path(GSE291468, "GSM8836088_cluster_0_peaks.narrowPeak"),
  file.path(GSE291468, "GSM8836088_cluster_1_peaks.narrowPeak")
)
BW_CLUSTERS <- c(
  file.path(GSE291468, "GSM8836088_cluster0_RPGC.bw"),
  file.path(GSE291468, "GSM8836088_cluster1_RPGC.bw")
)
stopifnot(all(file.exists(MIX_PEAKS)),
          all(file.exists(BW_CLUSTERS)),
          file.exists(REFGENE_TSS), file.exists(RNASEQ_XLSX))

# ----- 1. Consensus peaks (mixture clusters only) -----
cat("Building consensus peaks...\n")
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

# ----- 2. Score + classify diff categories -----
cat("Scoring cluster signal (bw_loci) and classifying...\n")
cov.pk_all <- bw_loci(BW_CLUSTERS, loci = pk_all, labels = c("cl0", "cl1"))
mcols(pk_all) <- cbind(mcols(pk_all), mcols(cov.pk_all))
pk_all$lfc  <- log2((pk_all$cl1 + 0.01) / (pk_all$cl0 + 0.01))
pk_all$mean <- (pk_all$cl1 + pk_all$cl0) / 2
pk_all$diff <- factor(ifelse(
  pk_all$mean < 0.5, "low",
  ifelse(pk_all$lfc > 2, "ESC",
         ifelse(pk_all$lfc < (-2), "MEF", "unchanged"))
), levels = c("low", "unchanged", "MEF", "ESC"))

# ----- 3. Gene annotation (RefSeq TSS ± 1 kb) - nearest feature -----
cat("Annotating peaks to genes (RefSeq TSS nearest feature)...\n")
pk_all_anno <- annotate_nearby_features(
  pk_all,
  rtracklayer::import(REFGENE_TSS),
  name_field = "name",
  distance_cutoff = 1000
)
pk_all$gene <- pk_all_anno$annotation
cat(sprintf("  Peaks with gene annotation: %s (%.1f%%)\n",
            format(sum(!is.na(pk_all$gene)), big.mark = ","),
            100 * sum(!is.na(pk_all$gene)) / length(pk_all)))

# ----- 4. Export per-peak table -----
cat("Exporting per-peak table...\n")
df <- as.data.frame(pk_all)
df$gn <- sub(";.*", "", df$gene)  # first gene from semicolon-separated list
write.csv(df, file.path(OUT_DIR, "F24_peak_categories.csv"), row.names = FALSE)

cat("\n=== Differential category counts ===\n")
diff_tab <- as.data.table(table(diff = df$diff))
diff_tab[, pct := round(100 * N / sum(N), 1)]
print(diff_tab)

# ----- 5. Join RNA-seq (GSE90894) -----
cat("\nJoining RNA-seq RPKM table...\n")
rnaseq <- read_xlsx(RNASEQ_XLSX, sheet = 1, skip = 4)
cat(sprintf("  RNA-seq table: %s genes x %s columns (columns: %s)\n",
            format(nrow(rnaseq), big.mark = ","), ncol(rnaseq),
            paste(head(colnames(rnaseq), 6), collapse = ", ")))

df2 <- df %>%
  left_join(rnaseq, by = c("gn" = "name"))
df2$lfc_rnaseq <- log2((df2$ESCs + 0.001) / (df2$MEFs + 0.001))
df2 <- df2[!is.na(df2$lfc_rnaseq), ]
df2 <- df2[df2$diff != "low", ]
df2$diff <- droplevels(df2$diff)
write.csv(df2, file.path(OUT_DIR, "F24_peak_categories_rnaseq.csv"), row.names = FALSE)
cat(sprintf("  Peaks with matching RNA-seq (excl. low): %s (%.1f%% of annotated)\n",
            format(nrow(df2), big.mark = ","),
            100 * nrow(df2) / sum(!is.na(df$gene))))

# ----- 6. Theme -----
theme_fig1 <- theme_classic() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black")
  )

diff_colors <- c(
  "unchanged" = "#bdbdbd",
  "MEF"       = "#1f77b4",
  "ESC"       = "#d62728"
)

# ----- 6b. Promoter overlap stacked bar (RNA-seq matched peaks, direct TSS) -----
cat("\nComputing promoter overlap by diff category (RNA-seq matched peaks, direct TSS)...\n")
tss_gr2 <- rtracklayer::import(REFGENE_TSS)
df2_gr <- GRanges(seqnames = df2$seqnames,
                  ranges = IRanges(start = df2$start, end = df2$end))
ov2 <- findOverlaps(df2_gr, tss_gr2, ignore.strand = TRUE)
df2$has_tss <- FALSE
df2$has_tss[queryHits(ov2)] <- TRUE

promoter_dt <- data.table(
  diff = rep(c("unchanged", "MEF", "ESC"), each = 2),
  status = rep(c("Promoter-proximal", "Non-promoter"), 3),
  n = c(
    sum(df2$diff == "unchanged" & df2$has_tss),
    sum(df2$diff == "unchanged" & !df2$has_tss),
    sum(df2$diff == "MEF" & df2$has_tss),
    sum(df2$diff == "MEF" & !df2$has_tss),
    sum(df2$diff == "ESC" & df2$has_tss),
    sum(df2$diff == "ESC" & !df2$has_tss)
  )
)
promoter_dt[, pct := 100 * n / sum(n), by = diff]
fwrite(promoter_dt, file.path(OUT_DIR, "F24_promoter_overlap.csv"))

cat("  Promoter overlap (RNA-seq matched peaks, direct TSS):\n")
print(promoter_dt)

promoter_dt$status <- factor(promoter_dt$status,
                             levels = c("Non-promoter", "Promoter-proximal"))
p_promoter <- ggplot(promoter_dt, aes(x = diff, y = pct, fill = status)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
            position = position_stack(vjust = 0.5),
            size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Promoter-proximal" = "#2ca02c",
                                "Non-promoter" = "#bdbdbd")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "G4 peak differential category",
       y = "Fraction of peaks (%)",
       title = "RefSeq promoter overlap",
       fill = NULL) +
  theme_fig1 +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))

ggsave(file.path(OUT_DIR, "panel_F24_promoter_overlap.pdf"),
       p_promoter, width = 7, height = 6, device = "pdf")
ggsave(file.path(OUT_DIR, "panel_F24_promoter_overlap.png"),
       p_promoter, width = 7, height = 6, dpi = 300)

# ----- 7. Violin plot (log2 ESC/MEF) -----
cat("Plotting violin...\n")
p <- ggviolin(
  df2, "diff", "lfc_rnaseq",
  fill = "diff",
  add = "mean_sd",
  add.params = list(size = 0.6)
) +
  scale_fill_manual(values = diff_colors) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  labs(
    x = "G4 peak differential category",
    y = "RNA-seq log2(ESC RPKM / MEF RPKM)",
    title = "RNA-seq expression of scG4 peaks by differential category"
  ) +
  theme_fig1 +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

ggsave(file.path(OUT_DIR, "panel_F24_rnaseq_violin.pdf"), p,
       width = 8, height = 6, device = "pdf")
ggsave(file.path(OUT_DIR, "panel_F24_rnaseq_violin.png"), p,
       width = 8, height = 6, dpi = 300)

# ----- 8. Split violin: MEF and ESC RPKM per category + paired t-test -----
cat("Plotting split violin (MEF/ESC RPKM by category, log10)...\n")
df2_long <- tidyr::pivot_longer(
  df2[, c("seqnames", "start", "end", "diff", "MEFs", "ESCs", "gn")],
  cols = c(MEFs, ESCs),
  names_to = "sample", values_to = "rpkm"
)
df2_long$sample <- factor(df2_long$sample,
                          levels = c("MEFs", "ESCs"),
                          labels = c("MEF", "ESC"))

# Paired t-tests (MEF vs ESC RPKM within each category)
ttests <- do.call(rbind, lapply(levels(df2$diff), function(d) {
  sub <- df2[df2$diff == d, ]
  tt <- t.test(sub$MEFs, sub$ESCs, paired = TRUE)
  p_lab <- ifelse(tt$p.value < 0.001, sprintf("p=%.1e", tt$p.value),
                  sprintf("p=%.3f", tt$p.value))
  data.frame(diff = d, p_value = tt$p.value, p_label = p_lab,
             mean_MEF = mean(sub$MEFs), mean_ESC = mean(sub$ESCs),
             stringsAsFactors = FALSE)
}))
cat("  Paired t-tests (MEF vs ESC RPKM):\n")
print(ttests[, c("diff", "p_label", "mean_MEF", "mean_ESC")])
fwrite(ttests, file.path(OUT_DIR, "F24_paired_ttest_rnaseq.csv"))

y_pos <- 10^(log10(max(df2_long$rpkm, na.rm = TRUE)) * 1.05)

p_split <- ggplot(df2_long, aes(x = diff, y = rpkm, fill = sample)) +
  geom_violin(position = position_dodge(width = 0.9), alpha = 0.75,
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_boxplot(width = 0.15, position = position_dodge(width = 0.9),
               alpha = 0.8, outlier.shape = NA) +
  scale_y_log10() +
  scale_fill_manual(values = c("MEF" = "#1f77b4", "ESC" = "#d62728")) +
  annotate("text", x = seq_along(levels(df2$diff)), y = y_pos,
           label = ttests$p_label, size = 4.5, fontface = "italic") +
  labs(x = "G4 peak differential category",
       y = "RNA-seq RPKM (log10)",
       title = "MEF and ESC RNA-seq RPKM by scG4 peak category",
       fill = "Cell type") +
  theme_fig1 +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top")

ggsave(file.path(OUT_DIR, "panel_F24_split_violin_rpkms.pdf"), p_split,
       width = 8, height = 6, device = "pdf")
ggsave(file.path(OUT_DIR, "panel_F24_split_violin_rpkms.png"), p_split,
       width = 8, height = 6, dpi = 300)

# ----- 9. Summary stats -----
stats <- df2 %>%
  group_by(diff) %>%
  summarise(
    n_peaks = n(),
    median_lfc_rnaseq = round(median(lfc_rnaseq), 3),
    mean_lfc_rnaseq   = round(mean(lfc_rnaseq), 3),
    q25_lfc_rnaseq    = round(as.numeric(quantile(lfc_rnaseq, 0.25)), 3),
    q75_lfc_rnaseq    = round(as.numeric(quantile(lfc_rnaseq, 0.75)), 3)
  )
write.csv(stats, file.path(OUT_DIR, "rnaseq_violin_stats.csv"), row.names = FALSE)

cat("\nSummary (RNA-seq log2 ESC/MEF by G4 peak category):\n")
print(as.data.frame(stats))
cat("\nWrote outputs to", OUT_DIR, "\n")
