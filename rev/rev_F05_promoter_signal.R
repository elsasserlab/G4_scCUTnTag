#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# G4 vs gene expression correlation per scCUT&Tag cluster.
#
# Per cluster (0 = MEF, 1 = mESC per Fig 1G):
#   - Compute mean G4 BigWig signal in the promoter (TSS ± 2 kb) of each gene.
#   - Pair with gene expression (cluster 0 -> MEF RNA-seq; cluster 1 -> mESC).
#
# Outputs:
#   - Gene-level data table (G4 signal x expression per cluster x cell type)
#   - Decile plot: G4 signal binned by expression decile, per cluster
#   - Cross-correlation table: cluster G4 x cell-type expression (4 cells)
#   - Heatmap of the cross-correlation matrix
#
# Inputs assumed in DATA_DIR:
#   - gene_level_TPM.tsv
#   - mm10_.annotation.gtf.gz
#   - scG4CnT_pseudobulk_cluster0_RPGC.bw
#   - scG4CnT_pseudobulk_cluster1_RPGC.bw
# ---------------------------------------------------------------------------


source("rev/paths.R")
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# ===== CONFIG =====
inp <- function(fn) {
  hits <- list.files(DATA, pattern = fn, recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("File not found in data/: ", fn)
  hits[1]
}
GTF_PATH <- file.path(GENOME_DIR, "mm10_.annotation.gtf.gz")
TPM_PATH <- file.path(GENOME_DIR, "gene_level_TPM.tsv")
BW_CLUSTER0 <- BW_CL0
BW_CLUSTER1 <- BW_CL1
OUT_DIR  <- file.path(OUT_ROOT, "G4_vs_expression")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CANONICAL   <- c(paste0("chr", 1:19), "chrX", "chrY")
PROMOTER_UP <- 2000
PROMOTER_DN <- 2000

# Cluster -> cell type mapping (from manuscript Fig 1G correlation matrix)
CLUSTER_TO_CELL <- c("cluster0" = "MEF", "cluster1" = "mESC")

# ===== 1. Build per-gene TSS GRanges =====
cat("Loading GTF and extracting per-gene TSS...\n")
gtf <- rtracklayer::import(GTF_PATH)
gtf <- gtf[as.character(seqnames(gtf)) %in% CANONICAL]
seqlevels(gtf) <- intersect(seqlevels(gtf), CANONICAL)

tx <- gtf[gtf$type == "transcript"]
tx_dt <- data.table(
  chrom     = as.character(seqnames(tx)),
  strand    = as.character(strand(tx)),
  start     = start(tx),
  end       = end(tx),
  gene_id   = tx$gene_id,
  gene_name = tx$gene_name
)
# Per-transcript TSS (+ strand: start, - strand: end)
tx_dt[, tss := ifelse(strand == "+", start, end)]

# Per gene: take the most upstream TSS as the gene's primary TSS
gene_tss <- tx_dt[, .(
  chrom = data.table::first(chrom),
  strand = data.table::first(strand),
  tss = ifelse(data.table::first(strand) == "+", min(tss), max(tss)),
  gene_name = data.table::first(gene_name)
), by = gene_id]
cat(sprintf("  %d genes with a TSS\n", nrow(gene_tss)))

# Promoter window
promoter_gr <- GRanges(
  seqnames = gene_tss$chrom,
  ranges   = IRanges(
    start = pmax(1L, gene_tss$tss - PROMOTER_UP),
    end   = gene_tss$tss + PROMOTER_DN
  ),
  strand   = gene_tss$strand,
  gene_id  = gene_tss$gene_id,
  gene_name = gene_tss$gene_name
)
cat(sprintf("  promoter ranges built: %d\n", length(promoter_gr)))

# ===== 2. Extract mean signal per promoter from each BigWig =====
get_mean_signal <- function(bw_path, gr) {
  cat(sprintf("    reading %s...\n", basename(bw_path)))
  # Load BigWig as RleList over the chromosomes we need
  bw <- import(bw_path, which = gr, as = "RleList")
  # Per-chromosome views -> mean
  signals <- numeric(length(gr))
  for (chrom in unique(as.character(seqnames(gr)))) {
    idx <- which(as.character(seqnames(gr)) == chrom)
    if (!(chrom %in% names(bw))) {
      signals[idx] <- 0
      next
    }
    sub_ranges <- ranges(gr)[idx]
    v <- Views(bw[[chrom]], sub_ranges)
    signals[idx] <- viewMeans(v)
  }
  signals
}

cat("\nExtracting G4 signal at promoters from BigWigs...\n")
cluster0_signal <- get_mean_signal(BW_CLUSTER0, promoter_gr)
cluster1_signal <- get_mean_signal(BW_CLUSTER1, promoter_gr)
cat(sprintf("  cluster0 mean: %.3f  median: %.3f  max: %.1f\n",
            mean(cluster0_signal), median(cluster0_signal), max(cluster0_signal)))
cat(sprintf("  cluster1 mean: %.3f  median: %.3f  max: %.1f\n",
            mean(cluster1_signal), median(cluster1_signal), max(cluster1_signal)))

# ===== 3. Read TPM and average replicates per cell type =====
cat("\nLoading TPM matrix and averaging replicates per cell type...\n")
tpm <- fread(TPM_PATH)
mESC_cols <- grep("mESC", colnames(tpm), value = TRUE)
MEF_cols  <- grep("MEF",  colnames(tpm), value = TRUE)
cat(sprintf("  mESC columns: %s\n", paste(mESC_cols, collapse = ", ")))
cat(sprintf("  MEF columns:  %s\n", paste(MEF_cols, collapse = ", ")))

tpm[, mESC_TPM := rowMeans(.SD, na.rm = TRUE), .SDcols = mESC_cols]
tpm[, MEF_TPM  := rowMeans(.SD, na.rm = TRUE), .SDcols = MEF_cols]

# Strip Ensembl version suffix from gene_id for joining (TPM gene_ids are like
# ENSMUSG00000000001.4; GTF gene_ids may or may not have version)
strip_ver <- function(x) sub("\\.[0-9]+$", "", x)

# ===== 4. Join G4 signal with expression =====
g4 <- data.table(
  gene_id      = promoter_gr$gene_id,
  gene_id_base = strip_ver(promoter_gr$gene_id),
  gene_name    = promoter_gr$gene_name,
  cluster0_G4  = cluster0_signal,
  cluster1_G4  = cluster1_signal
)
tpm_join <- tpm[, .(gene_id, gene_id_base = strip_ver(gene_id),
                    gene_name, mESC_TPM, MEF_TPM)]

joined <- merge(g4, tpm_join[, .(gene_id_base, mESC_TPM, MEF_TPM)],
                by = "gene_id_base", all.x = FALSE, all.y = FALSE)
cat(sprintf("\nJoined data: %d genes with both G4 signal and TPM\n", nrow(joined)))

fwrite(joined, file.path(OUT_DIR, "gene_level_G4_vs_expression.tsv"), sep = "\t")

# ===== 5. Decile-binned signal: G4 vs matched-cell-type expression =====
cat("\nComputing decile-binned G4 signal by expression rank...\n")
joined[, log_c0_G4   := log2(cluster0_G4 + 1)]
joined[, log_c1_G4   := log2(cluster1_G4 + 1)]
joined[, log_mESC_TPM := log2(mESC_TPM + 1)]
joined[, log_MEF_TPM  := log2(MEF_TPM + 1)]

# Rank-based decile assignment (handles ties / many-zero values gracefully).
ntile <- function(x, n = 10) {
  r <- rank(x, ties.method = "first")
  factor(pmin(n, ceiling(r / length(x) * n)), levels = 1:n)
}

decile_data <- rbindlist(list(
  # Cluster 0 (MEF) vs MEF expression deciles  — MATCHED
  joined[, .(cluster = "Cluster 0 (MEF)", match = "matched",
             expr_decile = ntile(MEF_TPM, 10),
             G4 = cluster0_G4)],
  # Cluster 1 (mESC) vs mESC expression deciles  — MATCHED
  joined[, .(cluster = "Cluster 1 (mESC)", match = "matched",
             expr_decile = ntile(mESC_TPM, 10),
             G4 = cluster1_G4)]
))
decile_data <- decile_data[!is.na(expr_decile)]
decile_summary <- decile_data[, .(
  n         = .N,
  mean_G4   = mean(G4),
  median_G4 = median(G4),
  sem_G4    = sd(G4) / sqrt(.N)
), by = .(cluster, expr_decile)]
decile_summary$expr_decile <- as.integer(as.character(decile_summary$expr_decile))
fwrite(decile_summary, file.path(OUT_DIR, "decile_summary.csv"))
print(decile_summary)

# ===== 6. Cross-correlation table =====
cat("\nComputing G4 x expression cross-correlation (Spearman, log-transformed)...\n")

corr_mat <- matrix(NA_real_, 2, 2,
                   dimnames = list(
                     c("Cluster 0 G4", "Cluster 1 G4"),
                     c("MEF expression", "mESC expression")))

corr_mat["Cluster 0 G4", "MEF expression"]  <- cor(joined$log_c0_G4, joined$log_MEF_TPM,  method = "spearman")
corr_mat["Cluster 0 G4", "mESC expression"] <- cor(joined$log_c0_G4, joined$log_mESC_TPM, method = "spearman")
corr_mat["Cluster 1 G4", "MEF expression"]  <- cor(joined$log_c1_G4, joined$log_MEF_TPM,  method = "spearman")
corr_mat["Cluster 1 G4", "mESC expression"] <- cor(joined$log_c1_G4, joined$log_mESC_TPM, method = "spearman")

# Also do Pearson on log-transformed values for completeness
pear_mat <- corr_mat
pear_mat["Cluster 0 G4", "MEF expression"]  <- cor(joined$log_c0_G4, joined$log_MEF_TPM,  method = "pearson")
pear_mat["Cluster 0 G4", "mESC expression"] <- cor(joined$log_c0_G4, joined$log_mESC_TPM, method = "pearson")
pear_mat["Cluster 1 G4", "MEF expression"]  <- cor(joined$log_c1_G4, joined$log_MEF_TPM,  method = "pearson")
pear_mat["Cluster 1 G4", "mESC expression"] <- cor(joined$log_c1_G4, joined$log_mESC_TPM, method = "pearson")

cat("\nSpearman correlation (matched pairs on diagonal of the meaningful pattern):\n")
print(round(corr_mat, 3))
cat("\nPearson correlation (log2):\n")
print(round(pear_mat, 3))

corr_dt <- data.table(
  comparison = c("Cluster 0 G4 vs MEF expr",  "Cluster 0 G4 vs mESC expr",
                 "Cluster 1 G4 vs MEF expr",  "Cluster 1 G4 vs mESC expr"),
  match_type = c("MATCHED", "mismatched", "mismatched", "MATCHED"),
  spearman_r = c(corr_mat["Cluster 0 G4", "MEF expression"],
                 corr_mat["Cluster 0 G4", "mESC expression"],
                 corr_mat["Cluster 1 G4", "MEF expression"],
                 corr_mat["Cluster 1 G4", "mESC expression"]),
  pearson_r  = c(pear_mat["Cluster 0 G4", "MEF expression"],
                 pear_mat["Cluster 0 G4", "mESC expression"],
                 pear_mat["Cluster 1 G4", "MEF expression"],
                 pear_mat["Cluster 1 G4", "mESC expression"])
)
fwrite(corr_dt, file.path(OUT_DIR, "cross_correlation_table.csv"))
print(corr_dt)

# ===== 7. PLOTS =====
cat("\nDrawing figures...\n")

# 7a. Decile plot — matched cluster vs cell-type expression
p_dec <- ggplot(decile_summary,
                aes(x = factor(expr_decile), y = mean_G4, fill = cluster)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3,
           width = 0.7, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_G4 - sem_G4, ymax = mean_G4 + sem_G4),
                width = 0.25, linewidth = 0.4) +
  geom_text(aes(label = scales::comma(n)),
            vjust = -2.0, size = 2.4, color = "grey40") +
  facet_wrap(~ cluster, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("Cluster 0 (MEF)" = "#2ca02c",
                                "Cluster 1 (mESC)" = "#1f77b4")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Matched-cell-type expression decile (1 = lowest, 10 = highest)",
       y = "Mean G4 signal at promoter (RPGC)",
       title = "G4 promoter signal tracks with gene expression, per cluster",
       subtitle = "Cluster 0 vs MEF RNA-seq; Cluster 1 vs mESC RNA-seq (matched)") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "decile_G4_by_expression.pdf"),
       p_dec, width = 11, height = 5)
cat("  saved decile_G4_by_expression.pdf\n")

# (Cross-correlation heatmap and combined decile+heatmap panel removed —
#  not used in the slides. Only the decile plot (F06) is generated here;
#  the cross_correlation_table.csv is still written as a supporting table.)

cat("\nDone. Outputs in:", OUT_DIR, "\n")
