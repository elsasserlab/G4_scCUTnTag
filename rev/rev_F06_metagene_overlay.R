#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Alternative visualizations of the G4-vs-expression relationship.
#
# Produces 5 plots, each emphasizing a different aspect of the same data:
#   1. Boxplot per expression decile (distribution, not just mean)
#   2. Hex-binned scatter of every gene with LOESS smoother
#   3. Both clusters overlaid on a single line+ribbon plot
#   4. Metagene profile around TSS, stratified by expression quintile
#   5. Heatmap of genes (sorted by expression) x TSS-relative position
#
# Reuses the same inputs as 07_G4_vs_expression.R.
# ---------------------------------------------------------------------------


source("rev/paths.R")
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(data.table)
  library(ggplot2)
  library(scales)
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

# For metagene / heatmap (visualizations 4 and 5)
META_UP    <- 3000        # window upstream of TSS
META_DN    <- 3000        # window downstream of TSS
META_BIN   <- 100         # bin size in bp (60 bins for a 6 kb window)
N_QUANTILE <- 5           # 5 expression bins for the metagene curves

# ===== 1. Build per-gene TSS GRanges (same as 07_) =====
cat("Loading GTF and building per-gene TSS...\n")
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
tx_dt[, tss := ifelse(strand == "+", start, end)]
gene_tss <- tx_dt[, .(
  chrom = data.table::first(chrom),
  strand = data.table::first(strand),
  tss = ifelse(data.table::first(strand) == "+", min(tss), max(tss)),
  gene_name = data.table::first(gene_name)
), by = gene_id]

tss_gr <- GRanges(
  seqnames = gene_tss$chrom,
  ranges   = IRanges(start = gene_tss$tss, end = gene_tss$tss),
  strand   = gene_tss$strand,
  gene_id  = gene_tss$gene_id
)
cat(sprintf("  %d TSS\n", length(tss_gr)))

# ===== 2. Reuse the previous joined table for whole-promoter means =====
# (We re-extract here in case the user wants this script to stand alone)
get_mean_signal <- function(bw_path, gr) {
  cat(sprintf("    %s...\n", basename(bw_path)))
  bw <- import(bw_path, which = gr, as = "RleList")
  signals <- numeric(length(gr))
  for (chrom in unique(as.character(seqnames(gr)))) {
    idx <- which(as.character(seqnames(gr)) == chrom)
    if (!(chrom %in% names(bw))) { signals[idx] <- 0; next }
    v <- Views(bw[[chrom]], ranges(gr)[idx])
    signals[idx] <- viewMeans(v)
  }
  signals
}

promoter_gr <- GRanges(
  seqnames = gene_tss$chrom,
  ranges   = IRanges(
    start = pmax(1L, gene_tss$tss - PROMOTER_UP),
    end   = gene_tss$tss + PROMOTER_DN
  ),
  strand   = gene_tss$strand,
  gene_id  = gene_tss$gene_id
)

cat("Computing mean G4 signal per promoter (TSS ±2kb)...\n")
c0_signal <- get_mean_signal(BW_CLUSTER0, promoter_gr)
c1_signal <- get_mean_signal(BW_CLUSTER1, promoter_gr)

# Load TPM
tpm <- fread(TPM_PATH)
mESC_cols <- grep("mESC", colnames(tpm), value = TRUE)
MEF_cols  <- grep("MEF",  colnames(tpm), value = TRUE)
tpm[, mESC_TPM := rowMeans(.SD), .SDcols = mESC_cols]
tpm[, MEF_TPM  := rowMeans(.SD), .SDcols = MEF_cols]
strip_ver <- function(x) sub("\\.[0-9]+$", "", x)
tpm[, gene_id_base := strip_ver(gene_id)]

g4 <- data.table(
  gene_id      = promoter_gr$gene_id,
  gene_id_base = strip_ver(promoter_gr$gene_id),
  c0_G4 = c0_signal, c1_G4 = c1_signal
)
joined <- merge(g4, tpm[, .(gene_id_base, mESC_TPM, MEF_TPM)],
                by = "gene_id_base")
cat(sprintf("  joined: %d genes\n", nrow(joined)))

# Rank-decile helper
ntile <- function(x, n = 10) {
  r <- rank(x, ties.method = "first")
  factor(pmin(n, ceiling(r / length(x) * n)), levels = 1:n)
}

# (Plot 1 boxplot-per-decile and Plot 2 hex-binned scatter removed —
#  not used in the slides. Only Plot 3 overlaid lines (F07) and Plot 4
#  metagene profile (F05) are generated. box_data is still computed below
#  because Plot 3 reuses it.)
box_data <- rbindlist(list(
  joined[, .(cluster = "Cluster 0 (MEF)",
             decile  = ntile(MEF_TPM, 10),
             G4      = c0_G4)],
  joined[, .(cluster = "Cluster 1 (mESC)",
             decile  = ntile(mESC_TPM, 10),
             G4      = c1_G4)]
))

# ===== Visualization 3: Both clusters overlaid (line + ribbon) =====
cat("\n[Plot 3] Both clusters overlaid as lines with SEM ribbons...\n")

dec_overlay <- box_data[, .(
  mean_G4   = mean(G4),
  median_G4 = median(G4),
  sem       = sd(G4) / sqrt(.N),
  n         = .N
), by = .(cluster, decile)]
dec_overlay$decile <- as.integer(as.character(dec_overlay$decile))

p3 <- ggplot(dec_overlay, aes(x = decile, y = mean_G4, color = cluster, fill = cluster)) +
  geom_ribbon(aes(ymin = mean_G4 - sem, ymax = mean_G4 + sem),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.4) +
  scale_color_manual(values = c("Cluster 0 (MEF)"  = "#2ca02c",
                                 "Cluster 1 (mESC)" = "#1f77b4"),
                     name = NULL) +
  scale_fill_manual(values = c("Cluster 0 (MEF)"  = "#2ca02c",
                                "Cluster 1 (mESC)" = "#1f77b4"),
                    name = NULL) +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "Matched-cell-type expression decile",
       y = "Mean G4 signal at promoter (RPGC)",
       title = "G4–expression relationship is identical in both clusters",
       subtitle = "Mean signal per decile, SEM ribbon, both clusters overlaid") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = c(0.2, 0.85),
        legend.background = element_rect(fill = "white", color = "grey80"))
ggsave(file.path(OUT_DIR, "plot3_overlaid_lines.pdf"), p3, width = 7, height = 5)
cat("  saved plot3_overlaid_lines\n")

# ===== Visualizations 4 and 5: metagene + heatmap around TSS =====
# Need a per-gene signal matrix (rows = genes, cols = bins around TSS).
cat("\nExtracting signal matrix around TSS for metagene/heatmap...\n")

extract_tss_matrix <- function(bw_path, tss_gr,
                                up = META_UP, dn = META_DN, bin = META_BIN) {
  cat(sprintf("  loading %s (this takes a minute)...\n", basename(bw_path)))
  bw <- import(bw_path, as = "RleList")
  n_genes <- length(tss_gr)
  n_bins  <- (up + dn) / bin
  mat <- matrix(0, nrow = n_genes, ncol = n_bins)
  bin_offsets_plus  <- seq(-up, dn - bin, by = bin)   # bin starts for + strand
  bin_offsets_minus <- seq(dn - bin, -up, by = -bin)  # bin starts for - strand (reversed)
  for (chrom in unique(as.character(seqnames(tss_gr)))) {
    if (!(chrom %in% names(bw))) next
    idx <- which(as.character(seqnames(tss_gr)) == chrom)
    if (length(idx) == 0) next
    chrom_cov <- bw[[chrom]]
    chrom_len <- length(chrom_cov)
    tss_v    <- start(tss_gr)[idx]
    strand_v <- as.character(strand(tss_gr)[idx])
    for (j in seq_along(idx)) {
      tss <- tss_v[j]
      if (strand_v[j] == "+") {
        s <- tss + bin_offsets_plus
      } else {
        s <- tss + bin_offsets_minus
      }
      e <- s + bin - 1
      valid <- s >= 1 & e <= chrom_len
      if (any(valid)) {
        v <- Views(chrom_cov, start = s[valid], end = e[valid])
        mat[idx[j], which(valid)] <- viewMeans(v)
      }
    }
  }
  mat
}

mat_c0 <- extract_tss_matrix(BW_CLUSTER0, tss_gr)
mat_c1 <- extract_tss_matrix(BW_CLUSTER1, tss_gr)

# Align signal matrices to the joined gene table (by gene_id)
align_matrix <- function(mat, tss_gr, joined_dt) {
  rownames(mat) <- tss_gr$gene_id
  m <- mat[joined_dt$gene_id, , drop = FALSE]
  m
}
mat_c0_j <- align_matrix(mat_c0, tss_gr, joined)
mat_c1_j <- align_matrix(mat_c1, tss_gr, joined)

# ===== Visualization 4: Metagene profile, stratified by expression quintile =====
cat("\n[Plot 4] Metagene profile, stratified by expression quintile...\n")

bin_centers <- seq(-META_UP + META_BIN/2, META_DN - META_BIN/2, by = META_BIN)

meta_for_cluster <- function(mat, expr_vec, cluster_label) {
  q <- ntile(expr_vec, N_QUANTILE)
  out <- rbindlist(lapply(seq_len(N_QUANTILE), function(k) {
    rows <- which(as.integer(as.character(q)) == k)
    if (length(rows) == 0) return(NULL)
    means <- colMeans(mat[rows, , drop = FALSE], na.rm = TRUE)
    data.table(
      cluster = cluster_label,
      quintile = paste0("Q", k),
      pos = bin_centers,
      mean_G4 = means
    )
  }))
  out
}

meta_data <- rbindlist(list(
  meta_for_cluster(mat_c0_j, joined$MEF_TPM,  "Cluster 0 (MEF)"),
  meta_for_cluster(mat_c1_j, joined$mESC_TPM, "Cluster 1 (mESC)")
))
fwrite(meta_data, file.path(OUT_DIR, "metagene_data.csv"))

p4 <- ggplot(meta_data, aes(x = pos / 1000, y = mean_G4, color = quintile)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1) +
  facet_wrap(~ cluster, ncol = 2, scales = "free_y") +
  scale_color_manual(
    values = c("Q1" = "#3b6fac", "Q2" = "#7cb3d2", "Q3" = "#cccccc",
               "Q4" = "#e89c81", "Q5" = "#c73d3d"),
    labels = c("Q1 (lowest expr)", "Q2", "Q3", "Q4", "Q5 (highest expr)"),
    name   = "Expression quintile"
  ) +
  labs(x = "Distance from TSS (kb)",
       y = "Mean G4 signal (RPGC)",
       title = "Metagene profile of G4 signal around TSS",
       subtitle = "Genes split into expression quintiles in matched cell type") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold"),
        legend.position = "top")
ggsave(file.path(OUT_DIR, "plot4_metagene_profile.pdf"), p4, width = 12, height = 5)
cat("  saved plot4_metagene_profile\n")

# (Plot 5 expression-sorted heatmap removed — not used in the slides.)

cat("\nDone. Outputs in:", OUT_DIR, "\n")
