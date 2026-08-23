#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Stratified recovery of bulk G4 peaks in scCUT&Tag pseudobulk.
#
# Pipeline:
#   1. Load peak files (2 bulk reps + 1 pseudobulk per cell type).
#   2. Build bulk consensus via bedscout::loci_consensus(..., min_consensus = 2).
#   3. Re-attach signalValue to consensus (loci_consensus does not preserve it).
#   4. For each consensus peak, check overlap with pseudobulk.
#   5. Bin consensus peaks by signal decile.
#   6. Plot recovery rate per decile, per cell type.
#
# Install bedscout: devtools::install_github("cnluzon/bedscout")
# ---------------------------------------------------------------------------


source("rev/paths.R")
suppressPackageStartupMessages({
  library(bedscout)
  library(GenomicRanges)
  library(rtracklayer)
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
OUT_DIR  <- file.path(OUT_ROOT, "stratified_recovery_bedscout")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CANONICAL <- c(paste0("chr", 1:19), "chrX", "chrY")

peak_files <- c(
  mESC_bulk_rep1   = file.path(BULK_PEAKS, "GSM8836082_bulkG4CnT_mESC_rep1.broadPeak"),
  mESC_bulk_rep2   = file.path(BULK_PEAKS, "GSM8836083_bulkG4CnT_mESC_rep2.broadPeak"),
  MEF_bulk_rep1    = file.path(BULK_PEAKS, "GSM8836084_bulkG4CnT_3T3_rep1.broadPeak"),
  MEF_bulk_rep2    = file.path(BULK_PEAKS, "GSM8836085_bulkG4CnT_3T3_rep2.broadPeak"),
  mESC_pseudobulk  = CLUSTER_PEAKS_1,
  MEF_pseudobulk   = CLUSTER_PEAKS_0
)

# ===== 1. LOAD PEAKS =====
# rtracklayer::import() auto-handles broadPeak/narrowPeak and preserves
# signalValue, which we need for the decile ranking.
cat("Loading peaks via rtracklayer::import()...\n")
peaks <- lapply(peak_files, function(f) {
  fmt <- if (grepl("broadPeak$", f)) "broadPeak" else "narrowPeak"
  gr  <- rtracklayer::import(f, format = fmt)
  gr  <- gr[as.character(seqnames(gr)) %in% CANONICAL]
  gr$signal <- gr$signalValue
  sort(gr)
})
for (nm in names(peaks)) cat(sprintf("  %-18s %d peaks\n", nm, length(peaks[[nm]])))

# ===== 2. BULK CONSENSUS (bedscout) =====
# loci_consensus returns merged consensus intervals but does NOT preserve signal,
# so we re-attach by taking the max signal of contributing peaks from either rep.
attach_signal <- function(consensus_gr, rep1, rep2) {
  combined <- c(rep1, rep2)
  ov <- findOverlaps(consensus_gr, combined, ignore.strand = TRUE)
  agg <- tapply(combined$signal[subjectHits(ov)],
                queryHits(ov), max, default = NA_real_)
  out_signal <- rep(NA_real_, length(consensus_gr))
  out_signal[as.integer(names(agg))] <- as.numeric(agg)
  consensus_gr$signal <- out_signal
  consensus_gr
}

cat("\nBuilding bulk consensus via bedscout::loci_consensus(..., min_consensus = 2)...\n")
consensus <- list()
for (ct in c("mESC", "MEF")) {
  r1 <- peaks[[paste0(ct, "_bulk_rep1")]]
  r2 <- peaks[[paste0(ct, "_bulk_rep2")]]
  cons <- loci_consensus(list(r1, r2), min_consensus = 2)
  cons <- attach_signal(cons, r1, r2)
  consensus[[ct]] <- cons
  cat(sprintf("  %s consensus: %d peaks (rep1=%d, rep2=%d)\n",
              ct, length(cons), length(r1), length(r2)))
}

# ===== 3. OVERLAP STATS =====
cat("\nOverlap statistics (bedscout::jaccard_index + total_bp_overlap)...\n")
stats <- data.table()
for (ct in c("mESC", "MEF")) {
  bulk <- consensus[[ct]]
  sc   <- peaks[[paste0(ct, "_pseudobulk")]]

  bulk_recov  <- overlapsAny(bulk, sc, ignore.strand = TRUE)
  sc_in_bulk  <- overlapsAny(sc, bulk, ignore.strand = TRUE)
  jacc        <- jaccard_index(bulk, sc, ignore.strand = TRUE)
  bp_inter    <- total_bp_overlap(bulk, sc, ignore.strand = TRUE)

  stats <- rbind(stats, data.table(
    cell_type        = ct,
    bulk_consensus   = length(bulk),
    pseudobulk       = length(sc),
    bulk_recovered   = sum(bulk_recov),
    bulk_recovery    = 100 * mean(bulk_recov),
    sc_in_bulk       = sum(sc_in_bulk),
    sc_precision     = 100 * mean(sc_in_bulk),
    bp_intersection  = bp_inter,
    jaccard          = jacc
  ))
}
fwrite(stats, file.path(OUT_DIR, "overlap_stats.csv"))
print(stats)

# ===== 4. STRATIFIED RECOVERY =====
cat("\nStratified recovery by bulk signal decile...\n")
strat <- data.table()
for (ct in c("mESC", "MEF")) {
  bulk <- consensus[[ct]]
  sc   <- peaks[[paste0(ct, "_pseudobulk")]]
  bulk_recov <- overlapsAny(bulk, sc, ignore.strand = TRUE)

  dt <- data.table(
    cell_type = ct,
    signal    = bulk$signal,
    recovered = bulk_recov
  )
  dt[, decile := cut(signal,
                     breaks = quantile(signal, probs = seq(0, 1, 0.1),
                                       na.rm = TRUE),
                     include.lowest = TRUE, labels = 1:10)]

  agg <- dt[!is.na(decile),
            .(n            = .N,
              recovered    = sum(recovered),
              recovery_pct = 100 * mean(recovered),
              mean_signal  = mean(signal),
              min_signal   = min(signal),
              max_signal   = max(signal)),
            by = .(cell_type, decile)]
  strat <- rbind(strat, agg)
}
strat$decile <- as.integer(as.character(strat$decile))
fwrite(strat, file.path(OUT_DIR, "stratified_recovery.csv"))
print(strat)

# ===== 5. PLOT: STRATIFIED RECOVERY =====
cat("\nDrawing stratified recovery plot...\n")
p_strat <- ggplot(strat, aes(x = factor(decile), y = recovery_pct,
                              fill = cell_type)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3,
           width = 0.75, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%d/%d", recovered, n)),
            vjust = -0.4, size = 2.6) +
  facet_wrap(~ cell_type, ncol = 2) +
  scale_fill_manual(values = c(mESC = "#1f77b4", MEF = "#2ca02c")) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Bulk signal decile (1 = lowest, 10 = highest)",
       y = "Bulk peaks recovered in pseudobulk (%)",
       title = "Recovery of bulk G4 peaks in scCUT&Tag pseudobulk",
       subtitle = "Stratified by bulk consensus signalValue (decile)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "stratified_recovery.pdf"),
       p_strat, width = 10, height = 4.5)
cat("  saved stratified_recovery.pdf\n")

# (Per-cell-type Venn plots removed — not used in the slides. Only the
#  stratified-recovery bar chart (F04) is generated.)

cat("\nDone. Outputs in:", OUT_DIR, "\n")
