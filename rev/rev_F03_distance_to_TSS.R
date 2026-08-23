#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Distance-to-TSS distribution: bulk consensus vs scCUT&Tag pseudobulk.
#
# For each peak, compute signed distance to the nearest TSS (TSS is defined
# per transcript; positive = downstream of TSS on the gene's strand, negative
# = upstream). Plot density curves overlaying bulk vs sc per cell type, and
# also a binned bar plot.
#
# If the curves match in shape, the spatial relationship of G4 peaks to
# genes is identical between bulk and sc — a continuous complement to the
# discrete genomic-feature distribution.
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
GTF_PATH <- file.path(GENOME_DIR, "mm10_.annotation.gtf.gz")
OUT_DIR  <- file.path(OUT_ROOT, "distance_to_TSS")
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

# ===== LOAD PEAKS =====
cat("Loading peaks...\n")
peaks <- lapply(peak_files, function(f) {
  fmt <- if (grepl("broadPeak$", f)) "broadPeak" else "narrowPeak"
  gr  <- rtracklayer::import(f, format = fmt)
  gr  <- gr[as.character(seqnames(gr)) %in% CANONICAL]
  gr$signal <- gr$signalValue
  sort(gr)
})
for (nm in names(peaks)) cat(sprintf("  %-18s %d peaks\n", nm, length(peaks[[nm]])))

# ===== BULK CONSENSUS =====
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

cat("\nBuilding bulk consensus...\n")
consensus <- list()
for (ct in c("mESC", "MEF")) {
  r1 <- peaks[[paste0(ct, "_bulk_rep1")]]
  r2 <- peaks[[paste0(ct, "_bulk_rep2")]]
  cons <- loci_consensus(list(r1, r2), min_consensus = 2)
  cons <- attach_signal(cons, r1, r2)
  consensus[[ct]] <- cons
  cat(sprintf("  %s consensus: %d peaks\n", ct, length(cons)))
}

# ===== TSS GRanges =====
# A TSS is the start of a transcript on the + strand, or the end on the -.
# Build a strand-aware GRanges of single-bp TSS positions, one per transcript.
cat("\nBuilding TSS GRanges from GTF...\n")
gtf <- rtracklayer::import(GTF_PATH)
gtf <- gtf[as.character(seqnames(gtf)) %in% CANONICAL]
seqlevels(gtf) <- intersect(seqlevels(gtf), CANONICAL)

tx <- gtf[gtf$type == "transcript"]
tss_pos <- ifelse(as.character(strand(tx)) == "+", start(tx), end(tx))
tss <- GRanges(seqnames(tx), IRanges(tss_pos, tss_pos), strand = strand(tx))
cat(sprintf("  %d TSS records (one per transcript)\n", length(tss)))

# ===== DISTANCE FROM EACH PEAK TO NEAREST TSS =====
# Sign convention: positive = peak is downstream of TSS on gene strand
# (i.e., inside or past the gene body), negative = upstream of TSS.
# We use peak midpoint vs TSS, and look up strand of the nearest TSS.
compute_signed_distance <- function(peak_gr, tss_gr) {
  mids <- GRanges(seqnames(peak_gr),
                  IRanges(start = (start(peak_gr) + end(peak_gr)) %/% 2, width = 1))
  dn  <- distanceToNearest(mids, tss_gr, ignore.strand = TRUE)
  q   <- queryHits(dn); s <- subjectHits(dn)
  raw <- mcols(dn)$distance
  # Signed: midpoint - TSS, flipped on minus strand
  delta <- start(mids)[q] - start(tss_gr)[s]
  flip  <- as.character(strand(tss_gr))[s] == "-"
  signed <- ifelse(flip, -delta, delta)
  # Make signed magnitude match raw distance (raw is unsigned bp); restore sign
  out <- rep(NA_real_, length(peak_gr))
  out[q] <- sign(signed) * raw
  out
}

cat("\nComputing distance-to-TSS for each peak set...\n")
dist_dt <- rbindlist(lapply(c("mESC", "MEF"), function(ct) {
  bulk <- consensus[[ct]]; sc <- peaks[[paste0(ct, "_pseudobulk")]]
  rbindlist(list(
    data.table(cell_type = ct, source = "Bulk consensus",
               distance = compute_signed_distance(bulk, tss)),
    data.table(cell_type = ct, source = "Pseudobulk (sc)",
               distance = compute_signed_distance(sc, tss))
  ))
}))
dist_dt <- dist_dt[!is.na(distance)]

# Summary stats
cat("\nDistance summary (kb):\n")
print(dist_dt[, .(n          = .N,
                  median_kb  = median(distance) / 1000,
                  mean_kb    = mean(distance) / 1000,
                  pct_within_2kb  = 100 * mean(abs(distance) <= 2000),
                  pct_within_5kb  = 100 * mean(abs(distance) <= 5000),
                  pct_within_10kb = 100 * mean(abs(distance) <= 10000)),
              by = .(cell_type, source)])
fwrite(dist_dt[, .(cell_type, source, distance)],
       file.path(OUT_DIR, "distance_to_TSS_per_peak.csv"))

# ===== PLOTS =====
cat("\nDrawing figures...\n")

dist_dt$cell_type <- factor(dist_dt$cell_type, levels = c("mESC", "MEF"))
dist_dt$source    <- factor(dist_dt$source,
                            levels = c("Bulk consensus", "Pseudobulk (sc)"))

# (Density overlay and zoomed-density plots removed — not used in the slides.
#  Only the binned bar plot (F03) is generated.)

# 3. Binned bar plot in log-distance space
breaks_kb <- c(-Inf, -100, -10, -2, -0.5, 0.5, 2, 10, 100, Inf)
labels_kb <- c("<-100","-100..-10","-10..-2","-2..-0.5","-0.5..0.5",
               "0.5..2","2..10","10..100",">100")
dist_dt[, bin := cut(distance / 1000, breaks = breaks_kb,
                     labels = labels_kb, include.lowest = TRUE)]
binned <- dist_dt[, .N, by = .(cell_type, source, bin)][
  , pct := 100 * N / sum(N), by = .(cell_type, source)]
fwrite(binned, file.path(OUT_DIR, "distance_to_TSS_binned.csv"))

p_bin <- ggplot(binned, aes(x = bin, y = pct, fill = source)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", pct)),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 2.6) +
  facet_wrap(~ cell_type, ncol = 2) +
  scale_fill_manual(values = c("Bulk consensus"  = "#4682b4",
                                "Pseudobulk (sc)" = "#ff9933")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Distance to nearest TSS (kb)",
       y = "Fraction of peaks (%)",
       title = "G4 peak distance to nearest TSS — binned",
       subtitle = "Bulk consensus vs scCUT&Tag pseudobulk",
       fill = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text  = element_text(face = "bold", size = 12),
        plot.title  = element_text(face = "bold"),
        legend.position = "top")
ggsave(file.path(OUT_DIR, "distance_to_TSS_binned.pdf"),
       p_bin, width = 11, height = 5.5)
cat("  saved distance_to_TSS_binned.pdf\n")

cat("\nDone. Outputs in:", OUT_DIR, "\n")
