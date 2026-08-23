# =============================================================
# PQS overlap analysis: PQS predicted G4 sites vs single-cell peak sets.
#
# Datasets analyzed:
#   - scG4 MEF/ESC: cluster 0 (MEF), cluster 1 (mESC)
#   - Mouse brain GFP+: clusters 0, 1, 2
#   - Mouse brain unsorted: clusters 0, 1
#   - scATAC brain: width-adjusted to 270 bp (matching GFP+ scG4 median)
#
# Output:
#   R_outputs/PQS_venn_diagrams/
#     - bar_pct_peaks_overlap_PQS_MEF_ESC_scG4.pdf
#     - bar_pct_peaks_overlap_PQS_mousebrain_GFPpos.pdf
#     - bar_pct_peaks_overlap_PQS_mousebrain_unsorted.pdf
#     - bar_pct_peaks_overlap_PQS_scATAC.pdf
#     - overlap_counts.csv  (all sc datasets)
#     - venn_*.pdf  (individual Venn diagrams)
#     - permTest_summary.csv
# =============================================================


source("rev/paths.R")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)           # overlapsAny lives here
  library(rtracklayer)
  library(VennDiagram)       # draw.pairwise.venn with scaled = FALSE
  library(gridExtra)
  library(grid)
  library(ggplot2)
})
# Silence VennDiagram's chatty log files
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# ----- Paths -----
inp <- function(fn) {
  hits <- list.files(DATA, pattern = fn, recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("File not found in data/: ", fn)
  hits[1]
}
OUT_DIR  <- file.path(OUT_ROOT, "PQS_venn_diagrams")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PQS_BED     <- file.path(DATA, "pqsfinder/PQS_scores.mm10.bed")
PEAK_CLUST0 <- CLUSTER_PEAKS_0
PEAK_CLUST1 <- CLUSTER_PEAKS_1

CANONICAL <- c(paste0("chr", 1:19), "chrX", "chrY")

# ----- Load -----
cat("Loading PQS BED...\n")
pqs <- import(PQS_BED, format = "BED")
pqs <- pqs[as.character(seqnames(pqs)) %in% CANONICAL]
cat(sprintf("  PQS sites: %s\n", format(length(pqs), big.mark = ",")))

cat("Loading cluster 0 (MEF) peaks...\n")
c0 <- import(PEAK_CLUST0, format = "narrowPeak")
c0 <- c0[as.character(seqnames(c0)) %in% CANONICAL]
cat(sprintf("  Cluster 0 peaks: %s\n", format(length(c0), big.mark = ",")))

cat("Loading cluster 1 (mESC) peaks...\n")
c1 <- import(PEAK_CLUST1, format = "narrowPeak")
c1 <- c1[as.character(seqnames(c1)) %in% CANONICAL]
cat(sprintf("  Cluster 1 peaks: %s\n", format(length(c1), big.mark = ",")))

# =============================================================
# Compute symmetric counts for each pair (peaks_only / pqs_only / overlap)
# Overlap direction: peaks that overlap >=1 PQS site.
# =============================================================
# mm10 canonical-chromosome effective genome size (chr1-19, X, Y)
MM10_GENOME_BP <- 2.65e9

overlap_counts <- function(peaks_gr, pqs_gr, peaks_label,
                            genome_bp = MM10_GENOME_BP) {
  n_peaks  <- length(peaks_gr)
  n_pqs    <- length(pqs_gr)
  n_int    <- sum(IRanges::overlapsAny(peaks_gr, pqs_gr,
                                        ignore.strand = TRUE))
  n_pqs_in <- sum(IRanges::overlapsAny(pqs_gr, peaks_gr,
                                        ignore.strand = TRUE))

  # ----- Significance: binomial test of overlap vs random expectation -----
  # Random expectation: prob that a randomly placed peak overlaps a PQS site
  # is the fraction of the genome covered by PQS (after collapsing overlaps).
  pqs_reduced  <- GenomicRanges::reduce(pqs_gr, ignore.strand = TRUE)
  pqs_bp_total <- sum(as.numeric(width(pqs_reduced)))
  pqs_frac     <- pqs_bp_total / genome_bp

  # Account for finite peak width: a peak of width w overlaps a PQS if the
  # peak's center lies within (PQS_start - w/2, PQS_end + w/2). Approximate
  # by inflating PQS coverage by the median peak width.
  med_peak_w   <- median(width(peaks_gr))
  pqs_frac_adj <- min(1, pqs_frac * (1 + med_peak_w / 30))  # PQS ~30 bp avg
  exp_overlap  <- n_peaks * pqs_frac_adj
  fold_enrich  <- n_int / exp_overlap
  bin_p        <- binom.test(n_int, n_peaks,
                              p = pqs_frac_adj,
                              alternative = "greater")$p.value

  list(
    peaks_only        = n_peaks - n_int,
    pqs_only          = n_pqs   - n_int,
    overlap           = n_int,
    n_peaks           = n_peaks,
    n_pqs             = n_pqs,
    pct_peaks_overlap = 100 * n_int    / n_peaks,
    pct_pqs_overlap   = 100 * n_pqs_in / n_pqs,
    expected_overlap  = exp_overlap,
    fold_enrichment   = fold_enrich,
    binom_pvalue      = bin_p,
    label             = peaks_label
  )
}

cat("\nComputing overlap counts...\n")
oc0 <- overlap_counts(c0, pqs, "scG4 MEF (cluster 0)")
oc1 <- overlap_counts(c1, pqs, "scG4 mESC (cluster 1)")

cat("\nLoading mouse brain GFP+ sorted clusters...\n")
MB_GFP_CLUST0 <- file.path(GSE291468, "GSM8836086_GFPpos_cluster_0_peaks.narrowPeak")
MB_GFP_CLUST1 <- file.path(GSE291468, "GSM8836086_GFPpos_cluster_1_peaks.narrowPeak")
MB_GFP_CLUST2 <- file.path(GSE291468, "GSM8836086_GFPpos_cluster_2_peaks.narrowPeak")
stopifnot(all(file.exists(c(MB_GFP_CLUST0, MB_GFP_CLUST1, MB_GFP_CLUST2))))
mb_gfp_c0 <- import(MB_GFP_CLUST0, format = "narrowPeak")
mb_gfp_c0 <- mb_gfp_c0[as.character(seqnames(mb_gfp_c0)) %in% CANONICAL]
cat(sprintf("  GFP+ cluster 0 peaks: %s\n", format(length(mb_gfp_c0), big.mark = ",")))
mb_gfp_c1 <- import(MB_GFP_CLUST1, format = "narrowPeak")
mb_gfp_c1 <- mb_gfp_c1[as.character(seqnames(mb_gfp_c1)) %in% CANONICAL]
cat(sprintf("  GFP+ cluster 1 peaks: %s\n", format(length(mb_gfp_c1), big.mark = ",")))
mb_gfp_c2 <- import(MB_GFP_CLUST2, format = "narrowPeak")
mb_gfp_c2 <- mb_gfp_c2[as.character(seqnames(mb_gfp_c2)) %in% CANONICAL]
cat(sprintf("  GFP+ cluster 2 peaks: %s\n", format(length(mb_gfp_c2), big.mark = ",")))

cat("\nLoading mouse brain unsorted clusters...\n")
MB_UNS_CLUST0 <- file.path(GSE291468, "GSM8836087_unsorted_cluster_0_peaks.narrowPeak")
MB_UNS_CLUST1 <- file.path(GSE291468, "GSM8836087_unsorted_cluster_1_peaks.narrowPeak")
stopifnot(all(file.exists(c(MB_UNS_CLUST0, MB_UNS_CLUST1))))
mb_uns_c0 <- import(MB_UNS_CLUST0, format = "narrowPeak")
mb_uns_c0 <- mb_uns_c0[as.character(seqnames(mb_uns_c0)) %in% CANONICAL]
cat(sprintf("  Unsorted cluster 0 peaks: %s\n", format(length(mb_uns_c0), big.mark = ",")))
mb_uns_c1 <- import(MB_UNS_CLUST1, format = "narrowPeak")
mb_uns_c1 <- mb_uns_c1[as.character(seqnames(mb_uns_c1)) %in% CANONICAL]
cat(sprintf("  Unsorted cluster 1 peaks: %s\n", format(length(mb_uns_c1), big.mark = ",")))

cat("\nComputing mouse brain overlap counts...\n")
oc_mb_gfp_c0 <- overlap_counts(mb_gfp_c0, pqs, "Mouse brain GFP+ cluster 0")
oc_mb_gfp_c1 <- overlap_counts(mb_gfp_c1, pqs, "Mouse brain GFP+ cluster 1")
oc_mb_gfp_c2 <- overlap_counts(mb_gfp_c2, pqs, "Mouse brain GFP+ cluster 2")
oc_mb_uns_c0 <- overlap_counts(mb_uns_c0, pqs, "Mouse brain unsorted cluster 0")
oc_mb_uns_c1 <- overlap_counts(mb_uns_c1, pqs, "Mouse brain unsorted cluster 1")

# ----- scATAC brain (width-adjusted to 270 bp) -----
cat("\nLoading mouse brain scATAC-seq peaks (width-adjusted)...\n")
SCATAC_PEAKS <- file.path(GSE291468, "GSM8836087_unsorted_mouse_brain/CellRanger/peaks.bed")
stopifnot(file.exists(SCATAC_PEAKS))
scatac_raw <- import(SCATAC_PEAKS, format = "BED")
scatac_raw <- scatac_raw[as.character(seqnames(scatac_raw)) %in% CANONICAL]
cat(sprintf("  scATAC raw peaks: %s (median width %d bp)\n",
            format(length(scatac_raw), big.mark = ","),
            median(width(scatac_raw))))

# Resize to 270 bp (matching GFP+ sorted brain scG4 median width)
TARGET_WIDTH <- 270L
scatac_resized <- resize(scatac_raw, width = TARGET_WIDTH, fix = "center")
cat(sprintf("  scATAC resized to %d bp: %s peaks\n",
            TARGET_WIDTH, format(length(scatac_resized), big.mark = ",")))

oc_scatac <- overlap_counts(scatac_resized, pqs, "scATAC brain (width-adjusted)")
cat(sprintf("  scATAC overlap: %.1f%% (%s / %s)  %.2fx  p=%g\n",
            oc_scatac$pct_peaks_overlap,
            format(oc_scatac$overlap, big.mark = ","),
            format(oc_scatac$n_peaks, big.mark = ","),
            oc_scatac$fold_enrichment,
            oc_scatac$binom_pvalue))

print_oc <- function(oc) {
  cat(sprintf("  %s\n", oc$label))
  cat(sprintf("    peaks only:        %s\n",  format(oc$peaks_only, big.mark = ",")))
  cat(sprintf("    PQS only:          %s\n",  format(oc$pqs_only,   big.mark = ",")))
  cat(sprintf("    overlap:           %s\n",  format(oc$overlap,    big.mark = ",")))
  cat(sprintf("    %% peaks ovlp PQS:  %.1f%%   (%s / %s)\n",
              oc$pct_peaks_overlap,
              format(oc$overlap, big.mark = ","),
              format(oc$n_peaks, big.mark = ",")))
  cat(sprintf("    expected overlap:  %s   (random baseline from PQS coverage)\n",
              format(round(oc$expected_overlap), big.mark = ",")))
  cat(sprintf("    fold enrichment:   %.2fx\n", oc$fold_enrichment))
  cat(sprintf("    binomial p-value:  %g\n",   oc$binom_pvalue))
}
cat("\n=== MEF/ESC scG4 ===\n")
print_oc(oc0); print_oc(oc1)
cat("\n=== Mouse brain GFP+ ===\n")
print_oc(oc_mb_gfp_c0); print_oc(oc_mb_gfp_c1); print_oc(oc_mb_gfp_c2)
cat("\n=== Mouse brain unsorted ===\n")
print_oc(oc_mb_uns_c0); print_oc(oc_mb_uns_c1)
cat("\n=== scATAC brain (width-adjusted) ===\n")
print_oc(oc_scatac)

# ----- regioneR analysis: per-dataset (not per-cluster) -----
cat("\n=== regioneR analysis (per dataset, chr1, 1000 permutations) ===\n")
suppressPackageStartupMessages(library(regioneR))
CHROM <- "chr1"; N_TIMES <- 1000; N_CORES <- 4; SEED <- 42
CHR1_LEN <- 195154279; GENOME_CHR1 <- GRanges(CHROM, IRanges(1, CHR1_LEN))

restrict_chr1 <- function(gr) {
  gr <- gr[as.character(seqnames(gr)) == CHROM]
  GenomeInfoDb::seqlevels(gr, pruning.mode = "coarse") <- CHROM
  GenomeInfoDb::seqlengths(gr) <- CHR1_LEN
  gr
}

# Combine peaks by dataset category
c_mesc_combined <- c(c0, c1)
mb_gfp_combined <- c(mb_gfp_c0, mb_gfp_c1, mb_gfp_c2)
mb_uns_combined <- c(mb_uns_c0, mb_uns_c1)

cat(sprintf("  MEF/ESC combined: %s peaks\n", format(length(c_mesc_combined), big.mark = ",")))
cat(sprintf("  GFP+ combined:    %s peaks\n", format(length(mb_gfp_combined), big.mark = ",")))
cat(sprintf("  Unsorted combined: %s peaks\n", format(length(mb_uns_combined), big.mark = ",")))
cat(sprintf("  scATAC (270bp):   %s peaks\n", format(length(scatac_resized), big.mark = ",")))

pqs_chr <- restrict_chr1(pqs)
c_mesc_chr <- restrict_chr1(c_mesc_combined)
mb_gfp_chr <- restrict_chr1(mb_gfp_combined)
mb_uns_chr <- restrict_chr1(mb_uns_combined)
scatac_chr <- restrict_chr1(scatac_resized)

run_perm <- function(peaks_chr, label) {
  set.seed(SEED)
  res <- regioneR::overlapPermTest(A = peaks_chr, B = pqs_chr, ntimes = N_TIMES,
           genome = GENOME_CHR1, alternative = "greater", mc.cores = N_CORES, verbose = FALSE)
  no <- res$numOverlaps
  cat(sprintf("  [%s] observed=%s  mean=%.1f  fold=%.2fx  z=%.2f  p < %g\n",
      label, format(no$observed, big.mark = ","), mean(no$permuted),
      no$observed / mean(no$permuted), no$zscore, 1 / N_TIMES))
  list(
    label = label,
    observed = no$observed,
    mean_expected = mean(no$permuted),
    fold = no$observed / mean(no$permuted),
    zscore = no$zscore
  )
}

cat("\nRunning permutation tests...\n")
res_mesc <- run_perm(c_mesc_chr, "MEF/ESC combined")
res_gfp  <- run_perm(mb_gfp_chr, "GFP+ combined")
res_uns  <- run_perm(mb_uns_chr, "Unsorted combined")
res_scatac <- run_perm(scatac_chr, "scATAC (270bp)")

# Store results for bar plot subtitles
regioneR_stats <- list(
  mesc_mef = res_mesc,
  gfp = res_gfp,
  uns = res_uns,
  scatac = res_scatac
)

# =============================================================
# Build one pairwise (2-set) Venn per cluster.
# Uses VennDiagram::draw.pairwise.venn(scaled = FALSE) so the
# two circles are equal sized regardless of set-size disparity
# (PQS is ~300x larger than scG4, so any proportional rendering
# would hide the overlap).
# =============================================================

make_venn_grob <- function(oc, fill_cols, title) {
  # Use bigger set (PQS) as area1 so it lands on the LEFT (peach).
  # VennDiagram auto-flips if area2 > area1; passing in the natural
  # bigger-first order avoids that flip and keeps colors aligned.
  g_list <- VennDiagram::draw.pairwise.venn(
              area1      = oc$n_pqs,                 # bigger set first
              area2      = oc$n_peaks,
              cross.area = oc$overlap,
              category   = c("", ""),                # no labels above circles
              scaled     = FALSE,                    # equal-sized circles
              fill       = c(fill_cols[2], fill_cols[1]),  # PQS=peach, scG4=blue
              alpha      = 0.55,
              lwd        = 1.5,
              col        = "grey30",
              cex        = 1.3,                      # numbers inside circles
              fontface   = "bold",
              ind        = FALSE                     # return gList
            )
  # draw.pairwise.venn returns a gList; collapse to a single grob
  g <- grid::grobTree(g_list)

  # Two-line title: main title + significance subtitle
  subtitle <- sprintf("overlap = %s   |   expected = %s   |   %.2fx enrichment   |   %s",
                      format(oc$overlap,                    big.mark = ","),
                      format(round(oc$expected_overlap),    big.mark = ","),
                      oc$fold_enrichment,
                      format_pvalue(oc$binom_pvalue))

  title_grob <- gridExtra::arrangeGrob(
    grid::textGrob(title,
                    gp = grid::gpar(fontsize = 13, fontface = "bold")),
    grid::textGrob(subtitle,
                    gp = grid::gpar(fontsize = 10, fontface = "plain",
                                     col = "grey20")),
    ncol = 1, heights = grid::unit(c(1.2, 1), "lines")
  )

  gridExtra::arrangeGrob(g, top = title_grob)
}

cat("\nGenerating summary bar plots...\n")

# =============================================================
# Bar chart 1: MEF/ESC scG4 clusters
# =============================================================
bar_mesc_mef <- data.frame(
  cluster = c("MEF (cluster 0)", "mESC (cluster 1)"),
  pct     = c(oc0$pct_peaks_overlap, oc1$pct_peaks_overlap),
  overlap = c(oc0$overlap, oc1$overlap),
  n_peaks = c(oc0$n_peaks, oc1$n_peaks),
  fold    = c(oc0$fold_enrichment, oc1$fold_enrichment),
  pval    = c(oc0$binom_pvalue, oc1$binom_pvalue)
)
bar_mesc_mef$label <- sprintf("%.1f%%\n(%s / %s)",
                              bar_mesc_mef$pct,
                              format(bar_mesc_mef$overlap, big.mark = ","),
                              format(bar_mesc_mef$n_peaks, big.mark = ","))

p_mesc_mef <- ggplot(bar_mesc_mef, aes(x = cluster, y = pct)) +
  geom_col(width = 0.6, fill = "#9ecae1", color = "grey30") +
  geom_text(aes(label = label), vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_y_continuous(limits = c(0, 105),
                     breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL,
       y = "% of scG4 peaks overlapping >=1 PQS site",
       title = "MEF/ESC scG4 clusters overlap with PQS",
       subtitle = sprintf("%.2fx enrichment, z=%.1f, p<0.001 (regioneR chr1)",
                          regioneR_stats$mesc_mef$fold, regioneR_stats$mesc_mef$zscore)) +
  theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file.path(OUT_DIR, "bar_pct_peaks_overlap_PQS_MEF_ESC_scG4.pdf"), p_mesc_mef,
       width = 6, height = 6)
cat("  Saved bar_pct_peaks_overlap_PQS_MEF_ESC_scG4.pdf\n")

# =============================================================
# Bar chart 2: Mouse brain GFP+ clusters
# =============================================================
bar_mb_gfp <- data.frame(
  cluster = c("GFP+ (cluster 0)", "GFP+ (cluster 1)", "GFP+ (cluster 2)"),
  pct     = c(oc_mb_gfp_c0$pct_peaks_overlap, oc_mb_gfp_c1$pct_peaks_overlap, oc_mb_gfp_c2$pct_peaks_overlap),
  overlap = c(oc_mb_gfp_c0$overlap, oc_mb_gfp_c1$overlap, oc_mb_gfp_c2$overlap),
  n_peaks = c(oc_mb_gfp_c0$n_peaks, oc_mb_gfp_c1$n_peaks, oc_mb_gfp_c2$n_peaks),
  fold    = c(oc_mb_gfp_c0$fold_enrichment, oc_mb_gfp_c1$fold_enrichment, oc_mb_gfp_c2$fold_enrichment),
  pval    = c(oc_mb_gfp_c0$binom_pvalue, oc_mb_gfp_c1$binom_pvalue, oc_mb_gfp_c2$binom_pvalue)
)
bar_mb_gfp$label <- sprintf("%.1f%%\n(%s / %s)",
                            bar_mb_gfp$pct,
                            format(bar_mb_gfp$overlap, big.mark = ","),
                            format(bar_mb_gfp$n_peaks, big.mark = ","))

p_mb_gfp <- ggplot(bar_mb_gfp, aes(x = cluster, y = pct)) +
  geom_col(width = 0.6, fill = "#fc9272", color = "grey30") +
  geom_text(aes(label = label), vjust = -0.3, size = 3.3, fontface = "bold") +
  scale_y_continuous(limits = c(0, 105),
                     breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL,
       y = "% of scG4 peaks overlapping >=1 PQS site",
       title = "Mouse brain GFP+ clusters overlap with PQS",
       subtitle = sprintf("%.2fx enrichment, z=%.1f, p<0.001 (regioneR chr1)",
                          regioneR_stats$gfp$fold, regioneR_stats$gfp$zscore)) +
  theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file.path(OUT_DIR, "bar_pct_peaks_overlap_PQS_mousebrain_GFPpos.pdf"), p_mb_gfp,
       width = 7, height = 6)
cat("  Saved bar_pct_peaks_overlap_PQS_mousebrain_GFPpos.pdf\n")

# =============================================================
# Bar chart 3: Mouse brain unsorted clusters
# =============================================================
bar_mb_uns <- data.frame(
  cluster = c("Unsorted (cluster 0)", "Unsorted (cluster 1)"),
  pct     = c(oc_mb_uns_c0$pct_peaks_overlap, oc_mb_uns_c1$pct_peaks_overlap),
  overlap = c(oc_mb_uns_c0$overlap, oc_mb_uns_c1$overlap),
  n_peaks = c(oc_mb_uns_c0$n_peaks, oc_mb_uns_c1$n_peaks),
  fold    = c(oc_mb_uns_c0$fold_enrichment, oc_mb_uns_c1$fold_enrichment),
  pval    = c(oc_mb_uns_c0$binom_pvalue, oc_mb_uns_c1$binom_pvalue)
)
bar_mb_uns$label <- sprintf("%.1f%%\n(%s / %s)",
                            bar_mb_uns$pct,
                            format(bar_mb_uns$overlap, big.mark = ","),
                            format(bar_mb_uns$n_peaks, big.mark = ","))

p_mb_uns <- ggplot(bar_mb_uns, aes(x = cluster, y = pct)) +
  geom_col(width = 0.6, fill = "#a1d99b", color = "grey30") +
  geom_text(aes(label = label), vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_y_continuous(limits = c(0, 105),
                     breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL,
       y = "% of scG4 peaks overlapping >=1 PQS site",
       title = "Mouse brain unsorted clusters overlap with PQS",
       subtitle = sprintf("%.2fx enrichment, z=%.1f, p<0.001 (regioneR chr1)",
                          regioneR_stats$uns$fold, regioneR_stats$uns$zscore)) +
  theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file.path(OUT_DIR, "bar_pct_peaks_overlap_PQS_mousebrain_unsorted.pdf"), p_mb_uns,
       width = 6, height = 6)
cat("  Saved bar_pct_peaks_overlap_PQS_mousebrain_unsorted.pdf\n")

# =============================================================
# Bar chart 4: scATAC brain (width-adjusted to 270 bp)
# =============================================================
bar_scatac <- data.frame(
  assay = "scATAC brain",
  pct     = oc_scatac$pct_peaks_overlap,
  overlap = oc_scatac$overlap,
  n_peaks = oc_scatac$n_peaks,
  fold    = oc_scatac$fold_enrichment,
  pval    = oc_scatac$binom_pvalue
)
bar_scatac$label <- sprintf("%.1f%%\n(%s / %s)",
                            bar_scatac$pct,
                            format(bar_scatac$overlap, big.mark = ","),
                            format(bar_scatac$n_peaks, big.mark = ","))

p_scatac <- ggplot(bar_scatac, aes(x = assay, y = pct)) +
  geom_col(width = 0.3, fill = "#b3de69", color = "grey30") +
  geom_text(aes(label = label), vjust = -0.3, size = 4, fontface = "bold") +
  scale_y_continuous(limits = c(0, 105),
                     breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL,
       y = "% of scATAC peaks overlapping >=1 PQS site",
       title = "Mouse brain scATAC-seq overlap with PQS (width-adjusted to 270 bp)",
       subtitle = sprintf("%.2fx enrichment, z=%.1f, p<0.001 (regioneR chr1)",
                          regioneR_stats$scatac$fold, regioneR_stats$scatac$zscore)) +
  theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 11))

ggsave(file.path(OUT_DIR, "bar_pct_peaks_overlap_PQS_scATAC.pdf"), p_scatac,
       width = 5, height = 6)
cat("  Saved bar_pct_peaks_overlap_PQS_scATAC.pdf\n")

# ----- Save overlap table (all sc datasets) -----
overlap_tbl <- data.frame(
  comparison           = c("Cluster 0 (MEF) vs PQS",
                           "Cluster 1 (mESC) vs PQS",
                           "Mouse brain GFP+ cluster 0 vs PQS",
                           "Mouse brain GFP+ cluster 1 vs PQS",
                           "Mouse brain GFP+ cluster 2 vs PQS",
                           "Mouse brain unsorted cluster 0 vs PQS",
                           "Mouse brain unsorted cluster 1 vs PQS",
                           "scATAC brain (width-adjusted to 270 bp) vs PQS"),
  n_peaks              = c(oc0$n_peaks, oc1$n_peaks,
                           oc_mb_gfp_c0$n_peaks, oc_mb_gfp_c1$n_peaks, oc_mb_gfp_c2$n_peaks,
                           oc_mb_uns_c0$n_peaks, oc_mb_uns_c1$n_peaks,
                           oc_scatac$n_peaks),
  n_PQS_sites          = rep(oc0$n_pqs, 8),
  n_overlap_PQS        = c(oc0$overlap, oc1$overlap,
                           oc_mb_gfp_c0$overlap, oc_mb_gfp_c1$overlap, oc_mb_gfp_c2$overlap,
                           oc_mb_uns_c0$overlap, oc_mb_uns_c1$overlap,
                           oc_scatac$overlap),
  pct_peaks_overlap_PQS = round(c(oc0$pct_peaks_overlap, oc1$pct_peaks_overlap,
                                   oc_mb_gfp_c0$pct_peaks_overlap, oc_mb_gfp_c1$pct_peaks_overlap, oc_mb_gfp_c2$pct_peaks_overlap,
                                   oc_mb_uns_c0$pct_peaks_overlap, oc_mb_uns_c1$pct_peaks_overlap,
                                   oc_scatac$pct_peaks_overlap), 1),
  pct_PQS_overlap_peaks = round(c(oc0$pct_pqs_overlap, oc1$pct_pqs_overlap,
                                   oc_mb_gfp_c0$pct_pqs_overlap, oc_mb_gfp_c1$pct_pqs_overlap, oc_mb_gfp_c2$pct_pqs_overlap,
                                   oc_mb_uns_c0$pct_pqs_overlap, oc_mb_uns_c1$pct_pqs_overlap,
                                   oc_scatac$pct_pqs_overlap), 2),
  expected_overlap     = round(c(oc0$expected_overlap, oc1$expected_overlap,
                                  oc_mb_gfp_c0$expected_overlap, oc_mb_gfp_c1$expected_overlap, oc_mb_gfp_c2$expected_overlap,
                                  oc_mb_uns_c0$expected_overlap, oc_mb_uns_c1$expected_overlap,
                                  oc_scatac$expected_overlap)),
  fold_enrichment      = round(c(oc0$fold_enrichment,  oc1$fold_enrichment,
                                  oc_mb_gfp_c0$fold_enrichment,  oc_mb_gfp_c1$fold_enrichment,  oc_mb_gfp_c2$fold_enrichment,
                                  oc_mb_uns_c0$fold_enrichment,  oc_mb_uns_c1$fold_enrichment,
                                  oc_scatac$fold_enrichment), 2),
  binomial_pvalue      = c(oc0$binom_pvalue, oc1$binom_pvalue,
                           oc_mb_gfp_c0$binom_pvalue, oc_mb_gfp_c1$binom_pvalue, oc_mb_gfp_c2$binom_pvalue,
                           oc_mb_uns_c0$binom_pvalue, oc_mb_uns_c1$binom_pvalue,
                           oc_scatac$binom_pvalue)
)
write.csv(overlap_tbl,
          file.path(OUT_DIR, "overlap_counts.csv"),
          row.names = FALSE)

cat("\nDone. Files in:", OUT_DIR, "\n")
cat("  Bar plots (PQS overlap):\n")
cat("    - bar_pct_peaks_overlap_PQS_MEF_ESC_scG4.pdf\n")
cat("    - bar_pct_peaks_overlap_PQS_mousebrain_GFPpos.pdf\n")
cat("    - bar_pct_peaks_overlap_PQS_mousebrain_unsorted.pdf\n")
cat("    - bar_pct_peaks_overlap_PQS_scATAC.pdf  (width-adjusted to 270 bp)\n")
cat("  Data:\n")
cat("    - overlap_counts.csv  (all sc datasets)\n")
cat("    - permTest_summary.csv  (permutation test on chr1)\n")

# Save regioneR summary
permTest_summary <- do.call(rbind, lapply(regioneR_stats, function(rs) {
  data.frame(comparison = rs$label, chromosome = CHROM, ntimes = N_TIMES,
             observed_overlap = rs$observed,
             mean_expected = round(rs$mean_expected, 1),
             fold_enrichment = round(rs$fold, 2),
             zscore = round(rs$zscore, 2),
             pvalue_reported = sprintf("p < %g", 1 / N_TIMES))
}))
write.csv(permTest_summary, file.path(OUT_DIR, "permTest_summary.csv"), row.names = FALSE)
cat("\nWrote permTest_summary.csv\n"); print(permTest_summary)
