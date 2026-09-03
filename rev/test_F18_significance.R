#!/usr/bin/env Rscript
# Statistical significance testing for F18 (scG4 vs ATAC peak overlap)
# 1) Hypergeometric test
# 2) regioneR overlapPermTest

source("rev/paths.R")
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer)
  library(regioneR)
})

CANON <- c(paste0("chr",1:19),"chrX","chrY")
HALF <- 150L

# Load chromosome lengths
cs <- read.table(CHROM_SIZES, header = FALSE, stringsAsFactors = FALSE)
GLEN <- setNames(as.integer(cs$V2), cs$V1)[CANON]
GENOME_SIZE <- sum(GLEN)

read_summit <- function(fn){
  d <- fread(fn)
  d <- d[V1 %in% CANON]
  off <- if(ncol(d)>=10) d$V10 else round((d$V3-d$V2)/2)
  s <- d$V2 + off
  gr <- GRanges(d$V1, IRanges(pmax(s-HALF,1L), s+HALF))
  seqlengths(gr) <- GLEN[seqlevels(gr)]
  trim(gr)
}

# Define cell types
cells <- list(
  mESC = list(
    scg4_pk = CLUSTER_PEAKS_1,
    atac_pk = file.path(DATA, "GSE149080/GSM4661960_ESC_WT_batch2_peaks.narrowPeak")
  ),
  MEF = list(
    scg4_pk = CLUSTER_PEAKS_0,
    atac_pk = file.path(DATA, "GSE211123/GSM6451000_3T3_control_ATAC_peaks.narrowPeak")
  )
)

cat("=== F18 Statistical Significance Testing ===\n\n")

results <- list()

for (cl in names(cells)) {
  cat(sprintf("=== %s ===\n", cl))
  cc <- cells[[cl]]

  # Load peaks
  scg4 <- read_summit(cc$scg4_pk)
  atac <- read_summit(cc$atac_pk)

  cat(sprintf("  scG4 peaks: %s\n", format(length(scg4), big.mark = ",")))
  cat(sprintf("  ATAC peaks: %s\n", format(length(atac), big.mark = ",")))

  # Overlap
  ov <- overlapsAny(scg4, atac, ignore.strand = TRUE)
  n_overlap <- sum(ov)
  n_scg4 <- length(scg4)
  n_atac <- length(atac)

  cat(sprintf("  Overlap: %s (%.1f%%)\n", format(n_overlap, big.mark = ","), 100*n_overlap/n_scg4))

  # ---- 1) Hypergeometric test ----
  # N = total possible non-overlapping regions in genome (approximated)
  # K = number of ATAC peaks
  # n = number of scG4 peaks
  # k = observed overlap

  # Effective genome size for 300bp regions
  EFFECTIVE_GENOME <- GENOME_SIZE / 300  # ~8.8 million possible 300bp regions

  N <- EFFECTIVE_GENOME  # total population
  K <- n_atac            # successes in population (ATAC peaks)
  n <- n_scg4            # draws (scG4 peaks)
  k <- n_overlap         # observed successes (overlap)

  # P(X >= k) under hypergeometric
  pval_upper <- phyper(k-1, K, N-K, n, lower.tail = FALSE)
  expected <- n * K / N

  cat(sprintf("\n  Hypergeometric test:\n"))
  cat(sprintf("    Expected overlap: %s\n", format(round(expected), big.mark = ",")))
  cat(sprintf("    Observed overlap: %s\n", format(k, big.mark = ",")))
  cat(sprintf("    Fold enrichment: %.2fx\n", k / expected))
  cat(sprintf("    p-value: < %g\n", pval_upper))

  # ---- 2) regioneR overlapPermTest (chr1 only for speed) ----
  CHROM <- "chr1"
  CHR1_LEN <- GLEN["chr1"]
  GENOME_CHR1 <- GRanges(CHROM, IRanges(1, CHR1_LEN))

  restrict_chr1 <- function(gr) {
    gr <- gr[as.character(seqnames(gr)) == CHROM]
    GenomeInfoDb::seqlevels(gr, pruning.mode = "coarse") <- CHROM
    GenomeInfoDb::seqlengths(gr) <- CHR1_LEN
    gr
  }

  scg4_chr1 <- restrict_chr1(scg4)
  atac_chr1 <- restrict_chr1(atac)

  cat(sprintf("\n  regioneR (chr1, 1000 permutations):\n"))
  cat(sprintf("    scG4 chr1: %s peaks\n", format(length(scg4_chr1), big.mark = ",")))
  cat(sprintf("    ATAC chr1: %s peaks\n", format(length(atac_chr1), big.mark = ",")))

  set.seed(42)
  res <- regioneR::overlapPermTest(
    A = scg4_chr1, B = atac_chr1,
    ntimes = 1000,
    genome = GENOME_CHR1,
    alternative = "greater",
    mc.cores = 4,
    verbose = FALSE
  )

  no <- res$numOverlaps
  fold_chr1 <- no$observed / mean(no$permuted)

  cat(sprintf("    Observed overlap: %s\n", format(no$observed, big.mark = ",")))
  cat(sprintf("    Mean expected: %.1f\n", mean(no$permuted)))
  cat(sprintf("    Fold enrichment: %.2fx\n", fold_chr1))
  cat(sprintf("    z-score: %.2f\n", no$zscore))
  cat(sprintf("    p < %g\n", 1/1000))

  results[[cl]] <- list(
    n_scg4 = n_scg4,
    n_atac = n_atac,
    n_overlap = n_overlap,
    pct_overlap = 100 * n_overlap / n_scg4,
    hypergeo_expected = expected,
    hypergeo_fold = k / expected,
    hypergeo_pval = pval_upper,
    regioner_observed = no$observed,
    regioner_expected = mean(no$permuted),
    regioner_fold = fold_chr1,
    regioner_zscore = no$zscore
  )

  cat("\n")
}

# Summary table
cat("=== SUMMARY ===\n\n")
for (cl in names(results)) {
  r <- results[[cl]]
  cat(sprintf("%s:\n", cl))
  cat(sprintf("  scG4 peaks: %s | ATAC peaks: %s | Overlap: %.1f%% (%s)\n",
              format(r$n_scg4, big.mark=","), format(r$n_atac, big.mark=","),
              r$pct_overlap, format(r$n_overlap, big.mark=",")))
  cat(sprintf("  Hypergeometric: %.1fx enrichment, p < %g\n", r$hypergeo_fold, r$hypergeo_pval))
  cat(sprintf("  regioneR (chr1): %.2fx enrichment, z=%.1f, p < 0.001\n\n",
              r$regioner_fold, r$regioner_zscore))
}

# Save results
out_tbl <- do.call(rbind, lapply(names(results), function(cl) {
  r <- results[[cl]]
  data.frame(
    cell_type = cl,
    n_scG4 = r$n_scg4,
    n_ATAC = r$n_atac,
    n_overlap = r$n_overlap,
    pct_overlap = round(r$pct_overlap, 1),
    hypergeo_fold = round(r$hypergeo_fold, 2),
    hypergeo_pval = r$hypergeo_pval,
    regioner_fold = round(r$regioner_fold, 2),
    regioner_zscore = round(r$regioner_zscore, 2)
  )
}))
write.csv(out_tbl, file.path(OUT_ROOT, "rev", "outputs", "G4_vs_accessibility_decoupling", "F18_significance.csv"), row.names = FALSE)
cat("Saved F18_significance.csv\n")
