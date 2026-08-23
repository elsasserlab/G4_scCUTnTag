#!/usr/bin/env Rscript
# ===========================================================================
# END-TO-END (single self-contained run from raw inputs):
#   1. Merge the two unsorted-brain scG4 MACS2 cluster peak files -> summits.
#   2. Derive per-cluster brain ATAC accessible peaks DIRECTLY from the
#      GSE198467 ATAC Seurat object (peak detected in >= 5% of a cluster's cells).
#   3. regioneR overlapPermTest (genome-wide) of merged unsorted G4 AND each
#      ATAC cluster, vs PQS and eG4  -> fold enrichment + bootstrap 95% CI.
#   4. Significance: per-bar enrichment (perm z) + G4 > each ATAC cluster
#      (bootstrap one-sided p on fold diff + two-proportion z-test on rate).
#   5. Two annotated fold-enrichment plots (PQS, eG4) with G4-vs-cluster stars.
#
# Inputs : GSM8836087_scG4CnT_MouseBrain_unsorted_cluster_{0,1}_peaks.narrowPeak
#          GSE198467_ATAC_Seurat_object_clustered_renamed.Rds  (raw brain ATAC object)
#          PQS_scores_mm10.bed , Mouse_eG4.txt
# Outputs: brain_unsortedG4_MACS2_enrichment/
#            scG4CnT_MouseBrain_unsorted_MERGED_{peaks,summits}.bed
#            unsortedG4_enrichment_full_results.csv
#            foldenrichment_{PQS,eG4}_barplot.pdf
# ===========================================================================


source("rev/paths.R")
suppressMessages({
  library(regioneR); library(GenomicRanges); library(rtracklayer)
  library(ggplot2)
  library(Seurat); library(Signac); library(Matrix)   # to derive ATAC peaks from the object
})
BSGENOME_AVAILABLE <- requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)
if (BSGENOME_AVAILABLE) suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
inp <- function(fn) {
  hits <- list.files(DATA, pattern = fn, recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("File not found in data/: ", fn)
  hits[1]
}

## ---- parameters -----------------------------------------------------------
NTIMES    <- 1000   # permutations (enrichment p)
NBOOT_CI  <- 2000   # bootstrap resamples for fold CI
NBOOT_DIF <- 5000   # bootstrap resamples for G4-vs-cluster difference test
MC        <- 4
PQS_MIN   <- 50     # PQS score cutoff
AHALF     <- 75     # re-centre every tested set to +/- this (150 bp)
NP0   <- file.path(GSE291468, "GSM8836087_unsorted_cluster_0_peaks.narrowPeak")
NP1   <- file.path(GSE291468, "GSM8836087_unsorted_cluster_1_peaks.narrowPeak")
ATAC_RDS    <- file.path(GSE198467, "GSE198467_ATAC_Seurat_object_clustered_renamed.Rds")  # raw brain ATAC object
CLUSTER_COL <- "idents_short"   # 15 named brain cell types
DETECT_FRAC <- 0.05             # peak "accessible" in a cluster if detected in >= 5% of its cells
OUTDIR<- file.path(OUT_ROOT, "brain_unsortedG4_MACS2_enrichment")
set.seed(1)
dir.create(OUTDIR, showWarnings = FALSE)

STD_CHR <- paste0("chr", c(1:19, "X", "Y"))
if (BSGENOME_AVAILABLE) {
  sl <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)[STD_CHR]
} else {
  if (file.exists(CHROM_SIZES)) {
    cs <- read.table(CHROM_SIZES, header = FALSE, stringsAsFactors = FALSE)
    sl <- setNames(as.integer(cs$V2), cs$V1)[STD_CHR]
    message("Loaded chromosome lengths from ", CHROM_SIZES)
  } else {
    stop("BSgenome not available and chrom.sizes file not found: ", CHROM_SIZES)
  }
}
genome  <- toGRanges(data.frame(chr = STD_CHR, start = 1, end = as.integer(sl)))
onstd   <- function(gr) GenomeInfoDb::keepStandardChromosomes(
             gr[as.character(seqnames(gr)) %in% STD_CHR], pruning.mode = "coarse")
prep    <- function(gr) GenomicRanges::reduce(resize(onstd(gr), width = 2 * AHALF, fix = "center"))
stars   <- function(p) ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**", ifelse(p < 0.05, "*", "ns")))
narrow_cols <- c("chrom","start","end","name","score","strand","signal","pval","qval","peak")

## ---- STEP 1: merge the two cluster peak files -----------------------------
read_np <- function(f){ x <- read.table(f, sep="\t", stringsAsFactors=FALSE); colnames(x) <- narrow_cols; x }
np0 <- read_np(NP0); np1 <- read_np(NP1)
cat(sprintf("cluster_0 peaks: %d | cluster_1 peaks: %d\n", nrow(np0), nrow(np1)))

ivl <- c(GRanges(np0$chrom, IRanges(np0$start+1, np0$end)),
         GRanges(np1$chrom, IRanges(np1$start+1, np1$end)))
ivl <- sort(onstd(GenomicRanges::reduce(ivl)))
export.bed(ivl, file.path(OUTDIR, "scG4CnT_MouseBrain_unsorted_MERGED_peaks.bed"))

sm  <- c(GRanges(np0$chrom, IRanges(np0$start + np0$peak + 1, width=1)),
         GRanges(np1$chrom, IRanges(np1$start + np1$peak + 1, width=1)))
sm  <- onstd(sm)
export.bed(sm, file.path(OUTDIR, "scG4CnT_MouseBrain_unsorted_MERGED_summits.bed"))
cat(sprintf("merged: %d peak intervals, %d summits\n", length(ivl), length(sm)))

## ---- references + tested peak sets ----------------------------------------
cat("Loading references (genome-wide)...\n")
pqs <- onstd(import(PQS_BED, format="BED")); pqs <- pqs[score(pqs) >= PQS_MIN]
eg4_tab <- read.table(EG4_PATH, sep="\t", header=TRUE, quote="", comment.char="")
eg4 <- onstd(GRanges(eg4_tab[[1]], IRanges(eg4_tab[[2]]+1, eg4_tab[[3]])))
refs <- list(PQS = pqs, eG4 = eg4)
cat(sprintf("  PQS(score>=%d)=%d | eG4=%d\n", PQS_MIN, length(pqs), length(eg4)))

# ---- STEP 2: derive per-cluster brain ATAC accessible peaks from the object --
# For each cell type (idents_short), keep peaks detected in >= DETECT_FRAC of its
# cells, then centre to AHALF (same prep() as every tested set). This replaces
# the former pre-made brain_atac_cluster_beds — provenance is now in this script.
cat("Deriving per-cluster brain ATAC accessible peaks from the Seurat object...\n")
BED_HALF <- 250   # original ATAC peak half-width (500 bp), reduced, before final centring
atac_obj <- readRDS(ATAC_RDS)
labs     <- as.character(atac_obj@meta.data[[CLUSTER_COL]])
counts   <- GetAssayData(atac_obj, assay = "peaks", layer = "counts")  # peaks x cells
bin      <- counts; bin@x[] <- 1                                      # binarize detected/not
feat     <- rownames(counts)
clusters <- sort(unique(labs[!is.na(labs)]))
atac_sets <- setNames(lapply(clusters, function(cl) {
  cells <- which(labs == cl)
  keep  <- feat[Matrix::rowMeans(bin[, cells, drop = FALSE]) >= DETECT_FRAC]
  gr500 <- sort(GenomicRanges::reduce(resize(onstd(StringToGRanges(keep, sep = c("-", "-"))),
                                             width = 2 * BED_HALF, fix = "center")))
  prep(gr500)   # then re-centre to AHALF (150 bp), exactly as the former bed pipeline
}), clusters)
rm(atac_obj, counts, bin); gc()
cat(sprintf("  %d ATAC clusters derived (accessible peaks, detect>=%.0f%%)\n",
            length(atac_sets), 100 * DETECT_FRAC))
sets <- c(list(G4_unsorted_MACS2 = prep(sm)), atac_sets)

## ---- STEP 2: permTest (fold + CI) and cache overlap indicators -------------
ind <- list(); base <- list()
for (rn in names(refs)) {
  cat(sprintf("\n=== %s : permTest + bootstrap CI ===\n", rn))
  B <- refs[[rn]]
  for (sn in names(sets)) {
    A  <- sets[[sn]]
    no <- permTest(A=A, B=B, randomize.function=circularRandomizeRegions,
                   evaluate.function=numOverlaps, genome=genome, count.once=TRUE,
                   ntimes=NTIMES, mc.cores=MC, force.parallel=TRUE)$numOverlaps
    em <- mean(no$permuted); ic <- as.integer(overlapsAny(A, B)); n <- length(ic)
    bci <- quantile(vapply(seq_len(NBOOT_CI), function(i) sum(ic[sample.int(n,n,TRUE)]), numeric(1))/em,
                    c(.025,.975), names=FALSE)
    ind[[paste(rn,sn)]]  <- ic
    base[[paste(rn,sn)]] <- list(n=n, obs=no$observed, exp=em, z=no$zscore, p=no$pval,
                                 fold=no$observed/em, lo=bci[1], hi=bci[2])
    cat(sprintf("  %-20s n=%5d fold=%6.2fx [%.2f-%.2f]\n", sn, n, no$observed/em, bci[1], bci[2]))
  }
}

## ---- STEP 3: significance + assemble results ------------------------------
out <- list()
for (rn in names(refs)) {
  g <- base[[paste(rn,"G4_unsorted_MACS2")]]; iG <- ind[[paste(rn,"G4_unsorted_MACS2")]]
  for (sn in names(sets)) {
    b <- base[[paste(rn,sn)]]; ic <- ind[[paste(rn,sn)]]
    logp <- pnorm(b$z, lower.tail=FALSE, log.p=TRUE)
    p_enr <- if (logp/log(10) < -300) "<1e-300" else formatC(exp(logp), format="e", digits=2)
    if (sn == "G4_unsorted_MACS2") { p_dif <- NA; p_prop <- NA; sg <- "" } else {
      bd <- vapply(seq_len(NBOOT_DIF), function(k)
        (sum(iG[sample.int(g$n,g$n,TRUE)])/g$exp) - (sum(ic[sample.int(b$n,b$n,TRUE)])/b$exp), numeric(1))
      p_dif  <- max(mean(bd <= 0), 1/NBOOT_DIF)
      p_prop <- prop.test(c(g$obs,b$obs), c(g$n,b$n), alternative="greater")$p.value
      sg <- stars(p_dif)
    }
    out[[length(out)+1]] <- data.frame(
      reference=rn, set=sn, type=ifelse(sn=="G4_unsorted_MACS2","G4","ATAC"),
      n_peaks=b$n, observed=b$obs, exp_mean=round(b$exp,2),
      fold=round(b$fold,3), fold_lo=round(b$lo,3), fold_hi=round(b$hi,3),
      log2fold=round(log2(b$fold),3), perm_z=round(b$z,1), perm_p=signif(b$p,3),
      p_enrichment=p_enr,
      p_G4_gt_cluster_boot = if (is.na(p_dif)) NA else signif(p_dif,3),
      p_G4_gt_cluster_prop = if (is.na(p_prop)) NA else signif(p_prop,3),
      sig_vs_G4 = sg)
  }
}
res <- do.call(rbind, out)
write.csv(res, file.path(OUTDIR, "unsortedG4_enrichment_full_results.csv"), row.names=FALSE)
cat("\n=== Results ===\n"); print(res[,c("reference","set","fold","fold_lo","fold_hi","p_enrichment","p_G4_gt_cluster_boot","sig_vs_G4")], row.names=FALSE)

## ---- STEP 4: annotated fold plots (PQS, eG4 separately) -------------------
res$label_set <- ifelse(res$set=="G4_unsorted_MACS2", "G4 (unsorted)", res$set)
rank <- tapply(res$fold[res$type=="ATAC"], res$label_set[res$type=="ATAC"], mean)
res$label_set <- factor(res$label_set, levels=c("G4 (unsorted)", names(sort(rank, decreasing=TRUE))))

make_plot <- function(rn) {
  d <- res[res$reference==rn,]
  ggplot(d, aes(label_set, fold, fill=type)) +
    geom_col(color="black", width=0.78) +
    geom_errorbar(aes(ymin=fold_lo, ymax=fold_hi), width=0.3, linewidth=0.4) +
    geom_text(aes(label=sprintf("%.2fx", fold)), vjust=-0.6, size=3.0) +
    geom_text(aes(y=fold_hi, label=sig_vs_G4), vjust=-1.9, size=3.4, fontface="bold") +
    geom_hline(yintercept=1, linetype="dashed") +
    scale_fill_manual(values=c(G4="#B2182B", ATAC="#2166AC")) +
    scale_y_continuous(expand=expansion(mult=c(0,0.14))) +
    labs(title=sprintf("Fold enrichment over random - %s", rn),
         subtitle=paste0("Unsorted brain G4 (merged MACS2 summits) vs per-cluster ATAC\n",
                         "genome-wide | 1000 perms, 95% bootstrap CI | stars: G4 > cluster (bootstrap, ***p<0.001)"),
         x=NULL, y="Fold enrichment over random") +
    theme_bw(base_size=13) +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=11),
          legend.position="top", legend.title=element_blank(),
          plot.title=element_text(face="bold"),
          plot.subtitle=element_text(size=9.5, lineheight=1.1),
          plot.margin=margin(10,14,10,10))
}
fig_map <- c(PQS="F14a_brain_unsortedG4_PQS_foldenrichment", eG4="F14b_brain_unsortedG4_eG4_foldenrichment")
for (rn in c("PQS","eG4")) {
  p <- make_plot(rn)
  ggsave(file.path(OUTDIR, sprintf("foldenrichment_%s_barplot.pdf", rn)), p, width=9, height=6.2)
  cat("published:", fig_map[[rn]], "\n")
}
cat("\nDone (end-to-end).\n")
