#!/usr/bin/env Rscript
# ===========================================================================
# scG4 pseudobulk clusters vs matched bulk ATAC, enrichment for PQS and eG4.
#   3T3  group: scG4 cluster 0  vs  bulk 3T3 ATAC
#   mESC group: scG4 cluster 1  vs  bulk mESC ATAC
# Same method as the brain analysis: MACS2 summits, re-centred +/-75bp,
# regioneR overlapPermTest genome-wide, fold + bootstrap 95% CI, and
# significance (per-bar enrichment p; G4 > paired ATAC bootstrap + prop test).
# Two separate annotated plots (PQS, eG4).
# ===========================================================================


source("rev/paths.R")
suppressMessages({
  library(regioneR); library(GenomicRanges); library(rtracklayer)
  library(ggplot2)
})
BSGENOME_AVAILABLE <- requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)
if (BSGENOME_AVAILABLE) suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
inp <- function(fn) {
  hits <- list.files(DATA, pattern = fn, recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("File not found in data/: ", fn)
  hits[1]
}

NTIMES <- 1000; NBOOT_CI <- 2000; NBOOT_DIF <- 5000; MC <- 4
PQS_MIN <- 50; AHALF <- 75
OUTDIR <- file.path(OUT_ROOT, "scG4cluster_vs_bulkATAC_PQS_eG4")
set.seed(1); dir.create(OUTDIR, showWarnings = FALSE)

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
genome <- toGRanges(data.frame(chr = STD_CHR, start = 1, end = as.integer(sl)))
onstd  <- function(gr) GenomeInfoDb::keepStandardChromosomes(
            gr[as.character(seqnames(gr)) %in% STD_CHR], pruning.mode = "coarse")
narrow_cols <- c("chrom","start","end","name","score","strand","signal","pval","qval","peak")

# narrowPeak -> summit GRanges (summit = chromStart + col10), re-centred +/-AHALF
np_summits <- function(f) {
  x <- read.table(f, sep = "\t", stringsAsFactors = FALSE); colnames(x) <- narrow_cols
  gr <- GRanges(x$chrom, IRanges(x$start + x$peak + 1, width = 1))
  GenomicRanges::reduce(resize(onstd(gr), width = 2 * AHALF, fix = "center"))
}
stars <- function(p) ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**", ifelse(p < 0.05, "*", "ns")))

## ---- peak sets (with cell-type group + pairing) ---------------------------
meta <- list(
  list(key="scG4_cl0",  file=CLUSTER_PEAKS_0,          type="G4",   group="3T3",  label="scG4 cl0"),
  list(key="ATAC_3T3",  file=file.path(GSE211123, "GSE211123_3T3_control_ATAC_peaks.narrowPeak"),   type="ATAC", group="3T3",  label="3T3 ATAC"),
  list(key="scG4_cl1",  file=CLUSTER_PEAKS_1,           type="G4",   group="mESC", label="scG4 cl1"),
  list(key="ATAC_mESC", file=file.path(GSE149080, "GSE149080_ESC_WT_batch2_peaks.narrowPeak"),      type="ATAC", group="mESC", label="mESC ATAC"))
sets <- setNames(lapply(meta, function(m) np_summits(m$file)), sapply(meta, `[[`, "key"))
for (m in meta) cat(sprintf("%-10s %-9s %5d peaks\n", m$key, m$group, length(sets[[m$key]])))
pairs <- list(c(g="scG4_cl0", a="ATAC_3T3"), c(g="scG4_cl1", a="ATAC_mESC"))  # G4 vs paired ATAC

## ---- references -----------------------------------------------------------
cat("Loading references (genome-wide)...\n")
pqs <- onstd(import(PQS_BED, format = "BED")); pqs <- pqs[score(pqs) >= PQS_MIN]
eg4_tab <- read.table(EG4_PATH, sep = "\t", header = TRUE, quote = "", comment.char = "")
eg4 <- onstd(GRanges(eg4_tab[[1]], IRanges(eg4_tab[[2]] + 1, eg4_tab[[3]])))
refs <- list(PQS = pqs, eG4 = eg4)
cat(sprintf("  PQS(score>=%d)=%d | eG4=%d\n", PQS_MIN, length(pqs), length(eg4)))

## ---- permTest (fold + CI) + cache indicators ------------------------------
ind <- list(); base <- list()
for (rn in names(refs)) {
  cat(sprintf("\n=== %s ===\n", rn)); B <- refs[[rn]]
  for (k in names(sets)) {
    A <- sets[[k]]
    no <- permTest(A=A, B=B, randomize.function=circularRandomizeRegions,
                   evaluate.function=numOverlaps, genome=genome, count.once=TRUE,
                   ntimes=NTIMES, mc.cores=MC, force.parallel=TRUE)$numOverlaps
    em <- mean(no$permuted); ic <- as.integer(overlapsAny(A,B)); n <- length(ic)
    bci <- quantile(vapply(seq_len(NBOOT_CI), function(i) sum(ic[sample.int(n,n,TRUE)]), numeric(1))/em,
                    c(.025,.975), names=FALSE)
    ind[[paste(rn,k)]]  <- ic
    base[[paste(rn,k)]] <- list(n=n, obs=no$observed, exp=em, z=no$zscore, p=no$pval,
                                fold=no$observed/em, lo=bci[1], hi=bci[2])
    cat(sprintf("  %-10s n=%5d fold=%6.2fx [%.2f-%.2f]\n", k, n, no$observed/em, bci[1], bci[2]))
  }
}

## ---- significance + results table -----------------------------------------
mkey <- setNames(meta, sapply(meta, `[[`, "key"))
pair_for <- setNames(sapply(pairs, function(p) p["a"]), sapply(pairs, function(p) p["g"]))
out <- list()
for (rn in names(refs)) for (k in names(sets)) {
  b <- base[[paste(rn,k)]]; m <- mkey[[k]]
  logp <- pnorm(b$z, lower.tail=FALSE, log.p=TRUE)
  p_enr <- if (logp/log(10) < -300) "<1e-300" else formatC(exp(logp), format="e", digits=2)
  p_dif <- NA; p_prop <- NA; sg <- ""
  if (m$type == "G4" && k %in% names(pair_for)) {           # G4 vs its paired ATAC
    ak <- pair_for[[k]]; iG <- ind[[paste(rn,k)]]; iA <- ind[[paste(rn,ak)]]
    ba <- base[[paste(rn,ak)]]
    bd <- vapply(seq_len(NBOOT_DIF), function(j)
      (sum(iG[sample.int(b$n,b$n,TRUE)])/b$exp) - (sum(iA[sample.int(ba$n,ba$n,TRUE)])/ba$exp), numeric(1))
    p_dif  <- max(mean(bd <= 0), 1/NBOOT_DIF)
    p_prop <- prop.test(c(b$obs, ba$obs), c(b$n, ba$n), alternative="greater")$p.value
    sg <- stars(p_dif)
  }
  out[[length(out)+1]] <- data.frame(
    reference=rn, key=k, label=m$label, group=m$group, type=m$type,
    n_peaks=b$n, observed=b$obs, exp_mean=round(b$exp,2),
    fold=round(b$fold,3), fold_lo=round(b$lo,3), fold_hi=round(b$hi,3),
    log2fold=round(log2(b$fold),3), perm_z=round(b$z,1), perm_p=signif(b$p,3),
    p_enrichment=p_enr,
    p_G4_gt_ATAC_boot = if (is.na(p_dif)) NA else signif(p_dif,3),
    p_G4_gt_ATAC_prop = if (is.na(p_prop)) NA else signif(p_prop,3),
    sig = sg)
}
res <- do.call(rbind, out)
write.csv(res, file.path(OUTDIR, "scG4cluster_vs_bulkATAC_PQS_eG4_results.csv"), row.names=FALSE)
cat("\n=== Results ===\n"); print(res[,c("reference","label","group","type","fold","fold_lo","fold_hi","p_enrichment","p_G4_gt_ATAC_boot","sig")], row.names=FALSE)

## ---- plots: grouped bars per cell type, G4 vs ATAC, sig bracket -----------
res$group <- factor(res$group, levels=c("3T3","mESC"))
res$type  <- factor(res$type,  levels=c("G4","ATAC"))
DODGE <- 0.45; dodge <- position_dodge(width=DODGE); HALF <- DODGE/4

make_plot <- function(rn) {
  d <- res[res$reference==rn, ]
  # one bracket+star per group, from the G4 (cluster) row
  br <- do.call(rbind, lapply(levels(d$group), function(g){
    gg <- d[d$group==g, ]; xi <- as.integer(factor(g, levels=levels(d$group)))
    y  <- max(gg$fold_hi)*1.06
    data.frame(group=g, x=xi-HALF, xend=xi+HALF, y=y, lab=gg$sig[gg$type=="G4"])
  }))
  ggplot(d, aes(group, fold, fill=type)) +
    geom_col(position=dodge, width=0.40, color="black") +
    geom_errorbar(aes(ymin=fold_lo, ymax=fold_hi), position=dodge, width=0.12, linewidth=0.4) +
    geom_text(aes(label=sprintf("%.2fx", fold)), position=dodge, vjust=-0.6, size=3.0) +
    geom_segment(data=br, aes(x=x, xend=xend, y=y, yend=y), inherit.aes=FALSE, linewidth=0.4) +
    geom_text(data=br, aes(x=(x+xend)/2, y=y, label=lab), inherit.aes=FALSE, vjust=-0.3, fontface="bold", size=4) +
    geom_hline(yintercept=1, linetype="dashed") +
    scale_fill_manual(values=c(G4="#B2182B", ATAC="#2166AC"), labels=c(G4="scG4", ATAC="ATAC")) +
    scale_y_continuous(expand=expansion(mult=c(0,0.16))) +
    labs(title=sprintf("Fold enrichment over random - %s", rn),
         subtitle=paste0("scG4 pseudobulk cluster vs matched bulk ATAC\n",
                         "genome-wide | 1000 perms, 95% bootstrap CI | bracket: scG4 > ATAC (bootstrap, ***p<0.001)"),
         x=NULL, y="Fold enrichment over random", fill=NULL) +
    theme_bw(base_size=13) +
    theme(legend.position="top", plot.title=element_text(face="bold"),
          plot.subtitle=element_text(size=9.5, lineheight=1.1),
          axis.text.x=element_text(size=12, face="bold"))
}
for (rn in c("PQS","eG4")) {
  f <- file.path(OUTDIR, sprintf("scG4_vs_bulkATAC_%s_barplot.pdf", rn))
  ggsave(f, make_plot(rn), width=6.5, height=6); cat("Wrote", f, "\n")
}
cat("\nDone.\n")
