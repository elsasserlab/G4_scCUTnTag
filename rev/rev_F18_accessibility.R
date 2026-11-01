#!/usr/bin/env Rscript
# =============================================================
# Reviewer Q6 - scG4 signal is distinct from chromatin accessibility (F18).
#
# Shows the scG4 CUT&Tag signal is not a Tn5 accessibility artifact:
#   F18  scG4 vs ATAC peak-overlap Venn (a=mESC, b=MEF): a large fraction of
#        scG4 peaks have NO ATAC peak (accessibility-independent).
#
# Peak classes (summit +/-150 bp).
# Cluster 1 = mESC, cluster 0 = MEF/3T3.
#
# Run from repo root:  Rscript rev/rev_F18_accessibility.R
# =============================================================

source("rev/paths.R")
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(eulerr)
})
OUT_DIR  <- file.path(ROOT, "rev", "outputs", "G4_vs_accessibility_decoupling")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
inp <- function(fn) {
  if (file.exists(fn)) return(fn)
  hits <- list.files(DATA, pattern = basename(fn), recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("File not found in data/: ", fn)
  hits[1]
}
CANON <- c(paste0("chr",1:19),"chrX","chrY"); HALF <- 150L
if (file.exists(CHROM_SIZES)) {
  cs <- read.table(CHROM_SIZES, header = FALSE, stringsAsFactors = FALSE)
  GLEN <- setNames(as.integer(cs$V2), cs$V1)[CANON]
} else {
  stop("chrom.sizes file not found: ", CHROM_SIZES)
}

read_summit <- function(fn){ d<-fread(inp(fn)); d<-d[V1 %in% CANON]
  off <- if(ncol(d)>=10) d$V10 else round((d$V3-d$V2)/2); s<-d$V2+off
  gr<-GRanges(d$V1, IRanges(pmax(s-HALF,1L), s+HALF)); seqlengths(gr)<-GLEN[seqlevels(gr)]; trim(gr) }

cells <- list(
  mESC=list(scg4_pk=CLUSTER_PEAKS_1, atac_pk="GSE149080_ESC_WT_batch2_peaks.narrowPeak"),
  MEF =list(scg4_pk=CLUSTER_PEAKS_0, atac_pk="GSE211123_3T3_control_ATAC_peaks.narrowPeak")
)

venn<-list()
for (cl in names(cells)) {
  cc<-cells[[cl]]
  scg4<-read_summit(cc$scg4_pk); atac<-read_summit(cc$atac_pk)
  ov<-overlapsAny(scg4,atac,ignore.strand=TRUE)
  venn[[cl]]<-c(scG4_only=sum(!ov), inter=sum(ov),
                ATAC_only=sum(!overlapsAny(atac,scg4,ignore.strand=TRUE)),
                scG4=length(scg4), ATAC=length(atac))
}

# ---- F18 Venns (eulerr) ----
for (cl in names(cells)){ v<-venn[[cl]]
  fit<-euler(c(scG4=unname(v["scG4_only"]), ATAC=unname(v["ATAC_only"]), "scG4&ATAC"=unname(v["inter"])))
  p<-plot(fit, quantities=list(cex=1.1), fills=list(fill=c("#8e44ad","#2980b9"),alpha=.55),
          labels=list(cex=1.2,font=2), main=list(label=sprintf("%s scG4 vs ATAC peaks", cl), cex=1.05))
  id<-sprintf("F18%s_scG4_vs_ATAC_venn_%s", if(cl=="mESC")"a" else "b", cl)
  fn<-file.path(OUT_DIR,paste0(id,".pdf")); pdf(fn,width=6.2,height=5); print(p); dev.off() }

cat("Done. F18 outputs in", OUT_DIR, "\n")
