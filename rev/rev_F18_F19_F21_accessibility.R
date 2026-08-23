#!/usr/bin/env Rscript
# =============================================================
# Reviewer Q6 — scG4 signal is distinct from chromatin accessibility (F18, F19, F21).
#
# Shows the scG4 CUT&Tag signal is not a Tn5 accessibility artifact:
#   F18  scG4 vs ATAC peak-overlap Venn (a=mESC, b=MEF): a large fraction of
#        scG4 peaks have NO ATAC peak (accessibility-independent).
#   F19a G4-reference (eG4/PQS) enrichment by peak class: scG4-marked sites are
#        G4-enriched independent of accessibility, and accessible peaks WITHOUT
#        scG4 (ATAC-only) are G4-poor -> scG4 reads the G4 motif, not accessibility.
#   F19b scG4 vs ATAC RPGC signal by peak class: scG4-only peaks carry strong
#        scG4 signal at low ATAC.
#   F21  Is the eG4 overlap of scG4-only peaks "low"? Benchmark vs bulk G4 CUT&Tag
#        (gold standard) and random: eG4 is a sparse, accessibility-biased
#        reference (even bulk G4 reaches only ~32-36%); scG4-only is 21-40x over random.
#
# Peak classes (summit +/-150 bp): scG4&ATAC | scG4-only | ATAC-only | random.
# Cluster 1 = mESC, cluster 0 = MEF/3T3.
#
# Self-contained; reads only input_files/. Requires R >= 4.2 with the bundle deps
# (GenomicRanges, rtracklayer, Biostrings, BSgenome.Mmusculus.UCSC.mm10, pqsfinder,
# eulerr, ggplot2). BigWig inputs are large (hosted with the bundle's other tracks).
# Run from repo root:  Rscript scripts/34_G4_vs_accessibility_decoupling.R
# =============================================================

source("rev/paths.R")
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer)
  library(Biostrings); library(eulerr); library(ggplot2)
})
BSGENOME_AVAILABLE  <- requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)
PQSFINDER_AVAILABLE <- requireNamespace("pqsfinder", quietly = TRUE)
if (BSGENOME_AVAILABLE)  suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
if (PQSFINDER_AVAILABLE) suppressPackageStartupMessages(library(pqsfinder))
set.seed(1)
OUT_DIR  <- file.path(ROOT, "rev", "outputs", "G4_vs_accessibility_decoupling")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
inp <- function(fn) {
  if (file.exists(fn)) return(fn)
  hits <- list.files(DATA, pattern = basename(fn), recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("File not found in data/: ", fn)
  hits[1]
}
CANON <- c(paste0("chr",1:19),"chrX","chrY"); HALF <- 150L
if (BSGENOME_AVAILABLE) {
  GN  <- BSgenome.Mmusculus.UCSC.mm10
  GLEN <- seqlengths(GN)[CANON]
} else {
  if (file.exists(CHROM_SIZES)) {
    cs <- read.table(CHROM_SIZES, header = FALSE, stringsAsFactors = FALSE)
    GLEN <- setNames(as.integer(cs$V2), cs$V1)[CANON]
    message("Loaded chromosome lengths from ", CHROM_SIZES)
  } else {
    stop("BSgenome not available and chrom.sizes file not found: ", CHROM_SIZES)
  }
  GN <- NULL
}

read_summit <- function(fn){ d<-fread(inp(fn)); d<-d[V1 %in% CANON]
  off <- if(ncol(d)>=10) d$V10 else round((d$V3-d$V2)/2); s<-d$V2+off
  gr<-GRanges(d$V1, IRanges(pmax(s-HALF,1L), s+HALF)); seqlengths(gr)<-GLEN[seqlevels(gr)]; trim(gr) }
read_center <- function(fn){ d<-fread(inp(fn)); d<-d[V1 %in% CANON]; c<-round((d$V2+d$V3)/2)
  gr<-GRanges(d$V1, IRanges(pmax(c-HALF,1L), c+HALF)); seqlengths(gr)<-GLEN[seqlevels(gr)]; trim(gr) }
eg4 <- fread(EG4_PATH, header=TRUE); eg4<-eg4[Chr %in% CANON]
eg4_gr <- GRanges(eg4$Chr, IRanges(eg4$Start+1L, eg4$End))
pqs <- fread(PQS_BED, header=FALSE); pqs<-pqs[V1 %in% CANON & V5>=50]
pqs_gr <- GRanges(pqs$V1, IRanges(pqs$V2+1L, pqs$V3))
bw_mean <- function(fn, gr) as.numeric(unlist(summary(BigWigFile(fn), gr, size=1L, type="mean", defaultValue=0))$score)
pct <- function(a,b) 100*mean(overlapsAny(a,b,ignore.strand=TRUE))
rand_regions <- function(n){ chr<-sample(names(GLEN),n,replace=TRUE,prob=as.numeric(GLEN))
  st<-floor(runif(n)*(as.numeric(GLEN[chr])-2L*HALF-2L))+HALF+1L; GRanges(chr, IRanges(st-HALF, st+HALF)) }

cells <- list(
  mESC=list(scg4_pk=CLUSTER_PEAKS_1, atac_pk="GSE149080_ESC_WT_batch2_peaks.narrowPeak",
            bulk_pk=file.path(BULK_PEAKS, "GSM8836082_bulkG4CnT_mESC_rep1.broadPeak"), scg4_bw=BW_CL1, atac_bw=file.path(GSE149080, "GSM4661960_ATAC_ESC_WT_batch2.rpgc.bw")),
  MEF =list(scg4_pk=CLUSTER_PEAKS_0, atac_pk="GSE211123_3T3_control_ATAC_peaks.narrowPeak",
            bulk_pk=file.path(BULK_PEAKS, "GSM8836084_bulkG4CnT_3T3_rep1.broadPeak"), scg4_bw=BW_CL0, atac_bw=file.path(GSE211123, "GSM6451000_ATAC_3T3.rpgc.bw"))
)
CL_LEV<-c("scG4&ATAC","scG4-only","ATAC-only","random")
CL_LAB<-c("scG4 & ATAC","scG4-only\n(G4, no access.)","ATAC-only\n(access., no G4)","random")

cls_tab<-list(); ctx<-list(); venn<-list()
for (cl in names(cells)) {
  cc<-cells[[cl]]
  scg4<-read_summit(cc$scg4_pk); atac<-read_summit(cc$atac_pk); bulk<-read_center(cc$bulk_pk); rnd<-rand_regions(20000)
  ov<-overlapsAny(scg4,atac,ignore.strand=TRUE)
  cls<-list(`scG4&ATAC`=scg4[ov], `scG4-only`=scg4[!ov], `ATAC-only`=atac[!overlapsAny(atac,scg4,ignore.strand=TRUE)], `random`=rnd)
  venn[[cl]]<-c(scG4_only=sum(!ov), inter=sum(ov), ATAC_only=length(cls$`ATAC-only`), scG4=length(scg4), ATAC=length(atac))
  for (nm in names(cls)){ g<-cls[[nm]]
    cls_tab[[length(cls_tab)+1]]<-data.table(cell=cl, class=nm, n=length(g),
      pct_eG4=pct(g,eg4_gr), pct_PQS=pct(g,pqs_gr),
      scg4_signal=mean(bw_mean(cc$scg4_bw,g)), atac_signal=mean(bw_mean(cc$atac_bw,g))) }
  ctx[[cl]]<-data.table(cell=cl, set=c("bulk G4 CUT&Tag peaks","all scG4 peaks","scG4-only (no ATAC)","random"),
    eG4_pct=c(pct(bulk,eg4_gr),pct(scg4,eg4_gr),pct(cls$`scG4-only`,eg4_gr),pct(rnd,eg4_gr)),
    PQS_pct=c(pct(bulk,pqs_gr),pct(scg4,pqs_gr),pct(cls$`scG4-only`,pqs_gr),pct(rnd,pqs_gr)))
}
D<-rbindlist(cls_tab); D[,class:=factor(class,levels=CL_LEV)]
CTX<-rbindlist(ctx)
fwrite(D, file.path(OUT_DIR,"peak_class_motif_accessibility.csv"))
fwrite(CTX, file.path(OUT_DIR,"eG4_overlap_context.csv"))

sp<-function(p,id,w,h){ ggsave(file.path(OUT_DIR,paste0(id,".pdf")),p,width=w,height=h) }

# ---- F18 Venns (eulerr) ----
for (cl in names(cells)){ v<-venn[[cl]]
  fit<-euler(c(scG4=unname(v["scG4_only"]), ATAC=unname(v["ATAC_only"]), "scG4&ATAC"=unname(v["inter"])))
  p<-plot(fit, quantities=list(cex=1.1), fills=list(fill=c("#8e44ad","#2980b9"),alpha=.55),
          labels=list(cex=1.2,font=2), main=list(label=sprintf("%s scG4 vs ATAC peaks", cl), cex=1.05))
  id<-sprintf("F18%s_scG4_vs_ATAC_venn_%s", if(cl=="mESC")"a" else "b", cl)
  fn<-file.path(OUT_DIR,paste0(id,".pdf")); pdf(fn,width=6.2,height=5); print(p); dev.off() }

# ---- F19a eG4/PQS by class ----
m<-melt(D, id.vars=c("cell","class"), measure.vars=c("pct_eG4","pct_PQS"), variable.name="ref", value.name="pct")
m[,ref:=factor(ref,labels=c("eG4 (validated G4)","PQS (G4 motif)"))]
sp(ggplot(m,aes(class,pct,fill=ref))+geom_col(position=position_dodge())+facet_wrap(~cell)+
   scale_x_discrete(labels=CL_LAB)+scale_fill_manual(values=c("eG4 (validated G4)"="#c0392b","PQS (G4 motif)"="#e59866"),name=NULL)+
   labs(x=NULL,y="% peaks overlapping G4 reference",title="scG4-marked sites are G4-enriched independent of accessibility",
        subtitle="Accessible peaks WITHOUT scG4 (ATAC-only) are G4-poor")+theme_bw(base_size=10)+
   theme(plot.title=element_text(face="bold"),axis.text.x=element_text(size=8),legend.position="top"),
   "F19a_accessibility_decoupling_G4enrichment_by_class",10,5.5)

# ---- F19b signal by class ----
s<-melt(D, id.vars=c("cell","class"), measure.vars=c("scg4_signal","atac_signal"), variable.name="assay", value.name="signal")
s[,assay:=factor(assay,labels=c("scG4 CUT&Tag","ATAC-seq"))]
sp(ggplot(s,aes(class,signal,fill=assay))+geom_col(position=position_dodge())+facet_wrap(~cell)+
   scale_x_discrete(labels=CL_LAB)+scale_fill_manual(values=c("scG4 CUT&Tag"="#8e44ad","ATAC-seq"="#2980b9"),name=NULL)+
   labs(x=NULL,y="mean RPGC signal",title="scG4-only peaks carry strong G4 signal at low accessibility")+
   theme_bw(base_size=10)+theme(plot.title=element_text(face="bold"),axis.text.x=element_text(size=8),legend.position="top"),
   "F19b_accessibility_decoupling_signal_by_class",10,5.5)

# ---- F21 eG4 context benchmark ----
CTX[,set:=factor(set,levels=c("bulk G4 CUT&Tag peaks","all scG4 peaks","scG4-only (no ATAC)","random"))]
mc<-melt(CTX,id.vars=c("cell","set"),measure.vars=c("eG4_pct","PQS_pct"),variable.name="ref",value.name="pct")
mc[,ref:=factor(ref,labels=c("eG4 (validated G4)","PQS (G4 motif)"))]
sp(ggplot(mc,aes(set,pct,fill=ref))+geom_col(position=position_dodge(width=.8),width=.72)+
   geom_text(aes(label=sprintf("%.0f%%",pct)),position=position_dodge(width=.8),vjust=-.3,size=2.9)+facet_wrap(~cell)+
   scale_fill_manual(values=c("eG4 (validated G4)"="#c0392b","PQS (G4 motif)"="#e59866"),name=NULL)+
   labs(x=NULL,y="% peaks overlapping G4 reference",title="Accessibility-independent scG4 peaks are bona-fide G4 sites",
        subtitle="Even bulk G4 CUT&Tag reaches only ~32-36% eG4; scG4-only is 21-40x over random (0.5%)")+
   theme_bw(base_size=10)+theme(plot.title=element_text(face="bold"),plot.subtitle=element_text(size=8.2),
        axis.text.x=element_text(angle=25,hjust=1,size=8),legend.position="top"),
   "F21_scG4only_eG4_overlap_benchmark",11,5.5)

cat("Done. F18/F19/F21 outputs in", OUT_DIR, "\n")
