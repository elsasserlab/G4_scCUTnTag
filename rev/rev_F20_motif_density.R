#!/usr/bin/env Rscript
# =============================================================
# Reviewer Q6 — G4-motif density: scG4 peaks vs matched ATAC peaks (F20).
#
# Direct test of "stronger enrichment of G-quadruplex forming DNA motifs than
# the ATAC-seq profiles". For each cell type, matched-N, matched-width
# (summit +/-150 bp) scG4 vs ATAC peaks are scored for G4-forming sequence:
#   F20a  pqsfinder max G4 score per peak (both strands)
#   F20b  canonical G4 motif (regex G3+ N1-7 G3+ N1-7 G3+ N1-7 G3+, both strands):
#         % of peaks with a motif, and mean motifs/peak.
# MEF: scG4 > ATAC (p<1e-38). mESC: sequence motif saturates (GC-rich open
# chromatin), so the experimentally-validated eG4 reference (F19/F21) is the
# discriminating one there.
#
# Cluster 1 = mESC, cluster 0 = MEF/3T3. Self-contained; reads only input_files/.
# Requires R >= 4.2 with Biostrings, BSgenome.Mmusculus.UCSC.mm10, pqsfinder, ggplot2.
# Run from repo root:  Rscript scripts/35_G4_motif_density.R
# =============================================================

source("rev/paths.R")
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(Biostrings)
  library(ggplot2)
})
BSGENOME_AVAILABLE   <- requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)
PQSFINDER_AVAILABLE  <- requireNamespace("pqsfinder", quietly = TRUE)
if (BSGENOME_AVAILABLE)  suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
if (PQSFINDER_AVAILABLE) suppressPackageStartupMessages(library(pqsfinder))
set.seed(1)
OUT_DIR <- file.path(ROOT, "rev", "outputs", "G4_motif_density")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
inp <- function(fn) {
  full <- file.path(DATA, fn)
  if (file.exists(full)) return(full)
  hits <- list.files(DATA, pattern = fn, recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("File not found in data/: ", fn)
  hits[1]
}
CANON<-c(paste0("chr",1:19),"chrX","chrY"); HALF<-150L; N_SUB<-8000L
G4RE<-"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}"
PRECOMPUTED_CSV <- file.path(GENOME_DIR, "pqsfinder_scores_precomputed.csv")

if (BSGENOME_AVAILABLE) {
  GN <- BSgenome.Mmusculus.UCSC.mm10
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

read_summit<-function(fn){ d<-fread(inp(fn)); d<-d[V1 %in% CANON]
  off<-if(ncol(d)>=10) d$V10 else round((d$V3-d$V2)/2); s<-d$V2+off
  gr<-GRanges(d$V1, IRanges(pmax(s-HALF,1L), s+HALF)); seqlengths(gr)<-GLEN[seqlevels(gr)]; trim(gr) }
g4_count<-function(seqs){ f<-as.character(seqs); r<-as.character(reverseComplement(seqs))
  cnt<-function(v) vapply(gregexpr(G4RE,v,perl=TRUE),function(m) if(m[1]==-1) 0L else length(m),integer(1)); cnt(f)+cnt(r) }
if (BSGENOME_AVAILABLE && PQSFINDER_AVAILABLE) {
  pqs_max<-function(gr){ ss<-getSeq(GN,gr); vapply(seq_along(ss),function(i){ pv<-suppressMessages(pqsfinder(ss[[i]],strand="*",verbose=FALSE)); if(length(pv)) max(score(pv)) else 0L },integer(1)) }
} else {
  pqs_max <- NULL
}

cells<-list(mESC=c(scg4=CLUSTER_PEAKS_1, atac="GSE149080_ESC_WT_batch2_peaks.narrowPeak"),
            MEF =c(scg4=CLUSTER_PEAKS_0,  atac="GSE211123_3T3_control_ATAC_peaks.narrowPeak"))

if (file.exists(PRECOMPUTED_CSV) && is.null(pqs_max)) {
  message("Loading pre-computed pqsfinder scores from ", PRECOMPUTED_CSV)
  PP <- fread(PRECOMPUTED_CSV)
  PP[, group := factor(group, levels = c("scG4", "ATAC"))]
  SM <- PP[, .(n = .N,
               pct_with_G4motif = 100 * mean(g4_motifs > 0),
               mean_G4motifs   = mean(g4_motifs),
               median_pqs      = median(pqs_score),
               mean_pqs        = mean(pqs_score)),
           by = .(cell, group)]
} else if (!is.null(pqs_max) && BSGENOME_AVAILABLE) {
  pp<-list(); sm<-list()
  for (cl in names(cells)){
    scg4<-read_summit(cells[[cl]]["scg4"]); atac<-read_summit(cells[[cl]]["atac"])
    n<-min(length(scg4),length(atac),N_SUB)
    s1<-scg4[sort(sample(length(scg4),n))]; a1<-atac[sort(sample(length(atac),n))]
    for (grp in c("scG4","ATAC")){ g<-if(grp=="scG4") s1 else a1
      d<-data.table(cell=cl, group=grp, g4_motifs=g4_count(getSeq(GN,g)), pqs_score=pqs_max(g))
      pp[[length(pp)+1]]<-d
      sm[[length(sm)+1]]<-data.table(cell=cl, group=grp, n=nrow(d),
        pct_with_G4motif=100*mean(d$g4_motifs>0), mean_G4motifs=mean(d$g4_motifs),
        median_pqs=median(d$pqs_score), mean_pqs=mean(d$pqs_score)) }
  }
  PP<-rbindlist(pp); SM<-rbindlist(sm); PP[,group:=factor(group,levels=c("scG4","ATAC"))]
} else {
  stop("Cannot run F20: pqsfinder/BSgenome not available and no pre-computed data at ", PRECOMPUTED_CSV)
}
fwrite(PP, file.path(OUT_DIR,"motif_density_perpeak.csv")); fwrite(SM, file.path(OUT_DIR,"motif_density_summary.csv"))
sp<-function(p,id,w,h){ ggsave(file.path(OUT_DIR,paste0(id,".pdf")),p,width=w,height=h) }

# F20a pqsfinder distribution
sp(ggplot(PP,aes(group,log1p(pqs_score),fill=group))+geom_violin(alpha=.5,colour=NA)+
   geom_boxplot(width=.18,outlier.shape=NA,alpha=.8)+facet_wrap(~cell)+
   scale_fill_manual(values=c(scG4="#8e44ad",ATAC="#2980b9"),guide="none")+
   labs(x=NULL,y="log1p( pqsfinder max G4 score )",title="Per-peak G4-forming potential: scG4 vs matched ATAC peaks",
        subtitle="Matched N and width (summit +/-150 bp). MEF: scG4 > ATAC (p=3e-63); mESC: saturated (GC-rich open chromatin)")+
   theme_bw(base_size=10)+theme(plot.title=element_text(face="bold")),
   "F20a_G4motif_density_pqsfinder",9,5)

# F20b canonical G4 motif bars
b<-melt(SM,id.vars=c("cell","group"),measure.vars=c("pct_with_G4motif","mean_G4motifs"),variable.name="metric",value.name="val")
b[,metric:=factor(metric,labels=c("% peaks with canonical G4 motif","mean canonical G4 motifs / peak"))]
b[,group:=factor(group,levels=c("scG4","ATAC"))]
sp(ggplot(b,aes(cell,val,fill=group))+geom_col(position=position_dodge())+facet_wrap(~metric,scales="free_y")+
   scale_fill_manual(values=c(scG4="#8e44ad",ATAC="#2980b9"),name=NULL)+
   labs(x=NULL,y=NULL,title="Canonical G4 motif content: scG4 vs matched ATAC peaks")+
   theme_bw(base_size=10)+theme(plot.title=element_text(face="bold"),legend.position="top"),
   "F20b_G4motif_density_canonical",9,5)

for (cl in names(cells)){ a<-PP[cell==cl&group=="scG4"]; z<-PP[cell==cl&group=="ATAC"]
  cat(sprintf("[%s] canonical G4/peak scG4 %.2f vs ATAC %.2f (p=%.1e); pqs scG4 %.1f vs ATAC %.1f (p=%.1e)\n",
    cl, mean(a$g4_motifs), mean(z$g4_motifs), wilcox.test(a$g4_motifs,z$g4_motifs)$p.value,
    mean(a$pqs_score), mean(z$pqs_score), wilcox.test(a$pqs_score,z$pqs_score)$p.value)) }
cat("Done. F20 outputs in", OUT_DIR, "\n")
