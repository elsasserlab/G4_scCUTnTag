#!/usr/bin/env Rscript
# =============================================================
# Reviewer QC — per-cell data quality of the primary NEURAL scG4 dataset (F23).
#
# Documents reads-per-cell and FRiP (fraction of reads in peaks) for the
# unsorted mouse-brain scG4 CUT&Tag data, matching the QC reported for the
# mESC-MEF cell-line data. Metrics are taken per cell from the Cell Ranger ATAC
# QC fields carried in the brain scG4 Seurat object.
#
#   F23a  per-cell QC violins for the unsorted brain (reads/cell, FRiP,
#         TSS fragments, peaks detected), split by scG4 cluster
#   F23b  cross-dataset comparison of median HQ fragments/cell and FRiP
#         (mESC-MEF cell line; sorted / unsorted brain) — Cell Ranger values
#
# NOTE: reads/cell = passed_filters (unique high-quality fragments per cell);
# FRiP = peak_region_fragments / passed_filters. The sorted (GFP+) values in
# F23b are the Cell Ranger numbers reported in the manuscript (no per-cell GFP+
# object is bundled); the unsorted per-cell distributions (F23a) are computed
# here. Over the QC-passing analysed cells the unsorted medians are higher than
# the Cell Ranger raw median (2,194 vs 1,511 fragments/cell) because filtering
# removes low-quality barcodes; FRiP agrees (17.1% vs 16.3%).
#
# Self-contained; reads only input_files/. Requires R >= 4.2 with Seurat,
# Signac, data.table, ggplot2. Run from repo root:
#   Rscript scripts/37_neural_scG4_QC.R
# =============================================================

source("rev/paths.R")
suppressPackageStartupMessages({ library(Seurat); library(Signac); library(data.table); library(ggplot2) })
DATA  <- file.path(ROOT, "input_files")
OUT_DIR <- file.path(ROOT, "rev", "outputs", "neural_scG4_QC")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
inp <- function(fn) {
  hits <- list.files(DATA, pattern = fn, recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("File not found in data/: ", fn)
  hits[1]
}

o  <- readRDS(file.path(DATA, "GSE291468/GSM8836087_unsorted_Seurat_object.Rds"))
md <- as.data.table(o@meta.data)
md[, reads_per_cell := passed_filters]
md[, FRiP := peak_region_fragments / passed_filters]
md[, cluster := as.character(seurat_clusters)]
fwrite(md[, .(cluster, total, passed_filters, peak_region_fragments,
              reads_per_cell, FRiP, TSS_fragments, nCount_peaks)],
       file.path(OUT_DIR, "neural_scG4_per_cell_QC.csv"))
cat(sprintf("Unsorted brain scG4: n=%d cells | median reads/cell=%.0f | median FRiP=%.1f%%\n",
            nrow(md), median(md$reads_per_cell), 100*median(md$FRiP)))

save_fig <- function(p, id, w, h){
  ggsave(file.path(OUT_DIR, paste0(id, ".pdf")), p, width=w, height=h, dpi=200)
}

# ---- F23a: per-cell QC violins ----
md[, cluster := factor(cluster)]
L <- rbindlist(list(
  md[, .(cluster, metric="Reads/cell (log10)",            value=log10(reads_per_cell))],
  md[, .(cluster, metric="FRiP (fraction reads in peaks)", value=FRiP)],
  md[, .(cluster, metric="TSS fragments/cell (log10)",     value=log10(pmax(TSS_fragments,1)))],
  md[, .(cluster, metric="Peaks detected/cell (log10)",    value=log10(pmax(nCount_peaks,1)))]))
L[, metric := factor(metric, levels=c("Reads/cell (log10)","FRiP (fraction reads in peaks)",
                                       "TSS fragments/cell (log10)","Peaks detected/cell (log10)"))]
p1 <- ggplot(L, aes(cluster, value, fill=cluster)) +
  geom_violin(alpha=.55, colour=NA, scale="width") +
  geom_boxplot(width=.16, outlier.shape=NA, alpha=.85) +
  facet_wrap(~metric, scales="free_y", nrow=1) +
  scale_fill_manual(values=c("0"="#4c78a8","1"="#f58518"), guide="none") +
  labs(x="scG4 cluster", y=NULL,
       title=sprintf("Unsorted mouse-brain scG4 CUT&Tag - per-cell data quality (n = %d cells)", nrow(md)),
       subtitle=sprintf("Median %.0f reads/cell and %.1f%% FRiP over analysed cells (Cell Ranger raw median 1,511 / 16.3%%)",
                        median(md$reads_per_cell), 100*median(md$FRiP))) +
  theme_bw(base_size=10) + theme(plot.title=element_text(face="bold"), plot.subtitle=element_text(size=8))
save_fig(p1, "F23a_neural_scG4_per_cell_QC", 12, 4)

# ---- F23b: cross-dataset QC (Cell Ranger ATAC values from the manuscript) ----
comp <- data.table(dataset=factor(c("mESC-MEF\n(cell line)","sorted brain\n(GFP+)","unsorted brain"),
                     levels=c("mESC-MEF\n(cell line)","sorted brain\n(GFP+)","unsorted brain")),
                   median_fragments=c(3434,1983,1511), FRiP_pct=c(28.3,13.5,16.3))
fwrite(comp, file.path(OUT_DIR, "scG4_QC_across_datasets.csv"))
cm <- melt(comp, id.vars="dataset"); cm[, variable := factor(variable, labels=c("median HQ fragments / cell","FRiP (%)"))]
p2 <- ggplot(cm, aes(dataset, value, fill=dataset)) + geom_col(width=.65) +
  geom_text(aes(label=ifelse(grepl("FRiP",variable), sprintf("%.1f%%",value), format(value,big.mark=","))), vjust=-.3, size=3) +
  facet_wrap(~variable, scales="free_y") +
  scale_fill_manual(values=c("#54a24b","#e45756","#b279a2"), guide="none") +
  labs(x=NULL, y=NULL, title="scG4 CUT&Tag data quality across datasets (Cell Ranger ATAC)",
       subtitle="Cell-line (mESC-MEF) and both primary neural datasets (sorted / unsorted brain)") +
  theme_bw(base_size=10) + theme(plot.title=element_text(face="bold"), axis.text.x=element_text(size=8))
save_fig(p2, "F23b_scG4_QC_across_datasets", 9, 4.2)
cat("Done. F23a/F23b outputs in", OUT_DIR, "\n")
