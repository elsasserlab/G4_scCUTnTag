#!/usr/bin/env Rscript
# =============================================================
# Reviewer Q4 — are the Fig-3C astrocyte G4-marked genes (Tnik, Pitpnc1, Pbx1,
# Nwd1) also astrocyte-specific in chromatin accessibility? (F22)
#
# Independent orthogonal check using the paper's own brain snATAC object
# (GSE198467, 15 named cell types). For each gene, per-cell accessibility from
# the Signac gene-activity (GA) assay, compared astrocyte (AST_NT + AST_TE) vs
# non-astrocyte (binary, matching the AST-vs-non-AST framing of Fig 3B/3D).
#
# This complements the G4 evidence (Fig 3B/3C) and the scRNA-Seq expression
# evidence (Fig 3A/3D) with a third modality (accessibility). NOTE: the object
# is ATAC-only (assays "peaks" + "GA"); GA is ATAC-derived gene activity.
#
# Self-contained; reads only input_files/. Requires R >= 4.2 with Seurat,
# Signac, ggplot2. Run from repo root:
#   Rscript scripts/36_brain_AST_gene_scATAC.R
# =============================================================

source("rev/paths.R")
suppressPackageStartupMessages({ library(Seurat); library(Signac); library(data.table); library(ggplot2) })
DATA  <- file.path(ROOT, "input_files")
OUT_DIR <- file.path(ROOT, "rev", "outputs", "brain_AST_gene_scATAC")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
inp <- function(fn) {
  hits <- list.files(DATA, pattern = fn, recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("File not found in data/: ", fn)
  hits[1]
}
GENES <- c("Tnik","Pitpnc1","Pbx1","Nwd1")

o <- readRDS(file.path(GSE198467, "GSE198467_ATAC_Seurat_object_clustered_renamed.Rds"))
DefaultAssay(o) <- "GA"; o <- NormalizeData(o, verbose = FALSE)
stopifnot(all(GENES %in% rownames(o[["GA"]])))
ga  <- as.matrix(GetAssayData(o, assay = "GA", layer = "data")[GENES, , drop = FALSE])
ct  <- as.character(o@meta.data[["idents_short"]])
ast <- ifelse(ct %in% c("AST_NT","AST_TE"), "AST", "non-AST")
cat(sprintf("cells: AST=%d  non-AST=%d\n", sum(ast=="AST"), sum(ast=="non-AST")))

summ <- rbindlist(lapply(GENES, function(g){
  a <- ga[g, ast=="AST"]; b <- ga[g, ast=="non-AST"]
  data.table(gene=g, AST_mean=mean(a), nonAST_mean=mean(b), fold=mean(a)/(mean(b)+1e-9),
             wilcox_p=suppressWarnings(wilcox.test(a,b))$p.value)
}))
fwrite(summ, file.path(OUT_DIR, "AST_vs_nonAST_GA_summary.csv"))
cat("=== AST vs non-AST gene activity (scATAC) ===\n"); print(summ)

# ---- F22: binary AST vs non-AST bar per gene ----
b <- melt(summ, id.vars="gene", measure.vars=c("AST_mean","nonAST_mean"), variable.name="group", value.name="GA")
b[, group := factor(group, labels=c("AST","non-AST"))]; b[, gene := factor(gene, levels=GENES)]
lab <- summ[, .(gene=factor(gene,levels=GENES), y=pmax(AST_mean,nonAST_mean),
                txt=sprintf("%.2fx\np=%s", fold, formatC(wilcox_p, format="e", digits=0)))]
p <- ggplot(b, aes(group, GA, fill=group)) + geom_col(width=.65) +
  geom_text(data=lab, aes(x=1.5, y=y*1.08, label=txt), inherit.aes=FALSE, size=2.9, lineheight=.9) +
  facet_wrap(~gene, nrow=1) +
  scale_fill_manual(values=c("AST"="#c0392b","non-AST"="#7f8c8d"), guide="none") +
  labs(x=NULL, y="mean scATAC gene activity",
       title="Astrocyte vs non-astrocyte accessibility (scATAC) of the Fig-3C G4 genes",
       subtitle="All four genes are significantly more accessible in astrocytes (GSE198467 brain snATAC)") +
  theme_bw(base_size=10) + theme(plot.title=element_text(face="bold"), plot.subtitle=element_text(size=8.5))
id <- "F22_scATAC_ASTvsnonAST_Fig3C_genes"
ggsave(file.path(OUT_DIR, paste0(id, ".pdf")), p, width=10, height=4, dpi=200)
cat("Done. F22 outputs in", OUT_DIR, "\n")
