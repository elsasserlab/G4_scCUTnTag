#!/usr/bin/env Rscript
# =============================================================
# Cell-specific G4 peaks -> g:Profiler functional enrichment (F16, F17).
#
# Independent, multi-source (GO / KEGG / Reactome / WikiPathways / TF) answer
# to the reviewer question: "what is the functionality of the mESC-unique /
# MEF-unique G4 peak genes; are they related to distinct features of the two
# cell types (e.g. mESC self-renewal)?"
#
# Self-contained: reads ONLY the original input peak files + GTF, redoes the
# Fig-1D overlap analysis to extract the cluster-discriminating (unique) peaks
# (identical logic to 13_cell_specific_GO_enrichment.R), annotates them to
# nearest genes, and queries g:Profiler (gprofiler2::gost, whole-genome
# background, g:SCS correction, evcodes = TRUE) for each gene set.
#
# NOTE on peak vs region counts: the published Fig 1D reports 6,543 mESC-unique
# and 3,321 MEF-unique peaks (region-level eulerr inclusion-exclusion counts).
# The per-peak extractable sets are 6,221 / 2,816 peaks -> 5,522 / 2,594 nearest
# genes, which are what an enrichment tool can consume. This script uses those
# gene sets (same as F09/F10).
#
# Inputs  (input_files/):
#   cluster_spec_peaks/0_peaks.narrowPeak   (scG4 cluster 0 = MEF/3T3 context)
#   cluster_spec_peaks/1_peaks.narrowPeak   (scG4 cluster 1 = mESC context)
#   bulkG4CnT_mESC_rep1.broadPeak          (bulk mESC G4)
#   bulkG4CnT_3T3_rep1.broadPeak           (bulk MEF/3T3 G4)
#   mm10_.annotation.gtf.gz                (gene models for nearest-TSS)
#
# Outputs (rev/outputs/gprofiler2_cell_specific/):
#   gprofiler2_wholegenome_evcodes.csv                     (full result, all sources)
#   gprofiler2_mESCunique_selfRenewal.csv                  (self-renewal/signalling subset)
#   gprofiler2_MEFunique_fibroblast.csv                    (fibroblast/cytoskeleton/ECM subset)
#   F16_mESCunique_gProfiler_enrichment.{pdf,png}          
#   F17_MEFunique_gProfiler_enrichment.{pdf,png}
#
# Requires internet (queries the g:Profiler API) and R >= 4.2 with a matching
# ggplot2 (the bundle dependency set). Run from the repo root:
#   Rscript scripts/33_gprofiler2_cell_specific_enrichment.R
# =============================================================

source("rev/paths.R")
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(gprofiler2)
  library(ggplot2)
})

# ----- Paths -----
inp <- function(fn) {
  hits <- list.files(DATA, pattern = fn, recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("File not found in data/: ", fn)
  hits[1]
}
GTF_PATH <- file.path(GENOME_DIR, "mm10_.annotation.gtf.gz")
OUT_DIR  <- file.path(OUT_ROOT, "gprofiler2_cell_specific")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
CANONICAL <- c(paste0("chr", 1:19), "chrX", "chrY")

# ============================================================
# 1. Extract cluster-discriminating (unique) peaks  [same as script 13]
# ============================================================
read_peak_gr <- function(path, label) {
  dt <- fread(path)
  GRanges(seqnames = dt$V1, ranges = IRanges(dt$V2, dt$V3, names = rep(label, nrow(dt))))
}
keep_can <- function(gr) gr[as.character(seqnames(gr)) %in% CANONICAL]
cluster0 <- keep_can(read_peak_gr(CLUSTER_PEAKS_0, "c0"))
cluster1 <- keep_can(read_peak_gr(CLUSTER_PEAKS_1, "c1"))
mesc     <- keep_can(read_peak_gr(file.path(BULK_PEAKS, "GSM8836082_bulkG4CnT_mESC_rep1.broadPeak"),        "mESC"))
mef      <- keep_can(read_peak_gr(file.path(BULK_PEAKS, "GSM8836084_bulkG4CnT_3T3_rep1.broadPeak"),         "MEF"))

mesc_specific <- mesc[ overlapsAny(mesc, cluster1, ignore.strand = TRUE) &
                      !overlapsAny(mesc, cluster0, ignore.strand = TRUE)]
mef_specific  <- mef[  overlapsAny(mef,  cluster0, ignore.strand = TRUE) &
                      !overlapsAny(mef,  cluster1, ignore.strand = TRUE)]
cat(sprintf("Unique peaks: mESC=%d  MEF=%d (per-peak; Fig1D region counts 6543/3321)\n",
            length(mesc_specific), length(mef_specific)))

# ============================================================
# 2. Nearest-TSS gene per unique peak
# ============================================================
gtf <- rtracklayer::import(GTF_PATH); gtf <- gtf[as.character(seqnames(gtf)) %in% CANONICAL]
tx  <- gtf[gtf$type == "transcript"]
tdt <- data.table(chrom = as.character(seqnames(tx)), strand = as.character(strand(tx)),
                  start = start(tx), end = end(tx), gene_name = tx$gene_name, gene_id = tx$gene_id)
tdt[, tss := ifelse(strand == "+", start, end)]
gts <- tdt[, .(chrom = chrom[1], strand = strand[1],
               tss = ifelse(strand[1] == "+", min(tss), max(tss)), gene_name = gene_name[1]), by = gene_id]
gts <- gts[!is.na(gene_name) & nzchar(gene_name)]
tss_gr <- GRanges(gts$chrom, IRanges(gts$tss, gts$tss), gene_name = gts$gene_name)
nearest_genes <- function(q) {
  dn <- distanceToNearest(q, tss_gr, ignore.strand = TRUE)
  g <- unique(as.character(tss_gr$gene_name[subjectHits(dn)])); g[!is.na(g) & nzchar(g)]
}
mesc_genes <- nearest_genes(mesc_specific)
mef_genes  <- nearest_genes(mef_specific)
cat(sprintf("Nearest genes: mESC-unique=%d  MEF-unique=%d\n", length(mesc_genes), length(mef_genes)))

# ============================================================
# 3. g:Profiler enrichment (whole-genome background, g:SCS, evcodes)
# ============================================================
SOURCES <- c("GO:BP","GO:MF","GO:CC","KEGG","REAC","WP","TF","MIRNA","CORUM","HP")
g <- gost(query = list(`mESC-unique` = mesc_genes, `MEF-unique` = mef_genes),
          organism = "mmusculus", ordered_query = FALSE, significant = TRUE,
          user_threshold = 0.05, correction_method = "g_SCS",
          sources = SOURCES, evcodes = TRUE, domain_scope = "annotated")
res <- as.data.table(g$result)
if ("parents" %in% names(res)) res[, parents := sapply(parents, paste, collapse = "|")]
cat("g:Profiler data version:", g$meta$version, "\n")
fwrite(res, file.path(OUT_DIR, "gprofiler2_wholegenome_evcodes.csv"))

sr_pat  <- "stem cell|pluripot|self.renew|Wnt|SMAD|Nodal|TGF|BMP|cell fate|renewal|maintenance of cell|proliferation"
fib_pat <- "actin|collagen|extracellular matrix|adhesion|migrat|cytoskeleton|focal|integrin|ECM|fibrobl|ossif|mesenchym"
mesc_hi <- res[query == "mESC-unique" & grepl(sr_pat,  term_name, ignore.case = TRUE)][order(p_value)]
mef_hi  <- res[query == "MEF-unique"  & grepl(fib_pat, term_name, ignore.case = TRUE)][order(p_value)]
fwrite(mesc_hi, file.path(OUT_DIR, "gprofiler2_mESCunique_selfRenewal.csv"))
fwrite(mef_hi,  file.path(OUT_DIR, "gprofiler2_MEFunique_fibroblast.csv"))

# ============================================================
# 4. Two standalone dotplots (F16 = mESC, F17 = MEF)
# ============================================================
SRC_KEEP <- c("GO:BP","GO:MF","GO:CC","KEGG","REAC","WP","TF")
src_cols <- c("GO:BP"="#1b9e77","GO:MF"="#d95f02","GO:CC"="#7570b3","KEGG"="#e7298a",
              "REAC"="#66a61e","WP"="#e6ab02","TF"="#a6761d")
dotplot_std <- function(d, title, subtitle, n_top = 14) {
  d <- d[source %in% SRC_KEEP][order(p_value)]
  d <- head(d[!duplicated(term_name)], n_top)
  d[, neglog10p := -log10(p_value)]
  d[, source := factor(source, levels = SRC_KEEP)]
  d[, term_name := factor(term_name, levels = d[order(neglog10p)]$term_name)]
  ggplot(d, aes(neglog10p, term_name, size = intersection_size, color = source)) +
    geom_point() +
    scale_color_manual(values = src_cols, name = "Source", drop = TRUE) +
    scale_size_continuous(range = c(3, 9), name = "Genes") +
    labs(x = expression(-log[10]~"(g:SCS adjusted p)"), y = NULL, title = title, subtitle = subtitle) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 10, color = "grey30"),
          axis.text.y = element_text(size = 9.5))
}
save_fig <- function(p, id) {
  ggsave(file.path(OUT_DIR, paste0(id, ".pdf")), p, width = 10, height = 6.5)
}
save_fig(dotplot_std(mesc_hi, "mESC-unique G4 peak genes - g:Profiler enrichment",
                     sprintf("%d nearest genes | whole-genome background | self-renewal / pluripotency / signalling", length(mesc_genes))),
         "F16_mESCunique_gProfiler_enrichment")
save_fig(dotplot_std(mef_hi, "MEF-unique G4 peak genes - g:Profiler enrichment",
                     sprintf("%d nearest genes | whole-genome background | fibroblast / cytoskeleton / ECM", length(mef_genes))),
         "F17_MEFunique_gProfiler_enrichment")

cat("Done. Outputs in", OUT_DIR, "\n")
