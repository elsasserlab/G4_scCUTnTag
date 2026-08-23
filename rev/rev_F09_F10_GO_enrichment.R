#!/usr/bin/env Rscript
# =============================================================
# Cell-specific G4 peaks -> GO Biological Process enrichment (F09, F10).
#
# Self-contained: reads ONLY the original input peak files + GTF, does the
# Fig-1D overlap analysis to extract the cluster-discriminating (unique)
# peaks, annotates them to nearest genes, and runs GO BP enrichment.
# No Venn plot is produced (removed by request); the overlap analysis and
# unique-peak extraction are kept.
#
# Inputs  (input_files/):
#   cluster_spec_peaks/0_peaks.narrowPeak   (scG4 cluster 0 = MEF/3T3 context)
#   cluster_spec_peaks/1_peaks.narrowPeak   (scG4 cluster 1 = mESC context)
#   bulkG4CnT_mESC_rep1.broadPeak          (bulk mESC G4)
#   bulkG4CnT_3T3_rep1.broadPeak           (bulk MEF/3T3 G4)
#   mm10_.annotation.gtf.gz                (gene models for nearest-TSS)
#
# Outputs (rev/outputs/cell_specific_peaks_functional/):
#   bed/mESC_specific_peaks.bed , bed/MEF_specific_peaks.bed   (unique peaks)
#   mESC_specific_nearest_genes.tsv , MEF_specific_nearest_genes.tsv
#   GO_mESC_BP.csv , GO_MEF_BP.csv
#   GO_BP_mESC_dotplot.{pdf,png} (F10) , GO_BP_MEF_dotplot.{pdf,png} (F09)
# =============================================================


source("rev/paths.R")
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
})

# ----- Paths -----
inp <- function(fn) {
  hits <- list.files(DATA, pattern = fn, recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("File not found in data/: ", fn)
  hits[1]
}

GTF_PATH <- file.path(GENOME_DIR, "mm10_.annotation.gtf.gz")
OUT_DIR  <- file.path(OUT_ROOT, "cell_specific_peaks_functional")
BED_DIR  <- file.path(OUT_DIR, "bed")
dir.create(BED_DIR, recursive = TRUE, showWarnings = FALSE)

CANONICAL <- c(paste0("chr", 1:19), "chrX", "chrY")

# ============================================================
# 1. Load the original peak files
# ============================================================
read_peak_gr <- function(path, label) {
  dt <- fread(path)
  GRanges(seqnames = dt$V1,
          ranges   = IRanges(start = dt$V2, end = dt$V3,
                             names = rep(label, nrow(dt))))
}
cat("Loading original peak files...\n")
cluster0 <- read_peak_gr(CLUSTER_PEAKS_0, "cluster0")
cluster1 <- read_peak_gr(CLUSTER_PEAKS_1, "cluster1")
mesc     <- read_peak_gr(file.path(BULK_PEAKS, "GSM8836082_bulkG4CnT_mESC_rep1.broadPeak"),        "mESC")
mef      <- read_peak_gr(file.path(BULK_PEAKS, "GSM8836084_bulkG4CnT_3T3_rep1.broadPeak"),         "MEF")

keep_canonical <- function(gr) gr[as.character(seqnames(gr)) %in% CANONICAL]
cluster0 <- keep_canonical(cluster0); cluster1 <- keep_canonical(cluster1)
mesc     <- keep_canonical(mesc);     mef      <- keep_canonical(mef)

# ============================================================
# 2. Overlap analysis -> cluster-discriminating unique peaks (Fig 1D)
#    mESC-specific = bulk mESC peaks overlapping cluster 1 but NOT cluster 0
#    MEF-specific  = bulk MEF  peaks overlapping cluster 0 but NOT cluster 1
# ============================================================
cat("\nExtracting cluster-specific (unique) peaks...\n")
mesc_specific <- mesc[ overlapsAny(mesc, cluster1, ignore.strand = TRUE) &
                      !overlapsAny(mesc, cluster0, ignore.strand = TRUE)]
mef_specific  <- mef [ overlapsAny(mef,  cluster0, ignore.strand = TRUE) &
                      !overlapsAny(mef,  cluster1, ignore.strand = TRUE)]
cat(sprintf("  mESC-specific peaks: %s  (Fig 1D reports ~6,543)\n",
            format(length(mesc_specific), big.mark = ",")))
cat(sprintf("  MEF-specific  peaks: %s  (Fig 1D reports ~3,321)\n",
            format(length(mef_specific),  big.mark = ",")))
export(mesc_specific, file.path(BED_DIR, "mESC_specific_peaks.bed"), format = "BED")
export(mef_specific,  file.path(BED_DIR, "MEF_specific_peaks.bed"),  format = "BED")

# ============================================================
# 3. Per-gene TSS from the GTF (one primary TSS per gene)
# ============================================================
cat("\nBuilding per-gene TSS from GTF...\n")
gtf <- rtracklayer::import(GTF_PATH)
gtf <- gtf[as.character(seqnames(gtf)) %in% CANONICAL]
seqlevels(gtf) <- intersect(seqlevels(gtf), CANONICAL)
tx  <- gtf[gtf$type == "transcript"]
tx_dt <- data.table(chrom = as.character(seqnames(tx)),
                    strand = as.character(strand(tx)),
                    start = start(tx), end = end(tx),
                    gene_id = tx$gene_id, gene_name = tx$gene_name)
tx_dt[, tss := ifelse(strand == "+", start, end)]
gene_tss <- tx_dt[, .(
  chrom     = data.table::first(chrom),
  strand    = data.table::first(strand),
  tss       = ifelse(data.table::first(strand) == "+", min(tss), max(tss)),
  gene_name = data.table::first(gene_name)
), by = gene_id]
gene_tss <- gene_tss[!is.na(gene_name) & nzchar(gene_name)]
tss_gr <- GRanges(seqnames = gene_tss$chrom,
                  ranges   = IRanges(start = gene_tss$tss, end = gene_tss$tss),
                  strand   = gene_tss$strand,
                  gene_id  = gene_tss$gene_id, gene_name = gene_tss$gene_name)

# ============================================================
# 4. Annotate each peak to its nearest TSS -> unique gene lists
# ============================================================
cat("\nAnnotating peaks to nearest TSS...\n")
annotate_nearest <- function(query_gr, feat_gr, label) {
  if (length(query_gr) == 0) return(character(0))
  dn  <- distanceToNearest(query_gr, feat_gr, ignore.strand = TRUE)
  out <- rep(NA_character_, length(query_gr))
  out[queryHits(dn)] <- as.character(feat_gr$gene_name[subjectHits(dn)])
  cat(sprintf("  %s: %s peaks -> %s gene assignments\n", label,
              format(length(query_gr), big.mark = ","),
              format(sum(!is.na(out) & nzchar(out)), big.mark = ",")))
  out
}
mesc_genes <- unique(annotate_nearest(mesc_specific, tss_gr, "mESC-specific"))
mesc_genes <- mesc_genes[!is.na(mesc_genes) & nzchar(mesc_genes)]
mef_genes  <- unique(annotate_nearest(mef_specific,  tss_gr, "MEF-specific"))
mef_genes  <- mef_genes[!is.na(mef_genes)  & nzchar(mef_genes)]
fwrite(data.table(gene = mesc_genes), file.path(OUT_DIR, "mESC_specific_nearest_genes.tsv"), col.names = FALSE)
fwrite(data.table(gene = mef_genes),  file.path(OUT_DIR, "MEF_specific_nearest_genes.tsv"),  col.names = FALSE)
cat(sprintf("  unique genes: %d (mESC-specific)  %d (MEF-specific)\n",
            length(mesc_genes), length(mef_genes)))

# ============================================================
# 5. GO Biological Process enrichment (clusterProfiler)
# ============================================================
cat("\nRunning GO BP enrichment...\n")
to_entrez <- function(symbols) {
  m <- suppressMessages(bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
  unique(m$ENTREZID)
}
mesc_entrez     <- to_entrez(mesc_genes)
mef_entrez      <- to_entrez(mef_genes)
universe_entrez <- unique(to_entrez(unique(gene_tss$gene_name)))
cat(sprintf("  ENTREZ: %d (mESC)  %d (MEF)  | universe: %d\n",
            length(mesc_entrez), length(mef_entrez), length(universe_entrez)))

run_go_bp <- function(entrez) {
  res <- enrichGO(gene = entrez, universe = universe_entrez, OrgDb = org.Mm.eg.db,
                  keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  minGSSize = 10, maxGSSize = 500, readable = TRUE)
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  simplify(res, cutoff = 0.7, by = "p.adjust", select_fun = min)
}
go_mesc <- run_go_bp(mesc_entrez)
go_mef  <- run_go_bp(mef_entrez)
if (!is.null(go_mesc)) fwrite(as.data.table(as.data.frame(go_mesc)), file.path(OUT_DIR, "GO_mESC_BP.csv"))
if (!is.null(go_mef))  fwrite(as.data.table(as.data.frame(go_mef)),  file.path(OUT_DIR, "GO_MEF_BP.csv"))

# ============================================================
# 6. GO BP dot plots (F09 = MEF, F10 = mESC)
# ============================================================
cat("\nDrawing GO BP dot plots...\n")
plot_go <- function(res, title, n_top = 15) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) { cat("  (no terms for", title, ")\n"); return(NULL) }
  dotplot(res, showCategory = n_top, font.size = 9) +
    ggtitle(title) + theme(plot.title = element_text(face = "bold"))
}
p_mesc <- plot_go(go_mesc, "mESC-specific peaks — GO Biological Process")
p_mef  <- plot_go(go_mef,  "MEF-specific peaks — GO Biological Process")
if (!is.null(p_mesc)) {
  ggsave(file.path(OUT_DIR, "GO_BP_mESC_dotplot.pdf"), p_mesc, width = 9, height = 7)
}
if (!is.null(p_mef)) {
  ggsave(file.path(OUT_DIR, "GO_BP_MEF_dotplot.pdf"),  p_mef,  width = 9, height = 7)
}
cat("\nDone. Outputs in:", OUT_DIR, "\n")
