# =============================================================
# F25: Enrichr GO Biological Process heatmap for the scG4 G4 peak
# differential categories (unchanged / MEF / ESC).
#
# Analogous to fig1 panel G: top-10 enrichR GO BP terms per set,
# heatmap of -log10 P.value, same colour scale as panel G.
#
# Pipeline:
#   1. Reads the per-peak differential categories exported by F24
#      (rev/outputs/RNAseq_lfc_violin/F24_peak_categories.csv):
#      same peak set, classification, and "low" exclusion as F24.
#   2. Per category (unchanged / MEF / ESC) the unique genes
#      annotated to category peaks (refGene TSS +/- 1 kb) form the
#      enrichR gene set.
#   3. enrichR "GO_Biological_Process_2023" per category (queried
#      separately); significant terms = P.value < 0.05 (same
#      convention as the fig1F enrichR run).
#   4. Heatmap of top-10 significant terms per category (by
#      P.value), -log10 P.value.
#
# Output (rev/outputs/enrichR_GO_BP/):
#     panel_F25_enrichR_heatmap.pdf / .png
#     GO_Biological_Process_2023_significant.tsv
#     GO_Biological_Process_2023_top10.tsv
#     gene_sets.tsv
# =============================================================


source("rev/paths.R")
suppressPackageStartupMessages({
  library(data.table)
  library(enrichR)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ----- Paths -----
OUT_DIR <- file.path(OUT_ROOT, "enrichR_GO_BP")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

F24_CATS <- file.path(OUT_ROOT, "RNAseq_lfc_violin", "F24_peak_categories.csv")
stopifnot(file.exists(F24_CATS))

ENRICHR_DB <- "GO_Biological_Process_2023"
N_TOP <- 10L

# ----- 1. F24 differential categories -----
cat("Reading F24 peak categories...\n")
df <- fread(F24_CATS)
required <- c("diff", "gene", "gn", "cl0", "cl1")
stopifnot(all(required %in% names(df)))
df <- df[diff != "low", ]
df$diff <- factor(df$diff, levels = c("unchanged", "MEF", "ESC"))
cat(sprintf("  Peaks per category: %s\n",
            paste(sprintf("%s=%s", names(table(df$diff)), table(df$diff)),
                  collapse = ", ")))

# ----- 2. Gene sets per category -----
cat("Building gene sets...\n")
gene_sets <- df[!is.na(gn) & gn != "", .(set = diff, gene = gn)][!duplicated(paste(set, gene))]
fwrite(gene_sets, file.path(OUT_DIR, "gene_sets.tsv"), sep = "\t")
cat("  Gene set sizes (unique genes):\n")
print(gene_sets[, .(n_genes = uniqueN(gene)), by = set])

# ----- 3. enrichR (one query per category) -----
cat("Running enrichR (", ENRICHR_DB, ")...\n", sep = "")
setEnrichrSite("Enrichr")
res <- lapply(levels(df$diff), function(s) {
  genes <- gene_sets[set == s, unique(gene)]
  cat(sprintf("  Querying %s (%d genes)...\n", s, length(genes)))
  out <- enrichr(genes, databases = ENRICHR_DB)[[1]]
  out$set <- s
  out
})
terms <- as.data.table(do.call(rbind, res))
required <- c("set", "Term", "P.value", "Adjusted.P.value", "Combined.Score", "Genes")
stopifnot(all(required %in% names(terms)))
terms <- terms[, ..required]
fwrite(terms, file.path(OUT_DIR, paste0(ENRICHR_DB, ".tsv")), sep = "\t")

sig <- terms[P.value < 0.05]
cat(sprintf("  Significant terms (P < 0.05): %s\n",
            paste(sprintf("%s=%d", names(table(sig$set)), table(sig$set)),
                  collapse = ", ")))
if (!nrow(sig)) {
  cat("  No significant terms; using top-10 by P.value per set\n")
  sig <- terms[order(P.value), head(.SD, N_TOP), by = set]
}
fwrite(sig, file.path(OUT_DIR, paste0(ENRICHR_DB, "_significant.tsv")), sep = "\t")

selected <- sig[order(P.value), head(.SD, N_TOP), by = set]
fwrite(selected, file.path(OUT_DIR, paste0(ENRICHR_DB, "_top10.tsv")), sep = "\t")
cat("  Top-10 terms per set:\n")
print(selected[, .(set, Term, P.value = signif(P.value, 3), Adjusted.P.value = signif(Adjusted.P.value, 3))])

# ----- 6. Heatmap (analogous to fig1 panel G) -----
cat("Plotting heatmap...\n")
sets <- c("unchanged", "MEF", "ESC")
wide <- dcast(selected, Term ~ set, value.var = "P.value", fill = 1)
for (set_name in sets) {
  if (!set_name %in% names(wide)) wide[, (set_name) := 1]
}
mat <- as.matrix(wide[, ..sets])
rownames(mat) <- wide$Term
logp <- -log10(pmax(mat, .Machine$double.xmin))
term_names <- sub("\\s*\\(GO:[0-9]+\\)", "", rownames(logp))
dimnames(logp) <- list(term_names, colnames(logp))

col_fun <- colorRamp2(c(0, 1, 2, 4, 8),
                      c("grey95", "#fee0d2", "#fc9272", "#de2d26", "#99000d"))
hm <- Heatmap(
  logp,
  name = "-log10 p",
  col = col_fun,
  column_title = "Enrichr GO BP terms: scG4 G4 peak differential categories",
  cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE,
  rect_gp = gpar(col = "black", lwd = 0.1),
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  row_names_max_width = unit(8, "cm"),
  column_names_gp = gpar(fontsize = 9), column_names_rot = 45,
  heatmap_width = unit(15, "cm"),
  heatmap_height = unit(max(5, nrow(logp) * 0.6 + 1), "cm")
)

pdf(file.path(OUT_DIR, "panel_F25_enrichR_heatmap.pdf"),
    width = 18, height = max(5, nrow(logp) * 0.25 + 1))
draw(hm)
dev.off()

png(file.path(OUT_DIR, "panel_F25_enrichR_heatmap.png"),
    width = 18, height = max(5, nrow(logp) * 0.25 + 1), units = "in", res = 300)
draw(hm)
dev.off()

cat("\nWrote outputs to", OUT_DIR, "\n")
