#!/usr/bin/env Rscript
# ===========================================================================
# Figure 3: scBridge annotation, astrocyte G4, Cicero coaccessibility
#
# Panels:
#   A - Coembedded UMAPs: cell types, modality, predictions, reliability
#   B - AST vs non-AST: 2x2 grid (UMAP, full-gene GA volcano, differential
#       G4 peaks at promoters, scRNA-seq expression volcano)
#   C - Coverage at astrocyte marker genes
#   D - Feature plots: Tnik, Pitpnc1, Pbx1, Nwd1 (normalized RNA, viridis)
#   E - Cicero browser tracks
#   S1 - Supplementary: AST-specific G4 peaks intersected with ENCODE4 cCRE
#       classes (data/genome/cCRE.mm10.bed); horizontal stacked bar with
#       a "no cCRE overlap" category (primary class = max bp overlap).
#
# Run from repo root:
#   Rscript fig3/fig3.R
#   Rscript fig3/fig3.R --only=A,B
#   Rscript fig3/fig3.R --skip=E
#   Rscript fig3/fig3.R --only=S1
# ===========================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(ggrepel)
  library(ggrastr)
  library(viridis)
  library(ggpubr)
  library(RColorBrewer)
  library(stringr)
  library(EnsDb.Mmusculus.v79)
  library(GenomicRanges)
  library(rtracklayer)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
ROOT       <- dirname(SCRIPT_DIR)
RESULTS    <- file.path(SCRIPT_DIR, "results")
OUT        <- file.path(SCRIPT_DIR, "outputs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SC_DIR  <- file.path(RESULTS, "scBridge")
INT_DIR <- file.path(ROOT, "fig2", "results", "integration", "outputs")

save_fig <- function(p, name, width = 10, height = 8) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = width, height = height, device = "pdf")
}

args <- commandArgs(trailingOnly = TRUE)
only <- NULL
skip <- NULL
for (i in seq_along(args)) {
  a <- args[i]
  if (startsWith(a, "--only=")) only <- strsplit(sub("--only=", "", a), ",")[[1]]
  if (startsWith(a, "--skip=")) skip <- strsplit(sub("--skip=", "", a), ",")[[1]]
  if (a == "--only" && !is.na(args[i + 1])) only <- strsplit(args[i + 1], ",")[[1]]
  if (a == "--skip" && !is.na(args[i + 1])) skip <- strsplit(args[i + 1], ",")[[1]]
}
run_panel <- function(p) {
  if (!is.null(skip) && p %in% skip) return(FALSE)
  if (!is.null(only) && !(p %in% only)) return(FALSE)
  TRUE
}

pastel1 <- brewer.pal(8, "Pastel1")
CT_COLS <- c(
  'AST' = pastel1[1], 'COP-NFOL' = pastel1[2], 'MOL' = pastel1[3],
  'OPC' = pastel1[4], 'OEC' = pastel1[5], 'VEC' = pastel1[6],
  'VLMC' = pastel1[7], 'Pericytes' = pastel1[8]
)
UMAP_THEME <- theme(
  plot.title = element_text(size = 16, face = "bold"),
  axis.text = element_text(size = 11, color = "black"),
  axis.title = element_text(size = 13),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 11)
)

# ===========================================================================
# Load shared data
# ===========================================================================
if (any(vapply(LETTERS[1:4], run_panel, logical(1)))) {
  cat("== Loading data ==\n")

  cat("  Coembedded UMAP...\n")
  umap_df <- fread(file.path(SC_DIR, "umap_coembedded.csv"))
  umap_df[, modality := ifelse(
    Domain == "Bartosovic_scRNA-Seq", "scRNA-Seq", "G4 scCUT&Tag"
  )]

  cat("  scBridge predictions + reliability...\n")
  pred <- fread(file.path(SC_DIR, "scbridge_predictions.csv"), header = TRUE)
  rel  <- fread(file.path(SC_DIR, "scbridge_reliability.csv"), header = TRUE)
  setnames(pred, "V1", "barcode")
  setnames(rel, "V1", "barcode")
  pred <- merge(pred, rel, by = "barcode")
  pred[, Prediction_clean := fifelse(
    str_detect(Prediction, "Novel"),
    "unreliable",
    str_replace_all(Prediction, "Astrocytes", "AST")
  )]
  pred[, Prediction_clean := str_replace_all(Prediction_clean, "Oligodendrocytes", "MOL")]

  cat("  Seurat objects...\n")
  rna <- readRDS(file.path(INT_DIR, "scRNA_Seq_Seurat_object.Rds"))
  rna@meta.data$cell_type <- str_replace_all(rna@meta.data$cell_type, "Astrocytes", "AST")
  rna@meta.data$cell_type <- str_replace_all(rna@meta.data$cell_type, "Oligodendrocytes", "MOL")

  g4 <- readRDS(file.path(INT_DIR, "G4_scRNA_integration.Rds"))
}

# ===========================================================================
# Panel A: Coembedded UMAPs
# ===========================================================================
if (run_panel("A")) {
  cat("== Panel A: Coembedded UMAPs ==\n")

  # A1: scRNA-seq reference cell types (scRNA-seq cells only, using coembedded UMAP)
  rna_barcodes <- umap_df[Domain == "Bartosovic_scRNA-Seq", barcode]
  rna_umap <- umap_df[Domain == "Bartosovic_scRNA-Seq"]
  celltype_vec <- setNames(rna@meta.data$cell_type, rownames(rna@meta.data))
  rna_umap[, cell_type := celltype_vec[barcode]]
  rna_umap[, cell_type := factor(cell_type, levels = rev(names(CT_COLS)))]

  p_A1 <- ggplot(rna_umap, aes(UMAP1, UMAP2, color = cell_type)) +
    geom_point_rast(size = 0.8) +
    scale_color_manual(values = CT_COLS, drop = FALSE) +
    labs(title = "Cell type (scRNA-Seq)", color = "Cell type") +
    theme_classic() + UMAP_THEME +
    guides(color = guide_legend(override.aes = list(size = 3)))

  # A2: Modality
  p_A2 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = modality)) +
    geom_point_rast(size = 0.8) +
    scale_color_manual(values = c("scRNA-Seq" = "#a6bddb", "G4 scCUT&Tag" = "#fc9272")) +
    labs(title = "Modality", color = "Modality") +
    theme_classic() + UMAP_THEME +
    guides(color = guide_legend(override.aes = list(size = 3)))

  # A3: scBridge predictions (G4 cells only)
  g4_umap <- umap_df[Domain == "scCutTag_gene_activity_scores"]
  pred_match <- pred[match(g4_umap$barcode, pred$barcode)]
  g4_umap[, prediction := pred_match$Prediction_clean]
  pred_levels <- c(names(CT_COLS), "unreliable")
  pred_cols <- c(CT_COLS, "unreliable" = "#f0f0f0")
  g4_umap[, prediction := factor(prediction, levels = rev(pred_levels))]

  p_A3 <- ggplot(g4_umap, aes(UMAP1, UMAP2, color = prediction)) +
    geom_point_rast(size = 1.0) +
    scale_color_manual(values = pred_cols, drop = FALSE) +
    labs(title = "Prediction", color = "Cell type") +
    theme_classic() + UMAP_THEME +
    guides(color = guide_legend(override.aes = list(size = 3)))

  # A4: Reliability (G4 cells only)
  g4_umap[, reliability := pred_match$Reliability]

  p_A4 <- ggplot(g4_umap, aes(UMAP1, UMAP2, color = reliability)) +
    geom_point_rast(size = 1.0) +
    scale_color_viridis(limits = c(0, 1), oob = scales::squish) +
    labs(title = "Reliability", color = "Reliability") +
    theme_classic() + UMAP_THEME

  p_A <- (p_A1 | p_A2) / (p_A3 | p_A4)
  save_fig(p_A, "panel_A_coembedded_UMAPs", 16, 16)
  cat("  Saved panel_A_coembedded_UMAPs.pdf\n")
}

# ===========================================================================
# Panel B: AST vs non-AST UMAP + volcano
# ===========================================================================
if (run_panel("B")) {
  cat("== Panel B: AST vs non-AST ==\n")

  # B1: UMAP colored by AST status
  g4_umap_b <- umap_df[Domain == "scCutTag_gene_activity_scores"]
  pred_match_b <- pred[match(g4_umap_b$barcode, pred$barcode)]
  g4_umap_b[, AST_status := fifelse(pred_match_b$Prediction_clean == "AST", "AST", "non-AST")]
  g4_umap_b[, AST_status := factor(AST_status, levels = c("AST", "non-AST"))]

  pastel1_9 <- brewer.pal(9, "Pastel1")
  ast_cols <- c("AST" = pastel1_9[1], "non-AST" = pastel1_9[9])
  p_B1 <- ggplot(g4_umap_b, aes(UMAP1, UMAP2, color = AST_status)) +
    geom_point_rast(size = 0.5) +
    scale_color_manual(values = ast_cols) +
    labs(title = "", color = "") +
    theme_classic(base_size = 16) + UMAP_THEME +
    theme(
      axis.text = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 17)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))

  # B2/B3: Volcano plots mirroring the original figure legend ("differential
  # G-quadruplex region analysis ... gene activity score was considered").
  # B2 uses the full-gene GA assay; B3 recomputes GA restricted to promoter
  # windows (TSS-2kb..TSS) for a promoter-specific comparison.
  cat("  Running FindMarkers AST vs non-AST (LR test on full-gene GA)...\n")
  barcodes_ast <- pred[Prediction_clean == "AST", barcode]
  barcodes_ast <- barcodes_ast[barcodes_ast %in% rownames(g4@meta.data)]
  barcodes_nonast <- pred[Prediction_clean != "AST" & Prediction_clean != "unreliable", barcode]
  barcodes_nonast <- barcodes_nonast[barcodes_nonast %in% rownames(g4@meta.data)]

  FOCUS_GENES <- c("Tnik", "Pitpnc1", "Pbx1", "Nwd1")
  volc_cols <- c("up in AST" = pastel1[1], "down in AST" = brewer.pal(9, "Set3")[9], "unaltered" = "grey80")
  volc_sizes <- c("up in AST" = 5, "down in AST" = 5, "unaltered" = 2)
  volc_alphas <- c("up in AST" = 1, "down in AST" = 1, "unaltered" = 0.5)

  prepare_volcano <- function(df) {
    df <- as.data.table(df)
    if (!"gene" %in% names(df)) stop("prepare_volcano: 'gene' column required")
    df[, group := fifelse(
      avg_log2FC > 0.5 & p_val_adj < 0.05, "up in AST",
      fifelse(avg_log2FC < -0.5 & p_val_adj < 0.05, "down in AST", "unaltered")
    )]
    df[, sign_label := fifelse(
      group != "unaltered" & p_val_adj < 0.05 & abs(avg_log2FC) > 0.5, gene, ""
    )]
    df[gene %in% FOCUS_GENES, sign_label := gene]
    df
  }

  make_volcano <- function(df, title_str) {
    rng <- max(6, ceiling(max(abs(df$avg_log2FC))))
    ggplot(df, aes(avg_log2FC, -log10(p_val_adj),
                   fill = group, size = group, alpha = group)) +
      geom_point_rast(shape = 21, colour = "black") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
      scale_fill_manual(values = volc_cols) +
      scale_size_manual(values = volc_sizes) +
      scale_alpha_manual(values = volc_alphas) +
      scale_x_continuous(breaks = seq(-10, 10, 2), limits = c(-rng, rng)) +
      labs(x = "log2 Fold Change", y = "-log10 adj. p-value",
           title = title_str, fill = " ") +
      guides(alpha = "none", size = "none",
             fill = guide_legend(override.aes = list(size = 7))) +
      theme_minimal(base_size = 18) +
      theme(
        plot.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 17)
      ) +
      geom_text_repel(aes(label = sign_label), size = 5, max.overlaps = 100,
                      segment.size = 0.35, min.segment.length = 0)
  }

  DefaultAssay(g4) <- "GA"
  diff_AST <- FindMarkers(
    g4, ident.1 = barcodes_ast, ident.2 = barcodes_nonast,
    only.pos = FALSE, assay = "GA", logfc.threshold = 0,
    test.use = "LR", latent.vars = "peak_region_fragments"
  )
  diff_GA_full <- data.table(
    gene = rownames(diff_AST), as.data.table(diff_AST), check.names = FALSE
  )
  fwrite(diff_GA_full, file.path(OUT, "panel_B_diff_GA_full_AST_vs_nonAST.csv"))
  volc_ga <- prepare_volcano(diff_GA_full)
  p_B2 <- make_volcano(volc_ga, "G4 gene activity (full gene)")

  # B3: differential G4 peaks at promoters. All CellRanger peaks are tested
  # (full peak set = proper background for BH correction), then results are
  # subset to peaks within +/-2kb of a TSS for display and gene lists.
  # Wilcoxon on log-normalized values (LR+depth-covariate loses power on
  # sparse single-cell G4 counts).
  cat("  Running FindMarkers on all G4 peaks (wilcox, log-normalized)...\n")
  DefaultAssay(g4) <- "peaks"
  g4 <- NormalizeData(g4, assay = "peaks",
                      normalization.method = "LogNormalize", verbose = FALSE)
  diff_peaks <- FindMarkers(
    g4, ident.1 = barcodes_ast, ident.2 = barcodes_nonast,
    only.pos = FALSE, assay = "peaks", logfc.threshold = 0,
    test.use = "wilcox"
  )
  edb_genes <- suppressWarnings(GenomicFeatures::genes(EnsDb.Mmusculus.v79))
  tss_pos <- ifelse(as.character(strand(edb_genes)) == "-",
                    end(edb_genes), start(edb_genes))
  peaks_gr <- g4[["peaks"]]@ranges
  tss_gr <- GenomicRanges::GRanges(
    paste0("chr", as.character(seqnames(edb_genes))),
    IRanges::IRanges(tss_pos, width = 1), strand = "*"
  )
  names(tss_gr) <- edb_genes$gene_name
  tss_gr <- tss_gr[as.character(seqnames(tss_gr)) %in%
                     paste0("chr", c(1:19, "X", "Y"))]
  # nearest() returns NA for peaks on contigs without TSS annotation
  # (chrM, scaffolds) - handle explicitly.
  nn <- GenomicRanges::nearest(peaks_gr, tss_gr, ignore.strand = TRUE)
  pk_gene <- rep(NA_character_, length(peaks_gr))
  pk_dist <- rep(NA_integer_, length(peaks_gr))
  has_nn <- !is.na(nn)
  pk_gene[has_nn] <- names(tss_gr)[nn[has_nn]]
  pk_dist[has_nn] <- as.integer(GenomicRanges::distance(
    peaks_gr[has_nn], tss_gr[nn[has_nn]], ignore.strand = TRUE
  ))
  pk_key <- paste(as.character(seqnames(peaks_gr)), start(peaks_gr),
                  end(peaks_gr), sep = "-")
  sel_row <- match(rownames(diff_peaks), pk_key)
  stopifnot(!any(is.na(sel_row)))
  diff_peaks_out <- data.table(
    peak = rownames(diff_peaks),
    gene = pk_gene[sel_row],
    distance_to_TSS = pk_dist[sel_row],
    as.data.table(diff_peaks), check.names = FALSE
  )
  fwrite(diff_peaks_out,
         file.path(OUT, "panel_B_diff_G4peaks_AST_vs_nonAST.csv"))
  sel_pkp <- diff_peaks_out[distance_to_TSS <= 2000]
  fwrite(sel_pkp,
         file.path(OUT, "panel_B_diff_G4peaks_promoters_AST_vs_nonAST.csv"))
  cat("  Peaks tested:", nrow(diff_peaks_out), "| promoter-proximal:",
      nrow(sel_pkp), "| up-hits at promoters:",
      sel_pkp[avg_log2FC > 0.5 & p_val_adj < 0.05, .N], "\n")
  p_B3 <- make_volcano(prepare_volcano(copy(sel_pkp)),
                       "G4 peaks at promoters")

  # Genes that are AST-upregulated hits in EITHER differential (full-gene GA
  # or promoter-proximal peaks): used to label the scRNA-seq volcano.
  up_labeled <- union(
    volc_ga[group == "up in AST" & sign_label != "", gene],
    sel_pkp[avg_log2FC > 0.5 & p_val_adj < 0.05, gene]
  )
  cat("  RNA label set:", length(up_labeled), "genes (GA-full or peak-promoter hits)\n")

  # B4: scRNA-seq DGE. The RNA data is the co-analyzed reference dataset (no
  # shared cells with the G4 experiment); groups are the reference's own
  # Astrocytes vs all other reference cells. Wilcoxon on log-normalized
  # values (the bundled FindAllMarkers table is pre-filtered to ~2.5k
  # genes/cluster and would hide most G4-hit labels).
  cat("  Running FindMarkers on scRNA-seq (wilcox, log-normalized):",
      sum(rna@meta.data$cell_type == "AST"), "AST vs",
      sum(rna@meta.data$cell_type != "AST"), "non-AST cells...\n")
  rna_cells_ast <- colnames(rna)[rna@meta.data$cell_type == "AST"]
  rna_cells_nonast <- setdiff(colnames(rna), rna_cells_ast)
  DefaultAssay(rna) <- "RNA"
  rna <- NormalizeData(rna, assay = "RNA",
                       normalization.method = "LogNormalize", verbose = FALSE)
  diff_RNA <- FindMarkers(
    rna, ident.1 = rna_cells_ast, ident.2 = rna_cells_nonast,
    only.pos = FALSE, assay = "RNA", logfc.threshold = 0,
    test.use = "wilcox"
  )
  diff_RNA_out <- data.table(
    gene = rownames(diff_RNA), as.data.table(diff_RNA), check.names = FALSE
  )
  fwrite(diff_RNA_out, file.path(OUT, "panel_B_diff_RNA_AST_vs_nonAST.csv"))
  volc_rna <- prepare_volcano(diff_RNA_out)
  volc_rna[, sign_label := fifelse(gene %in% up_labeled, gene, "")]
  cat("  RNA labels drawn:", sum(volc_rna$sign_label != ""), "genes\n")
  # Sparse display: every gene as a minimal dot; only the G4-hit label set as
  # prominent circles (light red = up in AST, grey = anything else).
  # p_val_adj can underflow to exactly 0 -> cap -log10 at y_cap.
  y_cap <- 250
  volc_rna[, y_show := pmin(-log10(pmax(p_val_adj, .Machine$double.xmin)), y_cap)]
  volc_rna_fg <- copy(volc_rna[sign_label != ""])
  volc_rna_fg[, fg_class := fifelse(group == "up in AST", "up in AST", "other")]
  fg_cols <- c("up in AST" = pastel1[1], "other" = "grey75")
  rng_rna <- max(6, ceiling(max(abs(volc_rna$avg_log2FC))))
  p_B4 <- ggplot(volc_rna, aes(avg_log2FC, y_show)) +
    geom_point_rast(size = 0.15, colour = "grey65", alpha = 0.5) +
    geom_point(data = volc_rna_fg,
               aes(avg_log2FC, y_show, fill = fg_class),
               shape = 21, size = 5, colour = "black", stroke = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = fg_cols, name = "AST-upregulated\nG4 hits") +
    scale_x_continuous(breaks = seq(-10, 10, 2),
                       limits = c(-rng_rna, rng_rna)) +
    labs(x = "log2 Fold Change", y = "-log10 adj. p-value",
         title = "scRNA-seq expression") +
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 17)
    ) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    geom_text_repel(data = volc_rna_fg, aes(label = sign_label), size = 5,
                    max.overlaps = 100, segment.size = 0.35,
                    min.segment.length = 0)

  p_B <- (p_B1 | p_B2) / (p_B3 | p_B4)
  save_fig(p_B, "panel_B_AST_vs_nonAST", 22, 18)
  cat("  Saved panel_B_AST_vs_nonAST.pdf\n")
}

# ===========================================================================
# Panel C: Coverage at differential-G4-promoter-peak genes
# ===========================================================================
if (run_panel("C")) {
  cat("== Panel C: Coverage at differential-G4-promoter-peak genes ==\n")

  # Export group-split CellRanger peaks as BED files (for genome-browser viewing).
  # Cell groups mirror the corrected RPGC bigwig tracks: predicted AST = 183
  # scBridge-Astrocytes, non-AST = authors' non_AST list (2521 cells).
  bc_dir_c <- file.path(ROOT, "data/GSE291468/scBridge_predictions/barcodes")
  ast_cells_c <- intersect(fread(file.path(bc_dir_c, "barcodes_Astrocytes.tsv"),
                                 header = FALSE)$V1,
                           colnames(g4[["peaks"]]@counts))
  nonast_cells_c <- intersect(fread(file.path(bc_dir_c, "barcodes_non_AST.tsv"),
                                    header = FALSE)$V1,
                              colnames(g4[["peaks"]]@counts))
  counts_c <- GetAssayData(g4, assay = "peaks", layer = "counts")
  pct_ast_c <- Matrix::rowMeans(counts_c[, ast_cells_c] > 0)
  pct_non_c <- Matrix::rowMeans(counts_c[, nonast_cells_c] > 0)
  mean_ast_c <- as.numeric(Matrix::rowMeans(counts_c[, ast_cells_c]))
  mean_non_c <- as.numeric(Matrix::rowMeans(counts_c[, nonast_cells_c]))
  log2fc_c <- log2((mean_ast_c + 0.01) / (mean_non_c + 0.01))
  peaks_gr_c <- g4[["peaks"]]@ranges
  n_pk_c <- length(peaks_gr_c)
  bed_rows_c <- function(sel, score) {
    data.frame(
      seqnames = as.character(seqnames(peaks_gr_c))[sel],
      start = start(peaks_gr_c)[sel] - 1L,
      end = end(peaks_gr_c)[sel],
      name = paste0(as.character(seqnames(peaks_gr_c))[sel], ":",
                    start(peaks_gr_c)[sel], "-", end(peaks_gr_c)[sel]),
      score = as.integer(pmin(1000, pmax(0, round(score * 1000)))),
      strand = "*"
    )
  }
  write_bed_c <- function(df, path) {
    fwrite(df, path, sep = "\t", col.names = FALSE, quote = FALSE)
  }
  write_bed_c(bed_rows_c(seq_len(n_pk_c), rep(0.5, n_pk_c)),
              file.path(OUT, "GSM8836086_all_peaks.bed"))
  # group-specific peaks: detected in >=2% of that group's cells and
  # >=4-fold higher detection rate than the other group (|log2 ratio| >= 2)
  log2det_c <- log2((pct_ast_c + 0.001) / (pct_non_c + 0.001))
  sel_a_c <- which(pct_ast_c >= 0.02 & log2det_c >= 2)
  sel_n_c <- which(pct_non_c >= 0.02 & log2det_c <= -2)
  write_bed_c(bed_rows_c(sel_a_c, pct_ast_c[sel_a_c]),
              file.path(OUT, "GSM8836086_AST_specific_peaks.bed"))
  write_bed_c(bed_rows_c(sel_n_c, pct_non_c[sel_n_c]),
              file.path(OUT, "GSM8836086_nonAST_specific_peaks.bed"))
  o_c <- order(-log2fc_c)
  write_bed_c(bed_rows_c(o_c, pmin(1, pmax(0, log2fc_c[o_c] / 4))),
              file.path(OUT, "GSM8836086_AST_vs_nonAST_log2fc.bed"))
  fwrite(
    data.table(
      chr = as.character(seqnames(peaks_gr_c)),
      start = start(peaks_gr_c) - 1L,
      end = end(peaks_gr_c),
      pct_AST = round(pct_ast_c, 3),
      pct_nonAST = round(pct_non_c, 3),
      mean_AST = round(mean_ast_c, 3),
      mean_nonAST = round(mean_non_c, 3),
      log2FC = round(log2fc_c, 3)
    ),
    file.path(OUT, "GSM8836086_AST_vs_nonAST_per_peak.tsv"),
    sep = "\t", quote = FALSE
  )
  cat("  Exported peak BEDs:", length(sel_a_c), "AST-specific,",
      length(sel_n_c), "non-AST-specific ->", OUT, "\n")

  # Panel C fragment-track grouping mirrors the corrected RPGC bigwigs
  # exactly: scBridge Astrocytes (183 cells) vs authors' non-AST list,
  # both restricted to cells present in the Seurat object.
  ast_barcodes <- intersect(fread(file.path(bc_dir_c, "barcodes_Astrocytes.tsv"),
                                  header = FALSE)$V1,
                            rownames(g4@meta.data))
  nonast_barcodes <- intersect(fread(file.path(bc_dir_c, "barcodes_non_AST.tsv"),
                                     header = FALSE)$V1,
                               rownames(g4@meta.data))
  cat("  Panel C groups:", length(ast_barcodes), "AST /",
      length(nonast_barcodes), "non-AST (scBridge)\n")

  g4$AST_status <- ifelse(rownames(g4@meta.data) %in% ast_barcodes,
                          "predicted AST", "predicted non-AST")
  Idents(g4) <- "AST_status"

  # ---- Final Panel C: curated showcase genes --------------------------------
  # Manually selected to showcase different scenarios (AST-up promoter peaks
  # with GA/RNA support, peak-only hits, and negative controls). Windows span
  # 10 kb upstream / 35 kb downstream (strand-aware) of each gene's strongest
  # promoter-proximal peak, falling back to the TSS when no peak exists.
  sel_genes_c <- c("Plcl1", "Nwd1", "Gli2", "Smad9", "Tmem74", "Gm25493",
                   "Scarna17", "Acaa2", "Stambpl1", "App", "Rsg1")
  pk_promo_c <- fread(file.path(OUT, "panel_B_diff_G4peaks_promoters_AST_vs_nonAST.csv"))
  ga_full_c <- fread(file.path(OUT, "panel_B_diff_GA_full_AST_vs_nonAST.csv"))
  rna_dge_c <- fread(file.path(OUT, "panel_B_diff_RNA_AST_vs_nonAST.csv"))
  hit_c <- function(d, dir = 1)
    d[avg_log2FC * dir > 0.5 & p_val_adj < 0.05, unique(gene)]
  pk_up_c <- hit_c(pk_promo_c); ga_up_c <- hit_c(ga_full_c)
  rna_up_c <- hit_c(rna_dge_c)
  ga_dn_c <- hit_c(ga_full_c, -1); rna_dn_c <- hit_c(rna_dge_c, -1)

  DefaultAssay(g4) <- "peaks"
  # CoveragePlot always derives a fragment-based coverage track internally,
  # so the fragment file must be attached. Collapse Idents to a single dummy
  # level so only one minimal built-in lane is drawn next to the bigwigs.
  frag_path_c <- file.path(ROOT, "data/GSE291468/GSM8836086_GFP_sorted_mouse_brain/CellRanger/fragments.tsv.gz")
  frg_c <- Signac::CreateFragmentObject(path = frag_path_c,
                                        cells = colnames(g4))
  slot(g4[["peaks"]], "fragments") <- list(frg_c)
  # Built-in fragment track pooled over all cells (single ident level).
  Idents(g4) <- setNames(rep("all GFP+ cells", ncol(g4)), colnames(g4))

  edb <- EnsDb.Mmusculus.v79
  gdb_c <- genes(edb)

  bw_ast_c <- file.path(ROOT, "data/GSE291468/GSM8836086_Predicted_Astrocytes_RPGC.smooth150.bw")
  bw_non_c <- file.path(ROOT, "data/GSE291468/GSM8836086_Predicted_non-Astrocytes_RPGC.smooth150.bw")
  pqs_bw_c <- file.path(ROOT, "data/pqsfinder/PQS_scores.mm10.bw")
  stopifnot(file.exists(bw_ast_c), file.exists(bw_non_c), file.exists(pqs_bw_c))

  plots <- list()
  rec_c <- list()
  for (gname in sel_genes_c) {
    gi <- which(gdb_c$gene_name == gname)[1]
    if (is.na(gi)) {
      warning("Gene ", gname, " not found in EnsDb; skipped")
      next
    }
    ggr <- gdb_c[gi]
    fwd <- as.character(strand(ggr)) != "-"
    bpk <- pk_promo_c[gene == gname][which.min(p_val_adj)]
    if (nrow(bpk)) {
      pk_chr <- sub("-.*", "", bpk$peak)
      pk_se <- as.integer(strsplit(sub("^[^-]+-", "", bpk$peak), "-")[[1]])
      anchor <- sprintf("peak %s:%d-%d", sub("^chr", "", pk_chr),
                        pk_se[1], pk_se[2])
    } else {
      pk_chr <- paste0("chr", as.character(seqnames(ggr)))
      tss <- if (fwd) start(ggr) else end(ggr)
      pk_se <- c(tss, tss)
      anchor <- "TSS (no promoter peak)"
    }
    up_ext_c <- 10000L   # upstream of the anchor region
    dn_ext_c <- 35000L   # into the gene body (~30-40 kb requested)
    win_s <- max(1L, if (fwd) pk_se[1] - up_ext_c else pk_se[1] - dn_ext_c)
    win_e <- if (fwd) pk_se[2] + dn_ext_c else pk_se[2] + up_ext_c
    reg_gr <- GenomicRanges::GRanges(pk_chr, IRanges(win_s, win_e))

    tags <- c(if (gname %in% pk_up_c) "peak-up",
              if (gname %in% ga_up_c) "GA-up",
              if (gname %in% rna_up_c) "RNA-up")
    if (!length(tags)) {
      tags <- if (gname %in% ga_dn_c || gname %in% rna_dn_c) "neg ctrl" else "no call"
    }
    title_c <- paste0(gname, " [", paste(tags, collapse = " + "), "]")

    cov_c <- tryCatch(
      CoveragePlot(
        g4, region = reg_gr, annotation = TRUE, peaks = TRUE,
        show.bulk = FALSE,
        bigwig = list("predicted AST (RPGC)" = bw_ast_c,
                      "predicted non-AST (RPGC)" = bw_non_c,
                      "PQS score" = pqs_bw_c),
        bigwig.type = c("coverage", "coverage", "line"),
        bigwig.scale = "common"
      ) & ggtitle(title_c),
      error = function(e) {
        warning(gname, ": CoveragePlot failed (", conditionMessage(e), ")")
        NULL
      }
    )
    if (is.null(cov_c)) next
    plots[[length(plots) + 1]] <- cov_c
    rec_c[[length(rec_c) + 1]] <- data.table(
      gene = gname, anchor = anchor, window_start = win_s, window_end = win_e,
      strand = if (fwd) "+" else "-", tags = paste(tags, collapse = "+")
    )
  }
  sel_out_c <- rbindlist(rec_c, fill = TRUE)
  fwrite(sel_out_c, file.path(OUT, "panel_C_selected_genes.csv"))
  print(sel_out_c)
  cat("  Coverage plots built:", length(plots), "of", length(sel_genes_c),
      "genes\n")
  p_C <- wrap_plots(plots, ncol = 2)
  save_fig(p_C, "panel_C_coverage", 16, 28)
  cat("  Saved panel_C_coverage.pdf\n")
}

# ===========================================================================
# Panel D: Feature plots (Tnik, Pitpnc1, Pbx1, Nwd1) - normalized scRNA-seq
# ===========================================================================
if (run_panel("D")) {
  cat("== Panel D: Feature plots ==\n")

  genes <- c("Tnik", "Pitpnc1", "Pbx1", "Nwd1")
  DefaultAssay(rna) <- "RNA"
  norm <- GetAssayData(rna, layer = "data")

  # Use the same saved co-embedding coordinates as panel A. Re-running
  # RunUMAP on the RNA object produces a different coordinate system.
  rna_umap_d <- umap_df[Domain == "Bartosovic_scRNA-Seq"]
  rna_umap_d[, barcode := as.character(barcode)]

  plots <- list()
  for (g in genes) {
    keep <- match(rna_umap_d$barcode, colnames(rna))
    keep_ok <- !is.na(keep)
    plot_data <- rna_umap_d[keep_ok, ]
    plot_data[, expr := as.numeric(norm[g, keep[keep_ok]])]
    plots[[g]] <- ggplot(plot_data, aes(UMAP1, UMAP2, color = expr)) +
      geom_point_rast(size = 0.8) +
      scale_color_viridis() +
      labs(title = g, color = "Expression") +
      theme_classic() + UMAP_THEME
  }
  # Keep the panel on the same co-embedded coordinate system as panel A.
  p_D <- wrap_plots(plots, ncol = 2)
  save_fig(p_D, "panel_D_featureplots", 14, 12)
  cat("  Saved panel_D_featureplots.pdf\n")
}

# ===========================================================================
# Panel E: Cicero browser tracks
# ===========================================================================
if (run_panel("E")) {
  cat("== Panel E: Cicero browser tracks ==\n")
  source(file.path(SCRIPT_DIR, "cog4_browser.R"))

  cicero_dir <- file.path(RESULTS, "cicero")
  ast_conns <- load_saved_r_object(
    file.path(cicero_dir, "cicero_GFPsorted-predAST.Rds"), "conns"
  )
  nonast_conns <- load_saved_r_object(
    file.path(cicero_dir, "cicero_GFPsorted-pred_nonAST.Rds"), "conns"
  )

  ast_bw <- file.path(
    ROOT, "data", "GSE291468", "GSM8836086_Predicted_Astrocytes_RPGC.bw"
  )
  nonast_bw <- file.path(
    ROOT, "data", "GSE291468", "GSM8836086_Predicted_non-Astrocytes_RPGC.bw"
  )
  pqs_bw <- file.path(ROOT, "data", "pqsfinder", "PQS_scores.mm10.bw")
  ccre_bed <- file.path(
    ROOT, "data", "cCRE", "Li_et_al-mousebrain_cCRE_with_K27ac.bed"
  )
  panel_e_inputs <- c(ast_bw, nonast_bw, pqs_bw, ccre_bed)
  if (any(!file.exists(panel_e_inputs))) {
    stop("Missing Panel E input(s): ", paste(panel_e_inputs[!file.exists(panel_e_inputs)], collapse = ", "))
  }

  ccre <- fread(ccre_bed, header = FALSE, select = 1:4)
  setnames(ccre, c("chr", "start", "end", "k27ac"))
  peaks <- cicero_peak_universe(ast_conns, nonast_conns)
  gene_gr <- GenomicFeatures::genes(EnsDb.Mmusculus.v79)
  gene_chr <- as.character(seqnames(gene_gr))
  gene_chr <- ifelse(startsWith(gene_chr, "chr"), gene_chr, paste0("chr", gene_chr))
  genes <- data.table(
    chr = gene_chr,
    start = start(gene_gr),
    end = end(gene_gr),
    strand = as.character(strand(gene_gr)),
    gene = gene_gr$gene_name
  )[!is.na(gene) & gene != ""]

  panel_e_loci <- c(
    Rsg1 = "chr4-141212918-141213521",
    Tmem74 = "chr15-43869537-43870384",
    Amot = "chrX-145505676-145506438",
    Sox10 = "chr15-79140838-79141742",
    Akt1s1 = "chr7-44848458-44849363"
  )
  panel_e_kinds <- c(
    rep("AST-up G4", 3),
    "oligodendrocyte-lineage constitutive G4",
    "balanced constitutive multi-link G4"
  )
  p_E <- make_cog4_browser_panel(
    loci = panel_e_loci,
    panel_kinds = panel_e_kinds,
    ast_conns = ast_conns,
    nonast_conns = nonast_conns,
    peaks = peaks,
    genes = genes,
    ccre = ccre,
    ast_bigwig = ast_bw,
    nonast_bigwig = nonast_bw,
    pqs_bigwig = pqs_bw
  )
  save_fig(p_E, "panel_E_cicero_tracks", width = 10, height = 27.5)
  cat("  Saved panel_E_cicero_tracks.pdf\n")
}

# ===========================================================================
# Panel S1: ENCODE4 cCRE classes intersected with AST-specific G4 peaks
# ===========================================================================
if (run_panel("S1")) {
  cat("== Panel S1: AST-specific peaks vs ENCODE4 cCRE classes ==\n")

  # AST-specific peak set is exported by Panel C
  # (outputs/GSM8836086_AST_specific_peaks.bed).
  ast_bed_s1 <- file.path(OUT, "GSM8836086_AST_specific_peaks.bed")
  ccre_bed_s1 <- file.path(ROOT, "data", "genome", "cCRE.mm10.bed")
  if (!file.exists(ast_bed_s1))
    stop("Panel S1 needs outputs/GSM8836086_AST_specific_peaks.bed; run Panel C first.")
  if (!file.exists(ccre_bed_s1))
    stop("Panel S1 needs data/genome/cCRE.mm10.bed (ENCODE4 Registry cCREs).")

  canonical_s1 <- c(paste0("chr", 1:19), "chrX", "chrY")
  # ENCODE4 Registry cCRE classes encoded in BED column 9 (itemRgb);
  # official ZLab palette (https://wiki.wenglab.org/references/color-mappings/).
  ccre_label_s1 <- c(
    "255,0,0" = "PLS", "255,167,0" = "pELS", "255,205,0" = "dELS",
    "255,170,170" = "CA-H3K4me3", "0,176,240" = "CA-CTCF",
    "6,218,147" = "CA-only", "190,40,229" = "CA-TF", "216,118,236" = "TF-only"
  )
  ccre_col_s1 <- c(
    PLS = "#FF0000", pELS = "#FFA700", dELS = "#FFCD00",
    "CA-H3K4me3" = "#FFAAAA", "CA-CTCF" = "#00B0F0",
    "CA-only" = "#06DA93", "CA-TF" = "#BE28E5", "TF-only" = "#D876EC"
  )
  ccre_full_s1 <- c(
    dELS = "dELS — distal enhancer-like signature",
    "CA-CTCF" = "CA-CTCF — chromatin accessibility + CTCF",
    "CA-only" = "CA-only — chromatin accessibility only",
    PLS = "PLS — promoter-like signature",
    pELS = "pELS — proximal enhancer-like signature",
    "CA-H3K4me3" = "CA-H3K4me3 — chromatin accessibility + H3K4me3",
    "TF-only" = "TF-only — transcription factor binding only",
    "CA-TF" = "CA-TF — chromatin accessibility + transcription factor"
  )

  rows_s1 <- fread(ast_bed_s1, header = FALSE, fill = TRUE)
  rows_s1 <- rows_s1[V1 %in% canonical_s1]
  ast_s1 <- GenomicRanges::GRanges(
    seqnames = rows_s1$V1, ranges = IRanges(rows_s1$V2 + 1L, rows_s1$V3)
  )
  n_s1 <- length(ast_s1)

  cref_s1 <- fread(ccre_bed_s1, header = FALSE, fill = TRUE)
  cref_s1 <- cref_s1[V1 %in% canonical_s1]
  cref_s1[, cls := unname(ccre_label_s1[V9])]
  cc_s1 <- GenomicRanges::GRanges(
    seqnames = cref_s1$V1, ranges = IRanges(cref_s1$V2 + 1L, cref_s1$V3),
    cls = cref_s1$cls
  )

  # Each AST peak is assigned a single class: the one it overlaps with the
  # greatest bp overlap (ties broken by genomic order); peaks without cCRE
  # overlap form a "No cCRE" category.
  hit_s1 <- as.data.table(GenomicRanges::findOverlaps(ast_s1, cc_s1))
  hit_s1[, cls := cc_s1$cls[subjectHits]]
  hit_s1[, ov_bp := width(pintersect(ast_s1[queryHits], cc_s1[subjectHits]))]
  best_s1 <- hit_s1[
    order(queryHits, -ov_bp, subjectHits),
    .(cls = cls[1]), by = queryHits
  ]
  peak_class_s1 <- rep("No cCRE", n_s1)
  peak_class_s1[best_s1$queryHits] <- best_s1$cls

  tab_s1 <- as.data.table(table(peak_class_s1))[, .(class = peak_class_s1, n = N)]
  tab_s1[, pct := n / sum(n)]
  setorder(tab_s1, -n)
  fwrite(tab_s1, file.path(OUT, "panel_S1_AST_cCRE_partition.csv"))

  lvl_s1 <- unique(c(tab_s1$class[tab_s1$class != "No cCRE"], "No cCRE"))
  tab_s1[, class := factor(class, levels = rev(lvl_s1))]
  tab_s1[, lab := sprintf("%d (%.0f%%)", n, 100 * pct)]
  cols_s1 <- c(ccre_col_s1, "No cCRE" = "grey45")
  leg_lbl_s1 <- function(cl) {
    v <- unname(sapply(as.character(cl),
      function(x) ifelse(x == "No cCRE", "No cCRE overlap", ccre_full_s1[x])))
    factor(v, levels = v)
  }

  p_S1 <- ggplot(tab_s1, aes(x = "AST-specific G4 peaks", y = n, fill = class)) +
    geom_col(width = 0.55, colour = "white", linewidth = 0.4) +
    geom_text(aes(label = lab), colour = "white",
              position = position_stack(vjust = 0.5),
              size = 3.2, fontface = "bold") +
    coord_flip() +
    scale_fill_manual(values = cols_s1, name = NULL, labels = leg_lbl_s1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    labs(x = NULL, y = "number of AST-specific G4 peaks") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = 10),
          axis.text.y = element_blank())
  save_fig(p_S1, "panel_S1_AST_cCRE_stacked", width = 8, height = 3.2)
  cat("  Saved panel_S1_AST_cCRE_stacked.pdf (", n_s1,
      "AST peaks; largest class:", as.character(tab_s1$class[1]),
      sprintf("%d (%.0f%%)", tab_s1$n[1], 100 * tab_s1$pct[1]), ")\n")
}

cat("\nAll Figure 3 outputs written to", OUT, "/\n")
