#!/usr/bin/env Rscript
# ===========================================================================
# Figure 1: G4 profiling by scCUT&Tag separates different cell types
#
# Produces panels B-H. Run from repo root:
#   Rscript fig1/fig1.R
#
# Requires (all from GEO GSE291468):
#   - data/GSE291468/GSM8836088_mESCMEF_Seurat_object.Rds
#   - data/GSE291468/GSM8836088_mESC_MEF/CellRanger/
#   - data/GSE291468/GSM8836088_cluster_0_peaks.narrowPeak
#   - data/GSE291468/GSM8836088_cluster_1_peaks.narrowPeak
#   - data/GSE291468/GSM8836088_cluster0_RPGC.bw
#   - data/GSE291468/GSM8836088_cluster1_RPGC.bw
# ===========================================================================

# --- Packages ---
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(rtracklayer)
  library(EnsDb.Mmusculus.v79)
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(ComplexHeatmap)
  library(circlize)
  library(patchwork)
  library(cowplot)
  library(glue)
  library(bedscout)
  library(wigglescout)
  library(ChIPseeker)
  library(org.Mm.eg.db)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(enrichR)
})

# --- Parse CLI flags ---
# Usage:
#   Rscript fig1/fig1.R                  # run all panels
#   Rscript fig1/fig1.R --only=B,C       # run panels B and C only
#   Rscript fig1/fig1.R --skip=F,G       # run all except F and G
args <- commandArgs(trailingOnly = TRUE)
panels_all <- c("B", "C", "D", "E", "F", "G", "S1", "S2", "S3", "S4", "H")
run_panels <- panels_all

if (any(grepl("^--only=", args))) {
  run_panels <- unlist(strsplit(sub("^--only=", "", args[grep("^--only=", args)]), ","))
  run_panels <- toupper(trimws(run_panels))
}
if (any(grepl("^--only-", args))) {
  run_panels <- toupper(sub("^--only-", "", args[grep("^--only-", args)]))
}
if (any(grepl("^--skip=", args))) {
  skip_panels <- unlist(strsplit(sub("^--skip=", "", args[grep("^--skip=", args)]), ","))
  skip_panels <- toupper(trimws(skip_panels))
  run_panels <- setdiff(run_panels, skip_panels)
}
run_panels <- intersect(run_panels, panels_all)
cat("Running panels:", paste(run_panels, collapse = ", "), "\n")

# --- Paths ---
ROOT     <- normalizePath(getwd())
DATA     <- file.path(ROOT, "data")
RESULTS  <- file.path(ROOT, "results")
OUT      <- file.path(ROOT, "fig1", "outputs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SEURAT_RDS      <- file.path(DATA, "GSE291468/GSM8836088_mESCMEF_Seurat_object.Rds")
FRAGMENTS_TSV   <- file.path(DATA, "GSE291468/GSM8836088_mESC_MEF/CellRanger/fragments.tsv.gz")
CLUSTER_PEAKS_DIR <- file.path(DATA, "GSE291468")
BULK_PEAKS_DIR  <- file.path(DATA, "GSE291468")

# --- Theme ---
theme_fig1 <- theme_classic() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    legend.text = element_text(size = 14)
  )

# --- Helpers ---
load_seurat <- function() readRDS(SEURAT_RDS)

get_markers <- function(seurat, force = FALSE) {
  cache_file <- file.path(OUT, "markers_mesc_mef.rds")
  if (!force && file.exists(cache_file)) return(readRDS(cache_file))
  markers <- FindAllMarkers(seurat, test.use = "LR", latent.vars = "peak_region_fragments", only.pos = TRUE, logfc.threshold = 0.5)
  markers <- markers[markers$p_val_adj < 0.05, ]
  if (nrow(markers) > 0) {
    markers$peak <- sub("\\.[0-9]+$", "", rownames(markers))
  }
  saveRDS(markers, cache_file)
  markers
}

# --- Shared helpers (used by multiple panels) ---
cluster0_file <- file.path(CLUSTER_PEAKS_DIR, "GSM8836088_cluster_0_peaks.narrowPeak")
cluster1_file <- file.path(CLUSTER_PEAKS_DIR, "GSM8836088_cluster_1_peaks.narrowPeak")

# --- Load seurat lazily if any seurat-dependent panel is requested ---
  seurat_needs <- intersect(c("B", "C", "E", "F", "S1", "H"), run_panels)
if (length(seurat_needs) > 0) {
  cat("Loading Seurat object...\n")
  seurat <- load_seurat()
}
cols <- c("0" = "#9ecae1", "1" = "#fc9272")

# ===========================================================================
# Panel B: UMAP
# ===========================================================================
if ("B" %in% run_panels) {
  cat("== Panel B: UMAP ==\n")
  p_B <- DimPlot(seurat, label = TRUE, label.size = 7, repel = TRUE, raster = FALSE) +
    xlim(-10, 10) + ylim(-10, 10) +
    scale_colour_manual(values = cols, breaks = c("0", "1"), labels = c("MEF", "mESC")) +
    ggtitle("") + theme_fig1 +
    theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25))
  ggsave(file.path(OUT, "panel_B_umap.pdf"), p_B, width = 10, height = 10, dpi = 300)
}

# ===========================================================================
# Panel C: QC violin plots
# ===========================================================================
if ("C" %in% run_panels) {
  cat("== Panel C: QC ==\n")
  if (!"FRiP" %in% colnames(seurat@meta.data)) {
    seurat$FRiP <- seurat$peak_region_fragments / seurat$passed_filters
  }
  if (!"TSS_enrichment" %in% colnames(seurat@meta.data) || !"nucleosome_signal" %in% colnames(seurat@meta.data)) {
    frag_list <- Fragments(seurat[["peaks"]])
    old_path <- if (length(frag_list) > 0) frag_list[[1]]@path else ""
    if (length(frag_list) == 0 || !file.exists(old_path)) {
      if (!file.exists(FRAGMENTS_TSV)) {
        stop("Fragment file not found at '", FRAGMENTS_TSV, "'.\n",
             "Download from GEO and place at:\n  ", FRAGMENTS_TSV)
      }
      message("Updating fragment path to: ", FRAGMENTS_TSV)
      for (i in seq_along(frag_list)) frag_list[[i]]@path <- FRAGMENTS_TSV
      seurat[["peaks"]]@fragments <- frag_list
    }
  }
  if (!"TSS_enrichment" %in% colnames(seurat@meta.data)) {
    cat("  Computing TSS enrichment from fragment file (may take a minute)...\n")
    seurat <- TSSEnrichment(seurat, fast = TRUE)
    if ("TSS.enrichment" %in% colnames(seurat@meta.data)) {
      seurat$TSS_enrichment <- seurat$TSS.enrichment
    }
  }
  if (!"nucleosome_signal" %in% colnames(seurat@meta.data)) {
    cat("  Computing nucleosome signal from fragment file (may take a minute)...\n")
    seurat <- NucleosomeSignal(seurat)
  }
  clust_colors <- c("#9ecae1", "#fc9272")
  make_vln <- function(feature, ylab_text, log10 = FALSE) {
    p <- VlnPlot(seurat, group.by = "seurat_clusters", features = feature, pt.size = 0.1) +
      scale_fill_manual(values = clust_colors) +
      ggtitle("") + xlab("cluster") + ylab(ylab_text) +
      stat_summary(fun = median, geom = "crossbar", width = 0.4, linewidth = 0.4,
                   color = "black") +
      stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
      theme_fig1 +
      theme(axis.text.x = element_text(size = 25, angle = 0),
            axis.text.y = element_text(size = 25)) +
      NoLegend()
    if (log10) p <- p + scale_y_log10()
    p
  }
  p_C <- ggarrange(
    make_vln("nFeature_peaks", "nFeature (peaks)", log10 = TRUE),
    make_vln("nCount_peaks", "nCount (peaks)", log10 = TRUE),
    make_vln("TSS_fragments", "TSS fragments", log10 = TRUE),
    make_vln("mitochondrial", "mitochondrial fragments", log10 = TRUE),
    make_vln("FRiP", "Fraction of reads in peaks"),
    make_vln("TSS_enrichment", "TSS enrichment"),
    make_vln("nucleosome_signal", "Nucleosome Signal"),
    ncol = 2, nrow = 4
  )
  ggsave(file.path(OUT, "panel_C_qc.pdf"), p_C, width = 12, height = 20, dpi = 300)
}

# ===========================================================================
# Panel D: Venn diagram
# ===========================================================================
  if ("D" %in% run_panels || "E" %in% run_panels || "G" %in% run_panels || "F" %in% run_panels || "S2" %in% run_panels || "S3" %in% run_panels) {
  # Load peak objects shared by D, E, F, G
  required <- c(cluster0_file, cluster1_file)
  missing <- required[!file.exists(required)]
  if (length(missing) > 0) {
    stop("Missing cluster-specific peaks:\n  ", paste(missing, collapse = "\n  "),
         "\n\nRun data_prep.R first to generate cluster_spec_peaks from pseudobulk MACS3 calls.")
  }
  cl0 <- rtracklayer::import(cluster0_file)
  cl1 <- rtracklayer::import(cluster1_file)
  mesc <- rtracklayer::import(file.path(BULK_PEAKS_DIR, "GSM8836082_bulkG4CnT_mESC_rep1.broadPeak"))
  mef  <- rtracklayer::import(file.path(BULK_PEAKS_DIR, "GSM8836084_bulkG4CnT_3T3_rep1.broadPeak"))
}

  if ("D" %in% run_panels) {
  cat("== Panel D: Venn ==\n")
  euler_from_gr <- function(s1, s2, s3, label1, label2, label3) {
    grl <- GRangesList(s1, s2, s3)
    p <- plot_euler(grl, names = c(label1, label2, label3),
                    fills = c("#9ecae1", "#fc9272", "#f0f0f0"))
    ov12 <- countOverlaps(s1, s2) > 0
    ov13 <- countOverlaps(s1, s3) > 0
    ov23 <- countOverlaps(s2, s3) > 0
    n123 <- sum(ov12 & ov13)
    raw_counts <- c(length(s1) - sum(ov12) - sum(ov13) + n123,
                    length(s2) - sum(ov12) - sum(ov23) + n123,
                    length(s3) - sum(ov13) - sum(ov23) + n123,
                    sum(ov12) - n123, sum(ov13) - n123, sum(ov23) - n123, n123)
    raw_names <- c(label1, label2, label3,
                   paste0(label1, "&", label2), paste0(label1, "&", label3),
                   paste0(label2, "&", label3), paste0(label1, "&", label2, "&", label3))
    keep <- raw_counts > 0
    label_df <- data.frame(region = raw_names[keep], count = raw_counts[keep])
    list(plot = p, counts = data.frame(comparison = NA, region = label_df$region, count = label_df$count))
  }
  mesc_result <- euler_from_gr(cl0, cl1, mesc, "cluster0", "cluster1", "mESC bulk")
  mesc_result$counts$comparison <- "mESC"
  mef_result <- euler_from_gr(cl0, cl1, mef, "cluster0", "cluster1", "MEF bulk")
  mef_result$counts$comparison <- "MEF"
  venn_counts <- rbind(mesc_result$counts, mef_result$counts)
  p_D <- ggarrange(mesc_result$plot, mef_result$plot, ncol = 2)
  ggsave(file.path(OUT, "panel_D_venn.pdf"), p_D, width = 10, height = 5, dpi = 300)
  fwrite(venn_counts, file.path(OUT, "panel_D_venn_counts.csv"))
}

# ===========================================================================
# Define cluster-specific peak sets (used by Panels E, F, S1)
# ===========================================================================
cl0_np <- rtracklayer::import(cluster0_file)
cl1_np <- rtracklayer::import(cluster1_file)

make_summit_gr <- function(np) {
  summit_pos <- start(np) + np$peak
  GRanges(seqnames = seqnames(np),
          ranges = IRanges(start = summit_pos - 250, end = summit_pos + 250))
}
gr0 <- make_summit_gr(cl0_np)
gr1 <- make_summit_gr(cl1_np)

hits <- findOverlaps(gr0, gr1)
shared0 <- unique(queryHits(hits))
shared1 <- unique(subjectHits(hits))

cl0_only <- gr0[-shared0]
cl1_only <- gr1[-shared1]
both <- GenomicRanges::reduce(sort(c(gr0[shared0], gr1[shared1])))

cat("Cluster-specific peak sets: cl0_only =", length(cl0_only),
    ", cl1_only =", length(cl1_only), ", shared (both) =", length(both), "\n")

# ===========================================================================
# Panel E: PCA of 6 samples (2 clusters + 4 bulk)
# ===========================================================================
if ("E" %in% run_panels) {
  cat("== Panel E: PCA ==\n")
  marker_gr <- c(cl0_only, cl1_only)
  cat("  Using", length(marker_gr), "cluster-specific peaks (cl0_only + cl1_only)\n")
  rtracklayer::export.bed(marker_gr, file.path(OUT, "panel_E_marker_regions.mm10.bed"))
  cat("  Saved panel_E_marker_regions.mm10.bed\n")

  # --- (a) Binary overlap PCA (diagnostic) ---
  mesc_r1 <- file.path(BULK_PEAKS_DIR, "GSM8836082_bulkG4CnT_mESC_rep1.broadPeak")
  mesc_r2 <- file.path(BULK_PEAKS_DIR, "GSM8836083_bulkG4CnT_mESC_rep2.broadPeak")
  mef_r1  <- file.path(BULK_PEAKS_DIR, "GSM8836084_bulkG4CnT_3T3_rep1.broadPeak")
  mef_r2  <- file.path(BULK_PEAKS_DIR, "GSM8836085_bulkG4CnT_3T3_rep2.broadPeak")
  stopifnot(all(file.exists(c(cluster0_file, cluster1_file, mesc_r1, mesc_r2, mef_r1, mef_r2))))
  sample_grs <- list(
    "cluster0"       = cl0, "cluster1"      = cl1,
    "mESC bulk rep1" = rtracklayer::import(mesc_r1), "mESC bulk rep2" = rtracklayer::import(mesc_r2),
    "MEF bulk rep1"  = rtracklayer::import(mef_r1),  "MEF bulk rep2"  = rtracklayer::import(mef_r2)
  )
  peak_mat <- t(sapply(sample_grs, function(gr) as.integer(countOverlaps(marker_gr, gr) > 0)))
  col_var <- apply(peak_mat, 2, var)
  peak_mat <- peak_mat[, col_var > 0]
  pca_res <- prcomp(peak_mat, scale. = TRUE, center = TRUE)
  var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
  pca_df <- data.frame(
    PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2],
    sample = rownames(pca_res$x),
    group = c("cluster", "cluster", "mESC", "mESC", "MEF", "MEF")
  )
  sample_colors <- c(
    "cluster0" = "#9ecae1", "mESC bulk rep1" = "#d62728", "mESC bulk rep2" = "#d62728",
    "cluster1" = "#fc9272", "MEF bulk rep1" = "#1f77b4", "MEF bulk rep2" = "#1f77b4"
  )
  p_E <- ggplot(pca_df, aes(x = PC1, y = PC2, color = sample)) +
    geom_point(size = 6) +
    scale_color_manual(values = sample_colors) +
    labs(title = "PCA of peak profiles (binary)",
         x = glue("PC1 ({var_exp[1]}%)"), y = glue("PC2 ({var_exp[2]}%)")) +
    theme_fig1 + geom_text(aes(label = sample), vjust = -1, size = 3) +
    NoLegend()
   ggsave(file.path(OUT, "panel_E_binary.pdf"), p_E, width = 8, height = 6, dpi = 300)
  cat("  Saved panel_E_binary.pdf\n")

  # --- (b) Signal-based PCA (RPGC bigwigs x marker regions) ---
  bw_6 <- c(
    "cluster0 (MEF)" = file.path(DATA, "GSE291468/GSM8836088_cluster0_RPGC.bw"),
    "cluster1 (mESC)" = file.path(DATA, "GSE291468/GSM8836088_cluster1_RPGC.bw"),
    "mESC rep1"       = file.path(DATA, "GSE291468/GSM8836082_bulkG4CnT_mESC_rep1.bw"),
    "mESC rep2"       = file.path(DATA, "GSE291468/GSM8836083_bulkG4CnT_mESC_rep2.bw"),
    "MEF rep1"        = file.path(DATA, "GSE291468/GSM8836084_bulkG4CnT_3T3_rep1.bw"),
    "MEF rep2"        = file.path(DATA, "GSE291468/GSM8836085_bulkG4CnT_3T3_rep2.bw")
  )
  missing_bw <- bw_6[!file.exists(unlist(bw_6))]
  if (length(missing_bw) > 0) {
    cat("    Skipping signal PCA: missing bigwigs:", paste(names(missing_bw), collapse = ", "), "\n")
  } else {
     sig <- bw_loci(unlist(bw_6), marker_gr, labels = names(bw_6), default_na = 0)
     sig_mat <- as.matrix(mcols(sig))
     if ("name" %in% colnames(sig_mat)) sig_mat <- sig_mat[, setdiff(colnames(sig_mat), "name"), drop = FALSE]
    row_var <- apply(sig_mat, 1, var)
    sig_mat <- sig_mat[row_var > 0, ]
     pca_sig <- prcomp(t(sig_mat), scale. = TRUE, center = TRUE)
    var_sig <- round(summary(pca_sig)$importance[2, 1:2] * 100, 1)
    pca_sig_df <- data.frame(
      PC1 = pca_sig$x[, 1], PC2 = pca_sig$x[, 2],
       sample = rownames(pca_sig$x),
      group = c("cluster", "cluster", "mESC", "mESC", "MEF", "MEF")
    )
    sample_colors_sig <- c(
      "cluster0 (MEF)" = "#9ecae1", "mESC rep1" = "#d62728", "mESC rep2" = "#d62728",
      "cluster1 (mESC)" = "#fc9272", "MEF rep1" = "#1f77b4", "MEF rep2" = "#1f77b4"
    )
    p_E_sig <- ggplot(pca_sig_df, aes(x = PC1, y = PC2, color = sample)) +
      geom_point(size = 6) +
      scale_color_manual(values = sample_colors_sig) +
      labs(title = "PCA of signal profiles",
           x = glue("PC1 ({var_sig[1]}%)"), y = glue("PC2 ({var_sig[2]}%)")) +
      theme_fig1 + geom_text(aes(label = sample), vjust = -1, size = 3) +
      NoLegend()
     # The manuscript panel is based on normalized signal, not binary peak
     # membership. Keep the binary version as a diagnostic alternative.
     ggsave(file.path(OUT, "panel_E_pca.pdf"), p_E_sig, width = 8, height = 6, dpi = 300)
     cat("  Saved panel_E_pca.pdf\n")
  }
}

# ===========================================================================
# Panel F: G4 occupancy heatmaps & average profiles (wigglescout)
# ===========================================================================
if ("F" %in% run_panels) {
  cat("== Panel F: wigglescout heatmaps & profiles ==\n")

  names(cl0_only) <- rep("cluster0", length(cl0_only))
  names(both) <- rep("both", length(both))
  names(cl1_only) <- rep("cluster1", length(cl1_only))
  merged_peaks <- c(cl0_only, both, cl1_only)
  rtracklayer::export.bed(merged_peaks, file.path(OUT, "panel_F_merged_peaks.bed"))
  cat("  Saved panel_F_merged_peaks.bed:", length(merged_peaks), "regions\n")

  peak_categories <- list(cl0 = cl0_only, both = both, cl1 = cl1_only)

  bw_files <- list(
    "bulk mESC"    = file.path(DATA, "GSE291468/GSM8836082_bulkG4CnT_mESC_rep1.bw"),
    "bulk MEF"     = file.path(DATA, "GSE291468/GSM8836084_bulkG4CnT_3T3_rep1.bw"),
    "mESC cl1"     = file.path(DATA, "GSE291468/GSM8836088_cluster1_RPGC.bw"),
    "MEF cl0"      = file.path(DATA, "GSE291468/GSM8836088_cluster0_RPGC.bw"),
    "mESC ATAC"    = file.path(DATA, "GSE149080/GSM4661960_ATAC_ESC_WT_batch2.rpgc.bw"),
    "MEF ATAC"     = file.path(DATA, "GSE211123/GSM6451000_ATAC_3T3.rpgc.bw"),
    "PQS score"    = file.path(DATA, "pqsfinder/PQS_scores.mm10.bw")
  )
  heatmap_cmaps <- c(
    "bulk mESC" = "Reds",
    "bulk MEF" = "Blues",
    "mESC cl1" = "Reds",
    "MEF cl0" = "Blues",
    "mESC ATAC" = "Reds",
    "MEF ATAC" = "Blues",
    "PQS score" = "Purples"
  )

  # --- Heatmaps: 3 peak categories x 6 bigwigs = 18 plots, one PDF per bigwig ---
  # Use a common zmax per bigwig so all 3 heatmaps share the same color scale
  cat("  Generating heatmaps...\n")
  compute_zmax <- function(bw_path, gr_list, upstream = 3000, downstream = 3000) {
    bw <- BigWigFile(bw_path)
    all_vals <- c()
    for (gr in gr_list) {
      centers <- start(gr) + (end(gr) - start(gr)) %/% 2
      query_gr <- GRanges(seqnames = seqnames(gr),
                          ranges = IRanges(start = centers - upstream, end = centers + downstream))
      sig <- import(bw, selection = query_gr, as = "NumericList")
      vals <- unlist(lapply(sig, function(x) x[!is.na(x)]))
      if (length(vals) > 0) all_vals <- c(all_vals, vals)
    }
    if (length(all_vals) == 0) return(10)
    quantile(all_vals, 0.99, na.rm = TRUE)
  }

  all_cat_grs <- unname(peak_categories)
  # Fix row ordering from the cluster tracks before rendering any heatmap.
  cluster_order_tracks <- c(
    cl0 = bw_files[["MEF cl0"]],
    cl1 = bw_files[["mESC cl1"]]
  )
  heatmap_orders <- lapply(names(peak_categories), function(cat_name) {
    gr <- peak_categories[[cat_name]]
    if (!length(gr)) return(integer())
    cluster_signal <- as.data.frame(mcols(bw_loci(
      cluster_order_tracks, gr, labels = names(cluster_order_tracks),
      default_na = 0
    )))
    score <- switch(cat_name,
      cl0 = cluster_signal$cl0,
      cl1 = cluster_signal$cl1,
      both = cluster_signal$cl0 * cluster_signal$cl1
    )
    order(score)
  })
  names(heatmap_orders) <- names(peak_categories)
  for (bw_name in names(bw_files)) {
    bw_path <- bw_files[[bw_name]]
    if (!file.exists(bw_path)) {
      message("  Skipping heatmap for '", bw_name, "': file not found: ", bw_path)
      next
    }
    zmax_val <- if (bw_name == "PQS score") 20 else compute_zmax(bw_path, all_cat_grs)
    message("  zmax for ", bw_name, ": ", round(zmax_val, 2))

    # Combine the average profile and the three category heatmaps into one
    # panel. Canvas heights are proportional to the original category sizes,
    # while wigglescout is capped at 100 displayed rows per heatmap.
    cat("  Generating composite profile and heatmap...\n")
    cat_colors <- c("#2166ac", "#333333", "#b2182b")
    profile <- plot_bw_profile(
      bw_path, peak_categories,
      upstream = 3000, downstream = 3000,
      mode = "center", colors = cat_colors,
      labels = names(peak_categories), default_na = 0,
      show_error = TRUE, verbose = FALSE
    )
    # Make the profile 10% wider than the previous size and align it left.
    profile <- cowplot::ggdraw() + cowplot::draw_plot(
      profile, x = 0, y = 0.1875, width = 0.825, height = 0.625
    )
    heatmaps <- lapply(names(peak_categories), function(cat_name) {
      gr <- peak_categories[[cat_name]]
      plot_bw_heatmap(
        bw_path, gr,
        upstream = 3000, downstream = 3000,
        mode = "center",
        zmin = if (bw_name == "PQS score") 5 else 0,
        zmax = zmax_val,
        # wigglescout accepts named palettes here, not a raw colour vector.
        # Purples provides the requested white-to-purple appearance.
        cmap = if (bw_name == "PQS score") "RdPu" else heatmap_cmaps[[bw_name]],
        max_rows_allowed = max(1, ceiling(length(gr) / 100)),
        order_by = heatmap_orders[[cat_name]],
        default_na = 0, verbose = FALSE
      ) + ggtitle(cat_name)
    })
    # Profile gets a compact header; each heatmap receives height according to
    # its source row count (10,000 rows is approximately six inches).
    panel <- ggarrange(
      profile,
      ggarrange(plotlist = heatmaps, ncol = 1,
                heights = pmax(1, vapply(peak_categories, length, numeric(1)) / 10000),
                align = "v"),
      ncol = 1,
      # Keep the profile square while preserving the enlarged heatmap section.
      heights = c(6, max(1, 3 * max(8, 3 + sum(pmax(1,
        vapply(peak_categories, length, numeric(1)) / 10000))) - 6)),
      align = "v"
    )
    composite_file <- file.path(OUT, paste0("panel_F_", make.names(bw_name), ".pdf"))
    composite_height <- 3 * max(8, 3 + sum(pmax(1, vapply(peak_categories, length, numeric(1)) / 10000)))
    ggsave(composite_file, panel, width = 6, height = composite_height, limitsize = FALSE)
    message("  Saved ", basename(composite_file), " (height ", round(composite_height, 2), ")")
  }

  # Figure 1F: Gene annotation and enrichR analysis for the three peak categories
  cat("  Annotating peaks to genes and running enrichR...\n")
  annotate_set <- function(gr) {
    if (!length(gr)) return(character(0))
    tss_gr <- rtracklayer::import(file.path(DATA, "genome/refGene.tss.1kb.mm10.bed"))
    anno_dt <- as.data.table(as.data.frame(annotate_nearby_features(
      gr, tss_gr, name_field = "name", distance_cutoff = 1000
    )))
    anno_dt <- anno_dt[!is.na(annotation) & annotation != ""]
    if (!nrow(anno_dt)) return(character(0))
    gene_symbols <- sub(";.*", "", anno_dt$annotation)
    unique(gene_symbols[gene_symbols != ""])
  }

  names(peak_categories) <- c("cl0_only", "both", "cl1_only")
  annotated <- lapply(peak_categories, annotate_set)
  genes <- lapply(annotated, function(x) {
    if (length(x) == 0) return(character(0))
    x
  })
  cat("    Annotated genes:", paste(names(genes), vapply(genes, length, integer(1)), collapse = "; "), "\n")

  gene_table <- rbindlist(lapply(names(genes), function(name) {
    data.table(set = name, gene = genes[[name]])
  }), fill = TRUE)
  fwrite(gene_table, file.path(OUT, "panel_F_gene_sets.tsv"), sep = "\t")

  dbs <- c("GO_Biological_Process_2023")
  enrich <- lapply(genes, function(g) {
    if (length(g) < 2) return(NULL)
    enrichr(g, dbs)
  })

  bp <- rbindlist(lapply(names(enrich), function(name) {
    if (is.null(enrich[[name]])) return(NULL)
    x <- as.data.table(enrich[[name]][["GO_Biological_Process_2023"]])
    if (!nrow(x)) return(NULL)
    x[, set := name]
    x
  }), fill = TRUE)

  if (nrow(bp)) {
    bp_sig <- bp[P.value < 0.05, .(set, Term, P.value, Adjusted.P.value, Combined.Score, Genes)]
    fwrite(bp, file.path(OUT, "panel_F_GO_Biological_Process_2023.tsv"), sep = "\t")
    fwrite(bp_sig, file.path(OUT, "panel_F_GO_Biological_Process_2023_significant.tsv"), sep = "\t")
    cat("    Significant GO BP terms:", nrow(bp_sig), "\n")
  } else {
    message("    No significant GO Biological Process terms found")
  }
}

if ("G" %in% run_panels) {
  enrich_terms_file <- file.path(OUT, "panel_F_GO_Biological_Process_2023_significant.tsv")
  if (file.exists(enrich_terms_file)) {
    enrich_terms <- fread(enrich_terms_file)
    enrich_terms <- enrich_terms[order(P.value), head(.SD, 10), by = set]
    enrich_wide <- dcast(enrich_terms, Term ~ set, value.var = "P.value", fill = 1)
    enrich_sets <- c("cl0_only", "both", "cl1_only")
    for (set_name in enrich_sets) {
      if (!set_name %in% names(enrich_wide)) enrich_wide[, (set_name) := 1]
    }
    enrich_mat <- as.matrix(enrich_wide[, ..enrich_sets])
    rownames(enrich_mat) <- enrich_wide$Term
    enrich_logp <- -log10(pmax(enrich_mat, .Machine$double.xmin))
    term_names <- sub("\\s*\\(GO:[0-9]+\\)", "", rownames(enrich_mat))
    dimnames(enrich_logp) <- list(term_names, colnames(enrich_mat))
    enrich_col_fun <- colorRamp2(c(0, 1, 2, 4, 8),
      c("grey95", "#fee0d2", "#fc9272", "#de2d26", "#99000d"))
    enrich_heatmap <- Heatmap(
      enrich_logp,
      name = "-log10 p",
      col = enrich_col_fun,
      column_title = "Figure 1F peak sets: top 10 enrichR GO BP terms",
      cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE,
      rect_gp = gpar(col = "black", lwd = 0.1),
      show_row_names = TRUE,
      row_names_gp = gpar(fontsize = 9),
      row_names_max_width = unit(8, "cm"),
      column_names_gp = gpar(fontsize = 9), column_names_rot = 45,
      # The width includes the term-label viewport; keep the matrix body
      # narrow so its three cells remain close to square.
      heatmap_width = unit(15, "cm"),
      heatmap_height = unit(max(5, nrow(enrich_logp) * 0.6 + 1), "cm")
    )
    enrich_file <- file.path(OUT, "panel_G_enrichR_heatmap.pdf")
    pdf(enrich_file, width = 18, height = max(5, nrow(enrich_logp) * 0.25 + 1))
    draw(enrich_heatmap)
    dev.off()
    fwrite(enrich_terms, file.path(OUT, "panel_G_enrichR_GO_terms.csv"))
    cat("  Saved panel_G_enrichR_heatmap.pdf and panel_G_enrichR_GO_terms.csv\n")
  } else {
    message("  Skipping panel_G: enrichR results not found (run Panel F first)")
  }
}

# ===========================================================================
# Panel S2/S3/S4: BigWig signal correlations (with optional PQS highlights)
# ===========================================================================
if (any(c("S2", "S3", "S4") %in% run_panels)) {
  # Reuse cluster-specific peak sets defined at top of script
  peak_sets <- list(
    cl0_only = cl0_only,
    cl1_only = cl1_only,
    shared = both
  )

  bw_esc_atac <- file.path(DATA, "GSE149080/GSM4661960_ATAC_ESC_WT_batch2.rpgc.bw")
  bw_mef_atac <- file.path(DATA, "GSE211123/GSM6451000_ATAC_3T3.rpgc.bw")
  bw_cl0 <- file.path(DATA, "GSE291468/GSM8836088_cluster0_RPGC.bw")
  bw_cl1 <- file.path(DATA, "GSE291468/GSM8836088_cluster1_RPGC.bw")
  stopifnot(all(file.exists(c(bw_esc_atac, bw_mef_atac, bw_cl0, bw_cl1))))

  mef_peaks <- rtracklayer::import(
    file.path(DATA, "GSE211123/GSM6451000_3T3_control_ATAC_peaks.narrowPeak")
  )
  esc_peaks <- rtracklayer::import(
    file.path(DATA, "GSE149080/GSM4661960_ESC_WT_batch2_peaks.narrowPeak")
  )

  correlation_rows <- list()
  make_scatter <- function(set_name, x, y, x_label, y_label, plot_id = set_name,
                           highlight = NULL, loci = NULL,
                           hl_colors = c("#d7301f"),
                           point_alpha = 0.2) {
    if (is.null(loci)) loci <- peak_sets[[set_name]]
    signal <- as.data.frame(mcols(bw_loci(
      c(x = x, y = y), loci, labels = c("x", "y"), default_na = NA
    )))
    x_values <- signal$x
    y_values <- signal$y
    keep <- is.finite(x_values) & is.finite(y_values)
    pearson <- if (sum(keep) >= 2) cor(x_values[keep], y_values[keep], method = "pearson") else NA_real_
    correlation_rows[[plot_id]] <<- data.table(
      peak_set = set_name,
      plot_id = plot_id,
      x_track = x_label,
      y_track = y_label,
      n_regions = sum(keep),
      pearson = pearson
    )

    scatter <- plot_bw_loci_scatter(
      x, y, loci,
      norm_mode_x = "none", norm_mode_y = "none",
      highlight = highlight,
      highlight_label = if (is.null(highlight)) NULL else "red2",
      highlight_colors = hl_colors,
      minoverlap = 1L,
      default_na = NA_real_, verbose = FALSE
    )

    for (layer in scatter$layers) {
      if (inherits(layer$geom, "GeomPoint")) {
        layer$aes_params$size <- 0.75
        layer$aes_params$alpha <- point_alpha
      }
    }
    # Rasterize the point layers (300 dpi) so the scatter plots stay manageable
    # when imported into Illustrator.
    ggrastr::rasterise(
      scatter +
        geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.4) +
        scale_x_continuous(limits = c(1, 1024), trans = "log2") +
        scale_y_continuous(limits = c(1, 1024), trans = "log2") +
        labs(
          title = plot_id,
          x = paste0(x_label, " (log2)"),
          y = paste0(y_label, " (log2)"),
          subtitle = paste0("Pearson r = ", formatC(pearson, format = "f", digits = 3),
                            " (n = ", sum(keep), ")")
        ),
      dpi = 150
    )
  }
}

# Panel S2: ATAC and pseudobulk signal correlations
# ===========================================================================
if ("S2" %in% run_panels) {
  cat("== Panel S2: ATAC and pseudobulk signal correlations ==\n")

  p_s2 <- wrap_plots(
    make_scatter("cl0_only", bw_cl0, bw_mef_atac, "cluster 0 RPGC", "MEF ATAC RPGC"),
    make_scatter("cl1_only", bw_cl1, bw_esc_atac, "cluster 1 RPGC", "mESC ATAC RPGC"),
    make_scatter("shared", bw_cl0, bw_cl1, "cluster 0 RPGC", "cluster 1 RPGC", "shared: cl0 vs cl1"),
    make_scatter("shared", bw_cl0, bw_mef_atac, "cluster 0 RPGC", "MEF ATAC RPGC", "shared: cl0 vs MEF ATAC"),
    make_scatter("shared", bw_cl1, bw_esc_atac, "cluster 1 RPGC", "mESC ATAC RPGC", "shared: cl1 vs mESC ATAC"),
    ncol = 3
  )
  ggsave(file.path(OUT, "panel_S2_G4_ATAC_correlation.pdf"), p_s2, width = 18, height = 12, device = "pdf")
  fwrite(rbindlist(correlation_rows), file.path(OUT, "panel_S2_G4_ATAC_correlation_pearson.csv"))
  cat("  Saved panel_S2_G4_ATAC_correlation.pdf and Pearson summary CSV\n")
}

# Panel S3: BigWig signal correlations with PQS-overlapping regions highlighted
# ===========================================================================
if ("S3" %in% run_panels) {
  cat("== Panel S3: ATAC and pseudobulk correlations with PQS highlights ==\n")
  pqs_bed <- file.path(DATA, "pqsfinder/PQS_scores.min50.mm10.bed") # using min 50 PQS score here
  if (!file.exists(pqs_bed)) stop("PQS BED not found: ", pqs_bed)
 # pqs_gr <- rtracklayer::import(pqs_bed)
  p_s3 <- wrap_plots(
    make_scatter("cl0_only", bw_cl0, bw_mef_atac, "cluster 0 RPGC", "MEF ATAC RPGC",
                 "cl0-only: cl0 vs MEF ATAC", highlight = pqs_bed),
    make_scatter("cl1_only", bw_cl1, bw_esc_atac, "cluster 1 RPGC", "mESC ATAC RPGC",
                 "cl1-only: cl1 vs mESC ATAC", highlight = pqs_bed),
    make_scatter("shared", bw_cl0, bw_cl1, "cluster 0 RPGC", "cluster 1 RPGC",
                 "shared: cl0 vs cl1", highlight = pqs_bed),
    make_scatter("shared", bw_cl0, bw_mef_atac, "cluster 0 RPGC", "MEF ATAC RPGC",
                 "shared: cl0 vs MEF ATAC", highlight = pqs_bed),
    make_scatter("shared", bw_cl1, bw_esc_atac, "cluster 1 RPGC", "mESC ATAC RPGC",
                 "shared: cl1 vs mESC ATAC", highlight = pqs_bed),
    ncol = 3
  )
  ggsave(file.path(OUT, "panel_S3_G4_ATAC_correlation_PQS.pdf"),
         p_s3, width = 18, height = 12, device = "pdf")
  cat("  Saved panel_S3_G4_ATAC_correlation_PQS.pdf\n")

  # Alternate S3 panels: pseudobulk peak sets merged with the matching ATAC-seq
  # peak set (concatenate + reduce), one scatter per set, saved as separate PDFs.
  alt_sets <- list(
    list(name = "MEF_ATAC",
         loci = mef_peaks,
         x = bw_cl0, y = bw_mef_atac,
         x_label = "cluster 0 RPGC", y_label = "MEF ATAC RPGC"),
    list(name = "ESC_ATAC",
         loci = esc_peaks,
         x = bw_cl1, y = bw_esc_atac,
         x_label = "cluster 1 RPGC", y_label = "mESC ATAC RPGC"),
    list(name = "shared_plus_MEF_ATAC",
         loci = GenomicRanges::reduce(c(peak_sets$shared, mef_peaks)),
         x = bw_cl0, y = bw_mef_atac,
         x_label = "cluster 0 RPGC_ATAC", y_label = "MEF ATAC RPGC"),
    list(name = "shared_plus_ESC",
         loci = GenomicRanges::reduce(c(peak_sets$shared, esc_peaks)),
         x = bw_cl1, y = bw_esc_atac,
         x_label = "cluster 1 RPGC", y_label = "mESC ATAC RPGC")
  )
  # make_scatter already rasterizes the point layers; combine the four merged
  # peak-set scatters into one 3x2 grid.
  alt_plots <- lapply(alt_sets, function(alt) {
    make_scatter(alt$name, alt$x, alt$y, alt$x_label, alt$y_label,
                 plot_id = sub("_plus_", " + ", alt$name),
                 loci = alt$loci,
                 highlight = pqs_bed,
                 point_alpha = 0.2)
  })
  p_alt <- wrap_plots(alt_plots, ncol = 3)
  ggsave(file.path(OUT, "panel_S3_alt_merged_peaks.pdf"),
         p_alt, width = 18, height = 12, device = "pdf")
  cat("  Saved panel_S3_alt_merged_peaks.pdf (3x2 rasterised grid)\n")
}


# Panel S4: ggscatterhist alternative (cluster peaks + matching ATAC-seq peaks)
# ==============================================================================
if ("S4" %in% run_panels) {
  cat("== Panel S4: ggscatterhist with shared log2 axes ==\n")

  # Signals are computed with bw_loci, log2-transformed and clipped to the same
  # [0, 10] axis limits (= log2(1) .. log2(1024)) used by the S2/S3 scatters so
  # all ggscatterhist panels share identical, log-scaled axes.
  get_sh_sig <- function(x_bw, y_bw, gr) {
    as.data.frame(mcols(bw_loci(
      c(x = x_bw, y = y_bw), gr, labels = c("x", "y"), default_na = NA
    ))[, c("x", "y")])
  }
  prep_sh <- function(sig) {
    sig$x <- log2(sig$x)
    sig$y <- log2(sig$y)
    keep <- is.finite(sig$x) & is.finite(sig$y) &
      sig$x >= 0 & sig$x <= 10 & sig$y >= 0 & sig$y <= 10
    sig[keep, , drop = FALSE]
  }
  # Fair comparison: narrow the ATAC-seq peaks to their summits +/-250 bp, the
  # same window used for the pseudobulk cluster peaks (make_summit_gr above).
  mef_peaks_s <- make_summit_gr(mef_peaks)
  esc_peaks_s <- make_summit_gr(esc_peaks)

  # Eight single-group panels: ATAC and cluster peaks are plotted separately to
  # avoid overplotting. Red = cluster peaks / PQS+ (alpha 0.2); black = ATAC
  # peaks / non-PQS (alpha 0.1). Point z-order follows data row order, but with
  # single-group panels there is no stacking issue.
  sh_cl0_red  <- prep_sh(transform(get_sh_sig(bw_cl0, bw_mef_atac, peak_sets$cl0_only), pt_alpha = 0.2))
  sh_mef_black <- prep_sh(transform(get_sh_sig(bw_cl0, bw_mef_atac, mef_peaks_s), pt_alpha = 0.1))
  sh_cl1_red  <- prep_sh(transform(get_sh_sig(bw_cl1, bw_esc_atac, peak_sets$cl1_only), pt_alpha = 0.2))
  sh_esc_black <- prep_sh(transform(get_sh_sig(bw_cl1, bw_esc_atac, esc_peaks_s), pt_alpha = 0.1))

  # Merged peak-set panels: union of the pseudobulk cluster peaks and the
  # matching ATAC-seq peaks (both summit +/-250), split into PQS+ (red,
  # overlapping PQS_scores.min50.mm10.bed) and non-PQS (black, the rest).
  pqs_gr <- rtracklayer::import(file.path(DATA, "pqsfinder/PQS_scores.min50.mm10.bed"))
  split_merged_sh <- function(x_bw, y_bw, merged_gr) {
    sig <- get_sh_sig(x_bw, y_bw, merged_gr)
    ov <- findOverlaps(merged_gr, pqs_gr)
    is_pqs <- seq_along(merged_gr) %in% unique(queryHits(ov))
    list(
      pqs = prep_sh(transform(sig[is_pqs, , drop = FALSE], pt_alpha = 0.2)),
      non_pqs = prep_sh(transform(sig[!is_pqs, , drop = FALSE], pt_alpha = 0.1))
    )
  }
  sh_merged0 <- split_merged_sh(
    bw_cl0, bw_mef_atac, GenomicRanges::reduce(c(peak_sets$cl0_only, mef_peaks_s)))
  sh_merged1 <- split_merged_sh(
    bw_cl1, bw_esc_atac, GenomicRanges::reduce(c(peak_sets$cl1_only, esc_peaks_s)))

  # ggscatterhist returns a list of ggplots; assemble each into a single
  # ggplot via its print method on a null device before combining.
  assemble_sh <- function(p) {
    cowplot::pdf_null_device(width = 7, height = 7)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
  }
  make_sh_panel <- function(dat, x_label, y_label, title, color) {
    p <- ggscatterhist(
      dat, x = "x", y = "y", color = color,
      margin.params = list(fill = color, color = "black", size = 0.2),
      xlab = paste0(x_label, " (log2)"), ylab = paste0(y_label, " (log2)"),
      title = title, shape = 16, alpha = "pt_alpha", print = FALSE)
    p$sp <- p$sp +
      scale_x_continuous(limits = c(0, 10)) +
      scale_y_continuous(limits = c(0, 10)) +
      scale_alpha_identity()
    p$sp <- ggrastr::rasterise(p$sp, layers = "Point", dpi = 300)
    assemble_sh(p)
  }
  # cl0 / MEF row
  p_s4_1 <- make_sh_panel(sh_cl0_red, "cluster 0 RPGC", "MEF ATAC RPGC",
                          "cl0 peaks", "red2")
  p_s4_2 <- make_sh_panel(sh_mef_black, "cluster 0 RPGC", "MEF ATAC RPGC",
                          "MEF ATAC peaks", "black")
  p_s4_3 <- make_sh_panel(sh_merged0$pqs, "cluster 0 RPGC", "MEF ATAC RPGC",
                          "merged cl0+MEF: PQS+ (min50)", "red2")
  p_s4_4 <- make_sh_panel(sh_merged0$non_pqs, "cluster 0 RPGC", "MEF ATAC RPGC",
                          "merged cl0+MEF: non-PQS", "black")
  # cl1 / ESC row
  p_s4_5 <- make_sh_panel(sh_cl1_red, "cluster 1 RPGC", "mESC ATAC RPGC",
                          "cl1 peaks", "red2")
  p_s4_6 <- make_sh_panel(sh_esc_black, "cluster 1 RPGC", "mESC ATAC RPGC",
                          "ESC ATAC peaks", "black")
  p_s4_7 <- make_sh_panel(sh_merged1$pqs, "cluster 1 RPGC", "mESC ATAC RPGC",
                          "merged cl1+ESC: PQS+ (min50)", "red2")
  p_s4_8 <- make_sh_panel(sh_merged1$non_pqs, "cluster 1 RPGC", "mESC ATAC RPGC",
                          "merged cl1+ESC: non-PQS", "black")
  ggsave(file.path(OUT, "panel_S4_ggscatterhist.pdf"),
         wrap_plots(p_s4_1, p_s4_2, p_s4_3, p_s4_4,
                    p_s4_5, p_s4_6, p_s4_7, p_s4_8, ncol = 4),
         width = 28, height = 14, device = "pdf")
  cat("  Saved panel_S4_ggscatterhist.pdf\n")
}


# Panel S1: Spearman similarity heatmap
# ===========================================================================
if ("S1" %in% run_panels) {
  cat("== Panel S1: Spearman heatmap ==\n")

  # Main analysis: cluster-specific peaks only (no overlap)
  marker_gr <- c(cl0_only, cl1_only)
  cat("  Using", length(marker_gr), "cluster-specific peaks (cl0_only + cl1_only)\n")
  rtracklayer::export.bed(marker_gr, file.path(OUT, "panel_S1_marker_regions.mm10.bed"))
  cat("  Saved panel_S1_marker_regions.mm10.bed\n")

  # For |lfc|>2 filtered analysis: include shared peaks since strong fold-change
  # ensures specificity regardless of peak-calling overlap
  marker_gr_all <- c(cl0_only, cl1_only, both)
  cat("  For |lfc|>2 analysis: ", length(marker_gr_all), " peaks (including shared)\n", sep="")

  bw_4 <- c(
    cl0  = file.path(DATA, "GSE291468/GSM8836088_cluster0_RPGC.bw"),
    cl1  = file.path(DATA, "GSE291468/GSM8836088_cluster1_RPGC.bw"),
    mESC = file.path(DATA, "GSE291468/GSM8836082_bulkG4CnT_mESC_rep1.bw"),
    MEF  = file.path(DATA, "GSE291468/GSM8836084_bulkG4CnT_3T3_rep1.bw")
  )
missing_bw <- bw_4[!file.exists(unlist(bw_4))]
  if (length(missing_bw) > 0) {
    cat("    Skipping correlation: missing bigwigs:", paste(names(missing_bw), collapse = ", "), "\n")
  } else {
    sig <- bw_loci(unlist(bw_4), marker_gr, labels = names(bw_4), default_na = 0)
    sig_mat <- as.matrix(mcols(sig))

    # Compute both Spearman and Pearson correlations for all cluster-specific peaks
    cor_spearman <- cor(sig_mat, method = "spearman", use = "complete.obs")
    cor_pearson <- cor(sig_mat, method = "pearson", use = "complete.obs")

    order_G <- c("MEF", "cl0", "cl1", "mESC")
    cor_spearman <- cor_spearman[order_G, order_G, drop = FALSE]
    cor_pearson <- cor_pearson[order_G, order_G, drop = FALSE]

    display_names <- c("MEF bulk", "cluster0", "cluster1", "mESC bulk")
    dimnames(cor_spearman) <- list(display_names, display_names)
    dimnames(cor_pearson) <- list(display_names, display_names)

    col_fun_G <- colorRamp2(c(0.3, 1), c("white", "red2"))

    heatmap_spearman <- Heatmap(cor_spearman, name = "Spearman rho", col = col_fun_G,
                column_title = "Spearman correlation", row_title = "",
            rect_gp = gpar(col = "black", lwd = 1),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.2f", cor_spearman[i, j]), x, y, gp = gpar(fontsize = 14))
            },
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_dend = FALSE, show_column_dend = FALSE,
            heatmap_width = unit(6, "cm"), heatmap_height = unit(6, "cm"),
            row_names_gp = gpar(fontsize = 12), column_names_gp = gpar(fontsize = 12))

    heatmap_pearson <- Heatmap(cor_pearson, name = "Pearson r", col = col_fun_G,
                column_title = "Pearson correlation", row_title = "",
            rect_gp = gpar(col = "black", lwd = 1),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.2f", cor_pearson[i, j]), x, y, gp = gpar(fontsize = 14))
            },
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_dend = FALSE, show_column_dend = FALSE,
            heatmap_width = unit(6, "cm"), heatmap_height = unit(6, "cm"),
            row_names_gp = gpar(fontsize = 12), column_names_gp = gpar(fontsize = 12))

    pdf(file.path(OUT, "panel_S1.pdf"), width = 15, height = 7.5)
    draw(heatmap_spearman + heatmap_pearson, ht_gap = unit(1, "cm"))
    dev.off()
    cat("  Saved panel_S1.pdf (Spearman + Pearson)\n")

    # Export correlation matrices to CSV
    fwrite(as.data.table(cor_spearman, keep.rownames = "track"),
           file.path(OUT, "panel_S1_spearman_correlation.csv"))
    fwrite(as.data.table(cor_pearson, keep.rownames = "track"),
           file.path(OUT, "panel_S1_pearson_correlation.csv"))
    cat("  Saved panel_S1_spearman_correlation.csv and panel_S1_pearson_correlation.csv\n")

# --- |lfc| > 2 filtered analysis (includes shared peaks) ---
    cat("\n  Computing |log2FC| > 2 filtered correlations (all peaks including shared)...\n")
    sig_all <- bw_loci(unlist(bw_4), marker_gr_all, labels = names(bw_4), default_na = 0)
    sig_mat_all <- as.matrix(mcols(sig_all))

    cl0_signal_all <- sig_mat_all[, "cl0"]
    cl1_signal_all <- sig_mat_all[, "cl1"]
    lfc_all <- log2((cl1_signal_all + 0.01) / (cl0_signal_all + 0.01))
    lfc2_mask_all <- abs(lfc_all) > 2
    cat("    Peaks with |lfc| > 2:", sum(lfc2_mask_all), "of", nrow(sig_mat_all), "\n")

    if (sum(lfc2_mask_all) > 0) {
      sig_mat_lfc2 <- sig_mat_all[lfc2_mask_all, , drop = FALSE]
      cor_spearman_lfc2 <- cor(sig_mat_lfc2, method = "spearman", use = "complete.obs")
      cor_pearson_lfc2 <- cor(sig_mat_lfc2, method = "pearson", use = "complete.obs")

      cor_spearman_lfc2 <- cor_spearman_lfc2[order_G, order_G, drop = FALSE]
      cor_pearson_lfc2 <- cor_pearson_lfc2[order_G, order_G, drop = FALSE]
      dimnames(cor_spearman_lfc2) <- list(display_names, display_names)
      dimnames(cor_pearson_lfc2) <- list(display_names, display_names)

      heatmap_spearman_lfc2 <- Heatmap(cor_spearman_lfc2, name = "Spearman rho", col = col_fun_G,
                  column_title = sprintf("Spearman (|lfc|>2, all peaks, n=%d)", sum(lfc2_mask_all)), row_title = "",
              rect_gp = gpar(col = "black", lwd = 1),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", cor_spearman_lfc2[i, j]), x, y, gp = gpar(fontsize = 14))
              },
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_dend = FALSE, show_column_dend = FALSE,
              heatmap_width = unit(6, "cm"), heatmap_height = unit(6, "cm"),
              row_names_gp = gpar(fontsize = 12), column_names_gp = gpar(fontsize = 12))

      heatmap_pearson_lfc2 <- Heatmap(cor_pearson_lfc2, name = "Pearson r", col = col_fun_G,
                  column_title = sprintf("Pearson (|lfc|>2, all peaks, n=%d)", sum(lfc2_mask_all)), row_title = "",
              rect_gp = gpar(col = "black", lwd = 1),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", cor_pearson_lfc2[i, j]), x, y, gp = gpar(fontsize = 14))
              },
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_dend = FALSE, show_column_dend = FALSE,
              heatmap_width = unit(6, "cm"), heatmap_height = unit(6, "cm"),
              row_names_gp = gpar(fontsize = 12), column_names_gp = gpar(fontsize = 12))

      pdf(file.path(OUT, "panel_S1_lfc2.pdf"), width = 15, height = 7.5)
      draw(heatmap_spearman_lfc2 + heatmap_pearson_lfc2, ht_gap = unit(1, "cm"))
      dev.off()
      cat("  Saved panel_S1_lfc2.pdf (Spearman + Pearson, |lfc|>2)\n")

      fwrite(as.data.table(cor_spearman_lfc2, keep.rownames = "track"),
             file.path(OUT, "panel_S1_lfc2_spearman_correlation.csv"))
      fwrite(as.data.table(cor_pearson_lfc2, keep.rownames = "track"),
             file.path(OUT, "panel_S1_lfc2_pearson_correlation.csv"))
      cat("  Saved panel_S1_lfc2_spearman_correlation.csv and panel_S1_lfc2_pearson_correlation.csv\n")
    } else {
      cat("    No peaks with |lfc| > 2, skipping _lfc2 outputs\n")
    }
  }
}

# ===========================================================================
# Panel H: Coverage tracks (Lin28a and Cdhr3 ±10kb)
# ===========================================================================
if ("H" %in% run_panels) {
  cat("== Panel H: Tracks ==\n")
  expected_path <- FRAGMENTS_TSV
  frag_list <- Fragments(seurat[["peaks"]])
  old_path <- if (length(frag_list) > 0) frag_list[[1]]@path else ""
  if (length(frag_list) == 0 || !file.exists(old_path)) {
    if (!file.exists(expected_path)) {
      stop("Fragment file not found at '", expected_path, "'.\n",
           "This file is required for CoveragePlot (panel H).\n",
           "Download from GEO and place at:\n  ", expected_path,
           "\nThen index with: tabix -p bed ", expected_path)
    }
    message("Updating fragment path to: ", expected_path)
    for (i in seq_along(frag_list)) frag_list[[i]]@path <- expected_path
    seurat[["peaks"]]@fragments <- frag_list
  }
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations) <- "UCSC"
  Annotation(seurat) <- annotations
  genes_to_plot <- list(Lin28a = list(flank = 10000), Cdhr3 = list(flank = 10000))
  for (gene_name in names(genes_to_plot)) {
    gene_gr <- annotations[annotations$gene_name == gene_name]
    if (length(gene_gr) == 0) {
      message("Panel H: gene ", gene_name, " not found in annotation, skipping")
      next
    }
    gene_gr <- gene_gr[1]
    flank <- genes_to_plot[[gene_name]]$flank
    region_gr <- GRanges(seqnames = seqnames(gene_gr),
                         ranges = IRanges(start(gene_gr) - flank, end(gene_gr) + flank))
    p <- CoveragePlot(object = seurat, region = region_gr, annotation = TRUE,
                      peaks = TRUE, show.bulk = TRUE) +
      ggtitle(gene_name) + theme_fig1
    out_file <- file.path(OUT, paste0("panel_H_", gene_name, ".pdf"))
    ggsave(out_file, p, width = 12, height = 6, dpi = 300)
    message("  Saved ", out_file)
  }
}

cat(glue("\nAll Figure 1 outputs written to {OUT}/\n"))
