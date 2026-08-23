#!/usr/bin/env Rscript
# ===========================================================================
# Figure 4: Unsorted brain G4 mapping + label transfer validation
#
# Panels:
#   A - Unsorted G4 CnT UMAP (cl0=green, cl1=gray) + GFP+ UMAP
#   B - Unsorted → GFP+ projection (labeled + pred.score), Venn, boxplots
#   C - Signal profiles + heatmaps (cl0-only/both/cl1-only peaks × bigwigs × PQS)
#   D - enrichR GO Biological Process heatmap (cl0 vs cl1, CNS terms)
#   E - Genome browser tracks at Cc2da, Lias, Prkacb
#   F - scBridge neuron integration UMAP (copied)
#   G - Neuron prediction treemap (copied)
#
# Run from repo root:
#   Rscript fig4/fig4.R
#   Rscript fig4/fig4.R --only=A,B
#   Rscript fig4/fig4.R --skip=E
# ===========================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Gviz)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(stringr)
  library(eulerr)
  library(bedscout)
  library(ggpubr)
  library(RColorBrewer)
    library(EnsDb.Mmusculus.v79)
    library(rtracklayer)
    library(ChIPseeker)
  library(org.Mm.eg.db)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(GenomicRanges)
  library(wigglescout)
  library(enrichR)
  library(ComplexHeatmap)
  library(circlize)
})

ROOT    <- normalizePath(getwd())
RESULTS <- file.path(ROOT, "results")
DATA    <- file.path(ROOT, "github/data")
OUT     <- file.path(ROOT, "github/fig4", "outputs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

UNSORTED_RDS <- file.path(DATA, "GSE291468/GSM8836087_unsorted_Seurat_object.Rds")
SORTED_RDS   <- file.path(DATA, "GSE291468/GSM8836086_GFPpos_Seurat_object.Rds")
SC_DIR       <- file.path(RESULTS, "scBridge/output/unsorted_cl1_Zeisel")

save_fig <- function(p, name, width = 10, height = 8) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = width, height = height, device = "pdf")
}

read_peak_gr <- function(path) {
  dt <- fread(path, select = 1:3)
  GRanges(seqnames = dt[[1]], ranges = IRanges(start = dt[[2]], end = dt[[3]]))
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
cat("== Loading data ==\n")

cat("  Unsorted Seurat object...\n")
unsorted <- readRDS(UNSORTED_RDS)

cat("  Sorted Seurat object...\n")
sorted <- readRDS(SORTED_RDS)

# Fix fragment paths (Windows → Linux)
cat("  Fixing fragment paths...\n")
fix_fragments <- function(obj, correct_path) {
  frags <- Fragments(obj[["peaks"]])
  if (length(frags) > 0 && !file.exists(frags[[1]]@path)) {
    frag <- CreateFragmentObject(correct_path, cells = colnames(obj), validate = FALSE)
    slot(obj[["peaks"]], "fragments") <- list(frag)
  }
  obj
}
  unsorted <- fix_fragments(unsorted, file.path(DATA, "GSE291468/GSM8836087_unsorted_mouse_brain/CellRanger/fragments.tsv.gz"))
  sorted   <- fix_fragments(sorted,   file.path(DATA, "GSE291468/GSM8836086_GFP_sorted_mouse_brain/CellRanger/fragments.tsv.gz"))

# ===========================================================================
# Panel A: Unsorted UMAP + GFP+ UMAP
# ===========================================================================
if (run_panel("A")) {
  cat("== Panel A: Unsorted + GFP+ UMAPs ==\n")

  cl0_cl1_cols <- c("0" = "#addd8e", "1" = "#bcbcbc")

  p_A1 <- DimPlot(unsorted, reduction = "umap", group.by = "seurat_clusters",
                  label = TRUE, label.size = 4, repel = TRUE, pt.size = 0.35) +
    scale_color_manual(values = cl0_cl1_cols) +
    labs(title = "Unsorted G4 CnT", x = "UMAP_1", y = "UMAP_2") +
    theme_classic() + UMAP_THEME +
    guides(color = guide_legend(override.aes = list(size = 3)))

  sorted_cl_cols <- setNames(brewer.pal(4, "Set3")[1:4], c("0", "1", "2", "3"))
  p_A2 <- DimPlot(sorted, reduction = "umap", group.by = "seurat_clusters",
                  label = TRUE, label.size = 4, repel = TRUE, pt.size = 0.35) +
    scale_color_manual(values = sorted_cl_cols) +
    labs(title = "GFP+ G4 CnT", x = "UMAP_1", y = "UMAP_2") +
    theme_classic() + UMAP_THEME +
    guides(color = guide_legend(override.aes = list(size = 3)))

  p_A <- p_A1 | p_A2
  save_fig(p_A, "panel_A_umaps", 16, 7)
  cat("  Saved panel_A_umaps.pdf\n")
}

# ===========================================================================
# Panel B: Unsorted → GFP+ projection, Venn, boxplots
# ===========================================================================
if (run_panel("B")) {
  cat("== Panel B: Unsorted → GFP+ projection ==\n")

   # Reproduce the successful GEO testbed workflow exactly: the unsorted
   # object is the reference and the GFP+ object is projected onto it.
   DefaultAssay(unsorted) <- "GA"
   DefaultAssay(sorted) <- "GA"
   unsorted[["GA_unsorted"]] <- unsorted[["GA"]]
   DefaultAssay(unsorted) <- "GA_unsorted"
   unsorted[["GA"]] <- NULL
   sorted[["GA_sorted"]] <- sorted[["GA"]]
   DefaultAssay(sorted) <- "GA_sorted"
   sorted[["GA"]] <- NULL

   unsorted <- RunTFIDF(unsorted)
   unsorted <- FindTopFeatures(unsorted, min.cutoff = "q0")
   unsorted <- RunSVD(unsorted)
   unsorted <- RunUMAP(unsorted, dims = 1:10, reduction = "lsi",
                       return.model = TRUE, verbose = FALSE)

   cat("  Finding transfer anchors (GFP+ query → unsorted reference)...\n")
   common_ga <- intersect(rownames(sorted), rownames(unsorted))
   cat("  Common GA features:", length(common_ga), "\n")

   anchors <- FindTransferAnchors(
     reference = unsorted,
     query = sorted,
     reduction = "cca",
     query.assay = "GA_sorted",
      reference.assay = "GA_unsorted",
      k.filter = NA,
      features = common_ga
   )

   sorted_proj <- MapQuery(
     anchorset = anchors,
     query = sorted,
     reference = unsorted,
     refdata = list(seurat_clusters = "seurat_clusters"),
     reference.reduction = "umap",
     reduction.model = "umap"
   )

   # B1: Projected unsorted cells labeled by cl0/cl1. Keep the archived
   # plotting settings because these match the manuscript panel.
   cl0_cl1_cols <- c("0" = "#addd8e", "1" = "#bcbcbc")
    p_B1 <- DimPlot(sorted_proj, reduction = "umap",
                   group.by = "predicted.seurat_clusters",
                   label = TRUE, label.size = 3, repel = TRUE) +
                   scale_color_manual(values = cl0_cl1_cols) +
                   labs(title = "predicted unsorted labels on GFP+ UMAP",
                        x = "UMAP_1", y = "UMAP_2", fill = NULL) +
                   xlim(-5, 8) + ylim(-5, 5)
   save_fig(p_B1, "panel_B_umap1", 7, 6)

   # B2: Prediction score UMAP, using the archived gradient and limits.
    p_B2 <- FeaturePlot(sorted_proj, reduction = "umap",
                       features = "predicted.seurat_clusters.score",
                       label = FALSE, label.size = 3, repel = TRUE) +
                       scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
                       labs(title = "predicted unsorted labels on GFP+ UMAP",
                            x = "UMAP_1", y = "UMAP_2", fill = NULL) +
                       xlim(-5, 8) + ylim(-5, 5)
   save_fig(p_B2, "panel_B_umap2", 7, 6)

   # --- B3: Venn diagram (peak overlaps) ---
   cat("  Computing peak overlaps for Venn...\n")
   sorted_gr <- read_peak_gr(file.path(
     DATA, "GSE291468/GSM8836086_GFP_sorted_mouse_brain/CellRanger/peaks.bed"
   ))
   cl0_gr <- read_peak_gr(file.path(
     DATA, "GSE291468/GSM8836087_unsorted_cluster_0_peaks.narrowPeak"
   ))
   cl1_gr <- read_peak_gr(file.path(
     DATA, "GSE291468/GSM8836087_unsorted_cluster_1_peaks.narrowPeak"
   ))

   # bedscout::plot_euler() calculates all intersections with GRanges and then
   # fits eulerr with input = "union", matching the original manuscript Venn.
   venn_sets <- list(`GFP+` = sorted_gr, `unsorted cl. 0` = cl0_gr, `cl. 1` = cl1_gr)
   venn_counts <- bedscout:::calculate_venn_intersections(
     venn_sets, names = names(venn_sets), ignore.strand = TRUE
   )
   write.csv(data.frame(region = names(venn_counts), count = unname(venn_counts)),
             file.path(OUT, "panel_B_venn_counts.csv"), row.names = FALSE)
   euler_obj <- eulerr::euler(venn_counts, input = "union")
   p_B3 <- plot(euler_obj,
                fills = c("#ffffff", "#addd8e", "#bcbcbc"),
                quantities = TRUE,
                main = "")
   save_fig(wrap_elements(p_B3), "panel_B_venn", 7, 6)

  # --- B4: Boxplots of prediction score ---
  cat("  Creating prediction score boxplots...\n")
   meta <- sorted_proj@meta.data
  meta$cell <- rownames(meta)
  meta$pred_score <- meta$predicted.seurat_clusters.score
   meta$pred_cluster <- meta$predicted.seurat_clusters
   meta$orig_cluster <- meta$seurat_clusters

  meta$pred_cluster <- factor(meta$pred_cluster, levels = c("0", "1"))
   meta$orig_cluster <- factor(meta$orig_cluster, levels = c("0", "1", "2", "3"))

   p_B4 <- ggplot(meta,
                  aes(x = seurat_clusters,
                      y = predicted.seurat_clusters.score,
                      fill = predicted.seurat_clusters)) +
     geom_boxplot() +
     scale_fill_manual(values = c("#addd8e", "#bcbcbc")) +
     ylim(0.5, 1) +
     labs(title = "", x = "Seurat cluster (GFP+)", y = "prediction score",
          fill = "predicted \nunsorted cluster") +
     theme_minimal() +
     theme(
       text = element_text(size = 10),
       plot.title = element_text(size = 10),
       axis.title.y = element_text(size = 20, color = "black"),
       axis.title.x = element_text(size = 20, color = "black"),
       axis.text.x = element_text(size = 20, color = "black", angle = 45,
                                  vjust = 0.75, hjust = 0.5),
       axis.text.y = element_text(size = 20, color = "black")
     ) +
     stat_compare_means(label.y = 1, label = "p.signif")
   save_fig(p_B4, "panel_B_predictionscore", 7, 6)

  # Compose B panel
  p_B_top <- p_B1 | p_B2
  p_B_bot <- wrap_plots(list(
    wrap_elements(p_B3 + ggtitle("") + theme(plot.title = element_text(size = 1))),
    p_B4
  ), widths = c(1, 1.5))
   p_B <- p_B_top / p_B_bot + plot_layout(heights = c(1, 1))
  save_fig(p_B, "panel_B_projection_venn_boxplots", 16, 14)
  cat("  Saved panel_B_projection_venn_boxplots.pdf\n")
}

# ===========================================================================
# Panel C: Signal profiles and heatmaps for cl0-only / both / cl1-only peaks
# ===========================================================================
if (run_panel("C")) {
  cat("== Panel C: Signal profiles + heatmaps (cl0-only/both/cl1-only) ==\n")

  bw_cl0 <- file.path(DATA, "GSE291468/GSM8836087_cluster0_RPGC.bw")
  bw_cl1 <- file.path(DATA, "GSE291468/GSM8836087_cluster1_RPGC.bw")
  bw_gfp <- file.path(DATA, "GSE291468/GSM8836086_GFP_sorted_mousebrain.rpgc.bw")
  bw_pqs <- file.path(DATA, "pqsfinder/PQS_scores.mm10.bw")

   cl0_np <- file.path(DATA, "GSE291468/GSM8836087_unsorted_cluster_0_peaks.narrowPeak")
   cl1_np <- file.path(DATA, "GSE291468/GSM8836087_unsorted_cluster_1_peaks.narrowPeak")

  if (!all(file.exists(bw_cl0, bw_cl1, bw_gfp, bw_pqs, cl0_np, cl1_np))) {
    missing <- c(bw_cl0, bw_cl1, bw_gfp, bw_pqs, cl0_np, cl1_np)[
      !file.exists(c(bw_cl0, bw_cl1, bw_gfp, bw_pqs, cl0_np, cl1_np))
    ]
    message("  Missing files, skipping Panel C: ", paste(basename(missing), collapse = ", "))
  } else {
    # Match Figure 1F exactly: classify fixed summit-centered +/-500 bp
    # windows rather than fragmented overlaps of the raw peak intervals.
    cl0_np_gr <- rtracklayer::import(cl0_np)
    cl1_np_gr <- rtracklayer::import(cl1_np)
    make_summit_gr <- function(np) {
      summit_pos <- start(np) + np$peak
      GRanges(
        seqnames = seqnames(np),
        ranges = IRanges(start = summit_pos - 500, end = summit_pos + 500)
      )
    }
    cl0_gr <- make_summit_gr(cl0_np_gr)
    cl1_gr <- make_summit_gr(cl1_np_gr)
    hits <- findOverlaps(cl0_gr, cl1_gr)
    shared0 <- unique(queryHits(hits))
    shared1 <- unique(subjectHits(hits))
    cl0_only <- cl0_gr[-shared0]
    cl1_only <- cl1_gr[-shared1]
    shared <- GenomicRanges::reduce(sort(c(cl0_gr[shared0], cl1_gr[shared1])))

    cat("  cl0-only:", length(cl0_only), " | both:", length(shared), " | cl1-only:", length(cl1_only), "\n")

    # Write temporary BED files for wigglescout
    dir.create(file.path(OUT, "tmp_bed"), showWarnings = FALSE, recursive = TRUE)
    write.bed <- function(gr, name) {
      df <- data.frame(seqnames = as.character(seqnames(gr)),
                       start = start(gr), end = end(gr))
      path <- file.path(OUT, "tmp_bed", paste0(name, ".bed"))
      fwrite(df, path, sep = "\t", col.names = FALSE)
      path
    }
    bed_cl0 <- write.bed(cl0_only, "cl0_only")
    bed_both <- write.bed(shared, "both")
    bed_cl1 <- write.bed(cl1_only, "cl1_only")

    beds <- list(cl0_only = bed_cl0, both = bed_both, cl1_only = bed_cl1)
    bws_signal <- c("cluster0" = bw_cl0, "cluster1" = bw_cl1, "GFP+" = bw_gfp)
    bw_pqs_vec <- c("PQS score" = bw_pqs)

    # --- Figure 1F-style composite profile and heatmaps ---
    cat("  Generating composite profiles and heatmaps...\n")
    profile_colors <- c(cl0_only = "#555555", both = "black", cl1_only = "#4daf4a")
    category_heights <- pmax(1, vapply(beds, function(x) {
      nrow(fread(x, header = FALSE))
    }, numeric(1)) / 10000)
    sort_tracks <- c(cl0 = bw_cl0, cl1 = bw_cl1)
    heatmap_orders <- lapply(names(beds), function(bed_name) {
      gr <- rtracklayer::import(beds[[bed_name]])
      if (!length(gr)) return(integer())
      cluster_signal <- as.data.frame(mcols(bw_loci(
        sort_tracks, gr, labels = names(sort_tracks), default_na = 0
      )))
      score <- switch(bed_name,
        cl0_only = cluster_signal$cl0,
        cl1_only = cluster_signal$cl1,
        both = cluster_signal$cl0 * cluster_signal$cl1
      )
      order(score)
    })
    names(heatmap_orders) <- names(beds)
    compute_zmax <- function(bw_path, gr_list, upstream = 3000, downstream = 3000) {
      bw <- BigWigFile(bw_path)
      all_vals <- c()
      for (gr in gr_list) {
        centers <- start(gr) + (end(gr) - start(gr)) %/% 2
        query_gr <- GRanges(
          seqnames = seqnames(gr),
          ranges = IRanges(start = centers - upstream, end = centers + downstream)
        )
        sig <- import(bw, selection = query_gr, as = "NumericList")
        vals <- unlist(lapply(sig, function(x) x[is.finite(x)]))
        if (length(vals)) all_vals <- c(all_vals, vals)
      }
      if (!length(all_vals)) return(10)
      as.numeric(quantile(all_vals, 0.99, na.rm = TRUE))
    }
    all_cat_grs <- lapply(beds, rtracklayer::import)
    composite_height <- 3 * max(8, 3 + sum(category_heights))

    make_composite <- function(track, track_name, cmap, zmin, zmax) {
      profile <- plot_bw_profile(
        track, beds, mode = "center", upstream = 3000, downstream = 3000,
        colors = profile_colors, show_error = TRUE, verbose = FALSE
      ) + labs(title = paste0("Signal: ", track_name)) +
        theme(plot.title = element_text(size = 12, face = "bold"))
      # Match Figure 1F: 82.5% width, left aligned, with the compact profile
      # occupying the upper slot of the enlarged composite canvas.
      profile <- cowplot::ggdraw() + cowplot::draw_plot(
        profile, x = 0, y = 0.1875, width = 0.825, height = 0.625
      )
      heatmaps <- lapply(names(beds), function(bed_name) {
        plot_bw_heatmap(
          track, beds[[bed_name]], mode = "center", upstream = 3000,
          downstream = 3000, cmap = cmap, zmin = zmin, zmax = zmax,
          max_rows_allowed = max(1, ceiling(nrow(fread(beds[[bed_name]], header = FALSE)) / 100)),
          order_by = heatmap_orders[[bed_name]],
          verbose = FALSE
        ) + ggtitle(bed_name)
      })
      panel <- ggarrange(
        profile,
        ggarrange(plotlist = heatmaps, ncol = 1, heights = category_heights, align = "v"),
        ncol = 1,
        heights = c(6, max(1, 3 * max(8, 3 + sum(category_heights)) - 6)),
        align = "v"
      )
      file <- file.path(OUT, paste0("panel_C_", make.names(track_name), ".pdf"))
      ggsave(file, panel, width = 6, height = composite_height, limitsize = FALSE)
      cat("  Saved", basename(file), "\n")
    }

    # wigglescout expects palette names; these provide the requested white-to-
    # dark-gray, white-to-green, and white-to-purple gradients.
    make_composite(
      bws_signal[["cluster0"]], "cluster0", "Greens", 0,
      compute_zmax(bws_signal[["cluster0"]], all_cat_grs)
    )
    make_composite(
      bws_signal[["cluster1"]], "cluster1", "Greys", 0,
      compute_zmax(bws_signal[["cluster1"]], all_cat_grs)
    )
    make_composite(
      bws_signal[["GFP+"]], "GFP+", "Greens", 0,
      compute_zmax(bws_signal[["GFP+"]], all_cat_grs)
    )
    make_composite(bw_pqs_vec, "PQS_score", "PuRd", 5, 20)

    # Clean up temp BEDs
    unlink(file.path(OUT, "tmp_bed"), recursive = TRUE)

    cat("  Saved Figure 1F-style Panel C composites\n")
  }
}

# ===========================================================================
# Panel D: enrichR GO Biological Process heatmap (unsorted cl0 vs cl1)
# ===========================================================================
if (run_panel("D")) {
  cat("== Panel D: enrichR GO term heatmap ==\n")

   cl0_np <- file.path(DATA, "GSE291468/GSM8836087_unsorted_cluster_0_peaks.narrowPeak")
   cl1_np <- file.path(DATA, "GSE291468/GSM8836087_unsorted_cluster_1_peaks.narrowPeak")

  if (!all(file.exists(cl0_np, cl1_np))) {
    message("  Missing narrowPeak files, skipping Panel D")
  } else {
    # Recreate the archived bedtools unique sets directly from GEO files: each
    # cluster peak must be absent from both the other unsorted cluster and the
    # GFP+ peak universe.
    gfp_gr <- read_peak_gr(file.path(
      DATA, "GSE291468/GSM8836086_GFP_sorted_mouse_brain/CellRanger/peaks.bed"
    ))
    cl0_raw <- fread(cl0_np)
    cl1_raw <- fread(cl1_np)
    cl0_raw_gr <- GRanges(seqnames = cl0_raw$V1,
                          ranges = IRanges(start = cl0_raw$V2 + 1L, end = cl0_raw$V3))
    cl1_raw_gr <- GRanges(seqnames = cl1_raw$V1,
                          ranges = IRanges(start = cl1_raw$V2 + 1L, end = cl1_raw$V3))
    cl0_dt <- cl0_raw[countOverlaps(cl0_raw_gr, c(gfp_gr, cl1_raw_gr)) == 0]
    cl1_dt <- cl1_raw[countOverlaps(cl1_raw_gr, c(gfp_gr, cl0_raw_gr)) == 0]

    annotate_to_nearest <- function(dt) {
      gr <- GRanges(seqnames = dt$V1, ranges = IRanges(start = dt$V2, end = dt$V3))
      as.data.table(as.data.frame(annotatePeak(
        gr,
        tssRegion = c(-3000, 3000),
        TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
        annoDb = "org.Mm.eg.db"
      )))
    }
    cat("  Deriving unique peaks and annotating nearest genes...\n")
    cl0_dt <- cl0_dt[V7 > quantile(V7, 0.90)]
    cl1_dt <- cl1_dt[V7 > quantile(V7, 0.90)]
    cl0_dt <- annotate_to_nearest(cl0_dt)
    cl0_dt <- cl0_dt[abs(distanceToTSS) < 3000 & !is.na(SYMBOL)]
    cl1_dt <- annotate_to_nearest(cl1_dt)
    cl1_dt <- cl1_dt[abs(distanceToTSS) < 3000 & !is.na(SYMBOL)]

    cat("  cl0 genes:", length(unique(cl0_dt$SYMBOL)),
        " | cl1 genes:", length(unique(cl1_dt$SYMBOL)), "\n")

    # Run enrichR
    cat("  Running enrichR (GO_Biological_Process_2023)...\n")
     dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023",
              "GO_Biological_Process_2023")

     enrich_cl0 <- enrichr(unique(cl0_dt$SYMBOL), dbs)
     enrich_cl1 <- enrichr(unique(cl1_dt$SYMBOL), dbs)

     bp_cl0 <- as.data.table(enrich_cl0[["GO_Biological_Process_2023"]])
     bp_cl1 <- as.data.table(enrich_cl1[["GO_Biological_Process_2023"]])

     bp_cl0 <- bp_cl0[P.value < 0.05][, .(Term, P.value)]
     bp_cl1 <- bp_cl1[P.value < 0.05][, .(Term, P.value)]

    bp_cl0$type <- "unique_cl0_terms"
    bp_cl1$type <- "unique_cl1_terms"

    # Combine and filter to CNS-related terms
    terms <- rbind(bp_cl0, bp_cl1)
    terms[, Term := str_to_lower(Term)]
    cns_pattern <- "nervous|neural|oligodendrocyte|microglia|neuron|axon|synaptic|synapsis|myelination|glial|glia"
     terms <- terms[str_detect(Term, cns_pattern)]

     if (nrow(terms) == 0) {
       stop("enrichR returned no CNS-related GO terms after the archived filtering workflow")
     }

    # Pivot to matrix
    terms_wide <- dcast(terms, Term ~ type, value.var = "P.value", fill = 1)
    terms_mat <- as.matrix(terms_wide[, -1, with = FALSE])
    rownames(terms_mat) <- terms_wide$Term

    # Rename columns for display
    colnames(terms_mat) <- gsub("unique_cl0_terms", "Unsorted cl0", colnames(terms_mat))
    colnames(terms_mat) <- gsub("unique_cl1_terms", "Unsorted cl1", colnames(terms_mat))

    cat("  CNS-related GO terms:", nrow(terms_mat), "\n")

    # Match the archived heatmap scale and layout.
    col_fun <- colorRamp2(c(1, 0.05, 0.025, 0.001, 0.002),
                          c("grey", "#fee0d2", "#fc9272", "#de2d26", "#99000d"))

    hm <- Heatmap(
      terms_mat,
      column_title = "enrichR",
      row_title = "",
      name = "p-value",
      col = col_fun,
      heatmap_legend_param = list(
        title = "p-value",
        at = c(1, 0.05, 0.001),
        labels = c("not enriched", "p < 0.05", "p < 0.001"),
        legend_height = unit(2, "cm")
      ),
      rect_gp = gpar(col = "black", lwd = 0.1),
      show_column_dend = FALSE,
      cluster_columns = FALSE,
      cluster_rows = TRUE,
      show_row_dend = FALSE,
      heatmap_width = unit(7, "cm"),
      heatmap_height = unit(10, "cm"),
      row_names_gp = gpar(fontsize = 7),
      column_names_gp = gpar(fontsize = 6),
      column_names_rot = 90
    )

    pdf(file.path(OUT, "panel_D_enrichR_heatmap.pdf"), width = 5, height = max(5, nrow(terms_mat) * 0.15 + 1))
    print(hm)
    dev.off()

    # Also save the GO terms table
    fwrite(terms, file.path(OUT, "panel_D_enrichR_GO_terms.csv"))
    cat("  Saved panel_D_enrichR_heatmap.pdf\n")
  }
}

# ===========================================================================
# Panel E: Genome browser tracks
# ===========================================================================
if (run_panel("E")) {
  cat("== Panel E: Genome browser tracks ==\n")

    cl1_bed <- file.path(DATA, "GSE291468/GSM8836087_unsorted_cluster_1_peaks.narrowPeak")
   fragment_file <- file.path(DATA, "GSE291468/GSM8836087_unsorted_mouse_brain/CellRanger/fragments.tsv.gz")

    if (!file.exists(cl1_bed)) {
      message("  Missing unsorted cluster peaks, skipping Panel E")
    } else {
       gene_names <- c("Cc2d1a", "Lias", "Prkacb")
       edb <- EnsDb.Mmusculus.v79
        DefaultAssay(unsorted) <- "peaks"
       # BigWig-only panel: prevent CoveragePlot from consulting the loaded
       # fragment object or its stale Tabix index.
       slot(unsorted[["peaks"]], "fragments") <- list()
      rpgc_tracks <- list(
        `cluster 0 RPGC` = file.path(DATA, "GSE291468/GSM8836087_cluster0_RPGC.bw"),
        `cluster 1 RPGC` = file.path(DATA, "GSE291468/GSM8836087_cluster1_RPGC.bw")
      )
       rpgc_tracks <- rpgc_tracks[vapply(rpgc_tracks, file.exists, logical(1))]
       tss_roi <- function(gene_name, flank = 5000) {
         gene_gr <- genes(edb)
         gene_gr <- gene_gr[gene_gr$gene_name == gene_name]
         if (!length(gene_gr)) return(NULL)
         gene_gr <- gene_gr[1]
         tss <- if (as.character(strand(gene_gr)) == "-") end(gene_gr) else start(gene_gr)
          seqname <- as.character(seqnames(gene_gr))
          if (!startsWith(seqname, "chr")) seqname <- paste0("chr", seqname)
          paste0(seqname, "-", max(1, tss - flank), "-", tss + flank)
       }
        for (g in gene_names) {
       cat("    Gene:", g, "\n")
       gene_info <- genes(edb, filter = ~ gene_name == g)
      if (length(gene_info) == 0) {
        warning("  Gene ", g, " not found in EnsDb, skipping")
        next
        }
         out_pdf <- file.path(OUT, paste0("panel_E_genome_track_", g, ".pdf"))
         roi <- tss_roi(g)
         if (is.null(roi)) {
           warning("TSS not found for ", g, "; skipping")
           next
         }
         # Plot only the precomputed RPGC BigWigs, without fragment coverage.
          p <- CoveragePlot(unsorted, region = roi, annotation = TRUE,
                           peaks = FALSE, show.bulk = FALSE, bigwig = rpgc_tracks,
                          bigwig.type = "coverage", bigwig.scale = "common") +
         scale_fill_manual(values = c("0" = "#addd8e", "1" = "#bcbcbc"), drop = FALSE) +
          ggtitle(paste0(g, " promoter +/- 5 kb; cluster-specific G4"))
       ggsave(out_pdf, p, width = 10, height = 6, device = "pdf")
     }
    cat("  Saved panel_E_genome_track_*.pdf\n")
  }
}

# ===========================================================================
# Panel F: scBridge neuron integration UMAP
# ===========================================================================
if (run_panel("F")) {
  cat("== Panel F: scBridge neuron UMAP ==\n")
  src <- file.path(SC_DIR, "Seurat_UMAPs-unsorted_cl1_int.pdf")
  if (file.exists(src)) {
    file.copy(src, file.path(OUT, "panel_F_neuron_integration.pdf"), overwrite = TRUE)
    cat("  Copied Seurat_UMAPs-unsorted_cl1_int.pdf\n")
  } else {
    message("  Not found: ", src)
  }
}

# ===========================================================================
# Panel G: Neuron prediction treemap
# ===========================================================================
if (run_panel("G")) {
  cat("== Panel G: Neuron treemap ==\n")
  src <- file.path(SC_DIR, "tree_plot-unsorted_cl1_predictions.pdf")
  if (file.exists(src)) {
    file.copy(src, file.path(OUT, "panel_G_neuron_treemap.pdf"), overwrite = TRUE)
    cat("  Copied tree_plot-unsorted_cl1_predictions.pdf\n")
  } else {
    message("  Not found: ", src)
  }
}

# ===========================================================================
# Panel S1: QC violin plots for unsorted mouse brain
# ===========================================================================
if (run_panel("S1")) {
  cat("== Panel S1: QC violin plots ==\n")
  
  # Load unsorted data
  SEURAT_UNSORTED <- file.path(DATA, "GSE291468/GSM8836087_unsorted_Seurat_object.Rds")
  FRAGMENTS_UNSORTED <- file.path(DATA, "GSE291468/GSM8836087_unsorted_mouse_brain/CellRanger/fragments.tsv.gz")
  
  if (!file.exists(SEURAT_UNSORTED) || !file.exists(FRAGMENTS_UNSORTED)) {
    stop("Required files not found:\n  ", SEURAT_UNSORTED, "\n  ", FRAGMENTS_UNSORTED)
  }
  
  seurat <- readRDS(SEURAT_UNSORTED)
  cat(sprintf("  Loaded Seurat object: %s cells, %s features\n", ncol(seurat), nrow(seurat)))
  
  # Compute QC metrics
  cat("  Computing QC metrics...\n")
  fragments <- Signac::CreateFragmentObject(FRAGMENTS_UNSORTED, cells = colnames(seurat))
  seurat[["peaks"]]@fragments <- list(fragments)
  
  if (!"TSS_enrichment" %in% colnames(seurat@meta.data)) {
    seurat <- TSSEnrichment(seurat, fast = TRUE)
    if ("TSS.enrichment" %in% colnames(seurat@meta.data)) {
      seurat$TSS_enrichment <- seurat$TSS.enrichment
    }
  }
  if (!"nucleosome_signal" %in% colnames(seurat@meta.data)) {
    seurat <- NucleosomeSignal(seurat)
  }
  if (!"FRiP" %in% colnames(seurat@meta.data)) {
    total_frags <- CountFragments(FRAGMENTS_UNSORTED, cells = colnames(seurat), verbose = FALSE)
    seurat <- AddMetaData(seurat, metadata = total_frags$reads_count, col.name = "total_fragments")
    seurat <- FRiP(seurat, assay = "peaks", total.fragments = "total_fragments", col.name = "FRiP", verbose = FALSE)
  }
  
  # Generate cluster colors
  n_clusters <- length(levels(seurat$seurat_clusters))
  colors <- c("#9ecae1", "#fc9272", "#a1d99b", "#fdb462", "#b3de69")[1:n_clusters]
  
  # Create QC plots
  make_vln_s1 <- function(feature, ylab_text, log10 = FALSE, ylim = NULL) {
    p <- VlnPlot(seurat, group.by = "seurat_clusters", features = feature, pt.size = 0.1) +
      scale_fill_manual(values = colors) +
      ggtitle("") + xlab("cluster") + ylab(ylab_text) +
      stat_summary(fun = median, geom = "crossbar", width = 0.4, linewidth = 0.4, color = "black") +
      stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
      theme_classic() +
      theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 12)) +
      NoLegend()
    if (log10) {
      p <- p + scale_y_log10()
    } else if (!is.null(ylim)) {
      p <- p + coord_cartesian(ylim = ylim)
    }
    p
  }
  
  plots <- list(
    make_vln_s1("nFeature_peaks", "nFeature (peaks)", log10 = TRUE),
    make_vln_s1("nCount_peaks", "nCount (peaks)", log10 = TRUE),
    make_vln_s1("TSS_fragments", "TSS fragments", log10 = TRUE),
    make_vln_s1("mitochondrial", "mitochondrial fragments", log10 = TRUE),
    make_vln_s1("FRiP", "Fraction of reads in peaks", ylim = c(0, 0.20)),
    make_vln_s1("TSS_enrichment", "TSS enrichment", ylim = c(0, 10)),
    make_vln_s1("nucleosome_signal", "Nucleosome signal", ylim = c(0, 3))
  )
  
  p_combined <- cowplot::plot_grid(plotlist = plots, ncol = 2, nrow = 4)
  ggsave(file.path(OUT, "panel_S1_qc_unsorted.pdf"), p_combined, width = 12, height = 16, dpi = 300)
  cat("  Saved panel_S1_qc_unsorted.pdf\n")
}

cat("\nAll Figure 4 outputs written to", OUT, "/\n")
