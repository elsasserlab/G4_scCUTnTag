# =============================================================
# Heatmap + correlation panel: published bulk 3T3 G4 CUT&Tag
# signal at scG4 MEF (Cluster 0) peaks, plus PQS score track.
#
# Outputs:
#   1. heatmap_clust0_published_G4_R.{pdf,png}
#      - 5-column heatmap: scG4 | bulk rep1-3 (red) | PQS (blue)
#      - rows sorted by scG4 signal
#   2. correlation_matrix_clust0.{pdf,png,csv}
#      - Spearman correlation between scG4, 3 bulk reps and PQS,
#        computed from per-peak mean signal in +/- 500 bp.
#
# Install once if missing:
#   BiocManager::install(c("EnrichedHeatmap", "rtracklayer",
#                          "circlize", "ComplexHeatmap"))
# =============================================================


source("rev/paths.R")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})
ENRICHEDHEATMAP_AVAILABLE <- requireNamespace("EnrichedHeatmap", quietly = TRUE) &&
                            requireNamespace("circlize", quietly = TRUE) &&
                            requireNamespace("ComplexHeatmap", quietly = TRUE)
if (ENRICHEDHEATMAP_AVAILABLE) {
  suppressPackageStartupMessages({
    library(EnrichedHeatmap)
    library(circlize)
    library(ComplexHeatmap)
  })
} else {
  message("EnrichedHeatmap/circlize/ComplexHeatmap not available — skipping heatmap plot.")
}

# ----- Paths -----
inp <- function(fn) {
  hits <- list.files(DATA, pattern = fn, recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("File not found in data/: ", fn)
  hits[1]
}
OUT_DIR  <- file.path(OUT_ROOT, "heatmap_published_G4")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PEAKS    <- CLUSTER_PEAKS_0
BIGWIGS  <- list(
  "scG4 MEF"   = BW_CL0,
  "Bulk rep 1" = file.path(GSE217860, "GSE217860_3T3-rep1_R1_val_1.bw"),
  "Bulk rep 2" = file.path(GSE217860, "GSM6729105_3T3-rep2_R1_val_1.bw"),
  "Bulk rep 3" = file.path(GSE217860, "GSM6729106_3T3-rep3_R1_val_1.bw"),
  "PQS score"  = file.path(DATA, "pqsfinder/PQS_scores.mm10.bw")
)

# Which tracks should use a blue (motif) palette vs red (signal)
BLUE_TRACKS <- c("PQS score")

# ----- Parameters -----
HALF_WIN  <- 2000   # +/- bp around peak center
BIN       <- 50     # bp per bin
CORE_WIN  <- 500    # +/- bp around center for per-peak signal aggregation
CANONICAL <- c(paste0("chr", 1:19), "chrX", "chrY")

if (!ENRICHEDHEATMAP_AVAILABLE) {
  cat("Skipping F13 heatmap and correlation: required packages not available.\n")
  quit(save = "no", status = 0)
}

# ----- Peaks: center points only, on canonical chroms -----
cat("Loading peaks...\n")
peaks <- import(PEAKS, format = "narrowPeak")
peaks <- peaks[as.character(seqnames(peaks)) %in% CANONICAL]
centers <- GRanges(seqnames(peaks),
                   IRanges(start = round((start(peaks) + end(peaks)) / 2),
                           width = 1))
cat(sprintf("  %d peaks (canonical chromosomes)\n", length(centers)))

# ----- Build normalized matrix per BigWig -----
cat("\nBuilding signal matrices (this is the slow step)...\n")
mats <- lapply(names(BIGWIGS), function(lbl) {
  cat(sprintf("  %s\n", lbl))
  bw  <- import(BIGWIGS[[lbl]], format = "BigWig")
  normalizeToMatrix(bw, centers,
                    value_column = "score",
                    extend       = HALF_WIN,
                    mean_mode    = "w0",
                    w            = BIN,
                    background   = 0,
                    smooth       = TRUE)
})
names(mats) <- names(BIGWIGS)

# ----- Sort peaks by mean scG4 signal (descending) -----
ord <- order(-rowMeans(mats[[1]]))
mats_sorted <- lapply(mats, function(m) m[ord, ])

# ----- Color palettes -----
make_col_red <- function(m) {
  q <- quantile(m, c(0.01, 0.99), na.rm = TRUE)
  colorRamp2(c(0, q[2]), c("white", "#cb181d"))
}
make_col_blue <- function(m) {
  q <- quantile(m, c(0.01, 0.99), na.rm = TRUE)
  colorRamp2(c(0, q[2]), c("white", "#2171b5"))
}

# ----- Build heatmap list -----
cat("\nDrawing heatmaps...\n")
ht_list <- NULL
for (i in seq_along(mats_sorted)) {
  lbl <- names(mats_sorted)[i]
  m   <- mats_sorted[[i]]
  col_fn <- if (lbl %in% BLUE_TRACKS) make_col_blue(m) else make_col_red(m)
  ht  <- EnrichedHeatmap(
           m,
           col           = col_fn,
           name          = lbl,
           column_title  = lbl,
           column_title_gp = gpar(fontsize = 11, fontface = "bold"),
           axis_name     = c(sprintf("-%d", HALF_WIN), "center",
                             sprintf("+%d", HALF_WIN)),
           axis_name_gp  = gpar(fontsize = 8),
           use_raster     = TRUE,
           raster_device  = if (requireNamespace("ragg", quietly = TRUE))
                              "agg_png" else "CairoPNG",
           raster_quality = 4,
           width         = unit(3.5, "cm"),
           top_annotation = HeatmapAnnotation(
             enriched = anno_enriched(
               gp = gpar(col = "black", lwd = 1.5),
               axis_param = list(side = "right")
             )
           )
         )
  ht_list <- if (is.null(ht_list)) ht else ht_list + ht
}

# ----- Save heatmap panel -----
pdf(file.path(OUT_DIR, "heatmap_clust0_published_G4_R.pdf"),
    width = 15, height = 9)
draw(ht_list,
     ht_gap            = unit(4, "mm"),
     row_title         = sprintf("scG4 MEF peaks (n = %d)", length(centers)),
     padding           = unit(c(4, 6, 4, 4), "mm"),
     newpage           = FALSE)
dev.off()

cat("  heatmap panel saved\n")

# =============================================================
# Correlation heatmap (per-peak signal, +/- 500 bp around center)
# =============================================================
cat("\nComputing per-peak signal & Spearman correlation matrix...\n")

# Identify the bin columns within +/- CORE_WIN of center.
# normalizeToMatrix returns a matrix with 2*HALF_WIN/BIN columns,
# center is between columns HALF_WIN/BIN and HALF_WIN/BIN+1.
n_bins  <- 2 * HALF_WIN / BIN                    # 80
n_side  <- CORE_WIN / BIN                        # 10
mid_lo  <- HALF_WIN / BIN - n_side + 1           # 31
mid_hi  <- HALF_WIN / BIN + n_side                # 50

# Per-peak mean signal in the central window for each track
peak_sig <- sapply(mats, function(m) rowMeans(m[, mid_lo:mid_hi, drop = FALSE]))
peak_sig_df <- as.data.frame(peak_sig)
colnames(peak_sig_df) <- names(BIGWIGS)

# Spearman correlation
cor_mat <- cor(peak_sig_df, method = "spearman", use = "pairwise.complete.obs")
write.csv(round(cor_mat, 3),
          file.path(OUT_DIR, "correlation_matrix_clust0.csv"),
          row.names = TRUE)

# (Correlation-matrix heatmap plot removed — not used in the slides.
#  The Spearman correlation table is still written as a supporting CSV.)

cat("\nCorrelation matrix (Spearman):\n")
print(round(cor_mat, 3))

cat("\nDone. Files in:", OUT_DIR, "\n")
cat("  - heatmap_clust0_published_G4_R.pdf (F13)\n")
cat("  - correlation_matrix_clust0.csv              (supporting table)\n")
