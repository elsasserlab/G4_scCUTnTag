#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Genomic feature distribution for scG4 peak categories:
# cl0_only, cl1_only, shared (both)
#
# Pipeline:
#   1. Load cluster-specific peak files and define categories
#   2. Build feature interval sets from GENCODE mm10 GTF
#   3. Classify each peak by midpoint using priority cascade
#   4. Plot grouped + stacked bars per category
# ---------------------------------------------------------------------------

source("rev/paths.R")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ===== CONFIG =====
inp <- function(fn) {
  hits <- list.files(DATA, pattern = fn, recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("File not found in data/: ", fn)
  hits[1]
}
GTF_PATH <- file.path(GENOME_DIR, "mm10_.annotation.gtf.gz")
OUT_DIR  <- file.path(getwd(), "fig1", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CANONICAL   <- c(paste0("chr", 1:19), "chrX", "chrY")
PROMOTER_UP <- 2000
PROMOTER_DN <- 2000

# ===== 1. LOAD SCG4 PEAKS AND DEFINE CATEGORIES =====
cat("Loading scG4 cluster peaks...\n")
cl0_np <- rtracklayer::import(CLUSTER_PEAKS_0, format = "narrowPeak")
cl1_np <- rtracklayer::import(CLUSTER_PEAKS_1, format = "narrowPeak")
cl0_np <- cl0_np[as.character(seqnames(cl0_np)) %in% CANONICAL]
cl1_np <- cl1_np[as.character(seqnames(cl1_np)) %in% CANONICAL]

# Create summit-centered regions (±250 bp)
make_summit_gr <- function(np) {
  summit_pos <- start(np) + np$peak
  GRanges(seqnames = seqnames(np),
          ranges = IRanges(start = summit_pos - 250, end = summit_pos + 250))
}
gr0 <- make_summit_gr(cl0_np)
gr1 <- make_summit_gr(cl1_np)

# Define categories by overlap
hits <- findOverlaps(gr0, gr1)
shared0 <- unique(queryHits(hits))
shared1 <- unique(subjectHits(hits))
cl0_only <- gr0[-shared0]
cl1_only <- gr1[-shared1]
both <- GenomicRanges::reduce(sort(c(gr0[shared0], gr1[shared1])))

cat(sprintf("  cl0_only: %d peaks\n", length(cl0_only)))
cat(sprintf("  cl1_only: %d peaks\n", length(cl1_only)))
cat(sprintf("  shared:   %d peaks\n", length(both)))

peak_sets <- list(
  "cl0_only (MEF)" = cl0_only,
  "cl1_only (mESC)" = cl1_only,
  "shared (both)" = both
)

# ===== 2. BUILD FEATURE INTERVAL SETS FROM GTF =====
cat("\nBuilding feature intervals from GTF...\n")
gtf <- rtracklayer::import(GTF_PATH)
gtf <- gtf[as.character(seqnames(gtf)) %in% CANONICAL]
seqlevels(gtf) <- intersect(seqlevels(gtf), CANONICAL)

# 2a. Promoters: TSS ± 2 kb on transcripts
tx <- gtf[gtf$type == "transcript"]
promoters_gr <- reduce(
  promoters(tx, upstream = PROMOTER_UP, downstream = PROMOTER_DN),
  ignore.strand = TRUE
)

# 2b. UTRs — GENCODE combines them; split by strand + start_codon position
utr     <- gtf[gtf$type == "UTR"]
start_c <- gtf[gtf$type == "start_codon"]

sc_dt <- as.data.table(start_c)[, .(transcript_id, sc_start = start,
                                    sc_end = end, sc_strand = strand)]
sc_dt <- sc_dt[, .(sc_start = min(sc_start),
                   sc_end   = max(sc_end),
                   sc_strand = data.table::first(sc_strand)),
               by = transcript_id]
utr_dt <- as.data.table(utr)[, .(seqnames, start, end, strand, transcript_id)]
utr_dt <- merge(utr_dt, sc_dt, by = "transcript_id")

utr_dt[, type := fcase(
  (strand == "+" & end   <= sc_start), "5'UTR",
  (strand == "-" & start >= sc_end),   "5'UTR",
  (strand == "+" & start >= sc_end),   "3'UTR",
  (strand == "-" & end   <= sc_start), "3'UTR",
  default = NA_character_
)]
utr_dt <- utr_dt[!is.na(type)]

utr5_gr <- reduce(GRanges(utr_dt[type == "5'UTR"]), ignore.strand = TRUE)
utr3_gr <- reduce(GRanges(utr_dt[type == "3'UTR"]), ignore.strand = TRUE)

# 2c. CDS exons = all exons minus UTRs
exons_all      <- reduce(gtf[gtf$type == "exon"], ignore.strand = TRUE)
exon_nonutr_gr <- GenomicRanges::setdiff(
  exons_all, reduce(c(utr5_gr, utr3_gr), ignore.strand = TRUE),
  ignore.strand = TRUE
)

# 2d. Introns = transcript span minus exons
tx_span   <- reduce(tx, ignore.strand = TRUE)
intron_gr <- GenomicRanges::setdiff(tx_span, exons_all, ignore.strand = TRUE)

# 2e. Stash in a named list (one GRanges per priority class)
feature_list <- list(
  Promoter     = promoters_gr,
  `5' UTR`     = utr5_gr,
  `3' UTR`     = utr3_gr,
  `Exon (CDS)` = exon_nonutr_gr,
  Intron       = intron_gr
)
cat("  feature interval counts:\n")
for (nm in names(feature_list))
  cat(sprintf("    %-12s %d intervals (%.1f Mb)\n",
              nm, length(feature_list[[nm]]),
              sum(width(feature_list[[nm]])) / 1e6))

# ===== 3. CLASSIFY PEAKS (priority cascade, single category per peak) =====
CATS <- c("Promoter", "5' UTR", "3' UTR", "Exon (CDS)", "Intron", "Intergenic")

classify_priority <- function(peaks_gr, feature_list) {
  mids <- GRanges(
    seqnames(peaks_gr),
    IRanges(start = (start(peaks_gr) + end(peaks_gr)) %/% 2, width = 1)
  )
  cat_assign <- rep("Intergenic", length(peaks_gr))
  for (nm in c("Promoter", "5' UTR", "3' UTR", "Exon (CDS)", "Intron")) {
    unassigned <- cat_assign == "Intergenic"
    if (!any(unassigned)) break
    hits <- overlapsAny(mids[unassigned], feature_list[[nm]],
                        ignore.strand = TRUE)
    cat_assign[unassigned][hits] <- nm
  }
  factor(cat_assign, levels = CATS)
}

cat("\nClassifying peaks (priority cascade by midpoint)...\n")

distros <- list()
for (label in names(peak_sets)) {
  gr  <- peak_sets[[label]]
  cls <- classify_priority(gr, feature_list)
  tab <- table(cls)
  distros[[label]] <- as.numeric(tab)[match(CATS, names(tab))]
  distros[[label]][is.na(distros[[label]])] <- 0
  total <- sum(distros[[label]])
  cat(sprintf("  %s (n = %d):\n", label, total))
  for (i in seq_along(CATS))
    cat(sprintf("    %-12s %6d  (%5.1f%%)\n",
                CATS[i], distros[[label]][i], 100 * distros[[label]][i] / total))
}

# Long table for plotting / CSV
dist_dt <- data.table(
  category = rep(names(peak_sets), each = length(CATS)),
  feature  = rep(CATS, times = length(peak_sets)),
  count    = unlist(distros)
)
dist_dt[, pct := 100 * count / sum(count), by = category]
fwrite(dist_dt, file.path(OUT_DIR, "panel_D_feature_distribution_scG4.csv"), sep = "\t")
cat("\nSaved panel_D_feature_distribution_scG4.csv\n")

# ===== 4. PLOTS =====
cat("\nDrawing figures...\n")
dist_dt$feature <- factor(dist_dt$feature, levels = CATS)
dist_dt$category <- factor(dist_dt$category,
                           levels = c("cl0_only (MEF)", "cl1_only (mESC)", "shared (both)"))

# 4a. Grouped bars (side by side per category)
p_grouped <- ggplot(dist_dt, aes(x = feature, y = pct, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", pct)),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("cl0_only (MEF)" = "#9ecae1",
                               "cl1_only (mESC)" = "#fc9272",
                               "shared (both)" = "#bdbdbd")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Fraction of peaks (%)",
       title = "Genomic feature distribution of scG4 peaks",
       subtitle = sprintf("cl0_only: %s | cl1_only: %s | shared: %s",
                          format(length(cl0_only), big.mark = ","),
                          format(length(cl1_only), big.mark = ","),
                          format(length(both), big.mark = ",")),
       fill = NULL) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        plot.title  = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        legend.position = "top",
        legend.text = element_text(size = 11))

ggsave(file.path(OUT_DIR, "panel_D_feature_distribution_scG4.pdf"),
       p_grouped, width = 12, height = 6)
cat("  Saved panel_D_feature_distribution_scG4.pdf\n")

# 4b. Stacked bars (composition view)
p_stacked <- ggplot(dist_dt, aes(x = category, y = pct, fill = feature)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.4, width = 0.65) +
  geom_text(aes(label = ifelse(pct >= 3.0, sprintf("%.1f%%", pct), "")),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.2) +
  scale_fill_manual(values = c(
    "Promoter"   = "#d62728",
    "5' UTR"     = "#ff7f0e",
    "3' UTR"     = "#bcbd22",
    "Exon (CDS)" = "#2ca02c",
    "Intron"     = "#1f77b4",
    "Intergenic" = "#9467bd"
  )) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = NULL, y = "Fraction of peaks (%)",
       fill = NULL,
       title = "Genomic feature composition (scG4 categories)",
       subtitle = "Peaks classified by midpoint overlap priority") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(size = 11, angle = 0),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        legend.position = "right",
        legend.text = element_text(size = 10))

ggsave(file.path(OUT_DIR, "panel_D_feature_distribution_scG4_stacked.pdf"),
       p_stacked, width = 14, height = 6)
cat("  Saved panel_D_feature_distribution_scG4_stacked.pdf\n")

cat("\nDone. Outputs in:", OUT_DIR, "\n")
