#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Genomic feature distribution comparison: bulk consensus vs scCUT&Tag
# pseudobulk, in mESC and MEF.
#
# Pipeline:
#   1. Load peaks (rtracklayer::import, preserves signalValue).
#   2. Build bulk consensus via bedscout::loci_consensus(..., min_consensus = 2)
#      and re-attach max signal.
#   3. Build feature interval sets from GENCODE mm10 GTF:
#        Promoter (TSS ± 2 kb), 5' UTR, 3' UTR, Exon (CDS), Intron.
#   4. Classify each peak (by midpoint) using a priority cascade with
#      overlapsAny() from GenomicRanges — this guarantees mutually exclusive
#      categories. (This replaces the comma-separated-string parsing approach
#      via bedscout::annotate_overlapping_features which was unreliable.)
#   5. Plot grouped + stacked bars: bulk consensus vs pseudobulk per cell type.
#
# A separate per-peak annotation file is also produced via
# bedscout::annotate_overlapping_features() for inspection.
#
# Install bedscout: devtools::install_github("cnluzon/bedscout")
# ---------------------------------------------------------------------------


source("rev/paths.R")
suppressPackageStartupMessages({
  library(bedscout)
  library(GenomicRanges)
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
OUT_DIR  <- file.path(OUT_ROOT, "feature_distribution_bedscout")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CANONICAL   <- c(paste0("chr", 1:19), "chrX", "chrY")
PROMOTER_UP <- 2000
PROMOTER_DN <- 2000

peak_files <- c(
  mESC_bulk_rep1   = file.path(BULK_PEAKS, "GSM8836082_bulkG4CnT_mESC_rep1.broadPeak"),
  mESC_bulk_rep2   = file.path(BULK_PEAKS, "GSM8836083_bulkG4CnT_mESC_rep2.broadPeak"),
  MEF_bulk_rep1    = file.path(BULK_PEAKS, "GSM8836084_bulkG4CnT_3T3_rep1.broadPeak"),
  MEF_bulk_rep2    = file.path(BULK_PEAKS, "GSM8836085_bulkG4CnT_3T3_rep2.broadPeak"),
  mESC_pseudobulk  = CLUSTER_PEAKS_1,
  MEF_pseudobulk   = CLUSTER_PEAKS_0
)

# ===== 1. LOAD PEAKS =====
cat("Loading peaks...\n")
peaks <- lapply(peak_files, function(f) {
  fmt <- if (grepl("broadPeak$", f)) "broadPeak" else "narrowPeak"
  gr  <- rtracklayer::import(f, format = fmt)
  gr  <- gr[as.character(seqnames(gr)) %in% CANONICAL]
  gr$signal <- gr$signalValue
  sort(gr)
})
for (nm in names(peaks)) cat(sprintf("  %-18s %d peaks\n", nm, length(peaks[[nm]])))

# ===== 2. BULK CONSENSUS =====
attach_signal <- function(consensus_gr, rep1, rep2) {
  combined <- c(rep1, rep2)
  ov <- findOverlaps(consensus_gr, combined, ignore.strand = TRUE)
  agg <- tapply(combined$signal[subjectHits(ov)],
                queryHits(ov), max, default = NA_real_)
  out_signal <- rep(NA_real_, length(consensus_gr))
  out_signal[as.integer(names(agg))] <- as.numeric(agg)
  consensus_gr$signal <- out_signal
  consensus_gr
}

cat("\nBuilding bulk consensus via bedscout::loci_consensus()...\n")
consensus <- list()
for (ct in c("mESC", "MEF")) {
  r1 <- peaks[[paste0(ct, "_bulk_rep1")]]
  r2 <- peaks[[paste0(ct, "_bulk_rep2")]]
  cons <- loci_consensus(list(r1, r2), min_consensus = 2)
  cons <- attach_signal(cons, r1, r2)
  consensus[[ct]] <- cons
  cat(sprintf("  %s consensus: %d peaks\n", ct, length(cons)))
}

# ===== 3. BUILD FEATURE INTERVAL SETS FROM GTF =====
cat("\nBuilding feature intervals from GTF...\n")
gtf <- rtracklayer::import(GTF_PATH)
gtf <- gtf[as.character(seqnames(gtf)) %in% CANONICAL]
seqlevels(gtf) <- intersect(seqlevels(gtf), CANONICAL)

# 3a. Promoters: TSS ± 2 kb on transcripts
tx <- gtf[gtf$type == "transcript"]
promoters_gr <- reduce(
  promoters(tx, upstream = PROMOTER_UP, downstream = PROMOTER_DN),
  ignore.strand = TRUE
)

# 3b. UTRs — GENCODE combines them; split by strand + start_codon position
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

# 3c. CDS exons = all exons minus UTRs
exons_all      <- reduce(gtf[gtf$type == "exon"], ignore.strand = TRUE)
exon_nonutr_gr <- GenomicRanges::setdiff(
  exons_all, reduce(c(utr5_gr, utr3_gr), ignore.strand = TRUE),
  ignore.strand = TRUE
)

# 3d. Introns = transcript span minus exons
tx_span   <- reduce(tx, ignore.strand = TRUE)
intron_gr <- GenomicRanges::setdiff(tx_span, exons_all, ignore.strand = TRUE)

# 3e. Stash in a named list (one GRanges per priority class)
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

# ===== 4. CLASSIFY PEAKS (priority cascade, single category per peak) =====
# Robust approach: take peak midpoints, then check overlap against each feature
# class in priority order. First match wins; otherwise label as "Intergenic".
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
peak_sets <- list(
  `mESC bulk consensus` = consensus$mESC,
  `mESC pseudobulk`     = peaks$mESC_pseudobulk,
  `MEF bulk consensus`  = consensus$MEF,
  `MEF pseudobulk`      = peaks$MEF_pseudobulk
)

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

# Long table for plotting / CSV.
# Use explicit source assignment via match() — avoids the previous grepl("bulk")
# bug where "pseudobulk" also matched "bulk".
peak_set_source <- c(
  `mESC bulk consensus` = "Bulk consensus",
  `mESC pseudobulk`     = "Pseudobulk (sc)",
  `MEF bulk consensus`  = "Bulk consensus",
  `MEF pseudobulk`      = "Pseudobulk (sc)"
)
peak_set_celltype <- c(
  `mESC bulk consensus` = "mESC",
  `mESC pseudobulk`     = "mESC",
  `MEF bulk consensus`  = "MEF",
  `MEF pseudobulk`      = "MEF"
)

dist_dt <- rbindlist(lapply(names(distros), function(label) {
  total <- sum(distros[[label]])
  data.table(
    peak_set  = label,
    feature   = CATS,
    n         = distros[[label]],
    pct       = 100 * distros[[label]] / total,
    source    = peak_set_source[[label]],
    cell_type = peak_set_celltype[[label]]
  )
}))
fwrite(dist_dt, file.path(OUT_DIR, "feature_distribution.csv"))

cat("\nSanity check: source label counts in dist_dt:\n")
print(table(dist_dt$source))

# ===== 5. (Optional) bedscout::annotate_overlapping_features() for inspection =====
# Build one combined features GRanges with a feature_class column, call once
# per peak set, and save the per-peak annotation (does NOT drive the bar
# plot — that's the priority cascade above).
build_named <- function(gr, label) {
  strand(gr) <- "*"
  gr$feature_class <- label
  gr
}
features_gr <- c(
  build_named(promoters_gr,   "Promoter"),
  build_named(utr5_gr,        "5' UTR"),
  build_named(utr3_gr,        "3' UTR"),
  build_named(exon_nonutr_gr, "Exon (CDS)"),
  build_named(intron_gr,      "Intron")
)
mcols(features_gr) <- mcols(features_gr)[, "feature_class", drop = FALSE]

cat("\nPer-peak annotation via bedscout::annotate_overlapping_features() (for inspection)...\n")
for (label in names(peak_sets)) {
  ann <- tryCatch(
    annotate_overlapping_features(
      peak_sets[[label]], features_gr,
      name_field = "feature_class", ignore.strand = TRUE
    ),
    error = function(e) { cat("  failed:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(ann)) {
    fn <- gsub(" ", "_", paste0(label, "_features.tsv"))
    fwrite(as.data.table(ann), file.path(OUT_DIR, fn), sep = "\t")
  }
}

# ===== 6. PLOTS =====
cat("\nDrawing figures...\n")
dist_dt$feature <- factor(dist_dt$feature, levels = CATS)
dist_dt$source  <- factor(dist_dt$source,
                          levels = c("Bulk consensus", "Pseudobulk (sc)"))
dist_dt$cell_type <- factor(dist_dt$cell_type, levels = c("mESC", "MEF"))

# 6a. Grouped bars (bulk vs sc, side by side per category)
p_grouped <- ggplot(dist_dt, aes(x = feature, y = pct, fill = source)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", pct)),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 2.8) +
  facet_wrap(~ cell_type, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("Bulk consensus"  = "#4682b4",
                                "Pseudobulk (sc)" = "#ff9933")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Fraction of peaks (%)",
       title = "Genomic feature distribution of G4 peaks",
       subtitle = "Bulk consensus vs scCUT&Tag pseudobulk",
       fill = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        strip.text  = element_text(face = "bold", size = 12),
        plot.title  = element_text(face = "bold"),
        legend.position = "top")
ggsave(file.path(OUT_DIR, "feature_distribution_grouped.pdf"),
       p_grouped, width = 11, height = 5.5)
cat("  saved feature_distribution_grouped.pdf\n")

# 6b. Stacked bars (composition view)
p_stacked <- ggplot(dist_dt, aes(x = source, y = pct, fill = feature)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.4, width = 0.65) +
  geom_text(aes(label = ifelse(pct >= 2.5, sprintf("%.1f%%", pct), "")),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3) +
  facet_wrap(~ cell_type, ncol = 2) +
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
       title = "Genomic feature composition",
       subtitle = "Bulk consensus vs scCUT&Tag pseudobulk") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold"),
        legend.position = "right")
ggsave(file.path(OUT_DIR, "feature_distribution_stacked.pdf"),
       p_stacked, width = 11, height = 5.5)
cat("  saved feature_distribution_stacked.pdf\n")

cat("\nDone. Outputs in:", OUT_DIR, "\n")
