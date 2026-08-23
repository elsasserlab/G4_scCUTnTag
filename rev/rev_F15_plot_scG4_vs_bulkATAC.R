#!/usr/bin/env Rscript
# Re-render the scG4-vs-bulkATAC fold plots from the results CSV.
# Tweaks: legend 'G4' -> 'scG4'; thinner bars, closer together within group.

source("rev/paths.R")
suppressMessages(library(ggplot2))
OUTDIR   <- file.path(OUT_ROOT, "scG4cluster_vs_bulkATAC_PQS_eG4")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# prefer the analysis script's CSV; fall back to the bundled tables copy (F15)
csv_out <- file.path(OUTDIR, "scG4cluster_vs_bulkATAC_PQS_eG4_results.csv")
csv_tab <- file.path(ROOT, "tables", "F15_scG4cluster_vs_bulkATAC_PQS_eG4_results.csv")
res <- read.csv(if (file.exists(csv_out)) csv_out else csv_tab)
res$group <- factor(res$group, levels = c("3T3", "mESC"))
res$type  <- factor(res$type,  levels = c("G4", "ATAC"))

DODGE <- 0.45                      # narrow spacing between the two bars in a group
dodge <- position_dodge(width = DODGE)
HALF  <- DODGE / 4                 # bar-centre offset from group centre (2 bars)

make_plot <- function(rn) {
  d <- res[res$reference == rn, ]
  br <- do.call(rbind, lapply(levels(d$group), function(g) {
    gg <- d[d$group == g, ]; xi <- as.integer(factor(g, levels = levels(d$group)))
    data.frame(x = xi - HALF, xend = xi + HALF, y = max(gg$fold_hi) * 1.06,
               lab = gg$sig[gg$type == "G4"])
  }))
  ggplot(d, aes(group, fold, fill = type)) +
    geom_col(position = dodge, width = 0.40, color = "black") +
    geom_errorbar(aes(ymin = fold_lo, ymax = fold_hi), position = dodge, width = 0.12, linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2fx", fold)), position = dodge, vjust = -0.6, size = 3.0) +
    geom_segment(data = br, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = br, aes(x = (x + xend) / 2, y = y, label = lab), inherit.aes = FALSE,
              vjust = -0.3, fontface = "bold", size = 4) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_fill_manual(values = c(G4 = "#B2182B", ATAC = "#2166AC"),
                      labels = c(G4 = "scG4", ATAC = "ATAC")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(title = sprintf("Fold enrichment over random - %s", rn),
         subtitle = paste0("scG4 pseudobulk cluster vs matched bulk ATAC\n",
                           "genome-wide | 1000 perms, 95% bootstrap CI | bracket: scG4 > ATAC (bootstrap, ***p<0.001)"),
         x = NULL, y = "Fold enrichment over random", fill = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9.5, lineheight = 1.1),
          axis.text.x = element_text(size = 12, face = "bold"))
}
fig_map <- c(eG4="F15a_scG4_vs_matched_bulkATAC_eG4_foldenrichment", PQS="F15b_scG4_vs_matched_bulkATAC_PQS_foldenrichment")
for (rn in c("PQS", "eG4")) {
  p <- make_plot(rn)
  ggsave(file.path(OUTDIR, sprintf("scG4_vs_bulkATAC_%s_barplot.pdf", rn)), p, width = 6.5, height = 6)
  cat("published:", fig_map[[rn]], "\n")
}
