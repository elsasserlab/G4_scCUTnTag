# Reusable browser tracks for direct Cicero co-G4 links.

load_saved_r_object <- function(path, object_name) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(object_name, envir = env, inherits = FALSE)) {
    stop("Expected object '", object_name, "' in ", path)
  }
  get(object_name, envir = env, inherits = FALSE)
}

parse_peak_id <- function(peak) {
  fields <- strsplit(peak, "-", fixed = TRUE)[[1]]
  if (length(fields) != 3L) stop("Invalid peak identifier: ", peak)
  list(chr = fields[1], start = as.numeric(fields[2]), end = as.numeric(fields[3]))
}

cicero_peak_universe <- function(...) {
  tables <- list(...)
  peaks <- unique(unlist(lapply(tables, function(x) {
    c(as.character(x$Peak1), as.character(x$Peak2))
  })))
  fields <- data.table::tstrsplit(peaks, "-", fixed = TRUE)
  unique(data.table::data.table(
    chr = fields[[1]], start = as.numeric(fields[[2]]), end = as.numeric(fields[[3]])
  ))
}

direct_cicero_links <- function(conns, focal_peak, group, threshold = 0.2) {
  links <- data.table::as.data.table(conns)[
    (Peak1 == focal_peak | Peak2 == focal_peak) &
      !is.na(coaccess) & coaccess > threshold
  ]
  if (!nrow(links)) return(data.table::data.table())

  links[, partner := ifelse(
    as.character(Peak1) == focal_peak, as.character(Peak2), as.character(Peak1)
  )]
  fields <- data.table::tstrsplit(links$partner, "-", fixed = TRUE)
  links[, `:=`(
    partner_chr = fields[[1]],
    partner_start = as.numeric(fields[[2]]),
    partner_end = as.numeric(fields[[3]]),
    group = group
  )]
  links[order(-coaccess), .SD[1], by = partner]
}

cog4_browser_window <- function(focal_peak, ast_conns, nonast_conns,
                                threshold = 0.2, flank = 100000L,
                                cap = 500000L) {
  focal <- parse_peak_id(focal_peak)
  midpoint <- (focal$start + focal$end) / 2
  links <- data.table::rbindlist(list(
    direct_cicero_links(ast_conns, focal_peak, "AST", threshold),
    direct_cicero_links(nonast_conns, focal_peak, "non-AST", threshold)
  ), fill = TRUE)
  links <- links[partner_chr == focal$chr]
  coordinates <- c(focal$start, focal$end)
  if (nrow(links)) {
    coordinates <- c(coordinates, links$partner_start, links$partner_end)
  }
  c(
    max(1, midpoint - cap, min(coordinates) - flank),
    min(midpoint + cap, max(coordinates) + flank)
  )
}

bigwig_summary_track <- function(path, region, bins = 1400L) {
  values <- summary(
    rtracklayer::BigWigFile(path), region, size = bins,
    type = "mean", defaultValue = 0
  )[[1]]
  data.table::data.table(
    x = (GenomicRanges::start(values) + GenomicRanges::end(values)) / 2,
    score = values$score
  )
}

assign_gene_lanes <- function(genes, span) {
  data.table::setorder(genes, plot_start, plot_end)
  lane_ends <- numeric()
  padding <- span * 0.05
  genes[, lane := 0L]
  for (i in seq_len(nrow(genes))) {
    lane <- which(lane_ends + padding < genes$plot_start[i])[1]
    if (is.na(lane)) {
      lane <- length(lane_ends) + 1L
      lane_ends <- c(lane_ends, -Inf)
    }
    genes$lane[i] <- lane
    lane_ends[lane] <- genes$plot_end[i]
  }
  genes
}

cog4_arc_points <- function(links, focal_mid, x_min, x_max) {
  if (!nrow(links)) return(data.table::data.table())
  span <- x_max - x_min
  data.table::rbindlist(lapply(seq_len(nrow(links)), function(i) {
    partner_mid <- (links$partner_start[i] + links$partner_end[i]) / 2
    left <- min(focal_mid, partner_mid)
    right <- max(focal_mid, partner_mid)
    theta <- seq(pi, 0, length.out = 100)
    direction <- if (links$group[i] == "AST") 1 else -1
    data.table::data.table(
      x = (left + right) / 2 + (right - left) / 2 * cos(theta),
      y = direction * (right - left) / span * sin(theta),
      link_id = paste(links$group[i], links$partner[i], sep = "_"),
      coaccess = links$coaccess[i]
    )
  }))
}

signal_browser_track <- function(values, label, color, ymax, focal, xlim) {
  ggplot2::ggplot(values, ggplot2::aes(x, score)) +
    ggplot2::annotate(
      "rect", xmin = focal$start, xmax = focal$end,
      ymin = -Inf, ymax = Inf, fill = "grey45", alpha = 0.14
    ) +
    ggplot2::geom_area(fill = color, alpha = 0.8) +
    ggplot2::geom_line(color = color, linewidth = 0.25) +
    ggplot2::scale_y_continuous(
      limits = c(0, ymax), expand = ggplot2::expansion(mult = c(0, 0.03))
    ) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::labs(x = NULL, y = label) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 10)
    )
}

make_cog4_browser_plot <- function(
    gene, focal_peak, ast_conns, nonast_conns, peaks, genes, ccre,
    ast_bigwig, nonast_bigwig, pqs_bigwig,
    panel_kind = "AST-up G4", threshold = 0.2,
    ast_color = "#fc5656", nonast_color = "#0bb8b8",
    pqs_color = "#6a3d9a") {
  focal <- parse_peak_id(focal_peak)
  focal_mid <- (focal$start + focal$end) / 2
  window <- cog4_browser_window(
    focal_peak, ast_conns, nonast_conns, threshold = threshold
  )
  x_min <- window[1]
  x_max <- window[2]
  xlim <- c(x_min, x_max)
  region <- GenomicRanges::GRanges(
    focal$chr, IRanges::IRanges(round(x_min), round(x_max))
  )

  ast_links <- direct_cicero_links(ast_conns, focal_peak, "AST", threshold)
  nonast_links <- direct_cicero_links(nonast_conns, focal_peak, "non-AST", threshold)
  links <- data.table::rbindlist(list(ast_links, nonast_links), fill = TRUE)
  links <- links[
    partner_chr == focal$chr & partner_start <= x_max & partner_end >= x_min
  ]

  ast_signal <- bigwig_summary_track(ast_bigwig, region)
  nonast_signal <- bigwig_summary_track(nonast_bigwig, region)
  signal_max <- max(c(ast_signal$score, nonast_signal$score), na.rm = TRUE)
  if (!is.finite(signal_max) || signal_max <= 0) signal_max <- 1

  pqs_gr <- rtracklayer::import(rtracklayer::BigWigFile(pqs_bigwig), which = region)
  pqs <- data.table::data.table(
    x = (GenomicRanges::start(pqs_gr) + GenomicRanges::end(pqs_gr)) / 2,
    score = as.numeric(S4Vectors::mcols(pqs_gr)$score)
  )
  pqs_max <- max(pqs$score, na.rm = TRUE)
  if (!is.finite(pqs_max) || pqs_max <= 0) pqs_max <- 1

  local_ccre <- ccre[chr == focal$chr & start <= x_max & end >= x_min]
  local_peaks <- peaks[chr == focal$chr & start <= x_max & end >= x_min]
  local_genes <- data.table::copy(
    genes[chr == focal$chr & start <= x_max & end >= x_min]
  )
  local_genes[, `:=`(
    plot_start = pmax(start, x_min),
    plot_end = pmin(end, x_max)
  )]
  local_genes[, label_x := (plot_start + plot_end) / 2]
  local_genes <- assign_gene_lanes(local_genes, x_max - x_min)
  n_lanes <- max(c(1L, local_genes$lane))

  p_ast <- signal_browser_track(
    ast_signal, "AST\nRPGC", ast_color, signal_max, focal, xlim
  )
  p_nonast <- signal_browser_track(
    nonast_signal, "non-AST\nRPGC", nonast_color, signal_max, focal, xlim
  )
  p_ccre <- ggplot2::ggplot(local_ccre) +
    ggplot2::annotate(
      "rect", xmin = focal$start, xmax = focal$end,
      ymin = -Inf, ymax = Inf, fill = "grey45", alpha = 0.14
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = k27ac),
      color = NA
    ) +
    ggplot2::scale_fill_gradient(low = "grey80", high = "#e34a33", guide = "none") +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "H3K27ac+\ncCRE") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 10)
    )

  p_genes <- ggplot2::ggplot(local_genes) +
    ggplot2::annotate(
      "rect", xmin = focal$start, xmax = focal$end,
      ymin = -Inf, ymax = Inf, fill = "grey45", alpha = 0.14
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = plot_start, xend = plot_end, y = lane, yend = lane),
      color = "#253494", linewidth = 3.2, lineend = "butt"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = ifelse(strand == "-", plot_start, plot_end), y = lane),
      shape = 18, color = "#253494", size = 1.8
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = label_x, y = lane + 0.28, label = gene), size = 3.2
    ) +
    ggplot2::coord_cartesian(
      xlim = xlim, ylim = c(0, n_lanes + 0.7), clip = "off"
    ) +
    ggplot2::labs(x = NULL, y = "genes") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 10)
    )

  p_peaks <- ggplot2::ggplot(local_peaks) +
    ggplot2::annotate(
      "rect", xmin = focal$start, xmax = focal$end,
      ymin = -Inf, ymax = Inf, fill = "grey45", alpha = 0.14
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = start, xmax = end, ymin = 0.2, ymax = 0.8),
      fill = "grey30", color = NA
    ) +
    ggplot2::annotate(
      "rect", xmin = focal$start, xmax = focal$end,
      ymin = 0.08, ymax = 0.92, fill = "black"
    ) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "G4 peaks") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 10)
    )

  arcs <- cog4_arc_points(links, focal_mid, x_min, x_max)
  arc_max <- max(c(threshold + 0.01, arcs$coaccess), na.rm = TRUE)
  p_arcs <- ggplot2::ggplot(
    arcs, ggplot2::aes(x, y, group = link_id, color = coaccess, linewidth = coaccess)
  ) +
    ggplot2::annotate(
      "rect", xmin = focal$start, xmax = focal$end,
      ymin = -Inf, ymax = Inf, fill = "grey45", alpha = 0.14
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3) +
    ggplot2::geom_path(alpha = 0.95) +
    ggplot2::annotate(
      "text", x = x_min, y = 0.92, label = "AST", hjust = 0,
      color = ast_color, fontface = "bold", size = 4
    ) +
    ggplot2::annotate(
      "text", x = x_min, y = -0.92, label = "non-AST", hjust = 0,
      color = nonast_color, fontface = "bold", size = 4
    ) +
    ggplot2::scale_color_gradient(
      low = "grey75", high = "#d7191c", limits = c(threshold, arc_max),
      name = "coaccess"
    ) +
    ggplot2::scale_linewidth_continuous(
      range = c(0.25, 0.8), limits = c(threshold, arc_max), guide = "none"
    ) +
    ggplot2::scale_x_continuous(labels = function(x) sprintf("%.2f Mb", x / 1e6)) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(-1, 1)) +
    ggplot2::labs(x = focal$chr, y = "direct co-G4\nlinks") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 10), legend.position = "right"
    )

  p_pqs <- ggplot2::ggplot(pqs, ggplot2::aes(x = x, ymin = 0, ymax = score)) +
    ggplot2::annotate(
      "rect", xmin = focal$start, xmax = focal$end,
      ymin = -Inf, ymax = Inf, fill = "grey45", alpha = 0.14
    ) +
    ggplot2::geom_linerange(color = pqs_color, linewidth = 0.22, alpha = 0.8) +
    ggplot2::scale_y_continuous(
      limits = c(0, pqs_max), expand = ggplot2::expansion(mult = c(0, 0.03))
    ) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::labs(x = NULL, y = "PQS\nscore") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 10)
    )

  title <- sprintf(
    "%s | %s | direct links shown: AST=%d, non-AST=%d",
    gene, panel_kind, sum(links$group == "AST"), sum(links$group == "non-AST")
  )
  patchwork::wrap_plots(
    list(p_ast, p_nonast, p_ccre, p_genes, p_peaks, p_arcs, p_pqs),
    ncol = 1,
    heights = c(1, 1, 0.38, max(0.9, n_lanes * 0.48), 0.38, 2.2, 0.9)
  ) +
    patchwork::plot_annotation(
      title = title,
      subtitle = paste(
        "Unsmoothed RPGC BigWigs; discrete PQS intervals;",
        "all annotated genes in view"
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 10, color = "grey30")
      )
    )
}

make_cog4_browser_panel <- function(loci, panel_kinds, ...) {
  if (is.null(names(loci)) || any(names(loci) == "")) {
    stop("loci must be a named character vector of peak identifiers")
  }
  if (length(panel_kinds) == 1L) panel_kinds <- rep(panel_kinds, length(loci))
  if (length(panel_kinds) != length(loci)) {
    stop("panel_kinds must have length 1 or length(loci)")
  }
  plots <- lapply(seq_along(loci), function(i) {
    make_cog4_browser_plot(
      gene = names(loci)[i], focal_peak = loci[[i]],
      panel_kind = panel_kinds[i], ...
    )
  })
  patchwork::wrap_plots(
    lapply(plots, function(x) patchwork::wrap_elements(full = x)), ncol = 1
  )
}
