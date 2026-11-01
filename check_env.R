# =============================================================
# check_env.R — verify that the R environment has all packages
# used by the figure and revision scripts in this repository.
#
# Run inside the `scg4` conda environment AFTER install_dependencies.R:
#   conda activate scg4
#   Rscript check_env.R
#
# Exits with status 1 and lists any missing packages, so that a
# wrong/incomplete R installation is caught before running the
# figure scripts (common cause of scripts "not running").
# =============================================================

required <- c(
  # Bioconductor
  "GenomicRanges", "IRanges", "rtracklayer", "regioneR",
  "ChIPseeker", "ComplexHeatmap", "Gviz", "circlize",
  "EnsDb.Mmusculus.v79", "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "org.Mm.eg.db", "monocle3",
  # CRAN
  "Seurat", "Signac", "bedscout", "wigglescout", "cicero",
  "data.table", "readxl", "dplyr", "tidyr", "ggplot2", "ggpubr",
  "ggrastr", "ggrepel", "patchwork", "cowplot", "eulerr",
  "VennDiagram", "RColorBrewer", "scales", "stringr", "glue",
  "viridis", "enrichR", "ggplotify", "gridExtra"
)

missing <- required[!vapply(required, function(p)
  requireNamespace(p, quietly = TRUE), logical(1))]

cat("R version :", R.version.string, "\n")
cat(".libPaths :", .libPaths(), sep = "\n  ")
cat("\nPackages checked :", length(required), "\n")

if (length(missing) > 0) {
  cat("\nMISSING packages:\n  ", paste(missing, collapse = ", "), "\n\n",
      "Fix: activate the correct environment and (re)run the installer:\n",
      "  conda activate scg4\n",
      "  Rscript install_dependencies.R\n", sep = "")
  quit(status = 1)
}

cat("\nAll required packages are installed.\n")
