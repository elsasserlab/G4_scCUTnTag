#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = max(1L, parallel::detectCores() - 1L))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran <- c(
  "Seurat", "Signac", "tidyverse", "data.table", "glue", "argparse",
  "cowplot", "ggpubr", "ggrepel", "ggrastr", "patchwork", "viridis",
  "scales", "RColorBrewer", "readxl", "gridExtra", "ggplotify", "eulerr",
  "VennDiagram", "hrbrthemes", "Matrix", "pacman", "enrichR", "circlize"
)
bioc <- c(
  "GenomicRanges", "IRanges", "GenomeInfoDb", "GenomicFeatures",
  "rtracklayer", "BSgenome", "BSgenome.Mmusculus.UCSC.mm10",
  "EnsDb.Mmusculus.v79", "org.Mm.eg.db",
  "TxDb.Mmusculus.UCSC.mm10.knownGene", "ChIPseeker", "pqsfinder",
  "ComplexHeatmap", "Gviz", "regioneR", "cicero", "monocle3"
)
github <- c(
  "cnluzon/wigglescout",
  "cnluzon/bedscout",
  "egeulgen/DoubletFinder"
)

installed <- rownames(installed.packages())
to_cran <- setdiff(cran, installed)
if (length(to_cran)) install.packages(to_cran)

installed <- rownames(installed.packages())
to_bioc <- setdiff(bioc, installed)
if (length(to_bioc)) BiocManager::install(to_bioc, update = FALSE, ask = FALSE)

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github(github, upgrade = "never")
