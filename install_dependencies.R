#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = max(1L, parallel::detectCores() - 1L))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran <- c(
  "Signac", "tidyverse", "data.table", "glue", "argparse",
  "cowplot", "ggpubr", "ggrepel", "ggrastr", "patchwork", "viridis",
  "scales", "RColorBrewer", "readxl", "gridExtra", "ggplotify", "eulerr",
  "VennDiagram", "hrbrthemes", "Matrix", "pacman", "circlize",
  "matrixStats", "fields"
)
bioc <- c(
  "GenomicRanges", "IRanges", "GenomeInfoDb", "GenomicFeatures",
  "rtracklayer", "BSgenome", "BSgenome.Mmusculus.UCSC.mm10",
  "EnsDb.Mmusculus.v79", "org.Mm.eg.db", "ensembldb", "glmGamPoi",
  "plyranges",
  "TxDb.Mmusculus.UCSC.mm10.knownGene", "ChIPseeker", "pqsfinder",
  "ComplexHeatmap", "Gviz", "regioneR", "cicero", "monocle3"
)
github <- c(
  "cnluzon/wigglescout",
  "cnluzon/bedscout",
  "egeulgen/DoubletFinder",
  "mojaveazure/loomR"
)

installed <- rownames(installed.packages())
to_cran <- setdiff(cran, installed)
if (length(to_cran)) install.packages(to_cran)

installed <- rownames(installed.packages())
to_bioc <- setdiff(bioc, installed)
if (length(to_bioc)) BiocManager::install(to_bioc, update = FALSE, ask = FALSE)

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# Seurat is version-pinned: the shipped Seurat objects and integration
# intermediates were produced with Seurat 5.1.0 / SeuratObject 5.0.2
# (see scripts/*/session_info.txt). Newer Seurat releases change the
# integration anchors and are NOT compatible with this repository.
if (!"SeuratObject" %in% installed ||
    compareVersion(as.character(packageVersion("SeuratObject")), "5.0.2") != 0) {
  remotes::install_version("SeuratObject", version = "5.0.2", upgrade = "never")
}
if (!"Seurat" %in% installed ||
    compareVersion(as.character(packageVersion("Seurat")), "5.1.0") != 0) {
  remotes::install_version("Seurat", version = "5.1.0", upgrade = "never")
}

remotes::install_github(github, upgrade = "never")

# enrichR was removed from CRAN; install the last available version from the
# CRAN archive (fig4/fig4_enrichment.R and rev/rev_F25_enrichR_heatmap.R)
if (!requireNamespace("enrichR", quietly = TRUE)) {
  tryCatch(
    remotes::install_version("enrichR", version = "3.2", upgrade = "never"),
    error = function(e) {
      remotes::install_url(
        "https://cran.r-project.org/src/contrib/Archive/enrichR/enrichR_3.2.tar.gz"
      )
    }
  )
}
