#!/usr/bin/env Rscript
# Regenerate the Cicero connection and CCAN intermediates used by Figure 3E.
# This compute-intensive step is optional because its outputs are shipped in
# fig3/results/cicero/.

suppressPackageStartupMessages({
  library(cicero)
  library(monocle3)
  library(Seurat)
  library(data.table)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
ROOT <- dirname(SCRIPT_DIR)
OUT <- file.path(SCRIPT_DIR, "results", "cicero")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

run_cicero_analysis <- function(seurat_object, connection_file, ccan_file) {
  set.seed(42)
  counts <- GetAssayData(seurat_object, assay = "peaks", layer = "counts")
  counts <- counts[Matrix::rowSums(counts) > 0, , drop = FALSE]
  cell_metadata <- seurat_object@meta.data[colnames(counts), , drop = FALSE]
  peak_metadata <- data.frame(
    gene_short_name = rownames(counts), row.names = rownames(counts)
  )

  input_cds <- monocle3::new_cell_data_set(
    expression_data = counts,
    cell_metadata = cell_metadata,
    gene_metadata = peak_metadata
  )
  input_cds <- monocle3::estimate_size_factors(input_cds)
  input_cds <- monocle3::preprocess_cds(input_cds, method = "LSI")
  input_cds <- monocle3::reduce_dimension(
    input_cds, reduction_method = "UMAP", preprocess_method = "LSI"
  )
  cicero_cds <- cicero::make_cicero_cds(
    input_cds,
    reduced_coordinates = SingleCellExperiment::reducedDims(input_cds)$UMAP
  )

  chrom_sizes <- fread(file.path(ROOT, "data", "genome", "mm10.chrom.sizes.txt"))
  conns <- cicero::run_cicero(cicero_cds, chrom_sizes)
  save(conns, file = connection_file)

  CCAN_assigns <- cicero::generate_ccans(
    conns, coaccess_cutoff_override = 0.1
  )
  save(CCAN_assigns, file = ccan_file)
}

g4_path <- file.path(
  ROOT, "data", "GSE291468", "GSM8836086_GFPpos_Seurat_object.Rds"
)
prediction_path <- file.path(
  SCRIPT_DIR, "results", "scBridge", "scbridge_predictions.csv"
)
if (!file.exists(g4_path) || !file.exists(prediction_path)) {
  stop("Missing Cicero input. Run download_data.sh --fig3 and verify scBridge outputs.")
}

cat("Loading GFP+ G4 object and scBridge predictions...\n")
g4 <- readRDS(g4_path)
predictions <- fread(prediction_path)
setnames(predictions, names(predictions)[1], "barcode")
ast_cells <- intersect(
  predictions[Prediction == "Astrocytes", barcode], colnames(g4)
)
nonast_cells <- setdiff(colnames(g4), ast_cells)
if (!length(ast_cells) || !length(nonast_cells)) {
  stop("Could not define both AST and non-AST cell groups from scBridge predictions")
}

analyses <- list(
  combined = list(
    object = g4,
    conns = "cicero_GFPsorted.Rds",
    ccans = "cicero_GFPsorted_coG4networks.Rds"
  ),
  predAST = list(
    object = subset(g4, cells = ast_cells),
    conns = "cicero_GFPsorted-predAST.Rds",
    ccans = "cicero_GFPsorted_coG4networks-predAST.Rds"
  ),
  pred_nonAST = list(
    object = subset(g4, cells = nonast_cells),
    conns = "cicero_GFPsorted-pred_nonAST.Rds",
    ccans = "cicero_GFPsorted_coG4networks-pred_nonAST.Rds"
  )
)

for (label in names(analyses)) {
  analysis <- analyses[[label]]
  cat("Running Cicero:", label, "\n")
  run_cicero_analysis(
    analysis$object,
    file.path(OUT, analysis$conns),
    file.path(OUT, analysis$ccans)
  )
}

cat("Cicero intermediates written to", OUT, "\n")
