#!/usr/bin/env Rscript
# Extract cluster barcodes from GSM8836086 Seurat object
# This script creates the Seurat_cluster_barcodes/ directory with barcode files

library(Seurat)

# Read Seurat object
seurat_file <- "GSM8836086_GFPpos_Seurat_object.Rds"
cat("Reading Seurat object:", seurat_file, "\n")
seurat <- readRDS(seurat_file)

# Get cluster assignments
clusters <- Idents(seurat)
cat("Found", length(clusters), "cells with", length(unique(clusters)), "clusters\n")

# Create output directory
outdir <- "Seurat_cluster_barcodes"
dir.create(outdir, showWarnings = FALSE)

# Extract barcodes per cluster
for (cluster in sort(unique(clusters))) {
  outfile <- file.path(outdir, paste0("barcodes_cluster_", cluster, ".tsv"))
  barcodes <- names(clusters[clusters == cluster])
  
  cat("Cluster", cluster, ":", length(barcodes), "cells ->", outfile, "\n")
  
  # Write as TSV (matching original format - barcodes only in first column)
  write.table(
    data.frame(barcode = barcodes, cluster = cluster),
    file = outfile,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  # Also create barcode-only file for compatibility
  write.table(
    barcodes,
    file = sub(".tsv", "_barcodes.txt", outfile),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

cat("\nDone! Created files in", outdir, "/\n")
cat("Total cells:", length(clusters), "\n")
cat("Cluster sizes:\n")
print(table(clusters))