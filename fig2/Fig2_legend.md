# Figure 2 — Legend and Derivation Notes

**Figure 2: Integration of G4 scCUT&Tag data with scRNA-seq using canonical correlation analysis (CCA).**

## Data sources

- **Sorted GFP+ brain G4 scCUT&Tag (GSE291468, GSM8836086)** — `GSM8836086_GFPpos_Seurat_object.Rds`, `GSM8836086_GFP_sorted_mouse_brain/CellRanger/fragments.tsv.gz`, peak assay, gene-activity (GA) assay
- **Brain scRNA-seq reference (GSE163484, Bartosovic et al.)** — processed Seurat object with cell_type annotations
- **scBridge co-embedding** — precomputed outputs from `results/scBridge/output/GFPsorted_Bartosovic/`: `umap_coembedded.csv`, `scbridge_predictions.csv`, `scbridge_reliability.csv`

Produced by `fig2/fig2.R`, with QC by `fig2/fig2_S1_QC.R` and preprocessing in `fig2/seurat_integration.R`.

## Panels

**Panel S1 — QC violin plots for GFP+ sorted brain.** Seven per-cluster metrics computed with Signac: nFeature_peaks, nCount_peaks, TSS_fragments (`TSSEnrichment(fast=TRUE)`), mitochondrial fragments, FRiP (`FRiP()` with total_fragments from `CountFragments()`), TSS_enrichment, nucleosome_signal (`NucleosomeSignal()`). Plotted as violin plots per Seurat cluster with median (crossbar) and mean (point) overlaid. Tool: Seurat `VlnPlot()`. Output: `panel_S1_qc_GFPpos.pdf`.

**Panel B — Co-embedded UMAPs.** Four-view grid showing joint embedding of GFP+ G4 cells (GSM8836086) and scRNA-seq reference (GSE163484). UMAP computed on merged Seurat object: VariableFeatures from RNA assay, ScaleData (do.scale=FALSE), RunPCA, RunUMAP (dims=1:15, reduction.name="fig2_umap"). Four views share same coordinates: (top-left) reference cell types (scRNA-seq cells only), (top-right) data modality (scRNA-seq vs G4 scCUT&Tag), (bottom-left) transferred predictions for G4 cells (from scbridge_predictions.csv, "Astrocytes"→"AST", "Oligodendrocytes"→"MOL", NA→"unreliable"), (bottom-right) integration reliability score (from scbridge_reliability.csv, viridis scale). Tools: Seurat `RunUMAP()`, ggplot2. Output: `panel_B_umap.pdf`.

**Panel C — Prediction-score heatmap.** Heatmap of Seurat label transfer prediction scores for GFP+ G4 cells across reference cell types. Matrix from `g4_cell_label_preds.Rds` (@data slot), "max" row removed. Rendered with ComplexHeatmap `Heatmap()` (col=viridis(100), show_column_names=FALSE). Tool: ComplexHeatmap. Output: `panel_C_prediction_scores.pdf`.

**Panel D — Cluster/prediction correspondence.** Two panels: **D1** (left): UMAP of GFP+ G4 cells colored by Seurat clusters (from `seurat_clusters` metadata). Tool: Seurat `DimPlot(label=TRUE, repel=TRUE, pt.size=1.2)`. **D2** (right): symmetric overlap-score heatmap quantifying agreement between Seurat clusters and predicted cell types. Overlap computed as min(cluster_fraction, type_fraction) from contingency table (clusters × predictions, with "Astrocytes"→"AST", "Oligodendrocytes"→"MOL"). Heatmap rendered with ComplexHeatmap `Heatmap()` (col=viridis(100, option="D"), row_title="Predicted cell type", column_title="Seurat cluster", cluster_rows=TRUE, cluster_columns=TRUE, cell values shown). Tools: Seurat, ComplexHeatmap. Outputs: `panel_D1_umap_clusters.pdf`, `panel_D2_overlap_heatmap.pdf`.

**Panel E — scRNA-seq marker analysis.** Two panels: (left) Reference UMAP colored by cell_type (scRNA-seq cells only from joint embedding). (right) Violin plots of top differentially expressed marker gene per cell type from `scRNA-Seq-FindAllMarkers_output.tsv` (Seurat `FindAllMarkers()` on GSE163484: filter p_val_adj<0.05, avg_log2FC>0, group_by cluster, slice_max avg_log2FC, n=1). Tools: Seurat, ggplot2. Outputs: `panel_E_umap.pdf`, `panel_E_markers.pdf`.

**Panel F — G4/RNA co-enrichment feature plots.** For each marker gene from Panel E, two UMAP feature plots on shared embedding coordinates: (left) normalized scRNA-seq expression (RNA assay, layer="data" from GSE163484), (right) G4 gene-activity score (GA assay, layer="data" from GSM8836086). Plots show co-enrichment of G4 activity and expression within corresponding cell types. RNA plotted with viridis scale (direction=-1), G4 plotted with white-to-red gradient. Tools: Seurat `GetAssayData()`, ggplot2. Outputs: `panel_F_RNA_featureplots.pdf`, `panel_F_G4_featureplots.pdf` (one PDF per modality with all marker genes arranged in 3-column grid).