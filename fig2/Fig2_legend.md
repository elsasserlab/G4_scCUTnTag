# Figure 2 — Legend and Derivation Notes

**Figure 2: Integration of G4 scCUT&Tag data with scRNA-seq using canonical correlation analysis (CCA).**

## Data sources

- **Sorted GFP+ brain G4 scCUT&Tag (GSE291468, GSM8836086)** — `GSM8836086_GFPpos_Seurat_object.Rds`, `GSM8836086_GFP_sorted_mouse_brain/CellRanger/fragments.tsv.gz`, peak assay, gene-activity (GA) assay
- **Brain scRNA-seq reference (GSE163484, Bartosovic et al.)** — processed Seurat object with cell_type annotations
- **Seurat CCA integration outputs** — `scRNA_Seq_Seurat_object.Rds`, `G4_scRNA_integration.Rds` and `g4_cell_label_preds.Rds` under `fig2/results/integration/outputs/`

Produced by `fig2/fig2.R`, with QC by `fig2/fig2_S1_QC.R` and preprocessing in `fig2/seurat_integration.R`.

## Panels

**Panel S1 — QC violin plots for GFP+ sorted brain.** The 2,402 cells comprised clusters 0-3 with 1,102, 513, 402 and 385 cells (an empty factor level 4 is unused). Median values by cluster 0/1/2/3 were 292.5/416/495.5/452 detected peaks, 601.5/862/1,029/955 peak counts, 246/387/415.5/381 TSS fragments, 43/38/29/69 mitochondrial fragments and FRiP 0.128/0.149/0.135/0.135. Crossbars denote medians and points means; no inferential test was performed. Output: `panel_S1_qc_GFPpos.pdf`.

**Panel B — Co-embedded UMAPs.** Four views share a Seurat CCA-derived embedding of 5,488 reference RNA cells and 2,402 G4 cells (7,890 total). G4 label counts were AST 1,476, MOL 389, OEC 177, pericytes 146, COP-NFOL 73, VLMC 60, OPC 41 and VEC 40. The fourth view shows Seurat `pred_max_score` (range 0.185-0.995), not scBridge reliability. Output: `panel_B_umap.pdf`.

**Panel C — Prediction-score heatmap.** Heatmap of Seurat label-transfer scores for eight reference cell types across all 2,402 GFP+ G4 cells (8 x 2,402 matrix after removing the `max` row). Output: `panel_C_prediction_scores.pdf`.

**Panel D — Cluster/prediction correspondence.** **D1:** UMAP of the 2,402 GFP+ G4 cells, colored by Seurat cluster (n = 1,102, 513, 402 and 385). **D2:** overlap score `min(cluster fraction, cell-type fraction)` from the cluster x prediction contingency table. The highest score for clusters 0-3 was AST (0.353), AST (0.282), AST (0.248) and MOL (0.252), respectively. No inferential test was performed. Outputs: `panel_D1_umap_clusters.pdf`, `panel_D2_overlap_heatmap.pdf`.

**Panel E — scRNA-seq marker analysis.** Reference composition was AST 2,632, MOL 1,191, OEC 636, pericytes 424, COP-NFOL 238, VLMC 153, VEC 126 and OPC 88. The top displayed markers were *Igsf1*, *Gpr17*, *Foxd3*, *S1pr5*, *Dcx*, *Aspn*, *AU021092* and *Fam180a*. Marker filtering used Seurat Bonferroni-adjusted `p_val_adj < 0.05` and `avg_log2FC > 0`; violin plots are descriptive and have no additional test. Outputs: `panel_E_umap.pdf`, `panel_E_markers.pdf`.

**Panel F — G4/RNA co-enrichment feature plots.** Nine genes are plotted: the eight Panel E markers plus *Mbp*. Normalized RNA expression and G4 gene activity are displayed on the shared 7,890-cell embedding; RNA uses a yellow-teal-blue gradient and G4 a pale-to-dark-red gradient. These plots are descriptive and no test is performed. Outputs: `panel_F_RNA_featureplots.pdf`, `panel_F_G4_featureplots.pdf`, `panel_F_combined_RNA_G4.pdf`.
