# Figure 4 — Legend and Derivation Notes

**Figure 4: Integration of the unsorted (GFP−/GFP+) brain G4 scCUT&Tag dataset with neuronal scRNA-seq aids identification of neuronal cell types.**

## Data sources

- **Unsorted mouse-brain G4 scCUT&Tag (GSE291468, GSM8836087)** — `GSM8836087_unsorted_Seurat_object.Rds`, `GSM8836087_unsorted_mouse_brain/CellRanger/fragments.tsv.gz`, per-cluster narrowPeak files (`unsorted_cluster_0_peaks.narrowPeak`, `unsorted_cluster_1_peaks.narrowPeak`), RPGC bigWigs (`cluster0_RPGC.bw`, `cluster1_RPGC.bw`)
- **Sorted GFP+ brain G4 scCUT&Tag (GSE291468, GSM8836086)** — `GSM8836086_GFPpos_Seurat_object.Rds`, `GSM8836086_GFP_sorted_mouse_brain/CellRanger/peaks.bed`, `GSM8836086_GFP_sorted_mousebrain.rpgc.bw`
- **Zeisel et al. neuronal scRNA-seq** — `scRNA_Seq-Zeisel_et_al-neuron.Rds`
- **PQS scores** — `PQS_scores.mm10.bw`
- **Mouse gene annotation** — EnsDb.Mmusculus.v79

Produced by `fig4/fig4.R`, with preprocessing in `fig4/map_unsorted.R` and `data/Zeisel2018/processing_neuron_data.R`.

## Panels

**Panel S1 — QC violin plots for unsorted brain.** The 3,636 cells comprised cluster 0 (n = 2,673) and cluster 1 (n = 963). Median values for clusters 0/1 were 326/383 detected peaks, 654/773 peak counts, 266/323 TSS fragments, 7/10 mitochondrial fragments and FRiP 0.161/0.208. Crossbars denote medians and points means; no inferential test was performed. Output: `panel_S1_qc_unsorted.pdf`.

**Panel A — Query and reference UMAPs.** The unsorted dataset contains 3,636 cells (cluster 0, 2,673; cluster 1, 963), and the GFP+ dataset contains 2,402 cells (clusters 0-3: 1,102, 513, 402 and 385). Output: `panel_A_umaps.pdf`.

**Panel B — Query-to-reference label transfer.** Mapping used 21,808 common GA genes. Predicted unsorted labels 0/1 by GFP+ cluster were 1,051/51, 369/144, 397/5 and 370/15 for clusters 0-3. Inclusive peak-set counts were GFP+ 50,399, unsorted cluster 0 26,360 and cluster 1 15,627; pairwise overlaps were 21,626, 10,423 and 9,273, with 9,010 shared by all three. Median prediction scores for GFP+ clusters 0-3 were 0.847, 0.793, 0.925 and 0.856. The boxplot annotation is one global Kruskal-Wallis test (chi-squared = 218.12, df = 3, nominal P < 2.2 x 10^-16), not pairwise Wilcoxon tests; no multiplicity adjustment was applied. Outputs: panel composite, component PDFs and `panel_B_venn_counts.csv`.

**Panel C — G4 occupancy heatmaps and average profiles.** Summit-centered windows are classified as 16,761 cluster0-only, 9,344 shared and 6,088 cluster1-only regions. Signal is shown +/-3 kb for cluster0, cluster1, GFP+ and PQS tracks. Heatmaps use the 99th percentile per track, at most 100 displayed rows per category and cluster-signal ordering; profiles use `show_error=TRUE`. No inferential test is performed. Outputs: `panel_C_*.pdf`.

**Panel D — enrichR GO analysis.** Excluding peaks found in the other unsorted cluster or GFP+ set leaves 4,509 cluster0 and 4,952 cluster1 peaks. Strict filtering above each set's 90th MACS2-score percentile retains 401 and 472 peaks; nearest-gene annotation within 3 kb of a TSS yields 93 unique genes per set. The heatmap contains all 24 CNS-related terms with nominal P < 0.05 (18 cluster0, 6 cluster1), rather than a top-ten or adjusted-P selection. No adjusted P-values are used in this panel. Outputs: `panel_D_enrichR_heatmap.pdf`, `panel_D_enrichR_GO_terms.csv`.

**Panel E — Genome browser tracks.** Three representative cluster1-specific G4 peaks at genes linked to GO term *primary neural tube formation*: *Cc2d1a*, *Lias*, *Prkacb*. Each window: gene TSS ± 5 kb. Tracks shown: cluster0 RPGC bigWig, cluster1 RPGC bigWig (common scale). Plotted with Signac `CoveragePlot()` (annotation=TRUE, peaks=FALSE, show.bulk=FALSE, bigwig.type="coverage"). Tools: Signac `CoveragePlot()`, rtracklayer. Outputs: `panel_E_genome_track_Cc2d1a.pdf`, `panel_E_genome_track_Lias.pdf`, `panel_E_genome_track_Prkacb.pdf`.

**Panel F — Neuron integration UMAP.** Joint scBridge embedding contains all 963 unsorted cluster-1 G4 cells plus the Zeisel neuronal reference. Four views show neuron type, modality, predictions and reliability. The reference-cell count is not recoverable from the bundled PDF alone. Output: `panel_F_neuron_integration.pdf`.

**Panel G — Neuron-type treemap.** The treemap contains all 963 G4 cells: 446 reliable assignments and 517 unreliable. Reliable counts were amygdala 54, cerebellum 50, olfactory 49, cortex 38, striatum-ventral 36, hippocampus 33, thalamus 31, midbrain-ventral 28, enteric 28, striatum-dorsal 24, hypothalamus 22, medulla 15, pons 13, midbrain-dorsal 10, dorsal-root ganglia 8 and sympathetic 7. Area and color encode log2(count). Output: `panel_G_neuron_treemap.pdf`.
