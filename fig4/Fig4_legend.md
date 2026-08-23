# Figure 4 — Legend and Derivation Notes

**Figure 4: Integration of the unsorted (GFP−/GFP+) brain G4 scCUT&Tag dataset with neuronal scRNA-seq aids identification of neuronal cell types.**

## Data sources

- **Unsorted mouse-brain G4 scCUT&Tag (GSE291468, GSM8836087)** — `GSM8836087_unsorted_Seurat_object.Rds`, `GSM8836087_unsorted_mouse_brain/CellRanger/fragments.tsv.gz`, per-cluster narrowPeak files (`unsorted_cluster_0_peaks.narrowPeak`, `unsorted_cluster_1_peaks.narrowPeak`), RPGC bigWigs (`cluster0_RPGC.bw`, `cluster1_RPGC.bw`)
- **Sorted GFP+ brain G4 scCUT&Tag (GSE291468, GSM8836086)** — `GSM8836086_GFPpos_Seurat_object.Rds`, `GSM8836086_GFP_sorted_mouse_brain/CellRanger/peaks.bed`, `GSM8836086_GFP_sorted_mousebrain.rpgc.bw`
- **Zeisel et al. neuronal scRNA-seq** — `scRNA_Seq-Zeisel_et_al-neuron.Rds`
- **PQS scores** — `PQS_scores.mm10.bw`
- **Mouse gene annotation** — EnsDb.Mmusculus.v79

Produced by `fig4/fig4.R`, with preprocessing in `fig4/map_unsorted.R` and `fig4/processing_neuron_data.R`.

## Panels

**Panel S1 — QC violin plots for unsorted brain.** Seven per-cluster metrics computed with Signac: nFeature_peaks, nCount_peaks, TSS_fragments (`TSSEnrichment(fast=TRUE)`), mitochondrial fragments, FRiP (`FRiP()` with total_fragments from `CountFragments()`), TSS_enrichment, nucleosome_signal (`NucleosomeSignal()`). Plotted as violin plots per Seurat cluster with median (crossbar) and mean (point) overlaid. Tool: Seurat `VlnPlot()`. Output: `panel_S1_qc_unsorted.pdf`.

**Panel A — Query and reference UMAPs.** Two UMAPs side-by-side: (left) unsorted brain G4 cells (GSM8836087) colored by two Seurat clusters (cluster 0, cluster 1); (right) sorted GFP+ G4 cells (GSM8836086) colored by four Seurat clusters. Tools: Seurat `DimPlot(label=TRUE, repel=TRUE, pt.size=0.35)`. Output: `panel_A_umaps.pdf`.

**Panel B — Query-to-reference label transfer.** Four-panel composite: **B1** (top-left): UMAP of GFP+ cells colored by predicted unsorted-cluster label from Seurat `MapQuery()` (predicted.seurat_clusters). **B2** (top-right): UMAP of GFP+ cells colored by prediction score (predicted.seurat_clusters.score, RdBu reversed gradient). **B3** (bottom-left): Venn/Euler diagram of peak overlap between GFP+ peak set (`peaks.bed`) and two unsorted cluster peak sets (narrowPeak files), computed with `bedscout::plot_euler()` (input="union"). **B4** (bottom-right): boxplots of prediction scores per GFP+ cluster with group-wise significance (Wilcoxon). Label transfer: Seurat `FindTransferAnchors(reference=unsorted, query=sorted, reduction="cca", features=common_ga_genes)`, then `MapQuery()` with reference.reduction="umap". Tools: Seurat, bedscout, eulerr, ggplot2. Outputs: `panel_B_projection_venn_boxplots.pdf`, `panel_B_umap1.pdf`, `panel_B_umap2.pdf`, `panel_B_venn.pdf`, `panel_B_venn_counts.csv`.

**Panel C — G4 occupancy heatmaps and average profiles.** Peak regions classified as cluster0-only, shared (both), cluster1-only from summit-centered narrowPeak windows (±500 bp from summit) using `GenomicRanges::findOverlaps()` and `reduce()`. For each class, RPGC signal (±3 kb around center) shown as average-profile curves and heatmaps for four bigWigs: GSM8836087_cluster0, GSM8836087_cluster1, GSM8836086_GFP+, PQS_scores. Heatmaps computed with `wigglescout::plot_bw_heatmap()` (zmax=99th percentile per track, capped at 100 rows per category, rows sorted by cluster track signal). Profiles computed with `wigglescout::plot_bw_profile()` (show_error=TRUE). Composite layout: profile on top (62.5% height), three heatmaps below (heights proportional to category sizes). Tools: wigglescout, GenomicRanges, cowplot. Outputs: `panel_C_cluster0.pdf`, `panel_C_cluster1.pdf`, `panel_C_GFP+.pdf`, `panel_C_PQS_score.pdf`, `panel_C_heatmap_*.pdf`, `panel_C_profile_*.pdf`.

**Panel D — enrichR GO analysis.** Unique cluster-specific peak sets defined as: cluster0 peaks absent from both cluster1 and GFP+ peak universe (bedtools-style unique), cluster1 peaks absent from both cluster0 and GFP+. Peaks filtered to top 10% by MACS2 score (V7 column), annotated to nearest genes with ChIPseeker `annotatePeak()` (tssRegion=±3kb, TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb=org.Mm.eg.db). Gene sets queried against Enrichr GO_Biological_Process_2023 library with `enrichR::enrichr()`. Heatmap shows CNS-related GO terms (top 10 by P-value per category) as −log10 P. Tools: ChIPseeker, enrichR, ComplexHeatmap. Outputs: `panel_D_enrichR_heatmap.pdf`, `panel_D_enrichR_GO_terms.csv`.

**Panel E — Genome browser tracks.** Three representative cluster1-specific G4 peaks at genes linked to GO term *primary neural tube formation*: *Cc2d1a*, *Lias*, *Prkacb*. Each window: gene TSS ± 5 kb. Tracks shown: cluster0 RPGC bigWig, cluster1 RPGC bigWig (common scale). Plotted with Signac `CoveragePlot()` (annotation=TRUE, peaks=FALSE, show.bulk=FALSE, bigwig.type="coverage"). Tools: Signac `CoveragePlot()`, rtracklayer. Outputs: `panel_E_genome_track_Cc2d1a.pdf`, `panel_E_genome_track_Lias.pdf`, `panel_E_genome_track_Prkacb.pdf`.

**Panel F — Neuron integration UMAP.** Joint embedding (scBridge) of unsorted cluster 1 G4 cells (GSM8836087) with Zeisel et al. neuronal scRNA-seq reference. Four views: neuron type, modality, predicted labels, integration reliability. Precomputed scBridge output copied from `results/scBridge/output/unsorted_cl1_Zeisel/Seurat_UMAPs-unsorted_cl1_int.pdf`. Output: `panel_F_neuron_integration.pdf`.

**Panel G — Neuron-type treemap.** Treemap of predicted neuron types from scBridge label imputation onto Zeisel reference (reliable cells only). Area and color represent log2-transformed predicted cell count per neuron type. Precomputed scBridge output copied from `results/scBridge/output/unsorted_cl1_Zeisel/tree_plot-unsorted_cl1_predictions.pdf`. Output: `panel_G_neuron_treemap.pdf`.