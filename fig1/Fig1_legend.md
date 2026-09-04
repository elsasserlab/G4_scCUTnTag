# Figure 1 — Legend and Derivation Notes

**Figure 1: G4 profiling by scCUT&Tag separates different cell types.**

## Data sources

All data from **GSE291468**:

- **mESC-MEF single-cell mixture (GSM8836088)** — processed Seurat object, Cell Ranger fragment file (`GSM8836088_fragments.tsv.gz`), per-cluster MACS2 narrowPeak sets, RPGC-normalized bigWig tracks
- **Bulk mESC G4 CUT&Tag (GSM8836082, GSM8836083)** — broadPeak calls and RPGC bigWigs
- **Bulk MEF/3T3 G4 CUT&Tag (GSM8836084, GSM8836085)** — broadPeak calls and RPGC bigWigs
- **mESC ATAC-seq (GSE149080, GSM4661960)** — RPGC bigWig and narrowPeak calls
- **3T3 ATAC-seq (GSE211123, GSM6451000)** — RPGC bigWig and narrowPeak calls
- **PQS scores** — pqsfinder predictions on mm10 (`PQS_scores.mm10.bed`, `PQS_scores.mm10.bw`)

Produced by `fig1/fig1.R` with preprocessing in `fig1/seurat_workflow.R`.

## Panels

**Panel B — UMAP of the mESC-MEF mixture.** LSI dimensionality reduction → UMAP (n = 10 dimensions, local connectivity = 15, min dist = 0.5) and SNN clustering (resolution 0.1) of GSM8836088 Seurat object; cells colored by two clusters (MEF = cluster 0, mESC = cluster 1). Tool: Seurat `DimPlot()`. Output: `panel_B_umap.pdf`.

**Panel C — QC violin plots and PQS overlap bar plot.** Seven per-cell metrics computed with Signac: nFeature_peaks, nCount_peaks, TSS_fragments (from `TSSEnrichment()`), mitochondrial fragments, FRiP (peak_region_fragments / passed_filters), TSS_enrichment (`TSSEnrichment(fast=TRUE)`), nucleosome_signal (`NucleosomeSignal()`). Plotted as violin plots per cluster with median (crossbar) and mean (point) overlaid. Bar plot shows % of cluster-specific peaks (cl0_only + cl1_only) overlapping ≥1 PQS site (score ≥50 from `PQS_scores.mm10.bed`). Tool: Seurat `VlnPlot()`, ggplot2. Outputs: `panel_C_qc.pdf` (7 plots, 2×4 grid), `panel_C_bar_pct_peaks_overlap_PQS_MEF_ESC_scG4.pdf`.

**Panel D — Peak-overlap Euler diagrams.** Two diagrams: (left) cluster0 vs cluster1 vs bulk mESC (GSM8836082 broadPeak); (right) cluster0 vs cluster1 vs bulk MEF (GSM8836084 broadPeak). Peaks imported with `rtracklayer::import()`, overlaps computed with `bedscout::plot_euler()`. Output: `panel_D_venn.pdf`, `panel_D_venn_counts.csv`.

**Panel E — |log2FC| > 2 filtered correlation heatmap.** Spearman and Pearson correlations among four bigWig tracks (cluster0, cluster1, GSM8836082, GSM8836084) at peaks with |log2FC(cluster1/cluster0)| > 2 (including shared peaks; 10,825 peaks total). Signal extracted with `wigglescout::bw_loci()`, correlations computed with `cor(method="spearman")` and `cor(method="pearson")`, heatmaps rendered with ComplexHeatmap `Heatmap()`. Tools: wigglescout, GenomicRanges, ComplexHeatmap. Outputs: `panel_E.pdf` (Spearman + Pearson), `panel_E_spearman_correlation.csv`, `panel_E_pearson_correlation.csv`.

**Panel F — G4 occupancy heatmaps and average profiles.** Summit-centered regions (±250 bp from narrowPeak summit) classified as cluster0-only, shared (both), cluster1-only using `GenomicRanges::findOverlaps()` and `reduce()`. For each class, RPGC signal (±3 kb around center) shown as average-profile curves and heatmaps for seven bigWigs: GSM8836082 (bulk mESC), GSM8836084 (bulk MEF), GSM8836088_cluster0, GSM8836088_cluster1, GSM4661960 (mESC ATAC), GSM6451000 (MEF ATAC), PQS_scores.mm10.bw. Heatmaps computed with `wigglescout::plot_bw_heatmap()` (zmax = 99th percentile per track, capped at 100 rows per category, rows sorted by cluster track signal). Profiles computed with `wigglescout::plot_bw_profile()` (show_error=TRUE). Tools: wigglescout, GenomicRanges, cowplot. Outputs: `panel_F_*.pdf` (7 composite PDFs), `panel_F_merged_peaks.bed`, `panel_F_gene_sets.tsv`, `panel_F_GO_Biological_Process_2023.tsv`.

**Panel G — enrichR GO Biological Process heatmap.** Peak-associated genes from Panel F (cluster0-only, shared, cluster1-only) queried against Enrichr GO_Biological_Process_2023 library with `enrichR::enrichr()`. Top 10 terms per category by adjusted P-value shown as heatmap (−log10 P). Tool: enrichR, ComplexHeatmap. Outputs: `panel_G_enrichR_heatmap.pdf`, `panel_G_enrichR_GO_terms.csv`.

**Panel H — Genome browser tracks.** Two loci: *Lin28a* (mESC marker) and *Cdhr3* (MEF marker). Fragment coverage from GSM8836088 Cell Ranger fragments (`GSM8836088_fragments.tsv.gz`) plotted with Signac `CoveragePlot()` (peaks=TRUE, features=TRUE, flank=10kb). Tool: Signac `CoveragePlot()` using EnsDb.Mmusculus.v79 annotation. Outputs: `panel_H_Lin28a.pdf`, `panel_H_Cdhr3.pdf`.

## Supplementary panels

**Panel S1 — Spearman and Pearson correlation heatmaps.** Correlations among four bigWig tracks (cluster0, cluster1, GSM8836082, GSM8836084) at cluster-specific peaks (cl0_only + cl1_only; 23,994 total). Signal extracted with `wigglescout::bw_loci()`, correlations computed with `cor(method="spearman")` and `cor(method="pearson")`, heatmaps rendered with ComplexHeatmap `Heatmap()`. Outputs: `panel_S1.pdf` (Spearman + Pearson), `panel_S1_spearman_correlation.csv`, `panel_S1_pearson_correlation.csv`, `panel_S1_marker_regions.mm10.bed`.

**Panel S2 — G4 vs ATAC scatter correlations.** Three peak categories (cluster0-only, cluster1-only, shared) from Panel F. For each category, log2 RPGC scatter plot comparing G4 signal (cluster0 or cluster1 bigWig) vs matched ATAC signal (GSM4661960 for cluster1/mESC, GSM6451000 for cluster0/MEF). Signal extracted with `wigglescout::bw_loci()`, scatter plotted with `wigglescout::plot_bw_loci_scatter()` (rasterized at 150 dpi with ggrastr). Pearson r computed per category. Tools: wigglescout, ggrastr. Outputs: `panel_S2_G4_ATAC_correlation.pdf`, `panel_S2_G4_ATAC_correlation_pearson.csv`.

**Panel S3 — PQS-highlighted correlations.** Same as S2, but with loci overlapping PQS sites (score ≥50 from `PQS_scores.mm10.bed`) highlighted in red. Additional panel shows merged peak sets (pseudobulk cluster peaks ∪ matching ATAC peaks) with PQS overlap highlighted. Tools: wigglescout, GenomicRanges, ggrastr. Outputs: `panel_S3_G4_ATAC_correlation_PQS.pdf`, `panel_S3_alt_merged_peaks.pdf`.

**Panel S4 — Density scatterplots.** Eight ggscatterhist plots comparing cluster peaks vs matching ATAC peaks (log2 axes) with marginal density distributions. Merged peaks colored by PQS overlap (score ≥50). Tool: ggpubr `ggscatterhist()`. Output: `panel_S4_ggscatterhist.pdf`.