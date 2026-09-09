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

**Panel B — UMAP of the mESC-MEF mixture.** LSI dimensionality reduction → UMAP (n = 10 dimensions, local connectivity = 15, min dist = 0.5) and SNN clustering (resolution 0.1) of 1,621 GSM8836088 cells; cells are colored as MEF/cluster 0 (n = 1,198) or mESC/cluster 1 (n = 423). Tool: Seurat `DimPlot()`. Output: `panel_B_umap.pdf`.

**Panel C — QC violin plots and PQS overlap bar plot.** Seven per-cell metrics computed with Signac: nFeature_peaks, nCount_peaks, TSS_fragments (from `TSSEnrichment()`), mitochondrial fragments, FRiP (peak_region_fragments / passed_filters), TSS_enrichment (`TSSEnrichment(fast=TRUE)`), nucleosome_signal (`NucleosomeSignal()`). Plotted as violin plots per cluster with median (crossbar) and mean (point) overlaid. Median values for MEF/mESC were 839/1,660 detected peaks, 1,660.5/3,600 peak counts, 701/1,872 TSS fragments, 369/1,366 mitochondrial fragments, and FRiP 0.300/0.276. The bar plot shows PQS overlap for 9,134/10,408 MEF-specific peaks (87.8%), 12,672/13,586 mESC-specific peaks (93.3%) and 19,579/20,480 shared peaks (95.6%), using all PQS sites in `data/pqsfinder/PQS_scores.mm10.bed` (scores 20–367). A regioneR permutation test (1,000 permutations, chr1) on the union of all 44,398 MEF/ESC peaks versus the same PQS set gave 11,942 observed overlaps versus 5,319 expected (2.3-fold enrichment, z = 66.1, Monte Carlo p = 0.001), shown as the plot subtitle. Tool: Seurat `VlnPlot()`, ggplot2. Outputs: `panel_C_qc.pdf` (7 plots, 2×4 grid), `panel_C_bar_pct_peaks_overlap_PQS_MEF_ESC_scG4.pdf`.

**Panel D — Peak-overlap Euler diagrams.** Input sets contained 31,019 cluster-0, 34,152 cluster-1, 34,674 bulk-mESC and 22,857 bulk-MEF peaks. For the mESC comparison, the displayed region counts were 8,833 cluster0-only, 6,680 cluster1-only, 14,407 bulk-only, 8,510 cluster0+cluster1, 1,305 cluster0+bulk, 6,591 cluster1+bulk and 12,371 shared by all three. For the MEF comparison they were 6,730, 12,069, 1,200, 3,834, 3,408, 1,202 and 17,047, respectively. These are interval-anchored overlap categories exported by the script; no significance test was performed. Output: `panel_D_venn.pdf`, `panel_D_venn_counts.csv`.

**Panel E — |log2FC| > 2 filtered correlation heatmap.** Of 44,474 candidate regions, 10,825 with |log2FC(cluster1/cluster0)| > 2 were retained. Spearman/Pearson correlations were 0.647/0.741 for cluster0 versus bulk MEF, 0.666/0.731 for cluster1 versus bulk mESC, and -0.242/-0.050 between cluster0 and cluster1. Correlations were computed with `cor()` without P-values. Tools: wigglescout, GenomicRanges, ComplexHeatmap. Outputs: `panel_E.pdf`, `panel_E_spearman_correlation.csv`, `panel_E_pearson_correlation.csv`.

**Panel F — G4 occupancy heatmaps and average profiles.** Summit-centered regions (±250 bp from narrowPeak summit) were classified as 10,408 cluster0-only, 20,480 shared and 13,586 cluster1-only regions. The corresponding exported sets contained 1,452, 8,953 and 2,318 unique nearest genes. For each class, RPGC signal (±3 kb around center) is shown for seven tracks. Heatmaps use the 99th percentile per track, at most 100 displayed rows per category, and rows sorted by cluster-track signal; profiles show the mean with the error representation produced by `show_error=TRUE`. No significance test was performed. Outputs: `panel_F_*.pdf`, `panel_F_merged_peaks.bed`, `panel_F_gene_sets.tsv`, `panel_F_GO_Biological_Process_2023.tsv`.

**Panel G — enrichR GO Biological Process heatmap.** Peak-associated genes from Panel F were queried against Enrichr GO_Biological_Process_2023. The heatmap displays the top 10 terms per category by nominal P-value as -log10(nominal P), not adjusted P. Across all returned terms, nominal/adjusted-P<0.05 counts were 306/0 (cluster0-only), 1,134/548 (shared) and 283/11 (cluster1-only); thus none of the ten displayed cluster0-only terms passed adjusted P<0.05. Outputs: `panel_G_enrichR_heatmap.pdf`, `panel_G_enrichR_GO_terms.csv`.

**Panel H — Genome browser tracks.** Two loci: *Lin28a* (mESC marker) and *Cdhr3* (MEF marker). Fragment coverage from GSM8836088 Cell Ranger fragments (`GSM8836088_fragments.tsv.gz`) plotted with Signac `CoveragePlot()` (peaks=TRUE, features=TRUE, flank=10kb). Tool: Signac `CoveragePlot()` using EnsDb.Mmusculus.v79 annotation. Outputs: `panel_H_Lin28a.pdf`, `panel_H_Cdhr3.pdf`.

## Supplementary panels

**Panel S1 — Spearman and Pearson correlation heatmaps.** Correlations among four bigWig tracks were calculated at 23,994 cluster-specific peaks. Spearman/Pearson correlations were 0.651/0.670 for cluster0 versus bulk MEF, 0.552/0.579 for cluster1 versus bulk mESC and 0.121/0.452 between cluster0 and cluster1. `cor()` was used without P-values. Outputs: `panel_S1.pdf`, correlation CSVs and `panel_S1_marker_regions.mm10.bed`.

**Panel S2 — G4 vs ATAC scatter correlations.** Pearson r (n regions) was 0.382 (10,407) for cluster0-only versus MEF ATAC, 0.260 (13,576) for cluster1-only versus mESC ATAC, 0.893 (20,480) for shared cluster0 versus cluster1, 0.205 and 0.112 (20,471 each) for shared cluster0/MEF-ATAC and cluster1/mESC-ATAC, and 0.143 (44,454) for all cluster1 regions versus mESC ATAC. No P-values were calculated. Outputs: `panel_S2_G4_ATAC_correlation.pdf`, `panel_S2_G4_ATAC_correlation_pearson.csv`.

**Panel S3 — PQS-highlighted correlations.** Same as S2, but with loci overlapping PQS sites (score ≥50 from `PQS_scores.mm10.bed`) highlighted in red. Additional panel shows merged peak sets (pseudobulk cluster peaks ∪ matching ATAC peaks) with PQS overlap highlighted. Tools: wigglescout, GenomicRanges, ggrastr. Outputs: `panel_S3_G4_ATAC_correlation_PQS.pdf`, `panel_S3_alt_merged_peaks.pdf`.

**Panel S4 — Density scatterplots.** Eight ggscatterhist plots comparing cluster peaks vs matching ATAC peaks (log2 axes) with marginal density distributions. Merged peaks colored by PQS overlap (score ≥50). Tool: ggpubr `ggscatterhist()`. Output: `panel_S4_ggscatterhist.pdf`.
