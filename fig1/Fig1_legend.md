# Figure 1 — Legend and Derivation Notes

**Figure 1: G4 profiling by scCUT&Tag separates different cell types.**

## Data sources used throughout Figure 1

All single-cell G4/CUT&Tag data come from **GSE291468**:

- **mESC-MEF single-cell mixture (GSM8836088)** — processed Seurat object, Cell Ranger fragment file, per-cluster MACS2 peak sets and RPGC-normalized cluster signal tracks.
- **Bulk mESC G4 CUT&Tag (GSM8836082, GSM8836083)** — two replicates; broadPeak peak calls and RPGC bigWigs.
- **Bulk MEF/3T3 G4 CUT&Tag (GSM8836084, GSM8836085)** — two replicates; broadPeak peak calls and RPGC bigWigs.
- **mESC ATAC-seq (GSE149080, GSM4661960)** — RPGC track and narrowPeak calls, used for cross-assay comparison.
- **3T3 ATAC-seq (GSE211123, GSM6451000)** — RPGC track and narrowPeak calls.
- **PQS (G-quadruplex-predicting sequence) scores** — genome-wide pqsfinder predictions on mm10 (in-house track; a minimum-score BED and a score bigWig are used).

Panels are produced by the R script `fig1/fig1.R` (panel-by-panel) with supporting preprocessing in `fig1/seurat_workflow.R`.

## Panels

**Panel B — UMAP of the mESC-MEF mixture.** Dimensionality reduction (LSI → UMAP) and SNN clustering (resolution 0.1) of the GSM8836088 single-cell object; cells are colored by the two Seurat clusters (MEF = cluster 0, mESC = cluster 1). Tool: Seurat `DimPlot`.

**Panel C — QC violin plots.** Distribution of per-cell quality metrics (number of peak features, peak counts, TSS fragments, mitochondrial fragments, fraction of reads in peaks [FRiP], TSS enrichment, nucleosome signal) split by cluster, computed with Seurat/Signac from the GSM8836088 peak matrix and fragment file. Medians and means are overlaid.

**Panel D — Peak-overlap Euler/Venn diagrams.** Pairwise overlaps of the two mESC-MEF cluster peak sets (GSM8836088 cluster-specific narrowPeak calls) with the bulk mESC (GSM8836082) and bulk MEF (GSM8836084) broadPeak sets, drawn as two Euler diagrams with overlap counts. Tools: `GenomicRanges` and `bedscout::plot_euler`.

**Panel E — PCA of marker regions.** Differentially enriched marker peaks are identified with Seurat `FindAllMarkers` (logistic-regression test with fragment count as latent variable) on the GSM8836088 object. Principal-component analysis is then performed on RPGC-normalized signal at these marker regions across six tracks: the two single-cell clusters (GSM8836088) and the two bulk replicates of mESC (GSM8836082, GSM8836083) and MEF (GSM8836084, GSM8836085). Cluster 0 and cluster 1 group with the MEF and mESC bulk replicates, respectively.

**Panel F — G4 occupancy heatmaps and average profiles.** Peak regions are classified into cluster-0-only, shared, and cluster-1-only sets from summit-centered narrowPeak windows (GSM8836088). For each class, G4 occupancy (±3 kb around peak center, RPGC normalized) is shown as average-profile curves and heatmaps for: bulk mESC/MEF G4 CUT&Tag (GSM8836082, GSM8836084), single-cell cluster tracks (GSM8836088), mESC/3T3 ATAC-seq (GSE149080 GSM4661960, GSE211123 GSM6451000), and PQS score. Heatmap rows are sorted by the matching cluster signal. Tools: `wigglescout::plot_bw_heatmap` and `wigglescout::plot_bw_profile`.

**Panel G — enrichR GO Biological Process heatmap.** Top significant GO-BP terms (top 10 per category) for the cluster-0-only, shared, and cluster-1-only peak sets of Panel F, queried against the Enrichr GO_Biological_Process_2023 library; heatmap shows −log10 P-values.

**Panel H — Genome browser tracks.** Representative cluster-specific G4 loci at the marker genes *Lin28a* (mESC) and *Cdhr3* (MEF), showing single-cell fragment coverage, the cluster-specific peak, and bulk tracks over ±10 kb. Tool: Signac `CoveragePlot` from the GSM8836088 fragment file.

## Supplementary panels produced alongside Figure 1

**Panel S1 — Spearman correlation heatmap.** Spearman correlation of RPGC signal among the two cluster tracks and the mESC/MEF bulk tracks at positive marker regions (Wilcoxon markers, log2FC ≥ 1.25; GSM8836088 and bulk replicates).

**Panel S2 — G4 vs ATAC signal scatter correlations.** Log2 RPGC scatter plots of single-cell G4 signal vs the matched ATAC-seq track (mESC/3T3) for each peak category (cluster-0-only, cluster-1-only, shared), with Pearson r values. Tools: `wigglescout::plot_bw_loci_scatter`.

**Panel S3 — PQS-highlighted correlations.** The same G4-vs-ATAC scatters with loci overlapping a PQS site (minimum score 50) highlighted, plus an alternative rendering on merged peak sets (pseudobulk cluster peaks ∪ matching ATAC peaks).

**Panel S4 — ggscatterhist.** Density scatterplots with marginal distributions (log2 axes) comparing cluster peaks vs matching ATAC peaks, and merged peaks colored by PQS overlap. Tool: `ggpubr::ggscatterhist`.
