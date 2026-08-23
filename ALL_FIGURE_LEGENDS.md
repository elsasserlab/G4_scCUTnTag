# All Figure Legends (combined)

Generated from fig1-4 and rev legend sources. Individual per-figure
copies are kept in their figure directories.


---

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

---

# Figure 2 — Legend and Derivation Notes

**Figure 2: Integration of G4 scCUT&Tag data with scRNA-seq using canonical correlation analysis (CCA).**

## Data sources used throughout Figure 2

- **Sorted GFP+ brain G4 scCUT&Tag (GSE291468, GSM8836086)** — single-cell peak/gene-activity matrices and processed Seurat integration object.
- **Brain scRNA-seq reference (GSE163484, Bartosovic et al.)** — processed reference object used for label transfer.
- **scBridge co-embedding and label-transfer outputs** — precomputed joint embedding (UMAP coordinates), per-cell predicted cell types and integration reliability for the GFP+ G4 cells. These de-novo outputs (not deposited on GEO) are read from `results/scBridge/output/GFPsorted_Bartosovic/`.

Panels are produced by the R script `fig2/fig2.R`, with QC panel generated by `fig2/fig2_S1_QC.R`, and preprocessing in `fig2/seurat_integration.R`.

## Panels

**Panel S1 — QC violin plots for GFP+ sorted mouse brain.** Quality control metrics for the GFP+ sorted brain G4 scCUT&Tag dataset (GSM8836086) shown as violin plots per Seurat cluster: nFeature (peaks), nCount (peaks), TSS fragments, mitochondrial fragments, fraction of reads in peaks (FRiP), TSS enrichment, and nucleosome signal. Median (black bar) and mean (black dot) are indicated. Y-axis limits: FRiP 0-0.15, TSS enrichment 0-10, nucleosome signal 0-3. Tools: Signac `TSSEnrichment()`, `NucleosomeSignal()`, `FRiP()`, and `VlnPlot()`.

**Panel B — Co-embedded UMAPs.** Joint embedding of the GFP+ G4 cells (GSM8836086) and the scRNA-seq reference (GSE163484) obtained with scBridge. Four views share the same coordinates: data modality (G4 scCUT&Tag vs scRNA-seq), reference cell-type labels, transferred (predicted) labels for the G4 cells, and integration reliability score. Tools: ggplot2 on the precomputed scBridge coordinates (`umap_coembedded.csv`, `scbridge_predictions.csv`, `scbridge_reliability.csv`).

**Panel C — Prediction-score heatmap.** Heatmap of the class-transfer prediction scores for the GFP+ G4 cells across the reference cell types (Seurat label transfer on gene-activity; the `max` row is removed), colored with a viridis scale. Tool: `ComplexHeatmap`.

**Panel D — Cluster/prediction correspondence.** D1: UMAP of the GFP+ G4 cells colored by Seurat cluster. D2: symmetric overlap-score heatmap quantifying the agreement between Seurat clusters and predicted cell types (fraction-of-cluster vs fraction-of-cell-type overlap, cell types abbreviated as in the manuscript: AST, MOL, OPC, COP-NFOL, OEC, VEC, VLMC, Pericytes).

**Panel E — scRNA-seq marker analysis.** Reference UMAP colored by cell type plus violin plots of the top differentially expressed marker gene per cell type, from `FindAllMarkers` on the GSE163484 scRNA-seq reference.

**Panel F — G4/RNA co-enrichment feature plots.** For the same marker genes, UMAP feature plots on the shared embedding showing the G4 gene-activity score (GA assay, GSM8836086) and the normalized scRNA-seq expression (GSE163484), demonstrating co-enrichment of G4 activity and expression within the corresponding cell types.

---

# Figure 3 — Legend and Derivation Notes

**Figure 3: Reliability-based annotation enhances cell-type assignment of G4 scCUT&Tag data and identifies astrocyte-specific G4 signals.**

## Data sources used throughout Figure 3

- **Sorted GFP+ brain G4 scCUT&Tag (GSE291468, GSM8836086)** — gene-activity (GA) and peak assays, predicted-astrocyte and predicted-non-astrocyte RPGC signal tracks, and label-transfer prediction scores.
- **Brain scRNA-seq reference (GSE163484, Bartosovic et al.)** — processed reference object for co-embedding and marker expression.
- **scBridge outputs** — precomputed joint embedding, predictions and reliability (`results/scBridge/output/GFPsorted_Bartosovic/`), used for the co-embedding and AST/non-AST definition.
- **Cicero co-accessibility** — co-accessibility model on the G4 peaks with candidate regulatory elements from a mouse-brain cCRE resource (Li et al., cCRE intervals with H3K27ac signal).
- **Ensembl mouse gene annotation (EnsDb, mm10)** — for TSS windows in coverage panels.

Panels are produced by the R script `fig3/fig3.R`, with supporting analyses in `fig3/cicero.R` and `fig3/scbridge_input.R`.

## Panels

**Panel A — scBridge joint-embedding UMAP.** Four views on a single co-embedding of the GFP+ G4 cells (GSM8836086) and the scRNA-seq reference (GSE163484): reference cell types, data modality, transferred (predicted) labels for the G4 cells, and integration reliability. Tools: ggplot2 on precomputed scBridge coordinates.

**Panel B — Astrocyte (AST) vs non-AST differential G4.** Left: UMAP of the GFP+ G4 cells colored by AST/non-AST status (prediction label "AST" vs all other reliable labels). Right: volcano plot of differential gene-activity scores between predicted AST and non-AST cells, computed with Seurat `FindMarkers` (logistic-regression test on the GA assay with fragment count as latent variable). Thresholds: |log2 fold change| > 0.5 and adjusted P < 0.05. Genes highlighted in the volcano (e.g. *Tnik*, *Pitpnc1*, *Pbx1*, *Nwd1*) are followed up in C and D.

**Panel C — Coverage at astrocyte marker genes.** Genome coverage (RPGC) at the TSS ± 5 kb of top astrocyte marker genes (plus the four highlighted genes), comparing the predicted-astrocyte and predicted-non-astrocyte RPGC tracks (GSM8836086). Tool: Signac `CoveragePlot` on bigWig tracks with a common scale.

**Panel D — Astrocyte-specific gene expression feature plots.** UMAP feature plots of normalized scRNA-seq expression for the four highlighted astrocyte genes on the scRNA-seq reference cells (GSE163484), on the same co-embedding coordinates as panel A.

**Panel E — Cicero co-accessibility browser tracks.** Candidate regulatory-element browser views connecting astrocyte-specific G4 peaks with mouse-brain cCREs (cCREs with H3K27ac signal above the 75th percentile), with G4 peaks, PQS sites and regulatory elements indicated; aggregated AST and non-AST G4 coverage are shown alongside. Derived from the Cicero co-accessibility analysis on the G4 peaks.

---

# Figure 4 — Legend and Derivation Notes

**Figure 4: Integration of the unsorted (GFP−/GFP+) brain G4 scCUT&Tag dataset with neuronal scRNA-seq aids identification of neuronal cell types.**

## Data sources used throughout Figure 4

- **Unsorted mouse-brain G4 scCUT&Tag (GSE291468, GSM8836087)** — Seurat object, Cell Ranger fragment file, per-cluster (cluster 0/1) MACS2 peak sets and cluster RPGC signal tracks.
- **Sorted GFP+ brain G4 scCUT&Tag (GSE291468, GSM8836086)** — Seurat object, fragment file, genome-wide peak set and RPGC signal track, used as the reference for label transfer.
- **Zeisel et al. neuronal scRNA-seq** — reference for the neuron-focused embedding (processed object).
- **PQS (G-quadruplex-predicting sequence) scores** — genome-wide pqsfinder predictions on mm10.
- **Ensembl mouse gene annotation (EnsDb, mm10)** — TSS windows for genome-browser panels.

Panels are produced by the R script `fig4/fig4.R`, with QC panel generated by `fig4/fig4_S1_QC.R`, and preprocessing in `fig4/map_unsorted.R` and `fig4/processing_neuron_data.R`.

## Panels

**Panel S1 — QC violin plots for unsorted mouse brain.** Quality control metrics for the unsorted brain G4 scCUT&Tag dataset (GSM8836087) shown as violin plots per Seurat cluster: nFeature (peaks), nCount (peaks), TSS fragments, mitochondrial fragments, fraction of reads in peaks (FRiP), TSS enrichment, and nucleosome signal. Median (black bar) and mean (black dot) are indicated. Y-axis limits: FRiP 0-0.20, TSS enrichment 0-10, nucleosome signal 0-3. Tools: Signac `TSSEnrichment()`, `NucleosomeSignal()`, `FRiP()`, and `VlnPlot()`.

**Panel A — Query and reference UMAPs.** UMAP of the unsorted brain G4 cells (GSM8836087) colored by the two Seurat clusters, and UMAP of the sorted GFP+ G4 cells (GSM8836086) colored by Seurat cluster.

**Panel B — Query-to-reference label transfer.** The unsorted dataset (query) is projected onto the GFP+ dataset (reference) using Seurat `FindTransferAnchors`/`MapQuery` on gene-activity scores. Shown are: UMAP of the GFP+ cells colored by the predicted unsorted-cluster label, UMAP of the prediction score, a Venn/Euler diagram of peak overlap between the GFP+ peak set and the two unsorted cluster peak sets (GSM8836086, GSM8836087), and boxplots of prediction scores per GFP+ cluster with group-wise significance (Wilcoxon).

**Panel C — G4 occupancy heatmaps and average profiles.** Peak regions are classified into cluster-0-only, shared, and cluster-1-only sets from summit-centered narrowPeak windows (GSM8836087). For each class, RPGC signal (±3 kb around peak center) is shown as average-profile curves and heatmaps for the unsorted cluster 0/1 tracks (GSM8836087), the GFP+ track (GSM8836086), and the PQS score. Tools: `wigglescout::plot_bw_heatmap` and `wigglescout::plot_bw_profile`.

**Panel D — enrichR GO analysis.** Unique cluster-specific peak sets (each cluster peak absent from the other cluster and from the GFP+ peak universe) are annotated to nearest genes (ChIPseeker, ±3 kb of TSS) and queried against Enrichr; the heatmap shows CNS-related GO Biological Process terms with their P-values for the two clusters.

**Panel E — Genome browser tracks.** Representative cluster-1-specific G4 peaks at genes linked to the GO term *primary neural tube formation* (*Cc2d1a*, *Lias*, *Prkacb*), showing the cluster 0/1 RPGC tracks (GSM8836087) over the TSS ± 5 kb. Tool: Signac `CoveragePlot` on bigWig tracks with a common scale.

**Panel F — Neuron integration UMAP.** Joint embedding (scBridge) of the unsorted cluster 1 G4 cells (GSM8836087) with the Zeisel et al. neuronal scRNA-seq reference, colored by neuron type, modality, predicted labels and integration reliability (precomputed output).

**Panel G — Neuron-type treemap.** Treemap of the predicted neuron types obtained by label imputation onto the Zeisel reference (reliable cells only, reliability > 0.90), with area/color representing the log2-transformed predicted cell count per neuron type.

---

# Revision Figures — Legend and Derivation Notes

The revision/rebuttal panels (F01–F23) answer reviewer questions about G4
signal reproducibility, functional annotation, independence from chromatin
accessibility, and QC of the neural datasets. Each panel was produced by an
R script under `rev/` (`rev_F*.R`), writing to `rev/outputs/`.

## Shared data sources

- **mESC-MEF single-cell G4 scCUT&Tag (GSE291468, GSM8836088)** — cluster-specific (cluster 0 = MEF, cluster 1 = mESC) MACS2 narrowPeak sets, pseudobulk cluster RPGC tracks, and the single-cell Seurat object.
- **Bulk G4 CUT&Tag (GSE291468)** — bulk mESC replicates (GSM8836082, GSM8836083) and bulk MEF/3T3 replicates (GSM8836084, GSM8836085) as broadPeak calls and RPGC bigWigs.
- **Sorted GFP+ brain G4 (GSE291468, GSM8836086)** — predicted-astrocyte / predicted-non-astrocyte RPGC tracks and label-transfer outputs.
- **Unsorted mouse-brain G4 (GSE291468, GSM8836087)** — Seurat object, cluster-specific peak sets and cluster RPGC tracks.
- **mESC ATAC-seq (GSE149080, GSM4661960)** and **3T3 ATAC-seq (GSE211123, GSM6451000)** — matched bulk accessibility peak sets and RPGC tracks.
- **Published 3T3 G4 CUT&Tag (GSE217860, GSM6729104–GSM6729106)** — three replicates used as an independent bulk benchmark.
- **Brain snATAC reference (GSE198467)** — clustered Seurat object with 15 named cell types.
- **Gene expression reference** — per-gene TPM table (mESC/MEF) and the mm10 GENCODE annotation GTF.
- **G4 references** — genome-wide pqsfinder PQS predictions (mm10; BED + score bigWig) and an experimentally validated eG4 catalog.

## Panels

**F01 — Genomic feature distribution (bulk vs scCUT&Tag).** Fraction of G4 peaks assigned to genomic features (promoter TSS ± 2 kb, 5′ UTR, 3′ UTR, CDS exon, intron, intergenic) for bulk-consensus vs single-cell pseudobulk peaks in mESC and MEF. Bulk consensus built from the two bulk replicates (GSM8836082/83 mESC; GSM8836084/85 MEF) with `bedscout::loci_consensus` (min 2 replicates); pseudobulk = GSM8836088 cluster peaks. Features derived from the mm10 GENCODE GTF; each peak classified by a midpoint priority cascade (GenomicRanges).

**F03 — Distance-to-TSS distribution.** Signed distance of each bulk-consensus and pseudobulk G4 peak to the nearest transcript TSS (mm10 GENCODE), shown as binned bars per cell type. Same peak inputs as F01; demonstrates that bulk and single-cell G4 peaks share the same spatial relationship to genes.

**F04 — Stratified recovery of bulk peaks.** For each cell type, bulk-consensus peaks are binned into signal deciles (signalValue) and the percentage recovered in the scCUT&Tag pseudobulk is plotted per decile, plus Jaccard/overlap statistics. Inputs as F01; tools `bedscout::jaccard_index` and `total_bp_overlap`.

**F05 — G4 promoter signal vs gene expression.** Mean G4 RPGC signal at gene promoters (TSS ± 2 kb) of the two scCUT&Tag clusters plotted against matched-cell-type RNA-seq expression deciles (mESC/MEF TPM reference). Cluster pseudobulk RPGC tracks from GSE291468; expression from the TPM reference.

**F06 — G4–expression relationship (overlay and metagene).** Left: both clusters overlaid as mean promoter-signal vs expression-decile lines with SEM ribbons. Right: metagene profile of G4 signal ± 3 kb around the TSS stratified by expression quintile. Inputs as F05.

**F09 / F10 — GO Biological Process of cell-specific G4 peaks.** Cluster-discriminating (unique) peaks are extracted from the bulk-vs-cluster overlap analysis (bulk G4 CUT&Tag GSM8836082/84 vs GSM8836088 cluster peaks), annotated to nearest genes, and tested for GO-BP enrichment (clusterProfiler `enrichGO`, BH-corrected, simplified). Dot plots show the top terms for MEF-specific (F09) and mESC-specific (F10) peaks.

**F11 — Fig 3C threshold validation.** For the four astrocyte G4-marked genes (*Tnik*, *Pitpnc1*, *Pbx1*, *Nwd1*), G4 peaks are re-called from the predicted-astrocyte vs predicted-non-astrocyte RPGC tracks (GSM8836086) by signal thresholding (≥ 30/40/50 RPGC, min width 100 bp) and the number of peaks overlapping a PQS site is plotted per gene and threshold. Supporting table records the AST/non-AST peak-count ratio across thresholds.

**F12 — PQS overlap of scG4 peaks.** Percentage of scCUT&Tag cluster peaks (GSM8836088) overlapping at least one PQS site (pqsfinder predictions on mm10, minimum score 20), pairwise Euler diagrams (cluster 0 / cluster 1 vs PQS), and significance assessment (binomial test against genomic expectation and a regioneR permutation test, 1,000 permutations on chr1).

**F13 — Published bulk 3T3 G4 signal at scG4 peaks.** Heatmap of the published bulk 3T3 G4 CUT&Tag replicates (GSE217860) and PQS score at scCUT&Tag MEF cluster-0 peaks (± 2 kb, rows sorted by scG4 signal). Tool: `EnrichedHeatmap::normalizeToMatrix` + `ComplexHeatmap`; supporting Spearman correlation matrix among scG4, bulk replicates and PQS.

**F14 — Unsorted brain G4 vs per-cluster brain ATAC enrichment.** Genome-wide fold enrichment over random (regioneR `overlapPermTest`, 1,000 permutations; 95% bootstrap CI) of merged unsorted brain G4 peaks (GSM8836087) and per-cell-type brain ATAC accessible peaks (derived from the GSE198467 snATAC object, ≥ 5% cells per cluster) against the PQS and eG4 references. Stars indicate bootstrap tests that G4 enrichment exceeds each ATAC cluster.

**F15 — scCUT&Tag pseudobulk vs matched bulk ATAC enrichment.** Same permTest design for the cell-line data: cluster 0 (MEF) vs 3T3 ATAC (GSE211123, GSM6451000) and cluster 1 (mESC) vs mESC ATAC (GSE149080, GSM4661960), for the PQS and eG4 references, with a bracket testing scG4 > matched ATAC.

**F16 / F17 — g:Profiler enrichment of cell-specific G4 peak genes.** g:Profiler (`gost`, multi-source: GO, KEGG, Reactome, WikiPathways, TF, etc.; whole-genome background, g:SCS correction) on the genes nearest to the mESC-unique (F16) and MEF-unique (F17) G4 peaks. Dot plots highlight self-renewal/pluripotency/signalling terms for mESC and fibroblast/cytoskeleton/ECM terms for MEF.

**F18 — scG4 vs ATAC peak-overlap Venn.** Euler diagrams of scG4 cluster peaks (GSM8836088) vs matched bulk ATAC peaks (GSE149080 mESC; GSE211123 3T3), showing a large fraction of scG4 peaks with no ATAC peak.

**F19 — Accessibility-independent G4 enrichment.** F19a: percentage of peaks overlapping the eG4 and PQS references per class (scG4 & ATAC, scG4-only, ATAC-only, random), showing scG4-marked sites are G4-enriched independent of accessibility. F19b: mean scG4 and ATAC RPGC signal by class, showing strong G4 signal at low accessibility.

**F21 — eG4-overlap benchmark.** Percentage of peaks overlapping the validated eG4 reference for bulk G4 CUT&Tag peaks, all scG4 peaks, scG4-only peaks and random regions per cell type, benchmarking that even bulk G4 reaches only ~32–36% and scG4-only peaks remain 21–40× over random.

**F20 — G4 motif density in scG4 vs ATAC peaks.** Matched-N, matched-width (summit ± 150 bp) scCUT&Tag and ATAC peaks (GSM8836088; GSE149080/GSE211123) are scored for G4-forming potential: pqsfinder maximum G4 score (F20a) and canonical G4 motif counts (regex, both strands; F20b). Demonstrates stronger G4-motif enrichment in scG4 than ATAC (MEF; p ≪ 0.001, Wilcoxon).

**F22 — Astrocyte genes in brain snATAC.** Independent accessibility check of the Fig 3C genes (*Tnik*, *Pitpnc1*, *Pbx1*, *Nwd1*): mean scATAC gene-activity in astrocyte (AST_NT + AST_TE) vs non-astrocyte cells of the GSE198467 brain snATAC object, with fold change and Wilcoxon P per gene.

**F23 — Neural scG4 data quality (QC).** F23a: per-cell QC violins (reads/cell, FRiP, TSS fragments, peaks detected) for the unsorted mouse-brain G4 dataset (GSM8836087 Seurat object) split by cluster. F23b: cross-dataset comparison of median high-quality fragments/cell and FRiP across the mESC-MEF cell line, sorted (GFP+) brain and unsorted brain datasets (Cell Ranger values).
