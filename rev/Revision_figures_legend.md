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

**F24 — RNA-seq expression of scG4 G4 peaks by differential category.** Consensus G4 peaks from the scG4 mixture clusters 0/1 (resized to ±500 bp) are scored with the cluster RPGC tracks and classified as unchanged/MEF/ESC by mean signal and log2(cluster1/cluster0) fold change (the low-signal class is excluded). Peaks are annotated to genes (RefSeq TSS ± 1 kb) and the violin shows log2(ESC RPKM / MEF RPKM) from the Chronis et al. mESC/MEF RNA-seq table (GSE90894) per category. A paired split violin (log10 RPKM) with paired t-tests (MEF vs ESC RPKM per category) is also shown. Figure colours match the Fig 1 scheme (MEF blue, mESC red).

**F25 — Enrichr GO biological process terms by G4 peak differential category.** Using the F24 per-peak categories, the unique genes annotated (RefSeq TSS ± 1 kb) to peaks in each category (unchanged/MEF/ESC) are tested with enrichR against GO_Biological_Process_2023 (queried per category). Heatmap of top-10 significant terms per category (P < 0.05) as −log10 P, colours and layout analogous to fig1 panel G. ESC peaks are enriched for synapse/axonogenesis terms, MEF peaks for cell-matrix adhesion/non-canonical Wnt signalling. The unchanged category has no significant enrichments (best raw P = 0.06) and is shown in grey.

**F26 — Spearman and Pearson correlation between scG4 and bulk tracks at marker G4 peaks.** Remake of fig1 panel S1 on the combined marker G4 peaks (union of the ESC and MEF differential categories from F24, RNA-seq matched; n = 4,420 peaks, 1,791 MEF + 2,629 ESC). Signal is scored with the cluster0/cluster1 RPGC tracks and the mESC/MEF bulk rep1 tracks; heatmaps show Spearman and Pearson correlation matrices (Spearman: cluster0↔MEF bulk ρ = 0.66, cluster1↔mESC bulk ρ = 0.67; Pearson: cluster0↔MEF r = 0.73, cluster1↔mESC r = 0.72; cross terms negative).

**F27 — G4 peak overlap and signal at differentially expressed genes.** Reverse of F24: genes with |log2(ESC RPKM / MEF RPKM)| > 4 from the Chronis et al. RNA-seq table (GSE90894) are classified as ESC-up (1,818 genes) or MEF-up (1,327 genes). Each gene set is intersected with the G4 peaks (cluster0+cluster1 union, RefSeq TSS ±1 kb overlap). Left: stacked bar graph (x-axis: ESC-up / MEF-up; fill: G4 overlapping / non-overlapping; ESC-up 372/1818 = 20.5%, MEF-up 394/1327 = 29.7%). Middle: violin of G4 signal lfc (log2 cluster1 / cluster0) per DE category. Right: paired dodge violins of MEF (cl0) and ESC (cl1) G4 RPGC signal at overlapping peaks (y-axis 0–100), with paired t-tests (ESC-up: mean MEF 15.6, ESC 22.2, P = 5.2×10⁻¹⁸; MEF-up: mean MEF 25.5, ESC 20.6, P = 4.1×10⁻⁹).
