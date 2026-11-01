# Additional Figures — Legend and Derivation Notes

Additional Figures. Each panel produced by R scripts under `rev/`.

## Shared data sources

- **mESC-MEF single-cell G4 scCUT&Tag (GSE291468, GSM8836088)** — cluster-specific MACS2 narrowPeak sets (`GSM8836088_cluster_0_peaks.narrowPeak`, `GSM8836088_cluster_1_peaks.narrowPeak`), pseudobulk RPGC bigWigs (`cluster0_RPGC.bw`, `cluster1_RPGC.bw`), Seurat object
- **Bulk G4 CUT&Tag (GSE291468)** — bulk mESC replicates (GSM8836082, GSM8836083) and bulk MEF/3T3 replicates (GSM8836084, GSM8836085) as broadPeak calls and RPGC bigWigs
- **Sorted GFP+ brain G4 (GSE291468, GSM8836086)** — predicted-astrocyte and predicted-non-astrocyte RPGC bigWigs (`GSM8836086_Predicted_Astrocytes_RPGC.smooth150.bw`, `GSM8836086_Predicted_non-Astrocytes_RPGC.smooth150.bw`), label-transfer outputs
- **Unsorted mouse-brain G4 (GSE291468, GSM8836087)** — Seurat object, cluster-specific narrowPeak sets and RPGC bigWigs
- **mESC ATAC-seq (GSE149080, GSM4661960)** and **3T3 ATAC-seq (GSE211123, GSM6451000)** — matched bulk accessibility narrowPeak sets and RPGC bigWigs
- **Published 3T3 G4 CUT&Tag (GSE217860, GSM6729104–06)** — three replicates (aligned bigWigs) used as independent bulk benchmark
- **Brain snATAC reference (GSE198467)** — `GSE198467_ATAC_Seurat_object_clustered_renamed.Rds`, clustered Seurat object with 15 named cell types
- **Gene expression reference** — `gene_level_TPM.tsv` (mESC/MEF TPM table), `gene_level_RPKM.xlsx` (Chronis et al. GSE90894), mm10 GENCODE annotation GTF (`mm10_.annotation.gtf.gz`)
- **G4 references** — pqsfinder PQS predictions (`PQS_scores.mm10.bed`, `PQS_scores.mm10.bw`), experimentally validated eG4 catalog (`eG4_catalog.bed`)
- **RefSeq annotation** — `refGene.tss.1kb.mm10.bed` (TSS ±1kb)

## Panels

**F01 — Genomic feature distribution (bulk vs scCUT&Tag).** Peak totals were 23,600 mESC bulk-consensus, 34,133 mESC pseudobulk, 19,866 MEF bulk-consensus and 31,009 MEF pseudobulk. Promoter assignments were 13,113 (55.6%), 15,004 (44.0%), 11,266 (56.7%) and 13,583 (43.8%), respectively. Peaks were assigned by midpoint using the priority promoter → 5′UTR → 3′UTR → CDS exon → intron → intergenic; this panel is descriptive and no test was performed. Outputs include `feature_distribution.csv`, two PDFs and per-peak feature TSVs.

**F12 — PQS overlap of scG4 peaks.** At PQS score ≥20, overlaps were 26,979/31,009 MEF peaks (87.0%), 30,244/34,133 mESC (88.6%), 16,145/18,652, 11,608/13,288 and 8,182/9,735 GFP+ clusters (86.6%, 87.4%, 84.0%), and 22,827/26,356 and 13,726/15,620 unsorted clusters (86.6%, 87.9%). Width-adjusted brain scATAC overlap was 34,798/46,204 (75.3%). Chr1 permutation tests used 1,000 permutations and gave observed/expected overlap 17,311/6,237.6 (2.78-fold), 9,426/3,044.6 (3.10-fold) and 9,602/3,240.1 (2.96-fold) for cell-line, GFP+ and unsorted sets. Permutation P-values are raw with resolution approximately P≤0.001; no multiplicity adjustment was applied. Outputs: bar/Euler PDFs and count/permutation CSVs.

**F18 — scG4 vs ATAC peak-overlap Venn.** Summit-centered 301-bp intervals gave mESC counts of 16,428 shared, 17,705 scG4-only and 46,311 ATAC-only, and MEF counts of 20,794, 10,215 and 34,263. Thus scG4-only fractions were 51.9% for mESC and 32.9% for MEF. The panel is descriptive; an optional post-hoc hypergeometric and regioneR permutation test is provided in `rev/test_F18_significance.R` but is not plotted. Outputs: two Euler PDFs.

**F24 — RNA-seq expression of scG4 G4 peaks by differential category.** Signal classification yielded 17 low, 33,249 unchanged, 4,183 MEF and 6,530 ESC peaks; RNA matching retained 10,300 unchanged, 269 MEF and 723 ESC peaks. Median/mean log2(ESC/MEF RPKM) was 0.011/−0.093, −0.677/−0.724 and 0.656/1.545, respectively. Mean MEF/ESC raw RPKM was 34.90/39.32 (unchanged), 19.64/23.57 (MEF) and 10.61/18.72 (ESC), with raw paired t-test P=3.06 x 10^-4, 0.589 and 4.44 x 10^-4. No multiplicity adjustment was applied. Outputs: two PDFs and peak/RNA/statistics CSVs.

**F25 — Enrichr GO biological process terms by G4 peak differential category.** Gene-set sizes were 10,520 unchanged, 316 MEF and 847 ESC genes. Numbers of nominal-P<0.05 terms were 1,173, 136 and 200; adjusted-P<0.05 counts were 446, 0 and 0. Best nominal/adjusted P-values were 3.11 x 10^-21/1.66 x 10^-17, 3.53 x 10^-4/0.422 and 3.31 x 10^-5/0.0894, respectively. The heatmap uses the top ten nominally significant terms and −log10 nominal P, not adjusted P. Outputs: PDF and enrichment/gene-set TSVs.

**F26 — Spearman and Pearson correlation between scG4 and bulk tracks at marker G4 peaks.** The RNA-matched analysis contains 269 MEF + 723 ESC regions (n=992), with matched scG4-bulk Spearman/Pearson values 0.660/0.688 and 0.597/0.649. The all-differential analysis contains 4,183 + 6,530 regions (n=10,713), with values 0.644/0.737 and 0.672/0.718. Called-peak membership was 3,396 cluster0-only, 6,370 cluster1-only and 947 shared. Correlations use untransformed `bw_loci()` values and no P-values. Outputs: correlation PDFs/CSVs and category-count CSV.

**F27 — G4 peak overlap and signal at differentially expressed genes.** At |log2(ESC/MEF RPKM)|>4, 372/1,818 ESC-up genes (20.5%) and 394/1,327 MEF-up genes (29.7%) overlapped a promoter G4 peak. One peak per gene was selected by maximum absolute G4 log2FC. The log2FC violin uses a raw Wilcoxon test. Paired mean MEF/ESC RPGC was 14.705/21.454 for ESC-up genes (raw paired t-test P=1.86 x 10^-15) and 25.732/20.427 for MEF-up genes (P=2.55 x 10^-8). Outputs: two PDFs and two CSVs.

**F28 — scG4 vs scATAC peak overlap (Venn).** scG4 peaks are the union of the four GFP+ cluster MACS2 narrowPeaks (26,982 regions after merging, median 285 bp), and scATAC peaks are the MACS2 re-called consensus narrowPeak (92,926 regions, median 466 bp). Both restricted to canonical chromosomes. Merging all peaks produces 70,364 unique regions, of which 26,656 carry scG4 signal, 64,251 carry scATAC signal, and 20,543 carry both (77.1% of scG4; 32.0% of scATAC), leaving 6,113 scG4-only (22.9%) and 43,708 scATAC-only (68.0%). No test was performed. Outputs: Euler PDF and statistics CSV.

**F30 — scG4 vs scATAC scatter plot (Union vs Intersection).** Uses the same consistent narrowPeak sets as F28: scG4 union of GFP+ cluster peaks (26,982) and scATAC re-called consensus (92,926); overlapping peak pairs reduce to 20,542 intersection intervals and all peaks to 70,364 union intervals. Continuous coverage (RPGC-normalized, 10 bp bins) for both assays is scored over the whole interval with `wigglescout::bw_loci` (default 0 when no read covers a bin), and plotted as log1p. Pearson/Spearman correlation was 0.166/0.135 for intersection intervals and 0.274/0.295 for the union; 6,113 union regions lack scG4 coverage and 168 lack scATAC coverage, giving no artificial zeros at single-modality regions. No P-values were calculated. Outputs: scatter PDF and intersection/union CSVs.

## Supplementary QC (not a numbered panel)

`rev/rev_QC_GSM5501245.R` and `rev/rev_QC_GSM5501245_comparison.R` generate
the QC panels for the third-party GSM5501245 (SINAG1, U2OS/MCF7) dataset:
the comparison panel `panel_QC_GSM5501245.pdf` (written by
`rev_QC_GSM5501245_comparison.R`; summarises peak-width, chromosome
distribution and published per-cell QC metrics versus the ESC-MEF data)
and the single-dataset panel `panel_QC_GSM5501245_single.pdf` (written by
`rev_QC_GSM5501245.R`).
These QC plots support the revision but are not assigned an F-number.

**Note on permutation z-scores:** the regioneR permutation tests (`permTest`) are run in parallel (`mc.cores > 1`). The random permutation stream is therefore not fully deterministic across runs/machines, and reported z-scores can vary slightly (observed/expected ratios and resolution-limited P-values, ~P<=0.001, are stable). The z-scores in the CSVs are approximate.