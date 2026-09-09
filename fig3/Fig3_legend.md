# Figure 3 — Legend and Derivation Notes

**Figure 3: Reliability-based annotation enhances cell-type assignment of G4 scCUT&Tag data and identifies astrocyte-specific G4 signals.**

## Data sources

- **Sorted GFP+ brain G4 scCUT&Tag (GSE291468, GSM8836086)** — `GSM8836086_GFPpos_Seurat_object.Rds`, gene-activity (GA) assay, peak assay, `GSM8836086_GFP_sorted_mouse_brain/CellRanger/fragments.tsv.gz`
- **Brain scRNA-seq reference (GSE163484, Bartosovic et al.)** — processed Seurat object with cell_type annotations
- **scBridge co-embedding** — precomputed outputs: `umap_coembedded.csv`, `scbridge_predictions.csv` (barcodes_Astrocytes.tsv, barcodes_non_AST.tsv), `scbridge_reliability.csv`
- **Cicero co-accessibility** — shipped AST and non-AST connection tables in `results/cicero/`
- **Mouse-brain cCREs with H3K27ac** — `Li_et_al-mousebrain_cCRE_with_K27ac.bed`
- **ENCODE4 Registry cCREs** — `data/genome/cCRE.mm10.bed` (9-column BED, itemRgb encodes the regulatory class)
- **Predicted AST/non-AST G4 signal** — unsmoothed RPGC bigWigs for GSM8836086
- **Mouse gene annotation** — EnsDb.Mmusculus.v79
- **PQS scores** — `PQS_scores.mm10.bw`

Produced by `fig3/fig3.R`, with Cicero analysis in `fig3/cicero.R`.

## Panels

**Panel A — scBridge joint-embedding UMAP.** Four views share a co-embedding of 5,488 scRNA-seq and 2,767 G4 cells (8,255 total). G4 predictions comprised AST 183, COP-NFOL 199, MOL 126, OEC 746, OPC 599, pericytes 50, VEC 99, VLMC 101 and 664 `Novel (Most Unreliable)` cells. Median reliability was 0.996 among reliably assigned cells and 1.28 x 10^-40 among unreliable cells. Output: `panel_A_coembedded_UMAPs.pdf`.

**Panel B — Astrocyte (AST) vs non-AST differential G4.** **B1:** UMAP contains 183 AST and 2,584 other G4 cells; the latter includes the 664 unreliable cells. **B2/B3:** after matching scBridge barcodes to the Seurat object, differential tests use 159 AST and 1,679 reliably assigned non-AST cells. B2 tests 17,120 GA genes by logistic regression with peak-region fragments as a latent variable and identifies 18 AST-up and 8 AST-down genes. B3 tests peaks by Wilcoxon rank-sum, retaining 8,833 promoter-proximal peaks and identifying 23 AST-up and no AST-down promoter peaks. **B4:** Wilcoxon testing of reference AST versus non-AST expression returns 14,344 genes, of which 1,763 are AST-up and 5,068 AST-down. All calls require |avg_log2FC| > 0.5 and Seurat Bonferroni-adjusted `p_val_adj < 0.05`; `FindMarkers()` itself uses `logfc.threshold=0`. Outputs: `panel_B_AST_vs_nonAST.pdf` and four differential-result CSVs.

**Panel C — Coverage at differential-G4-promoter-peak genes.** Coverage is shown for *Scarna17*, *Tmem74*, *Nwd1*, *Gli2* and *Plcl1*. Group tracks represent 183 predicted AST and the authors' broader 2,521-cell non-AST list, which includes unreliable cells and therefore differs from Panel B. Independently of the significance tests in Panel B, peak detection-rate classification across all 47,126 peaks identified 994 AST-specific and 76 non-AST-specific peaks, defined as detection in at least 2% of the group and at least four-fold higher detection frequency than the other group. This classification has no P-value. Outputs include `panel_C_coverage.pdf`, selected-gene and per-peak tables, and group-specific BED files.

**Panel D — Astrocyte-specific gene expression feature plots.** Normalized expression is shown for *Tmem74*, *Nwd1*, *Gli2* and *Plcl1* across the 5,488 reference cells on the Panel A coordinates. *Scarna17* is absent from the RNA assay and is skipped. Expression color and alpha are capped at each gene's 99th percentile; alpha ranges from 0.2 to 1. These are descriptive ggplot2 feature plots with no test. Output: `panel_D_featureplots.pdf`.

**Panel E — Cicero co-accessibility browser tracks.** Five views show *Rsg1*, *Tmem74*, *Amot*, *Sox10* and *Akt1s1*. At coaccessibility >0.2, the numbers of direct AST/non-AST links from the focal peak were 3/0, 2/0, 2/1, 4/1 and 2/4, respectively. Arc direction, color and width encode group and coaccessibility; no hypothesis test is performed. Output: `panel_E_cicero_tracks.pdf`.

**Panel S1 — AST-specific G4 peaks vs ENCODE4 cCRE classes.** Of 994 AST-specific peaks, 993 lie on canonical chromosomes and are partitioned as dELS 315 (31.7%), CA-CTCF 248 (25.0%), CA-only 120 (12.1%), no cCRE 107 (10.8%), PLS 92 (9.3%), pELS 73 (7.4%), CA-H3K4me3 28 (2.8%), TF-only 6 (0.6%) and CA-TF 4 (0.4%). Each peak receives the class with greatest bp overlap; no test is performed. Outputs: `panel_S1_AST_cCRE_stacked.pdf`, partition CSV and `panel_S1_AST_specific_peaks_annotated.bed`.
