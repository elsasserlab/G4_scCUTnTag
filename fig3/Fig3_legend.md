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

**Panel C — Coverage at differential-G4-promoter-peak genes.** Genome coverage (fragment tracks: predicted AST vs predicted non-AST, GSM8836086) at genes with significant AST-up G4 peaks at their promoter (Panel B), prioritised by cross-assay support (tier 1 = also GA up-hit and scRNA-seq upregulated). Each window spans the hit peak ±5 kb plus ~30 kb into the gene body. Below each track: called peaks (Signac `CoveragePlot`, `peaks = TRUE`) and experimental G4 sites (eG4 catalogue, colored by confidence level).

**Panel D — Astrocyte-specific gene expression feature plots.** UMAP feature plots of normalized scRNA-seq expression for the four highlighted astrocyte genes on the scRNA-seq reference cells (GSE163484), on the same co-embedding coordinates as panel A.

**Panel E — Cicero co-accessibility browser tracks.** Candidate regulatory-element browser views connecting astrocyte-specific G4 peaks with mouse-brain cCREs (cCREs with H3K27ac signal above the 75th percentile), with G4 peaks, PQS sites and regulatory elements indicated; aggregated AST and non-AST G4 coverage are shown alongside. Derived from the Cicero co-accessibility analysis on the G4 peaks.
