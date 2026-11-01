# scBridge Output Files

These files are shipped precomputed. scBridge is **not** a dependency of this
repository and does not need to be installed or run to reproduce any figure.

## Contents

- `umap_coembedded.csv` (516KB) - UMAP coordinates for co-embedded cells
- `scbridge_predictions.csv` (83KB) - Cell type predictions for G4 cells
- `scbridge_reliability.csv` (108KB) - Reliability scores for predictions
- `predictions_reliability.csv` (133KB) - Combined predictions and reliability

**Total size:** ~840KB

## Provenance

The files were generated once with scBridge
(Li et al., 2023, *Nature Communications*;
<https://github.com/XLearning-SCU/scBridge>) from:

1. G4 scCUT&Tag gene activity scores:
   `data/GSE291468/GSM8836086_GFPpos_Seurat_object.Rds`
2. Bartosovic scRNA-seq reference:
   `data/GSE163484/scRNA_Seq-mouse_brain.Rds`

`fig3/scbridge_input.R` documents how the two input tables were prepared. The
CSV files in this directory were copied from the scBridge output directory after
that single run.

### scBridge parameters used

- Integration method: multi-omics CCA
- UMAP dimensions: 2
- Prediction method: label transfer from scRNA-seq to G4

## Figure usage

`fig3/fig3.R` reads `umap_coembedded.csv`, `scbridge_predictions.csv` and
`scbridge_reliability.csv`; `fig3/cicero.R` reads `scbridge_predictions.csv`.
No scBridge execution is required.
