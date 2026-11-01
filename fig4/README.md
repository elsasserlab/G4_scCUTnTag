# Figure 4: Unsorted Brain Neuronal Analysis

## Overview
This figure shows analysis of unsorted G4 CUT&Tag data projected onto neuronal reference.

## Prerequisites
Panels F and G use two precomputed scBridge PDFs shipped in
`fig4/results/scBridge/`:

- `Seurat_UMAPs-unsorted_cl1_int.pdf` (Panel F)
- `tree_plot-unsorted_cl1_predictions.pdf` (Panel G)

scBridge is not a dependency and does not need to be installed or run; the
PDFs are shipped so that Panels F and G are reproduced directly.

## Inputs (from data/)
- `data/GSE291468/GSM8836087_unsorted_Seurat_object.Rds` - Unsorted G4 data
- `data/GSE291468/GSM8836086_GFPpos_Seurat_object.Rds` - GFP+ sorted G4 data
- `data/Zeisel2018/scRNA_Seq-Zeisel_et_al-neuron.Rds` - Zeisel neuron reference
- `data/GSE163484/scRNA_Seq-mouse_brain.Rds` - Bartosovic reference
- `data/pqsfinder/PQS_scores.mm10.bw` - PQS predictions

## Production Scripts
```bash
Rscript fig4/fig4.R              # Main figure
Rscript fig4/map_unsorted.R      # Unsorted mapping (optional)
Rscript data/Zeisel2018/processing_neuron_data.R  # Neuron processing (optional)
```

## Outputs
- `outputs/panel_*.pdf` - Figure panels

## Regeneration Time
- `fig4.R`: ~20-30 minutes
- `map_unsorted.R`: ~30-45 minutes (integration)
