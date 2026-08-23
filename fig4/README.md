# Figure 4: Unsorted Brain Neuronal Analysis

## Overview
This figure shows analysis of unsorted G4 CUT&Tag data projected onto neuronal reference.

## Prerequisites
None - all inputs are from data/ directory.

## Inputs (from data/)
- `../data/GSE291468/GSM8836087_unsorted_Seurat_object.Rds` - Unsorted G4 data
- `../data/GSE291468/GSM8836086_GFPpos_Seurat_object.Rds` - GFP+ sorted G4 data
- `../data/scRNA-Seq/scRNA_Seq-Zeisel_et_al-neuron.Rds` - Zeisel neuron reference
- `../data/scRNA-Seq/scRNA_Seq-mouse_brain.Rds` - Bartosovic reference
- `../data/pqsfinder/PQS_scores.mm10.bw` - PQS predictions

## Production Scripts
```bash
cd fig4
Rscript fig4.R              # Main figure
Rscript map_unsorted.R      # Unsorted mapping (optional)
Rscript processing_neuron_data.R  # Neuron processing (optional)
```

## Outputs
- `outputs/panel_*.pdf` - Figure panels

## Regeneration Time
- `fig4.R`: ~20-30 minutes
- `map_unsorted.R`: ~30-45 minutes (integration)