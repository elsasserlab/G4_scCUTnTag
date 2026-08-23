# Figure 3: Astrocyte-specific G4 Analysis

## Overview
This figure shows astrocyte-specific G4 patterns and Cicero coaccessibility analysis.

## Prerequisites
Before running `fig3.R` and `cicero.R`, you must:

### 1. Generate integration objects (from fig2)
```bash
cd ../fig2
Rscript seurat_integration.R
```

This creates integration outputs that fig3 depends on.

**Note:** scBridge CSV files are shipped with this directory (GPU required to regenerate).

## Inputs
- From `../fig2/results/integration/outputs/`:
  - `scRNA_Seq_Seurat_object.Rds`
  - `G4_scRNA_integration.Rds`
- From `results/scBridge/`:
  - `umap_coembedded.csv` (shipped)
  - `scbridge_predictions.csv` (shipped)
  - `scbridge_reliability.csv` (shipped)
- From `../data/`:
  - `GSE291468/GSM8836086_GFPpos_Seurat_object.Rds`
  - `pqsfinder/PQS_scores.mm10.bw`
  - `cCRE/Li_et_al-mousebrain.union.cCRE_with_K27ac.Rds`

## Production Scripts
```bash
cd fig3
Rscript fig3.R      # Main figure
Rscript cicero.R    # Cicero coaccessibility (optional)
```

## Outputs
- `outputs/panel_*.pdf` - Figure panels

## Regeneration Time
- `seurat_integration.R` (fig2): ~30-60 minutes
- `fig3.R`: ~20-30 minutes
- `cicero.R`: ~1-2 hours (coaccessibility computation)

## Note on scBridge Files

The scBridge CSV files in `results/scBridge/` are **shipped with the repository** because:
- scBridge requires GPU for training
- Files are small (~840KB total)
- See `results/scBridge/README.md` for regeneration instructions if needed