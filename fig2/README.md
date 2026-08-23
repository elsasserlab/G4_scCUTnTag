# Figure 2: CCA Integration with scRNA-seq

## Overview
This figure shows integration of scG4 CUT&Tag data with scRNA-seq reference using Seurat CCA.

## Prerequisites
Before running `fig2.R`, you must generate the integrated Seurat objects:

```bash
cd fig2
Rscript seurat_integration.R
```

This creates:
- `results/integration/outputs/scRNA_Seq_Seurat_object.Rds` (scRNA-seq object)
- `results/integration/outputs/G4_scRNA_integration.Rds` (integrated G4+scRNA object)
- `results/integration/outputs/g4_cell_label_preds.Rds` (cell type predictions)
- `results/integration/outputs/scRNA-Seq-FindAllMarkers_output.tsv` (marker genes)

**Note:** These files are gitignored (~1GB total) and must be regenerated locally.

## Inputs (from data/)
- `../data/GSE291468/GSM8836086_GFPpos_Seurat_object.Rds` - G4 scCUT&Tag data
- `../data/scRNA-Seq/scRNA_Seq-mouse_brain.Rds` - Bartosovic scRNA-seq reference

## Production Script
After generating intermediates:
```bash
Rscript fig2.R
```

## Outputs
- `outputs/panel_*.pdf` - Figure panels

## Regeneration Time
- `seurat_integration.R`: ~30-60 minutes (CCA integration)
- `fig2.R`: ~10-20 minutes (plotting)