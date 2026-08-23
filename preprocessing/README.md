# Preprocessing Scripts (OPTIONAL)

This directory contains scripts that regenerate intermediate files (Seurat objects, 
QC plots, etc.). **You do NOT need to run these** if you have downloaded the 
precomputed results from `results/`.

## When to Run Preprocessing

Run these scripts only if:
- You want to modify preprocessing parameters (clustering resolution, QC thresholds)
- Precomputed files in `results/` are missing or corrupted
- You want to reproduce the entire pipeline from raw data
- You're adding new samples or analyses

## Scripts

### Legacy Bulk Processing (`pipeline_*`)

These scripts reference old cluster paths and are kept for reference only:
- `bulk_CnT_workflow.sh` - Bulk CUT&Tag processing (alignment, dedup, bigwig)
- `qc_plots.R` - QC violin plot generation
- `macs3.sh` - Peak calling
- `featureCounts.sh` - Gene counting

**Note**: These are NOT used for Figures 1-4. The figures use preprocessed data
from GEO (GSE291468).

### Annotation Utilities

- `annotation.R` - Gene annotation helper functions

## Output

Preprocessing scripts write to `../results/`. See `../results/README.md` for
documentation of what files are precomputed vs regenerable.

## Relationship to Figure Generation

```
Raw Data (GEO) → Preprocessing → Intermediate Files (results/) → Figure Scripts → Outputs (fig*/outputs/)
```

For Figures 1-4:
1. Download data via `../data/download_data.sh`
2. (Optional) Run preprocessing if needed
3. Run figure scripts (`fig1/fig1.R`, etc.)

## Seurat Processing

The main Seurat processing workflows are in the figure directories:
- `../fig1/seurat_workflow.R` - mESC-MEF mixture processing
- `../fig2/seurat_integration.R` - GFP+ brain integration
- `../fig4/map_unsorted.R` - Unsorted brain mapping
- `../fig4/processing_neuron_data.R` - Neuron data processing

These are the authoritative preprocessing scripts for the figures.