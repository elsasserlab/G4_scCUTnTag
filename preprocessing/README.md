# Preprocessing Scripts (OPTIONAL)

This directory contains scripts that regenerate intermediate files (Seurat objects,
QC plots, etc.). **You do NOT need to run these** if you have downloaded the
precomputed results from the figure result directories.

## When to Run Preprocessing

Run these scripts only if:
- You want to modify preprocessing parameters (clustering resolution, QC thresholds)
- Precomputed files are missing or corrupted
- You want to reproduce the entire pipeline from raw data
- You're adding new samples or analyses

## Annotation Utilities

- `annotation.R` - Gene annotation helper functions

## Output

Preprocessing scripts write to `results/` at the repository root. Figure scripts
document which intermediate files they read.

## Relationship to Figure Generation

```
Raw Data (GEO) -> Preprocessing -> Intermediate Files -> Figure Scripts -> Outputs (fig*/outputs/)
```

For Figures 1-4:
1. Download data with `bash data/download_data.sh`
2. (Optional) Run preprocessing if needed
3. Run figure scripts (`Rscript fig1/fig1.R`, etc.)

## Seurat Processing

The main Seurat processing workflows are in the figure directories:
- `fig1/seurat_workflow.R` - mESC-MEF mixture processing
- `fig2/seurat_integration.R` - GFP+ brain integration
- `fig4/map_unsorted.R` - Unsorted brain mapping
- `data/Zeisel2018/processing_neuron_data.R` - Neuron data processing

These are the authoritative preprocessing scripts for the figures. Bulk
alignment and peak-calling scripts from the original cluster environment are
not shipped because they contain machine-specific paths; the corresponding
processed files are available from GEO or shipped with the repository.
