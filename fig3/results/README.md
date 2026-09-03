# Figure 3: Results Directory

## Contents

This directory contains intermediate files required by fig3 production scripts:

### scBridge/ (shipped - GPU required)
- `umap_coembedded.csv` - UMAP coordinates
- `scbridge_predictions.csv` - Cell type predictions
- `scbridge_reliability.csv` - Reliability scores
- `predictions_reliability.csv` - Combined data

**Why shipped:** scBridge requires GPU for training (~1 hour computation)

### cicero/ (shipped - computationally intensive)
- `cicero_GFPsorted*.Rds` - combined, AST, and non-AST Cicero connections
- `cicero_GFPsorted_coG4networks*.Rds` - corresponding CCAN assignments

**Why shipped:** Cicero computation takes 1-2 hours. Figure 3E reads the AST
and non-AST connection tables directly, so this computation is not required
for routine figure generation.

## Total Size: ~14MB

## Regeneration

If you need to regenerate these files:

### scBridge (requires GPU)
```bash
Rscript scbridge_input.R
# Then run scBridge with GPU...
```

### cicero
```bash
# From the repository root; requires cicero and monocle3
Rscript fig3/cicero.R
# Regenerates all six shipped Cicero/CCAN files
```

## Note

These files are kept in the repository because:
1. They are required by production scripts (`fig3.R`, `cicero.R`)
2. Regeneration requires significant compute resources (GPU, 1-2 hours)
3. Total size is small (~14MB)
