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
- `cicero_GFPsorted*.Rds` - Cicero coaccessibility results

**Why shipped:** Cicero computation takes 1-2 hours

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
Rscript cicero.R
# This will regenerate all cicero outputs
```

## Note

These files are kept in the repository because:
1. They are required by production scripts (`fig3.R`, `cicero.R`)
2. Regeneration requires significant compute resources (GPU, 1-2 hours)
3. Total size is small (~14MB)