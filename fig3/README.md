# Figure 3: Astrocyte-specific G4 Analysis

## Overview
This figure shows astrocyte-specific G4 patterns and Cicero coaccessibility analysis.

## Prerequisites
Before running all panels in `fig3.R`, you must:

### 1. Generate integration objects (from fig2)
```bash
cd ../fig2
Rscript seurat_integration.R
```

This creates integration outputs that fig3 depends on.

**Note:** scBridge CSV files are shipped with this directory (GPU required to regenerate).

### 2. Regenerate Smooth BigWig Files (Required for Figure 3C)

For Figure 3C, you need to generate smoothed bigWig files (150bp moving average):

```bash
cd ../../data/GSE291468
bash regenerate_bigwigs.sh
# Takes ~2-4 hours total (2 files)
cd ../..
```

This generates:
- `GSM8836086_Predicted_Astrocytes_RPGC.smooth150.bw`
- `GSM8836086_Predicted_non-Astrocytes_RPGC.smooth150.bw`

**Requirements:**
- Compute cluster with SLURM (or run manually with deepTools)
- `bedtools`, `samtools`, and `deepTools` installed
- ~2 hours per file

**Note:** The unsmoothed corrected bigWigs are already included in the repository.

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
  - `GSE291468/GSM8836086_Predicted_Astrocytes_RPGC.bw`
  - `GSE291468/GSM8836086_Predicted_non-Astrocytes_RPGC.bw`
  - `pqsfinder/PQS_scores.mm10.bw`
  - `cCRE/Li_et_al-mousebrain_cCRE_with_K27ac.bed`
- From `results/cicero/` (shipped):
  - `cicero_GFPsorted-predAST.Rds`
  - `cicero_GFPsorted-pred_nonAST.Rds`

## Production Scripts
```bash
# From the repository root
Rscript fig3/fig3.R
Rscript fig3/fig3.R --only=E
```

Panel E is generated directly by `fig3.R` using `cog4_browser.R`. It produces
one long multiplot containing Rsg1, Tmem74, Amot, Sox10, and Akt1s1. The saved
Cicero connection tables are shipped, so running `cicero.R` first is not
required.

To regenerate the shipped Cicero connection and CCAN intermediates instead:

```bash
Rscript fig3/cicero.R  # Optional; requires cicero + monocle3, ~1-2 hours
```

## Outputs
- `outputs/panel_*.pdf` - Figure panels
- `outputs/panel_E_cicero_tracks.pdf` - Directly generated five-locus browser multiplot

## Regeneration Time
- `seurat_integration.R` (fig2): ~30-60 minutes
- `fig3.R`: ~20-30 minutes
- `cicero.R`: ~1-2 hours (coaccessibility computation)

## Note on scBridge Files

The scBridge CSV files in `results/scBridge/` are **shipped with the repository** because:
- scBridge requires GPU for training
- Files are small (~840KB total)
- See `results/scBridge/README.md` for regeneration instructions if needed
