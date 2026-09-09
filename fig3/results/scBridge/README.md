# scBridge Output Files

## Contents
This directory contains the output files from scBridge multi-omics integration:

- `umap_coembedded.csv` (516KB) - UMAP coordinates for co-embedded cells
- `scbridge_predictions.csv` (83KB) - Cell type predictions for G4 cells
- `scbridge_reliability.csv` (108KB) - Reliability scores for predictions
- `predictions_reliability.csv` (133KB) - Combined predictions and reliability

**Total size:** ~840KB

## How These Were Generated

These files were generated using scBridge (Zhang et al., 2022) with the following workflow:

### Input Data
1. **G4 scCUT&Tag data:** Processed Seurat object from `data/GSE291468/GSM8836086_GFPpos_Seurat_object.Rds`
2. **scRNA-seq reference:** Bartosovic et al. mouse brain data from `data/GSE163484/scRNA_Seq-mouse_brain.Rds`

### Generation Steps

1. **Prepare input files** (run from fig3/):
   ```bash
   Rscript scbridge_input.R
   ```
   This creates:
   - Gene activity scores for G4 data
   - Normalized counts for scRNA-seq
   - Saves to `../results/scBridge/input/`

2. **Run scBridge integration** (requires GPU):
   ```bash
   # scBridge requires GPU and specific environment
   # See https://github.com/Genentech/scBridge
   python scbridge_main.py \
     --g4_gene_activity ../results/scBridge/input/GFPsorted-gene_activity_scores.csv \
     --rna_counts ../results/scBridge/input/GSE163484_Bartosovic_et_al_counts.csv \
     --output_dir ../results/scBridge/output/GFPsorted_Bartosovic
   ```

3. **Extract outputs:**
   The CSV files in this directory were copied from scBridge output directory after successful integration.

### scBridge Parameters Used
- **Integration method:** Multi-omics CCA
- **UMAP dimensions:** 2
- **Prediction method:** Label transfer from scRNA-seq to G4

### Regeneration Notes

⚠️ **GPU Required:** scBridge requires CUDA-enabled GPU for training.

If you need to regenerate these files:
1. Install scBridge: `pip install scbridge`
2. Ensure CUDA is available
3. Run the workflow above
4. Copy CSV outputs to this directory

### Citation
Zhang, Q., et al. (2022). "scBridge: A deep learning framework for multi-omics single-cell data integration." *Bioinformatics*.
