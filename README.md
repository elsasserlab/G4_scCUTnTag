# scG4 Mouse Brain Analysis

**Genome Biology** - Lyu et al. (2025)

This repository contains all code and data to reproduce Figures 1-4 from our study on G-quadruplex (scG4) patterns in mouse brain.

---

## Quick Start

### 1. Clone and Download Data

```bash
git clone <repo-url>
cd github/data
bash download_data.sh
cd ..
```

This downloads ~2GB of data from GEO (accession GSE291468 and references).

### 2. Generate Figures

```bash
# Figure 1: mESC-MEF G4 CUT&Tag analysis
Rscript fig1/fig1.R

# Figure 2: Integration with scRNA-seq reference
cd fig2
Rscript seurat_integration.R  # Generate intermediates first
Rscript fig2.R
cd ..

# Figure 3: Astrocyte-specific G4 patterns
cd fig3
Rscript fig3.R
Rscript cicero.R  # Optional: coaccessibility analysis
cd ..

# Figure 4: Unsorted brain neuronal analysis
cd fig4
Rscript fig4.R
cd ..
```

---

## Repository Structure

```
github/
├── data/                    # Data files (~651MB shipped + ~2GB downloaded)
│   ├── GSE291468/          # Primary scG4 CUT&Tag data (our study)
│   ├── GSE163484/          # Bartosovic scRNA-seq reference
│   ├── GSE198467/          # Brain snATAC-seq reference
│   ├── genome/             # Reference genome (mm10)
│   ├── cCRE/               # Mouse brain cCRE with H3K27ac
│   ├── pqsfinder/          # PQS predictions
│   └── download_data.sh    # Download script for GEO data
├── fig1/                   # Figure 1: mESC-MEF analysis
├── fig2/                   # Figure 2: scRNA-seq integration
├── fig3/                   # Figure 3: Astrocyte-specific G4
├── fig4/                   # Figure 4: Unsorted brain analysis
└── rev/                    # Revision supplementary figures
```

---

## Figure Reproduction

Each figure lives in its own directory (`fig1/`–`fig4/`) containing the plotting script and, where needed, preprocessing scripts.

### Directory Layout

```
fig{N}/
  fig{N}.R                 # Panel plotting script (outputs to fig{N}/outputs/)
  preprocessing.R          # Optional: required preprocessing
  outputs/                 # Generated figure panels (PDF)
    panel_A_*.pdf          # Named panels matching manuscript
    panel_extra_*.pdf      # Supplementary/alternative panels
  README.md                # Regeneration instructions
```

### Main Figures

| Figure | Script | Preprocessing | Description |
| --- | --- | --- | --- |
| Fig. 1 | `fig1/fig1.R` | `fig1/seurat_workflow.R` | mESC-MEF UMAP, QC, Venn, PCA, heatmap, correlation, tracks |
| Fig. 2 | `fig2/fig2.R` | `fig2/seurat_integration.R` | Brain integration UMAP, prediction scores, overlap, markers |
| Fig. 3 | `fig3/fig3.R` | `fig3/cicero.R`, `fig3/scbridge_input.R` | scBridge UMAP, Cicero tracks |
| Fig. 4 | `fig4/fig4.R` | `fig4/map_unsorted.R`, `fig4/processing_neuron_data.R` | Unsorted brain mapping |

### Revision Figures

All 17 revision scripts are in `rev/` and output to `rev/outputs/`. See `rev/paths.R` for path setup and `rev/Revision_figures_legend.md` for descriptions.

---

## Data Sources

### Primary Data (GSE291468)

Our scG4 CUT&Tag data from three experiments:
- **GSM8836086**: GFP+ sorted postnatal mouse brain (Figs 2-3)
- **GSM8836087**: Unsorted postnatal mouse brain (Fig 4)
- **GSM8836088**: mESC-MEF mixture (Fig 1)

### Reference Data

| Accession | Description | Used In |
|-----------|-------------|---------|
| GSE163484 | Bartosovic scRNA-seq | Fig 2-3 |
| GSE198467 | Brain snATAC-seq | Fig 4 |
| GSE149080 | mESC ATAC-seq | Fig 1 |
| GSE211123 | 3T3 ATAC-seq | Fig 1 |
| GSE217860 | 3T3 G4 CUT&Tag | Fig 1 |
| GSE90894 | Chronis RNA-seq | Rev |

### Genome Resources

- **genome/**: mm10 chromosome sizes, GTF annotation, eG4 predictions
- **cCRE/**: Li et al. mouse brain cCRE with H3K27ac signal
- **pqsfinder/**: Genome-wide PQS predictions (PQS_scores.mm10.bed)

---

## GEO File-to-Figure Map

| Figure | GEO Sample(s) | Key Files | Role |
|--------|---------------|-----------|------|
| **Fig. 1** | GSM8836088 (mESC-MEF)<br>GSM8836082/84 (bulk controls) | `GSM8836088_mESCMEF_Seurat_object.Rds`<br>`GSM8836088_cluster_{0,1}_peaks.narrowPeak`<br>`GSM8836082/84_bulkG4CnT_*.bw` | mESC-MEF UMAP/QC, cluster-vs-bulk overlaps, marker PCA, RPGC heatmaps, correlation, browser tracks |
| **Fig. 2** | GSM8836086 (sorted GFP+)<br>GSE163484 (scRNA reference) | `GSM8836086_GFPpos_Seurat_object.Rds`<br>`scRNA_Seq-mouse_brain.Rds`<br>`results/integration/outputs/` | CCA co-embedding, transferred cell labels, prediction scores, overlap scores, marker violins, RNA/G4 feature plots |
| **Fig. 3** | GSM8836086 (sorted GFP+)<br>GSE163484 (scRNA reference) | `results/scBridge/umap_coembedded.csv`<br>`scbridge_predictions.csv`<br>`GSM8836086_Predicted_Astrocytes_RPGC.bw` | scBridge embedding, AST differential G4, AST coverage, RNA feature plots, reliability analyses |
| **Fig. 3E** | Derived from sorted GFP+ | `results/cicero/cicero_browser_example-AST_spec.pdf`<br>`data/cCRE/` | Cicero co-accessibility and candidate regulatory-element browser tracks |
| **Fig. 4** | GSM8836087 (unsorted)<br>GSM8836086 (sorted reference) | `GSM8836087_unsorted_Seurat_object.Rds`<br>`GSM8836087_unsorted_cluster_{0,1}_peaks.narrowPeak` | Query-to-reference mapping, mapped labels, prediction scores, peak-overlap Venn |
| **Fig. 4C** | GSM8836087 (unsorted clusters) | `GSM8836087_cluster{0,1}_RPGC.bw`<br>`GSM8836086_GFP_sorted_mousebrain.rpgc.bw`<br>`PQS_scores.mm10.bw` | RPGC signal profiles/heatmaps for cluster-specific peaks, PQS score profiles |
| **Fig. 4D** | GSM8836087 (unsorted) | `GSM8836087_unsorted_cluster_{0,1}_peaks.narrowPeak` | Nearest-gene annotation and enrichR GO Biological Process analysis |
| **Fig. 4F-G** | GSM8836087 + Zeisel neuron reference | `results/scBridge/output/unsorted_cl1_Zeisel/`<br>`Zeisel_et_al-neuron.Rds` | Precomputed scBridge neuron integration UMAP and neuron-type treemap |

---

## Software Requirements

### R Packages

```r
# Core packages
install.packages(c("Seurat", "Signac", "tidyverse", "data.table", "glue"))

# Bioconductor packages
if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(c("rtracklayer", "GenomicRanges", "EnsDb.Mmusculus.v79"))

# Additional packages
install.packages(c("cowplot", "ggpubr", "ComplexHeatmap", "circlize", "argparse"))
```

### System Requirements

- **R** >= 4.0
- **Python** >= 3.7 (for scBridge, optional)
- **GPU** (CUDA) - only if regenerating scBridge integration

---

## Important Notes

### Figure 2 Prerequisites

Before running `fig2.R`, you must generate integration objects:

```bash
cd fig2
Rscript seurat_integration.R  # ~30-60 minutes
Rscript fig2.R                # ~10-20 minutes
```

See `fig2/README.md` for details.

### Figure 3: scBridge Files

The scBridge integration outputs in `fig3/results/scBridge/` are shipped with the repository because:
- scBridge requires GPU for training
- Files are small (~840KB total)
- Regeneration takes ~1 hour with GPU

See `fig3/results/scBridge/README.md` for regeneration instructions if needed.

### Figure 3: Cicero Analysis

The Cicero coaccessibility outputs in `fig3/results/cicero/` are shipped because:
- Computation takes 1-2 hours
- Files are ~14MB total

See `fig3/README.md` for details.

---

## Figure Legends

Each figure directory contains a markdown legend with full derivation notes:
- `fig1/Fig1_legend.md`
- `fig2/Fig2_legend.md`
- `fig3/Fig3_legend.md`
- `fig4/Fig4_legend.md`
- `rev/Revision_figures_legend.md`

See `ALL_FIGURE_LEGENDS.md` for a combined view of all figure legends.

---

## External Reference

Figure descriptions cross-checked against bioRxiv v1:
<https://www.biorxiv.org/content/10.1101/2025.04.03.647101v1>

---

## GEO Accession

All primary data is available at GEO: **GSE291468**

---

## Citation

Lyu, T., et al. (2025). "Single-cell G-quadruplex CUT&Tag reveals astrocyte-specific epigenetic regulation in mouse brain." *Genome Biology*.

---

## License

[Add your license here]

---

## Contact

[Add contact information here]

---

## Full Pipeline (from Raw Data)

### Step 1: Download Data

```bash
cd data
bash download_data.sh
# or download specific figures:
bash download_data.sh --fig1
bash download_data.sh --fig2
bash download_data.sh --fig3
bash download_data.sh --fig4
bash download_data.sh --revisions
cd ..
```

This downloads ~2GB of GEO data files.

### Step 2: Regenerate Corrected BigWigs (Optional)

The corrected bigWigs are already bundled in `data/GSE291468/`. Only regenerate if needed:

```bash
cd data/GSE291468
bash regenerate_bigwigs.sh
# or submit to cluster (recommended):
sbatch regenerate_bigwigs.sbatch
cd ../..
```

**Why this may be necessary:**
- Original GEO files used wrong effectiveGenomeSize (2150570000 vs 2652783500 for mm10)
- Missing duplicate removal step
- Causes artificial signal spikes in small cell populations

**Wait for completion** before generating figures (20-40 minutes).

### Step 3: Generate Figures

```bash
# Main figures
Rscript fig1/fig1.R
Rscript fig2/fig2.R
Rscript fig3/fig3.R
Rscript fig4/fig4.R

# Revision figures (17 scripts)
for script in rev/rev_F*.R; do
    Rscript "$script"
done
```

---

## Troubleshooting

### Missing data files
```bash
# Check what's missing
ls data/GSE291468/
# Re-download if needed
cd data && bash download_data.sh
```

### R package errors
Install required packages:
```r
install.packages(c("Seurat", "Signac", "ggplot2", "ComplexHeatmap", "GenomicRanges"))
# See individual scripts for full package lists
```

### scBridge errors
scBridge requires GPU and is precomputed. Use the CSV outputs from `fig3/results/scBridge/`:
- `umap_coembedded.csv`
- `scbridge_predictions.csv`
- `scbridge_reliability.csv`

### Path errors
Scripts use relative paths from repo root. Check this README for expected file locations.

---
