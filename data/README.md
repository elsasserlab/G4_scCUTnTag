# Data Directory

This directory contains all data files required to reproduce Figures 1-4 and revision supplementary figures.

## Quick Start

### Step 1: Download all data

Most Source Data files are too large to be shipped with the repo. All necessary source data can be downloaded from GEO with a simple script:

```bash
cd data
bash download_data.sh
# or download specific figures:
bash download_data.sh --fig1
bash download_data.sh --revisions
```

This downloads all GEO data files to their proper subdirectories. Compressed files (`.gz`) are automatically decompressed during download, so you will have ready-to-use uncompressed files (e.g., `.narrowPeak`, `.broadPeak`, `.Rds`).

**IMPORTANT**

Several GEO bigwig files for GSM8836086 (predicted cell types) are replaced in this repo with
updated versions that have duplicates removed. This is important for downstream analysis. 
Following bigwig files are shipped with this repo and hence should not be overwritten with their GEO
versions (they are commented out in the download script):

GSM8836086_Predicted_Astrocytes_RPGC.bw
GSM8836086_Predicted_non-Astrocytes_RPGC.bw
GSM8836086_Predicted_Oligodendrocyte_progenitor_cells_RPGC.bw
GSM8836086_Predicted_Committed_and_Newly_formed_oligodendr_RPGC.bw
GSM8836086_Predicted_Mature_oligodendrocytes_RPGC.bw
GSM8836086_Predicted_Olfactory_ensheathing_cells_RPGC.bw
GSM8836086_Predicted_Pericytes_RPGC.bw
GSM8836086_Predicted_Vascular_endothelial_cells_RPGC.bw
GSM8836086_highly_predicted_Astrocytes_RPGC.bw
GSM8836086_Unreliable_predictions_RPGC.bw

### Step 2: Generate additional bigwigs

Generation of smoothed bigwigs (`data/GSE291468/regenerate_bigwigs.sh`) is necessary for Figure 3.
The script also can be used to regenerate the deduplicated bigwigs, but these commands are commented out by default

```bash
cd GSE291468
bash regenerate_bigwigs.sh
# or submit to cluster:
sbatch regenerate_bigwigs.sbatch LABEL OUTPUT_BW [SMOOTH_BP]
```

The regeneration script needs:
1. `GSM8836086_GFP_sorted_mouse_brain/CellRanger/fragments.tsv.gz` - Downloaded from GEO
2. `scBridge_predictions/barcodes/barcodes_*.tsv` - Authors' original cell-group barcode lists, bundled (~225KB)
3. `Seurat_cluster_barcodes/barcodes_cluster_*.tsv` - Bundled in repo (small, ~50KB total)

All auxiliary files are included in the repository and tracked by git (not gitignored).

### Step 3: Generate figures

```bash
cd ..
Rscript fig1/fig1.R
Rscript fig2/fig2.R
Rscript fig3/fig3.R
Rscript fig4/fig4.R
```

## Directory Structure

```
data/
├── download_data.sh          # Download all GEO data (run first)
├── GSE291468/regenerate_bigwigs.sh     # Optional regeneration (smooth150 pair)
├── GSE291468/regenerate_bigwigs.sbatch # SLURM submission for regeneration
├── README.md                 # This file
├── .gitignore_data           # Git ignore rules for data files
├── genome/                   # Reference genome (shipped: ~42M)
├── cCRE/                     # Mouse brain cCRE with H3K27ac (shipped: ~1M)
├── pqsfinder/                # PQS predictions (shipped: ~407M)
├── GSE291468/                # Primary scG4 CUT&Tag data (partial: metadata + corrected BWs)
├── GSE163484/                # Bartosovic scRNA-seq (downloaded from GEO)
├── GSE198467/                # Brain snATAC-seq (downloaded from GEO)
├── GSE149080/                # mESC ATAC-seq (downloaded from GEO)
├── GSE211123/                # 3T3 ATAC-seq (downloaded from GEO)
├── GSE217860/                # 3T3 G4 CUT&Tag (downloaded from GEO)
└── GSE90894/                 # Chronis mESC/MEF RNA-seq table (downloaded from GEO)
```

## Data Sources

### Primary Data (GSE291468)

Our scG4 CUT&Tag data from three experiments:
- **GSM8836086**: GFP+ sorted postnatal mouse brain
- **GSM8836087**: Unsorted postnatal mouse brain
- **GSM8836088**: mESC-MEF mixture

### Reference Data

- **GSE163484**: Bartosovic et al. mouse brain scRNA-seq (Fig 2-3 reference)
- **GSE198467**: Brain snATAC-seq (Fig 4 reference)
- **GSE149080**: mESC ATAC-seq (bulk reference)
- **GSE211123**: 3T3 ATAC-seq (bulk reference)
- **GSE217860**: 3T3 G4 CUT&Tag (bulk reference)
- **GSE90894**: Chronis mESC/MEF RNA-seq table

### Genome Resources

- **genome/**: mm10 chromosome sizes, GTF annotation
- **cCRE/**: Mouse brain cCRE with H3K27ac signal
- **pqsfinder/**: G-quadruplex predicting sequences (PQS)

## File Naming

Downloaded files retain their GEO accession numbers (GSM prefix) for traceability.
Some files are renamed for convenience (e.g., `GSM8836088_cluster0_RPGC.bw`).

See `../FIGURE_REPRODUCIBILITY.md` for detailed file-to-figure mapping.

## BigWig Regeneration Details

The bundled GSM8836086 bigwigs were regenerated with `GSE291468/regenerate_bigwigs.sh`, fixing three issues in the GEO originals:

1. **Wrong effectiveGenomeSize**: 
   - Original: 2,150,570,000 (incorrect)
   - Corrected: 2,652,783,500 (mm10 effective genome size from ENCODE)

2. **Missing duplicate removal**:
   - Original: No deduplication
   - Corrected: `samtools markdup -r` before bigwig generation

This is critical for the predicted cell type bigwigs (Astrocytes, non-Astrocytes, etc.)
which show artificial spikes due to PCR duplicates in small cell populations.

The Seurat cluster bigwigs (cluster_0-3) are also regenerated for consistency.

## What's Shipped vs Downloaded

### Shipped with the repository (~651 MB)
These files are tracked by git and synced when you clone the repo:

- **Reference genome** (`genome/`, ~42M): mm10 chrom.sizes, GTF, eG4, TSS regions, TPM table
- **cCRE** (`cCRE/`, ~1M): Li et al. mouse brain cCRE with H3K27ac
- **PQS predictions** (`pqsfinder/`, ~407M): Genome-wide PQS scores (bed + bigWig)
- **GSE291468 metadata** (~54M): CellRanger CSV/TSV files, Seurat cluster barcodes
- **Corrected bigWigs** (~158M): 9 GSM8836086 bigWigs with duplicates removed
- **Scripts**: `download_data.sh`, `regenerate_bigwigs.sh`, `filter_pqs.sh`, etc.

### Downloaded from GEO (~2+ GB)
Run `bash download_data.sh` to download these before generating figures:

- **GSE163484**: Bartosovic scRNA-seq
- **GSE198467**: Brain snATAC-seq
- **GSE149080**: mESC ATAC-seq
- **GSE211123**: 3T3 ATAC-seq
- **GSE217860**: 3T3 G4 CUT&Tag
- **GSE90894**: Chronis RNA-seq table
- **GSE291468**: Most bigWigs, Seurat objects, narrowPeak files

See `.gitignore_data` for the complete list of gitignored files.

## Troubleshooting

### Missing files after download
Check `download_data.sh` output for failed downloads. Some GEO URLs may require manual download.

### BigWig spikes still visible
The supplied `GSE291468/GSM8836086_*_RPGC.bw` files are already fixed; do not re-download the GEO originals (masked in download_data.sh).

### Path errors in figure scripts
Scripts should use relative paths from repo root. Check `../FIGURE_REPRODUCIBILITY.md` for expected file locations.