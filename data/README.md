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
GSM8836086_Unreliable_predictions_RPGC.bw

### Step 2: Generate additional bigwigs

Generation of smoothed bigwigs (`data/GSE291468/regenerate_bigwigs.sh`) is necessary for Figure 3.
The script also can be used to regenerate the deduplicated bigwigs, but these commands are commented out by default

```bash
# From the repository root
bash data/GSE291468/regenerate_bigwigs.sh
```

The regeneration script needs:
1. `GSM8836086_GFP_sorted_mouse_brain/CellRanger/fragments.tsv.gz` - Downloaded from GEO
2. `scBridge_predictions/barcodes/barcodes_*.tsv` - Authors' original cell-group barcode lists, bundled (~225KB)
3. `Seurat_cluster_barcodes/barcodes_cluster_*.tsv` - Bundled in repo (small, ~50KB total)

All auxiliary files are included in the repository and tracked by git (not gitignored).

### Step 3: Generate figures

```bash
Rscript fig1/fig1.R
Rscript fig2/seurat_integration.R
Rscript fig2/fig2.R
Rscript fig3/fig3.R
Rscript fig4/fig4.R
```

## Directory Structure

```
data/
├── download_data.sh          # Download all GEO data (run first)
├── GSE291468/regenerate_bigwigs.sh     # Optional regeneration (smooth150 pair)
├── README.md                 # This file
├── ../.gitignore             # Git ignore rules for data files
├── genome/                   # Reference genome (shipped: ~42M)
├── cCRE/                     # Mouse brain cCRE with H3K27ac (shipped: ~1M)
├── pqsfinder/                # PQS predictions (shipped: ~407M)
├── GSE291468/                # Primary scG4 CUT&Tag data (partial: metadata + corrected BWs)
├── GSE163484/                # Bartosovic scRNA-seq (downloaded from GEO)
├── Zeisel2018/               # Zeisel neuron scRNA-seq Rds + processing script (author-provided)
├── GSE198467/                # Brain snATAC-seq (downloaded from GEO)
├── GSE149080/                # mESC ATAC-seq (bigWig downloaded; narrowPeak shipped)
├── GSE211123/                # 3T3 ATAC-seq (bigWig downloaded; narrowPeak shipped)
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
- **Zeisel2018**: Zeisel et al. neuron scRNA-seq (Fig 4 reference)
- **GSE198467**: Brain snATAC-seq (Fig 4 reference). The GEO Seurat object
  (`GSE198467_ATAC_Seurat_object_clustered_renamed.Rds`) contains only broad
  peak windows (median ~3.7 kb) with no summit offsets. Consensus peaks with
  real summits are re-called from the raw fragments file
  (`GSE198467_ATAC_fragments.tsv.gz`, downloaded from GEO, gitignored) by
  `data/GSE198467/reprocess_atac_peaks.sh` (MACS2, q=1e-6, ~93k peaks):
  `GSE198467_atac_consensus_peaks.narrowPeak`.
- **GSE149080**: mESC ATAC-seq (bulk reference; RPGC bigWig downloaded, narrowPeak shipped because GEO has no peak file)
- **GSE211123**: 3T3 ATAC-seq (bulk reference; RPGC bigWig downloaded, narrowPeak shipped because GEO has no peak file)
- **GSE217860**: 3T3 G4 CUT&Tag (bulk reference)
- **GSE90894**: Chronis mESC/MEF RNA-seq table

### Genome Resources

- **genome/**: mm10 chromosome sizes (the GENCODE M25 GTF is downloaded by `download_data.sh`)
- **cCRE/**: Mouse brain cCRE with H3K27ac signal
- **pqsfinder/**: G-quadruplex predicting sequences (PQS)

## File Naming

Downloaded files retain their GEO accession numbers (GSM prefix) for traceability.
Some files are renamed for convenience (e.g., `GSM8836088_cluster0_RPGC.bw`).

See the "GEO File-to-Figure Map" in `../README.md` for detailed file-to-figure mapping.

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

- **Reference genome** (`genome/`, ~42M tracked): mm10 chrom.sizes, eG4, TSS regions, TPM table (the GENCODE GTF is downloaded by `download_data.sh`)
- **cCRE** (`cCRE/`, ~1M): Li et al. mouse brain cCRE with H3K27ac
- **PQS predictions** (`pqsfinder/`, ~407M): Genome-wide PQS scores (bed + bigWig)
- **GSE291468 auxiliary tables** (~0.25M): scBridge cell-group barcodes and Seurat cluster barcodes
- **Bulk ATAC narrowPeaks** (~11M): GSE149080 and GSE211123 peak calls used by Fig 1 and revision F18 (not distributed by GEO)
- **Corrected bigWigs** (~158M): 9 GSM8836086 bigWigs with duplicates removed
- **Scripts**: `download_data.sh`, `regenerate_bigwigs.sh`, `filter_pqs.sh`, etc.

### Downloaded from GEO (~2+ GB)
Run `bash download_data.sh` to download these before generating figures:

- **GSE163484**: Bartosovic scRNA-seq (counts CSV + annotation CSV)
- **GSE198467**: Brain snATAC-seq
- **GSE149080**: mESC ATAC-seq RPGC bigWig (narrowPeak is shipped)
- **GSE211123**: 3T3 ATAC-seq RPGC bigWig (narrowPeak is shipped)
- **GSE217860**: 3T3 G4 CUT&Tag
- **GSE90894**: Chronis RNA-seq table
- **GSE291468**: Most bigWigs, Seurat objects, narrowPeak files
- **genome/**: GENCODE M25 mm10 GTF annotation

### Author-provided Seurat objects (not on GEO)

Two processed Seurat reference objects used in this study are **not available
on GEO** and therefore cannot be downloaded by `download_data.sh`. They were
provided directly by the study authors and must be obtained by contacting them.
Place them at the locations indicated below (both are gitignored due to file
size). A direct download link (e.g. data-repository/Zenodo) will be added to
`download_data.sh` once available.

| File | Source study | Place at |
|------|--------------|----------|
| `scRNA_Seq-mouse_brain.Rds` | Bartosovic et al., *Nat Biotechnol* 2021 (GSE163484); author-provided Seurat object | `GSE163484/scRNA_Seq-mouse_brain.Rds` |
| `scRNA_Seq-Zeisel_et_al-neuron.Rds` | Zeisel et al., *Cell* 2018, 174:999-1014, DOI: 10.1016/j.cell.2018.06.021 (mousebrain.org); processed Seurat object of the neuron population, built by `Zeisel2018/processing_neuron_data.R` from the (now offline) level-2 neuron loom files | `Zeisel2018/scRNA_Seq-Zeisel_et_al-neuron.Rds` |

See `../.gitignore` for the complete list of gitignored files.

## Troubleshooting

### Missing files after download
Check `download_data.sh` output for failed downloads. Some GEO URLs may require manual download.

### BigWig spikes still visible
The supplied `GSE291468/GSM8836086_*_RPGC.bw` files are already fixed; do not re-download the GEO originals (masked in download_data.sh).

### Path errors in figure scripts
Scripts should use paths derived from the repository root. See `../README.md` for expected file locations.