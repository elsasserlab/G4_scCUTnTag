# G-quadruplex Profiling in Complex Tissues Using Single-cell CUT&Tag

Lyu et al. (2025)

This repository contains the analysis scripts and data for our study on G-quadruplex (scG4) patterns in mouse brain using single-cell CUT&Tag. All figure scripts, generated figure outputs, and the data required to reproduce Figures 1-4 and Supplementary Figures are included. Upstream preprocessing that cannot be re-run from the repository (Cell Ranger runs, scBridge model training) is described in the respective READMEs; see "Known non-reproducible steps" below.

---

## Quick Start

### 1. Clone and Download Data

```bash
git clone <repo-url>
cd G4_scCUTnTag
git lfs install          # requires git-lfs; several data files are LFS-tracked
git lfs pull             # materialize LFS files (bigWigs, author-provided Rds)
bash data/download_data.sh
```

`git lfs pull` fetches files that are not on GEO (whole-population and
cell-type-predicted bigWigs, author-provided Seurat reference objects). The
download script warns if LFS files are still un-materialized pointers.

### 2. Generate Figures

Run all commands from the repository root (`fig1/fig1.R` and `fig2/fig2.R` resolve paths from the working directory).

```bash
# Figure 1: mESC-MEF G4 CUT&Tag analysis
Rscript fig1/fig1.R

# Figure 2: Integration with scRNA-seq reference
Rscript fig2/seurat_integration.R  # Generate intermediates first
Rscript fig2/fig2.R

# Figure 3: Astrocyte-specific G4 patterns
Rscript fig3/fig3.R
Rscript fig3/cicero.R  # Optional: coaccessibility analysis

# Figure 4: Unsorted brain neuronal analysis
Rscript fig4/fig4.R
```

---

## Repository Structure

```
github/
├── data/                    # Data files (~750MB shipped + ~2GB downloaded)
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

## Analysis Pipeline

### Main Steps

1. **Seurat processing** of G4 scCUT&Tag data from CellRanger outputs
2. **Exploratory analysis** and performance validation
3. **Integration with scRNA-Seq** data (Bartosovic et al., GSE163484)
4. **Differential G4 analysis** for cell-type-specific patterns
5. **Cicero coaccessibility** analysis for regulatory element prediction

### Main Figures

| Figure | Script | Preprocessing | Description |
| --- | --- | --- | --- |
| Fig. 1 | `fig1/fig1.R` | `fig1/seurat_workflow.R` | mESC-MEF UMAP, QC, Venn, PCA, heatmap, correlation, tracks |
| Fig. 2 | `fig2/fig2.R` | `fig2/seurat_integration.R` | Brain integration UMAP, prediction scores, overlap, markers |
| Fig. 3 | `fig3/fig3.R` | `fig3/cicero.R` | scBridge UMAP (shipped), Cicero tracks |
| Fig. 4 | `fig4/fig4.R` | `fig4/map_unsorted.R`, `data/Zeisel2018/processing_neuron_data.R` | Unsorted brain mapping |

### Additional Figures

Further scripts are in `rev/` and output to `rev/outputs/`. See `rev/paths.R` for path setup and `rev/Revision_figures_legend.md` for descriptions.

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
| GSE217860 | 3T3 G4 CUT&Tag | (not used by any figure script) |
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
| **Fig. 2** | GSM8836086 (sorted GFP+)<br>GSE163484 (scRNA reference) | `GSM8836086_GFPpos_Seurat_object.Rds`<br>`scRNA_Seq-mouse_brain.Rds`<br>`fig2/results/integration/outputs/` | CCA co-embedding, transferred cell labels, prediction scores, overlap scores, marker violins, RNA/G4 feature plots |
| **Fig. 3** | GSM8836086 (sorted GFP+)<br>GSE163484 (scRNA reference) | `fig3/results/scBridge/umap_coembedded.csv`<br>`fig3/results/scBridge/scbridge_predictions.csv`<br>`GSM8836086_Predicted_Astrocytes_RPGC.bw` | scBridge embedding, AST differential G4, AST coverage, RNA feature plots, reliability analyses |
| **Fig. 3E** | Derived from sorted GFP+ | `fig3/results/cicero/cicero_GFPsorted-{predAST,pred_nonAST}.Rds`<br>`GSM8836086_Predicted_{Astrocytes,non-Astrocytes}_RPGC.bw`<br>`PQS_scores.mm10.bw`<br>`data/cCRE/` | Directly generated Cicero co-accessibility browser tracks for five G4 loci |
| **Fig. 4** | GSM8836087 (unsorted)<br>GSM8836086 (sorted reference) | `GSM8836087_unsorted_Seurat_object.Rds`<br>`GSM8836087_unsorted_cluster_{0,1}_peaks.narrowPeak` | Query-to-reference mapping, mapped labels, prediction scores, peak-overlap Venn |
| **Fig. 4C** | GSM8836087 (unsorted clusters) | `GSM8836087_cluster{0,1}_RPGC.bw`<br>`GSM8836086_GFP_sorted_mousebrain.rpgc.bw`<br>`PQS_scores.mm10.bw` | RPGC signal profiles/heatmaps for cluster-specific peaks, PQS score profiles |
| **Fig. 4D** | GSM8836087 (unsorted) | `GSM8836087_unsorted_cluster_{0,1}_peaks.narrowPeak` | Nearest-gene annotation and enrichR GO Biological Process analysis |
| **Fig. 4F-G** | GSM8836087 + Zeisel neuron reference | `fig4/results/scBridge/`<br>`data/Zeisel2018/scRNA_Seq-Zeisel_et_al-neuron.Rds` | Precomputed scBridge neuron integration UMAP and neuron-type treemap (Rds only needed to re-run scBridge itself) |

---

## Software Requirements

### Recommended: conda environment + R package installer

```bash
conda env create -f environment.yml
conda activate scg4
Rscript install_dependencies.R
```

`environment.yml` provides R and the command-line tools used by the pipeline
(deepTools, bedtools, samtools, MACS2, tabix). `install_dependencies.R` installs
the remaining CRAN, Bioconductor and GitHub R packages.

`data/GSE291468/regenerate_bigwigs.sh` and `data/GSE198467/make_atac_bigwig.sh`
activate the conda environment named `scg4` by default. Override with
`CONDA_ENV=<name>` if you use a different environment. `fig1/fig1_deeptools.R`
uses deepTools from `PATH`; set `DEEPTOOLS_BIN=/path/to/deeptools/bin` if they are
not on `PATH`.

### R packages installed by `install_dependencies.R`

- CRAN: Seurat, Signac, tidyverse, data.table, glue, argparse, cowplot, ggpubr,
  ggrepel, ggrastr, patchwork, viridis, scales, RColorBrewer, readxl, gridExtra,
  ggplotify, eulerr, VennDiagram, hrbrthemes, Matrix, pacman, matrixStats, fields
- Bioconductor: GenomicRanges, IRanges, GenomeInfoDb, GenomicFeatures, rtracklayer,
  BSgenome, BSgenome.Mmusculus.UCSC.mm10, EnsDb.Mmusculus.v79, org.Mm.eg.db,
  ensembldb, glmGamPoi, plyranges,
  TxDb.Mmusculus.UCSC.mm10.knownGene, ChIPseeker, pqsfinder, ComplexHeatmap,
  circlize, Gviz, regioneR, cicero, monocle3
- GitHub: cnluzon/wigglescout, cnluzon/bedscout, egeulgen/DoubletFinder,
  mojaveazure/loomR
- `enrichR` was archived on CRAN and is installed from the CRAN archive
  (version 3.2) by `install_dependencies.R`

### Command-line tools

- deepTools (`computeMatrix`, `plotHeatmap`, `bamCoverage`)
- bedtools, samtools, MACS2, tabix/htslib
- `curl` and `gzip` (for `data/download_data.sh`)

### Other requirements

- **R** >= 4.3 (tested with R 4.4)
- No scBridge installation, Python environment or GPU is required; all
  scBridge-derived inputs used by Figures 3 and 4 are shipped precomputed
- Internet access for the Enrichr web API (Fig. 1G, 4D and rev F25)

---

## Known non-reproducible steps

The following upstream steps cannot be re-run from this repository; their
outputs are either shipped with the repo or regenerated with the provided
scripts:

- **Cell Ranger runs** (raw FASTQ -> fragments/peak matrices): not included;
  the Cell Ranger outputs used by the figures are downloaded from GEO or
  shipped with the repository (`data/GSE291468/GSM*/CellRanger/`).
- **scBridge model training** (GPU required): not re-runnable; all scBridge
  outputs used by Figures 3-4 are shipped in `fig3/results/scBridge/` and
  `fig4/results/scBridge/`.
- **150 bp-smoothed bigWigs** (`*.smooth150.bw`, used by Fig 3C/S1): too large
  to ship; regenerate with `bash data/GSE291468/regenerate_bigwigs.sh` before
  running Figure 3 (see `fig3/README.md`).
- **Author-provided Seurat reference objects** (Bartosovic, Zeisel): shipped
  via Git LFS; not available on GEO.

---

## Important Notes

### Figure 2 Prerequisites

Before running `fig2.R`, you must generate integration objects:

```bash
Rscript fig2/seurat_integration.R  # ~30-60 minutes
Rscript fig2/fig2.R                # ~10-20 minutes
```

See `fig2/README.md` for details.

### Figure 3: scBridge Integration

scBridge is a neural network-based single-cell multi-omics integration tool that leverages data heterogeneity ([Yunfan Li et al., 2023](https://www.nature.com/articles/s41467-023-41795-5)).
scBridge is **not** a dependency of this repository and does not need to be
installed or run: every scBridge-derived input is shipped precomputed.

- Fig. 3: `fig3/results/scBridge/*.csv`
- Fig. 4F-G: `fig4/results/scBridge/*.pdf`

No GPU, Python scBridge environment or model weights are required.

### Figure 3: Cicero Analysis

The Cicero coaccessibility outputs in `fig3/results/cicero/` are shipped because:
- Computation takes 1-2 hours
- Files are ~14MB total
- `fig3/fig3.R --only=E` uses the shipped connection tables directly

See `fig3/README.md` for details.

### Shipped intermediate results

Only the intermediate files that figure scripts read as inputs are shipped
(~22MB, tracked); everything else under `results/` is regenerable and ignored:

- `fig2/results/integration/outputs/`: prediction scores, marker table,
  anchor matrix. The two large integration objects (`G4_scRNA_integration.Rds`,
  `scRNA_Seq_Seurat_object.Rds`) are **not** shipped — generate them with
  `Rscript fig2/seurat_integration.R` (required before `fig2/fig2.R`).
- `fig3/results/scBridge/*.csv`: scBridge embedding/predictions (Fig 3B-D).
- `fig3/results/cicero/*.Rds`: coaccessibility tables (Fig 3E).
- `fig4/results/scBridge/*.pdf`: scBridge neuron integration (Fig 4F-G).

---

## Figure Legends

Each figure directory contains a markdown legend with full derivation notes:
- `fig1/Fig1_legend.md`
- `fig2/Fig2_legend.md`
- `fig3/Fig3_legend.md`
- `fig4/Fig4_legend.md`
---

## Preprint

This manuscript is available as a preprint on bioRxiv:
<https://www.biorxiv.org/content/10.1101/2025.04.03.647101>

---

## GEO Accession

All primary data is available at GEO: **GSE291468**

---

## Citation

Lyu, T., et al. (2025). "Single-cell G-quadruplex CUT&Tag reveals astrocyte-specific epigenetic regulation in mouse brain." *bioRxiv*. https://doi.org/10.1101/2025.04.03.647101


---

## Full Pipeline (from Raw Data)

### Step 1: Download Data

```bash
bash data/download_data.sh
# or download specific figures:
bash data/download_data.sh --fig1
bash data/download_data.sh --fig2
bash data/download_data.sh --fig3
bash data/download_data.sh --fig4
bash data/download_data.sh --revisions
```

This downloads ~2GB of GEO data files.

### Step 2: Regenerate Corrected BigWigs (Optional)

Updated bigWigs (include PCR duplicate removal) are already bundled in `data/GSE291468/`. Only regenerate if needed:

```bash
bash data/GSE291468/regenerate_bigwigs.sh
```

### Step 3: Generate Figures

```bash
# Main figures
Rscript fig1/fig1.R
Rscript fig2/seurat_integration.R
Rscript fig2/fig2.R
Rscript fig3/fig3.R
Rscript fig4/fig4.R
```

---

## Troubleshooting

### R package errors
Install all required packages with:
```bash
Rscript install_dependencies.R
```
See `README.md#software-requirements` for the full package list.

### scBridge errors
scBridge is not required. The inputs used by the figures are shipped in
`fig3/results/scBridge/` and `fig4/results/scBridge/`; if they are missing,
re-clone or re-download the repository (they are regular Git files, not LFS).

### LFS files are tiny "pointer" files
If bigWigs or `.Rds` files are only ~130 bytes, git-lfs did not materialize
them. Run:
```bash
git lfs install
git lfs pull
```
