# scG4 vs scATAC Accessibility Analysis

**Question from reviewer**: In how far does G4 signal merely reflect an accessibility artifact?

**Answer**: **A large fraction of G4 signal is independent of chromatin accessibility.**

## Figures

### F28: Peak Overlap Venn Diagram (Euler)
- **Script**: `rev_F28_scG4_vs_scATAC_venn.R`
- **Output**: `F28_peak_overlap_venn.pdf`, `F28_peak_overlap_stats.csv`
- **Key finding**: 23% of scG4 peaks occur at sites NOT accessible by ATAC-seq; 77% of scG4 peaks overlap ATAC peaks while only 32% of the much larger ATAC peak set overlaps scG4.
- **Note**: Diagram areas are proportional to reduced genomic regions (not raw peak counts). Original peak totals shown outside circles.

### F30: Scatter Plot (Union vs Intersection)
- **Script**: `rev_F30_scG4_vs_scATAC_scatter.R`
- **Output**: `F30_scatter_union_intersection.pdf`, `F30_intersection_data.csv`, `F30_union_data.csv`
- **Key finding**: Weak overall correlation (UNION r = 0.274); correlation at overlapping regions is even lower (INTERSECTION r = 0.166), i.e. even where scG4 and scATAC peaks coincide, G4 levels poorly track accessibility.

## Key Results

### Peak Counts (consistent MACS2 narrowPeak sets)
| Metric | Value | Interpretation |
|--------|-------|----------------|
| scG4 peaks (merged cluster-narrowPeak union) | 26,982 (median 285 bp) | Canonical-chr union of GFP+ cluster peaks |
| scATAC peaks (re-called consensus) | 92,926 (median 466 bp) | MACS2 narrowPeak re-call |
| scG4-only | 6,113 (23%) | **G4 without accessibility** |
| scATAC-only | 43,708 (68%) | Accessibility without G4 |

### Reduced Regions (for Venn/Euler diagram)
| Metric | Value | Interpretation |
|--------|-------|----------------|
| Total unique regions | 70,364 | After merging overlapping peaks |
| scG4 signal | 26,656 | Regions with scG4 signal |
| scATAC signal | 64,251 | Regions with scATAC signal |
| Overlap | 20,543 (77%/32%) | Regions with both signals |

### Correlation Results (continuous coverage, log scale)
| Metric | Value | Interpretation |
|--------|-------|----------------|
| **F30 INTERSECTION r** | **0.166** | Weak correlation at shared peaks |
| **F30 UNION r** | **0.274** | Weak correlation across all regions |

## Conclusion

Using consistent MACS2 narrowPeak definitions for both assays:

1. **23% of G4 peaks are at sites with no ATAC signal** - a substantial G4 fraction occurs without chromatin accessibility
2. **F30: only weak correlations** - union and intersection correlations are both low (r <= 0.27); even where peaks overlap, G4 levels do not track accessibility levels

This demonstrates that scG4 CUT&Tag cannot be explained as a mere accessibility artifact.

## Data Sources

- **scG4 peaks**: per-cluster MACS2 narrowPeaks `GSM8836086_GFPpos_cluster_{0..3}_peaks.narrowPeak`, merged into a union (canonical chr only)
- **scG4 coverage**: `GSM8836086_GFP_sorted_mousebrain.rpgc.bw` (RPGC, 10 bp bins)
- **scATAC peaks**: re-called consensus `GSE198467_atac_consensus_peaks.narrowPeak`
- **scATAC coverage**: `GSE198467_ATAC_fragments.RPGC.bw` (RPGC, 10 bp bins)
- scG4: GSM8836086 (GFP+ sorted mouse brain), 2,402 cells; scATAC: GSE198467 (mouse brain), 4,960 cells