# MYC Mouse RNAseq Analysis

## Project Overview

This project investigates the effect of the **Myc oncogene** in early breast tumourigenesis using a mouse model. Myc is selectively and constitutively expressed in breast epithelial cells using the MMTV promoter.

### Experimental Design

- **Model**: MMTV-Myc transgenic mice
- **Groups**: 4 experimental groups (n=6 per group)
  - `6W_neg`: 6 weeks, Myc negative (control)
  - `6W_pos`: 6 weeks, Myc positive
  - `12W_neg`: 12 weeks, Myc negative (control)
  - `12W_pos`: 12 weeks, Myc positive
- **Total samples**: 24
- **Data type**: Bulk RNAseq

### Research Questions

1. What is the transcriptional effect of Myc expression in breast epithelial cells?
2. How does the Myc effect evolve between 6 and 12 weeks?
3. Which pathways and gene sets are affected by Myc expression?

---

## Analysis Workflow

### Scripts

| Script | Description |
|--------|-------------|
| `00_setup_packages.R` | Load/install required R packages |
| `01_load_data.R` | Load count data, create DESeq2 object, load gene sets |
| `02_qc.R` | Quality control: PCA, sample distances, variance transforms |

### Directory Structure

```
myc_mouse/
├── data/
│   ├── coldata.csv                    # Sample metadata
│   ├── FULL.DAT.csv                   # Raw count matrix
│   ├── mitocarta_pathways.csv         # Mitocarta gene sets
│   ├── myc_signature_genesets.gmx     # MYC signature gene sets (Felsher)
│   └── felsher_integrative_signature.csv
├── scripts/
│   ├── 00_setup_packages.R
│   ├── 01_load_data.R
│   └── 02_qc.R
├── results/
│   ├── dds_int.rds                    # DESeq2 object (interaction design)
│   ├── count_matrix.rds
│   ├── coldata.rds
│   ├── gene_sets_list.rds
│   └── ortholog_table.rds             # Cached human-mouse orthologs
├── outputs/
│   └── qc/                            # QC plots (PDF)
└── README.md
```

---

## QC Summary

### Library Size (Total Counts)

| Group | Mean Counts | SD | Min | Max |
|-------|-------------|-----|-----|-----|
| 6W_neg | 22.6M | 2.5M | 18.8M | 25.3M |
| 6W_pos | 23.8M | 3.3M | 19.7M | 29.0M |
| 12W_neg | 12.9M | 1.8M | 10.2M | 14.7M |
| 12W_pos | 13.5M | 2.6M | 9.9M | 16.1M |

**Note**: 6W samples have approximately 2x more reads than 12W samples. This is a technical/batch effect that DESeq2's size factor normalization addresses.

### Size Factors

| Group | Mean | SD |
|-------|------|-----|
| 6W_neg | 1.30 | 0.25 |
| 6W_pos | 1.29 | 0.24 |
| 12W_neg | 0.81 | 0.21 |
| 12W_pos | 0.83 | 0.20 |

Size factors reflect the library size differences between timepoints. Within-group variation is modest, indicating consistent library preparation.

### PCA Analysis

- **PC1 (36%)**: Captures timepoint effect (6W vs 12W)
- **PC2 (14%)**: Captures Myc status effect, particularly at 12W

Key observations:
- **12W_pos** separates clearly from other groups on PC2, suggesting a stronger Myc transcriptional effect at 12 weeks
- **6W_pos** clusters tightly, with modest separation from 6W_neg
- **6W_neg and 12W_neg** (controls) show considerable overlap across timepoints
- No major outliers requiring removal

### Sample Distance Heatmap

- Samples cluster primarily by timepoint
- 12W_pos samples show more internal variability
- No obvious outliers

---

## QC Conclusions

1. **Data quality is acceptable** — no major outliers, size factors are reasonable
2. **Strong timepoint effect** — 6W vs 12W is the dominant source of variation (driven partly by library size differences)
3. **Myc effect visible at 12W** — 12W_pos separates from 12W_neg on PC2, suggesting Myc-driven transcriptional changes become more pronounced over time
4. **Interaction model is appropriate** — the Myc effect appears to differ between timepoints, supporting the `~ timepoint * myc_status` design
5. **Group-based comparisons also warranted** — direct 12W_pos vs 6W_pos comparison will capture progressive Myc effects

---

## Statistical Design

### Interaction Model

```r
design = ~ timepoint * myc_status
```

This model tests:
- Main effect of timepoint (12W vs 6W)
- Main effect of Myc status (pos vs neg)
- **Interaction**: whether the Myc effect differs between timepoints

### Group-Based Model (planned)

```r
design = ~ group
```

Allows direct pairwise comparisons (e.g., 12W_pos vs 6W_pos) to capture progressive Myc effects.

---

## Gene Sets

The analysis includes curated gene sets for pathway-level interpretation:

- **Mitocarta 3.0**: Mitochondrial pathways (OXPHOS, TCA cycle, FAO, etc.)
- **MYC signatures**: Multiple gene sets from Felsher et al. (2022)
- **Apoptosis**: Pro- and anti-apoptotic genes

Human gene symbols are mapped to mouse orthologs via biomaRt (cached in `results/ortholog_table.rds`).

---

## Requirements

### R Packages

**CRAN:**
- here, dplyr, tibble, readr, stringr, purrr, magrittr
- ggplot2, ggstatsplot, pheatmap, RColorBrewer
- grid, gridExtra, reshape2, ggrepel

**Bioconductor:**
- DESeq2, biomaRt, org.Mm.eg.db, AnnotationDbi
- ComplexHeatmap, EnhancedVolcano, fgsea, msigdbr
- PoiClaClu, vsn, sva

---

## Usage

```r
# From the project root directory
library(here)

# Run scripts in order
source(here("scripts", "00_setup_packages.R"))
source(here("scripts", "01_load_data.R"))
source(here("scripts", "02_qc.R"))
```

---

## Authors

[Add collaborator information]

## Date

Analysis updated: February 2026
