RNAseq 6w vs 12w +/- MMTV-MYC mouse merged data June-September 2024

# First attempt
k693_rerun_analysis_GS.R --> analysis with interaction myc-status vs timepoint

  - coldata setup: 

from: data/FULL.DAT.COL.DATA.txt
columns:
- sample
- group
- myc_status
- timepoint
- prefix
- suffix
- numeric_part
- alphabetical_part

- counts: cts
- issue?: 12W has about half total counts comapred to 6W

- dds: design = ~ timepoint * myc_status) #this expands to ~ timepoint + myc_status + timepoint:myc_status
  
- coverage: 19K genes >10 counts
  
- issue?: data are noisy - none of the variance stabilising transformations, followed by distance measures or PCA ("log2(x + 1)", "vst", "rlog") reveal grouping by the actual experimental groups
PC1 and 2 top genes hint to immune cells...
  
- further exploration with variancePartition: only small part of the variation is explained by any of the parameters: see coldata: (1|myc_status) + (1|timepoint) + (1|numeric_part) + (1|alphabetical_part) + (1|myc_status) + (1|timepoint) (1|prefix) (1|suffix)

- moving on for DGE anyway: 

- altogether this analysis provides a list, using the interaction terms to create lists of genes whihc are significant at 6W and 12W. The genes are used in Gprofiler, but it is is difficult to interpret the changes.

- I have also tried a Cytoscpae analysis to visualise the gene lists, but it is not yet completed and/or informative.

- Conclusion: this file is used as the base of the second analysis where both the interactiion and group-based design is used to create gene lists, to understand the difference in the Myc effect between 6W and 12W timepoints.

# Second attempt
Myc_timecourse_analysis_GS.R and the sandbox version, later updated to:

# MYC-Dependent Temporal Transcriptome Analysis in Mouse Tissues

This repository contains the analysis pipeline and results for a transcriptomic study comparing MYC-positive and MYC-negative samples across two timepoints (6 weeks and 12 weeks). It uses RNA-seq data processed with DESeq2 to identify MYC-dependent expression programs and their temporal changes, with visualizations based on gene set heatmaps.

---

## 📁 Project Structure

<pre>
myc-temporal-analysis/
├── data/                      # Raw input files (counts, metadata, gene sets)
├── results/                  # RDS results files (LFCs, annotations)
├── outputs/heatmaps_int/     # Heatmaps generated from interaction design
├── scripts/                  # Modular R scripts used for pipeline
├── functions/                # Utility functions (heatmap generation, etc.)
├── 00_setup_packages.R       # Package loading and environment setup
├── run_all_scripts.R         # Wrapper to run full pipeline
└── README.md                 # This file
</pre>

---

## 🔬 Analysis Overview

### Experimental Design
- 4 experimental conditions:
  - Timepoints: 6W vs 12W
  - Genotype: MYC-positive vs MYC-negative
- 6 biological replicates per group

### DESeq2 Designs
1. **Interaction model**: `~ timepoint * myc_status`
2. **Group-based design**: 4-level `group` factor with contrast `"12W_pos vs 6W_pos"`

### Gene Sets Used
- 🧬 [MitoCarta 3.0](https://www.broadinstitute.org/mitocarta) mitochondrial pathways
- 🔬 MYC signature gene sets (`.gmx` format, converted to mouse orthologs)
- 📄 Felsher integrative MYC target list

---

## ✅ Implemented Features

- 🧪 DESeq2 pipeline with `ashr` shrinkage and `IHW` p-value filtering
- 🧬 Classification of MYC effect over time:
  - `direct_Myc_reduction/increase`
  - `baseline_driven_reduction/increase`
  - `no_change`
- 📊 Heatmap generation:
  - Log2FC matrix with annotations
  - Z-scaled group expression heatmaps
  - Per-gene annotations:
    - MYC temporal classification
    - Group-based significance
- 🧭 Modular structure with individual scripts for:
  - Data loading
  - Model fitting
  - LFC classification
  - Heatmap generation
  - QC exploration

---

## 📈 Output Preview

Each gene set generates:
- One PDF per heatmap type (shrunk and raw)
- Row-annotated by MYC classification and group-level DE
- OXPHOS subunits additionally annotated by complex and mtDNA encoding

---

## 📦 Setup

Required R packages (auto-installed in `00_setup_packages.R`):
- `DESeq2`, `apeglm`, `ashr`, `IHW`
- `ComplexHeatmap`, `circlize`, `ggplot2`, `biomaRt`, `edgeR`, etc.

Run full pipeline:

```r
source("run_all_scripts.R")

NOTE: this has not been fully implemented yet. Testing on OXPHOS (scripts/07_heatmap_oxphos_annotated.R) currently, not satisfying heatmap structure and annotations
aims to develop:

Create: Fig 1 Myc-mito paper: mitochondrial adaptation defines early tumourigenesis and late progression


- To understand early evolution, looked at timecourse of early myc regulated genes. 
- Histology of tumours - 6w and 12W, apoptosis and proliferation
- What do we see at the gene level?
- What are the gene sets changing most? - pathway analysis
- gene sets from hocklebbery? - better to define it ourselves, eventually comapre to that
- Clusters of gene sets - mitochondrial most affected
- While Myc goes on, some mitochondrial genes are reduced - Complex I and Complex IV
- What is special about Complex I - a lot.

- where to get the myc-ER 6 vs 12 in?
- 
go to in vitro: MYAZ has the same
---
What to do:
- Find best way to cluster genes: try clustering to enhance 12to6 downreg (12_6_d geneset) - use old PGT4o or a new Claude4?
- apoptotic pathway?
- p19 pathway?
- separate mtDNA - why LFC does not fit?
- Pathway analysis













  

