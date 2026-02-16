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
| `03_deseq_results_qc.R` | DESeq2 results extraction (raw + shrunken LFCs) with IHW, MA plots |
| `04_fgsea_pathway_analysis.R` | fGSEA pathway analysis using Wald statistic ranking |
| `05_fgsea_visualisation.R` | fGSEA visualisation: dot plots, bar charts, enrichment plots |

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
│   ├── 02_qc.R
│   └── 03_deseq_results_qc.R          # Includes extended interaction QC
├── results/
│   ├── dds_int.rds                    # DESeq2 object (interaction design)
│   ├── dds_int_run.rds                # DESeq2 object after DESeq() run
│   ├── dds_group_run.rds              # DESeq2 object (group design) after run
│   ├── interaction_results.rds        # All interaction model results
│   ├── group_results.rds              # All group model results
│   ├── extended_qc_summary.rds        # Extended QC statistics
│   ├── lfc_comparison_timepoint.rds   # LFC comparison data for timepoint effects
│   ├── count_matrix.rds
│   ├── coldata.rds
│   ├── gene_sets_list.rds
│   └── ortholog_table.rds             # Cached human-mouse orthologs
│   └── fgsea_results.rds              # fGSEA results (Wald statistic ranking)
├── outputs/
│   ├── qc/                            # Initial QC plots (PDF)
│   ├── deseq_qc/                      # MA plots, extended QC, and results summary
│   │   ├── MA_*.pdf                   # MA plots for each contrast
│   │   ├── extended_qc_interaction_analysis.pdf
│   │   └── results_summary.csv
│   └── fgsea/                         # fGSEA visualisations
│       ├── dotplot_*.pdf              # Dot plots by pathway category
│       ├── barplot_*_comparison.pdf   # Comparative bar plots (Myc+ vs Myc-)
│       ├── enrichment_*.pdf           # Running enrichment score plots
│       ├── category_summary.pdf       # Pathway classification summary
│       └── heatmap_top_pathways.pdf   # NES heatmap
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

## DESeq2 Results and QC

### Analysis Scripts

| Script | Description |
|--------|-------------|
| `03_deseq_results_qc.R` | Generate DESeq2 results (raw + shrunken LFCs) with IHW filtering, produce MA plots |

### Results Summary

Results were extracted from two model types:

1. **Interaction model** (`~ timepoint * myc_status`): Tests main effects and interaction
2. **Group model** (`~ group`): Direct pairwise comparisons between groups

**IHW** (Independent Hypothesis Weighting) was used for p-value adjustment, using mean expression as the covariate.

### Summary Statistics (padj < 0.1)

| Contrast | Significant | Up | Down |
|----------|-------------|-----|------|
| Myc effect at 6W | 2,777 | 1,955 | 822 |
| Myc effect at 12W | 239 | 193 | 46 |
| Timepoint effect (Myc−) | 1,867 | 950 | 917 |
| Timepoint effect (Myc+) | 2,623 | 1,078 | 1,545 |
| Interaction (Myc × Time) | 0 | 0 | 0 |
| Group: 12W_neg vs 6W_neg | 1,869 | 952 | 917 |
| Group: 12W_pos vs 6W_pos | 2,619 | 1,081 | 1,538 |

### Key Observations

1. **Strong Myc effect at 6W** (2,777 genes), but **weak at 12W** (239 genes) in the interaction model
2. **Interaction term yields 0 significant genes** — no detectable differential Myc effect between timepoints
3. **Timepoint effect (Myc−) ≈ Group 12W_neg vs 6W_neg** — as expected, these are equivalent comparisons (1,867 vs 1,869 genes, near-identical)
4. **Timepoint effect (Myc+) ≈ Group 12W_pos vs 6W_pos** — also equivalent (2,623 vs 2,619 genes)

### MA Plot QC Findings

- **Raw vs shrunken LFCs**: `ashr` shrinkage is well-behaved for main effects
- **Interaction term**: Raw results show variance, but shrinkage collapses LFCs to near-zero, explaining the lack of significant hits
- **Low-count genes**: Appropriately shrunk toward zero, especially at 12W where library sizes are smaller

---

## Extended QC: Understanding the Interaction Term

### The Puzzle

A key puzzle emerges from the results:

- **No significant interactions** (0 genes at padj < 0.1)
- **Yet ~2,600 genes differ between 12W_pos and 6W_pos**

If so many genes change over time in Myc+ samples, why doesn't the interaction term detect differential Myc effects?

### Clarification: DESeq2 Normalization Handles Library Size

DESeq2's size factors correctly account for the ~2x library size difference between 6W and 12W samples. The group comparisons (e.g., 12W_pos vs 6W_pos) are **not confounded** by library size—they capture true biological differences (developmental changes, tumour progression) plus any residual technical variation beyond what size factors address.

### Resolution: Timepoint Effects Are Highly Correlated

We investigated the overlap between timepoint effects in Myc+ and Myc− samples:

| Metric | Value |
|--------|-------|
| Significant in Myc+ (timepoint effect) | 2,623 |
| Significant in Myc− (timepoint effect) | 1,867 |
| **Overlap** | **961** (51% of Myc−) |
| Unique to Myc+ | 1,662 |
| Unique to Myc− | 906 |

**LFC correlations:**
- Genes significant in both: r = **0.96** (near-perfect)
- All genes significant in either: r = **0.74**

This explains the lack of interaction: genes changing over time follow **the same direction and similar magnitude** in both Myc+ and Myc− samples. The interaction term (difference-of-differences) is therefore close to zero.

### Interaction P-value Distribution

The interaction term p-value histogram is **near-uniform**, indicating no enrichment of true signal:

- Minimum padj: **0.16** (no genes below 0.1)
- Genes with p < 0.05: **~900** (expected by chance: ~925)

However, for the 1,662 genes "unique to Myc+ timepoint effect":
- **357 have interaction p < 0.05** (expected by chance: ~83)
- **0 survive FDR correction** (padj < 0.1)

This suggests a **weak signal** exists (more than expected by chance), but effect sizes are too small and/or variance too high for reliable detection after multiple testing correction.

---

## Interpretation: Interaction Term vs Group Comparisons

Understanding the difference between the interaction term and the group comparisons is critical for this analysis.

### What the Interaction Term Tests

The interaction term (`timepoint12W.myc_statuspos`) asks:

> **"Does the Myc effect differ between 6W and 12W?"**

Mathematically:
```
Interaction = (Myc effect at 12W) − (Myc effect at 6W)
            = (12W_pos − 12W_neg) − (6W_pos − 6W_neg)
```

A significant interaction would mean that Myc activates (or represses) a gene **differently** depending on the timepoint. For example:
- A gene upregulated 4-fold by Myc at 6W but only 1.5-fold at 12W → negative interaction
- A gene unchanged by Myc at 6W but strongly repressed at 12W → negative interaction

**Result**: No significant interactions detected. This means the **Myc transcriptional program is stable between timepoints** — genes affected by Myc at 6W are affected similarly at 12W.

### What the Group Comparison Tests

The comparison `12W_pos vs 6W_pos` asks:

> **"What genes differ between Myc+ samples at 12W versus 6W?"**

This is a **simple difference**, capturing:
1. **Developmental/tissue changes** that occur from 6W to 12W
2. **Tumour progression** effects in Myc+ samples
3. Any combination of the above

Critically, this comparison does **not** tell us whether observed changes are Myc-specific or would occur anyway in Myc− controls.

**Result**: 2,619 significant genes. This is similar to the timepoint effect in Myc+ samples (2,623 genes), confirming equivalence.

### Key Insight: Most Timepoint Changes Are Shared

The group comparisons reveal that:
- **12W_neg vs 6W_neg**: 1,869 genes (developmental baseline)
- **12W_pos vs 6W_pos**: 2,619 genes

The ~750 additional genes in the Myc+ comparison could reflect:
1. Myc-specific temporal dynamics (true interaction signal, too weak to detect)
2. Power differences (Myc+ comparison may have slightly different variance)
3. Borderline genes just crossing significance in one comparison

The **high correlation** (r = 0.74–0.96) between Myc+ and Myc− timepoint effects confirms that most of these changes are **shared developmental effects**, not Myc-specific.

### Summary Table: Which Comparison to Use

| Biological Question | Comparison to Use |
|---------------------|-------------------|
| Does Myc change gene X expression? | Myc effect at 6W or 12W |
| Does Myc affect gene X differently over time? | Interaction term |
| How do Myc+ tumours change from 6W to 12W? | Group: 12W_pos vs 6W_pos |
| Is a change Myc-specific vs developmental? | Compare 12W_pos/6W_pos to 12W_neg/6W_neg |
| What's the baseline developmental effect? | Timepoint effect (Myc−) or Group: 12W_neg vs 6W_neg |

### Biological Conclusions

1. **The Myc transcriptional program is largely stable between 6W and 12W**
   - No significant interactions = the Myc effect doesn't change dramatically over time
   
2. **Timepoint effects are predominantly developmental**
   - Most genes changing from 6W to 12W do so similarly in both Myc+ and Myc− samples
   
3. **Genes "unique" to Myc+ timepoint are likely borderline cases**
   - Not true Myc-specific temporal dynamics, but genes near the significance threshold

### Practical Recommendations

1. **Use the GROUP MODEL for clear pairwise comparisons**
2. **Reserve the interaction term for hypothesis-driven checks** on specific candidate genes
3. **Consider relaxed thresholds (padj < 0.2)** for exploratory interaction analysis if needed
4. **Interpret timepoint differences cautiously** — most are shared between Myc+ and Myc−

---

## fGSEA Pathway Analysis

### Script

| Script | Description |
|--------|-------------|
| `04_fgsea_pathway_analysis.R` | Gene set enrichment analysis for Myc+ vs Myc- progression |

### Strategy

Given the lack of significant interaction term hits, we use a comparative fGSEA approach:

1. **Run fGSEA on 12W_pos vs 6W_pos** (Q1: How do Myc+ tumours change over time?)
2. **Run fGSEA on 12W_neg vs 6W_neg** (Q3: What's the baseline developmental effect?)
3. **Compare enrichment between comparisons** (Q2: Which changes are Myc-specific?)

### Ranking Metric: Wald Statistic

We use the **Wald statistic** (`stat` column from DESeq2) for ranking genes in fGSEA:

```
Wald = log2FoldChange / lfcSE
```

**Why Wald over sign(LFC) × -log10(p)?**

| Metric | Wald statistic | sign(LFC) × -log10(p) |
|--------|----------------|------------------------|
| Source | Single model quantity | Derived from LFC + pvalue |
| Properties | Already signed, incorporates effect size and precision | Amplifies extreme values |
| Correlation | — | r = 0.95 with Wald |

Empirical comparison showed:
- **Wald detected 11 additional pathways** (44 vs 33 significant in Myc+)
- **No pathways lost** — all sign(LFC) × -log10(p) hits also significant with Wald
- Higher sensitivity for pathways with moderate but consistent effects

### Gene Sets Analysed (89 total)

| Source | Sets | Description |
|--------|------|-------------|
| MitoCarta 3.0 | 22 | Mitochondrial pathways |
| MYC signatures | 17 | Felsher et al. (2022) + others |
| Apoptosis | 2 | Pro/anti-apoptotic genes |
| MSigDB Hallmark | 50 | Canonical pathway collection |

### Output

Results saved to `results/fgsea_results.rds` containing:
- Individual fGSEA results for Myc+ and Myc- comparisons
- Combined comparison table with pathway classification
- Ranked gene lists

### Interpretation Notes

#### Temporal vs cross-sectional comparisons

The fGSEA results reflect **changes over time within each genotype** (12W vs 6W), 
not differences between Myc+ and Myc- at any given timepoint. A pathway classified 
as "Myc+ specific" shows significant temporal change only in Myc+ tumours, but 
this doesn't confirm the pathway differs between genotypes at 12W.

| Observation | What it means | What it doesn't mean |
|-------------|---------------|----------------------|
| Negative NES in both genotypes | Both decrease over time | Nothing about absolute levels at 12W |
| "Myc+ specific" | Significant change only in Myc+ | Not necessarily direct Myc regulation |
| "Opposite effects" | Genotypes diverge over time | Could be direct or indirect |

#### MYC target signatures decrease despite stable Myc mRNA

All MYC target gene sets show significant negative enrichment in Myc+ tumours 
(12W vs 6W), suggesting reduced transcriptional output from Myc over time. 
However, **Myc mRNA itself shows no significant change** (log2FC = -0.10, 
padj = 0.76), implying:

- Post-transcriptional regulation (protein stability, localisation, modification)
- Cofactor limitation (Max, Miz1, etc.)
- Chromatin accessibility changes
- Negative feedback from Myc targets

#### Limitations

These results identify pathways with differential temporal dynamics between 
genotypes but cannot determine:

- Absolute pathway activity levels at either timepoint
- Whether "Myc+ specific" effects reflect direct Myc regulation
- Causal relationships between Myc and pathway changes

To address these limitations, additional analyses could include:
- Myc+ vs Myc- fGSEA at each timepoint (requires genotype contrasts)
- Leading edge analysis of MYC signatures
- Expression analysis of Myc cofactors and regulators

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
