│   ├── felsher_integrative_signature.csv
│   └── Mouse.MitoCarta3.0.xls        # MitoCarta3.0 mouse annotation (Broad Institute)
├── scripts/s project investigates the effect of the **Myc oncogene** in early breast tumourigenesis using a mouse model. Myc is selectively and constitutively expressed in breast epithelial cells using the MMTV promoter.

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
| `06_fgsea_cross_sectional.R` | Cross-sectional fGSEA: Myc+ vs Myc- at each timepoint |
| `07_fgsea_xs_visualisation.R` | Cross-sectional fGSEA visualisation and NES correlation plots |
| `08_mitoPPS_analysis.R` | MitoPPS: mitochondrial pathway prioritisation scores (Monzel et al. 2025) |

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
│   └── mitopps_scores.rds             # mitoPPS scores, statistics, pathway annotations
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
│   └── fgsea_cross_sectional/         # Cross-sectional fGSEA visualisations
│       ├── xs_category_summary.pdf    # Category bar chart
│       ├── dotplot_*_6W.pdf           # Dot plots at 6W
│       ├── dotplot_*_12W.pdf          # Dot plots at 12W
│       ├── xs_nes_correlation.pdf     # Plot 1: NES correlation (maintenance vs baseline)
│       ├── xs_dev_contribution.pdf    # Plot 2: Developmental contribution
│       ├── xs_nes_temporal_comparison.pdf     # Plot 3: NES Myc- vs Myc+ temporal
│       └── xs_nes_crosssectional_comparison.pdf  # Plot 4: NES 6W vs 12W cross-sectional
│   └── mitopps/                           # mitoPPS visualisations
│       ├── pca_raw_pathway_scores.pdf     # PCA of raw MitoPathway scores
│       ├── pca_mitopps.pdf                # PCA of normalised mitoPPS
│       ├── heatmap_mitopps.pdf            # Group mean heatmap (row z-scored)
│       ├── heatmap_raw_pathway_scores.pdf # Raw scores heatmap (row z-scored)
│       ├── heatmap_mitopps_col_zscore.pdf # Column z-scored (prioritisation within conditions)
│       ├── heatmap_mitopps_unscaled.pdf   # Unscaled log10(mitoPPS) group means
│       ├── dotplot_mitopps_myc_effect.pdf # Myc-driven reprioritisation at 6W and 12W
│       ├── scatter_mitopps_6W_vs_12W_myc_effect.pdf  # Consistency of Myc effect across timepoints
│       ├── dotplot_mitopps_temporal_effect.pdf         # Temporal reprioritisation in Myc- and Myc+
│       ├── scatter_mitopps_temporal_mycneg_vs_mycpos.pdf  # Correlation of temporal trajectories
│       ├── boxplots_top_mitopps_myc.pdf   # Top 12 Myc-affected pathways
│       └── boxplot_total_mito_expression.pdf  # Total mito expression by group
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

## Cross-Sectional fGSEA Analysis

### Scripts

| Script | Description |
|--------|-------------|
| `06_fgsea_cross_sectional.R` | fGSEA on Myc+ vs Myc- contrasts at 6W and 12W separately |
| `07_fgsea_xs_visualisation.R` | Visualisation of cross-sectional results and combined NES correlation plots |

### Strategy

The cross-sectional analysis complements the temporal analysis by asking a different question:

- **Temporal** (scripts 04–05): How do pathways change *over time* within each genotype?
- **Cross-sectional** (scripts 06–07): How do pathways differ *between genotypes* at each timepoint?

fGSEA is run on two contrasts:
1. **Myc+ vs Myc- at 6W** (early Myc effect)
2. **Myc+ vs Myc- at 12W** (late Myc effect)

Results are classified into categories: Stable (significant at both timepoints), 6W only (lost by 12W), 12W only (gained by 12W), Reversed, or Not significant.

### Output

Results saved to `results/fgsea_xs_results.rds` containing:
- Individual fGSEA results for 6W and 12W cross-sectional contrasts
- Combined comparison table with pathway classification
- Pathway subsets (MYC, MitoCarta, Hallmark)

### NES Correlation Plots

Four plots combine information from the temporal and cross-sectional analyses to characterise the relationship between developmental and Myc-driven pathway changes.

#### Plot 1: NES Correlation — Maintenance of the Myc Effect (`xs_nes_correlation.pdf`)

- **X-axis**: Developmental baseline effect (NES: Myc- 12W vs 6W)
- **Y-axis**: Change in cross-sectional Myc effect over time (NES₁₂W − NES₆W from Myc+ vs Myc-)
- **Size**: -log10(padj) of the most significant comparison across all four contrasts
- **Colour**: Developmental contribution (1/ΔNES where ΔNES = NES_pos − NES_neg from temporal)

Interpretation: Pathways in the upper half *gain* Myc enrichment between 6W and 12W; pathways in the lower half *lose* it. The colour indicates whether the temporal difference between Myc+ and Myc- is large (near zero, grey) or small (saturated colour, meaning development drives a similar trajectory in both genotypes).

#### Plot 2: Developmental Contribution (`xs_dev_contribution.pdf`)

- **X-axis**: Developmental baseline effect (NES: Myc- 12W vs 6W)
- **Y-axis**: Developmental contribution (1/[NES_pos − NES_neg])
- **Size**: -log10(padj)
- **Colour**: Myc+ temporal trajectory (NES: Myc+ 12W vs 6W)

Interpretation: Pathways far from y = 0 have similar temporal trajectories in both genotypes (small ΔNES → large 1/ΔNES), suggesting development rather than Myc drives the change. Pathways near y = 0 have large Myc-specific temporal effects. The colour shows the direction of change in Myc+ cells.

#### Plot 3: Temporal NES — Myc- vs Myc+ (`xs_nes_temporal_comparison.pdf`)

- **X-axis**: NES (Myc- 12W vs 6W) — developmental trajectory
- **Y-axis**: NES (Myc+ 12W vs 6W) — Myc+ trajectory
- **Size**: -log10(padj)
- **Colour**: Pathway category (MitoCarta, MYC Signature, Hallmark)

Interpretation: Points on the diagonal have identical temporal trajectories in both genotypes. The deviation from the diagonal represents the Myc-specific component of temporal change. Key observations:

- Most pathways cluster near the diagonal, confirming that developmental effects dominate temporal changes (consistent with the non-significant interaction term from DESeq2)
- MitoCarta pathways (red) and MYC signatures (blue) tend to fall below the diagonal — both genotypes decline over time, but Myc+ declines more steeply
- Pathways above the diagonal (e.g., estrogen response, androgen response) are enhanced in Myc+ relative to the developmental trend

#### Plot 4: Cross-Sectional NES — 6W vs 12W (`xs_nes_crosssectional_comparison.pdf`)

- **X-axis**: NES (Myc+ vs Myc- at 6W) — early Myc effect
- **Y-axis**: NES (Myc+ vs Myc- at 12W) — late Myc effect
- **Size**: -log10(padj)
- **Colour**: Pathway category (MitoCarta, MYC Signature, Hallmark)

Interpretation: Points on the diagonal have a stable Myc effect across timepoints. Deviations indicate pathways where the Myc effect changes over time. Key observations:

- The majority of pathways are positive at both timepoints, clustering in the upper-right quadrant — Myc+ cells are enriched for these pathways relative to Myc- at both 6W and 12W
- MYC signatures (blue) show slightly higher NES at 12W than at 6W (above the diagonal), despite declining in absolute terms over time (Plot 3) — because Myc- cells decline faster
- Translation, OXPHOS, and mTORC1 signaling are strongly positive at both timepoints, confirming these as robust Myc-driven pathways

### Reconciling Plots 3 and 4: A Key Insight

MYC signature pathways show **negative NES in Plot 3** (both genotypes decline over time) but **positive NES in Plot 4** (Myc+ remains enriched relative to Myc- at both timepoints), with a slight increase at 12W (above the diagonal in Plot 4).

This means: **MYC signatures decline in both genotypes over time, but they decline faster in Myc- than in Myc+.** The Myc transgene does not prevent the developmental decline in MYC pathway activity — it buffers against it. As the Myc- baseline drops further at 12W, the relative enrichment in Myc+ cells actually increases.

This is consistent with the non-significant interaction term: the *absolute* Myc effect (difference-of-differences) is small, but the *relative* Myc effect (Myc+ vs Myc- at each timepoint) is maintained or even slightly enhanced.

---

## Mitochondrial Pathway Prioritisation Score (mitoPPS) Analysis

### Script

| Script | Description |
|--------|-------------|
| `08_mitoPPS_analysis.R` | Compute raw MitoPathway scores and mitoPPS, statistical analysis, visualisation |

### Method

MitoPPS was computed following Monzel et al. (2025) (*bioRxiv* 2025.02.03.635951). For each sample, raw MitoPathway scores were calculated as the mean DESeq2-normalised count of member genes for each MitoCarta3.0 pathway (Sheet 4, mouse). MitoPPS was then derived via a three-step pairwise ratio normalisation: (1) all pairwise pathway ratios computed per sample, (2) each ratio corrected by its global mean across samples, (3) corrected ratios averaged per pathway per sample. This normalisation removes both total mitochondrial content and intrinsic scale differences between pathways, so mitoPPS reflects relative mitochondrial resource allocation — which pathways are selectively prioritised — independent of overall mitochondrial abundance. Values centre around 1.0 (dataset average prioritisation), with >1.0 indicating relative up-prioritisation and <1.0 down-prioritisation.

### mtDNA Gene Handling

Mouse mtDNA-encoded genes (13 protein-coding: mt-Nd1–6, mt-Co1–3, mt-Cytb, mt-Atp6, mt-Atp8) are transcribed at orders-of-magnitude higher levels than nuclear-encoded mitochondrial genes. Leaving them in their canonical MitoCarta3.0 pathways (OXPHOS complexes, mitochondrial central dogma) would dominate and distort those pathway scores. Therefore:

1. All mt-* genes are **removed** from their original MitoCarta3.0 pathways
2. A synthetic pathway **"mtDNA-encoded OXPHOS subunits"** is created containing all detected mt-* genes

This preserves the interpretability of nuclear-encoded pathway scores while retaining the mtDNA-encoded contribution as its own pathway in the mitoPPS framework.

### Additional Curated Pathways

The parent MitoCarta3.0 "Apoptosis" pathway is retained, and two child pathways are added:

- **Apoptosis-PRO** (25 genes): pro-apoptotic factors (Bax, Bak1, Bad, Bid, Casp3/8/9, Cycs, etc.)
- **Apoptosis-ANTI** (9 genes): anti-apoptotic factors (Bcl2, Bcl2l1, Mcl1, etc.)

### Results

**Gene coverage:** 952 / 1037 MitoCarta genes found in expression data (91.8%), yielding 142 scoreable pathways.

#### Myc Effect

Myc+ samples show significant mitoPPS reprioritisation compared to Myc− at both 6W and 12W timepoints. The effect is highly consistent across timepoints (scatter plot: 6W vs 12W Myc effect), indicating that Myc imposes a stable mitochondrial prioritisation signature rather than a time-dependent one.

| Effect (ANOVA) | Pathways padj < 0.05 | Pathways padj < 0.10 |
|---|---|---|
| myc_status | 34 | 49 |
| timepoint | 7 | 21 |
| interaction | 0 | 0 |

#### Temporal Effect

No pathways reach significance in the Myc− timecourse (6W→12W), but the Myc+ timecourse shows significant temporal reprioritisation. Critically, the two temporal trajectories are strongly correlated across all pathways (Pearson r = 0.75, p < 2.2e-16, slope = 0.76), indicating that the direction of age-associated mitochondrial reprioritisation is largely shared between Myc− and Myc+. The Myc− effect does not reach significance likely due to higher within-group variance rather than a genuinely absent effect. Myc therefore appears to amplify and stabilise an underlying physiological ageing-associated mitochondrial remodelling programme rather than inducing a novel one.

| Pairwise contrast | Pathways padj < 0.05 | Pathways padj < 0.10 |
|---|---|---|
| Myc effect at 12W | 4 | 16 |
| Myc effect at 6W | 0 | 6 |
| Temporal in Myc+ | 0 | 9 |
| Temporal in Myc− | 0 | 0 |

### Heatmap Scaling

Three heatmap variants are produced for mitoPPS, each answering a different question:

| Heatmap | Scaling | Question answered |
|---|---|---|
| `heatmap_mitopps.pdf` | Row z-score | How does each pathway's prioritisation change across conditions? |
| `heatmap_mitopps_col_zscore.pdf` | Column z-score | How are pathways ranked within each condition? (Closest to the mitoPPS concept) |
| `heatmap_mitopps_unscaled.pdf` | None (log10) | Absolute prioritisation — both between-pathway and between-condition differences preserved |

### Output

Results saved to `results/mitopps_scores.rds` containing per-sample raw and mitoPPS scores, group means, ANOVA and pairwise statistics, pathway annotations, and PCA objects. All figures saved to `outputs/mitopps/`.

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
