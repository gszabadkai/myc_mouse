
# MYC Mouse Time-Course Analysis Summary

_Last updated: 2025-10-02_

## Overview

This project investigates MYC-driven transcriptional changes in a time-course mouse model using interaction and group-based DESeq2 designs. It includes:
- Differential expression analysis (shrunken and raw)
- Classification of genes based on temporal MYC effects
- ComplexHeatmap-based visualisation for pathway gene sets

---

## Scripts and Structure

```
scripts/
├── 00_setup_packages.R           # Load all required packages
├── 01_load_data.R                # Load count matrix, coldata, and gene sets (mitocarta, MYC, felsher)
├── 02_deseq_interaction_model.R  # DESeq2 interaction model + shrinkage
├── 03_lfc_classification.R       # Classify LFC trends across timepoints
├── 04_group_comparison.R         # DESeq2 group comparison (12W_pos vs 6W_pos)
├── 05_heatmap_utils.R            # Core ComplexHeatmap function
├── 06_diagnostics_QC.R           # QC: size factors, PCA, clustering
├── 07_heatmap_oxphos_annotated.R # Special annotated heatmap for OXPHOS subunits
├── rebuild_all.R                 # Re-run scripts in order
```

Data & results are saved under:
```
data/         # Input gene sets and full data matrix
results/      # Saved RDS files: DESeq2 objects, LFC tables, gene sets
outputs/      # Heatmaps and plots
```

---

## Gene Sets

- **Mitocarta v3** gene sets parsed from `mitocarta_pathways.csv`
- **MYC signature gene sets** from `.gmx` format (converted to mouse using Ensembl orthologs)
- **Felsher integrative signature** (human → mouse mapped)

All are merged into `gene_sets_list.rds`

---

## Differential Expression

### Interaction Model
- Design: `~ myc_status * timepoint`
- Extracted contrasts:
  - `myc_6W_log2FC`, `myc_12W_log2FC`
  - `timepoint_neg_log2FC`, `timepoint_pos_log2FC`
- Shrinkage: **applied with `ashr`**
- P-values adjusted with **IHW**

### Group Model
- Design: `~ group` (factor with 4 levels)
- Contrast: `12W_pos vs 6W_pos`
- Shrinkage and IHW applied
- Significant group-based changes annotated as: `group_sig_status`

---

## LFC Classification Logic (Revised)

```r
combined_df_annotated <- combined_df_annotated %>%
  mutate(time_effect_category_shrunk = case_when(
    (myc_12W_log2FC_raw - myc_6W_log2FC_raw) < -0.5 & abs(timepoint_pos_log2FC_raw - timepoint_neg_log2FC_raw) > 1 ~ "direct_Myc_reduction",
    (myc_12W_log2FC_raw - myc_6W_log2FC_raw) > 0.5 & abs(timepoint_pos_log2FC_raw - timepoint_neg_log2FC_raw) > 1 ~ "direct_Myc_increase",
    (myc_12W_log2FC_raw - myc_6W_log2FC_raw) < -0.5 & abs(timepoint_pos_log2FC_raw - timepoint_neg_log2FC_raw) < 1 ~ "baseline_driven_reduction",
    (myc_12W_log2FC_raw - myc_6W_log2FC_raw) > 0.5 & abs(timepoint_pos_log2FC_raw - timepoint_neg_log2FC_raw) < 1 ~ "baseline_driven_increase",
    TRUE ~ "no_change"
  ))
```

- Classification aims to disentangle whether change in MYC effect over time is due to:
  - Decreased MYC levels/interactions → `direct_Myc_reduction`
  - Baseline shift in MYC-neg → `baseline_driven_reduction`

Group comparison significance (`group_sig_status`) is overlaid to reinforce interpretation.

---

## Heatmap Generation

- Function: `generate_heatmaps_for_gene_set()`
- Two panels:
  - **Left**: LFC heatmap (columns: `myc_6W`, `myc_12W`, `timepoint_neg`)
  - **Right**: Sample group expression (log2norm, z-scored)
- Heatmaps are annotated by:
  - `time_effect_category`
  - `group_sig_status`

Other features:
- Row name font size auto-scales by number of genes
- Column clustering disabled (fixed order: 6W_neg, 6W_pos, 12W_neg, 12W_pos)
- Row clustering: `ward.D2` on Pearson correlation

---

## Special OXPHOS Heatmap

- Added custom annotation:
  - **Complex (CI–CV)**
  - **mtDNA-encoded genes** (genes with prefix `mt-`)
- Built with extra row annotation (`rowAnnotation()`)
- Saved separately in `scripts/07_heatmap_oxphos_annotated.R`

---

## Git & Documentation

- Project under version control (GitHub)
- `.gitignore` excludes `results/`, `outputs/`, intermediate RDS
- Markdown files:
  - `README.md`: analysis structure & scripts
  - `NOTES.md`: technical notes (e.g., commit shortcuts, workflow logic)

---

## Outstanding / To-do

- 🔁 Revisit LFC classification logic (e.g. thresholds, significance rules)
- 📊 Consider alternative heatmap clustering methods (e.g. Spearman, complete)
- 📄 Create R Markdown summaries for collaborators
- 🔍 Explore PCA, variancePartition, or batch effect analysis further

---

