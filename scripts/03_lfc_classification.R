# scripts/03_lfc_classification.R

source("scripts/00_setup_packages.R")

# === Load annotated LFC tables ===
combined_df_annotated <- readRDS("results/combined_df_annotated.rds")
combined_df_annotated_raw <- readRDS("results/combined_df_annotated_raw.rds")
dds_int <- readRDS("results/dds_int_run.rds")

# === Get interaction results (raw and IHW-filtered) ===
interaction_res_raw <- results(dds_int, name = "timepoint12W.myc_statuspos", filterFun = ihw) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, interaction_log2FC_raw = log2FoldChange, interaction_padj_raw = padj)

# === Shrunk interaction result (ashr + IHW) ===
interaction_res_shrunk <- lfcShrink(dds_int, coef = "timepoint12W.myc_statuspos", type = "ashr", res = results(dds_int, name = "timepoint12W.myc_statuspos", filterFun = ihw)) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, interaction_log2FC_shrunk = log2FoldChange, interaction_padj_shrunk = padj)

# === Merge with LFC tables ===
combined_df_annotated <- combined_df_annotated %>%
  left_join(interaction_res_shrunk, by = "gene")

combined_df_annotated_raw <- combined_df_annotated_raw %>%
  left_join(interaction_res_raw, by = "gene")

# === Classification logic ===
# there are no significant interaction terms, so this classification indicates trends behind changes in log2FC at 6W and 12W. Combined with group based significance it will indicate the number of directly reduced or increased myc affected genes.
combined_df_annotated <- combined_df_annotated %>%
  mutate(time_effect_category_shrunk = case_when(
    (myc_12W_log2FC - myc_6W_log2FC) < -0.3 & abs(timepoint_pos_log2FC - timepoint_neg_log2FC) > 0.4 ~ "direct_Myc_reduction",
    (myc_12W_log2FC - myc_6W_log2FC) > 0.3 & abs(timepoint_pos_log2FC - timepoint_neg_log2FC) > 0.4 ~ "direct_Myc_increase",
    (myc_12W_log2FC - myc_6W_log2FC) < -0.3 & abs(timepoint_pos_log2FC - timepoint_neg_log2FC) < 0.4 ~ "baseline_driven_reduction",
    (myc_12W_log2FC - myc_6W_log2FC) > 0.3 & abs(timepoint_pos_log2FC - timepoint_neg_log2FC) < 0.4 ~ "baseline_driven_increase",
    TRUE ~ "no_change"
  ))

combined_df_annotated_raw <- combined_df_annotated_raw %>%
  mutate(time_effect_category_raw = case_when(
    (myc_12W_log2FC_raw - myc_6W_log2FC_raw) < -0.3 & abs(timepoint_pos_log2FC_raw - timepoint_neg_log2FC_raw) > 0.4 ~ "direct_Myc_reduction",
    (myc_12W_log2FC_raw - myc_6W_log2FC_raw) > 0.3 & abs(timepoint_pos_log2FC_raw - timepoint_neg_log2FC_raw) > 0.4 ~ "direct_Myc_increase",
    (myc_12W_log2FC_raw - myc_6W_log2FC_raw) < -0.3 & abs(timepoint_pos_log2FC_raw - timepoint_neg_log2FC_raw) < 0.4 ~ "baseline_driven_reduction",
    (myc_12W_log2FC_raw - myc_6W_log2FC_raw) > 0.3 & abs(timepoint_pos_log2FC_raw - timepoint_neg_log2FC_raw) < 0.4 ~ "baseline_driven_increase",
    TRUE ~ "no_change"
  ))

# === Save updated annotated tables ===
saveRDS(combined_df_annotated, "results/combined_df_annotated.rds")
saveRDS(combined_df_annotated_raw, "results/combined_df_annotated_raw.rds")
