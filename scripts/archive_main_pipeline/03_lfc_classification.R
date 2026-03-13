# scripts/03_lfc_classification.R

source("scripts/00_setup_packages.R")

# === Load annotated LFC tables ===
combined_df_annotated <- readRDS("results/combined_df_annotated.rds")
dds_int <- readRDS("results/dds_int_run.rds")

# # === Get interaction results (raw and IHW-filtered) ===
# interaction_res_raw <- results(dds_int, name = "timepoint12W.myc_statuspos", filterFun = ihw) %>%
#   as.data.frame() %>%
#   rownames_to_column("gene") %>%
#   dplyr::select(gene, interaction_log2FC_raw = log2FoldChange, interaction_padj_raw = padj)
# 
# # === Shrunk interaction result (ashr + IHW) ===
# interaction_res_shrunk <- lfcShrink(dds_int, coef = "timepoint12W.myc_statuspos", type = "ashr", res = results(dds_int, name = "timepoint12W.myc_statuspos", filterFun = ihw)) %>%
#   as.data.frame() %>%
#   rownames_to_column("gene") %>%
#   dplyr::select(gene, interaction_log2FC_shrunk = log2FoldChange, interaction_padj_shrunk = padj)
# 
# # === Merge with LFC tables ===
# combined_df_annotated <- combined_df_annotated %>%
#   left_join(interaction_res_shrunk, by = "gene")
# 
# combined_df_annotated_raw <- combined_df_annotated_raw %>%
#   left_join(interaction_res_raw, by = "gene")

# # === Classification logic ===
# # there are no significant interaction terms, so this classification indicates trends behind changes in log2FC at 6W and 12W. Combined with group based significance it will indicate the number of directly reduced or increased myc affected genes.
# combined_df_annotated <- combined_df_annotated %>%
#   mutate(time_effect_category_shrunk = case_when(
#     (myc_12W_log2FC - myc_6W_log2FC) < -0.3 & abs(timepoint_pos_log2FC - timepoint_neg_log2FC) > 0.4 ~ "direct_Myc_reduction",
#     (myc_12W_log2FC - myc_6W_log2FC) > 0.3 & abs(timepoint_pos_log2FC - timepoint_neg_log2FC) > 0.4 ~ "direct_Myc_increase",
#     (myc_12W_log2FC - myc_6W_log2FC) < -0.3 & abs(timepoint_pos_log2FC - timepoint_neg_log2FC) < 0.4 ~ "baseline_driven_reduction",
#     (myc_12W_log2FC - myc_6W_log2FC) > 0.3 & abs(timepoint_pos_log2FC - timepoint_neg_log2FC) < 0.4 ~ "baseline_driven_increase",
#     TRUE ~ "no_change"
#   ))
# 
# combined_df_annotated_raw <- combined_df_annotated_raw %>%
#   mutate(time_effect_category_raw = case_when(
#     (myc_12W_log2FC_raw - myc_6W_log2FC_raw) < -0.3 & abs(timepoint_pos_log2FC_raw - timepoint_neg_log2FC_raw) > 0.4 ~ "direct_Myc_reduction",
#     (myc_12W_log2FC_raw - myc_6W_log2FC_raw) > 0.3 & abs(timepoint_pos_log2FC_raw - timepoint_neg_log2FC_raw) > 0.4 ~ "direct_Myc_increase",
#     (myc_12W_log2FC_raw - myc_6W_log2FC_raw) < -0.3 & abs(timepoint_pos_log2FC_raw - timepoint_neg_log2FC_raw) < 0.4 ~ "baseline_driven_reduction",
#     (myc_12W_log2FC_raw - myc_6W_log2FC_raw) > 0.3 & abs(timepoint_pos_log2FC_raw - timepoint_neg_log2FC_raw) < 0.4 ~ "baseline_driven_increase",
#     TRUE ~ "no_change"
#   ))

alpha <- 0.10
trend_z_buffer <- 1.00   # how much smaller |stat_12W| must be than |stat_6W| to call "suspected attenuation"

combined_df_annotated <- combined_df_annotated %>%
  mutate(
    # Evidence decisions (minimum-effect tests)
    mu6_GE   = !is.na(mu6_padj)   & mu6_padj   < alpha,
    mu12_GE  = !is.na(mu12_padj)  & mu12_padj  < alpha,
    tau_GE   = !is.na(tau_padj)   & tau_padj   < alpha,
    delta_GE = !is.na(delta_padj) & delta_padj < alpha,
    
    class_evidence = dplyr::case_when(
      tau_GE & !mu6_GE & !mu12_GE ~ "baseline_driven",
      !tau_GE & (mu6_GE | mu12_GE) ~ "direct_MYC_effect_stable",
      tau_GE & (mu6_GE | mu12_GE)  ~ "mixed_baseline_and_MYC",
      TRUE ~ "ambiguous"
    ),
    
    # Trend flags (use raw signs + Wald stats to avoid differential shrinkage)
    same_sign_raw = !is.na(mu6_stat) & !is.na(mu12_stat) &
      sign(mu6_stat) == sign(mu12_stat),
    suspected_attenuation = same_sign_raw &
      !is.na(mu6_stat) & !is.na(mu12_stat) &
      (abs(mu12_stat) + 1e-9) < (abs(mu6_stat) - trend_z_buffer),
    
    suspected_gain = same_sign_raw &
      !is.na(mu6_stat) & !is.na(mu12_stat) &
      (abs(mu12_stat) - 1e-9) > (abs(mu6_stat) + trend_z_buffer),
    
    class_display = dplyr::case_when(
      class_evidence == "mixed_baseline_and_MYC"  & suspected_attenuation ~ "mixed + suspected_MYC_attenuation",
      class_evidence == "direct_MYC_effect_stable" & suspected_attenuation ~ "direct_MYC (suspected attenuation)",
      TRUE ~ class_evidence
    )
  )

# === Save updated annotated tables ===
saveRDS(combined_df_annotated,     "results/combined_df_annotated.rds")

# QC (optional)
write.csv(as.data.frame(table(combined_df_annotated$class_evidence)), "results/class_counts_evidence.csv", row.names=FALSE)
write.csv(as.data.frame(table(combined_df_annotated$class_display)),  "results/class_counts_display.csv",  row.names=FALSE)


