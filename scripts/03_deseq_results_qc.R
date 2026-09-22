# scripts/03_deseq_results_qc.R
#
# Purpose: Run DESeq2, extract contrasts, and QC with MA plots
#
# Background:
# - QC revealed a batch effect: 6W samples have ~2x more reads than 12W
# - DESeq2 size factors normalize for this, but shrinkage (ashr/apeglm) can
#   behave differently when variance differs between groups
# - We generate BOTH raw and shrunken LFCs to compare and detect artifacts
# - IHW (Independent Hypothesis Weighting) is used for p-value adjustment:
#   it stratifies genes by baseMean and applies different FDR thresholds,
#   which can increase power when effect sizes correlate with expression level

source(here::here("scripts", "00_setup_packages.R"))

# === Load data ===
dds_int <- readRDS(here("results", "dds_int.rds"))

# Create output directory for QC plots
results_qc_dir <- here("outputs", "deseq_qc")
dir.create(results_qc_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PART 1: INTERACTION MODEL
# =============================================================================
# Design: ~ timepoint * myc_status
# Expands to: ~ timepoint + myc_status + timepoint:myc_status
#
# Coefficients:
#   - Intercept: baseline (6W_neg)
#   - timepoint_12W_vs_6W: effect of time in myc_neg samples
#   - myc_status_pos_vs_neg: effect of Myc at 6W
#   - timepoint12W.myc_statuspos: INTERACTION - how Myc effect changes from 6W to 12W

message("Running DESeq2 on interaction model...")
dds_int <- DESeq(dds_int)

# Check available coefficients
message("Available coefficients:")
print(resultsNames(dds_int))

# -----------------------------------------------------------------------------
# Extract contrasts - RAW (unshrunken) with IHW filtering
# -----------------------------------------------------------------------------
# IHW uses baseMean as a covariate to weight hypotheses. Genes with higher
# expression (and thus lower variance in LFC estimates) get more "statistical
# budget" for discovery.

message("Extracting raw (unshrunken) results with IHW...")

# 1. Myc effect at 6W (main effect of myc_status)
myc_6W_raw <- results(dds_int, 
                      name = "myc_status_pos_vs_neg", 
                      filterFun = ihw)

# 2. Timepoint effect in myc_neg (main effect of timepoint)
timepoint_neg_raw <- results(dds_int, 
                              name = "timepoint_12W_vs_6W", 
                              filterFun = ihw)

# 3. Myc effect at 12W (main effect + interaction)
#    This is: myc_status_pos_vs_neg + timepoint12W.myc_statuspos
myc_12W_raw <- results(dds_int, 
                        contrast = list(c("myc_status_pos_vs_neg", 
                                         "timepoint12W.myc_statuspos")), 
                        filterFun = ihw)

# 4. Timepoint effect in myc_pos (main effect + interaction)
#    This is: timepoint_12W_vs_6W + timepoint12W.myc_statuspos
timepoint_pos_raw <- results(dds_int, 
                              contrast = list(c("timepoint_12W_vs_6W", 
                                               "timepoint12W.myc_statuspos")), 
                              filterFun = ihw)

# 5. Interaction term: difference in Myc effect between timepoints
#    Positive = Myc effect is stronger at 12W; Negative = Myc effect weakens
interaction_raw <- results(dds_int, 
                           name = "timepoint12W.myc_statuspos", 
                           filterFun = ihw)

# -----------------------------------------------------------------------------
# Extract contrasts - SHRUNKEN (ashr) with IHW filtering
# -----------------------------------------------------------------------------
# ashr shrinkage is adaptive: it shrinks noisy estimates more than precise ones.
# This is generally preferred for visualization and ranking, but can behave
# differently when variance structure differs (as with our batch effect).

message("Extracting shrunken (ashr) results...")

myc_6W_shrunk <- lfcShrink(dds_int, 
                           coef = "myc_status_pos_vs_neg", 
                           type = "ashr", 
                           res = myc_6W_raw)

timepoint_neg_shrunk <- lfcShrink(dds_int, 
                                   coef = "timepoint_12W_vs_6W", 
                                   type = "ashr", 
                                   res = timepoint_neg_raw)

myc_12W_shrunk <- lfcShrink(dds_int, 
                            contrast = list(c("myc_status_pos_vs_neg", 
                                             "timepoint12W.myc_statuspos")), 
                            type = "ashr", 
                            res = myc_12W_raw)

timepoint_pos_shrunk <- lfcShrink(dds_int, 
                                   contrast = list(c("timepoint_12W_vs_6W", 
                                                    "timepoint12W.myc_statuspos")), 
                                   type = "ashr", 
                                   res = timepoint_pos_raw)

interaction_shrunk <- lfcShrink(dds_int, 
                                coef = "timepoint12W.myc_statuspos", 
                                type = "ashr", 
                                res = interaction_raw)

# =============================================================================
# PART 2: GROUP-BASED MODEL
# =============================================================================
# Design: ~ group (where group = 6W_neg, 6W_pos, 12W_neg, 12W_pos)
#
# This allows direct pairwise comparisons without interaction terms.
# Useful when the interaction model lacks power (noisy data) or when you
# want to directly compare specific groups.

message("Running DESeq2 on group model...")

cts <- readRDS(here("results", "count_matrix.rds"))
coldata <- readRDS(here("results", "coldata.rds"))

dds_group <- DESeqDataSetFromMatrix(countData = cts,
                                    colData = coldata,
                                    design = ~ group)

# Apply same filtering as interaction model
keep <- rowSums(counts(dds_group) >= 10) >= 4
dds_group <- dds_group[keep, ]
dds_group <- DESeq(dds_group)

# Key comparison: 12W_pos vs 6W_pos (progressive Myc effect)
# This captures genes that change in Myc+ samples over time
group_12Wpos_vs_6Wpos_raw <- results(dds_group, 
                                     contrast = c("group", "12W_pos", "6W_pos"), 
                                     filterFun = ihw)

group_12Wpos_vs_6Wpos_shrunk <- lfcShrink(dds_group, 
                                          contrast = c("group", "12W_pos", "6W_pos"), 
                                          type = "ashr", 
                                          res = group_12Wpos_vs_6Wpos_raw)

# Additional comparisons for completeness
group_12Wneg_vs_6Wneg_raw <- results(dds_group, 
                                     contrast = c("group", "12W_neg", "6W_neg"), 
                                     filterFun = ihw)

group_12Wneg_vs_6Wneg_shrunk <- lfcShrink(dds_group, 
                                          contrast = c("group", "12W_neg", "6W_neg"), 
                                          type = "ashr", 
                                          res = group_12Wneg_vs_6Wneg_raw)

# =============================================================================
# PART 3: MA PLOT QC
# =============================================================================
# MA plots show log2FC vs mean expression (baseMean).
# Key things to look for:
# - Symmetry around y=0 (no global bias)
# - Spread should be wider at low baseMean (expected)
# - Shrinkage should reduce spread, especially at low expression
# - Compare raw vs shrunk to see if shrinkage is appropriate

message("Generating MA plots for QC...")

# Helper function to generate paired MA plots (raw vs shrunk)
plot_ma_comparison <- function(res_raw, res_shrunk, title_prefix, filename) {
  pdf(file.path(results_qc_dir, filename), width = 10, height = 5)
  par(mfrow = c(1, 2))
  
  # Raw
  DESeq2::plotMA(res_raw, ylim = c(-4, 4), main = paste(title_prefix, "- Raw"))
  abline(h = c(-1, 1), col = "dodgerblue", lty = 2)
  

  # Shrunk
  DESeq2::plotMA(res_shrunk, ylim = c(-4, 4), main = paste(title_prefix, "- Shrunk (ashr)"))
  abline(h = c(-1, 1), col = "dodgerblue", lty = 2)
  
  dev.off()
}

# Generate MA plots for all contrasts
plot_ma_comparison(myc_6W_raw, myc_6W_shrunk, 
                   "Myc effect at 6W", "MA_myc_6W.pdf")

plot_ma_comparison(myc_12W_raw, myc_12W_shrunk, 
                   "Myc effect at 12W", "MA_myc_12W.pdf")

plot_ma_comparison(timepoint_neg_raw, timepoint_neg_shrunk, 
                   "Timepoint effect (Myc-)", "MA_timepoint_neg.pdf")

plot_ma_comparison(timepoint_pos_raw, timepoint_pos_shrunk, 
                   "Timepoint effect (Myc+)", "MA_timepoint_pos.pdf")

plot_ma_comparison(interaction_raw, interaction_shrunk, 
                   "Interaction (Myc x Time)", "MA_interaction.pdf")

plot_ma_comparison(group_12Wpos_vs_6Wpos_raw, group_12Wpos_vs_6Wpos_shrunk, 
                   "Group: 12W_pos vs 6W_pos", "MA_group_12Wpos_vs_6Wpos.pdf")

plot_ma_comparison(group_12Wneg_vs_6Wneg_raw, group_12Wneg_vs_6Wneg_shrunk, 
                   "Group: 12W_neg vs 6W_neg", "MA_group_12Wneg_vs_6Wneg.pdf")

# =============================================================================
# PART 4: SUMMARY STATISTICS
# =============================================================================

message("Generating summary statistics...")

summarize_results <- function(res, name) {
  data.frame(
    contrast = name,
    total_genes = nrow(res),
    sig_01 = sum(res$padj < 0.1, na.rm = TRUE),
    sig_005 = sum(res$padj < 0.05, na.rm = TRUE),
    up_01 = sum(res$padj < 0.1 & res$log2FoldChange > 0, na.rm = TRUE),
    down_01 = sum(res$padj < 0.1 & res$log2FoldChange < 0, na.rm = TRUE)
  )
}

summary_table <- bind_rows(
  summarize_results(myc_6W_raw, "myc_6W"),
  summarize_results(myc_12W_raw, "myc_12W"),
  summarize_results(timepoint_neg_raw, "timepoint_neg"),
  summarize_results(timepoint_pos_raw, "timepoint_pos"),
  summarize_results(interaction_raw, "interaction"),
  summarize_results(group_12Wpos_vs_6Wpos_raw, "group_12Wpos_vs_6Wpos"),
  summarize_results(group_12Wneg_vs_6Wneg_raw, "group_12Wneg_vs_6Wneg")
)

write.csv(summary_table, file.path(results_qc_dir, "results_summary.csv"), row.names = FALSE)

message("Summary of significant genes (padj < 0.1):")
print(summary_table)

# =============================================================================
# PART 5: SAVE RESULTS
# =============================================================================

message("Saving results...")

# Save DESeq2 objects
saveRDS(dds_int, here("results", "dds_int_run.rds"))
saveRDS(dds_group, here("results", "dds_group_run.rds"))

# Save all results in a list for easy access
interaction_results <- list(
  myc_6W_raw = myc_6W_raw,
  myc_6W_shrunk = myc_6W_shrunk,
  myc_12W_raw = myc_12W_raw,
  myc_12W_shrunk = myc_12W_shrunk,
  timepoint_neg_raw = timepoint_neg_raw,
  timepoint_neg_shrunk = timepoint_neg_shrunk,
  timepoint_pos_raw = timepoint_pos_raw,
  timepoint_pos_shrunk = timepoint_pos_shrunk,
  interaction_raw = interaction_raw,
  interaction_shrunk = interaction_shrunk
)

group_results <- list(
  pos_12W_vs_6W_raw = group_12Wpos_vs_6Wpos_raw,
  pos_12W_vs_6W_shrunk = group_12Wpos_vs_6Wpos_shrunk,
  neg_12W_vs_6W_raw = group_12Wneg_vs_6Wneg_raw,
  neg_12W_vs_6W_shrunk = group_12Wneg_vs_6Wneg_shrunk
)

saveRDS(interaction_results, here("results", "interaction_results.rds"))
saveRDS(group_results, here("results", "group_results.rds"))

message("Done! QC plots saved to: ", results_qc_dir)
message("Results saved to: results/")

# =============================================================================
# PART 6: EXTENDED QC - UNDERSTANDING THE INTERACTION TERM
# =============================================================================
# 
# Key puzzle: No significant interactions, but ~2,600 genes differ between
# 12W_pos and 6W_pos. This section investigates why.
#
# Note on batch effects and normalization:
# DESeq2's size factors correctly account for library size differences 
# (~2x between 6W and 12W). The group comparisons are NOT confounded by this.
# The timepoint effects capture true biological differences (developmental,
# tumour progression) plus any residual technical variation beyond library size.
# =============================================================================

message("\n", strrep("=", 70))
message("EXTENDED QC: INTERACTION TERM ANALYSIS")
message(strrep("=", 70))

# -----------------------------------------------------------------------------
# 6.1 Extract significant gene sets from each comparison
# -----------------------------------------------------------------------------

sig_timepoint_pos <- interaction_results$timepoint_pos_raw |>
  as.data.frame() |>
  filter(padj < 0.1) |>
  rownames()

sig_timepoint_neg <- interaction_results$timepoint_neg_raw |>
  as.data.frame() |>
  filter(padj < 0.1) |>
  rownames()

sig_myc_6W <- interaction_results$myc_6W_raw |>
  as.data.frame() |>
  filter(padj < 0.1) |>
  rownames()

sig_myc_12W <- interaction_results$myc_12W_raw |>
  as.data.frame() |>
  filter(padj < 0.1) |>
  rownames()

# -----------------------------------------------------------------------------
# 6.2 Overlap analysis: timepoint effects in Myc+ vs Myc- samples
# -----------------------------------------------------------------------------

overlap_timepoint <- intersect(sig_timepoint_pos, sig_timepoint_neg)

timepoint_overlap_summary <- list(
  timepoint_pos_sig = length(sig_timepoint_pos),
  timepoint_neg_sig = length(sig_timepoint_neg),
  overlap = length(overlap_timepoint),
  unique_to_pos = length(setdiff(sig_timepoint_pos, sig_timepoint_neg)),
  unique_to_neg = length(setdiff(sig_timepoint_neg, sig_timepoint_pos)),
  percent_overlap_of_neg = round(100 * length(overlap_timepoint) / length(sig_timepoint_neg), 1)
)

message("\nTimepoint effect overlap (Myc+ vs Myc-):")
message(sprintf("  Significant in Myc+: %d", timepoint_overlap_summary$timepoint_pos_sig))
message(sprintf("  Significant in Myc-: %d", timepoint_overlap_summary$timepoint_neg_sig))
message(sprintf("  Overlap: %d (%.1f%% of Myc-)", 
                timepoint_overlap_summary$overlap, 
                timepoint_overlap_summary$percent_overlap_of_neg))
message(sprintf("  Unique to Myc+: %d", timepoint_overlap_summary$unique_to_pos))
message(sprintf("  Unique to Myc-: %d", timepoint_overlap_summary$unique_to_neg))

# -----------------------------------------------------------------------------
# 6.3 LFC correlation analysis
# -----------------------------------------------------------------------------
# If genes change similarly in both Myc+ and Myc-, the interaction would be ~0

all_timepoint_genes <- union(sig_timepoint_pos, sig_timepoint_neg)

lfc_comparison <- data.frame(
  gene = all_timepoint_genes,
  lfc_pos = interaction_results$timepoint_pos_raw[all_timepoint_genes, "log2FoldChange"],
  lfc_neg = interaction_results$timepoint_neg_raw[all_timepoint_genes, "log2FoldChange"],
  sig_in_pos = all_timepoint_genes %in% sig_timepoint_pos,
  sig_in_neg = all_timepoint_genes %in% sig_timepoint_neg,
  sig_in_both = all_timepoint_genes %in% overlap_timepoint
)

# Correlations
cor_all_sig <- cor(lfc_comparison$lfc_pos, lfc_comparison$lfc_neg, use = "complete.obs")
cor_overlap_only <- cor(
  lfc_comparison$lfc_pos[lfc_comparison$sig_in_both], 
  lfc_comparison$lfc_neg[lfc_comparison$sig_in_both], 
  use = "complete.obs"
)

message("\nLFC correlation (timepoint effects):")
message(sprintf("  All sig genes (n=%d): r = %.3f", nrow(lfc_comparison), cor_all_sig))
message(sprintf("  Overlap genes (n=%d): r = %.3f", sum(lfc_comparison$sig_in_both), cor_overlap_only))

# -----------------------------------------------------------------------------
# 6.4 Interaction p-value analysis for "unique to Myc+" genes
# -----------------------------------------------------------------------------
# These genes SHOULD show interaction if the Myc effect truly differs

unique_to_pos <- setdiff(sig_timepoint_pos, sig_timepoint_neg)

interaction_for_unique <- interaction_results$interaction_raw[unique_to_pos, ] |>
  as.data.frame() |>
  mutate(gene = unique_to_pos)

# Expected number of p < 0.05 by chance
expected_by_chance <- length(unique_to_pos) * 0.05

interaction_unique_summary <- list(
  n_unique_to_pos = length(unique_to_pos),
  pval_below_0.05 = sum(interaction_for_unique$pvalue < 0.05, na.rm = TRUE),
  pval_below_0.01 = sum(interaction_for_unique$pvalue < 0.01, na.rm = TRUE),
  expected_by_chance_0.05 = round(expected_by_chance),
  padj_below_0.1 = sum(interaction_for_unique$padj < 0.1, na.rm = TRUE),
  padj_below_0.2 = sum(interaction_for_unique$padj < 0.2, na.rm = TRUE)
)

message("\nInteraction p-values for genes unique to Myc+ timepoint effect:")
message(sprintf("  Total genes: %d", interaction_unique_summary$n_unique_to_pos))
message(sprintf("  p < 0.05: %d (expected by chance: %d)", 
                interaction_unique_summary$pval_below_0.05,
                interaction_unique_summary$expected_by_chance_0.05))
message(sprintf("  p < 0.01: %d", interaction_unique_summary$pval_below_0.01))
message(sprintf("  padj < 0.1: %d", interaction_unique_summary$padj_below_0.1))
message(sprintf("  padj < 0.2: %d", interaction_unique_summary$padj_below_0.2))

# -----------------------------------------------------------------------------
# 6.5 Overall interaction term statistics
# -----------------------------------------------------------------------------

interaction_all <- interaction_results$interaction_raw |> as.data.frame()

interaction_overall_summary <- list(
  min_padj = min(interaction_all$padj, na.rm = TRUE),
  n_padj_below_0.5 = sum(interaction_all$padj < 0.5, na.rm = TRUE),
  n_padj_below_0.2 = sum(interaction_all$padj < 0.2, na.rm = TRUE),
  n_pval_below_0.05 = sum(interaction_all$pvalue < 0.05, na.rm = TRUE),
  n_pval_below_0.01 = sum(interaction_all$pvalue < 0.01, na.rm = TRUE)
)

message("\nOverall interaction term statistics:")
message(sprintf("  Minimum padj: %.3f", interaction_overall_summary$min_padj))
message(sprintf("  Genes with padj < 0.2: %d", interaction_overall_summary$n_padj_below_0.2))
message(sprintf("  Genes with p < 0.05: %d (expected: %d)", 
                interaction_overall_summary$n_pval_below_0.05,
                round(nrow(interaction_all) * 0.05)))

# -----------------------------------------------------------------------------
# 6.6 Generate extended QC plots
# -----------------------------------------------------------------------------

# Plot 1: LFC correlation between timepoint effects
p_lfc_correlation <- ggplot(lfc_comparison, aes(x = lfc_neg, y = lfc_pos)) +
  geom_point(aes(color = sig_in_both), alpha = 0.3, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.8) +
  scale_color_manual(
    values = c("TRUE" = "darkgreen", "FALSE" = "grey50"),
    labels = c("TRUE" = "Sig in both", "FALSE" = "Sig in one only"),
    name = "Significance"
  ) +
  labs(
    title = "Timepoint effect: Myc+ vs Myc- samples",
    subtitle = sprintf("Genes sig in either comparison (n=%d), r = %.2f", 
                       nrow(lfc_comparison), cor_all_sig),
    x = "log2FC (12W vs 6W) in Myc- samples",
    y = "log2FC (12W vs 6W) in Myc+ samples"
  ) +
  coord_fixed(xlim = c(-6, 6), ylim = c(-6, 6)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

# Plot 2: Interaction p-value histogram (overall)
p_interaction_pval_hist <- ggplot(interaction_all, aes(x = pvalue)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", boundary = 0) +
  geom_hline(yintercept = nrow(interaction_all) / 50, color = "red", linetype = "dashed") +
  labs(
    title = "Interaction term p-value distribution",
    subtitle = "Red line = expected under null (uniform)",
    x = "p-value",
    y = "Count"
  ) +
  theme_bw(base_size = 10)

# Plot 3: Interaction p-values for genes unique to Myc+ timepoint
p_interaction_unique_hist <- ggplot(interaction_for_unique, aes(x = pvalue)) +
  geom_histogram(bins = 30, fill = "darkorange", color = "white", boundary = 0) +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
  labs(
    title = "Interaction p-values: genes unique to Myc+ timepoint effect",
    subtitle = sprintf("n = %d genes; %d have p < 0.05 (expected: %d)", 
                       length(unique_to_pos),
                       sum(interaction_for_unique$pvalue < 0.05, na.rm = TRUE),
                       round(length(unique_to_pos) * 0.05)),
    x = "Interaction p-value",
    y = "Count"
  ) +
  theme_bw(base_size = 10)

# Save extended QC plots
pdf(here(results_qc_dir, "extended_qc_interaction_analysis.pdf"), width = 10, height = 10)

# Page 1: LFC correlation
print(p_lfc_correlation)

# Page 2: P-value histograms
gridExtra::grid.arrange(
  p_interaction_pval_hist,
  p_interaction_unique_hist,
  ncol = 1
)

dev.off()

message(sprintf("\nSaved: %s", here(results_qc_dir, "extended_qc_interaction_analysis.pdf")))

# -----------------------------------------------------------------------------
# 6.7 Save extended QC summaries
# -----------------------------------------------------------------------------

extended_qc_summary <- list(
  timepoint_overlap = timepoint_overlap_summary,
  lfc_correlation = list(
    all_sig_genes = cor_all_sig,
    overlap_genes = cor_overlap_only,
    n_all = nrow(lfc_comparison),
    n_overlap = sum(lfc_comparison$sig_in_both)
  ),
  interaction_unique_to_pos = interaction_unique_summary,
  interaction_overall = interaction_overall_summary
)

saveRDS(extended_qc_summary, here("results", "extended_qc_summary.rds"))
saveRDS(lfc_comparison, here("results", "lfc_comparison_timepoint.rds"))

message(sprintf("Saved: %s", here("results", "extended_qc_summary.rds")))
message(sprintf("Saved: %s", here("results", "lfc_comparison_timepoint.rds")))

# =============================================================================
# SUMMARY: INTERPRETATION
# =============================================================================
message("\n", strrep("=", 70))
message("INTERPRETATION SUMMARY")
message(strrep("=", 70))
message("
1. INTERACTION TERM IS UNDERPOWERED
   - No genes reach padj < 0.1 (minimum padj ~ 0.16)
   - P-value histogram is near-uniform: no enrichment of true signal
   - The ~357 genes with p < 0.05 among 'unique to Myc+' exceed chance (~83),
     suggesting weak signal, but insufficient for FDR correction

2. TIMEPOINT EFFECTS ARE HIGHLY CORRELATED
   - Genes changing over time do so similarly in Myc+ and Myc- (r = 0.74-0.96)
   - This explains why interaction ~ 0: the difference-of-differences is small
   - The ~961 overlapping genes (51% of Myc- sig) change in the same direction

3. BIOLOGICAL CONCLUSION
   - The Myc transcriptional program is largely STABLE between 6W and 12W
   - Timepoint effects are predominantly developmental/tissue changes
   - Genes 'unique' to Myc+ timepoint are likely borderline cases, not true
     Myc-specific temporal dynamics

4. PRACTICAL RECOMMENDATION
   - Use GROUP MODEL for clear pairwise comparisons
   - Reserve interaction term for hypothesis-driven checks on specific genes
   - Consider relaxed thresholds (padj < 0.2) for exploratory interaction analysis
")

message("\n", strrep("=", 70))
message("DESeq2 RESULTS QC COMPLETE")
message(strrep("=", 70))
