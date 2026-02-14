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
