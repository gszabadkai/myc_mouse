# scripts/14_gate2_apoptosis_readout.R
# =============================================================================
# GATE 2 (Block A decision layer): Apoptosis readout -- decision point vs readout
# =============================================================================
#
# Purpose (Gate 2 in docs/2026-07-04_block_A_build_spec.md):
#   Decide whether the mitochondrion acts as a DECISION POINT (principle 5:
#   Myc's fading apoptotic drive reflects a coordinated, module-wide shift in
#   the pro-apoptotic programme -> selection / active decision) or merely a
#   READOUT (isolated / absent group-level shift -> the apoptosis signal is a
#   passive by-product; lean on the metabolic-readout interpretation).
#
#   The test is selection vs remodelling on the interaction LFC:
#     - SELECTION  predicts a COORDINATED shift across the whole pro-apoptotic
#       module (Apoptosis-PRO mean interaction LFC != 0 as a group).
#     - REMODELLING predicts ISOLATED gene-level changes (Bbc3/PUMA an outlier;
#       the PRO group as a whole not shifted).
#
#   Gate logic:
#     - If the group-level shift is NULL (current padj ~0.91 on the interaction
#       term for these genes tilts readout/metabolic), AP4 + the apoptosis
#       framing shrink; cell death is held at Supp 1.
#     - If the PRO module IS coordinately shifted, the decision-point framing
#       (selection) survives and the apoptosis arm stays in the main figures.
#
#   Sign convention (matches script 11):
#     interaction log2FC negative = Myc effect WEAKENS at 12W.
#
# CRITICAL (tension T3, resolved): the PRO/ANTI stats read here were computed by
#   script 11 on interaction_raw_df = interaction_results$interaction_raw, i.e.
#   the UNSHRUNKEN MLE fit. IHW reweights only padj, never log2FoldChange/stat.
#   So Gate 2 already runs on raw/Wald values -- no shrinkage recompute needed.
#
# WRINKLE handled here: script 11 persists the per-gene apoptosis_pro /
#   apoptosis_anti tables (with log2FoldChange) but NOT its one-sample t-tests or
#   the Bbc3 outlier stats. This script RE-DERIVES those from the saved LFC
#   columns, reproducing script 11 exactly (t.test on the same values).
#
# This is a REFRAME script: reads existing outputs, re-derives the group tests,
# collates the gate decision. No model fitting.
#
# Input:
#   - results/interaction_gene_characterisation.rds (script 11)
#       $interaction_by_geneset$apoptosis_pro   (25 genes: gene_symbol,
#                                                 log2FoldChange, lfcSE, pvalue,
#                                                 stat, sig, direction)
#       $interaction_by_geneset$apoptosis_anti  (9 genes, same columns)
#
# Output:
#   - results/gate2_apoptosis_readout.rds
#   - outputs/gates/gate2_summary.csv
#   - outputs/gates/gate2_apoptosis_forest.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

gate_dir <- here("outputs", "gates")
dir.create(gate_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PART 1: LOAD PERSISTED PRO / ANTI INTERACTION STATS (raw / Wald-based)
# =============================================================================

message("Loading interaction gene characterisation (script 11 outputs)...")

characterisation <- readRDS(here("results", "interaction_gene_characterisation.rds"))

stats_pro  <- characterisation$interaction_by_geneset$apoptosis_pro
stats_anti <- characterisation$interaction_by_geneset$apoptosis_anti

stopifnot(all(c("gene_symbol", "log2FoldChange", "lfcSE", "pvalue") %in% names(stats_pro)))
stopifnot(all(c("gene_symbol", "log2FoldChange", "lfcSE", "pvalue") %in% names(stats_anti)))

message(sprintf("  Apoptosis-PRO:  %d genes (%d with interaction p<0.05)",
                nrow(stats_pro), sum(stats_pro$pvalue < 0.05, na.rm = TRUE)))
message(sprintf("  Apoptosis-ANTI: %d genes (%d with interaction p<0.05)",
                nrow(stats_anti), sum(stats_anti$pvalue < 0.05, na.rm = TRUE)))

# =============================================================================
# PART 2: RE-DERIVE THE GROUP-LEVEL SHIFT TESTS (reproduce script 11)
# =============================================================================
# One-sample t-tests: is each module's mean interaction LFC != 0?
#   - PRO shifted negative  = pro-apoptotic Myc drive coordinately fades at 12W
#   - ANTI shifted positive = anti-apoptotic Myc drive coordinately rises at 12W
# Either coordinated shift is evidence for a module-wide (selection-like) change.

message("\n", strrep("=", 70))
message("PART 2: GROUP-LEVEL SHIFT (one-sample t-tests on interaction LFC)")
message(strrep("=", 70))

pro_ttest  <- t.test(stats_pro$log2FoldChange,  mu = 0)
anti_ttest <- t.test(stats_anti$log2FoldChange, mu = 0)

pro_median  <- median(stats_pro$log2FoldChange,  na.rm = TRUE)
anti_median <- median(stats_anti$log2FoldChange, na.rm = TRUE)

message(sprintf("  Apoptosis-PRO : mean LFC = %+.4f, median = %+.4f, one-sample t p = %.4f",
                pro_ttest$estimate, pro_median, pro_ttest$p.value))
message(sprintf("  Apoptosis-ANTI: mean LFC = %+.4f, median = %+.4f, one-sample t p = %.4f",
                anti_ttest$estimate, anti_median, anti_ttest$p.value))

# =============================================================================
# PART 3: Bbc3 / PUMA OUTLIER STATUS
# =============================================================================
# If Bbc3 carries the apoptosis signal alone (extreme vs the PRO median), the
# fading is an isolated remodelling event, not a module-wide shift.

message("\n", strrep("=", 70))
message("PART 3: Bbc3 OUTLIER STATUS")
message(strrep("=", 70))

bbc3_row <- stats_pro |> dplyr::filter(gene_symbol == "Bbc3")
if (nrow(bbc3_row) == 1) {
  bbc3_lfc  <- bbc3_row$log2FoldChange
  bbc3_rank <- which(sort(stats_pro$log2FoldChange) == bbc3_lfc)[1]
  bbc3_is_outlier <- abs(bbc3_lfc) > abs(pro_median) * 2
  message(sprintf("  Bbc3 interaction LFC = %+.4f (rank %d/%d in PRO; pathway median = %+.4f)",
                  bbc3_lfc, bbc3_rank, nrow(stats_pro), pro_median))
  message(sprintf("  Bbc3 |LFC| > 2x |PRO median|: %s", bbc3_is_outlier))
} else {
  bbc3_lfc <- NA_real_; bbc3_rank <- NA_integer_; bbc3_is_outlier <- NA
  message("  Bbc3 not found in Apoptosis-PRO set.")
}

# =============================================================================
# PART 4: GATE 2 CALL -- DECISION POINT vs READOUT
# =============================================================================
# Decision-point (selection) evidence: a coordinated module shift, i.e. either
# t-test significant in the expected direction. Readout: neither module shifted
# and the PRO signal (if any) isolated to Bbc3.

pro_shifted  <- pro_ttest$p.value  < 0.05
anti_shifted <- anti_ttest$p.value < 0.05
any_coordinated_shift <- pro_shifted || anti_shifted

gate2_readout <- if (any_coordinated_shift) {
  "DECISION_POINT"   # module-wide coordinated shift -> selection framing survives
} else if (isTRUE(bbc3_is_outlier)) {
  "READOUT_ISOLATED" # no group shift; Bbc3 an outlier -> remodelling / readout
} else {
  "READOUT_ABSENT"   # no group shift, no dominant driver -> metabolic-readout lean
}

gate2_call <- dplyr::case_when(
  gate2_readout == "DECISION_POINT" ~ paste(
    "Coordinated module shift detected -> mitochondrion as DECISION POINT",
    "(selection); apoptosis arm stays in main figures, AP4 retained."),
  gate2_readout == "READOUT_ISOLATED" ~ paste(
    "No module-wide shift; Bbc3 an isolated outlier -> READOUT (remodelling).",
    "Apoptosis framing shrinks; cell death held at Supp 1; lean metabolic."),
  gate2_readout == "READOUT_ABSENT" ~ paste(
    "No module-wide shift and no dominant driver -> READOUT.",
    "Apoptosis framing shrinks; cell death held at Supp 1; lean metabolic."),
  TRUE ~ "Indeterminate."
)

message("\n", strrep("=", 70))
message("GATE 2 CALL")
message(strrep("=", 70))
message(sprintf("  PRO shifted: %s  |  ANTI shifted: %s  |  Bbc3 outlier: %s",
                pro_shifted, anti_shifted, bbc3_is_outlier))
message(sprintf("  Readout class: %s", gate2_readout))
message(sprintf("  %s", gate2_call))

# =============================================================================
# PART 5: FOREST PLOT (PRO + ANTI, interaction LFC +- 95% CI)
# =============================================================================

forest_df <- dplyr::bind_rows(
  stats_pro  |> dplyr::mutate(module = "Apoptosis-PRO"),
  stats_anti |> dplyr::mutate(module = "Apoptosis-ANTI")
) |>
  dplyr::mutate(
    module = factor(module, levels = c("Apoptosis-PRO", "Apoptosis-ANTI")),
    sig    = !is.na(pvalue) & pvalue < 0.05
  )

p_forest <- ggplot(forest_df,
                   aes(x = log2FoldChange,
                       y = reorder(gene_symbol, log2FoldChange))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = log2FoldChange - 1.96 * lfcSE,
                    xmax = log2FoldChange + 1.96 * lfcSE,
                    colour = sig),
                height = 0.3, linewidth = 0.5, orientation = "y") +
  geom_point(aes(colour = sig, size = sig)) +
  facet_wrap(~ module, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey60"),
                      labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."),
                      name = "Interaction") +
  scale_size_manual(values = c("TRUE" = 3, "FALSE" = 2), guide = "none") +
  labs(
    title = "Gate 2: apoptosis module interaction LFC (Myc effect change 6W -> 12W)",
    subtitle = sprintf(
      "Negative = Myc effect weakens at 12W  |  PRO group t p = %.3f  |  ANTI group t p = %.3f  |  class: %s",
      pro_ttest$p.value, anti_ttest$p.value, gate2_readout
    ),
    x = "Interaction log2FC (timepoint12W:myc_statuspos)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave(file.path(gate_dir, "gate2_apoptosis_forest.pdf"),
       p_forest, width = 8, height = 9)

# =============================================================================
# PART 6: WRITE OUTPUTS
# =============================================================================

gate2_summary <- tibble::tibble(
  module        = c("Apoptosis-PRO", "Apoptosis-ANTI"),
  n_genes       = c(nrow(stats_pro), nrow(stats_anti)),
  n_sig_p05     = c(sum(stats_pro$pvalue < 0.05, na.rm = TRUE),
                    sum(stats_anti$pvalue < 0.05, na.rm = TRUE)),
  mean_lfc      = c(unname(pro_ttest$estimate), unname(anti_ttest$estimate)),
  median_lfc    = c(pro_median, anti_median),
  ttest_p       = c(pro_ttest$p.value, anti_ttest$p.value),
  group_shifted = c(pro_shifted, anti_shifted),
  gate2_readout = gate2_readout
)

write.csv(gate2_summary, file.path(gate_dir, "gate2_summary.csv"), row.names = FALSE)

gate2_results <- list(
  stats_pro       = stats_pro,
  stats_anti      = stats_anti,
  pro_ttest       = pro_ttest,
  anti_ttest      = anti_ttest,
  bbc3 = list(lfc = bbc3_lfc, rank = bbc3_rank, is_outlier = bbc3_is_outlier,
              pro_median = pro_median),
  gate2_readout   = gate2_readout,
  gate2_call      = gate2_call,
  summary         = gate2_summary,
  analysis_date   = Sys.Date()
)

saveRDS(gate2_results, here("results", "gate2_apoptosis_readout.rds"))

message("\nSaved:")
message("  results/gate2_apoptosis_readout.rds")
message("  outputs/gates/gate2_summary.csv")
message("  outputs/gates/gate2_apoptosis_forest.pdf")

# =============================================================================
# SANDBOX (skipped by source()/Rscript; run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  # Confirm the persisted PRO/ANTI tables carry raw interaction LFCs
  head(stats_pro) |> print()
  head(stats_anti) |> print()

  # Reproduce script 11's t-tests independently and confirm identical values
  t.test(stats_pro$log2FoldChange,  mu = 0)
  t.test(stats_anti$log2FoldChange, mu = 0)

  # Where does Bbc3 sit relative to the PRO distribution?
  stats_pro |>
    dplyr::arrange(log2FoldChange) |>
    dplyr::select(gene_symbol, log2FoldChange, lfcSE, pvalue) |>
    print(n = 30)

  # Eyeball the forest
  print(p_forest)
}
