# scripts/13_gate1_divergence_timing.R
# =============================================================================
# GATE 1 (Block A decision layer): Divergence timing across the 6W -> 12W window
# =============================================================================
#
# Purpose (AP5 / Gate 1 in docs/2026-07-04_block_A_build_spec.md):
#   Resolve the "permission vs divergence" fork in the plan's framing. The
#   question is WHEN the Myc+ / Myc- genotypes diverge:
#
#     (i)  How much is already divergent at the EARLIEST sampled timepoint (6W)?
#          -> genotype (myc_status) main effect at 6W (the model reference).
#     (ii) Does that genotype divergence GROW, stay FLAT, or SHRINK from
#          6W -> 12W?
#          -> compare genotype effect at 6W vs 12W, and read the interaction
#             direction (how the Myc effect itself changes over the window).
#
#   Gate logic:
#     - If divergence GROWS across the window, "developmental change licenses
#       Myc" weakens -> the title verb moves off "licenses".
#     - If divergence is ALREADY-LARGE-AND-FLAT at 6W, the soft "already
#       diverged" claim holds.
#   This feeds the Day-1 figure lock and the dated gates note.
#
#   Sign convention (matches script 11):
#     interaction log2FC (timepoint12W.myc_statuspos)
#       negative = Myc effect WEAKENS at 12W vs 6W
#       positive = Myc effect STRENGTHENS at 12W
#
# This is a REFRAME script: it reads existing DESeq2 / characterisation outputs
# and collates them. No new model fitting, no shrinkage. All contrasts read are
# the RAW (unshrunken MLE) results; shrinkage does not enter divergence-timing.
#
# Input:
#   - results/interaction_results.rds            (script 03; raw contrasts)
#       myc_6W_raw    = genotype effect at 6W  (myc_status_pos_vs_neg)
#       myc_12W_raw   = genotype effect at 12W (composite contrast)
#       timepoint_pos_raw / timepoint_neg_raw = temporal effect within genotype
#       interaction_raw = timepoint12W.myc_statuspos
#   - results/interaction_gene_characterisation.rds (script 11; $direction)
#   - results/ortholog_table.rds                 (ensembl -> symbol labels)
#
# Output:
#   - results/gate1_divergence_timing.rds
#   - outputs/gates/gate1_summary.csv
#   - outputs/gates/gate1_divergence_timing.pdf
#
# NOTE (Day-3 hook): the stage-annotated developmental-set projection depends on
#   the GSVA build (AP2, script 15) and dev composition (script 18). It is left
#   here as an explicit placeholder (`dev_projection_hook`) to be filled from
#   those outputs on Day 3; it is NOT computed in this script.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

gate_dir <- here("outputs", "gates")
dir.create(gate_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PART 1: LOAD
# =============================================================================

message("Loading interaction results and characterisation outputs...")

interaction_results <- readRDS(here("results", "interaction_results.rds"))
characterisation    <- readRDS(here("results", "interaction_gene_characterisation.rds"))
ortholog_table      <- readRDS(here("results", "ortholog_table.rds"))

ensembl_to_symbol <- setNames(
  ortholog_table$external_gene_name,
  ortholog_table$ensembl_gene_id
)

# Contrasts used here (all RAW / unshrunken)
myc_6W    <- as.data.frame(interaction_results$myc_6W_raw)     # genotype @ 6W
myc_12W   <- as.data.frame(interaction_results$myc_12W_raw)    # genotype @ 12W
int_raw   <- as.data.frame(interaction_results$interaction_raw)

# =============================================================================
# PART 2 (i): HOW MUCH IS ALREADY DIVERGENT AT 6W?
# =============================================================================
# Genotype (Myc+ vs Myc-) main effect at 6W = myc_status_pos_vs_neg, the
# reference-timepoint contrast. Count genes already separating the genotypes
# at the earliest sample, at two FDR thresholds.

message("\n", strrep("=", 70))
message("PART 2 (i): GENOTYPE DIVERGENCE AT 6W")
message(strrep("=", 70))

count_sig <- function(df, thr) sum(df$padj < thr, na.rm = TRUE)

divergence_6W <- tibble::tibble(
  contrast     = "genotype @ 6W (myc_status_pos_vs_neg)",
  n_tested     = sum(!is.na(myc_6W$padj)),
  n_padj_0.1   = count_sig(myc_6W, 0.1),
  n_padj_0.05  = count_sig(myc_6W, 0.05),
  median_abs_lfc_sig = median(abs(myc_6W$log2FoldChange[which(myc_6W$padj < 0.1)]), na.rm = TRUE)
)

message(sprintf("  Genes divergent at 6W (padj<0.1):  %d", divergence_6W$n_padj_0.1))
message(sprintf("  Genes divergent at 6W (padj<0.05): %d", divergence_6W$n_padj_0.05))

# =============================================================================
# PART 2 (ii): DOES DIVERGENCE GROW, STAY FLAT, OR SHRINK 6W -> 12W?
# =============================================================================
# Two complementary reads:
#   (a) Genotype effect at 6W vs 12W: counts + magnitude side by side.
#   (b) Interaction direction: fraction of interaction-signal genes whose Myc
#       effect weakens (negative) vs strengthens (positive) at 12W.

message("\n", strrep("=", 70))
message("PART 2 (ii): DIVERGENCE TREND ACROSS THE WINDOW")
message(strrep("=", 70))

divergence_by_timepoint <- tibble::tibble(
  timepoint    = c("6W", "12W"),
  contrast     = c("myc_status_pos_vs_neg", "myc @ 12W (composite)"),
  n_padj_0.1   = c(count_sig(myc_6W, 0.1), count_sig(myc_12W, 0.1)),
  n_padj_0.05  = c(count_sig(myc_6W, 0.05), count_sig(myc_12W, 0.05)),
  median_abs_lfc_sig = c(
    median(abs(myc_6W$log2FoldChange[which(myc_6W$padj < 0.1)]), na.rm = TRUE),
    median(abs(myc_12W$log2FoldChange[which(myc_12W$padj < 0.1)]), na.rm = TRUE)
  )
)

# Trend classification on the genotype-divergent gene count
delta_n <- divergence_by_timepoint$n_padj_0.1[2] - divergence_by_timepoint$n_padj_0.1[1]
rel_change <- delta_n / max(divergence_by_timepoint$n_padj_0.1[1], 1)
divergence_trend <- dplyr::case_when(
  rel_change >  0.15 ~ "GROWS",
  rel_change < -0.15 ~ "SHRINKS",
  TRUE               ~ "FLAT"
)

message(sprintf("  Genotype-divergent genes (padj<0.1): 6W=%d, 12W=%d (delta=%+d, %.0f%%) -> %s",
                divergence_by_timepoint$n_padj_0.1[1],
                divergence_by_timepoint$n_padj_0.1[2],
                delta_n, 100 * rel_change, divergence_trend))

# --- Interaction direction (reuse script 11's persisted $direction summary) ---
dir_summary <- characterisation$direction
message(sprintf("\n  Interaction-signal genes (p<0.05, unique-to-Myc+, from script 11):"))
message(sprintf("    n = %d  |  negative (Myc weakens): %d (%.1f%%)  |  positive: %d",
                dir_summary$n_negative + dir_summary$n_positive,
                dir_summary$n_negative, dir_summary$pct_negative,
                dir_summary$n_positive))
message(sprintf("    median interaction LFC = %.3f  |  binomial p = %.2e",
                dir_summary$median_lfc, dir_summary$binom_p))

# Direction on ALL genes with a nominal interaction signal (context, not the 357)
int_all_sig <- int_raw[which(int_raw$pvalue < 0.05), , drop = FALSE]
n_int_neg <- sum(int_all_sig$log2FoldChange < 0, na.rm = TRUE)
n_int_pos <- sum(int_all_sig$log2FoldChange > 0, na.rm = TRUE)
pct_int_neg_all <- 100 * n_int_neg / (n_int_neg + n_int_pos)
message(sprintf("  Context - all genes with interaction p<0.05: %d neg (%.1f%%), %d pos",
                n_int_neg, pct_int_neg_all, n_int_pos))

# =============================================================================
# PART 2b: COMMON-GENE-SET EFFECT-SIZE CHECK (resolve the power caveat / H2)
# =============================================================================
# The 2777 -> 239 drop in PART 2 is a significant-gene COUNT. Because myc_12W is a
# composite (list) contrast with a larger SE than the 6W reference-level contrast,
# the count drop is confounded with lower 12W power. A LFC point estimate is
# UNAFFECTED by SE, so comparing |LFC| across a FIXED gene set (not re-thresholded
# at 12W) discriminates:
#   - |LFC| magnitude PRESERVED 6W -> 12W  -> count drop is a POWER artifact (H2)
#   - |LFC| magnitude COLLAPSES 6W -> 12W  -> BIOLOGICAL convergence (H1/H3/H4)
#
# Two complementary sets:
#   (A) 6W-divergent set (padj<0.1 at 6W): "did the divergent genes stay
#       divergent?". CAVEAT: selecting on 6W significance inflates the 6W |LFC|
#       (winner's curse / regression-to-mean), biasing this comparison TOWARD
#       apparent attenuation. So ratio >= ~0.9 here is STRONG evidence for a power
#       artifact; ratio < 1 is only SUGGESTIVE of real convergence.
#   (B) expression-selected set (baseMean >= median), independent of the genotype
#       LFC: an UNBIASED read of the 12W-vs-6W LFC relationship (slope, corr),
#       free of the selection-on-effect bias in (A).

message("\n", strrep("=", 70))
message("PART 2b: COMMON-GENE-SET EFFECT-SIZE CHECK (power caveat / H2)")
message(strrep("=", 70))

# Joined genotype LFCs across the two contrasts (shared gene universe from dds_int)
lfc_join <- data.frame(
  ensembl_id = rownames(myc_6W),
  baseMean   = myc_6W$baseMean,
  lfc_6W     = myc_6W$log2FoldChange,
  padj_6W    = myc_6W$padj,
  lfc_12W    = myc_12W$log2FoldChange,
  row.names  = NULL
) |>
  dplyr::filter(!is.na(lfc_6W), !is.na(lfc_12W))

summarise_effize <- function(df, set_label) {
  # paired Wilcoxon signed-rank on |LFC| across the SAME genes
  wp <- suppressWarnings(
    wilcox.test(abs(df$lfc_12W), abs(df$lfc_6W), paired = TRUE)
  )$p.value
  tibble::tibble(
    gene_set           = set_label,
    n                  = nrow(df),
    median_abs_lfc_6W  = median(abs(df$lfc_6W)),
    median_abs_lfc_12W = median(abs(df$lfc_12W)),
    ratio_12W_over_6W  = median(abs(df$lfc_12W)) / median(abs(df$lfc_6W)),
    slope_12W_on_6W    = unname(coef(lm(lfc_12W ~ 0 + lfc_6W, data = df))[1]),
    pearson_r          = cor(df$lfc_6W, df$lfc_12W),
    frac_attenuated    = mean(abs(df$lfc_12W) < abs(df$lfc_6W)),
    paired_wilcox_p    = wp
  )
}

set_A <- lfc_join |> dplyr::filter(padj_6W < 0.1)                 # 6W-divergent (biased toward attenuation)
set_B <- lfc_join |> dplyr::filter(baseMean >= median(baseMean)) # expressed (LFC-independent, unbiased)

effize_A     <- summarise_effize(set_A, "6W-divergent (padj<0.1)")
effize_B     <- summarise_effize(set_B, "expressed (baseMean>=median)")
effize_check <- dplyr::bind_rows(effize_A, effize_B)

effize_check |> print()

# --- Verdict: magnitude collapse (biology) vs preserved magnitude (power) ---
ratioA <- effize_A$ratio_12W_over_6W
effize_class <- dplyr::case_when(
  ratioA >= 0.9 ~ "PRESERVED_POWER_ARTIFACT",
  ratioA <= 0.6 ~ "COLLAPSED_BIOLOGICAL",
  TRUE          ~ "PARTIAL"
)

effize_verdict <- switch(
  effize_class,
  PRESERVED_POWER_ARTIFACT = sprintf(
    paste("PRESERVED: genotype |LFC| on the 6W-divergent set holds at 12W",
          "(ratio %.2f, slope %.2f, r %.2f). The count drop is largely a POWER",
          "artifact (H2); effect-size 'convergence' NOT supported."),
    effize_A$ratio_12W_over_6W, effize_A$slope_12W_on_6W, effize_A$pearson_r),
  COLLAPSED_BIOLOGICAL = sprintf(
    paste("COLLAPSED: genotype |LFC| falls at 12W (ratio %.2f, slope %.2f, r %.2f).",
          "Supports BIOLOGICAL convergence (H1/H3/H4). 6W-selection biases toward",
          "attenuation - corroborated by the expressed-set slope %.2f."),
    effize_A$ratio_12W_over_6W, effize_A$slope_12W_on_6W, effize_A$pearson_r,
    effize_B$slope_12W_on_6W),
  PARTIAL = sprintf(
    paste("PARTIAL: |LFC| ratio %.2f intermediate (slope %.2f, r %.2f); the count",
          "drop mixes power + magnitude. Expressed-set slope %.2f is the unbiased",
          "tiebreak."),
    effize_A$ratio_12W_over_6W, effize_A$slope_12W_on_6W, effize_A$pearson_r,
    effize_B$slope_12W_on_6W)
)

message(sprintf("\n  Class: %s", effize_class))
message(sprintf("  %s", effize_verdict))

# --- Money plot: genotype LFC 6W vs 12W on the 6W-divergent set ---
p_effize <- ggplot(set_A, aes(x = lfc_6W, y = lfc_12W)) +
  geom_point(alpha = 0.25, size = 0.7, colour = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  geom_abline(slope = effize_A$slope_12W_on_6W, intercept = 0,
              colour = "darkblue", linewidth = 0.7) +
  coord_equal() +
  labs(
    title = "Gate 1 effect-size check: genotype LFC, 6W-divergent genes",
    subtitle = sprintf(
      paste0("n=%d  |  median |LFC| 6W=%.2f, 12W=%.2f (ratio %.2f)  |  ",
             "slope(12W~6W)=%.2f, r=%.2f\n",
             "Red dashed = y=x (magnitude preserved)  |  Blue = fitted through-origin slope"),
      effize_A$n, effize_A$median_abs_lfc_6W, effize_A$median_abs_lfc_12W,
      effize_A$ratio_12W_over_6W, effize_A$slope_12W_on_6W, effize_A$pearson_r),
    x = "Genotype LFC at 6W (myc_status_pos_vs_neg)",
    y = "Genotype LFC at 12W (composite)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(gate_dir, "gate1_effect_size_check.pdf"), p_effize, width = 7, height = 7)
write.csv(effize_check, file.path(gate_dir, "gate1_effect_size_check.csv"),
          row.names = FALSE)

# =============================================================================
# PART 3: DAY-3 HOOK (developmental-set stage projection) -- placeholder only
# =============================================================================
# Filled on Day 3 from results/gsva_scores.rds (script 15) + results/
# dev_composition.rds (script 18). Left as a labelled stub so the gate object
# has a stable slot; NOT computed here.

dev_projection_hook <- list(
  status = "PENDING_DAY3",
  fill_from = c("results/gsva_scores.rds", "results/dev_composition.rds"),
  note = paste(
    "Stage-annotated developmental-set GSVA projection: place WT 6W/12W and",
    "Myc+/- samples on the developmental axis to test whether genotype",
    "separation tracks developmental stage. Depends on AP2 (script 18)."
  )
)

# =============================================================================
# PART 4: GATE 1 CALL
# =============================================================================

gate1_call <- dplyr::case_when(
  divergence_trend == "GROWS" ~ paste(
    "Divergence GROWS across 6W->12W: 'developmental change licenses Myc'",
    "weakens; move the title verb off 'licenses'."),
  divergence_trend == "FLAT" ~ paste(
    "Divergence ALREADY-LARGE-AND-FLAT at 6W: the soft 'already diverged'",
    "claim holds."),
  divergence_trend == "SHRINKS" ~ paste(
    "Divergence SHRINKS across the window: genotypes converge by 12W;",
    "reconsider the trajectory framing."),
  TRUE ~ "Indeterminate."
)

message("\n", strrep("=", 70))
message("GATE 1 CALL")
message(strrep("=", 70))
message(sprintf("  Trend: %s", divergence_trend))
message(sprintf("  %s", gate1_call))

# =============================================================================
# PART 5: PLOT
# =============================================================================

plot_df <- divergence_by_timepoint |>
  dplyr::mutate(timepoint = factor(timepoint, levels = c("6W", "12W")))

p_divergence <- ggplot(plot_df, aes(x = timepoint, y = n_padj_0.1)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_text(aes(label = n_padj_0.1), vjust = -0.4, size = 4) +
  labs(
    title = "Gate 1: genotype divergence across the timecourse",
    subtitle = sprintf(
      "Genes with genotype (Myc+/-) effect padj<0.1  |  trend: %s  |  interaction: %.0f%% weaken at 12W",
      divergence_trend, dir_summary$pct_negative
    ),
    x = "Timepoint (genotype contrast)",
    y = "Divergent genes (padj < 0.1)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(gate_dir, "gate1_divergence_timing.pdf"),
       p_divergence, width = 7, height = 5)

# =============================================================================
# PART 6: WRITE OUTPUTS
# =============================================================================

gate1_summary <- divergence_by_timepoint |>
  dplyr::mutate(
    divergence_trend = divergence_trend,
    pct_interaction_negative = dir_summary$pct_negative,
    interaction_binom_p = dir_summary$binom_p
  )

write.csv(gate1_summary, file.path(gate_dir, "gate1_summary.csv"), row.names = FALSE)

gate1_results <- list(
  divergence_6W          = divergence_6W,
  divergence_by_timepoint = divergence_by_timepoint,
  divergence_trend       = divergence_trend,
  interaction_direction  = dir_summary,
  interaction_direction_all_context = list(
    n_negative = n_int_neg, n_positive = n_int_pos, pct_negative = pct_int_neg_all
  ),
  effect_size_check      = effize_check,
  effect_size_class      = effize_class,
  effect_size_verdict    = effize_verdict,
  dev_projection_hook    = dev_projection_hook,
  gate1_call             = gate1_call,
  analysis_date          = Sys.Date()
)

saveRDS(gate1_results, here("results", "gate1_divergence_timing.rds"))

message("\nSaved:")
message("  results/gate1_divergence_timing.rds")
message("  outputs/gates/gate1_summary.csv")
message("  outputs/gates/gate1_divergence_timing.pdf")
message("  outputs/gates/gate1_effect_size_check.csv")
message("  outputs/gates/gate1_effect_size_check.pdf")

# =============================================================================
# SANDBOX (skipped by source()/Rscript; run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  # Inspect the raw contrasts directly
  head(myc_6W) |> print()
  head(int_raw) |> print()

  # Confirm the genotype-@-6W contrast really is the reference-timepoint effect
  # (should equal script 03's myc_status_pos_vs_neg 'name' extraction)
  summary(interaction_results$myc_6W_raw)

  # Top genes already divergent at 6W (label with symbols)
  myc_6W |>
    tibble::rownames_to_column("ensembl_id") |>
    dplyr::mutate(symbol = ensembl_to_symbol[ensembl_id]) |>
    dplyr::filter(padj < 0.05) |>
    dplyr::arrange(dplyr::desc(abs(log2FoldChange))) |>
    head(20) |>
    print()

  # Sanity: interaction direction split reproduced from script 11's persisted set
  characterisation$direction

  # --- Effect-size check (PART 2b) inspection ---
  # Two summary rows: 6W-divergent (biased toward attenuation) + expressed (unbiased).
  #  - set_A ratio near/above 1 or slope near 1  => power artifact (H2)
  #  - set_A ratio well below 1 AND set_B slope below 1 => real convergence
  effize_check |> print()

  # Eyeball both plots
  print(p_divergence)
  print(p_effize)
}
