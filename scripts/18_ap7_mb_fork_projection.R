# scripts/18_ap7_mb_fork_projection.R
# =============================================================================
# AP7 -- METABRIC mito-fork projection (Block A, Day 3, CENTREPIECE)
# =============================================================================
#
# The direct cross-species validation of the analytical companion paper
# (Menegollo, Bentham et al., Cancer Res 2024, CAN-23-3172; docs/). MCbiclust
# found large gene-cluster SWITCHES that stratify breast cancer:
#
#   - MB1 & MB2 are a shared MULTISTATE switch: their lower forks overlap almost
#     completely (MB12_LF, the common base -- low mito biogenesis, glycolytic,
#     stemness/MaSC, stress/hypoxia), then diverge into two distinct UPPER forks:
#       * MB1_UF = mito-biogenesis-high, proliferative, mature-luminal (mL),
#                  WITHOUT Myc/miRNA.
#       * MB2_UF = mito-biogenesis-high, proliferative, luminal-progenitor (LP),
#                  ER-NEGATIVE, WITH Myc/miRNA activation.
#   - MB3 = an INDEPENDENT bistate switch ('switch in function') -> specificity
#     control.
#
# Falsifiable hypothesis: Myc+ mouse tissue shifts toward MB2_UF SPECIFICALLY.
#   Because MB1_UF and MB2_UF are BOTH biogenesis-high proliferative upper forks,
#   the only feature distinguishing MB2_UF is Myc -- so the decisive test is
#   MB2_UF vs MB1_UF, which isolates 'Myc drives the Myc fork' from 'Myc drives
#   biogenesis generically'. MMTV-Myc (Myc-driven, ER-negative) landing in the
#   MB2_UF (LP / Myc / biogenesis-high) state is the point.
#
# Gene-set treatment (fixed by the paper + author):
#   - UF = group1, LF = group2, ANTICORRELATED poles -> a per-sample fork score is
#     the CONTRAST UF - LF (not an average).
#   - HI_CV lists are the fork-position definition (highest/lowest whole-
#     transcriptome CV genes) -> primary. BICLUSTER (core seed) and METAB
#     (metabolic subset) are robustness. Sets are already mouse-mapped and
#     GSVA-scored in results/gsva_scores.rds (consume as-is; no re-biclustering,
#     no MCbiclust needed).
#
# Statistics: the primary test is the genotype MAIN EFFECT (Myc+ vs Myc-,
#   timepoint-adjusted) -- the powered regime (cf. Gate 1), a single pre-specified
#   externally-anchored hypothesis (MB2 specifically), not the underpowered
#   interaction. Trajectory (6W vs 12W) is secondary/descriptive; biogenesis is
#   front-loaded at 6W so MB2_UF resemblance is predicted to peak at 6W.
#
# Input:
#   - results/gsva_scores.rds : scores [set x 24] + sample_meta
# Output:
#   - results/ap7_mb_fork.rds : fork scores, tests, robustness, notes
#   - outputs/ap7/ : mb2_over_mb1_boxplot.pdf, fork_space_projection.pdf,
#                    mb2_trajectory.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# =============================================================================
# PART 1: LOAD + SAMPLE METADATA
# =============================================================================

gsva_out    <- readRDS(here::here("results", "gsva_scores.rds"))
scores      <- gsva_out$scores
sample_meta <- as.data.frame(gsva_out$sample_meta)
sample_meta <- sample_meta[colnames(scores), , drop = FALSE]
stopifnot(identical(rownames(sample_meta), colnames(scores)))

sample_meta$timepoint  <- stats::relevel(as.factor(sample_meta$timepoint),  "6W")
sample_meta$myc_status <- stats::relevel(as.factor(sample_meta$myc_status), "neg")
sample_meta$group      <- factor(sample_meta$group,
                                 levels = c("6W_neg", "6W_pos", "12W_neg", "12W_pos"))

# =============================================================================
# PART 2: MB SET NAMES + PRESENCE CHECK
# =============================================================================

mb_set <- function(subtype, type, group) {
  sprintf("METABRIC_%s_%s_GROUP%d", subtype, type, group)
}
subtypes <- c("MB1", "MB2", "MB3")
types    <- c("HI_CV", "BICLUSTER", "METAB")

needed <- unlist(lapply(subtypes, function(s)
  unlist(lapply(types, function(t) c(mb_set(s, t, 1), mb_set(s, t, 2))))))
missing_sets <- setdiff(needed, rownames(scores))
if (length(missing_sets) > 0) {
  stop("MB fork sets missing from gsva_scores.rds: ",
       paste(missing_sets, collapse = ", "),
       " -- re-run scripts/15_gsva_scoring.R.")
}
message(sprintf("All %d MB fork sets present in GSVA scores.", length(needed)))

sc_row <- function(name) scores[name, ]

# =============================================================================
# PART 3: PER-SAMPLE FORK SCORES (UF - LF contrasts)
# =============================================================================
# fork(subtype, type) = GSVA(UF = group1) - GSVA(LF = group2): position on the
# switch (positive = toward the UF).

fork <- function(subtype, type) sc_row(mb_set(subtype, type, 1)) -
                                 sc_row(mb_set(subtype, type, 2))

# Primary (HI_CV) fork-position scores
MB1_score <- fork("MB1", "HI_CV")
MB2_score <- fork("MB2", "HI_CV")
MB3_score <- fork("MB3", "HI_CV")

# The decisive specificity contrast: MB2_UF vs MB1_UF (both share ~MB12_LF, so
# this is pure upper-fork divergence -- the Myc fork vs the non-Myc mito fork).
MB2_over_MB1_UF <- sc_row(mb_set("MB2", "HI_CV", 1)) -
                   sc_row(mb_set("MB1", "HI_CV", 1))

fork_df <- tibble::tibble(
  sample          = colnames(scores),
  group           = sample_meta$group,
  myc_status      = sample_meta$myc_status,
  timepoint       = sample_meta$timepoint,
  MB1_score       = MB1_score,
  MB2_score       = MB2_score,
  MB3_score       = MB3_score,
  MB2_over_MB1_UF = MB2_over_MB1_UF
)

# =============================================================================
# PART 4: TESTS -- genotype main effect (powered) + trajectory (secondary)
# =============================================================================
# Force treatment contrasts so coefficient names are stable regardless of the
# session's global options(contrasts) (a real trap here; see script 17).

tc <- list(myc_status = "contr.treatment", timepoint = "contr.treatment")

test_metric <- function(y) {
  d     <- data.frame(y = y, myc_status = sample_meta$myc_status,
                      timepoint = sample_meta$timepoint)
  m_add <- stats::lm(y ~ myc_status + timepoint, data = d, contrasts = tc)
  m_int <- stats::lm(y ~ timepoint * myc_status, data = d, contrasts = tc)
  sa    <- summary(m_add)$coefficients
  si    <- summary(m_int)$coefficients
  tibble::tibble(
    geno_pooled_beta = sa["myc_statuspos", "Estimate"],   # PRIMARY (powered)
    geno_pooled_p    = sa["myc_statuspos", "Pr(>|t|)"],
    myc_6W_beta      = si["myc_statuspos", "Estimate"],    # Myc effect at 6W
    myc_6W_p         = si["myc_statuspos", "Pr(>|t|)"],
    int_beta         = si["timepoint12W:myc_statuspos", "Estimate"],  # trajectory
    int_p            = si["timepoint12W:myc_statuspos", "Pr(>|t|)"]
  )
}

metric_list <- list(
  MB2_over_MB1_UF = MB2_over_MB1_UF,   # decisive specificity contrast
  MB2_score       = MB2_score,         # Myc-fork switch position
  MB1_score       = MB1_score,         # specificity control (non-Myc mito fork)
  MB3_score       = MB3_score          # specificity control (independent switch)
)
tests <- do.call(rbind, lapply(names(metric_list), function(nm) {
  cbind(metric = nm, test_metric(metric_list[[nm]]))
}))
tests <- tibble::as_tibble(tests)

# Simple robust headline: Myc+ vs Myc- on the decisive contrast (Wilcoxon).
wilcox_mb2 <- suppressWarnings(stats::wilcox.test(
  MB2_over_MB1_UF ~ sample_meta$myc_status))

# =============================================================================
# PART 5: ROBUSTNESS ACROSS SUB-SIGNATURE TYPES (HI_CV / BICLUSTER / METAB)
# =============================================================================
# Does the MB2 genotype effect hold when the fork is defined by the core
# bicluster or the metabolic subset, not just HI_CV?

robustness <- do.call(rbind, lapply(types, function(t) {
  mb2_t <- fork("MB2", t)
  mb2_over_t <- sc_row(mb_set("MB2", t, 1)) - sc_row(mb_set("MB1", t, 1))
  data.frame(
    type              = t,
    MB2_score_geno_p  = test_metric(mb2_t)$geno_pooled_p,
    MB2_score_geno_b  = test_metric(mb2_t)$geno_pooled_beta,
    MB2_over_MB1_geno_p = test_metric(mb2_over_t)$geno_pooled_p,
    MB2_over_MB1_geno_b = test_metric(mb2_over_t)$geno_pooled_beta
  )
}))
robustness <- tibble::as_tibble(robustness)

# =============================================================================
# PART 6: FIGURES
# =============================================================================

out_dir <- here::here("outputs", "ap7")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

geno_cols <- c(neg = "#4575B4", pos = "#D73027")
mb2_geno_p <- tests$geno_pooled_p[tests$metric == "MB2_over_MB1_UF"]

# Fig A: the decisive contrast by the 4 groups (named figure)
p_box <- ggplot2::ggplot(fork_df,
    ggplot2::aes(x = group, y = MB2_over_MB1_UF, fill = myc_status)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  ggplot2::geom_jitter(width = 0.12, size = 2, ggplot2::aes(colour = myc_status)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::scale_fill_manual(values = geno_cols) +
  ggplot2::scale_colour_manual(values = geno_cols) +
  ggplot2::labs(
    title = "AP7: shift toward MB2_UF (Myc fork) vs MB1_UF (non-Myc mito fork)",
    subtitle = sprintf("MB2_UF - MB1_UF (HI_CV); genotype main effect p = %.3g", mb2_geno_p),
    x = NULL, y = "MB2_UF - MB1_UF  (positive = toward the Myc fork)") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "mb2_over_mb1_boxplot.pdf"), p_box,
                width = 7, height = 5)

# Fig B: fork-space projection (MB1 vs MB2 switch position), coloured by group
p_space <- ggplot2::ggplot(fork_df,
    ggplot2::aes(x = MB1_score, y = MB2_score, colour = group)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(
    title = "Fork-space projection (UF - LF switch positions)",
    subtitle = "Myc+ predicted high on MB2 (toward MB2_UF), not specifically MB1",
    x = "MB1 switch position (UF - LF)", y = "MB2 switch position (UF - LF)") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "fork_space_projection.pdf"), p_space,
                width = 6.5, height = 5.5)

# Fig C: MB2 trajectory (6W->12W) by genotype -- the 6W-peak prediction
traj <- fork_df |>
  dplyr::group_by(timepoint, myc_status) |>
  dplyr::summarise(MB2 = mean(MB2_score),
                   MB2_over = mean(MB2_over_MB1_UF), .groups = "drop")
p_traj <- ggplot2::ggplot(traj,
    ggplot2::aes(x = timepoint, y = MB2_over, colour = myc_status,
                 group = myc_status)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::scale_colour_manual(values = geno_cols, name = "Myc") +
  ggplot2::labs(
    title = "MB2_UF - MB1_UF trajectory (prediction: Myc+ peaks at 6W)",
    x = NULL, y = "mean MB2_UF - MB1_UF") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "mb2_trajectory.pdf"), p_traj,
                width = 6, height = 4.5)

# =============================================================================
# PART 7: SAVE
# =============================================================================

group_means <- fork_df |>
  dplyr::group_by(group) |>
  dplyr::summarise(dplyr::across(c(MB1_score, MB2_score, MB3_score, MB2_over_MB1_UF),
                                 mean), .groups = "drop")

ap7_out <- list(
  fork_df      = fork_df,
  tests        = tests,
  robustness   = robustness,
  wilcox_mb2   = wilcox_mb2,
  group_means  = group_means,
  notes        = paste(
    "Fork scores = GSVA(UF=group1) - GSVA(LF=group2), HI_CV primary.",
    "MB2_over_MB1_UF = MB2_UF - MB1_UF = the decisive Myc-fork specificity",
    "contrast (both UFs share ~MB12_LF). PRIMARY test = genotype main effect",
    "(geno_pooled_*, timepoint-adjusted, powered); myc_6W_* = Myc effect at 6W;",
    "int_* = trajectory (secondary/underpowered). Falsifier: MB2 metrics positive",
    "and MB1_score/MB3_score not. Sets pre-mapped + GSVA-scored (script 15).")
)
saveRDS(ap7_out, here::here("results", "ap7_mb_fork.rds"))
message("Saved results/ap7_mb_fork.rds")

message(sprintf("AP7 primary: MB2_UF-vs-MB1_UF genotype effect beta = %.3f, p = %.3g",
                tests$geno_pooled_beta[tests$metric == "MB2_over_MB1_UF"],
                tests$geno_pooled_p[tests$metric == "MB2_over_MB1_UF"]))

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  ap7 <- readRDS(here::here("results", "ap7_mb_fork.rds"))

  # --- Validate the gene-set treatment against the paper's topology ---
  # (1) UF and LF are anticorrelated poles -> group1 vs group2 should be NEGATIVELY
  #     correlated across samples (justifies the UF - LF contrast).
  for (s in c("MB1", "MB2", "MB3")) {
    r <- stats::cor(scores[mb_set(s, "HI_CV", 1), ],
                    scores[mb_set(s, "HI_CV", 2), ])
    cat(sprintf("%s UF vs LF (HI_CV) correlation: %.2f (expect negative)\n", s, r))
  }
  # (2) MB1_LF and MB2_LF overlap (MB12_LF) -> their LF scores should be strongly
  #     POSITIVELY correlated across samples.
  cat(sprintf("MB1_LF vs MB2_LF correlation: %.2f (expect high positive = MB12_LF)\n",
              stats::cor(scores[mb_set("MB1", "HI_CV", 2), ],
                         scores[mb_set("MB2", "HI_CV", 2), ])))

  # --- The result: specificity table (MB2 positive + significant; MB1/MB3 not) ---
  ap7$tests |> print()
  ap7$group_means |> print()

  # Robustness: does MB2 hold across HI_CV / BICLUSTER / METAB?
  ap7$robustness |> print()

  # Headline Wilcoxon (Myc+ vs Myc-) on the decisive contrast
  ap7$wilcox_mb2

  # Trajectory: is MB2_UF resemblance strongest at 6W in Myc+? (int_beta < 0)
  ap7$tests |> dplyr::filter(metric %in% c("MB2_over_MB1_UF", "MB2_score")) |>
    dplyr::select(metric, myc_6W_beta, int_beta, int_p) |> print()

  # Optional: confirm the raw fork sets' structure / gene lists
  # mbs <- readRDS(here::here("data", "genesets_from_library", "metabric_sets.rds"))
  # str(mbs, max.level = 2)

  list.files(here::here("outputs", "ap7"), pattern = "\\.pdf$")
}
