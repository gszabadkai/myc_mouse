# scripts/30_attenuation_moderation.R
# =============================================================================
# Block A revision -- Issue #5: can MYC-independent factors EXPLAIN the OXPHOS attenuation?
# =============================================================================
#
# Issue #4 (script 29) resolved WHAT the attenuation is (magnitude shrinks, not rank) and
# verified in absolute mRNA that Myc RAISES nuclear OXPHOS (geno d~1.27) -- so mitoPPS
# "OXPHOS down" is reprioritisation, not an absolute drop. Its Part D then found the WT
# OXPHOS trajectory co-varies within WT with the developmental luminal axis (partial 0.71)
# and ER/PGC1a TF activity (0.70-0.83). The author then asked the SHARPER question: can we
# "blame" MYC-INDEPENDENT factors (development, other TFs) to FULLY EXPLAIN the reduction?
#
# WHY PART D WAS THE WRONG INSTRUMENT: it tested LEVEL association (does the axis track
# OXPHOS?), not whether the axis ABSORBS the attenuation. The attenuation is a MODERATION
# phenomenon -- the genotype (Myc) effect DEPENDS on timepoint = the genotype x time
# INTERACTION (b_int). "Explaining the reduction" = does conditioning on a candidate
# MYC-independent axis SHRINK b_int toward 0? That is the test Issue #5 runs.
#
# THE HONEST ANSWER, up front: NO on current evidence, and this script BOUNDS rather than
# settles. Three reasons it can only bound:
#   (1) ENDOGENEITY -- Myc drives the candidate axes too (Issue #2 endogenous Myc gates the
#       pubertal program; script 24 Myc co-opts ESRRA/NRF1/GABPA). Conditioning on a
#       Myc-driven mediator BIASES the absorption estimate (over-controls).
#   (2) NEAR-TAUTOLOGY -- the two timepoints ARE two developmental stages, so "developmental
#       stage moderates Myc" is nearly definitional; WHICH molecular factor cannot be
#       separated here (n=6/group; dev axis, TF activity, proliferation, Myc's own targets
#       collinear -- Felsher~biogenesis rho 0.90).
#   (3) POWER FLOOR -- the per-sample composite interaction is itself DIRECTIONAL not
#       significant (VST b_int -0.279, int_p 0.61, n=6/tp; script 29 absolute_stats). A small,
#       noisy baseline makes the absorption FRACTION (1 - b_int_adj/b_int_base) numerically
#       UNSTABLE (small denominator). So the fraction is reported WITH its bootstrap CI, but
#       the PRIMARY metrics are the raw change Delta-b_int and the nested-model LRT, which do
#       not divide by a near-zero baseline.
#
# PRIOR THAT POINTS AGAINST ABSORPTION: script 17 PART 6 adjusted the Myc+ biogenesis
# 6W->12W decline for the developmental axis and it did NOT shrink -- it GREW (-0.075 ->
# -0.143). Development CO-DRIVES the program (positively correlated), it is not an
# antagonist that absorbs the Myc effect. Part D here extends that from biogenesis to OXPHOS
# and from the Myc+ slope to the full interaction.
#
# Reframe on already-fitted data (script 29 saved set-lists in $defs; GSVA/mitoPPS scores;
# dds for VST). NO DESeq/GSVA re-run. Construction MATCHES script 29 exactly (same
# composite()/comp_gsva() helpers, same set-lists) so the baseline b_int reproduces
# absolute_stats.
#
# Input:  results/attenuation_decomposition.rds  ($defs: nuclear_oxphos_ens, luminal_sets,
#                                                  tf_sets, mito_ox_gsva_sets)
#         results/dds_int_run.rds                (VST source for the absolute OXPHOS composite)
#         results/gsva_scores.rds                (dev / TF / proliferation / OXPHOS composites)
#         results/mitopps_scores.rds             (mtDNA reprioritisation axis)
# Output: results/attenuation_moderation.rds
#         outputs/attenuation_moderation/*.pdf
#
# CEILING (state prominently; do NOT let a large absorption fraction be read as "development
# explains it"): endogenous mediators -> BIASED; n=6/group; collinear axes; baseline
# interaction ns. This test BOUNDS the claim (expected: the axes do NOT fully absorb the
# attenuation, echoing script 17). Association not causation; the identifying evidence is
# experimental (Part E). Report the CI and the endogeneity caveat together.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# =============================================================================
# PART 1: LOAD + REBUILD THE SCRIPT-29 CONSTRUCTION (exact-match set-lists)
# =============================================================================
ad <- readRDS(here::here("results", "attenuation_decomposition.rds"))
defs <- ad$defs

dds <- readRDS(here::here("results", "dds_int_run.rds"))
sm  <- as.data.frame(SummarizedExperiment::colData(dds))
sm$timepoint  <- stats::relevel(as.factor(sm$timepoint),  "6W")
sm$myc_status <- stats::relevel(as.factor(sm$myc_status), "neg")
sm$group      <- factor(sm$group, levels = c("6W_neg", "6W_pos", "12W_neg", "12W_pos"))
samples       <- colnames(dds)

gsva_out <- readRDS(here::here("results", "gsva_scores.rds"))
scores   <- gsva_out$scores[, samples, drop = FALSE]
set_meta <- gsva_out$set_meta
mp       <- readRDS(here::here("results", "mitopps_scores.rds"))

# VST matrix (deterministic transform of the fitted dds), aligned to master sample order
vst_mat <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))[, samples, drop = FALSE]

# script-29 helpers, verbatim behaviour
zrow      <- function(m) t(scale(t(m)))
composite <- function(m, genes) colMeans(zrow(m[intersect(genes, rownames(m)), , drop = FALSE]))
comp_gsva <- function(sets) {
  s <- intersect(sets, rownames(scores)); stopifnot(length(s) > 0)
  colMeans(scores[s, , drop = FALSE])
}

# proliferation composite = the 14 Proliferation-category sets (matches script 27 prol_sets)
prol_sets <- set_meta$set_name[set_meta$category_primary == "Proliferation"]

# =============================================================================
# PART 2: PER-SAMPLE MODERATION FRAME (all 24 samples)
# =============================================================================
# OUTCOMES:
#   oxphos_abs  = VST nuclear-OXPHOS subunit composite (mtDNA-clean; the powered level)
#   oxphos_gsva = MitoCarta OXPHOS-family GSVA composite (secondary; the 19-set family in
#                 $defs includes MITOCARTA_OXPHOS_MT, so a small mtDNA contribution rides in
#                 -- kept for exact-match to script 29 Part D; oxphos_abs is the clean primary)
# CANDIDATE MYC-INDEPENDENT AXES (endogenous -- see ceiling):
#   dev_luminal   = LASP/LHS luminal state composite (developmental maturation axis)
#   tf_biogenesis = ESRRA/GABPA/NRF1 (ER/PGC1a mito-biogenesis TFs)
# GENERIC COVARIATE:
#   prolif        = Proliferation composite ("is the attenuation just proliferation?")
# SECONDARY AXIS:
#   mtdna_reprior = mtDNA-encoded OXPHOS mitoPPS (the reprioritisation shift)
fr <- data.frame(
  sample        = samples,
  group         = sm$group,
  tp            = sm$timepoint,
  myc           = sm$myc_status,
  oxphos_abs    = composite(vst_mat, defs$nuclear_oxphos_ens),
  oxphos_gsva   = comp_gsva(defs$mito_ox_gsva_sets),
  dev_luminal   = as.numeric(scale(comp_gsva(defs$luminal_sets))),
  tf_biogenesis = as.numeric(scale(comp_gsva(defs$tf_sets))),
  prolif        = as.numeric(scale(comp_gsva(prol_sets))),
  mtdna_reprior = as.numeric(scale(
    mp$mitopps_scores$`mtDNA-encoded OXPHOS subunits`[match(samples, mp$mitopps_scores$sample)])),
  stringsAsFactors = FALSE)
stopifnot(!anyNA(fr[, c("oxphos_abs","oxphos_gsva","dev_luminal","tf_biogenesis",
                        "prolif","mtdna_reprior")]))

outcomes  <- c("oxphos_abs", "oxphos_gsva")
# each axis-definition is a vector of frame columns entered as main + timepoint interaction
axis_defs <- list(
  dev_luminal   = "dev_luminal",
  tf_biogenesis = "tf_biogenesis",
  prolif        = "prolif",
  mtdna_reprior = "mtdna_reprior",
  dev_plus_tf   = c("dev_luminal", "tf_biogenesis"))

INT <- "tp12W:mycpos"        # the attenuation coefficient (genotype x time interaction)

# adjusted formula: y ~ tp*myc + <axes> + <tp:axes>  (axis moderated by timepoint)
adj_formula <- function(axes) {
  stats::as.formula(paste(
    "y ~ tp*myc +", paste(axes, collapse = " + "), "+",
    paste0("tp:", axes, collapse = " + ")))
}
b_int_of <- function(fit) {
  co <- stats::coef(fit); stopifnot(INT %in% names(co)); unname(co[INT])
}

# =============================================================================
# PART A: BASELINE ATTENUATION (the quantity every adjusted model must shrink)
# =============================================================================
baseline <- dplyr::bind_rows(lapply(outcomes, function(oc) {
  d  <- data.frame(y = fr[[oc]], tp = fr$tp, myc = fr$myc, grp = fr$group)
  m  <- summary(stats::lm(y ~ tp * myc, d))$coefficients
  gm <- tapply(d$y, d$grp, mean)
  tibble::tibble(
    outcome     = oc,
    b_int_base  = m[INT, "Estimate"],
    b_int_p     = m[INT, "Pr(>|t|)"],
    r2_base     = summary(stats::lm(y ~ tp * myc, d))$r.squared,
    m6_neg = unname(gm["6W_neg"]),  m6_pos  = unname(gm["6W_pos"]),
    m12_neg = unname(gm["12W_neg"]), m12_pos = unname(gm["12W_pos"]))
}))
# NOTE reproduced from Issue #4: on oxphos_abs this interaction is DIRECTIONAL not
# significant (b_int ~ -0.28, p ~ 0.61) -- this caps the test and makes the absorption
# FRACTION unstable; Delta-b_int + LRT are the primary read.

# =============================================================================
# PART B: COVARIATE ABSORPTION (the core test)
# =============================================================================
# For each outcome x candidate axis: refit with the axis (main + tp:axis) and ask whether
# b_int moves toward 0. absorption_frac = 1 - b_int_adj/b_int_base (toward 1 = fully
# absorbed; ~0 = untouched; <0 or >1 = grew / overshoot, common when the baseline is small).
# delta_b_int = b_int_adj - b_int_base is the STABLE primary metric (no division). LRT_F/p =
# nested anova of the added axis terms.
absorb_one <- function(oc, axname) {
  axes <- axis_defs[[axname]]
  d    <- cbind(data.frame(y = fr[[oc]], tp = fr$tp, myc = fr$myc), fr[axes])
  base <- stats::lm(y ~ tp * myc, d)
  adj  <- stats::lm(adj_formula(axes), d)
  b0   <- b_int_of(base); b1 <- b_int_of(adj)
  av   <- stats::anova(base, adj)               # nested F for the added axis terms
  ip   <- summary(adj)$coefficients[INT, "Pr(>|t|)"]
  tibble::tibble(
    outcome = oc, axis = axname,
    b_int_base = b0, b_int_adj = b1,
    delta_b_int = b1 - b0,                       # STABLE primary metric
    absorption_frac = 1 - b1 / b0,               # unstable when b0 small -- read with CI
    int_survives_p = ip,                         # does the interaction term still stand?
    delta_r2 = summary(adj)$r.squared - summary(base)$r.squared,
    lrt_F = av$F[2], lrt_p = av$`Pr(>F)`[2])
}
absorption <- dplyr::bind_rows(
  lapply(outcomes, function(oc)
    dplyr::bind_rows(lapply(names(axis_defs), function(ax) absorb_one(oc, ax)))))

# =============================================================================
# PART C: BOOTSTRAP CI on the absorption metrics (honest uncertainty at n=24)
# =============================================================================
# Case bootstrap, STRATIFIED by the 4 design groups (preserves n=6/group so tp*myc stays
# estimable). `mediation` is not installed -> manual via boot::boot. Reports CIs for
# delta_b_int (stable) AND absorption_frac (unstable) so the reader sees whether "fully
# explains" (frac ~ 1) is even inside the interval -- and how wide the fraction gets when
# the baseline interaction is near zero. ENDOGENEITY caveat applies to every number here.
grid <- expand.grid(outcome = outcomes, axis = names(axis_defs),
                    stringsAsFactors = FALSE)
metric_vec <- function(d) {
  unlist(lapply(seq_len(nrow(grid)), function(i) {
    oc <- grid$outcome[i]; axes <- axis_defs[[grid$axis[i]]]
    dd   <- cbind(data.frame(y = d[[oc]], tp = d$tp, myc = d$myc), d[axes])
    b0   <- b_int_of(stats::lm(y ~ tp * myc, dd))
    b1   <- b_int_of(stats::lm(adj_formula(axes), dd))
    stats::setNames(c(b1 - b0, 1 - b1 / b0),
                    paste(oc, grid$axis[i], c("delta", "frac"), sep = "__"))
  }))
}
set.seed(1)
boot_out <- boot::boot(fr, statistic = function(data, idx) metric_vec(data[idx, ]),
                       R = 2000, strata = fr$group)
# manual percentile CIs from the bootstrap replicate matrix (robust to the heavy tail)
perc_ci <- function(j) stats::quantile(boot_out$t[, j], c(0.025, 0.975),
                                       na.rm = TRUE, names = FALSE)
boot_ci <- tibble::tibble(
  key   = names(boot_out$t0),
  t0    = unname(boot_out$t0),
  lo    = vapply(seq_along(boot_out$t0), function(j) perc_ci(j)[1], numeric(1)),
  hi    = vapply(seq_along(boot_out$t0), function(j) perc_ci(j)[2], numeric(1))) |>
  tidyr::separate(key, into = c("outcome", "axis", "metric"), sep = "__")

# =============================================================================
# PART D: WITHIN-WT SLOPE ABSORPTION (the script-17 extension, OXPHOS analogue)
# =============================================================================
# Within WT (n=12) the transgene is fixed. Does the WT OXPHOS temporal slope VANISH when the
# developmental + TF axes are held? Script 17 found the biogenesis Myc+ decline did NOT
# shrink under the dev axis (grew -0.075 -> -0.143). This is the OXPHOS test of the same
# question. HYPOTHESIS-GENERATING (n=12, correlated covariates).
wt <- fr[fr$myc == "neg", ]
withinwt_one <- function(oc) {
  d    <- data.frame(y = wt[[oc]], tp = wt$tp,
                     dev = wt$dev_luminal, tf = wt$tf_biogenesis)
  base <- stats::lm(y ~ tp, d)
  adj  <- stats::lm(y ~ tp + dev + tf, d)
  s0   <- summary(base)$coefficients["tp12W", ]
  s1   <- summary(adj)$coefficients["tp12W", ]
  tibble::tibble(
    outcome = oc,
    slope_base = unname(s0["Estimate"]), slope_base_p = unname(s0["Pr(>|t|)"]),
    slope_adj  = unname(s1["Estimate"]), slope_adj_p  = unname(s1["Pr(>|t|)"]),
    slope_absorption_frac = 1 - unname(s1["Estimate"]) / unname(s0["Estimate"]))
}
withinwt <- dplyr::bind_rows(lapply(outcomes, withinwt_one))

# =============================================================================
# PART E: "WHAT WOULD SETTLE IT" (documented; not a statistical claim)
# =============================================================================
# These are the identifying experiments the bulk data cannot substitute for. Recorded so the
# manuscript states the ceiling honestly rather than over-reading a bootstrap fraction.
settle_it <- tibble::tribble(
  ~experiment,                       ~what_it_would_show,
  "Inducible Myc off/on (tet)",      "does the Myc+ OXPHOS arm re-DIVERGE on re-induction at 12W? separates a fading transgene effect from a fixed developmental ceiling",
  "Perturb ESRRA / NRF1 / GABPA",    "genetic necessity of the ER/PGC1a TF axis for the OXPHOS level -- turns co-variation into a causal test",
  "Single-cell / deconvolution",     "is the bulk change a CELL-COMPOSITION shift (luminal fraction) or a per-cell program change? the two are confounded at bulk resolution",
  "ATAC / ChIP at OXPHOS promoters", "is the TF axis actually bound / active at OXPHOS genes (Issue #3 left the non-MYC OXPHOS axis a candidate, not proof)")

# =============================================================================
# PART F: FIGURES
# =============================================================================
out_dir <- here::here("outputs", "attenuation_moderation")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ax_levels <- names(axis_defs)

# F1 -- b_int before vs after each adjustment (does it move to 0?). Dumbbell base->adj.
p1 <- absorption |>
  dplyr::mutate(axis = factor(axis, levels = ax_levels)) |>
  ggplot2::ggplot() +
  ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
  ggplot2::geom_segment(ggplot2::aes(x = b_int_base, xend = b_int_adj,
                                     y = axis, yend = axis), colour = "grey70",
                        arrow = grid::arrow(length = grid::unit(0.10, "cm"))) +
  ggplot2::geom_point(ggplot2::aes(x = b_int_base, y = axis, shape = "baseline"),
                      size = 2.6, colour = "#4575B4") +
  ggplot2::geom_point(ggplot2::aes(x = b_int_adj, y = axis, shape = "adjusted"),
                      size = 2.6, colour = "#D73027") +
  ggplot2::facet_wrap(~ outcome, scales = "free_x") +
  ggplot2::scale_shape_manual(values = c(baseline = 16, adjusted = 17)) +
  ggplot2::labs(
    title = "Part B: does conditioning on a MYC-independent axis shrink the attenuation?",
    subtitle = "b_int (genotype x time) before -> after adjustment. Toward 0 = absorbs; away = grows (echoes script 17)",
    x = "interaction coefficient b_int", y = NULL, shape = NULL) +
  ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "B_bint_before_after.pdf"), p1, width = 10, height = 4.5)

# F2 -- absorption fraction with bootstrap CI (ref 0 = no absorption, 1 = full). The wide
# CIs are the message: at a directional baseline "fully explains" is not distinguishable.
p2 <- boot_ci |>
  dplyr::filter(metric == "frac") |>
  dplyr::mutate(axis = factor(axis, levels = ax_levels)) |>
  ggplot2::ggplot(ggplot2::aes(x = t0, y = axis)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "solid",  colour = "grey40") +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi), orientation = "y",
                         width = 0.25, colour = "grey55") +
  ggplot2::geom_point(size = 2.6, colour = "#D73027") +
  ggplot2::facet_wrap(~ outcome) +
  ggplot2::coord_cartesian(xlim = c(-2, 2)) +
  ggplot2::labs(
    title = "Part C: absorption fraction with bootstrap 95% CI (0 = none, 1 = fully explains)",
    subtitle = "clipped to [-2,2]; wide/heavy-tailed because the baseline interaction is directional (int_p~0.61). ENDOGENOUS axes -> biased",
    x = "absorption fraction  (1 - b_int_adj / b_int_base)", y = NULL) +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "C_absorption_fraction_bootCI.pdf"), p2, width = 10, height = 4.5)

# F2b -- the STABLE companion: delta_b_int with bootstrap CI (no division). CI crossing 0 =
# the adjustment does not reliably move the attenuation.
p2b <- boot_ci |>
  dplyr::filter(metric == "delta") |>
  dplyr::mutate(axis = factor(axis, levels = ax_levels)) |>
  ggplot2::ggplot(ggplot2::aes(x = t0, y = axis)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi), orientation = "y",
                         width = 0.25, colour = "grey55") +
  ggplot2::geom_point(size = 2.6, colour = "#1A9850") +
  ggplot2::facet_wrap(~ outcome, scales = "free_x") +
  ggplot2::labs(
    title = "Part C (stable metric): change in the attenuation, Delta b_int, with bootstrap 95% CI",
    subtitle = "b_int_adj - b_int_base; no division by a small baseline. CI crossing 0 = adjustment does not reliably absorb",
    x = "Delta b_int (adjusted - baseline)", y = NULL) +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "C_delta_bint_bootCI.pdf"), p2b, width = 10, height = 4.5)

# F3 -- within-WT OXPHOS slope before/after adjusting for dev + TF (the script-17 analogue)
p3 <- withinwt |>
  tidyr::pivot_longer(c(slope_base, slope_adj), names_to = "model", values_to = "slope") |>
  dplyr::mutate(model = factor(ifelse(model == "slope_base", "WT slope (unadjusted)",
                                      "WT slope | dev + TF held"),
                               levels = c("WT slope (unadjusted)", "WT slope | dev + TF held"))) |>
  ggplot2::ggplot(ggplot2::aes(x = slope, y = outcome, fill = model)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
  ggplot2::scale_fill_manual(values = c("WT slope (unadjusted)" = "#4575B4",
                                        "WT slope | dev + TF held" = "#D73027")) +
  ggplot2::labs(
    title = "Part D (within WT, n=12): does the WT OXPHOS 6W->12W slope vanish when dev + TF are held?",
    subtitle = "script-17 analogue for OXPHOS. If the adjusted slope does not shrink toward 0, the axes do not absorb the WT trajectory",
    x = "WT temporal slope (6W -> 12W)", y = NULL, fill = NULL) +
  ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "D_withinwt_slope_absorption.pdf"), p3, width = 9, height = 4)

# =============================================================================
# PART G: SAVE
# =============================================================================
moderation <- list(
  baseline    = baseline,
  absorption  = absorption,
  boot_ci     = boot_ci,
  withinwt    = withinwt,
  settle_it   = settle_it,
  meta = list(
    boot_R          = boot_out$R,
    outcomes        = outcomes,
    axis_defs       = axis_defs,
    interaction_term = INT),
  notes = paste(
    "Issue #5: can MYC-independent factors EXPLAIN the OXPHOS attenuation? The attenuation is",
    "the genotype x time INTERACTION (b_int); 'explaining' it = does conditioning on a",
    "candidate axis SHRINK b_int? baseline = per-outcome b_int (oxphos_abs directional,",
    "int_p~0.61 -- reproduces script 29 absolute_stats). absorption = per axis: b_int_adj,",
    "delta_b_int (STABLE primary metric, no division), absorption_frac (=1-b1/b0, UNSTABLE at",
    "a small baseline -> read with boot_ci), int_survives_p, delta_r2, nested LRT (lrt_F/",
    "lrt_p). boot_ci = stratified (by the 4 groups) case bootstrap, R=2000, percentile CI for",
    "delta_b_int and absorption_frac (mediation pkg not installed -> manual boot::boot).",
    "withinwt = script-17 analogue: WT (n=12) OXPHOS 6W->12W slope before/after holding",
    "dev_luminal + tf_biogenesis. settle_it = the identifying experiments (inducible Myc,",
    "TF perturbation, single-cell, ATAC/ChIP). CEILING: candidate axes are ENDOGENOUS (Myc",
    "drives them) -> absorption estimates BIASED; n=6/group; collinear axes; baseline",
    "interaction ns. This BOUNDS the claim -- it cannot prove a MYC-independent cause.",
    "Expected (echoing script 17 biogenesis non-absorption): the axes do NOT fully absorb the",
    "attenuation. Do NOT read a large absorption_frac as 'development explains it' -- report",
    "the CI and the endogeneity caveat together. See docs/2026-07-08_BlockA_revision_plan.md."))
saveRDS(moderation, here::here("results", "attenuation_moderation.rds"))
message("Saved results/attenuation_moderation.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  am <- readRDS(here::here("results", "attenuation_moderation.rds"))

  # --- Part A: baseline attenuation to be shrunk (b_int; note oxphos_abs int_p ~0.61) ---
  am$baseline |> print()

  # --- Part B: does any axis absorb b_int? delta_b_int is the stable read; absorption_frac
  #     toward 1 = "fully explains" (but read WITH the CI); lrt_p = do the axis terms add?
  am$absorption |>
    dplyr::select(outcome, axis, b_int_base, b_int_adj, delta_b_int,
                  absorption_frac, int_survives_p, lrt_p) |> print(n = Inf)

  # --- Part C: bootstrap CIs. frac CIs should be WIDE (unstable baseline); delta CIs the
  #     honest primary -- if they straddle 0, no reliable absorption. ---
  am$boot_ci |> dplyr::filter(metric == "delta") |> print(n = Inf)
  am$boot_ci |> dplyr::filter(metric == "frac")  |> print(n = Inf)

  # --- Part D: within-WT OXPHOS slope -- does it vanish when dev + TF held? (script-17) ---
  am$withinwt |> print()

  # --- Part E: what would actually settle it ---
  am$settle_it |> print()

  list.files(here::here("outputs", "attenuation_moderation"), pattern = "\\.pdf$")
}
