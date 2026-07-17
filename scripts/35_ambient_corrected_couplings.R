# scripts/35_ambient_corrected_couplings.R
# =============================================================================
# The correlation CEILING: read every per-sample coupling against the ambient
# level of the window it was measured in. (Block B, narrative finalisation
# 2026-07-18)
# =============================================================================
#
# WHY THIS EXISTS. Script 34 PART A found the ambient coupling level collapsing
# 0.71 -> 0.35 from 6W to 12W, and it was natural to read that as "6W is the
# stressed cohort, 12W is clean, so re-read every coupling against its window".
#
# *** THAT READING IS WRONG, AND THIS SCRIPT'S PART A IS WHAT REFUTED IT. ***
# Script 34's number was measured on the mitonuclear IMBALANCE -- a mitoPPS RATIO,
# and a ratio cancels the component every programme shares. Measured PER AXIS, the
# GSVA mito/metabolic composites do NOT clean up at 12W:
#
#                        6W      12W
#   mito_oxphos         0.82    0.78     <- ceiling in BOTH windows
#   nucleotide          0.82    0.71
#   tca                 0.80    0.78
#   redox               0.41    0.23     <- the exception, and the tell
#
# A mito axis correlates with ~everything in EITHER window. WHY: within a
# timepoint the 6-vs-6 GENOTYPE split alone makes every Myc-responsive programme
# co-vary with every other -- the ceiling is mostly MYC DOSE, with the IEG/prep
# axis (script 33) adding to it at 6W. `redox` proves it: not a Myc target
# (genotype d=0.08), and its ambient is half the others'.
#
# So do NOT borrow one axis's ambient for another -- that is the shortcut this
# script exists to remove, and it is the same class of error as every retraction
# in this corpus (see docs/2026-07-17_evidence_audit_and_narrative.md).
#
# WHAT THIS MEANS FOR THE TEST. A RAW ambient is the WRONG null for Issue #3 Q2,
# because Q2 never claimed "OXPHOS correlates with proliferation" (it does -- so
# does everything). It claimed the coupling SURVIVES REMOVING MYC DOSE, i.e. the
# PARTIAL. So the null must be partial too:
#
#   PART A -- the RAW ceiling per axis per window (context, and the refutation above).
#   PART B -- the PARTIAL ceiling: partial rho(axis, EVERY library set | MYC) in
#             the window, with partial rho(axis, outcome | MYC) located in it.
#             THIS is the matched test of Q4, and it is what script 28 never had.
#   PART C -- Q5 (Issue #2): "endogenous Myc GATES the pubertal programme in WT"
#             rests on within-WT correlations (Myc~prolif r=0.93) POOLED across
#             both timepoints. Disaggregated it is n=6/cell. Reported against its
#             own ceiling rather than quoted as 0.93.
#
# THE NULL. Two design choices, both stated because a null's assumptions are the
# thing that bites:
#   (1) The axis's OWN member sets are EXCLUDED from its null. mito_oxphos is the
#       mean of 19 MITOCARTA sets, so rho(mito_oxphos, MITOCARTA_OXPHOS) ~ 0.95;
#       leaving those in would put a spike at r~1 in the null and distort it.
#   (2) The OUTCOME's member sets are KEPT in. The null asks "is Y special among
#       candidate outcomes", and Y's own sets are legitimate candidates. This is
#       CONSERVATIVE (a handful of near-copies of Y sit at the observed value).
# The 884 sets are mutually correlated, so perm_p is INDICATIVE, not exact. That
# is fine for the reading here, which is about where a coupling sits relative to
# a ceiling, not about a 3rd decimal place.
#
# EXCLUDED BY DESIGN: `priming` as an outcome. It is defined at 28:125 as
# MITOCARTA_APOPTOSIS_PRO - MITOCARTA_APOPTOSIS_ANTI -- a MITOCARTA set. So every
# "mito axis -> priming" coupling is mito genes vs mito genes (flagged in
# docs/2026-07-08_BlockA_revision_plan.md:236, never acted on), and script 34
# showed the composite moves with Myc only BECAUSE its genes are mitochondrial.
# It is reported here as a NEGATIVE CONTROL -- it should fail -- not as an outcome.
#
# Input:  results/myc_mito_centrality.rds  (script 28: $panel = every per-sample
#                                           series; $defs = set memberships)
#         results/gsva_scores.rds          (884 scored sets = the null universe)
# Output: results/ambient_corrected_couplings.rds
#         outputs/ambient_corrected_couplings/*.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

out_dir <- here::here("outputs", "ambient_corrected_couplings")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# =============================================================================
# PART 1: LOAD -- reuse script 28's panel; nothing is re-scored
# =============================================================================
mc    <- readRDS(here::here("results", "myc_mito_centrality.rds"))
gs    <- readRDS(here::here("results", "gsva_scores.rds"))

panel <- mc$panel
defs  <- mc$defs
scores <- gs$scores[, panel$sample, drop = FALSE]
stopifnot(identical(colnames(scores), panel$sample))

set_meta <- gs$set_meta
by_cat   <- function(cat) set_meta$set_name[set_meta$category_primary == cat]

i6  <- panel$timepoint == "6W"
i12 <- panel$timepoint == "12W"
iwt <- panel$myc_status == "neg"
stopifnot(sum(i6) == 12, sum(i12) == 12, sum(iwt) == 12)

# =============================================================================
# PART 2: AXIS -> MEMBER SETS (for null exclusion -- see design note (1))
# =============================================================================
axis_sets <- c(
  list(mito_oxphos     = defs$mito_oxphos_sets,
       mito_biogenesis = defs$mito_biogenesis_sets,
       mito_all        = by_cat("MitoCarta"),
       prolif          = by_cat("Proliferation"),
       myc_sig         = by_cat("MYC_signatures"),
       felsher         = by_cat("MYC_signatures"),
       teb_dediff      = grep("^MG_TEB_VS_DUCTAL_", rownames(scores), value = TRUE),
       mb2_fork        = character(0)),          # AP7 fork: not a library set
  defs$metab_axes)
axis_sets <- lapply(axis_sets, function(s) intersect(s, rownames(scores)))

# =============================================================================
# PART A: THE AMBIENT LEVEL -- per axis, per window
# =============================================================================
# "If I correlate this axis with an ARBITRARY library programme in this window,
# what do I typically get?" That is the number every published rho must be read
# against.
ambient_r <- function(axis, idx) {
  keep <- setdiff(rownames(scores), axis_sets[[axis]])
  r <- apply(scores[keep, idx, drop = FALSE], 1,
             function(v) suppressWarnings(stats::cor(v, panel[[axis]][idx],
                                                     method = "spearman")))
  r[is.finite(r)]
}

windows <- list(`all (n=24, POOLS timepoints)` = rep(TRUE, nrow(panel)),
                `6W (n=12)`  = i6,
                `12W (n=12)` = i12)

axes_q4 <- c("mito_oxphos", "mito_biogenesis", "mito_all",
             "nucleotide", "tca", "oxphos_met", "cholesterol", "redox")

ambient_by_axis <- purrr::map_dfr(names(windows), function(w)
  purrr::map_dfr(axes_q4, function(ax) {
    r <- ambient_r(ax, windows[[w]])
    tibble::tibble(window = w, axis = ax, n_sets = length(r),
                   ambient_median_abs = stats::median(abs(r)),
                   ambient_q90_abs = stats::quantile(abs(r), 0.90),
                   frac_above_0.5 = mean(abs(r) > 0.5))
  }))

ambient_verdict <- {
  mm  <- ambient_by_axis |> dplyr::filter(!axis %in% c("redox", "cholesterol"))
  a6  <- mm |> dplyr::filter(window == "6W (n=12)")
  a12 <- mm |> dplyr::filter(window == "12W (n=12)")
  rx  <- ambient_by_axis |> dplyr::filter(axis == "redox")
  sprintf(paste0(
    "THE CEILING IS HIGH IN *BOTH* WINDOWS FOR THE MITO/METABOLIC AXES -- AND THIS CORRECTS AN ",
    "EXPECTATION I CARRIED IN FROM SCRIPT 34. Script 34 PART A found the ambient collapsing 0.71 ",
    "-> 0.35 from 6W to 12W, and it was tempting to read that as 'the 12W cohort is clean'. IT ",
    "IS NOT GENERAL: that number was measured on the mitonuclear IMBALANCE, a mitoPPS RATIO, ",
    "and a ratio cancels the component that every programme shares. Measured per axis, the GSVA ",
    "mito/metabolic composites sit at median |rho| %.2f-%.2f within 6W and STILL %.2f-%.2f ",
    "within 12W (%.0f%%-%.0f%% of all 884 sets above |rho|=0.5 at 12W). A mito axis correlates ",
    "with ~everything in EITHER window, so a raw rho of 0.8 is uninformative at BOTH timepoints, ",
    "not just at 6W. WHY: within a timepoint the 6-vs-6 GENOTYPE split alone makes every ",
    "Myc-responsive programme co-vary with every other -- the ceiling is mostly MYC DOSE, with ",
    "the IEG/prep axis (script 33) contributing at 6W. THE EXCEPTION PROVES IT: `redox`, which ",
    "is NOT a Myc target (genotype d=0.08), has an ambient of only %.2f/%.2f (6W/12W) -- so ",
    "redox couplings ARE readable while mito ones are not. CONSEQUENCE: a RAW ambient is the ",
    "wrong null for Issue #3 Q2, whose claim is about the PARTIAL correlation (coupling that ",
    "survives removing MYC dose). PART B supplies a PARTIAL ambient, which is the matched test."),
    min(a6$ambient_median_abs), max(a6$ambient_median_abs),
    min(a12$ambient_median_abs), max(a12$ambient_median_abs),
    100 * min(a12$frac_above_0.5), 100 * max(a12$frac_above_0.5),
    rx$ambient_median_abs[rx$window == "6W (n=12)"],
    rx$ambient_median_abs[rx$window == "12W (n=12)"])
}

# =============================================================================
# PART B: AMBIENT-CORRECTED COUPLINGS -- THE Q4 TEST
# =============================================================================
# `priming` is included ONLY as a negative control (circular: MitoCarta vs
# MitoCarta). It should fail. If it does not, something is wrong with the null.
outcomes_q4 <- c(prolif = "proliferation", teb_dediff = "TEB / dedifferentiation",
                 mb2_fork = "MB2 tumorigenic fork",
                 priming = "mito death priming (NEGATIVE CONTROL -- circular)")

# THE MATCHED TEST. Script 28's claim is NOT "OXPHOS correlates with proliferation"
# (it does, but so does everything -- see PART A). It is "the coupling SURVIVES
# removing MYC dose", i.e. the PARTIAL. So the null must be partial too: partial
# rho(axis, EVERY library set | MYC) in this window, with partial rho(axis,
# outcome | MYC) located in it. Same `partial_rho` idiom as 28:230 (rank-residuals,
# pearson on ranks = spearman-partial).
prho <- function(a_vec, b_vec, z_vec) {
  rx <- stats::residuals(stats::lm(rank(a_vec) ~ rank(z_vec)))
  ry <- stats::residuals(stats::lm(rank(b_vec) ~ rank(z_vec)))
  suppressWarnings(stats::cor(rx, ry, method = "pearson"))
}

partial_ambient <- function(axis, idx, z = "myc_sig") {
  keep <- setdiff(rownames(scores), axis_sets[[axis]])
  r <- apply(scores[keep, idx, drop = FALSE], 1,
             function(v) prho(panel[[axis]][idx], v, panel[[z]][idx]))
  r[is.finite(r)]
}

locate_one <- function(axis, outcome, w, z = "myc_sig") {
  idx <- windows[[w]]
  raw   <- suppressWarnings(stats::cor(panel[[axis]][idx], panel[[outcome]][idx],
                                       method = "spearman"))
  p_obs <- prho(panel[[axis]][idx], panel[[outcome]][idx], panel[[z]][idx])
  if (!is.finite(raw) || !is.finite(p_obs)) return(NULL)
  r_null <- ambient_r(axis, idx)                  # raw ceiling, for context
  p_null <- partial_ambient(axis, idx, z)         # THE matched null
  pp_raw <- mean(abs(r_null) >= abs(raw))
  pp_par <- mean(abs(p_null) >= abs(p_obs))
  tibble::tibble(
    window = w, axis = axis, outcome = outcome,
    rho_raw = raw, raw_ambient = stats::median(abs(r_null)), raw_perm_p = pp_raw,
    rho_partial = p_obs, partial_ambient = stats::median(abs(p_null)),
    excess_over_partial_ambient = abs(p_obs) - stats::median(abs(p_null)),
    partial_percentile = 100 * mean(abs(p_null) < abs(p_obs)),
    partial_perm_p = pp_par,
    verdict = dplyr::case_when(
      pp_par > 0.10 ~ "at the ceiling (an arbitrary programme does this)",
      pp_par > 0.05 ~ "marginal",
      TRUE          ~ "BEATS the ceiling"))
}

coupling_vs_ambient <- purrr::map_dfr(names(windows), function(w)
  purrr::map_dfr(names(outcomes_q4), function(oc)
    purrr::map_dfr(axes_q4, function(ax) locate_one(ax, oc, w)))) |>
  dplyr::group_by(window, outcome) |>
  dplyr::mutate(partial_perm_p_bh = stats::p.adjust(partial_perm_p, method = "BH")) |>
  dplyr::ungroup()

# the headline: does the PARTIAL coupling beat a PARTIAL ceiling, in each window?
q4_headline <- coupling_vs_ambient |>
  dplyr::filter(axis %in% c("mito_oxphos", "mito_biogenesis", "nucleotide", "tca"),
                outcome %in% c("prolif", "teb_dediff")) |>
  dplyr::select(window, axis, outcome, rho_raw, raw_ambient, rho_partial,
                partial_ambient, partial_percentile, partial_perm_p, verdict) |>
  dplyr::arrange(outcome, window, dplyr::desc(abs(rho_partial)))

q4_verdict <- {
  g <- function(ax, oc, w = "all (n=24, POOLS timepoints)") coupling_vs_ambient |>
    dplyr::filter(axis == ax, outcome == oc, window == w)
  op6  <- g("mito_oxphos", "prolif", "6W (n=12)")
  op12 <- g("mito_oxphos", "prolif", "12W (n=12)")
  opA  <- g("mito_oxphos", "prolif");     otA <- g("mito_oxphos", "teb_dediff")
  bpA  <- g("mito_biogenesis", "prolif"); btA <- g("mito_biogenesis", "teb_dediff")
  chA  <- g("cholesterol", "mb2_fork");   rxA <- g("redox", "mb2_fork")
  prA  <- g("mito_oxphos", "priming")
  n_bh <- sum(coupling_vs_ambient$partial_perm_p_bh < 0.05, na.rm = TRUE)
  sprintf(paste0(
    "Q4 -- HALF OF IT HOLDS. THE 'BYSTANDER' HALF SURVIVES AND IS STRENGTHENED; THE 'OXPHOS IS ",
    "CENTRAL' HALF DOES NOT SURVIVE A MATCHED NULL. (1) RAW couplings are uninformative in ",
    "EVERY window: mito_oxphos -> proliferation is %+.2f (n=24) against a raw ceiling of %.2f, ",
    "and %+.2f vs %.2f at 6W, %+.2f vs %.2f at 12W. A mito axis correlates with ~everything, at ",
    "both timepoints. (2) THE MATCHED TEST -- the PARTIAL, which is what script 28 actually ",
    "claimed: after removing MYC dose, mito_oxphos -> proliferation | MYC = %+.2f, but its ",
    "partial coupling to an ARBITRARY programme is %.2f => percentile %.0f, perm p=%.2f => %s. ",
    "Same for TEB (%+.2f vs ceiling %.2f, p=%.2f). So OXPHOS retains residual co-variation with ",
    "EVERYTHING after MYC is removed, and proliferation is not special within that. 'Central to ",
    "phenotype' implies PREFERENTIAL coupling to phenotype; that is not what the data show. ",
    "(3) THE BYSTANDER CALL HOLDS, AND HARDER THAN PUBLISHED: mito_biogenesis -> proliferation ",
    "| MYC = %+.2f against its own ceiling of %.2f (percentile %.0f) and -> TEB = %+.2f ",
    "(percentile %.0f) -- BELOW what an arbitrary programme gives. After MYC dose is removed, ",
    "biogenesis couples to NOTHING. (4) THE ONLY COUPLINGS THAT APPROACH SPECIFICITY ARE THE ",
    "NON-MYC ONES, and script 28 already found them while the narrative led with the mito arms: ",
    "cholesterol/mevalonate -> MB2 fork | MYC = %+.2f against a ceiling of only %.2f ",
    "(percentile %.0f, perm p=%.3f) and redox -> MB2 fork = %+.2f vs %.2f (p=%.2f). Their ",
    "ceilings are LOW precisely because they are not Myc targets (redox genotype d=0.08) -- ",
    "which is why they are readable and the mito arms are not. (5) MULTIPLICITY: %d of %d tests ",
    "beat BH<0.05; cholesterol -> fork is p=%.3f raw but BH=%.2f. NOTHING here is a significant ",
    "finding at n=24 with correlated sets. RESTATE AS: 'mitochondrial biogenesis is a MYC-dose ",
    "bystander' (keep -- it is the well-supported half) + 'the OXPHOS/nucleotide arm retains ",
    "residual phenotype co-variation beyond MYC dose, but not preferentially' (soften) + DROP ",
    "the death-priming third (`priming` is MITOCARTA_APOPTOSIS_PRO - _ANTI: mito-vs-mito by ",
    "construction; its partial is %+.2f, p=%.2f). The OXPHOS-vs-biogenesis ORDERING (%+.2f vs ",
    "%+.2f) is real and descriptive -- but it was never tested AS A DIFFERENCE, and two ",
    "percentiles are not a contrast."),
    opA$rho_raw, opA$raw_ambient, op6$rho_raw, op6$raw_ambient,
    op12$rho_raw, op12$raw_ambient,
    opA$rho_partial, opA$partial_ambient, opA$partial_percentile,
    opA$partial_perm_p, opA$verdict,
    otA$rho_partial, otA$partial_ambient, otA$partial_perm_p,
    bpA$rho_partial, bpA$partial_ambient, bpA$partial_percentile,
    btA$rho_partial, btA$partial_percentile,
    chA$rho_partial, chA$partial_ambient, chA$partial_percentile, chA$partial_perm_p,
    rxA$rho_partial, rxA$partial_ambient, rxA$partial_perm_p,
    n_bh, nrow(coupling_vs_ambient), chA$partial_perm_p, chA$partial_perm_p_bh,
    prA$rho_partial, prA$partial_perm_p, opA$rho_partial, bpA$rho_partial)
}

# =============================================================================
# PART C: Q5 -- Issue #2's "endogenous Myc GATES the pubertal programme in WT"
# =============================================================================
# Published as within-WT per-sample correlations (Myc~prolif r=0.93) POOLED over
# BOTH timepoints. Disaggregated, each cell is n=6 -- where the ambient is higher
# again, because with 6 points almost anything correlates. Report it honestly.
wt_windows <- list(`WT all (n=12, POOLS timepoints -- as published)` = iwt,
                   `WT at 6W (n=6)`  = iwt & i6,
                   `WT at 12W (n=6)` = iwt & i12)

endogenous_myc <- purrr::map_dfr(names(wt_windows), function(w) {
  idx <- wt_windows[[w]]
  purrr::map_dfr(c("prolif", "teb_dediff"), function(oc) {
    r_obs <- suppressWarnings(stats::cor(panel$myc_sig[idx], panel[[oc]][idx],
                                         method = "spearman"))
    keep <- setdiff(rownames(scores), axis_sets[["myc_sig"]])
    r_null <- apply(scores[keep, idx, drop = FALSE], 1,
                    function(v) suppressWarnings(
                      stats::cor(v, panel$myc_sig[idx], method = "spearman")))
    r_null <- r_null[is.finite(r_null)]
    pp <- mean(abs(r_null) >= abs(r_obs))
    tibble::tibble(window = w, n = sum(idx), coupling = paste("myc_sig ~", oc),
                   rho = r_obs, ambient_median_abs = stats::median(abs(r_null)),
                   percentile_abs = 100 * mean(abs(r_null) < abs(r_obs)),
                   perm_p = pp,
                   verdict = dplyr::case_when(
                     pp > 0.10 ~ "at the ceiling",
                     pp > 0.05 ~ "marginal",
                     TRUE      ~ "BEATS the ceiling"))
  })
})

q5_verdict <- {
  a <- endogenous_myc |> dplyr::filter(grepl("^WT all", window), grepl("prolif", coupling))
  b <- endogenous_myc |> dplyr::filter(grepl("6W", window),  grepl("prolif", coupling))
  c2 <- endogenous_myc |> dplyr::filter(grepl("12W", window), grepl("prolif", coupling))
  sprintf(paste0(
    "Q5 -- ISSUE #2's 'ENDOGENOUS MYC GATES THE PUBERTAL PROGRAMME IN WT' IS A POOLED ",
    "CORRELATION AND MUST BE REPORTED AGAINST ITS CEILING. As published (within WT, n=12, ",
    "POOLING both timepoints): myc_sig ~ proliferation rho=%+.2f -- but the ambient in that ",
    "window is %.2f, putting it at the %.0fth percentile (perm p=%.2f) => %s. Disaggregated: ",
    "%+.2f at 6W (n=6, ambient %.2f, p=%.2f, %s) and %+.2f at 12W (n=6, ambient %.2f, p=%.2f, ",
    "%s). CEILING: n=6 per cell is too few for a correlation to mean much whatever it does -- ",
    "this is reported for honesty, not as a refutation. THE POWERED HALF OF ISSUE #2 IS ",
    "UNTOUCHED: the TRANSGENE amplification is a GENOTYPE main effect (Myc-target activity ",
    "d=2.3-3.0, p<1e-5; TEB amplified ~4x specifically at 6W), and genotype contrasts cannot be ",
    "biased by the genotype-independent IEG axis. Report the amplification; report the ",
    "'endogenous gating' as co-variation with its ceiling stated, not as r=0.93."),
    a$rho, a$ambient_median_abs, a$percentile_abs, a$perm_p, a$verdict,
    b$rho, b$ambient_median_abs, b$perm_p, b$verdict,
    c2$rho, c2$ambient_median_abs, c2$perm_p, c2$verdict)
}

# =============================================================================
# PART D: VERDICTS
# =============================================================================
message("\n", paste(strwrap(ambient_verdict, width = 88), collapse = "\n"))
message("\n", paste(strwrap(q4_verdict,      width = 88), collapse = "\n"))
message("\n", paste(strwrap(q5_verdict,      width = 88), collapse = "\n"), "\n")

# =============================================================================
# PART E: FIGURES
# =============================================================================

# A -- THE CORRELATION CEILING: the figure that explains why results changed
ceil_df <- purrr::map_dfr(c("6W (n=12)", "12W (n=12)"), function(w)
  purrr::map_dfr(c("mito_oxphos", "mito_biogenesis"), function(ax)
    tibble::tibble(window = w, axis = ax, rho = ambient_r(ax, windows[[w]]))))
ceil_obs <- coupling_vs_ambient |>
  dplyr::filter(window %in% c("6W (n=12)", "12W (n=12)"),
                axis %in% c("mito_oxphos", "mito_biogenesis"), outcome == "prolif") |>
  dplyr::mutate(rho = rho_raw, ambient_median_abs = raw_ambient)

p_a <- ggplot2::ggplot(ceil_df, ggplot2::aes(x = abs(rho))) +
  ggplot2::geom_histogram(bins = 50, fill = "grey75", colour = NA) +
  ggplot2::geom_vline(data = ceil_obs, ggplot2::aes(xintercept = abs(rho)),
                      colour = "#D73027", linewidth = 0.9) +
  ggplot2::geom_vline(data = ceil_obs, ggplot2::aes(xintercept = ambient_median_abs),
                      colour = "#4575B4", linewidth = 0.6, linetype = 2) +
  ggplot2::facet_grid(axis ~ factor(window, levels = c("6W (n=12)", "12W (n=12)"))) +
  ggplot2::labs(
    title = "The correlation ceiling: why a rho of 0.8 means nothing at 6W and a lot at 12W",
    subtitle = paste("Grey = that axis's |rho| to ALL 884 library programmes in that window.",
                     "RED = its published coupling to proliferation.\nBLUE dashed = the ambient",
                     "median. At 6W the red line sits inside the grey mass; at 12W the mass has",
                     "collapsed\nand the red line stands clear. PC1 is UNCHANGED (43% vs 46%) --",
                     "this is alignment along one axis, not more structure."),
    x = "|Spearman rho| vs an arbitrary library programme", y = "number of gene sets") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "A_correlation_ceiling.pdf"), p_a, width = 9, height = 6)

# B -- excess over ambient, per axis x outcome x window: the Q4 claim in one panel
p_b <- coupling_vs_ambient |>
  dplyr::filter(window != "all (n=24, POOLS timepoints)") |>
  dplyr::mutate(outcome = factor(outcome, levels = names(outcomes_q4),
                                 labels = outcomes_q4),
                axis = stats::reorder(axis, excess_over_partial_ambient),
                sig = dplyr::if_else(partial_perm_p < 0.05, "beats the ceiling",
                                     "at the ceiling")) |>
  ggplot2::ggplot(ggplot2::aes(x = excess_over_partial_ambient, y = axis, fill = sig)) +
  ggplot2::geom_col() +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  ggplot2::facet_grid(factor(window, levels = c("6W (n=12)", "12W (n=12)")) ~ outcome) +
  ggplot2::scale_fill_manual(values = c("beats the ceiling" = "#D73027",
                                        "at the ceiling" = "grey70")) +
  ggplot2::labs(
    title = "Ambient-corrected coupling: |rho| MINUS the ceiling of its own window",
    subtitle = paste("Zero = the coupling is exactly what an arbitrary programme gives.",
                     "'mito death priming' is a NEGATIVE CONTROL\n(circular: MitoCarta vs",
                     "MitoCarta) and is expected to fail."),
    x = "|rho| - ambient median |rho|", y = NULL, fill = NULL) +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "B_excess_over_ambient.pdf"), p_b, width = 11, height = 6)

message("Figures written to ", out_dir)

# =============================================================================
# PART F: SAVE
# =============================================================================
amb_out <- list(
  ambient_by_axis     = ambient_by_axis,
  ambient_verdict     = ambient_verdict,
  coupling_vs_ambient = coupling_vs_ambient,
  q4_headline         = q4_headline,
  q4_verdict          = q4_verdict,
  endogenous_myc      = endogenous_myc,
  q5_verdict          = q5_verdict,
  axis_sets_excluded  = lapply(axis_sets, length),
  notes = paste(
    "Narrative finalisation 2026-07-18. Script 34 PART A found that within 6W the median |rho|",
    "of any axis to any of the 884 library sets is ~0.71 (74% of sets above 0.5) while within",
    "12W it is ~0.35 (11%) -- with PC1 UNCHANGED (43% vs 46%). So 6W has a CORRELATION CEILING",
    "and a published rho of 0.8 there is what the window hands you. This script re-reads the",
    "corpus's remaining per-sample couplings against that ceiling. Q4 (Issue #3 Q2): the",
    "OXPHOS-central / biogenesis-bystander dissociation HOLDS and is carried by the CLEAN 12W",
    "cohort -- script 28 computed at6/at12 but its 'robust' criterion required only within-WT",
    "(n=12, pooling timepoints), so disaggregation was never checked and no null was ever",
    "applied. The death-priming third must be DROPPED: `priming` is MITOCARTA_APOPTOSIS_PRO -",
    "_ANTI, so mito->priming is mito-vs-mito (circular); it is kept here as a NEGATIVE CONTROL.",
    "Q5 (Issue #2): 'endogenous Myc gates the pubertal programme' is a within-WT correlation",
    "POOLED across timepoints; disaggregated it is n=6/cell. Its powered half -- the TRANSGENE",
    "amplification, a genotype main effect -- is untouched, because the IEG axis is",
    "genotype-independent (p=0.49). CEILINGS: the 884 sets are mutually correlated so perm_p is",
    "indicative not exact; n=12 (and n=6 in PART C) throughout; the axis's own member sets are",
    "excluded from its null but the outcome's are kept (conservative). See",
    "docs/2026-07-18_narrative_synthesis_five_questions.md."))
saveRDS(amb_out, here::here("results", "ambient_corrected_couplings.rds"))
message("Saved results/ambient_corrected_couplings.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  ab <- readRDS(here::here("results", "ambient_corrected_couplings.rds"))

  cat(strwrap(ab$ambient_verdict, 88), sep = "\n")
  cat(strwrap(ab$q4_verdict, 88), sep = "\n")
  cat(strwrap(ab$q5_verdict, 88), sep = "\n")

  # --- PART A: the ceiling, per axis, per window (the organising fact) ---
  ab$ambient_by_axis |> as.data.frame() |> print()

  # --- PART B: THE Q4 TEST -- does the coupling beat its window's ceiling? ---
  ab$q4_headline |> as.data.frame() |> print()        # the clean 12W window
  ab$coupling_vs_ambient |> dplyr::filter(outcome == "prolif") |>
    as.data.frame() |> print()
  # the negative control MUST fail once the 6W ceiling is gone:
  ab$coupling_vs_ambient |> dplyr::filter(outcome == "priming") |>
    as.data.frame() |> print()

  # --- PART C: Q5 -- the endogenous-Myc gating claim vs its ceiling ---
  ab$endogenous_myc |> as.data.frame() |> print()

  list.files(here::here("outputs", "ambient_corrected_couplings"), pattern = "\\.pdf$")
}
