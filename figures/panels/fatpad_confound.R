# =============================================================================
# fatpad_confound.R -- why the fat-pad series cannot carry a respiratory claim
# -----------------------------------------------------------------------------
# DISCUSSION PANEL, not a manuscript slot. The `fatpad_` prefix keeps it outside
# rebuild_panels.R and panels_to_pdf.R, both of which glob `^fig.*\.R$`: it takes
# no slot, does not enter the 28-panel count, and does not appear in
# paper/analysis_record.qmd. It exists to be shown to collaborators.
#
# THE ARGUMENT, IN THE ORDER THE PANEL DRAWS IT.
#
# (a) The respiratory score tracks the fat. Within the 12-week-to-tumour limb the
#     endpoint that states the reprioritisation claim -- nuclear OXPHOS subunits
#     relative to the mitoribosome -- correlates with adipocyte content at rho
#     +0.58. It is a tissue-composition axis wearing a mitochondrial label.
#
# (b) And that is not bad luck, it is the construction. The two halves of the
#     endpoint load on the adipocyte fraction with OPPOSITE SIGNS, so subtracting
#     one from the other ADDS their adipose components instead of cancelling them.
#     The difference is more adipose-loaded than either part it is made of. This is
#     the panel's point: the construction built to isolate the biology is the one
#     the tissue gradient most contaminates.
#
# (c) And it cannot be adjusted away. Proliferation genuinely rises across this
#     limb. Covarying out the four adipose markers destroys that rise as well,
#     because on this limb adipose depletion and tumour progression are the same
#     variable. An adjustment that removes the confound removes the biology, so
#     there is no version of this analysis that works.
#
# WHY THE 6-WEEK GROUPS ARE NOT DRAWN: they are not part of the limb. Composition
# moves 1.87x over 6W->12W against 1.20x across the limb, so the developmental
# step is confound-aligned and was excluded before any test was run (script 49).
#
# WHAT IS RECOMPUTED HERE AND WHY: panel (c) needs the per-sample ADJUSTED Mki67,
# which the object stores only as a summary. It is refitted from `per_sample`,
# which carries all four markers, and asserted against the stored `adj_control`
# tau so the panel cannot drift from the script that made the claim.
#
# Reads : results/fatpad_tumour_limb_trend.rds (script 49 PARTS C, D, G)
# Output: outputs/figures/panels/fatpad_confound.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

fp_path <- here::here("results", "fatpad_tumour_limb_trend.rds")
if (!file.exists(fp_path)) stop("run scripts/49_fatpad_tumour_limb_oxphos_trend.R first")
fp <- readRDS(fp_path)

ADJ <- c("Adipoq", "Cidec", "Fabp4", "Plin1")
ps  <- as.data.frame(fp$per_sample)
ps  <- ps[ps$in_limb, ]
ps$group <- factor(ps$group, levels = names(fatpad_group_cols))
stopifnot(nrow(ps) == 20L, !anyNA(ps$group))

gi <- as.integer(ps$group)
sp <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman"))

# ---------------------------------------------------------------------------
# (a) the score tracks the fat
# ---------------------------------------------------------------------------
rho_a <- sp(ps$ox_nuc_mtrib, ps$Adipoq)
stopifnot(abs(rho_a - fp$load_direction$rho_primary[fp$load_direction$marker == "Adipoq"]) < 1e-8)

# No fitted line. The statistic quoted is a RANK correlation, and a least-squares
# line through these points is carried by the single adipose-depleted animal at
# Adipoq 10.7 -- it would draw a precision the data do not have. The leave-one-out
# range is the honest robustness statement and is annotated instead.
pa <- ggplot2::ggplot(ps, ggplot2::aes(Adipoq, ox_nuc_mtrib)) +
  ggplot2::geom_point(ggplot2::aes(fill = group), shape = 21, size = 2,
                      colour = "white", stroke = 0.3) +
  ggplot2::annotate("text", x = min(ps$Adipoq), y = max(ps$ox_nuc_mtrib),
                    hjust = 0, vjust = 1, size = 2.5, colour = ms_diverging[["neg"]],
                    label = sprintf("rho = %+.2f\n(%.2f to %.2f\nleaving any one out)",
                                    rho_a, fp$loo_range$loo_min[2], fp$loo_range$loo_max[2])) +
  ggplot2::scale_fill_manual(values = fatpad_group_cols, labels = fatpad_group_labs,
                             name = NULL) +
  ggplot2::labs(x = "adipocyte content (Adipoq, log2)",
                y = "OXPHOS subunits\nminus mitoribosome") +
  theme_panel() +
  ggplot2::theme(legend.position = "top",
                 legend.margin = ggplot2::margin(0, 0, -3, 0))

# ---------------------------------------------------------------------------
# (b) the two halves load with opposite signs, so the difference amplifies
# ---------------------------------------------------------------------------
MEAS <- c(rho_ox_sub  = "nuclear\nOXPHOS",
          rho_mtrib   = "mito-\nribosome",
          rho_primary = "their difference\n= the endpoint")
ld <- as.data.frame(fp$load_direction)
ld <- ld[ld$marker %in% ADJ, c("marker", names(MEAS))]
lb <- stats::reshape(ld, direction = "long", varying = names(MEAS),
                     v.names = "rho", timevar = "measure", times = names(MEAS),
                     idvar = "marker")
lb$measure <- factor(MEAS[lb$measure], levels = MEAS)
lb$marker  <- factor(lb$marker, levels = ADJ)
lb$sign    <- ifelse(lb$rho >= 0, "up", "down")

pb <- ggplot2::ggplot(lb, ggplot2::aes(marker, rho, fill = sign)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey40") +
  ggplot2::geom_col(width = 0.68) +
  ggplot2::facet_wrap(~ measure, nrow = 1) +
  ggplot2::scale_fill_manual(values = direction_cols, guide = "none") +
  # the negative side needs a labelled tick or the reader cannot judge how far
  # the mitoribosome bars fall; the data stop at about -0.33, so the limit is set
  ggplot2::scale_y_continuous(labels = lab_signed, breaks = c(-0.4, 0, 0.4, 0.8),
                              limits = c(-0.45, 0.82)) +
  ggplot2::labs(x = NULL, y = "correlation with each\nadipocyte marker") +
  theme_panel() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 strip.text = ggplot2::element_text(size = 6.4, lineheight = 0.95),
                 panel.grid.major.x = ggplot2::element_blank())

# ---------------------------------------------------------------------------
# (c) the negative control: adjustment destroys a real trend
# ---------------------------------------------------------------------------
fit_adj  <- stats::lm(stats::as.formula(paste("Mki67 ~", paste(ADJ, collapse = " + "))), ps)
ps$Mki67_adj <- stats::resid(fit_adj)
tau_raw <- stats::cor(ps$Mki67, gi, method = "kendall")
tau_adj <- stats::cor(ps$Mki67_adj, gi, method = "kendall")
# the panel must not be able to drift from the script that made the claim
stopifnot(abs(tau_raw - fp$adj_control$tau_adj[1]) < 1e-8,
          abs(tau_adj - fp$adj_control$tau_adj[2]) < 1e-8)

STATE <- c(raw = "proliferation,\nas measured",
           adj = "proliferation, after\nadjusting for fat")
pc <- rbind(
  data.frame(group = ps$group, state = "raw",
             z = as.numeric(scale(ps$Mki67))),
  data.frame(group = ps$group, state = "adj",
             z = as.numeric(scale(ps$Mki67_adj))))
pc$state <- factor(STATE[pc$state], levels = STATE)
tau_lab <- data.frame(state = factor(STATE, levels = STATE),
                      tau = c(tau_raw, tau_adj),
                      col = c(ms_diverging[["pos"]], ms_diverging[["neg"]]))

pcp <- ggplot2::ggplot(pc, ggplot2::aes(group, z)) +
  ggplot2::stat_summary(fun = stats::median, geom = "crossbar",
                        width = 0.55, linewidth = 0.28, colour = "grey45") +
  ggplot2::geom_point(ggplot2::aes(fill = group), shape = 21, size = 1.8,
                      colour = "white", stroke = 0.3,
                      position = ggplot2::position_jitter(width = 0.11, height = 0,
                                                          seed = 1)) +
  ggplot2::geom_text(data = tau_lab, ggplot2::aes(label = sprintf("tau = %+.2f", tau),
                                                  colour = I(col)),
                     x = 0.55, y = max(pc$z) * 1.05, hjust = 0, vjust = 1, size = 2.5,
                     inherit.aes = FALSE) +
  ggplot2::facet_wrap(~ state, nrow = 1) +
  ggplot2::scale_fill_manual(values = fatpad_group_cols, guide = "none") +
  ggplot2::scale_x_discrete(labels = fatpad_group_labs) +
  ggplot2::labs(x = NULL, y = "Mki67 (z within panel)") +
  theme_panel() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 panel.grid.major.x = ggplot2::element_blank())

p <- patchwork::wrap_plots(pa, pb, pcp, nrow = 1, widths = c(1, 1.35, 1.1)) +
  patchwork::plot_annotation(tag_levels = "a") &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 8, face = "bold"))

LEGEND <- panel_legend(
  slot = "Discussion FP1",
  what = paste(
    "Why the whole-tissue fat-pad progression series cannot test whether respiratory",
    "priority recovers as tumours establish. (a) On the 12-week-to-tumour limb the",
    "endpoint tracks adipocyte content. (b) Its two halves load on that content with",
    "opposite signs, so their difference amplifies the tissue signal rather than",
    "cancelling it. (c) Adjusting the tissue signal away also removes a proliferation",
    "increase that is certainly real."),
  detail = c(
    sprintf("n = 20 (12WK_POS 5, SMALL_TUMOUR 6, LARGE_TUMOUR 9). Spearman rho in (a) and (b); Kendall tau-b against group order in (c). Scores are z-composites on log2(DESeq-normalised counts + 1)."),
    sprintf("(a) rho(endpoint, Adipoq) = %+.3f. The line is an unweighted least-squares fit, drawn to show the direction only.", rho_a),
    sprintf("(b) nuclear OXPHOS subunits load %+.2f to %+.2f on the four adipocyte markers; the mitoribosome loads %+.2f to %+.2f; the difference loads %+.2f to %+.2f. The four adipocyte markers explain %.0f%% of the endpoint's variance and %.0f%% of ox_rel's.",
            min(ld$rho_ox_sub), max(ld$rho_ox_sub), min(ld$rho_mtrib), max(ld$rho_mtrib),
            min(ld$rho_primary), max(ld$rho_primary),
            100 * fp$adj_joint$adj_r2[1], 100 * fp$adj_joint$adj_r2[2]),
    sprintf("(c) Mki67 rises across the limb at tau %+.3f (one-sided p = %.3f) and goes to tau %+.3f (p = %.3f) after covarying out Adipoq, Cidec, Fabp4 and Plin1 -- while those four explain only %.0f%% of its variance. Points are z-scored WITHIN each facet so the two states share an axis; the tau values are computed on the unscaled data.",
            tau_raw, fp$adj_control$p_one_sided[1], tau_adj, fp$adj_control$p_one_sided[2],
            100 * fp$adj_control$adj_r2[2]),
    sprintf("The primary trend itself is tau %+.3f [%.3f, %.3f], one-sided p = %.3f -- no rise, and uninterpretable in either direction because the confound pushes the same way.",
            fp$trend$tau[fp$trend$endpoint == "ox_nuc_mtrib"],
            fp$trend$tau_lo[fp$trend$endpoint == "ox_nuc_mtrib"],
            fp$trend$tau_hi[fp$trend$endpoint == "ox_nuc_mtrib"],
            fp$trend$p_one_sided[fp$trend$endpoint == "ox_nuc_mtrib"])),
  bounds = c(
    "WHOLE FAT PAD, not purified MECs. Adipose is the majority tissue throughout and only partly dilutes; this is the opposite situation to the 6W/12W MEC cohort, where an adipocyte signal is contamination.",
    "The 6-week groups are excluded and are not drawn: composition moves 1.87x over 6W->12W against 1.20x across this limb, so the developmental step is confound-aligned. That step is already clean in the MEC cohort.",
    "n = 20, four near-collinear covariates (VIF 41-64). Panel (c) is a demonstration that adjustment fails here, not an estimate of anything.",
    "Cohort-relative scores. Nothing on this panel may be compared numerically with the MEC cohort or the orthotopic series -- directions and orderings only.",
    "This panel shows why a respiratory claim is unavailable. It says nothing against the dataset's two contributions, which need no respiratory score: the 6-week genotype panel and the flat BCL2L1 progression."),
  source = c("results/fatpad_tumour_limb_trend.rds (script 49 PARTS C, D, G)",
             "docs/2026-09-07_fatpad_tumour_limb_oxphos_trend.md sections 3, 7b"))

save_panel_p(p, "fatpad_confound", width = fig_w[["double"]], height = 66)

if (FALSE) {
  print(p)
  ## the three numbers the panel is built on
  fp$load_direction |> print()
  fp$adj_control |> print()
  fp$adj_joint |> print()
  LEGEND |> print()
}
