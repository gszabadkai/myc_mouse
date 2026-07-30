# =============================================================================
# fig1B_cell_state_composition.R -- MEC state composition and the TEB-ductal axis
# -----------------------------------------------------------------------------
# SUPPORTS (Results, "Myc drives early breast tumourigenesis by inducing OXPHOS
# and biosynthetic pathways", paragraph 1):
#   "GSVA scoring showed no major changes in the overall composition of major
#    mammary cell states (BMYO - basal-myoepithelial, LHS - luminal hormone
#    sensitive and luminal alveolar secretory precursor - LASP), but a clear
#    TEB-ductal shift away from the pubertal proliferative state at 12W"
#
# FORM (author, 2026-07-30): a square heatmap laid out AS THE DESIGN SCHEMATIC IS.
# Each of the four programmes gets its own 2x2 tile grid with age on x and
# genotype on y, exactly the geometry of Fig. S1B, so a trajectory is read the
# same way in both panels: left to right within a row is the 6->12W change,
# top to bottom within a column is the Myc effect. The previous version plotted
# 24 points against a group axis and the trajectories were hard to see; the
# per-sample distribution is retained in the sandbox at the end of this file.
#
# WHAT THE DATA SAY, AND WHERE THE SENTENCE NEEDS ATTENTION. In within-group SD
# units the wild-type 6->12W shift is BMYO +0.10, LASP -0.70, LHS -0.60, and
# TEB-ductal -0.91. So BMYO is flat, but the two LUMINAL states decline by two
# thirds of what the TEB axis does -- "no major changes in the overall
# composition" is not what the object says, and script 26 carries the same note
# in its own PART C2 header ("the data do NOT support a BMYO expansion (flat);
# the WT change is luminal (LASP/LHS) DECLINE"). Every one of those four shifts
# is non-significant at n=6 (p 0.24-0.86), so the honest statement is a shape,
# not a set of tests: one state flat, two down, the TEB axis furthest down.
#
# WHY EVERYTHING IS IN WITHIN-GROUP SD UNITS. The state composites are means of
# GSVA scores; TEB-ductal is a DIFFERENCE of two such means, so its raw spread is
# about 1.7x theirs (within_sd 0.36 vs 0.18-0.22) for arithmetic reasons alone.
# A shared raw fill scale would manufacture the contrast the sentence claims.
# Dividing each programme by its own pooled within-group SD is the project's
# existing convention for exactly this (scripts 26:335, 26:351-352, 27).
#
# Input:  results/dev_program_myc_integration.rds   (script 26 -- $annot, $state_profile,
#                                                    $state_stats)
#         results/myc_endogenous_amplification.rds  (script 27 -- $prog_wide, $prog_stats)
#         results/gsva_scores.rds                   (script 15 -- per-sample scores)
# Output: outputs/figures/panels/fig1B_cell_state_composition.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

dev_path  <- here::here("results", "dev_program_myc_integration.rds")
prog_path <- here::here("results", "myc_endogenous_amplification.rds")
require_fresher_than(dev_path)
require_fresher_than(prog_path)

dev   <- readRDS(dev_path)
prog  <- readRDS(prog_path)
gsva  <- readRDS(here::here("results", "gsva_scores.rds"))

stopifnot(all(c("annot", "state_profile", "state_stats") %in% names(dev)),
          all(c("prog_wide", "prog_stats") %in% names(prog)))

# --- per-sample state composites, rebuilt exactly as scripts/26:333-336 -------
# (26 saves only the group means; the panel needs the 24 values, so the composite
# is recomputed here from the same inputs and then ASSERTED against 26's own
# state_profile -- if the two ever diverge the panel stops.)
scores      <- gsva$scores
sample_meta <- as.data.frame(gsva$sample_meta)[colnames(scores), , drop = FALSE]

annot <- dev$annot
stopifnot(all(grep("^MG_", rownames(scores), value = TRUE) %in% annot$set))

main_sets <- annot |>
  dplyr::filter(as.character(state) %in% c("BMYO", "LASP", "LHS")) |>
  dplyr::select(set, state)

state_comp <- as.data.frame(scores[main_sets$set, , drop = FALSE]) |>
  tibble::rownames_to_column("set") |>
  tidyr::pivot_longer(-set, names_to = "sample", values_to = "gsva") |>
  dplyr::left_join(main_sets, by = "set") |>
  dplyr::group_by(sample, state) |>
  dplyr::summarise(value = mean(gsva), .groups = "drop") |>
  dplyr::mutate(group      = sample_meta[sample, "group"],
                timepoint  = sample_meta[sample, "timepoint"],
                myc_status = sample_meta[sample, "myc_status"],
                programme  = as.character(state)) |>
  dplyr::select(sample, group, timepoint, myc_status, programme, value)

# assertion 1: group means reproduce script 26's saved state_profile
chk <- state_comp |>
  dplyr::group_by(group, programme) |>
  dplyr::summarise(mean_gsva = mean(value), .groups = "drop") |>
  dplyr::inner_join(dplyr::transmute(dev$state_profile,
                                     group = as.character(group),
                                     programme = as.character(state),
                                     ref = mean_gsva),
                    by = c("group", "programme"))
stopifnot(nrow(chk) == 12, max(abs(chk$mean_gsva - chk$ref)) < 1e-10)

# --- the TEB-ductal axis, per sample, from script 27 --------------------------
teb <- prog$prog_wide |>
  dplyr::transmute(sample, group, timepoint, myc_status,
                   programme = "TEB - ductal", value = teb_ductal)

dat <- dplyr::bind_rows(state_comp, teb)

# --- standardise each programme by its own pooled within-group SD ------------
dat <- dat |>
  dplyr::group_by(programme) |>
  dplyr::mutate(wsd = wsd_of(value, group),
                z   = (value - mean(value)) / wsd) |>
  dplyr::ungroup()

# assertion 2: the recomputed SDs reproduce the saved ones
ref_wsd <- c(
  stats::setNames(dev$state_stats$within_sd, as.character(dev$state_stats$state)),
  "TEB - ductal" = prog$prog_stats$within_sd[prog$prog_stats$program == "teb_ductal"])
got_wsd <- dat |> dplyr::distinct(programme, wsd)
stopifnot(all(abs(got_wsd$wsd - ref_wsd[got_wsd$programme]) < 1e-8))

# --- the 2x2 grid, one per programme -----------------------------------------
prog_levels <- c("BMYO", "LASP", "LHS", "TEB - ductal")

cells <- dat |>
  dplyr::group_by(programme, timepoint, myc_status) |>
  dplyr::summarise(mean_z = mean(z), .groups = "drop") |>
  dplyr::mutate(
    programme = factor(programme, levels = prog_levels),
    age       = factor(as.character(timepoint), levels = c("6W", "12W"),
                       labels = c("6 weeks", "12 weeks")),
    # genotype descending so WT sits on top, as in Fig. S1B
    geno      = factor(as.character(myc_status), levels = c("pos", "neg"),
                       labels = c(unname(geno_labels[["pos"]]),
                                  unname(geno_labels[["neg"]]))))
stopifnot(nrow(cells) == 16)

# assertion 3: the four differences a reader takes off the tiles are the numbers
# script 26 and 27 report (d6/d12 for the states; group means for TEB-ductal)
gap_from_tiles <- cells |>
  dplyr::select(programme, age, geno, mean_z) |>
  tidyr::pivot_wider(names_from = geno, values_from = mean_z) |>
  dplyr::mutate(gap = `Myc+` - WT)
st6  <- gap_from_tiles$gap[gap_from_tiles$age == "6 weeks" &
                             gap_from_tiles$programme == "BMYO"]
stopifnot(abs(st6 - dev$state_stats$d6[as.character(dev$state_stats$state) == "BMYO"]) < 1e-8)

# symmetric fill limit, rounded up to a clean tick
LIM <- ceiling(max(abs(cells$mean_z)) * 10) / 10

p <- ggplot2::ggplot(cells, ggplot2::aes(x = age, y = geno, fill = mean_z)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
  # white ink only where the fill is genuinely dark; the PRGn ramp is still light
  # at two thirds of the limit, so the switch sits high
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.2f", mean_z),
                                  colour = abs(mean_z) > 0.82 * LIM),
                     size = 2.1, show.legend = FALSE) +
  ggplot2::facet_wrap(~ programme, nrow = 2) +
  heat_fill(LIM, name = "mean z\n(within-group SD)") +
  ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  ggplot2::scale_x_discrete(expand = c(0, 0)) +
  ggplot2::scale_y_discrete(expand = c(0, 0)) +
  ggplot2::labs(x = NULL, y = NULL) +
  theme_panel() +
  ggplot2::theme(
    axis.line       = ggplot2::element_blank(),
    axis.ticks      = ggplot2::element_blank(),
    axis.text       = ggplot2::element_text(colour = "black"),
    legend.position = "right",
    legend.key.width  = ggplot2::unit(2.4, "mm"),
    legend.key.height = ggplot2::unit(7, "mm"),
    panel.spacing   = ggplot2::unit(2.0, "mm"))

# --- the legend text (never drawn) -------------------------------------------
ss <- dev$state_stats
ps <- prog$prog_stats
f  <- function(x) sprintf("%+.2f", x)
teb_row <- ps[ps$program == "teb_ductal", ]

LEGEND <- panel_legend(
  slot = "Fig. 1B",
  what = paste0(
    "Group-mean scores for the three consensus mammary epithelial states (BMYO ",
    "basal-myoepithelial, LASP luminal alveolar secretory precursor, LHS luminal ",
    "hormone-sensing) and for the TEB-versus-ductal axis. Each programme is drawn ",
    "on the 2x2 design of Fig. S1B: age across, genotype down."),
  detail = c(
    "n = 6 animals per group, n = 24 total. Each tile is the mean of six animals; the number in the tile is that mean.",
    sprintf("State composites are the mean GSVA score over the sets assigned to that state by the curated annotation (BMYO %d sets, LASP %d, LHS %d of 179 developmental sets); TEB-ductal is mean(TEB_VS_DUCTAL UP) minus mean(DN).",
            sum(annot$state == "BMYO"), sum(annot$state == "LASP"),
            sum(annot$state == "LHS")),
    sprintf("Scores are centred on the 24-sample mean and divided by each programme's own pooled within-group SD (BMYO %.2f, LASP %.2f, LHS %.2f, TEB-ductal %.2f), so the four grids share one fill scale: TEB-ductal is a difference of two GSVA means and its raw spread is ~1.7x the others for arithmetic reasons alone.",
            ss$within_sd[ss$state == "BMYO"], ss$within_sd[ss$state == "LASP"],
            ss$within_sd[ss$state == "LHS"], teb_row$within_sd),
    sprintf("Wild-type 6 to 12 weeks (top row, left to right), in within-group SD: BMYO %s (p = %.2f), LASP %s (p = %.2f), LHS %s (p = %.2f), TEB-ductal %s (p = %.2f).",
            f(ss$wt_d[ss$state == "BMYO"]),  ss$wt_p[ss$state == "BMYO"],
            f(ss$wt_d[ss$state == "LASP"]),  ss$wt_p[ss$state == "LASP"],
            f(ss$wt_d[ss$state == "LHS"]),   ss$wt_p[ss$state == "LHS"],
            f(teb_row$wt_d), teb_row$wt_p),
    sprintf("Myc+ 6 to 12 weeks (bottom row): TEB-ductal %s (p = %.4f), the one temporal shift that reaches significance at this n.",
            f(teb_row$mycpos_shift / teb_row$within_sd), teb_row$mycpos_p),
    sprintf("Genotype effect (top minus bottom within a column), in within-group SD: BMYO %s at 6 weeks and %s at 12; LASP %s and %s; LHS %s and %s. BMYO is the only state with a powered pooled effect (%s, p = %.3f); LHS has none pooled (p = %.2f) but flips sign with age (interaction p = %.2f).",
            f(ss$d6[ss$state == "BMYO"]), f(ss$d12[ss$state == "BMYO"]),
            f(ss$d6[ss$state == "LASP"]), f(ss$d12[ss$state == "LASP"]),
            f(ss$d6[ss$state == "LHS"]),  f(ss$d12[ss$state == "LHS"]),
            f(ss$geno_beta[ss$state == "BMYO"] / ss$within_sd[ss$state == "BMYO"]),
            ss$geno_p[ss$state == "BMYO"], ss$geno_p[ss$state == "LHS"],
            ss$int_p[ss$state == "LHS"])),
  bounds = c(
    "Every 6-to-12-week statement here is exposed to batch = timepoint and is described, not claimed. The genotype comparisons are clean.",
    "No wild-type temporal shift is significant at n = 6 (p 0.24 to 0.86), so the panel shows a shape rather than a set of tests. Per-set consistency agrees with the shape: the fraction of sets moving down over the wild-type timeline is 0.46 for BMYO, 0.86 for LASP, 0.79 for LHS.",
    "The word 'clear' for the TEB shift is licensed by a different ruler, not by this panel: against expression-matched random gene sets the TEB-versus-ductal (HS) arm moves -0.422 at the 0th percentile (p_emp 0.0000), the largest mover in that comparison, and MG_HS_GRAY sits at the 0th percentile too (docs/2026-07-26_introduction_alignment_and_the_question.md section 0.7).",
    "A tile is a mean of six animals and the panel does not show their spread; the per-animal distributions are in the source object and in this script's sandbox block.",
    "The state composites are means over correlated gene sets (within-state r ~ 0.25-0.36, measured in script 17), so they are not independent replication and the p-values are indicative.",
    "GSVA is cohort-relative: a score is a position within these 24 samples, not an absolute abundance, and the fill is centred on the 24-sample mean by construction. These are transcriptional programme scores, not cell counts - a shift can be a change in state abundance or a change in programme activity within a state, and bulk RNA-seq cannot separate the two (see docs/deconvolution_subproject_plan.md)."),
  source = c(
    "results/dev_program_myc_integration.rds (scripts/26, PART C/C2) -- $annot, $state_profile, $state_stats",
    "results/myc_endogenous_amplification.rds (scripts/27) -- $prog_wide, $prog_stats",
    "results/gsva_scores.rds (scripts/15) -- per-sample scores; state annotation from data/dev_mec_annotation.csv"))

save_panel_p(p, "fig1B_cell_state_composition",
             width = fig_w[["single"]], height = 62)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the sixteen tile values
  cells |> dplyr::select(programme, age, geno, mean_z) |>
    tidyr::pivot_wider(names_from = c(geno, age), values_from = mean_z) |> print()

  ## the per-animal spread the tiles average over -- the previous form of this
  ## panel, kept because a mean of six is worth checking against its points
  ggplot2::ggplot(dplyr::mutate(dat, programme = factor(programme, levels = prog_levels),
                                group = factor(as.character(group),
                                               levels = rev(names(group_cols)))),
                  ggplot2::aes(x = z, y = group)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
    ggplot2::geom_point(ggplot2::aes(fill = group), shape = 21, size = 1.4,
                        colour = "grey25", stroke = 0.2) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 124, size = 2.4) +
    ggplot2::facet_wrap(~ programme, nrow = 1) +
    ggplot2::scale_fill_manual(values = group_cols, guide = "none") +
    ggplot2::scale_y_discrete(labels = group_labels) +
    ggplot2::labs(x = "GSVA composite (within-group SD)", y = NULL) +
    theme_panel()

  ## the shape the sentence is about: WT 6W -> 12W, in SD units
  dev$state_stats |> dplyr::select(state, wt_d, wt_p, wt_setfrac_down) |> print()
}
