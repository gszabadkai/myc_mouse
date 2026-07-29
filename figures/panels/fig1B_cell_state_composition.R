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
# WHAT THE DATA SAY, AND WHERE THE SENTENCE NEEDS ATTENTION. In within-group SD
# units the wild-type 6->12W shift is BMYO +0.10, LASP -0.70, LHS -0.60, and
# TEB-ductal -0.91. So BMYO is flat, but the two LUMINAL states decline by two
# thirds of what the TEB axis does -- "no major changes in the overall
# composition" is not what the object says, and script 26 carries the same note
# in its own PART C2 header ("the data do NOT support a BMYO expansion (flat);
# the WT change is luminal (LASP/LHS) DECLINE"). Every one of those four shifts
# is non-significant at n=6 (p 0.24-0.86), so the honest statement is a shape,
# not a set of tests: one state flat, two down, the TEB axis furthest down.
# The panel plots that shape and lets it be read; see LEGEND for the numbers and
# for the rulers on which the TEB shift does clear a null.
#
# WHY EVERYTHING IS IN WITHIN-GROUP SD UNITS. The state composites are means of
# GSVA scores; TEB-ductal is a DIFFERENCE of two such means, so its raw spread is
# about 1.7x theirs (within_sd 0.36 vs 0.18-0.22) for arithmetic reasons alone.
# A shared raw axis would manufacture the contrast the sentence claims. Dividing
# each programme by its own pooled within-group SD is the project's existing
# convention for exactly this (scripts 26:351-352, 27), and makes the comparison
# the sentence is actually making.
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
# (26 saves only the group means; the panel needs the 24 points, so the composite
# is recomputed here from the same inputs and then ASSERTED against 26's own
# state_profile -- if the two ever diverge the panel stops.)
scores      <- gsva$scores
sample_meta <- as.data.frame(gsva$sample_meta)[colnames(scores), , drop = FALSE]
dev_sets    <- grep("^MG_", rownames(scores), value = TRUE)

annot <- dev$annot
stopifnot(all(dev_sets %in% annot$set))

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
wsd_of <- function(x, g) sqrt(mean(tapply(x, g, stats::var)))   # scripts/26:335

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

prog_levels <- c("BMYO", "LASP", "LHS", "TEB - ductal")
dat <- dplyr::mutate(
  dat,
  programme = factor(programme, levels = prog_levels),
  group     = factor(as.character(group),
                     levels = rev(c("6W_neg", "6W_pos", "12W_neg", "12W_pos"))))

# --- panel -------------------------------------------------------------------
# ggbeeswarm where available (n=6 per row, so overplotting is real), plain jitter
# otherwise -- the same fallback figures/fig01_mito_content.R uses.
pts_layer <- if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
  ggbeeswarm::geom_quasirandom(ggplot2::aes(fill = group), orientation = "y",
                               width = 0.26, shape = 21, size = 1.25,
                               colour = "grey25", stroke = 0.18)
} else {
  ggplot2::geom_jitter(ggplot2::aes(fill = group), height = 0.16, width = 0,
                       shape = 21, size = 1.25, colour = "grey25", stroke = 0.18)
}

p <- ggplot2::ggplot(dat, ggplot2::aes(x = z, y = group)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
  pts_layer +
  ggplot2::stat_summary(fun = mean,
                        fun.min = function(v) mean(v) - stats::sd(v) / sqrt(length(v)),
                        fun.max = function(v) mean(v) + stats::sd(v) / sqrt(length(v)),
                        geom = "errorbar", width = 0, linewidth = 0.45,
                        colour = "grey15") +
  ggplot2::stat_summary(fun = mean, geom = "point", shape = 124, size = 2.4,
                        colour = "grey15") +
  ggplot2::facet_wrap(~ programme, nrow = 1) +
  ggplot2::scale_fill_manual(values = group_cols, labels = group_labels,
                             breaks = c("6W_neg", "6W_pos", "12W_neg", "12W_pos"),
                             name = NULL) +
  ggplot2::scale_y_discrete(labels = group_labels) +
  ggplot2::labs(x = "GSVA composite (within-group SD)", y = NULL) +
  theme_panel() +
  ggplot2::theme(legend.position = "none",
                 axis.line.y  = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 panel.spacing.x = ggplot2::unit(2.4, "mm"))

# --- the legend text (never drawn) -------------------------------------------
ss <- dev$state_stats
ps <- prog$prog_stats
f  <- function(x) sprintf("%+.2f", x)
teb_row <- ps[ps$program == "teb_ductal", ]

LEGEND <- panel_legend(
  slot = "Fig. 1B",
  what = paste0(
    "Per-sample GSVA composites for the three consensus mammary epithelial ",
    "states (BMYO basal-myoepithelial, LASP luminal alveolar secretory ",
    "precursor, LHS luminal hormone-sensing) and for the TEB-versus-ductal ",
    "axis, in all four groups."),
  detail = c(
    "n = 6 animals per group, n = 24 total; each point is one animal. Tick and bar are the group mean and its standard error.",
    sprintf("State composites are the mean GSVA score over the sets assigned to that state by the curated annotation (BMYO %d sets, LASP %d, LHS %d of 179 developmental sets); TEB-ductal is mean(TEB_VS_DUCTAL UP) minus mean(DN).",
            sum(annot$state == "BMYO"), sum(annot$state == "LASP"),
            sum(annot$state == "LHS")),
    sprintf("Scores are divided by each programme's own pooled within-group SD (BMYO %.2f, LASP %.2f, LHS %.2f, TEB-ductal %.2f) so the four are on one scale: TEB-ductal is a difference of two GSVA means and its raw spread is ~1.7x the others for arithmetic reasons alone.",
            ss$within_sd[ss$state == "BMYO"], ss$within_sd[ss$state == "LASP"],
            ss$within_sd[ss$state == "LHS"], teb_row$within_sd),
    sprintf("Wild-type 6 to 12 weeks, in within-group SD: BMYO %s (p = %.2f), LASP %s (p = %.2f), LHS %s (p = %.2f), TEB-ductal %s (p = %.2f).",
            f(ss$wt_d[ss$state == "BMYO"]),  ss$wt_p[ss$state == "BMYO"],
            f(ss$wt_d[ss$state == "LASP"]),  ss$wt_p[ss$state == "LASP"],
            f(ss$wt_d[ss$state == "LHS"]),   ss$wt_p[ss$state == "LHS"],
            f(teb_row$wt_d), teb_row$wt_p),
    sprintf("Myc+ 6 to 12 weeks: TEB-ductal %s (p = %.4f), the one temporal shift that reaches significance at this n.",
            f(teb_row$mycpos_shift / teb_row$within_sd), teb_row$mycpos_p),
    sprintf("Genotype main effect (Myc+ minus wild type, pooled over age): BMYO %s (p = %.3f) is the only state with a powered effect; LASP %s (p = %.2f) and LHS %s (p = %.2f) have none, and LHS instead shows a time-dependent flip (interaction p = %.2f; Myc+ is below wild type at 6 weeks and above it at 12).",
            f(ss$geno_beta[ss$state == "BMYO"] / ss$within_sd[ss$state == "BMYO"]),
            ss$geno_p[ss$state == "BMYO"],
            f(ss$geno_beta[ss$state == "LASP"] / ss$within_sd[ss$state == "LASP"]),
            ss$geno_p[ss$state == "LASP"],
            f(ss$geno_beta[ss$state == "LHS"] / ss$within_sd[ss$state == "LHS"]),
            ss$geno_p[ss$state == "LHS"], ss$int_p[ss$state == "LHS"])),
  bounds = c(
    "Every 6-to-12-week statement here is exposed to batch = timepoint and is described, not claimed. The genotype comparisons are clean.",
    "No wild-type temporal shift is significant at n = 6 (p 0.24 to 0.86), so the panel shows a shape rather than a set of tests. Per-set consistency agrees with the shape: the fraction of sets moving down over the wild-type timeline is 0.46 for BMYO, 0.86 for LASP, 0.79 for LHS.",
    "The word 'clear' for the TEB shift is licensed by a different ruler, not by this panel: against expression-matched random gene sets the TEB-versus-ductal (HS) arm moves -0.422 at the 0th percentile (p_emp 0.0000), the largest mover in that comparison, and MG_HS_GRAY sits at the 0th percentile too (docs/2026-07-26_introduction_alignment_and_the_question.md section 0.7).",
    "The state composites are means over correlated gene sets (within-state r ~ 0.25-0.36, measured in script 17), so they are not independent replication and the p-values are indicative.",
    "GSVA is cohort-relative: a score is a position within these 24 samples, not an absolute abundance. These are transcriptional programme scores, not cell counts - a shift can be a change in state abundance or a change in programme activity within a state, and bulk RNA-seq cannot separate the two (see docs/deconvolution_subproject_plan.md)."),
  source = c(
    "results/dev_program_myc_integration.rds (scripts/26, PART C/C2) -- $annot, $state_profile, $state_stats",
    "results/myc_endogenous_amplification.rds (scripts/27) -- $prog_wide, $prog_stats",
    "results/gsva_scores.rds (scripts/15) -- per-sample scores; state annotation from data/dev_mec_annotation.csv"))

save_panel_p(p, "fig1B_cell_state_composition",
             width = fig_w[["onehalf"]], height = 46)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the four programmes side by side, in SD units
  dat |>
    dplyr::group_by(programme, group) |>
    dplyr::summarise(mean_z = mean(z), .groups = "drop") |>
    tidyr::pivot_wider(names_from = group, values_from = mean_z) |>
    print()

  ## the shape the sentence is about: WT 6W -> 12W, in SD units
  dev$state_stats |> dplyr::select(state, wt_d, wt_p, wt_setfrac_down) |> print()
}
