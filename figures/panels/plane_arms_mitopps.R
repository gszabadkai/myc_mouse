# =============================================================================
# plane_arms_mitopps.R -- the mitochondrial arms on the mitoPPS ruler, both
# timelines: every arm lies on the diagonal
# -----------------------------------------------------------------------------
# SLOT: not currently cited. Built 2026-09-21 for the Figure-1 verification pass
# (script 54); `plane_` keeps it outside the runners' `^fig` glob until the
# allocation note gives it a slot.
#
# WHY A SECOND ARMS PANEL (author's ruling 5, 2026-09-21: "report the diagonal
# drop on both rulers, draw mitoPPS"). mitoPPS is a pairwise ratio WITHIN the
# mitochondrial compartment, so a shift the whole compartment makes cancels, and
# what is left is how the compartment reallocates its budget. On this ruler the
# respiratory chain sits close to the diagonal; on the content ruler
# (plane_arms_content.R) it sits nearer the declared boundary. Both are drawn, so
# neither is shown alone.
#
# WHAT THIS RULER CANNOT CARRY, said on the page's behalf in the legend: mitoPPS
# exists only for MitoCarta pathways, so there is no proliferative or TEB arm here;
# no expression-matched null exists on it, so there are no bars; and its units are
# not a gene's log2 fold change, so it does NOT share axes with the four-gene panel.
# The declared 0.20 is applied here as the author declared it (ruling 5), in this
# ruler's own units.
#
# Reads (read-only, no re-run):
#   results/two_timeline_verification.rds (script 54) -- $arms_mitopps (script 40's
#       ruler, read and checked by 54), $arm_diagonal, $params
# Output: outputs/figures/panels/plane_arms_mitopps.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

tv_path <- here::here("results", "two_timeline_verification.rds")
require_fresher_than(tv_path)
tv <- readRDS(tv_path)

am <- as.data.frame(tv$arms_mitopps)
ad <- as.data.frame(tv$arm_diagonal)
DIAG <- tv$params$diagonal_threshold

d <- merge(am, ad[, c("arm", "content_int", "draw")], by = "arm")
stopifnot(nrow(d) == nrow(am), identical(DIAG, 0.20),
          # the vertical distance IS the interaction on this ruler too
          max(abs(d$p_myc - (d$p_wt + d$p_int))) < 1e-9)
d <- d[d$draw, ]
# This panel shares its axes with no other, so it takes its own frame: script 54's
# `plane_limits_mitopps` leaves 8% headroom, and the first render put the
# respiratory chain in the corner with its label printed over it. Asserted no
# smaller than 54's.
LIM <- c(-1, 1) * max(abs(c(d$p_wt, d$p_myc))) * 1.30
stopifnot(nrow(d) == 7L, !"OXPHOS (all)" %in% d$arm,
          LIM[2] >= max(tv$plane_limits_mitopps),
          all(abs(c(d$p_wt, d$p_myc)) < LIM[2]))
d$kind <- ifelse(d$arm == "OXPHOS subunits", "respiratory", "mitochondrial")

fx <- (d$p_wt - LIM[1]) / diff(LIM); fy <- (d$p_myc - LIM[1]) / diff(LIM)
DIAG_AT <- 0.79
stopifnot(!any(fx > DIAG_AT - 0.06 & fy > DIAG_AT - 0.06))   # the label's run is empty

p <- ggplot2::ggplot(d, ggplot2::aes(p_wt, p_myc)) +
  two_timeline_base(LIM, diag_at = DIAG_AT, band = DIAG) +
  ggplot2::geom_point(ggplot2::aes(fill = kind), shape = 21, size = 1.9,
                      stroke = 0.4, colour = "grey20") +
  ggrepel::geom_text_repel(ggplot2::aes(label = arm), size = 1.7, colour = "grey15",
                           seed = 12, max.overlaps = Inf, min.segment.length = 0,
                           segment.size = 0.2, segment.colour = "grey60",
                           box.padding = 0.35, point.padding = 0.2) +
  ggplot2::scale_fill_manual(values = c(respiratory = unname(pole_cols[["down"]]),
                                        mitochondrial = "grey62"), guide = "none") +
  ggplot2::labs(x = "wild-type 6>12W  (mitoPPS)", y = "Myc+ 6>12W  (mitoPPS)") +
  theme_panel(base_size = 6) +
  ggplot2::theme(plot.margin = ggplot2::margin(1.5, 2.5, 1, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
ox <- d[d$arm == "OXPHOS subunits", ]
far <- d[which.max(abs(d$p_int)), ]

LEGEND <- panel_legend(
  slot = "not currently cited",
  what = paste0(
    "The seven mitochondrial arms on the mitoPPS ruler, on both developmental ",
    "timelines: each arm's change in within-compartment priority across the wild-type ",
    "window (horizontal) and across the same window in the Myc+ gland (vertical). The ",
    "dashed diagonal is development alone and the dotted lines are the declared ",
    "threshold for lying on it."),
  detail = c(
    "n = 6 animals per group. mitoPPS (Monzel et al. 2025) from linear-scale DESeq2 normalised counts, script 08; the per-pathway contrasts are script 40's ruler, read and checked by script 54. The vertical distance from the diagonal is the interaction, exactly.",
    sprintf("THE RESPIRATORY CHAIN (OXPHOS subunits, filled dark): %+.3f across the wild-type window and %+.3f across the Myc+ one, %+.3f from the diagonal -- inside the declared %.2f by %.3f. The furthest of the seven is %s at %+.3f.",
            ox$p_wt, ox$p_myc, ox$p_int, DIAG, DIAG - abs(ox$p_int), far$arm, far$p_int),
    sprintf("ON THE CONTENT RULER the same arm is %+.3f from the diagonal (plane_arms_content.R): inside the declared %.2f by only %.3f, and beyond its expression-matched null. This ruler sits further from the declared boundary; the content ruler sits nearer it.",
            ox$content_int, DIAG, DIAG - abs(ox$content_int))),
  bounds = c(
    "BATCH = TIMEPOINT, on BOTH axes: each coordinate is a temporal contrast and is DESCRIBED, not claimed; the vertical distance from the diagonal is the batch-clean quantity. Read the panel down from the line, never along an axis.",
    "mitoPPS is RELATIVE: a pairwise ratio within the compartment, so a shift the whole compartment makes cancels. A point on the diagonal here means the arm kept its SHARE of the budget in both genotypes, not that its transcripts did not move.",
    "WHAT IS NOT HERE: no proliferative or TEB arm (neither is a MitoCarta pathway, and mitoPPS exists only for those); no expression-matched null (none exists on this ruler); and not the four-gene panel's axes (mitoPPS units are not a gene's log2 fold change). OXPHOS (all), the union of the two OXPHOS arms drawn, is not drawn.",
    "THE 0.20 was declared for log2 fold changes and is applied here as declared, in this ruler's own units (the author's ruling 5); it is not recalibrated for mitoPPS."),
  source = c(
    "results/two_timeline_verification.rds (scripts/54_two_timeline_verification.R) -- $arms_mitopps (script 40's ruler, via results/background_vs_myc.rds), $arm_diagonal, $plane_limits_mitopps",
    "The grammar: two_timeline_base() in figures/panels/_panel_common.R, with its declared band"))

save_panel_p(p, "plane_arms_mitopps", height = 84)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  d[order(d$p_int), c("arm", "pathway", "n_genes", "p_wt", "p_myc", "p_int",
                      "content_int")] |>
    print(row.names = FALSE, digits = 3)
}
