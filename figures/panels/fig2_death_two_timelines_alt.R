# =============================================================================
# fig2_death_two_timelines_alt.R -- one BH3-only sensor is flat in development
# and lost under Myc, and it is PUMA
# -----------------------------------------------------------------------------
# SLOT: Fig. 2G (alt). An ALTERNATIVE to fig2_priming_ratios.R, not a
# replacement -- both are built and the author picks.
#
#   "... the PUMA/BCL-XL ratio showed a striking reversal ... This reduction was
#    HIGHLY SPECIFIC TO THE Bbc3 TRANSCRIPT, as other BH3-only sensors did not
#    show similar trends. Notably, Bbc3 mRNA levels remained constant in the WT
#    gland but decreased significantly within the 6>12W_myc context (-0.48,
#    padj = 0.016), demonstrating a specific interaction between genotype and
#    timeline effects."
#
# WHY A SECOND VERSION. fig2_priming_ratios.R draws the RATIOS against the global
# rescaling line, which is the first half of the sentence and is the right picture
# for it. It cannot show the second half -- the specificity -- because a ratio
# hides which of its two members moved, and it cannot show the third, that Bbc3 is
# flat in the wild-type gland, because the wild-type timeline is not on it.
#
# THE PLANE IS FIG. 2F (alt)'s, one level down: x = what the wild-type gland does
# across the window, y = what the Myc+ gland does, dashed diagonal = development
# alone. Read the four regions:
#
#   on the diagonal, upper right   up in both -- purely developmental. Bmf sits
#                                  here (+1.02 / +0.84, both significant): the
#                                  one death gene the window moves, and Myc does
#                                  nothing to it.
#   at the origin                  Bcl2l1 (-0.11 / +0.07, neither significant).
#                                  BCL-XL DOES NOT MOVE, so the ratio's reversal
#                                  is entirely its numerator.
#   below the line, below zero     lost under Myc and not in development. Bbc3
#                                  (+0.06 / -0.48, padj 0.016) is alone among the
#                                  BH3-only sensors.
#
# THE ONE HONEST QUALIFIER, and it is on the panel rather than buried: Bax also
# falls significantly on the Myc timeline (-0.38, padj 0.008). Bax is an EFFECTOR,
# not a BH3-only sensor, so "specific to Bbc3" holds as the sentence writes it --
# but the sentence should say "among the BH3-only sensors", and a reader can see
# Bax sitting there.
#
# Reads (read-only, no re-run):
#   results/priming_arm_teb.rds (script 42) -- $machinery, the curated death
#                                  roster with all four contrasts and padj
# Output: outputs/figures/panels/fig2_death_two_timelines_alt.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

pa_path <- here::here("results", "priming_arm_teb.rds")
require_fresher_than(pa_path)
m <- as.data.frame(readRDS(pa_path)$machinery)
stopifnot(all(c("gene", "arm", "lfc_wt_time", "padj_wt_time", "lfc_myc_time",
                "padj_myc_time", "lfc_interaction") %in% names(m)))

# =============================================================================
# the roster
# =============================================================================
# Script 42's `machinery` is a mixed table -- death machinery, biogenesis TFs and
# OXPHOS blot targets in one object -- so the death arm is selected by its own
# `arm` labels rather than by a list re-typed here. Everything that is a trigger,
# an effector, a brake or an execution step; nothing that is biogenesis or OXPHOS.
DEATH <- "BH3-only|effector|brake|execution|apoptosome|IAP"
d <- m[grepl(DEATH, m$arm), ]
d$sig <- !is.na(d$padj_myc_time) & d$padj_myc_time < 0.05
d$sensor <- grepl("BH3-only", d$arm)
stopifnot(nrow(d) >= 12L, "Bbc3" %in% d$gene, "Bcl2l1" %in% d$gene)

g <- function(x, col) d[[col]][d$gene == x]
# The claim, asserted so a re-run cannot flip it silently: Bbc3 is flat in the
# wild-type window and significantly down under Myc, and it is the only BH3-only
# sensor that is.
sens_sig <- d$gene[d$sensor & d$sig & d$lfc_myc_time < 0]
stopifnot(identical(sens_sig, "Bbc3"),
          g("Bbc3", "padj_wt_time") > 0.5, g("Bbc3", "padj_myc_time") < 0.05,
          g("Bcl2l1", "padj_wt_time") > 0.5, g("Bcl2l1", "padj_myc_time") > 0.5,
          # and Bmf is the developmental one: significant on BOTH timelines
          g("Bmf", "padj_wt_time") < 0.05, g("Bmf", "padj_myc_time") < 0.05)

# =============================================================================
# the panel
# =============================================================================
LIM <- c(-1, 1) * max(abs(c(d$lfc_wt_time, d$lfc_myc_time))) * 1.08

p <- ggplot2::ggplot(d, ggplot2::aes(lfc_wt_time, lfc_myc_time)) +
  two_timeline_base(LIM, diag_at = 0.62,
                    quadrant = "below the line = lost under Myc as well",
                    quadrant_at = c(0.50, 0.02), quadrant_hjust = 0.5) +
  ggplot2::geom_point(ggplot2::aes(shape = sig), size = 1.7, stroke = 0.35,
                      colour = "grey15", fill = "grey15") +
  ggrepel::geom_text_repel(ggplot2::aes(label = gene), size = 1.7,
                           colour = "grey15", fontface = "italic", seed = 6,
                           max.overlaps = Inf, min.segment.length = 0,
                           segment.size = 0.2, segment.colour = "grey60",
                           box.padding = 0.30, point.padding = 0.20) +
  ggplot2::scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1),
                              breaks = c(TRUE, FALSE),
                              labels = c("padj < 0.05 under Myc", "n.s."),
                              name = NULL) +
  ggplot2::labs(x = "wild-type 6>12W  (raw log2FC)", y = "Myc+ 6>12W") +
  ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 1.6))) +
  theme_panel(base_size = 6) +
  # Key inside, top left: nothing rises under Myc while falling in development,
  # so the wedge above the diagonal on the left cannot be occupied.
  ggplot2::theme(
    legend.position        = "inside",
    legend.position.inside = c(0.01, 0.99),
    legend.justification   = c(0, 1),
    legend.background      = ggplot2::element_blank(),
    legend.margin          = ggplot2::margin(0, 0, 0, 0),
    legend.key.size        = ggplot2::unit(2.4, "mm"),
    legend.spacing.y       = ggplot2::unit(0.3, "mm"),
    plot.margin            = ggplot2::margin(1.5, 2.5, 1, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
say <- function(x, nm = x)
  sprintf("%s %+.3f in development (padj %.2f) and %+.3f under Myc (padj %.3f)",
          nm, g(x, "lfc_wt_time"), g(x, "padj_wt_time"),
          g(x, "lfc_myc_time"), g(x, "padj_myc_time"))

LEGEND <- panel_legend(
  slot = "Fig. 2G (alt)",
  what = paste0(
    "The death machinery on both timelines at once: how each transcript changes ",
    "across the wild-type window on the horizontal axis, and across the same ",
    "window in the Myc+ gland on the vertical. The dashed diagonal is what a ",
    "gene would do under development alone. Filled points are significant on ",
    "the Myc timeline."),
  detail = c(
    sprintf("n = 6 per group; %d transcripts, selected from script 42's curated roster by its own arm labels (BH3-only triggers, effectors, brakes, execution steps). Raw (unshrunken) DESeq2 log2 fold changes, Benjamini-Hochberg adjusted genome-wide.",
            nrow(d)),
    sprintf("THE CLAIM, AND THE SCRIPT ASSERTS IT: %s -- and it is the ONLY BH3-only sensor that is flat in development and significantly down under Myc. The other sensors do not do it: %s; %s; %s; %s.",
            say("Bbc3"), say("Bcl2l11"), say("Bid"), say("Pmaip1"), say("Bik")),
    sprintf("BCL-XL DOES NOT MOVE ON EITHER TIMELINE: %s. So the reversal of the PUMA:Bcl-xL ratio that Fig. 2G reports is entirely its numerator, which is what makes the ratio worth quoting as a PUMA result rather than a balance result.",
            say("Bcl2l1")),
    sprintf("AND THE PANEL CARRIES ITS OWN NEGATIVE CONTROL: %s. Bmf is the one death transcript the window itself moves, and it moves in BOTH genotypes by nearly the same amount -- it sits on the diagonal. A gene on the diagonal is developmental; Bbc3 is as far off it as anything here.",
            say("Bmf")),
    sprintf("ONE HONEST QUALIFIER, VISIBLE ON THE PANEL: %s. Bax is an EFFECTOR, not a BH3-only sensor, so \"specific to Bbc3\" is true as the sentence writes it -- but the sentence should say \"among the BH3-only sensors\", because a reader can see Bax in the same quadrant.",
            say("Bax")),
    sprintf("The interaction terms rank the same way: Bbc3 %+.3f, Bax %+.3f, Bid %+.3f, Bak1 %+.3f, against Bcl2l1 %+.3f and Bmf %+.3f.",
            g("Bbc3", "lfc_interaction"), g("Bax", "lfc_interaction"),
            g("Bid", "lfc_interaction"), g("Bak1", "lfc_interaction"),
            g("Bcl2l1", "lfc_interaction"), g("Bmf", "lfc_interaction"))),
  bounds = c(
    "BATCH = TIMEPOINT, and it bites BOTH axes: each coordinate is a temporal contrast and is DESCRIBED, not claimed. What is batch-CLEAN is the vertical distance from the diagonal, which is the interaction, because genotype is balanced within each extraction batch. Read the panel down from the line, not along the axes.",
    "The adjusted p-values are genome-wide, so a gene that is significant here is significant against the whole transcriptome and not against this roster. Only three of the transcripts drawn reach padj < 0.05 on the Myc timeline.",
    "A TRANSCRIPT IS NOT A PRIMING STATE. Apoptotic priming is a property of the protein complement and of how close the mitochondrion sits to the threshold; BH3 profiling is the measurement, and this panel is a reason to do it rather than a substitute.",
    "Several of these transcripts are lowly expressed (Bid 148, Bbc3 152, Birc5 150 mean counts), so a real half-log2 effect could be missed in either timeline.",
    "This panel and fig2_priming_ratios.R are two readings of the same event, not two results: one draws the ratios against the global rescaling, this one draws the members against development."),
  source = c(
    "results/priming_arm_teb.rds (scripts/42_priming_arm_and_teb_substrate.R) -- $machinery, the curated death roster with all four contrasts and their adjusted p-values",
    "The grammar: two_timeline_base() in figures/panels/_panel_common.R, shared with Figs. 2F (alt), 2H+I (alt) and S2D (alt)"))

save_panel_p(p, "fig2_death_two_timelines_alt", height = 76)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the drawn roster, ranked by what the Myc timeline does
  d[order(d$lfc_myc_time),
    c("gene", "arm", "baseMean", "lfc_wt_time", "padj_wt_time", "lfc_myc_time",
      "padj_myc_time", "lfc_interaction")] |>
    print(row.names = FALSE, digits = 3)

  ## the BH3-only sensors alone -- the sentence's actual scope
  d[d$sensor, c("gene", "lfc_wt_time", "lfc_myc_time", "padj_myc_time")] |>
    print(row.names = FALSE, digits = 3)

  ## what the same roster does on the GENOTYPE axis, which is Fig. 2G's ruler
  d[, c("gene", "lfc_myc_6W", "padj_myc_6W", "lfc_myc_12W")] |>
    (\(x) x[order(-x$lfc_myc_6W), ])() |> print(row.names = FALSE, digits = 3)
}
