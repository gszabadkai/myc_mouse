# =============================================================================
# plane_four_genes.R -- four transcripts on the two-timeline plane: two lost under
# MYC, one that did not follow the wild-type rise, one on the diagonal
# -----------------------------------------------------------------------------
# SLOT: not currently cited. Built 2026-09-21 for the Figure-1 verification pass
# (script 54); `plane_` keeps it outside the runners' `^fig` glob until the
# allocation note gives it a slot.
#
# A THREE-WAY CONTRAST, NOT A TWO-WAY ONE (author's ruling, 2026-09-21). The draft
# read this plane as Bbc3 off the diagonal against Bax and Bcl2l1 on it. Check 1,
# read against the threshold declared before the numbers were retrieved, says
# otherwise: Bax is below the diagonal beyond the declared magnitude, with Bbc3's
# sign, and is MYC-specific by the declared rule. The positions this panel draws:
#
#   lost under MYC                    Bbc3, Bax -- displaced DOWN the Myc+ axis
#   did not follow the wild-type rise Foxo3 -- displaced ALONG the wild-type axis;
#                                     MYC-specific by the same rule, by the other
#                                     route (it rose in the normal gland and did
#                                     not rise under MYC; the Myc+ arm is flat)
#   on the diagonal                   Bcl2l1
#
# The classes are DERIVED from the declared rule and the direction of
# displacement, and then asserted against the ruling -- so a change in the object
# stops the panel rather than redrawing it silently. Bax is not labelled as a
# control: it is not one.
#
# THE LICENCE IS NOT THE SAME FOR THE FOUR (ruling 2). Bbc3 was pre-specified;
# Bcl2l1 is the fixed denominator of the pre-specified pair; Bax and Foxo3 are
# exploratory. Ruling 2 governs what the manuscript TEXT may claim, and it quotes
# no p-value for either. The legend SHOWS their interaction p, labelled
# exploratory (author's correction, 2026-09-21): pre-specification governs what
# may be claimed, not what may be shown, and a value withheld from a legend cannot
# be calibrated. For these two the legend prints the interaction -- the contrast
# the rule reads -- and not the arms' p-values.
#
# SAME AXES AS plane_arms_content.R, by construction (plane_content_lim()), so the
# genes can be read against the arms: the respiratory chain inside the declared
# band there, Bbc3 and Bax outside it here.
#
# Reads (read-only, no re-run):
#   results/two_timeline_verification.rds (script 54) -- $four_genes, $rules,
#       $licence, $params
# Output: outputs/figures/panels/plane_four_genes.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

tv_path <- here::here("results", "two_timeline_verification.rds")
require_fresher_than(tv_path)
tv <- readRDS(tv_path)

fg   <- as.data.frame(tv$four_genes)
DIAG <- tv$params$diagonal_threshold
stopifnot(nrow(fg) == 4L, setequal(fg$gene, c("Bbc3", "Foxo3", "Bax", "Bcl2l1")),
          identical(DIAG, 0.20),
          # the plane's geometry, gene by gene
          max(abs(fg$identity_dev)) < 1e-12)

LOST <- "lost under MYC"; NORISE <- "did not follow the wild-type rise"; ON <- "on the diagonal"
fg$cls <- ifelse(fg$on_diagonal, ON,
          ifelse(fg$myc_specific & abs(fg$myc_lfc) >= abs(fg$wt_lfc), LOST,
          ifelse(fg$myc_specific, NORISE, "off the diagonal, not MYC-specific")))
cls_of <- function(g) fg$cls[fg$gene == g]
# the ruling, asserted: a re-run cannot quietly turn this back into a two-way panel
stopifnot(cls_of("Bbc3") == LOST, cls_of("Bax") == LOST,
          cls_of("Bcl2l1") == ON, cls_of("Foxo3") == NORISE)

LIM <- plane_content_lim(tv)
fx <- (fg$wt_lfc - LIM[1]) / diff(LIM); fy <- (fg$myc_lfc - LIM[1]) / diff(LIM)
DIAG_AT <- 0.79
stopifnot(all(abs(c(fg$wt_lfc, fg$myc_lfc)) < LIM[2]),
          !any(fx > DIAG_AT - 0.06 & fy > DIAG_AT - 0.06))   # the label's run is empty

fill_of <- c(unname(pole_cols[["down"]]), unname(pole_cols[["up"]]),
             unname(pole_cols[["other"]]))
names(fill_of) <- c(LOST, NORISE, ON)

p <- ggplot2::ggplot(fg, ggplot2::aes(wt_lfc, myc_lfc)) +
  two_timeline_base(LIM, diag_at = DIAG_AT, band = DIAG) +
  ggplot2::geom_point(ggplot2::aes(fill = cls), shape = 21, size = 2.3,
                      stroke = 0.45, colour = "grey15") +
  ggrepel::geom_text_repel(ggplot2::aes(label = gene), size = 2.0,
                           fontface = "italic", colour = "grey10", seed = 4,
                           max.overlaps = Inf, min.segment.length = 0,
                           segment.size = 0.2, segment.colour = "grey60",
                           box.padding = 0.45, point.padding = 0.25) +
  ggplot2::scale_fill_manual(values = fill_of, guide = "none") +
  ggplot2::labs(x = "wild-type 6>12W  (log2FC)", y = "Myc+ 6>12W  (log2FC)") +
  theme_panel(base_size = 6) +
  ggplot2::theme(plot.margin = ggplot2::margin(1.5, 2.5, 1, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
g <- function(x) fg[fg$gene == x, ]
genes_in <- function(k) paste(fg$gene[fg$cls == k], collapse = " and ")
# p-values only where the licence allows them (ruling 2)
with_p <- function(x) {
  r <- g(x)
  sprintf("%s %+.3f across the wild-type window and %+.3f across the Myc+ one (SE %.3f, raw p %.4f); interaction %+.3f (SE %.3f, raw p %.4f)",
          x, r$wt_lfc, r$myc_lfc, r$myc_se, r$myc_p, r$int_lfc, r$int_se, r$int_p)
}
with_int_p <- function(x) {
  r <- g(x)
  sprintf("%s %+.3f across the wild-type window and %+.3f across the Myc+ one; interaction %+.3f (SE %.3f, raw p %.4f)",
          x, r$wt_lfc, r$myc_lfc, r$int_lfc, r$int_se, r$int_p)
}
# Bax is off the diagonal by the magnitude half of the rule alone; the legend says
# so, so it is asserted
ALPHA <- tv$params$alpha
stopifnot(identical(ALPHA, 0.05), g("Bax")$int_p > ALPHA, abs(g("Bax")$int_lfc) >= DIAG)

LEGEND <- panel_legend(
  slot = "not currently cited",
  what = paste0(
    "Four transcripts on both developmental timelines, on the same axes as the ",
    "gene-set arms: the change across the wild-type window (horizontal) and across ",
    "the same window in the Myc+ gland (vertical). The dashed diagonal is development ",
    "alone and the dotted lines are the declared threshold for lying on it. Dark: lost ",
    sprintf("under MYC (%s). Green: did not follow the wild-type rise (%s). Grey: on ",
            genes_in(LOST), genes_in(NORISE)),
    sprintf("the diagonal (%s).", genes_in(ON))),
  detail = c(
    sprintf("n = 6 animals per group. Raw (unshrunken) DESeq2 log2 fold changes from the genotype x timepoint fit. The vertical distance from the diagonal is the interaction term, exactly: Myc+ = wild type + interaction holds to %.0e for these four and to %.0e across all %d genes.",
            max(abs(fg$identity_dev)), max(tv$identity$max_abs_dev), tv$identity$n_genes[1]),
    sprintf("THE RULE, fixed before the numbers were retrieved: on the diagonal if %s; MYC-specific if %s.",
            tv$rules[["on_diagonal"]], tv$rules[["myc_specific"]]),
    sprintf("PRE-SPECIFIED: %s.", with_p("Bbc3")),
    sprintf("THE DENOMINATOR OF THE PRE-SPECIFIED PAIR, on the diagonal on both halves of the rule: %s.",
            with_p("Bcl2l1")),
    sprintf("EXPLORATORY, NOT PRE-SPECIFIED: %s -- beyond the declared magnitude, with Bbc3's sign, and off the diagonal by the magnitude half of the rule alone (its p is above %.2f); and %s -- the same rule, by the other route: it rises in the normal gland and is flat under MYC. Neither was named in advance, so both p-values are reported as exploratory.",
            with_int_p("Bax"), ALPHA, with_int_p("Foxo3"))),
  bounds = c(
    "BATCH = TIMEPOINT, on BOTH axes: each coordinate is a temporal contrast and is DESCRIBED, not claimed. What is batch-clean is the vertical distance from the diagonal, because genotype is balanced within each extraction batch. Read the panel down from the line, never along an axis.",
    "THE LICENCE, GENE BY GENE. Bbc3 was named in advance from the PGC1a westerns; Bcl2l1 is the fixed denominator of that pair. Bax and Foxo3 were not pre-specified: their p-values are shown and reported as exploratory, and the manuscript text quotes none for either.",
    "THE DOTTED BAND IS HALF THE RULE: the declared magnitude. The rule's other half, a raw interaction p above 0.05, is in the detail above rather than on the page.",
    "A transcript's position is not a protein level and not a measurement of how close a cell sits to the apoptotic threshold."),
  source = c(
    "results/two_timeline_verification.rds (scripts/54_two_timeline_verification.R) -- $four_genes (from results/interaction_results.rds, script 03), $rules, $licence",
    "The grammar: two_timeline_base() in figures/panels/_panel_common.R, with its declared band; the frame: plane_content_lim(), shared with plane_arms_content.R"))

save_panel_p(p, "plane_four_genes", height = 84)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  fg[, c("gene", "cls", "wt_lfc", "myc_lfc", "int_lfc", "int_se", "int_p",
         "on_diagonal", "myc_specific", "licence")] |>
    print(row.names = FALSE, digits = 3)
}
