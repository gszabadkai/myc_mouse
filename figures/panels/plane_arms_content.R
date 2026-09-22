# =============================================================================
# plane_arms_content.R -- the gene-set arms on both timelines, each against the
# vertical distance from the diagonal that comparably expressed genes reach
# -----------------------------------------------------------------------------
# SLOT: not currently cited. Built 2026-09-21 for the Figure-1 verification pass
# (script 54); the slot is assigned in the allocation note, and the `plane_`
# prefix keeps it outside the runners' `^fig` glob until then -- the same device
# the `biogax_` discussion panels use.
#
# THE PLANE (two_timeline_base, the figure's shared grammar): x is what the
# wild-type gland does across the window, y is what the Myc+ gland does, and the
# dashed diagonal is development alone. The vertical distance from it is the
# interaction -- exact, not approximate (script 54 asserts it gene by gene). The
# dotted lines are the DECLARED threshold for "on the diagonal", +/- 0.20, fixed
# before the numbers were read.
#
# THE NULL, per arm. Each grey bar is the 95% range of the vertical distance from
# the diagonal that 2,000 expression-matched random sets of the arm's own size
# reach, placed at the arm's own position. A point outside its bar has moved off
# the diagonal further than comparably expressed genes do. Script 43 nulled the
# wild-type axis only; the bar is the new null, and it is the batch-clean one.
#
# WHAT THE PANEL SHOWS, stated before anyone reads it off the page. The
# respiratory chain (OXPHOS subunits) is inside the declared 0.20 -- and outside
# its own bar. Both are true, both are drawn, and the legend says which rule each
# one is. The same arm on the mitoPPS ruler is plane_arms_mitopps.R; on this
# ruler it sits nearer the declared boundary (author's ruling 5).
#
# SAME AXES AS plane_four_genes.R, by construction: plane_content_lim() in
# _panel_common.R computes one frame for both.
#
# Reads (read-only, no re-run):
#   results/two_timeline_verification.rds (script 54) -- $arms_content,
#       $arm_null, $arm_diagonal, $params
# Output: outputs/figures/panels/plane_arms_content.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

tv_path <- here::here("results", "two_timeline_verification.rds")
ss_path <- here::here("results", "substrate_specificity_tradeoff.rds")
require_fresher_than(tv_path)
require_fresher_than(ss_path)
tv <- readRDS(tv_path)
w43 <- as.data.frame(readRDS(ss_path)$wt_null)   # script 43's wild-type null, the text's

ac <- as.data.frame(tv$arms_content)
an <- as.data.frame(tv$arm_null)
ad <- as.data.frame(tv$arm_diagonal)
DIAG <- tv$params$diagonal_threshold

d <- merge(ac, an[, c("arm", "n_matched", "pct_wt", "pct_myc", "pct_int",
                      "null_med_int", "null_int_lo", "null_int_hi")], by = "arm")
d <- merge(d, ad[, c("arm", "mitopps_int", "content_on_diagonal", "draw")], by = "arm")
stopifnot(nrow(d) == nrow(ac),
          # the vertical distance IS the interaction
          max(abs(d$c_int - (d$c_myc - d$c_wt))) < 1e-9,
          identical(DIAG, 0.20))
d <- d[d$draw, ]
stopifnot(nrow(d) == 9L, !"OXPHOS (all)" %in% d$arm)

# three kinds of arm, and the only encoding is the fill: the respiratory chain the
# sentence is about, the other mitochondrial arms, and the two arms that are not
# mitochondrial at all (drawn open, so a reader cannot take them for compartment)
NONMITO <- c("TEB vs ductal (HS)", "PROLIF_* pooled")
d$kind <- ifelse(d$arm == "OXPHOS subunits", "respiratory",
          ifelse(d$arm %in% NONMITO, "non-mitochondrial", "mitochondrial"))
stopifnot(sum(d$kind == "respiratory") == 1L, sum(d$kind == "non-mitochondrial") == 2L)

# the null bar: the 95% null range of the diagonal distance, at the arm's own x
d$bar_lo <- d$c_wt + d$null_int_lo
d$bar_hi <- d$c_wt + d$null_int_hi
d$outside_bar <- d$c_int < d$null_int_lo | d$c_int > d$null_int_hi

LIM <- plane_content_lim(tv)
stopifnot(all(abs(c(d$c_wt, d$c_myc, d$bar_lo, d$bar_hi)) < LIM[2]))

# Furniture placed where nothing is, checked in panel fractions: the diagonal's
# label on the upper line, and the quadrant note in the bottom-right corner.
# The null bars are furniture too: the first render put the label across two of
# them, so the check covers bars as well as points.
fx  <- (d$c_wt - LIM[1]) / diff(LIM); fy <- (d$c_myc - LIM[1]) / diff(LIM)
bhi <- (d$bar_hi - LIM[1]) / diff(LIM)
DIAG_AT <- 0.79
stopifnot(!any(fx > DIAG_AT - 0.06 & fy  > DIAG_AT - 0.06),   # points clear of the label
          !any(fx > DIAG_AT - 0.06 & bhi > DIAG_AT - 0.06),   # bars clear of the label
          !any(fx > 0.60 & fy < 0.12))                          # the quadrant note's strip

LAB <- c(`TEB vs ductal (HS)` = "TEB vs ductal", `PROLIF_* pooled` = "proliferation")
d$label <- ifelse(d$arm %in% names(LAB), LAB[d$arm], d$arm)

p <- ggplot2::ggplot(d, ggplot2::aes(c_wt, c_myc)) +
  two_timeline_base(LIM, diag_at = DIAG_AT, band = DIAG,
                    quadrant = "below the line = lost under MYC",
                    quadrant_at = c(0.99, 0.02), quadrant_hjust = 1) +
  ggplot2::geom_segment(ggplot2::aes(x = c_wt, xend = c_wt, y = bar_lo, yend = bar_hi),
                        linewidth = 1.1, colour = "grey86", lineend = "butt") +
  ggplot2::geom_point(ggplot2::aes(fill = kind), shape = 21, size = 1.9,
                      stroke = 0.4, colour = "grey20") +
  ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 1.7,
                           colour = "grey15", seed = 21, max.overlaps = Inf,
                           min.segment.length = 0, segment.size = 0.2,
                           segment.colour = "grey60", box.padding = 0.35,
                           point.padding = 0.2) +
  ggplot2::scale_fill_manual(values = c(respiratory = unname(pole_cols[["down"]]),
                                        mitochondrial = "grey62",
                                        `non-mitochondrial` = "white"),
                             guide = "none") +
  ggplot2::labs(x = "wild-type 6>12W  (set average log2FC)",
                y = "Myc+ 6>12W  (set average log2FC)") +
  theme_panel(base_size = 6) +
  ggplot2::theme(plot.margin = ggplot2::margin(1.5, 2.5, 1, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
row_of <- function(a) d[d$arm == a, ]
ox <- row_of("OXPHOS subunits"); pf <- row_of("PROLIF_* pooled"); tb <- row_of("TEB vs ductal (HS)")
off <- d[!d$content_on_diagonal & d$kind != "non-mitochondrial", ]
off <- off[order(off$c_int), ]
mito <- d[d$kind != "non-mitochondrial", ]
out_bar <- mito$label[mito$outside_bar]
pct43 <- function(a) w43$percentile[w43$arm == a]
# the legend says the arms beyond 0.20 are the ones MYC raised most at six weeks;
# asserted, so a change in the object cannot leave that sentence behind
stopifnot(setequal(off$arm, head(mito$arm[order(-mito$c_myc_6W)], nrow(off))),
          # and that the matched sets themselves sit below the diagonal
          all(d$null_med_int < 0), !anyNA(mito$mitopps_int))

LEGEND <- panel_legend(
  slot = "not currently cited",
  what = paste0(
    "Nine gene-set arms on both developmental timelines at once: each arm's average ",
    "log2 fold change across the wild-type window (horizontal) and across the same ",
    "window in the Myc+ gland (vertical). The dashed diagonal is development alone ",
    "and the dotted lines are the declared threshold for lying on it. Each grey bar ",
    "is the range of vertical distances from the diagonal that expression-matched ",
    "random sets of the arm's size reach."),
  detail = c(
    sprintf("n = 6 animals per group. Raw (unshrunken) DESeq2 log2 fold changes, averaged with equal weight over each arm's genes -- the content ruler, unweighted. Arms and membership are script 43's (%d to %d genes), resolved through the symbol reconciler and reproduced against script 43 to 1e-9. OXPHOS (all), the union of the two OXPHOS arms drawn, is not drawn.",
            min(d$n_genes), max(d$n_genes)),
    sprintf("THE BARS are the 2.5th to 97.5th percentiles of the interaction over 2,000 random sets matched to each arm on expression (20 baseMean bins), placed at the arm's own horizontal position. The same draws give the wild-type and Myc+ percentiles below, so an arm's three percentiles describe one set of draws."),
    sprintf("THE RESPIRATORY CHAIN (OXPHOS subunits, filled dark): %+.3f in the wild-type window and %+.3f in the Myc+ one, %+.3f below the diagonal. That is INSIDE the declared %.2f, by %.3f -- and OUTSIDE its bar: comparably expressed sets fall %+.3f below the diagonal at the median (95%% range %+.3f to %+.3f), and the chain sits at percentile %.1f of theirs.",
            ox$c_wt, ox$c_myc, ox$c_int, DIAG, DIAG - abs(ox$c_int),
            ox$null_med_int, ox$null_int_lo, ox$null_int_hi, ox$pct_int),
    sprintf("BEYOND THE DECLARED 0.20, among the mitochondrial arms: %s -- the arms MYC raised most at six weeks (%s). Outside their own bars, below the diagonal: %d of the %d mitochondrial arms (%s).",
            paste(sprintf("%s %+.3f", off$label, off$c_int), collapse = ", "),
            paste(sprintf("%+.2f", off$c_myc_6W), collapse = ", "),
            length(out_bar), nrow(mito), paste(out_bar, collapse = ", ")),
    sprintf("THE TWO ARMS THAT ARE NOT MITOCHONDRIAL (open): proliferation %+.3f / %+.3f, %+.3f from the diagonal (percentile %.1f of its null); the TEB programme %+.3f / %+.3f, %+.3f ABOVE it -- the one arm the Myc+ gland loses less of than the wild-type gland does.",
            pf$c_wt, pf$c_myc, pf$c_int, pf$pct_int, tb$c_wt, tb$c_myc, tb$c_int),
    sprintf("ON THE WILD-TYPE AXIS ALONE, which is what the text's percentiles are: the respiratory chain at percentile %.1f of its matched null and proliferation at percentile %.1f in these draws, against %.1f and %.1f in script 43's -- the same null redrawn, and the text quotes script 43's.",
            ox$pct_wt, pf$pct_wt, pct43("OXPHOS subunits"), pct43("PROLIF_* pooled"))),
  bounds = c(
    "BATCH = TIMEPOINT, on BOTH axes: each coordinate is a temporal contrast and is DESCRIBED, not claimed. What is batch-clean is the vertical distance from the diagonal, because genotype is balanced within each extraction batch. Read the panel down from the line, never along an axis.",
    sprintf("TWO RULES, AND THEY DISAGREE ABOUT THE RESPIRATORY CHAIN. The declared magnitude (%.2f, fixed before the numbers were read and applied as declared) puts it on the diagonal; its matched null puts it below. The magnitude rule has no arm-level p-value to pair with, so the bar is reported beside it rather than used to override it.",
            DIAG),
    sprintf("THIS RULER SITS NEARER THE DECLARED BOUNDARY THAN mitoPPS DOES: the chain is %+.3f from the diagonal here and %+.3f on the mitoPPS ruler (plane_arms_mitopps.R), where every mitochondrial arm lies within %.3f of it. The two rulers answer different questions -- this one what the compartment's transcripts do, mitoPPS what the compartment spends its budget on -- and neither is shown alone.",
            ox$c_int, ox$mitopps_int, max(abs(mito$mitopps_int))),
    "A matched null answers \"more than comparably expressed genes\", never \"more than chance\", and the matched sets are themselves below the diagonal: MYC's effect fades across the window for ordinary genes too.",
    "An unweighted set average counts a rare transcript as much as an abundant one; the expression-weighted reading of the same arm differs (script 45), so a number quoted from this panel must name this ruler."),
  source = c(
    "results/two_timeline_verification.rds (scripts/54_two_timeline_verification.R) -- $arms_content (coordinates), $arm_null (the bars and percentiles), $arm_diagonal (the declared threshold on both rulers)",
    "results/substrate_specificity_tradeoff.rds (scripts/43_substrate_specificity_and_tradeoff.R) -- $wt_null, the wild-type percentiles the text quotes",
    "The grammar: two_timeline_base() in figures/panels/_panel_common.R, with its declared band; the frame: plane_content_lim(), shared with plane_four_genes.R"))

save_panel_p(p, "plane_arms_content", height = 84)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  d[order(d$c_int), c("arm", "n_genes", "c_wt", "c_myc", "c_int", "pct_wt", "pct_myc",
                      "pct_int", "null_int_lo", "null_int_hi", "mitopps_int")] |>
    print(row.names = FALSE, digits = 3)
}
