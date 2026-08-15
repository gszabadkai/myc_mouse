# =============================================================================
# figS2_oxphos_subunit_heatmap.R -- why the two abundance rulers disagree about
# the respiratory chain, drawn gene by gene
# -----------------------------------------------------------------------------
# SLOT: not currently cited.
#
# THE PROBLEM THIS PANEL EXISTS TO SHOW. On the 6W_wt -> 12W_myc diagonal the two
# abundance rulers give different answers about the SAME 87 nuclear OXPHOS
# subunits:
#
#     unweighted per-gene mean log2FC   +0.061   ("returns to baseline")
#     expression-WEIGHTED mean          +0.201
#     log2(sum of normalised counts)    +0.226   (the compartment-level ruler)
#
# Both are correct. They weight differently -- a mean of log ratios against a log
# of the sum -- and the weighted mean landing between the other two is what proves
# the gap IS the weighting rather than a membership or reconciliation difference
# (script 45 PART H asserts that ordering). What the numbers cannot show is WHERE
# the weight sits, and that is what this panel draws.
#
# THE THREE REGIONS, all sharing the row axis:
#   LEFT    the WEIGHT. Each gene's absolute expression in the 6W wild-type gland,
#           log10 normalised counts, as a bar whose length IS its contribution to
#           the summed ruler. Genes run highest-first inside each complex, so the
#           gradient runs top to bottom in every block.
#   MIDDLE  the four states, as asked: the median of the six animals' normalised
#           counts per group, row-centred against that gene's own 24-sample mean.
#   RIGHT   the CHANGE. Each gene's diagonal log2FC on a continuous axis, with the
#           two summary lines drawn on it -- so the reader sees the unweighted mean
#           sitting in the middle of the cloud and the weighted mean pulled right,
#           towards where the long bars are.
#
# THE READING: long bars at the top of each block sit to the RIGHT on the change
# axis; short bars at the bottom sit to the LEFT, several of them below zero. That
# juxtaposition IS the +0.061-versus-+0.226 gap.
#
# AND IT IS NOT AN ARTEFACT OF NOISY LOW-EXPRESSED GENES. That was the obvious
# alternative -- lowly expressed genes carry noisier fold changes, so some
# expression gradient is expected generically. Script 45 PART H tests it against
# 2000 expression-matched random gene sets and the observed gradient sits at the
# 100th percentile on both the wild-type window and the diagonal, against a null
# median of the OPPOSITE sign (-0.17 and -0.22). Matched random sets show high
# expressers moving DOWN relative to low ones; the respiratory chain does the
# reverse. The legend block carries the numbers.
#
# WHERE THE GRADIENT COMES FROM, and a correction to an earlier reading of it: it
# is present in the wild-type window (Spearman +0.33) AND in the Myc effect at six
# weeks (+0.39), weak in the Myc effect at twelve weeks (+0.20) and absent in the
# Myc+ timeline (+0.07). It is NOT a wild-type-only phenomenon -- the oncogene's
# induction is itself expression-graded, preferentially amplifying the subunits
# that are already abundant. Script 45's null covers all four (an earlier version
# covered only two, which would have left the strongest one untested); the legend
# reads the roster off the object rather than naming it, so it stays true either
# way.
#
# Reads (read-only, no re-run):
#   results/state_readings.rds  (script 45 PART H) -- $oxphos_genes (per gene: the
#                                  four group medians, the 6W_wt level, the complex,
#                                  and the log2FC on every contrast), $weighting,
#                                  $gradient, $gradient_quartiles, $gradient_null
# Output: outputs/figures/panels/figS2_oxphos_subunit_heatmap.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

sr_path <- here::here("results", "state_readings.rds")
require_fresher_than(sr_path)
sr <- readRDS(sr_path)

g    <- as.data.frame(sr$oxphos_genes)
wt_  <- as.data.frame(sr$weighting)
grd  <- as.data.frame(sr$gradient)
gq   <- as.data.frame(sr$gradient_quartiles)
gnul <- as.data.frame(sr$gradient_null)

# --- guards ------------------------------------------------------------------
stopifnot(nrow(g) == 87L, !anyNA(g$level_6W_wt), !anyNA(g$lfc_cross))
# THE ORDERING IS THE PANEL'S PREMISE: the gap is the weighting only if the
# weighted mean lands between the unweighted mean and the sum ratio.
UW <- wt_$value[wt_$summary == "unweighted per-gene mean"]
WW <- wt_$value[wt_$summary == "expression-weighted mean"]
SR <- wt_$value[wt_$summary == "log2(sum ratio)"]
stopifnot(length(UW) == 1L, UW < WW, abs(WW - SR) < abs(UW - SR))

# =============================================================================
# rows: six complex blocks, highest expresser first inside each
# =============================================================================
# Complex membership is MitoCarta's own CI-CV subunit sets. The two genes that
# belong to no complex (cytochrome c and its synthase) get a labelled block rather
# than being dropped, so the 87 on this panel are the same 87 every number above
# is computed on.
CX <- c("CI", "CII", "CIII", "CIV", "CV", "cytochrome c / other")
stopifnot(setequal(g$complex, CX))
# The strip column is as wide as its longest label, and "cytochrome c / other"
# alone cost more width than the five complex names put together. Shortened for
# the strip only; the legend block spells out what the block holds.
CX_SHORT <- c(CX[1:5], "cyt c")
g$complex <- factor(g$complex, levels = CX, labels = CX_SHORT)

# Three of 87 have no symbol in the cached ortholog table; the ENSMUSG id is drawn
# rather than a blank row, so the panel still accounts for all 87.
g$lab <- ifelse(is.na(g$symbol) | g$symbol == "", g$gene, g$symbol)
stopifnot(!anyNA(g$lab), !anyDuplicated(g$lab))

# THE COMPETING EXPLANATION, AND IT HAD TO BE TESTED BEFORE THIS PANEL COULD BE
# DRAWN. The lowest-expressed subunits of complex IV are not minor members of the
# mammary respiratory chain -- they are the TISSUE-RESTRICTED PARALOGS, whose
# canonical expression is heart and skeletal muscle (Cox6a2, Cox7a1, Cox8b), lung
# (Cox4i2) and testis (Cox6b2). In an enzymatically dissociated MEC prep their
# signal is plausibly residual non-epithelial tissue, so a fall in them could be a
# drop in CONTAMINATION rather than a withdrawal from respiration -- and they sit
# at expression ranks 83 to 87 of 87, exactly where the gradient is steepest.
#
# The sensitivity is computed below and reported in the legend rather than being
# left as a worry: dropping all five moves the gradient from +0.39 to +0.30 and the
# unweighted mean from +0.061 to +0.105. So they INFLATE the effect and do not
# create it, and the weighted mean does not move at all (they carry no weight,
# which is the point of the panel). They are marked on the page as triangles.
PARALOG <- c("Cox6a2", "Cox7a1", "Cox8b", "Cox4i2", "Cox6b2")
g$paralog <- !is.na(g$symbol) & g$symbol %in% PARALOG
stopifnot(sum(g$paralog) == length(PARALOG))
gx <- g[!g$paralog, ]
wx <- gx$level_6W_wt / sum(gx$level_6W_wt)
SENS <- list(
  n         = nrow(gx),
  rho_all   = stats::cor(log10(g$level_6W_wt),  g$lfc_cross,  method = "spearman"),
  rho_drop  = stats::cor(log10(gx$level_6W_wt), gx$lfc_cross, method = "spearman"),
  uw_drop   = mean(gx$lfc_cross),
  ww_drop   = sum(wx * gx$lfc_cross),
  ranks     = sort(match(PARALOG, g$symbol[order(-g$level_6W_wt)])))
# the sensitivity must point the way the comment says, or the comment is wrong
stopifnot(SENS$rho_drop > 0.2, SENS$rho_drop < SENS$rho_all,
          SENS$uw_drop > UW, abs(SENS$ww_drop - WW) < 0.02)

g <- g[order(g$complex, -g$level_6W_wt), ]
g$row <- stats::ave(seq_len(nrow(g)), g$complex, FUN = function(i) -seq_along(i))
g$lab <- factor(g$lab, levels = g$lab[order(g$complex, g$row)])

blank_y <- ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                          axis.ticks.y = ggplot2::element_blank())
blocks <- ggplot2::facet_grid(rows = ggplot2::vars(complex),
                              scales = "free_y", space = "free_y", switch = "y")

# =============================================================================
# LEFT -- the weight
# =============================================================================
# Bar length is log10 expression, and the fill repeats it on the declared
# sequential ramp: the "encoded twice" idiom of Figs. 1B and 2F, so no colour key
# is needed and the axis is the key. A sequential ramp, NOT ms_diverging -- an
# expression level has a low end and a high end but no meaningful zero, and on a
# page where brown and mint mean down and up a diverging fill would read as sign.
g$lev10 <- log10(g$level_6W_wt)
# The bar runs from the AXIS FLOOR, not from zero: on a log10 axis spanning 1.25
# to 3.82 a bar anchored at zero would be three-quarters constant and would carry
# almost no length information (and geom_col, which does anchor at zero, silently
# drops every bar when the scale excludes it). The floor is padded just below the
# smallest value so the shortest bar is still visible, and the axis is labelled in
# counts so the anchor is inspectable.
LLIM <- c(min(g$lev10) - 0.07 * diff(range(g$lev10)), max(g$lev10))

pW <- ggplot2::ggplot(g, ggplot2::aes(x = lev10, y = lab, fill = lev10)) +
  ggplot2::geom_segment(ggplot2::aes(x = LLIM[1], xend = lev10,
                                     y = lab, yend = lab, colour = lev10),
                        linewidth = 0.9, lineend = "butt") +
  level_fill(LLIM) +
  ggplot2::scale_colour_gradientn(
    colours = unname(ms_sequential[c("low", "mid", "high")]),
    limits = LLIM, oob = scales::squish, guide = "none") +
  ggplot2::guides(fill = "none") +
  blocks +
  ggplot2::scale_x_continuous(
    limits = LLIM, expand = ggplot2::expansion(mult = 0),
    breaks = c(2, 3), labels = c("100", "1k")) +
  ggplot2::labs(x = "level, 6W_wt", y = NULL) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    axis.text.y     = ggplot2::element_text(size = 3.7, margin = ggplot2::margin(r = 0.4)),
    axis.ticks.y    = ggplot2::element_blank(),
    axis.title.x    = ggplot2::element_text(size = 5),
    strip.placement = "outside",
    strip.text.y.left = ggplot2::element_text(angle = 0, size = 5, face = "bold"),
    panel.grid      = ggplot2::element_blank(),
    panel.spacing.y = ggplot2::unit(0.6, "mm"),
    plot.margin     = ggplot2::margin(1.5, 0.6, 1, 1.5, "mm"))

# =============================================================================
# MIDDLE -- the four states, as asked
# =============================================================================
# Per gene, the MEDIAN of its six animals' normalised counts in each group,
# row-centred as log2 against that gene's own 24-sample mean. Row-centring is what
# makes all four columns informative; against the 6W_wt column the first column
# would be zero by construction and would carry nothing.
STATES <- c(med_6W_neg = "6W_wt", med_12W_neg = "12W_wt",
            med_6W_pos = "6W_myc", med_12W_pos = "12W_myc")
bod <- do.call(rbind, lapply(names(STATES), function(k) {
  data.frame(lab = g$lab, complex = g$complex,
             state = factor(unname(STATES[k]), levels = unname(STATES)),
             v = log2(g[[k]] + 1) - log2(g$mean_all + 1))
}))
BLIM <- range(bod$v)
stopifnot(BLIM[1] < 0, BLIM[2] > 0)

pB <- ggplot2::ggplot(bod, ggplot2::aes(x = state, y = lab, fill = v)) +
  ggplot2::geom_tile(colour = NA) +
  heat_fill(BLIM, name = "log2 vs\ngene mean") +
  blocks +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0)) +
  ggplot2::labs(x = NULL, y = NULL) +
  theme_panel(base_size = 6) +
  blank_y +
  ggplot2::theme(
    axis.text.x       = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5,
                                              size = 4.6),
    legend.position   = "bottom",
    legend.key.width  = ggplot2::unit(5, "mm"),
    legend.key.height = ggplot2::unit(1.4, "mm"),
    legend.title      = ggplot2::element_text(size = 4.2),
    legend.text       = ggplot2::element_text(size = 4.2),
    legend.margin     = ggplot2::margin(-1, 0, 0, 0, "mm"),
    strip.text.y      = ggplot2::element_blank(),
    panel.grid        = ggplot2::element_blank(),
    panel.spacing.y   = ggplot2::unit(0.6, "mm"),
    plot.margin       = ggplot2::margin(1.5, 0.6, 1, 0.6, "mm"))

# =============================================================================
# RIGHT -- the change, on a continuous axis, with the two summaries drawn on it
# =============================================================================
# A continuous axis rather than a fourth heat column, because the payoff has to be
# DRAWN: the two summary lines can only be placed against a real scale. The
# unweighted mean sits in the middle of the cloud; the weighted mean sits to its
# right, pulled towards the long bars at the top of each block.
CLIM <- range(g$lfc_cross) + c(-1, 1) * diff(range(g$lfc_cross)) * 0.04

pC <- ggplot2::ggplot(g, ggplot2::aes(x = lfc_cross, y = lab)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
  ggplot2::geom_vline(xintercept = UW, linewidth = 0.35, linetype = "22",
                      colour = "grey30") +
  ggplot2::geom_vline(xintercept = WW, linewidth = 0.35, linetype = "solid",
                      colour = unname(contrast_cols[[contrast_net]])) +
  # triangles are the five tissue-restricted paralogs -- the competing explanation
  # for the gradient, marked on the page rather than left to the legend
  ggplot2::geom_point(ggplot2::aes(fill = lfc_cross,
                                   shape = ifelse(paralog, "paralog", "subunit")),
                      size = 0.9, stroke = 0.15, colour = "grey40") +
  ggplot2::scale_shape_manual(values = c(subunit = 21, paralog = 24),
                              guide = "none") +
  heat_fill(range(g$lfc_cross)) +
  ggplot2::guides(fill = "none") +
  blocks +
  ggplot2::scale_x_continuous(limits = CLIM, labels = lab_signed,
                              breaks = c(-1, -0.5, 0, 0.5),
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::labs(x = "diagonal (log2FC)", y = NULL) +
  theme_panel(base_size = 6) +
  blank_y +
  ggplot2::theme(
    axis.title.x    = ggplot2::element_text(size = 5),
    axis.text.x     = ggplot2::element_text(size = 4.6),
    strip.text.y    = ggplot2::element_blank(),
    panel.grid      = ggplot2::element_blank(),
    panel.spacing.y = ggplot2::unit(0.6, "mm"),
    plot.margin     = ggplot2::margin(1.5, 1.5, 1, 0.6, "mm"))

# The key for the two lines, named on the page because the panel must talk by
# naming its drawn elements. Placed as its own strip under the three regions so it
# cannot collide with a facet.
keydf <- data.frame(
  x = c(0.06, 0.06), y = c(1, 0),
  lab = c(sprintf("unweighted mean of the 87 genes  %s", lab_signed(round(UW, 3))),
          sprintf("expression-weighted mean  %s", lab_signed(round(WW, 3)))),
  col = c("grey30", unname(contrast_cols[[contrast_net]])),
  lty = c("22", "solid"))
pK <- ggplot2::ggplot(keydf) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 0.045, y = y, yend = y,
                                     colour = col, linetype = lty),
                        linewidth = 0.35) +
  ggplot2::geom_text(ggplot2::aes(x = x, y = y, label = lab), hjust = 0,
                     size = 1.6, colour = "grey20") +
  # the third drawn element the reader has to be told the name of
  ggplot2::annotate("point", x = 0.545, y = 0.5, shape = 24, size = 1.0,
                    stroke = 0.15, colour = "grey40", fill = "grey85") +
  ggplot2::annotate("text", x = 0.572, y = 0.5, hjust = 0, size = 1.6,
                    colour = "grey20",
                    label = sprintf("tissue-restricted paralog (%d)", sum(g$paralog))) +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_linetype_identity() +
  ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = 0)) +
  ggplot2::scale_y_continuous(limits = c(-0.6, 1.6), expand = ggplot2::expansion(mult = 0)) +
  ggplot2::labs(x = NULL, y = NULL) +
  theme_panel(base_size = 6) +
  ggplot2::theme(axis.text = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_blank(),
                 axis.line = ggplot2::element_blank(),
                 panel.grid = ggplot2::element_blank(),
                 plot.margin = ggplot2::margin(0, 1.5, 0.5, 1.5, "mm"))

# =============================================================================
# assembly
# =============================================================================
rows <- patchwork::wrap_plots(pW, pB, pC, nrow = 1, widths = c(0.40, 0.26, 0.34))
p    <- patchwork::wrap_plots(rows, pK, ncol = 1, heights = c(24, 1))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
GR <- function(k, col) grd[[col]][grd$contrast == k]
NU <- function(k, st, col) gnul[[col]][gnul$contrast == k & gnul$statistic == st]
# Which contrasts script 45 actually nulled, read from the object rather than
# assumed. The first version of PART H nulled two; it now nulls four, and this
# legend must not claim more or fewer than the object in front of it holds.
NULLED <- unique(gnul$contrast)
stopifnot(length(NULLED) >= 2L, all(c("wt_time", "cross") %in% NULLED),
          all(vapply(NULLED, NU, numeric(1), "spearman rho", "null_median") < 0),
          all(vapply(NULLED, NU, numeric(1), "spearman rho", "percentile") > 90))
q  <- gq[order(gq$median_level), ]
top5 <- utils::head(g[order(-g$level_6W_wt), ], 5)
bot5 <- utils::head(g[order(g$level_6W_wt), ], 5)

LEGEND <- panel_legend(
  slot = "not currently cited",
  what = paste0(
    "The 87 nuclear-encoded OXPHOS subunits, one row each, in six complex blocks, ",
    "ordered by expression within every block. Left, the WEIGHT: each gene's ",
    "absolute level in the six-week wild-type gland, log10 normalised counts, bar ",
    "length and fill carrying the same number. Middle, the four states: the median ",
    "of the six animals' normalised counts per group, row-centred against that ",
    "gene's own 24-sample mean. Right, the CHANGE: that gene's log2 fold change on ",
    "the 6W_wt to 12W_myc diagonal, with the unweighted and expression-weighted ",
    "means of all 87 drawn as vertical lines and the five tissue-restricted ",
    "paralogs drawn as triangles."),
  detail = c(
    "n = 6 animals per group. Medians are drawn; the fitted log2 fold changes on the right are DESeq2 raw (unshrunken) values, and the two agree at Spearman 0.98, so the choice of median over mean changes nothing.",
    sprintf("THE THREE SUMMARIES OF THE SAME 87 GENES: unweighted per-gene mean %+.3f, expression-weighted mean %+.3f, log2 of the summed normalised counts %+.3f. The weighted mean landing BETWEEN the other two is what shows the gap is the WEIGHTING and not a difference of membership or symbol reconciliation -- script 45 asserts that ordering rather than assuming it.",
            UW, WW, SR),
    sprintf("THE GRADIENT, WHICH IS WHAT THE LEFT AND RIGHT REGIONS SHOW TOGETHER: expression correlates with the diagonal fold change at Spearman %+.2f (Pearson %+.2f). By expression quartile the diagonal runs %s -- the lowest-expressed quarter of the respiratory chain FALLS over this window while the highest-expressed quarter rises.",
            GR("cross", "spearman"), GR("cross", "pearson"),
            paste(sprintf("%s %+.3f", q$quartile, q$cross), collapse = ", ")),
    sprintf("AND IT IS NOT THE NOISE OF LOW-EXPRESSED GENES, which is the obvious alternative and the reason this needed a null. Against 2000 expression-matched random gene sets the observed gradient sits at %s -- and in every case the null MEDIAN has the OPPOSITE SIGN. Matched random genes show high expressers moving DOWN relative to low ones; the respiratory chain does the reverse. %s",
            paste(sprintf("the %.0fth percentile on %s (null median %+.2f)",
                          vapply(NULLED, NU, numeric(1), "spearman rho", "percentile"),
                          NULLED,
                          vapply(NULLED, NU, numeric(1), "spearman rho", "null_median")),
                  collapse = ", "),
            sprintf("The quartile GAP gives %s.",
                    paste(sprintf("%s %.0fth", NULLED,
                                  vapply(NULLED, NU, numeric(1), "Q4 - Q1 gap",
                                         "percentile")), collapse = ", "))),
    sprintf("WHERE THE GRADIENT SITS AMONG THE CONTRASTS: wild-type 6 to 12 weeks %+.2f and the Myc effect at six weeks %+.2f both carry it; the Myc effect at twelve weeks is weaker (%+.2f) and the Myc+ timeline is flat (%+.2f). It is therefore NOT a wild-type-only phenomenon: the oncogene's induction is itself expression-graded, preferentially amplifying the subunits that are already abundant. %s",
            GR("wt_time", "spearman"), GR("myc_6W", "spearman"),
            GR("myc_12W", "spearman"), GR("myc_time", "spearman"),
            if (all(c("wt_time", "cross", "myc_6W", "myc_12W") %in% NULLED))
              "All four have been tested against the matched null."
            else sprintf("Only %s have been tested against the matched null.",
                         paste(NULLED, collapse = " and "))),
    sprintf("BY EXPRESSION QUARTILE ON THE WILD-TYPE WINDOW the fall is monotone -- %s -- while the Myc effect at twelve weeks is not (%s). The maturing gland withdraws disproportionately from the minor subunits and holds the high-abundance core. Note that the quartile means for the Myc effect AT TWELVE WEEKS are roughly even; that is specific to twelve weeks and is NOT true of the Myc effect at six, which is graded as steeply as the wild-type window.",
            paste(sprintf("%s %+.3f", q$quartile, q$wt_time), collapse = ", "),
            paste(sprintf("%s %+.3f", q$quartile, q$myc_12W), collapse = ", ")),
    sprintf("THE FIVE HIGHEST-EXPRESSED SUBUNITS, which carry most of the summed ruler: %s. And the five lowest: %s.",
            paste(sprintf("%s (%s, %s counts, %+.2f)", top5$lab, top5$complex,
                          format(round(top5$level_6W_wt), big.mark = ","),
                          top5$lfc_cross), collapse = "; "),
            paste(sprintf("%s (%s, %s counts, %+.2f)", bot5$lab, bot5$complex,
                          format(round(bot5$level_6W_wt), big.mark = ","),
                          bot5$lfc_cross), collapse = "; ")),
    sprintf("COMPLEX MEMBERSHIP is MitoCarta's own CI-CV subunit sets (CI %d, CII %d, CIII %d, CIV %d, CV %d). The %d genes belonging to no complex -- cytochrome c and its synthase -- are drawn in their own block rather than dropped, so the rows on this panel are exactly the 87 every number above is computed on.",
            sum(g$complex == "CI"), sum(g$complex == "CII"), sum(g$complex == "CIII"),
            sum(g$complex == "CIV"), sum(g$complex == "CV"),
            sum(g$complex == "cyt c"))),
  bounds = c(
    "BATCH = TIMEPOINT, AND EVERY CONTRAST ON THIS PANEL SPANS IT. The diagonal and the wild-type window both cross the two extraction batches, so the gradient is DESCRIBED, not claimed. What the matched null rules out is that the gradient is a generic property of the expression range; it cannot rule out that a batch offset is itself expression-dependent.",
    sprintf("THE FIVE LOWEST-EXPRESSED SUBUNITS ARE TISSUE-RESTRICTED PARALOGS, and they are the competing explanation for the gradient. %s are canonically heart and skeletal muscle, %s is lung and %s is testis; in an enzymatically dissociated MEC prep their signal is plausibly residual non-epithelial tissue, so a fall in them could be a drop in CONTAMINATION rather than a withdrawal from respiration. They occupy expression ranks %s of 87 -- exactly where the gradient is steepest -- and they are drawn as TRIANGLES. The sensitivity: dropping all five moves the gradient from %+.2f to %+.2f and the unweighted mean from %+.3f to %+.3f, while the expression-weighted mean does not move (%+.3f, they carry no weight). So they INFLATE the discrepancy and do not create it -- but the headline unweighted value is partly them, and a sentence quoting %+.3f should quote %+.3f beside it.",
            "Cox6a2, Cox7a1 and Cox8b", "Cox4i2", "Cox6b2",
            paste(SENS$ranks, collapse = ", "),
            SENS$rho_all, SENS$rho_drop, UW, SENS$uw_drop, SENS$ww_drop,
            UW, SENS$uw_drop),
    "THE GRADIENT IS A CORRELATION AMONG 87 GENES, not a per-gene test. No individual subunit on this panel is claimed to move, and none is marked as significant.",
    "A SHARED-DENOMINATOR BIAS RUNS AGAINST THE OBSERVED SIGN, which is why the null median is negative. The level is measured in the 6W wild-type group and every drawn contrast has that group in its denominator, so noise pushes level down and fold change up -- inducing a NEGATIVE correlation. The observed positive gradient is therefore conservative, and the matched draws carry the same structure so the comparison is like for like.",
    "ROW-CENTRING IS AGAINST EACH GENE'S OWN 24-SAMPLE MEAN, so the middle region shows each subunit's four states relative to itself and NOT relative to other subunits. Nothing about the relative abundance of two genes can be read from the middle region -- that is the left region's job.",
    "THE TWO SUMMARY LINES ARE MEANS OF THE 87 DRAWN VALUES. The third summary quoted in the text, the log2 of the summed counts, is not a mean of these points and so is not drawn as a line on their axis; it appears in this legend only.",
    "Three of the 87 have no symbol in the cached ortholog table and are labelled with their ENSMUSG identifier. They are included in every computation.",
    "n = 6 per cell. This panel EXPLAINS a discrepancy between two summaries; it is not a test of either."),
  source = c(
    "results/state_readings.rds (scripts/45_state_readings.R PART H) -- $oxphos_genes (per gene: the four group medians, the 6W_wt level, complex membership, and the log2 fold change on all five contrasts), $weighting (the three summaries), $gradient and $gradient_quartiles, $gradient_null (2000 expression-matched draws)",
    "expression levels: DESeq2::sizeFactors(results/dds_group_run.rds) applied to results/count_matrix.rds (unfiltered)",
    "the diagonal contrast: results(dds_group, contrast = c(\"group\", \"12W_pos\", \"6W_neg\"), filterFun = ihw) -- no refit",
    "complex membership: results/mitopps_scores.rds (scripts/08_mitoPPS_analysis.R), the CI-CV subunit sets from Mouse.MitoCarta3.0.xls Sheet 4"))

save_panel_p(p, "figS2_oxphos_subunit_heatmap", height = 185)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the three summaries, and the ordering the panel's premise rests on
  wt_ |> print(row.names = FALSE, digits = 4)

  ## the gradient on every contrast, and the quartile tables behind it
  grd |> print(row.names = FALSE, digits = 3)
  gq  |> print(row.names = FALSE, digits = 3)
  gnul |> print(row.names = FALSE, digits = 3)

  ## the genes, highest expresser first -- the top of every block
  g[order(-g$level_6W_wt),
    c("lab", "complex", "level_6W_wt", "lfc_cross", "lfc_wt_time", "lfc_myc_12W")] |>
    head(15) |> print(row.names = FALSE, digits = 3)

  ## and the bottom, which is where the developmental withdrawal sits
  g[order(g$level_6W_wt),
    c("lab", "complex", "level_6W_wt", "lfc_cross", "lfc_wt_time", "lfc_myc_12W")] |>
    head(15) |> print(row.names = FALSE, digits = 3)

  ## per complex: does the gradient hold inside each one, or is it between them?
  do.call(rbind, lapply(levels(g$complex), function(k) {
    d <- g[g$complex == k, ]
    data.frame(complex = k, n = nrow(d),
               rho = if (nrow(d) >= 5)
                 stats::cor(log10(d$level_6W_wt), d$lfc_cross,
                            method = "spearman") else NA_real_,
               mean_lfc = mean(d$lfc_cross))
  })) |> print(row.names = FALSE, digits = 3)
}
