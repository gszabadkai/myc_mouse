# =============================================================================
# fig2_priming_ratios.R -- nine pro:anti transcript ratios against a line fitted
# inside the set, and PUMA:Bcl-xL read against that line without a verdict
# -----------------------------------------------------------------------------
# SLOT: Fig. 2G. Rebuilt 2026-09-21 on the author's ruling 4, from script 54's
# object and nothing else (docs/2026-09-21_two_timeline_verification.md,
# sections 1.2, 1.3, 1.5 and 3). The slug is unchanged: a filename is a filename
# (the third-round ruling on the N3 residue), and PANELS.md, rebuild_panels.R
# and the analysis record all key on it.
#
# WHAT RULING 4 RETIRED. The 2026-08-05 panel drew seven ratios against a line
# IMPORTED from another estimator -- script 44's gene-level rescaling rate -- and
# read each ratio's retention against it. A ratio's retention is a quotient of
# two per-animal OLS coefficients; the rate is a slope over the DESeq2 fold
# changes of the genes Myc moves. Setting one against the other compares two
# estimators, so "fades at the programme-wide rate" was never a like-for-like
# statement. Ruling 4: the line is fitted INSIDE the ratio set -- through the
# origin, to the eight ratios that are not PUMA:Bcl-xL -- with no imported rate
# and no hard-coded 0.55. Script 54 PART B fits it; this panel draws it and
# asserts that what it draws IS that fit.
#
# NO VERDICT. No criterion for "off the line" was declared before the numbers
# were read, so none is applied here. The legend reports PUMA:Bcl-xL's residual,
# its size in residual SDs, its rank among the nine, and whether it lies inside
# the line's 95% prediction interval. It is the largest departure of the nine,
# and it lies inside the interval. Each of those is asserted, so a re-run of
# script 54 that changes one stops the panel instead of leaving the legend behind.
#
# WHAT IS DRAWN
#   the band    the line's 95% prediction interval, script 54's `$ratio_band`,
#               drawn IN FULL. It is what PUMA:Bcl-xL is read against, and a band
#               cropped to the points would hide how wide it is.
#   the line    the through-origin fit to the other eight. Its slope is printed on
#               the page, as Fig. 1G prints its own.
#   the points  all NINE ratios (ruling 3). Shape is the DENOMINATOR -- circles
#               Bcl-xL, triangles Mcl-1 -- the second encoding the 2026-08-05
#               version said a second denominator needs, and its reason for
#               leaving the two Mcl-1 ratios out. PUMA:Bcl-xL is the one OPEN
#               point: it was not used to fit the line, so it is read against a
#               line it did not help to draw.
#
# THE KEY SAYS "not used to fit the line", NOT "not in the fit" (author's
# question, 2026-09-21). The shorter form reads as a verdict -- "does not fit
# the line" -- which is the opposite of what the band shows: PUMA:Bcl-xL lies
# INSIDE the 95% prediction interval. The key names the method, not a result.
#
# THE MEMBERS ARE NOT HERE. The 2026-08-05 legend used the members' own
# retentions to explain why the ratios held; ruling 4 retired that comparison,
# and the manuscript's clause about the pro- and anti-apoptotic members now cites
# the four-gene plane (plane_four_genes.R). The legend says so in its first
# sentence and in a bound, so it cannot be read the other way.
#
# TWO ELEMENTS OF THE 2026-08-05 PANEL ARE GONE, AND WHY
#   * The dashed residual drops. PUMA:Bcl-xL and Bak1:Bcl-xL sit 0.014 apart on
#     the six-week axis, so a drop from PUMA:Bcl-xL up to the line would pass
#     about half a millimetre from Bak1's point, which sits halfway up it, and
#     read as Bak1's. The band carries the reading, and the legend carries the
#     residual.
#   * The filled/open "significantly induced at 6W" shape. It scoped a RETENTION,
#     which is a quotient and meaningless where the six-week effect is near zero.
#     A residual from a fitted line is defined for every ratio, and the fit takes
#     all eight comparators whatever their six-week p. Marking significance on a
#     contrast no claim rests on is what the third round removed from 2G (alt).
#     The nine six-week p-values are printed in the legend rather than dropped.
#     Whether the encoding comes back is the author's call.
#
# AS OF SCRIPT 54's RUN (2026-09-21, 11:53 UTC). The legend formats every number
# from the object and types none:
#
#   line over the other eight   slope 0.429 (SE 0.137), residual SD 0.221, 7 df
#   PUMA:Bcl-xL                 +0.674 -> -0.061, predicted +0.289
#                               residual -0.350 = -1.58 SD, largest of the nine
#                               95% prediction interval -0.276 to +0.855: inside
#   check 3                     IQR of the six-week effect 0.444, floor 0.30
#
# R2 IS THE SQUARED PEARSON CORRELATION (0.43), Fig. 1G's definition. A
# through-origin lm reports an uncentred R2 (0.58 here), which is not comparable
# with anything else in the figure set and is never the number quoted. It stays
# off the page: slope and R2 both round to 0.43, and printed together they read
# as a typing error.
#
# THE WORD (N3). The legend block never gives a transcript ratio the name of a
# cellular state: how close a mitochondrion sits to the apoptotic threshold is a
# property of its proteins. Script 42's object is cited as "Script 42's ..."
# rather than by its filename, and a check before the save asserts the block is
# clear.
#
# Reads (read-only, no re-run):
#   results/two_timeline_verification.rds (script 54) -- $ratios, $ratio_line,
#       $ratio_band, $ratio_resid, $ratio_range, $retention_provenance,
#       $four_genes, $licence, $params
# Output: outputs/figures/panels/fig2_priming_ratios.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

tv_path <- here::here("results", "two_timeline_verification.rds")
require_fresher_than(tv_path)
tv <- readRDS(tv_path)

pr <- as.data.frame(tv$ratios)
ln <- as.data.frame(tv$ratio_line)
rs <- as.data.frame(tv$ratio_resid)
bd <- as.data.frame(tv$ratio_band)
rg <- as.data.frame(tv$ratio_range)
rp <- as.data.frame(tv$retention_provenance)
fg <- as.data.frame(tv$four_genes)
TARGET <- tv$params$target_ratio

stopifnot(nrow(pr) == 9L, !anyDuplicated(pr$pair), identical(TARGET, "Bbc3:Bcl2l1"),
          TARGET %in% pr$pair, setequal(pr$anti, c("Bcl2l1", "Mcl1")),
          identical(rs$pair, pr$pair), identical(rp$pair, pr$pair),
          nrow(ln) == 2L, nrow(bd) >= 50L)

L8 <- ln[ln$comparator == "the other eight ratios", ]
L7 <- ln[ln$comparator == "the seven that are not a PUMA ratio", ]
stopifnot(nrow(L8) == 1L, nrow(L7) == 1L, L8$n_ratios == 8L, L7$n_ratios == 7L)

# Ratio names for the legend: the pro member's gene symbol, except that Bbc3's
# ratios take PUMA, the name the text uses for PUMA:Bcl-xL; the denominator by
# its protein name. Set here, before anything is ranked, so every table below
# carries them. On the page each point carries its pro gene's symbol only,
# because the shape already says which denominator it has.
ANTI_PROT <- c(Bcl2l1 = "Bcl-xL", Mcl1 = "Mcl-1")
pr$name   <- paste0(ifelse(pr$pro == "Bbc3", "PUMA", pr$pro), ":",
                    unname(ANTI_PROT[pr$anti]))
rs$name   <- pr$name[match(rs$pair, pr$pair)]
stopifnot(!anyNA(pr$name), !anyNA(rs$name), !anyDuplicated(pr$name))

# =============================================================================
# the checks: what is drawn IS script 54's fit, and what the legend says holds
# =============================================================================

# --- CHECK 3 (ruling 3): why this is still a scatter --------------------------
# The panel stays a scatter only because the six-week effect has the spread the
# rule asked for. Recomputed over all nine, so the verdict is read from the
# ratios themselves and not only from the saved row.
stopifnot(rg$n_ratios == nrow(pr),
          abs(stats::median(pr$d6) - rg$median_d6) < 1e-12,
          abs(stats::IQR(pr$d6) - rg$iqr_d6) < 1e-12,
          identical(rg$threshold, tv$params$iqr_min),
          rg$iqr_d6 >= rg$threshold, identical(rg$verdict, "keep the scatter"))

# --- the line: refitted from the saved ratios and matched ---------------------
others <- pr[pr$pair != TARGET, ]
tgt    <- pr[pr$pair == TARGET, ]
m8     <- stats::lm(d12 ~ 0 + d6, data = others)
s8     <- summary(m8)
SLOPE  <- unname(stats::coef(m8)[1])
stopifnot(abs(SLOPE - L8$slope) < 1e-12,
          abs(s8$coefficients[1, 2] - L8$slope_se) < 1e-12,
          abs(s8$sigma - L8$sigma) < 1e-12,
          m8$df.residual == L8$df_resid,
          abs(stats::cor(others$d6, others$d12)^2 - L8$r2_pearson) < 1e-12)

# --- the band IS the drawn line's 95% prediction interval ---------------------
# Its centre is the line, and its edges are the closed form for a through-origin
# fit: t on the residual df, times sqrt(sigma^2 + x^2 SE^2). Script 54 computed
# it on a grid, next to the fit, so the panel takes it as saved and checks it.
TQ <- stats::qt(0.975, L8$df_resid)
half_width <- function(x) TQ * sqrt(L8$sigma^2 + x^2 * L8$slope_se^2)
stopifnot(max(abs(bd$fit - SLOPE * bd$d6)) < 1e-12,
          max(abs(bd$lwr - (bd$fit - half_width(bd$d6)))) < 1e-10,
          max(abs(bd$upr - (bd$fit + half_width(bd$d6)))) < 1e-10,
          min(bd$d6) == 0, max(bd$d6) > max(pr$d6))

# --- PUMA:Bcl-xL against the line ---------------------------------------------
pred <- SLOPE * tgt$d6
stopifnot(abs(pred - L8$target_pred) < 1e-12,
          abs((tgt$d12 - pred) - L8$target_resid) < 1e-12,
          abs(L8$target_resid / L8$sigma - L8$target_resid_in_sigma) < 1e-12,
          abs((pred - half_width(tgt$d6)) - L8$pi_lwr) < 1e-10,
          abs((pred + half_width(tgt$d6)) - L8$pi_upr) < 1e-10)

# --- every ratio's residual, and the three statements the legend makes --------
stopifnot(max(abs(rs$residual - (rs$d12 - SLOPE * rs$d6))) < 1e-12,
          identical(as.logical(rs$in_fit), rs$pair != TARGET))
rk <- rs[order(-abs(rs$residual)), ]
# ASSERTED so a re-run of script 54 cannot leave the legend behind: PUMA:Bcl-xL
# is the largest departure of the nine, it lies BELOW the line, and it lies
# INSIDE the line's 95% prediction interval -- on the drawn fit and on the
# sensitivity fit the legend also reports.
stopifnot(rk$pair[1] == TARGET, rs$abs_rank[rs$pair == TARGET] == 1,
          L8$target_resid < 0, isTRUE(L8$target_inside_pi),
          tgt$d12 > L8$pi_lwr, tgt$d12 < L8$pi_upr,
          isTRUE(L7$target_inside_pi))

# --- a retention is computed, never compared into place (the note, 1.3) --------
# Script 54 checked script 42's saved d12/d6 against the recomputed quotient for
# all nine; this re-reads that check. The panel draws no retention and no rate;
# the one retention the legend prints (Bax:Bcl-xL) is formatted from here.
stopifnot(max(abs(rp$retention_saved - rp$d12 / rp$d6)) < 1e-12,
          max(abs(rp$retention_recomputed - rp$d12 / rp$d6)) < 1e-12,
          all(rp$computed_not_compared),
          identical(rp$d6, pr$d6), identical(rp$d12, pr$d12))

# =============================================================================
# the panel
# =============================================================================
KEY <- c(bcl = "ratio to Bcl-xL", mcl = "ratio to Mcl-1", out = "not used to fit the line")
pr$cls <- factor(ifelse(pr$pair == TARGET, KEY[["out"]],
                        ifelse(pr$anti == "Bcl2l1", KEY[["bcl"]], KEY[["mcl"]])),
                 levels = unname(KEY))
stopifnot(sum(pr$cls == KEY[["out"]]) == 1L, sum(pr$cls == KEY[["mcl"]]) == 2L)

# THE FRAME IS THE BAND's. x runs from the origin, where the through-origin line
# starts, to the end of script 54's grid; y takes the whole band.
XR <- range(bd$d6)
YR <- range(c(bd$lwr, bd$upr, pr$d12))
YR <- YR + c(-1, 1) * diff(YR) * 0.03

# Labels are placed, not repelled, as in Figs. 2F and the 2026-08-05 version:
# nine points on a mostly empty plane, each label beside its own point. Offsets
# are in mm, converted with the plotting area's size inside the axes, which is
# MEASURED for this panel at 89 x 64 mm (the sandbox re-measures it) and must be
# re-checked if the size changes. The check below is that no label box contains
# another point and every label stays inside the frame.
PLOT_W <- 78.1; PLOT_H <- 56.4                      # mm, inside the axes
MMX <- diff(XR) / PLOT_W; MMY <- diff(YR) / PLOT_H  # data units per mm
CHW <- 0.95; LNH <- 1.9                             # a character and a line, mm, at size 1.7
GAP <- 1.4                                          # point centre to text, mm

LAB <- data.frame(
  pair = c("Bax:Bcl2l1", "Bid:Bcl2l1", "Bak1:Bcl2l1", "Bbc3:Bcl2l1", "Pmaip1:Bcl2l1",
           "Bcl2l11:Bcl2l1", "Bmf:Bcl2l1", "Bax:Mcl1", "Bbc3:Mcl1"),
  side = c("left", "below", "right", "left", "right",
           "left", "left", "right", "right"),
  stringsAsFactors = FALSE)
stopifnot(setequal(LAB$pair, pr$pair))
LAB <- merge(LAB, pr[, c("pair", "pro", "d6", "d12")], by = "pair")
LAB$hj <- c(left = 1, right = 0, below = 0.5)[LAB$side]
LAB$vj <- c(left = 0.5, right = 0.5, below = 1)[LAB$side]
LAB$x  <- LAB$d6 + c(left = -GAP, right = GAP, below = 0)[LAB$side] * MMX
LAB$y  <- LAB$d12 - c(left = 0, right = 0, below = GAP)[LAB$side] * MMY

# each label's box in data units, and the two checks on it
w <- nchar(LAB$pro) * CHW * MMX
LAB$x1 <- LAB$x - LAB$hj * w;           LAB$x2 <- LAB$x1 + w
LAB$y1 <- LAB$y - LAB$vj * LNH * MMY;   LAB$y2 <- LAB$y1 + LNH * MMY
clash <- vapply(seq_len(nrow(LAB)), function(i)
  any(pr$pair != LAB$pair[i] & pr$d6 > LAB$x1[i] & pr$d6 < LAB$x2[i] &
      pr$d12 > LAB$y1[i] & pr$d12 < LAB$y2[i]), logical(1))
stopifnot(!any(clash), min(LAB$x1) > XR[1], max(LAB$x2) < XR[2],
          min(LAB$y1) > YR[1], max(LAB$y2) < YR[2])

# The slope, top left, in Fig. 1G's form: plotmath, the number QUOTED so the
# trailing digit survives. The corner above the band is empty by construction --
# the band's upper edge starts at +sigma * t at the origin and rises from there.
ANN <- data.frame(x = XR[1] + diff(XR) * 0.025, y = YR[2] - diff(YR) * 0.02,
                  lab = sprintf('slope~"%.2f"', SLOPE))
stopifnot(ANN$y - 2 * LNH * MMY > max(bd$upr[bd$d6 < XR[1] + diff(XR) * 0.2]))

p <- ggplot2::ggplot(pr, ggplot2::aes(d6, d12)) +
  ggplot2::geom_ribbon(data = bd, inherit.aes = FALSE,
                       ggplot2::aes(x = d6, ymin = lwr, ymax = upr),
                       fill = "grey91", colour = NA) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75") +
  # the line over the band's own grid, so the two are drawn from one object
  ggplot2::geom_line(data = bd, inherit.aes = FALSE, ggplot2::aes(x = d6, y = fit),
                     linewidth = 0.45, colour = "grey20") +
  ggplot2::geom_point(ggplot2::aes(shape = cls, fill = cls), size = 1.6,
                      stroke = 0.4, colour = "grey15") +
  ggplot2::geom_text(data = LAB, inherit.aes = FALSE,
                     ggplot2::aes(x = x, y = y, label = pro, hjust = hj, vjust = vj),
                     size = 1.7, colour = "grey15", fontface = "italic") +
  ggplot2::geom_text(data = ANN, inherit.aes = FALSE,
                     ggplot2::aes(x = x, y = y, label = lab), parse = TRUE,
                     hjust = 0, vjust = 1, size = 1.8, colour = "grey25") +
  # one key for two aesthetics: shape is the denominator, the open fill is the
  # ratio the line was not fitted to
  ggplot2::scale_shape_manual(values = stats::setNames(c(21, 24, 21), KEY),
                              breaks = unname(KEY), name = NULL) +
  ggplot2::scale_fill_manual(values = stats::setNames(c("grey15", "grey15", "white"), KEY),
                             breaks = unname(KEY), name = NULL) +
  ggplot2::scale_x_continuous(limits = XR, labels = lab_signed,
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::scale_y_continuous(limits = YR, labels = lab_signed,
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::labs(x = "Myc effect at 6W  (log2 pro:anti ratio)",
                y = "Myc effect at 12W") +
  ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 1.5))) +
  theme_panel(base_size = 6) +
  # Key inside, BOTTOM RIGHT: the band's lower edge rises to the right, and the
  # wedge beneath it in that corner is the one region no point, no band and no
  # line reaches.
  ggplot2::theme(
    legend.position        = "inside",
    legend.position.inside = c(0.995, 0.01),
    legend.justification   = c(1, 0),
    legend.background      = ggplot2::element_blank(),
    legend.margin          = ggplot2::margin(0, 0, 0, 0),
    legend.key.size        = ggplot2::unit(2.4, "mm"),
    legend.spacing.y       = ggplot2::unit(0.3, "mm"),
    plot.margin            = ggplot2::margin(1.5, 2.5, 1, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
row_of <- function(pair) pr[pr$pair == pair, ]
bax    <- row_of("Bax:Bcl2l1")

# THE MEMBERS ARE THE PLANE's, NOT THIS PANEL's (author, 2026-09-21): the
# manuscript's clause about the pro- and anti-apoptotic members cites the
# four-gene plane. That panel has no slot yet, so the legend names its script,
# and the name is asserted: the rename that comes with a slot stops this panel
# rather than leaving its legend pointing at a file that is gone.
PLANE_MEMBERS <- "plane_four_genes.R"
stopifnot(file.exists(here::here("figures", "panels", PLANE_MEMBERS)))
bbc3_g <- fg[fg$gene == "Bbc3", ]
stopifnot(nrow(bax) == 1L, nrow(bbc3_g) == 1L,
          # the licence the legend states, read from the object rather than typed
          grepl("^pre-specified", tv$licence[["Bbc3"]]),
          grepl("denominator", tv$licence[["Bcl2l1"]]))

# p to two significant figures, never in scientific notation; and counts as words,
# so a computed count reads like the typed "nine" beside it
fmt_p <- function(p) ifelse(p < 0.001, sprintf("%.5f", p),
                            ifelse(p < 0.01, sprintf("%.4f", p), sprintf("%.2f", p)))
as_word <- function(n) {
  stopifnot(n == round(n), n >= 1, n <= 9)
  c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine")[n]
}
by_d6 <- pr[order(-pr$d6), ]
n_smaller <- sum(pr$d12 < pr$d6)                    # a smaller Myc effect at 12W
stopifnot(n_smaller == sum(pr$int < 0))             # the same thing, as the interaction

LEGEND <- panel_legend(
  slot = "Fig. 2G",
  what = paste0(
    "Nine pro- to anti-apoptotic transcript ratios, each drawn as the Myc effect on ",
    "its log2 ratio at twelve weeks against the same effect at six. The line is ",
    "fitted through the origin to the eight ratios other than PUMA:Bcl-xL. ",
    "PUMA:Bcl-xL (open circle) was not used to fit it, so it is read against a line ",
    "it did not help to draw. The shaded band is the line's 95% prediction interval. ",
    "Circles are ratios to Bcl-xL and triangles ratios to Mcl-1. Each point is a ",
    "ratio, labelled with its pro-apoptotic member (Bbc3 encodes PUMA); the members ",
    "themselves are not drawn here."),
  detail = c(
    "n = 6 animals per group. Each ratio is log2(pro) - log2(anti) per animal, on log2(normalised count + 1), and its Myc effect at each age is the genotype coefficient at that age in lm(ratio ~ timepoint x genotype) over the 24 animals: Script 42's fits, carried unchanged in script 54's object.",
    sprintf("THE LINE IS FITTED INSIDE THE RATIO SET: through the origin, over the %s ratios other than PUMA:Bcl-xL. Slope %.3f (SE %.3f); residual SD %.3f on %d df; R2 %.2f as the squared Pearson correlation, the definition Fig. 1G uses. No rate from outside the ratio set is drawn or compared.",
            as_word(L8$n_ratios), L8$slope, L8$slope_se, L8$sigma, L8$df_resid,
            L8$r2_pearson),
    sprintf("PUMA:Bcl-xL, NOT USED TO FIT THE LINE: Myc effect %+.3f at six weeks and %+.3f at twelve, against %+.3f predicted by the line. Its residual is %+.3f, or %.2f residual SDs, the largest of the nine in absolute size; the next three are %s %+.3f, %s %+.3f and %s %+.3f. The line's 95%% prediction interval at its six-week effect runs from %+.3f to %+.3f, and PUMA:Bcl-xL lies inside it.",
            tgt$d6, tgt$d12, L8$target_pred, L8$target_resid, L8$target_resid_in_sigma,
            rk$name[2], rk$residual[2], rk$name[3], rk$residual[3],
            rk$name[4], rk$residual[4], L8$pi_lwr, L8$pi_upr),
    sprintf("SENSITIVITY, REPORTED AND NOT DRAWN. The comparison set holds one other PUMA ratio, PUMA:Mcl-1. Fitted to the %s ratios that are not a PUMA ratio, the slope is %.3f and PUMA:Bcl-xL's residual %+.3f (%.2f residual SDs), inside an interval of %+.3f to %+.3f. The drawn line keeps the comparison set as specified: all eight ratios other than PUMA:Bcl-xL.",
            as_word(L7$n_ratios), L7$slope, L7$target_resid, L7$target_resid_in_sigma,
            L7$pi_lwr, L7$pi_upr),
    sprintf("WHY A SCATTER: across all nine ratios the six-week Myc effect has median %+.3f and quartiles %+.3f and %+.3f, an interquartile range of %.3f against the %.2f fixed in advance as the floor for a scatter (below it the panel would have become a retention plot). All nine are drawn, including the two ratios to Mcl-1.",
            rg$median_d6, rg$q1_d6, rg$q3_d6, rg$iqr_d6, rg$threshold),
    sprintf("THE SIX-WEEK EFFECTS, with the raw p of each genotype coefficient: %s. The line is fitted to all eight comparators, whatever their six-week p.",
            paste(sprintf("%s %+.3f (p %s)", by_d6$name, by_d6$d6, fmt_p(by_d6$p6)),
                  collapse = "; "))),
  bounds = c(
    "NO CRITERION FOR \"OFF THE LINE\" WAS DECLARED before the numbers were read, so none is applied, and PUMA:Bcl-xL's residual is reported without a verdict: it is the largest departure of the nine, and it lies inside the line's 95% prediction interval. The panel does not support \"significantly\" for that departure.",
    sprintf("THE RATIO'S OWN INTERACTION IS A DIFFERENT COMPARISON, and not this panel's. PUMA:Bcl-xL's twelve-week Myc effect minus its six-week one is %+.3f (raw p %.3f; Benjamini-Hochberg across the nine ratios %.2f; %+.3f, raw p %.3f, with the epithelial and immune composites as covariates). Its null is no attenuation at all, and %s of the nine ratios have a smaller Myc effect at twelve weeks than at six, so that p asks whether PUMA:Bcl-xL attenuates, not whether it attenuates more than the set does; the line fitted inside the set is the comparison this panel draws. The gene-level Bbc3 interaction, %+.3f (raw p %.4f), is a third statistic and belongs to the MEMBER, which the four-gene plane draws; none of the three stands in for another.",
            tgt$int, tgt$int_p, tgt$int_p_bh, tgt$int_adj, tgt$int_p_adj,
            as_word(n_smaller), bbc3_g$int_lfc, bbc3_g$int_p),
    "PUMA:Bcl-xL IS THE PRE-SPECIFIED PAIR: Bbc3 was named in advance from the PGC1a westerns, and Bcl-xL is its fixed denominator. The other eight ratios are its comparison set, and nothing is claimed for any of them.",
    sprintf("A RETENTION IS COMPUTED, NOT COMPARED. Bax:Bcl-xL's twelve-week effect over its six-week one is %.6f / %.6f = %.4f, taken from its own two coordinates. It is not evidence that the ratio fades at a programme-wide rate: that would set a ratio-level estimate against a gene- or pathway-level slope, two different estimators, and this panel's only comparator is the line fitted inside the ratio set.",
            bax$d12, bax$d6, bax$d12 / bax$d6),
    sprintf("THIS PANEL DRAWS RATIOS, NOT THEIR MEMBERS. A ratio's place against the line says how the balance between its two members moved; it does not say whether either member fell, or whether the two fell together. How the members move is drawn on the four-gene plane (%s), and a sentence about the members cites that panel, not this one.",
            PLANE_MEMBERS),
    "The nine ratios are not independent: seven share the Bcl-xL denominator, two share Mcl-1, and Bax and Bbc3 each enter two ratios. The prediction interval treats the eight comparators' scatter about the line as independent, an assumption the shared members break.",
    "Both axes are the Myc genotype contrast at one age, which is clean: genotype is balanced within each extraction batch. The panel therefore does not carry the batch = timepoint caveat of the temporal panels. It does carry n = 6 per cell.",
    "These are TRANSCRIPT ratios. How close a mitochondrion sits to the apoptotic threshold is a property of its protein complement; BH3 profiling is the measurement, and this panel is a reason to do it rather than a substitute."),
  source = c(
    "results/two_timeline_verification.rds (scripts/54_two_timeline_verification.R) -- $ratios (Script 42's nine ratio fits, carried unchanged; drawn), $ratio_line (the drawn line and the sensitivity fit), $ratio_band (the drawn 95% prediction interval), $ratio_resid (every ratio's residual and its rank), $ratio_range (the scatter check), $retention_provenance, $four_genes (the gene-level Bbc3 interaction), $licence"))

# (That every item arrived is checked by panel_legend() itself, for every panel.
# This panel's first draft is why: its PUMA:Bcl-xL item evaluated to character(0)
# and dropped out of the block without an error.)

# THE N3 CLEARANCE, asserted rather than proof-read: no field of the block may
# carry the word, so a later edit that brings it back stops the panel.
stopifnot(!any(grepl("priming", unlist(LEGEND[c("what", "detail", "bounds", "source")]),
                     ignore.case = TRUE)))

save_panel_p(p, "fig2_priming_ratios", height = 64)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## all nine ratios against the line, largest departure first
  rs[order(-abs(rs$residual)), c("name", "d6", "p6", "d12", "fitted", "residual",
                                 "in_fit", "abs_rank")] |>
    print(row.names = FALSE, digits = 3)

  ## the drawn fit and the sensitivity fit, side by side
  ln[, c("comparator", "n_ratios", "slope", "slope_se", "sigma", "df_resid",
         "r2_pearson", "target_pred", "target_resid", "target_resid_in_sigma",
         "pi_lwr", "pi_upr", "target_inside_pi")] |>
    print(row.names = FALSE, digits = 4)

  ## check 3
  rg |> print(row.names = FALSE, digits = 4)

  ## re-measure the plotting area if the panel size changes: PLOT_W and PLOT_H
  ## above must match these (mm, inside the axes, at the saved size)
  {
    grDevices::pdf(NULL, width = fig_w[["single"]] / 25.4, height = 64 / 25.4)
    g <- ggplot2::ggplotGrob(p)
    fixed <- function(u) grid::unitType(u) != "null"
    c(width  = fig_w[["single"]] -
        sum(grid::convertWidth(g$widths[fixed(g$widths)], "mm", valueOnly = TRUE)),
      height = 64 -
        sum(grid::convertHeight(g$heights[fixed(g$heights)], "mm", valueOnly = TRUE))) |>
      print()
    grDevices::dev.off()
  }
}
