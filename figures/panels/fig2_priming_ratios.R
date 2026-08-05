# =============================================================================
# fig2_priming_ratios.R -- the priming ratios follow the global rescaling, and
# one of them does not
# -----------------------------------------------------------------------------
# SLOT: Fig. 2G.
#
#   "While most apoptotic priming ratios established by MYC at 6W remained
#    stable, since both pro- and anti-apoptotic proteins diminished in accordance
#    with the global rescaling of the MYC effect (Fig. 2G), the PUMA/BCL-XL ratio
#    showed a striking reversal, deviating significantly from the expected
#    pattern."
#
# (Author's correction, 2026-08-05, and it is the right one: what declines
# together is the MEMBERS of each ratio, which is why the ratios themselves are
# left where Myc put them.)
#
# THE PANEL IS FIG. 1G's PLOT, ONE LEVEL UP. There the whole transcriptome sat on
# a line of slope ~0.49 -- the Myc effect at twelve weeks against the Myc effect
# at six. Here the same axes carry seven pro-apoptotic:BCL-XL log ratios, and the
# same line is drawn as the EXPECTED pattern. A ratio that merely follows the
# global rescaling lands on it. A ratio that is selectively dismantled does not.
#
#   Bax    0.901 -> 0.495   retention 0.55   on the line
#   Bid    0.946 -> 0.396              0.42   on the line
#   Bak1   0.688 -> 0.185              0.27   a little under it
#   Bbc3   0.674 -> -0.061            -0.09   CROSSES ZERO
#
# ONE SHARED DENOMINATOR, which is what makes the comparison internal (script 42's
# own design note): every ratio on this panel is against BCL-XL, so "it is just
# the global attenuation" is answered from inside the panel rather than against an
# outside null. If BAX priming retains the global rate while PUMA priming reverses
# against the SAME denominator, the attenuation cannot be the explanation.
#
# WHY THE MEMBERS ARE IN THE LEGEND AND NOT ON THE PANEL. The sentence's clause
# about pro- and anti-apoptotic transcripts declining together is the REASON the
# ratios sit on the line, and it is exact for the two large members: Bax retains
# 0.52 and Bcl-xL 0.57, both at the global 0.49, so their ratio retains 0.55.
# Bbc3 retains -1.10 against that same denominator. Drawing the members as well
# would double the panel to make a point the ratios already carry.
#
# SHAPE IS "ESTABLISHED AT SIX WEEKS", and it is load-bearing rather than
# decorative: retention is a RATIO, so it is meaningless where the six-week effect
# is not distinguishable from zero. Three of the seven pairs are in that state
# (p6 = 0.40 to 0.71) and their positions must not be read as retention.
#
# Reads (read-only, no re-run):
#   results/priming_arm_teb.rds           (script 42) -- $priming, $pair_null
#   results/collapse_module_ownership.rds (script 44) -- $defs$global_rate_fitted,
#                                            $wt_genes for the member retentions
# Output: outputs/figures/panels/fig2_priming_ratios.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

pa_path <- here::here("results", "priming_arm_teb.rds")
cm_path <- here::here("results", "collapse_module_ownership.rds")
require_fresher_than(pa_path)
require_fresher_than(cm_path)
pa  <- readRDS(pa_path)
cmo <- readRDS(cm_path)

pr  <- as.data.frame(pa$priming)
pn  <- as.data.frame(pa$pair_null)
wg  <- as.data.frame(cmo$wt_genes)
RATE <- as.numeric(cmo$defs$global_rate_fitted)     # 0.487, script 44 / Fig. 1G

stopifnot(nrow(pr) == 9L,
          all(c("pair", "pro", "anti", "d6", "p6", "d12", "retention",
                "int", "int_p", "int_p_adj") %in% names(pr)),
          abs(RATE - 0.4872) < 1e-3)

# =============================================================================
# the seven ratios that share BCL-XL
# =============================================================================
# The two Mcl-1 pairs are quoted in the legend, not drawn: a panel with two
# denominators needs a second encoding to say which is which, and the comparison
# the sentence makes only works when the denominator is held fixed.
d <- pr[pr$anti == "Bcl2l1", ]
d$established <- d$p6 < 0.05
stopifnot(nrow(d) == 7L, sum(d$established) == 4L,
          identical(sort(d$pro[d$established]), sort(c("Bax", "Bak1", "Bid", "Bbc3"))),
          # the reversal, asserted so a re-run cannot flip it silently
          d$d12[d$pro == "Bbc3"] < 0, all(d$d12[d$established & d$pro != "Bbc3"] > 0))

# Cross-instrument check: script 42 fits the log ratio per sample, script 44
# reports each member's DESeq2 fold change. The ratio must be the difference of
# its members. They are different fits, so this is a sanity bound, not identity.
# `which()`, not a bare logical: two rows of $wt_genes have an NA symbol, and a
# logical index containing NA returns extra NA elements rather than dropping them.
mem <- function(g, col) wg[[col]][which(wg$gene == g)]
stopifnot(abs((mem("Bax", "myc_6W") - mem("Bcl2l1", "myc_6W")) -
              d$d6[d$pro == "Bax"]) < 0.1)

# =============================================================================
# the panel
# =============================================================================
XR <- range(c(d$d6, 0)) + c(-1, 1) * diff(range(c(d$d6, 0))) * 0.06
YR <- range(c(d$d12, RATE * d$d6)) +
      c(-1, 1) * diff(range(c(d$d12, RATE * d$d6))) * 0.10

# The residual from the expected line, drawn ONLY where retention is interpretable
# -- i.e. only for the four ratios Myc established at six weeks. For the other
# three the six-week effect is not distinguishable from zero, so the distance from
# the line is not a deviation from anything.
res <- d[d$established, ]
res$y_exp <- RATE * res$d6

# Labels are placed, not repelled, with the leader leaving the text on its own
# line where one is needed (Fig. 2F's convention). Seven points on an otherwise
# empty plane, so each label sits beside its own point and the check below is that
# no label box contains another point.
LAB <- data.frame(
  pro   = c("Bax", "Bid", "Bak1",  "Bbc3", "Pmaip1", "Bcl2l11", "Bmf"),
  dx    = c(-0.022, 0,     0.022,   0,      0.022,    0.022,     0.022),
  dy    = c(0,     -0.055, 0,      -0.055,  0,        0,         0),
  hj    = c(1,      0.5,   0,       0.5,    0,        0,         0),
  vj    = c(0.5,    1,     0.5,     1,      0.5,      0.5,       0.5),
  stringsAsFactors = FALSE)
LAB <- merge(LAB, d[, c("pro", "d6", "d12")], by = "pro")
LAB$x <- LAB$d6 + LAB$dx
LAB$y <- LAB$d12 + LAB$dy

# a character's width and a line's height, in data units, from the drawn text size
CH <- diff(XR) / 72
LH <- diff(YR) / 26
LAB$x1 <- LAB$x - ifelse(LAB$hj == 0, 0, ifelse(LAB$hj == 1, nchar(LAB$pro) * CH,
                                                nchar(LAB$pro) * CH / 2))
LAB$x2 <- LAB$x1 + nchar(LAB$pro) * CH
LAB$y1 <- LAB$y - ifelse(LAB$vj == 1, LH, LH / 2)
LAB$y2 <- LAB$y1 + LH
clash <- vapply(seq_len(nrow(LAB)), function(i)
  any(d$pro != LAB$pro[i] & d$d6 > LAB$x1[i] & d$d6 < LAB$x2[i] &
      d$d12 > LAB$y1[i] & d$d12 < LAB$y2[i]), logical(1))
stopifnot(!any(clash), min(LAB$x1) > XR[1], max(LAB$x2) < XR[2])

p <- ggplot2::ggplot(d, ggplot2::aes(d6, d12)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey80") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
  # THE EXPECTED PATTERN: the Myc effect rescaled by the global rate of Fig. 1G.
  ggplot2::geom_abline(slope = RATE, intercept = 0, linewidth = 0.35,
                       colour = "grey35") +
  ggplot2::geom_segment(data = res, inherit.aes = FALSE,
                        ggplot2::aes(x = d6, xend = d6, y = d12, yend = y_exp),
                        linetype = "22", linewidth = 0.3, colour = "grey45") +
  ggplot2::geom_point(ggplot2::aes(shape = established), size = 1.6,
                      stroke = 0.35, colour = "grey15", fill = "grey15") +
  # The line is IMPOSED, not fitted, and a bare line through a scatter reads as a
  # regression -- so it carries its factor. 0.49 is Fig. 1G's global rescaling
  # rate, not something estimated from these seven points.
  ggplot2::annotate("text", x = 0.40, y = RATE * 0.40 + 0.045, label = "x 0.49",
                    hjust = 0, vjust = 0, size = 1.7, colour = "grey35") +
  ggplot2::geom_text(data = LAB, inherit.aes = FALSE,
                     ggplot2::aes(x = x, y = y, label = pro, hjust = hj, vjust = vj),
                     size = 1.7, colour = "grey15", fontface = "italic") +
  ggplot2::scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1),
                              breaks = c(TRUE, FALSE),
                              labels = c("established at 6W", "not established"),
                              name = NULL) +
  ggplot2::scale_x_continuous(limits = XR, labels = lab_signed,
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::scale_y_continuous(limits = YR, labels = lab_signed,
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::labs(x = "Myc effect at 6W  (log2 ratio to Bcl-xL)",
                y = "Myc effect at 12W") +
  ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 1.5))) +
  theme_panel(base_size = 6) +
  # Key inside, top left: the expected line runs up to the right, so the wedge
  # above it on the left is where nothing can sit.
  ggplot2::theme(
    legend.position        = "inside",
    legend.position.inside = c(0.005, 0.99),
    legend.justification   = c(0, 1),
    legend.background      = ggplot2::element_blank(),
    legend.margin          = ggplot2::margin(0, 0, 0, 0),
    legend.key.size        = ggplot2::unit(2.4, "mm"),
    legend.spacing.y       = ggplot2::unit(0.3, "mm"),
    plot.margin            = ggplot2::margin(1.5, 2.5, 1, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
row_of  <- function(x) d[d$pro == x, ]
mcl     <- pr[pr$anti == "Mcl1", ]
null_of <- function(x) pn[pn$pair == paste0(x, ":Bcl2l1"), ]
ret     <- function(g) mem(g, "myc_12W") / mem(g, "myc_6W")

LEGEND <- panel_legend(
  slot = "Fig. 2G",
  what = paste0(
    "Seven pro-apoptotic transcripts as a log2 ratio to Bcl-xL: the Myc effect ",
    "on each ratio at twelve weeks against the same effect at six. The line is ",
    "the expected pattern -- the Myc effect rescaled by the global rate of Fig. ",
    "1G -- and the dashed drops are each ratio's distance from it. Filled points ",
    "are the ratios Myc established at six weeks."),
  detail = c(
    sprintf("n = 6 per group. Each ratio is fitted per sample as log2(pro) - log2(Bcl-xL) and the Myc effect is the genotype coefficient of `ratio ~ timepoint * genotype` within each age. The expected line has slope %.3f, the global rescaling rate of the Myc effect fitted over the 2,648 genes Myc moves (script 44; Fig. 1G quotes it as 0.487, against 0.552 for the mitochondrial compartment).",
            RATE),
    sprintf("THE THREE THAT FOLLOW THE RESCALING: %s %+.3f to %+.3f (retention %.2f), %s %+.3f to %+.3f (%.2f), %s %+.3f to %+.3f (%.2f). All three stay positive -- Myc still raises them at twelve weeks, by about half as much, which is what every other Myc effect in this dataset does.",
            "Bax",  row_of("Bax")$d6,  row_of("Bax")$d12,  row_of("Bax")$retention,
            "Bid",  row_of("Bid")$d6,  row_of("Bid")$d12,  row_of("Bid")$retention,
            "Bak1", row_of("Bak1")$d6, row_of("Bak1")$d12, row_of("Bak1")$retention),
    sprintf("AND THE ONE THAT DOES NOT: %s %+.3f to %+.3f, retention %.2f -- it crosses zero, so at twelve weeks Myc no longer raises the PUMA:Bcl-xL ratio at all. The interaction is %+.3f (p = %.3f), and adjusted for epithelial and immune content %+.3f (p = %.3f).",
            "Bbc3", row_of("Bbc3")$d6, row_of("Bbc3")$d12, row_of("Bbc3")$retention,
            row_of("Bbc3")$int, row_of("Bbc3")$int_p,
            pr$int_adj[pr$pair == "Bbc3:Bcl2l1"], pr$int_p_adj[pr$pair == "Bbc3:Bcl2l1"]),
    sprintf("WHY THE OTHER RATIOS HOLD, WHICH IS THE SENTENCE'S CLAUSE ABOUT THE MEMBERS: the two large members decline together at the global rate, so their ratio is left where Myc put it. Bax retains %.2f of its six-week effect and Bcl-xL %.2f, against the global %.2f, and their ratio retains %.2f. Bbc3 retains %.2f against the same denominator, which is the whole of the difference.",
            ret("Bax"), ret("Bcl2l1"), RATE, row_of("Bax")$retention, ret("Bbc3")),
    sprintf("THE DENOMINATOR IS SHARED, and that is what makes the comparison internal rather than a claim against an outside null: every ratio here is against Bcl-xL, so \"it is just the global attenuation\" is refuted from inside the panel. Two further pairs use Mcl-1 and are not drawn: %s %+.3f to %+.3f (retention %.2f, established at six weeks, p = %.4f) and %s %+.3f to %+.3f (%.2f, NOT established, p = %.2f). PUMA reverses against both denominators, but only the Bcl-xL one was there to begin with.",
            mcl$pair[1], mcl$d6[1], mcl$d12[1], mcl$retention[1], mcl$p6[1],
            mcl$pair[2], mcl$d6[2], mcl$d12[2], mcl$retention[2], mcl$p6[2]),
    sprintf("THE OPEN POINTS ARE NOT WEAK RESULTS, THEY ARE ABSENT ONES. Myc did not establish those three ratios at six weeks (%s p = %.2f, %s p = %.2f, %s p = %.2f), and retention is a quotient, so their positions carry no information about loss -- Pmaip1's retention is %.1f and Bmf's %.2f purely because the denominators are near zero. No residual is drawn for them.",
            "Bcl2l11", row_of("Bcl2l11")$p6, "Bmf", row_of("Bmf")$p6,
            "Pmaip1", row_of("Pmaip1")$p6,
            row_of("Pmaip1")$retention, row_of("Bmf")$retention)),
  bounds = c(
    sprintf("\"SIGNIFICANTLY\" NEEDS ITS SCOPE. The PUMA:Bcl-xL interaction is nominally significant (p = %.3f raw, %.3f purity-adjusted), but across the nine pairs its Benjamini-Hochberg value is %.2f, and against script 42's matched-pair null -- random pro-like/anti-like pairs matched on expression, conditioned on a six-week effect at least as large -- the empirical p is %.3f (the %.0fth percentile of %d matched pairs). It is a nominal result with a pre-specified licence, not a multiplicity-surviving one.",
            row_of("Bbc3")$int_p, pr$int_p_adj[pr$pair == "Bbc3:Bcl2l1"],
            pr$int_p_bh[pr$pair == "Bbc3:Bcl2l1"], null_of("Bbc3")$p_emp_cond,
            null_of("Bbc3")$pct_retention_cond, null_of("Bbc3")$n_conditional),
    sprintf("THE LICENCE IS PRE-SPECIFICATION, and the text should say so once: PUMA was named in advance from the PGC1a cell experiments, not selected from this panel. The honest counterweight is that on the empirical null the most extreme pair is not PUMA but %s (p = %.4f) -- and Myc never established that ratio (p = %.2f), which is exactly why the conditional null is the one to read.",
            "Bmf", null_of("Bmf")$p_emp_cond, row_of("Bmf")$p6),
    "PRIMING IS NOT DEATH. These are transcript ratios; apoptotic priming is a property of the protein complement and of how close the mitochondrion sits to the threshold. BH3 profiling is the measurement, and Fig. S2C is the reason to do it.",
    "The gene-level interaction for `Bbc3` itself (p = 0.0081) is a DIFFERENT statistic from the ratio interaction quoted here (p = 0.037 raw). Fig. 2H is where the gene-level one belongs; the two must not be interchanged.",
    "Genotype contrasts are clean -- genotype is balanced within each extraction batch -- so this panel does not carry the batch caveat that the wild-type temporal panels do. What it does carry is n = 6 per cell: every interaction on it is a difference of two 6-versus-6 contrasts, which is the least powered thing in the design (median lfcSE 0.333 against 0.233 for a genotype effect)."),
  source = c(
    "results/priming_arm_teb.rds (scripts/42_priming_arm_and_teb_substrate.R) -- $priming for the nine pro:anti ratios and their interactions, $pair_null for the matched-pair empirical null",
    "results/collapse_module_ownership.rds (scripts/44_collapse_module_and_ownership.R) -- $defs$global_rate_fitted, the expected line; $wt_genes for the member-level retentions"))

save_panel_p(p, "fig2_priming_ratios", height = 58)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## all nine pairs, drawn or not, with both interaction fits
  pr[, c("pair", "d6", "p6", "d12", "retention", "int", "int_p", "int_p_adj",
         "int_p_bh")] |> print(row.names = FALSE, digits = 3)

  ## the matched-pair null, which is the honest test and does not reach 0.05
  pn |> print(row.names = FALSE, digits = 3)

  ## the members behind each ratio -- the sentence's clause, gene by gene. Note
  ## that a member retention is only interpretable where the six-week effect is
  ## large: Bak1 and Pmaip1 have small ones and their quotients are unstable.
  g <- c("Bax", "Bak1", "Bid", "Bbc3", "Bcl2l11", "Pmaip1", "Bcl2l1", "Mcl1")
  data.frame(gene = g, baseMean = vapply(g, mem, numeric(1), "baseMean"),
             myc_6W = vapply(g, mem, numeric(1), "myc_6W"),
             myc_12W = vapply(g, mem, numeric(1), "myc_12W"),
             retention = vapply(g, ret, numeric(1))) |>
    print(row.names = FALSE, digits = 3)

  ## the drawn ratios against the expected line, as residuals
  data.frame(pro = d$pro, d6 = d$d6, d12 = d$d12, expected = RATE * d$d6,
             residual = d$d12 - RATE * d$d6, established = d$established) |>
    (\(x) x[order(x$residual), ])() |> print(row.names = FALSE, digits = 3)
}
