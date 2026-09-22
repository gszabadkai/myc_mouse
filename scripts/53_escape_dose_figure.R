# =============================================================================
# 53 -- THE ESCAPE SERIES AS A DOSE FIGURE
# -----------------------------------------------------------------------------
# Reads results/orthotopic_escape_series.rds and draws it. NO ANALYSIS: every
# value is re-read from that object and nothing is transcribed or recomputed. If
# a number here disagrees with `docs/2026-09-10_orthotopic_escape_series.md`, the
# object is right and both of them are wrong.
#
# THE READING THE FIGURE HAS TO CARRY, AND THE ONE A READER WILL GET WRONG.
# The x-axis is `Ppargc1a` RETAINED, which is an OUTCOME OF SELECTION and not an
# assigned dose. Nobody titrated PGC1a. PGC1a kills MYAZ cells, so each arm's
# tumours are the ones that survived expressing it, and how much they kept is set
# by how much of the death pathway they had already broken -- transgene lost
# entirely (KOPgc1a), a partial dose (NTPgc1a), the full dose under a downstream
# buffer (Pgc1a_BclxL). A reader will assume a dose-response experiment. It is a
# survivorship series read backwards, and the caption says so.
#
# EVERY POINT IS AN EFFECT AGAINST ITS OWN MATCHED CONTROL, WITHIN SERIES:
# Pgc1a_BclxL vs BclxL, NTPgc1a vs NTEV, KOPgc1a vs KOEV. The vector arms are a
# polyclonal pool and the CRISPR arms are clones, so nothing is ever compared
# across series -- three backgrounds, one manipulation, each with its own control.
#
# PANELS
#   A  the three respiratory rulers as a percentage of the full-dose arm's effect,
#      one line each, over a strip carrying the MitoCarta SHARE difference -- the
#      mass readout, and the thing that separates the full-dose arm from the rest.
#   B  `Bbc3`, log2 difference with its bootstrap interval.
#   C  `Myc` 80th-percentile contrast on BOTH normalisations side by side. CPM
#      deflation is the obvious objection to C2 and the answer is that both
#      normalisations agree, so both are drawn rather than one.
#
# Reads : results/orthotopic_escape_series.rds (script 52)
# Output: outputs/orthotopic/53_escape_dose.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

res_path <- here::here("results", "orthotopic_escape_series.rds")
if (!file.exists(res_path)) stop("run scripts/52_orthotopic_escape_series.R first")
res <- readRDS(res_path)

# --- the three PGC1a arms, ordered by the dose they RETAINED ------------------
W  <- as.data.frame(res$within_series)          # hi / lo / series / status
ET <- as.data.frame(res$escape_table)
dose <- stats::setNames(ET$Ppargc1a_CPM, as.character(ET$arm))[W$hi]
W$dose <- as.numeric(dose)
W <- W[order(W$dose), ]
W$lab <- sprintf("%s\n%.3g CPM", W$hi, W$dose)
W$lab <- factor(W$lab, levels = W$lab)          # ordered by retained dose
FULL  <- W$hi[which.max(W$dose)]                # the full-dose arm, from the data
stopifnot(nrow(W) == 3L, !anyNA(W$dose))

lab_of <- function(hi) W$lab[match(hi, W$hi)]

# =============================================================================
# PANEL A -- the respiratory rulers, as a percentage of the full-dose arm
# -----------------------------------------------------------------------------
# The three normalisations of the same respiratory question: the absolute level,
# the within-compartment share, and mitoPPS's pairwise ratio. Each arm's effect is
# expressed against the FULL-DOSE arm's effect on the same ruler, so the three
# rulers are comparable on one axis; the full-dose arm is 100% by construction.
# =============================================================================
RUL <- c(ox_lvl = "absolute level", ox_rel = "compartment share",
         mitopps_oxphos = "mitoPPS (pairwise ratio)")
RU  <- as.data.frame(res$rulers)
ra  <- do.call(rbind, lapply(names(RUL), function(rn) {
  d <- RU[RU$ruler == rn & RU$hi %in% W$hi & RU$lo %in% W$lo, ]
  d <- d[match(W$hi, d$hi), ]
  full <- d$diff_median[d$hi == FULL]
  data.frame(ruler = RUL[[rn]], arm = d$hi, lab = lab_of(d$hi),
             pct = 100 * d$diff_median / full,
             raw = d$diff_median, stringsAsFactors = FALSE) }))
ra$ruler <- factor(ra$ruler, levels = unname(RUL))

pA <- ggplot2::ggplot(ra, ggplot2::aes(lab, pct, colour = ruler, group = ruler)) +
  ggplot2::geom_hline(yintercept = c(0, 100), linewidth = 0.25,
                      colour = "grey70", linetype = c("solid", "dotted")) +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::scale_colour_manual(values = c("#1B7837", "#4393C3", "#B2182B"), name = NULL) +
  ggplot2::labs(x = NULL, y = sprintf("effect as %% of the %s arm", FULL)) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(legend.position = "top",
                 legend.margin = ggplot2::margin(0, 0, -4, 0),
                 panel.grid.minor = ggplot2::element_blank())

# --- the strip: the MASS readout, as a share difference against the control ---
DG <- as.data.frame(res$c2_diagnostic)
sh <- stats::setNames(DG$mitocarta_share, DG$arm)
ms <- data.frame(lab = lab_of(W$hi),
                 d_share = 100 * (sh[W$hi] - sh[W$lo]), stringsAsFactors = FALSE)
pAs <- ggplot2::ggplot(ms, ggplot2::aes(lab, d_share)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  ggplot2::geom_col(width = 0.45, fill = "grey35") +
  ggplot2::labs(x = NULL, y = "MitoCarta share\n(pp vs control)") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

# =============================================================================
# PANEL B -- Bbc3, log2 difference with its bootstrap interval
# =============================================================================
P5 <- as.data.frame(res$p53_panel)
bb <- P5[P5$gene == "Bbc3" & P5$hi %in% W$hi, ]
bb <- bb[match(W$hi, bb$hi), ]
bb$lab <- lab_of(bb$hi)
bb$excl <- bb$ci_hi < 0 | bb$ci_lo > 0

pB <- ggplot2::ggplot(bb, ggplot2::aes(lab, diff_median)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
                         width = 0.12, linewidth = 0.4) +
  ggplot2::geom_point(ggplot2::aes(shape = excl), size = 2.8, fill = "white") +
  ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
  ggplot2::labs(x = NULL, y = "Bbc3, log2 difference\nvs matched control") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

# =============================================================================
# PANEL C -- Myc 80th percentile, BOTH normalisations
# -----------------------------------------------------------------------------
# The units differ between scales by three orders of magnitude, so the facets
# carry free axes; what the panel is for is the SIGN and the p, side by side.
# =============================================================================
C2 <- as.data.frame(res$c2)
mc <- C2[C2$statistic == "pct80" & C2$hi %in% W$hi, ]
mc$lab   <- lab_of(mc$hi)
mc$scale <- factor(mc$scale, levels = c("CPM", "median-of-ratios"))
mc$sig   <- mc$p_one_sided_lower < 0.05

pC <- ggplot2::ggplot(mc, ggplot2::aes(lab, stat_obs)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  ggplot2::geom_col(ggplot2::aes(fill = sig), width = 0.5) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("p %.3f", p_one_sided_lower),
                                  vjust = ifelse(stat_obs < 0, 1.4, -0.6)),
                     size = 2.5) +
  ggplot2::scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey65"),
                             guide = "none") +
  ggplot2::facet_wrap(~ scale, nrow = 1, scales = "free_y") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.18)) +
  ggplot2::labs(x = NULL, y = "Myc 80th-percentile contrast\nvs matched control") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

# =============================================================================
# ASSEMBLE
# =============================================================================
ko_max <- ET$Myc_max[ET$arm == "KOPgc1a"]
CAPTION <- paste0(
  "Every point is an effect against its own MATCHED CONTROL, within series: ",
  paste(sprintf("%s vs %s", W$hi, W$lo), collapse = "; "),
  ". The vector arms are a polyclonal pool and the CRISPR arms are clones, so ",
  "nothing is compared across series. x is Ppargc1a RETAINED -- an outcome of ",
  "selection, not an assigned dose: PGC1a kills MYAZ cells, so each arm is the ",
  "population that survived expressing it. (A) three normalisations of the same ",
  "respiratory question, each as a percentage of the ", FULL, " arm's effect on ",
  "that ruler; the strip is the MitoCarta share against control, in percentage ",
  "points, and is the mass readout. (B) filled points exclude zero. (C) both ",
  "normalisations, free axes -- the units differ, the sign and the p do not. ",
  "Myc MAX is not drawn: KOPgc1a carries an outlier at ", format(round(ko_max)),
  " CPM which makes its max contrast positive and uninterpretable. n = 7 per arm; ",
  "intervals are 5,000-resample bootstraps; the Myc null is exact over all 3432 ",
  "relabellings. Transcript associations throughout.")

# tags are assigned EXPLICITLY, with an empty one for the mass strip: it belongs
# to A and must not take "B" and push the caption's letters out of step
fig <- patchwork::wrap_plots(
  patchwork::wrap_plots(pA, pAs, ncol = 1, heights = c(3, 1)),
  patchwork::wrap_plots(pB, pC, ncol = 2, widths = c(1, 1.45)),
  ncol = 1, heights = c(1.35, 1)) +
  patchwork::plot_annotation(
    title = "The orthotopic escape series: what tumours that survived PGC1a expression look like",
    subtitle = "ordered by the PGC1a dose each arm RETAINED, which is an outcome of selection",
    caption = paste(strwrap(CAPTION, width = 150), collapse = "\n"),
    tag_levels = list(c("A", "", "B", "C"))) &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 10, face = "bold"),
                 plot.caption = ggplot2::element_text(size = 6.5, hjust = 0,
                                                      colour = "grey25"))

dir.create(here::here("outputs", "orthotopic"), showWarnings = FALSE, recursive = TRUE)
ggplot2::ggsave(here::here("outputs", "orthotopic", "53_escape_dose.pdf"),
                fig, width = 10, height = 8.5)
message("53: wrote outputs/orthotopic/53_escape_dose.pdf")

cat("\n=== what is drawn, re-read from the object ===\n")
print(ra, digits = 4, row.names = FALSE)
cat("\n"); print(ms, digits = 3, row.names = FALSE)
cat("\n"); print(bb[, c("lab", "diff_median", "ci_lo", "ci_hi")], digits = 3, row.names = FALSE)
cat("\n"); print(mc[, c("lab", "scale", "stat_obs", "p_one_sided_lower")], digits = 4, row.names = FALSE)

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(fig)

  ## the reading the figure exists to protect: x is RETAINED dose, an outcome of
  ## selection. Nobody titrated PGC1a.
  cat(CAPTION, "\n")

  ## panel A's numbers, and the raw effects behind the percentages
  ra |> print()
  ## the mass readout -- only the full-dose arm expands the compartment
  ms |> print()

  ## panel B -- Bbc3 in dose order
  bb[, c("hi", "lo", "diff_median", "ci_lo", "ci_hi")] |> print()

  ## panel C -- both normalisations agree in sign and in significance
  mc[, c("hi", "lo", "scale", "stat_obs", "p_one_sided_lower")] |> print()

  ## nothing here is computed: everything came from
  res$notes |> cat(sep = "\n")
}
