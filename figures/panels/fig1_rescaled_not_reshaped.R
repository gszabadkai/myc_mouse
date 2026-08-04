# =============================================================================
# fig1_rescaled_not_reshaped.R -- the twelve-week Myc effect IS the six-week Myc
# effect times a constant
# -----------------------------------------------------------------------------
# SLOT: Fig. 1G.
#
#   "At 12W weeks the entire structure is present at half amplitude and unchanged
#    in shape both in the whole and mitochondrial transcriptome (overall R2 = ...,
#    OXPHOS R2 = ..., p < ..., Fig. 1G)"
#
#   LEFT    per gene, the 2,648 genes Myc moves at six weeks
#   RIGHT   per pathway, the 143 nuclear-encoded MitoPathways, OXPHOS marked by
#           colour only, with the key inside the frame
#
# Each facet carries the slope and R2 of its own line and nothing else.
#
# WHAT IS DELIBERATELY NOT DRAWN (author, 2026-08-04): the OXPHOS fit line and its
# numbers, and the null. All three go in the text. They are still computed here and reported in
# the legend block, because the sentence needs them and because one of them is the
# distinction the sentence turns on:
#
#     slope 0.552   64 of 500 expression-matched shuffles reach it -> p = 0.13
#     R2    0.801    0 of 500 reach it (null max 0.719)            -> p < 0.002
#
# HALF AMPLITUDE is the slope; UNCHANGED IN SHAPE is the R2. A shuffled set of the
# same size and expression can halve an effect. What it cannot do is halve it
# COHERENTLY -- so the p belongs to the shape claim and not to the amplitude one.
#
# WHY THE LEFT FACET IS 2,648 GENES AND NOT 8,774. The rate of record, script 44's
# `global_rate_fitted` = 0.487, is fitted on exactly this set (padj < 0.1 at six
# weeks, |LFC| >= 0.2, baseMean >= 20) and the filter is reproduced and asserted
# here. Over all 8,774 reported genes the slope barely moves (0.450) but the R2
# falls to 0.381, because six thousand genes Myc does not move contribute noise to
# both axes -- regression dilution, not a change in the result. The panel draws the
# set the number belongs to and the strip says so; the all-gene figures are in the
# legend block. The manuscript's "whole transcriptome" should say which.
#
# Reads (read-only, no re-run):
#   results/collapse_module_ownership.rds (script 44) -- $collapse_genes, $defs
#   results/background_vs_myc.rds (script 40) -- $ruler, $regressions,
#       $regression_boot, $regression_null, $null_draws
# Output: outputs/figures/panels/fig1_rescaled_not_reshaped.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

cm  <- readRDS(here::here("results", "collapse_module_ownership.rds"))
bg  <- readRDS(here::here("results", "background_vs_myc.rds"))

cg  <- as.data.frame(cm$collapse_genes)          # scripts 40 and 44 save tibbles
dfs <- cm$defs
ruler <- as.data.frame(bg$ruler)
reg   <- as.data.frame(bg$regressions)
boot  <- as.data.frame(bg$regression_boot)
rnull <- as.data.frame(bg$regression_null)

# --- the gene-level universe, reproduced from its own definition --------------
in_fit <- !is.na(cg$padj_6W) & cg$padj_6W < dfs$padj6 &
  abs(cg$lfc_6W) >= dfs$lfc6_floor & cg$baseMean >= dfs$basemean_floor
stopifnot(identical(as.logical(in_fit), as.logical(cg$in_ranking)),
          sum(in_fit) == dfs$n_ranking_genes,
          nrow(cg) == dfs$n_reported_genes)

gsub_ <- cg[in_fit, c("gene", "lfc_6W", "lfc_12W")]
# THROUGH THE ORIGIN, because the quantity is a multiplier: "the 12W effect is
# 0.49 times the 6W effect" is a statement with no intercept in it, and this is
# the fit script 44's rate of record comes from. Freeing the intercept gives 0.086
# and a slope of 0.455 -- both in the legend block.
g_fit <- stats::lm(lfc_12W ~ 0 + lfc_6W, data = gsub_)
g_all <- stats::lm(lfc_12W ~ 0 + lfc_6W, data = cg)
# R2 as the SQUARED PEARSON CORRELATION throughout, which is what
# summary(lm)$r.squared gives for a fit WITH an intercept and is therefore the
# same definition the pathway facet uses. A through-origin lm reports an uncentred
# R2 instead (0.565 rather than 0.513 here), which is not comparable to anything
# else on the panel and must not be the number quoted.
g_r2  <- function(d) stats::cor(d$lfc_6W, d$lfc_12W)^2

# ASSERTION: the drawn slope IS script 44's rate of record.
stopifnot(abs(unname(stats::coef(g_fit)) - dfs$global_rate_fitted) < 1e-6)

# --- the pathway-level fits, taken from script 40 and reproduced --------------
r143 <- ruler[!ruler$is_mtdna, ]
ox   <- r143[r143$tier == "OXPHOS", ]
p_all <- stats::lm(c_m12 ~ c_m6, data = r143)
p_ox  <- stats::lm(c_m12 ~ c_m6, data = ox)
reg_row <- reg[reg$model == "rescale (Myc@12W ~ Myc@6W)" &
                 reg$scope == "content, all", ]

stopifnot(nrow(r143) == 143L, nrow(ox) == 19L, nrow(reg_row) == 1L,
          abs(unname(stats::coef(p_all)[2]) - reg_row$slope) < 1e-6,
          abs(summary(p_all)$r.squared      - reg_row$r2)    < 1e-6)

# =============================================================================
# the two scatters
# =============================================================================
# One ggplot with two facets rather than two plots: the facet strip is where the
# universe of each panel is named, which is what removes the need for a key on the
# grey points. The only key is the OXPHOS highlight.
OX_COL <- unname(ms_diverging[["pos"]])          # declared; not a sample colour

# DISPLAY WINDOWS. A handful of genes reach |LFC| 7-9 and two pathways reach +1.9
# (glycine cleavage, four genes); left in frame they compress the cloud that
# carries the fit into a smear. Each facet is drawn over a window and the points
# outside it are not drawn -- but EVERY NUMBER ON THE PANEL, and every fitted line,
# is computed on the complete set. The counts left out are in the legend block, as
# Fig. 1E states its own capped point.
GLIM <- 3.0                                       # genes, symmetric
PHI  <- 1.05                                      # pathways, upper
PLO  <- min(r143$c_m6, r143$c_m12)

fac <- c(gene = sprintf("Myc-responsive genes (%s)", format(nrow(gsub_), big.mark = ",")),
         path = "MitoPathways (143)")

pts <- rbind(
  data.frame(facet = fac[["gene"]], x = gsub_$lfc_6W, y = gsub_$lfc_12W,  cls = "bg"),
  data.frame(facet = fac[["path"]], x = r143$c_m6,    y = r143$c_m12,     cls = "bg"),
  data.frame(facet = fac[["path"]], x = ox$c_m6,      y = ox$c_m12,       cls = "ox"))
pts$facet <- factor(pts$facet, levels = fac)

# per-facet drawing frames
lims <- data.frame(
  facet = factor(fac, levels = fac),
  lo = c(-GLIM, PLO),
  hi = c( GLIM, PHI))

n_out_g <- sum(pmax(abs(gsub_$lfc_6W), abs(gsub_$lfc_12W)) > GLIM)
n_out_p <- sum(pmax(r143$c_m6, r143$c_m12) > PHI)
inside  <- (pts$facet == fac[["gene"]] & pmax(abs(pts$x), abs(pts$y)) <= GLIM) |
           (pts$facet == fac[["path"]] & pmax(pts$x, pts$y) <= PHI)
pts <- pts[inside, ]
# ONE fitted line per facet (author, 2026-08-04): the whole set. The OXPHOS
# pathways are marked by colour and nothing else -- a second line at slope 0.525
# beside one at 0.552 is two lines saying they are the same, which the points
# already say more directly, and their numbers go in the text.
fits <- rbind(
  data.frame(facet = fac[["gene"]], slope = unname(stats::coef(g_fit)), int = 0),
  data.frame(facet = fac[["path"]], slope = unname(stats::coef(p_all)[2]),
             int = unname(stats::coef(p_all)[1])))
fits$facet <- factor(fits$facet, levels = fac)

# Each facet carries the slope and R2 of its own line. Two annotations rather than
# one two-line string, because the R2 is written with a real superscript: plotmath
# keeps the source ASCII (`R^2`, per the coding rules) and renders it properly,
# which a literal character would not do reliably through the base pdf() device.
# The numbers are QUOTED inside the plotmath string. Unquoted, `R^2~0.80` is
# parsed as a number and printed as 0.8 -- the trailing zero, which is the
# significant figure the reader is being given, silently disappears.
ann <- rbind(
  data.frame(facet = fac[["gene"]], x = -GLIM,
             y = GLIM - c(0, 1) * 2 * GLIM * 0.075,
             lab = c(sprintf('slope~"%.2f"', stats::coef(g_fit)),
                     sprintf('R^2~"%.2f"',   g_r2(gsub_)))),
  data.frame(facet = fac[["path"]], x = PLO,
             y = PHI - c(0, 1) * (PHI - PLO) * 0.075,
             lab = c(sprintf('slope~"%.2f"', stats::coef(p_all)[2]),
                     sprintf('R^2~"%.2f"',   summary(p_all)$r.squared))),
  stringsAsFactors = FALSE)
ann$facet <- factor(ann$facet, levels = fac)

cls_cols <- c(bg = "grey35", ox = OX_COL)

p <- ggplot2::ggplot(pts, ggplot2::aes(x, y)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "22",
                       linewidth = 0.25, colour = "grey65") +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey85") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") +
  ggplot2::geom_point(data = pts[pts$cls == "bg", ], colour = "grey55",
                      size = 0.35, alpha = 0.35, stroke = 0) +
  ggplot2::geom_point(data = pts[pts$cls == "ox", ], ggplot2::aes(colour = cls),
                      size = 0.9, alpha = 0.95, stroke = 0) +
  ggplot2::geom_abline(data = fits,
                       ggplot2::aes(slope = slope, intercept = int),
                       linewidth = 0.45, colour = "grey20") +
  ggplot2::geom_text(data = ann, ggplot2::aes(x = x, y = y, label = lab),
                     parse = TRUE, hjust = 0, vjust = 1, size = 1.8,
                     colour = "grey25", inherit.aes = FALSE) +
  ggplot2::geom_blank(data = rbind(
    data.frame(facet = lims$facet, x = lims$lo, y = lims$lo),
    data.frame(facet = lims$facet, x = lims$hi, y = lims$hi))) +
  ggplot2::facet_wrap(~ facet, nrow = 1, scales = "free") +
  ggplot2::scale_colour_manual(values = cls_cols, breaks = "ox",
                               labels = c(ox = "OXPHOS"), name = NULL) +
  ggplot2::scale_x_continuous(labels = lab_signed) +
  ggplot2::scale_y_continuous(labels = lab_signed) +
  ggplot2::labs(x = "Myc effect at 6W (log2)", y = "Myc effect at 12W (log2)") +
  ggplot2::guides(colour = ggplot2::guide_legend(
    override.aes = list(size = 1.5, alpha = 1))) +
  theme_panel(base_size = 6) +
  # The one-item key sits INSIDE the right facet, bottom right (author,
  # 2026-08-04), which is the wedge the diagonal cloud leaves empty -- high six-week
  # effect with a low twelve-week one is exactly what does not happen. Under the
  # plot it cost a whole row of height for one dot.
  ggplot2::theme(
    strip.text      = ggplot2::element_text(face = "plain", size = 5.6,
                                            margin = ggplot2::margin(0, 0, 0.6, 0, "mm")),
    panel.spacing.x = ggplot2::unit(2.4, "mm"),
    legend.position        = "inside",
    legend.position.inside = c(0.995, 0.02),
    legend.justification   = c(1, 0),
    legend.background      = ggplot2::element_blank(),
    legend.margin   = ggplot2::margin(0, 0, 0, 0),
    legend.key.size = ggplot2::unit(2.4, "mm"),
    plot.margin     = ggplot2::margin(1, 1.5, 0.5, 1.5, "mm"))

# =============================================================================
# the null -- COMPUTED, NOT DRAWN
# =============================================================================
# It is the statistic that makes "unchanged in shape" a claim rather than a
# description, but the author's call (2026-08-04) is that a null distribution on
# the page is one significance level's worth of ink for a number the text can
# carry in five words. So it is computed, asserted, and reported in the legend
# block, and the panel stays a scatter. The drawn form is kept in the sandbox.
#
# 500 label shuffles WITHIN expression deciles, which preserve set size, expression
# and the pathway OVERLAP structure (MitoPathways nest, and the nesting alone
# creates cross-pathway correlation). Script 40's PART B; the resampled secondary
# null is in the `_rs` columns and agrees.
nd     <- bg$null_draws
r2_obs <- reg_row$r2
r2_nul <- as.numeric(nd["rescale_r2", ])
sl_nul <- as.numeric(nd["rescale_slope", ])
p_r2   <- (sum(r2_nul >= r2_obs) + 1) / (length(r2_nul) + 1)
p_sl   <- (sum(sl_nul >= reg_row$slope) + 1) / (length(sl_nul) + 1)

stopifnot(abs(stats::median(r2_nul) - rnull$null_median[rnull$statistic == "rescale_r2"]) < 1e-6,
          sum(r2_nul >= r2_obs) == 0L)

# =============================================================================
# the legend text (never drawn)
# =============================================================================
bt <- boot[boot$model == "rescale content", ]

LEGEND <- panel_legend(
  slot = "Fig. 1G",
  what = paste0(
    "The twelve-week Myc effect is the six-week Myc effect at about half ",
    "amplitude and in the same shape. LEFT: per gene, over the genes Myc moves ",
    "at six weeks. RIGHT: per MitoPathway, over the 143 nuclear-encoded ",
    "MitoPathways, with the nineteen OXPHOS pathways marked and fitted ",
    "separately. BOTTOM: the R2 of the right-hand regression against 500 ",
    "expression-matched shuffled sets."),
  detail = c(
    "Both axes are RAW, unshrunken DESeq2 log2 fold changes of the Myc genotype contrast (Myc+ minus wild type) at each age, which is what CLAUDE.md requires of a delta-LFC visual. Both contrasts are clean: genotype is balanced within each extraction batch. The dashed line is y = x, so distance below it is the amplitude that was lost.",
    sprintf("LEFT: %s genes, being those with padj < %.1f at six weeks, |log2FC| >= %.1f and baseMean >= %d -- script 44's own filter, reproduced and asserted here. The fitted slope through the origin is %.3f, which is script 44's global rate of record, and R2 is %.3f.",
            format(nrow(gsub_), big.mark = ","), dfs$padj6, dfs$lfc6_floor,
            dfs$basemean_floor, stats::coef(g_fit), g_r2(gsub_)),
    sprintf("Over ALL %s reported genes the slope barely moves, %.3f, but R2 falls to %.3f: the six thousand genes Myc does not move add noise to both axes and dilute the correlation without changing the rate. The panel draws the set the rate belongs to.",
            format(nrow(cg), big.mark = ","), stats::coef(g_all), g_r2(cg)),
    sprintf("RIGHT: 143 pathways, slope %.3f with an intercept of %.4f -- indistinguishable from zero, so this is a pure rescaling and not a shift -- and R2 %.3f. Bootstrap 95%% interval on the slope %.3f to %.3f. Within the OXPHOS tier alone (19 pathways) the slope is %.3f and R2 %.3f.",
            reg_row$slope, reg_row$intercept, reg_row$r2, bt$lo, bt$hi,
            stats::coef(p_ox)[2], summary(p_ox)$r.squared),
    sprintf("BOTTOM: the null is 500 label shuffles WITHIN expression deciles, which preserve set size, expression level and the nesting structure of the MitoPathways. Its median R2 is %.3f and its maximum is %.3f; the observed %.3f is above every draw, giving an empirical p < %.3f.",
            stats::median(r2_nul), max(r2_nul), r2_obs, p_r2),
    sprintf("The SLOPE is not the extreme statistic and should not be quoted as though it were: %d of the 500 shuffles reach 0.552 or more (p = %.2f). A random collection of genes of the same size and expression can halve an effect. What it cannot do is halve it coherently, which is the R2.",
            sum(sl_nul >= reg_row$slope), p_sl),
    "On the content-blind priority ruler the same regression gives slope 0.644 and R2 0.789, and within the OXPHOS tier R2 0.945 -- drawn in Fig. 1F as the connector on every row, which is what the sentence's back-reference points at."),
  bounds = c(
    sprintf("BOTH facets are WINDOWED for display and the points outside are not drawn: +/-%.1f log2 on the left, which leaves out %d of the %s genes, and an upper bound of %+.2f on the right, which leaves out %d of the 143 pathways (Glycine metabolism at %+.2f and the four-gene Glycine cleavage system at %+.2f). Every fitted line and every number on the panel is computed on the COMPLETE set.",
            GLIM, n_out_g, format(nrow(gsub_), big.mark = ","), PHI, n_out_p,
            max(r143$c_m6[r143$pathway == "Glycine metabolism"]),
            max(r143$c_m6[r143$pathway == "Glycine cleavage system"])),
    "The empty band down the middle of the left facet is the filter, not a gap in the data: genes with |log2FC| below 0.2 at six weeks are not in the fitted set.",
    "The two facets are not on the same footing and must not be read as a replication. The left is per gene and the right is per pathway; a pathway average of 100 genes is far less noisy than one gene, which is most of why its R2 is higher. What the two share is the SLOPE, and that is the comparison the panel supports.",
    "\"Unchanged in shape\" is an R2 statement about a linear rescaling. It does not say that no pathway changed rank -- Fig. 1F shows several that did, most sharply pyruvate metabolism (+0.328 to +0.006) -- only that a single multiplier explains four fifths of the twelve-week profile.",
    "The MitoPathways NEST inside one another, so the 143 points are not independent. That is precisely what the null is built to absorb: the shuffle is within expression deciles and preserves the overlap structure, so the 100th-percentile result is not an artefact of the nesting.",
    "The genotype contrast at each age is clean; the comparison BETWEEN the two ages made here is a comparison of two clean contrasts, so it is not exposed to batch = timepoint. The interpretation of the slope as a Myc DOSE effect rests on the western blot (Fig. S1D), not on this panel.",
    "The synthetic mtDNA-encoded pathway is excluded, as it is from every fit and null in script 40: it is a 4 to 6.5 SD outlier on the temporal contrasts and alone moves the shared-vector slope from 0.75 to 0.84."),
  source = c(
    "results/background_vs_myc.rds (scripts/40_background_vs_myc_decomposition.R) -- $ruler, $regressions, $regression_boot, $regression_null, $null_draws",
    "results/collapse_module_ownership.rds (scripts/44) -- $collapse_genes, $defs$global_rate_fitted",
    "Earlier double-column form: figures/fig03_background_vs_myc.R panel B"))

save_panel_p(p, "fig1_rescaled_not_reshaped", height = 46)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## every rescale fit script 40 reports, by scope
  reg[reg$model == "rescale (Myc@12W ~ Myc@6W)", ] |> print(row.names = FALSE, digits = 3)

  ## the gene-level fit on both universes, and with an intercept
  rbind(
    data.frame(universe = "Myc-responsive (2,648)", n = nrow(gsub_),
               slope = stats::coef(g_fit), r2 = g_r2(gsub_)),
    data.frame(universe = "all reported (8,774)",  n = nrow(cg),
               slope = stats::coef(g_all), r2 = g_r2(cg))) |>
    print(row.names = FALSE, digits = 3)

  ## the two nulls side by side -- only the R2 is beyond its own
  data.frame(statistic = c("slope", "R2"),
             observed  = c(reg_row$slope, r2_obs),
             null_med  = c(stats::median(sl_nul), stats::median(r2_nul)),
             null_max  = c(max(sl_nul), max(r2_nul)),
             p_emp     = c(p_sl, p_r2)) |> print(row.names = FALSE, digits = 3)

  ## the null, drawn -- the form that was on the panel until 2026-08-04
  {
    dn2 <- stats::density(r2_nul, adjust = 0.9)
    ggplot2::ggplot(data.frame(r2 = r2_nul), ggplot2::aes(x = r2)) +
      ggplot2::geom_density(adjust = 0.9, fill = "grey88", colour = "grey40",
                            linewidth = 0.3) +
      ggplot2::annotate("segment", x = r2_obs, xend = r2_obs, y = 0,
                        yend = max(dn2$y) * 1.10, linewidth = 0.5, colour = "grey10") +
      ggplot2::annotate("point", x = r2_obs, y = max(dn2$y) * 1.10, size = 1.1,
                        shape = 25, fill = "grey10", colour = "grey10") +
      ggplot2::annotate("text", x = r2_obs, y = max(dn2$y) * 1.10,
                        label = sprintf("observed %.2f  ", r2_obs),
                        hjust = 1, vjust = 0.35, size = 1.7, colour = "grey10") +
      ggplot2::scale_x_continuous(limits = c(0, r2_obs * 1.06)) +
      ggplot2::labs(x = "R2 on expression-matched shuffled sets", y = NULL) +
      theme_panel(base_size = 6) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.line.y = ggplot2::element_blank())
  }

  ## which pathways sit furthest from the fitted line -- the reshaping that IS there
  r143$resid <- stats::resid(p_all)
  r143[order(r143$resid), c("pathway", "tier", "c_m6", "c_m12", "resid")] |>
    head(6) |> print(row.names = FALSE, digits = 3)
  r143[order(-r143$resid), c("pathway", "tier", "c_m6", "c_m12", "resid")] |>
    head(6) |> print(row.names = FALSE, digits = 3)
}
