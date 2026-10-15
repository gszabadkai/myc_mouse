# =============================================================================
# fig1d_mito_content.R -- Fig. 1d: MYC raises mitochondrial content
# -----------------------------------------------------------------------------
# The text: "MYC raised mitochondrial content by 20-33%, depending on the module
# (Fig. 1d)", and, across the window, "the mitochondrial share of the transcriptome
# did not fall (Fig. 1d,f)". One PDF, two stacked parts, one shared
# group key:
#
#   TOP     each module's share of the nuclear transcriptome, every animal, four
#           groups, one bracket per module for the genotype main effect. The
#           mtDNA-encoded module is left out (author, 2026-10-10: not informative
#           here).
#   BOTTOM  the MYC effect on each of the 143 nuclear MitoPathways, one point per
#           pathway, at 6W and at 12W, with a thin tick at each age's median. The
#           points are in the MYC+ 6W and 12W colours of the shared key.
#
# FORM. Top: one value per animal, four groups of six -> every animal over a box
#   (median and interquartile range; whiskers minimum to maximum since 2026-10-14, as Prism)
#   (data-to-viz, "Do boxplots hide information?"). Each module has its own y axis, from
#   the tick closest below its lowest animal to the first tick above its bracket's P
#   label (author, 2026-10-15: a box plot's axis need not start at 0, "in order to
#   ensure resolution"). The x axis is not drawn, so the key names the groups. Groups
#   run genotype-major because the bracket is the genotype main effect.
#   Bottom: 143 values per age -> the points themselves (author, 2026-10-10: a
#   histogram added nothing the points do not show). Both rows have a single
#   peak, so a median is a fair summary here (FIGURE_RULES.md, Showing data).
#
# FORMAT (author, 2026-10-15): 100 x 60 mm (was 100 x 66); tick labels, the key and P values
#   5 pt, strip and axis titles 6 pt, as the Prism-format panel (Fig. 1b); the exponent
#   of a P value reduced (fmt_p's scriptstyle: 4 pt inside the 5-pt label).
#
# THE Y TITLE names the ruler: script 32's share_nomt leaves the 13 mtDNA-encoded
#   transcripts (7-16% of each library) out of the denominator, so the values are a
#   fraction of the nuclear-encoded transcriptome, not of the total. Script 32's
#   share_all is the fraction of the total; it moves the text's 20-33% to 21-33%.
#
# STATISTICS. Top: script 32's model, lm(log2 share ~ timepoint + myc_status); the
#   bracket P is its genotype main effect, read from results/mito_content_proxies.rds
#   $share_stats and asserted against a refit here. Bottom: descriptive; values are
#   script 40's set-average raw (unshrunken) log2 fold changes.
#
# Reads:  results/mito_content_proxies.rds (script 32: $shares, $share_stats)
#         results/background_vs_myc.rds    (script 40: $ruler)
# Output: outputs/natmetab/Fig1/Fig1d_mito_content.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

content <- readRDS(here::here("results", "mito_content_proxies.rds"))
bg      <- readRDS(here::here("results", "background_vs_myc.rds"))

TXT   <- 5      # pt; tick labels, the key and P values
TITLE <- 6      # pt; strip titles and axis titles

# =============================================================================
# TOP: the module shares
# =============================================================================
arms <- c("MASS_MARKERS_NOCHAP", "MITOCARTA_NUCLEAR_ENCODED",
          "MITOCARTA_OXPHOS_NU", "BIOGENESIS_FULL")
arm_name <- c(MASS_MARKERS_NOCHAP       = "Mass\nmarkers",
              MITOCARTA_NUCLEAR_ENCODED = "Nuclear\nMitoCarta",
              MITOCARTA_OXPHOS_NU       = "Nuclear\nOXPHOS",
              BIOGENESIS_FULL           = "Mitochondrial\nbiogenesis")

shares <- as.data.frame(content$shares)
stats  <- as.data.frame(content$share_stats)
stats  <- stats[stats$denominator == "share_nomt", ]
rownames(stats) <- stats$panel

df <- shares[shares$panel %in% arms, ]
df$panel <- factor(df$panel, levels = arms, labels = arm_name[arms])
df$group <- factor(df$group, levels = order_by_geno)
stopifnot(nrow(df) == length(arms) * 24L, !anyNA(df$group),
          all(table(df$panel, df$group) == 6L))

# the drawn contrast IS script 32's: refit its model and assert
refit <- vapply(arms, function(a) {
  d <- shares[shares$panel == a, ]
  d$tp <- factor(d$timepoint, levels = c("6W", "12W"))
  d$mc <- factor(d$myc_status, levels = c("neg", "pos"))
  co <- summary(stats::lm(log2(share_nomt) ~ tp + mc, data = d))$coefficients
  c(beta = co["mcpos", 1], p = co["mcpos", 4])
}, numeric(2))
stopifnot(max(abs(refit["beta", ] - stats[arms, "geno_beta"])) < 1e-8,
          max(abs(refit["p", ] - stats[arms, "geno_p"])) < 1e-8)

# THE TEXT'S NUMBERS (author, 2026-10-10): "MYC raised mitochondrial content by
# 20-33%, depending on the module". The range is the genotype effect across the
# four modules drawn, as whole percents; asserted here, where it is computed.
TEXT_CONTENT_PCT <- c(20, 33)
pct <- 100 * (2^stats[arms, "geno_beta"] - 1)
stopifnot(identical(round(range(pct)), TEXT_CONTENT_PCT))

# THE TEXT (author, 2026-10-11): in the WT gland "the mitochondrial share of the
# transcriptome did not fall" between 6 and 12 weeks -- the nuclear MitoCarta
# module's mean log2 share at 12W is not below that at 6W.
wt_mc <- shares[shares$panel == "MITOCARTA_NUCLEAR_ENCODED" & shares$myc_status == "neg", ]
stopifnot(nrow(wt_mc) == 12L)
wt_window <- stats::t.test(log2(share_nomt) ~ factor(timepoint, levels = c("12W", "6W")),
                           data = wt_mc)
stopifnot(diff(rev(wt_window$estimate)) >= 0)

# Per-module axis (author, 2026-10-15): a box plot's axis need not start at 0, so for
# resolution each runs from the tick closest below the lowest animal (the lower whisker:
# min to max) to the first tick above the bracket's P label, the highest mark
# (box_axis()). Worked in printed mm on graphs GRAPH_H mm tall: the bracket GAP mm
# above the top animal; the label 1% of the graph above the bracket (ggpubr) at vjust
# -0.15, so its top 1.15 label heights higher; then at least CLEAR mm to the axis end.
# A lowest animal close above its tick overhangs the axis end (Nuclear OXPHOS, 0.31 on
# an axis from 0.3), so the graphs do not clip.
GRAPH_H <- 28     # mm; the upper graphs' height, fixed in the stack below
GAP     <- 1.0    # mm; top animal to bracket
CLEAR   <- 0.2    # mm; top of the label to the axis end, at least
STEP    <- c(MASS_MARKERS_NOCHAP = 0.05, MITOCARTA_NUCLEAR_ENCODED = 1,
             MITOCARTA_OXPHOS_NU = 0.2, BIOGENESIS_FULL = 0.2)   # labelled ticks; minor ones halfway

lo_animal <- tapply(df$share_nomt, df$panel, min)
fac_max   <- tapply(df$share_nomt, df$panel, max)
stopifnot(identical(names(fac_max), unname(arm_name[arms])))
comps <- data.frame(panel = factor(names(fac_max), levels = levels(df$panel)),
                    group1 = 1.5, group2 = 3.5, p = stats[arms, "geno_p"])
need  <- GAP + 0.01 * GRAPH_H + 1.15 * label_mm(fmt_p(comps$p, exponent = "reduced"), TXT) + CLEAR
axes  <- lapply(seq_along(arms), function(i)
  box_axis(lo_animal[[i]], fac_max[[i]], STEP[[arms[i]]], GRAPH_H, above = need[i]))
lims  <- t(vapply(axes, function(a) a$limits, numeric(2)))   # one row per module: bottom, top
stopifnot(!anyDuplicated(lims[, 2]))                         # each facet's ticks are found by its top
pick   <- function(l) axes[[which.min(abs(lims[, 2] - l[2]))]]
per_mm <- (lims[, 2] - lims[, 1]) / GRAPH_H                  # data units per printed mm

br    <- brackets(comps, base = unname(fac_max) + GAP * per_mm,
                  step = (need - GAP) * per_mm, txt = TXT, exponent = "reduced")
frame <- data.frame(panel = factor(rep(names(fac_max), 2), levels = levels(df$panel)),
                    group = factor(order_by_geno[1], levels = order_by_geno),
                    y = c(lims[, 1], lims[, 2]))

p_top <- ggplot(df, aes(group, share_nomt)) +
  geom_blank(data = frame, aes(x = group, y = y), inherit.aes = FALSE) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.62, coef = Inf,   # min-max whiskers
               colour = box_line, linewidth = NM_LINE) +
  ggbeeswarm::geom_quasirandom(aes(colour = group), width = 0.22, size = 0.7, shape = 16) +
  br +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = group_cols, labels = group_keys, name = NULL) +
  scale_fill_manual(values = group_fill, guide = "none") +
  scale_y_continuous(breaks = function(l) pick(l)$breaks,
                     minor_breaks = function(l) pick(l)$minor_breaks,
                     guide = guide_axis(minor.ticks = TRUE), labels = lab_num, expand = c(0, 0)) +
  labs(x = NULL, y = "Fraction of nuclear­encoded\ntranscriptome (%)") +   # soft hyphen: pdf() prints "-" as a minus
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.5))) +
  coord_cartesian(clip = "off") +              # a point on the axis end is drawn whole
  theme_nm(legend = "bottom", txt = TXT, txt_title = TITLE) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x  = element_blank(),
        strip.text   = element_text(size = TITLE, lineheight = 0.95),
        panel.spacing.x = unit(2.5, "mm"),
        legend.margin = margin(0, 0, 0, 0))

# =============================================================================
# BOTTOM: the MYC effect on every MitoPathway, at both ages
# =============================================================================
ruler <- as.data.frame(bg$ruler)
r143  <- ruler[!ruler$is_mtdna, ]
stopifnot(nrow(r143) == 143L)

eff <- rbind(data.frame(contrast = "myc_6W",  lfc = r143$c_m6),
             data.frame(contrast = "myc_12W", lfc = r143$c_m12))
eff$contrast <- factor(eff$contrast, levels = names(contrast_cols))
eff$row      <- factor(eff$contrast, levels = rev(names(contrast_cols)))   # 6W on top

meds <- data.frame(contrast = factor(names(contrast_cols), levels = names(contrast_cols)),
                   med = c(stats::median(r143$c_m6), stats::median(r143$c_m12)))
meds$yi <- as.integer(factor(meds$contrast, levels = levels(eff$row)))

p_bot <- ggplot(eff, aes(lfc, row)) +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
  ggbeeswarm::geom_quasirandom(aes(colour = contrast), orientation = "y", width = 0.32,
                               size = 0.8, alpha = 0.65, stroke = 0) +
  geom_segment(data = meds, aes(x = med, xend = med, y = yi - 0.4, yend = yi + 0.4),
               inherit.aes = FALSE, linewidth = 0.4, colour = "grey15") +
  scale_colour_manual(values = contrast_cols, guide = "none") +
  scale_y_discrete(labels = c(myc_6W = "6W", myc_12W = "12W")) +
  scale_x_nm(step = 0.5, top = max(eff$lfc), bottom = min(eff$lfc), minor = TRUE,
             labels = lab_signed) +
  labs(x = "MYC+ vs WT, per MitoPathway (log2 fold change)", y = NULL) +
  theme_nm(txt = TXT, txt_title = TITLE) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# The upper graphs are GRAPH_H mm tall, the height their axis ends are worked on; the
# lower graph takes what is left of the page (author, 2026-10-15: 100 x 60 mm).
PAGE_W  <- 100
PAGE_H  <- 60
LOWER_H <- PAGE_H - natural_size(stack_parts(p_top, p_bot, heights = grid::unit(c(GRAPH_H, 0), "mm"),
                                             collect_guides = TRUE))[["height"]]
stopifnot(LOWER_H >= 8)
panel <- stack_parts(p_top, p_bot, heights = grid::unit(c(GRAPH_H, LOWER_H), "mm"), collect_guides = TRUE)
save_panel(panel, fig = "Fig1", panel = "d", name = "mito_content", width = PAGE_W, height = PAGE_H)

# --- numbers for the legend (printed, never drawn) --------------------------------
cat("\nFig. 1d -- for the legend\n")
cat("  TOP: n = 6 mice per group; share of the nuclear (non-mtDNA) transcriptome;\n")
cat("  P: genotype main effect of lm(log2 share ~ age + genotype), script 32\n")
for (a in arms) {
  cat(sprintf("    %-26s %4d genes  MYC+ vs WT %+.1f%%  P = %.2g  (interaction P = %.2f)\n",
              sub("\n", " ", arm_name[[a]]), stats[a, "n_genes"],
              100 * (2^stats[a, "geno_beta"] - 1), stats[a, "geno_p"], stats[a, "int_p"]))
}
cat(sprintf("  WT window, nuclear MitoCarta share 12W vs 6W: %+.1f%% (Welch on log2 share, P = %.2g)\n",
            100 * (2^diff(rev(wt_window$estimate)) - 1), wt_window$p.value))
cat(sprintf("  BOTTOM: %d nuclear MitoPathways, one point each; ticks are medians\n", nrow(r143)))
cat(sprintf("    6W:  median %+.3f, %.1f%% above zero\n", stats::median(r143$c_m6), 100 * mean(r143$c_m6 > 0)))
cat(sprintf("    12W: median %+.3f, %.1f%% above zero\n", stats::median(r143$c_m12), 100 * mean(r143$c_m12 > 0)))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(.nm_drawable(panel))
  data.frame(module = names(fac_max), lowest = lo_animal, highest = fac_max,
             axis_from = lims[, 1], axis_to = lims[, 2], need_mm = need)
  stats[arms, c("n_genes", "geno_beta", "geno_p", "int_p")]
  r143[order(-r143$c_m6), c("pathway", "c_m6", "c_m12")][1:5, ]
}
