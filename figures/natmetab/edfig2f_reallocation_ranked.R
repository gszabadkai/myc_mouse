# =============================================================================
# edfig2f_reallocation_ranked.R -- Extended Data Fig. 2f: MYC re-weights the
# mitochondrial compartment -- every MitoPathway ranked by MYC's effect on its priority
# -----------------------------------------------------------------------------
# NEW 2026-10-16 (author): the exploratory fig02b (figures/fig02_reallocation_ranked.R,
# outputs/figures/fig02b_reallocation_ranked.pdf) "reproduced basically as it is",
# about 45 mm tall across the page, above ED 2g. The old ED 2f, the window scatter,
# moved to Fig. 1f (retired/edfig2f_window_two_rulers.R).
#
# FORM: the nuclear MitoPathways ranked by MYC's effect on their priority at 6W
#   (mitoPPS, MYC+ vs WT), promoted (left) to demoted (right). Each large dot is the
#   6W effect, joined by a thin line to a small dot, the 12W effect: the fade. Colour
#   is the MitoPathway tier, as in the original (mito_tier_cols, _style.R).
# SIZE (author, later on 2026-10-16): "Make the actual plot shorter - around 120 mm,
#   move the legend key to the right side, so the plot can be taller." The data area
#   is 120 mm wide and as tall as the 45-mm row allows; the key stands to its right.
# LABELS (author, the same note): "as in fig02_reallocation_ranked.pdf -- the tier
#   labels larger, circled points, any significant point marked with an asterisk
#   (Nature doesn't want asterisk, but in this case it is an exception). Point
#   legends only for the tier labels, not individual pathways. The tier labels
#   aligned on the bottom or the top." So the seven top-level MitoPathways, the tiers
#   themselves, are ringed and named in 7-pt type (black: the text rule, where the
#   original coloured them; not bold, author, the evening of 2026-10-16); promoted ones along the foot of the empty lower left,
#   demoted ones along the head of the empty upper right, which keeps the axis to the
#   data's height (ranked_tier_marks(), _style.R). Each 6W dot at BH-adjusted
#   P < 0.05 carries an asterisk, on the side away from its 12W dot. The 33 pathway
#   names of the morning's version are gone; the legend lists them.
# THE EVENING'S NOTE (author): labels not bold; no x-axis title ("that goes in the
#   legend"); the key's title "MitoPathways (MitoCarta3.0, Hierarchy top level)", on
#   two lines so that the panel stays within 180 mm.
# STATISTICS: script 08's Welch t-test on per-animal mitoPPS scores, BH across the
#   MitoPathways within each contrast (results/mitopps_scores.rds, $mitopps_pairwise),
#   asserted equal to script 40's copy, the one Fig. 1e reads.
# 143 PATHWAYS: the nuclear MitoPathways of Fig. 1e and ED 2e. The mtDNA-encoded
#   module is left out (the original drew 144).
#
# Reads:  results/mitopps_scores.rds   (script 08: $mitopps_pairwise, $pathway_tier1_map,
#                                       $pathway_levels)
#         results/background_vs_myc.rds (script 40: $ruler, the 143 and the check)
# Output: outputs/natmetab/EDFig2/EDFig2f_reallocation_ranked.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

DATA_W  <- 120    # mm, the data area's width (author: "around 120 mm")
PANEL_H <- 45     # mm, the panel's row on the page

mp    <- readRDS(here::here("results", "mitopps_scores.rds"))
pw    <- as.data.frame(mp$mitopps_pairwise)
tier  <- mp$pathway_tier1_map
ruler <- as.data.frame(readRDS(here::here("results", "background_vs_myc.rds"))$ruler)
r143  <- ruler[!ruler$is_mtdna, ]
stopifnot(nrow(r143) == 143L)

e6  <- pw[pw$contrast == "Myc_effect_6W", ]
e12 <- pw[pw$contrast == "Myc_effect_12W", ]
d <- data.frame(pathway = r143$pathway,
                eff6  = e6$diff[match(r143$pathway, e6$pathway)],
                padj6 = e6$padj[match(r143$pathway, e6$pathway)],
                eff12 = e12$diff[match(r143$pathway, e12$pathway)])
d$tier <- factor(unname(mito_tier_labels[tier[d$pathway]]), levels = mito_tier_labels)
# the values of record: script 40's copy, which Fig. 1e draws
stopifnot(!anyNA(d), max(abs(d$eff6 - r143$p_m6)) < 1e-12,
          max(abs(d$padj6 - r143$p_m6_padj)) < 1e-12, max(abs(d$eff12 - r143$p_m12)) < 1e-12)

d <- d[order(-d$eff6), ]
d$rank <- seq_len(nrow(d))
d$sig  <- d$padj6 < 0.05
stopifnot(sum(d$sig & d$eff6 > 0) == 15L, sum(d$sig & d$eff6 < 0) == 18L)

pts <- rbind(data.frame(rank = d$rank, y = d$eff12, tier = d$tier, age = "12W"),
             data.frame(rank = d$rank, y = d$eff6,  tier = d$tier, age = "6W"))
pts$age <- factor(pts$age, levels = c("6W", "12W"))
yl  <- range(pts$y) + c(-1, 1) * 0.04 * diff(range(pts$y))
DOT <- c(`6W` = 1.1, `12W` = 0.5)

# the seven top-level MitoPathways, each its own tier: ringed and named
l1  <- mp$pathway_levels$Pathway[mp$pathway_levels$Level == "Pathway_Level1"]
anc <- d[d$pathway %in% l1, ]
stopifnot(nrow(anc) == 7L, identical(as.character(anc$tier), unname(mito_tier_labels[anc$pathway])))
anchors <- data.frame(rank = anc$rank, y = anc$eff6, tier = anc$tier,
                      label = unname(mito_tier_text[as.character(anc$tier)]))
# an asterisk on each 6W dot at adjusted P < 0.05, on the side away from its 12W dot
stars <- data.frame(rank = d$rank[d$sig], y = d$eff6[d$sig], up = d$eff12[d$sig] < d$eff6[d$sig],
                    r_mm = point_r_mm(DOT[["6W"]]))

base <- ggplot(d, aes(rank, eff6)) +
  geom_hline(yintercept = 0, linewidth = NM_LINE, colour = "grey55") +
  geom_segment(aes(xend = rank, yend = eff12), colour = "grey80", linewidth = 0.2) +
  geom_point(data = pts, aes(rank, y, colour = tier, size = age), stroke = 0) +
  scale_colour_manual(values = mito_tier_cols, name = "MitoPathways\n(MitoCarta3.0, Hierarchy top level)", drop = FALSE) +
  scale_fill_manual(values = mito_tier_cols, guide = "none") +
  scale_size_manual(values = DOT, name = NULL,
                    labels = c(`6W` = "MYC+ vs WT, 6W", `12W` = "MYC+ vs WT, 12W")) +
  scale_x_continuous(expand = expansion(add = 1.5)) +
  scale_y_nm(step = 0.2, top = yl[2], bottom = yl[1], minor = TRUE, labels = lab_signed) +
  labs(x = NULL, y = "MYC+ vs WT\n(mitoPPS)") +    # the ranking is the legend's
  guides(colour = guide_legend(order = 1, ncol = 1, override.aes = list(size = 1.5)),
         size   = guide_legend(order = 2, ncol = 1, override.aes = list(colour = "grey45"))) +
  theme_nm(legend = "right") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), legend.key.size = unit(2.6, "mm"),
        legend.spacing.y = unit(2, "mm"), legend.justification = c(0, 0.5))

# The data area is DATA_W wide and as tall as PANEL_H allows; the marks are worked
# in printed mm on it.
rest   <- natural_size(fixed_size(base, DATA_W, 30))[["height"]] - 30
DATA_H <- PANEL_H - rest
bp <- ggplot2::ggplot_build(base)$layout$panel_params[[1]]
p  <- base + ranked_tier_marks(anchors, stars, pts[, c("rank", "y")], bp$x.range, bp$y.range,
                               DATA_W, DATA_H)
panel   <- fixed_size(p, DATA_W, DATA_H)
PANEL_W <- ceiling(natural_size(panel)[["width"]])
stopifnot(PANEL_W <= NM_ED_MAX_W)

save_panel(panel, fig = "EDFig2", panel = "f", name = "reallocation_ranked",
           width = PANEL_W, height = PANEL_H)

# --- numbers for the legend (printed, never drawn) --------------------------------
cat(sprintf("\nED Fig. 2f -- for the legend: %d nuclear MitoPathways ranked by MYC's effect on priority at 6W\n", nrow(d)))
cat("  (mitoPPS, MYC+ vs WT); large dot 6W, small dot 12W; colour = MitoPathway tier\n")
cat(sprintf("  data area %d x %.1f mm; panel %d x %d mm\n", DATA_W, DATA_H, PANEL_W, PANEL_H))
cat("  ringed and named: the seven top-level MitoPathways (the tiers themselves):\n")
for (i in order(anc$rank)) cat(sprintf("    %-44s rank %3d  6W %+.3f (adj P %.3f)  12W %+.3f\n",
  anc$pathway[i], anc$rank[i], anc$eff6[i], anc$padj6[i], anc$eff12[i]))
cat(sprintf("  asterisk: BH-adjusted P < 0.05 at 6W (script 08, Welch): %d (%d promoted, %d demoted)\n",
            sum(d$sig), sum(d$sig & d$eff6 > 0), sum(d$sig & d$eff6 < 0)))
o <- d[d$sig, ]
for (i in seq_len(nrow(o))) cat(sprintf("    %-44s %-30s 6W %+.3f (adj P %.3f)  12W %+.3f\n",
  o$pathway[i], as.character(o$tier[i]), o$eff6[i], o$padj6[i], o$eff12[i]))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p)
  table(d$tier, d$sig)
}
