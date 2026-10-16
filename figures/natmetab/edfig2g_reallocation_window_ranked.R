# =============================================================================
# edfig2g_reallocation_window_ranked.R -- Extended Data Fig. 2g: across the window,
# every MitoPathway ranked by the WT gland's change in priority, with the MYC+ change
# -----------------------------------------------------------------------------
# NEW 2026-10-16 (author): the exploratory figS6 (figures/figS6_reallocation_temporal_
# ranked.R, outputs/figures/figS6_reallocation_temporal_ranked.pdf) "reproduced
# basically as it is", about 45 mm tall across the page, under ED 2f.
#
# FORM: the nuclear MitoPathways ranked by the WT change in priority from 6W to 12W
#   (mitoPPS), promoted (left) to demoted (right). The filled dot is the WT change; a
#   line drops to the open dot, the MYC+ change: where MYC bends the WT trajectory.
#   Colour is the MitoPathway tier, as in the original (mito_tier_cols, _style.R).
# SIZE AND LABELS: as ED 2f (author, later on 2026-10-16; see that script's header).
#   The data area is 120 mm wide, the key to its right. The seven top-level
#   MitoPathways are ringed on their WT dots and named in 7-pt type: those the WT
#   window promotes along the foot of the empty lower left, those it demotes along
#   the head of the empty upper right (ranked_tier_marks(), _style.R). An asterisk
#   marks each change at BH-adjusted P < 0.05. No WT change passes (0 of 143), so the
#   asterisks are the ten MYC+ changes that do, each on its open dot, on the side
#   away from its WT dot. The original drew no significance on purpose: the window
#   coincides with an extraction batch (6W and 12W were extracted separately), so a
#   single genotype's change is descriptive, and the legend must say so.
#   The evening's note (as ED 2f): labels not bold, no x-axis title, the key titled
#   "MitoPathways (MitoCarta3.0, Hierarchy top level)".
# STATISTICS: script 08's Welch t-test on per-animal mitoPPS scores, BH across the
#   MitoPathways within each contrast (results/mitopps_scores.rds, $mitopps_pairwise),
#   asserted equal to script 40's copy (Fig. 1f's scatter).
# 143 PATHWAYS: the nuclear MitoPathways, as ED 2f; the mtDNA-encoded module is left
#   out (the original drew 144).
#
# Reads:  results/mitopps_scores.rds   (script 08: $mitopps_pairwise, $pathway_tier1_map,
#                                       $pathway_levels)
#         results/background_vs_myc.rds (script 40: $ruler, the 143 and the check)
# Output: outputs/natmetab/EDFig2/EDFig2g_reallocation_window_ranked.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

DATA_W  <- 120    # mm, the data area's width, as ED 2f
PANEL_H <- 45     # mm, the panel's row on the page

mp    <- readRDS(here::here("results", "mitopps_scores.rds"))
pw    <- as.data.frame(mp$mitopps_pairwise)
tier  <- mp$pathway_tier1_map
ruler <- as.data.frame(readRDS(here::here("results", "background_vs_myc.rds"))$ruler)
r143  <- ruler[!ruler$is_mtdna, ]
stopifnot(nrow(r143) == 143L)

w <- pw[pw$contrast == "Temporal_Myc-", ]
m <- pw[pw$contrast == "Temporal_Myc+", ]
d <- data.frame(pathway = r143$pathway,
                wt  = w$diff[match(r143$pathway, w$pathway)], wt_padj  = w$padj[match(r143$pathway, w$pathway)],
                myc = m$diff[match(r143$pathway, m$pathway)], myc_padj = m$padj[match(r143$pathway, m$pathway)])
d$tier <- factor(unname(mito_tier_labels[tier[d$pathway]]), levels = mito_tier_labels)
# the values of record: script 40's copy, which Fig. 1f's scatter draws
stopifnot(!anyNA(d), max(abs(d$wt - r143$p_tn)) < 1e-12, max(abs(d$myc - r143$p_tp)) < 1e-12)

d <- d[order(-d$wt), ]
d$rank <- seq_len(nrow(d))
d$sig  <- d$myc_padj < 0.05
stopifnot(!any(d$wt_padj < 0.05), sum(d$sig) == 10L)     # no WT change passes; ten MYC+ ones do

pts <- rbind(data.frame(rank = d$rank, y = d$wt,  tier = d$tier, geno = "neg"),
             data.frame(rank = d$rank, y = d$myc, tier = d$tier, geno = "pos"))
pts$geno <- factor(pts$geno, levels = c("neg", "pos"))
yl  <- range(pts$y) + c(-1, 1) * 0.04 * diff(range(pts$y))
DOT <- c(neg = 1.1, pos = 0.9)
OPEN_STROKE <- 0.4

# the seven top-level MitoPathways, each its own tier: ringed on the WT dot and named
l1  <- mp$pathway_levels$Pathway[mp$pathway_levels$Level == "Pathway_Level1"]
anc <- d[d$pathway %in% l1, ]
stopifnot(nrow(anc) == 7L, identical(as.character(anc$tier), unname(mito_tier_labels[anc$pathway])))
anchors <- data.frame(rank = anc$rank, y = anc$wt, tier = anc$tier,
                      label = unname(mito_tier_text[as.character(anc$tier)]))
# an asterisk on each MYC+ change at adjusted P < 0.05, on the side away from its WT dot
stars <- data.frame(rank = d$rank[d$sig], y = d$myc[d$sig], up = d$myc[d$sig] > d$wt[d$sig],
                    r_mm = point_r_mm(DOT[["pos"]], OPEN_STROKE))

base <- ggplot(d, aes(rank, wt)) +
  geom_hline(yintercept = 0, linewidth = NM_LINE, colour = "grey55") +
  geom_segment(aes(xend = rank, yend = myc), colour = "grey80", linewidth = 0.25) +
  geom_point(data = pts, aes(rank, y, colour = tier, shape = geno, size = geno), stroke = OPEN_STROKE) +
  scale_colour_manual(values = mito_tier_cols, name = "MitoPathways\n(MitoCarta3.0, Hierarchy top level)", drop = FALSE) +
  scale_fill_manual(values = mito_tier_cols, guide = "none") +
  scale_shape_manual(values = c(neg = 16, pos = 1), labels = dev_labels, name = NULL) +
  scale_size_manual(values = DOT, guide = "none") +
  scale_x_continuous(expand = expansion(add = 1.5)) +
  scale_y_nm(step = 0.2, top = yl[2], bottom = yl[1], minor = TRUE, labels = lab_signed) +
  labs(x = NULL, y = "12W vs 6W\n(mitoPPS)") +     # the ranking is the legend's
  guides(colour = guide_legend(order = 1, ncol = 1, override.aes = list(size = 1.5, shape = 16)),
         shape  = guide_legend(order = 2, ncol = 1, override.aes = list(colour = "grey35", size = 1.2))) +
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

save_panel(panel, fig = "EDFig2", panel = "g", name = "reallocation_window_ranked",
           width = PANEL_W, height = PANEL_H)

# --- numbers for the legend (printed, never drawn) --------------------------------
cat(sprintf("\nED Fig. 2g -- for the legend: %d nuclear MitoPathways ranked by the WT change in priority, 12W vs 6W\n", nrow(d)))
cat("  (mitoPPS); filled dot WT, open dot MYC+, joined; colour = MitoPathway tier\n")
cat(sprintf("  data area %d x %.1f mm; panel %d x %d mm\n", DATA_W, DATA_H, PANEL_W, PANEL_H))
cat("  ringed (on the WT dot) and named: the seven top-level MitoPathways (the tiers themselves):\n")
for (i in order(anc$rank)) cat(sprintf("    %-44s rank %3d  WT %+.3f (adj P %.2f)  MYC+ %+.3f (adj P %.2f)\n",
  anc$pathway[i], anc$rank[i], anc$wt[i], anc$wt_padj[i], anc$myc[i], anc$myc_padj[i]))
cat(sprintf("  asterisk: BH-adjusted P < 0.05 (script 08, Welch): WT 0 of %d; MYC+ %d, on their open dots\n",
            nrow(d), sum(d$sig)))
cat("  BATCH: 6W and 12W were extracted as separate batches, so each window coincides with a batch;\n")
cat("  a single genotype's change is descriptive\n")
o <- d[d$sig, ][order(-d$myc[d$sig]), ]
for (i in seq_len(nrow(o))) cat(sprintf("    %-44s %-30s MYC+ %+.3f (adj P %.3f)  WT %+.3f (adj P %.2f)\n",
  o$pathway[i], as.character(o$tier[i]), o$myc[i], o$myc_padj[i], o$wt[i], o$wt_padj[i]))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p)
  stats::cor(d$wt, d$myc)
}
