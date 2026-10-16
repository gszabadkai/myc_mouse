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
# LABELS (author): the pathways at BH-adjusted P < 0.05 are named, in black. No WT
#   change passes (0 of 143), so the names are the ten MYC+ changes that pass, each on
#   its open dot. The original drew no significance on purpose: the window coincides
#   with an extraction batch (6W and 12W were extracted separately), so a single
#   genotype's change is descriptive, and the legend must say so. Names keep
#   MitoCarta's wording; a hyphen is the soft hyphen.
# STATISTICS: script 08's Welch t-test on per-animal mitoPPS scores, BH across the
#   MitoPathways within each contrast (results/mitopps_scores.rds, $mitopps_pairwise),
#   asserted equal to script 40's copy (Fig. 1f's scatter).
# 143 PATHWAYS: the nuclear MitoPathways, as ED 2f; the mtDNA-encoded module is left
#   out (the original drew 144).
#
# Reads:  results/mitopps_scores.rds   (script 08: $mitopps_pairwise, $pathway_tier1_map)
#         results/background_vs_myc.rds (script 40: $ruler, the 143 and the check)
# Output: outputs/natmetab/EDFig2/EDFig2g_reallocation_window_ranked.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

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
d$name <- gsub("-", "­", d$pathway)                     # soft hyphen
stopifnot(!any(d$wt_padj < 0.05), sum(d$sig) == 10L)     # no WT change passes; ten MYC+ ones do

pts <- rbind(data.frame(rank = d$rank, y = d$wt,  tier = d$tier, geno = "neg"),
             data.frame(rank = d$rank, y = d$myc, tier = d$tier, geno = "pos"))
pts$geno <- factor(pts$geno, levels = c("neg", "pos"))
yl <- range(pts$y) + c(-1, 1) * 0.04 * diff(range(pts$y))

# The label layer sees every point, WT and MYC+, unnamed ones with an empty label, so
# the names keep off all points; names of gains start above their dot, of losses below.
rd <- rbind(data.frame(rank = d$rank, y = d$myc, label = ifelse(d$sig, paste0("  ", d$name, "  "), ""),
                       ny = ifelse(d$sig, sign(d$myc) * 0.1, 0)),
            data.frame(rank = d$rank, y = d$wt, label = "", ny = 0))

p <- ggplot(d, aes(rank, wt)) +
  geom_hline(yintercept = 0, linewidth = NM_LINE, colour = "grey55") +
  geom_segment(aes(xend = rank, yend = myc), colour = "grey80", linewidth = 0.25) +
  geom_point(data = pts, aes(rank, y, colour = tier, shape = geno, size = geno), stroke = 0.4) +
  ggrepel::geom_text_repel(data = rd, aes(rank, y, label = label), inherit.aes = FALSE,
                           nudge_y = rd$ny, size = pt(NM_TXT_MIN), family = NM_FONT, colour = "black",
                           segment.colour = "grey55", segment.size = 0.2, min.segment.length = 0,
                           box.padding = 0.08, point.padding = 0.05, force = 3, force_pull = 0.2,
                           max.overlaps = Inf, max.iter = 100000, max.time = 10, seed = 1) +
  scale_colour_manual(values = mito_tier_cols, name = "MitoPathway tier", drop = FALSE) +
  scale_shape_manual(values = c(neg = 16, pos = 1), labels = dev_labels, name = NULL) +
  scale_size_manual(values = c(neg = 1.1, pos = 0.9), guide = "none") +
  scale_x_continuous(expand = expansion(add = 1.5)) +
  scale_y_nm(step = 0.2, top = yl[2], bottom = yl[1], minor = TRUE, labels = lab_signed) +
  labs(x = "MitoPathways, ranked by the WT change (promoted to demoted)",
       y = "12W vs 6W\n(mitoPPS)") +
  guides(colour = guide_legend(order = 1, nrow = 2, override.aes = list(size = 1.5, shape = 16)),
         shape  = guide_legend(order = 2, nrow = 2, override.aes = list(colour = "grey35", size = 1.2))) +
  theme_nm(legend = "bottom") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), legend.key.size = unit(2.6, "mm"),
        legend.spacing.y = unit(1, "mm"))

save_panel(p, fig = "EDFig2", panel = "g", name = "reallocation_window_ranked", width = 180, height = 45)

# --- numbers for the legend (printed, never drawn) --------------------------------
cat(sprintf("\nED Fig. 2g -- for the legend: %d nuclear MitoPathways ranked by the WT change in priority, 12W vs 6W\n", nrow(d)))
cat("  (mitoPPS); filled dot WT, open dot MYC+, joined; colour = MitoPathway tier\n")
cat(sprintf("  named: BH-adjusted P < 0.05 (script 08, Welch): WT 0 of %d; MYC+ %d, on their open dots\n", nrow(d), sum(d$sig)))
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
