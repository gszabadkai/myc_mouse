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
# LABELS (author): the pathways at BH-adjusted P < 0.05 at 6W are named, in black,
#   and outlined: they replace the original's asterisks. The promoted names sit in
#   the empty space above the curve, the demoted ones below it. The original's seven
#   ringed top-level categories are not drawn: with every significant pathway named
#   they would add seven more labels, and the three of them that pass are named.
#   Names keep MitoCarta's wording; a hyphen is the soft hyphen (base pdf() prints "-"
#   as a minus).
# STATISTICS: script 08's Welch t-test on per-animal mitoPPS scores, BH across the
#   MitoPathways within each contrast (results/mitopps_scores.rds, $mitopps_pairwise),
#   asserted equal to script 40's copy, the one Fig. 1e reads.
# 143 PATHWAYS: the nuclear MitoPathways of Fig. 1e and ED 2e. The mtDNA-encoded
#   module is left out (the original drew 144).
#
# Reads:  results/mitopps_scores.rds   (script 08: $mitopps_pairwise, $pathway_tier1_map)
#         results/background_vs_myc.rds (script 40: $ruler, the 143 and the check)
# Output: outputs/natmetab/EDFig2/EDFig2f_reallocation_ranked.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

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
d$name <- gsub("-", "­", d$pathway)                     # soft hyphen
up <- d[d$sig & d$eff6 > 0, ]; dn <- d[d$sig & d$eff6 < 0, ]
stopifnot(nrow(up) == 15L, nrow(dn) == 18L)

pts <- rbind(data.frame(rank = d$rank, y = d$eff12, tier = d$tier, age = "12W"),
             data.frame(rank = d$rank, y = d$eff6,  tier = d$tier, age = "6W"))
pts$age <- factor(pts$age, levels = c("6W", "12W"))
yl <- range(pts$y) + c(-1, 1) * 0.04 * diff(range(pts$y))

# Each label layer sees every point, 6W and 12W, unnamed ones with an empty label, so
# the names keep off all points; each set is held to the empty side of the curve.
# Two spaces either side of a name keep names that share a row visibly apart.
rep_data <- function(keep) rbind(
  data.frame(rank = d$rank, y = d$eff6,  label = ifelse(d$pathway %in% keep$pathway,
                                                        paste0("  ", d$name, "  "), "")),
  data.frame(rank = d$rank, y = d$eff12, label = ""))
# Names start spread over their empty triangle: those of promoted pathways right of
# and above their points, those of demoted ones left of and below theirs.
repel <- function(keep, xlim, ylim, hjust, nx, ny) {
  rd <- rep_data(keep)
  on <- rd$label != ""
  ggrepel::geom_text_repel(
    data = rd, aes(rank, y, label = label), inherit.aes = FALSE,
    nudge_x = ifelse(on, nx, 0), nudge_y = ifelse(on, ny, 0),
    xlim = xlim, ylim = ylim, hjust = hjust, size = pt(NM_TXT_MIN), family = NM_FONT,
    colour = "black", segment.colour = "grey55", segment.size = 0.2, min.segment.length = 0,
    box.padding = 0.08, point.padding = 0.05, force = 3, force_pull = 0.2,
    max.overlaps = Inf, max.iter = 100000, max.time = 10, seed = 1)
}

p <- ggplot(d, aes(rank, eff6)) +
  geom_hline(yintercept = 0, linewidth = NM_LINE, colour = "grey55") +
  geom_segment(aes(xend = rank, yend = eff12), colour = "grey80", linewidth = 0.2) +
  geom_point(data = pts, aes(rank, y, colour = tier, size = age), stroke = 0) +
  geom_point(data = d[d$sig, ], aes(fill = tier), shape = 21, size = 1.15, colour = "black",
             stroke = 0.3, show.legend = FALSE) +
  repel(up, xlim = c(max(up$rank) + 3, NA), ylim = c(0.05, NA), hjust = 0, nx = 45, ny = 0.06) +
  repel(dn, xlim = c(NA, min(dn$rank) - 3), ylim = c(NA, -0.05), hjust = 1, nx = -45, ny = -0.06) +
  scale_colour_manual(values = mito_tier_cols, name = "MitoPathway tier", drop = FALSE) +
  scale_fill_manual(values = mito_tier_cols, guide = "none") +
  scale_size_manual(values = c(`6W` = 1.1, `12W` = 0.5), name = NULL,
                    labels = c(`6W` = "MYC+ vs WT, 6W", `12W` = "MYC+ vs WT, 12W")) +
  scale_x_continuous(expand = expansion(add = 1.5)) +
  scale_y_nm(step = 0.2, top = yl[2], bottom = yl[1], minor = TRUE, labels = lab_signed) +
  labs(x = "MitoPathways, ranked by MYC's effect at 6W (promoted to demoted)",
       y = "MYC+ vs WT\n(mitoPPS)") +
  guides(colour = guide_legend(order = 1, nrow = 2, override.aes = list(size = 1.5)),
         size   = guide_legend(order = 2, nrow = 2, override.aes = list(colour = "grey45"))) +
  theme_nm(legend = "bottom") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), legend.key.size = unit(2.6, "mm"),
        legend.spacing.y = unit(1, "mm"))

save_panel(p, fig = "EDFig2", panel = "f", name = "reallocation_ranked", width = 180, height = 45)

# --- numbers for the legend (printed, never drawn) --------------------------------
cat(sprintf("\nED Fig. 2f -- for the legend: %d nuclear MitoPathways ranked by MYC's effect on priority at 6W\n", nrow(d)))
cat("  (mitoPPS, MYC+ vs WT); large dot 6W, small dot 12W; colour = MitoPathway tier\n")
cat(sprintf("  named and outlined: BH-adjusted P < 0.05 at 6W (script 08, Welch): %d (%d promoted, %d demoted)\n",
            sum(d$sig), nrow(up), nrow(dn)))
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
