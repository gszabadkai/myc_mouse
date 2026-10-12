# =============================================================================
# edfig2a_design.R -- Extended Data Fig. 2a: the design, four groups and the
# contrasts drawn between them
# -----------------------------------------------------------------------------
# The text: "To understand the permissive mechanism, we profiled MECs purified
# from 6- and 12-week glands ... (ED Fig. 2a-c)." From the old Fig. S1A
# (figures/panels/figS1_design_contrasts.R), redrawn in this layer.
#
# A 2 x 2 schematic: genotype in rows, age in columns, six mice per group. The
# arrows are the contrasts every effect panel draws, in the colours those panels
# use, so this panel doubles as the key to them:
#   down   "MYC+ vs WT, 6W" / "MYC+ vs WT, 12W"   contrast_cols (MYC+ hue of the age)
#   across "12W vs 6W, WT" / "12W vs 6W, MYC+"    geno_cols (the genotype's colour)
#   middle the interaction: does MYC's effect differ between the ages?
# Boxes take the group tints with black text (no text in colour, FIGURE_RULES.md).
#
# NOT ON THE PANEL: each age was extracted as one batch, so batch = age. The
# genotype contrasts and the interaction are within batch; the window contrasts
# coincide with it. That goes in the legend (printed below), as it did for S1A.
#
# Canvas in millimetres: the panel is the page (no margins, no axes), so every
# coordinate below is a position on the printed panel.
#
# Reads:  nothing (the group sizes are asserted against results/mitopps_scores.rds)
# Output: outputs/natmetab/EDFig2/EDFig2a_design.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

W <- 46; H <- 40                                    # mm, the printed panel

# six mice per group, from the sample sheet the analyses ran on
sm <- readRDS(here::here("results", "mitopps_scores.rds"))$mitopps_scores
n_per <- table(sm$group)
stopifnot(all(n_per[order_by_age] == 6L))

# 20-mm gaps between the boxes, each way: the longest arrow label is about 17.5 mm
cx <- c(`6W` = 13, `12W` = 41)                      # column centres
ry <- c(neg = 32, pos = 6.5)                        # row centres
bw <- 8; bh <- 5.5                                  # box size
boxes <- data.frame(
  group = c("6W_neg", "12W_neg", "6W_pos", "12W_pos"),
  x = cx[c("6W", "12W", "6W", "12W")], y = ry[c("neg", "neg", "pos", "pos")])
boxes$group <- factor(boxes$group, levels = order_by_age)

g <- 0.8                                            # gap between arrow tip and box
arrows <- data.frame(
  key  = c("myc_6W", "myc_12W", "neg", "pos"),
  x    = c(cx[["6W"]], cx[["12W"]], cx[["6W"]] + bw / 2 + g, cx[["6W"]] + bw / 2 + g),
  xend = c(cx[["6W"]], cx[["12W"]], cx[["12W"]] - bw / 2 - g, cx[["12W"]] - bw / 2 - g),
  y    = c(ry[["neg"]] - bh / 2 - g, ry[["neg"]] - bh / 2 - g, ry[["neg"]], ry[["pos"]]),
  yend = c(ry[["pos"]] + bh / 2 + g, ry[["pos"]] + bh / 2 + g, ry[["neg"]], ry[["pos"]]),
  col  = c(contrast_cols[["myc_6W"]], contrast_cols[["myc_12W"]], geno_cols[["neg"]],
           geno_cols[["pos"]]))

mid_y <- mean(ry)
labs_df <- data.frame(
  x     = c(cx[["6W"]] - 2.2, cx[["12W"]] + 2.2, mean(cx), mean(cx), mean(cx)),
  y     = c(mid_y, mid_y, ry[["neg"]] + 2.0, ry[["pos"]] - 2.2, mid_y + 1.7),
  label = c(contrast_labels[["myc_6W"]], contrast_labels[["myc_12W"]], dev_labels[["neg"]],
            dev_labels[["pos"]], "Interaction"),
  angle = c(90, 90, 0, 0, 0))

heads <- data.frame(x = c(cx, 0, 0), y = c(H - 1.4, H - 1.4, ry[["neg"]], ry[["pos"]]),
                    label = c("6W", "12W", "WT", "MYC+"), hjust = c(0.5, 0.5, 0, 0))

arw <- grid::arrow(length = grid::unit(1.1, "mm"), angle = 25, type = "closed")
p <- ggplot() +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend, colour = I(col)),
               linewidth = 0.5, arrow = arw, linejoin = "mitre") +
  # the interaction: a dotted line with solid heads (a dotted head draws as dots)
  geom_segment(aes(x = cx[["6W"]] + 4.4, xend = cx[["12W"]] - 4.4, y = mid_y, yend = mid_y),
               linewidth = 0.35, linetype = "11", colour = "grey45") +
  geom_segment(data = data.frame(x = cx[["6W"]] + c(4.6, 0) + c(0, cx[["12W"]] - cx[["6W"]] - 4.6),
                                 xend = c(cx[["6W"]] + 3.4, cx[["12W"]] - 3.4)),
               aes(x = x, xend = xend, y = mid_y, yend = mid_y), linewidth = 0.35,
               colour = "grey45",
               arrow = grid::arrow(length = grid::unit(0.9, "mm"), angle = 25, type = "closed")) +
  geom_tile(data = boxes, aes(x, y, fill = group, colour = group), width = bw, height = bh,
            linewidth = 0.5) +
  geom_text(data = boxes, aes(x, y), label = "n = 6", size = pt(NM_TXT), family = NM_FONT) +
  geom_text(data = labs_df, aes(x, y, label = label, angle = angle), size = pt(NM_TXT),
            family = NM_FONT) +
  geom_text(data = heads, aes(x, y, label = label, hjust = hjust), size = pt(NM_TXT_TITLE),
            family = NM_FONT) +
  scale_fill_manual(values = group_fill, guide = "none") +
  scale_colour_manual(values = group_cols, guide = "none") +
  scale_x_continuous(limits = c(0, W), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, H), expand = c(0, 0)) +
  coord_fixed(clip = "off") +
  theme_void(base_size = NM_TXT_TITLE, base_family = NM_FONT) +
  theme(plot.margin = margin(0, 0, 0, 0))

save_panel(p, fig = "EDFig2", panel = "a", name = "design", width = W, height = H)

# --- numbers for the legend (printed, never drawn) --------------------------------
cat("\nED Fig. 2a -- for the legend\n")
cat("  MMTV-Myc (MYC+) and wild-type (WT) littermates at 6 and 12 weeks; n = 6 mice per group;\n")
cat("  mammary epithelial cells purified by enzymatic dissociation, bulk RNA-seq;\n")
cat("  DESeq2 ~ age * genotype. Arrows: the genotype effect at each age, the 12W vs 6W\n")
cat("  change in each genotype, and their interaction.\n")
cat("  BATCH: each age was extracted as one batch, so the 12W vs 6W changes coincide with\n")
cat("  batch; the genotype effects and the interaction are within batch.\n")

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p)
}
