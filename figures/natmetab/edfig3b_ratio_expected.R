# =============================================================================
# edfig3b_ratio_expected.R -- Extended Data Fig. 3b: each BCL-2 ratio's MYC effect
# at 6 weeks, scaled by the common shrink, against its 12-week effect
# -----------------------------------------------------------------------------
# MOVED FROM FIG. 1h (author, 2026-10-16: "Figure 1h I think is not the one which
# should be on the main figure"); the file is fig1h_ratio_expected.R, renamed with
# git mv. Fig. 1h is now the BCL-2 family's departure from MYC's fade
# (fig1h_bcl2_fade.R). Two changes (author, 2026-10-16):
#   ROWS run by the size of the 6W effect, largest on top -- not by the gap between
#     the two points, as in Fig. 1h.
#   SIGNIFICANCE of the 6W effect (MYC+ vs WT at 6W): filled = P < 0.05, open = not.
#     The test is script 42's: lm(log2 ratio ~ age * genotype) over the 24 animals
#     (log2 of normalised counts + 1), the genotype coefficient at 6W, two-sided,
#     unadjusted -- read from script 54's copy ($ratios: d6, p6), never refitted.
#     The 12W points carry no test and stay filled; the key draws the significance
#     symbols in the 6W colour.
# THE SENTENCE (author's reading, 2026-10-16, asserted below): Bbc3:Bcl2l1 is the one
# ratio that rises significantly at 6W and ends below zero at 12W.
#
# THE POINTS, as in Fig. 1h (2026-10-14): the two points of each row are aligned with
# the common shrink, so the gap between them is only the added effect.
#   large dot (#D55E00)  the 6W MYC effect x 0.43: where the ratio would be at 12W
#                        if it shrank with the other eight (the common shrink)
#   small dot (#E69F00)  the 12W MYC effect: where it is
# The common shrink is the line through the origin fitted to the eight ratios other
# than Bbc3:Bcl2l1 (script 54, ruling 4); Bbc3:Bcl2l1 (bold) was not used to fit it.
# Its 95% prediction interval is in the legend, not drawn: Bbc3:Bcl2l1 lies inside it.
# The x axis reads "MYC+ vs WT (normalised to overall MYC activity)" and the key
# reads 6W and 12W, without the factor -- the legend carries it. The factor is the
# RATIO SET's common shrink (0.43), not the gene-level fade of script 55 (0.49); the
# legend must name which one it is. Scaling does not change the 6W test: the drawn
# 6W point is 0 exactly when the 6W effect is.
#
# TYPE AND SIZE: Extended Data's 6/7 pt (author, 2026-10-16: the ED figures keep
# their type); 86 x 45 mm, under ED 3a (86 x 44), the two as tall as ED 3c (92).
#
# Reads: results/two_timeline_verification.rds (script 54: $ratios, $ratio_line)
# Output: outputs/natmetab/EDFig3/EDFig3b_ratio_expected.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

tv <- readRDS(here::here("results", "two_timeline_verification.rds"))
pr <- as.data.frame(tv$ratios)
ln <- as.data.frame(tv$ratio_line)
TARGET <- tv$params$target_ratio
L8 <- ln[ln$comparator == "the other eight ratios", ]
stopifnot(nrow(pr) == 9L, identical(TARGET, "Bbc3:Bcl2l1"), nrow(L8) == 1L, !anyNA(pr$p6))

SLOPE <- unname(stats::coef(stats::lm(d12 ~ 0 + d6, data = pr[pr$pair != TARGET, ]))[1])
stopifnot(abs(SLOPE - L8$slope) < 1e-12)

pr$expected <- SLOPE * pr$d6
pr$gap      <- pr$expected - pr$d12                 # the added loss beyond the shrink
pr$sig6     <- pr$p6 < 0.05                         # the 6W effect, script 42's test
pr$row      <- factor(pr$pair, levels = pr$pair[order(pr$d6)])    # largest 6W effect on top

# --- the sentences, asserted -------------------------------------------------------------
tg <- pr[pr$pair == TARGET, ]
stopifnot(identical(pr$pair[pr$sig6 & pr$d12 < 0], TARGET),   # the one that rises at 6W, ends below 0
          pr$pair[which.max(pr$gap)] == TARGET,           # and the largest added loss
          abs(tg$d12) < 0.1,                              # PUMA:BCL-XL ends at about zero
          abs(L8$target_pred - tg$expected) < 1e-12,
          isTRUE(L8$target_inside_pi))                    # inside the band: legend says so

long <- rbind(data.frame(row = pr$row, x = pr$expected, what = "expected",
                         mark = ifelse(pr$sig6, "sig", "ns")),
              data.frame(row = pr$row, x = pr$d12, what = "observed", mark = "none"))
long$what <- factor(long$what, levels = c("expected", "observed"))
long$mark <- factor(long$mark, levels = c("sig", "ns", "none"))
ylab <- stats::setNames(ifelse(levels(pr$row) == TARGET,
                               sprintf('bolditalic("%s")', levels(pr$row)),
                               sprintf('italic("%s")', levels(pr$row))), levels(pr$row))
lo <- min(long$x); hi <- max(long$x); pad <- 0.04 * (hi - lo)

p <- ggplot(long, aes(x, row)) +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
  geom_segment(data = pr, aes(x = expected, xend = d12, y = row, yend = row),
               linewidth = 0.4, colour = "grey60") +
  geom_point(aes(colour = what, size = what, shape = mark), stroke = 0.5) +
  scale_colour_manual(values = c(expected = contrast_cols[["myc_6W"]],
                                 observed = contrast_cols[["myc_12W"]]),
                      labels = c(expected = "6W", observed = "12W"), name = NULL) +
  scale_size_manual(values = c(expected = 1.6, observed = 1.1), guide = "none") +
  scale_shape_manual(values = c(sig = 16, ns = 1, none = 16), breaks = c("sig", "ns"), name = NULL,
                     labels = c(sig = expression(italic(P) < "0.05"), ns = "Not significant")) +
  scale_y_discrete(labels = parse(text = ylab)) +
  scale_x_nm(step = 0.2, top = hi + pad, bottom = lo - pad, minor = TRUE, labels = lab_signed) +
  labs(x = "MYC+ vs WT (normalised to overall MYC activity)", y = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = c(1.6, 1.1), shape = 16)),
         shape  = guide_legend(order = 2, override.aes = list(size = 1.6,
                                                              colour = contrast_cols[["myc_6W"]]))) +
  theme_nm(legend = "bottom") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.key.size = unit(2.4, "mm"), legend.box = "horizontal",
        legend.spacing.x = unit(4, "mm"))

save_panel(p, fig = "EDFig3", panel = "b", name = "ratio_expected", width = 86, height = 45)

cat(sprintf("\nED 3b -- for the legend: the 6W points are scaled by the common shrink %.3f (SE %.3f), the line\n",
            SLOPE, L8$slope_se))
cat("  through the other eight ratios, so both ages read at the 12W level of overall MYC activity\n")
cat("  filled 6W point: the 6W MYC effect on the ratio P < 0.05 (script 42: lm(log2 ratio ~ age * genotype),\n")
cat("  24 animals, the genotype coefficient at 6W, two-sided t, unadjusted); 12W points carry no test\n")
cat(sprintf("  Bbc3:Bcl2l1 expected %+.3f, observed %+.3f; inside the 95%% prediction interval [%+.3f, %+.3f]\n",
            tg$expected, tg$d12, L8$pi_lwr, L8$pi_upr))
o <- pr[order(-pr$d6), ]
for (i in seq_len(nrow(o))) cat(sprintf("  %-15s 6W %+.3f (P %.2g%s)  expected 12W %+.3f  observed %+.3f  gap %+.3f\n",
  o$pair[i], o$d6[i], o$p6[i], if (o$sig6[i]) ", filled" else "", o$expected[i], o$d12[i], o$gap[i]))
