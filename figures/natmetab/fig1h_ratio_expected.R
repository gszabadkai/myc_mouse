# =============================================================================
# fig1h_ratio_expected.R -- Fig. 1h: each BCL-2 ratio's 12-week MYC effect,
# expected from the common shrink against observed
# -----------------------------------------------------------------------------
# PROMOTED from proposals/ on 2026-10-14 (author: "promote them to final panels").
# The author's plan (2026-10-14): "Notably, the steepest drop was in Bbc3,
# decreasing the PUMA/BCL-XL ratio" -- the simplified form of the three-step rows:
# the two points of each row are aligned with the common shrink, so the gap
# between them is only the added effect.
#   large dot (#D55E00)  the 6W MYC effect x 0.43: where the ratio would be at 12W
#                        if it shrank with the other eight (the common shrink)
#   small dot (#E69F00)  the 12W MYC effect: where it is
# Rows run by the gap, largest first. The common shrink is the line through the
# origin fitted to the eight ratios other than Bbc3:Bcl2l1 (script 54, ruling 4);
# Bbc3:Bcl2l1 (bold) was not used to fit it. Its 95% prediction interval is in the
# legend, not drawn: Bbc3:Bcl2l1 lies inside it, so no row is "significant".
#
# LABELS (author, 2026-10-14): the x axis reads "MYC+ vs WT (normalised to overall
# MYC activity)" and the key reads 6W and 12W, without the factor -- the legend
# carries it. The factor is the RATIO SET's common shrink (0.43), not the
# gene-level fade of script 55 (0.49); the legend must name which one it is.
#
# Reads: results/two_timeline_verification.rds (script 54: $ratios, $ratio_line)
# Output: outputs/natmetab/Fig1/Fig1h_ratio_expected.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

tv <- readRDS(here::here("results", "two_timeline_verification.rds"))
pr <- as.data.frame(tv$ratios)
ln <- as.data.frame(tv$ratio_line)
TARGET <- tv$params$target_ratio
L8 <- ln[ln$comparator == "the other eight ratios", ]
stopifnot(nrow(pr) == 9L, identical(TARGET, "Bbc3:Bcl2l1"), nrow(L8) == 1L)

SLOPE <- unname(stats::coef(stats::lm(d12 ~ 0 + d6, data = pr[pr$pair != TARGET, ]))[1])
stopifnot(abs(SLOPE - L8$slope) < 1e-12)

pr$expected <- SLOPE * pr$d6
pr$gap      <- pr$expected - pr$d12                 # the added loss beyond the shrink
ord <- pr$pair[order(pr$gap)]                       # largest gap on top
pr$row <- factor(pr$pair, levels = ord)

# --- the sentence, asserted -------------------------------------------------------------
tg <- pr[pr$pair == TARGET, ]
stopifnot(pr$pair[which.max(pr$gap)] == TARGET,           # the largest added loss
          abs(tg$d12) < 0.1,                              # PUMA:BCL-XL ends at about zero
          abs(L8$target_pred - tg$expected) < 1e-12,
          isTRUE(L8$target_inside_pi))                    # inside the band: legend says so

long <- rbind(data.frame(row = pr$row, x = pr$expected, what = "expected"),
              data.frame(row = pr$row, x = pr$d12,      what = "observed"))
long$what <- factor(long$what, levels = c("expected", "observed"))
ylab <- stats::setNames(ifelse(levels(pr$row) == TARGET,
                               sprintf('bolditalic("%s")', levels(pr$row)),
                               sprintf('italic("%s")', levels(pr$row))), levels(pr$row))
lo <- min(long$x); hi <- max(long$x); pad <- 0.04 * (hi - lo)

p <- ggplot(long, aes(x, row)) +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
  geom_segment(data = pr, aes(x = expected, xend = d12, y = row, yend = row),
               linewidth = 0.4, colour = "grey60") +
  geom_point(aes(colour = what, size = what)) +
  scale_colour_manual(values = c(expected = contrast_cols[["myc_6W"]],
                                 observed = contrast_cols[["myc_12W"]]),
                      labels = c(expected = "6W", observed = "12W"), name = NULL) +
  scale_size_manual(values = c(expected = 1.6, observed = 1.1), guide = "none") +
  scale_y_discrete(labels = parse(text = ylab)) +
  scale_x_nm(step = 0.2, top = hi + pad, bottom = lo - pad, minor = TRUE, labels = lab_signed) +
  labs(x = "MYC+ vs WT (normalised to overall MYC activity)", y = NULL) +
  guides(colour = guide_legend(override.aes = list(size = c(1.6, 1.1)), nrow = 1)) +
  theme_nm(legend = "bottom") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.key.size = unit(2.4, "mm"))

save_panel(p, fig = "Fig1", panel = "h", name = "ratio_expected", width = 100, height = 60)

cat(sprintf("\n1h -- for the legend: the 6W points are scaled by the common shrink %.3f (SE %.3f), the line\n",
            SLOPE, L8$slope_se))
cat("  through the other eight ratios, so both ages read at the 12W level of overall MYC activity\n")
cat(sprintf("  Bbc3:Bcl2l1 expected %+.3f, observed %+.3f; inside the 95%% prediction interval [%+.3f, %+.3f]\n",
            tg$expected, tg$d12, L8$pi_lwr, L8$pi_upr))
o <- pr[order(-pr$gap), ]
for (i in seq_len(nrow(o))) cat(sprintf("  %-15s 6W %+.3f  expected 12W %+.3f  observed %+.3f  gap %+.3f\n",
  o$pair[i], o$d6[i], o$expected[i], o$d12[i], o$gap[i]))
