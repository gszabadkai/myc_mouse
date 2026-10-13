# =============================================================================
# PROPOSAL (alternative to p1g) -- the BCL-2 balance in three steps, one row per
# ratio: the 6W tilt, where the common shrink puts it at 12W, and where it is
# -----------------------------------------------------------------------------
# The author's argument (2026-10-13): MYC builds a pro-death tilt at 6W; MYC's
# overall fade shrinks that tilt (the dose part); and on top of the general shrink,
# PUMA:BCL-XL loses all of it. The same numbers as p1g's plane, in Fig. 1e's row
# grammar, so each of the three steps has its own mark:
#   large dot (#D55E00)  MYC+ vs WT at 6W: the tilt (all nine to the right of zero)
#   black tick           0.43 x the 6W tilt: the line fitted through the origin to the
#                        eight ratios other than Bbc3:Bcl2l1 (the common shrink)
#   grey bar             that line's 95% prediction interval at this 6W tilt (as the
#                        band in p1g)
#   small dot (#E69F00)  MYC+ vs WT at 12W: where the ratio is
# Bbc3:Bcl2l1 (bold) was not used to fit the line, so its tick is a prediction; for
# the other eight the tick is their fitted value. Its 12W dot is the largest
# departure of the nine, and it lies inside the interval, as in p1g. No verdict.
#
# Reads: results/two_timeline_verification.rds (script 54)
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

tv <- readRDS(here::here("results", "two_timeline_verification.rds"))
pr <- as.data.frame(tv$ratios)
ln <- as.data.frame(tv$ratio_line)
rs <- as.data.frame(tv$ratio_resid)
TARGET <- tv$params$target_ratio
L8 <- ln[ln$comparator == "the other eight ratios", ]
stopifnot(nrow(pr) == 9L, identical(TARGET, "Bbc3:Bcl2l1"), nrow(L8) == 1L,
          identical(rs$pair, pr$pair))

# --- the line and its interval ARE script 54's (as in p1g) -----------------------
m8    <- stats::lm(d12 ~ 0 + d6, data = pr[pr$pair != TARGET, ])
SLOPE <- unname(stats::coef(m8)[1])
TQ    <- stats::qt(0.975, L8$df_resid)
hw    <- function(x) TQ * sqrt(L8$sigma^2 + x^2 * L8$slope_se^2)
tg    <- pr[pr$pair == TARGET, ]
stopifnot(abs(SLOPE - L8$slope) < 1e-12,
          abs(SLOPE * tg$d6 - L8$target_pred) < 1e-12,
          abs(SLOPE * tg$d6 - hw(tg$d6) - L8$pi_lwr) < 1e-10,
          abs(SLOPE * tg$d6 + hw(tg$d6) - L8$pi_upr) < 1e-10)

# --- the three steps, asserted -------------------------------------------------------
pr$fit <- SLOPE * pr$d6
pr$res <- pr$d12 - pr$fit
stopifnot(all(pr$d6 > 0),                                       # 1. the tilt: all nine
          SLOPE < 0.5,                                          # 2. the common shrink
          pr$pair[which.max(abs(pr$res))] == TARGET,            # 3. the added effect ...
          tg$d12 < tg$fit, abs(tg$d12) < 0.1,                   #    ... to about zero
          tg$d12 > tg$fit - hw(tg$d6))                          #    inside the interval

# --- the panel -------------------------------------------------------------------------
pr$row <- factor(pr$pair, levels = pr$pair[order(pr$d6)])        # largest tilt on top
pr$y   <- as.numeric(pr$row)
pr$lwr <- pr$fit - hw(pr$d6); pr$upr <- pr$fit + hw(pr$d6)
ylab <- stats::setNames(
  ifelse(levels(pr$row) == TARGET, sprintf('bolditalic("%s")', levels(pr$row)),
         sprintf('italic("%s")', levels(pr$row))), levels(pr$row))
long <- rbind(data.frame(row = pr$row, x = pr$d6,  contrast = "myc_6W"),
              data.frame(row = pr$row, x = pr$d12, contrast = "myc_12W"))
long$contrast <- factor(long$contrast, levels = names(contrast_cols))
H <- 0.32                                                        # half the bar height, rows
lo <- min(pr$lwr, pr$d12); hi <- max(pr$upr, pr$d6)
pad <- 0.04 * (hi - lo)

p <- ggplot() +
  geom_blank(data = pr, aes(x = d6, y = row)) +
  geom_rect(data = pr, aes(xmin = lwr, xmax = upr, ymin = y - H, ymax = y + H, fill = "pi")) +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
  geom_segment(data = pr, aes(x = d6, xend = d12, y = row, yend = row),
               linewidth = 0.4, colour = "grey60") +
  geom_segment(data = pr, aes(x = fit, xend = fit, y = y - H, yend = y + H, linetype = "fit"),
               linewidth = 0.55, colour = "black", key_glyph = "vpath") +
  geom_point(data = long, aes(x, row, colour = contrast, size = contrast)) +
  scale_y_discrete(labels = parse(text = ylab)) +
  scale_colour_manual(values = contrast_cols, labels = contrast_labels, name = NULL) +
  scale_size_manual(values = c(myc_6W = 1.6, myc_12W = 1.1), guide = "none") +
  scale_linetype_manual(values = c(fit = "solid"), name = NULL,
                        labels = c(fit = sprintf("Common shrink (%.2f × 6W)", SLOPE))) +
  scale_fill_manual(values = c(pi = "grey90"), name = NULL,
                    labels = c(pi = "95% prediction interval")) +
  scale_x_nm(step = 0.2, top = hi + pad, bottom = lo - pad, minor = TRUE, labels = lab_signed) +
  labs(x = "MYC+ vs WT (log2 ratio)", y = NULL) +
  guides(colour   = guide_legend(order = 1, override.aes = list(size = c(1.6, 1.1))),
         linetype = guide_legend(order = 2),
         fill     = guide_legend(order = 3)) +
  theme_nm(legend = "bottom") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.box = "vertical", legend.spacing.y = unit(0, "mm"),
        legend.key.size = unit(2.4, "mm"))

save_panel(p, fig = "Fig1", panel = "g", name = "alt_ratio_rows", width = 62, height = 58,
           dir = here::here("outputs", "natmetab", "proposals"))

cat(sprintf("\n1g rows: common shrink %.3f (SE %.3f) over the other eight\n", SLOPE, L8$slope_se))
o <- pr[order(-pr$d6), ]
for (i in seq_len(nrow(o))) cat(sprintf("  %-15s 6W %+.3f  expected 12W %+.3f [%+.3f, %+.3f]  12W %+.3f  departure %+.3f\n",
  o$pair[i], o$d6[i], o$fit[i], o$lwr[i], o$upr[i], o$d12[i], o$res[i]))
