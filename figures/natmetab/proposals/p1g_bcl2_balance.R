# =============================================================================
# PROPOSAL -- Fig. 1g: the BCL-2 balance at 6 and 12 weeks
# -----------------------------------------------------------------------------
# The ruling-4 Fig. 2G (figures/panels/fig2_priming_ratios.R, signed off
# 2026-09-21) redrawn in the Nat Metab layer. Nine pro:anti transcript ratios of
# the BCL-2 family: the MYC effect at 6W (x) against the MYC effect at 12W (y).
# Both axes compare genotypes within one extraction batch, so the panel is free of
# the batch = age confound. The line is fitted through the origin to the eight
# ratios that are not Bbc3:Bcl2l1; the band is its 95% prediction interval, drawn
# in full. Bbc3:Bcl2l1 is read against a line it did not help to draw.
#
# Change from the old panel: Bbc3:Bcl2l1 was the OPEN point. In this figure open
# means "not significant" (Fig. 1e), so here it is black and the eight are grey.
#
# Text it carries (draft paragraph): "At 6 weeks, MYC raised all nine pro- to
# anti-apoptotic transcript ratios ... By 12 weeks, MYC shifted the other eight by
# less than half as much, but no longer shifted the PUMA:BCL-XL ratio at all."
#
# Reads: results/two_timeline_verification.rds (script 54)
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

tv <- readRDS(here::here("results", "two_timeline_verification.rds"))
pr <- as.data.frame(tv$ratios)
ln <- as.data.frame(tv$ratio_line)
rs <- as.data.frame(tv$ratio_resid)
bd <- as.data.frame(tv$ratio_band)
TARGET <- tv$params$target_ratio
stopifnot(nrow(pr) == 9L, identical(TARGET, "Bbc3:Bcl2l1"), identical(rs$pair, pr$pair),
          setequal(pr$anti, c("Bcl2l1", "Mcl1")))

# --- what is drawn IS script 54's fit (ported from the old 2G) -----------------
L8 <- ln[ln$comparator == "the other eight ratios", ]
others <- pr[pr$pair != TARGET, ]
tgt    <- pr[pr$pair == TARGET, ]
m8     <- stats::lm(d12 ~ 0 + d6, data = others)
SLOPE  <- unname(stats::coef(m8)[1])
TQ     <- stats::qt(0.975, L8$df_resid)
hw     <- function(x) TQ * sqrt(L8$sigma^2 + x^2 * L8$slope_se^2)
stopifnot(nrow(L8) == 1L, abs(SLOPE - L8$slope) < 1e-12,
          max(abs(bd$fit - SLOPE * bd$d6)) < 1e-12,
          max(abs(bd$lwr - (bd$fit - hw(bd$d6)))) < 1e-10,
          max(abs(bd$upr - (bd$fit + hw(bd$d6)))) < 1e-10)

# --- the sentences, asserted -----------------------------------------------------
stopifnot(all(pr$d6 > 0),                                  # MYC raised all nine at 6W
          SLOPE < 0.5,                                     # the other eight: less than half
          abs(tgt$d12) < 0.1,                              # PUMA:BCL-XL: no longer shifted
          rs$abs_rank[rs$pair == TARGET] == 1,             # the largest departure of nine
          L8$target_resid < 0, isTRUE(L8$target_inside_pi))  # below the line, inside the band

# --- the panel -----------------------------------------------------------------
pr$den  <- factor(pr$anti, levels = c("Bcl2l1", "Mcl1"))
pr$role <- factor(ifelse(pr$pair == TARGET, "target", "fit"), levels = c("fit", "target"))
pr$lab  <- ifelse(pr$pair == TARGET, "Bbc3:Bcl2l1", pr$pro)   # plain strings, set italic below

XT <- 1.0
YB <- min(bd$lwr); YT <- max(bd$upr)
ann <- data.frame(x = 0.02, y = 1.06, lab = sprintf('slope~"%.2f"', SLOPE))

p <- ggplot(pr, aes(d6, d12)) +
  geom_ribbon(data = bd, aes(x = d6, ymin = lwr, ymax = upr), inherit.aes = FALSE,
              fill = "grey91") +
  geom_hline(yintercept = 0, linewidth = NM_LINE, colour = "grey70") +
  geom_line(data = bd, aes(d6, fit), linewidth = 0.45, colour = "black") +
  geom_point(aes(shape = den, colour = role), size = 1.5) +
  ggrepel::geom_text_repel(aes(label = lab), size = pt(NM_TXT), fontface = "italic",
                           family = NM_FONT, colour = "black", seed = 3,
                           min.segment.length = 0.2, segment.size = 0.2,
                           segment.colour = "grey50", box.padding = 0.18,
                           point.padding = 0.12, max.overlaps = Inf) +
  geom_text(data = ann, aes(x, y, label = lab), parse = TRUE, hjust = 0,
            size = pt(NM_TXT), family = NM_FONT) +
  scale_shape_manual(values = c(Bcl2l1 = 16, Mcl1 = 17), name = NULL,
                     labels = c(Bcl2l1 = expression("Ratio to" ~ italic("Bcl2l1")),
                                Mcl1   = expression("Ratio to" ~ italic("Mcl1")))) +
  # grey: the eight ratios the line is fitted to; black: Bbc3:Bcl2l1, read against it
  # (named on the panel; the legend says it). No key entry: it would not fit.
  scale_colour_manual(values = c(fit = "grey55", target = "black"), guide = "none") +
  scale_x_nm(step = 0.2, top = XT, minor = TRUE, labels = lab_signed) +
  scale_y_nm(step = 0.2, top = YT, bottom = YB, minor = TRUE, labels = lab_signed) +
  labs(x = "MYC+ vs WT, 6W (log2 ratio)", y = "MYC+ vs WT, 12W (log2 ratio)") +
  guides(shape = guide_legend(override.aes = list(colour = "grey55", size = 1.3))) +
  theme_nm(legend = "inside") +
  theme(legend.position.inside = c(0.03, 0.93), legend.justification.inside = c(0, 1),
        legend.key.size = unit(2.4, "mm"), legend.background = element_blank(),
        legend.key = element_blank())

PANEL_W <- 56; PANEL_H <- 56
save_panel(p, fig = "Fig1", panel = "g", name = "bcl2_balance", width = PANEL_W, height = PANEL_H,
           dir = here::here("outputs", "natmetab", "proposals"))

cat(sprintf("\n1g: slope %.3f (SE %.3f); Bbc3:Bcl2l1 %+.3f -> %+.3f; residual %.3f (%.2f SD), PI %.3f to %.3f\n",
            SLOPE, L8$slope_se, tgt$d6, tgt$d12, L8$target_resid, L8$target_resid_in_sigma,
            L8$pi_lwr, L8$pi_upr))
cat("6W P (per-animal OLS, script 42):", paste(sprintf("%s %.2g", pr$pair, pr$p6), collapse = "; "), "\n")
