# =============================================================================
# PROPOSAL -- Fig. 1i: OXPHOS priority against the PUMA:BCL-XL balance, animal
# by animal (old Fig. 2I / 2H+I (alt) part B, in the Nat Metab layer)
# -----------------------------------------------------------------------------
# Both variables z-scored over the 24 animals; one line per genotype (12 animals),
# drawn over that genotype's own range. Left: no covariates. Right: adjusted for
# the epithelial and immune RNA composites (the pre-specified model). Each half
# prints its slope difference with its parametric and its permutation P
# (within-timepoint permutation, 5,000 draws): the impasse of ruling 7 -- the
# unadjusted interaction clears its permutation null and the adjusted one does
# not. Hypothesis-generating; it asserts nothing (verification note, 1.6).
#
# Reads: results/two_timeline_verification.rds (script 54)
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

tv <- readRDS(here::here("results", "two_timeline_verification.rds"))
pd <- as.data.frame(tv$coupling_panel)
ln <- as.data.frame(tv$coupling_lines)
cf <- as.data.frame(tv$coupling_fits)
ix <- as.data.frame(tv$coupling_interaction)
pm <- as.data.frame(tv$coupling_perm)
stopifnot(nrow(pd) == 24L, nrow(ln) == 4L, all(table(pd$genotype) == 12L))

# --- the drawn lines ARE the reported fits (ported from coupling_two_fits()) ------
sl  <- function(fit, g) cf$slope[cf$axis == "ox_ppd" & cf$fit == fit & cf$genotype == g]
dr  <- function(fit, g) ln$slope[ln$fit == fit & ln$genotype == g]
i_u <- ix[ix$axis == "ox_ppd" & ix$covariates == "none", ]
i_a <- ix[ix$axis == "ox_ppd" & ix$covariates == "epi + imm", ]
pp  <- function(cv) pm$p_emp[pm$axis == "ox_ppd" & pm$covariates == cv]
FU  <- "U unadjusted, within genotype"; FP <- "P pooled, shared covariates: epi + imm"
stopifnot(abs(dr("unadjusted", "neg") - sl(FU, "neg")) < 1e-10,
          abs(dr("unadjusted", "pos") - sl(FU, "pos")) < 1e-10,
          abs(dr("adjusted", "neg") - sl(FP, "neg")) < 1e-10,
          abs(dr("adjusted", "pos") - sl(FP, "pos")) < 1e-10,
          abs((dr("unadjusted", "pos") - dr("unadjusted", "neg")) - i_u$interaction) < 1e-10,
          abs((dr("adjusted", "pos") - dr("adjusted", "neg")) - i_a$interaction) < 1e-10,
          # the impasse, as ruled: unadjusted clears its permutation null, adjusted does not
          pp("none") < 0.05, pp("epi + imm") >= 0.05)

FITS <- c(unadjusted = "Unadjusted", adjusted = "Adjusted for epithelial\nand immune content")
long <- rbind(data.frame(pd[, c("group", "genotype", "x")], y = pd$y_unadjusted, fit = "unadjusted"),
              data.frame(pd[, c("group", "genotype", "x")], y = pd$y_adjusted,   fit = "adjusted"))
long$fit   <- factor(FITS[long$fit], levels = FITS)
long$group <- factor(long$group, levels = order_by_geno)
xr  <- tapply(pd$x, pd$genotype, range)
seg <- do.call(rbind, lapply(seq_len(nrow(ln)), function(i) {
  x <- xr[[ln$genotype[i]]]
  data.frame(fit = FITS[[ln$fit[i]]], genotype = ln$genotype[i], x = x[1], xend = x[2],
             y = ln$intercept[i] + ln$slope[i] * x[1], yend = ln$intercept[i] + ln$slope[i] * x[2])
}))
seg$fit <- factor(seg$fit, levels = FITS)

# the statistics, three lines at the top of each half (a band above the data)
XL <- c(-1.5, 3); YB <- -3
ytop <- max(c(long$y, seg$y, seg$yend))
YT <- ytop + 1.9
fmt <- function(v) formatC(signif(v, 2), format = "fg", digits = 2, flag = "#")
lab <- do.call(rbind, lapply(c("unadjusted", "adjusted"), function(f) {
  i  <- if (f == "unadjusted") i_u else i_a
  pe <- pp(if (f == "unadjusted") "none" else "epi + imm")
  data.frame(fit = FITS[[f]], x = XL[1] + 0.12, y = YT - c(0.12, 0.72, 1.32),
             lab = c(sprintf('"Slope difference"~"%+.2f"', i$interaction),
                     sprintf("italic(P)=='%s'", fmt(i$p)),
                     sprintf("'Permutation'~italic(P)=='%s'", fmt(pe))))
}))
lab$fit <- factor(lab$fit, levels = FITS)

p <- ggplot(long, aes(x, y)) +
  geom_hline(yintercept = 0, linewidth = NM_LINE, colour = "grey85") +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = "grey85") +
  geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend, colour = genotype),
               inherit.aes = FALSE, linewidth = 0.5) +
  geom_point(aes(fill = group), shape = 21, size = 1.3, stroke = 0.2, colour = "grey20") +
  geom_text(data = lab, aes(x, y, label = lab), parse = TRUE, hjust = 0, vjust = 1,
            size = pt(NM_TXT), family = NM_FONT, inherit.aes = FALSE) +
  facet_wrap(~ fit, nrow = 1) +
  scale_colour_manual(values = geno_cols, guide = "none") +
  scale_fill_manual(values = group_cols, labels = group_keys, name = NULL) +
  scale_x_nm(step = 1, top = XL[2], bottom = XL[1], minor = TRUE, labels = lab_signed) +
  scale_y_nm(step = 1, top = YT, bottom = YB, minor = TRUE, labels = lab_signed) +
  coord_equal() +
  labs(x = "OXPHOS priority (mitoPPS, z score)",
       y = expression(italic("Bbc3:Bcl2l1") ~ "(log2 ratio, z score)")) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 1.5))) +
  theme_nm(legend = "bottom") +
  theme(strip.text = element_text(size = NM_TXT, hjust = 0, lineheight = 0.9,
                                  margin = margin(0, 0, 0.8, 0, "mm")),
        panel.spacing.x = unit(2.5, "mm"), legend.key.size = unit(2.4, "mm"),
        legend.margin = margin(-1, 0, 0, 0, "mm"))

save_panel(p, fig = "Fig1", panel = "i", name = "coupling", width = 62, height = 56,
           dir = here::here("outputs", "natmetab", "proposals"))
cat(sprintf("\n1i: unadjusted slopes WT %+.2f, MYC+ %+.2f; difference %+.2f, P %.2g, permutation P %.3f\n",
            dr("unadjusted", "neg"), dr("unadjusted", "pos"), i_u$interaction, i_u$p, pp("none")))
cat(sprintf("    adjusted   slopes WT %+.2f, MYC+ %+.2f; difference %+.2f, P %.2g, permutation P %.3f\n",
            dr("adjusted", "neg"), dr("adjusted", "pos"), i_a$interaction, i_a$p, pp("epi + imm")))
