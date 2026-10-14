# =============================================================================
# fig1i_coupling.R -- Fig. 1i: OXPHOS priority against the PUMA:BCL-XL balance,
# animal by animal, adjusted for immune content
# -----------------------------------------------------------------------------
# PROMOTED from proposals/ on 2026-10-14 (author: "promote them to final panels").
# The author's plan (2026-10-14): "along with OXPHOS deprioritisation" needs a
# panel -- the adjusted half of the coupling, the WT line dashed and thinner, the
# MYC+ line's R2 on the panel. P values in the legend (author): the slope
# difference +0.78, P 0.0052 parametric and 0.085 against the within-age
# permutation null. The unadjusted fit, and both fits side by side with all four P
# values, are ED Fig. 3d (edfig3d_coupling_both_fits.R): the impasse of the 2026-09-21 ruling.
#
# WHAT THE R2 IS: the squared correlation of the drawn MYC+ points (n = 12), whose
# y values are partial residuals -- the ratio with the shared epithelial and immune
# terms removed. Without the adjustment it is 0.35. The MYC+ line spans both ages,
# so part of it is the 6W-to-12W difference (the verification note, section 5).
#
# Reads: results/two_timeline_verification.rds (script 54)
# Output: outputs/natmetab/Fig1/Fig1i_coupling.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

tv <- readRDS(here::here("results", "two_timeline_verification.rds"))
pd <- as.data.frame(tv$coupling_panel)
ln <- as.data.frame(tv$coupling_lines)
cf <- as.data.frame(tv$coupling_fits)
ix <- as.data.frame(tv$coupling_interaction)
pm <- as.data.frame(tv$coupling_perm)
stopifnot(nrow(pd) == 24L, all(table(pd$genotype) == 12L))

# --- the drawn lines ARE the reported pooled adjusted fit -----------------------------
FP  <- "P pooled, shared covariates: epi + imm"
la  <- ln[ln$fit == "adjusted", ]
sl  <- function(g) cf$slope[cf$axis == "ox_ppd" & cf$fit == FP & cf$genotype == g]
i_a <- ix[ix$axis == "ox_ppd" & ix$covariates == "epi + imm", ]
p_perm <- pm$p_emp[pm$axis == "ox_ppd" & pm$covariates == "epi + imm"]
stopifnot(nrow(la) == 2L,
          abs(la$slope[la$genotype == "neg"] - sl("neg")) < 1e-10,
          abs(la$slope[la$genotype == "pos"] - sl("pos")) < 1e-10,
          abs((sl("pos") - sl("neg")) - i_a$interaction) < 1e-10,
          i_a$p < 0.01, p_perm >= 0.05)                    # the impasse, as ruled

pos <- pd$genotype == "pos"
R2  <- stats::cor(pd$x[pos], pd$y_adjusted[pos])^2
stopifnot(round(R2, 2) == 0.73)

xr  <- tapply(pd$x, pd$genotype, range)
seg <- do.call(rbind, lapply(c("neg", "pos"), function(g) {
  x <- xr[[g]]; l <- la[la$genotype == g, ]
  data.frame(genotype = g, x = x[1], xend = x[2], y = l$intercept + l$slope * x[1],
             yend = l$intercept + l$slope * x[2])
}))
pd$group <- factor(pd$group, levels = order_by_geno)
XL <- c(-1.5, 3); YL <- c(-2, 2.5)
stopifnot(all(pd$x > XL[1], pd$x < XL[2], pd$y_adjusted > YL[1], pd$y_adjusted < YL[2]))

p <- ggplot(pd, aes(x, y_adjusted)) +
  geom_hline(yintercept = 0, linewidth = NM_LINE, colour = "grey85") +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = "grey85") +
  geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend, colour = genotype,
                               linetype = genotype, linewidth = genotype), inherit.aes = FALSE) +
  geom_point(aes(fill = group), shape = 21, size = 1.4, stroke = 0.2, colour = "grey20") +
  annotate("text", x = XL[1] + 0.15, y = YL[2] - 0.1, hjust = 0, vjust = 1, parse = TRUE,
           label = sprintf('"MYC+"~italic(R)^displaystyle(2)=="%.2f"', R2),
           size = pt(NM_TXT), family = NM_FONT) +
  scale_colour_manual(values = geno_cols, labels = geno_labels, name = NULL) +
  scale_linetype_manual(values = c(neg = "22", pos = "solid"), labels = geno_labels, name = NULL) +
  scale_linewidth_manual(values = c(neg = 0.35, pos = 0.55), labels = geno_labels, name = NULL) +
  scale_fill_manual(values = group_cols, labels = group_keys, name = NULL) +
  scale_x_nm(step = 1, top = XL[2], bottom = XL[1], minor = TRUE, labels = lab_signed) +
  scale_y_nm(step = 1, top = YL[2], bottom = YL[1], minor = TRUE, labels = lab_signed) +
  coord_equal() +
  labs(x = "OXPHOS priority (mitoPPS, z score)",
       # "adjusted for immune content", as ED 3d (author, 2026-10-14): "adjusted" alone, and
       # "epithelial content", read as adjusting for the MECs themselves
       y = expression(atop(italic("Bbc3:Bcl2l1") ~ "(z score),", "adjusted for immune content"))) +
  guides(fill = guide_legend(order = 1, nrow = 2, override.aes = list(size = 1.5)),
         colour = guide_legend(order = 2, nrow = 2), linetype = guide_legend(order = 2, nrow = 2),
         linewidth = guide_legend(order = 2, nrow = 2)) +
  theme_nm(legend = "bottom") +
  theme(legend.key.size = unit(2.4, "mm"), legend.key.width = unit(4, "mm"),
        legend.box = "horizontal")

save_panel(p, fig = "Fig1", panel = "i", name = "coupling", width = 80, height = 60)

cat(sprintf("\n1i -- for the legend: adjusted for immune content -- the pooled model's shared covariates are\n"))
cat(sprintf("  script 42's immune and epithelial marker scores, which estimate each unsorted prep's\n"))
cat(sprintf("  non-epithelial share (the immune score carries most of the change);\n"))
cat(sprintf("  slopes WT %+.2f, MYC+ %+.2f; difference %+.2f (SE %.2f), P = %.2g parametric, %.3f against the\n",
            sl("neg"), sl("pos"), i_a$interaction, i_a$se, i_a$p, p_perm))
cat(sprintf("  within-age permutation null (5,000 draws; null median %+.2f); MYC+ R2 of the drawn points %.2f\n",
            pm$null_median[pm$covariates == "epi + imm"], R2))
