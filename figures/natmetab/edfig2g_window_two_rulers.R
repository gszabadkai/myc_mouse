# =============================================================================
# edfig2g_window_two_rulers.R -- Extended Data Fig. 2g: across the window, every
# MitoPathway's change in expression against its change in priority
# -----------------------------------------------------------------------------
# The text: "both the WT and the MYC+ gland demoted OXPHOS in priority and in
# expression, while fatty acid oxidation and the carnitine shuttle ... gained
# priority (Fig. 1f, ED Fig. 2g)." Fig. 1f draws six pathways; this is all 143,
# one facet per genotype. Lower left: falls on both rulers; upper right: gains on
# both.
#
# FORM: two signed quantities per pathway, two windows -> two scatters sharing
#   both axes (data-to-viz, "Consistency between charts": the genotype colours of
#   Fig. 1f). Named: the three pathways the sentence names.
#
# SIZE: 84 x 64 mm, beside ED 2f (92 mm) in one 180-mm Extended Data row.
#
# Reads:  results/background_vs_myc.rds (script 40: $ruler)
# Output: outputs/natmetab/EDFig2/EDFig2g_window_two_rulers.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

ruler <- as.data.frame(readRDS(here::here("results", "background_vs_myc.rds"))$ruler)
r143  <- ruler[!ruler$is_mtdna, ]
stopifnot(nrow(r143) == 143L)

d <- rbind(data.frame(pathway = r143$pathway, geno = "neg", x = r143$c_tn, y = r143$p_tn),
           data.frame(pathway = r143$pathway, geno = "pos", x = r143$c_tp, y = r143$p_tp))
d$geno <- factor(d$geno, levels = c("neg", "pos"))

# the three pathways the sentence names (Fig. 1f draws three more), on two lines
named <- c("OXPHOS subunits" = "OXPHOS\nsubunits", "Carnitine shuttle" = "Carnitine\nshuttle",
           "Fatty acid oxidation" = "Fatty acid\noxidation")
stopifnot(all(names(named) %in% r143$pathway))
hl <- d[d$pathway %in% names(named), ]
# every pathway goes to the label layer, unnamed ones with an empty label, so the
# labels keep off all points, not just the named ones
d$label <- ifelse(d$pathway %in% names(named), unname(named[d$pathway]), "")

rr <- vapply(c("neg", "pos"), function(g) stats::cor(d$x[d$geno == g], d$y[d$geno == g]),
             numeric(1))
# small pads: the x axis then ends on the minor tick at +0.6, so the last label of the
# left facet (+0.4) cannot run into the first of the right one (-0.8)
xpad <- 0.03 * diff(range(d$x)); ypad <- 0.04 * diff(range(d$y))
xl <- range(d$x) + c(-xpad, xpad); yl <- range(d$y) + c(-ypad, ypad)
ann <- data.frame(geno = factor(c("neg", "pos"), levels = c("neg", "pos")),
                  x = xl[1] + 0.03, y = yl[2] - 0.02,
                  lab = sprintf("italic(r)=='%.2f'", rr))

p <- ggplot(d, aes(x, y)) +
  geom_hline(yintercept = 0, linewidth = NM_LINE, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = "grey80") +
  geom_point(aes(colour = geno), size = 0.8, alpha = 0.5, stroke = 0) +
  geom_point(data = hl, aes(fill = geno), shape = 21, size = 1.3, colour = "black", stroke = 0.35) +
  ggrepel::geom_text_repel(aes(label = label), size = pt(NM_TXT), family = NM_FONT,
                           colour = "black", lineheight = 0.9, segment.colour = "grey50",
                           segment.size = 0.25, min.segment.length = 0, box.padding = 0.2,
                           point.padding = 0.05, max.overlaps = Inf, seed = 1) +
  geom_text(data = ann, aes(x, y, label = lab), parse = TRUE, hjust = 0, size = pt(NM_TXT)) +
  facet_wrap(~ geno, nrow = 1, labeller = as_labeller(dev_labels)) +
  scale_colour_manual(values = geno_cols, guide = "none") +
  scale_fill_manual(values = geno_cols, guide = "none") +
  scale_x_nm(step = 0.4, top = xl[2], bottom = xl[1], minor = TRUE, labels = lab_signed) +
  scale_y_nm(step = 0.2, top = yl[2], bottom = yl[1], minor = TRUE, labels = lab_signed) +
  labs(x = "Expression, 12W vs 6W (log2 fold change)", y = "Priority, 12W vs 6W (mitoPPS)") +
  theme_nm() +
  theme(strip.text = element_text(size = NM_TXT), panel.spacing.x = unit(4, "mm"))

save_panel(p, fig = "EDFig2", panel = "g", name = "window_two_rulers", width = 84, height = 64)

# --- numbers for the legend (printed, never drawn) --------------------------------
cat("\nED Fig. 2g -- for the legend\n")
cat(sprintf("  %d nuclear MitoPathways; expression = set-average raw log2 fold change; priority = mitoPPS\n",
            nrow(r143)))
cat(sprintf("  r(expression, priority): WT window %.2f, MYC+ window %.2f\n", rr[1], rr[2]))
cat("  outlined and named: the three pathways the sentence names\n")

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p)
  d[d$geno == "neg", ][order(d$x[d$geno == "neg"]), ][1:10, ]
}
