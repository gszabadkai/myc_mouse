# =============================================================================
# fig1e_myc_reallocation.R -- Fig. 1e: MYC re-weights the mitochondrial
# compartment, the same way at both ages
# -----------------------------------------------------------------------------
# The text: "Within the nuclear-encoded mitochondrial transcriptome, scored as
# pairwise ratios of pathways independent of its overall size (mitoPPS), MYC
# promoted biogenesis and biosynthetic pathways and demoted mitochondrial
# dynamics, surveillance and signalling, in the same pattern at both ages
# (Fig. 1e, ED Fig. 2e)."
#
# THE ROWS are chosen by a rule, not by eye: every row except OXPHOS subunits
# is significant at 6W (BH-adjusted P < 0.05) and keeps its sign at 12W (both
# asserted). OXPHOS subunits is drawn as the reference for Fig. 1f: MYC does
# not demote OXPHOS -- the window does. ED Fig. 2e draws all 143 pathways.
#
# FORM: ranked effects, two conditions -> a dot plot, 6W and 12W per row joined
#   by a line (data-to-viz, "Order your data"; "The Moire effect": dots, not
#   bars). Colour is the age, as in ED Fig. 1 and Fig. 1d; a filled point is
#   significant, an open one is not. The 12W point is smaller: 6W is the effect
#   the rows are chosen on, 12W the same pattern at lower amplitude.
#
# SIZE: 78 x 55 mm (author, 2026-10-15; 80 x 66 until then), beside Fig. 1d (100 x 60)
#   on the 183-mm page (83 wide until 2026-10-14, which left no gap). Type as Fig. 1d
#   (author, 2026-10-15): row labels, tick labels and the key 5 pt, the group names and
#   the axis title 6 pt.
#
# STATISTICS: script 08's Welch t-test on per-animal mitoPPS scores, 6 vs 6 mice,
#   BH-adjusted across the MitoPathways within each contrast; carried into
#   results/background_vs_myc.rds$ruler by script 40 (p_m6/p_m12 and their padj).
#
# Reads:  results/background_vs_myc.rds (script 40: $ruler)
# Output: outputs/natmetab/Fig1/Fig1e_myc_reallocation.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

ruler <- as.data.frame(readRDS(here::here("results", "background_vs_myc.rds"))$ruler)

TXT   <- 5      # pt; row labels, tick labels and the key
TITLE <- 6      # pt; the group names and the axis title

GROUPS <- c(bio = "Biogenesis and\nbiosynthesis", ox = "OXPHOS",
            dyn = "Dynamics,\nsurveillance and\nsignalling")
rows <- data.frame(
  pathway = c("Translation", "Transcription", "Chaperones", "Serine metabolism",
              "Folate and 1-C metabolism",
              "OXPHOS subunits",
              "Mitochondrial dynamics and surveillance", "Fission", "Apoptosis",
              "Calcium homeostasis", "SLC25A family"),
  # labels avoid hyphens: base pdf() draws '-' as a minus sign (FIGURE_RULES.md)
  label   = c("Mitochondrial translation", "Mitochondrial transcription", "Chaperones",
              "Serine metabolism", "Folate and 1C metabolism",
              "OXPHOS subunits",
              "Dynamics and surveillance", "Fission", "Apoptosis",
              "Calcium homeostasis", "SLC25A carriers"),
  group   = c(rep("bio", 5), "ox", rep("dyn", 5)),
  stringsAsFactors = FALSE)
k <- match(rows$pathway, ruler$pathway)
stopifnot(!anyNA(k))
rows <- cbind(rows, ruler[k, c("n_genes", "p_m6", "p_m6_padj", "p_m12", "p_m12_padj")])

# the rule the rows were chosen by
named <- rows$group != "ox"
stopifnot(all(rows$p_m6_padj[named] < 0.05),
          all(sign(rows$p_m6) == sign(rows$p_m12)),
          all(rows[!named, c("p_m6_padj", "p_m12_padj")] >= 0.05))

rows$label_n <- sprintf("%s (%d)", rows$label, rows$n_genes)
ord <- order(match(rows$group, names(GROUPS)), -rows$p_m6)        # top row first
rows$label_n <- factor(rows$label_n, levels = rev(rows$label_n[ord]))
rows$group   <- factor(GROUPS[rows$group], levels = GROUPS)

long <- rbind(
  data.frame(rows[, c("label_n", "group")], contrast = "myc_6W",  eff = rows$p_m6,
             sig = rows$p_m6_padj < 0.05),
  data.frame(rows[, c("label_n", "group")], contrast = "myc_12W", eff = rows$p_m12,
             sig = rows$p_m12_padj < 0.05))
long$contrast <- factor(long$contrast, levels = names(contrast_cols))
# 6W is the effect the row rule is built on; 12W is drawn smaller (author, 2026-10-11)
PT_SIZE <- c(myc_6W = 1.6, myc_12W = 1.1)
long$sig <- factor(ifelse(long$sig, "sig", "ns"), levels = c("sig", "ns"))

p <- ggplot(long, aes(eff, label_n)) +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
  geom_segment(data = rows, aes(x = p_m6, xend = p_m12, y = label_n, yend = label_n),
               inherit.aes = FALSE, linewidth = 0.4, colour = "grey70") +
  geom_point(aes(colour = contrast, shape = sig, size = contrast), stroke = 0.5) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = contrast_cols, labels = contrast_labels, name = NULL) +
  scale_size_manual(values = PT_SIZE, guide = "none") +
  scale_shape_manual(values = c(sig = 16, ns = 1), name = NULL,
                     labels = c(sig = expression("adjusted" ~ italic(P) < "0.05"),
                                ns  = "not significant")) +
  scale_x_nm(step = 0.2, top = max(long$eff), bottom = min(long$eff), minor = TRUE,
             labels = lab_signed) +
  labs(x = "MYC+ vs WT (mitoPPS)", y = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = unname(PT_SIZE))),
         shape  = guide_legend(order = 2, override.aes = list(size = 1.4))) +
  theme_nm(legend = "bottom", txt = TXT, txt_title = TITLE) +
  theme(axis.line.y  = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(size = TITLE, angle = 0, hjust = 0, lineheight = 0.95),
        panel.spacing.y = unit(1.5, "mm"),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "mm"),
        legend.margin = margin(0, 0, 0, 0))

save_panel(p, fig = "Fig1", panel = "e", name = "myc_reallocation", width = 78, height = 55)

# --- numbers for the legend (printed, never drawn) --------------------------------
r143 <- ruler[!ruler$is_mtdna, ]
cat("\nFig. 1e -- for the legend\n")
cat("  mitoPPS (Monzel 2025), MYC+ vs WT at each age; n = 6 mice per group;\n")
cat("  Welch t-test on per-animal scores, BH-adjusted across MitoPathways (script 08)\n")
for (i in ord) cat(sprintf("    %-38s 6W %+.2f (P %.2g)   12W %+.2f (P %.2g)\n", rows$label_n[i],
                           rows$p_m6[i], rows$p_m6_padj[i], rows$p_m12[i], rows$p_m12_padj[i]))
cat(sprintf("  across all %d MitoPathways: r(6W, 12W) = %.2f; %d significant at 6W, all %d with the same sign at 12W\n",
            nrow(r143), stats::cor(r143$p_m6, r143$p_m12), sum(r143$p_m6_padj < 0.05),
            sum(r143$p_m6_padj < 0.05 & sign(r143$p_m6) == sign(r143$p_m12))))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p)
  r143[order(-r143$p_m6), c("pathway", "tier", "p_m6", "p_m6_padj", "p_m12", "p_m12_padj")][1:20, ]
}
