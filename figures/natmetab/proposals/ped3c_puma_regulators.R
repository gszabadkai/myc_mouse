# =============================================================================
# PROPOSAL -- Extended Data Fig. 3c: the PUMA regulators, ranked by how far their
# 12-week MYC effect departs from MYC's overall fade
# -----------------------------------------------------------------------------
# The author's plan (2026-10-14): the FOXO3 argument goes to Extended Data, drawn
# over the old Fig. S2D (alt) roster -- script 42's p53 arm (9) and other PUMA
# inducers (12: the ISR, E2f1, Trp73 and the FOXO family) -- with Bbc3 added, as a
# dot plot like ED 3b, the genes at P < 0.05 labelled. Trp73 and Cdkn2a are below
# the expression floor (script 55), so 20 genes are drawn.
#
# THE AXIS is script 55's departure z, UNSIGNED: negative = the 12W MYC effect is
# lower than MYC's overall fade predicts. No function sign: Mdm2 restrains p53 and
# the others induce PUMA.
#
# THE LABELS (author: label all four, 2026-10-14): raw interaction P < 0.05, the
# DESeq2 test of record (the departure's own nominal P agrees, script 55). None
# survives correction across the 20. Bbc3 and Foxo3 depart BELOW the fade; the p53
# targets Zmat3 and Eda2r depart ABOVE it, through a fall in the WT gland while
# the MYC+ gland is flat -- so on the MYC+ timeline the p53 arm does not move.
# Bbc3's licence is pre-specification; Foxo3, Zmat3 and Eda2r are exploratory.
#
# Reads: results/departure_from_myc_fade.rds (script 55)
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

dep  <- readRDS(here::here("results", "departure_from_myc_fade.rds"))
BAND <- dep$band$q_abs_z[dep$band$reference == "responsive"]
r <- as.data.frame(dep$rosters)
r <- r[r$roster == "PUMA regulators" & is.na(r$dropped), ]
r$sig <- r$int_p < 0.05
stopifnot(nrow(r) == 20L,
          identical(sort(r$gene[r$sig]), sort(dep$defs$text$regulators_p05)),
          all(r$int_p_bh >= 0.05),
          all(r$z[r$gene %in% c("Bbc3", "Foxo3")] < 0), all(r$z[r$gene %in% c("Zmat3", "Eda2r")] > 0))

ROLE <- c(PUMA = "PUMA", `p53 arm` = "p53 arm", `PUMA inducer` = "Other PUMA inducers")
r$role <- factor(ROLE[r$role], levels = ROLE)
r$gene_f <- factor(r$gene, levels = r$gene[order(-r$z)])        # most negative on top
r$plab <- ifelse(r$sig, fmt_p(pmax(r$int_p, 0)), NA_character_)
r$lx <- ifelse(r$z < 0, r$z + 0.25, r$z - 0.25)
r$hj <- ifelse(r$z < 0, 0, 1)
lo <- min(r$z, -BAND); hi <- max(r$z, BAND); pad <- 0.04 * (hi - lo)

p <- ggplot(r, aes(z, gene_f)) +
  annotate("rect", xmin = -BAND, xmax = BAND, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
  geom_point(aes(shape = role, fill = sig), size = 1.8, stroke = 0.3, colour = "grey15") +
  geom_text(data = r[r$sig, ], aes(x = lx, label = plab, hjust = hj), parse = TRUE,
            size = pt(NM_TXT_MIN), family = NM_FONT) +
  scale_shape_manual(values = c(PUMA = 23, `p53 arm` = 21, `Other PUMA inducers` = 24), name = NULL) +
  scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "grey80"), breaks = c("TRUE", "FALSE"),
                    labels = c(`TRUE` = expression("Interaction" ~ italic(P) < "0.05 (raw)"),
                               `FALSE` = expression(italic(P) >= "0.05")), name = NULL) +
  scale_x_nm(step = 1, top = hi + pad, bottom = lo - pad, minor = TRUE, labels = lab_signed) +
  labs(x = expression(atop("Difference from the overall MYC effect (" * italic(z) * ")",
                           "" %<-% "lower than the fade" ~~~~~~ "higher" %->% "")),
       y = NULL) +
  guides(shape = guide_legend(order = 1, override.aes = list(size = 1.6, fill = "grey80")),
         fill  = guide_legend(order = 2, override.aes = list(shape = 21, size = 1.6))) +
  theme_nm(legend = "bottom") +
  theme(axis.text.y = element_text(face = "italic"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.box = "vertical", legend.spacing.y = unit(0.3, "mm"),
        legend.key.size = unit(2.4, "mm"))

save_panel(p, fig = "EDFig3", panel = "c", name = "puma_regulators", width = 72, height = 86,
           dir = here::here("outputs", "natmetab", "proposals"))

cat("\nED 3c -- for the legend (departure z, script 55; raw interaction P, DESeq2; BH within the 20)\n")
o <- r[order(r$z), ]
for (i in seq_len(nrow(o))) cat(sprintf("  %-10s %-20s z %+.2f  interaction P %.3g  BH %.2f\n",
  o$gene[i], as.character(o$role[i]), o$z[i], o$int_p[i], o$int_p_bh[i]))
