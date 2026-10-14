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
# THE LINES, as ED 3b (author, 2026-10-14): two dashed pairs, labelled P = 0.05
# (|z| = 1.73, the 95th percentile over the 8,774 genes MYC moved at 6W) and P = 0.01
# (|z| = 2.38, the 99th), so that the 16 genes inside are seen not to depart
# (every one P >= 0.10). No gene falls between the pairs.
#
# THE LABELS (author: label all four, 2026-10-14). They now print the SAME P as the
# lines: the departure's empirical P against the MYC-responsive genes (script 55,
# p_emp_resp; all four < 0.01). They printed the raw interaction P before, and a
# P = 0.01 line beside those would put Zmat3 beyond it with "P = 0.046" printed
# next to it. The interaction P (Foxo3 0.0074, Bbc3 0.0081, Eda2r 0.014, Zmat3
# 0.046) goes to the legend; the four genes are the same either way. The nominal
# P of the same departure is larger (Bbc3 0.013, Foxo3 0.011): departures among
# MYC-responsive genes are narrower than a standard normal, and the line is drawn
# against them. Bbc3 and Foxo3 depart BELOW the fade; the p53
# targets Zmat3 and Eda2r depart ABOVE it, through a fall in the WT gland while
# the MYC+ gland is flat -- so on the MYC+ timeline the p53 arm does not move.
# Bbc3's licence is pre-specification; Foxo3, Zmat3 and Eda2r are exploratory.
#
# Reads: results/departure_from_myc_fade.rds (script 55)
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

dep  <- readRDS(here::here("results", "departure_from_myc_fade.rds"))
BAND   <- dep$band$q_abs_z[dep$band$reference == "responsive"]     # P = 0.05
BAND99 <- dep$band$q99_abs_z[dep$band$reference == "responsive"]   # P = 0.01
LINES  <- c(BAND, BAND99)
r <- as.data.frame(dep$rosters)
r <- r[r$roster == "PUMA regulators" & is.na(r$dropped), ]
r$sig <- r$outside_99                                   # beyond the lines: empirical P < 0.01
stopifnot(nrow(r) == 20L, round(BAND99, 2) == 2.38,
          identical(sort(r$gene[r$sig]), sort(dep$defs$text$outside99_regulators)),
          all(r$p_emp_resp[r$sig] < 0.01),
          identical(sort(r$gene[r$int_p < 0.05]), sort(r$gene[r$sig])),     # the same four genes
          # two classes are enough: beyond both pairs (P < 0.01) or inside both (P > 0.05)
          !any(abs(r$z) > BAND & abs(r$z) <= BAND99), all(r$p_emp_resp[!r$sig] > 0.05),
          all(r$int_p_bh >= 0.05),
          all(r$z[r$gene %in% c("Bbc3", "Foxo3")] < 0), all(r$z[r$gene %in% c("Zmat3", "Eda2r")] > 0))

ROLE <- c(PUMA = "PUMA", `p53 arm` = "p53 arm", `PUMA inducer` = "Other PUMA inducers")
r$role <- factor(ROLE[r$role], levels = ROLE)
r$gene_f <- factor(r$gene, levels = r$gene[order(-r$z)])        # most negative on top
r$plab <- ifelse(r$sig, fmt_p(pmax(r$p_emp_resp, 0)), NA_character_)
r$lx <- ifelse(r$z < 0, r$z + 0.25, r$z - 0.25)
r$hj <- ifelse(r$z < 0, 0, 1)
lo <- min(r$z, -BAND99); hi <- max(r$z, BAND99); pad <- 0.04 * (hi - lo)
line_lab <- data.frame(x = c(-LINES - 0.13, LINES + 0.13), y = (nrow(r) + 1) / 2,
                       label = rep(c("italic(P)=='0.05'", "italic(P)=='0.01'"), 2))

p <- ggplot(r, aes(z, gene_f)) +
  geom_vline(xintercept = c(-LINES, LINES), linetype = "22", linewidth = NM_LINE, colour = "grey45") +
  geom_text(data = line_lab, aes(x, y, label = label), inherit.aes = FALSE, parse = TRUE,
            angle = 90, size = pt(NM_TXT_MIN), family = NM_FONT) +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
  geom_point(aes(shape = role, fill = sig), size = 1.8, stroke = 0.3, colour = "grey15") +
  # on a white backing: the band's dashed line runs behind Foxo3's and Bbc3's labels
  geom_label(data = r[r$sig, ], aes(x = lx, label = plab, hjust = hj), parse = TRUE,
             size = pt(NM_TXT_MIN), family = NM_FONT, fill = "white", border.colour = NA,
             label.padding = unit(0.3, "mm"), label.r = unit(0, "mm")) +
  scale_shape_manual(values = c(PUMA = 23, `p53 arm` = 21, `Other PUMA inducers` = 24), name = NULL) +
  scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "grey80"), breaks = c("TRUE", "FALSE"),
                    labels = c(`TRUE` = expression(italic(P) < "0.01"),
                               `FALSE` = expression(italic(P) > "0.05")), name = NULL) +
  scale_x_nm(step = 1, top = hi + pad, bottom = lo - pad, minor = TRUE, labels = lab_signed) +
  scale_y_discrete() +                                    # rows; the line labels sit between them
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

save_panel(p, fig = "EDFig3", panel = "c", name = "puma_regulators", width = 86, height = 92,
           dir = here::here("outputs", "natmetab", "proposals"))

cat(sprintf("\nED 3c -- for the legend: departure z (script 55); dashed lines at |z| = %.2f and %.2f, the 95th and\n", BAND, BAND99))
cat("  99th percentiles over the 8,774 MYC-responsive genes (empirical two-sided P = 0.05 and 0.01);\n")
cat("  printed P = the same empirical P (p_emp_resp);\n")
cat("  interaction P (DESeq2, raw; BH within the 20) for the legend text:\n")
o <- r[order(r$z), ]
for (i in seq_len(nrow(o))) cat(sprintf("  %-10s %-20s z %+.2f  empirical P %.4f  nominal P %.3f  interaction P %.3g  BH %.2f\n",
  o$gene[i], as.character(o$role[i]), o$z[i], o$p_emp_resp[i], o$p_norm[i], o$int_p[i], o$int_p_bh[i]))
