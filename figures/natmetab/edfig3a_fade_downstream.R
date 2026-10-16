# =============================================================================
# edfig3a_fade_downstream.R -- Extended Data Fig. 3a: the downstream genes, ranked
# by how far their 12-week MYC effect departs from MYC's overall fade
# -----------------------------------------------------------------------------
# SPLIT FROM THE OLD ED FIG. 3b (author, 2026-10-16): its BCL-2 family half moved to
# Fig. 1h (fig1h_bcl2_fade.R); "the downstream effectors remain separate as an ED
# Fig. panel", this one. The old ED 3b (edfig3b_fade_ranking.R) was promoted from
# proposals/ on 2026-10-14; this file is it, renamed with git mv and cut to the
# downstream genes. The Extended Data type stays 6/7 pt (author, 2026-10-16).
#
# THE AXIS: script 55's departure z, signed by function -- positive = toward
# survival: a pro-apoptotic transcript lower, or an anti-apoptotic one higher, at
# 12W than MYC's overall fade (0.487 x the 6W effect) predicts.
# THE LINES, as Fig. 1h (author, 2026-10-14): two dashed pairs, each labelled along
# its length, for the empirical two-sided P against the 8,774 genes MYC moved at 6W
# (script 55's p_emp_resp): P = 0.05 at |z| = 1.73 and P = 0.01 at |z| = 2.38. With
# Apaf1 out, every downstream gene is inside both (P >= 0.28).
# SHAPE: circle pro-apoptotic, square anti-apoptotic (Diablo, SMAC, is
# pro-apoptotic: the IAP antagonist; corrected in script 55 on 2026-10-14).
# FILL LIMITS +/-0.5 (author, 2026-10-14): Birc5 (+1.16) and Htra2 (+1.49) are drawn
# at the end colour, and the legend names them.
#
# APAF1 IS LEFT OUT (author, 2026-10-14). It sits outside the band (z -2.14 on
# the raw axis), but MYC had already lowered it at 6W (-0.24 log2), the change
# across the window does not differ by genotype (interaction P 0.39), its
# apoptosome module relaxes rather than tightens under MYC (script 17), and it acts
# after the BCL-2 decision this figure is about. The legend must say it is left out
# and give its value, so the omission is visible (printed below).
#
# Reads: results/departure_from_myc_fade.rds (script 55)
# Output: outputs/natmetab/EDFig3/EDFig3a_fade_downstream.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

dep    <- readRDS(here::here("results", "departure_from_myc_fade.rds"))
BAND   <- dep$band$q_abs_z[dep$band$reference == "responsive"]     # P = 0.05
BAND99 <- dep$band$q99_abs_z[dep$band$reference == "responsive"]   # P = 0.01
LINES  <- c(BAND, BAND99)
r <- as.data.frame(dep$rosters)
r <- r[r$roster == "downstream" & is.na(r$dropped), ]
r$fn <- ifelse(r$fn_sign > 0, "anti", "pro")
OMIT <- "Apaf1"
om <- r[r$gene %in% OMIT, ]
stopifnot(nrow(om) == 1L, om$outside_band, om$int_p > 0.05)
r <- r[!r$gene %in% OMIT, ]

# --- the record is the corrected one, and the claims hold --------------------------------
stopifnot(nrow(r) == 7L, round(BAND, 2) == 1.73, round(BAND99, 2) == 2.38,
          r$fn[r$gene == "Diablo"] == "pro",
          setequal(r$gene[r$fn == "anti"], c("Xiap", "Birc5")),
          !any(r$outside_band), all(r$p_emp_resp >= 0.05))   # with Apaf1 out, none departs

r$gene_f <- factor(r$gene, levels = r$gene[order(r$z_function)])     # highest on top
lo <- min(r$z_function, -BAND99); hi <- max(r$z_function, BAND99); pad <- 0.04 * (hi - lo)
# the line labels: along each line, just outside it, half-way down
line_lab <- data.frame(x = c(-LINES - 0.15, LINES + 0.15), y = (nrow(r) + 1) / 2,
                       label = rep(c("italic(P)=='0.05'", "italic(P)=='0.01'"), 2))
FILL_LIM <- 0.5
beyond <- r[abs(r$lfc6) > FILL_LIM, ]                    # drawn at the end colour (squished)
stopifnot(setequal(beyond$gene, c("Birc5", "Htra2")))

p <- ggplot(r, aes(z_function, gene_f)) +
  geom_vline(xintercept = c(-LINES, LINES), linetype = "22", linewidth = NM_LINE, colour = "grey45") +
  geom_text(data = line_lab, aes(x, y, label = label), inherit.aes = FALSE, parse = TRUE,
            angle = 90, size = pt(NM_TXT_MIN), family = NM_FONT) +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
  geom_point(aes(shape = fn, fill = lfc6), size = 1.8, stroke = 0.3, colour = "grey15") +
  scale_shape_manual(values = c(pro = 21, anti = 22), breaks = c("pro", "anti"), name = NULL,
                     labels = c(pro = "Pro­apoptotic", anti = "Anti­apoptotic")) +
  scale_fill_div(c(-FILL_LIM, FILL_LIM), name = "MYC+ vs WT, 6W\n(log2 fold change)",
                 breaks = c(-FILL_LIM, 0, FILL_LIM)) +
  scale_x_nm(step = 1, top = hi + pad, bottom = lo - pad, minor = TRUE, labels = lab_signed) +
  scale_y_discrete() +                                    # rows; the line labels sit between them
  labs(x = expression(atop("Difference from the overall MYC effect (" * italic(z) * ")",
                           "" %<-% "toward death" ~~~~~~ "toward survival" %->% "")),
       y = NULL) +
  guides(shape = guide_legend(order = 1, ncol = 1, override.aes = list(size = 1.6, fill = "grey80")),
         fill  = guide_colourbar(order = 2, direction = "horizontal", barwidth = unit(16, "mm"),
                                 barheight = unit(2, "mm"), title.position = "top",
                                 ticks.colour = "grey30", frame.colour = "grey30")) +
  theme_nm(legend = "right") +
  theme(axis.text.y = element_text(face = "italic"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.box = "vertical", legend.spacing.y = unit(1.5, "mm"),
        legend.key.size = unit(2.4, "mm"),
        legend.title = element_text(size = NM_TXT, lineheight = 0.9))

save_panel(p, fig = "EDFig3", panel = "a", name = "fade_downstream", width = 86, height = 44)

cat(sprintf("\nED 3a -- for the legend: the downstream genes (%d drawn); fade rate %.4f; dashed lines at |z| = %.2f\n",
            nrow(r), dep$defs$rate, BAND))
cat(sprintf("  and %.2f, the 95th and 99th percentiles over the 8,774 genes MYC moved at 6W (empirical two-sided\n", BAND99))
cat(sprintf("  P = 0.05 and 0.01); every drawn gene P >= %.2f\n", min(r$p_emp_resp)))
cat(sprintf("  fill limits +/-%.1f; beyond them, at the end colour: %s\n", FILL_LIM,
            paste(sprintf("%s %+.2f", beyond$gene, beyond$lfc6), collapse = ", ")))
cat(sprintf("  NOT DRAWN: %s (%s), 6W %+.2f, 12W %+.2f, z %+.2f (empirical P %.3f: between the P = 0.05\n    and P = 0.01 lines), interaction P %.2f\n",
            om$gene, om$role, om$lfc6, om$lfc12, om$z_function, om$p_emp_resp, om$int_p))
o <- r[order(-r$z_function), ]
for (i in seq_len(nrow(o))) cat(sprintf("  %-8s %-28s %-4s 6W %+.2f  12W %+.2f  z %+.2f  P %.3f\n",
  o$gene[i], o$role[i], o$fn[i], o$lfc6[i], o$lfc12[i], o$z_function[i], o$p_emp_resp[i]))
