# =============================================================================
# fig1h_bcl2_fade.R -- Fig. 1h: the BCL-2 family, ranked by how far each gene's
# 12-week MYC effect departs from MYC's overall fade
# -----------------------------------------------------------------------------
# MOVED FROM ED FIG. 3b (author, 2026-10-16): "I would move current Ed Fig. 3b in
# the place of Fig. 1h ... without the downstream effectors, only the Bcl-2
# family", at the old Fig. 1h's size, 52 x 56 mm (the row with g and i), in Fig.
# 1's 5/6-pt type. The downstream genes stay in Extended Data as ED Fig. 3a
# (edfig3a_fade_downstream.R, the old ED 3b's script); the old Fig. 1h, the
# ratios, is ED Fig. 3b (edfig3b_ratio_expected.R).
#
# THE AXIS: script 55's departure z, signed by function -- positive = toward
# survival: a pro-apoptotic transcript lower, or an anti-apoptotic one higher, at
# 12W than MYC's overall fade (0.487 x the 6W effect) predicts.
# THE LINES (author, 2026-10-14): two dashed pairs, each labelled along its length,
# for the empirical two-sided P against the 8,774 genes MYC moved at 6W (script
# 55's p_emp_resp: the fraction of them that depart at least as far):
#   P = 0.05 at |z| = 1.73 (their 95th percentile): every other gene is inside;
#   P = 0.01 at |z| = 2.38 (their 99th percentile): Bbc3 alone is beyond it.
# SHAPE: circle pro-apoptotic, square anti-apoptotic.
# FILL: MYC's 6-week effect, espresso-white-mint, limits +/-0.5 (author,
# 2026-10-14); Bid (+0.61) is drawn at the end colour, and the legend names it. The
# bar has no frame, as Fig. 1g's key (author, 2026-10-15).
# One group, so no strip: the legend says these are the BCL-2 family.
#
# Reads: results/departure_from_myc_fade.rds (script 55)
# Output: outputs/natmetab/Fig1/Fig1h_bcl2_fade.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

TXT   <- 5      # pt; gene and tick labels, the line labels, the key
TITLE <- 6      # pt; the axis title and the colour bar's title

dep    <- readRDS(here::here("results", "departure_from_myc_fade.rds"))
BAND   <- dep$band$q_abs_z[dep$band$reference == "responsive"]     # P = 0.05
BAND99 <- dep$band$q99_abs_z[dep$band$reference == "responsive"]   # P = 0.01
LINES  <- c(BAND, BAND99)
r <- as.data.frame(dep$rosters)
r <- r[r$roster == "BCL-2 family" & is.na(r$dropped), ]
r$fn <- ifelse(r$fn_sign > 0, "anti", "pro")

# --- the record, and the claims ------------------------------------------------------------
stopifnot(nrow(r) == 16L, round(BAND, 2) == 1.73, round(BAND99, 2) == 2.38,
          identical(r$gene[r$outside_band], "Bbc3"), identical(r$gene[r$outside_99], "Bbc3"),
          r$gene[which.max(r$z_function)] == "Bbc3",
          # every gene is beyond both pairs or inside both: the labels need no third class
          !any(abs(r$z) > BAND & abs(r$z) <= BAND99), all(r$p_emp_resp[abs(r$z) <= BAND] >= 0.05))

r$gene_f <- factor(r$gene, levels = r$gene[order(r$z_function)])     # Bbc3 on top
lo <- min(r$z_function, -BAND99); hi <- max(r$z_function, BAND99); pad <- 0.04 * (hi - lo)
# the line labels: along each line, just outside it, half-way down
line_lab <- data.frame(x = c(-LINES - 0.2, LINES + 0.2), y = (nrow(r) + 1) / 2,
                       label = rep(c("italic(P)=='0.05'", "italic(P)=='0.01'"), 2))
FILL_LIM <- 0.5
beyond <- r[abs(r$lfc6) > FILL_LIM, ]                    # drawn at the end colour (squished)
stopifnot(identical(beyond$gene, "Bid"))

p <- ggplot(r, aes(z_function, gene_f)) +
  geom_vline(xintercept = c(-LINES, LINES), linetype = "22", linewidth = NM_LINE, colour = "grey45") +
  geom_text(data = line_lab, aes(x, y, label = label), inherit.aes = FALSE, parse = TRUE,
            angle = 90, size = pt(TXT), family = NM_FONT) +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
  geom_point(aes(shape = fn, fill = lfc6), size = 1.5, stroke = 0.3, colour = "grey15") +
  scale_shape_manual(values = c(pro = 21, anti = 22), breaks = c("pro", "anti"), name = NULL,
                     labels = c(pro = "Pro­apoptotic", anti = "Anti­apoptotic")) +
  scale_fill_div(c(-FILL_LIM, FILL_LIM), name = "MYC+ vs WT, 6W\n(log2 fold change)",
                 breaks = c(-FILL_LIM, 0, FILL_LIM)) +
  scale_x_nm(step = 1, top = hi + pad, bottom = lo - pad, minor = TRUE, labels = lab_signed) +
  scale_y_discrete() +                                    # rows; the line labels sit between them
  labs(x = NULL, y = NULL) +                              # the x title is a strip of its own, below
  guides(shape = guide_legend(order = 1, ncol = 1, override.aes = list(size = 1.4, fill = "grey80")),
         fill  = guide_colourbar(order = 2, direction = "horizontal", barwidth = unit(14, "mm"),
                                 barheight = unit(1.6, "mm"), title.position = "top",
                                 ticks.colour = "grey30")) +            # no frame, as Fig. 1g
  theme_nm(legend = "bottom", txt = TXT, txt_title = TITLE) +
  theme(axis.text.y = element_text(face = "italic"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.box = "horizontal", legend.spacing.x = unit(3, "mm"),
        legend.key.size = unit(2.2, "mm"),
        legend.title = element_text(size = TITLE, lineheight = 0.9))

# The x title: two lines LEAD mm apart, centre to centre (6 pt at 1.15), closer than
# plotmath's atop() sets them, so it is a strip of its own under the graph, as Fig.
# 1b's y title. Patchwork aligns the strip's panel with the graph's, so x = 0.5 centres
# both lines under the axis (38.9 and 35.5 mm against about 40). The key goes below.
LEAD <- 2.4; TITLE_H <- 5                                 # mm
xtitle <- ggplot() +
  annotate("text", x = 0.5, y = 1 - (1.2 + c(0, LEAD)) / TITLE_H, parse = TRUE,
           label = c('"Difference from the overall MYC effect ("*italic(z)*")"',
                     '""%<-%"toward death"~~~~~~"toward survival"%->%""'),
           size = pt(TITLE), family = NM_FONT, colour = "black") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = NM_FONT) + theme(plot.margin = margin(0, 0, 0, 0))

panel <- stack_parts(p, xtitle, heights = grid::unit(c(1, TITLE_H), c("null", "mm")),
                     collect_guides = TRUE)
save_panel(panel, fig = "Fig1", panel = "h", name = "bcl2_fade", width = 52, height = 56)

cat(sprintf("\nFig. 1h -- for the legend: the BCL-2 family (%d genes); fade rate %.4f; dashed lines at |z| = %.2f and\n",
            nrow(r), dep$defs$rate, BAND))
cat(sprintf("  %.2f, the 95th and 99th percentiles over the 8,774 genes MYC moved at 6W (empirical two-sided\n", BAND99))
cat(sprintf("  P = 0.05 and 0.01); Bbc3 beyond both (P %.4f); every other gene P >= %.2f\n",
            r$p_emp_resp[r$gene == "Bbc3"], min(r$p_emp_resp[r$gene != "Bbc3"])))
cat(sprintf("  fill limits +/-%.1f; beyond them, at the end colour: %s\n", FILL_LIM,
            paste(sprintf("%s %+.2f", beyond$gene, beyond$lfc6), collapse = ", ")))
o <- r[order(-r$z_function), ]
for (i in seq_len(nrow(o))) cat(sprintf("  %-8s %-15s %-4s 6W %+.2f  12W %+.2f  z %+.2f  P %.3f%s\n",
  o$gene[i], o$role[i], o$fn[i], o$lfc6[i], o$lfc12[i], o$z_function[i], o$p_emp_resp[i],
  if (o$outside_band[i]) "  OUTSIDE" else ""))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p)
  r[order(-r$z_function), c("gene", "role", "lfc6", "lfc12", "z_function", "p_emp_resp")]
}
