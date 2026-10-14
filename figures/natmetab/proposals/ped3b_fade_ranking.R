# =============================================================================
# PROPOSAL -- Extended Data Fig. 3b: the BCL-2 family and the downstream genes,
# ranked by how far their 12-week MYC effect departs from MYC's overall fade
# -----------------------------------------------------------------------------
# The 2026-10-13 ranking (p_fade_ranking.R, version B) as the author specified on
# 2026-10-14: the points are coloured by MYC's 6-week effect on a graded scale
# (the style's diverging espresso-white-mint) instead of the black/grey split, and
# the axis reads "Difference from the overall MYC effect". It now reads script 55's
# record instead of computing anything.
#
# THE AXIS: script 55's departure z, signed by function -- positive = toward
# survival: a pro-apoptotic transcript lower, or an anti-apoptotic one higher, at
# 12W than MYC's overall fade (0.487 x the 6W effect) predicts.
# THE LINES (author, 2026-10-14): two dashed pairs, each labelled along its length,
# for the empirical two-sided P against the 8,774 genes MYC moved at 6W (script
# 55's p_emp_resp: the fraction of them that depart at least as far):
#   P = 0.05 at |z| = 1.73 (their 95th percentile) -- "more important to show that
#            the others are not significant" (author): every gene inside has P >= 0.10;
#   P = 0.01 at |z| = 2.38 (their 99th percentile) -- Bbc3 alone is beyond it.
# They replace the grey 95% band, whose label read as a panel title. No drawn gene
# falls between the pairs; Apaf1 (P 0.018) would, and it stays out on the author's
# other grounds below.
# SHAPE: circle pro-apoptotic, square anti-apoptotic (Diablo, SMAC, is
# pro-apoptotic: the IAP antagonist; corrected in script 55 on 2026-10-14).
# FILL LIMITS +/-0.5 (author, 2026-10-14: at +/-1.5 "the differences are not
# visible"): the 6W effects have a median |LFC| of 0.26, and only three exceed
# 0.5 (Bid +0.61, Birc5 +1.16, Htra2 +1.49), which are drawn at the end colour.
# The legend names them (printed below).
#
# APAF1 IS LEFT OUT (author, 2026-10-14). It sits outside the band (z -2.14 on
# the raw axis), but MYC had already lowered it at 6W (-0.24 log2), the change
# across the window does not differ by genotype (interaction P 0.39), its
# apoptosome module relaxes rather than tightens under MYC (script 17), and it acts
# after the BCL-2 decision this figure is about. The legend must say it is left out
# and give its value, so the omission is visible (printed below).
#
# Reads: results/departure_from_myc_fade.rds (script 55)
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

dep  <- readRDS(here::here("results", "departure_from_myc_fade.rds"))
BAND <- dep$band$q_abs_z[dep$band$reference == "responsive"]       # the 95% band (legend)
BAND99 <- dep$band$q99_abs_z[dep$band$reference == "responsive"]   # the drawn lines
r <- as.data.frame(dep$rosters)
r <- r[r$roster %in% c("BCL-2 family", "downstream") & is.na(r$dropped), ]
r$fn <- ifelse(r$fn_sign > 0, "anti", "pro")
OMIT <- "Apaf1"
om <- r[r$gene %in% OMIT, ]
stopifnot(nrow(om) == 1L, om$roster == "downstream", om$outside_band, om$int_p > 0.05)
r <- r[!r$gene %in% OMIT, ]

# --- the record is the corrected one, and the claims hold --------------------------------
stopifnot(nrow(r) == 23L, round(BAND, 2) == 1.73,
          r$fn[r$gene == "Diablo"] == "pro",
          setequal(r$gene[r$roster == "downstream" & r$fn == "anti"], c("Xiap", "Birc5")),
          identical(r$gene[r$roster == "BCL-2 family" & r$outside_band], "Bbc3"),
          round(BAND99, 2) == 2.38, identical(r$gene[r$outside_99], "Bbc3"),
          !any(r$outside_band[r$roster == "downstream"]),      # with Apaf1 out, none
          r$gene[r$roster == "BCL-2 family"][which.max(r$z_function[r$roster == "BCL-2 family"])] == "Bbc3")

r$group <- factor(ifelse(r$roster == "BCL-2 family", "BCL­2 family", "Downstream"),
                  levels = c("BCL­2 family", "Downstream"))
r$row <- factor(paste(r$group, r$gene), levels = paste(r$group, r$gene)[order(r$group, r$z_function)])
lab <- stats::setNames(r$gene, paste(r$group, r$gene))
lo <- min(r$z_function, -BAND99); hi <- max(r$z_function, BAND99); pad <- 0.04 * (hi - lo)
# the line labels: along each line, just outside it, half-way down the top facet
n_top <- sum(r$group == levels(r$group)[1])
LINES <- c(BAND, BAND99)
line_lab <- data.frame(group = factor(levels(r$group)[1], levels = levels(r$group)),
                       x = c(-LINES - 0.15, LINES + 0.15), y = (n_top + 1) / 2,
                       label = rep(c("italic(P)=='0.05'", "italic(P)=='0.01'"), 2))
# every drawn gene is beyond both pairs or inside both: the labels need no third class
stopifnot(!any(abs(r$z) > BAND & abs(r$z) <= BAND99), all(r$p_emp_resp[abs(r$z) <= BAND] >= 0.05))
FILL_LIM <- 0.5
beyond <- r[abs(r$lfc6) > FILL_LIM, ]                    # drawn at the end colour (squished)
stopifnot(nrow(beyond) <= 3L)

p <- ggplot(r, aes(z_function, row)) +
  geom_vline(xintercept = c(-LINES, LINES), linetype = "22", linewidth = NM_LINE, colour = "grey45") +
  geom_text(data = line_lab, aes(x, y, label = label), inherit.aes = FALSE, parse = TRUE,
            angle = 90, size = pt(NM_TXT_MIN), family = NM_FONT) +   # top facet only
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
  geom_point(aes(shape = fn, fill = lfc6), size = 1.8, stroke = 0.3, colour = "grey15") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = lab) +
  scale_shape_manual(values = c(pro = 21, anti = 22), breaks = c("pro", "anti"), name = NULL,
                     labels = c(pro = "Pro­apoptotic", anti = "Anti­apoptotic")) +
  scale_fill_div(c(-FILL_LIM, FILL_LIM), name = "MYC+ vs WT, 6W\n(log2 fold change)",
                 breaks = c(-FILL_LIM, 0, FILL_LIM)) +
  scale_x_nm(step = 1, top = hi + pad, bottom = lo - pad, minor = TRUE, labels = lab_signed) +
  labs(x = expression(atop("Difference from the overall MYC effect (" * italic(z) * ")",
                           "" %<-% "toward death" ~~~~~~ "toward survival" %->% "")),
       y = NULL) +
  guides(shape = guide_legend(order = 1, override.aes = list(size = 1.6, fill = "grey80")),
         fill  = guide_colourbar(order = 2, direction = "horizontal", barwidth = unit(16, "mm"),
                                 barheight = unit(2, "mm"), title.position = "top",
                                 ticks.colour = "grey30", frame.colour = "grey30")) +
  theme_nm(legend = "bottom") +
  theme(axis.text.y = element_text(face = "italic"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_text(size = NM_TXT, angle = 0, hjust = 0),
        panel.spacing.y = unit(1.5, "mm"), legend.box = "vertical",
        legend.spacing.y = unit(0.5, "mm"), legend.key.size = unit(2.4, "mm"),
        legend.title = element_text(size = NM_TXT, lineheight = 0.9))

save_panel(p, fig = "EDFig3", panel = "b", name = "fade_ranking", width = 86, height = 100,
           dir = here::here("outputs", "natmetab", "proposals"))

cat(sprintf("\nED 3b -- for the legend: fade rate %.4f; dashed lines at |z| = %.2f and %.2f, the 95th and 99th\n",
            dep$defs$rate, BAND, BAND99))
cat(sprintf("  percentiles over the 8,774 genes MYC moved at 6W (empirical two-sided P = 0.05 and 0.01); Bbc3 beyond\n"))
cat(sprintf("  both (P %.4f); every other drawn gene P >= %.2f\n", r$p_emp_resp[r$gene == "Bbc3"],
            min(r$p_emp_resp[r$gene != "Bbc3"])))
cat(sprintf("  fill limits +/-%.1f; beyond them, at the end colour: %s\n", FILL_LIM,
            paste(sprintf("%s %+.2f", beyond$gene, beyond$lfc6), collapse = ", ")))
cat(sprintf("  NOT DRAWN: %s (%s), 6W %+.2f, 12W %+.2f, z %+.2f (empirical P %.3f: between the P = 0.05\n    and P = 0.01 lines), interaction P %.2f\n",
            om$gene, om$role, om$lfc6, om$lfc12, om$z_function, om$p_emp_resp, om$int_p))
o <- r[order(r$group, -r$z_function), ]
for (i in seq_len(nrow(o))) cat(sprintf("  %-12s %-8s %-4s 6W %+.2f  12W %+.2f  z %+.2f%s\n",
  o$roster[i], o$gene[i], o$fn[i], o$lfc6[i], o$lfc12[i], o$z_function[i],
  if (o$outside_band[i]) "  OUTSIDE" else ""))
