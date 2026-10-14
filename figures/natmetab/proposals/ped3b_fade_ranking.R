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
# THE BAND: the central 95% of the same quantity over the 8,774 genes MYC moved at
# 6W (|z| < 1.73; 1.71 for the genes it did not move, script 55).
# SHAPE: circle pro-apoptotic, square anti-apoptotic (Diablo, SMAC, is
# pro-apoptotic: the IAP antagonist; corrected in script 55 on 2026-10-14).
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
BAND <- dep$band$q_abs_z[dep$band$reference == "responsive"]
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
          !any(r$outside_band[r$roster == "downstream"]),      # with Apaf1 out, none
          r$gene[r$roster == "BCL-2 family"][which.max(r$z_function[r$roster == "BCL-2 family"])] == "Bbc3")

r$group <- factor(ifelse(r$roster == "BCL-2 family", "BCL­2 family", "Downstream"),
                  levels = c("BCL­2 family", "Downstream"))
r$row <- factor(paste(r$group, r$gene), levels = paste(r$group, r$gene)[order(r$group, r$z_function)])
lab <- stats::setNames(r$gene, paste(r$group, r$gene))
lo <- min(r$z_function, -BAND); hi <- max(r$z_function, BAND); pad <- 0.04 * (hi - lo)
FILL_LIM <- 1.5
stopifnot(max(abs(r$lfc6)) < FILL_LIM)

p <- ggplot(r, aes(z_function, row)) +
  annotate("rect", xmin = -BAND, xmax = BAND, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
  geom_point(aes(shape = fn, fill = lfc6), size = 1.8, stroke = 0.3, colour = "grey15") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = lab) +
  scale_shape_manual(values = c(pro = 21, anti = 22), breaks = c("pro", "anti"), name = NULL,
                     labels = c(pro = "Pro­apoptotic", anti = "Anti­apoptotic")) +
  scale_fill_div(c(-FILL_LIM, FILL_LIM), name = "MYC+ vs WT, 6W\n(log2 fold change)",
                 breaks = c(-1, 0, 1)) +
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

save_panel(p, fig = "EDFig3", panel = "b", name = "fade_ranking", width = 72, height = 96,
           dir = here::here("outputs", "natmetab", "proposals"))

cat(sprintf("\nED 3b -- for the legend: fade rate %.4f; band |z| < %.2f (genes MYC moved at 6W)\n",
            dep$defs$rate, BAND))
cat(sprintf("  NOT DRAWN: %s (%s), 6W %+.2f, 12W %+.2f, z %+.2f (outside the band), interaction P %.2f\n",
            om$gene, om$role, om$lfc6, om$lfc12, om$z_function, om$int_p))
o <- r[order(r$group, -r$z_function), ]
for (i in seq_len(nrow(o))) cat(sprintf("  %-12s %-8s %-4s 6W %+.2f  12W %+.2f  z %+.2f%s\n",
  o$roster[i], o$gene[i], o$fn[i], o$lfc6[i], o$lfc12[i], o$z_function[i],
  if (o$outside_band[i]) "  OUTSIDE" else ""))
