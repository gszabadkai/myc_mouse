# =============================================================================
# edfig2d_window_independent.R -- Extended Data Fig. 2d: the window's
# re-weighting is the gland's own; the MYC+ gland follows it
# -----------------------------------------------------------------------------
# The text: "The pubertal-adult transition, by contrast, reshaped the compartment
# independently of MYC (ED Fig. 2d)". Two scatters side by side, one PDF, all 143
# MitoPathways on the priority ruler (mitoPPS):
#
#   LEFT   the WT window against MYC's effect at 6W: unrelated (r near 0), so
#          the window does not undo or extend MYC's re-weighting
#   RIGHT  the MYC+ window against the WT window: the same direction, so the
#          MYC+ gland follows the gland's own programme. Dashed: equal change.
#
# FORM: two signed quantities per pathway -> scatters on equal axes, r printed.
# SIZE (author, 2026-10-14: "make them the same size to align"): every data area in
# ED 2d and 2e is a fixed ED2_SCATTER_MM square on one shared range (ED2_RANGE,
# both scripts), so the three scatters have the same size and the same units, and
# align along the row whatever their axis labels need.
#
# Reads:  results/background_vs_myc.rds (script 40: $ruler)
# Output: outputs/natmetab/EDFig2/EDFig2d_window_independent.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

ruler <- as.data.frame(readRDS(here::here("results", "background_vs_myc.rds"))$ruler)
r143  <- ruler[!ruler$is_mtdna, ]
stopifnot(nrow(r143) == 143L)

# one range for all three scatters of ED 2d-e (the same line is in edfig2e)
ED2_SCATTER_MM <- 44
v   <- c(r143$p_m6, r143$p_tn, r143$p_tp, r143$p_m12)
pad <- 0.04 * diff(range(v)); lo <- min(v) - pad; hi <- max(v) + pad

panel_xy <- function(x, y, xlab, ylab, identity) {
  d  <- data.frame(x = x, y = y)
  rv <- stats::cor(x, y)
  g <- ggplot(d, aes(x, y)) +
    geom_hline(yintercept = 0, linewidth = NM_LINE, colour = "grey80") +
    geom_vline(xintercept = 0, linewidth = NM_LINE, colour = "grey80")
  if (identity) g <- g + geom_abline(slope = 1, intercept = 0, linetype = "22",
                                     linewidth = NM_LINE, colour = "grey45")
  g + geom_point(size = 0.9, colour = "grey30", alpha = 0.7, stroke = 0) +
    annotate("text", x = lo + 0.02, y = hi - 0.03, label = sprintf("italic(r)=='%.2f'", rv),
             parse = TRUE, hjust = 0, size = pt(NM_TXT)) +
    scale_x_nm(step = 0.2, top = hi, bottom = lo, minor = TRUE, labels = lab_signed) +
    scale_y_nm(step = 0.2, top = hi, bottom = lo, minor = TRUE, labels = lab_signed) +
    coord_equal() +
    labs(x = xlab, y = ylab) +
    theme_nm()
}

# the text's "independently of MYC" rests on this correlation being near zero
stopifnot(abs(stats::cor(r143$p_m6, r143$p_tn)) < 0.1)

p_left  <- panel_xy(r143$p_m6, r143$p_tn, "MYC+ vs WT, 6W (mitoPPS)", "12W vs 6W, WT (mitoPPS)",
                    identity = FALSE)
p_right <- panel_xy(r143$p_tn, r143$p_tp, "12W vs 6W, WT (mitoPPS)", "12W vs 6W, MYC+ (mitoPPS)",
                    identity = TRUE)

save_panel(stack_parts(p_left, p_right, ncol = 2, widths = grid::unit(rep(ED2_SCATTER_MM, 2), "mm"),
                       heights = grid::unit(ED2_SCATTER_MM, "mm")),
           fig = "EDFig2", panel = "d", name = "window_independent", width = 114, height = 56)

# --- numbers for the legend (printed, never drawn) --------------------------------
cat("\nED Fig. 2d -- for the legend\n")
cat(sprintf("  %d nuclear MitoPathways, mitoPPS; n = 6 mice per group\n", nrow(r143)))
cat(sprintf("  left:  r(MYC+ vs WT at 6W, WT window) = %.2f\n", stats::cor(r143$p_m6, r143$p_tn)))
cat(sprintf("  right: r(WT window, MYC+ window) = %.2f; dashed: equal change\n",
            stats::cor(r143$p_tn, r143$p_tp)))
cat("  BATCH: each window coincides with an extraction batch\n")

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(.nm_drawable(stack_parts(p_left, p_right, ncol = 2,
                                 widths = grid::unit(rep(ED2_SCATTER_MM, 2), "mm"),
                                 heights = grid::unit(ED2_SCATTER_MM, "mm"))))
}
