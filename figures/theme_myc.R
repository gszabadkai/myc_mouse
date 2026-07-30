# =============================================================================
# theme_myc.R -- shared publication backbone for the MMTV-Myc paper figures
# -----------------------------------------------------------------------------
# Reusable, LIGHT-WEIGHT layer for the figure scripts under figures/. It does
# NOT source 00_setup_packages.R (that pulls DESeq2 etc.): the figure layer only
# reads results/*.rds and renders, so it needs ggplot2 + a couple of helpers.
#
# Provides:
#   geno_cols / geno_labels / group_labels  -- the project genotype palette+labels
#   theme_myc()                              -- theme_classic-based publication theme
#   save_panel()                             -- cairo_pdf export (Type-42 fonts, mm)
#   fig_w                                    -- Nature column widths (mm)
#
# Palette: Okabe-Ito, fixed by the author 2026-07-30 as the project-wide sample
# encoding. Hue = genotype (blue WT / orange-red Myc+), lightness = age
# (6W saturated, 12W light) -- note this INVERTS the earlier 6W-light convention,
# so the four assembled manuscript figures change colour on their next rebuild.
# That is the intent: one sample palette across every figure.
# =============================================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("theme_myc.R needs ggplot2")
}

# --- four-group palette: THE project sample encoding -------------------------
# Author's specification (2026-07-30), Okabe-Ito:
#   6W_wt  #0072B2 dark blue    12W_wt  #56B4E9 sky blue
#   6W_myc #D55E00 vermilion    12W_myc #E69F00 orange
# Keys use the factor levels on disk (wt = neg, myc = pos). Order is genotype-
# major so a legend built with breaks = names(group_cols) reads WT pair then
# Myc+ pair. Use these WHEREVER INDIVIDUAL SAMPLE VALUES ARE SHOWN.
group_cols <- c("6W_neg"  = "#0072B2", "12W_neg" = "#56B4E9",
                "6W_pos"  = "#D55E00", "12W_pos" = "#E69F00")

# --- genotype encoding (identical across all panels) -------------------------
# The two-level palette is the 6W anchor of the four-group one, so a genotype
# key and a group key cannot disagree about which hue means Myc+.
geno_cols   <- c(neg = "#0072B2", pos = "#D55E00")   # WT blue, Myc+ vermilion
geno_labels <- c(neg = "WT", pos = "Myc+")
group_labels <- c("6W_neg" = "6W WT", "12W_neg" = "12W WT",
                  "6W_pos" = "6W Myc+", "12W_pos" = "12W Myc+")

# --- Nature figure widths (mm); the real size control is physical size -------
fig_w <- c(single = 89, onehalf = 120, double = 183)

# --- publication theme -------------------------------------------------------
theme_myc <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid        = ggplot2::element_blank(),
      axis.text         = ggplot2::element_text(colour = "black"),
      axis.ticks        = ggplot2::element_line(colour = "black", linewidth = 0.3),
      axis.line         = ggplot2::element_line(colour = "black", linewidth = 0.3),
      strip.background  = ggplot2::element_blank(),
      strip.text        = ggplot2::element_text(face = "bold", size = base_size,
                                                lineheight = 1.05),
      plot.title        = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle     = ggplot2::element_text(size = base_size, colour = "grey20"),
      plot.caption      = ggplot2::element_text(size = base_size - 2, hjust = 0,
                                                colour = "grey30", lineheight = 1.1),
      plot.tag          = ggplot2::element_text(face = "bold", size = base_size + 4),
      legend.position   = "none")
}

# --- export: cairo_pdf embeds TrueType (editable text), size in mm ------------
save_panel <- function(plot, filename, width = fig_w[["double"]], height = 100,
                       units = "mm") {
  ggplot2::ggsave(filename = filename, plot = plot,
                  width = width, height = height, units = units,
                  device = grDevices::cairo_pdf)
  message("wrote ", filename, "  (", width, "x", height, " ", units, ")")
  invisible(filename)
}
