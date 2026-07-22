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
# Palette: the project's CVD-safe genotype pair (blue WT / red Myc+), already
# validated this project (CVD dE 17.4). Kept identical across every panel.
# =============================================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("theme_myc.R needs ggplot2")
}

# --- genotype encoding (identical across all panels) -------------------------
geno_cols   <- c(neg = "#4575B4", pos = "#D73027")   # WT blue, Myc+ red
geno_labels <- c(neg = "WT", pos = "Myc+")
group_labels <- c("6W_neg" = "6W WT", "6W_pos" = "6W Myc+",
                  "12W_neg" = "12W WT", "12W_pos" = "12W Myc+")

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
