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
# Drawn labels are the author's own names for the four groups (2026-07-31), not a
# prettified version of them: the Results text, the analysis objects and the
# figures now all say 6W_wt / 12W_wt / 6W_myc / 12W_myc. The keys stay the
# on-disk factor levels (neg/pos), because that is what the data carry.
group_labels <- c("6W_neg" = "6W_wt", "12W_neg" = "12W_wt",
                  "6W_pos" = "6W_myc", "12W_pos" = "12W_myc")

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

# --- export: base pdf(), Helvetica, size in mm -------------------------------
# WHY NOT cairo_pdf (changed 2026-07-31, author reported uneven character
# spacing). cairo writes text as individually positioned glyphs, and at the 6 pt
# type these panels use the inter-word advance rounds away: on this system
# "Hallmark comparator" renders as "Hallmarkcomparator" and "nucleotide" as
# "nudeotide" in Preview and in every Core Graphics viewer, at every font family
# tried. The PDF's text content is correct -- pdftotext reads it back with the
# spaces -- so it is a rendering artefact, but it is the artefact the reader
# sees. The base pdf() device writes real text strings with font metrics instead
# and is clean at the same size.
#
# THE TRADE. Base pdf() uses Helvetica, one of the 14 standard PDF fonts, which
# is referenced rather than embedded (there is no Ghostscript on this machine, so
# embedFonts() is not available). That is fine for working figures and for
# journals that accept the standard 14; if a submission demands full embedding,
# set options(myc.fig.cairo = TRUE) to switch back and check the spacing at final
# size, or embed once in Illustrator on the way out. Every string in this project
# is ASCII by the coding rules, so nothing needs cairo's UTF-8 handling.
save_panel <- function(plot, filename, width = fig_w[["double"]], height = 100,
                       units = "mm") {
  dev <- if (isTRUE(getOption("myc.fig.cairo"))) grDevices::cairo_pdf else
    function(filename, width, height, ...)
      grDevices::pdf(file = filename, width = width, height = height,
                     family = "Helvetica", useDingbats = FALSE)
  ggplot2::ggsave(filename = filename, plot = plot,
                  width = width, height = height, units = units, device = dev)
  message("wrote ", filename, "  (", width, "x", height, " ", units, ")")
  invisible(filename)
}
