# =============================================================================
# _panel_common.R -- shared spine for the panel-by-panel figure layer
# -----------------------------------------------------------------------------
# figures/panels/ holds ONE SCRIPT PER PANEL, each built to support a specific
# sentence of the Results section as it is written. Assembly into numbered
# manuscript figures happens last, once the narrative fixes the numbering --
# see figures/panels/PANELS.md for the slug -> slot map.
#
# THE PUBLICATION RULE THIS FILE ENFORCES
#
# Nature-family figures carry no explanatory information: the panel shows data,
# axes, units and a minimal key, and everything else belongs in the figure
# legend. theme_panel() therefore BLANKS title, subtitle and caption, so a panel
# cannot accidentally acquire explanatory text. The explanation is written into
# each script as a panel_legend() block, which rebuild_panels.R collects into
# outputs/figures/panels/legends.md -- that file is the raw material for the
# legends, not the figure.
#
# RELATIONSHIP TO figures/theme_myc.R
#
# This file SOURCES theme_myc.R and never edits it, so the four assembled
# manuscript figures (figure1, figure2, figureS1, figureS2) keep rendering
# byte-identically. What it adds is (a) the panel theme, (b) the palettes that
# are currently copy-pasted across a dozen scripts, and (c) a save wrapper that
# honours myc.fig.nosave and writes to outputs/figures/panels/.
#
# CONVENTION every panel script follows
#
#   source(here::here("figures", "panels", "_panel_common.R"))
#   LEGEND <- panel_legend(slot = "Fig. 1B", what = ..., detail = ..., bounds = ...)
#   ... build `p` ...
#   save_panel_p(p, "fig1B_cell_state_composition", height = 70)
#   if (FALSE) { print(p) }        # sandbox, skipped by source()
#
# The composite object is ALWAYS named `p`: both rebuild_panels.R and
# figures/rebuild_manuscript_figures.R locate it by that name to render a dry
# run through a null device.
# =============================================================================

source(here::here("figures", "theme_myc.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("panels need ggplot2")

panel_dir <- here::here("outputs", "figures", "panels")

# --- the panel theme ---------------------------------------------------------
# base_size 7: these are single- or half-column panels that will later be
# composed by patchwork, so type is set for the FINAL size, not the preview.
theme_panel <- function(base_size = 7) {
  theme_myc(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_blank(),   # the publication rule,
      plot.subtitle = ggplot2::element_blank(),   # made mechanical
      plot.caption  = ggplot2::element_blank(),
      legend.title  = ggplot2::element_text(size = base_size),
      legend.text   = ggplot2::element_text(size = base_size),
      legend.key.size = ggplot2::unit(3.2, "mm"),
      plot.margin   = ggplot2::margin(2, 2, 2, 2, "mm"))
}

# --- palettes, declared once -------------------------------------------------
# Every one of these is currently re-declared verbatim in several scripts under
# figures/ (contrast: fig01b:100, figS1b:109, figS2b:73; verdict: fig03:202,
# fig04:67, figure2:76; tier: fig02:40, fig03:43, figure1:60 and five more).
# New panels take them from here. The existing scripts are deliberately left
# untouched.

contrast_cols <- c("Myc@6W"   = "#F4A582", "Myc@12W"    = "#B2182B",
                   "Time WT"  = "#BDBDBD", "Time Myc+"  = "#7B7B7B")

verdict_cols  <- c("withdraws" = "#762A83", "at chance" = "grey55",
                   "rises"     = "#1B7837")

direction_cols <- c(up = "#D6604D", down = "#4393C3")

# MitoPathway Level-1 tiers: Okabe-Ito, with grey72 for the Metabolism catch-all.
tier_cols <- c("OXPHOS"                  = "#E69F00",
               "Mitochondrial dynamics"  = "#56B4E9",
               "Metabolism"              = "grey72",
               "Protein import, sorting" = "#009E73",
               "Signaling"               = "#F0E442",
               "Small molecule transport" = "#0072B2",
               "Central dogma"           = "#D55E00")

# Gene-set quantification methods, as tagged in the library provenance table.
# A single-hue purple ramp, deliberately NOT the genotype blue/red: across the
# figure set a reader learns blue = WT and red = Myc+, and reusing that pair for
# a methods key would spend it on something that is not a genotype. Sequential
# because the routing is ordinal (fGSEA only -> either -> GSVA only).
method_cols <- c("fgsea" = "#54278F", "both" = "#9E9AC8", "gsva" = "#DADAEB")

# --- export ------------------------------------------------------------------
# Wraps theme_myc.R's save_panel() so a panel script never needs to know the
# output directory, and so myc.fig.nosave is honoured in ONE place rather than
# being re-implemented as an `if` in every script.
save_panel_p <- function(plot, slug,
                         width = fig_w[["single"]], height = 70, units = "mm") {
  stopifnot(is.character(slug), length(slug) == 1L, nzchar(slug))
  if (isTRUE(getOption("myc.fig.nosave"))) {
    message("myc.fig.nosave = TRUE -- not writing ", slug, ".pdf")
    return(invisible(NULL))
  }
  if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)
  save_panel(plot, file.path(panel_dir, paste0(slug, ".pdf")),
             width = width, height = height, units = units)
}

# --- the legend block --------------------------------------------------------
# Everything the panel is NOT allowed to say on the page. Four fields, because
# four things have to survive the trip to the legend:
#   what   -- one sentence naming what is plotted (the legend's first line)
#   detail -- n, units, the test, thresholds, what any error bar is
#   bounds -- the caveat that must travel with the claim (batch = timepoint,
#             n=6 power floor, a circularity, whatever applies). Not optional in
#             this project: most of its numbers carry one.
#   source -- the results/*.rds objects and script the numbers came from
panel_legend <- function(slot, what, detail, bounds, source = NULL) {
  stopifnot(is.character(slot),   length(slot)   == 1L,
            is.character(what),   length(what)   == 1L,
            is.character(detail), length(detail) >= 1L,
            is.character(bounds), length(bounds) >= 1L)
  structure(list(slot = slot, what = what, detail = detail,
                 bounds = bounds, source = source),
            class = "panel_legend")
}

# Markdown rendering, used by rebuild_panels.R to build legends.md.
legend_md <- function(x, slug = NULL) {
  stopifnot(inherits(x, "panel_legend"))
  bul <- function(v) paste0("- ", v, collapse = "\n")
  c(paste0("### ", x$slot, if (!is.null(slug)) paste0("  (`", slug, "`)") else ""),
    "", x$what, "",
    "**Detail.**", bul(x$detail), "",
    "**Bounds.**", bul(x$bounds), "",
    if (!is.null(x$source)) c("**Source.**", bul(x$source), "") else NULL)
}

print.panel_legend <- function(x, ...) {
  cat("\n", strrep("-", 74), "\n", sep = "")
  cat(paste(legend_md(x), collapse = "\n"), "\n", sep = "")
  cat(strrep("-", 74), "\n\n", sep = "")
  invisible(x)
}

# --- freshness guard ---------------------------------------------------------
# results/gsva_scores.rds was rebuilt 2026-07-24 with the gene-symbol reconciler
# (docs/2026-07-24_symbol_reconciliation.md). Anything downstream of it that
# predates that rebuild is stale, and a stale panel is worse than a missing one
# because it looks finished. Panels that read a derived GSVA object call this.
require_fresher_than <- function(path, reference = here::here("results", "gsva_scores.rds")) {
  if (!file.exists(path)) {
    stop(basename(path), " is missing -- re-source its producing script (see PANELS.md)",
         call. = FALSE)
  }
  if (file.exists(reference) && file.mtime(path) < file.mtime(reference)) {
    stop(basename(path), " (", format(file.mtime(path), "%Y-%m-%d"), ") is OLDER than ",
         basename(reference), " (", format(file.mtime(reference), "%Y-%m-%d"), ").\n",
         "  It was built before the gene-symbol reconciliation. Re-source, in order:\n",
         "    scripts/17_gsva_overview.R -> scripts/26_dev_program_myc_integration.R\n",
         "    scripts/27_myc_endogenous_amplification.R",
         call. = FALSE)
  }
  invisible(TRUE)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  ## what the legend block looks like
  panel_legend(slot = "Fig. 9Z", what = "A demonstration panel.",
               detail = c("n = 24.", "Bars are SE."),
               bounds = "Nothing is claimed here.") |> print()

  ## dry run of every panel
  options(myc.fig.nosave = TRUE)
  source(here::here("figures", "panels", "rebuild_panels.R"))
}
