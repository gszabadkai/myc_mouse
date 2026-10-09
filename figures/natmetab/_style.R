# =============================================================================
# _style.R -- the shared style for the Nature Metabolism panels
# -----------------------------------------------------------------------------
# figures/natmetab/ holds ONE SCRIPT PER PANEL. Each panel is drawn at its final
# printed size and assembled by the author in Illustrator, together with blots,
# images and Prism graphs. Panels are placed at 100%: scaling a panel changes the
# size of its text.
#
# The rules are written out in FIGURE_RULES.md. This file implements them, and
# save_panel() enforces the checkable ones: a panel that breaks a rule is refused,
# with the reason, rather than written. figures/natmetab/_style_selftest.R proves
# that each check refuses what it should.
#
# This layer is new (2026-10-09) and shares no code with figures/panels/, the
# August long-form version, which stays untouched as the record.
#
# USE, in a panel script:
#   source(here::here("figures", "natmetab", "_style.R"))
#   ... build `p` with theme_nm(), the palettes, brackets() and scale_y_nm() ...
#   save_panel(p, fig = "Fig1", panel = "b", name = "myc_transcript",
#              width = 48, height = 52)
#   -> outputs/natmetab/Fig1/Fig1b_myc_transcript.pdf
#
#   options(natmetab.dry_run = TRUE)   # run every check, draw to a null device,
#                                      # write nothing
# =============================================================================

for (pkg in c("ggplot2", "ggpubr", "here", "scales")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("figures/natmetab needs the '", pkg, "' package", call. = FALSE)
  }
}
# Attached here on purpose: the figure layer does not source 00_setup_packages.R
# (it only reads results/*.rds and draws), and panel scripts read better with
# ggplot2 on the search path.
suppressPackageStartupMessages(library(ggplot2))

# --- type, lines, page -------------------------------------------------------
NM_FONT      <- "Helvetica"  # throughout: R panels, Prism graphs, Illustrator labels
NM_TXT_MIN   <- 5            # pt; the floor -- nothing on a panel may be smaller
NM_TXT       <- 6            # pt; tick labels, keys, on-panel labels, P values
NM_TXT_TITLE <- 7            # pt; axis titles and strips (Nature's maximum)
NM_LINE      <- 0.25         # ggplot linewidth: axes, ticks, brackets, box outlines
NM_MAX_W     <- 183          # mm; double column
NM_MAX_H     <- 170          # mm; tallest Nature figure that leaves room for a legend
nm_widths    <- c(single = 89, double = 183)

pt <- function(x) x / ggplot2::.pt   # points -> ggplot text size units

# --- the sample palette (Okabe-Ito): hue = genotype, lightness = age ----------
# Keys are the factor levels on disk (neg/pos); labels are what gets printed.
group_cols <- c("6W_neg" = "#0072B2", "12W_neg" = "#56B4E9",
                "6W_pos" = "#D55E00", "12W_pos" = "#E69F00")

# Box and bar fills: the same hue at 28% over white, as SOLID colours, so no
# transparency reaches the PDF and Prism can use identical values.
FILL_ALPHA <- 0.28
group_fill <- c("6W_neg" = "#B8D8E9", "12W_neg" = "#D0EAF9",
                "6W_pos" = "#F3D2B8", "12W_pos" = "#F8E4B8")
local({
  tint <- function(h) {
    v <- grDevices::col2rgb(h)[, 1]
    w <- round(255 - FILL_ALPHA * (255 - v))
    grDevices::rgb(w[1], w[2], w[3], maxColorValue = 255)
  }
  stopifnot(identical(vapply(group_cols, tint, character(1)), group_fill))
})

group_labels <- c("6W_neg" = "WT\n6W",   "12W_neg" = "WT\n12W",
                  "6W_pos" = "MYC+\n6W", "12W_pos" = "MYC+\n12W")
geno_cols    <- c(neg = "#0072B2", pos = "#D55E00")   # the 6W pair of group_cols
geno_labels  <- c(neg = "WT", pos = "MYC+")
box_line     <- "#595959"                             # box outlines and medians

# Effects (contrasts) are labelled by what is compared, then where. The genotype
# effects take the MYC+ hues of their age, because that is what they measure.
# Development contrasts get their labels and colours when the first panel needs them.
contrast_cols   <- c(myc_6W = "#D55E00", myc_12W = "#E69F00")
contrast_labels <- c(myc_6W = "MYC+ vs WT, 6W", myc_12W = "MYC+ vs WT, 12W")
stopifnot(identical(unname(contrast_cols), unname(group_cols[c("6W_pos", "12W_pos")])))

# The order follows the test: age-major when the comparisons are within-age gaps,
# genotype-major when the comparison is a genotype main effect.
order_by_age  <- c("6W_neg", "6W_pos", "12W_neg", "12W_pos")
order_by_geno <- c("6W_neg", "12W_neg", "6W_pos", "12W_pos")
stopifnot(setequal(names(group_labels), names(group_cols)),
          setequal(names(group_fill), names(group_cols)))

# --- diverging and sequential fills (author's specification, 2026-07-30) ------
# Espresso -> white -> mint, white pinned at zero. Not blue-red: blue and orange
# already mean genotype. Each arm is scaled to its own end of `limits`, so on a
# lopsided range equal ink is not equal magnitude across zero -- use the
# asymmetric form only where an axis carries the magnitude, and show the bar.
div_cols <- c(neg = "#4A3525", zero = "#FAFAFA", pos = "#2A8A6D")
seq_cols <- c(low = "#F2F4F6", mid = "#8FA3B0", high = "#243642")  # unsigned levels

scale_fill_div <- function(limits, name = NULL, ...) {
  if (length(limits) == 1L) limits <- c(-abs(limits), abs(limits))
  stopifnot(length(limits) == 2L, limits[1] < 0, limits[2] > 0)
  ggplot2::scale_fill_gradientn(
    colours = unname(div_cols[c("neg", "zero", "pos")]),
    values = scales::rescale(c(limits[1], 0, limits[2])),
    limits = limits, oob = scales::squish, name = name, labels = lab_signed, ...)
}

scale_fill_seq <- function(limits, name = NULL, ...) {
  stopifnot(length(limits) == 2L, limits[1] < limits[2])
  ggplot2::scale_fill_gradientn(colours = unname(seq_cols), limits = limits,
                                oob = scales::squish, name = name, ...)
}

# --- numbers on axes -----------------------------------------------------------
# Thousands take commas (Nature: 1,000), zero is "0", no trailing zeros.
lab_num <- function(x) {
  out <- format(x, big.mark = ",", scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
  out[is.na(x)] <- NA_character_
  out
}
# Signed axes print "+" on the positive side. The ASCII hyphen prints as a true
# minus sign: R's PDF encoding maps it to the minus glyph.
lab_signed <- function(x) {
  out <- lab_num(x)
  pos <- !is.na(x) & x > 0
  out[pos] <- paste0("+", out[pos])
  out
}

# --- P values: formatted here, never computed here -----------------------------
# Exact, two significant figures, italic capital P (the journal's house style),
# x 10^k below 0.001. The exponent is a quoted string inside displaystyle():
# plotmath otherwise prints superscripts at ~70% size (4.2 pt inside a 6 pt label)
# and puts a gap after an unquoted minus. Plotmath keeps the source ASCII.
fmt_p <- function(p) {
  stopifnot(is.numeric(p), !anyNA(p), all(p >= 0 & p <= 1))
  vapply(p, function(v) {
    if (v >= 0.001) {
      return(sprintf("italic(P)=='%s'", formatC(signif(v, 2), format = "fg", digits = 2, flag = "#")))
    }
    if (v < 1e-300) return("italic(P)<10^displaystyle('-300')")
    e <- floor(log10(v))
    m <- signif(v / 10^e, 2)
    if (m >= 10) { m <- m / 10; e <- e + 1 }
    sprintf("italic(P)=='%s'%%*%%10^displaystyle('%d')", format(m), e)
  }, character(1))
}

# --- comparison brackets --------------------------------------------------------
# P values arrive computed by the analysis of record (DESeq2, the models of
# record) and are only drawn here. There is deliberately no way to run a test:
# ggpubr::stat_compare_means() and ggstatsplot test the plotted points and would
# print values that disagree with the text.
#
#   comps  data.frame with group1, group2 (x values as drawn) and p; optional
#          `level` stacks brackets (1 = lowest; default: row order)
#   base   y of the lowest bracket, on the drawn scale (just above the data)
#   step   vertical distance between stacked brackets; also the room left above
#          the top label
#
# Returns the layer, with attr "top": the y the axis must reach so the top label
# is not clipped. Pass it on: scale_y_nm(step = ..., top = attr(br, "top")).
brackets <- function(comps, base, step) {
  stopifnot(is.data.frame(comps), all(c("group1", "group2", "p") %in% names(comps)),
            nrow(comps) >= 1L, is.numeric(base), length(base) == 1L,
            is.numeric(step), length(step) == 1L, step > 0)
  if (is.null(comps$level)) comps$level <- seq_len(nrow(comps))
  comps$y.position <- base + (comps$level - 1) * step
  comps$label <- fmt_p(comps$p)
  layer <- ggpubr::geom_bracket(
    data = comps, inherit.aes = FALSE,
    mapping = ggplot2::aes(xmin = group1, xmax = group2,
                           y.position = y.position, label = label),
    type = "expression", label.size = pt(NM_TXT), linewidth = NM_LINE,
    tip.length = 0.012, vjust = -0.15, colour = "black", family = NM_FONT)
  structure(list(layer), top = max(comps$y.position) + 0.9 * step)
}

# --- axes that end on ticks -----------------------------------------------------
# Both ends sit on a tick and there is no padding. `bottom` is 0 for magnitudes;
# a signed axis passes its own (negative) bottom and labels = lab_signed.
# minor = TRUE draws minor ticks halfway between the labelled ones and lets an end
# fall on one of them, which trims the empty band rounding up to a major tick can
# leave (e.g. 18,000 rather than 20,000 on a 4,000 step).
.nm_axis <- function(step, top, bottom, minor) {
  stopifnot(is.numeric(step), length(step) == 1L, step > 0, top > bottom)
  unit <- if (minor) step / 2 else step
  lo <- floor(bottom / unit + 1e-9) * unit
  hi <- ceiling(top / unit - 1e-9) * unit
  list(limits = c(lo, hi),
       breaks = seq(ceiling(lo / step - 1e-9) * step, hi + 1e-9 * step, by = step),
       minor_breaks = if (minor) seq(lo, hi + 1e-9 * unit, by = unit) else NULL,
       guide = ggplot2::guide_axis(minor.ticks = minor))
}

scale_y_nm <- function(step, top, bottom = 0, minor = FALSE, labels = lab_num, ...) {
  a <- .nm_axis(step, top, bottom, minor)
  ggplot2::scale_y_continuous(limits = a$limits, breaks = a$breaks, minor_breaks = a$minor_breaks,
                              guide = a$guide, expand = c(0, 0), labels = labels, ...)
}

scale_x_nm <- function(step, top, bottom = 0, minor = FALSE, labels = lab_num, ...) {
  a <- .nm_axis(step, top, bottom, minor)
  ggplot2::scale_x_continuous(limits = a$limits, breaks = a$breaks, minor_breaks = a$minor_breaks,
                              guide = a$guide, expand = c(0, 0), labels = labels, ...)
}

# --- the theme -------------------------------------------------------------------
# Titles, subtitles and captions are blanked: on a Nature panel the explanation
# belongs in the legend, so a panel cannot acquire explanatory text by accident.
theme_nm <- function(legend = "none") {
  ggplot2::theme_classic(base_size = NM_TXT_TITLE, base_family = NM_FONT) +
    ggplot2::theme(
      text              = ggplot2::element_text(colour = "black"),
      axis.text         = ggplot2::element_text(size = NM_TXT, colour = "black"),
      axis.title        = ggplot2::element_text(size = NM_TXT_TITLE, colour = "black"),
      axis.line         = ggplot2::element_line(linewidth = NM_LINE, colour = "black"),
      axis.ticks        = ggplot2::element_line(linewidth = NM_LINE, colour = "black"),
      axis.ticks.length = ggplot2::unit(1, "mm"),
      axis.minor.ticks.length = ggplot2::rel(0.5),
      legend.text       = ggplot2::element_text(size = NM_TXT),
      legend.title      = ggplot2::element_text(size = NM_TXT),
      legend.key.size   = ggplot2::unit(3, "mm"),
      legend.position   = legend,
      strip.text        = ggplot2::element_text(size = NM_TXT_TITLE),
      strip.background  = ggplot2::element_blank(),
      plot.title        = ggplot2::element_blank(),
      plot.subtitle     = ggplot2::element_blank(),
      plot.caption      = ggplot2::element_blank(),
      plot.margin       = ggplot2::margin(1, 2.5, 1, 1, "mm"))  # right: end labels overhang the axis
}

# --- the checks ------------------------------------------------------------------
# Run by save_panel() on every panel. Each rule names its source in
# FIGURE_RULES.md. All problems are collected and reported together.
.nm_text_geoms <- c("GeomText", "GeomLabel", "GeomBracket", "GeomTextRepel", "GeomLabelRepel")

# Are minor ticks DRAWN on this axis? ggplot computes minor breaks whether or not
# they are drawn, so an axis ending on an undrawn one would look tickless. Drawn
# means minor.ticks = TRUE on the scale's guide or on guides(<axis> = ...).
.nm_minor_drawn <- function(p, v, ax) {
  on_scale  <- inherits(v$scale$guide, "Guide") && isTRUE(v$scale$guide$params$minor.ticks)
  g <- tryCatch(p$guides$guides[[ax]], error = function(e) NULL)
  on_guides <- inherits(g, "Guide") && isTRUE(g$params$minor.ticks)
  on_scale || on_guides
}

.nm_achromatic <- function(cols) {
  cols <- unique(stats::na.omit(as.character(cols)))
  if (!length(cols)) return(TRUE)
  rgb <- grDevices::col2rgb(cols)
  all(apply(rgb, 2, function(v) diff(range(v)) == 0))
}

check_panel <- function(p, width, height) {
  errs <- character(0)
  add <- function(...) errs <<- c(errs, sprintf(...))

  # page size
  if (width > NM_MAX_W)  add("width %g mm is over the %g mm page", width, NM_MAX_W)
  if (height > NM_MAX_H) add("height %g mm is over the %g mm page", height, NM_MAX_H)

  # theme text: size, colour, font
  thm <- ggplot2::complete_theme(p$theme)
  for (e in c("axis.text.x.bottom", "axis.text.y.left", "axis.title.x.bottom",
              "axis.title.y.left", "legend.text", "legend.title",
              "strip.text.x.top", "strip.text.y.right", "plot.tag")) {
    el <- ggplot2::calc_element(e, thm)
    if (!inherits(el, "element_text")) next
    if (!is.null(el$size) && el$size < NM_TXT_MIN) add("theme %s is %.1f pt", e, el$size)
    if (!is.null(el$colour) && !.nm_achromatic(el$colour)) add("theme %s is coloured", e)
    if (!is.null(el$family) && !el$family %in% c("", NM_FONT)) add("theme %s uses '%s'", e, el$family)
  }

  # text layers: size, colour, font, superscripts
  b <- ggplot2::ggplot_build(p)
  for (i in seq_along(p$layers)) {
    g <- class(p$layers[[i]]$geom)[1]
    if (!g %in% .nm_text_geoms) next
    d <- b$data[[i]]
    sz <- if (g == "GeomBracket") d$label.size else d$size
    if (length(sz) && any(sz * ggplot2::.pt < NM_TXT_MIN - 1e-9)) {
      add("layer %d (%s): text at %.2f pt", i, g, min(sz) * ggplot2::.pt)
    }
    if (!.nm_achromatic(d$colour)) add("layer %d (%s): coloured text", i, g)
    if (!is.null(d$family) && any(!d$family %in% c("", NM_FONT))) {
      add("layer %d (%s): font is not %s", i, g, NM_FONT)
    }
    if (!is.null(d$label) && any(grepl("\\^(?!displaystyle\\()", as.character(d$label), perl = TRUE))) {
      add("layer %d (%s): a superscript would print at ~70%% size", i, g)
    }
  }

  # axes: both ends on a tick; y from 0 (or spanning 0) whenever an x axis is drawn
  x_drawn <- !all(vapply(c("axis.text.x.bottom", "axis.line.x.bottom"), function(e)
    inherits(ggplot2::calc_element(e, thm), "element_blank"), logical(1)))
  yv <- unlist(lapply(b$data, function(d)
    unlist(d[intersect(names(d), c("y", "ymin", "ymax", "lower", "upper", "middle"))],
           use.names = FALSE)), use.names = FALSE)
  for (k in seq_along(b$layout$panel_params)) {
    pp <- b$layout$panel_params[[k]]
    for (ax in c("x", "y")) {
      v <- pp[[ax]]
      if (is.null(v) || v$is_discrete()) next
      r <- v$continuous_range
      ticks <- stats::na.omit(c(v$get_breaks(),
                                if (.nm_minor_drawn(p, v, ax)) v$get_breaks_minor()))
      tol <- 1e-6 * max(1, abs(diff(r)))
      on_tick <- function(z) length(ticks) > 0 && any(abs(ticks - z) <= tol)
      if (!on_tick(r[1]) || !on_tick(r[2])) {
        add("panel %d: the %s axis does not end on a tick (it runs %s to %s)",
            k, ax, format(signif(r[1], 4)), format(signif(r[2], 4)))
      }
      if (ax == "y" && x_drawn) {
        tr <- v$scale$get_transformation()$name
        if (tr != "identity") {
          add("panel %d: the y axis is %s, which has no 0; with an x axis drawn it must start at 0", k, tr)
        } else if (r[1] > tol) {
          add("panel %d: the y axis starts at %s, not 0", k, format(signif(r[1], 4)))
        } else if (r[2] < -tol) {
          add("panel %d: the y axis does not reach 0", k)
        } else if (r[1] < -tol && length(yv) && all(yv >= 0, na.rm = TRUE)) {
          add("panel %d: the data are non-negative, so the y axis must start at 0, not %s",
              k, format(signif(r[1], 4)))
        }
      }
    }
  }

  if (length(errs)) stop("panel refused:\n  - ", paste(errs, collapse = "\n  - "), call. = FALSE)
  invisible(TRUE)
}

# --- export ----------------------------------------------------------------------
# Base pdf() with Helvetica: one of the standard PDF fonts, referenced rather than
# embedded. Illustrator supplies the Mac's Helvetica, keeps the text editable and
# embeds the font when the assembled figure is saved. NOT cairo_pdf: it embeds the
# font but places every letter separately, and at 5-7 pt Apple's viewers collapse
# the spaces between words (commit a58fa1c; re-checked 2026-10-09).
#
#   fig    "Fig1".."Fig4" or "EDFig1".."EDFig10"
#   panel  the lowercase letter; two R pieces of one panel share it and differ by name
#   name   lowercase, digits and underscores
#   dir    override the output folder (default outputs/natmetab/<fig>/)
save_panel <- function(p, fig, panel, name, width, height, dir = NULL) {
  stopifnot(grepl("^(Fig[1-4]|EDFig([1-9]|10))$", fig), grepl("^[a-z]$", panel),
            grepl("^[a-z0-9_]+$", name), is.numeric(width), is.numeric(height))
  check_panel(p, width, height)
  if (is.null(dir)) dir <- here::here("outputs", "natmetab", fig)
  file <- file.path(dir, sprintf("%s%s_%s.pdf", fig, panel, name))
  if (isTRUE(getOption("natmetab.dry_run"))) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
    message("dry run: checked and drew ", basename(file), " (", width, " x ", height,
            " mm); nothing written")
    return(invisible(file))
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file, p, width = width, height = height, units = "mm",
                  device = function(filename, width, height, ...)
                    grDevices::pdf(filename, width = width, height = height,
                                   family = NM_FONT, useDingbats = FALSE))
  message("wrote ", file, " (", width, " x ", height, " mm)")
  invisible(file)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  ## what a P label looks like before it is drawn
  fmt_p(c(0.48, 0.0074, 2.65e-17))

  ## the checks, on their own
  source(here::here("figures", "natmetab", "_style_selftest.R"))

  ## a dry run of one panel: every check, drawn to a null device, nothing written
  options(natmetab.dry_run = TRUE)
  # source(here::here("figures", "natmetab", "<panel script>.R"))
  options(natmetab.dry_run = NULL)
}
