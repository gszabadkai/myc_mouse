# =============================================================================
# p3_pgc1a_bclxl_colours.R -- PROPOSAL (2026-10-23): a colour for PGC-1a + BCL-xL
# (Figs 3-4), drawn with the whole construct palette and the PGC-1a no-Dox tint.
# -----------------------------------------------------------------------------
# docs/2026-10-17_fig2_construct_colours.md left this open: "the combination needs its
# own colour (a darker green is one candidate)". The author (2026-10-23) asked for the
# darker greens, together with the light green #B8E4D8 they use for PGC-1a in the
# no-Dox condition on some bar graphs (it is Option A's 28% tint of #009E73).
#
# Page 1: the palette together, then five candidates for PGC-1a + BCL-xL, each in a
#   Prism-format bar graph (EV, BCL-xL, PGC-1a no Dox, PGC-1a, the combination) under
#   normal vision, deuteranopia, protanopia and tritanopia. Four are the green
#   darkened (the hue times 0.60 / 0.55 / 0.45) or Paul Tol's green; the fifth is the
#   PGC-1a green with a black hatch, the combination drawn as its two colours.
# Page 2: the no-Dox tint against EV. Under deuteranopia and protanopia only
#   lightness is left, and #B8E4D8 has nearly EV's lightness; no tint of the green
#   escapes, so the page shows the tint in use and the tint with a green hatch.
# Prism format as p2_bar_fills.R: graph 20 mm tall, 0.9-pt lines, 5/6-pt type, white
# dots with a black edge, error bar s.e.m. upward. Illustrative values, not data.
# A demo page, not a panel: cairo_pdf (for the Greek alpha), never save_panel().
# Output: outputs/natmetab/proposals/p3_pgc1a_bclxl_colours.pdf
#         (options(natmetab.proposal_dir = ...) writes it elsewhere)
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))
out_dir <- getOption("natmetab.proposal_dir", here::here("outputs", "natmetab", "proposals"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

LINE <- lw_pt(0.9)                           # the Prism format's every line
HATCH <- lw_pt(0.5)                          # hatch lines, thinner than the border
GRAPH_H <- 20                                # mm, the graph (data area) height
GRAPH_W <- 40                                # mm, five bars
alpha_hex <- function(h, a) {                # the hue at a over white, as a solid colour
  v <- grDevices::col2rgb(h)
  w <- round(255 - a * (255 - v))
  grDevices::rgb(w[1, ], w[2, ], w[3, ], maxColorValue = 255)
}
shade_hex <- function(h, k) {                # the hue times k (towards black)
  v <- grDevices::col2rgb(h)
  w <- round(k * v)
  grDevices::rgb(w[1, ], w[2, ], w[3, ], maxColorValue = 255)
}
vision <- list(`Normal vision` = identity,
               Deuteranopia = function(x) colorspace::deutan(x, severity = 1),
               Protanopia   = function(x) colorspace::protan(x, severity = 1),
               Tritanopia   = function(x) colorspace::tritan(x, severity = 1))
de76 <- function(a, b) {                     # CIELAB distance between two colours
  lab <- colorspace::coords(as(colorspace::hex2RGB(c(a, b)), "LAB"))
  sqrt(sum((lab[1, ] - lab[2, ])^2))
}

GREEN <- unname(construct_cols["PGC1A"])
NODOX <- alpha_hex(GREEN, FILL_ALPHA)        # #B8E4D8, the author's no-Dox PGC-1a
stopifnot(NODOX == "#B8E4D8")

cands <- list(
  a = list(fill = shade_hex(GREEN, 0.60), hatch = NA, name = "green x 0.60"),
  b = list(fill = shade_hex(GREEN, 0.55), hatch = NA, name = "green x 0.55"),
  c = list(fill = shade_hex(GREEN, 0.45), hatch = NA, name = "green x 0.45"),
  d = list(fill = "#117733",              hatch = NA, name = "Paul Tol's green"),
  e = list(fill = GREEN,                  hatch = "#000000", name = "PGC-1α green, black hatch")
)

# ---- illustrative values -------------------------------------------------------
set.seed(5)
draw_n <- function(m, sd, n = 4) pmax(0, stats::rnorm(n, m, sd))
bars <- c("EV", "BCL-xL", "PGC-1α\nno Dox", "PGC-1α", "PGC-1α +\nBCL-xL")
vals <- list(draw_n(8, 1.5), draw_n(5, 1), draw_n(9, 1.5), draw_n(25, 3), draw_n(12, 2))

# ---- hatching: 45-degree lines (on the page) clipped to a rectangle ---------------
hatch_segs <- function(x0, x1, y0, y1, mm_x, mm_y, spacing = 0.9) {
  W <- (x1 - x0) * mm_x; H <- (y1 - y0) * mm_y
  cs <- seq(-W, H, by = spacing * sqrt(2))   # lines y = x + c, in mm from the corner
  xa <- pmax(0, -cs); xb <- pmin(W, H - cs)
  ok <- xa < xb
  data.frame(x = x0 + xa[ok] / mm_x, y = y0 + (xa[ok] + cs[ok]) / mm_y,
             xend = x0 + xb[ok] / mm_x, yend = y0 + (xb[ok] + cs[ok]) / mm_y)
}

# fills, hatches: one per bar (NA hatch = none); f is the vision transform
bar_graph <- function(vals, labels, fills, hatches, f, title = NULL, ylab = "Cell death (%)") {
  n <- length(vals)
  pts <- do.call(rbind, lapply(seq_len(n), function(i) data.frame(x = i, y = vals[[i]])))
  sm <- data.frame(x = seq_len(n), m = sapply(vals, mean),
                   se = sapply(vals, function(y) stats::sd(y) / sqrt(length(y))))
  sm$fill <- f(fills)
  top <- 35
  mm_x <- GRAPH_W / (n + 0.2); mm_y <- GRAPH_H / top
  hs <- do.call(rbind, lapply(which(!is.na(hatches)), function(i) {
    h <- hatch_segs(i - 0.35, i + 0.35, 0, sm$m[i], mm_x, mm_y)
    if (nrow(h)) h$col <- f(hatches[i])
    h
  }))
  g <- ggplot() +
    geom_rect(data = sm, aes(xmin = x - 0.35, xmax = x + 0.35, ymin = 0, ymax = m, fill = I(fill)),
              colour = NA)
  if (!is.null(hs) && nrow(hs))
    g <- g + geom_segment(data = hs, aes(x, y, xend = xend, yend = yend, colour = I(col)),
                          linewidth = HATCH, lineend = "butt")
  g +
    geom_rect(data = sm, aes(xmin = x - 0.35, xmax = x + 0.35, ymin = 0, ymax = m),
              fill = NA, colour = "black", linewidth = LINE) +
    geom_errorbar(data = sm, aes(x, ymin = m, ymax = m + se), width = 0.3, linewidth = LINE) +
    geom_point(data = pts, aes(x, y), shape = 21, fill = "white", colour = "black",
               stroke = 0.4, size = 1.1, position = position_jitter(width = 0.13, height = 0, seed = 4)) +
    scale_x_continuous(breaks = seq_len(n), labels = labels, limits = c(0.4, n + 0.6),
                       expand = expansion(0)) +
    scale_y_continuous(limits = c(0, top), breaks = c(0, 10, 20, 30), expand = expansion(0)) +
    labs(x = NULL, y = ylab, title = title) +
    theme_nm(txt = 5, txt_title = 6, line = LINE) +
    theme(plot.title = element_text(size = 5, face = "plain", hjust = 0, margin = margin(0, 0, 1, 0, "mm")),
          axis.text.x = element_text(size = 5, lineheight = 0.9),
          panel.heights = unit(GRAPH_H, "mm"), panel.widths = unit(GRAPH_W, "mm"))
}
row_text <- function(txt) {
  ggplot() + annotate("text", x = 0, y = 0.5, label = txt, hjust = 0, vjust = 0.5,
                      size = pt(6), family = NM_FONT, lineheight = 0.95) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(0)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(0)) +
    coord_cartesian(clip = "off") + theme_void()
}

# ---- page 1 ------------------------------------------------------------------------
swatch <- function(cols, hatches, labels, title) {
  d <- data.frame(col = cols, lab = labels, x = seq_along(cols) - 1)
  hs <- do.call(rbind, lapply(which(!is.na(hatches)), function(i) {
    h <- hatch_segs(d$x[i] - 0.43, d$x[i] + 0.43, -0.45, 0.45, 13, 9.5)
    h$col <- hatches[i]; h
  }))
  ggplot(d) +
    geom_tile(aes(x, 0, fill = I(col)), width = 0.86, height = 0.9, colour = NA) +
    geom_segment(data = hs, aes(x, y, xend = xend, yend = yend, colour = I(col)), linewidth = HATCH) +
    geom_tile(aes(x, 0), fill = NA, width = 0.86, height = 0.9, colour = "black", linewidth = LINE) +
    geom_text(aes(x, -0.58, label = lab), size = pt(5), family = NM_FONT, vjust = 1, lineheight = 0.95) +
    scale_x_continuous(limits = c(-0.5, length(cols) - 0.5), expand = expansion(0)) +
    scale_y_continuous(limits = c(-1.9, 0.5), expand = expansion(0)) +
    coord_cartesian(clip = "off") +
    labs(title = title) + theme_void(base_family = NM_FONT) +
    theme(plot.title = element_text(size = 6, colour = "black", margin = margin(0, 0, 1.5, 0, "mm")))
}
in_use <- c(EV = "#CCCCCC", `BCL-xL` = "#000000", `PGC-1α,\nNRF1` = GREEN,
            `PGC-1α\nno Dox` = NODOX, `Ndi1, LbNOX\nmatrix` = "#CC79A7", `LbNOX\ncytosol` = "#F0E442")
cand_fill <- sapply(cands, `[[`, "fill"); cand_hatch <- sapply(cands, `[[`, "hatch")
sw <- swatch(c(in_use, cand_fill), c(rep(NA, length(in_use)), cand_hatch),
             c(paste0(names(in_use), "\n", toupper(in_use)),
               paste0("(", names(cands), ")\n", toupper(cand_fill),
                      ifelse(is.na(cand_hatch), "", "\n+ black hatch"))),
             paste("The palette together: Option A in use, the PGC-1α no-Dox tint (in use),",
                   "and (a)-(e), the candidates for PGC-1α + BCL-xL"))

ref <- c(EV = "#CCCCCC", `BCL-xL` = "#000000", `PGC-1α` = GREEN, `no Dox` = NODOX)
cand_row <- function(k) {
  cd <- cands[[k]]
  fills <- c("#CCCCCC", "#000000", NODOX, GREEN, cd$fill)
  hatches <- c(NA, NA, NA, NA, cd$hatch)
  if (is.na(cd$hatch)) {
    worst <- sapply(names(ref), function(r) min(sapply(vision, function(f) de76(f(cd$fill), f(ref[[r]])))))
    note <- sprintf("closest under any of the four: PGC-1α %.0f, BCL-xL %.0f, no Dox %.0f, EV %.0f",
                    worst["PGC-1α"], worst["BCL-xL"], worst["no Dox"], worst["EV"])
  } else note <- "a pattern: apart from PGC-1α and BCL-xL by shape, under any vision"
  graphs <- lapply(names(vision), function(v)
    bar_graph(vals, bars, fills, hatches, vision[[v]], if (k == "a") v else NULL))
  hdr <- row_text(sprintf("(%s) %s, %s.  CIELAB ΔE %s", k, toupper(cd$fill), cd$name, note))
  patchwork::wrap_plots(hdr, patchwork::wrap_plots(graphs, nrow = 1), ncol = 1,
                        heights = c(5, 31))
}
page1 <- patchwork::wrap_plots(c(list(sw), lapply(names(cands), cand_row), list(patchwork::plot_spacer())),
                               ncol = 1, heights = c(24, rep(44, length(cands)), 30)) +   # mm, roughly; the spacer takes the rest
  patchwork::plot_annotation(
    title = "Figs 3-4: a colour for PGC-1α + BCL-xL, with the whole construct palette (print size; graphs 20 mm tall)",
    subtitle = paste("Bars mean, error bar s.e.m. (upward), dots white with a black edge; 0.9-pt lines; 5/6-pt type.",
                     "ΔE: CIELAB distance; above ~20 reads as different.\nRecommended: (b), at least 25 from PGC-1α and 31 from BCL-xL under every simulation.",
                     "Illustrative values, not data."),
    theme = theme(plot.title = element_text(size = 8, family = NM_FONT),
                  plot.subtitle = element_text(size = 6, family = NM_FONT)))

# ---- page 2: the no-Dox tint against EV ------------------------------------------------
nd_bars <- c("EV", "PGC-1α\nno Dox", "PGC-1α")
nd_vals <- list(draw_n(8, 1.5), draw_n(10, 1.5), draw_n(25, 3))
nd_row <- function(hatch, label) {
  graphs <- lapply(names(vision), function(v) {
    f <- vision[[v]]
    bar_graph(nd_vals, nd_bars, c("#CCCCCC", NODOX, GREEN), c(NA, hatch, NA), f,
              sprintf("%s: EV / no Dox ΔE %.0f", v, de76(f("#CCCCCC"), f(NODOX))))
  })
  graphs <- lapply(graphs, function(g) g + theme(panel.widths = unit(24, "mm")))
  patchwork::wrap_plots(row_text(label), patchwork::wrap_plots(graphs, nrow = 1),
                        ncol = 1, heights = c(5, 31))
}
tints <- c(0.28, 0.40, 0.50, 0.60)
tint_tab <- sapply(tints, function(a) {
  t <- alpha_hex(GREEN, a)
  sprintf("%2.0f%% %s: EV %.0f / %.0f, PGC-1α %.0f / %.0f", 100 * a, t,
          de76("#CCCCCC", t), de76(vision$Deuteranopia("#CCCCCC"), vision$Deuteranopia(t)),
          de76(GREEN, t), de76(vision$Deuteranopia(GREEN), vision$Deuteranopia(t)))
})
page2 <- patchwork::wrap_plots(
  nd_row(NA, sprintf("In use: PGC-1α no Dox %s, the green's 28%% tint", NODOX)),
  nd_row(GREEN, "The same tint with a hatch in the PGC-1α green (0.5 pt): the pattern carries what the lightness cannot"),
  row_text(paste(c("Other tints of the green, ΔE normal / deuteranopia: no tint is 20 from both EV and PGC-1α",
                   "under deuteranopia, because EV's lightness lies between the two.", tint_tab), collapse = "\n")),
  patchwork::plot_spacer(),
  ncol = 1, heights = c(44, 44, 26, 160)) +             # mm, roughly; the spacer takes the rest
  patchwork::plot_annotation(
    title = "The PGC-1α no-Dox tint against EV",
    subtitle = paste("Under deuteranopia and protanopia only lightness is left: #B8E4D8 and EV's #CCCCCC differ",
                     "by little more than their lightness.\nIt matters only where EV and PGC-1α no Dox share a graph.",
                     "Illustrative values, not data."),
    theme = theme(plot.title = element_text(size = 8, family = NM_FONT),
                  plot.subtitle = element_text(size = 6, family = NM_FONT)))

f <- file.path(out_dir, "p3_pgc1a_bclxl_colours.pdf")
grDevices::cairo_pdf(f, width = 210 / 25.4, height = 297 / 25.4, family = NM_FONT, onefile = TRUE)  # demo page only: cairo draws the Greek alpha
print(page1)
print(page2)
grDevices::dev.off()

for (k in names(cands)) if (is.na(cands[[k]]$hatch)) for (v in names(vision)) {
  fv <- vision[[v]]
  cat(sprintf("(%s) %s %-14s PGC-1a %5.1f  BCL-xL %5.1f  noDox %5.1f  EV %5.1f\n", k, cands[[k]]$fill, v,
              de76(fv(cands[[k]]$fill), fv(GREEN)), de76(fv(cands[[k]]$fill), fv("#000000")),
              de76(fv(cands[[k]]$fill), fv(NODOX)), de76(fv(cands[[k]]$fill), fv("#CCCCCC"))))
}
cat(tint_tab, sep = "\n")
cat("written:", f, "\n")
