# =============================================================================
# p23_nodox_tints.R -- the construct palette as recorded on 2026-10-23, with the
# no-Dox tints, drawn plain (page 1) and hatched in the full colour (page 2).
# -----------------------------------------------------------------------------
# The author (2026-10-23) chose #00573F for PGC-1a + BCL-xL (p3_pgc1a_bclxl_colours.R,
# candidate (b)) and uses the 28% tint for PGC-1a's no-Dox bar on some bar graphs;
# the same rule gives Ndi1 / LbNOX matrix and LbNOX cytosol theirs. All values are
# read from _style.R (construct_cols, construct_nodox), so this page shows what is
# recorded. Whether the tints stay plain or take a hatch is the open decision.
#
# Each page: the palette, then three graphs under normal vision, deuteranopia,
# protanopia and tritanopia -- the three tints together; as Fig. 3 (EV, BCL-xL,
# PGC-1a no Dox, PGC-1a, PGC-1a + BCL-xL); as Fig. 2k (EV, LbNOX matrix and cytosol,
# each no Dox and Dox) -- and the CIELAB distances of the pairs that matter.
# Prism format as p2_bar_fills.R: graph 20 mm tall, 0.9-pt lines, 5/6-pt type, white
# dots with a black edge, error bar s.e.m. upward. Illustrative values, not data.
# A demo page, not a panel: cairo_pdf (for the Greek alpha), never save_panel().
# Output: outputs/natmetab/proposals/p23_nodox_tints.pdf
#         (options(natmetab.proposal_dir = ...) writes it elsewhere)
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))
out_dir <- getOption("natmetab.proposal_dir", here::here("outputs", "natmetab", "proposals"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

LINE <- lw_pt(0.9)                           # the Prism format's every line
HATCH <- lw_pt(0.5)                          # hatch lines, thinner than the border
GRAPH_H <- 20                                # mm, the graph (data area) height
BAR_W <- 8                                   # mm per bar slot
vision <- list(`Normal vision` = identity,
               Deuteranopia = function(x) colorspace::deutan(x, severity = 1),
               Protanopia   = function(x) colorspace::protan(x, severity = 1),
               Tritanopia   = function(x) colorspace::tritan(x, severity = 1))
de76 <- function(a, b) {                     # CIELAB distance between two colours
  lab <- colorspace::coords(as(colorspace::hex2RGB(c(a, b)), "LAB"))
  sqrt(sum((lab[1, ] - lab[2, ])^2))
}
cc <- construct_cols; nd <- construct_nodox

# ---- hatching: 45-degree lines (on the page) clipped to a rectangle ---------------
hatch_segs <- function(x0, x1, y0, y1, mm_x, mm_y, spacing = 0.9) {
  W <- (x1 - x0) * mm_x; H <- (y1 - y0) * mm_y
  cs <- seq(-W, H, by = spacing * sqrt(2))   # lines y = x + c, in mm from the corner
  xa <- pmax(0, -cs); xb <- pmin(W, H - cs)
  ok <- xa < xb
  data.frame(x = x0 + xa[ok] / mm_x, y = y0 + (xa[ok] + cs[ok]) / mm_y,
             xend = x0 + xb[ok] / mm_x, yend = y0 + (xb[ok] + cs[ok]) / mm_y)
}

# ---- one bar graph: a list of bars, each list(label, fill, hatch, values) ---------
bar_graph <- function(bars, f, title = NULL, ylab = "Cell death (%)") {
  n <- length(bars)
  pts <- do.call(rbind, lapply(seq_len(n), function(i) data.frame(x = i, y = bars[[i]]$values)))
  sm <- data.frame(x = seq_len(n), m = sapply(bars, function(b) mean(b$values)),
                   se = sapply(bars, function(b) stats::sd(b$values) / sqrt(length(b$values))),
                   fill = f(sapply(bars, `[[`, "fill")))
  top <- 35; gw <- BAR_W * n
  mm_x <- gw / (n + 0.2); mm_y <- GRAPH_H / top
  hs <- do.call(rbind, lapply(which(!is.na(sapply(bars, `[[`, "hatch"))), function(i) {
    h <- hatch_segs(i - 0.35, i + 0.35, 0, sm$m[i], mm_x, mm_y)
    if (nrow(h)) h$col <- f(bars[[i]]$hatch)
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
    scale_x_continuous(breaks = seq_len(n), labels = sapply(bars, `[[`, "label"),
                       limits = c(0.4, n + 0.6), expand = expansion(0)) +
    scale_y_continuous(limits = c(0, top), breaks = c(0, 10, 20, 30), expand = expansion(0)) +
    labs(x = NULL, y = ylab, title = title) +
    theme_nm(txt = 5, txt_title = 6, line = LINE) +
    theme(plot.title = element_text(size = 5, face = "plain", hjust = 0, margin = margin(0, 0, 1, 0, "mm")),
          axis.text.x = element_text(size = 5, lineheight = 0.9),
          panel.heights = unit(GRAPH_H, "mm"), panel.widths = unit(gw, "mm"))
}
row_text <- function(txt, size = 6) {
  ggplot() + annotate("text", x = 0, y = 0.5, label = txt, hjust = 0, vjust = 0.5,
                      size = pt(size), family = NM_FONT, lineheight = 0.95) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(0)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(0)) +
    coord_cartesian(clip = "off") + theme_void()
}

# ---- the palette: full colour above, no-Dox tint below (where one is set) ----------
swatch <- function(hatched) {
  keys <- c("EV", "BCL_xL", "PGC1A", "PGC1A_BCL_xL", "Ndi1", "LbNOX_cytosol")
  labs_ <- c("EV", "BCL-xL", "PGC-1α, NRF1", "PGC-1α +\nBCL-xL", "Ndi1, LbNOX\nmatrix", "LbNOX\ncytosol")
  d <- data.frame(x = rep(seq_along(keys) - 1, 2), y = rep(c(0, -1.25), each = length(keys)),
                  col = c(cc[keys], ifelse(keys %in% names(nd), nd[keys], NA)),
                  hatch = c(rep(NA, length(keys)), ifelse(keys %in% names(nd) & hatched, cc[keys], NA)))
  d$lab <- c(paste0(labs_, "\n", cc[keys]),
             ifelse(keys %in% names(nd), paste0("no Dox ", nd[keys]), "no tint set"))
  hs <- do.call(rbind, lapply(which(!is.na(d$hatch)), function(i) {
    h <- hatch_segs(d$x[i] - 0.43, d$x[i] + 0.43, d$y[i] - 0.4, d$y[i] + 0.4, 16, 10)
    h$col <- d$hatch[i]; h
  }))
  sw <- d[!is.na(d$col), ]
  g <- ggplot(sw) +
    geom_tile(aes(x, y, fill = I(col)), width = 0.86, height = 0.8, colour = NA)
  if (!is.null(hs)) g <- g + geom_segment(data = hs, aes(x, y, xend = xend, yend = yend, colour = I(col)),
                                          linewidth = HATCH)
  g +
    geom_tile(aes(x, y), fill = NA, width = 0.86, height = 0.8, colour = "black", linewidth = LINE) +
    geom_text(data = d[d$y == 0, ], aes(x, y + 0.48, label = lab), size = pt(5), family = NM_FONT,
              vjust = 0, lineheight = 0.95) +
    geom_text(data = d[d$y < 0, ], aes(x, y - 0.46, label = lab), size = pt(5), family = NM_FONT,
              vjust = 1, colour = ifelse(is.na(d$col[d$y < 0]), "grey45", "black")) +
    scale_x_continuous(limits = c(-0.5, length(keys) - 0.5), expand = expansion(0)) +
    scale_y_continuous(limits = c(-1.85, 1.15), expand = expansion(0)) +
    coord_cartesian(clip = "off") + theme_void(base_family = NM_FONT)
}

# ---- illustrative values -----------------------------------------------------------
set.seed(7)
draw_n <- function(m, sd, n = 4) pmax(0, stats::rnorm(n, m, sd))
V <- list(EV = draw_n(8, 1.5), BCL = draw_n(5, 1), PGC_n = draw_n(9, 1.5), PGC = draw_n(25, 3),
          PB = draw_n(12, 2), Ndi_n = draw_n(10, 1.5), Ndi = draw_n(28, 3),
          Mat_n = draw_n(10, 1.5), Mat = draw_n(30, 3), Cyt_n = draw_n(9, 1.5), Cyt = draw_n(18, 2.5))
bar <- function(label, key, nodox, values, hatched)
  list(label = label, fill = if (nodox) nd[[key]] else cc[[key]],
       hatch = if (nodox && hatched) cc[[key]] else NA, values = values)
graphs_for <- function(hatched) list(
  `The three no-Dox tints together` = list(
    bar("EV", "EV", FALSE, V$EV, hatched),
    bar("PGC-1α\nno Dox", "PGC1A", TRUE, V$PGC_n, hatched), bar("PGC-1α", "PGC1A", FALSE, V$PGC, hatched),
    bar("Ndi1\nno Dox", "Ndi1", TRUE, V$Ndi_n, hatched), bar("Ndi1", "Ndi1", FALSE, V$Ndi, hatched),
    bar("LbNOX\ncytosol\nno Dox", "LbNOX_cytosol", TRUE, V$Cyt_n, hatched),
    bar("LbNOX\ncytosol", "LbNOX_cytosol", FALSE, V$Cyt, hatched)),
  `As Fig. 3: PGC-1α + BCL-xL` = list(
    bar("EV", "EV", FALSE, V$EV, hatched), bar("BCL-xL", "BCL_xL", FALSE, V$BCL, hatched),
    bar("PGC-1α\nno Dox", "PGC1A", TRUE, V$PGC_n, hatched), bar("PGC-1α", "PGC1A", FALSE, V$PGC, hatched),
    bar("PGC-1α +\nBCL-xL", "PGC1A_BCL_xL", FALSE, V$PB, hatched)),
  `As Fig. 2k: the two LbNOX` = list(
    bar("EV", "EV", FALSE, V$EV, hatched),
    bar("LbNOX\nmatrix\nno Dox", "LbNOX_matrix", TRUE, V$Mat_n, hatched),
    bar("LbNOX\nmatrix", "LbNOX_matrix", FALSE, V$Mat, hatched),
    bar("LbNOX\ncytosol\nno Dox", "LbNOX_cytosol", TRUE, V$Cyt_n, hatched),
    bar("LbNOX\ncytosol", "LbNOX_cytosol", FALSE, V$Cyt, hatched))
)

# ---- the distances that matter -----------------------------------------------------------
pairs <- list(
  c("PGC-1α no Dox", nd[["PGC1A"]], "EV", cc[["EV"]]),
  c("Ndi1 no Dox", nd[["Ndi1"]], "EV", cc[["EV"]]),
  c("LbNOX cytosol no Dox", nd[["LbNOX_cytosol"]], "EV", cc[["EV"]]),
  c("PGC-1α no Dox", nd[["PGC1A"]], "Ndi1 no Dox", nd[["Ndi1"]]),
  c("Ndi1 no Dox", nd[["Ndi1"]], "LbNOX cytosol no Dox", nd[["LbNOX_cytosol"]]),
  c("PGC-1α no Dox", nd[["PGC1A"]], "LbNOX cytosol no Dox", nd[["LbNOX_cytosol"]]),
  c("LbNOX cytosol no Dox", nd[["LbNOX_cytosol"]], "LbNOX cytosol", cc[["LbNOX_cytosol"]]),
  c("LbNOX cytosol no Dox", nd[["LbNOX_cytosol"]], "the white page", "#FFFFFF"),
  c("PGC-1α + BCL-xL", cc[["PGC1A_BCL_xL"]], "PGC-1α", cc[["PGC1A"]]),
  c("PGC-1α + BCL-xL", cc[["PGC1A_BCL_xL"]], "BCL-xL", cc[["BCL_xL"]])
)
de_tab <- do.call(rbind, lapply(pairs, function(p)
  data.frame(pair = paste(p[1], "vs", p[3]),
             t(sapply(vision, function(f) de76(f(p[2]), f(p[4])))), check.names = FALSE)))
de_plot <- function() {
  long <- data.frame(r = rep(seq_len(nrow(de_tab)), length(vision)),
                     c = rep(seq_along(vision), each = nrow(de_tab)),
                     v = unlist(de_tab[names(vision)]))
  ggplot() +
    annotate("text", x = 0, y = 0, label = "CIELAB ΔE (above ~20 reads as different; bold under 15)",
             hjust = 0, size = pt(6), family = NM_FONT) +
    annotate("text", x = 2.6 + seq_along(vision) * 0.75, y = -0.9, label = names(vision),
             hjust = 1, size = pt(5), family = NM_FONT) +
    annotate("text", x = 0, y = -1.7 - (seq_len(nrow(de_tab)) - 1) * 0.8, label = de_tab$pair,
             hjust = 0, size = pt(5), family = NM_FONT) +
    geom_text(data = long, aes(2.6 + c * 0.75, -1.7 - (r - 1) * 0.8, label = sprintf("%.0f", v),
                               fontface = ifelse(v < 15, "bold", "plain")),
              hjust = 1, size = pt(5), family = NM_FONT) +
    scale_x_continuous(limits = c(0, 6), expand = expansion(0)) +
    coord_cartesian(clip = "off") + theme_void()
}

page <- function(hatched) {
  gl <- graphs_for(hatched)
  rows <- lapply(names(vision), function(v) {
    gs <- lapply(names(gl), function(k) bar_graph(gl[[k]], vision[[v]], if (v == names(vision)[1]) k else NULL))
    patchwork::wrap_plots(row_text(v), patchwork::wrap_plots(gs, nrow = 1), ncol = 1,
                          heights = c(4, if (v == names(vision)[1]) 35 else 31))
  })
  patchwork::wrap_plots(c(list(swatch(hatched)), rows, list(de_plot(), patchwork::plot_spacer())),
                        ncol = 1, heights = c(36, 41, 37, 37, 37, 45, 30)) +   # mm, roughly
    patchwork::plot_annotation(
      title = if (hatched) "No-Dox tints HATCHED in the full colour (0.5 pt), as recorded 2026-10-23 (print size)"
              else "No-Dox tints PLAIN, as recorded 2026-10-23 (print size; graphs 20 mm tall)",
      subtitle = paste("Values from _style.R: construct_cols (PGC-1α + BCL-xL #00573F, author 2026-10-23) and",
                       "construct_nodox (the colour at 28% over white).",
                       if (hatched) "\nThe pattern separates the tints from EV and from each other where lightness cannot."
                       else "\nUnder deuteranopia and protanopia only lightness is left, and the tints sit near EV and near each other.",
                       "Illustrative values, not data."),
      theme = theme(plot.title = element_text(size = 8, family = NM_FONT),
                    plot.subtitle = element_text(size = 6, family = NM_FONT)))
}

f <- file.path(out_dir, "p23_nodox_tints.pdf")
grDevices::cairo_pdf(f, width = 210 / 25.4, height = 297 / 25.4, family = NM_FONT, onefile = TRUE)  # demo page only: cairo draws the Greek alpha
print(page(FALSE))
print(page(TRUE))
grDevices::dev.off()

print(cbind(de_tab[1], round(de_tab[-1], 1)), row.names = FALSE)
cat("written:", f, "\n")
