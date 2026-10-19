# =============================================================================
# p2_bar_fills.R -- PROPOSAL (2026-10-19): how strong Fig. 2's bar fills should be.
# -----------------------------------------------------------------------------
# The author (2026-10-19), making the MYAZ PGC-1a bar charts with Option A
# (docs/2026-10-17_fig2_construct_colours.md): "the bar comes almost invisible, the
# fill colour is too light". Option A borrowed Fig. 1's fill rule (the hue at 28% over
# white), which was made for BOXES, where the dots carry the hue. In a bar the fill is
# the only mark, so the rule does not carry over.
#
# The page, at print size, in the Prism format (graph 20 mm tall, 0.9-pt lines, 5/6-pt
# type): three MYAZ EV-vs-PGC-1a graphs and one graph of the five Option A hues, with
# the fill at 28%, 60% and 100%; then 100% under deuteranopia, protanopia and tritanopia.
# Bars are the mean, the error bar s.e.m. (upward only, as Prism draws it on a bar),
# dots white with a black edge, so they read on every fill, black included.
# Illustrative values only -- this page shows fills, not data.
# A demo page, not a panel: cairo_pdf (for the Greek alpha), never save_panel().
# Output: outputs/natmetab/proposals/p2_bar_fills.pdf
#         (options(natmetab.proposal_dir = ...) writes it elsewhere)
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))
out_dir <- getOption("natmetab.proposal_dir", here::here("outputs", "natmetab", "proposals"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

LINE <- lw_pt(0.9)                           # the Prism format's every line
GRAPH_H <- 20                                # mm, the graph (data area) height
alpha_hex <- function(h, a) {                # the hue at a over white, as a solid colour
  v <- grDevices::col2rgb(h)
  w <- round(255 - a * (255 - v))
  grDevices::rgb(w[1, ], w[2, ], w[3, ], maxColorValue = 255)
}
hues <- c(EV = "#CCCCCC", `BCL-xL` = "#000000", `PGC-1α` = "#009E73",
          Ndi1 = "#CC79A7", `LbNOX\ncytosol` = "#F0E442")       # Option A, one per hue

# illustrative values
set.seed(3)
draw_n <- function(m, sd, n = 4) pmax(0, stats::rnorm(n, m, sd))
myaz <- list(
  "OCR (relative)"  = list(EV = draw_n(1, 0.08),  `PGC-1α` = draw_n(1.6, 0.12)),
  "PUMA (relative)" = list(EV = draw_n(1, 0.10),  `PGC-1α` = draw_n(2.2, 0.25)),
  "Cell death (%)"  = list(EV = draw_n(8, 1.5),   `PGC-1α` = draw_n(25, 3))
)
allc <- list(EV = draw_n(8, 1.5), `BCL-xL` = draw_n(5, 1), `PGC-1α` = draw_n(25, 3),
             Ndi1 = draw_n(30, 3), `LbNOX\ncytosol` = draw_n(18, 2.5))

bar_graph <- function(vals, ylab, fills, title = NULL) {
  pts <- do.call(rbind, lapply(names(vals), function(k) data.frame(g = k, y = vals[[k]])))
  pts$g <- factor(pts$g, levels = names(vals))
  sm <- do.call(rbind, lapply(split(pts$y, pts$g), function(y)
    data.frame(m = mean(y), se = stats::sd(y) / sqrt(length(y)))))
  sm$g <- factor(rownames(sm), levels = names(vals))
  top <- max(c(pts$y, sm$m + sm$se)) * 1.08
  ggplot() +
    geom_col(data = sm, aes(g, m, fill = I(unname(fills[as.character(g)]))),
             colour = "black", linewidth = LINE, width = 0.7) +
    geom_errorbar(data = sm, aes(g, ymin = m, ymax = m + se), width = 0.3, linewidth = LINE) +
    geom_point(data = pts, aes(g, y), shape = 21, fill = "white", colour = "black",
               stroke = 0.4, size = 1.1,
               position = position_jitter(width = 0.13, seed = 4)) +
    scale_y_continuous(limits = c(0, top), expand = expansion(0)) +
    labs(x = NULL, y = ylab, title = title) +
    theme_nm(txt = 5, txt_title = 6, line = LINE) +
    theme(plot.title = element_text(size = 6, face = "plain", hjust = 0),
          axis.text.x = element_text(size = 5, lineheight = 0.9),
          panel.heights = unit(GRAPH_H, "mm"))
}

row_of <- function(fills, label) {
  g <- c(lapply(names(myaz), function(k) bar_graph(myaz[[k]], k, fills,
                                                   if (k == names(myaz)[1]) label else NULL)),
         list(bar_graph(allc, "Cell death (%)", fills, " ")))
  g[[1]] <- g[[1]] + theme(panel.widths = unit(12, "mm"))
  g[[2]] <- g[[2]] + theme(panel.widths = unit(12, "mm"))
  g[[3]] <- g[[3]] + theme(panel.widths = unit(12, "mm"))
  g[[4]] <- g[[4]] + theme(panel.widths = unit(36, "mm"))
  patchwork::wrap_plots(g, nrow = 1)
}
vision <- list(Deuteranopia = function(x) colorspace::deutan(x, severity = 1),
               Protanopia   = function(x) colorspace::protan(x, severity = 1),
               Tritanopia   = function(x) colorspace::tritan(x, severity = 1))
fill_at <- function(a, f = identity) stats::setNames(f(alpha_hex(hues, a)), names(hues))

rows <- list(
  row_of(fill_at(0.28), "28%: Fig. 1's box rule (Option A as proposed)"),
  row_of(fill_at(0.60), "60% over white"),
  row_of(fill_at(1.00), "100%: the hue itself"),
  row_of(fill_at(1.00, vision$Deuteranopia), "100%, deuteranopia"),
  row_of(fill_at(1.00, vision$Protanopia),   "100%, protanopia"),
  row_of(fill_at(1.00, vision$Tritanopia),   "100%, tritanopia")
)
page <- patchwork::wrap_plots(rows, ncol = 1) +
  patchwork::plot_annotation(
    title = "Fig. 2 bar fills: Option A's hues at 28%, 60% and 100% (print size; graphs 20 mm tall)",
    subtitle = paste("Bars mean, error bar s.e.m. (upward), dots white with a black edge; 0.9-pt lines;",
                     "5/6-pt type. Illustrative values, not data."),
    theme = theme(plot.title = element_text(size = 8, family = NM_FONT),
                  plot.subtitle = element_text(size = 6, family = NM_FONT)))

f <- file.path(out_dir, "p2_bar_fills.pdf")
grDevices::cairo_pdf(f, width = 210 / 25.4, height = 297 / 25.4, family = NM_FONT)
print(page)
grDevices::dev.off()

# the fills for Prism (hex), at 60% and 100%
tab <- data.frame(construct = sub("\n", " ", names(hues)), hue = hues,
                  fill_60 = alpha_hex(hues, 0.60), fill_28 = alpha_hex(hues, 0.28), row.names = NULL)
print(tab)
cat("written:", f, "\n")
