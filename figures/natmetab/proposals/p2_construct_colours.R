# =============================================================================
# p2_construct_colours.R -- PROPOSAL (2026-10-17): Fig. 2's colours for the cell-line
# constructs. Page 1, the palettes in use (Fig. 1, ED 1-3); pages 2-3, Option A
# (recommended: Nature's example palette, the half Fig. 1 does not use, one hue per
# mechanism) and Option B (Paul Tol's muted, one hue per construct), each under normal
# vision, deuteranopia, protanopia, tritanopia and greyscale. The decision is the
# author's; the note is docs/2026-10-17_fig2_construct_colours.md.
# Illustrative values only -- this page shows colours, not data.
# A demo page, not a panel: cairo_pdf (for the Greek alpha), never save_panel().
# Output: outputs/natmetab/proposals/p2_construct_colours.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))
out_dir <- here::here("outputs", "natmetab", "proposals")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tint <- function(h) {                       # the Fig. 1 rule: hue at 28% over white
  v <- grDevices::col2rgb(h)
  w <- round(255 - FILL_ALPHA * (255 - v))
  grDevices::rgb(w[1, ], w[2, ], w[3, ], maxColorValue = 255)
}
vision <- list(
  "Normal vision"      = identity,
  "Deuteranopia"       = function(x) colorspace::deutan(x, severity = 1),
  "Protanopia"         = function(x) colorspace::protan(x, severity = 1),
  "Tritanopia"         = function(x) colorspace::tritan(x, severity = 1),
  "Greyscale"          = function(x) colorspace::desaturate(x)
)
# smallest CIELAB distance (dE76) between any two colours of a set
min_de <- function(cols) {
  lab <- colorspace::coords(as(colorspace::hex2RGB(cols), "LAB"))
  d <- as.matrix(stats::dist(lab)); diag(d) <- Inf
  i <- which(d == min(d), arr.ind = TRUE)[1, ]
  list(de = min(d), pair = names(cols)[i])
}

constructs <- c("EV", "BCL-xL", "PGC-1α", "NRF1", "Ndi1", "LbNOX matrix", "LbNOX cytosol")
opts <- list(
  A = stats::setNames(c("#CCCCCC", "#000000", "#009E73", "#009E73", "#CC79A7", "#CC79A7", "#F0E442"),
                      constructs),
  B = stats::setNames(c("#CCCCCC", "#332288", "#44AA99", "#117733", "#882255", "#AA4499", "#DDCC77"),
                      constructs)
)
opt_title <- c(
  A = "Option A -- Nature's palette, the half Fig. 1 does not use: one hue per mechanism",
  B = "Option B -- Paul Tol's 'muted' palette: one hue per construct"
)
opt_note <- c(
  A = paste("light grey EV; black BCL-xL (anti-apoptotic); bluish green, mitochondrial biogenesis (PGC-1α, NRF1);",
            "reddish purple, matrix NADH oxidation (Ndi1, LbNOX matrix); yellow, cytosolic NADH oxidation."),
  B = paste("Seven distinct hues, none of them Fig. 1's blue or orange. Families kept: teal and green for",
            "biogenesis, wine and purple for matrix NADH oxidation; sand for the cytosol.")
)

# illustrative values: cell death (%), n = 4 per bar
set.seed(1)
mk <- function(con, myc) {
  base <- c(EV = 8, `BCL-xL` = 5, `PGC-1α` = 10, NRF1 = 10, Ndi1 = 9,
            `LbNOX matrix` = 9, `LbNOX cytosol` = 8)[con]
  bump <- c(EV = 12, `BCL-xL` = 1, `PGC-1α` = 32, NRF1 = 28, Ndi1 = 30,
            `LbNOX matrix` = 38, `LbNOX cytosol` = 18)[con]
  data.frame(con = con, myc = myc, y = base + (myc == "on") * bump + stats::rnorm(4, 0, 3))
}
all_d <- do.call(rbind, lapply(constructs, function(k) rbind(mk(k, "off"), mk(k, "on"))))
all_d$con <- factor(all_d$con, levels = constructs)
all_d$myc <- factor(all_d$myc, levels = c("off", "on"), labels = c("−", "+"))
k_d <- droplevels(all_d[all_d$con %in% c("EV", "LbNOX matrix", "LbNOX cytosol"), ])

# boxes take the tint (Fig. 1's rule), dots the hue
panel2 <- function(d, pal, title) {
  d$hue  <- unname(pal[as.character(d$con)])
  d$fill <- tint(d$hue)
  ggplot(d, aes(myc, y)) +
    geom_boxplot(aes(fill = fill), colour = box_line, linewidth = NM_LINE, width = 0.7,
                 outlier.shape = NA, coef = Inf) +
    geom_point(aes(fill = hue), shape = 21, colour = "grey20", stroke = 0.2, size = 0.9,
               position = position_jitter(width = 0.12, seed = 2)) +
    scale_fill_identity() +
    facet_grid(. ~ con, switch = "x", labeller = labeller(con = label_wrap_gen(7))) +
    scale_y_continuous(limits = c(0, 55), breaks = seq(0, 50, 25), expand = expansion(0)) +
    labs(x = NULL, y = "Cell death (%)", title = title) +
    theme_nm() +
    theme(strip.placement = "outside", strip.background = element_blank(),
          strip.text = element_text(size = NM_TXT), panel.spacing.x = unit(1.2, "mm"),
          plot.title = element_text(size = NM_TXT, face = "plain", hjust = 0))
}
swatch <- function(cols, labels, title, ncol = length(cols)) {
  nrow <- ceiling(length(cols) / ncol)
  d <- data.frame(col = cols, lab = labels,
                  x = (seq_along(cols) - 1) %% ncol, y = -((seq_along(cols) - 1) %/% ncol) * 2.4)
  ggplot(d) +
    geom_tile(aes(x, y, fill = I(col)), width = 0.86, height = 0.9, colour = "grey60", linewidth = 0.15) +
    geom_text(aes(x, y - 0.55, label = lab), size = pt(5), family = NM_FONT, vjust = 1, lineheight = 0.95) +
    scale_x_continuous(limits = c(-0.5, ncol - 0.5), expand = expansion(0)) +
    scale_y_continuous(limits = c(-(nrow - 1) * 2.4 - 1.9, 0.5), expand = expansion(0)) +
    coord_cartesian(clip = "off") +
    labs(title = title) + theme_void(base_family = NM_FONT) +
    theme(plot.title = element_text(size = NM_TXT, colour = "black"))
}

# ---- page 1: what is in use --------------------------------------------------
lab_hex <- function(nm, h) paste0(nm, "\n", h)
nat  <- c(Black = "#000000", Orange = "#E69F00", `Sky blue` = "#56B4E9", `Bluish green` = "#009E73",
          Yellow = "#F0E442", Blue = "#0072B2", Vermillion = "#D55E00", `Reddish purple` = "#CC79A7")
nat_use <- c("not in Fig. 1", "MYC+ 12W", "WT 12W", "not in Fig. 1\n(ED 2f-g tier)",
             "not used", "WT 6W", "MYC+ 6W", "not in Fig. 1\n(ED 2b-c, 2f-g)")
p1 <- list(
  swatch(nat, paste0(names(nat), "\n", nat, "\n", nat_use),
         "Nature's example colour-blind palette (research figure guide) = Okabe-Ito; what each colour is in Fig. 1"),
  swatch(c(group_cols, group_fill),
         lab_hex(c(sub("\n", " ", group_labels), paste(sub("\n", " ", group_labels), "fill")),
                 c(group_cols, group_fill)),
         "Fig. 1 and ED: the four groups (hue = genotype, lightness = age); box fills are the hue at 28% over white", ncol = 4),
  swatch(c(div_cols, seq_cols, box_line),
         lab_hex(c("diverging −", "zero", "diverging +", "levels low", "mid", "high", "box lines"),
                 c(div_cols, seq_cols, box_line)),
         "Heatmap scales (espresso-white-mint diverging; slate sequential) and the box line"),
  swatch(c(mito_class_cols, mito_tier_cols),
         lab_hex(c(mito_class_labels, names(mito_tier_cols)), c(mito_class_cols, mito_tier_cols)),
         "ED 2b-c mitochondrial class; ED 2f-g MitoPathway tiers (Okabe-Ito, Metabolism grey)", ncol = 5)
)
page1 <- patchwork::wrap_plots(p1, ncol = 1, heights = c(1, 2, 1, 2)) +
  patchwork::plot_annotation(title = "Colours in use: Fig. 1 and Extended Data Figs 1-3 (figures/natmetab/_style.R)",
                             theme = theme(plot.title = element_text(size = 8, family = NM_FONT)))

# ---- pages 2-3: the options --------------------------------------------------
option_page <- function(o) {
  pal <- opts[[o]]
  rows <- lapply(names(vision), function(v) {
    f <- vision[[v]]
    pv <- stats::setNames(f(pal), names(pal))
    de <- min_de(pv[!duplicated(pv)])
    patchwork::wrap_plots(
      swatch(pv, paste0(sub(" ", "\n", names(pal)), "\n", toupper(pal)), v, ncol = 7),
      panel2(all_d, pv, sprintf("all constructs; closest pair %s / %s, ΔE %.0f", de$pair[1], de$pair[2], de$de)),
      panel2(k_d, pv, "as Fig. 2k (MYAZ)"),
      nrow = 1, widths = c(0.9, 1.6, 0.65))
  })
  patchwork::wrap_plots(rows, ncol = 1) +
    patchwork::plot_annotation(title = opt_title[[o]], subtitle = paste0(opt_note[[o]],
      "\nMYC state (off or basal −, on or raised +) is on the axis, so it takes no colour. Illustrative values, not data.",
      " ΔE: CIELAB distance of the closest two distinct colours (above ~20 reads as different)."),
      theme = theme(plot.title = element_text(size = 8, family = NM_FONT),
                    plot.subtitle = element_text(size = 6, family = NM_FONT)))
}

f <- file.path(out_dir, "p2_construct_colours.pdf")
grDevices::cairo_pdf(f, width = 210 / 25.4, height = 297 / 25.4, family = NM_FONT, onefile = TRUE)  # demo page only: cairo draws the Greek alpha
print(page1)
for (o in names(opts)) print(option_page(o))
grDevices::dev.off()

for (o in names(opts)) for (v in names(vision)) {
  pv <- vision[[v]](opts[[o]]); names(pv) <- names(opts[[o]])
  de <- min_de(pv[!duplicated(pv)])
  cat(sprintf("Option %s  %-14s closest %-14s / %-14s dE %5.1f\n", o, v, de$pair[1], de$pair[2], de$de))
}
cat("written:", f, "\n")
