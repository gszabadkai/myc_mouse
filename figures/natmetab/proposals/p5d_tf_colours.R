# =============================================================================
# p5d_tf_colours.R -- PROPOSAL: three colour schemes for ED Fig. 5d, the TF
#                     activity charts, at printed size and on the ED 5 page
# -----------------------------------------------------------------------------
# The author (2026-11-07): "As for the colours, I need to see which is most
# fitting to the general scheme." Three options, each as the two charts
# (Bcl-xL + Pgc-1a against Bcl-xL + EV, and against EV + EV):
#   A  every bar in the Bcl-xL + Pgc-1a green (#00573F, construct_cols), filled
#      where BH-adjusted P < 0.05 and white otherwise -- Fig. 1e's filled/open rule,
#      colour naming the arm whose effect is drawn ("genotype effects take the
#      MYC+ hue").
#   B  each bar in the colour of the arm where the TF is more active: the green,
#      or the comparator's construct colour -- black (Bcl-xL + EV) or light grey
#      #CCCCCC (EV + EV, FIGURE_RULES; the assembled ED 5a draws EV + EV
#      near-white).
#   C  the paper's diverging fill (espresso - white - mint, div_cols) by score:
#      the literal swap for the original red-blue gradient. It encodes what the
#      bar's length already shows.
#
# The format is ED 5's Prism graphs': lines 0.72 pt (the assembled ED 5 measures
# 0.715), tick labels and TF names 5 pt, strip and axis titles 6 pt. TF names are
# roman: they name an inferred protein activity, not a transcript. The 25 TFs per
# contrast are script 58's panel table (the 25 largest |score|), the same x scale
# on both charts and all three options.
#
# Reads : results/orthotopic_tf_activity.rds (script 58), or the file named by
#         options(natmetab.s58_rds = ...) -- e.g. a dry run's
#         options(natmetab.ed5_png = ...): the assembled ED 5 (Drive) rendered to
#         an A4 PNG, for the page mock-ups; without it those pages are skipped
# Writes: outputs/natmetab/proposals/p5d/EDFig5d_tf_<a|b|c>_<vs_bclxl|vs_ev>.pdf
#         outputs/natmetab/proposals/p5d_tf_colours.pdf (the sheet and the pages)
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))
suppressPackageStartupMessages(library(grid))

S58  <- getOption("natmetab.s58_rds", here::here("results", "orthotopic_tf_activity.rds"))
ED5  <- getOption("natmetab.ed5_png", NA_character_)
OUTD <- here::here("outputs", "natmetab", "proposals")
o  <- readRDS(S58)
pn <- as.data.frame(o$panel)
stopifnot(all(table(pn$contrast) == 25L), setequal(unique(pn$contrast), c("vs_bclxl", "vs_ev")))

TXT    <- 5
TITLE  <- 6
LINE   <- lw_pt(NM_LINE_PRISM[["Fig3"]])   # 0.72 pt
ROW    <- 2.0                              # mm per TF
DATA_W <- 30                               # mm, the data area's width
SH     <- "\u00ad"                         # soft hyphen: base pdf() prints "-" as a minus
GREEN  <- construct_cols[["PGC1A_BCL_xL"]]
COMPARATOR <- c(vs_bclxl = construct_cols[["BCL_xL"]], vs_ev = construct_cols[["EV"]])
STRIP <- c(vs_bclxl = sprintf('paste("Bcl%sxL + Pgc%s1", alpha, " vs Bcl%sxL + EV")', SH, SH, SH),
           vs_ev    = sprintf('paste("Bcl%sxL + Pgc%s1", alpha, " vs EV + EV")', SH, SH))

lim <- range(c(0, pn$score))
pad <- 0.04 * diff(lim)                    # FIGURE_RULES: pad 3-4% before rounding to a tick
ax  <- .nm_axis(2, lim[2] + pad, lim[1] - pad, TRUE)

tf_chart <- function(k, opt) {
  d <- pn[pn$contrast == k, ]
  d$tf <- factor(d$source, levels = d$source[order(d$score)])
  d$strip <- STRIP[[k]]
  p <- ggplot(d, aes(x = score, y = tf))
  if (opt == "C") {
    p <- p + geom_col(aes(fill = score), orientation = "y", width = 0.72,
                      colour = "black", linewidth = LINE) +
      scale_fill_div(limits = ax$limits)
  } else {
    d$fill <- if (opt == "A") ifelse(d$bh_sig, GREEN, "white")
              else ifelse(d$score > 0, GREEN, COMPARATOR[[k]])
    p <- ggplot(d, aes(x = score, y = tf)) +
      geom_col(aes(fill = fill), orientation = "y", width = 0.72,
               colour = "black", linewidth = LINE) +
      scale_fill_identity()
  }
  p + geom_vline(xintercept = 0, linewidth = LINE, colour = "black") +
    facet_wrap(~ strip, labeller = label_parsed) +
    scale_x_nm(step = 2, top = lim[2] + pad, bottom = lim[1] - pad, minor = TRUE,
               labels = lab_signed) +
    labs(x = "TF activity (t)", y = NULL) +
    theme_nm(txt = TXT, txt_title = TITLE, line = LINE) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          strip.text = element_text(size = TITLE, colour = "black"),
          plot.margin = margin(0.5, 2.5, 0.5, 0.5, "mm"))
}

# --- each chart through save_panel(), so every rule is checked -------------------
op <- options(natmetab.capture = TRUE)
charts <- list(); sizes <- list()
for (opt in c("A", "B", "C")) for (k in c("vs_bclxl", "vs_ev")) {
  st <- fixed_size(tf_chart(k, opt), DATA_W, 25 * ROW)
  sz <- ceiling(natural_size(st) * 2) / 2
  stopifnot(!anyNA(sz))
  save_panel(st, fig = "EDFig5", panel = "d", name = sprintf("tf_%s_%s", tolower(opt), k),
             width = sz[["width"]], height = sz[["height"]], dir = file.path(OUTD, "p5d"))
  charts[[paste(opt, k)]] <- get(".nm_captured", envir = globalenv())
  sizes[[paste(opt, k)]]  <- sz
}
options(op)

# --- the sheet: the three options at printed size, then each on the ED 5 page ----
mm <- function(v) unit(v, "mm")
hy <- function(x) gsub("(?<=[A-Za-z0-9])-(?=[A-Za-z0-9])", SH, x, perl = TRUE)
PW <- 210; PH <- 297
place <- function(p, x, y, sz) {
  print(p, vp = viewport(x = mm(x), y = mm(PH - y), width = mm(sz[["width"]]),
                         height = mm(sz[["height"]]), just = c("left", "top")), newpage = FALSE)
}
caption <- c(
  A = "A  The Bcl-xL + Pgc-1a green, filled where BH-adjusted P < 0.05 and white otherwise (Fig. 1e's filled/open rule).",
  B = "B  The colour of the arm where the TF is more active: green against black (Bcl-xL + EV) or light grey (EV + EV).",
  C = "C  The paper's diverging fill by score (espresso to white to mint), the direct swap for the red-blue gradient.")
txt <- function(label, x, y, size = 7, face = "plain") {
  grid.text(hy(label), x = mm(x), y = mm(PH - y), just = c("left", "top"),
            gp = gpar(fontsize = size, fontface = face, fontfamily = NM_FONT))
}

pdf(file.path(OUTD, "p5d_tf_colours.pdf"), width = PW / 25.4, height = PH / 25.4,
    family = NM_FONT, useDingbats = FALSE)
grid.newpage()
txt("ED Fig. 5d, TF activity: three colour schemes, at printed size", 15, 10, 10, "bold")
txt(sprintf("Script 58 (%s): 21 vector-series tumours, construct transcripts left out, pinned CollecTRI; the 25 largest |t| per contrast.",
            if (identical(S58, here::here("results", "orthotopic_tf_activity.rds"))) "the record" else "a dry run"),
    15, 16, 6.5)
y <- 26
for (opt in c("A", "B", "C")) {
  txt(caption[[opt]], 15, y, 7)
  s1 <- sizes[[paste(opt, "vs_bclxl")]]; s2 <- sizes[[paste(opt, "vs_ev")]]
  place(charts[[paste(opt, "vs_bclxl")]], 15, y + 5, s1)
  place(charts[[paste(opt, "vs_ev")]], 15 + s1[["width"]] + 6, y + 5, s2)
  y <- y + 5 + max(s1[["height"]], s2[["height"]]) + 9
}

if (!is.na(ED5) && file.exists(ED5)) {
  bg <- png::readPNG(ED5, native = TRUE)
  for (opt in c("A", "B", "C")) {
    grid.newpage()
    grid.raster(bg, x = 0, y = 0, width = 1, height = 1, just = c("left", "bottom"))
    # the old panel d (two PNG bar charts), covered; the letter d above it stays
    grid.rect(x = mm(24.5), y = mm(PH - 75), width = mm(98), height = mm(38.5),
              just = c("left", "top"), gp = gpar(fill = "white", col = NA))
    s1 <- sizes[[paste(opt, "vs_bclxl")]]; s2 <- sizes[[paste(opt, "vs_ev")]]
    place(charts[[paste(opt, "vs_bclxl")]], 26, 76, s1)
    place(charts[[paste(opt, "vs_ev")]], 26 + s1[["width"]] + 4, 76, s2)
    txt(sprintf("Option %s on the assembled ED 5 (2026-11-06); panel d is %.0f mm wide and %.0f mm tall",
                opt, s1[["width"]] + 4 + s2[["width"]], max(s1[["height"]], s2[["height"]])), 15, 8, 8, "bold")
  }
} else {
  message("no natmetab.ed5_png given: the ED 5 page mock-ups are skipped")
}
invisible(dev.off())
message("wrote ", file.path(OUTD, "p5d_tf_colours.pdf"))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(tf_chart("vs_bclxl", "A"))
  sizes
  pn[pn$bh_sig, c("contrast", "source", "score", "padj_bh")]
}
