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

for (pkg in c("ggplot2", "ggpubr", "here", "scales", "patchwork")) {
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
# The two working sizes can be set for a whole figure before its panels are sourced,
# e.g. options(natmetab.txt = 5, natmetab.txt_title = 6) to try a figure one point
# smaller (2026-10-12, Fig. 1 page mock-ups). The floor above does not move.
NM_TXT       <- getOption("natmetab.txt", 6)        # pt; tick labels, keys, on-panel labels, P values
NM_TXT_TITLE <- getOption("natmetab.txt_title", 7)  # pt; axis titles and strips (Nature's maximum)
stopifnot(NM_TXT >= NM_TXT_MIN, NM_TXT_TITLE >= NM_TXT, NM_TXT_TITLE <= 7)
NM_LINE      <- 0.25         # ggplot linewidth: axes, ticks, brackets, box outlines (0.53 pt printed)
# The Prism graphs' stroke weight, in POINTS, per figure -- MEASURED from the built
# figure, not assumed, because the two figures' Prism graphs are set differently.
# Fig. 1b was set to 0.9 pt by hand (author, 2026-10-15: 1 pt read heavier than its
# Prism neighbours). Figure 3's panels a-h measure 0.72 pt: the dominant stroke in
# `outputs/natmetab/proposals/NM_Fig3.pdf` and the weight of its symbol outlines
# (2026-10-28). An R panel drawn at NM_LINE (0.53 pt) reads visibly thinner beside
# them. Pass it through lw_pt(): theme_nm(line = lw_pt(NM_LINE_PRISM[["Fig3"]])).
NM_LINE_PRISM <- c(Fig1 = 0.9, Fig3 = 0.72)
NM_MAX_W     <- 183          # mm; double column
NM_ED_MAX_W  <- 180          # mm; Extended Data page width (Nature research figure guide)
NM_MAX_H     <- 170          # mm; tallest Nature figure that leaves room for a legend
nm_widths    <- c(single = 89, double = 183)

pt <- function(x) x / ggplot2::.pt   # points -> ggplot text size units

# Points -> ggplot linewidth, for a line of a stated printed weight. ggplot draws
# linewidth x .pt as R's lwd, and base pdf() writes lwd 1 as 0.75 pt (R's lwd unit is
# 1/96 inch), so linewidth 1 prints 2.13 pt. lw_pt(1) prints 1 pt, as Prism and
# Illustrator measure it (read from the PDF's own line-width operators, 2026-10-15).
lw_pt <- function(x) x * 96 / 72 / ggplot2::.pt

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
group_keys   <- sub("\n", " ", group_labels)          # one line, for keys and legends
geno_cols    <- c(neg = "#0072B2", pos = "#D55E00")   # the 6W pair of group_cols
geno_labels  <- c(neg = "WT", pos = "MYC+")
box_line     <- "#595959"                             # box outlines and medians

# Effects (contrasts) are labelled by what is compared, then where. The genotype
# effects take the MYC+ hues of their age, because that is what they measure.
# Development contrasts get their labels and colours when the first panel needs them.
contrast_cols   <- c(myc_6W = "#D55E00", myc_12W = "#E69F00")
contrast_labels <- c(myc_6W = "MYC+ vs WT, 6W", myc_12W = "MYC+ vs WT, 12W")
stopifnot(identical(unname(contrast_cols), unname(group_cols[c("6W_pos", "12W_pos")])))

# Development contrasts -- 12W vs 6W within one genotype -- take that genotype's
# colour (geno_cols, keyed neg/pos) and these labels (author, 2026-10-11).
dev_labels <- c(neg = "12W vs 6W, WT", pos = "12W vs 6W, MYC+")

# The order follows the test: age-major when the comparisons are within-age gaps,
# genotype-major when the comparison is a genotype main effect.
order_by_age  <- c("6W_neg", "6W_pos", "12W_neg", "12W_pos")
order_by_geno <- c("6W_neg", "12W_neg", "6W_pos", "12W_pos")
stopifnot(setequal(names(group_labels), names(group_cols)),
          setequal(names(group_fill), names(group_cols)))

# --- the cell-line constructs (Fig. 2 on; author, 2026-10-19) ----------------
# Nature's example palette (Okabe-Ito), the half the mouse groups do not use, one hue
# per mechanism: EV light grey, BCL-xL black, biogenesis (PGC-1a, NRF1) bluish green,
# matrix NADH oxidation (Ndi1, LbNOX in the matrix) reddish purple, cytosolic NADH
# oxidation (LbNOX in the cytosol) yellow. Two constructs share a hue only where they
# share a mechanism, and never in one panel; the axis names the construct. MYC state
# (off/on, basal/raised) is on the axis and takes no colour. EV is light grey, not
# mid-grey: #999999 has reddish purple's lightness and merges with it under
# deuteranopia. docs/2026-10-17_fig2_construct_colours.md has the comparison.
# BARS take the hue itself (100%), not the 28% tint: in a bar the fill is the only
# mark (author, 2026-10-19: at 28% "the bar comes almost invisible"). Border black,
# 0.9 pt; symbols white with a black border; error bars upward only.
# Keys are ASCII; labels are set in Prism and Illustrator ("PGC-1" + alpha, "BCL-xL").
# PGC-1a + BCL-xL (Figs 3-4; author, 2026-10-23): the biogenesis green times 0.55,
# CIELAB dE at least 25 from PGC-1a and 31 from BCL-xL under deuteranopia, protanopia
# and tritanopia (proposals/p3_pgc1a_bclxl_colours.R).
construct_cols <- c(EV = "#CCCCCC", BCL_xL = "#000000", PGC1A = "#009E73", NRF1 = "#009E73",
                    PGC1A_BCL_xL = "#00573F",
                    Ndi1 = "#CC79A7", LbNOX_matrix = "#CC79A7", LbNOX_cytosol = "#F0E442")
# No Dox (author, 2026-10-23): on bar graphs with a no-Dox bar, the construct's colour
# at 28% over white (FILL_ALPHA). Set for three hues only. Under deuteranopia these
# tints sit near EV (dE 4 for PGC-1a, 8 for Ndi1) and the PGC-1a and Ndi1 tints near
# each other (dE 4); hatching them in the full colour is open (proposals/p23_nodox_tints.R).
# Where Dox switches MYC (iMMEC rtTA-MYC) the tint marks MYC off, an exception to
# "MYC state takes no colour" above.
construct_nodox <- c(PGC1A = "#B8E4D8", NRF1 = "#B8E4D8", Ndi1 = "#F1D9E6",
                     LbNOX_matrix = "#F1D9E6", LbNOX_cytosol = "#FBF7CA")
local({
  tint <- function(h) {
    v <- grDevices::col2rgb(h)[, 1]
    w <- round(255 - FILL_ALPHA * (255 - v))
    grDevices::rgb(w[1], w[2], w[3], maxColorValue = 255)
  }
  stopifnot(identical(vapply(construct_cols[names(construct_nodox)], tint, character(1)),
                      construct_nodox))
})
stopifnot(!any(c(construct_cols, construct_nodox) %in% c(group_cols, group_fill)))  # never a mouse group's colour

# --- the genetic background (Fig. 3 on; author, 2026-10-28) -------------------
# A SCOPED RULE, NOT AN EXCEPTION. Where a panel holds the CONSTRUCT CONSTANT and
# varies the BACKGROUND, colour encodes the BACKGROUND. The construct-hue rule
# above is untouched everywhere it actually applies. Fig. 3i and its Extended Data
# are the case: all three informative arms carry PGC-1a, so the construct palette
# would give two of them one hue and say nothing about what differs.
# Keyed by ARM, so a panel indexes it directly with the arm names on disk.
# Three of the four ARE their construct's colour and are asserted to stay so; only
# the p19ARF-disrupted background needs a colour of its own.
#   #AA4499 (Paul Tol's muted purple) rather than Okabe-Ito's orange #E69F00: all
#   EIGHT Okabe-Ito colours are already assigned, and #E69F00 is the MYC+ 12W mouse
#   group. Orange would then mean one biological group in Fig. 1 and a different one
#   in Fig. 3. #AA4499 stays at least CIELAB dE 32 from black, #00573F and #009E73
#   under deuteranopia, protanopia and tritanopia, and is free paper-wide.
#   NOT grey: grey reads as "not significant" rather than "different background".
background_cols <- c(BclxL = "#000000", Pgc1a_BclxL = "#00573F",
                     NTPgc1a = "#009E73", KOPgc1a = "#AA4499")
background_labels <- c(BclxL = "no transgene", Pgc1a_BclxL = "effector blocked",
                       NTPgc1a = "pathway intact", KOPgc1a = "checkpoint disrupted")
stopifnot(
  identical(names(background_cols), names(background_labels)),
  # the three that reuse a construct colour must keep tracking it
  identical(unname(background_cols[c("BclxL", "Pgc1a_BclxL", "NTPgc1a")]),
            unname(construct_cols[c("BCL_xL", "PGC1A_BCL_xL", "PGC1A")])),
  !any(background_cols %in% c(group_cols, group_fill)))   # never a mouse group's colour

# --- diverging and sequential fills (author's specification, 2026-07-30) ------
# Espresso -> white -> mint, white pinned at zero. Not blue-red: blue and orange
# already mean genotype. Each arm is scaled to its own end of `limits`, so on a
# lopsided range equal ink is not equal magnitude across zero -- use the
# asymmetric form only where an axis carries the magnitude, and show the bar.
div_cols <- c(neg = "#4A3525", zero = "#FAFAFA", pos = "#2A8A6D")
seq_cols <- c(low = "#F2F4F6", mid = "#8FA3B0", high = "#243642")  # unsigned levels

# Mitochondrial class of a gene set (ED 2b-c), script 37's three classes. The two
# mitochondrial classes take one hue, reddish purple (Okabe-Ito), darker = defined
# on MitoCarta; everything else is neutral grey, so mitochondrial against not reads
# first (author, 2026-10-14: the slate shades did not separate them). Not blue or
# orange, which mean genotype, and not green, which is the diverging scale's end.
# The middle class is mitochondrial BY CONSTRUCTION: transcription-factor target
# sets intersected with MitoCarta. Labels avoid hyphens (base pdf() draws a minus).
# options(natmetab.mito_class_cols = c(...)) tries another set (same three names).
mito_class_cols   <- getOption("natmetab.mito_class_cols",
                               c(mitocarta_proper = "#8E2C74", construction_MITO = "#CC79A7",
                                 non_mito = "#D9D9D9"))
stopifnot(identical(names(mito_class_cols), c("mitocarta_proper", "construction_MITO", "non_mito")))
mito_class_labels <- c(mitocarta_proper = "MitoCarta",
                       construction_MITO = "Curated (with mitochondria)",
                       non_mito = "Not mitochondrial")

# MitoPathway tiers (ED 2f-g), the seven top-level MitoCarta categories, coloured "as in
# the original panels" (author, 2026-10-16: fig02b and figS6): Okabe-Ito, Metabolism grey.
# Here blue and orange are tiers, not genotypes; these panels carry no genotype colour.
# Keys are script 08's tier names (results/mitopps_scores.rds $pathway_tier1_map).
mito_tier_labels <- c("Protein import, sorting and homeostasis" = "Protein import and homeostasis",
                      "Mitochondrial central dogma"             = "Central dogma",
                      "OXPHOS"                                  = "OXPHOS",
                      "Metabolism"                              = "Metabolism",
                      "Signaling"                               = "Signalling",
                      "Mitochondrial dynamics and surveillance" = "Dynamics and surveillance",
                      "Small molecule transport"                = "Small molecule transport")
mito_tier_cols <- stats::setNames(c("#D55E00", "#E69F00", "#009E73", "#B8B8B8", "#56B4E9",
                                    "#0072B2", "#CC79A7"), mito_tier_labels)
# The same names as point labels, the long ones on two lines so that four fit along
# the empty corner of a 120-mm data area (ranked_tier_marks()).
mito_tier_text <- stats::setNames(c("Protein import\nand homeostasis", "Central dogma", "OXPHOS",
                                    "Metabolism", "Signalling", "Dynamics and\nsurveillance",
                                    "Small molecule\ntransport"), mito_tier_labels)

# ED 2f-g's marks on a curve of MitoPathways ranked from promoted to demoted (author,
# 2026-10-16: "as in fig02_reallocation_ranked.pdf"). In data units, on a data area
# data_w x data_h mm whose ranges xr, yr are read from the built plot:
#   anchors  (rank, y, tier, label): the seven top-level MitoPathways, the tiers
#            themselves. Each is ringed and named in black, plain, at txt pt (larger
#            than the panel's 6 pt), its name on one of two lines of labels: promoted
#            pathways (y > 0) along the foot of the empty lower left, bottom-aligned,
#            demoted ones along the head of the empty upper right, top-aligned, so
#            that the axis stays the data's height. Each line runs only where no
#            point comes within `clear` mm of it; spread_on_line() sets the names on
#            it, `gap` mm apart; a leader joins each name to its ring.
#   stars    (rank, y, up, r_mm): an asterisk above (up) or below each significant
#            point, clear of its dot (radius r_mm; a ringed one, of its ring). Nature
#            asks for exact P values, not asterisks; the author makes these two
#            panels the exception, as the legend must say.
#   pts      (rank, y): every point drawn, which the lines of labels keep off.
# Returns the layers; the caller maps the ring fills with scale_fill_manual().
ranked_tier_marks <- function(anchors, stars, pts, xr, yr, data_w, data_h, txt = NM_TXT_TITLE,
                              ring = 2, clear = 1, gap = 2.5, edge = 1, lineheight = 0.9) {
  stopifnot(all(c("rank", "y", "tier", "label") %in% names(anchors)),
            all(c("rank", "y", "up", "r_mm") %in% names(stars)),
            all(c("rank", "y") %in% names(pts)), txt >= NM_TXT_MIN, txt <= 7)
  ux <- diff(xr) / data_w; uy <- diff(yr) / data_h           # data units per mm
  em <- txt * 25.4 / 72                                        # one em at txt pt, in mm
  r_ring <- point_r_mm(ring, stroke = 0.4)
  # Helvetica: ascenders and capitals 0.718 em above the baseline, descenders 0.21 below
  nl  <- lengths(strsplit(anchors$label, "\n", fixed = TRUE))
  asc <- ((nl - 1) * lineheight + 0.718) * em                  # above the last baseline
  dsc <- 0.21 * em
  w   <- text_mm(anchors$label, txt, "plain", lineheight)[, "w"]
  bot <- anchors$y > 0
  y_bot <- yr[1] + (edge + dsc) * uy                          # the bottom line's last baseline
  y_top <- yr[2] - edge * uy                                   # the top line's top
  lo_y  <- tapply(pts$y, pts$rank, min); hi_y <- tapply(pts$y, pts$rank, max)
  rk    <- as.numeric(names(lo_y))
  cx    <- numeric(nrow(anchors))
  for (b in c(TRUE, FALSE)) {
    k <- bot == b
    if (!any(k)) next
    blocked <- if (b) rk[lo_y < y_bot + (max(asc[k]) + clear) * uy]
               else   rk[hi_y > y_top - (max(asc[k]) + dsc + clear) * uy]
    a <- anchors$rank[k]
    stopifnot(!any(blocked >= min(a) & blocked <= max(a)))   # the line is free around its anchors
    lo <- max(c(xr[1] + edge * ux, blocked[blocked < min(a)] + 0.5 + clear * ux))
    hi <- min(c(xr[2] - edge * ux, blocked[blocked > max(a)] - 0.5 - clear * ux))
    cx[k] <- spread_on_line(a, w[k] * ux, lo, hi, gap * ux)
  }
  lab <- data.frame(x = cx, y = ifelse(bot, y_bot, y_top), vjust = ifelse(bot, 0, 1),
                    label = anchors$label)
  # leaders: from just off each name to just off its ring
  sy <- ifelse(bot, y_bot + (asc + 0.4) * uy, y_top - (asc + dsc + 0.4) * uy)
  dx <- (anchors$rank - cx) / ux; dy <- (anchors$y - sy) / uy
  f  <- pmax(0, 1 - (r_ring + 0.3) / sqrt(dx^2 + dy^2))
  lead <- data.frame(x = cx, y = sy, xend = cx + f * dx * ux, yend = sy + f * dy * uy)
  # asterisks: the glyph (0.431-0.718 em above its baseline) centred 0.25 mm off the dot
  ringed <- paste(stars$rank, stars$y) %in% paste(anchors$rank, anchors$y)
  r  <- ifelse(ringed, r_ring, stars$r_mm)
  ct <- stars$y + ifelse(stars$up, 1, -1) * (r + 0.25 + 0.1435 * em) * uy
  st <- data.frame(x = stars$rank, y = ct - 0.5745 * em * uy)
  list(
    ggplot2::geom_segment(data = lead, ggplot2::aes(x, y, xend = xend, yend = yend),
                          inherit.aes = FALSE, linewidth = 0.2, colour = "grey55"),
    ggplot2::geom_point(data = anchors, ggplot2::aes(rank, y, fill = tier), inherit.aes = FALSE,
                        shape = 21, size = ring, colour = "black", stroke = 0.4, show.legend = FALSE),
    ggplot2::geom_text(data = st, ggplot2::aes(x, y), label = "*", inherit.aes = FALSE, vjust = 0,
                       size = pt(txt), family = NM_FONT, colour = "black"),
    ggplot2::geom_text(data = lab, ggplot2::aes(x, y, label = label, vjust = vjust), inherit.aes = FALSE,
                       size = pt(txt), family = NM_FONT, lineheight = lineheight,
                       colour = "black"))
}

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
# x 10^k below 0.001. The exponent is a quoted string, because plotmath puts a gap
# after an unquoted minus, inside a style named on purpose:
#   exponent = "full"     displaystyle(): the exponent at the label's size (the default)
#   exponent = "reduced"  scriptstyle(): at 70% of it, which base pdf() rounds to a
#                         whole point: 4 pt inside a 5 pt label, read from the PDF
#                         (author, 2026-10-15, Fig. 1d: "superscript size should be reduced")
# check_panel() refuses a superscript in neither style. Plotmath keeps the source ASCII.
fmt_p <- function(p, exponent = c("full", "reduced")) {
  style <- c(full = "displaystyle", reduced = "scriptstyle")[[match.arg(exponent)]]
  stopifnot(is.numeric(p), !anyNA(p), all(p >= 0 & p <= 1))
  vapply(p, function(v) {
    if (v >= 0.001) {
      return(sprintf("italic(P)=='%s'", formatC(signif(v, 2), format = "fg", digits = 2, flag = "#")))
    }
    if (v < 1e-300) return(sprintf("italic(P)<10^%s('-300')", style))
    e <- floor(log10(v))
    m <- signif(v / 10^e, 2)
    if (m >= 10) { m <- m / 10; e <- e + 1 }
    sprintf("italic(P)=='%s'%%*%%10^%s('%d')", format(m), style, e)
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
#   base and step may be one value, or one per bracket -- for facets that each
#   have their own y scale (comps then carries the facet column).
#   group1/group2 may also be numeric x positions, e.g. 1.5 and 3.5 to span the
#   two WT boxes against the two MYC+ boxes.
#   txt, line  the P label's size (pt) and the line's linewidth
#   tips       FALSE draws the bar alone, without the downward ends (the Prism format,
#              Fig. 1b). ggpubr would still draw each tip as a zero-length segment, a
#              stray point in Illustrator, so the bar is all that is drawn.
#   exponent   "full" or "reduced": the size of a P value's exponent (fmt_p)
#
# Returns the layer, with attr "top": the y the axis must reach so the top label
# is not clipped. Pass it on: scale_y_nm(step = ..., top = attr(br, "top")).
# attr "tops" gives the same per bracket, for per-facet scales.
brackets <- function(comps, base, step, txt = NM_TXT, line = NM_LINE, tips = TRUE,
                     exponent = "full") {
  stopifnot(is.data.frame(comps), all(c("group1", "group2", "p") %in% names(comps)),
            nrow(comps) >= 1L, is.numeric(base), length(base) %in% c(1L, nrow(comps)),
            is.numeric(step), length(step) %in% c(1L, nrow(comps)), all(step > 0),
            is.numeric(txt), is.numeric(line), is.logical(tips))
  # Default stacking: by row order when all brackets share one base; level 1 each
  # when every bracket has its own base (one per facet), or they would be stacked
  # across facets as if they shared a panel.
  if (is.null(comps$level)) {
    comps$level <- if (length(base) > 1L) 1 else seq_len(nrow(comps))
  }
  comps$y.position <- base + (comps$level - 1) * step
  comps$label <- fmt_p(comps$p, exponent)
  layer <- ggpubr::geom_bracket(
    data = comps, inherit.aes = FALSE,
    mapping = ggplot2::aes(xmin = group1, xmax = group2,
                           y.position = y.position, label = label),
    type = "expression", label.size = pt(txt), linewidth = line,
    tip.length = if (tips) 0.012 else 0, vjust = -0.15, colour = "black", family = NM_FONT)
  if (!tips) {
    # ggpubr's stat gives each bracket three segments (tip, bar, tip); with no tips the
    # tips have x == xend. Drawing the bar alone puts the label where ggpubr would: the
    # midpoint of the first and last segments, which are now the bar's two ends.
    with_tips <- layer$geom
    layer$geom <- ggplot2::ggproto("GeomBracketBar", with_tips,
      draw_group = function(data, panel_params, coord, type = "text", coord.flip = FALSE,
                            orientation = "horizontal") {
        with_tips$draw_group(data[data$x != data$xend, , drop = FALSE], panel_params, coord,
                             type = type, coord.flip = coord.flip, orientation = orientation)
      })
  }
  tops <- comps$y.position + 0.9 * step
  structure(list(layer), top = max(tops), tops = tops)
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
  # Everything is rounded to 10 decimals: in floating point -0.2 + 3 * 0.2 is
  # 0.4000000000000001, a hair past an axis that ends at 0.4, and ggplot then drops
  # that end tick as out of range (caught by check_panel on Fig. 1f, 2026-10-11).
  r  <- function(x) round(x, 10)
  lo <- r(floor(bottom / unit + 1e-9) * unit)
  hi <- r(ceiling(top / unit - 1e-9) * unit)
  list(limits = c(lo, hi),
       breaks = r(seq(ceiling(lo / step - 1e-9) * step, hi + 1e-9 * step, by = step)),
       minor_breaks = if (minor) r(seq(lo, hi + 1e-9 * unit, by = unit)) else NULL,
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

# A box plot's y axis need not start at 0 (author, 2026-10-15: "for box and whiskers plots
# it is acceptable that the y axis does not start at 0"). For resolution it runs from the
# tick closest below the lowest whisker to the first tick above the highest mark, which
# may be a bracket's P label; whiskers run min to max, so the whiskers end at the lowest
# and highest animals (lo, hi). Ticks are labelled (step) or minor (halfway, drawn).
# Worked in printed mm on a graph `graph` mm tall: `above` mm are kept over the highest
# animal (a bracket and its label; label_mm() measures a label) and `below` mm under the
# lowest. With below = 0, the author's rule, the lowest point can sit on the axis end:
# draw with coord_cartesian(clip = "off") so it is not cut. Returns .nm_axis()'s list,
# for scale_y_continuous(), one panel or one facet at a time.
box_axis <- function(lo, hi, step, graph, below = 0, above = 0) {
  stopifnot(is.numeric(lo), is.numeric(hi), hi > lo, below >= 0, above >= 0, graph > below + above)
  s <- (hi - lo) / (graph - below - above)       # data units per printed mm, at the tightest
  repeat {
    a <- .nm_axis(step, top = hi + above * s, bottom = lo - below * s, minor = TRUE)
    # rounding out to ticks lengthens the axis, so each mm is more data: work it again
    # on the axis just found until that axis keeps its own margins
    s_axis <- diff(a$limits) / graph
    if (s_axis <= s * (1 + 1e-9)) return(a)
    s <- s_axis
  }
}

# --- the theme -------------------------------------------------------------------
# Titles, subtitles and captions are blanked: on a Nature panel the explanation
# belongs in the legend, so a panel cannot acquire explanatory text by accident.
# txt, txt_title and line default to the R panels' sizes; a panel drawn in another
# format passes its own (Fig. 1b, the Prism format: 5, 6 and lw_pt(1)).
theme_nm <- function(legend = "none", txt = NM_TXT, txt_title = NM_TXT_TITLE, line = NM_LINE) {
  stopifnot(txt >= NM_TXT_MIN, txt_title >= txt, txt_title <= 7)
  ggplot2::theme_classic(base_size = txt_title, base_family = NM_FONT) +
    ggplot2::theme(
      text              = ggplot2::element_text(colour = "black"),
      axis.text         = ggplot2::element_text(size = txt, colour = "black"),
      axis.title        = ggplot2::element_text(size = txt_title, colour = "black"),
      axis.line         = ggplot2::element_line(linewidth = line, colour = "black"),
      axis.ticks        = ggplot2::element_line(linewidth = line, colour = "black"),
      axis.ticks.length = ggplot2::unit(1, "mm"),
      axis.minor.ticks.length = ggplot2::rel(0.5),
      legend.text       = ggplot2::element_text(size = txt),
      legend.title      = ggplot2::element_text(size = txt),
      legend.key.size   = ggplot2::unit(3, "mm"),
      legend.key.spacing.y = ggplot2::unit(0.4, "mm"),
      legend.margin     = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = ggplot2::unit(2, "mm"),     # key to plot; the 7-pt default is 2.5 mm
      legend.position   = legend,
      strip.text        = ggplot2::element_text(size = txt_title),
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

check_panel <- function(p, width, height, max_w = NM_MAX_W) {
  errs <- character(0)
  add <- function(...) errs <<- c(errs, sprintf(...))

  # page size (Extended Data pages are 180 mm wide, main figures 183)
  if (width > max_w)  add("width %g mm is over the %g mm page", width, max_w)
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
    bracket <- inherits(p$layers[[i]]$geom, "GeomBracket")   # tipless brackets too
    if (!bracket && !g %in% .nm_text_geoms) next
    d <- b$data[[i]]
    sz <- if (bracket) d$label.size else d$size
    if (length(sz) && any(sz * ggplot2::.pt < NM_TXT_MIN - 1e-9)) {
      add("layer %d (%s): text at %.2f pt", i, g, min(sz) * ggplot2::.pt)
    }
    if (!.nm_achromatic(d$colour)) add("layer %d (%s): coloured text", i, g)
    if (!is.null(d$family) && any(!d$family %in% c("", NM_FONT))) {
      add("layer %d (%s): font is not %s", i, g, NM_FONT)
    }
    # a superscript prints at the size fmt_p() names; a bare one comes out at ~70% by
    # default, with a gap after its minus
    if (!is.null(d$label) && any(grepl("\\^(?!(displaystyle|scriptstyle)\\()", as.character(d$label), perl = TRUE))) {
      add("layer %d (%s): a superscript in no named style would print at ~70%% size", i, g)
    }
  }

  # boxes: whiskers run minimum to maximum (author, 2026-10-14). Prism draws every point
  # only over min-max whiskers, so R matches it and one legend sentence covers both:
  # "boxes show the median and interquartile range, whiskers the minimum and maximum".
  for (i in seq_along(p$layers)) {
    if (!inherits(p$layers[[i]]$geom, "GeomBoxplot")) next
    cf <- p$layers[[i]]$stat_params$coef
    if (is.null(cf) || !is.infinite(cf)) {
      add("layer %d (box): whiskers must run min to max -- geom_boxplot(coef = Inf)", i)
    }
  }

  # axes: both ends on a tick; y from 0 (or spanning 0) whenever an x axis is drawn,
  # except under a box plot, whose axis may start above 0 (author, 2026-10-15; box_axis())
  drawn <- function(ax) !all(vapply(sprintf(c("axis.text.%s", "axis.line.%s"),
                                            if (ax == "x") "x.bottom" else "y.left"),
    function(e) inherits(ggplot2::calc_element(e, thm), "element_blank"), logical(1)))
  x_drawn <- drawn("x")
  boxes   <- any(vapply(p$layers, function(l) inherits(l$geom, "GeomBoxplot"), logical(1)))
  yv <- unlist(lapply(b$data, function(d)
    unlist(d[intersect(names(d), c("y", "ymin", "ymax", "lower", "upper", "middle"))],
           use.names = FALSE)), use.names = FALSE)
  for (k in seq_along(b$layout$panel_params)) {
    pp <- b$layout$panel_params[[k]]
    for (ax in c("x", "y")) {
      v <- pp[[ax]]
      if (is.null(v) || v$is_discrete()) next
      if (!drawn(ax)) next          # an undrawn axis (a schematic) has no ticks to end on
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
        } else if (r[1] > tol && !boxes) {
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

# --- panels made of stacked parts -------------------------------------------------
# A panel can be two or more plots stacked vertically, e.g. a histogram over a
# strip of points that shares its x axis. Pass save_panel() stack_parts(...):
# every part is checked on its own, then patchwork stacks them and aligns their
# panels. Give the parts the same x scale (scale_x_nm with the same arguments)
# when they share an x axis; patchwork aligns panel edges, not scales.
# collect_guides = TRUE gathers the parts' keys into one shared key at the bottom.
# ncol = 2 (with widths) puts parts side by side instead, e.g. two rulers sharing
# the same rows.
# widths and heights may be grid units: unit(44, "mm") fixes a part's DATA AREA
# at 44 mm whatever its axis labels need, so panels in different scripts come out
# the same size and align on the page (ED 2d and 2e). A single part is allowed
# for that purpose: fixed_size(p, 44, 44).
stack_parts <- function(..., heights = NULL, widths = NULL, ncol = 1L,
                        collect_guides = FALSE) {
  parts <- list(...)
  stopifnot(length(parts) >= 1L, all(vapply(parts, inherits, logical(1), "ggplot")),
            is.null(heights) || length(heights) == ceiling(length(parts) / ncol),
            is.null(widths) || length(widths) == ncol)
  structure(list(parts = parts, heights = heights, widths = widths, ncol = ncol,
                 collect = isTRUE(collect_guides)),
            class = "nm_stack")
}
fixed_size <- function(p, width_mm, height_mm) {
  stack_parts(p, widths = grid::unit(width_mm, "mm"), heights = grid::unit(height_mm, "mm"))
}

# The printed size a fixed-size panel needs: its data areas plus its axes, titles and
# margins, in mm, measured with the PDF device's Helvetica metrics. NA where a row or
# column stretches (a relative width or height), which takes whatever room is left.
# A null inside a sum does not stretch its row: patchwork centres a collected key that
# way, and grid's layout reads the null as 0, as convertHeight() does here (Fig. 1d).
# A smaller page cuts the panel's edges; save_panel() refuses one. A larger page
# centres the panel.
natural_size <- function(p) {
  stopifnot(inherits(p, "nm_stack"))
  grDevices::pdf(NULL, family = NM_FONT)
  on.exit(grDevices::dev.off())
  g <- patchwork::patchworkGrob(.nm_drawable(p))
  fixed <- function(u) !any(grid::unitType(u) == "null")
  c(width  = if (fixed(g$widths))  sum(grid::convertWidth(g$widths, "mm", valueOnly = TRUE))  else NA,
    height = if (fixed(g$heights)) sum(grid::convertHeight(g$heights, "mm", valueOnly = TRUE)) else NA)
}

# The printed height of a plotmath label (a P value from fmt_p(), say) at txt pt, in mm,
# measured with the PDF device's Helvetica metrics: for an axis that ends just above a
# bracket's label (Fig. 1d). 5-pt P labels measure 1.30 mm, and 1.90 mm with a reduced
# exponent.
label_mm <- function(label, txt = NM_TXT) {
  stopifnot(is.character(label), is.numeric(txt))
  grDevices::pdf(NULL, family = NM_FONT)
  on.exit(grDevices::dev.off())
  vapply(label, function(l) grid::convertHeight(grid::grobHeight(grid::textGrob(
    parse(text = l), gp = grid::gpar(fontsize = txt, fontfamily = NM_FONT))), "mm", valueOnly = TRUE),
    numeric(1), USE.NAMES = FALSE)
}

# The printed size of plain text at txt pt, in mm, measured with the PDF device's
# Helvetica metrics: w, the widest line; h, the block (for a single line, without
# descenders). One row per label.
text_mm <- function(label, txt = NM_TXT, fontface = "plain", lineheight = 0.9) {
  stopifnot(is.character(label), is.numeric(txt))
  grDevices::pdf(NULL, family = NM_FONT)
  on.exit(grDevices::dev.off())
  gp <- grid::gpar(fontsize = txt, fontfamily = NM_FONT, fontface = fontface, lineheight = lineheight)
  m <- vapply(label, function(l) {
    g <- grid::textGrob(l, gp = gp)
    c(grid::convertWidth(grid::grobWidth(g), "mm", valueOnly = TRUE),
      grid::convertHeight(grid::grobHeight(g), "mm", valueOnly = TRUE))
  }, numeric(2), USE.NAMES = FALSE)
  matrix(m, ncol = 2L, byrow = TRUE, dimnames = list(NULL, c("w", "h")))
}

# The printed radius of a ggplot point of the given size and stroke, its outline
# included, in mm: R draws the circle at 0.375 of the point's font size, which
# ggplot sets to size * .pt + stroke * .stroke / 2, and the outline, as wide as that
# stroke term, straddles it.
point_r_mm <- function(size, stroke = 0) {
  (0.375 * (size * ggplot2::.pt + stroke * ggplot2::.stroke / 2) + stroke * ggplot2::.stroke / 4) * 25.4 / 72
}

# Labels set along one line, each as near its anchor as the others allow (ED 2f-g):
# the centres of labels `width` wide that would each like to be centred on `anchor`,
# at least `gap` apart and inside [lo, hi], all in one unit; in the input order.
# The left edges, less the room the labels before each one take, must not decrease,
# so the least-squares answer is an isotonic regression; the bounds then clip it.
spread_on_line <- function(anchor, width, lo, hi, gap = 0) {
  n <- length(anchor)
  stopifnot(n >= 1L, length(width) == n, all(width > 0), gap >= 0,
            sum(width) + gap * (n - 1) <= hi - lo)
  o   <- order(anchor)
  w   <- width[o]
  off <- c(0, cumsum(w + gap))[seq_len(n)]
  q   <- anchor[o] - w / 2 - off
  if (n > 1L) q <- stats::isoreg(q)$yf
  q   <- pmin(pmax(q, lo), hi - w[n] - off[n])
  out <- numeric(n)
  out[o] <- q + off + w / 2
  out
}

.nm_drawable <- function(p) {
  if (!inherits(p, "nm_stack")) return(p)
  w <- patchwork::wrap_plots(p$parts, ncol = p$ncol, heights = p$heights, widths = p$widths,
                             guides = if (p$collect) "collect" else "auto")
  # patchwork lays out the collected key, and a margin around the whole stack, with
  # its own theme: theme_get() at 11 pt unless one is given. Give it the panel theme
  # and no outer margin, so the key sits as close as a single plot's would.
  w <- w + patchwork::plot_annotation(
    theme = theme_nm(legend = if (p$collect) "bottom" else "none") +
      ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0)))
  if (p$collect) w <- w & ggplot2::theme(legend.position = "bottom")
  w
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
  max_w <- if (grepl("^EDFig", fig)) NM_ED_MAX_W else NM_MAX_W
  if (inherits(p, "nm_stack")) {
    for (part in p$parts) check_panel(part, width, height, max_w)
    need <- natural_size(p)                    # NA where a row or column stretches
    if (isTRUE(any(c(width, height) < need - 0.05))) {
      stop(sprintf("panel refused:\n  - saved at %g x %g mm, it needs %.1f x %.1f mm: its edges would be cut",
                   width, height, need[["width"]], need[["height"]]), call. = FALSE)
    }
  } else {
    check_panel(p, width, height, max_w)
  }
  p <- .nm_drawable(p)
  # options(natmetab.capture = TRUE) keeps the checked, drawable panel in
  # .nm_captured (global), for assembling page mock-ups from the panel scripts.
  if (isTRUE(getOption("natmetab.capture"))) assign(".nm_captured", p, envir = globalenv())
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
