# =============================================================================
# p3i_xy_forms.R -- PROPOSAL (2026-10-28): should Fig. 3i be an XY plot?
# -----------------------------------------------------------------------------
# The author (2026-10-28): "I think these figures would be easier to follow if
# presented in XY. Median +/- SD or SEM, or any range which fits in both
# directions ... groups labelled either as a legend or by thin lines pointing to
# the data point."
#
# WHAT THE PAGE IS FOR. Four forms of the same four numbers, at print size, so the
# choice is made from pictures. A is the built categorical panel; B-D are XY.
#
# TWO THINGS THE PAGE IS MEANT TO MAKE VISIBLE, both of which argue about the
# SCIENCE rather than the drawing:
#
#   1. THE FOUR POINTS ARE NOT MONOTONE IN DOSE. KOPgc1a kept 0.17 CPM and sits
#      at +0.27; NTPgc1a kept 128.6 CPM and sits at +0.18, BELOW it. A dose axis
#      puts that inversion on the page. The categorical form, ordered by route,
#      does not show it. Neither form is wrong; they make different things easy.
#
#   2. ox_rel DOES NOT TRACK RETAINED DOSE WITHIN ANY ARM. The analogue of script
#      57 PART 5, computed on what this panel plots (n = 7 per arm, log2 dose):
#        Pgc1a_BclxL  r -0.046 [-0.773, +0.732]
#        NTPgc1a      r +0.172 [-0.667, +0.819]
#        KOPgc1a      r +0.729 [-0.052, +0.957]   (the transgene is LOST here, so
#                                                   its "dose" is 0.10-0.32 CPM)
#        BclxL        r +0.465 [-0.443, +0.902]
#      Every interval spans zero. So a continuous x axis can only be read ACROSS
#      arms, never as a titration within tumours -- the same objection the author
#      raised against a continuous axis for Extended Data a. It does not forbid
#      an XY form; it forbids drawing a fitted line through the four points, and
#      none is drawn here.
#
# THE TWO ERROR BARS DO NOT MEAN THE SAME THING, and that is unavoidable:
#   x  the spread of RETAINED DOSE across the seven animals of an arm (IQR)
#   y  the uncertainty of a BETWEEN-GROUP CONTRAST -- script 52's 5,000-resample
#      bootstrap interval on the difference of medians against the matched control
# y has no per-animal value, because it is a difference against another arm. D
# shows the one form where both axes mean the same thing, at the cost of dropping
# the matched-control framing the section rests on.
# Neither is "median +/- SD or SEM": the saved interval is a bootstrap of the
# difference, and it is used as saved rather than replaced by a fresh statistic.
#
# LOG x. Retained dose runs 0.17 to 364 CPM, so a linear axis stacks three of the
# four points on the left margin (drawn as B1 to show it). A log axis is what
# panel 3g already uses for the same quantity, so the figure stays consistent.
# Every one of the 28 animals is above zero, so nothing is dropped.
#
# A demo page, not a panel: cairo_pdf (for the Greek alpha), never save_panel().
# Base pdf() would draw the hyphens in "Bcl-xL" as long dashes -- the open issue
# on the built panel -- so this page shows the LAYOUT, not that detail.
# Output: outputs/natmetab/proposals/p3i_xy_forms.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))
out_dir <- getOption("natmetab.proposal_dir",
                     here::here("outputs", "natmetab", "proposals"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

o  <- readRDS(here::here("results", "orthotopic_section3_panels.rds"))
pa <- as.data.frame(o$per_animal)
sn <- as.data.frame(o$summary_numbers)

ROUTE <- c("BclxL", "Pgc1a_BclxL", "NTPgc1a", "KOPgc1a")
LINE  <- lw_pt(NM_LINE_PRISM[["Fig3"]])
TXT   <- 5; TITLE <- 6
SYM   <- 1.9
STROKE <- NM_LINE_PRISM[["Fig3"]] * 96 / 72 / (ggplot2::.stroke / 2)
SHORT <- c(BclxL = "no transgene", Pgc1a_BclxL = "effector blocked",
           NTPgc1a = "pathway intact", KOPgc1a = "checkpoint disrupted")

d <- sn[match(ROUTE, sn$arm), ]
d$arm  <- factor(d$arm, levels = ROUTE)
d$excl <- d$ox_rel_hi < 0 | d$ox_rel_lo > 0
d$fill_col <- ifelse(d$excl, background_cols[as.character(d$arm)], "white")
d$lab  <- SHORT[as.character(d$arm)]
# the x spread: IQR of retained dose across the seven animals of each arm
q <- t(vapply(ROUTE, function(a) stats::quantile(pa$Ppargc1a_cpm[pa$arm == a], c(.25, .75)),
              numeric(2)))
d$x_lo <- q[, 1]; d$x_hi <- q[, 2]

base_xy <- function(p, logx = TRUE) {
  p <- p +
    geom_hline(yintercept = 0, linewidth = LINE, colour = box_line) +
    geom_errorbar(aes(ymin = ox_rel_lo, ymax = ox_rel_hi), width = 0, linewidth = LINE) +
    geom_errorbarh(aes(xmin = x_lo, xmax = x_hi), height = 0, linewidth = LINE) +
    geom_point(aes(fill = fill_col), shape = 21, colour = "black",
               size = SYM, stroke = STROKE) +
    scale_fill_identity() +
    scale_y_nm(step = 0.5, top = max(d$ox_rel_hi), bottom = min(d$ox_rel_lo),
               minor = TRUE, labels = lab_signed) +
    labs(y = expression(atop("Relative OXPHOS share", "(vs matched control)"))) +
    theme_nm(txt = TXT, txt_title = TITLE, line = LINE) +
    theme(panel.heights = unit(26, "mm"), panel.widths = unit(46, "mm"),
          plot.title = element_text(size = TITLE, face = "plain", hjust = 0),
          plot.margin = margin(1, 2, 0.5, 0, "mm"))
  if (logx) p + scale_x_log10(limits = c(0.05, 1000),
                              breaks = c(0.1, 1, 10, 100, 1000),
                              labels = c("0.1", "1", "10", "100", "1,000"),
                              expand = c(0, 0)) +
                labs(x = expression(italic(Ppargc1a) ~ "retained (CPM)"))
  else p + scale_x_continuous(limits = c(0, 400), breaks = seq(0, 400, 100),
                              expand = c(0, 0)) +
           labs(x = expression(italic(Ppargc1a) ~ "retained (CPM)"))
}

# --- A: the built categorical panel, for comparison ---------------------------
pA <- ggplot(d, aes(arm, ox_rel_diff)) +
  geom_hline(yintercept = 0, linewidth = LINE, colour = box_line) +
  geom_errorbar(aes(ymin = ox_rel_lo, ymax = ox_rel_hi), width = 0.18, linewidth = LINE) +
  geom_point(aes(fill = fill_col), shape = 21, colour = "black", size = SYM, stroke = STROKE) +
  scale_fill_identity() +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", SHORT[x])) +
  scale_y_nm(step = 0.5, top = max(d$ox_rel_hi), bottom = min(d$ox_rel_lo),
             minor = TRUE, labels = lab_signed) +
  labs(x = NULL, y = expression(atop("Relative OXPHOS share", "(vs matched control)")),
       title = "A  categorical, ordered by route (as built)") +
  theme_nm(txt = TXT, txt_title = TITLE, line = LINE) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = TXT, lineheight = 0.95),
        panel.heights = unit(26, "mm"), panel.widths = unit(46, "mm"),
        plot.title = element_text(size = TITLE, face = "plain", hjust = 0),
        plot.margin = margin(1, 2, 0.5, 0, "mm"))

# --- B1: XY on a LINEAR x, to show why it does not work -----------------------
pB1 <- base_xy(ggplot(d, aes(Ppargc1a_cpm_arm, ox_rel_diff)), logx = FALSE) +
  labs(title = "B1  XY, linear x -- three of four points stack on the margin")

# --- B2: XY on a log x, groups named by thin leader lines ---------------------
# labels placed by hand, clear of the data, each with a hairline to its point
lead <- data.frame(
  arm  = ROUTE,                            # BclxL, Pgc1a_BclxL, NTPgc1a, KOPgc1a
  lx   = c(0.55,  12,   5,    0.055),      # label anchor, data coords
  ly   = c(-0.44, 1.30, 0.62, 0.86),
  hjust = c(0,    0,    0,    0))
lead <- merge(lead, d[, c("arm", "Ppargc1a_cpm_arm", "ox_rel_diff", "lab")], by = "arm")
pB2 <- base_xy(ggplot(d, aes(Ppargc1a_cpm_arm, ox_rel_diff))) +
  geom_segment(data = lead, inherit.aes = FALSE,
               aes(x = lx, y = ly, xend = Ppargc1a_cpm_arm, yend = ox_rel_diff),
               linewidth = LINE * 0.6, colour = "grey45") +
  geom_text(data = lead, inherit.aes = FALSE,
            aes(lx, ly, label = lab, hjust = hjust),
            size = pt(TXT), family = NM_FONT, vjust = -0.4) +
  labs(title = "B2  XY, log x, groups named by leader lines")

# --- C: XY on a log x, groups in a key ----------------------------------------
dC <- d; dC$fill_col <- background_cols[as.character(dC$arm)]   # every point filled
pC <- ggplot(dC, aes(Ppargc1a_cpm_arm, ox_rel_diff)) +
  geom_hline(yintercept = 0, linewidth = LINE, colour = box_line) +
  geom_errorbar(aes(ymin = ox_rel_lo, ymax = ox_rel_hi), width = 0, linewidth = LINE) +
  geom_errorbarh(aes(xmin = x_lo, xmax = x_hi), height = 0, linewidth = LINE) +
  geom_point(aes(fill = lab), shape = 21, colour = "black", size = SYM, stroke = STROKE) +
  scale_fill_manual(values = stats::setNames(unname(background_cols[ROUTE]), unname(SHORT[ROUTE])),
                    name = NULL, limits = unname(SHORT[ROUTE])) +
  scale_x_log10(limits = c(0.05, 1000), breaks = c(0.1, 1, 10, 100, 1000),
                labels = c("0.1", "1", "10", "100", "1,000"), expand = c(0, 0)) +
  scale_y_nm(step = 0.5, top = max(d$ox_rel_hi), bottom = min(d$ox_rel_lo),
             minor = TRUE, labels = lab_signed) +
  labs(x = expression(italic(Ppargc1a) ~ "retained (CPM)"),
       y = expression(atop("Relative OXPHOS share", "(vs matched control)")),
       title = "C  XY, log x, groups in a key (every point filled; see note)") +
  theme_nm(txt = TXT, txt_title = TITLE, line = LINE) +
  theme(panel.heights = unit(26, "mm"), panel.widths = unit(46, "mm"),
        plot.title = element_text(size = TITLE, face = "plain", hjust = 0),
        plot.margin = margin(1, 2, 0.5, 0, "mm"),
        legend.position = "right", legend.text = element_text(size = TXT),
        legend.key.size = unit(2.6, "mm"))

# --- D: both axes mean the same thing -- per-animal, no matched control -------
# every animal of the four arms AND of their four controls, so the pairing is on
# the page; the cost is that the y value is no longer a difference against control
PAIR <- c(BclxL = "EV", Pgc1a_BclxL = "BclxL", NTPgc1a = "NTEV", KOPgc1a = "KOEV")
keep <- unique(c(ROUTE, unname(PAIR)))
dd <- pa[pa$arm %in% keep, ]
med <- do.call(rbind, lapply(keep, function(a) {
  k <- dd$arm == a
  data.frame(arm = a,
             x = stats::median(dd$Ppargc1a_cpm[k]),
             xl = stats::quantile(dd$Ppargc1a_cpm[k], .25),
             xh = stats::quantile(dd$Ppargc1a_cpm[k], .75),
             y = stats::median(dd$ox_rel[k]),
             yl = stats::quantile(dd$ox_rel[k], .25),
             yh = stats::quantile(dd$ox_rel[k], .75),
             is_arm = a %in% ROUTE) }))
seg <- data.frame(a = names(PAIR), b = unname(PAIR))
seg <- merge(merge(seg, med[, c("arm","x","y")], by.x = "b", by.y = "arm"),
             med[, c("arm","x","y")], by.x = "a", by.y = "arm", suffixes = c("_c", "_a"))
pD <- ggplot(med, aes(x, y)) +
  geom_segment(data = seg, inherit.aes = FALSE,
               aes(x = x_c, y = y_c, xend = x_a, yend = y_a),
               linewidth = LINE * 0.7, colour = "grey55") +
  geom_errorbar(aes(ymin = yl, ymax = yh), width = 0, linewidth = LINE) +
  geom_errorbarh(aes(xmin = xl, xmax = xh), height = 0, linewidth = LINE) +
  geom_point(aes(fill = ifelse(is_arm, "black", "white")), shape = 21,
             colour = "black", size = SYM, stroke = STROKE) +
  scale_fill_identity() +
  scale_x_log10(limits = c(0.02, 1000), breaks = c(0.1, 1, 10, 100, 1000),
                labels = c("0.1", "1", "10", "100", "1,000"), expand = c(0, 0)) +
  scale_y_nm(step = 0.5, top = max(med$yh), bottom = min(med$yl), minor = TRUE,
             labels = lab_signed) +
  labs(x = expression(italic(Ppargc1a) ~ "retained (CPM)"),
       y = expression(atop("Relative OXPHOS share", "(arm level, not a difference)")),
       title = "D  both axes per animal (median, IQR); line joins each arm to its control") +
  theme_nm(txt = TXT, txt_title = TITLE, line = LINE) +
  theme(panel.heights = unit(26, "mm"), panel.widths = unit(46, "mm"),
        plot.title = element_text(size = TITLE, face = "plain", hjust = 0),
        plot.margin = margin(1, 2, 0.5, 0, "mm"))

NOTE <- paste(
  "The same four contrasts in four forms, at print size. No fitted line is drawn in any of them:",
  "ox_rel does not track retained dose WITHIN any arm (every interval spans zero, n = 7),",
  "so a dose axis is readable only ACROSS arms.\nWhat an XY form puts on the page and the",
  "categorical form does not: the four points are NOT monotone in dose. KOPgc1a kept 0.17 CPM and",
  "sits at +0.27, ABOVE NTPgc1a, which kept 128.6 CPM and sits at +0.18.\nIn A-C the y bar is",
  "script 52's bootstrap interval on the difference of medians against the matched control, and",
  "the x bar is the IQR of retained dose across seven animals: the two mean different things.",
  "\nOnly D has one meaning on both axes, and it pays for that by plotting arm-level share rather",
  "than a difference against control.\nFilled = the interval excludes zero in A, B1 and B2;",
  "filled = a PGC-1a arm and open = its control in D. C is the exception: a KEY CANNOT SHOW A",
  "COLOUR AND AN OPEN SYMBOL AT ONCE, so every point there is filled by arm and whether the bar",
  "crosses the zero line carries what open/filled carries elsewhere.")

page <- patchwork::wrap_plots(pA, pB1, pB2, pC, pD, patchwork::plot_spacer(), ncol = 2) +
  patchwork::plot_annotation(
    title = "Fig. 3i: categorical against XY (print size; data areas 46 x 26 mm)",
    caption = paste(strwrap(NOTE, width = 135), collapse = "\n"),
    theme = theme(plot.title = element_text(size = 8, family = NM_FONT),
                  plot.caption = element_text(size = 6, hjust = 0, colour = "grey25",
                                              family = NM_FONT)))

f <- file.path(out_dir, "p3i_xy_forms.pdf")
grDevices::cairo_pdf(f, width = 210 / 25.4, height = 200 / 25.4, family = NM_FONT)
print(page)
grDevices::dev.off()

cat("\nwhat each form draws (x = retained CPM, y = ox_rel vs matched control)\n")
print(d[, c("route", "arm", "Ppargc1a_cpm_arm", "x_lo", "x_hi",
            "ox_rel_diff", "ox_rel_lo", "ox_rel_hi", "excl")], row.names = FALSE, digits = 4)
cat("\nwritten:", f, "\n")
