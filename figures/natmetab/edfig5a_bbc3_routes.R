# =============================================================================
# edfig5a_bbc3_routes.R -- Extended Data, Bbc3 across the three escape routes,
#                          with the p53 specificity strip
# -----------------------------------------------------------------------------
#   main   Bbc3, log2 difference against each arm's own MATCHED control, three
#          backgrounds including the one where the transgene was lost.
#   strip  the specificity control: Cdkn2a and the four canonical p53 targets,
#          same arms, same symbol conventions.
#
# THREE CATEGORIES, NOT A CONTINUOUS DOSE AXIS (author, 2026-10-28). Script 57
#   PART 5 found Bbc3 does NOT track Ppargc1a within Pgc1a_BclxL across its seven
#   animals (Pearson r -0.165, interval spanning zero), so a continuous dose axis
#   would assert a dose-response the data do not show. These are three
#   backgrounds that differ in genotype as well as in retained dose. The retained
#   dose is a SUB-LABEL beneath each route, not a coordinate.
#   This REVERSES the draft, which ran left to right in ascending dose.
#
# ROUTE ORDER, THE SAME AS Fig. 3i: effector blocked, pathway intact, checkpoint
#   disrupted. The zero arm is not here -- it carries no transgene, so it has no
#   retained dose and script 52 fitted the p53 panel only on the three clean
#   within-series contrasts.
#
# COLOUR ENCODES THE BACKGROUND (`background_cols`), as in Fig. 3i.
#
# ALL FOUR p53 TARGETS ARE LABELLED, Cdkn1a included. It is also in Fig. 3f,
#   falling in the in vitro surviving fraction: down in culture, flat in tumours,
#   which is what the narrative predicts because those arms differ in buffer
#   status. It is not dropped to avoid the collision.
#
# Cdkn2a IS SHOWN FOR THE VECTOR AND NT ARMS ONLY (author, 2026-10-28). The locus
#   is CRISPR-targeted in the p19ARF arm, so its value there is not a read-out of
#   p53 output. The legend says so; the panel carries no note, which leaves that
#   slot empty by design.
#
# NO ANALYSIS AND NO CAPTION. Every value is re-read from the saved object;
#   nothing is re-fitted and no interval is recomputed. Legends describe panels
#   only, so the interpretation is not on the page.
#
# THE FIGURE NUMBER IS PROVISIONAL. The spec carries Extended Data Fig. 5 as a
#   placeholder (section 4, item 8): the manuscript has ED 1-4 plus two unnumbered
#   figures in the human section. When the numbering settles this is a `git mv`
#   and a one-line change to the save_panel() call.
#
# PRISM FORMAT, as Fig. 3i: data area 20 mm tall, axis title 6 pt, tick labels
#   5 pt, every line 0.72 pt (NM_LINE_PRISM[["Fig3"]], measured from the built
#   figure).
#
# Reads : results/orthotopic_section3_panels.rds  (script 57, which reads script 52)
# Output: outputs/natmetab/EDFig5/EDFig5a_bbc3_routes.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

o  <- readRDS(here::here("results", "orthotopic_section3_panels.rds"))
bb <- as.data.frame(o$eda_bbc3)
st <- as.data.frame(o$eda_strip)
sn <- as.data.frame(o$summary_numbers)

ROUTE <- c("Pgc1a_BclxL", "NTPgc1a", "KOPgc1a")       # by route, not by dose
stopifnot(setequal(ROUTE, bb$arm), all(ROUTE %in% names(background_cols)))

# the retained dose beneath each route: 301.6, 128.6, 0.17 -- the arm medians as
# script 52 saved them, formatted so the two large ones keep the decimal the
# record quotes and the lost one does not gain spurious digits
dose <- stats::setNames(sn$Ppargc1a_cpm_arm[match(ROUTE, sn$arm)], ROUTE)
dose_txt <- ifelse(dose >= 1, sprintf("%.1f CPM", dose), sprintf("%.2f CPM", dose))
x_lab <- parse(text = sprintf('atop("%s", "%s")', background_labels[ROUTE], dose_txt))

GRAPH_H   <- 20
CAT_W     <- 20                                 # mm per route on the main plot, as Fig. 3i
STRIP_W   <- 14.4                               # mm per gene in the strip (three short labels)
STRIP_GAP <- 1.6                                # mm between the strip's facets
LINE    <- lw_pt(NM_LINE_PRISM[["Fig3"]])
TXT     <- 5
TITLE   <- 6
SYM     <- 1.6
STROKE  <- NM_LINE_PRISM[["Fig3"]] * 96 / 72 / (ggplot2::.stroke / 2)

# One symbol convention for both plots: filled where the interval excludes zero,
# open where it does not; the fill is the background, the outline always black.
# The fill is set per row and mapped with scale_fill_identity(), NOT with alpha:
# ggplot applies alpha to a shape-21 point's OUTLINE as well as its fill, so
# alpha = 0 renders the whole symbol transparent and an open point disappears
# (checked, 2026-10-28).
fill_of <- function(d) ifelse(d$hi < 0 | d$lo > 0,
                              background_cols[as.character(d$arm)], "white")
sym_layers <- function(size = SYM) list(
  geom_hline(yintercept = 0, linewidth = LINE, colour = box_line),
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = LINE, colour = "black"),
  geom_point(aes(fill = fill_col), shape = 21, colour = "black",
             size = size, stroke = STROKE),
  scale_fill_identity())

# =============================================================================
# MAIN -- Bbc3
# =============================================================================
bb$arm  <- factor(bb$arm, levels = ROUTE)
bb      <- bb[order(bb$arm), ]
bb$excl <- bb$hi < 0 | bb$lo > 0
bb$fill_col <- fill_of(bb)
stopifnot(nrow(bb) == 3L, !anyNA(bb$diff))

p_main <- ggplot(bb, aes(arm, diff)) +
  sym_layers() +
  scale_x_discrete(labels = x_lab) +
  scale_y_nm(step = 0.5, top = max(bb$hi), bottom = min(bb$lo), minor = TRUE,
             labels = lab_signed) +
  labs(x = NULL, y = expression(atop(italic(Bbc3) ~ "(log"[2] ~ "difference",
                                     "vs matched control)"))) +
  theme_nm(txt = TXT, txt_title = TITLE, line = LINE) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = TXT, lineheight = 0.95, vjust = 1),
        panel.heights = unit(GRAPH_H, "mm"),
        panel.widths  = unit(CAT_W * length(ROUTE), "mm"),
        plot.margin = margin(1, 1.5, 0.5, 0, "mm"))

# =============================================================================
# STRIP -- Cdkn2a and the four canonical p53 targets
# -----------------------------------------------------------------------------
# Cdkn2a leads, because it is the one that moves; the four that do not follow it
# come after, in the order script 52 ran them. The KO arm is absent from Cdkn2a
# only, which leaves that slot empty -- the legend accounts for it.
# Short route words here: the full two-line labels will not fit five times over.
# =============================================================================
GENES <- c("Cdkn2a", "Cdkn1a", "Mdm2", "Ccng1", "Zmat3")
SHORT <- c(Pgc1a_BclxL = "blocked", NTPgc1a = "intact", KOPgc1a = "disrupted")
st$arm  <- factor(st$arm, levels = ROUTE)
st$gene <- factor(as.character(st$gene), levels = GENES)
st$excl <- st$hi < 0 | st$lo > 0
st$fill_col <- fill_of(st)
stopifnot(!anyNA(st$gene), !anyNA(st$arm),
          nrow(st) == 14L,                                   # 5 x 3, less Cdkn2a in the KO arm
          sum(st$gene == "Cdkn2a") == 2L,
          !any(st$gene == "Cdkn2a" & st$arm == "KOPgc1a"))

p_strip <- ggplot(st, aes(arm, diff)) +
  sym_layers(size = SYM - 0.3) +
  facet_wrap(~ gene, nrow = 1) +
  scale_x_discrete(labels = SHORT[ROUTE], drop = FALSE) +
  scale_y_nm(step = 1, top = max(st$hi), bottom = min(st$lo), minor = TRUE,
             labels = lab_signed) +
  labs(x = NULL, y = expression(atop("log"[2] ~ "difference",
                                     "vs matched control"))) +
  theme_nm(txt = TXT, txt_title = TITLE, line = LINE) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = TXT, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = TITLE, face = "italic"),
        panel.spacing.x = unit(STRIP_GAP, "mm"),
        # NO panel.widths/panel.heights here. A FACETTED plot keeps a null width in
        # its layout whatever panel.widths says, so natural_size() returns NA and
        # save_panel() is handed a missing width (2026-10-28). The facet region is
        # fixed at the stack level instead, by fixed_size() below.
        plot.margin = margin(1, 1.5, 0.5, 0, "mm"))

# =============================================================================
# THE PANELS -- saved separately, so each can be placed on its own
# =============================================================================
main_panel <- stack_parts(p_main)
m_size <- ceiling(natural_size(main_panel) * 2) / 2
save_panel(main_panel,
           fig = "EDFig5", panel = "a", name = "bbc3_routes",
           width = m_size[["width"]], height = m_size[["height"]])

strip_panel <- fixed_size(p_strip, STRIP_W * length(GENES) + STRIP_GAP * (length(GENES) - 1),
                          GRAPH_H)
s_size <- ceiling(natural_size(strip_panel) * 2) / 2
stopifnot(!anyNA(s_size))
save_panel(strip_panel,
           fig = "EDFig5", panel = "a", name = "p53_specificity",
           width = s_size[["width"]], height = s_size[["height"]])

# --- numbers for the legend (printed, never drawn) ----------------------------
cat(sprintf("\nExtended Data a -- for the legend (main %g x %g mm, strip %g x %g mm)\n",
            m_size[["width"]], m_size[["height"]], s_size[["width"]], s_size[["height"]]))
cat("  n = 7 tumours per arm; difference in medians against each arm's own matched\n")
cat("  control, with a 5,000-resample bootstrap interval; filled where the interval\n")
cat("  excludes zero. Cdkn2a is not shown for the p19ARF-disrupted arm: the locus is\n")
cat("  CRISPR-targeted there.\n")
for (i in seq_len(nrow(bb))) {
  cat(sprintf("    %-22s %-13s vs %-6s  retained %8s   Bbc3 %+.3f [%+.3f, %+.3f]%s\n",
              background_labels[[as.character(bb$arm[i])]], as.character(bb$arm[i]),
              bb$control[i], dose_txt[as.character(bb$arm[i])],
              bb$diff[i], bb$lo[i], bb$hi[i], if (bb$excl[i]) "  (filled)" else ""))
}
cat("  strip, the only interval excluding zero:\n")
mv <- st[st$excl, ]
if (nrow(mv)) {
  for (i in seq_len(nrow(mv)))
    cat(sprintf("    %-7s in %-13s %+.3f [%+.3f, %+.3f]\n", as.character(mv$gene[i]),
                as.character(mv$arm[i]), mv$diff[i], mv$lo[i], mv$hi[i]))
} else cat("    none\n")

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p_main)
  print(p_strip)
  ## the empty slot: Cdkn2a has no KO arm by design, and nothing on the panel says
  ## so -- check the legend carries it before this goes out
  st[st$gene == "Cdkn2a", c("gene", "arm", "diff", "lo", "hi")]
  bb[, c("arm", "control", "diff", "lo", "hi", "excl")]
}
