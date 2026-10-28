# =============================================================================
# fig3i_route_comparison.R -- Fig. 3i, the route comparison: TWO panels
# -----------------------------------------------------------------------------
# Two plots placed side by side in Illustrator, both panel i (author, 2026-10-28):
#
#   Fig3i_retained_dose.pdf   Ppargc1a RETAINED, per animal, CPM. The SPREAD is
#                             the point, so every animal is drawn over a median
#                             line: NTPgc1a runs 52 to 272 CPM against
#                             Pgc1a_BclxL's 242 to 364.
#   Fig3i_mitopps.pdf         OXPHOS priority (mitoPPS), difference against each
#                             arm's own MATCHED control, with script 52's
#                             5,000-resample bootstrap interval. Filled where the
#                             interval excludes zero, open where it does not.
#
# THE ORDER IS AN ESCALATION (author, 2026-10-28): no transgene, checkpoint
#   disrupted, pathway intact, effector blocked -- what has to be done to the
#   death pathway before a tumour can hold a high respiratory state. This
#   REPLACES the order in the spec's drawing requirement 1, which ran the other
#   way (effector blocked first); the spec and FIGURE_RULES are updated to match.
#   It is deliberately NOT dose order: the dose panel beside it carries that, and
#   the two together are the point.
#
# WHY mitoPPS AND NOT ox_rel (author, 2026-10-28). Both are WITHIN-COMPARTMENT
#   rulers -- ox_rel is the 86 OXPHOS subunits against the other 880 nuclear
#   MitoCarta genes, not against the transcriptome -- and they agree on every arm:
#   the same two exclude zero, the same two span it. Nothing in the conclusion
#   turns on the choice. mitoPPS is on the main figure because Fig. 1e and 1f
#   already use it ("Priority (mitoPPS)"), and the section's argument is that the
#   two systems give the same answer when buffer status is matched: one ruler
#   across both makes that legible. No number is set beside a Fig. 1 number, so
#   the scoring-run rule is untouched.
#   ox_rel goes to Extended Data as a second ruler (edfig5b_ox_rel_routes.R).
#   ox_lvl, the absolute level, is NOT used: pathway intact spans zero on it,
#   which is the escape note's "tilt, not mass expansion" reading and would put a
#   null under the arm whose sentence the section depends on.
#
# "RETAINED", NEVER "ASSIGNED". Retention is an outcome of selection, not a
#   designed variable. Nobody titrated PGC-1a.
#
# NO FITTED LINE ANYWHERE, AND NO DOSE AXIS. The respiratory ruler does not track
#   retained dose WITHIN any arm (n = 7; Pgc1a_BclxL r -0.046, NTPgc1a +0.172,
#   KOPgc1a +0.729, BclxL +0.465, every interval spanning zero, 2026-10-28), so
#   the relation is readable only ACROSS arms and the x axis stays categorical.
#
# COLOUR ENCODES THE BACKGROUND, not the construct (FIGURE_RULES.md, "The
#   cell-line constructs"; `background_cols`). All three informative arms carry
#   PGC-1a, so the construct palette would give two of them one hue.
#
# PRISM FORMAT, to sit beside panels a-h: data areas 20 mm tall, axis title 6 pt,
#   tick labels 5 pt, every line 0.72 pt -- the weight MEASURED from the built
#   figure (NM_LINE_PRISM[["Fig3"]]).
#
# NO ANALYSIS AND NO CAPTION. Every value is re-read from the saved object;
#   nothing is re-fitted and no interval is recomputed. Legends describe panels
#   only, so no result and no interpretation is drawn.
#
# Reads : results/orthotopic_section3_panels.rds  (script 57, which reads script 52)
# Output: outputs/natmetab/Fig3/Fig3i_retained_dose.pdf
#         outputs/natmetab/Fig3/Fig3i_mitopps.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

o  <- readRDS(here::here("results", "orthotopic_section3_panels.rds"))
pa <- as.data.frame(o$per_animal)
sn <- as.data.frame(o$summary_numbers)

# --- the escalation order, fixed here and used by both panels -----------------
ROUTE <- c("BclxL", "KOPgc1a", "NTPgc1a", "Pgc1a_BclxL")
stopifnot(setequal(ROUTE, sn$arm), all(ROUTE %in% names(background_cols)))

# The construct beneath each route. Plotmath, because base pdf() has no Greek in
# the text font but reaches the alpha through the Symbol font.
# KNOWN AND OPEN: base pdf() draws the ASCII hyphen as a long dash, so "Bcl-xL"
# prints as a dash (checked under Quick Look, not a preview artefact; the trap is
# recorded in FIGURE_RULES, "Labels avoid hyphens"). Panels a-h carry
# Illustrator-set labels with true hyphens. Either correct these four in
# Illustrator after placing, or drop this line and set it there, as Fig. 1b does.
CONSTRUCT <- c(BclxL       = 'paste("Bcl-xL")',
               KOPgc1a     = 'paste("PGC-1", alpha, ", p19ARF")',
               NTPgc1a     = 'paste("PGC-1", alpha, ", NT")',
               Pgc1a_BclxL = 'paste("PGC-1", alpha, " + Bcl-xL")')
x_lab <- parse(text = sprintf('atop("%s", %s)',
                              background_labels[ROUTE], CONSTRUCT[ROUTE]))

# --- the Prism format ----------------------------------------------------------
GRAPH_H <- 20                                   # mm, each data area
CAT_W   <- 20                                   # mm per route; "checkpoint disrupted" at 5 pt needs ~18
LINE    <- lw_pt(NM_LINE_PRISM[["Fig3"]])       # 0.72 pt, measured from the built figure
TXT     <- 5
TITLE   <- 6
# A ROTATED AXIS TITLE'S LINE LENGTH RUNS ALONG THE PANEL HEIGHT. At 6 pt about 19
# characters fit in 20 mm, so every title is two short lines set with atop().
SYM     <- 1.3
STROKE  <- NM_LINE_PRISM[["Fig3"]] * 96 / 72 / (ggplot2::.stroke / 2)

# the shared x treatment: no axis line, two-line labels drawn once per panel
x_axis <- list(
  scale_x_discrete(labels = x_lab),
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = TXT, lineheight = 0.95, vjust = 1),
        panel.heights = unit(GRAPH_H, "mm"),
        panel.widths  = unit(CAT_W * length(ROUTE), "mm"),
        plot.margin = margin(1, 1.5, 0.5, 0, "mm")))

# =============================================================================
# PANEL 1 -- Ppargc1a retained, per animal
# -----------------------------------------------------------------------------
# Linear y to 400 (author, 2026-10-28). The top animal is 364.0 CPM, so an axis
# to 320 would cut two of Pgc1a_BclxL's seven.
# =============================================================================
up <- pa[pa$arm %in% ROUTE, c("animal_id", "arm", "Ppargc1a_cpm")]
up$arm <- factor(up$arm, levels = ROUTE)
stopifnot(nrow(up) == 28L, all(table(up$arm) == 7L))

TOP <- 400
stopifnot(max(up$Ppargc1a_cpm) < TOP)           # nothing is clipped

p_dose <- ggplot(up, aes(arm, Ppargc1a_cpm)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "errorbar", width = 0.5, linewidth = LINE, colour = "black") +
  geom_point(aes(fill = arm), shape = 21, colour = "black", size = SYM, stroke = STROKE,
             position = position_jitter(width = 0.13, height = 0, seed = 1)) +
  scale_fill_manual(values = background_cols) +
  scale_y_nm(step = 100, top = TOP) +
  labs(x = NULL, y = expression(atop(italic(Ppargc1a) ~ "retained", "(CPM)"))) +
  theme_nm(txt = TXT, txt_title = TITLE, line = LINE) +
  x_axis

# =============================================================================
# PANEL 2 -- OXPHOS priority (mitoPPS) against the matched control
# -----------------------------------------------------------------------------
# A signed difference, so the axis spans zero and carries a reference line there.
# The intervals are script 52's, read as saved and never recomputed.
# Open where the interval includes zero. The fill is set per row and mapped with
# scale_fill_identity(), NOT with alpha: ggplot applies alpha to a shape-21
# point's OUTLINE as well as its fill, so alpha = 0 renders the whole symbol
# transparent and an open point disappears (checked, 2026-10-28).
# =============================================================================
lo <- sn[match(ROUTE, sn$arm), ]
lo$arm  <- factor(lo$arm, levels = ROUTE)
lo$excl <- lo$mitopps_hi < 0 | lo$mitopps_lo > 0
lo$fill_col <- ifelse(lo$excl, background_cols[as.character(lo$arm)], "white")
stopifnot(nrow(lo) == 4L, !anyNA(lo$mitopps_diff))

p_pps <- ggplot(lo, aes(arm, mitopps_diff)) +
  geom_hline(yintercept = 0, linewidth = LINE, colour = box_line) +
  geom_errorbar(aes(ymin = mitopps_lo, ymax = mitopps_hi), width = 0.18, linewidth = LINE,
                colour = "black") +
  geom_point(aes(fill = fill_col), shape = 21, colour = "black",
             size = SYM + 0.3, stroke = STROKE) +
  scale_fill_identity() +
  scale_y_nm(step = 0.2, top = max(lo$mitopps_hi), bottom = min(lo$mitopps_lo),
             minor = TRUE, labels = lab_signed) +
  labs(x = NULL, y = expression(atop("Priority (mitoPPS)", "vs matched control"))) +
  theme_nm(txt = TXT, txt_title = TITLE, line = LINE) +
  x_axis

# =============================================================================
# THE TWO PANELS -- saved separately, both letter i, placed side by side
# =============================================================================
emit <- function(p, name) {
  st <- stack_parts(p)
  sz <- ceiling(natural_size(st) * 2) / 2
  stopifnot(!anyNA(sz))
  save_panel(st, fig = "Fig3", panel = "i", name = name,
             width = sz[["width"]], height = sz[["height"]])
  sz
}
sz_dose <- emit(p_dose, "retained_dose")
sz_pps  <- emit(p_pps,  "mitopps")

# --- numbers for the legend (printed, never drawn) ----------------------------
cat(sprintf("\nFig. 3i -- for the legend (data areas %g x %g mm; dose %g x %g mm, mitoPPS %g x %g mm)\n",
            CAT_W * length(ROUTE), GRAPH_H, sz_dose[["width"]], sz_dose[["height"]],
            sz_pps[["width"]], sz_pps[["height"]]))
cat("  n = 7 tumours per arm. Left: points are animals, the line is the median.\n")
cat("  Right: OXPHOS priority (mitoPPS, Monzel 2025), difference in medians against\n")
cat("  each arm's own matched control, with a 5,000-resample bootstrap interval;\n")
cat("  filled where the interval excludes zero.\n")
for (i in seq_len(nrow(lo))) {
  cat(sprintf("    %-22s %-13s vs %-6s  retained %8.4g CPM   mitoPPS %+.3f [%+.3f, %+.3f]%s\n",
              background_labels[[as.character(lo$arm[i])]], as.character(lo$arm[i]),
              lo$control[i], lo$Ppargc1a_cpm_arm[i],
              lo$mitopps_diff[i], lo$mitopps_lo[i], lo$mitopps_hi[i],
              if (lo$excl[i]) "  (filled)" else ""))
}
cat(sprintf("  retained dose per animal: %s\n",
            paste(sprintf("%s %.4g-%.4g", ROUTE,
                          tapply(up$Ppargc1a_cpm, up$arm, min),
                          tapply(up$Ppargc1a_cpm, up$arm, max)), collapse = "; ")))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p_dose)
  print(p_pps)
  ## the open item: the hyphens in the construct sub-labels print as long dashes
  text_mm(background_labels, txt = 5)          # does "checkpoint disrupted" fit 20 mm?
  lo[, c("arm", "control", "mitopps_diff", "mitopps_lo", "mitopps_hi", "excl")]
  tapply(up$Ppargc1a_cpm, up$arm, summary)
}
