# =============================================================================
# edfig5b_ox_rel_routes.R -- Extended Data, the same route comparison on the
#                            second within-compartment ruler
# -----------------------------------------------------------------------------
# Fig. 3i draws OXPHOS priority (mitoPPS). This panel draws `ox_rel` over the
# same four arms in the same order, as a SECOND RULER: two different constructions
# of "OXPHOS relative to the rest of the mitochondrial compartment", agreeing on
# every arm (author, 2026-10-28).
#
# WHAT ox_rel IS, because the name invites the wrong reading. It is NOT OXPHOS
#   against the nuclear transcriptome. Script 52:
#       rest_mito <- setdiff(mito_nuc, ox_sub)              # 880 genes
#       ox_rel    <- comp_e(ox_sub) - comp_e(rest_mito)     # 86 subunits minus those
#   so the denominator is the REST OF THE NUCLEAR MitoCarta COMPARTMENT, and the
#   quantity is a difference of mean per-gene z-scores, not a proportion. It is
#   therefore already normalised within the compartment, the same family as
#   mitoPPS -- which is why the two agree, and why neither reads as "overall
#   mitochondrial level". The quantity that does read that way is the MitoCarta
#   share of the library (14.7% in Pgc1a_BclxL against 7.8-8.6% elsewhere), which
#   is in the saved object as `mito_share` and is not drawn here.
#
# ORDER, COLOUR, FORMAT and the open hyphen issue: all as Fig. 3i. Read that
#   script's header first; nothing is repeated here.
#
# THE FIGURE NUMBER IS PROVISIONAL, as for EDFig5a: the spec carries Extended
#   Data Fig. 5 as a placeholder (section 4, item 8). When the numbering settles
#   this is a `git mv` and a one-line change to the save_panel() call.
#
# NO ANALYSIS AND NO CAPTION. Every value is re-read from the saved object.
#
# Reads : results/orthotopic_section3_panels.rds  (script 57, which reads script 52)
# Output: outputs/natmetab/EDFig5/EDFig5b_ox_rel_routes.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

o  <- readRDS(here::here("results", "orthotopic_section3_panels.rds"))
sn <- as.data.frame(o$summary_numbers)

ROUTE <- c("BclxL", "KOPgc1a", "NTPgc1a", "Pgc1a_BclxL")     # the escalation, as Fig. 3i
stopifnot(setequal(ROUTE, sn$arm), all(ROUTE %in% names(background_cols)))

CONSTRUCT <- c(BclxL       = 'paste("Bcl-xL")',
               KOPgc1a     = 'paste("PGC-1", alpha, ", p19ARF")',
               NTPgc1a     = 'paste("PGC-1", alpha, ", NT")',
               Pgc1a_BclxL = 'paste("PGC-1", alpha, " + Bcl-xL")')
x_lab <- parse(text = sprintf('atop("%s", %s)',
                              background_labels[ROUTE], CONSTRUCT[ROUTE]))

GRAPH_H <- 20
CAT_W   <- 20
LINE    <- lw_pt(NM_LINE_PRISM[["Fig3"]])
TXT     <- 5
TITLE   <- 6
SYM     <- 1.6
STROKE  <- NM_LINE_PRISM[["Fig3"]] * 96 / 72 / (ggplot2::.stroke / 2)

d <- sn[match(ROUTE, sn$arm), ]
d$arm  <- factor(d$arm, levels = ROUTE)
d$excl <- d$ox_rel_hi < 0 | d$ox_rel_lo > 0
d$fill_col <- ifelse(d$excl, background_cols[as.character(d$arm)], "white")
stopifnot(nrow(d) == 4L, !anyNA(d$ox_rel_diff))

# The axis title is short because it is rotated: at 6 pt about 19 characters fit
# in a 20 mm panel. The full definition goes in the legend, printed below.
p <- ggplot(d, aes(arm, ox_rel_diff)) +
  geom_hline(yintercept = 0, linewidth = LINE, colour = box_line) +
  geom_errorbar(aes(ymin = ox_rel_lo, ymax = ox_rel_hi), width = 0.18, linewidth = LINE,
                colour = "black") +
  geom_point(aes(fill = fill_col), shape = 21, colour = "black",
             size = SYM, stroke = STROKE) +
  scale_fill_identity() +
  scale_x_discrete(labels = x_lab) +
  scale_y_nm(step = 0.5, top = max(d$ox_rel_hi), bottom = min(d$ox_rel_lo),
             minor = TRUE, labels = lab_signed) +
  labs(x = NULL, y = expression(atop("Relative OXPHOS (z)", "vs matched control"))) +
  theme_nm(txt = TXT, txt_title = TITLE, line = LINE) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = TXT, lineheight = 0.95, vjust = 1),
        panel.heights = unit(GRAPH_H, "mm"),
        panel.widths  = unit(CAT_W * length(ROUTE), "mm"),
        plot.margin = margin(1, 1.5, 0.5, 0, "mm"))

panel <- stack_parts(p)
size  <- ceiling(natural_size(panel) * 2) / 2
stopifnot(!anyNA(size))
save_panel(panel, fig = "EDFig5", panel = "b", name = "ox_rel_routes",
           width = size[["width"]], height = size[["height"]])

# --- numbers for the legend (printed, never drawn) ----------------------------
cat(sprintf("\nExtended Data b -- for the legend (data area %g x %g mm, saved %g x %g mm)\n",
            CAT_W * length(ROUTE), GRAPH_H, size[["width"]], size[["height"]]))
cat("  n = 7 tumours per arm. The 86 nuclear OXPHOS subunits against the other 880\n")
cat("  nuclear MitoCarta genes, as a difference of mean per-gene z-scores; difference\n")
cat("  in medians against each arm's own matched control, with a 5,000-resample\n")
cat("  bootstrap interval; filled where the interval excludes zero. The same four\n")
cat("  arms and order as Fig. 3i, on a second within-compartment ruler.\n")
for (i in seq_len(nrow(d))) {
  cat(sprintf("    %-22s %-13s vs %-6s  ox_rel %+.3f [%+.3f, %+.3f]%s\n",
              background_labels[[as.character(d$arm[i])]], as.character(d$arm[i]),
              d$control[i], d$ox_rel_diff[i], d$ox_rel_lo[i], d$ox_rel_hi[i],
              if (d$excl[i]) "  (filled)" else ""))
}
cat("  mitoPPS agrees on every arm (Fig. 3i): same two exclude zero, same two span it.\n")

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p)
  d[, c("arm", "control", "ox_rel_diff", "ox_rel_lo", "ox_rel_hi", "excl")]
  ## the mass readout this panel is NOT: MitoCarta share of the library, per arm
  as.data.frame(o$summary_numbers)[, c("arm", "mito_share_arm", "mito_share_ctrl")]
}
