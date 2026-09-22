# =============================================================================
# biogax_tf_context.R -- why the question stayed open for a whole analysis block
# -----------------------------------------------------------------------------
# DISCUSSION PANEL (see biogax_regulon_split.R for why the `biogax_` prefix).
#
# THE POINT. The obvious way to ask "which transcription factor turned the
# respiratory chain down" is to run the TF-target sets and read the ranking. In
# this contrast that does not work, and this panel is the demonstration rather
# than the assertion.
#
# The mammary-context TF sets (Gray CHEA lane) come as TF x cell-state context.
# Plotted by context, the lanes separate almost completely BY CONTEXT: the spread
# between context medians is an order of magnitude larger than the spread of
# different transcription factors within one context. The layer is reporting
# which cell state the gland is in, not which factor is active.
#
# THE FOUR NAMED POINTS ARE THE ARGUMENT. Inside AP_LE -- the context where the
# mitochondrial TF programmes were detected in the first place -- a MICOS
# structural protein (CHCHD3) and general transcription factor IIIA (GTF3A) both
# rank BELOW ERRa and GABP. Neither is a mammary transcription factor in any
# useful sense. A layer that ranks them above the biogenesis factors cannot be
# used to nominate a biogenesis factor.
#
# Reads : results/biogenesis_axis_developmental.rds (script 47 PART C)
# Output: outputs/figures/panels/biogax_tf_context.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

ba_path <- here::here("results", "biogenesis_axis_developmental.rds")
if (!file.exists(ba_path)) stop("run scripts/47_... first")
ba <- readRDS(ba_path)

al <- as.data.frame(ba$tf_all_lanes)
tl <- as.data.frame(ba$tf_layer)
stopifnot(all(c("pathway", "NES", "ctx", "is_gray", "is_mito_lane") %in% names(al)))

g <- al[al$is_gray & !is.na(al$ctx), ]
g$tf <- vapply(strsplit(sub("^TFT_", "", sub("_MITO$", "", g$pathway)), "_GRAY_",
                        fixed = TRUE), `[`, character(1), 1L)

ctx_med <- tapply(g$NES, g$ctx, stats::median)
CTX_LEV <- names(sort(ctx_med))
g$ctx  <- factor(g$ctx, levels = CTX_LEV)
# The y axis is built as a CONTINUOUS one with named breaks rather than a discrete
# scale, because the four callouts have to be placed in empty space away from
# their own row and joined to it by a leader. A discrete scale will not take the
# numeric positions that needs.
g$ypos <- as.numeric(g$ctx)

# The four named lanes, inside ONE context so the comparison is like for like.
NAMED <- c("TFT_CHCHD3_GRAY_AP_LE_MITO", "TFT_GTF3A_GRAY_AP_LE_MITO",
           "TFT_ESRRA_GRAY_AP_LE_MITO",  "TFT_GABPA_GRAY_AP_LE_MITO")
stopifnot(all(NAMED %in% g$pathway))
lab <- g[g$pathway %in% NAMED, ]
lab$short <- c(CHCHD3 = "CHCHD3, a MICOS subunit", GTF3A = "GTF3A, general TF IIIA",
               ESRRA = "ERRa", GABPA = "GABP")[lab$tf]
lab$is_axis <- lab$tf %in% c("ESRRA", "GABPA")
lab <- lab[order(lab$NES), ]

# The callout stack sits in the empty upper-left region -- rows AP_DUCT through
# HS_HE hold nothing below NES -2.2 -- and each label is joined to its own point.
AP_LE_Y <- which(CTX_LEV == "AP_LE")
stopifnot(length(AP_LE_Y) == 1L)
lab$ty <- AP_LE_Y + c(1.15, 2.15, 3.15, 4.15)
lab$tx <- min(g$NES) - 0.06

# ASSERTION: the demonstration only works if the two non-factors really do
# outrank the two biogenesis factors.
stopifnot(
  lab$NES[lab$tf == "CHCHD3"] < lab$NES[lab$tf == "ESRRA"],
  lab$NES[lab$tf == "GTF3A"]  < lab$NES[lab$tf == "ESRRA"],
  tl$adj_r2_context > tl$adj_r2_tf)

AMB <- tl$ambient_NES

set.seed(11)
p <- ggplot2::ggplot(g, ggplot2::aes(x = NES, y = ypos)) +
  ggplot2::geom_vline(xintercept = AMB, linewidth = 0.3, linetype = "22",
                      colour = "grey45") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey80") +
  ggplot2::geom_jitter(height = 0.2, width = 0, size = 0.62, stroke = 0,
                       colour = "grey62", alpha = 0.8) +
  ggplot2::geom_segment(data = lab,
    ggplot2::aes(x = tx + 0.02, xend = NES, y = ty, yend = ypos),
    inherit.aes = FALSE, linewidth = 0.18, colour = "grey60") +
  ggplot2::geom_point(data = lab, ggplot2::aes(x = NES, y = ypos, fill = is_axis),
                      inherit.aes = FALSE, size = 1.5, shape = 21, stroke = 0.3,
                      colour = "white") +
  ggplot2::geom_text(data = lab,
    ggplot2::aes(x = tx, y = ty, label = short, colour = is_axis),
    inherit.aes = FALSE, size = 1.8, hjust = 0, vjust = 0.5) +
  ggplot2::scale_fill_manual(
    values = c(`TRUE` = unname(verdict_cols[["withdraws"]]), `FALSE` = "grey20"),
    guide = "none") +
  ggplot2::scale_colour_manual(
    values = c(`TRUE` = unname(verdict_cols[["withdraws"]]), `FALSE` = "grey20"),
    guide = "none") +
  ggplot2::annotate("text", x = AMB + 0.07, y = length(CTX_LEV) + 0.75,
                    label = "median of all 421 TF lanes",
                    size = 1.75, hjust = 0, colour = "grey45") +
  ggplot2::scale_x_continuous(name = "fGSEA NES, wild-type 6 to 12 weeks",
                              breaks = seq(-3, 2, 1), labels = lab_signed) +
  ggplot2::scale_y_continuous(name = "Gray cell-state context",
                              breaks = seq_along(CTX_LEV), labels = CTX_LEV,
                              limits = c(0.4, length(CTX_LEV) + 1.1),
                              expand = c(0, 0)) +
  theme_panel() +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

LEGEND <- panel_legend(
  slot = "Discussion D3",
  what = paste(
    "Every mammary-context TF-target lane in the wild-type 6-to-12-week contrast,",
    "one point per lane, arranged by the cell-state context the lane was defined",
    "in. Dashed line is the median of all 421 TF lanes in this contrast."),
  detail = c(
    sprintf("Context explains %.1f%% of the variance in lane score (adjusted R2) and the identity of the transcription factor %.1f%%, over %d lanes covering %d factors in %d contexts.",
            100 * tl$adj_r2_context, 100 * tl$adj_r2_tf, tl$n_lanes_gray, tl$n_tf, tl$n_context),
    sprintf("The spread BETWEEN context medians is %.2f NES; within AP_LE the interquartile spread across 87 different factors is about 0.50.",
            tl$ctx_median_range),
    sprintf("Inside AP_LE, ERRa ranks 21st and GABP 25th of 87 lanes. CHCHD3 (%.2f) and GTF3A (%.2f) both rank above them, at 6th and 9th.",
            lab$NES[lab$tf == "CHCHD3"], lab$NES[lab$tf == "GTF3A"]),
    sprintf("The whole layer is negative-going in this contrast: ambient NES %.3f, so a lane must be read against that and never against zero. Script 24's PGC1a-axis lane median is -1.14, i.e. ABOVE ambient.",
            AMB)),
  bounds = c(
    "`_MITO` lanes are built as TF programme INTERSECTED WITH MitoCarta, so a lane with an OXPHOS-heavy overlap falls by construction. They are drawn here as the demonstration and are never read as factor activity.",
    "One lane does clear its own context: ERRa's AP_LE programme, by 0.54 NES. It still fails the content filter -- that programme is exactly where the ERRa-mito overlap was detected (47 genes, p = 1.1e-22) -- so both filters, not one, have to pass.",
    "Gray CHEA contexts are enrichment of a TF's ChEA3 targets in a cell-state signature, not binding measured in these mice."),
  source = "results/biogenesis_axis_developmental.rds (script 47 PART C)")

save_panel_p(p, "biogax_tf_context", height = 72)

if (FALSE) {
  print(p)
  ba$tf_layer |> print()
  ba$tf_context |> print(n = 20)
  ba$tf_demo |> head(12) |> print()
}
