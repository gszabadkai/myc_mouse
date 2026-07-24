# =============================================================================
# fig03_background_vs_myc.R -- is Myc's mitochondrial effect shaped by a changing
# background?  Four panels, one answer.
# -----------------------------------------------------------------------------
# The genotype figures (fig01/fig02) read Myc WITHIN a timepoint; the temporal ones
# (figS6/figS7) read each genotype ACROSS the timeline. This figure confronts them.
#
#   A  STATE SPACE. Each Level-1 MitoPathway tier as four states in (priority,
#      content). Red arrows = the Myc effect at each age (6W solid, 12W dashed);
#      grey arrows = the temporal move of each genotype. The Myc arrow keeps its
#      direction and shortens; the background arrow points somewhere else.
#   B  RESCALED, NOT RESHAPED. The 12W Myc effect IS the 6W effect times ~0.55 --
#      the same pathways, a smaller amplitude. Its coherence (R2) is at the 100th
#      percentile of expression-matched shuffled sets.
#   C  SHARED VECTOR + MYC OFFSET. Myc+ 6->12W against WT 6->12W. The SLOPE (the
#      shared temporal component) is NOT beyond the null; the INTERCEPT is, at the
#      0th percentile -- the Myc-specific part is a uniform downward offset, not a
#      different set of pathways moving.
#   D  DOES WT-CONVERGENCE SURVIVE ITS CONTROL? Issue #6's genotype contrast and WT
#      temporal contrast share the 6W WT baseline, which manufactures convergence.
#      Splitting the six 6W_neg mice (open -> filled) collapses the GLOBAL
#      convergence to chance; only the biosynthetic arms survive, and OXPHOS /
#      MYC-core divergence survives.
#
# BATCH = TIMEPOINT: both temporal contrasts carry the same batch offset. It cancels
# in the interaction (panel C's intercept) but NOT in panel D's split, and
# "developmental" stays an interpretation of the shared vector. Exploratory (n=6).
#
# Reads (read-only; author runs scripts/40_background_vs_myc_decomposition.R first):
#   results/background_vs_myc.rds -- $state_table (A), $ruler + $regressions +
#       $regression_null (B, C), $split_summary (D), $defs$tier_levels.
# =============================================================================

source(here::here("figures", "theme_myc.R"))
if (!requireNamespace("patchwork", quietly = TRUE)) stop("fig03 needs patchwork")
if (!requireNamespace("ggrepel", quietly = TRUE))   stop("fig03 needs ggrepel")

out_dir <- here::here("outputs", "figures")

bg <- readRDS(here::here("results", "background_vs_myc.rds"))

# --- shared tier labels + palette (identical to fig02 / figS5 / figS7) --------
tier_lv  <- bg$defs$tier_levels
tier_lab <- c("Protein import / homeostasis", "Central dogma", "OXPHOS", "Metabolism",
              "Signaling", "Dynamics & surveillance", "SM transport")
names(tier_lab) <- tier_lv
tier_col <- c("Protein import / homeostasis" = "#D55E00", "Central dogma" = "#E69F00",
              "OXPHOS"                  = "#009E73", "Metabolism"    = "grey72",
              "Signaling"               = "#56B4E9",
              "Dynamics & surveillance" = "#0072B2", "SM transport"  = "#CC79A7")

# =============================================================================
# PANEL A -- state space: four states per tier, with the Myc and the time arrows
# =============================================================================
# short strip labels: the facets are narrow, the full tier names clip
tier_short <- c("Protein import", "Central dogma", "OXPHOS", "Metabolism",
                "Signaling", "Dynamics", "SM transport")
names(tier_short) <- tier_lv

st <- bg$state_table
st$Tier <- factor(unname(tier_short[st$tier]), levels = unname(tier_short[tier_lv]))
wide <- stats::reshape(
  as.data.frame(st[, c("Tier", "group", "priority", "content_share")]),
  idvar = "Tier", timevar = "group", direction = "wide")
names(wide) <- sub("^priority\\.", "x_", sub("^content_share\\.", "y_", names(wide)))

arr <- function(g1, g2, kind) data.frame(
  Tier = wide$Tier, kind = kind,
  x = wide[[paste0("x_", g1)]], y = wide[[paste0("y_", g1)]],
  xend = wide[[paste0("x_", g2)]], yend = wide[[paste0("y_", g2)]])
arrows_df <- rbind(
  arr("6W_neg",  "6W_pos",  "Myc effect, 6W"),
  arr("12W_neg", "12W_pos", "Myc effect, 12W"),
  arr("6W_neg",  "12W_neg", "time, WT"),
  arr("6W_pos",  "12W_pos", "time, Myc+"))
arrows_df$kind <- factor(arrows_df$kind,
  levels = c("Myc effect, 6W", "Myc effect, 12W", "time, WT", "time, Myc+"))

arrow_col <- c("Myc effect, 6W" = "#D73027", "Myc effect, 12W" = "#D73027",
               "time, WT" = "grey45", "time, Myc+" = "grey45")
arrow_lty <- c("Myc effect, 6W" = 1, "Myc effect, 12W" = 2,
               "time, WT" = 1, "time, Myc+" = 2)

pA <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = arrows_df,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                 colour = kind, linetype = kind), linewidth = 0.4,
    arrow = ggplot2::arrow(length = ggplot2::unit(1.6, "mm"), type = "closed")) +
  ggplot2::geom_point(data = st,
                      ggplot2::aes(priority, content_share, fill = group),
                      shape = 21, size = 1.9, colour = "grey20", stroke = 0.3) +
  ggplot2::facet_wrap(~ Tier, nrow = 2, scales = "free") +
  ggplot2::scale_colour_manual(values = arrow_col, name = NULL) +
  ggplot2::scale_linetype_manual(values = arrow_lty, name = NULL) +
  ggplot2::scale_fill_manual(values = group_cols, labels = group_labels, name = NULL) +
  ggplot2::labs(
    x = "Priority  (mitoPPS tier score; 1 = compartment average)",
    y = "Content  (% of the nuclear transcriptome)",
    title = "A   The Myc arrow keeps its direction and shortens; the background arrow points elsewhere") +
  theme_myc(base_size = 8) +
  ggplot2::theme(
    legend.position  = "right",
    legend.spacing.y = ggplot2::unit(0.5, "mm"),
    legend.key.size  = ggplot2::unit(3.2, "mm"),
    legend.text      = ggplot2::element_text(size = 6.2),
    strip.text       = ggplot2::element_text(size = 6.8, face = "bold"),
    axis.text        = ggplot2::element_text(size = 5.8),
    plot.title       = ggplot2::element_text(face = "bold", size = 9)) +
  # colour and linetype must share name+labels+order or ggplot draws the arrow key twice
  ggplot2::guides(
    colour   = ggplot2::guide_legend(order = 1),
    linetype = ggplot2::guide_legend(order = 1),
    fill     = ggplot2::guide_legend(order = 2, override.aes = list(size = 2)))

# =============================================================================
# PANELS B and C -- the two regressions (CONTENT ruler; mtDNA pathway excluded)
# =============================================================================
d   <- bg$ruler[!bg$ruler$is_mtdna, ]
d$Tier <- factor(unname(tier_lab[d$tier]), levels = unname(tier_lab[tier_lv]))
reg <- bg$regressions
nul <- bg$regression_null
gv  <- function(model, scope, field) reg[[field]][grepl(model, reg$model) & reg$scope == scope]
nv  <- function(stat, field) nul[[field]][nul$statistic == stat]

# The 7 Level-1 tiers are also pathways in their own right -> ringed anchor points.
# No text labels here: the tier colour legend already names them and repel labels
# swamp 143 points in a 60 mm panel.
anc <- d[d$pathway %in% tier_lv, ]

scatter_base <- function(dat, xv, yv, xlab, ylab, ttl, note, slope, intercept) {
  rng <- range(c(dat[[xv]], dat[[yv]]), na.rm = TRUE)
  pad <- diff(rng) * 0.04
  ggplot2::ggplot(dat, ggplot2::aes(.data[[xv]], .data[[yv]])) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey55",
                         linetype = 2, linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(colour = Tier), size = 1.3, alpha = 0.85) +
    ggplot2::geom_point(data = anc, ggplot2::aes(fill = Tier), shape = 21, size = 2.4,
                        colour = "black", stroke = 0.4, show.legend = FALSE) +
    ggplot2::geom_abline(slope = slope, intercept = intercept,
                         colour = "black", linewidth = 0.5) +
    ggplot2::annotate("text", x = rng[1], y = rng[2], hjust = 0, vjust = 1,
                      size = 1.9, colour = "grey20", lineheight = 1.2, label = note) +
    ggplot2::scale_colour_manual(values = tier_col, name = NULL, drop = FALSE) +
    ggplot2::scale_fill_manual(values = tier_col, guide = "none") +
    ggplot2::coord_cartesian(xlim = rng + c(-pad, pad), ylim = rng + c(-pad, pad)) +
    ggplot2::labs(x = xlab, y = ylab, title = ttl) +
    theme_myc(base_size = 8) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 8.5))
}

pB <- scatter_base(
  d, "c_m6", "c_m12",
  "Myc effect at 6W  (set log2FC)", "Myc effect at 12W  (set log2FC)",
  "B   Rescaled, not reshaped",
  sprintf("slope %.2f (boot %.2f-%.2f)\nR2 %.2f, rho %.2f\nR2: pct %.0f of the null\n\nsame pathways,\n%.0f%% of the amplitude",
          gv("^rescale", "content, all", "slope"),
          bg$regression_boot$lo[bg$regression_boot$model == "rescale content"],
          bg$regression_boot$hi[bg$regression_boot$model == "rescale content"],
          gv("^rescale", "content, all", "r2"), gv("^rescale", "content, all", "rho"),
          nv("rescale_r2", "percentile"),
          100 * gv("^rescale", "content, all", "slope")),
  gv("^rescale", "content, all", "slope"), gv("^rescale", "content, all", "intercept"))

pC <- scatter_base(
  d, "c_tn", "c_tp",
  "WT 6->12W  (set log2FC)", "Myc+ 6->12W  (set log2FC)",
  "C   Shared vector, Myc-specific offset",
  sprintf("slope %.2f = shared move\n  null %.2f, pct %.0f\n  (not distinguishable)\nintercept %.2f = Myc-specific\n  null %.2f, pct %.0f\n\nthe intercept is the\nclean interaction",
          gv("^shared", "content, all", "slope"), nv("shared_slope", "null_median"),
          nv("shared_slope", "percentile"),
          gv("^shared", "content, all", "intercept"), nv("shared_int", "null_median"),
          nv("shared_int", "percentile")),
  gv("^shared", "content, all", "slope"), gv("^shared", "content, all", "intercept"))

# =============================================================================
# PANEL D -- does WT-convergence survive the sample-split control?
# =============================================================================
ss <- bg$split_summary
ss <- ss[ss$program != "ALL" | TRUE, ]
lab_of <- c(ALL = "ALL genes (global)",
            MITOCARTA_AMINO_ACID_METABOLISM   = "Amino acid metab.",
            MITOCARTA_LIPID_METABOLISM        = "Lipid metab.",
            MITOCARTA_TCA_CYCLE               = "TCA cycle",
            MITOCARTA_NUCLEOTIDE_METABOLISM   = "Nucleotide metab.",
            MITOCARTA_MITOCHONDRIAL_RIBOSOME  = "Mitoribosome",
            MITOCARTA_OXPHOS                  = "OXPHOS",
            MITOCARTA_OXPHOS_SUBUNITS         = "OXPHOS subunits",
            MYC_HALLMARK_MYC_TARGETS_V2       = "MYC targets (Hallmark)",
            MYC_felsher_integrative_signature = "MYC signature (Felsher)",
            PROLIFERATION_pooled              = "Proliferation",
            MAMMARY_LUMINAL_pooled            = "Mammary luminal")
ss$lab <- unname(lab_of[ss$program])
ss <- ss[!is.na(ss$lab), ]
ss <- ss[order(ss$frac_toward_split), ]
ss$lab <- factor(ss$lab, levels = ss$lab)
ss$verdict <- ifelse(ss$frac_toward_split > 0.55, "converges",
                     ifelse(ss$frac_toward_split < 0.45, "diverges", "chance"))
verdict_col <- c(converges = "#1B7837", chance = "grey55", diverges = "#762A83")

pD <- ggplot2::ggplot(ss, ggplot2::aes(y = lab)) +
  ggplot2::annotate("rect", xmin = 0.45, xmax = 0.55, ymin = -Inf, ymax = Inf,
                    fill = "grey92") +
  ggplot2::geom_vline(xintercept = 0.5, colour = "grey40", linewidth = 0.35, linetype = 2) +
  ggplot2::geom_segment(ggplot2::aes(x = frac_toward_shared, xend = frac_toward_split,
                                     yend = lab), colour = "grey55", linewidth = 0.35,
                        arrow = ggplot2::arrow(length = ggplot2::unit(1.3, "mm"),
                                               type = "closed")) +
  ggplot2::geom_point(ggplot2::aes(x = frac_toward_shared), shape = 21, size = 1.9,
                      fill = "white", colour = "grey40", stroke = 0.4) +
  ggplot2::geom_point(ggplot2::aes(x = frac_toward_split, colour = verdict), size = 2.2,
                      show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = verdict_col, guide = "none") +
  ggplot2::annotate("text", x = 0.265, y = nrow(ss) + 1.1, label = "WT moves away",
                    size = 1.9, colour = verdict_col[["diverges"]], hjust = 0, vjust = 0) +
  ggplot2::annotate("text", x = 0.785, y = nrow(ss) + 1.1, label = "WT converges",
                    size = 1.9, colour = verdict_col[["converges"]], hjust = 1, vjust = 0) +
  ggplot2::coord_cartesian(xlim = c(0.25, 0.80), ylim = c(0.4, nrow(ss) + 2.1),
                           expand = FALSE) +
  ggplot2::labs(
    x = "fraction of genes whose WT change points toward the Myc state",
    y = NULL, title = "D   WT-convergence vs its control",
    subtitle = "open = shared 6W WT baseline (as published)\nfilled = baseline split, 20 three-vs-three partitions") +
  theme_myc(base_size = 8) +
  ggplot2::theme(axis.text.y   = ggplot2::element_text(size = 6.4),
                 plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                       lineheight = 1.15))

# =============================================================================
# ASSEMBLY
# =============================================================================
bottom <- patchwork::wrap_plots(pB, pC, pD, nrow = 1, widths = c(1, 1, 1.25)) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom",
                 legend.key.size = ggplot2::unit(3, "mm"),
                 legend.text     = ggplot2::element_text(size = 6.2))

p <- patchwork::wrap_plots(pA, bottom, ncol = 1, heights = c(1, 1.12)) +
  patchwork::plot_annotation(
    caption = paste(
      "Content = DESeq2 set-average RAW log2FC; priority = mitoPPS (Monzel 2025), content-blind; share excludes the 13 mtDNA-encoded genes from numerator and denominator.",
      "B/C exclude the synthetic mtDNA-encoded pathway (a 4-6.5 SD outlier on every temporal contrast). Nulls = within-decile gene-label shuffles preserving set size, expression and pathway overlap.",
      "BATCH = TIMEPOINT: the batch offset cancels in the interaction (panel C's intercept) but not in panel D's split, so 'developmental' is an interpretation of the shared vector. Exploratory, n=6/group.",
      sep = "\n"),
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(size = 5.8, hjust = 0, colour = "grey30",
                                           lineheight = 1.15)))

if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "fig03_background_vs_myc.pdf"),
             width = fig_w[["double"]], height = 175)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(wide)
  bg$regressions |> dplyr::filter(scope %in% c("content, all", "priority, all")) |> print()
  bg$regression_null |> print()
  ss[, c("lab", "frac_toward_shared", "frac_toward_split", "verdict")] |> print(n = 20)
  print(pA); print(pB); print(pC); print(pD)
  print(p)
}
