# =============================================================================
# fig02b_reallocation_ranked.R -- compact horizontal ranked dot chart
# -----------------------------------------------------------------------------
# The compact alternative to fig02_reallocation.R (which labels all ~144 pathways,
# tall). Same data, stripped for a MAIN panel: horizontal, NO per-pathway labels,
# colour = tier, asterisks mark padj<0.05, and only the 7 TOP-OF-HIERARCHY
# (Level-1) MitoPathways are labelled -- pushed OFF the data into a top strip
# (promoted) / bottom strip (demoted) with leader lines. Ranked DESCENDING:
# promoted (left) -> demoted (right). 12W is shown as small points joined to the
# 6W dot by a thin connector, so the 6W->12W FADE (attenuation of the
# reprioritisation) is visible without clutter.
#
# Myc mitoPPS effect (Myc+ - WT); each MitoPathway ratio-normalized to 1
# (content-blind), genotype is the clean axis. Descriptive/exploratory (n=6/group).
#
# Reads (read-only; NO re-run):
#   results/mitopps_scores.rds -- $mitopps_pairwise (diff/padj, Myc_effect_6W/_12W),
#       $pathway_levels (Level-1 = top hierarchy), $pathway_tier1_map.
# =============================================================================

source(here::here("figures", "theme_myc.R"))
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  stop("fig02b needs ggrepel for the top-hierarchy leader labels")
}

out_dir <- here::here("outputs", "figures")

mp   <- readRDS(here::here("results", "mitopps_scores.rds"))
pw   <- mp$mitopps_pairwise
tier <- mp$pathway_tier1_map

e6  <- pw[pw$contrast == "Myc_effect_6W",  c("pathway", "diff", "padj")]
e12 <- pw[pw$contrast == "Myc_effect_12W", c("pathway", "diff")]
names(e6)[2:3] <- c("eff6", "padj6"); names(e12)[2] <- "eff12"
d <- merge(e6, e12, by = "pathway")
d$tier <- unname(tier[d$pathway])
d <- d[!is.na(d$tier) & !is.na(d$eff6), ]

# --- tier labels + palette (Metabolism greyed) -------------------------------
tier_lv  <- c("Protein import, sorting and homeostasis", "Mitochondrial central dogma",
              "OXPHOS", "Metabolism", "Signaling",
              "Mitochondrial dynamics and surveillance", "Small molecule transport")
tier_lab <- c("Protein import / homeostasis", "Central dogma", "OXPHOS", "Metabolism",
              "Signaling", "Dynamics & surveillance", "SM transport")
names(tier_lab) <- tier_lv
tier_col <- c("Protein import / homeostasis" = "#D55E00", "Central dogma" = "#E69F00",
              "OXPHOS"                  = "#009E73", "Metabolism"    = "grey72",
              "Signaling"               = "#56B4E9",
              "Dynamics & surveillance" = "#0072B2", "SM transport"  = "#CC79A7")
d$tier <- factor(unname(tier_lab[d$tier]), levels = unname(tier_lab[tier_lv]))

# --- rank DESCENDING on the 6W effect (promoted left -> demoted right) --------
d <- d[order(-d$eff6), ]
d$rank <- seq_len(nrow(d))
d$sig  <- !is.na(d$padj6) & d$padj6 < 0.05

# --- the 7 top-of-hierarchy (Level-1) pathways = the anchor labels ------------
l1  <- mp$pathway_levels$Pathway[mp$pathway_levels$Level == "Pathway_Level1"]
anc <- d[d$pathway %in% l1, ]
anc$lab <- unname(tier_lab[names(tier)[match(anc$pathway, names(tier))]])
TOP <- 0.55; BOT <- -0.50                       # label strips (off the data)
anc_up <- anc[anc$eff6 > 0, ]; anc_dn <- anc[anc$eff6 <= 0, ]

p <- ggplot2::ggplot(d, ggplot2::aes(rank, eff6)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  # 6W -> 12W fade connector (thin) + small 12W point
  ggplot2::geom_segment(ggplot2::aes(xend = rank, y = eff6, yend = eff12),
                        colour = "grey80", linewidth = 0.2) +
  ggplot2::geom_point(ggplot2::aes(y = eff12, colour = tier), size = 0.5, alpha = 0.55) +
  # 6W point (the ranked value)
  ggplot2::geom_point(ggplot2::aes(colour = tier), size = 1.3, alpha = 0.9) +
  ggplot2::geom_text(data = d[d$sig, ], ggplot2::aes(y = eff6 + 0.02),
                     label = "*", size = 2.6, colour = "grey25") +
  # top-of-hierarchy anchors: ringed 6W dot + leader label in the top/bottom strip
  ggplot2::geom_point(data = anc, ggplot2::aes(fill = tier), shape = 21, size = 2.8,
                      colour = "black", stroke = 0.5, show.legend = FALSE) +
  ggrepel::geom_text_repel(
    data = anc_up, ggplot2::aes(label = lab, colour = tier),
    ylim = c(0.30, TOP + 0.04), nudge_y = TOP - anc_up$eff6, size = 2.4, fontface = "bold",
    seed = 1, min.segment.length = 0, segment.size = 0.3, segment.colour = "grey60",
    box.padding = 0.5, point.padding = 0.3, force = 6, max.overlaps = Inf, show.legend = FALSE) +
  ggrepel::geom_text_repel(
    data = anc_dn, ggplot2::aes(label = lab, colour = tier),
    ylim = c(BOT - 0.04, -0.30), nudge_y = BOT - anc_dn$eff6, size = 2.4, fontface = "bold",
    seed = 1, min.segment.length = 0, segment.size = 0.3, segment.colour = "grey60",
    box.padding = 0.5, point.padding = 0.3, force = 6, max.overlaps = Inf, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = tier_col, name = "MitoPathway tier") +
  ggplot2::scale_fill_manual(values = tier_col, guide = "none") +
  ggplot2::coord_cartesian(ylim = c(BOT - 0.04, TOP + 0.04)) +
  ggplot2::labs(
    x = "MitoPathways, ranked by Myc mitoPPS effect  (promoted -> demoted)",
    y = "Myc mitoPPS effect\n(Myc+ - WT)",
    title = "Myc reprioritises the mitochondrial compartment",
    subtitle = "Each dot = a MitoPathway (colour = tier); * = padj<0.05; 7 top-level categories labelled; small point/line = 12W (the fade).",
    caption = paste(
      "mitoPPS (Monzel 2025), content-blind (each pathway ratio-normalized to 1): a pure relative-priority shift, Myc+ vs WT.",
      "Promoted (left) = import + biosynthesis; demoted (right) = dynamics/apoptosis/signalling; OXPHOS spans the middle (split, fig02). Exploratory (n=6).",
      sep = "\n")) +
  theme_myc(base_size = 9) +
  ggplot2::theme(
    axis.text.x     = ggplot2::element_blank(),
    axis.ticks.x    = ggplot2::element_blank(),
    legend.position = "right",
    legend.key.size = ggplot2::unit(3.5, "mm"),
    legend.title    = ggplot2::element_text(size = 7)) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.6, alpha = 1)))

# Guard: sourced only to obtain `p` (e.g. Quarto) when myc.fig.nosave = TRUE.
if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "fig02b_reallocation_ranked.pdf"),
             width = fig_w[["double"]], height = 100)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  cat("pathways:", nrow(d), " significant:", sum(d$sig), "\n")
  print(anc[order(-anc$eff6), c("lab", "eff6", "eff12", "padj6")], row.names = FALSE, digits = 2)
  print(p)
}
