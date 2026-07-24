# =============================================================================
# fig02b_reallocation_ranked.R -- compact horizontal ranked dot chart
# -----------------------------------------------------------------------------
# The compact alternative to fig02_reallocation.R (which labels all ~144 pathways,
# tall). Same data, same ordering, stripped for a MAIN panel: horizontal, NO
# per-pathway labels, colour = tier, asterisks mark padj<0.05, and only the 7
# TOP-OF-HIERARCHY (Level-1) MitoPathways are labelled with leader lines -- they
# are themselves scored composites, so they anchor where each top-level category
# sits on the ranking. The message: reprioritisation is real and structured --
# the demoted (left) tail is dynamics/apoptosis/signalling, the promoted (right)
# tail is import + biosynthesis, OXPHOS spans the middle (its subunit/assembly
# split is fig02 / figS3).
#
# Myc mitoPPS effect (Myc+ - WT) at 6W; each MitoPathway ratio-normalized to 1
# (content-blind), genotype is the clean axis. Descriptive/exploratory (n=6/group).
#
# Reads (read-only; NO re-run):
#   results/mitopps_scores.rds -- $mitopps_pairwise (diff/padj, Myc_effect_6W),
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

e6 <- pw[pw$contrast == "Myc_effect_6W", c("pathway", "diff", "padj")]
e6$tier <- unname(tier[e6$pathway])
e6 <- e6[!is.na(e6$tier) & !is.na(e6$diff), ]

# --- tier labels + palette (same as fig02; Metabolism greyed) ----------------
tier_lv  <- c("Protein import, sorting and homeostasis", "Mitochondrial central dogma",
              "OXPHOS", "Metabolism", "Signaling",
              "Mitochondrial dynamics and surveillance", "Small molecule transport")
tier_lab <- c("Import / homeostasis", "Central dogma", "OXPHOS", "Metabolism",
              "Signaling", "Dynamics & surveillance", "SM transport")
names(tier_lab) <- tier_lv
tier_col <- c("Import / homeostasis"    = "#D55E00", "Central dogma" = "#E69F00",
              "OXPHOS"                  = "#009E73", "Metabolism"    = "grey72",
              "Signaling"               = "#56B4E9",
              "Dynamics & surveillance" = "#0072B2", "SM transport"  = "#CC79A7")
e6$tier <- factor(unname(tier_lab[e6$tier]), levels = unname(tier_lab[tier_lv]))

# --- rank on the effect (demoted left -> promoted right) ----------------------
e6 <- e6[order(e6$diff), ]
e6$rank <- seq_len(nrow(e6))
e6$sig  <- !is.na(e6$padj) & e6$padj < 0.05

# --- the 7 top-of-hierarchy (Level-1) pathways = the anchor labels ------------
l1 <- mp$pathway_levels$Pathway[mp$pathway_levels$Level == "Pathway_Level1"]
anc <- e6[e6$pathway %in% l1, ]
anc$lab <- unname(tier_lab[names(tier)[match(anc$pathway, names(tier))]])   # short tier name
# (Level-1 pathway == its tier, so the short tier label names the anchor)

p <- ggplot2::ggplot(e6, ggplot2::aes(rank, diff)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  ggplot2::geom_point(ggplot2::aes(colour = tier), size = 1.3, alpha = 0.9) +
  # significant pathways: a small asterisk just above the dot
  ggplot2::geom_text(data = e6[e6$sig, ], ggplot2::aes(y = diff + 0.018),
                     label = "*", size = 2.6, colour = "grey25") +
  # top-of-hierarchy anchors: larger ringed dot + leader-line label
  ggplot2::geom_point(data = anc, ggplot2::aes(fill = tier), shape = 21, size = 2.8,
                      colour = "black", stroke = 0.5, show.legend = FALSE) +
  ggrepel::geom_text_repel(
    data = anc, ggplot2::aes(label = lab, colour = tier), size = 2.7, fontface = "bold",
    seed = 1, min.segment.length = 0, box.padding = 0.7, max.overlaps = Inf,
    segment.size = 0.3, segment.colour = "grey55", show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = tier_col, name = "MitoPathway tier") +
  ggplot2::scale_fill_manual(values = tier_col, guide = "none") +
  ggplot2::labs(
    x = "MitoPathways, ranked by Myc mitoPPS effect  (demoted -> promoted)",
    y = "Myc mitoPPS effect\n(Myc+ - WT)",
    title = "Myc reprioritises the mitochondrial compartment",
    subtitle = "Each dot = one MitoPathway (colour = tier); * = padj<0.05; labels mark the 7 top-level categories.",
    caption = paste(
      "mitoPPS (Monzel 2025), content-blind (each pathway ratio-normalized to 1): a pure relative-priority shift, Myc+ vs WT at 6W.",
      "Demoted tail = dynamics/apoptosis/signalling; promoted = import + biosynthesis; OXPHOS spans the middle (split, fig02). Exploratory (n=6).",
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
             width = fig_w[["double"]], height = 95)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  cat("pathways:", nrow(e6), " significant:", sum(e6$sig), "\n")
  print(anc[order(anc$diff), c("lab", "diff", "padj", "rank")], row.names = FALSE, digits = 2)
  print(p)
}
