# =============================================================================
# figS6_reallocation_temporal_ranked.R -- mitoPPS reallocation over the timeline
# -----------------------------------------------------------------------------
# The temporal counterpart of fig02 (which reads the Myc effect WITHIN a
# timepoint). Here the contrast is the within-genotype change over the timeline:
# WT 6->12W (Temporal_Myc-) and Myc+ 6->12W (Temporal_Myc+). Pathways are ranked
# by the WT change (the developmental reallocation of the WT background), with a
# DROPLINE to the Myc+ change -- the dropline is where Myc BENDS the WT trajectory.
#
# The point (distinct from the earlier Myc genotype effect): the WT background
# itself reprioritises over time; WT-time and Myc+-time mitoPPS correlate ~0.75, so
# most of the temporal reallocation is shared (developmental), and Myc bends it --
# the biogenesis/biosynthesis arms it promoted at 6W fade hardest by 12W.
#
# BATCH-CONFOUND: the time axis is batch-confounded (batch = timepoint), so the
# individual WT/Myc+ changes are DESCRIPTIVE, not claims; the clean quantity is
# their DIVERGENCE (the dropline; the y=x deviation in figS7). No significance
# asterisks are drawn (individual temporal padj would mislead). Exploratory (n=6).
#
# Reads (read-only; NO re-run -- mitoPPS-only):
#   results/mitopps_scores.rds -- $mitopps_pairwise (Temporal_Myc- / Temporal_Myc+
#       diff), $pathway_levels (Level-1 top hierarchy), $pathway_tier1_map.
# =============================================================================

source(here::here("figures", "theme_myc.R"))
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  stop("figS6 needs ggrepel for the top-hierarchy leader labels")
}

out_dir <- here::here("outputs", "figures")

mp   <- readRDS(here::here("results", "mitopps_scores.rds"))
pw   <- mp$mitopps_pairwise
tier <- mp$pathway_tier1_map

w <- pw[pw$contrast == "Temporal_Myc-", c("pathway", "diff")]   # WT  6->12W
m <- pw[pw$contrast == "Temporal_Myc+", c("pathway", "diff")]   # Myc+ 6->12W
names(w)[2] <- "wt"; names(m)[2] <- "myc"
d <- merge(w, m, by = "pathway")
d$tier <- unname(tier[d$pathway])
d <- d[!is.na(d$tier) & !is.na(d$wt), ]

# --- tier labels + palette (same as fig02; Metabolism greyed) ----------------
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

# --- rank DESCENDING by the WT change (promoted left -> demoted right) --------
d <- d[order(-d$wt), ]
d$rank <- seq_len(nrow(d))

# --- 7 top-of-hierarchy anchors (labelled), split by WT sign into strips ------
l1  <- mp$pathway_levels$Pathway[mp$pathway_levels$Level == "Pathway_Level1"]
anc <- d[d$pathway %in% l1, ]
anc$lab <- unname(tier_lab[names(tier)[match(anc$pathway, names(tier))]])
TOP <- 0.30; BOT <- -0.28
anc_up <- anc[anc$wt > 0, ]; anc_dn <- anc[anc$wt <= 0, ]

p <- ggplot2::ggplot(d, ggplot2::aes(rank, wt)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  # WT -> Myc+ dropline (the Myc-specific bend)
  ggplot2::geom_segment(ggplot2::aes(xend = rank, y = wt, yend = myc),
                        colour = "grey78", linewidth = 0.25) +
  ggplot2::geom_point(ggplot2::aes(y = myc, colour = tier), shape = 1, size = 0.9, stroke = 0.4) +
  ggplot2::geom_point(ggplot2::aes(colour = tier), size = 1.3, alpha = 0.9) +
  # anchors: ringed WT dot + leader label in the top/bottom strip
  ggplot2::geom_point(data = anc, ggplot2::aes(fill = tier), shape = 21, size = 2.8,
                      colour = "black", stroke = 0.5, show.legend = FALSE) +
  ggrepel::geom_text_repel(
    data = anc_up, ggplot2::aes(label = lab, colour = tier),
    nudge_y = TOP - anc_up$wt, direction = "x", size = 2.4, fontface = "bold",
    seed = 1, min.segment.length = 0, segment.size = 0.3, segment.colour = "grey60",
    box.padding = 0.45, force = 9, force_pull = 0.15, max.overlaps = Inf, show.legend = FALSE) +
  ggrepel::geom_text_repel(
    data = anc_dn, ggplot2::aes(label = lab, colour = tier),
    nudge_y = BOT - anc_dn$wt, direction = "x", size = 2.4, fontface = "bold",
    seed = 1, min.segment.length = 0, segment.size = 0.3, segment.colour = "grey60",
    box.padding = 0.45, force = 9, force_pull = 0.15, max.overlaps = Inf, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = tier_col, name = "MitoPathway tier") +
  ggplot2::scale_fill_manual(values = tier_col, guide = "none") +
  ggplot2::coord_cartesian(ylim = c(BOT - 0.04, TOP + 0.04)) +
  ggplot2::labs(
    x = "MitoPathways, ranked by WT 6->12W mitoPPS change  (promoted -> demoted)",
    y = "mitoPPS change\n(6W -> 12W)",
    title = "The WT background reprioritises over the timeline; Myc bends it",
    subtitle = "Ranked by the WT change (filled); dropline to Myc+ (open) = the Myc-specific bend; labels = 7 top-level categories.",
    caption = paste(
      "mitoPPS (Monzel 2025), content-blind. TIME axis batch-confounded: the individual WT/Myc+ changes are descriptive, not claims.",
      "The clean quantity is their divergence (the dropline; figS7). Exploratory (n=6).",
      sep = "\n")) +
  theme_myc(base_size = 9) +
  ggplot2::theme(
    axis.text.x     = ggplot2::element_blank(),
    axis.ticks.x    = ggplot2::element_blank(),
    legend.position = "right",
    legend.key.size = ggplot2::unit(3.5, "mm"),
    legend.title    = ggplot2::element_text(size = 7)) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.6, alpha = 1)))

if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "figS6_reallocation_temporal_ranked.pdf"),
             width = fig_w[["double"]], height = 100)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  cat("pathways:", nrow(d), "\n")
  print(anc[order(-anc$wt), c("lab", "wt", "myc")], row.names = FALSE, digits = 2)
  print(p)
}
