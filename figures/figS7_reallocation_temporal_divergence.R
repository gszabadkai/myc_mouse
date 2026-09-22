# =============================================================================
# figS7_reallocation_temporal_divergence.R -- WT vs Myc+ temporal mitoPPS
# -----------------------------------------------------------------------------
# The synthesis of the timeline reallocation: each MitoPathway's mitoPPS change
# over 6W->12W in WT (x) vs in Myc+ (y). The y=x diagonal is "Myc+ changes exactly
# as WT does" -- purely developmental / shared; OFF the diagonal is where Myc BENDS
# the WT trajectory (the interaction). Most pathways sit near the diagonal
# (cor ~0.75, 70% same direction): the WT background reprioritises over the
# timeline on its own, and Myc modifies it -- a DIFFERENT phenomenon from the Myc
# genotype effect (fig02). Below the diagonal = Myc+ fades relative to WT (the
# biogenesis/biosynthesis arms it promoted at 6W).
#
# BATCH-CONFOUND: the two axes are each batch-confounded (batch = timepoint) and so
# DESCRIPTIVE; the clean quantity is the DEVIATION from y=x (the interaction, where
# the batch effect cancels). Exploratory (n=6).
#
# Reads (read-only; NO re-run -- mitoPPS-only):
#   results/mitopps_scores.rds -- $mitopps_pairwise (Temporal_Myc- / Temporal_Myc+
#       diff), $pathway_levels (Level-1 top hierarchy), $pathway_tier1_map.
# =============================================================================

source(here::here("figures", "theme_myc.R"))
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("figS7 needs ggrepel")

out_dir <- here::here("outputs", "figures")

mp   <- readRDS(here::here("results", "mitopps_scores.rds"))
pw   <- mp$mitopps_pairwise
tier <- mp$pathway_tier1_map

w <- pw[pw$contrast == "Temporal_Myc-", c("pathway", "diff")]
m <- pw[pw$contrast == "Temporal_Myc+", c("pathway", "diff")]
names(w)[2] <- "wt"; names(m)[2] <- "myc"
d <- merge(w, m, by = "pathway")
d$tier <- unname(tier[d$pathway])
d <- d[!is.na(d$tier) & !is.na(d$wt), ]

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
d$Tier <- factor(unname(tier_lab[d$tier]), levels = unname(tier_lab[tier_lv]))

rP  <- stats::cor(d$wt, d$myc)
l1  <- mp$pathway_levels$Pathway[mp$pathway_levels$Level == "Pathway_Level1"]
anc <- d[d$pathway %in% l1, ]
anc$lab <- unname(tier_lab[names(tier)[match(anc$pathway, names(tier))]])
LIM <- c(min(c(d$wt, d$myc)) - 0.03, max(c(d$wt, d$myc)) + 0.03)   # square, symmetric

p <- ggplot2::ggplot(d, ggplot2::aes(wt, myc)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey45", linewidth = 0.4,
                       linetype = 2) +
  ggplot2::annotate("text", x = LIM[2], y = LIM[2], hjust = 1, vjust = -0.4, angle = 45,
                    size = 2.3, colour = "grey45", label = "Myc+ = WT (developmental)") +
  ggplot2::annotate("text", x = LIM[2], y = LIM[1], hjust = 1, vjust = 0, size = 2.3,
                    colour = "grey40", fontface = "italic",
                    label = "below the line:\nMyc+ fades vs WT") +
  ggplot2::geom_point(ggplot2::aes(colour = Tier), size = 1.6, alpha = 0.85) +
  ggplot2::geom_point(data = anc, ggplot2::aes(fill = Tier), shape = 21, size = 2.8,
                      colour = "black", stroke = 0.5, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = anc, ggplot2::aes(label = lab, colour = Tier),
                           size = 2.5, fontface = "bold", seed = 1, max.overlaps = Inf,
                           min.segment.length = 0, segment.colour = "grey55",
                           box.padding = 0.6, point.padding = 0.4, show.legend = FALSE) +
  ggplot2::annotate("text", x = LIM[1], y = LIM[2], hjust = 0, vjust = 1, size = 2.6,
                    colour = "grey35", label = sprintf("r = %.2f", rP)) +
  ggplot2::scale_colour_manual(values = tier_col, name = NULL, drop = FALSE) +
  ggplot2::scale_fill_manual(values = tier_col, guide = "none") +
  ggplot2::coord_equal(xlim = LIM, ylim = LIM, expand = FALSE) +
  ggplot2::labs(
    x = "WT 6->12W mitoPPS change", y = "Myc+ 6->12W mitoPPS change",
    title = "Timeline reallocation: development vs Myc-specific change",
    subtitle = "WT vs Myc+ 6->12W mitoPPS per pathway. On y=x = developmental; off it = Myc bends it.",
    caption = paste(
      "mitoPPS content-blind. Both axes batch-confounded, descriptive; the clean quantity is the y=x deviation.",
      "Most pathways track the diagonal (r=0.75) = shared development; below-diagonal = the biogenesis arms Myc fades. Exploratory (n=6).",
      sep = "\n")) +
  theme_myc(base_size = 9) +
  ggplot2::theme(legend.position = "bottom", legend.key.size = ggplot2::unit(3.2, "mm"),
                 legend.text = ggplot2::element_text(size = 7)) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, override.aes = list(size = 2.4)))

if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "figS7_reallocation_temporal_divergence.pdf"),
             width = 150, height = 155)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  cat(sprintf("r = %.2f ; %% same direction = %.0f\n", rP, 100 * mean(sign(d$wt) == sign(d$myc))))
  print(anc[order(anc$myc - anc$wt), c("lab", "wt", "myc")], row.names = FALSE, digits = 2)
  print(p)
}
