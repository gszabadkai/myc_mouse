# =============================================================================
# figS5_reallocation_scatter.R -- priority vs content, one point per pathway
# -----------------------------------------------------------------------------
# The synthesis of the two rulers (a collapse of the figS3 heatmap): each of the
# ~144 MitoPathways plotted as PRIORITY (Myc mitoPPS effect, x) vs CONTENT (Myc
# set-average raw log2FC, y), at 6W, with MARGINAL DENSITIES on both axes. The
# message the margins carry:
#   * CONTENT (y) is almost all POSITIVE (~94% of pathways up) -- a narrow,
#     one-sided distribution: Myc raises nearly the whole compartment.
#   * PRIORITY (x) spans a WIDE range in BOTH directions (~56% up), centred on 0.
# So the reprioritisation is that WIDE PRIORITY SPREAD sitting on top of a
# near-uniform content increase; the demoted pathways are not switched off, they
# rise slower than the biogenesis machinery and lose relative share. No trend line
# is drawn (mitoPPS is a ratio ~content-relative-to-the-compartment-mean, so a
# content/priority slope is partly built in and would mislead). Only the 7 TOP-OF-
# HIERARCHY (Level-1) categories are labelled, as anchors.
#
# Points ringed = padj<0.05. Genotype is the clean axis; descriptive/exploratory.
#
# Reads (read-only; NO re-run):
#   results/mitopps_scores.rds      -- $mitopps_pairwise (priority, Myc_effect_6W),
#       $gene_to_pathway, $pathway_levels (Level-1), $pathway_tier1_map.
#   results/interaction_results.rds -- raw DESeqResults (content = set log2FC).
#   functions/reconcile_gene_symbols.R -- pathway genes -> Ensembl.
# =============================================================================

source(here::here("figures", "theme_myc.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))
if (!requireNamespace("ggrepel", quietly = TRUE))   stop("figS5 needs ggrepel")
if (!requireNamespace("patchwork", quietly = TRUE)) stop("figS5 needs patchwork (marginals)")
if (!requireNamespace("DESeq2", quietly = TRUE))    stop("figS5 needs DESeq2 (coerce DESeqResults)")

out_dir <- here::here("outputs", "figures")

mp   <- readRDS(here::here("results", "mitopps_scores.rds"))
ir   <- readRDS(here::here("results", "interaction_results.rds"))
tier <- mp$pathway_tier1_map
gp   <- mp$gene_to_pathway; pcol <- names(gp)[1]; gcol <- names(gp)[2]

L6  <- { r <- as.data.frame(ir$myc_6W_raw); stats::setNames(r$log2FoldChange, rownames(r)) }
set_lfc <- function(pw) {
  e <- recon_to_ensembl(unique(gp[[gcol]][gp[[pcol]] == pw]), names(L6))
  v <- L6[e]; if (length(v)) mean(v, na.rm = TRUE) else NA_real_
}
rel <- mp$mitopps_pairwise[mp$mitopps_pairwise$contrast == "Myc_effect_6W", ]

pw <- intersect(unique(rel$pathway), names(tier))
d <- data.frame(
  pathway = pw, tier = unname(tier[pw]),
  content  = vapply(pw, set_lfc, numeric(1)),
  priority = rel$diff[match(pw, rel$pathway)],
  padj     = rel$padj[match(pw, rel$pathway)],
  stringsAsFactors = FALSE)
d <- d[!is.na(d$content) & !is.na(d$priority), ]

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
d$Tier <- factor(unname(tier_lab[d$tier]), levels = unname(tier_lab[tier_lv]))
d$sig  <- !is.na(d$padj) & d$padj < 0.05

pct_c <- round(100 * mean(d$content  > 0))
pct_p <- round(100 * mean(d$priority > 0))

# --- the 7 top-of-hierarchy (Level-1) anchors --------------------------------
l1  <- mp$pathway_levels$Pathway[mp$pathway_levels$Level == "Pathway_Level1"]
anc <- d[d$pathway %in% l1, ]
anc$lab <- unname(tier_lab[names(tier)[match(anc$pathway, names(tier))]])

XLIM <- c(min(d$priority) - 0.03, max(d$priority) + 0.05)     # x = priority
YLIM <- c(min(d$content)  - 0.05, max(d$content)  + 0.08)     # y = content

# --- main scatter (priority x, content y; no trend line) ---------------------
main <- ggplot2::ggplot(d, ggplot2::aes(priority, content)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  ggplot2::geom_point(data = d[!d$sig, ], ggplot2::aes(colour = Tier), size = 1.5, alpha = 0.8) +
  ggplot2::geom_point(data = d[d$sig, ], ggplot2::aes(fill = Tier), shape = 21,
                      size = 1.9, colour = "black", stroke = 0.4) +
  ggplot2::geom_point(data = anc, ggplot2::aes(fill = Tier), shape = 21, size = 2.8,
                      colour = "black", stroke = 0.5, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = anc, ggplot2::aes(label = lab, colour = Tier),
                           size = 2.5, fontface = "bold", seed = 1, max.overlaps = Inf,
                           min.segment.length = 0, segment.colour = "grey55",
                           box.padding = 0.6, point.padding = 0.4, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = tier_col, name = NULL, drop = FALSE) +
  ggplot2::scale_fill_manual(values = tier_col, guide = "none") +
  ggplot2::coord_cartesian(xlim = XLIM, ylim = YLIM, expand = FALSE) +
  ggplot2::labs(x = "Priority:  Myc mitoPPS effect  (Myc+ - WT, 6W)",
                y = "Content:  Myc set-average log2FC  (Myc+ - WT, 6W)") +
  theme_myc(base_size = 9) +
  ggplot2::theme(legend.position = "bottom", legend.key.size = ggplot2::unit(3.2, "mm"),
                 legend.text = ggplot2::element_text(size = 7),
                 plot.margin = ggplot2::margin(2, 2, 2, 2)) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 2.2)))

# --- marginal densities (x priority = wide; y content = one-sided positive) ---
mtheme <- ggplot2::theme_void() + ggplot2::theme(plot.margin = ggplot2::margin(2, 2, 2, 2))
topd <- ggplot2::ggplot(d, ggplot2::aes(priority)) +
  ggplot2::geom_density(fill = "grey80", colour = "grey45", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  ggplot2::annotate("text", x = XLIM[1], y = Inf, hjust = 0, vjust = 1.3, size = 2.3,
                    colour = "grey35", label = sprintf("priority: %d%% up (wide)", pct_p)) +
  ggplot2::coord_cartesian(xlim = XLIM, expand = FALSE) + mtheme
rightd <- ggplot2::ggplot(d, ggplot2::aes(content)) +
  ggplot2::geom_density(fill = "grey80", colour = "grey45", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  ggplot2::annotate("text", x = YLIM[2], y = Inf, hjust = 1, vjust = 1.3, size = 2.3,
                    colour = "grey35", label = sprintf("content: %d%% up", pct_c)) +
  ggplot2::coord_flip(xlim = YLIM, expand = FALSE) + mtheme

p <- patchwork::wrap_plots(topd, patchwork::plot_spacer(), main, rightd,
                           ncol = 2, nrow = 2, widths = c(5, 1), heights = c(1, 4.6)) +
  patchwork::plot_annotation(
    title = "Content is almost all up; priority spans a wide range -- that spread is the reallocation",
    subtitle = sprintf("Per pathway (6W): content (log2FC) is one-sided positive (%d%% up); priority (mitoPPS) is wide, both ways (%d%% up). Ringed = padj<0.05.", pct_c, pct_p),
    caption = paste(
      "Content = DESeq set-average raw log2FC; priority = relative mitoPPS effect (Monzel 2025), content-blind (a ratio ~content vs the compartment mean).",
      "No trend line is drawn (a content/priority slope is partly built in). The reallocation is the wide priority spread over a near-uniform content rise. Genotype clean; exploratory (n=6).",
      sep = "\n"),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 10.5),
      plot.subtitle = ggplot2::element_text(size = 8.2, colour = "grey20"),
      plot.caption  = ggplot2::element_text(size = 6.2, hjust = 0, colour = "grey30",
                                            lineheight = 1.1)))

if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "figS5_reallocation_scatter.pdf"),
             width = fig_w[["double"]], height = 135)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  cat(sprintf("content %% up = %d ; priority %% up = %d\n", pct_c, pct_p))
  print(table(content_up = d$content > 0, priority_up = d$priority > 0))
  print(anc[order(anc$priority), c("lab", "priority", "content")], row.names = FALSE, digits = 2)
  print(p)
}
