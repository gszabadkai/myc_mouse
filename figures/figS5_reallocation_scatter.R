# =============================================================================
# figS5_reallocation_scatter.R -- content vs priority, one point per pathway
# -----------------------------------------------------------------------------
# The synthesis of the two rulers (a collapse of the figS3 heatmap): each of the
# ~144 MitoPathways plotted as CONTENT (Myc set-average raw log2FC) vs PRIORITY
# (Myc mitoPPS effect), at 6W. It makes three things visible at once:
#   1. The two rulers correlate (r~0.76) and the cloud sits almost entirely in the
#      content-UP half -- Myc raises nearly the whole compartment; priority mostly
#      tracks content. So reprioritisation is the SECOND-ORDER spread about the
#      trend, not the main axis.
#   2. The reallocation is DIRECTIONAL: the content-down / priority-up quadrant is
#      EMPTY -- nothing that falls in content gains priority. The demoted pathways
#      are not switched off; they simply rise slower than the biogenesis machinery
#      and lose relative share.
#   3. The residuals off the trend name the candidates: ABOVE (extra priority) =
#      chaperones, immune response, anaplerotic metabolism; BELOW (raised yet
#      demoted) = OXPHOS assembly, cAMP-PKA signalling, dynamics/apoptosis.
#
# CAVEAT (in the caption): mitoPPS is a ratio (~content relative to the compartment
# mean), so part of the correlation is built in -- the SIGNAL is the residual, not
# the trend. Genotype is the clean axis; descriptive/exploratory (n=6/group).
#
# Reads (read-only; NO re-run):
#   results/mitopps_scores.rds      -- $mitopps_pairwise (priority, Myc_effect_6W),
#       $gene_to_pathway, $pathway_tier1_map.
#   results/interaction_results.rds -- raw DESeqResults (content = set log2FC).
#   functions/reconcile_gene_symbols.R -- pathway genes -> Ensembl.
# =============================================================================

source(here::here("figures", "theme_myc.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("figS5 needs ggrepel")
if (!requireNamespace("DESeq2", quietly = TRUE))  stop("figS5 needs DESeq2 (coerce DESeqResults)")

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

# --- trend + residual outliers to label --------------------------------------
fit  <- stats::lm(priority ~ content, d)
rP   <- stats::cor(d$content, d$priority)
d$resid <- stats::resid(fit)
d$sig   <- !is.na(d$padj) & d$padj < 0.05
lab  <- d[order(-abs(d$resid)), ][seq_len(12), ]

p <- ggplot2::ggplot(d, ggplot2::aes(content, priority)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  ggplot2::geom_abline(slope = stats::coef(fit)[2], intercept = stats::coef(fit)[1],
                       colour = "grey45", linewidth = 0.4, linetype = 2) +
  # non-significant points (filled), then significant (ringed) on top
  ggplot2::geom_point(data = d[!d$sig, ], ggplot2::aes(colour = Tier), size = 1.5, alpha = 0.8) +
  ggplot2::geom_point(data = d[d$sig, ], ggplot2::aes(fill = Tier), shape = 21,
                      size = 1.9, colour = "black", stroke = 0.4) +
  ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = pathway, colour = Tier),
                           size = 2.3, fontface = "plain", seed = 1, max.overlaps = Inf,
                           min.segment.length = 0, segment.colour = "grey70",
                           box.padding = 0.35, show.legend = FALSE) +
  # annotate the empty quadrant (content down, priority up)
  ggplot2::annotate("text", x = -0.28, y = 0.42, hjust = 0, size = 2.4, colour = "grey45",
                    fontface = "italic",
                    label = "empty quadrant:\nnothing falls in content\nyet gains priority") +
  ggplot2::annotate("text", x = max(d$content) * 0.62, y = min(d$priority) * 0.98,
                    hjust = 0, size = 2.5, colour = "grey35",
                    label = sprintf("r = %.2f", rP)) +
  ggplot2::scale_colour_manual(values = tier_col, name = "MitoPathway tier", drop = FALSE) +
  ggplot2::scale_fill_manual(values = tier_col, guide = "none") +
  ggplot2::labs(
    x = "Content:  Myc set-average log2FC  (Myc+ - WT, 6W)",
    y = "Priority:  Myc mitoPPS effect  (Myc+ - WT, 6W)",
    title = "Content vs priority per pathway: reprioritisation is a skew on a broad rise",
    subtitle = "The two rulers correlate; the reallocation is the spread below the trend (raised yet demoted). Ringed = padj<0.05.",
    caption = paste(
      "Content = DESeq set-average raw log2FC; priority = relative mitoPPS effect (Monzel 2025). Dashed = linear trend.",
      "mitoPPS is a ratio (~content relative to the compartment mean), so part of the correlation is built in -- the SIGNAL is the residual, not the trend.",
      "Genotype is the clean axis. Descriptive/exploratory (n=6/group).",
      sep = "\n")) +
  theme_myc(base_size = 9) +
  ggplot2::theme(legend.position = "right", legend.key.size = ggplot2::unit(3.5, "mm"),
                 legend.title = ggplot2::element_text(size = 7))

# Guard: sourced only to obtain `p` (e.g. Quarto) when myc.fig.nosave = TRUE.
if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "figS5_reallocation_scatter.pdf"),
             width = fig_w[["double"]], height = 130)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  cat(sprintf("r = %.2f ; quadrants:\n", rP))
  print(table(content_up = d$content > 0, priority_up = d$priority > 0))
  print(lab[order(-lab$resid), c("pathway", "tier", "content", "priority", "resid")], digits = 2)
  print(p)
}
