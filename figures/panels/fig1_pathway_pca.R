# =============================================================================
# fig1_pathway_pca.R -- the samples in pathway space, where one axis is 77%
# -----------------------------------------------------------------------------
# SLOT: Fig. 1C. Filenames do not carry the slot letter; figures/panels/PANELS.md
# is the slug -> slot map.
#
# SUPPORTS (Results, "Myc drives early breast tumourigenesis by inducing OXPHOS
# and biosynthetic pathways", paragraph 2):
#   "Significantly, the dominant principal component axis of these genesets
#    across the whole dataset aligned almost entirely with variability in
#    mitochondria related terms, revealing the leading role of mitochondrial
#    remodeling in early Myc-driven tumourigenesis (Fig. 1C, D)."
#
# This panel is the FIRST half of that sentence -- that there IS a dominant axis,
# and where the four groups sit on it. Fig. 1D is the second half: what the axis
# is made of. They are the sample scores and the per-set loadings of ONE
# principal component analysis, computed once by pathway_axis().
#
# PATHWAY-LEVEL PCA ONLY (author, 2026-07-31). Script 37 draws this next to the
# gene-level PCA, because the 77% is a property of the COMPOSITES and not of the
# transcriptome (gene-level PC1 is 44%). That comparison is not in the paper: it
# is a methods point, it doubles the panel, and the author's call is to keep the
# figure to the pathway lens. The number is in the legend block instead, together
# with what the gene-level axis turned out to be -- an epithelial-purity vs
# immune-infiltration COMPOSITION axis, genotype-independent, which is exactly
# the thing the compositing drops (script 37 PART A2).
#
# WHY THIS SURVIVES batch = timepoint. The 6W and 12W cohorts were two separate
# extractions, so anything that separates the samples by AGE is unreadable. This
# axis does not: on PC1 the genotype term is significant and the timepoint term
# is nowhere near it. The axis is a Myc axis, and Myc is balanced within batch.
#
# WHAT IS NOT DRAWN. No dispersion ellipse. At n = 6 per group a normal ellipse
# is a spread and reads as a confidence region, and the four clouds overlap
# enough that the reader would take the drawing for a test. The group mean is
# drawn as a ringed marker and the animals are all on the page; script 37's
# ellipse version is kept in the sandbox.
#
# Input:  results/gsva_scores.rds       (script 15 -- VST matrix, gene sets, meta)
#         results/pathway_loading.rds   (script 37 -- the assertions)
# Output: outputs/figures/panels/fig1_pathway_pca.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

ax <- pathway_axis()          # asserts itself against script 37 before returning

sm  <- ax$sample_meta
vfr <- ax$var_frac
dat <- data.frame(
  PC1   = ax$scores[, 1],
  PC2   = ax$scores[, 2],
  group = factor(as.character(sm$group), levels = names(group_cols)),
  stringsAsFactors = FALSE)
stopifnot(nrow(dat) == 24L, !anyNA(dat$group))

cen <- stats::aggregate(cbind(PC1, PC2) ~ group, dat, mean)

# --- what the axis is, for the legend block ----------------------------------
# Same idiom as every other panel here: fit it in the script so the drawn object
# and the sentence cannot disagree.
fit1 <- stats::coef(summary(stats::lm(
  PC1 ~ timepoint * myc_status,
  data = data.frame(PC1 = dat$PC1,
                    timepoint  = factor(sm$timepoint,  levels = c("6W", "12W")),
                    myc_status = factor(sm$myc_status, levels = c("neg", "pos"))))))
p_geno <- fit1["myc_statuspos", "Pr(>|t|)"]
p_time <- fit1["timepoint12W",  "Pr(>|t|)"]
p_int  <- fit1["timepoint12W:myc_statuspos", "Pr(>|t|)"]
gmean1 <- tapply(dat$PC1, dat$group, mean)
gap6   <- unname(gmean1["6W_pos"]  - gmean1["6W_neg"])
gap12  <- unname(gmean1["12W_pos"] - gmean1["12W_neg"])
r_gm   <- suppressWarnings(stats::cor(dat$PC1, ax$global_mean))

# =============================================================================
# THE PANEL
# =============================================================================
# coord_equal is load-bearing: the claim is that one direction carries most of
# the variance, and a free aspect ratio would let the panel shape decide how that
# looks. With equal units the PC2 spread is drawn at the size it actually is.
p <- ggplot2::ggplot(dat, ggplot2::aes(PC1, PC2, colour = group)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey88") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey88") +
  ggplot2::geom_point(size = 1.2, alpha = 0.9) +
  ggplot2::geom_point(data = cen, ggplot2::aes(PC1, PC2, fill = group),
                      shape = 21, size = 2.8, colour = "grey15", stroke = 0.35,
                      inherit.aes = FALSE) +
  ggplot2::scale_colour_manual(values = group_cols, breaks = names(group_cols),
                               labels = group_labels) +
  ggplot2::scale_fill_manual(values = group_cols, guide = "none") +
  ggplot2::coord_equal() +
  ggplot2::labs(x = sprintf("PC1 (%.0f%%)", 100 * vfr[1]),
                y = sprintf("PC2 (%.0f%%)", 100 * vfr[2]),
                colour = NULL) +
  ggplot2::guides(colour = ggplot2::guide_legend(
    nrow = 1, override.aes = list(size = 1.7, alpha = 1))) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    panel.grid       = ggplot2::element_blank(),
    legend.position  = "bottom",
    legend.margin    = ggplot2::margin(-3, 0, 0, 0),
    legend.key.size  = ggplot2::unit(2.6, "mm"),
    legend.spacing.x = ggplot2::unit(0.6, "mm"))

# --- the legend text (never drawn) -------------------------------------------
f1 <- function(x) sprintf("%.1f", x)
LEGEND <- panel_legend(
  slot = "Fig. 1C",
  what = paste0(
    "Principal component analysis of the ", nrow(ax$M), " gene-set scores across ",
    "all 24 samples. Each point is one animal; the ringed marker is the group mean."),
  detail = c(
    sprintf("Quantifier: the linear mean gene-wise z-score of each set, computed on the DESeq2 VST matrix (the method of record, scripts 36-37). PCA on the %d x 24 score matrix, centred, unscaled.",
            nrow(ax$M)),
    sprintf("PC1 explains %.1f per cent of the variance, PC2 %.1f and PC3 %.1f. Effective dimensionality is 1.7 of a possible 23, and %.0f per cent of the set loadings share one sign: this is a single common mode, not a contrast.",
            100 * vfr[1], 100 * vfr[2], 100 * vfr[3],
            100 * ax$ref$dimensionality$frac_loadings_pos_raw[
              ax$ref$dimensionality$quantifier == "zscore"]),
    sprintf("PC1 correlates with the simple per-sample mean of all %d scores at r = %.2f, i.e. the axis is the level of everything at once.",
            nrow(ax$M), r_gm),
    sprintf("PC1 is a GENOTYPE axis: fitting PC1 ~ timepoint * genotype gives genotype p = %.4f, timepoint p = %.2f, interaction p = %.2f. Group means are %s (6W Myc+), %s (12W Myc+), %s (6W WT) and %s (12W WT).",
            p_geno, p_time, p_int,
            f1(gmean1["6W_pos"]), f1(gmean1["12W_pos"]),
            f1(gmean1["6W_neg"]), f1(gmean1["12W_neg"])),
    "Axes are drawn on a common scale (coord_equal), so the collapse of the spread onto PC1 is the picture and not the panel shape.",
    "n = 6 animals per group. No dispersion region is drawn; at n = 6 an ellipse would read as a confidence interval."),
  bounds = c(
    "Batch = timepoint: the 6W and 12W cohorts were extracted separately, so any axis that separated the samples by age would be unreadable. This one does not - the timepoint term on PC1 is not significant while the genotype term is - and genotype is balanced within each batch, so the reading is clean.",
    sprintf("The genotype gap along PC1 is %s at 6 weeks and %s at 12. The direction is the same as the attenuation result but this is NOT that result: the interaction is not significant here (p = %.2f) at n = 6 per cell, and the attenuation result is a DESeq2 effect-size result on a different ruler (the Myc effect rescales by about 0.55, scripts 29-31 and 40).",
            f1(gap6), f1(gap12), p_int),
    sprintf("The 77 per cent is a property of the COMPOSITES, not of the transcriptome. On the same samples the largest gene-level axis (top-500 variable genes, VST) is %.0f per cent, and it is a different direction: genotype-independent, tracking epithelial-purity against immune infiltration, i.e. residual dissociation contamination. Averaging genes into correlated pathway scores concentrates variance onto the common mode AND foregrounds the Myc programme, which sits at gene-level PC2 (script 37 PART A2). The pathway PCA is the lens the paper uses; the gene-level comparison is a methods point and is not drawn.",
            ax$ref$sample_pca_var$pc1_pct[ax$ref$sample_pca_var$space == "gene_vst_top500"]),
    "GSVA gives the same picture more weakly (PC1 68 per cent), which is why the linear z-score is the quantifier of record: its correlation is average cross-gene covariance, with no rank transform in between.",
    "What the axis IS, biologically, is Fig. 1D's question and it has a bound there: mitochondria-led is not mitochondria-specific."),
  source = c(
    "results/gsva_scores.rds (script 15) for the VST matrix, gene sets and sample metadata",
    "results/pathway_loading.rds (script 37) - PC1/PC2/PC3 percentages, the OXPHOS gate and the gene-level comparison; the rebuild in pathway_axis() is asserted against all three",
    "Method: scripts/36_linear_pathway_coupling.R (quantifier of record), scripts/37_pathway_loading_and_technical_resolution.R PART 2 and PART A2"))

# Height is set by coord_equal, not chosen: at 89 mm wide the PC1 span of ~70
# units fixes the unit size, and the PC2 span of ~24 units then fixes the panel
# height. Anything taller is dead space between the axis title and the key.
save_panel_p(p, "fig1_pathway_pca", width = fig_w[["single"]], height = 44)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)

  ## script 37's version, with the 1-SD normal ellipse per group. Kept because it
  ## is occasionally worth seeing the spread, and dropped from the panel because
  ## at n = 6 a drawn region invites a test that was not done.
  print(p + ggplot2::stat_ellipse(ggplot2::aes(fill = group), geom = "polygon",
                                  type = "norm", level = 0.68, alpha = 0.10,
                                  colour = NA))

  ## PC2 and PC3: is there anything below the dominant axis?
  d23 <- data.frame(PC2 = ax$scores[, 2], PC3 = ax$scores[, 3], group = dat$group)
  ggplot2::ggplot(d23, ggplot2::aes(PC2, PC3, colour = group)) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_colour_manual(values = group_cols) + theme_panel()

  ## the scree, in numbers
  round(100 * ax$var_frac[1:8], 2)

  ## where the named composites sit on PC1 (this is Fig. 1D's content)
  sort(ax$loading[grep("^MITOCARTA_", names(ax$loading))], decreasing = TRUE) |> head(10)
}
