# =============================================================================
# figS1_axis_loadings.R -- what the dominant axis is made of
# -----------------------------------------------------------------------------
# SLOT: Fig. S1C (was Fig. 1D for a day; the author's call on 2026-07-31 is that
# the enrichment ranking is the more relevant main-figure panel, so the two
# swapped). Filenames do not carry the slot letter; figures/panels/PANELS.md is
# the slug -> slot map.
#
# SUPPORTS (Results, "Myc drives early breast tumourigenesis by inducing OXPHOS
# and biosynthetic pathways", paragraph 2), the second half of the sentence
# Fig. 1C opens -- and the main figure cites it alongside 1C:
#   "the dominant principal component axis of these genesets across the whole
#    dataset aligned almost entirely with variability in mitochondria related
#    terms, revealing the leading role of mitochondrial remodeling in early
#    Myc-driven tumourigenesis (Fig. 1C, D)."
#
# Every scored gene set, ranked by how strongly it loads on the PC1 of Fig. 1C.
# Same PCA, computed once by pathway_axis(); this panel is its loading vector.
#
# THE PANEL HAS TO CARRY THE CLAIM AND ITS BOUND AT THE SAME TIME, which is why
# it is the whole ranking rather than a chosen handful of programmes, and why the
# fill is the three-class mitochondrial definition of Fig. S1B.
#
#   THE CLAIM. Mitochondrial sets are 44 per cent of the library and 100 per cent
#   of the top 50. Median |loading| is 0.90 for MitoCarta sets against 0.61 for
#   non-mitochondrial ones. On this axis the mitochondrial block really is the
#   strongest-loading block, and OXPHOS as a composite tracks it at 0.98 -- it
#   does not merely participate in the axis, it very nearly IS the axis.
#
#   THE BOUND, three parts, all visible on the page.
#   (1) The very top is a BUILD TAUTOLOGY. The Gray _MITO transcription-factor
#       lanes are that factor's target programme intersected with MitoCarta, so
#       they load by mitochondrial gene CONTENT, not by regulator identity; 94 of
#       the top 100 are these. They are drawn in their own colour rather than
#       hidden, because deleting them would flatter the result and merging them
#       into the mitochondrial block would inflate it.
#   (2) Non-mitochondrial GROWTH programmes reach 0.98 -- the human breast-cancer
#       biclusters, nucleotide metabolism, the pentose-phosphate pathway, the Myc
#       and E2F regulons. The axis is not choosing mitochondria over them.
#   (3) redox is the control that fails. It is a mitochondrial-function gene set
#       and it loads 0.67, in the bottom half -- and it is the one axis Myc does
#       not drive (script 35). So what the axis tracks is not "mitochondria": it
#       is a coordinated Myc anabolic-proliferative-mitochondrial GROWTH state.
#
# LOADING IS COVARIATION, NOT PRIMACY. Nothing here says the mitochondrial change
# causes the rest; the analysis of record says so in as many words (script 37,
# enrichment_verdict). A primacy claim needs perturbation, which is what the cell
# experiments are for.
#
# WHICH AXIS. The raw PC1 of Fig. 1C, so that the two panels are the same object.
# Script 37 ranks sets on the DESIGN-RESIDUAL axis instead (the within-group
# phenotype it dissects); the two rankings correlate 0.968, and pathway_axis()
# asserts that before returning, so this panel is not a different result drawn on
# a friendlier axis.
#
# Input:  results/gsva_scores.rds       (script 15)
#         results/pathway_loading.rds   (script 37 -- classification + assertions)
# Output: outputs/figures/panels/figS1_axis_loadings.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

if (!requireNamespace("ggrepel", quietly = TRUE)) stop("figS1_axis_loadings needs ggrepel")

ax  <- pathway_axis()
cls <- ax$ref$mito_classification

dat <- data.frame(set = names(ax$loading), loading = as.numeric(ax$loading),
                  stringsAsFactors = FALSE) |>
  dplyr::inner_join(dplyr::select(cls, set, class3), by = "set") |>
  dplyr::mutate(abs_load = abs(loading)) |>
  dplyr::arrange(dplyr::desc(abs_load)) |>
  dplyr::mutate(rank = dplyr::row_number(),
                class3 = factor(class3, levels = names(mito_class_cols)))
stopifnot(nrow(dat) == nrow(ax$M), !anyNA(dat$class3))

# --- the two reference lines --------------------------------------------------
# median = where an arbitrary library set sits, so "high loading" has a referent.
# redox = the mitochondrial-function control; its three sets, not the composite.
redox_sets <- c("GS_METAB_REDOX", "GS_METAB_GLUTATHIONE", "GS_METAB_REACTIVE_OXYGEN")
stopifnot(all(redox_sets %in% dat$set))
med_ref   <- stats::median(dat$abs_load)
redox_ref <- stats::median(dat$abs_load[dat$set %in% redox_sets])

# --- the labelled exemplars ---------------------------------------------------
# Named explicitly, not taken off the top of the ranking: the point of the panel
# is that the top of the ranking is not the interesting part. Three mitochondrial
# (one of them the build tautology) and four non-mitochondrial, spread down the
# curve to the redox control. Short labels because the panel is 89 mm; the exact
# set names are in the legend block.
#
# WHERE THE LABELS GO, and why it is forced rather than chosen. Five of the seven
# anchors sit in the left third of the ranking, and a label is ~200 rank units
# wide, so any left-hand column puts one label across the next label's leader --
# which is what the first version did. The one region that is provably free is
# the wedge ABOVE the curve on the right: the curve falls monotonically, so a
# straight line from a high-rank anchor to a point up and to the right can never
# re-cross it. So all seven labels are right-aligned at the right edge in RANK
# ORDER, and because both the anchors and the label slots are monotone in y the
# leaders fan out without crossing each other or any label box. Positions are in
# data units and repel is used only to draw the leaders (force = 0 leaves the
# given position alone).
lab_spec <- tibble::tribble(
  ~set,                                    ~label,                          ~y,
  "TFT_MYC_GRAY_BA_LE_MITO",               "Myc regulon (MitoCarta subset)", 1.05,
  "MITOCARTA_OXPHOS",                      "OXPHOS",                         0.98,
  "METABRIC_MB2_HI_CV_GROUP1",             "human BRCA-Myc (MB2)",           0.91,
  "MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA", "mitochondrial biogenesis",       0.84,
  "MYC_felsher_integrative_signature",     "Myc signature (Felsher)",        0.77,
  "PROLIF_E2F_HALLMARK",                   "E2F targets",                    0.70,
  "GS_METAB_GLUTATHIONE",                  "glutathione (redox)",            0.63)
stopifnot(all(lab_spec$set %in% dat$set))
lab_df <- dplyr::inner_join(lab_spec, dat, by = "set") |>
  dplyr::arrange(rank) |>
  dplyr::mutate(x = nrow(dat) - 5, nx = x - rank, ny = y - abs_load)
# rank order must equal label order, or the fan crosses
stopifnot(!is.unsorted(lab_df$rank), !is.unsorted(rev(lab_df$y)))

# =============================================================================
# THE PANEL
# =============================================================================
# The curve is concave, so the two empty corners are upper-right and lower-left:
# labels go in the first, the key in the second. Nothing else is placed by hand.
p <- ggplot2::ggplot(dat, ggplot2::aes(rank, abs_load, colour = class3)) +
  ggplot2::geom_hline(yintercept = med_ref, linewidth = 0.22,
                      linetype = "22", colour = "grey45") +
  ggplot2::geom_hline(yintercept = redox_ref, linewidth = 0.22,
                      linetype = "13", colour = "grey25") +
  ggplot2::geom_point(size = 0.5, alpha = 0.7, stroke = 0) +
  ggrepel::geom_text_repel(
    data = lab_df, ggplot2::aes(rank, abs_load, label = label),
    inherit.aes = FALSE, size = 1.75, colour = "grey15", segment.colour = "grey55",
    segment.size = 0.15, min.segment.length = 0, box.padding = 0.10,
    point.padding = 0.10, nudge_x = lab_df$nx, nudge_y = lab_df$ny, hjust = 1,
    force = 0, force_pull = 0, max.overlaps = Inf, seed = 1) +
  ggplot2::scale_colour_manual(values = mito_class_cols, labels = mito_class_labels,
                               breaks = names(mito_class_cols)) +
  ggplot2::scale_x_continuous(
    breaks = c(1, seq(200, nrow(dat), by = 200), nrow(dat)),
    expand = ggplot2::expansion(mult = c(0.02, 0.03))) +
  ggplot2::scale_y_continuous(limits = c(0, 1.09), breaks = seq(0, 1, by = 0.25),
                              expand = ggplot2::expansion(mult = c(0.01, 0))) +
  ggplot2::labs(x = "gene sets, ranked by loading",
                y = "absolute loading on PC1", colour = NULL) +
  ggplot2::guides(colour = ggplot2::guide_legend(
    ncol = 1, override.aes = list(size = 1.5, alpha = 1))) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    panel.grid            = ggplot2::element_blank(),
    legend.position       = "inside",
    legend.position.inside = c(0.015, 0.02),
    legend.justification  = c(0, 0),
    legend.background     = ggplot2::element_blank(),
    legend.key            = ggplot2::element_blank(),
    legend.key.size       = ggplot2::unit(2.4, "mm"),
    legend.spacing.y      = ggplot2::unit(0.2, "mm"))

# --- the legend text (never drawn) -------------------------------------------
me <- ax$ref$mito_enrichment
n_med <- function(x) sprintf("%.2f", x)
LEGEND <- panel_legend(
  slot = "Fig. S1C",
  what = paste0(
    "All ", nrow(dat), " scored gene sets ranked by the absolute value of their ",
    "loading on the dominant axis of Fig. 1C, coloured by how the set is defined ",
    "on mitochondrial genes."),
  detail = c(
    "Loading = the correlation of a set's per-sample score with the PC1 sample scores of Fig. 1C, so it is on a correlation scale and 1 means the set moves exactly with the axis.",
    sprintf("Dashed line = the median loading over all sets (%s), i.e. where an arbitrary library set sits. Dotted line = the redox control (%s), the median of GS_METAB_REDOX, GS_METAB_GLUTATHIONE and GS_METAB_REACTIVE_OXYGEN.",
            n_med(med_ref), n_med(redox_ref)),
    sprintf("Mitochondrially defined sets are %.0f per cent of the library, %.0f per cent of the top 50, %.0f per cent of the top 100 and %.0f per cent of the top 200. Median |loading| is %s for MitoCarta and respiratory sets against %s for non-mitochondrial ones (Wilcoxon p < 1e-300).",
            100 * me$base_rate_mito, 100 * me$frac_mito_top50, 100 * me$frac_mito_top100,
            100 * me$frac_mito_top200,
            n_med(stats::median(dat$abs_load[dat$class3 == "mitocarta_proper"])),
            n_med(stats::median(dat$abs_load[dat$class3 == "non_mito"]))),
    "As composites, OXPHOS loads +0.98, Myc signatures +0.97, mitochondrial biogenesis +0.95, proliferation +0.90 and redox +0.21. OXPHOS does not merely participate in the axis; on these samples it is very nearly the axis.",
    paste0("Labelled sets, in rank order, with their exact names: ",
           paste(sprintf("%s = %s (rank %d, %s)",
                         lab_df$label[order(lab_df$rank)],
                         lab_df$set[order(lab_df$rank)],
                         lab_df$rank[order(lab_df$rank)],
                         n_med(lab_df$loading[order(lab_df$rank)])),
                 collapse = "; "),
           ". They are named exemplars chosen to span the curve, not the top of the ranking."),
    sprintf("The classification is the same three classes as Fig. S1B: %d MitoCarta or respiratory sets, %d mitochondrial by construction, %d non-mitochondrial.",
            sum(dat$class3 == "mitocarta_proper"),
            sum(dat$class3 == "construction_MITO"),
            sum(dat$class3 == "non_mito"))),
  bounds = c(
    sprintf("THE TOP OF THE RANKING IS A BUILD TAUTOLOGY. The top 100 splits %d mitochondrial by construction, %d genuine MitoCarta and %d non-mitochondrial. The Gray _MITO transcription-factor lanes are that factor's target programme intersected with MitoCarta, so they load by mitochondrial gene content rather than by regulator identity, and the same transcription factor's non-mitochondrial lane loads far lower. They are drawn rather than deleted, because removing them would flatter the result.",
            me$top100_construction_mito, me$top100_mitocarta_proper, me$top100_non_mito),
    "MITO-LED IS NOT MITO-SPECIFIC. Non-mitochondrial growth programmes reach 0.98 - the human breast-cancer biclusters, nucleotide metabolism, the pentose-phosphate pathway, and the Myc and E2F regulons. And redox, a mitochondrial-function gene set, sits in the bottom half at 0.67; it is also the one axis Myc does not drive. What the axis tracks is a coordinated Myc anabolic-proliferative-mitochondrial growth state, not mitochondria as such.",
    "LOADING IS COVARIATION, NOT PRIMACY. A high loading says a programme moves with the dominant axis; it cannot say the mitochondrial change drives the others, and at n = 24 nothing here can. The primacy claim belongs to the perturbation experiments.",
    "This is the RAW axis of Fig. 1C, which contains the four-group design. Script 37 ranks the same sets on the design-residual axis - the within-group phenotype - and the two rankings correlate 0.968, asserted in pathway_axis().",
    "The classification is a naming and category rule, so the fresh MSigDB Hallmark comparators used elsewhere in the section are not part of this ranking at all; every set here is a library set scored by GSVA."),
  source = c(
    "results/gsva_scores.rds (script 15); results/pathway_loading.rds (script 37) - mito_classification, mito_enrichment",
    "The classification rule: scripts/37_pathway_loading_and_technical_resolution.R:445-462, shared as mito_class3() in figures/panels/_panel_common.R",
    "The construction caveat: docs/library_reference/Gray_et_al_developmental_TFS_selection.md sec 2c; the verdict it comes from: pathway_loading.rds$enrichment_verdict",
    "Composite loadings and the redox control: scripts/35_ambient_corrected_couplings.R and pathway_loading.rds$axis_loadings"))

save_panel_p(p, "figS1_axis_loadings", width = fig_w[["single"]], height = 56)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)

  ## the decile composition, script 37's second view: the enrichment is a
  ## gradient, not a cliff
  dec <- dat |>
    dplyr::mutate(decile = dplyr::ntile(rank, 10)) |>
    dplyr::count(decile, class3) |>
    dplyr::group_by(decile) |> dplyr::mutate(frac = n / sum(n)) |> dplyr::ungroup()
  ggplot2::ggplot(dec, ggplot2::aes(factor(decile), frac, fill = class3)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = mito_class_cols, labels = mito_class_labels) +
    theme_panel()

  ## the top non-mitochondrial loaders, which is where the bound lives
  dat |> dplyr::filter(class3 == "non_mito") |> head(20) |> print()

  ## the construction gap for one transcription factor: same regulator, two lanes
  dat[grep("^TFT_MYC_GRAY_", dat$set), c("rank", "set", "loading", "class3")]

  ## does the ranking change on script 37's design-residual axis?
  ref <- ax$ref$loading_movement
  plot(abs(ref$lambda_design), abs(ax$loading[ref$set]), pch = 16, cex = 0.4,
       xlab = "design-residual", ylab = "raw axis")
}
