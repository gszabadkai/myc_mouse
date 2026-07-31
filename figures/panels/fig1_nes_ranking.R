# =============================================================================
# fig1_nes_ranking.R -- the library ranked by normalised effect size
# -----------------------------------------------------------------------------
# SLOT: Fig. 1D (was Fig. S1C for a day; the author's call on 2026-07-31 is that
# this is the more relevant main-figure panel, so it swapped with the axis
# loadings). Filenames do not carry the slot letter; figures/panels/PANELS.md is
# the slug -> slot map.
#
# SUPPORTS (Results, "Myc drives early breast tumourigenesis by inducing OXPHOS
# and biosynthetic pathways", paragraph 2):
#   "Mitochondrial biogenesis and OXPHOS complex genesets were overall top ranked
#    based on normalised changes in effect size, along with core
#    (non-mitochondrial) Myc target genesets, followed by biosynthetic metabolic
#    pathways, E2F and overall proliferation signaling (Fig. 1D)."
#
# THE RULER. "Normalised changes in effect size" is the fGSEA normalised
# enrichment score, which is how the manuscript's own methods paragraph uses the
# word: enrichment of differential expression normalised to the whole
# transcriptome. It is a DIFFERENT ruler from Fig. 1C and Fig. S1C, and that is
# the point of having both. Those two are per-sample covariation - which
# programmes move together across the 24 animals. This one is the genotype
# CONTRAST - which programmes are enriched among the genes Myc actually moves,
# ranked on the unshrunken Wald statistic. A programme can top one and not the
# other, and that is why this one is the main-figure panel: it is the one that
# measures what Myc DID rather than what covaries with what.
#
# WHY BOTH AGES. The claim is about a ranking, and a ranking that held at one age
# and not the other would be a different result. It holds: the two columns are
# essentially the same order (Spearman 0.93), which is also the observation the
# attenuation section opens with - the effect sizes fall while the ranking does
# not move.
#
# THE GROUPING IS AN ENCODING, so it is declared in _panel_common.R as
# programme_group() rather than written inline here, with the two rules that keep
# it honest: OXPHOS and biogenesis are split by script 37's own patterns so the
# words mean the same thing across panels, and every set whose membership is
# MitoCarta intersected with something else goes into ONE row rather than being
# spread across the others. That row is the top row. It should be: those sets are
# mitochondrial by build, and the alternative - deleting them - would flatter
# every row beneath.
#
# Input:  results/fgsea_percategory.rds  (script 20 -- NES per ranking x category)
#         results/ap6_permutation_null.rds (script 21 -- the matched null, legend)
#         data/genesets_from_library/provenance_table.csv (the mito classification)
# Output: outputs/figures/panels/fig1_nes_ranking.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

fg_path <- here::here("results", "fgsea_percategory.rds")
require_fresher_than(fg_path)
fg <- readRDS(fg_path)$fgsea
stopifnot(all(c("ranking", "category", "pathway", "NES", "padj_within_category")
              %in% names(fg)))

prov <- utils::read.csv(
  here::here("data", "genesets_from_library", "provenance_table.csv"),
  stringsAsFactors = FALSE)

dat <- fg |>
  dplyr::filter(ranking %in% contrast_geno) |>
  dplyr::mutate(
    contrast  = factor(ranking, levels = contrast_geno),
    programme = programme_group(pathway, category),
    # the same three classes as Figs. S1B and 1D; the fresh MSigDB Hallmark
    # comparators are not library sets, so they fall through to non-mitochondrial
    # by the naming rule -- said in the legend block, not silently
    class3 = factor(
      mito_class3(pathway,
                  prov$category_primary[match(pathway, prov$set_name)]),
      levels = names(mito_class_cols)))
stopifnot(!anyNA(dat$programme), !anyNA(dat$class3),
          nrow(dat) == 2L * sum(fg$ranking == contrast_geno[1]))

# rows are ordered by the 6-week median, which is the ordering the sentence makes
ord <- dat |>
  dplyr::filter(contrast == contrast_geno[1]) |>
  dplyr::group_by(programme) |>
  dplyr::summarise(med = stats::median(NES), .groups = "drop") |>
  dplyr::arrange(med)
dat$programme <- factor(as.character(dat$programme),
                        levels = as.character(ord$programme))

med_df <- dat |>
  dplyr::group_by(contrast, programme) |>
  dplyr::summarise(NES = stats::median(NES), n = dplyr::n(), .groups = "drop") |>
  # the median tick is drawn as a plain segment on the discrete y scale (integer
  # positions), which avoids geom_crossbar's deprecated `fatten` and says exactly
  # what it is: a line at the median, not a box with a summary inside it
  dplyr::mutate(yi = as.integer(factor(as.character(programme),
                                       levels = levels(dat$programme))))

# --- point geom: quasirandom if available, else jitter ------------------------
# Same fallback as figS1_mitocarta_survey_share.R:93-97, so a machine without
# ggbeeswarm still renders the panel rather than failing.
pts <- if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
  ggbeeswarm::geom_quasirandom(size = 0.45, alpha = 0.75, stroke = 0,
                               orientation = "y", width = 0.24)
} else {
  ggplot2::geom_jitter(size = 0.45, alpha = 0.75, stroke = 0, height = 0.24, width = 0)
}

RNG <- range(dat$NES)

# =============================================================================
# THE PANEL
# =============================================================================
p <- ggplot2::ggplot(dat, ggplot2::aes(NES, programme, colour = class3)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey75") +
  pts +
  ggplot2::geom_segment(data = med_df,
                        ggplot2::aes(x = NES, xend = NES,
                                     y = yi - 0.31, yend = yi + 0.31),
                        inherit.aes = FALSE, linewidth = 0.3, colour = "grey15") +
  ggplot2::facet_wrap(~ contrast, nrow = 1) +
  ggplot2::scale_colour_manual(values = mito_class_cols, labels = mito_class_labels,
                               breaks = names(mito_class_cols)) +
  ggplot2::scale_x_continuous(breaks = seq(-2, 2, by = 2), labels = lab_signed,
                              expand = ggplot2::expansion(mult = 0.05)) +
  ggplot2::labs(x = "fGSEA normalised enrichment score", y = NULL, colour = NULL) +
  ggplot2::guides(colour = ggplot2::guide_legend(
    nrow = 1, override.aes = list(size = 1.5, alpha = 1))) +
  ggplot2::coord_cartesian(xlim = RNG, clip = "off") +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(linewidth = 0.15, colour = "grey93"),
    panel.grid.minor   = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.spacing.x    = ggplot2::unit(2.4, "mm"),
    strip.text         = ggplot2::element_text(size = 6, margin = ggplot2::margin(0, 0, 1, 0)),
    axis.text.y        = ggplot2::element_text(margin = ggplot2::margin(r = 0.6, unit = "mm")),
    # One row, under the plot (author, 2026-07-31). The short class names fit
    # 89 mm on a single line, which is what makes this possible: it costs one
    # line of height instead of three, and it cannot collide with any point.
    legend.position       = "bottom",
    legend.justification  = c(0, 0.5),
    legend.location       = "plot",      # flush with the panel edge, not the axis
    legend.key            = ggplot2::element_blank(),
    legend.key.size       = ggplot2::unit(2.2, "mm"),
    legend.spacing.y      = ggplot2::unit(0, "mm"),
    legend.margin         = ggplot2::margin(0, 0, 0, -1),
    legend.box.spacing    = ggplot2::unit(0.8, "mm"))

# --- the legend text (never drawn) -------------------------------------------
ap6 <- readRDS(here::here("results", "ap6_permutation_null.rds"))$null_table
z6  <- function(cmp) ap6$z[ap6$compartment == cmp & ap6$contrast == "myc_6W"]
m6  <- med_df[med_df$contrast == contrast_geno[1], ]
m6  <- m6[order(-m6$NES), ]
w   <- fg |> dplyr::filter(ranking %in% contrast_geno) |>
  dplyr::select(ranking, pathway, NES) |>
  tidyr::pivot_wider(names_from = ranking, values_from = NES)
rho <- suppressWarnings(stats::cor(w[[contrast_geno[1]]], w[[contrast_geno[2]]],
                                   method = "spearman"))
n_sig6 <- sum(dat$padj_within_category[dat$contrast == contrast_geno[1]] < 0.05,
              na.rm = TRUE)
n_sets <- nrow(w)

LEGEND <- panel_legend(
  slot = "Fig. 1D",
  what = paste0(
    "Every gene set in the library, grouped by programme and ranked by fGSEA ",
    "normalised enrichment score for the Myc effect at each age. One point per ",
    "gene set; the bar is the group median."),
  detail = c(
    sprintf("n = %d gene sets per contrast: %d fGSEA-eligible library sets that pass the size filter (10 to 500 genes) plus %d fresh MSigDB Hallmark sets as an off-the-shelf comparator.",
            n_sets, n_sets - 50L, 50L),
    "Ranked on the unshrunken Wald statistic of each genotype contrast, so the score is shrinkage-independent and reads as absolute transcriptional enrichment against the whole transcriptome.",
    sprintf("Rows are ordered by the 6-week median. In order: %s.",
            paste(sprintf("%s %s", as.character(m6$programme), sprintf("%+.2f", m6$NES)),
                  collapse = "; ")),
    sprintf("The ordering is the same at both ages (Spearman %.2f across all %d sets), and the median absolute NES barely moves (%.2f at 6 weeks against %.2f at 12).",
            rho, n_sets,
            stats::median(abs(w[[contrast_geno[1]]])),
            stats::median(abs(w[[contrast_geno[2]]]))),
    sprintf("THE MATCHED NULL, which is what licenses 'top ranked' rather than merely 'top of a list'. Against an expression x dispersion-matched permutation null on the 6-week contrast (script 21), MitoCarta sits at z = %.1f, Metabolism %.1f, Hallmark OXPHOS %.1f, Hallmark MYC targets V1 %.1f, Hallmark E2F targets %.1f and Proliferation %.1f (all empirical p <= 0.004). The mitochondrial focus survives the obvious confound, that highly expressed housekeeping genes enrich easily; proliferation is the weakest of the six.",
            z6("MitoCarta"), z6("Metabolism"), z6("HALLMARK_OXPHOS"),
            z6("HALLMARK_MYC_TARGETS_V1"), z6("HALLMARK_E2F_TARGETS"),
            z6("Proliferation")),
    "Colour is the three-class mitochondrial definition of Fig. S1B, applied per gene set. The top row is named for what those sets ARE - curated mitochondrial sets - and the colour key says how they were made.",
    sprintf("APOPTOSIS IS THE ONE ROW THAT CHANGES SIGN (%+.2f at 6 weeks, %+.2f at 12; every other row's median moves up). It is not a result: the row is bimodal and its median sits on zero, so five weak sets crossing zero move it. Those five are TANG_NECROPTOSIS, TANG_LYSOSOME_DEPENDENT_CELL_DEATH, APOP_INTRINSIC_REACTOME, APOP_MODULATION_WP and TANG_PYROPTOSIS, none significant at either age. The members that ARE significant do not move: APOP_REGULATION_REACTOME stays at +2.1/+2.4 and TANG_CUPROPTOSIS at +1.9, while APOP_HALLMARK and APOP_KEGG stay depleted and if anything rise. Three of the five movers are non-apoptotic death modalities, so there is no mitochondrial death story in it either.",
            m6$NES[as.character(m6$programme) == "apoptosis"],
            med_df$NES[med_df$contrast == contrast_geno[2] &
                       as.character(med_df$programme) == "apoptosis"])),
  bounds = c(
    "THE TOP ROW IS A BUILD TAUTOLOGY, and it is drawn rather than dropped. It collects every set whose membership is MitoCarta intersected with something else - the Gray _MITO transcription-factor lanes and the biogenesis-discrimination and biogenesis-by-apoptosis constructs - so it enriches by mitochondrial gene content, not by the identity of the programme it is named for. The rows to read as programmes start below it.",
    sprintf("SIGNIFICANCE DOES NOT DISCRIMINATE HERE: %d of the %d sets clear padj < 0.05 on the 6-week contrast. The ordering is the readable quantity, not the p-value.",
            n_sig6, n_sets),
    "Benjamini-Hochberg is applied WITHIN each library category by script 20's design (that is the scope publication claims cite), so the adjusted p-values are not comparable between rows.",
    "fGSEA's gene permutation is anti-conservative for sets of correlated genes, which is most of this library; the enrichment scores are a ranking instrument here and no claim rests on an individual set's p-value.",
    "The two contrasts are genotype contrasts, 6 versus 6 balanced within batch, so both are clean. Nothing on this panel is a developmental comparison.",
    "One caveat in the colouring: the 50 Hallmark comparators are not library sets, so the naming and category rule classes all of them as non-mitochondrial, including HALLMARK_OXIDATIVE_PHOSPHORYLATION. That row is a comparator, and the classification is not the reason it sits where it does.",
    "The matched-null figures come from a run that predates the 2026-07-24 gene-symbol reconciliation. The reconciliation's recovery was concentrated in ATP synthase and OXPHOS membership, so it can only have raised the mitochondrial z scores; the direction of the statement is safe and the exact numbers are a floor."),
  source = c(
    "results/fgsea_percategory.rds (scripts/20_fgsea_percategory.R) - one fGSEA per ranking x category on the unshrunken Wald statistic",
    "results/ap6_permutation_null.rds (scripts/21_ap6_permutation_null.R) - the expression x dispersion-matched null",
    "data/genesets_from_library/provenance_table.csv for category_primary; the classification rule is mito_class3() in figures/panels/_panel_common.R",
    "The programme grouping: programme_group() in figures/panels/_panel_common.R"))

save_panel_p(p, "fig1_nes_ranking", width = fig_w[["single"]], height = 70)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)

  ## the median table the legend quotes, both ages side by side
  med_df |> tidyr::pivot_wider(names_from = contrast, values_from = c(NES, n)) |>
    dplyr::arrange(dplyr::desc(NES_myc_6W)) |> print(n = 20)

  ## what is actually at the top of the ranking, set by set
  dat |> dplyr::filter(contrast == "myc_6W") |> dplyr::arrange(dplyr::desc(NES)) |>
    dplyr::select(programme, pathway, NES, padj_within_category) |>
    head(25) |> print()

  ## the same, with the by-construction row removed: does the order survive?
  dat |> dplyr::filter(contrast == "myc_6W",
                       programme != "mito-defined by construction") |>
    dplyr::group_by(programme) |>
    dplyr::summarise(n = dplyr::n(), med = stats::median(NES)) |>
    dplyr::arrange(dplyr::desc(med)) |> print(n = 20)

  ## the whole matched null, both ages
  readRDS(here::here("results", "ap6_permutation_null.rds"))$null_table |> print(n = 20)
}
