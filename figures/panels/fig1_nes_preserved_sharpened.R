# =============================================================================
# fig1_nes_preserved_sharpened.R -- the ranking at twelve weeks against the
# ranking at six
# -----------------------------------------------------------------------------
# SLOT: Fig. 1H.
#
#   "Interestingly, this ranking was not just maintained but enhanced: the
#    normalized enrichment score for the mitochondrial OXPHOS set and MYC
#    integrative signatures increased substantially (Fig. 1H)."
#
# MAINTAINED IS STRONGLY TRUE and is the whole panel: 866 gene sets, Spearman
# 0.933. The effect size halves (Fig. 1G) and the order does not move.
#
# ENHANCED NEEDS THE SENTENCE CHANGED, and the panel is drawn so the reader can
# see why. Every point sits above the identity line, not just the named ones:
#
#     all 866 sets                 81% rise, median +0.312
#     already-enriched (NES6 >= 2) 93% rise, median +0.319
#     OXPHOS programme (17 sets)  100% rise, median +0.308 -> the 49th percentile
#     Myc signatures  (16 sets)    94% rise, median +0.386 -> the 60th percentile
#
# So as FAMILIES the two the sentence names rise by exactly the median amount.
# The individual sets the earlier draft quoted are the best members of those
# families (MITOCARTA_OXPHOS +0.52, 76th percentile; MYC_felsher +0.55, 78th) --
# upper quartile, not exceptional, and not what "increased substantially" claims.
#
# AND THERE IS A MECHANICAL REASON TO EXPECT A GLOBAL NES RISE HERE. fGSEA
# normalises each enrichment score against a permutation null built from THAT
# ranking. The twelve-week Wald list is much flatter than the six-week one --
# SD 1.221 against 1.744, IQR 1.503 against 2.176, 9.5% of genes past |stat| = 2
# against 22.5% -- so random sets reach smaller enrichment scores, the normaliser
# shrinks, and the same relative enrichment scores a higher NES. A rise in NES is
# what a WEAKER ranking mechanically produces, which is the opposite of what
# "enhanced" is meant to convey.
#
# WHAT SURVIVES, AND IT IS WORTH SAYING: NES is enrichment RELATIVE to the rest of
# the transcriptome, so a mitochondrial or Myc set sitting higher at twelve weeks
# means the residual Myc signal is more concentrated on those programmes than it
# was. The programme holds its position while the amplitude halves -- which is the
# same statement as the retention asymmetry (Myc core 0.74-0.77 against the
# mitochondrial arms' 0.50-0.58), read on a different instrument. That is a
# statement about POSITION, not about magnitude, and the panel supports it.
#
# THE TWO FAMILIES ARE PICKED BY programme_group(), the declared encoding Fig. 1D
# groups its rows by, so "OXPHOS" and "Myc signatures" mean the same thing on both
# panels and neither is a list hand-assembled here.
#
# Reads (read-only, no re-run):
#   results/fgsea_percategory.rds (script 20) -- $fgsea, NES per ranking x category
# Output: outputs/figures/panels/fig1_nes_preserved_sharpened.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

fg_path <- here::here("results", "fgsea_percategory.rds")
require_fresher_than(fg_path)
fg <- readRDS(fg_path)$fgsea
stopifnot(all(c("ranking", "category", "pathway", "NES") %in% names(fg)))

w <- fg |>
  dplyr::filter(ranking %in% contrast_geno) |>
  dplyr::select(ranking, category, pathway, NES) |>
  tidyr::pivot_wider(names_from = ranking, values_from = NES) |>
  as.data.frame()
names(w)[match(contrast_geno, names(w))] <- c("n6", "n12")
w$d    <- w$n12 - w$n6
w$prog <- as.character(programme_group(w$pathway, w$category))

# The two families the sentence names, by the declared grouping.
FAM <- c(OXPHOS = "OXPHOS", `Myc signatures` = "Myc signatures")
w$fam <- ifelse(w$prog %in% FAM, w$prog, "other")
w$fam <- factor(w$fam, levels = c("other", unname(FAM)))

rho    <- stats::cor(w$n6, w$n12, method = "spearman")
n_sets <- nrow(w)
stopifnot(n_sets == 866L, !anyNA(w$n6), !anyNA(w$n12),
          sum(w$fam == "OXPHOS") == 17L, sum(w$fam == "Myc signatures") == 16L)

# --- the family colours ------------------------------------------------------
# The two poles of the declared manuscript ramp, used here as two categories.
# Neither is a sample colour, and the mint is the same mint Fig. 1G spends on
# OXPHOS, so a reader who has met one panel has met the other. `ms_diverging` is
# a ramp by design and this is the one place it is used categorically -- said in
# the legend block rather than left to be noticed.
fam_cols <- c("other" = "grey72",
              "OXPHOS" = unname(ms_diverging[["pos"]]),
              "Myc signatures" = unname(ms_diverging[["neg"]]))

RNG <- range(c(w$n6, w$n12))
RNG <- RNG + c(-1, 1) * diff(RNG) * 0.04

p <- ggplot2::ggplot(w, ggplot2::aes(n6, n12)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey88") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey88") +
  # The identity line, which is the whole reference: above it the set is more
  # enriched at twelve weeks than at six.
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "22",
                       linewidth = 0.3, colour = "grey45") +
  ggplot2::geom_point(data = w[w$fam == "other", ], colour = fam_cols[["other"]],
                      size = 0.4, alpha = 0.5, stroke = 0) +
  ggplot2::geom_point(data = w[w$fam != "other", ],
                      ggplot2::aes(colour = fam), size = 1.0, alpha = 0.95,
                      stroke = 0) +
  ggplot2::annotate("text",
                    x = RNG[1] + diff(RNG) * 0.025, y = RNG[2] - diff(RNG) * 0.02,
                    label = sprintf('rho~"%.2f"', rho),
                    parse = TRUE, hjust = 0, vjust = 1, size = 1.8, colour = "grey25") +
  ggplot2::scale_colour_manual(values = fam_cols, breaks = unname(FAM), name = NULL) +
  ggplot2::scale_x_continuous(limits = RNG, labels = lab_signed,
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::scale_y_continuous(limits = RNG, labels = lab_signed,
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::labs(x = "NES, Myc effect at 6W", y = "NES, Myc effect at 12W") +
  ggplot2::guides(colour = ggplot2::guide_legend(
    override.aes = list(size = 1.5, alpha = 1))) +
  theme_panel(base_size = 6) +
  # Key inside, bottom right -- the wedge nothing occupies, because a set strongly
  # enriched at six weeks and depleted at twelve is exactly what does not happen.
  ggplot2::theme(
    legend.position        = "inside",
    legend.position.inside = c(0.99, 0.02),
    legend.justification   = c(1, 0),
    legend.background      = ggplot2::element_blank(),
    legend.margin          = ggplot2::margin(0, 0, 0, 0),
    legend.key.size        = ggplot2::unit(2.4, "mm"),
    plot.margin            = ggplot2::margin(1.5, 2, 1, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
fam_line <- function(g) {
  s <- w[w$fam == g, ]
  sprintf("%s (%d sets): median NES %+.2f at six weeks and %+.2f at twelve, a median rise of %+.3f, which is the %.0fth percentile of the rises of all %d sets; %.0f%% of them rise.",
          g, nrow(s), stats::median(s$n6), stats::median(s$n12), stats::median(s$d),
          100 * mean(w$d <= stats::median(s$d)), n_sets, 100 * mean(s$d > 0))
}
top <- function(pw) {
  r <- w[w$pathway == pw, ]
  sprintf("%s %+.2f to %+.2f (rise %+.2f, %.0fth percentile)", pw, r$n6, r$n12, r$d,
          100 * mean(w$d <= r$d))
}
enr <- w[w$n6 >= 2, ]
# The two off-diagonal clouds, defined by crossing sign with room to spare, and
# what the Spearman becomes once they are set aside.
cross    <- w[(w$n6 < -1 & w$n12 > 0.5) | (w$n6 > 0.5 & w$n12 < -0.5), ]
core     <- w[!(w$pathway %in% cross$pathway), ]
rho_core <- stats::cor(core$n6, core$n12, method = "spearman")

LEGEND <- panel_legend(
  slot = "Fig. 1H",
  what = paste0(
    "The fGSEA normalised enrichment score of every gene set in the library for ",
    "the Myc effect at twelve weeks against the same score at six weeks. One ",
    "point per set; the dashed line is equality. The two programmes the text ",
    "names are picked out by the same grouping Fig. 1D uses for its rows."),
  detail = c(
    sprintf("n = %d gene sets, ranked on the unshrunken Wald statistic of each genotype contrast. Both contrasts are clean: genotype is balanced within each extraction batch.",
            n_sets),
    sprintf("THE ORDER IS PRESERVED, which is the panel's strong statement: Spearman %.3f (Pearson %.3f) across all %d sets, while the underlying effect size falls to about half (Fig. 1G).",
            rho, stats::cor(w$n6, w$n12), n_sets),
    fam_line("OXPHOS"), fam_line("Myc signatures"),
    sprintf("The individual sets an earlier draft quoted are the best members of those families, not typical of them: %s; %s; %s.",
            top("MITOCARTA_OXPHOS"), top("MYC_felsher_integrative_signature"),
            top("HALLMARK_MYC_TARGETS_V1")),
    sprintf("Enrichment scores rise almost everywhere, which is why the whole cloud sits above the line: %.0f%% of all %d sets rise, by a median of %+.3f, and among the %d sets already enriched at six weeks (NES >= 2) %.0f%% rise by a median of %+.3f.",
            100 * mean(w$d > 0), n_sets, stats::median(w$d), nrow(enr),
            100 * mean(enr$d > 0), stats::median(enr$d)),
    "The largest movements on the panel are in the upper left, and they are not a Myc result: several mammary-development sets go from strongly depleted at six weeks to enriched at twelve (MG_LHS_CONSENSUS -3.20 to +2.05, MG_LHOR_SAEKI -3.45 to +2.00). That is the release of Myc's suppression of luminal hormone-sensing identity, the largest single departure from the dose line in the transcriptome, and it is a beat the Results has not yet placed.",
    sprintf("THE RESHUFFLING IS NOT SPREAD EVENLY, and knowing where it sits sharpens the claim. The %d sets that cross sign between the two ages -- the two off-diagonal clouds -- are %d%% mammary-development and transcription-factor-target sets. Set them aside and the Spearman over the remaining %d rises from %.3f to %.3f. What holds its order is the metabolic, mitochondrial and Myc core; what reshuffles is lineage identity, which is the subject of section 2.",
            nrow(cross), round(100 * mean(cross$category %in%
                                            c("03_mammary_development", "06_tf_targets"))),
            n_sets - nrow(cross), rho, rho_core)),
  bounds = c(
    "\"ENHANCED\" IS NOT SPECIFIC, AND THE SENTENCE SHOULD SAY SO. The two families named rise by the MEDIAN amount for a gene set in this library -- OXPHOS at the 49th percentile of all rises, the Myc signatures at the 60th. What is true and worth writing is that they hold their position at the top of the ranking while the amplitude halves.",
    "AND A GLOBAL RISE IN NES IS WHAT A WEAKER RANKING MECHANICALLY PRODUCES. fGSEA normalises each enrichment score against a permutation null built from the same ranked list. The twelve-week Wald list is much flatter than the six-week one -- SD 1.221 against 1.744, IQR 1.503 against 2.176, and 9.5% of genes past |stat| = 2 against 22.5% -- so random sets reach smaller enrichment scores, the normaliser shrinks and the same relative enrichment scores a higher NES. No claim of increased Myc activity can rest on this panel.",
    "What the panel does support is a statement about POSITION rather than magnitude. NES is enrichment relative to the rest of the transcriptome, so a set sitting higher at twelve weeks means the residual Myc signal is more concentrated on it. That is the retention asymmetry of scripts 29 to 31 (Myc core retains 0.74 to 0.77 of its effect against the mitochondrial arms' 0.50 to 0.58) read on a different instrument.",
    "Benjamini-Hochberg is applied WITHIN each library category by script 20's design, so adjusted p-values are not comparable across the panel and none is drawn. fGSEA's gene permutation is also anti-conservative for sets of correlated genes, which is most of this library.",
    "The 50 MSigDB Hallmark comparators are in the cloud but are not library sets; HALLMARK_MYC_TARGETS_V1 and V2 are therefore NOT in the Myc-signatures family drawn here, which is the library's own 02_myc_signatures category. Both are quoted in the detail above.",
    "`ms_diverging` is a diverging RAMP everywhere else in this manuscript; this is the one panel that spends its two poles as categorical colours. They were chosen because neither is a sample colour and because the mint is the same mint Fig. 1G gives OXPHOS."),
  source = c(
    "results/fgsea_percategory.rds (scripts/20_fgsea_percategory.R) -- one fGSEA per ranking x category on the unshrunken Wald statistic",
    "The family grouping: programme_group() in figures/panels/_panel_common.R, the same encoding Fig. 1D groups its rows by",
    "The Wald-statistic spreads quoted in the bounds: results/interaction_results.rds, myc_6W_raw and myc_12W_raw"))

save_panel_p(p, "fig1_nes_preserved_sharpened", height = 62)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the two families, set by set
  w[w$fam != "other", c("pathway", "fam", "n6", "n12", "d")] |>
    (\(x) x[order(x$fam, -x$d), ])() |> print(row.names = FALSE, digits = 3)

  ## where every programme sits, and how far it moves
  aggregate(cbind(n6, n12, d) ~ prog, data = w, FUN = stats::median) |>
    (\(x) x[order(-x$n6), ])() |> print(row.names = FALSE, digits = 3)

  ## the rise, as a distribution -- the thing the sentence has to reckon with
  quantile(w$d, c(0.1, 0.25, 0.5, 0.75, 0.9)) |> round(3) |> print()

  ## the flatter twelve-week ranking, which is the mechanical explanation
  ir <- readRDS(here::here("results", "interaction_results.rds"))
  vapply(c("myc_6W_raw", "myc_12W_raw"), function(n) {
    s <- as.data.frame(ir[[n]])$stat; s <- s[is.finite(s)]
    c(sd = stats::sd(s), IQR = stats::IQR(s), frac_gt2 = mean(abs(s) > 2))
  }, numeric(3)) |> round(3) |> print()

  ## the upper-left movers -- the luminal release, unplaced in the narrative
  w[order(-w$d), c("pathway", "category", "n6", "n12", "d")] |> head(10) |>
    print(row.names = FALSE, digits = 3)
}
