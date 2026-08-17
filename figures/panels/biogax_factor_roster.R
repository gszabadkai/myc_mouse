# =============================================================================
# biogax_factor_roster.R -- nothing that could have lowered the chain went down
# -----------------------------------------------------------------------------
# DISCUSSION PANEL (see biogax_regulon_split.R for why the `biogax_` prefix).
#
# THE ARGUMENT. If the maturing gland turned the respiratory chain down through a
# transcription factor, the obvious place to look is the factor. Fifty-seven of
# them -- the PGC-1 family, the ERR/NRF1/GABP axis, the mtDNA machinery, the
# coregulators, the growth arm and the quiescence arm -- are plotted against the
# only null that means anything for a single gene at this n: where the gene's own
# fold change falls among expression-matched genes.
#
# A raw log2 fold change cannot carry this. `Pprc1` moves by -0.47 with an
# adjusted p of 1.9e-4, which looks decisive until you notice it sits at
# baseMean 2102, where a modest change is measured precisely enough to clear FDR.
# Against genes of its own abundance it is at the 7th percentile -- inside the
# band. The percentile is the honest axis; the p-value is drawn as a ring.
#
# WHAT THE PANEL SHOWS: on the DOWN side the band is empty. Nothing in the roster
# falls more than expression-matched genes fall. Everything that clears the band
# goes UP, and the two that clear it with an adjusted p below 0.05 are FOXO3 and
# SIRT1 -- the quiescence arm, not the biogenesis arm.
#
# Reads : results/biogenesis_axis_developmental.rds (script 47 PART D)
# Output: outputs/figures/panels/biogax_factor_roster.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

ba_path <- here::here("results", "biogenesis_axis_developmental.rds")
if (!file.exists(ba_path)) stop("run scripts/47_... first")
ba <- readRDS(ba_path)

ro <- as.data.frame(ba$roster)
stopifnot(all(c("sym", "baseMean", "wt_lfc", "wt_padj", "wt_pct") %in% names(ro)))
ro <- ro[!is.na(ro$wt_pct) & !is.na(ro$wt_lfc), ]

BAND <- c(5, 95)
ro$clears <- ro$wt_pct < BAND[1] | ro$wt_pct > BAND[2]
ro$sig    <- !is.na(ro$wt_padj) & ro$wt_padj < 0.05

# The two families that have to be findable on the page whatever they did: the
# biogenesis axis (the hypothesis under test) and whatever clears the band.
AXIS <- c("Ppargc1a", "Ppargc1b", "Pprc1", "Esrra", "Nrf1", "Gabpa", "Tfam")
ro$family <- ifelse(ro$sym %in% AXIS, "biogenesis axis",
                    ifelse(ro$clears, "clears the band", "other factors"))
lab <- ro[ro$family != "other factors", ]

# ASSERTIONS -- the panel's two readings, so a re-run cannot soften them silently.
stopifnot(
  sum(ro$wt_pct < BAND[1]) == 0 || all(!ro$sig[ro$wt_pct < BAND[1]]),
  ro$wt_pct[ro$sym == "Foxo3"] > BAND[2],
  ro$wt_pct[ro$sym == "Sirt1"] > BAND[2],
  ro$wt_pct[ro$sym == "Pprc1"] > BAND[1],           # inside the band, the point
  ro$baseMean[ro$sym == "Ppargc1a"] < 50)

fam_cols <- c("biogenesis axis" = unname(verdict_cols[["withdraws"]]),
              "clears the band" = unname(verdict_cols[["rises"]]),
              "other factors"   = "grey70")

set.seed(3)
p <- ggplot2::ggplot(ro, ggplot2::aes(x = wt_lfc, y = wt_pct)) +
  ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                    ymin = BAND[1], ymax = BAND[2], fill = "grey96", colour = NA) +
  ggplot2::geom_hline(yintercept = BAND, linewidth = 0.22, colour = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey80") +
  ggplot2::geom_point(ggplot2::aes(fill = family, size = log10(baseMean)),
                      shape = 21, stroke = 0.2, colour = "white", alpha = 0.95) +
  ggplot2::geom_point(data = ro[ro$sig, ], shape = 21, size = 2.6, fill = NA,
                      stroke = 0.35, colour = unname(sig_cols[["sig"]])) +
  # ALWAYS draw the leader (min.segment.length = 0). With repulsion this strong a
  # label can end up nearer a neighbour's point than its own, and an unlabelled
  # significant point sitting between the two makes the misreading easy.
  ggrepel::geom_text_repel(
    data = lab, ggplot2::aes(label = sym, colour = sig),
    size = 1.9, seed = 3, min.segment.length = 0, segment.size = 0.18,
    segment.colour = "grey55", box.padding = 0.2, max.overlaps = 30,
    show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = c(`TRUE` = unname(sig_cols[["sig"]]),
                                          `FALSE` = "grey25"), guide = "none") +
  ggplot2::scale_fill_manual(values = fam_cols, name = NULL) +
  ggplot2::scale_size_continuous(range = c(0.9, 2.4), guide = "none") +
  ggplot2::scale_y_continuous(
    name = "percentile among expression-matched genes",
    limits = c(0, 100), breaks = c(0, 5, 25, 50, 75, 95)) +
  ggplot2::scale_x_continuous(name = "log2 fold change, wild type 6 to 12 weeks",
                              labels = lab_signed) +
  theme_panel() +
  ggplot2::theme(legend.position = "top",
                 legend.margin = ggplot2::margin(0, 0, -3, 0))

get1 <- function(s, col) ro[[col]][ro$sym == s]
LEGEND <- panel_legend(
  slot = "Discussion D4",
  what = paste(
    "Fifty-seven transcription factors, coactivators and corepressors: their",
    "wild-type 6-to-12-week log2 fold change against the percentile that change",
    "occupies among genes of the same expression. The shaded band is the middle",
    "90 per cent. Red rings mark an adjusted p below 0.05; point area is",
    "log10 mean expression."),
  detail = c(
    sprintf("PGC-1a is not expressed in these cells: baseMean %.0f, against PRC's %.0f. Whatever the gland uses, it is not PGC-1a.",
            get1("Ppargc1a", "baseMean"), get1("Pprc1", "baseMean")),
    sprintf("PRC (Pprc1) is the only PGC-1 family member that moves -- %+.3f log2, adjusted p %.1e -- but it sits at the %.1f percentile, INSIDE the band. The small p is precision at high expression, not a large effect.",
            get1("Pprc1", "wt_lfc"), get1("Pprc1", "wt_padj"), get1("Pprc1", "wt_pct")),
    sprintf("The biogenesis axis is unremarkable throughout: Esrra %.1f pct, Nrf1 %.1f, Gabpa %.1f, Tfam %.1f.",
            get1("Esrra", "wt_pct"), get1("Nrf1", "wt_pct"),
            get1("Gabpa", "wt_pct"), get1("Tfam", "wt_pct")),
    sprintf("Exactly one gene clears the band downward -- %s at the %.1f percentile -- and it is not significant (padj %.2f) and not a mitochondrial factor. Upward and significant: Foxo3 %+.3f (%.1f pct, padj %.4f) and Sirt1 %+.3f (%.1f pct, padj %.3f).",
            paste(ro$sym[ro$wt_pct < BAND[1]], collapse = ", "),
            min(ro$wt_pct), ro$wt_padj[which.min(ro$wt_pct)],
            get1("Foxo3", "wt_lfc"), get1("Foxo3", "wt_pct"), get1("Foxo3", "wt_padj"),
            get1("Sirt1", "wt_lfc"), get1("Sirt1", "wt_pct"), get1("Sirt1", "wt_padj"))),
  bounds = c(
    "mRNA is not activity. FOXO3 in particular is controlled by nuclear exclusion, so a transcript change is concordant evidence at best -- which is why the target-set panel is drawn beside this one.",
    "BATCH = TIMEPOINT: the between-age contrast is confounded with cohort. The expression-matched null controls for abundance, not for batch.",
    "The band is the pre-set 5th-to-95th percentile bar; it is not a significance threshold and no multiplicity correction is applied across the 57 genes.",
    "n = 12 wild-type animals."),
  source = "results/biogenesis_axis_developmental.rds (script 47 PART D)")

save_panel_p(p, "biogax_factor_roster", height = 76)

if (FALSE) {
  print(p)
  ba$pgc1_family |> print()
  ro[order(ro$wt_pct), c("sym", "baseMean", "wt_lfc", "wt_padj", "wt_pct")] |>
    head(12) |> print()
  ro[order(-ro$wt_pct), c("sym", "baseMean", "wt_lfc", "wt_padj", "wt_pct")] |>
    head(12) |> print()
}
