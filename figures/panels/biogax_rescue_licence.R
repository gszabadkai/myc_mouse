# =============================================================================
# biogax_rescue_licence.R -- what licenses reverting the decline with PGC1a/NRF1
# -----------------------------------------------------------------------------
# DISCUSSION PANEL (see biogax_regulon_split.R for why the `biogax_` prefix).
#
# THE PROBLEM THIS ANSWERS. If the maturing gland did NOT lower its respiratory
# chain by lowering PGC1a-axis activity, why is the experimental part allowed to
# raise the chain back with PGC-1a and NRF1? The answer is not that the axis is a
# plausible mechanism -- it is not the mechanism. It is that the axis is the right
# INSTRUMENT, and the three panels here are the three properties an instrument
# needs, each measured rather than argued.
#
#   TOP     REACH. The axis did not lower these genes but it reaches them: two
#           thirds of the 89 subunits that fell are inside it. Reverting a deficit
#           needs reach over the affected genes, not authorship of the deficit.
#
#   MIDDLE  CONVERGENCE. ERRa and NRF1 regulons overlap in only 22 genes, and half
#           of that overlap is respiratory subunits -- roughly twice the density
#           in either regulon alone. Two interventions whose reach barely overlaps,
#           producing one phenotype, localise the effect to what they share. This
#           is why using BOTH factors is a control and not a redundancy.
#
#   BOTTOM  NON-CIRCULARITY, and this is the plank a reviewer will actually test.
#           If the death effectors sat inside the PGC1a regulon, a restored death
#           phenotype could be direct transcription of the machinery and the
#           experiment would prove nothing. They do not: the regulon holds five of
#           the 25 pro-apoptotic MitoCarta genes, below its expected share, and
#           none of them is PUMA, BAX or a caspase. They sit in MYC's regulon
#           instead, at 2.7 times expectation.
#
# Reads : results/biogenesis_axis_developmental.rds (script 47 PART H)
# Output: outputs/figures/panels/biogax_rescue_licence.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

ba_path <- here::here("results", "biogenesis_axis_developmental.rds")
if (!file.exists(ba_path)) stop("run scripts/47_... first")
ba <- readRDS(ba_path)

rc <- as.data.frame(ba$rescue$reach)
cv <- as.data.frame(ba$rescue$convergence)
ct <- as.data.frame(ba$rescue$convergence_test)
ow <- as.data.frame(ba$rescue$ownership)

# --- TOP: reach over the genes that fell --------------------------------------
rc$label <- sub("_MITO$", "", rc$regulon)
rc$label <- c(CORE = "ERRa + NRF1 + GABP", ESRRA = "ERRa", NRF1 = "NRF1",
              GABPA = "GABP", MYC = "MYC", DEVELOPMENTAL = "ER",
              E2F1 = "E2F1")[rc$label]
rc <- rc[order(rc$frac_of_declining), ]
rc$label <- factor(rc$label, levels = rc$label)
rc$is_axis <- rc$regulon %in% c("CORE_MITO", "ESRRA_MITO", "NRF1_MITO", "GABPA_MITO")
stopifnot(rc$frac_of_declining[rc$regulon == "CORE_MITO"] > 0.6)

p1 <- ggplot2::ggplot(rc, ggplot2::aes(x = 100 * frac_of_declining, y = label)) +
  ggplot2::geom_col(ggplot2::aes(fill = is_axis), width = 0.66) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%d of %d", covers_ox, n_ox_total)),
                     hjust = -0.14, size = 1.8, colour = "grey20") +
  ggplot2::scale_fill_manual(
    values = c(`TRUE` = unname(verdict_cols[["withdraws"]]), `FALSE` = "grey78"),
    guide = "none") +
  ggplot2::scale_x_continuous(name = "per cent of the 89 declining subunits reached",
                              limits = c(0, 88), expand = c(0, 0)) +
  ggplot2::scale_y_discrete(name = NULL) +
  theme_panel() +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

# --- MIDDLE: where the two regulons converge ----------------------------------
cvd <- cv[cv$part != "union", ]
cvd$part <- factor(cvd$part, levels = c("NRF1 only", "ESRRA only", "shared"),
                   labels = c("NRF1 alone", "ERRa alone", "the 22 genes\nthey share"))
cvd$hi <- grepl("share", as.character(cvd$part))
stopifnot(cvd$ox_fraction[cvd$hi] > 2 * min(cvd$ox_fraction), ct$p_enrich < 0.05)

p2 <- ggplot2::ggplot(cvd, ggplot2::aes(x = 100 * ox_fraction, y = part)) +
  ggplot2::geom_col(ggplot2::aes(fill = hi), width = 0.62) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%d of %d genes", n_oxphos, n_genes)),
                     hjust = -0.14, size = 1.8, colour = "grey20") +
  ggplot2::scale_fill_manual(
    values = c(`TRUE` = unname(verdict_cols[["withdraws"]]), `FALSE` = "grey78"),
    guide = "none") +
  ggplot2::scale_x_continuous(name = "per cent of the part that is respiratory subunits",
                              limits = c(0, 74), expand = c(0, 0)) +
  ggplot2::scale_y_discrete(name = NULL) +
  theme_panel() +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(lineheight = 0.95))

# --- BOTTOM: who owns the death effectors -------------------------------------
ows <- ow[ow$apoptosis_set == "MITOCARTA_APOPTOSIS_PRO", ]
ows$label <- c(CORE_MITO = "ERRa + NRF1 + GABP", ESRRA_MITO = "ERRa",
               NRF1_MITO = "NRF1", MYC_MITO = "MYC")[ows$programme]
ows <- ows[order(ows$fold), ]
ows$label <- factor(ows$label, levels = ows$label)
stopifnot(ows$fold[ows$programme == "CORE_MITO"] < 1,
          ows$fold[ows$programme == "MYC_MITO"]  > 2)

p3 <- ggplot2::ggplot(ows, ggplot2::aes(x = fold, y = label)) +
  ggplot2::geom_vline(xintercept = 1, linewidth = 0.25, colour = "grey55") +
  ggplot2::geom_segment(ggplot2::aes(x = 1, xend = fold, yend = label),
                        linewidth = 0.4, colour = "grey70") +
  # ONE INK, deliberately. Elsewhere in this set green means "rises with age";
  # using it here for "enriched over expectation" would put two meanings on one
  # colour. Position carries the effect and the p-values go in the legend, which
  # is _panel_common.R's rule for significance.
  ggplot2::geom_point(shape = 21, size = 2.1, stroke = 0.3,
                      fill = "grey25", colour = "white") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%d of 25", overlap)),
                     hjust = -0.4, size = 1.8, colour = "grey20") +
  ggplot2::scale_x_continuous(
    name = "pro-apoptotic genes, observed / expected",
    limits = c(0, 3.4), breaks = c(0, 1, 2, 3), expand = c(0, 0)) +
  ggplot2::scale_y_discrete(name = NULL) +
  theme_panel() +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

p <- patchwork::wrap_plots(p1, p2, p3, ncol = 1, heights = c(1.15, 0.85, 0.85))

LEGEND <- panel_legend(
  slot = "Discussion D6",
  what = paste(
    "Three properties an intervention needs in order to test a state it did not",
    "create. Top: how much of the declining respiratory chain each regulon",
    "reaches. Middle: where the ERRa and NRF1 regulons overlap, and how",
    "respiratory that overlap is. Bottom: which regulon owns the pro-apoptotic",
    "effectors, as observed over expected against the MitoCarta background."),
  detail = c(
    sprintf("REACH: the ERRa/NRF1/GABP regulon contains %d of the %d subunits that fell (%.0f per cent); ERRa alone contains %d.",
            rc$covers_ox[rc$regulon == "CORE_MITO"], rc$n_ox_total[1],
            100 * rc$frac_of_declining[rc$regulon == "CORE_MITO"],
            rc$covers_ox[rc$regulon == "ESRRA_MITO"]),
    sprintf("CONVERGENCE: the two regulons share only %d genes, of which %d are respiratory subunits (%.0f per cent), against %.0f per cent for ERRa alone and %.0f per cent for NRF1 alone -- a %.2f-fold enrichment over their union, hypergeometric p = %.4f.",
            ct$n_shared, ct$n_ox_in_shared, 100 * cv$ox_fraction[cv$part == "shared"],
            100 * cv$ox_fraction[cv$part == "ESRRA only"],
            100 * cv$ox_fraction[cv$part == "NRF1 only"], ct$fold, ct$p_enrich),
    sprintf("NON-CIRCULARITY: the PGC1a-axis regulon holds %d of the 25 pro-apoptotic MitoCarta genes against %.1f expected (%.2f-fold) -- Aifm1, Aifm2, Cycs, Endog, Ifi27, with no PUMA, no BAX and no caspase. MYC's regulon holds %d (%.2f-fold, p = %.3f) including Bax, Casp9, Diablo and Pmaip1.",
            ows$overlap[ows$programme == "CORE_MITO"], ows$expected[ows$programme == "CORE_MITO"],
            ows$fold[ows$programme == "CORE_MITO"], ows$overlap[ows$programme == "MYC_MITO"],
            ows$fold[ows$programme == "MYC_MITO"], ows$p_enrich[ows$programme == "MYC_MITO"]),
    "CALIBRATION, not drawn: the deficit is x1.19 in the wild type and x1.32 in the Myc+ gland at OXPHOS-subunit level. Since PGC-1a is not expressed in these cells, calibrate on the output, not on the transgene."),
  bounds = c(
    "State the fold honestly: 0.76 and 0.56 are below expectation but neither is significantly DEPLETED -- with 25 genes the test has no power for that. What carries is the membership fact, the specific absence of PUMA, BAX and the caspases, which needs no p-value.",
    "The one death gene the axis does reach is cytochrome c, which is simultaneously a respiratory carrier. It should be named rather than passed over.",
    "PRE-SPECIFY THE OVERSHOOT: the developmental change was arm-selective and the rescue will not be. Assembly factors, mitoribosome and central dogma never fell yet are 29 to 35 per cent inside the regulon, so they will rise. The supported claim is then 'raising respiratory capacity restores death competence', not 'reverting the developmental change restores it'.",
    "Regulons are in-silico memberships from curated resources, not binding measured in these mice."),
  source = "results/biogenesis_axis_developmental.rds (script 47 PART H)")

save_panel_p(p, "biogax_rescue_licence", height = 104)

if (FALSE) {
  print(p)
  ba$rescue$reach |> print()
  ba$rescue$convergence |> print()
  ba$rescue$convergence_test |> print()
  ba$rescue$ownership |> dplyr::filter(apoptosis_set == "MITOCARTA_APOPTOSIS_PRO") |>
    print()
  ba$rescue$calibration |> print()
  ba$rescue$overshoot |> print()
}
