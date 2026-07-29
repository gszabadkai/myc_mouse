# =============================================================================
# figS1A_geneset_library.R -- composition of the custom gene-set library
# -----------------------------------------------------------------------------
# SUPPORTS (Results, "Myc drives early breast tumourigenesis by inducing OXPHOS
# and biosynthetic pathways", paragraph 1):
#   "we used fGSEA ranking and GSVA scoring (Fig. S1A)"
#   "we constructed and quantified (by effect size, ssGSVA and fGSEA) a large
#    custom library of ~900 genesets, covering MYC identity (signatures),
#    oncogenic and cell fate related signalling, bioenergetic and biosynthetic
#    metabolic pathways (Fig. S1A)"
#
# The panel is the library's own provenance record, read straight off the v1.0
# snapshot -- CLAUDE.md: consume data/genesets_from_library/, never rebuild the
# sets. Each set carries a `method` tag saying which quantifier it is routed to;
# that tag is the fill, because "quantified by ssGSVA and fGSEA" is a claim
# about routing and the reader should be able to see the split.
#
# A NUMBER THE TEXT SHOULD BE MADE EXACT ON. There are four defensible counts
# and the manuscript currently says "~900":
#   986  sets in the library snapshot (this panel's total)
#   902  scored by GSVA          (method 'gsva' or 'both'; scripts/15:194)
#   817  fGSEA-eligible          (provenance column fgsea_eligible)
#   884  used in the per-sample coupling analyses (scripts 35-38)
# All four are in the legend block below. 986 is the right number for "the
# library"; 902 is the right number for "scored by GSVA".
#
# Input:  data/genesets_from_library/provenance_table.csv   (v1.0 snapshot)
# Output: outputs/figures/panels/figS1A_geneset_library.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

prov_path <- here::here("data", "genesets_from_library", "provenance_table.csv")
stopifnot(file.exists(prov_path))

# utils::read.csv keeps the figure layer free of a readr dependency; the two
# columns used here are ASCII (the latin1 risk in this project lives in the
# free-text description/notes columns, which are not read).
prov <- utils::read.csv(prov_path, stringsAsFactors = FALSE)
stopifnot(all(c("set_name", "category_primary", "method",
                "size_mouse", "fgsea_eligible") %in% names(prov)))

# --- counts that go in the legend, computed here so figure and text agree ----
n_total  <- nrow(prov)
n_gsva   <- sum(prov$method %in% c("gsva", "both"))
n_fgsea  <- sum(prov$fgsea_eligible)
size_med <- stats::median(prov$size_mouse)
size_rng <- range(prov$size_mouse)

# --- human-readable category labels (labelling, not explanation) -------------
cat_lab <- c(
  "TF_targets"                         = "TF regulons",
  "Mammary_development"                = "Mammary development",
  "MitoCarta"                          = "MitoCarta",
  "Metabolism"                         = "Metabolism",
  "Apoptosis"                          = "Apoptosis",
  "Biogenesis_apoptosis_intersections" = "Biogenesis x apoptosis",
  "MYC_signatures"                     = "MYC signatures",
  "Biogenesis_discrimination"          = "Biogenesis discrimination",
  "Proliferation"                      = "Proliferation")
stopifnot(setequal(names(cat_lab), unique(prov$category_primary)))

method_lab <- c(fgsea = "fGSEA", gsva = "GSVA", both = "both")
stopifnot(setequal(names(method_lab), unique(prov$method)))

# --- tabulate ----------------------------------------------------------------
tab <- prov |>
  dplyr::mutate(category = unname(cat_lab[category_primary]),
                quant    = factor(unname(method_lab[method]),
                                  levels = c("fGSEA", "both", "GSVA"))) |>
  dplyr::count(category, quant, name = "n")

totals <- tab |>
  dplyr::group_by(category) |>
  dplyr::summarise(n_cat = sum(n), .groups = "drop") |>
  dplyr::arrange(n_cat)
stopifnot(sum(totals$n_cat) == n_total)

cat_levels <- totals$category
tab    <- dplyr::mutate(tab,    category = factor(category, levels = cat_levels))
totals <- dplyr::mutate(totals, category = factor(category, levels = cat_levels))

quant_cols <- c("fGSEA" = unname(method_cols[["fgsea"]]),
                "both"  = unname(method_cols[["both"]]),
                "GSVA"  = unname(method_cols[["gsva"]]))

# --- panel -------------------------------------------------------------------
XMAX <- max(totals$n_cat) * 1.14

# position_stack(reverse = TRUE) so the segments read left-to-right in the same
# order as the key; ggplot's default stacks the reverse of the factor levels.
# The hairline white border keeps the palest segment delimited on white.
p <- ggplot2::ggplot(tab, ggplot2::aes(x = n, y = category, fill = quant)) +
  ggplot2::geom_col(width = 0.68, colour = "white", linewidth = 0.2,
                    position = ggplot2::position_stack(reverse = TRUE)) +
  ggplot2::geom_text(data = totals,
                     ggplot2::aes(x = n_cat, y = category, label = n_cat),
                     inherit.aes = FALSE, hjust = -0.28, size = 2.1) +
  ggplot2::scale_fill_manual(values = quant_cols, name = NULL) +
  ggplot2::scale_x_continuous(limits = c(0, XMAX), expand = c(0, 0)) +
  ggplot2::labs(x = "gene sets (n)", y = NULL) +
  theme_panel() +
  ggplot2::theme(legend.position = "bottom",
                 legend.margin = ggplot2::margin(-3, 0, 0, 0),
                 axis.line.y  = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank())

# --- the legend text (never drawn) -------------------------------------------
LEGEND <- panel_legend(
  slot = "Fig. S1A",
  what = paste0(
    "Composition of the custom gene-set library used throughout: ", n_total,
    " sets in nine primary categories, coloured by the quantifier each set is ",
    "routed to."),
  detail = c(
    sprintf("Library snapshot: mammary_geneset_library v1.0, %d sets.", n_total),
    sprintf("Quantifier routing is a per-set tag in the library's provenance table: fGSEA only (%d), GSVA only (%d), or either (%d).",
            sum(prov$method == "fgsea"), sum(prov$method == "gsva"),
            sum(prov$method == "both")),
    sprintf("%d sets are GSVA-scored in one cohort-relative run; %d are fGSEA-eligible.",
            n_gsva, n_fgsea),
    sprintf("Set size (mouse symbols): median %d genes, range %d-%d.",
            size_med, size_rng[1], size_rng[2]),
    "fGSEA reads absolute transcriptional enrichment for a given contrast, ranking on the unshrunken Wald statistic; GSVA gives a per-sample, cohort-relative score for each programme. The two answer different questions and are not interchangeable."),
  bounds = c(
    "Categories overlap by construction (a MitoCarta OXPHOS set and a TF regulon can share genes), so the counts are a description of the library, not of independent tests.",
    "TF regulons dominate the count (444 of 986) but are the least specific category; set number is not evidence weight.",
    paste0("The manuscript text currently says '~900 genesets'. The exact figure depends on which count is meant: ",
           n_total, " in the library, ", n_gsva, " GSVA-scored, ", n_fgsea,
           " fGSEA-eligible, 884 used in the per-sample coupling analyses.")),
  source = c(
    "data/genesets_from_library/provenance_table.csv (v1.0 snapshot; provenance in that directory's README.md)",
    "GSVA scoring: scripts/15_gsva_scoring.R; fGSEA per category: scripts/20_fgsea_percategory.R"))

save_panel_p(p, "figS1A_geneset_library",
             width = fig_w[["single"]], height = 58)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the counts behind the bars
  tab |> tidyr::pivot_wider(names_from = quant, values_from = n, values_fill = 0) |>
    print()

  ## the four candidate totals, side by side
  c(library = n_total, gsva_scored = n_gsva, fgsea_eligible = n_fgsea,
    coupling_analyses = 884L) |> print()
}
