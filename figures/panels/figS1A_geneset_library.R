# =============================================================================
# figS1A_geneset_library.R -- composition of the custom gene-set library
# -----------------------------------------------------------------------------
# SUPPORTS (Results, "Myc drives early breast tumourigenesis by inducing OXPHOS
# and biosynthetic pathways", paragraph 1):
#   "we constructed and quantified (by effect size, ssGSVA and fGSEA) a large
#    custom library of ~900 genesets, covering MYC identity (signatures),
#    oncogenic and cell fate related signalling, bioenergetic and biosynthetic
#    metabolic pathways (Fig. S1A)"
#
# The panel is the library's own provenance record, read straight off the v1.0
# snapshot -- CLAUDE.md: consume data/genesets_from_library/, never rebuild the
# sets.
#
# WHAT THE FILL IS, AND WHY IT CHANGED (author, 2026-07-30). It used to be the
# quantifier tag (fGSEA / GSVA / both). That is minor routing information and
# belongs in the internal write-up, not in a manuscript supplement. The fill is
# now the MITOCHONDRIAL definition of each set, in the three classes script 37
# uses and the reference figure outputs/pathway_loading/E_library_coverage.pdf
# draws. The middle class carries the argument: the Gray _MITO / _LE_MITO TF
# lanes are MitoCarta SUBSETS BY CONSTRUCTION, so they are mitochondrial by
# build, not by biology, and counting them as mitochondrial coverage would
# inflate it. 260 of the 393 mito-labelled sets are of that kind.
#
# THE CLASSIFICATION IS SCRIPT 37'S, NOT A SECOND OPINION. The rule (37:445-462)
# is a function of set name and category only, so it is re-derived here over the
# WHOLE library rather than over the 885 sets that reached the loading analysis --
# and then asserted identical to script 37's saved mito_classification on the 885
# they share. If the two ever diverge the panel stops.
#
# A NUMBER THE TEXT SHOULD BE MADE EXACT ON. There are four defensible counts and
# the manuscript currently says "~900":
#   986  sets in the library snapshot (this panel's total)
#   885  actually scored by GSVA      (gsva_scores.rds$n_sets_scored; 902 were
#                                      routed to GSVA, 17 fell to the size filter)
#   817  fGSEA-eligible               (provenance column fgsea_eligible)
#   885  used in the per-sample coupling and loading analyses (scripts 35-38)
# All four are in the legend block. 986 is the right number for "the library";
# 885 is the right number for "scored".
#
# Input:  data/genesets_from_library/provenance_table.csv   (v1.0 snapshot)
#         results/pathway_loading.rds                       (script 37 -- assertion only)
# Output: outputs/figures/panels/figS1A_geneset_library.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

prov_path <- here::here("data", "genesets_from_library", "provenance_table.csv")
stopifnot(file.exists(prov_path))

# utils::read.csv keeps the figure layer free of a readr dependency; the columns
# used here are ASCII (the latin1 risk in this project lives in the free-text
# description/notes columns, which are not read).
prov <- utils::read.csv(prov_path, stringsAsFactors = FALSE)
stopifnot(all(c("set_name", "category_primary", "method",
                "size_mouse", "fgsea_eligible") %in% names(prov)))

# --- the mitochondrial classification, script 37's rule ----------------------
lib <- prov |>
  dplyr::transmute(
    set       = set_name,
    category  = dplyr::coalesce(category_primary, ""),
    name_mito = grepl("_MITO$|_MITO_|MITO_NU|^MITO_|CORE_MITO", set_name)) |>
  dplyr::mutate(
    mito_defined = category == "MitoCarta" | name_mito |
                   (category == "Metabolism" &
                    grepl("OXPHOS|KREBS|TCA|ELECTRON|RESPIRAT", set)),
    class3 = dplyr::case_when(
      category == "MitoCarta" ~ "mitocarta_proper",
      name_mito               ~ "construction_MITO",
      mito_defined            ~ "mitocarta_proper",   # Metabolism OXPHOS/TCA sets
      TRUE                    ~ "non_mito"))

# assertion: identical to script 37 on every set the two share
pl_path <- here::here("results", "pathway_loading.rds")
stopifnot(file.exists(pl_path))
ref <- readRDS(pl_path)$mito_classification
chk <- dplyr::inner_join(dplyr::select(lib, set, mine = class3),
                         dplyr::select(ref, set, theirs = class3), by = "set")
stopifnot(nrow(chk) == nrow(ref), identical(chk$mine, chk$theirs))

# --- counts that go in the legend, computed here so figure and text agree ----
n_total  <- nrow(prov)
n_gsva   <- readRDS(here::here("results", "gsva_scores.rds"))$n_sets_scored
n_fgsea  <- sum(prov$fgsea_eligible)
size_med <- stats::median(prov$size_mouse)
size_rng <- range(prov$size_mouse)
n_class  <- table(lib$class3)

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
stopifnot(setequal(names(cat_lab), unique(lib$category)))

# --- tabulate ----------------------------------------------------------------
class_levels <- c("mitocarta_proper", "construction_MITO", "non_mito")

tab <- lib |>
  dplyr::mutate(category = unname(cat_lab[category]),
                class3   = factor(class3, levels = class_levels)) |>
  dplyr::count(category, class3, name = "n")

totals <- tab |>
  dplyr::group_by(category) |>
  dplyr::summarise(n_cat = sum(n), .groups = "drop") |>
  dplyr::arrange(n_cat)
stopifnot(sum(totals$n_cat) == n_total)

cat_levels <- totals$category
tab    <- dplyr::mutate(tab,    category = factor(category, levels = cat_levels))
totals <- dplyr::mutate(totals, category = factor(category, levels = cat_levels))

# --- panel -------------------------------------------------------------------
XMAX <- max(totals$n_cat) * 1.14

# position_stack(reverse = TRUE) so the segments read left-to-right in the same
# order as the key; ggplot's default stacks the reverse of the factor levels.
# The hairline white border keeps adjacent segments delimited.
p <- ggplot2::ggplot(tab, ggplot2::aes(x = n, y = category, fill = class3)) +
  ggplot2::geom_col(width = 0.68, colour = "white", linewidth = 0.2,
                    position = ggplot2::position_stack(reverse = TRUE)) +
  ggplot2::geom_text(data = totals,
                     ggplot2::aes(x = n_cat, y = category, label = n_cat),
                     inherit.aes = FALSE, hjust = -0.28, size = 2.1) +
  ggplot2::scale_fill_manual(values = mito_class_cols,
                             labels = mito_class_labels,
                             breaks = class_levels, name = NULL) +
  ggplot2::scale_x_continuous(limits = c(0, XMAX), expand = c(0, 0)) +
  ggplot2::labs(x = "gene sets (n)", y = NULL) +
  # one class per row: "mitochondrial by construction" alone is wider than half of
  # an 89 mm column at this type size, and the key must stay in stacking order
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)) +
  theme_panel() +
  ggplot2::theme(legend.position = "bottom",
                 legend.margin = ggplot2::margin(-3, 0, 0, 0),
                 legend.key.size = ggplot2::unit(2.8, "mm"),
                 axis.line.y  = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank())

# --- the legend text (never drawn) -------------------------------------------
LEGEND <- panel_legend(
  slot = "Fig. S1A",
  what = paste0(
    "Composition of the custom gene-set library used throughout: ", n_total,
    " sets in nine primary categories, coloured by whether the set is defined ",
    "on mitochondrial genes."),
  detail = c(
    sprintf("Library snapshot: mammary_geneset_library v1.0, %d sets. Set size (mouse symbols): median %d genes, range %d-%d.",
            n_total, size_med, size_rng[1], size_rng[2]),
    sprintf("Mitochondrial classification, three classes: MitoCarta or Metabolism respiratory sets (%d), sets that are MitoCarta subsets by construction (%d), and non-mitochondrial sets (%d).",
            n_class[["mitocarta_proper"]], n_class[["construction_MITO"]],
            n_class[["non_mito"]]),
    sprintf("'By construction' means the Gray transcription-factor lanes named _MITO or _LE_MITO, whose members are that factor's target programme intersected with MitoCarta. They are mitochondrial by build, so they cannot be evidence that a mitochondrial programme was independently detected; %d of the %d mito-labelled sets are of this kind, and all but %d of them sit in the TF-regulon category.",
            n_class[["construction_MITO"]],
            n_class[["construction_MITO"]] + n_class[["mitocarta_proper"]],
            n_class[["construction_MITO"]] -
              sum(lib$class3 == "construction_MITO" & lib$category == "TF_targets")),
    sprintf("%d of the %d sets were scored per sample by GSVA in one cohort-relative run (%d were routed to GSVA; the remainder fell to the set-size filter); %d are fGSEA-eligible and are ranked on the unshrunken Wald statistic per contrast.",
            n_gsva, n_total, sum(prov$method %in% c("gsva", "both")), n_fgsea),
    "fGSEA reads absolute transcriptional enrichment for a given contrast; GSVA gives a per-sample, cohort-relative score for each programme. The two answer different questions and are not interchangeable."),
  bounds = c(
    "Categories overlap by construction (a MitoCarta OXPHOS set and a TF regulon can share genes), so the counts describe the library, not a set of independent tests.",
    sprintf("TF regulons dominate the count (%d of %d) and are the least specific category; set number is not evidence weight.",
            sum(lib$category == "TF_targets"), n_total),
    "The mitochondrial share of the library is 40 per cent by label but 13 per cent once the construction lanes are set aside. Any statement of the form 'the mitochondrial programmes dominate' has to say which of the two it means.",
    paste0("The manuscript text currently says '~900 genesets'. The exact figure depends on which count is meant: ",
           n_total, " in the library, ", n_gsva, " scored per sample, ", n_fgsea,
           " fGSEA-eligible.")),
  source = c(
    "data/genesets_from_library/provenance_table.csv (v1.0 snapshot; provenance in that directory's README.md)",
    "Mitochondrial classification: scripts/37_pathway_loading_and_technical_resolution.R:445-462, asserted identical here; reference figure outputs/pathway_loading/E_library_coverage.pdf",
    "GSVA scoring: scripts/15_gsva_scoring.R; fGSEA per category: scripts/20_fgsea_percategory.R"))

save_panel_p(p, "figS1A_geneset_library",
             width = fig_w[["single"]], height = 64)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the counts behind the bars
  tab |> tidyr::pivot_wider(names_from = class3, values_from = n, values_fill = 0) |>
    print()

  ## mito share, with and without the construction lanes
  c(labelled_mito = unname(n_class[["mitocarta_proper"]] + n_class[["construction_MITO"]]),
    genuine_mito  = unname(n_class[["mitocarta_proper"]]),
    total         = n_total) |> print()

  ## the four candidate totals, side by side
  c(library = n_total, scored = n_gsva, fgsea_eligible = n_fgsea) |> print()
}
