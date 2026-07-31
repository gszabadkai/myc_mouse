# =============================================================================
# _panel_common.R -- shared spine for the panel-by-panel figure layer
# -----------------------------------------------------------------------------
# figures/panels/ holds ONE SCRIPT PER PANEL, each built to support a specific
# sentence of the Results section as it is written. Assembly into numbered
# manuscript figures happens last, once the narrative fixes the numbering --
# see figures/panels/PANELS.md for the slug -> slot map.
#
# THE PUBLICATION RULE THIS FILE ENFORCES
#
# Nature-family figures carry no explanatory information: the panel shows data,
# axes, units and a minimal key, and everything else belongs in the figure
# legend. theme_panel() therefore BLANKS title, subtitle and caption, so a panel
# cannot accidentally acquire explanatory text. The explanation is written into
# each script as a panel_legend() block, which rebuild_panels.R collects into
# outputs/figures/panels/legends.md -- that file is the raw material for the
# legends, not the figure.
#
# RELATIONSHIP TO figures/theme_myc.R
#
# This file SOURCES theme_myc.R and never edits it, so the four assembled
# manuscript figures (figure1, figure2, figureS1, figureS2) keep rendering
# byte-identically. What it adds is (a) the panel theme, (b) the palettes that
# are currently copy-pasted across a dozen scripts, and (c) a save wrapper that
# honours myc.fig.nosave and writes to outputs/figures/panels/.
#
# CONVENTION every panel script follows
#
#   source(here::here("figures", "panels", "_panel_common.R"))
#   LEGEND <- panel_legend(slot = "Fig. 1B", what = ..., detail = ..., bounds = ...)
#   ... build `p` ...
#   save_panel_p(p, "fig1B_myc_teb_proliferation", height = 80)
#   if (FALSE) { print(p) }        # sandbox, skipped by source()
#
# The composite object is ALWAYS named `p`: both rebuild_panels.R and
# figures/rebuild_manuscript_figures.R locate it by that name to render a dry
# run through a null device.
# =============================================================================

source(here::here("figures", "theme_myc.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("panels need ggplot2")

panel_dir <- here::here("outputs", "figures", "panels")

# --- the panel theme ---------------------------------------------------------
# base_size 7: these are single- or half-column panels that will later be
# composed by patchwork, so type is set for the FINAL size, not the preview.
theme_panel <- function(base_size = 7) {
  theme_myc(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_blank(),   # the publication rule,
      plot.subtitle = ggplot2::element_blank(),   # made mechanical
      plot.caption  = ggplot2::element_blank(),
      legend.title  = ggplot2::element_text(size = base_size),
      legend.text   = ggplot2::element_text(size = base_size),
      legend.key.size = ggplot2::unit(3.2, "mm"),
      plot.margin   = ggplot2::margin(2, 2, 2, 2, "mm"))
}

# --- palettes, declared once -------------------------------------------------
# Every one of these is currently re-declared verbatim in several scripts under
# figures/ (contrast: fig01b:100, figS1b:109, figS2b:73; verdict: fig03:202,
# fig04:67, figure2:76; tier: fig02:40, fig03:43, figure1:60 and five more).
# New panels take them from here. The existing scripts are deliberately left
# untouched.

# --- the contrast vocabulary (author's naming, 2026-07-30) -------------------
# Two families, and the figures say which is which by NAME, not by a note on the
# page: genotype contrasts are the Myc effect measured at one age, development
# contrasts are the 6->12W trajectory within one genotype. These strings are what
# gets drawn; the mapping to the DESeq2 slot names lives in the legend blocks.
contrast_geno <- c("myc_6W", "myc_12W")
contrast_dev  <- c("6>12W_wt", "6>12W_myc")
contrast_levels <- c(contrast_geno, contrast_dev)

# Colours: the genotype contrasts take the Myc+ hues, because that is what they
# measure; the development contrasts take greys, because the trajectory is the
# background against which the Myc effect is read. Deliberately NOT the WT blue
# for 6>12W_wt -- that would put a sample colour on a contrast.
contrast_cols <- c("myc_6W"    = "#D55E00", "myc_12W"   = "#E69F00",
                   "6>12W_wt"  = "#BDBDBD", "6>12W_myc" = "#7B7B7B")

verdict_cols  <- c("withdraws" = "#762A83", "at chance" = "grey55",
                   "rises"     = "#1B7837")

direction_cols <- c(up = "#D6604D", down = "#4393C3")

# MitoPathway Level-1 tiers: Okabe-Ito, with grey72 for the Metabolism catch-all.
tier_cols <- c("OXPHOS"                  = "#E69F00",
               "Mitochondrial dynamics"  = "#56B4E9",
               "Metabolism"              = "grey72",
               "Protein import, sorting" = "#009E73",
               "Signaling"               = "#F0E442",
               "Small molecule transport" = "#0072B2",
               "Central dogma"           = "#D55E00")

# Gene-set quantification methods, as tagged in the library provenance table.
# A single-hue purple ramp, deliberately NOT the genotype palette: across the
# figure set a reader learns blue = WT and orange-red = Myc+, and reusing those
# for a methods key would spend them on something that is not a genotype.
# NOT used in the paper figures -- the author's call (2026-07-30) is that the
# fGSEA/GSVA routing is minor and belongs in the internal write-up. Kept here for
# the paper/myc_mito.qmd version of the library panel.
method_cols <- c("fgsea" = "#54278F", "both" = "#9E9AC8", "gsva" = "#DADAEB")

# Mitochondrial definition of a gene set, three classes, as script 37 defines
# them (37:445-462) and as the reference figure outputs/pathway_loading/
# E_library_coverage.pdf draws them. The middle class is the one that matters:
# the Gray _MITO / _LE_MITO TF lanes ARE MitoCarta subsets by construction, so
# they are mitochondrial by build rather than by biology and must not be counted
# as independent mitochondrial coverage.
mito_class_cols <- c("mitocarta_proper"  = "#D73027",
                     "construction_MITO" = "#FC8D59",
                     "non_mito"          = "#4575B4")
mito_class_labels <- c("mitocarta_proper"  = "MitoCarta / OXPHOS",
                       "construction_MITO" = "mitochondrial by construction",
                       "non_mito"          = "non-mitochondrial")

# The rule itself, so the three panels that need it (the library composition, the
# axis loadings and the enrichment ranking) cannot each grow their own copy.
# Verbatim script 37 (37:445-462); `category` is the library's category_primary.
# Sets outside the library -- the fresh MSigDB Hallmark comparators in the fGSEA
# ranking -- fall through to non_mito, which is a naming rule, not a judgement:
# say so wherever they are drawn.
mito_class3 <- function(set, category) {
  category  <- ifelse(is.na(category), "", category)
  name_mito <- grepl("_MITO$|_MITO_|MITO_NU|^MITO_|CORE_MITO", set)
  mito_defined <- category == "MitoCarta" | name_mito |
    (category == "Metabolism" & grepl("OXPHOS|KREBS|TCA|ELECTRON|RESPIRAT", set))
  ifelse(category == "MitoCarta", "mitocarta_proper",
    ifelse(name_mito, "construction_MITO",
      ifelse(mito_defined, "mitocarta_proper", "non_mito")))
}

# --- THE manuscript diverging fill -------------------------------------------
# Author's specification (2026-07-30), to be used for every diverging quantity in
# the manuscript: deep espresso brown at the negative extreme, stark white at
# zero, crisp mint green at the positive extreme. Deliberately not blue-red --
# blue and orange-red carry genotype meaning everywhere else in the figure set,
# so a blue-red fill would invite the reader to see genotype in it.
ms_diverging <- c(neg = "#4A3525", zero = "#FAFAFA", pos = "#2A8A6D")

# Signed tick labels, with the author's rule that zero is written "0" and never
# "0.0" or "+0.0". %+g drops trailing zeros, so 3 -> "+3" and 1.5 -> "+1.5".
lab_signed <- function(x) ifelse(x == 0, "0", sprintf("%+g", x))

# WHITE IS PINNED TO ZERO, AND THE TWO SIDES ARE SCALED SEPARATELY. `limits` may
# be one number (symmetric, +/- that) or two (the observed range). The asymmetric
# form matters whenever the data are lopsided: on a +3.0 / -1.7 range a symmetric
# ramp leaves every negative value in the first third of the brown, where it is
# indistinguishable from zero, so the fill stops carrying magnitude on that side.
# Pinning white to zero and rescaling each arm keeps the sign unambiguous and
# spends the whole ramp -- at the cost that equal ink no longer means equal
# magnitude ACROSS the sign change. That trade is only acceptable where a
# quantitative axis carries the magnitude anyway; say so in the legend block.
heat_fill <- function(limits, name = NULL, breaks = ggplot2::waiver()) {
  if (length(limits) == 1L) limits <- c(-abs(limits), abs(limits))
  stopifnot(length(limits) == 2L, limits[1] < 0, limits[2] > 0)
  ggplot2::scale_fill_gradientn(
    colours = unname(ms_diverging[c("neg", "zero", "pos")]),
    values  = scales::rescale(c(limits[1], 0, limits[2])),
    limits  = limits, oob = scales::squish, name = name,
    breaks  = breaks, labels = lab_signed)
}

# Ink that stays legible on that ramp. The brown end goes dark fast and the mint
# end stays mid-tone, so the switch is asymmetric: white ink earlier on negative
# fills than on positive ones.
ink_on_fill <- function(x, limit) {
  ifelse(x < -0.45 * limit | x > 0.75 * limit, "white", "grey10")
}

# --- per-sample composites and the design contrasts --------------------------
# Both live here because fig1B and fig1C plot the same composites through
# different lenses and MUST NOT be allowed to drift apart.

# Pooled within-group SD: the project's standardisation convention (scripts/26:335).
# Not cosmetic -- a difference-of-two-means axis such as TEB minus ductal has ~1.7x
# the raw spread of a single composite for arithmetic reasons alone, so a shared
# raw axis would manufacture contrast.
wsd_of <- function(x, g) sqrt(mean(tapply(x, g, stats::var)))

# Mean GSVA score over a set of rows = the composite. One line, but naming it
# keeps every panel using the same definition.
composite_of <- function(scores, sets) {
  sets <- intersect(sets, rownames(scores))
  stopifnot(length(sets) > 0L)
  colMeans(scores[sets, , drop = FALSE])
}

# The four drawn contrasts plus the interaction, standardised by the programme's
# own within-group SD. Fitting idiom is script 27's (27:103-106) and the figure
# layer's own (fig01_mito_content.R:82-84): each contrast is an ordinary least
# squares fit on the relevant subset, so the two genotype gaps and the two
# trajectories are estimated the same way rather than read off one pooled model.
contrast_table <- function(y, timepoint, myc_status, group) {
  stopifnot(length(y) == length(timepoint), length(y) == length(myc_status),
            length(y) == length(group))
  d <- data.frame(y = as.numeric(y),
                  tp = factor(as.character(timepoint), levels = c("6W", "12W")),
                  myc = factor(as.character(myc_status), levels = c("neg", "pos")))
  cf <- function(form, sub, term) {
    stats::coef(summary(stats::lm(form, data = d[sub, , drop = FALSE])))[term, ]
  }
  g6  <- cf(y ~ myc, d$tp  == "6W",  "mycpos")
  g12 <- cf(y ~ myc, d$tp  == "12W", "mycpos")
  twt <- cf(y ~ tp,  d$myc == "neg", "tp12W")
  tmp <- cf(y ~ tp,  d$myc == "pos", "tp12W")
  int <- cf(y ~ tp * myc, rep(TRUE, nrow(d)), "tp12W:mycpos")
  w   <- wsd_of(d$y, group)
  data.frame(
    contrast = c(contrast_levels, "interaction"),
    effect = c(g6[1], g12[1], twt[1], tmp[1], int[1]) / w,
    p      = c(g6[4], g12[4], twt[4], tmp[4], int[4]),
    within_sd = w, row.names = NULL, stringsAsFactors = FALSE)
}

# --- the dominant pathway axis -----------------------------------------------
# Figs. 1C and 1D are the sample scores and the per-set loadings of ONE principal
# component analysis, so they must not each compute it. results/pathway_loading.rds
# saves the per-set loadings and the variance percentages but NOT the score matrix
# or the sample scores, so this rebuilds them -- verbatim script 37 PART 2 and
# PART A2 (37:103-122, 37:221-222) -- and then proves the rebuild against the
# analysis of record before returning anything.
#
# The quantifier is the LINEAR mean gene-wise z-score, not GSVA: script 36 made it
# the method of record because its correlation IS average cross-gene covariance.
# GSVA gives the same picture more weakly (PC1 68% against 77%).
#
# Deterministic -- no RNG, no permutation -- so a re-run cannot drift; what could
# drift is the input, and that is what the assertions catch.
pathway_axis <- function(gsva_path = here::here("results", "gsva_scores.rds"),
                         ref_path  = here::here("results", "pathway_loading.rds")) {
  require_fresher_than(ref_path)
  gs <- readRDS(gsva_path)
  pl <- readRDS(ref_path)

  M_gsva <- gs$scores
  expr   <- gs$expr_mat[, colnames(M_gsva), drop = FALSE]      # VST, genes x samples
  sets   <- rownames(M_gsva)
  pw     <- gs$pathways[intersect(names(gs$pathways), sets)]
  stopifnot(identical(as.character(gs$sample_meta$sample), colnames(M_gsva)))

  zg <- t(scale(t(expr)))
  zg <- zg[is.finite(rowSums(zg)), , drop = FALSE]
  M  <- t(vapply(sets, function(s) {
    g <- intersect(pw[[s]], rownames(zg))
    if (length(g) < 5L) rep(NA_real_, ncol(zg)) else colMeans(zg[g, , drop = FALSE])
  }, numeric(ncol(zg))))
  # a set counts only if it is scored in BOTH matrices, so the universe matches
  keep <- intersect(rownames(M)[stats::complete.cases(M)],
                    rownames(M_gsva)[stats::complete.cases(M_gsva)])
  M <- M[keep, , drop = FALSE]

  pc  <- stats::prcomp(t(M), center = TRUE, scale. = FALSE)
  vfr <- pc$sdev^2 / sum(pc$sdev^2)
  gm  <- colMeans(M)                                   # per-sample global mean
  sc  <- pc$x
  # orient PC1 to track the global mean, so "high" means "high everywhere"
  if (suppressWarnings(stats::cor(sc[, 1], gm)) < 0) sc[, 1] <- -sc[, 1]
  lam <- stats::setNames(as.numeric(suppressWarnings(stats::cor(t(M), sc[, 1]))),
                         rownames(M))

  # --- prove the rebuild against script 37 ------------------------------------
  ref_v <- pl$sample_pca_var[pl$sample_pca_var$space == "pathway_scores_884", ]
  mito  <- intersect(gs$set_meta$set_name[gs$set_meta$category_primary == "MitoCarta"], keep)
  ox    <- grep("OXPHOS|COMPLEX_[IV]|_SUBUNITS|ASSEMBLY_FACTORS|ELECTRON_CARRIERS|CRISTAE",
                mito, value = TRUE)
  ox_gate <- suppressWarnings(stats::cor(colMeans(M[ox, , drop = FALSE]), gm))
  lm_ref  <- pl$loading_movement
  conc    <- suppressWarnings(stats::cor(lam[lm_ref$set], lm_ref$lambda_design))
  stopifnot(
    length(keep) == pl$n_sets,                                    # 885 sets
    nrow(ref_v) == 1L,
    max(abs(100 * vfr[1:3] -
            c(ref_v$pc1_pct, ref_v$pc2_pct, ref_v$pc3_pct))) < 1e-6,
    abs(ox_gate - pl$ox_gate) < 1e-6,
    # the raw axis and script 37's design-RESIDUAL axis are different objects;
    # this is a concordance floor, not an identity (observed 0.968)
    conc >= 0.9)

  list(M = M, scores = sc, var_frac = vfr, loading = lam, global_mean = gm,
       ox_gate = ox_gate, concordance_with_residual_axis = conc,
       sample_meta = gs$sample_meta, set_meta = gs$set_meta, ref = pl)
}

# --- programme grouping for the enrichment ranking ---------------------------
# Fig. S1C's y axis. Declared here rather than inline because it IS an encoding:
# it decides what the reader sees as "a programme", and the sentence it supports
# names the tiers. Two rules keep it honest.
#
# (1) OXPHOS and mitochondrial biogenesis are split by script 37's OWN patterns
#     (37:137-138), so the panel and the loading analysis mean the same thing by
#     the words.
# (2) Every set whose membership is MitoCarta intersected with something else --
#     the Gray _MITO TF lanes, and the biogenesis-discrimination and
#     biogenesis x apoptosis constructs -- goes into ONE row, because they are
#     mitochondrial by build and cannot be evidence that a mitochondrial
#     programme was independently detected. Hiding them would flatter the result;
#     spreading them across the other rows would inflate every one of them.
#
# `category` is the fGSEA table's category column (the GMT file stem, e.g.
# "01_mitocarta", plus "hallmark_msigdb"), which is what the consumer has.
# Row names are the author's (2026-07-31). Two of them lean on the colour key
# rather than saying it twice: "curated mitochondrial" is the by-construction
# block, and the key that colours it says "mitochondrial by construction".
programme_levels <- c(
  "curated mitochondrial", "mitochondrial biogenesis", "OXPHOS", "TCA cycle",
  "mitochondrial metabolism & dynamics", "Myc signatures", "E2F / cell cycle",
  "biosynthetic metabolism", "metabolism, other", "apoptosis",
  "TF target sets", "mammary development", "MSigDB Hallmarks")

programme_group <- function(set, category) {
  OX   <- "OXPHOS|COMPLEX_[IV]|_SUBUNITS|ASSEMBLY_FACTORS|ELECTRON_CARRIERS|CRISTAE"
  BIOG <- "RIBOSOME|CENTRAL_DOGMA|MT_TRNA|MT_RRNA|MTRNA|MTDNA|IMPORT|TRANSLATION"
  BIOSYN <- paste0("NUCLEOTIDE|PURINE|PYRIMIDINE|AMINO|SER_GLY|BCAA|ONE_CARBON|",
                   "PPP|PENTOSE|POLYAMINE|CHOLESTEROL|MEVALONATE|LIPID|FATTY|GLYCOLYSIS")
  constructed <- grepl("_MITO$|_MITO_|MITO_NU|^MITO_|CORE_MITO", set) |
    category %in% c("07_biogenesis_discrimination",
                    "09_biogenesis_apoptosis_intersections")
  g <- ifelse(
    constructed, "curated mitochondrial",
    ifelse(category == "01_mitocarta",
           ifelse(grepl(BIOG, set), "mitochondrial biogenesis",
                  ifelse(grepl(OX, set), "OXPHOS",
                         "mitochondrial metabolism & dynamics")),
    ifelse(category == "04_metabolism",
           ifelse(grepl("OXPHOS|ELECTRON|RESPIRAT", set), "OXPHOS",
                  ifelse(grepl("KREBS|_TCA", set), "TCA cycle",
                         ifelse(grepl(BIOSYN, set), "biosynthetic metabolism",
                                "metabolism, other"))),
    ifelse(category == "02_myc_signatures",      "Myc signatures",
    ifelse(category == "05_proliferation",       "E2F / cell cycle",
    ifelse(category == "06_tf_targets",          "TF target sets",
    ifelse(category == "03_mammary_development", "mammary development",
    ifelse(category == "08_apoptosis",           "apoptosis",
    ifelse(category == "hallmark_msigdb",        "MSigDB Hallmarks",
           NA_character_)))))))))
  factor(g, levels = programme_levels)
}

# --- export ------------------------------------------------------------------
# Wraps theme_myc.R's save_panel() so a panel script never needs to know the
# output directory, and so myc.fig.nosave is honoured in ONE place rather than
# being re-implemented as an `if` in every script.
save_panel_p <- function(plot, slug,
                         width = fig_w[["single"]], height = 70, units = "mm") {
  stopifnot(is.character(slug), length(slug) == 1L, nzchar(slug))
  if (isTRUE(getOption("myc.fig.nosave"))) {
    message("myc.fig.nosave = TRUE -- not writing ", slug, ".pdf")
    return(invisible(NULL))
  }
  if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)
  save_panel(plot, file.path(panel_dir, paste0(slug, ".pdf")),
             width = width, height = height, units = units)
}

# --- the legend block --------------------------------------------------------
# Everything the panel is NOT allowed to say on the page. Four fields, because
# four things have to survive the trip to the legend:
#   what   -- one sentence naming what is plotted (the legend's first line)
#   detail -- n, units, the test, thresholds, what any error bar is
#   bounds -- the caveat that must travel with the claim (batch = timepoint,
#             n=6 power floor, a circularity, whatever applies). Not optional in
#             this project: most of its numbers carry one.
#   source -- the results/*.rds objects and script the numbers came from
panel_legend <- function(slot, what, detail, bounds, source = NULL) {
  stopifnot(is.character(slot),   length(slot)   == 1L,
            is.character(what),   length(what)   == 1L,
            is.character(detail), length(detail) >= 1L,
            is.character(bounds), length(bounds) >= 1L)
  structure(list(slot = slot, what = what, detail = detail,
                 bounds = bounds, source = source),
            class = "panel_legend")
}

# Markdown rendering, used by rebuild_panels.R to build legends.md.
legend_md <- function(x, slug = NULL) {
  stopifnot(inherits(x, "panel_legend"))
  bul <- function(v) paste0("- ", v, collapse = "\n")
  c(paste0("### ", x$slot, if (!is.null(slug)) paste0("  (`", slug, "`)") else ""),
    "", x$what, "",
    "**Detail.**", bul(x$detail), "",
    "**Bounds.**", bul(x$bounds), "",
    if (!is.null(x$source)) c("**Source.**", bul(x$source), "") else NULL)
}

# Panel filenames stopped carrying the slot letter on 2026-07-31, so name order is
# no longer figure order; legends.md is ordered by the slot each block declares
# instead. Main figures before supplementary, then number, then letter -- the
# order a reader meets them in.
slug_slot_order <- function(slugs, legends) {
  slot <- vapply(slugs, function(s) legends[[s]]$slot, character(1))
  key  <- gsub("[^A-Z0-9]", "", toupper(slot))          # "Fig. S1A" -> "FIGS1A"
  supp <- grepl("^FIGS", key)
  rest <- sub("^FIGS?", "", key)
  num  <- suppressWarnings(as.integer(sub("^([0-9]+).*$", "\\1", rest)))
  let  <- sub("^[0-9]+", "", rest)
  slugs[order(supp, num, let, slugs)]
}

print.panel_legend <- function(x, ...) {
  cat("\n", strrep("-", 74), "\n", sep = "")
  cat(paste(legend_md(x), collapse = "\n"), "\n", sep = "")
  cat(strrep("-", 74), "\n\n", sep = "")
  invisible(x)
}

# --- freshness guard ---------------------------------------------------------
# results/gsva_scores.rds was rebuilt 2026-07-24 with the gene-symbol reconciler
# (docs/2026-07-24_symbol_reconciliation.md). Anything downstream of it that
# predates that rebuild is stale, and a stale panel is worse than a missing one
# because it looks finished. Panels that read a derived GSVA object call this.
require_fresher_than <- function(path, reference = here::here("results", "gsva_scores.rds")) {
  if (!file.exists(path)) {
    stop(basename(path), " is missing -- re-source its producing script (see PANELS.md)",
         call. = FALSE)
  }
  if (file.exists(reference) && file.mtime(path) < file.mtime(reference)) {
    stop(basename(path), " (", format(file.mtime(path), "%Y-%m-%d"), ") is OLDER than ",
         basename(reference), " (", format(file.mtime(reference), "%Y-%m-%d"), ").\n",
         "  It was built before the gene-symbol reconciliation. Re-source, in order:\n",
         "    scripts/17_gsva_overview.R -> scripts/26_dev_program_myc_integration.R\n",
         "    scripts/27_myc_endogenous_amplification.R",
         call. = FALSE)
  }
  invisible(TRUE)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  ## what the legend block looks like
  panel_legend(slot = "Fig. 9Z", what = "A demonstration panel.",
               detail = c("n = 24.", "Bars are SE."),
               bounds = "Nothing is claimed here.") |> print()

  ## dry run of every panel
  options(myc.fig.nosave = TRUE)
  source(here::here("figures", "panels", "rebuild_panels.R"))
}
