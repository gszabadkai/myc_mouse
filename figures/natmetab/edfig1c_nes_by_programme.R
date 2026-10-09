# =============================================================================
# edfig1c_nes_by_programme.R -- Extended Data Fig. 1c: every gene set's MYC
# enrichment at 6 and 12 weeks, by programme
# -----------------------------------------------------------------------------
# The text: "... but it retained its overall composition without selective
# repression of specific target gene sets." The amplitude halves (ED Fig. 1b); the
# ranking of programmes does not move. Read each row across the two columns.
#
# THE RULER: fGSEA normalized enrichment score on the unshrunken Wald statistic of
#   each genotype contrast -- enrichment relative to the whole transcriptome, so a
#   uniform halving of the effect does not lower it. That is why it shows
#   composition rather than amplitude.
# THE SETS: the library's fGSEA-eligible sets plus 50 MSigDB Hallmark comparators.
# THE ROWS: the programme grouping of the August layer (old Fig. 1D), ported
#   verbatim: OXPHOS and biogenesis split by script 37's patterns, and every set
#   whose membership is MitoCarta crossed with something else in one row,
#   "Curated mitochondrial", because those sets are mitochondrial by construction.
#
# FORM: many values in a few groups, two conditions -> one point per set
#   (data-to-viz, "Too many distributions": one row per programme, not per set;
#   "Order your data": rows by the 6-week median). Colour is the age, in ED 1a's
#   colours (MYC+ vs WT at 6W and at 12W), so the two panels read the same way;
#   the facet titles name the ages, so no key is drawn. The POINTS carry the
#   result -- the shape of each row at the two ages -- so they are drawn large
#   enough to read.
# NO CENTRAL VALUE IS DRAWN BY DEFAULT (SHOW_MEDIAN). The lower five rows are
#   bimodal, with sets enriched in both directions, so a median or mean sits inside
#   one mode or between them and moves with the balance of membership, not with the
#   shape: Apoptosis's median goes +0.74 -> -0.94 and its mean +0.17 -> -0.35, and
#   the means of three more rows change sign. Checked 2026-10-10.
#
# Reads:  results/fgsea_percategory.rds    (script 20)
#         results/ap6_permutation_null.rds (script 21; legend numbers only)
# Output: outputs/natmetab/EDFig1/EDFig1c_nes_by_programme.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

fg <- as.data.frame(readRDS(here::here("results", "fgsea_percategory.rds"))$fgsea)
stopifnot(all(c("ranking", "category", "pathway", "NES", "padj_within_category") %in% names(fg)))

# --- the programme grouping (ported from figures/panels/_panel_common.R) --------
programme_levels <- c(
  "Curated mitochondrial", "Mitochondrial biogenesis", "OXPHOS", "TCA cycle",
  "Mitochondrial metabolism & dynamics", "MYC signatures", "E2F / cell cycle",
  "Biosynthetic metabolism", "Metabolism, other", "Apoptosis",
  "TF target sets", "Mammary development", "MSigDB Hallmarks")

programme_group <- function(set, category) {
  OX     <- "OXPHOS|COMPLEX_[IV]|_SUBUNITS|ASSEMBLY_FACTORS|ELECTRON_CARRIERS|CRISTAE"
  BIOG   <- "RIBOSOME|CENTRAL_DOGMA|MT_TRNA|MT_RRNA|MTRNA|MTDNA|IMPORT|TRANSLATION"
  BIOSYN <- paste0("NUCLEOTIDE|PURINE|PYRIMIDINE|AMINO|SER_GLY|BCAA|ONE_CARBON|",
                   "PPP|PENTOSE|POLYAMINE|CHOLESTEROL|MEVALONATE|LIPID|FATTY|GLYCOLYSIS")
  constructed <- grepl("_MITO$|_MITO_|MITO_NU|^MITO_|CORE_MITO", set) |
    category %in% c("07_biogenesis_discrimination", "09_biogenesis_apoptosis_intersections")
  g <- ifelse(
    constructed, "Curated mitochondrial",
    ifelse(category == "01_mitocarta",
           ifelse(grepl(BIOG, set), "Mitochondrial biogenesis",
                  ifelse(grepl(OX, set), "OXPHOS", "Mitochondrial metabolism & dynamics")),
    ifelse(category == "04_metabolism",
           ifelse(grepl("OXPHOS|ELECTRON|RESPIRAT", set), "OXPHOS",
                  ifelse(grepl("KREBS|_TCA", set), "TCA cycle",
                         ifelse(grepl(BIOSYN, set), "Biosynthetic metabolism",
                                "Metabolism, other"))),
    ifelse(category == "02_myc_signatures",      "MYC signatures",
    ifelse(category == "05_proliferation",       "E2F / cell cycle",
    ifelse(category == "06_tf_targets",          "TF target sets",
    ifelse(category == "03_mammary_development", "Mammary development",
    ifelse(category == "08_apoptosis",           "Apoptosis",
    ifelse(category == "hallmark_msigdb",        "MSigDB Hallmarks", NA_character_)))))))))
  factor(g, levels = programme_levels)
}

dat <- fg[fg$ranking %in% names(contrast_cols), ]
dat$contrast  <- factor(dat$ranking, levels = names(contrast_cols))
dat$programme <- programme_group(dat$pathway, dat$category)
stopifnot(!anyNA(dat$programme),
          nrow(dat) == 2L * sum(fg$ranking == "myc_6W"))

# rows ordered by the 6-week median, highest at the top
med <- stats::aggregate(NES ~ contrast + programme, data = dat, FUN = stats::median)
m6  <- med[med$contrast == "myc_6W", ]
dat$programme <- factor(dat$programme, levels = as.character(m6$programme[order(m6$NES)]))
med$programme <- factor(med$programme, levels = levels(dat$programme))
med$yi <- as.integer(med$programme)

SHOW_MEDIAN <- FALSE   # TRUE: a thin grey tick at each row's median

p <- ggplot(dat, aes(NES, programme)) +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
  ggbeeswarm::geom_quasirandom(aes(colour = contrast), orientation = "y", width = 0.33,
                               size = 0.85, alpha = 0.55, stroke = 0) +
  (if (SHOW_MEDIAN)
    geom_segment(data = med, aes(x = NES, xend = NES, y = yi - 0.38, yend = yi + 0.38),
                 inherit.aes = FALSE, linewidth = 0.3, colour = "grey45")) +
  facet_wrap(~ contrast, nrow = 1, labeller = as_labeller(contrast_labels)) +
  scale_colour_manual(values = contrast_cols, guide = "none") +
  scale_x_nm(step = 2, top = max(dat$NES), bottom = min(dat$NES), minor = TRUE,
             labels = lab_signed) +
  labs(x = "Normalized enrichment score", y = NULL) +
  theme_nm() +
  theme(axis.line.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing.x = unit(6, "mm"))

save_panel(p, fig = "EDFig1", panel = "c", name = "nes_by_programme", width = 112, height = 80)

# --- numbers for the legend (printed, never drawn) --------------------------------
w <- stats::reshape(dat[, c("pathway", "category", "ranking", "NES")],
                    idvar = c("pathway", "category"), timevar = "ranking", direction = "wide")
rho <- stats::cor(w$NES.myc_6W, w$NES.myc_12W, method = "spearman")
ap6 <- readRDS(here::here("results", "ap6_permutation_null.rds"))$null_table
z6  <- function(cmp) ap6$z[ap6$compartment == cmp & ap6$contrast == "myc_6W"]
med_w <- stats::reshape(med[, c("programme", "contrast", "NES")], idvar = "programme",
                        timevar = "contrast", direction = "wide")
med_w <- med_w[order(-med_w$NES.myc_6W), ]

cat("\nED Fig. 1c -- for the legend\n")
cat(sprintf("  %d gene sets per contrast (%d library sets passing 10-500 genes, plus 50 MSigDB Hallmarks);",
            nrow(w), nrow(w) - 50L),
    if (SHOW_MEDIAN) "one point per set, grey tick = row median\n" else "one point per set\n")
cat("  fGSEA on the unshrunken DESeq2 Wald statistic of each genotype contrast\n")
cat(sprintf("  rank agreement between ages: Spearman %.3f across all %d sets\n", rho, nrow(w)))
cat(sprintf("  median |NES|: 6W %.2f, 12W %.2f\n",
            stats::median(abs(w$NES.myc_6W)), stats::median(abs(w$NES.myc_12W))))
cat("  row medians, 6W / 12W:\n")
for (i in seq_len(nrow(med_w))) {
  cat(sprintf("    %-36s %+.2f / %+.2f\n", as.character(med_w$programme[i]),
              med_w$NES.myc_6W[i], med_w$NES.myc_12W[i]))
}
neg <- stats::aggregate(NES ~ ranking + programme, data = dat, FUN = function(x) sum(x < 0))
cat("  sets with a negative NES, 6W -> 12W (rows with both signs):\n")
for (r in levels(dat$programme)) {
  n6 <- neg$NES[neg$programme == r & neg$ranking == "myc_6W"]
  n12 <- neg$NES[neg$programme == r & neg$ranking == "myc_12W"]
  tot <- sum(dat$programme == r & dat$ranking == "myc_6W")
  if (n6 + n12 > 0) cat(sprintf("    %-36s %d -> %d of %d\n", r, n6, n12, tot))
}
ap <- dat[dat$programme == "Apoptosis", c("pathway", "ranking", "NES", "padj_within_category")]
ap <- stats::reshape(ap, idvar = "pathway", timevar = "ranking", direction = "wide")
fl <- ap[sign(ap$NES.myc_6W) != sign(ap$NES.myc_12W), ]
cat("  Apoptosis sets that change sign (NES 6W -> 12W, within-category BH P at 12W):\n")
for (i in seq_len(nrow(fl))) cat(sprintf("    %-36s %+.2f -> %+.2f  P = %.2g\n", fl$pathway[i],
                                         fl$NES.myc_6W[i], fl$NES.myc_12W[i],
                                         fl$padj_within_category.myc_12W[i]))
cat(sprintf("  sets with BH-adjusted P < 0.05 (within category) at 6W: %d of %d\n",
            sum(dat$padj_within_category[dat$contrast == "myc_6W"] < 0.05, na.rm = TRUE), nrow(w)))
cat(sprintf("  matched null at 6W (script 21), z: MitoCarta %.1f, Metabolism %.1f, Hallmark OXPHOS %.1f,",
            z6("MitoCarta"), z6("Metabolism"), z6("HALLMARK_OXPHOS")),
    sprintf("Hallmark MYC targets V1 %.1f, Hallmark E2F %.1f, Proliferation %.1f\n",
            z6("HALLMARK_MYC_TARGETS_V1"), z6("HALLMARK_E2F_TARGETS"), z6("Proliferation")))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p)
  med_w
  table(dat$programme[dat$contrast == "myc_6W"])
}
