# =============================================================================
# fig2_wt_mito_contraction.R -- the normal gland withdraws from the respiratory
# chain, on both rulers, while the rest of the compartment rises
# -----------------------------------------------------------------------------
# SLOT: Fig. 2F.
#
#   "MECs withdraw from the respiratory chain; OXPHOS subunit LFC AND MitoPPS
#    drop across all complexes while biogenesis pathways remain relatively
#    stable ... the maturing gland upregulates amino-acid and lipid catabolism."
#
# THREE CLAUSES AND TWO RULERS, so the panel is eleven MitoPathways ranked on the
# wild-type 6->12W contrast, drawn twice:
#
#   CONTENT   set-average raw log2FC. What the compartment HAS.
#   PRIORITY  mitoPPS (Monzel 2025), pairwise-ratio and content-blind. What the
#             compartment SPENDS ITS BUDGET ON. A uniform scaling cancels.
#
# They are separate facets on separate scales because they are separate units:
# putting them on one axis would invite a reader to compare a log2 fold change
# with a ratio score. What is comparable is the ORDER, and the order is the
# result -- across all 143 non-mtDNA MitoPathways the two rulers agree at
# Spearman 0.82, so the withdrawal is not a normalisation artifact of either.
#
# THE PICTURE, top to bottom:
#   CIV, CI, CIII, CV subunits    down on both rulers, the bottom 4% of the
#                                 compartment on both
#   mitoribosome, central dogma   at the middle of the compartment: biogenesis
#   OXPHOS assembly factors       AT ZERO -- the internal control, because these
#                                 are the assembly factors OF THE SAME COMPLEXES
#   CII subunits                  UP, and it is the one complex with no
#                                 mtDNA-encoded subunit (n = 4; see the legend)
#   lipid, amino acid, fatty acid oxidation   up, the top 7-24%
#
# AND THE BACKGROUND IS RISING, which is what makes "withdraw" the right verb:
# the median MitoPathway GAINS content over this window (+0.041, 72% above zero).
# The respiratory arm is not falling with the compartment, it is falling against
# it.
#
# THE mtDNA-ENCODED OXPHOS SUBUNITS ARE EXCLUDED, exactly as Figs. 1E and 1F
# exclude them and as script 40 does (`is_mtdna`). They are the largest single
# movement here (+0.63 content, +0.53 priority) and they are the one quantity in
# this project that cannot be read on a temporal axis -- see the legend block.
#
# BATCH = TIMEPOINT (CLAUDE.md): the 6W and 12W cohorts were extracted as two
# batches, so every value on this panel is DESCRIBED, not claimed.
#
# Reads (read-only, no re-run):
#   results/background_vs_myc.rds              (script 40) -- $ruler, both rulers
#                                                 x four contrasts x 144 pathways
#   results/substrate_specificity_tradeoff.rds (script 43) -- $comparator,
#                                                 $comparator_priority, $wt_null,
#                                                 $defs (arm map + stale guard)
# Output: outputs/figures/panels/fig2_wt_mito_contraction.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

bv_path <- here::here("results", "background_vs_myc.rds")
ss_path <- here::here("results", "substrate_specificity_tradeoff.rds")
require_fresher_than(bv_path)
require_fresher_than(ss_path)

bv <- readRDS(bv_path)
ss <- readRDS(ss_path)

# Scripts 40 and 43 save TIBBLES, and a tibble's `[` returns a tibble rather than
# a scalar -- which silently poisons every sprintf downstream (Fig. 1E's trap).
r   <- as.data.frame(bv$ruler)
cmp <- as.data.frame(ss$comparator)
cpr <- as.data.frame(ss$comparator_priority)
wtn <- as.data.frame(ss$wt_null)
arm <- as.data.frame(ss$defs$arms)

# `n_genes` and the content columns carry names in the saved object; strip them or
# a data.frame row lookup returns a named vector and sprintf prints the name too.
r$n <- as.integer(unname(r$n_genes))
r$c <- as.numeric(unname(r$c_tn))       # content, wild-type 6->12W
r$p <- as.numeric(unname(r$p_tn))       # priority, wild-type 6->12W
stopifnot(all(c("pathway", "tier", "is_mtdna") %in% names(r)), nrow(r) == 144L)

# --- guard: a stale results object must fail here, not at review --------------
# Script 43's own positive control, as figures/fig04 uses it: the wild-type
# OXPHOS-subunit value reproduces Issue #4's -0.2548.
ox_obs <- cmp$c_wt_time[cmp$arm == "OXPHOS subunits"]
if (!isTRUE(abs(ox_obs - ss$defs$oxphos_wt_reference) < 0.02))
  stop(sprintf("fig2F: stale results object -- OXPHOS-subunit wild-type %.4f against the reference %.4f. Re-run script 43.",
               ox_obs, ss$defs$oxphos_wt_reference))

# --- the compartment, and the two scripts agreeing on it ----------------------
# 143, not 144: the synthetic mtDNA-encoded pathway is dropped exactly as script
# 40 drops it, so this panel, Fig. 1E's lower strip and every regression in the
# corpus describe the same set of pathways.
r143 <- r[!r$is_mtdna, ]
stopifnot(nrow(r143) == 143L)

# Scripts 40 and 43 computed the wild-type arms independently; they must agree.
# `defs$arms` is script 43's own arm -> ruler-pathway map, so this is an identity
# check and not a re-derivation.
chk <- merge(cmp[, c("arm", "c_wt_time")], arm[, c("arm", "ruler_pathway")], by = "arm")
chk <- chk[!is.na(chk$ruler_pathway), ]
chk$c_ruler <- r$c[match(chk$ruler_pathway, r$pathway)]
cpr$p_ruler <- r$p[match(cpr$ruler_pathway, r$pathway)]
stopifnot(nrow(chk) >= 7L,
          max(abs(chk$c_wt_time - chk$c_ruler)) < 1e-6,
          max(abs(cpr$prio_wt_time - cpr$p_ruler)) < 1e-6)

# =============================================================================
# the eleven rows
# =============================================================================
# One row per clause of the sentence, plus the control the OXPHOS clause needs.
# Named here rather than derived, because the sentence names them -- and each is a
# MitoCarta pathway of the ruler, so nothing is re-aggregated in the figure layer.
#
#   respiratory  the five complexes' SUBUNITS (the sentence says "subunit"), and
#                the assembly factors of the same complexes as the control
#   biogenesis   the central-dogma tier and the mitoribosome
#   catabolic    fatty acid oxidation, amino acid, lipid
ROWS <- c("CIV subunits", "CI subunits", "CIII subunits", "CV subunits",
          "Mitochondrial ribosome", "OXPHOS assembly factors",
          "Mitochondrial central dogma", "CII subunits",
          "Lipid metabolism", "Amino acid metabolism", "Fatty acid oxidation")

w <- r143[match(ROWS, r143$pathway), c("pathway", "tier", "n", "c", "p")]
stopifnot(!anyNA(w$c), !anyNA(w$p), nrow(w) == length(ROWS))

# The four proton-circuit complexes fall on BOTH rulers and CII does not. Asserted
# so that a re-run cannot flip the exception without stopping the panel.
FOUR <- c("CI subunits", "CIII subunits", "CIV subunits", "CV subunits")
stopifnot(all(w$c[w$pathway %in% FOUR] < 0), all(w$p[w$pathway %in% FOUR] < 0),
          w$c[w$pathway == "CII subunits"] > 0, w$p[w$pathway == "CII subunits"] > 0)

# Display names: `tier_labels` already shortens the central-dogma tier, so the
# shortening is taken from the declared vector rather than typed again here.
disp <- ifelse(w$pathway %in% names(tier_labels), tier_labels[w$pathway], w$pathway)
w$row <- sprintf("%s (%d)", disp, w$n)
# Discrete y draws level 1 at the BOTTOM, so ordering by DESCENDING content puts
# the deepest withdrawal at the top -- the order the sentence reads in, and the
# order figures/fig04 panel A uses for the same quantity.
w$row <- factor(w$row, levels = w$row[order(-w$c)])

# Strip labels are short because a facet strip clips at the panel edge (the trap
# Fig. 1E hit with its arm names); "set-average" lives in the legend block.
RULERS <- c("content  (log2FC)", "priority  (mitoPPS)")
long <- rbind(data.frame(row = w$row, ruler = RULERS[1], value = w$c),
              data.frame(row = w$row, ruler = RULERS[2], value = w$p))
long$ruler <- factor(long$ruler, levels = RULERS)

# =============================================================================
# the panel
# =============================================================================
p <- ggplot2::ggplot(long, ggplot2::aes(value, row)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey45") +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = value, yend = row),
                        linewidth = 0.4, colour = "grey55") +
  ggplot2::geom_point(size = 1.4, colour = "grey15") +
  ggplot2::facet_wrap(~ ruler, nrow = 1, scales = "free_x") +
  ggplot2::scale_x_continuous(labels = lab_signed,
                              expand = ggplot2::expansion(mult = 0.10)) +
  ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.7)) +
  ggplot2::labs(x = "6>12W_wt", y = NULL) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    strip.text  = ggplot2::element_text(face = "plain", size = 6, hjust = 0,
                                        margin = ggplot2::margin(0, 0, 1, 0, "mm")),
    strip.clip  = "off",
    axis.text.y = ggplot2::element_text(size = 6),
    # No y axis line or ticks: with two facets theme_classic draws one down the
    # left of EACH, and the second reads as an unlabelled second axis. The zero
    # line inside each panel is the reference the rows are read against.
    axis.line.y  = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.spacing.x = ggplot2::unit(3, "mm"),
    plot.margin = ggplot2::margin(1.5, 2.5, 1, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
val <- function(pw, col) w[[col]][w$pathway == pw]
pct <- function(pw, col) 100 * mean(r143[[col]] < val(pw, col))
# "1th" and "3th" are what sprintf's %.0f plus a hard-coded "th" produces; the
# suffix has to follow the number.
ordinal <- function(x) {
  n <- round(x)
  suf <- if (n %% 100 %in% 11:13) "th" else
    switch(as.character(n %% 10), "1" = "st", "2" = "nd", "3" = "rd", "th")
  paste0(n, suf)
}
line <- function(pw, label = pw)
  sprintf("%s %+.3f content (%s percentile of the compartment) and %+.3f priority (%s)",
          label, val(pw, "c"), ordinal(pct(pw, "c")),
          val(pw, "p"), ordinal(pct(pw, "p")))
null_of <- function(a) wtn$percentile[wtn$arm == a]
mt <- r[r$is_mtdna, ]

LEGEND <- panel_legend(
  slot = "Fig. 2F",
  what = paste0(
    "Eleven MitoPathways over the wild-type 6 to 12 week window, on two rulers. ",
    "Left, content: the set-average raw log2 fold change, what the compartment ",
    "has. Right, priority: the mitoPPS pairwise-ratio score, which is blind to ",
    "content and reads how the compartment divides its budget. Rows are ranked ",
    "by content, deepest withdrawal at the top; the number after each name is ",
    "the pathway's gene count. The two panels are on separate scales because ",
    "they are separate units."),
  detail = c(
    sprintf("n = 6 wild-type animals per timepoint. Both rulers are computed over %d MitoPathways (the synthetic mtDNA-encoded pathway excluded, as in script 40 and Figs. 1E and 1F), and the percentiles quoted below are positions within that compartment.",
            nrow(r143)),
    sprintf("THE TWO RULERS AGREE, which is the reason both are drawn: across all %d pathways they correlate at Spearman %.2f (Pearson %.2f). A drop in content alone could be a normalisation effect; a drop in a content-blind ratio alone could be a reshuffle inside a growing compartment. Together they are a withdrawal.",
            nrow(r143), stats::cor(r143$c, r143$p, method = "spearman"),
            stats::cor(r143$c, r143$p)),
    sprintf("THE BACKGROUND IS RISING, which is what makes \"withdraw\" the right verb: the median MitoPathway GAINS content over this window (%+.3f, %.0f%% of the %d above zero), while on the content-blind ruler the compartment is by construction near zero (median %+.3f). The respiratory arm is not falling with the compartment, it is falling against it.",
            stats::median(r143$c), 100 * mean(r143$c > 0), nrow(r143),
            stats::median(r143$p)),
    sprintf("THE RESPIRATORY ARM. %s; %s; %s; %s. Pooled over all 87 nuclear-encoded OXPHOS subunits the content effect is %+.4f, which is percentile 0.0 of 2000 expression-matched random gene sets (script 43's empirical null, p < 0.0005), and the priority effect is %+.4f.",
            line("CIV subunits"), line("CI subunits"), line("CIII subunits"),
            line("CV subunits"),
            cmp$c_wt_time[cmp$arm == "OXPHOS subunits"],
            cpr$prio_wt_time[cpr$arm == "OXPHOS subunits"]),
    sprintf("THE INTERNAL CONTROL IS THE ASSEMBLY FACTORS OF THE SAME COMPLEXES, and they do not move: %s. Against the same expression-matched null they sit at percentile %.1f -- the middle of the distribution -- where the subunits sit at 0.0. What the gland withdraws is the structural stoichiometry of the chain, not the machinery that builds it.",
            line("OXPHOS assembly factors"), null_of("OXPHOS assembly")),
    sprintf("BIOGENESIS IS STABLE, on both rulers and against the matched null: %s; %s (the mitoribosome sits at percentile %.1f of the matched null). This is the clause that separates the result from a general shrinkage of the mitochondrial programme.",
            line("Mitochondrial central dogma", "central dogma"),
            line("Mitochondrial ribosome", "mitoribosome"), null_of("mitoribosome")),
    sprintf("CATABOLISM RISES: %s; %s; %s. Amino-acid and lipid metabolism are the two arms that beat their expression-matched nulls in the UPWARD direction, at percentiles %.1f and %.1f of 2000 draws.",
            line("Fatty acid oxidation"), line("Amino acid metabolism"),
            line("Lipid metabolism"), null_of("amino-acid metabolism"),
            null_of("lipid metabolism")),
    sprintf("\"ACROSS ALL COMPLEXES\" HAS ONE EXCEPTION AND IT IS INFORMATIVE: %s. Complex II is the only respiratory complex with no mtDNA-encoded subunit -- it is not in the proton circuit and it is also a TCA enzyme -- and it is the only one that does not fall. But the set is FOUR genes, and MitoCarta sets are membership-loose, so this is a direction to note and not a mechanism to claim.",
            line("CII subunits")),
    sprintf("EXCLUDED, AND IT IS THE LARGEST MOVEMENT IN THE COMPARTMENT: the 13 mtDNA-encoded OXPHOS subunits rise %+.3f on content and %+.3f on priority. They are dropped by the same rule Figs. 1E and 1F use. The reason is not tidiness: the mtDNA-encoded read fraction is confounded three ways in these data -- real content, the proliferation denominator, and dissociation leak -- and it is the one quantity that is time-associated rather than genotype-associated, so a temporal contrast is exactly where it cannot be read. The nuclear-encoded arm falling while the mtDNA arm rises is a mitonuclear discordance if it is real; on this axis it is not adjudicable.",
            as.numeric(unname(mt$c_tn)), as.numeric(unname(mt$p_tn)))),
  bounds = c(
    "BATCH = TIMEPOINT. The 6W and 12W cohorts were extracted as two separate batches, so every value on this panel is DESCRIBED, not claimed. Two things mitigate it and neither removes it: the withdrawal is SPECIFIC within the compartment (the assembly factors of the same complexes, measured in the same libraries on the same batches, do not move), and it appears on a content-blind ratio ruler as well as on the content one.",
    "No per-pathway p-value is drawn, and there is none in the saved object for this contrast: script 40 carries adjusted p-values for the genotype contrasts only. The uncertainty that IS available is script 43's expression-matched null, which exists for five of the eleven rows and is quoted above.",
    "mitoPPS is RELATIVE BY CONSTRUCTION: a pathway can be demoted while its absolute expression rises. On this panel both rulers point the same way for the respiratory arm, so no such reading is needed -- but a priority value must never be reported as a fall in expression.",
    "MitoCarta sets are membership-loose and several rows here are small (CII subunits 4 genes, CIII subunits 9). A large effect on a four-gene set is one gene, not a module, which is why the gene count is on the face of the panel.",
    "Within-compartment percentiles are a RANK AMONG MITOPATHWAYS, not a significance statement; the expression-matched null percentiles quoted for five rows are the empirical test and come from script 43."),
  source = c(
    "results/background_vs_myc.rds (scripts/40_background_vs_myc_decomposition.R) -- $ruler, content (c_tn) and priority (p_tn) on the wild-type temporal contrast for 144 MitoPathways",
    "results/substrate_specificity_tradeoff.rds (scripts/43_substrate_specificity_and_tradeoff.R) -- $comparator, $comparator_priority and $wt_null for the arm-level values and their expression-matched nulls; $defs$arms is the arm-to-pathway map the identity check uses",
    "mitoPPS: Monzel et al. 2025, external/mitotyping/; pairwise-ratio scores on linear-scale DESeq2 normalised counts (scripts/08_mitoPPS_analysis.R)"))

save_panel_p(p, "fig2_wt_mito_contraction", height = 56)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the eleven rows with both rulers and their compartment percentiles
  data.frame(pathway = w$pathway, n = w$n, content = w$c, priority = w$p,
             c_pct = vapply(w$pathway, pct, numeric(1), "c"),
             p_pct = vapply(w$pathway, pct, numeric(1), "p")) |>
    (\(x) x[order(x$content), ])() |> print(row.names = FALSE, digits = 3)

  ## the whole OXPHOS tier, which is where the subunit / assembly split shows
  r[r$tier == "OXPHOS", c("pathway", "n", "c", "p")] |>
    (\(x) x[order(x$c), ])() |> print(row.names = FALSE, digits = 3)

  ## the central-dogma tier -- nothing in it moves
  r[r$tier == "Mitochondrial central dogma", c("pathway", "n", "c", "p")] |>
    (\(x) x[order(x$c), ])() |> print(row.names = FALSE, digits = 3)

  ## the top and bottom of the compartment on the content ruler
  r143[order(-r143$c), c("pathway", "tier", "n", "c", "p")] |> head(12) |>
    print(row.names = FALSE, digits = 3)
  r143[order(r143$c), c("pathway", "tier", "n", "c", "p")] |> head(12) |>
    print(row.names = FALSE, digits = 3)

  ## the arm-level view with the expression-matched nulls, which the panel does
  ## not draw -- script 43's own table, and figures/fig04 panel A draws it
  merge(cmp[, c("arm", "n_genes", "c_wt_time")],
        wtn[, c("arm", "percentile", "p_emp_lower")], by = "arm") |>
    (\(x) x[order(x$c_wt_time), ])() |> print(row.names = FALSE, digits = 3)

  ## the same eleven rows on the GENOTYPE contrast at six weeks, which is Fig. 1F's
  ## question -- the two are near mirror images and that is Fig. 2G's subject
  data.frame(pathway = w$pathway,
             myc_6W_content = r$c_m6[match(w$pathway, r$pathway)],
             wt_time_content = w$c) |> print(row.names = FALSE, digits = 3)
}
