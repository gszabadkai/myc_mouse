# =============================================================================
# fig1_compartment_levels.R -- the four states in absolute units: how much
# transcript each mitochondrial compartment actually holds, per animal
# -----------------------------------------------------------------------------
# SLOT: not currently cited. Built to answer a question the corpus never asked,
# so its sentence comes after it. `figS1_mb_fork_specificity.R` carries the same
# string.
#
# WHAT WAS MISSING. Abundance has only ever been read two ways in this project:
# as a % SHARE of the transcriptome (script 32, Fig. 1E) or as a LOG FOLD CHANGE
# (script 40's content ruler, Figs. 1E lower strip, 1F, 2F). Both are ratios.
# Nothing has ever shown the ABSOLUTE LEVEL -- how large each compartment is, in
# DESeq2-normalised units, in each of the four states, with the per-animal spread.
# `figures/fig03_background_vs_myc.R` panel A comes closest and facets by these
# same seven tiers, but it draws `bg$state_table`, which holds GROUP MEANS: four
# points per tier, no spread and no magnitude.
#
# WHY THAT IS WORTH A PANEL AND NOT JUST A TABLE:
#   * THE MAGNITUDE GAP IS ENORMOUS AND HAS NEVER BEEN DRAWN. Thirteen mtDNA-
#     encoded genes carry ~1.38 million normalised counts against ~65,000 for the
#     152 nuclear OXPHOS genes -- roughly 21x the transcript from a twelfth of the
#     genes. On a shared axis nothing else would be visible, which is exactly the
#     point, and is why the facets are free.
#   * AND THAT SAME facet IS THE ONE WITH NO GENOTYPE EFFECT (p = 0.75) and the
#     largest developmental change of anything here. Myc builds the nuclear
#     compartment and does nothing to the mtDNA-encoded output.
#   * THE CONTENT CLAIM REPRODUCES ON A RULER WITH NO DENOMINATOR IN IT. Every
#     nuclear compartment clears BH < 0.05 on the genotype main effect. That
#     matters because a share divides by a denominator that itself rises with Myc.
#
# THE BRACKET IS THE GENOTYPE MAIN EFFECT, and that is what fixes the group order.
# Groups are drawn GENOTYPE-MAJOR (both wild-type boxes, then both Myc+, which is
# `names(group_cols)` order), so a single bracket spanning the two halves IS the
# additive contrast -- Fig. 1E's construction exactly. The additive model is only
# licensed if the interaction is far from significant, so the interaction is
# fitted first and the panel STOPS if any facet's interaction reaches p < 0.05.
#
# LOG2 Y AXIS, LABELLED IN NORMALISED COUNTS. The model is fitted on log2 (script
# 32's `share_stat_one` idiom), so the drawn scale is the fitted scale and the
# bracket's beta is a distance on the page. Ticks are printed as counts, because
# "37,000 normalised counts" is the readable quantity and "15.2" is not.
#
# Reads (read-only, no re-run):
#   results/state_readings.rds  (script 45 PART G) -- $levels (per sample x set),
#                                  $level_stats (the additive model), $share_agreement
# Output: outputs/figures/panels/fig1_compartment_levels.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

sr_path <- here::here("results", "state_readings.rds")
require_fresher_than(sr_path)
sr <- readRDS(sr_path)

lv <- as.data.frame(sr$levels)
ls_ <- as.data.frame(sr$level_stats)
sa  <- as.data.frame(sr$share_agreement)
ss_ <- as.data.frame(sr$share_summary)

# =============================================================================
# the roster: the seven Level-1 tiers, plus the mtDNA-encoded set on its own
# =============================================================================
# The tier roster is the one fig03 panel A facets by. The mtDNA-encoded pathway is
# held OUT of every tier (it sits under OXPHOS in MitoCarta's own hierarchy) and
# given its own facet, because it is the comparison the panel exists to make and
# because it is the one quantity in this project that cannot be read on a temporal
# axis. The 13 mt- genes appear in NO other MitoPathway, so the nuclear tiers are
# clean by construction, not by subtraction.
TIERS <- names(tier_short)
MT    <- "mtDNA-encoded"
FACETS <- c(TIERS, MT)

lv <- lv[lv$group_set %in% FACETS, ]
ls_ <- ls_[ls_$group_set %in% FACETS, ]
stopifnot(nrow(ls_) == 8L, setequal(ls_$group_set, FACETS),
          all(lv$roster == "tier"))

# --- guards: a stale or half-run object must fail here, not at review ---------
# (1) THE ADDITIVE MODEL NEEDS ITS GATE. A bracket that reports a genotype main
# effect is only meaningful where the genotype x timepoint interaction is not.
# Script 32 checked this before quoting its own additive betas (its arms came out
# at p 0.35-0.81); the same gate is enforced here rather than assumed.
stopifnot(all(ls_$int_p > 0.05))
# (2) the mtDNA facet is the null one, and the panel's whole contrast rests on it
stopifnot(ls_$geno_p[ls_$group_set == MT] > 0.5,
          all(ls_$geno_padj[ls_$group_set != MT] < 0.05))
# (3) six animals per group in every facet
stopifnot(all(table(lv$group_set, lv$group) == 6L))

# Short strip labels: the declared `tier_short`, plus the mtDNA set. `tier_labels`
# is the axis form and is too long for a ~40 mm strip.
facet_lab <- c(tier_short, stats::setNames("mtDNA-encoded (13)", MT))
facet_lab[TIERS] <- sprintf("%s (%d)", tier_short[TIERS],
                            ls_$n_genes[match(TIERS, ls_$group_set)])

# Facet order: the seven nuclear tiers by descending 6W_wt level, then the mtDNA
# set LAST, so the one facet that behaves differently is the one the eye ends on.
ord <- ls_$group_set[ls_$group_set != MT]
ord <- ord[order(-ls_$level_6W_wt[match(ord, ls_$group_set)])]
FACET_ORDER <- c(ord, MT)

lv$facet <- factor(unname(facet_lab[lv$group_set]),
                   levels = unname(facet_lab[FACET_ORDER]))
lv$group <- factor(as.character(lv$group), levels = names(group_cols))
lv$y     <- log2(lv$norm_sum)
stopifnot(!anyNA(lv$facet), !anyNA(lv$group))

# =============================================================================
# the genotype bracket, one per facet
# =============================================================================
# Genotype-major x order means groups 1-2 are the wild-type pair and 3-4 the Myc+
# pair, so a bracket from 1.5 to 3.5 spans the two halves and IS the main effect.
brk <- do.call(rbind, lapply(FACET_ORDER, function(g) {
  s <- ls_[ls_$group_set == g, ]
  d <- lv[lv$group_set == g, ]
  # pad 0.05, not the 0.09 default: eight facets in a 108 mm column cannot spend
  # a fifth of each one on the gap under a bracket, and there is only one bracket
  # per facet so nothing has to stack above it.
  b <- bracket_frame(
    data.frame(x1 = 1.5, x2 = 3.5, level = 1L,
               lab = sprintf("%+.2f, p %s", s$geno_beta,
                             if (s$geno_p < 0.001) "< 0.001"
                             else sprintf("= %.3f", s$geno_p))),
    range(d$y), pad = 0.05)
  b$facet <- factor(unname(facet_lab[g]), levels = levels(lv$facet))
  b$col   <- unname(sig_cols[if (s$geno_padj < 0.05) "sig" else "ns"])
  b$headroom <- attr(b, "headroom")
  b
}))
# a geom_text has no data extent, so a free scale would clip the topmost label
hr <- data.frame(facet = brk$facet, group = factor(names(group_cols)[1],
                                                   levels = names(group_cols)),
                 y = brk$headroom)

# n = 6 per box: every animal is drawn, never summarised away
pts_layer <- function(dat) {
  if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
    ggbeeswarm::geom_quasirandom(data = dat, width = 0.20, size = 0.62,
                                 stroke = 0.15, shape = 21, colour = "grey30")
  } else {
    ggplot2::geom_jitter(data = dat, width = 0.14, height = 0, size = 0.62,
                         stroke = 0.15, shape = 21, colour = "grey30")
  }
}

# Ticks in normalised counts, not in log2. The axis is log2 because the model is,
# but "37,000" is the quantity a reader can hold and "15.2" is not.
count_lab <- function(x) {
  v <- 2^x
  ifelse(v >= 1e6, sprintf("%.1fM", v / 1e6),
         ifelse(v >= 1e3, sprintf("%.0fk", v / 1e3), sprintf("%.0f", v)))
}

p <- ggplot2::ggplot(lv, ggplot2::aes(group, y)) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = group), outlier.shape = NA,
                        width = 0.62, alpha = 0.30, colour = "grey35",
                        linewidth = 0.22) +
  pts_layer(lv) +
  bracket_layers(brk, size = 1.5, linewidth = 0.2) +
  ggplot2::geom_blank(data = hr, ggplot2::aes(x = group, y = y)) +
  ggplot2::facet_wrap(~ facet, ncol = 2, scales = "free_y") +
  ggplot2::scale_fill_manual(values = group_cols, labels = group_labels,
                             name = NULL) +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_y_continuous(labels = count_lab,
                              breaks = scales::breaks_pretty(3)) +
  ggplot2::labs(x = NULL, y = "compartment total (DESeq2 normalised counts)") +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    axis.text.x       = ggplot2::element_blank(),
    axis.ticks.x      = ggplot2::element_blank(),
    legend.position   = "bottom",
    legend.key.size   = ggplot2::unit(2.6, "mm"),
    legend.margin     = ggplot2::margin(-1.5, 0, 0, 0, "mm"),
    strip.text        = ggplot2::element_text(size = 5.4, face = "bold",
                                              margin = ggplot2::margin(0.6, 0, 0.9, 0, "mm")),
    panel.spacing.x   = ggplot2::unit(2.2, "mm"),
    panel.spacing.y   = ggplot2::unit(1.4, "mm"),
    plot.margin       = ggplot2::margin(1.5, 2, 1, 1.5, "mm")) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
S  <- function(g, col) ls_[[col]][ls_$group_set == g]
fmtp <- function(x) if (x < 0.001) "< 0.001" else sprintf("= %.3f", x)
lo <- ls_[ls_$group_set != MT, ]
lo <- lo[order(lo$geno_beta), ]

LEGEND <- panel_legend(
  slot = "not currently cited",
  what = paste0(
    "How much transcript each mitochondrial compartment holds, in DESeq2-",
    "normalised units, in each of the four states. One facet per MitoCarta ",
    "Level-1 tier with the thirteen mtDNA-encoded genes on their own, six animals ",
    "per box, boxes ordered genotype-major so the bracket spanning the two halves ",
    "is the genotype main effect. The y axis is log2 and its ticks are printed as ",
    "counts; facets are free because the compartments span two orders of ",
    "magnitude."),
  detail = c(
    sprintf("n = 6 animals per group, all drawn. Values are the SUM of DESeq2-normalised counts over the genes of each tier: size factors from the fitted group model applied to the FULL count matrix, not to the fitted object, whose expression filter would have silently under-counted every compartment. Membership resolved through functions/reconcile_gene_symbols.R."),
    sprintf("THE BRACKET IS THE GENOTYPE MAIN EFFECT from lm(log2(total) ~ timepoint + myc_status) -- script 32's model, column for column, so this table and the published share table read side by side. Its label is the beta (a log2 fold change) and its unadjusted p; red where BH-adjusted p < 0.05 across the eight facets."),
    sprintf("THE ADDITIVE MODEL IS LICENSED, AND THE PANEL CHECKS IT. Every facet's genotype x timepoint interaction is far from significance (p %.2f to %.2f), so the genotype effect does not differ between the ages and one bracket per facet is the right summary. The panel stops if any interaction reaches p < 0.05.",
            min(ls_$int_p), max(ls_$int_p)),
    sprintf("EVERY NUCLEAR COMPARTMENT RISES WITH MYC, on a ruler with no compartment denominator in it: %s. All seven clear BH < 0.05.",
            paste(sprintf("%s %+.2f (p %s)", unname(tier_short[lo$group_set]),
                          lo$geno_beta, vapply(lo$geno_p, fmtp, character(1))),
                  collapse = "; ")),
    sprintf("AND THE mtDNA-ENCODED SET DOES NOT: %+.2f, p %s -- the only facet whose bracket is not significant, and the comparison this panel exists to draw. Thirteen genes carry %s normalised counts against %s for the %d nuclear OXPHOS genes, about %.0fx the transcript from a twelfth of the genes. Myc builds the nuclear compartment and leaves the mtDNA-encoded output alone.",
            S(MT, "geno_beta"), fmtp(S(MT, "geno_p")),
            format(round(S(MT, "level_6W_wt")), big.mark = ","),
            format(round(S("OXPHOS", "level_6W_wt")), big.mark = ","),
            S("OXPHOS", "n_genes"),
            S(MT, "level_6W_wt") / S("OXPHOS", "level_6W_wt")),
    sprintf("THE ONE COMPARTMENT THE WILD-TYPE GLAND WITHDRAWS FROM IS THE RESPIRATORY ONE, and it reproduces here on this third ruler: over the wild-type 6 to 12 week window the OXPHOS tier moves %+.3f, the largest negative of the seven, against %+.3f for metabolism and %+.3f for the central dogma. The mtDNA facet moves %+.3f, the largest change of anything on the panel and in the opposite direction.",
            S("OXPHOS", "wt_temporal_beta"), S("Metabolism", "wt_temporal_beta"),
            S("Mitochondrial central dogma", "wt_temporal_beta"),
            S(MT, "wt_temporal_beta")),
    sprintf("HOW MUCH OF THIS IS NEW, MEASURED RATHER THAN ASSERTED. A summed normalised count and a %% share are close relatives -- they differ only by each sample's total. On the SAME genes the two correlate at median r = %.2f (range %.2f to %.2f across the eight facets). What is new is the absolute magnitude, the per-animal spread, and a genotype test with no denominator in it.",
            ss_$median_r_own, min(sa$r_own[sa$group_set %in% FACETS], na.rm = TRUE),
            max(sa$r_own[sa$group_set %in% FACETS], na.rm = TRUE)),
    sprintf("AND IT CORROBORATES THE PUBLISHED CONTENT CLAIM RATHER THAN COMPETING WITH IT. Against script 32's share model on the same panels the genotype betas agree at r = %.3f and are LARGER here on %d of %d matched panels (median offset %+.3f). That is the expected direction: the share divides by the nuclear transcriptome, which itself rises with Myc, so it subtracts part of the effect. Script 32 already states that its share effect is a LOWER BOUND on content; this puts a number on the gap.",
            ss_$beta_r, ss_$n_levels_larger, ss_$n_matched,
            ss_$median_beta_offset)),
  bounds = c(
    "BATCH = TIMEPOINT. The 6W and 12W cohorts were extracted as two separate batches, so every WITHIN-GENOTYPE comparison across the two ages on this panel is DESCRIBED, not claimed. The genotype bracket is the clean contrast: genotype is balanced within each batch, which is precisely why it, and not the age difference, is the thing the panel marks.",
    "THE mtDNA FACET IS THREE-WAY CONFOUNDED and is drawn as a comparison, not as a measurement. The mtDNA-encoded read fraction reflects real content, the proliferation denominator and dissociation leak at once, and it is time-associated rather than genotype-associated. Its FLAT GENOTYPE BRACKET is the readable half; its large developmental change is the half that is not.",
    "A TRANSCRIPT TOTAL IS NOT AN ORGANELLE COUNT. These are normalised read sums, not protein and not mitochondrial volume. Under global RNA amplification by Myc a constant total would already imply more mitochondria per cell, so as a per-cell content statement this is a LOWER BOUND, exactly as script 32's share is.",
    "ABSOLUTE LEVELS ARE NOT COMPARABLE BETWEEN COMPARTMENTS as biology: gene length is constant across samples but not across sets, so a set's total is comparable ACROSS the four groups and not against another set's total. The 21x mtDNA gap is quoted as what the sequencer sees, which is the claim being made about it.",
    "MitoCarta sets are membership-loose. A tier's shift can be carried by a few highly expressed genes -- the summed ruler weights by expression by construction, which is exactly the effect the companion OXPHOS-subunit panel dissects. The per-gene unweighted ruler gives systematically different values on the respiratory arm.",
    "n = 6 per cell. The brackets are the powered contrast in this design; nothing about the between-age differences on this panel is a confirmatory test."),
  source = c(
    "results/state_readings.rds (scripts/45_state_readings.R PART G) -- $levels (per sample x compartment totals), $level_stats (the additive genotype model, script 32's share_stat_one idiom), $share_agreement and $share_summary (the measured overlap with the published share ruler)",
    "size factors: DESeq2::sizeFactors(results/dds_group_run.rds), applied to results/count_matrix.rds (unfiltered)",
    "the tier partition and the mtDNA-encoded set: results/mitopps_scores.rds (scripts/08_mitoPPS_analysis.R), from Mouse.MitoCarta3.0.xls Sheet 4",
    "the published share comparator: results/mito_content_proxies.rds (scripts/32_mito_content_proxies.R), $shares and $share_stats on the share_nomt denominator"))

save_panel_p(p, "fig1_compartment_levels", height = 100)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the model behind every bracket, with the interaction gate
  ls_[, c("group_set", "n_genes", "level_6W_wt", "geno_beta", "geno_p",
          "geno_padj", "int_p", "wt_temporal_beta", "cross_log2")] |>
    (\(x) x[order(-x$geno_beta), ])() |> print(row.names = FALSE, digits = 3)

  ## the named arms, which script 45 also computes but this panel does not draw
  as.data.frame(sr$level_stats) |>
    (\(x) x[x$roster == "arm", c("group_set", "n_genes", "level_6W_wt",
                                 "geno_beta", "geno_padj", "cross_log2")])() |>
    print(row.names = FALSE, digits = 3)

  ## how much the levels ruler and the published share ruler overlap
  sa |> print(row.names = FALSE, digits = 3)
  ss_ |> print(row.names = FALSE, digits = 3)

  ## per-animal totals for one compartment, to see the spread the boxes summarise
  lv[lv$group_set == "OXPHOS", c("sample", "group", "norm_sum")] |>
    (\(x) x[order(x$group), ])() |> print(row.names = FALSE, digits = 6)
}
