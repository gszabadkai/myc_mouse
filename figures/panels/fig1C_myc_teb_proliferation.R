# =============================================================================
# fig1C_myc_teb_proliferation.R -- what Myc does at each age, and what the
# gland does over time
# -----------------------------------------------------------------------------
# SUPPORTS (Results, "Myc drives early breast tumourigenesis by inducing OXPHOS
# and biosynthetic pathways", paragraph 1):
#   "an endogenous Myc program contributed to the pubertal TEB state in WT
#    animals, promoting proliferation (Fig. 1C), evidenced by the reduction of
#    the canonical Myc, TEB and proliferation signatures in the adult stage
#    (WT 6W -> 12W). The addition of the Myc transgene amplified the TEB and
#    proliferation effects and suppressed the BMYO lineage in favor of
#    differentiation to LHS. However, the Myc+ 6W -> 12W TEB-proliferation
#    trajectory also showed a reduction similar to the WT trend"
#
# WHAT CHANGED, AND WHY (author, 2026-07-30)
#
# 1. THE POOLED GENOTYPE COLUMN IS GONE. "Myc+ vs WT" was the genotype main
#    effect, averaged over both ages, and it was hiding the thing the sentence is
#    about. The two genotype columns are now myc_6W and myc_12W, and they differ a
#    lot: TEB-ductal is +1.02 SD at 6 weeks and +0.23 SD at 12, so a Myc effect on
#    the TEB axis exists exactly where the phenotype is and has gone by 12 weeks.
#    The pooled effect (+0.63 SD, p = 0.14) is the average of those two and reads
#    as a null. THIS RESOLVES the objection recorded in PANELS.md that "the
#    transgene amplified the TEB effect" had no support: it has 6-week support.
#
# 2. THE TWO TIMELINES ARE STILL SIDE BY SIDE, AND THEY REALLY DO LOOK ALIKE on
#    the MYC signatures (-0.60 SD in wild type, -0.68 SD in Myc+). That is not a
#    drawing problem, it is the result: in cohort-relative GSVA space the genotype
#    gap on the MYC programme is the same at both ages (2.36 vs 2.28 SD), i.e.
#    THERE IS NO ATTENUATION HERE. The attenuation result lives on a different
#    ruler -- DESeq2 effect sizes and their fGSEA enrichment, where the Myc effect
#    rescales x0.55 (scripts 29-31, 40) -- and this panel must not be read as
#    evidence for or against it. What DOES attenuate here is the developmental
#    arm: TEB +1.02 -> +0.23, and every lineage-identity axis in block 3 weakens.
#    That contrast (MYC programme constant, developmental identity fading) is the
#    honest reading and it is the one that sets up the attenuation section.
#
# 3. BLOCK 3 IS NEW: the Gray 2023 lineage-identity axes. For each major type
#    (AP alveolar progenitor, BA basal, HS hormone-sensing) Gray's CHEA3 analysis
#    gives the regulator programmes of the HIGH-expressing and the LOW-expressing
#    cells of that type; HE is the differentiated end of the axis and LE the
#    lineage-suppressed end (docs/library_reference/
#    Gray_et_al_developmental_TFS_selection.md; LE = low-EXPRESSING, see the
#    project memory note -- the provenance table's gloss is wrong). Myc lowers
#    every HE composite and raises every LE composite, in all three lineages, at
#    both ages and most strongly at 6 weeks: this is the de-differentiation axis
#    the narrative needs. The signs are opposite within a lineage, which is the
#    one thing a single common-mode axis cannot produce -- see the bounds, because
#    the LE half sits almost on top of that axis (r 0.85-0.90).
#
# Input:  results/myc_endogenous_amplification.rds  (script 27 -- $prog_stats, $prog_wide)
#         results/dev_program_myc_integration.rds   (script 26 -- $annot, $state_stats)
#         results/gsva_scores.rds                   (script 15 -- per-sample scores)
# Output: outputs/figures/panels/fig1C_myc_teb_proliferation.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

prog_path <- here::here("results", "myc_endogenous_amplification.rds")
dev_path  <- here::here("results", "dev_program_myc_integration.rds")
require_fresher_than(prog_path)
require_fresher_than(dev_path)

prog <- readRDS(prog_path)
dev  <- readRDS(dev_path)
gsva <- readRDS(here::here("results", "gsva_scores.rds"))
stopifnot(all(c("prog_stats", "prog_wide") %in% names(prog)),
          all(c("annot", "state_stats") %in% names(dev)))

ps <- prog$prog_stats
ss <- dev$state_stats

scores      <- gsva$scores
sample_meta <- as.data.frame(gsva$sample_meta)[colnames(scores), , drop = FALSE]

# --- the per-sample matrix: 15 programmes on one footing ----------------------
# script 27's six composites come ready-made per sample; the three MEC states and
# the six Gray lineage lanes are rebuilt here from the same GSVA scores, using the
# same definition (mean over the member sets).
pw <- as.data.frame(prog$prog_wide)
rownames(pw) <- pw$sample
pw <- pw[colnames(scores), , drop = FALSE]
stopifnot(identical(rownames(pw), colnames(scores)))

main_sets <- dev$annot |>
  dplyr::filter(as.character(state) %in% c("BMYO", "LASP", "LHS"))
state_mat <- vapply(c("BMYO", "LASP", "LHS"),
                    function(s) composite_of(scores,
                                             main_sets$set[as.character(main_sets$state) == s]),
                    numeric(ncol(scores)))

# Gray HE / LE lanes: TFT_<TF>_GRAY_<lineage>_<HE|LE>. The anchored "$" is what
# excludes the _MITO promotions (TFT_..._AP_LE_MITO), which are that factor's
# programme INTERSECTED WITH MitoCarta and so are mitochondrial by construction --
# the project note is to always read the mito-removed version of this axis.
lanes <- list()
for (lin in c("AP", "BA", "HS")) for (end in c("HE", "LE")) {
  lanes[[paste(lin, end)]] <-
    grep(paste0("^TFT_.*_GRAY_", lin, "_", end, "$"), rownames(scores), value = TRUE)
}
stopifnot(lengths(lanes) >= 10, !any(grepl("_MITO$", unlist(lanes))))
lane_mat <- vapply(lanes, function(s) composite_of(scores, s), numeric(ncol(scores)))

Y <- cbind(pw[, c("myc", "felsher", "hallmark_v2", "myc_in_teb", "prolif",
                  "teb_ductal")],
           state_mat, lane_mat)

# --- the contrasts, one table -------------------------------------------------
tab <- do.call(rbind, lapply(names(Y), function(nm) {
  out <- contrast_table(Y[[nm]], sample_meta$timepoint, sample_meta$myc_status,
                        sample_meta$group)
  out$key <- nm
  out
}))

# --- assertions: the numbers drawn are script 26's and 27's -------------------
# 27's six composites: both timelines, the interaction and the SD must reproduce.
get1 <- function(k, cn) tab$effect[tab$key == k & tab$contrast == cn]
for (k in ps$program) {
  stopifnot(
    abs(get1(k, "6>12W_wt")  - ps$wt_d[ps$program == k])                          < 1e-8,
    abs(get1(k, "6>12W_myc") - ps$mycpos_shift[ps$program == k] /
                               ps$within_sd[ps$program == k])                     < 1e-8,
    abs(get1(k, "interaction") - ps$int_beta[ps$program == k] /
                                 ps$within_sd[ps$program == k])                   < 1e-8,
    abs(unique(tab$within_sd[tab$key == k]) - ps$within_sd[ps$program == k])       < 1e-8)
}
# 26's three states: d6 / d12 ARE the per-timepoint genotype gaps in SD units.
for (k in as.character(ss$state)) {
  stopifnot(abs(get1(k, "myc_6W")  - ss$d6[as.character(ss$state) == k])  < 1e-8,
            abs(get1(k, "myc_12W") - ss$d12[as.character(ss$state) == k]) < 1e-8)
}

# --- labels and blocks -------------------------------------------------------
prog_meta <- data.frame(
  key = c("myc", "felsher", "hallmark_v2", "myc_in_teb", "prolif",
          "teb_ductal", "BMYO", "LASP", "LHS",
          "AP HE", "AP LE", "BA HE", "BA LE", "HS HE", "HS LE"),
  label = c("MYC signatures", "Felsher", "Hallmark MYC V2", "MYC-in-TEB",
            "Proliferation",
            "TEB - ductal", "BMYO", "LASP", "LHS",
            "AP HE", "AP LE", "BA HE", "BA LE", "HS HE", "HS LE"),
  block = c(rep("MYC, proliferation", 5),
            rep("developmental state", 4),
            rep("lineage identity", 6)),
  stringsAsFactors = FALSE)
stopifnot(setequal(prog_meta$key, unique(tab$key)))

block_levels <- c("MYC, proliferation", "developmental state", "lineage identity")

dat <- tab |>
  dplyr::filter(contrast %in% contrast_levels) |>
  dplyr::inner_join(prog_meta, by = "key") |>
  dplyr::mutate(
    contrast = factor(contrast, levels = contrast_levels),
    block    = factor(block, levels = block_levels),
    # ggplot draws a discrete axis bottom-up, so reverse for top-to-bottom reading
    label    = factor(label, levels = rev(prog_meta$label)),
    sig      = p < 0.05)
stopifnot(nrow(dat) == 60)

# --- panel -------------------------------------------------------------------
sig_fill <- c(`TRUE` = "#D73027", `FALSE` = "grey85")

p <- ggplot2::ggplot(dat, ggplot2::aes(x = effect, y = label)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = effect,
                                     y = label, yend = label),
                        linewidth = 0.32, colour = "grey55") +
  ggplot2::geom_point(ggplot2::aes(fill = sig), shape = 21, size = 1.7,
                      colour = "grey25", stroke = 0.25) +
  ggplot2::facet_grid(block ~ contrast, scales = "free_y", space = "free_y") +
  ggplot2::scale_fill_manual(values = sig_fill, name = NULL,
                             breaks = "TRUE", labels = "p < 0.05") +
  ggplot2::labs(x = "effect (within-group SD)", y = NULL) +
  theme_panel() +
  ggplot2::theme(legend.position = "bottom",
                 legend.margin = ggplot2::margin(-3, 0, 0, 0),
                 axis.line.y   = ggplot2::element_blank(),
                 axis.ticks.y  = ggplot2::element_blank(),
                 strip.text.y  = ggplot2::element_text(angle = -90, size = 6),
                 panel.spacing.x = ggplot2::unit(2.0, "mm"),
                 panel.spacing.y = ggplot2::unit(1.6, "mm"))

# --- the legend text (never drawn) -------------------------------------------
e <- function(k, cn) tab$effect[tab$key == k & tab$contrast == cn]
q <- function(k, cn) tab$p[tab$key == k & tab$contrast == cn]
f <- function(x) sprintf("%+.2f", x)

LEGEND <- panel_legend(
  slot = "Fig. 1C",
  what = paste0(
    "Effect of each of the four contrasts of Fig. S1B on the canonical MYC and ",
    "proliferation signatures, on the developmental state axes, and on the Gray ",
    "lineage-identity axes, in units of each programme's own within-group ",
    "standard deviation."),
  detail = c(
    "n = 6 per group, n = 24. Points filled where p < 0.05. Each effect is an ordinary least-squares contrast on the per-sample GSVA composite, fitted on the relevant subset: myc_6W and myc_12W are the genotype difference within one age (6 versus 6), 6>12W_wt and 6>12W_myc the age difference within one genotype (6 versus 6).",
    "Composites are means of GSVA scores over member sets: MYC signatures 17 sets, Felsher 1, Hallmark MYC V2 1, MYC-in-TEB 4, proliferation 14, TEB-ductal = mean(UP) - mean(DN) over 3 + 3, the three MEC states 13 to 29 sets each from the curated annotation.",
    sprintf("Lineage identity, from Gray et al. 2023: for each major type (AP alveolar progenitor, BA basal, HS hormone-sensing) the composite of the regulator programmes detected in the HIGH-expressing (HE) and the LOW-expressing (LE) cells of that type - AP %d and %d sets, BA %d and %d, HS %d and %d. HE is the differentiated end of the axis, LE the lineage-suppressed end.",
            length(lanes[["AP HE"]]), length(lanes[["AP LE"]]),
            length(lanes[["BA HE"]]), length(lanes[["BA LE"]]),
            length(lanes[["HS HE"]]), length(lanes[["HS LE"]])),
    sprintf("MYC AXIS - the genotype effect is large and essentially unchanged with age: MYC signatures %s at 6 weeks (p = %.3f) and %s at 12 (p = %.4f); Hallmark MYC V2 %s and %s; MYC-in-TEB %s (p = %.3f) and %s (p = %.3f). Proliferation is smaller and does not reach significance at either age (%s, p = %.2f; %s, p = %.2f).",
            f(e("myc","myc_6W")), q("myc","myc_6W"),
            f(e("myc","myc_12W")), q("myc","myc_12W"),
            f(e("hallmark_v2","myc_6W")), f(e("hallmark_v2","myc_12W")),
            f(e("myc_in_teb","myc_6W")), q("myc_in_teb","myc_6W"),
            f(e("myc_in_teb","myc_12W")), q("myc_in_teb","myc_12W"),
            f(e("prolif","myc_6W")), q("prolif","myc_6W"),
            f(e("prolif","myc_12W")), q("prolif","myc_12W")),
    sprintf("TEB AXIS - the Myc effect is a 6-week effect: %s at 6 weeks (p = %.2f) against %s at 12 (p = %.2f). Both genotypes then travel down the axis, the Myc+ gland significantly so (%s, p = %.4f) and the wild type in the same direction (%s, p = %.2f).",
            f(e("teb_ductal","myc_6W")), q("teb_ductal","myc_6W"),
            f(e("teb_ductal","myc_12W")), q("teb_ductal","myc_12W"),
            f(e("teb_ductal","6>12W_myc")), q("teb_ductal","6>12W_myc"),
            f(e("teb_ductal","6>12W_wt")), q("teb_ductal","6>12W_wt")),
    sprintf("MEC STATES - BMYO is suppressed by Myc at both ages and more deeply at 12 weeks (%s, p = %.2f; %s, p = %.3f). LHS reverses: %s at 6 weeks, %s at 12 (interaction p = %.2f). LASP does the same, weakly (%s then %s).",
            f(e("BMYO","myc_6W")), q("BMYO","myc_6W"),
            f(e("BMYO","myc_12W")), q("BMYO","myc_12W"),
            f(e("LHS","myc_6W")), f(e("LHS","myc_12W")),
            ss$int_p[as.character(ss$state) == "LHS"],
            f(e("LASP","myc_6W")), f(e("LASP","myc_12W"))),
    sprintf("LINEAGE IDENTITY - Myc lowers the differentiated (HE) composite and raises the lineage-suppressed (LE) composite in all three lineages, at 6 weeks: AP %s / %s (p = %.2f / %.3f), BA %s / %s (p = %.2f / %.3f), HS %s / %s (p = %.2f / %.3f). At 12 weeks the same pattern is present but weaker (AP %s / %s, BA %s / %s, HS %s / %s).",
            f(e("AP HE","myc_6W")), f(e("AP LE","myc_6W")),
            q("AP HE","myc_6W"), q("AP LE","myc_6W"),
            f(e("BA HE","myc_6W")), f(e("BA LE","myc_6W")),
            q("BA HE","myc_6W"), q("BA LE","myc_6W"),
            f(e("HS HE","myc_6W")), f(e("HS LE","myc_6W")),
            q("HS HE","myc_6W"), q("HS LE","myc_6W"),
            f(e("AP HE","myc_12W")), f(e("AP LE","myc_12W")),
            f(e("BA HE","myc_12W")), f(e("BA LE","myc_12W")),
            f(e("HS HE","myc_12W")), f(e("HS LE","myc_12W"))),
    "The interaction is not drawn (it is the difference between the two genotype columns, and equivalently between the two development columns). No programme here has a significant interaction: p = 0.26 to 0.93 for the six signatures and states of script 27, and 0.26 to 0.86 for the six lineage axes."),
  bounds = c(
    "The two contrast families are not equally powered and must not be reported as if they were. Each genotype contrast is 6 versus 6 balanced within batch and is clean. Each development contrast is 6 versus 6 confounded with batch (batch = timepoint), and none of them is significant except the Myc+ TEB trajectory - they are a consistent direction across programmes, described and not claimed.",
    "THIS PANEL IS NOT THE ATTENUATION RESULT. On the MYC programme the genotype effect is the same size at both ages, so in cohort-relative GSVA space there is no attenuation to see; the attenuation result is a DESeq2 effect-size result (the Myc effect rescales by about 0.55 between the two ages, scripts 29-31 and 40) and is measured on a different ruler. What attenuates in this panel is the developmental arm, not the MYC arm.",
    "GSVA is cohort-relative, so an effect is a difference in position within these 24 samples. It cannot be read as an absolute change in programme activity, and the four contrasts are internally comparable but not comparable with the DESeq2 log fold changes.",
    "The LE composites sit almost on top of the global common-mode axis that every per-sample composite in this dataset shares (correlation with the mean of all 885 scores: AP LE 0.88, BA LE 0.90, HS LE 0.85), so on their own they are not independent evidence - see scripts 36 and 37 and the interpretability trap recorded there. The HE composites do NOT (AP -0.37, BA -0.57, HS +0.07), and HE and LE move in OPPOSITE directions within a lineage, which a single common-mode axis cannot produce. The de-differentiation reading rests on that opposition and on the HE half.",
    "The HE and LE lanes are not symmetric objects: the LE lanes are larger (165 to 788 genes against 26 to 159) and carry more mitochondrial genes (median 5 to 9 per cent against 0.7 to 3), both of which make an LE composite track the global axis more closely. The HE-minus-LE difference is therefore reported as two arms rather than as one difference score.",
    "These are regulator programmes detected in a cell state, not the state's own expression signature, so a change is evidence about the regulatory programme and only indirectly about lineage identity. Gray's HE/LE are human breast atlas contrasts mapped to mouse orthologs.",
    "The signatures are heavily inter-correlated (MYC signatures with Felsher r = 0.99, with Hallmark V2 r = 0.98, with proliferation r = 0.77), so the rows are not independent observations and the composite p-values are indicative.",
    "Endogenous Myc establishing the pubertal TEB phenotype is a literature claim; what these data show is wild-type co-variation consistent with it, plus that the transgene moves the same axis at 6 weeks. Not causal."),
  source = c(
    "results/myc_endogenous_amplification.rds (scripts/27) -- $prog_stats, $prog_wide",
    "results/dev_program_myc_integration.rds (scripts/26, PART C2) -- $annot, $state_stats",
    "results/gsva_scores.rds (scripts/15) -- per-sample scores; Gray lineage lanes built here as the mean over TFT_<TF>_GRAY_<lineage>_<HE|LE>, excluding the _MITO promotions",
    "Gray lane construction: docs/library_reference/Gray_et_al_developmental_TFS_selection.md sections 2a and 2c"))

# FULL WIDTH, and the reason is the shared x axis. All four columns share one
# scale in within-group SD, because the comparison the sentence makes is between
# a genotype effect (up to +3.0 SD) and a development effect (down to -1.7 SD) --
# free scales would hide exactly that. Spanning 4.7 SD legibly needs a wide
# column: at 120 mm a -0.6 SD point sits 3 mm off the zero line and reads as
# nothing, which is what made the previous version unreadable.
save_panel_p(p, "fig1C_myc_teb_proliferation",
             width = fig_w[["double"]], height = 92)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the whole table behind the panel, one row per programme
  tab |>
    dplyr::select(key, contrast, effect, p) |>
    tidyr::pivot_wider(names_from = contrast, values_from = c(effect, p)) |>
    as.data.frame() |> print()

  ## how close each composite sits to the global common-mode axis
  gm <- colMeans(scores)
  vapply(Y, function(y) stats::cor(y, gm), numeric(1)) |> round(3) |> print()

  ## the Gray lanes that went into block 3
  lanes |> lapply(length) |> unlist() |> print()
  lanes[["AP LE"]] |> print()
}
