# =============================================================================
# fig1B_myc_teb_proliferation.R -- what Myc does at each age, and what the
# gland does over time
# -----------------------------------------------------------------------------
# WAS Fig. 1C, IS NOW Fig. 1B (author, 2026-07-30). The cell-state panel that
# held the 1B slot plotted the same GSVA z-scores as this one, one lens further
# back, and was dropped as redundant; its 2x2 group-mean view is kept in this
# script's sandbox block so nothing is lost.
#
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
# 0. FORM. A Cleveland dot plot, single column (89 mm), with the effect encoded
#    TWICE: as the dot's position and as its fill on the manuscript diverging
#    scale (deep espresso brown / stark white / crisp mint green, declared once in
#    _panel_common.R as ms_diverging). Significance is a small asterisk rather
#    than a second fill colour, because at n=6 per cell the effect size is the
#    interesting quantity and the p-value is the footnote.
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
# 1b. AND THE BLOCKS ARE MERGED. Developmental state and lineage identity were
#    two blocks asking one question, so they are one block now, with each MEC
#    state followed by the LE composite of its own lineage (BMYO with BA, LASP
#    with AP, LHS with HS) labelled "-LE". Only the LE arm is drawn: HE and LE are
#    near-mirror images, so the HE rows spent six lines saying one thing. Their
#    numbers stay in the legend block, because the OPPOSITION of the two signs is
#    the argument (see bounds).
#
# 2. THE TWO TIMELINES ARE STILL SIDE BY SIDE, AND THEY REALLY DO LOOK ALIKE on
#    the Myc signatures (-0.60 SD in wild type, -0.68 SD in Myc+). That is not a
#    drawing problem, it is the result: in cohort-relative GSVA space the genotype
#    gap on the Myc programme is the same at both ages (2.36 vs 2.28 SD), i.e.
#    THERE IS NO ATTENUATION HERE. The attenuation result lives on a different
#    ruler -- DESeq2 effect sizes and their fGSEA enrichment, where the Myc effect
#    rescales x0.55 (scripts 29-31, 40) -- and this panel must not be read as
#    evidence for or against it. What DOES attenuate here is the developmental
#    arm: TEB +1.02 -> +0.23, and every lineage-identity axis weakens with it.
#    That contrast (Myc programme constant, developmental identity fading) is the
#    honest reading and it is the one that sets up the attenuation section.
#
# 3. THE LINEAGE-IDENTITY BLOCK: the Gray 2023 axes. For each major type
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
# 4. THE HUMAN ANCHOR (author, 2026-07-30): the METABRIC MB2 fork joins the top
#    block, which is renamed "Myc-tumourigenesis" because it is no longer only
#    signatures -- it now runs from the Myc programme through proliferation to a
#    human breast-cancer state. Drawn as "Human BRCA-MYC", the one label on the
#    panel that keeps the human MYC spelling; everything else is mouse Myc.
#    It does not yet appear in the narrative.
#
# Input:  results/myc_endogenous_amplification.rds  (script 27 -- $prog_stats, $prog_wide)
#         results/dev_program_myc_integration.rds   (script 26 -- $annot, $state_stats)
#         results/gsva_scores.rds                   (script 15 -- per-sample scores)
# Output: outputs/figures/panels/fig1B_myc_teb_proliferation.pdf
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

# --- the human anchor: the METABRIC MB2 fork ---------------------------------
# MB2 is the Myc arm of the MCbiclust multistate switch found in human breast
# cancer (Menegollo, Bentham et al., Cancer Res 2024, the analytical companion
# paper): biogenesis-high, proliferative, luminal-progenitor, ER-negative, WITH
# Myc activation -- against MB1, the same biogenesis-high upper fork WITHOUT Myc.
#
# WHICH SCORE: THE RAW UPPER FORK (author, 2026-07-30). MB2_UF is
# METABRIC_MB2_HI_CV_GROUP1, scored directly. The alternative is the switch
# POSITION, UF minus the lower fork GROUP2, which is script 18's convention for
# these anticorrelated poles (18:96-101) -- it is a weaker number (+1.46 p = 0.053
# at 6 weeks, +1.09 p = 0.045 at 12, against the upper fork's +1.63 p = 0.039 and
# +1.59 p = 0.005) and is reported in the legend block. Drawing the upper fork
# alone is the right call HERE because this row is a resemblance statement -- how
# far the mouse tissue sits toward a human MYC-driven state -- not a claim about
# where a switch is thrown. Fig. S1C carries the switch, and the specificity test
# against MB1 that goes with it.
mb2_uf <- "METABRIC_MB2_HI_CV_GROUP1"
mb2_lf <- "METABRIC_MB2_HI_CV_GROUP2"
stopifnot(all(c(mb2_uf, mb2_lf) %in% rownames(scores)))
mb2_uf_score <- scores[mb2_uf, ]
mb2_fork     <- scores[mb2_uf, ] - scores[mb2_lf, ]   # for the assertion + bounds

Y <- cbind(pw[, c("myc", "felsher", "hallmark_v2", "myc_in_teb", "prolif",
                  "teb_ductal")],
           state_mat, lane_mat, mb2_uf = mb2_uf_score)

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
# 18's MB2 fork score, per sample. results/ap7_mb_fork.rds is a Jul-6 object and
# so predates the gene-symbol reconciliation, which is why it is NOT read through
# require_fresher_than() -- but the MB2 HI_CV sets were untouched by it and the
# two agree exactly, so it is still a valid guard on the fork definition.
ap7_path <- here::here("results", "ap7_mb_fork.rds")
if (file.exists(ap7_path)) {
  ap7 <- as.data.frame(readRDS(ap7_path)$fork_df)
  rownames(ap7) <- ap7$sample
  stopifnot(max(abs(mb2_fork[rownames(ap7)] - ap7$MB2_score)) < 1e-8)
}

# --- labels and blocks -------------------------------------------------------
# Two blocks, not three: the developmental-state axes and the lineage-identity
# axes are the same question asked twice, so they are merged under "lineage
# identity" with each state followed by the LE composite of its own lineage --
# BMYO with BA, LASP with AP, LHS with HS. The paired row is labelled "-LE" so it
# reads as a subordinate of the state above it.
#
# ONLY THE LE ARM IS DRAWN. HE and LE are near-mirror images (they anticorrelate
# at -0.71 / -0.76 / -0.22 within a lineage), so drawing both cost six rows to say
# one thing. The HE numbers stay in the legend block, because the OPPOSITION of
# the two signs is what licenses the de-differentiation reading -- see bounds.
# MOUSE NOMENCLATURE (author, 2026-07-30): the gene and every programme scored on
# these mouse samples is "Myc", not "MYC". The one exception on the panel is the
# human anchor, which keeps the human convention -- "Human BRCA-MYC".
prog_meta <- data.frame(
  key = c("myc", "felsher", "hallmark_v2", "myc_in_teb", "prolif", "mb2_uf",
          "teb_ductal", "BMYO", "BA LE", "LASP", "AP LE", "LHS", "HS LE"),
  label = c("Myc signatures", "Myc-Felsher", "Hallmark Myc V2", "Myc-in-TEB",
            "Proliferation", "Human BRCA-MYC",
            "TEB - ductal", "BMYO", "-LE", "LASP", "-LE", "LHS", "-LE"),
  block = c(rep("Myc-tumourigenesis", 6),
            rep("lineage identity", 7)),
  stringsAsFactors = FALSE)
stopifnot(all(prog_meta$key %in% unique(tab$key)))

block_levels <- c("Myc-tumourigenesis", "lineage identity")

dat <- tab |>
  dplyr::filter(contrast %in% contrast_levels, key %in% prog_meta$key) |>
  dplyr::inner_join(prog_meta, by = "key") |>
  dplyr::mutate(
    contrast = factor(contrast, levels = contrast_levels),
    block    = factor(block, levels = block_levels),
    # the KEY carries the axis order (three rows share the label "-LE"); ggplot
    # draws a discrete axis bottom-up, hence rev()
    key      = factor(key, levels = rev(prog_meta$key)),
    sig      = p < 0.05)
stopifnot(nrow(dat) == 52)

# --- panel -------------------------------------------------------------------
# The fill is the SAME quantity as the x position, on the manuscript diverging
# scale. That redundancy is the point: at single-column width the four columns get
# about 15 mm each, so a -0.8 SD dot sits 2.5 mm off the zero line and its size
# does not read from position alone.
#
# WHITE IS AT ZERO AND THE RANGE IS ASYMMETRIC (-1.7 to +3.0), because these data
# are: a symmetric +/-3.1 ramp put every negative effect in the first third of the
# brown, where BMYO at -1.16 was indistinguishable from LASP at -0.20 and the fill
# had stopped doing its job on that side. The colour bar is therefore KEPT, so the
# asymmetry is inspectable rather than implied; equal ink does not mean equal
# magnitude across the sign change, and the x axis is what carries magnitude.
RNG <- c(floor(min(dat$effect) * 10) / 10, ceiling(max(dat$effect) * 10) / 10)

y_labels <- stats::setNames(prog_meta$label, prog_meta$key)

p <- ggplot2::ggplot(dat, ggplot2::aes(x = effect, y = key)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey80") +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = effect,
                                     y = key, yend = key),
                        linewidth = 0.25, colour = "grey60") +
  ggplot2::geom_point(ggplot2::aes(fill = effect), shape = 21, size = 1.9,
                      colour = "grey30", stroke = 0.22) +
  # significance is a small asterisk above the dot, not a second colour
  ggplot2::geom_text(data = dat[dat$sig, , drop = FALSE],
                     ggplot2::aes(x = effect, y = key), label = "*",
                     size = 2.0, colour = "grey15", vjust = -0.45) +
  ggplot2::facet_grid(block ~ contrast, scales = "free_y", space = "free_y") +
  heat_fill(RNG, breaks = c(RNG[1], 0, RNG[2])) +
  # The default 5 per cent expansion is not enough at this width: the dot at the
  # extreme of the range (Hallmark Myc V2 at myc_6W, TEB-ductal at 6>12W_myc) is
  # drawn half outside the panel and clipped. The expansion has to clear the dot
  # RADIUS plus the asterisk, so it is set explicitly rather than left to default.
  # check.overlap drops any tick label that would still collide in a ~15 mm column.
  ggplot2::scale_x_continuous(
    breaks = seq(-1, 3, by = 1),
    expand = ggplot2::expansion(mult = 0.11),
    guide = ggplot2::guide_axis(check.overlap = TRUE)) +
  # and the same for the top and bottom rows, where the asterisk sits above the dot
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::scale_y_discrete(labels = y_labels) +
  ggplot2::labs(x = "effect (within-group SD)", y = NULL) +
  ggplot2::guides(fill = ggplot2::guide_colourbar(
    direction = "horizontal", title = NULL,
    barwidth = ggplot2::unit(26, "mm"), barheight = ggplot2::unit(1.5, "mm"),
    ticks.colour = "grey30", frame.colour = "grey30",
    frame.linewidth = 0.15)) +
  theme_panel(base_size = 6) +
  ggplot2::theme(axis.line.y   = ggplot2::element_blank(),
                 axis.ticks.y  = ggplot2::element_blank(),
                 strip.text.x  = ggplot2::element_text(face = "bold", size = 5.8),
                 strip.text.y  = ggplot2::element_text(angle = -90, size = 5.4,
                                                       face = "plain"),
                 panel.spacing.x = ggplot2::unit(1.2, "mm"),
                 panel.spacing.y = ggplot2::unit(1.4, "mm"),
                 legend.position = "bottom",
                 legend.margin = ggplot2::margin(-2, 0, 0, 0),
                 plot.margin   = ggplot2::margin(1.5, 1.5, 1, 1, "mm"))

# --- the legend text (never drawn) -------------------------------------------
e <- function(k, cn) tab$effect[tab$key == k & tab$contrast == cn]
q <- function(k, cn) tab$p[tab$key == k & tab$contrast == cn]
f <- function(x) sprintf("%+.2f", x)

# the alternative score: the switch POSITION (upper fork minus lower), quoted in
# the bounds so the choice of the raw upper fork is visible rather than silent
fork_alt <- contrast_table(mb2_fork, sample_meta$timepoint,
                           sample_meta$myc_status, sample_meta$group)

LEGEND <- panel_legend(
  slot = "Fig. 1B",
  what = paste0(
    "Effect of each of the four contrasts of Fig. S1B on the canonical Myc and ",
    "proliferation signatures and on the lineage-identity axes, in units of each ",
    "programme's own within-group standard deviation. Dot position and dot fill ",
    "are the same quantity."),
  detail = c(
    "n = 6 per group, n = 24. Each effect is an ordinary least-squares contrast on the per-sample GSVA composite, fitted on the relevant subset: myc_6W and myc_12W are the genotype difference within one age (6 versus 6), 6>12W_wt and 6>12W_myc the age difference within one genotype (6 versus 6). An asterisk marks p < 0.05, uncorrected.",
    "Composites are means of GSVA scores over member sets: Myc signatures 17 sets, Myc-Felsher 1, Hallmark Myc V2 1, Myc-in-TEB 4, proliferation 14, TEB-ductal = mean(UP) - mean(DN) over 3 + 3, the three MEC states 13 to 29 sets each from the curated annotation.",
    sprintf("Human BRCA-MYC is the upper fork of the METABRIC MB2 switch (METABRIC_MB2_HI_CV_GROUP1, %d mouse-mapped genes): one arm of an MCbiclust multistate switch that stratifies human breast cancer, biogenesis-high, proliferative, luminal-progenitor, ER-negative and MYC-activated, against a lower fork that is glycolytic and stem-like (Menegollo, Bentham et al., Cancer Res 2024). Myc+ tissue resembles it at both ages, %s at 6 weeks (p = %.3f) and %s at 12 (p = %.3f).",
            length(gsva$pathways[[mb2_uf]]),
            f(e("mb2_uf","myc_6W")), q("mb2_uf","myc_6W"),
            f(e("mb2_uf","myc_12W")), q("mb2_uf","myc_12W")),
    sprintf("Lineage identity pairs each mammary epithelial state with the LE composite of its own lineage - BMYO with basal (BA), LASP with alveolar progenitor (AP), LHS with hormone-sensing (HS) - and the paired row is labelled '-LE'. From Gray et al. 2023: for each major type the regulator programmes detected in the LOW-expressing cells of that type, i.e. the lineage-suppressed end of the axis (AP %d sets, BA %d, HS %d). The matching HIGH-expressing (HE, differentiated) composites are not drawn (AP %d sets, BA %d, HS %d) but are reported below, because the two arms move in opposite directions and that is the argument.",
            length(lanes[["AP LE"]]), length(lanes[["BA LE"]]),
            length(lanes[["HS LE"]]), length(lanes[["AP HE"]]),
            length(lanes[["BA HE"]]), length(lanes[["HS HE"]])),
    sprintf("MYC AXIS - the genotype effect is large and essentially unchanged with age: Myc signatures %s at 6 weeks (p = %.3f) and %s at 12 (p = %.4f); Hallmark Myc V2 %s and %s; Myc-in-TEB %s (p = %.3f) and %s (p = %.3f). Proliferation is smaller and does not reach significance at either age (%s, p = %.2f; %s, p = %.2f).",
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
    sprintf("MAMMARY EPITHELIAL STATES - BMYO is suppressed by Myc at both ages and more deeply at 12 weeks (%s, p = %.2f; %s, p = %.3f). LHS reverses: %s at 6 weeks, %s at 12 (interaction p = %.2f). LASP does the same, weakly (%s then %s).",
            f(e("BMYO","myc_6W")), q("BMYO","myc_6W"),
            f(e("BMYO","myc_12W")), q("BMYO","myc_12W"),
            f(e("LHS","myc_6W")), f(e("LHS","myc_12W")),
            ss$int_p[as.character(ss$state) == "LHS"],
            f(e("LASP","myc_6W")), f(e("LASP","myc_12W"))),
    sprintf("LINEAGE IDENTITY - Myc raises every lineage-suppressed (LE) composite and lowers every differentiated (HE) composite, in all three lineages, at 6 weeks; drawn arm first, undrawn arm second: AP %s LE against %s HE (p = %.3f / %.2f), BA %s against %s (p = %.3f / %.2f), HS %s against %s (p = %.3f / %.2f). At 12 weeks the same pattern is present but weaker (AP %s / %s, BA %s / %s, HS %s / %s).",
            f(e("AP LE","myc_6W")), f(e("AP HE","myc_6W")),
            q("AP LE","myc_6W"), q("AP HE","myc_6W"),
            f(e("BA LE","myc_6W")), f(e("BA HE","myc_6W")),
            q("BA LE","myc_6W"), q("BA HE","myc_6W"),
            f(e("HS LE","myc_6W")), f(e("HS HE","myc_6W")),
            q("HS LE","myc_6W"), q("HS HE","myc_6W"),
            f(e("AP LE","myc_12W")), f(e("AP HE","myc_12W")),
            f(e("BA LE","myc_12W")), f(e("BA HE","myc_12W")),
            f(e("HS LE","myc_12W")), f(e("HS HE","myc_12W"))),
    "The interaction is not drawn (it is the difference between the two genotype columns, and equivalently between the two development columns). No programme here has a significant interaction: p = 0.26 to 0.93 for the six signatures and states of script 27, and 0.26 to 0.86 for the six lineage axes."),
  bounds = c(
    "The two contrast families are not equally powered and must not be reported as if they were. Each genotype contrast is 6 versus 6 balanced within batch and is clean. Each development contrast is 6 versus 6 confounded with batch (batch = timepoint), and none of them is significant except the Myc+ TEB trajectory - they are a consistent direction across programmes, described and not claimed.",
    "THIS PANEL IS NOT THE ATTENUATION RESULT. On the Myc programme the genotype effect is the same size at both ages, so in cohort-relative GSVA space there is no attenuation to see; the attenuation result is a DESeq2 effect-size result (the Myc effect rescales by about 0.55 between the two ages, scripts 29-31 and 40) and is measured on a different ruler. What attenuates in this panel is the developmental arm, not the Myc arm.",
    "GSVA is cohort-relative, so an effect is a difference in position within these 24 samples. It cannot be read as an absolute change in programme activity, and the four contrasts are internally comparable but not comparable with the DESeq2 log fold changes.",
    sprintf("The Human BRCA-MYC row is a RESEMBLANCE, not a claim that this tissue is that tumour: it says where the mouse samples sit on a switch defined in human data, mapped to mouse orthologs. It is a single GSVA composite over %d mouse-mapped genes (within-group SD %.2f). The upper fork scored on its own gives a stronger number (%s at 6 weeks, p = %.3f; %s at 12, p = %.4f), but the project convention for these anticorrelated poles is the contrast (scripts/18) and that is what is drawn. Its correlation with the global common-mode axis is 0.94, so the same bound as the LE rows applies, harder. And this row does NOT by itself separate 'Myc drives the Myc fork' from 'Myc drives mitochondrial biogenesis generically': the non-Myc upper fork MB1 rises just as much (+1.73 at 6 weeks, +1.34 at 12), and the two scores correlate 0.98. Fig. S1C carries that specificity test.",
            length(gsva$pathways[[mb2_uf]]),
            unique(tab$within_sd[tab$key == "mb2_uf"]),
            f(fork_alt$effect[fork_alt$contrast == "myc_6W"]),
            fork_alt$p[fork_alt$contrast == "myc_6W"],
            f(fork_alt$effect[fork_alt$contrast == "myc_12W"]),
            fork_alt$p[fork_alt$contrast == "myc_12W"]),
    "THE DRAWN LE ARM IS THE ONE MOST EXPOSED TO THE COMMON-MODE AXIS. Every per-sample composite in this dataset shares one dominant axis, and the LE composites sit almost on top of it (correlation with the mean of all 885 scores: AP LE 0.88, BA LE 0.90, HS LE 0.85), so an LE rise is not on its own independent evidence - see scripts 36 and 37 and the interpretability trap recorded there. The undrawn HE composites do NOT sit on it (AP -0.37, BA -0.57, HS +0.07), and HE and LE move in OPPOSITE directions within a lineage (they anticorrelate at -0.71, -0.76 and -0.22), which a single common-mode axis cannot produce. The de-differentiation reading rests on that opposition, so the HE arm has to be reported in the text even though it is not on the panel.",
    "The HE and LE lanes are not symmetric objects: the LE lanes are larger (165 to 788 genes against 26 to 159) and carry more mitochondrial genes (median 5 to 9 per cent against 0.7 to 3), both of which make an LE composite track the global axis more closely. They are therefore reported as two arms and never combined into one HE-minus-LE difference score.",
    "These are regulator programmes detected in a cell state, not the state's own expression signature, so a change is evidence about the regulatory programme and only indirectly about lineage identity. Gray's HE/LE are human breast atlas contrasts mapped to mouse orthologs.",
    "The signatures are heavily inter-correlated (Myc signatures with Myc-Felsher r = 0.99, with Hallmark V2 r = 0.98, with proliferation r = 0.77), so the rows are not independent observations and the composite p-values are indicative. The asterisks are uncorrected and there are 52 of them on the panel; they mark which effects clear a nominal threshold, not which survive multiplicity.",
    "Endogenous Myc establishing the pubertal TEB phenotype is a literature claim; what these data show is wild-type co-variation consistent with it, plus that the transgene moves the same axis at 6 weeks. Not causal."),
  source = c(
    "results/myc_endogenous_amplification.rds (scripts/27) -- $prog_stats, $prog_wide",
    "results/dev_program_myc_integration.rds (scripts/26, PART C2) -- $annot, $state_stats",
    "results/gsva_scores.rds (scripts/15) -- per-sample scores; Gray lineage lanes built here as the mean over TFT_<TF>_GRAY_<lineage>_<HE|LE>, excluding the _MITO promotions",
    "Gray lane construction: docs/library_reference/Gray_et_al_developmental_TFS_selection.md sections 2a and 2c"))

# SINGLE COLUMN (89 mm). All four columns share one x scale in within-group SD,
# because the comparison the sentence makes is between a genotype effect (up to
# +3.0 SD) and a development effect (down to -1.7 SD), and free scales would hide
# exactly that. At this width 4.7 SD spans about 15 mm per column, so a -0.6 SD
# dot sits only 2 mm off the zero line -- which is why the fill carries the
# magnitude as well. Type is set for the final size, not the preview.
save_panel_p(p, "fig1B_myc_teb_proliferation",
             width = fig_w[["single"]], height = 84)

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

  ## the six undrawn HE / LE rows, which carry the opposite-signs argument
  tab |>
    dplyr::filter(key %in% names(lanes), contrast %in% contrast_geno) |>
    dplyr::select(key, contrast, effect, p) |>
    tidyr::pivot_wider(names_from = contrast, values_from = c(effect, p)) |>
    as.data.frame() |> print()

  ## how close each composite sits to the global common-mode axis. The LE arms are
  ## on it (0.85-0.90); the HE arms are not. This is the bound that matters.
  gm <- colMeans(scores)
  vapply(Y, function(y) stats::cor(y, gm), numeric(1)) |> round(3) |> print()

  ## the Gray lanes, and the mitochondrial content that makes LE and HE asymmetric
  lanes |> lapply(length) |> unlist() |> print()
  mito_u <- unique(unlist(gsva$pathways[grep("^MITOCARTA_", names(gsva$pathways))]))
  vapply(lanes, function(ss)
    stats::median(vapply(ss, function(s) mean(gsva$pathways[[s]] %in% mito_u),
                         numeric(1))), numeric(1)) |> round(3) |> print()

  ## THE PANEL THAT USED TO HOLD THE 1B SLOT, kept because a contrast plot cannot
  ## show a LEVEL: the four groups as the 2x2 design, per programme, group-mean
  ## z-score in within-group SD. Dropped from the figure as redundant with the
  ## contrasts above, not because it was wrong.
  {
    lv <- c("BMYO", "LASP", "LHS", "TEB - ductal")
    zz <- data.frame(
      programme = rep(lv, each = 24),
      timepoint = rep(sample_meta$timepoint, times = 4),
      myc_status = rep(sample_meta$myc_status, times = 4),
      group = rep(sample_meta$group, times = 4),
      value = c(state_mat[, "BMYO"], state_mat[, "LASP"], state_mat[, "LHS"],
                pw$teb_ductal))
    zz <- zz |>
      dplyr::group_by(programme) |>
      dplyr::mutate(z = (value - mean(value)) / wsd_of(value, group)) |>
      dplyr::ungroup() |>
      dplyr::group_by(programme, timepoint, myc_status) |>
      dplyr::summarise(mean_z = mean(z), .groups = "drop") |>
      dplyr::mutate(programme = factor(programme, levels = lv),
                    age  = factor(as.character(timepoint), levels = c("6W", "12W"),
                                  labels = c("6 weeks", "12 weeks")),
                    geno = factor(as.character(myc_status), levels = c("pos", "neg"),
                                  labels = c("Myc+", "WT")))
    L <- ceiling(max(abs(zz$mean_z)) * 10) / 10
    ggplot2::ggplot(zz, ggplot2::aes(age, geno, fill = mean_z)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.2f", mean_z),
                                      colour = ink_on_fill(mean_z, L)),
                         size = 2.1, show.legend = FALSE) +
      ggplot2::facet_wrap(~ programme, nrow = 2) +
      heat_fill(L, name = "mean z") +
      ggplot2::scale_colour_identity() +
      ggplot2::labs(x = NULL, y = NULL) +
      theme_panel()
  }
}
