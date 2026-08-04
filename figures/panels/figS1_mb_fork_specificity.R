# =============================================================================
# figS1_mb_fork_specificity.R -- is the human resemblance MYC-specific, or just
# mitochondrial biogenesis?
# -----------------------------------------------------------------------------
# SLOT: NONE (2026-08-04). It held S1C, then S1D, and is now displaced out of the
# figure altogether: the written Fig. S1D is the MYC western blot, and this panel
# is still cited nowhere. It keeps building and keeps its legend block, and it
# gets its letter back the day a sentence asks for it -- renumbering it a fourth
# time before then would be three renames for nothing. figures/panels/PANELS.md
# is the slug -> slot map; anything that lays panels out by slot (panels_to_pdf.R)
# skips this one because the slot below is not a Figure 1 or S1 letter.
#
# SUPPORTS the Human BRCA-MYC row of Fig. 1B. That row says the Myc+ gland
# resembles the MB2 upper fork of the human breast-cancer switch. This panel is
# the control that says whether the resemblance is to the MYC arm SPECIFICALLY.
#
# THE DESIGN OF THE TEST (AP7, scripts/18, and the analytical companion paper
# Menegollo, Bentham et al., Cancer Res 2024). MCbiclust found that MB1 and MB2
# share a lower fork almost completely and then diverge into two distinct UPPER
# forks, both mitochondrial-biogenesis-high and proliferative:
#   MB1_UF -- mature-luminal, WITHOUT MYC
#   MB2_UF -- luminal-progenitor, ER-negative, WITH MYC
# The only feature distinguishing them is MYC. So MB2_UF minus MB1_UF isolates
# "Myc drives the Myc fork" from "Myc drives biogenesis generically", and it is
# the decisive contrast -- MB2_UF on its own cannot make that distinction.
#
# FORM (author, 2026-07-30): just the difference, as a small boxplot, because the
# difference IS the test. The earlier three-facet version also drew MB2_UF and
# MB1_UF separately to show that Myc raises both by almost the same amount; those
# numbers are in the legend block instead, and the two component scores stay in
# the sandbox at the end of this file.
#
# WHAT IT SAYS, AND IT IS NOT WHAT SCRIPT 18 PREDICTED. Split by age, the
# specificity contrast is absent at 6 weeks (+0.39 SD, p = 0.59) and present at 12
# (+1.62 SD, p = 0.003). Script 18's header predicted the opposite ("biogenesis is
# front-loaded at 6W so MB2_UF resemblance is predicted to peak at 6W"), and its
# headline was the POOLED genotype main effect (p = 0.026), which averages the two.
# The resemblance to MB2_UF itself is flat with age (+1.63 then +1.59), so nothing
# about this axis peaks at 6 weeks. Flagged, not resolved: the interaction is not
# significant either (p = 0.15) at n = 6 per cell.
#
# Input:  results/gsva_scores.rds    (script 15 -- per-sample scores)
#         results/ap7_mb_fork.rds    (script 18 -- assertion only; see below)
# Output: outputs/figures/panels/figS1_mb_fork_specificity.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

gsva        <- readRDS(here::here("results", "gsva_scores.rds"))
scores      <- gsva$scores
sample_meta <- as.data.frame(gsva$sample_meta)[colnames(scores), , drop = FALSE]

mb <- c(MB2_UF = "METABRIC_MB2_HI_CV_GROUP1",
        MB1_UF = "METABRIC_MB1_HI_CV_GROUP1")
stopifnot(all(mb %in% rownames(scores)))

mb2  <- scores[mb[["MB2_UF"]], ]
mb1  <- scores[mb[["MB1_UF"]], ]
d_spec <- mb2 - mb1   # NOT `diff`: base::diff is used below

# results/ap7_mb_fork.rds is a Jul-6 object and predates the gene-symbol
# reconciliation, so it is deliberately NOT gated by require_fresher_than() -- but
# the MB HI_CV sets were untouched by that rebuild and the two agree exactly, so
# it is still a valid guard on the contrast's definition (scripts/18:105).
ap7_path <- here::here("results", "ap7_mb_fork.rds")
if (file.exists(ap7_path)) {
  ap7 <- as.data.frame(readRDS(ap7_path)$fork_df)
  rownames(ap7) <- ap7$sample
  stopifnot(max(abs(d_spec[rownames(ap7)] - ap7$MB2_over_MB1_UF)) < 1e-8)
}

# AGE-MAJOR ORDER, so the two boxes a reader compares are adjacent: the test is
# the genotype gap WITHIN an age, and the design's ages are two batches.
grp_levels <- c("6W_neg", "6W_pos", "12W_neg", "12W_pos")
dat <- data.frame(
  group = factor(as.character(sample_meta$group), levels = grp_levels),
  value = as.numeric(d_spec))
stopifnot(nrow(dat) == 24, !any(is.na(dat$group)))

ct <- function(y) contrast_table(y, sample_meta$timepoint, sample_meta$myc_status,
                                 sample_meta$group)
td <- ct(d_spec)

# --- the within-age genotype brackets ----------------------------------------
# One bracket per age, over the pair it compares, carrying an asterisk only where
# the contrast clears p < 0.05. A bracket with a bare asterisk is unambiguous in a
# way a floating asterisk over one box is not.
YR   <- range(dat$value)
PAD  <- diff(YR) * 0.12
brk <- data.frame(
  age  = c("6W", "12W"),
  x    = c(1, 3), xend = c(2, 4),
  y    = YR[2] + PAD * c(0.7, 0.7),
  p    = c(td$p[td$contrast == "myc_6W"], td$p[td$contrast == "myc_12W"]),
  stringsAsFactors = FALSE)
brk$lab <- ifelse(brk$p < 0.05, "*", "")

# Layers and encoding are the figure layer's established idiom for a per-group
# distribution (figures/fig01_mito_content.R:118-170): colour AND fill both mapped
# to group, box at alpha 0.28 with a grey outline, solid points over it, the x
# axis blanked because the group key does that job, and group_cols / group_labels
# carried by the colour scale so the legend IS the project's group key.
pts_layer <- if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
  ggbeeswarm::geom_quasirandom(width = 0.22, size = 1.1, alpha = 0.95)
} else {
  ggplot2::geom_jitter(width = 0.16, height = 0, size = 1.1, alpha = 0.95)
}

p <- ggplot2::ggplot(dat, ggplot2::aes(x = group, y = value,
                                       colour = group, fill = group)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey80") +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.28,
                        colour = "grey35", linewidth = 0.3) +
  pts_layer +
  ggplot2::geom_segment(data = brk,
                        ggplot2::aes(x = x, xend = xend, y = y, yend = y),
                        inherit.aes = FALSE, linewidth = 0.25, colour = "grey35") +
  ggplot2::geom_text(data = brk,
                     ggplot2::aes(x = (x + xend) / 2, y = y, label = lab),
                     inherit.aes = FALSE, size = 2.2, colour = "grey15",
                     vjust = -0.15) +
  ggplot2::scale_colour_manual(values = group_cols, breaks = names(group_cols),
                               labels = group_labels) +
  ggplot2::scale_fill_manual(values = group_cols, guide = "none") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.16))) +
  ggplot2::labs(x = NULL, y = "MB2_UF - MB1_UF  (GSVA)", colour = NULL) +
  ggplot2::guides(colour = ggplot2::guide_legend(
    nrow = 2, override.aes = list(size = 1.6, alpha = 1, shape = 16))) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    axis.text.x     = ggplot2::element_blank(),
    axis.ticks.x    = ggplot2::element_blank(),
    axis.line.x     = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.margin   = ggplot2::margin(-2, 0, 0, 0),
    legend.key.size = ggplot2::unit(2.6, "mm"),
    plot.margin     = ggplot2::margin(1.5, 1.5, 1, 1, "mm"))

# --- the legend text (never drawn) -------------------------------------------
t2 <- ct(mb2); t1 <- ct(mb1)
g  <- function(t, cn) t$effect[t$contrast == cn]
gp <- function(t, cn) t$p[t$contrast == cn]
f  <- function(x) sprintf("%+.2f", x)
gm <- colMeans(scores)

LEGEND <- panel_legend(
  slot = "not currently cited",
  what = paste0(
    "Specificity control for the human breast-cancer resemblance of Fig. 1B. ",
    "Difference between the per-sample GSVA scores for the two upper forks of the ",
    "MCbiclust multistate switch: MB2, the MYC arm, minus MB1, the same ",
    "mitochondrial-biogenesis-high proliferative state without MYC."),
  detail = c(
    "n = 6 animals per group, n = 24; each point is one animal. Boxes are median and quartiles, whiskers 1.5x the interquartile range. Groups are ordered by age, and the key rather than an x axis carries the group labels. Positive means closer to the MYC arm than to the non-MYC one.",
    sprintf("MB2_UF is METABRIC_MB2_HI_CV_GROUP1 (%d mouse-mapped genes) and MB1_UF is METABRIC_MB1_HI_CV_GROUP1 (%d); they share %d genes. Both are the highest-variance gene lists defining the upper fork of their switch (Menegollo, Bentham et al., Cancer Res 2024; scripts/18).",
            length(gsva$pathways[[mb[["MB2_UF"]]]]),
            length(gsva$pathways[[mb[["MB1_UF"]]]]),
            length(intersect(gsva$pathways[[mb[["MB2_UF"]]]],
                             gsva$pathways[[mb[["MB1_UF"]]]]))),
    sprintf("The bracket is the genotype contrast within one age and the asterisk marks p < 0.05: %s within-group SD at 6 weeks (p = %.2f) and %s at 12 (p = %.4f). The interaction is p = %.2f and neither genotype moves along the contrast with age (%s SD in wild type, p = %.2f; %s SD in Myc+, p = %.2f).",
            f(g(td,"myc_6W")), gp(td,"myc_6W"),
            f(g(td,"myc_12W")), gp(td,"myc_12W"), gp(td,"interaction"),
            f(g(td,"6>12W_wt")), gp(td,"6>12W_wt"),
            f(g(td,"6>12W_myc")), gp(td,"6>12W_myc")),
    sprintf("The two component scores are not drawn, and they are the reason this contrast is the test rather than a detail: Myc raises BOTH forks, and by almost the same amount. In within-group SD, MB2_UF %s at 6 weeks (p = %.3f) and %s at 12 (p = %.4f); MB1_UF %s (p = %.3f) and %s (p = %.3f).",
            f(g(t2,"myc_6W")), gp(t2,"myc_6W"), f(g(t2,"myc_12W")), gp(t2,"myc_12W"),
            f(g(t1,"myc_6W")), gp(t1,"myc_6W"), f(g(t1,"myc_12W")), gp(t1,"myc_12W"))),
  bounds = c(
    sprintf("The two scores correlate %.3f across the 24 samples, so their difference is a small residual: within-group SD %.3f against %.3f and %.3f, about a quarter. Almost everything the two forks have in common cancels, which is the design of the test but also means the contrast is estimated on a fraction of the signal, and the y axis here spans a fraction of the range either score does.",
            stats::cor(mb2, mb1), unique(td$within_sd), unique(t2$within_sd),
            unique(t1$within_sd)),
    sprintf("That cancellation is what makes the contrast worth having. MB2_UF and MB1_UF each sit almost on the global common-mode axis that every per-sample composite in this dataset shares (correlation with the mean of all 885 scores %.2f and %.2f); their difference falls to %.2f - see scripts 36 and 37.",
            stats::cor(mb2, gm), stats::cor(mb1, gm), stats::cor(d_spec, gm)),
    "THE TIMING CONTRADICTS THE PREDICTION IT WAS BUILT TO TEST. Script 18 predicted MB2 resemblance would peak at 6 weeks, biogenesis being front-loaded, and reported a pooled genotype main effect (p = 0.026) that averages the two ages. Split by age the specificity is a 12-week effect and there is none at 6, while the resemblance to MB2_UF itself is flat with age. The interaction is not significant (p = 0.15) at n = 6 per cell, so this is a discrepancy to report, not a reversal to claim. It also runs opposite to the TEB result of Fig. 1B, which is a 6-week effect.",
    "Every age comparison is exposed to batch = timepoint. The genotype comparisons within an age - the two brackets - are clean.",
    "GSVA is cohort-relative: a score is a position within these 24 samples. A resemblance to a human tumour state is a statement about relative transcriptional programme scores, not about the tissue being that tumour.",
    "The MB gene lists are human METABRIC signatures mapped to mouse orthologs; the switch was defined by MCbiclust on human breast cancer, not on mouse mammary gland."),
  source = c(
    "results/gsva_scores.rds (scripts/15) -- per-sample scores",
    "Contrast definition and the AP7 hypothesis: scripts/18_ap7_mb_fork_projection.R:96-106; earlier form of this panel, outputs/ap7/mb2_over_mb1_boxplot.pdf",
    "Menegollo, Bentham et al., Cancer Res 2024 (CAN-23-3172), the analytical companion paper"))

save_panel_p(p, "figS1_mb_fork_specificity", width = 50, height = 52)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the contrast and its two components, all four contrasts, in within-group SD
  rbind(cbind(quantity = "MB2_UF", t2), cbind(quantity = "MB1_UF", t1),
        cbind(quantity = "difference", td)) |> print()

  ## the two component scores per sample -- the earlier form of this panel, which
  ## showed WHY the difference is small: Myc moves both forks together
  {
    lv <- c("MB2_UF", "MB1_UF", "MB2_UF - MB1_UF")
    dd <- data.frame(
      group = factor(rep(as.character(sample_meta$group), 3), levels = grp_levels),
      quantity = factor(rep(lv, each = ncol(scores)), levels = lv),
      value = c(mb2, mb1, d_spec))
    ggplot2::ggplot(dd, ggplot2::aes(x = group, y = value)) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey80") +
      ggplot2::geom_boxplot(ggplot2::aes(fill = group), width = 0.6,
                            outlier.shape = NA, linewidth = 0.3,
                            colour = "grey35", alpha = 0.28) +
      ggplot2::geom_jitter(ggplot2::aes(colour = group), width = 0.12, height = 0,
                           size = 0.8, alpha = 0.95) +
      ggplot2::facet_wrap(~ quantity, nrow = 1) +
      ggplot2::scale_fill_manual(values = group_cols, guide = "none") +
      ggplot2::scale_colour_manual(values = group_cols, breaks = names(group_cols),
                                   labels = group_labels) +
      ggplot2::labs(x = NULL, y = "GSVA score", colour = NULL) +
      theme_panel(base_size = 6) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.line.x = ggplot2::element_blank(),
                     legend.position = "bottom")
  }

  ## MB3 is the independent bistate switch, script 18's own specificity control
  ## on the specificity control -- unchanged by genotype would be the clean result
  ct(scores["METABRIC_MB3_HI_CV_GROUP1", ] - scores["METABRIC_MB3_HI_CV_GROUP2", ]) |>
    print()
}
