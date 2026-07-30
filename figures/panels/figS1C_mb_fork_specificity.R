# =============================================================================
# figS1C_mb_fork_specificity.R -- is the human resemblance MYC-specific, or just
# mitochondrial biogenesis?
# -----------------------------------------------------------------------------
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
# The only feature distinguishing them is MYC. So MB2_UF against MB1_UF isolates
# "Myc drives the Myc fork" from "Myc drives biogenesis generically", and it is
# the decisive contrast -- MB2_UF on its own cannot make that distinction.
#
# WHAT THE PANEL SHOWS, AND WHY IT IS DRAWN IN RAW GSVA UNITS. The two upper-fork
# scores correlate 0.982 (they share 230 genes of 626 and 419), and the Myc effect
# on each is large and almost identical. Their difference is therefore a small
# residual: within-group SD 0.082 against 0.335 and 0.297 for the scores
# themselves, about a quarter. Standardising each row by its own within-group SD
# -- the convention everywhere else in this figure set -- would put all three on
# one axis and make the residual look as large as the thing it is a residual of,
# which is exactly the contrast this panel exists to NOT manufacture. So the three
# rows share ONE RAW GSVA AXIS and the difference is drawn small, because it is.
# The SD-standardised numbers are in the legend block for comparability with 1B.
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
# Output: outputs/figures/panels/figS1C_mb_fork_specificity.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

gsva        <- readRDS(here::here("results", "gsva_scores.rds"))
scores      <- gsva$scores
sample_meta <- as.data.frame(gsva$sample_meta)[colnames(scores), , drop = FALSE]

mb <- c(MB2_UF = "METABRIC_MB2_HI_CV_GROUP1",
        MB1_UF = "METABRIC_MB1_HI_CV_GROUP1")
stopifnot(all(mb %in% rownames(scores)))

mb2 <- scores[mb[["MB2_UF"]], ]
mb1 <- scores[mb[["MB1_UF"]], ]

# results/ap7_mb_fork.rds is a Jul-6 object and predates the gene-symbol
# reconciliation, so it is deliberately NOT gated by require_fresher_than() -- but
# the MB HI_CV sets were untouched by that rebuild and the two agree exactly, so
# it is still a valid guard on the contrast's definition (scripts/18:105).
ap7_path <- here::here("results", "ap7_mb_fork.rds")
if (file.exists(ap7_path)) {
  ap7 <- as.data.frame(readRDS(ap7_path)$fork_df)
  rownames(ap7) <- ap7$sample
  stopifnot(max(abs((mb2 - mb1)[rownames(ap7)] - ap7$MB2_over_MB1_UF)) < 1e-8)
}

row_levels <- c("MB2_UF", "MB1_UF", "MB2_UF - MB1_UF")

dat <- data.frame(
  sample   = rep(colnames(scores), 3),
  group    = rep(as.character(sample_meta$group), 3),
  quantity = factor(rep(row_levels, each = ncol(scores)), levels = row_levels),
  value    = c(mb2, mb1, mb2 - mb1),
  stringsAsFactors = FALSE)
dat$group <- factor(dat$group, levels = rev(names(group_cols)))
stopifnot(nrow(dat) == 72, !any(is.na(dat$group)))

# --- panel -------------------------------------------------------------------
pts_layer <- if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
  ggbeeswarm::geom_quasirandom(ggplot2::aes(fill = group), orientation = "y",
                               width = 0.24, shape = 21, size = 1.2,
                               colour = "grey25", stroke = 0.18)
} else {
  ggplot2::geom_jitter(ggplot2::aes(fill = group), height = 0.15, width = 0,
                       shape = 21, size = 1.2, colour = "grey25", stroke = 0.18)
}

p <- ggplot2::ggplot(dat, ggplot2::aes(x = value, y = group)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey80") +
  pts_layer +
  ggplot2::stat_summary(fun = mean,
                        fun.min = function(v) mean(v) - stats::sd(v) / sqrt(length(v)),
                        fun.max = function(v) mean(v) + stats::sd(v) / sqrt(length(v)),
                        geom = "errorbar", width = 0, linewidth = 0.4,
                        colour = "grey15") +
  ggplot2::stat_summary(fun = mean, geom = "point", shape = 124, size = 2.2,
                        colour = "grey15") +
  ggplot2::facet_wrap(~ quantity, nrow = 1) +
  ggplot2::scale_fill_manual(values = group_cols, guide = "none") +
  ggplot2::scale_y_discrete(labels = group_labels) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.09),
                              guide = ggplot2::guide_axis(check.overlap = TRUE)) +
  ggplot2::labs(x = "GSVA score", y = NULL) +
  theme_panel(base_size = 6) +
  ggplot2::theme(axis.line.y  = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 strip.text.x = ggplot2::element_text(face = "bold", size = 5.8),
                 panel.spacing.x = ggplot2::unit(1.6, "mm"),
                 plot.margin  = ggplot2::margin(1.5, 1.5, 1, 1, "mm"))

# --- the legend text (never drawn) -------------------------------------------
ct <- function(y) contrast_table(y, sample_meta$timepoint, sample_meta$myc_status,
                                 sample_meta$group)
t2 <- ct(mb2); t1 <- ct(mb1); td <- ct(mb2 - mb1)
g  <- function(t, cn) t$effect[t$contrast == cn]
gp <- function(t, cn) t$p[t$contrast == cn]
f  <- function(x) sprintf("%+.2f", x)
gm <- colMeans(scores)

LEGEND <- panel_legend(
  slot = "Fig. S1C",
  what = paste0(
    "Specificity control for the human breast-cancer resemblance of Fig. 1B. ",
    "Per-sample GSVA scores for the two upper forks of the MCbiclust multistate ",
    "switch - MB2, the MYC arm, and MB1, the same mitochondrial-biogenesis-high ",
    "proliferative state without MYC - and their difference."),
  detail = c(
    "n = 6 animals per group, n = 24; each point is one animal, tick and bar are the group mean and its standard error.",
    sprintf("MB2_UF is METABRIC_MB2_HI_CV_GROUP1 (%d mouse-mapped genes) and MB1_UF is METABRIC_MB1_HI_CV_GROUP1 (%d); they share %d genes. Both are the highest-variance gene lists defining the upper fork of their switch (Menegollo, Bentham et al., Cancer Res 2024; scripts/18).",
            length(gsva$pathways[[mb[["MB2_UF"]]]]),
            length(gsva$pathways[[mb[["MB1_UF"]]]]),
            length(intersect(gsva$pathways[[mb[["MB2_UF"]]]],
                             gsva$pathways[[mb[["MB1_UF"]]]]))),
    "All three rows share one RAW GSVA axis, not the within-group SD scaling used in Fig. 1B, so that the difference is drawn at its true size relative to the scores it is a difference of.",
    sprintf("Myc raises BOTH forks, and by almost the same amount. In within-group SD: MB2_UF %s at 6 weeks (p = %.3f) and %s at 12 (p = %.4f); MB1_UF %s (p = %.3f) and %s (p = %.3f).",
            f(g(t2,"myc_6W")), gp(t2,"myc_6W"), f(g(t2,"myc_12W")), gp(t2,"myc_12W"),
            f(g(t1,"myc_6W")), gp(t1,"myc_6W"), f(g(t1,"myc_12W")), gp(t1,"myc_12W")),
    sprintf("The specificity contrast, MB2_UF minus MB1_UF, is absent at 6 weeks (%s SD, p = %.2f) and present at 12 (%s SD, p = %.4f); the interaction is p = %.2f.",
            f(g(td,"myc_6W")), gp(td,"myc_6W"),
            f(g(td,"myc_12W")), gp(td,"myc_12W"), gp(td,"interaction")),
    sprintf("Neither genotype moves along the contrast with age (%s SD in wild type, p = %.2f; %s SD in Myc+, p = %.2f).",
            f(g(td,"6>12W_wt")), gp(td,"6>12W_wt"),
            f(g(td,"6>12W_myc")), gp(td,"6>12W_myc"))),
  bounds = c(
    sprintf("The two scores correlate %.3f across the 24 samples, so their difference is a small residual: within-group SD %.3f against %.3f and %.3f. Almost everything the two forks have in common cancels, which is the design of the test but also means the contrast is estimated on a fraction of the signal.",
            stats::cor(mb2, mb1), unique(td$within_sd), unique(t2$within_sd),
            unique(t1$within_sd)),
    sprintf("That cancellation is what makes the contrast worth having. MB2_UF and MB1_UF each sit almost on the global common-mode axis that every per-sample composite in this dataset shares (correlation with the mean of all 885 scores %.2f and %.2f); their difference falls to %.2f. The difference is the only one of the three rows that is not largely the common-mode axis - see scripts 36 and 37.",
            stats::cor(mb2, gm), stats::cor(mb1, gm), stats::cor(mb2 - mb1, gm)),
    "THE TIMING CONTRADICTS THE PREDICTION IT WAS BUILT TO TEST. Script 18 predicted MB2 resemblance would peak at 6 weeks, biogenesis being front-loaded, and reported a pooled genotype main effect (p = 0.026) that averages the two ages. Split by age the specificity is a 12-week effect and there is none at 6, while the resemblance to MB2_UF itself is flat with age. The interaction is not significant (p = 0.15) at n = 6 per cell, so this is a discrepancy to report, not a reversal to claim.",
    "Every age comparison is exposed to batch = timepoint. The genotype comparisons within an age are clean.",
    "GSVA is cohort-relative: a score is a position within these 24 samples. A resemblance to a human tumour state is a statement about relative transcriptional programme scores, not about the tissue being that tumour.",
    "The MB gene lists are human METABRIC signatures mapped to mouse orthologs; the switch was defined by MCbiclust on human breast cancer, not on mouse mammary gland."),
  source = c(
    "results/gsva_scores.rds (scripts/15) -- per-sample scores",
    "Contrast definition and the AP7 hypothesis: scripts/18_ap7_mb_fork_projection.R:96-106; earlier form of this panel, outputs/ap7/mb2_over_mb1_boxplot.pdf",
    "Menegollo, Bentham et al., Cancer Res 2024 (CAN-23-3172), the analytical companion paper"))

save_panel_p(p, "figS1C_mb_fork_specificity",
             width = fig_w[["single"]], height = 46)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the three rows, all four contrasts, in within-group SD
  rbind(cbind(quantity = "MB2_UF", t2), cbind(quantity = "MB1_UF", t1),
        cbind(quantity = "difference", td)) |> print()

  ## the raw group means, which is what the panel draws
  dat |>
    dplyr::group_by(quantity, group) |>
    dplyr::summarise(mean = mean(value), .groups = "drop") |>
    tidyr::pivot_wider(names_from = group, values_from = mean) |> print()

  ## MB3 is the independent bistate switch, script 18's own specificity control
  ## on the specificity control -- unchanged by genotype would be the clean result
  ct(scores["METABRIC_MB3_HI_CV_GROUP1", ] - scores["METABRIC_MB3_HI_CV_GROUP2", ]) |>
    print()
}
