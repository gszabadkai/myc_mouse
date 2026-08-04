# =============================================================================
# figS2_reallocation_independence.R -- the normal gland's temporal programme and
# Myc's reallocation are unrelated
# -----------------------------------------------------------------------------
# SLOT: Fig. S2B.
#
#   "The WT temporal program operated independently of the MYC-driven
#    reallocation, showing a negligible correlation across the mitochondrial and
#    the whole transcriptome."
#
# TWO SCOPES, ONE CLAIM. Each point is a pathway (left) or a gene (right); the
# horizontal axis is what Myc does at six weeks, the vertical axis is what the
# wild-type gland does between six and twelve. If the two programmes were the
# same programme, or opposites, the cloud would have a slope.
#
#   mitochondrial       143 MitoPathways on the mitoPPS priority ruler, which is
#                       the ruler the word "reallocation" refers to (Fig. 1F).
#                       r = -0.02.
#   whole transcriptome 15,191 genes on the DESeq2 log2 fold change. r = +0.23,
#                       AND THAT NUMBER IS AN ARTEFACT OF THE DESIGN -- see below.
#
# THE SHARED BASELINE, WHICH IS THE WHOLE DIFFICULTY OF THIS PANEL. The two
# contrasts are
#
#     myc_6W        = 6W_myc  - 6W_wt
#     6>12W_wt      = 12W_wt  - 6W_wt
#
# so THE SAME SIX WILD-TYPE ANIMALS ARE SUBTRACTED IN BOTH. A fluctuation in that
# baseline pushes both contrasts the same way, which manufactures a positive
# correlation out of nothing. Two independent estimates of how much:
#
#   (1) FROM THE STANDARD ERRORS. If the four group means are independent with
#       equal variance v, Cov(myc_6W, 6>12W_wt) = Var(6W_wt) = v = SE^2/2, which
#       predicts r = +0.198 of the observed +0.228.
#   (2) FROM THE DATA, and this is the one the panel draws. SPLIT THE SIX
#       WILD-TYPE ANIMALS: estimate the Myc effect against one half and the
#       temporal effect against the other, so no animal is subtracted twice.
#       There are exactly 20 such assignments, so the control is EXHAUSTIVE and
#       needs no seed. Over all 20 the correlation is -0.025 (median -0.021,
#       range -0.44 to +0.16, 60% negative).
#
# So the whole-transcriptome correlation is +0.23 with the baseline shared and
# -0.03 without it -- the same as the mitochondrial value, and the sentence is
# right on both scopes. The dashed line on the right-hand panel is that
# corrected relationship; the solid line is the observed one.
#
# Script 40 uses the same split-baseline technique on a different statistic
# (Issue #6's wild-type convergence, $split_runs), and its $artifact_ledger
# records this class of artefact as "STRUCTURAL -- do not cite". This panel is
# the constructive version of that warning.
#
# BATCH = TIMEPOINT (CLAUDE.md): the wild-type temporal axis is confounded with
# extraction batch, so it is DESCRIBED, not claimed. A batch effect would not
# create independence, though -- it would have to be orthogonal to the Myc
# programme to leave this picture, so the caveat does not undo the negative.
#
# Reads (read-only, no re-run):
#   results/background_vs_myc.rds    (script 40) -- $ruler (both rulers x four
#                                       contrasts), $geometry (identity check)
#   results/interaction_results.rds  (script 03) -- raw DESeqResults, gene level
#   results/gsva_scores.rds          (script 15) -- $expr_mat, the VST matrix the
#                                       split-baseline control is computed on
# Output: outputs/figures/panels/figS2_reallocation_independence.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))
if (!requireNamespace("DESeq2", quietly = TRUE))
  stop("figS2B needs DESeq2 to coerce the DESeqResults in interaction_results.rds")

bv_path <- here::here("results", "background_vs_myc.rds")
gs_path <- here::here("results", "gsva_scores.rds")
ir_path <- here::here("results", "interaction_results.rds")
require_fresher_than(bv_path)
# interaction_results.rds predates the 2026-07-24 gene-symbol reconciliation on
# purpose: that was a set-MEMBERSHIP fix and never touched counts or per-gene
# results, so the Step 1 objects are not stale (the same exemption Fig. S1E takes).

bv <- readRDS(bv_path)
gs <- readRDS(gs_path)
ir <- readRDS(ir_path)

# =============================================================================
# PART 1 -- the mitochondrial scope
# =============================================================================
# Script 40 saves tibbles and its content columns carry names; strip them or a
# row lookup returns a named vector (Fig. 1E's trap).
r <- as.data.frame(bv$ruler)
r$c_m6 <- as.numeric(unname(r$c_m6)); r$c_tn <- as.numeric(unname(r$c_tn))
r143 <- r[!r$is_mtdna, ]                    # the 143 every regression in the corpus uses
stopifnot(nrow(r143) == 143L)

mito_r_prio <- stats::cor(r143$p_m6, r143$p_tn)
mito_r_cont <- stats::cor(r143$c_m6, r143$c_tn)

# Identity check against the analysis of record. Script 40 reports the UNCENTRED
# cosine of the same two vectors; recomputing it here proves the panel is drawing
# script 40's quantities and not a lookalike. Note the two differ on the content
# ruler (cosine +0.210, correlation -0.038) because both content vectors have a
# large positive mean -- the cosine is picking up the shared offset, not a
# relationship, which is why the sentence's word is CORRELATION.
cosine <- function(x, y) sum(x * y) / sqrt(sum(x^2) * sum(y^2))
geo <- as.data.frame(bv$geometry)
geo_all <- geo[geo$scope == "all pathways", ]
stopifnot(nrow(geo_all) == 1L,
          abs(cosine(r143$c_m6, r143$c_tn) - geo_all$c_cos_wt_myc) < 1e-6,
          abs(cosine(r143$p_m6, r143$p_tn) - geo_all$p_cos_wt_myc) < 1e-6)

# =============================================================================
# PART 2 -- the whole transcriptome
# =============================================================================
d6 <- as.data.frame(ir$myc_6W_raw)
dt <- as.data.frame(ir$timepoint_neg_raw)
stopifnot(identical(rownames(d6), rownames(dt)))
ok <- is.finite(d6$log2FoldChange) & is.finite(dt$log2FoldChange) &
      is.finite(d6$lfcSE) & d6$baseMean >= 20
g <- data.frame(a = d6$log2FoldChange[ok], b = dt$log2FoldChange[ok],
                se = d6$lfcSE[ok])
gene_r <- stats::cor(g$a, g$b)

# (1) the standard-error prediction. Var(contrast) = 2v and Cov = v, so v = SE^2/2.
v_hat      <- mean(g$se^2) / 2
r_from_se  <- v_hat / (stats::sd(g$a) * stats::sd(g$b))

# =============================================================================
# PART 3 -- the split-baseline control, exhaustive
# =============================================================================
# Estimate the Myc effect against one half of the six wild-type animals and the
# temporal effect against the other, so no animal is subtracted twice. There are
# choose(6,3) = 20 assignments and all 20 are used, so this is EXHAUSTIVE and has
# no RNG in it -- nothing to seed and nothing to drift.
#
# It runs on the VST matrix rather than on DESeq2 fits because re-fitting is an
# analysis, not a figure. VST is log2-scale, so a difference of group means is
# the same kind of quantity as a log2 fold change -- and the assertion below is
# what licenses the surrogate: with the baseline SHARED, the VST version must
# reproduce the DESeq2 correlation the panel draws.
E  <- gs$expr_mat
sm <- as.data.frame(gs$sample_meta)
stopifnot(identical(as.character(sm$sample), colnames(E)))
E  <- E[rowMeans(E) >= stats::quantile(rowMeans(E), 0.20), , drop = FALSE]
i6n  <- which(sm$group == "6W_neg")
i6p  <- which(sm$group == "6W_pos")
i12n <- which(sm$group == "12W_neg")
stopifnot(length(i6n) == 6L, length(i6p) == 6L, length(i12n) == 6L)

gm <- function(i) rowMeans(E[, i, drop = FALSE])
a_shared <- gm(i6p) - gm(i6n)
b_shared <- gm(i12n) - gm(i6n)
r_shared_vst <- stats::cor(a_shared, b_shared)

assign_h <- utils::combn(6, 3)
split_tab <- as.data.frame(t(apply(assign_h, 2, function(h) {
  A <- i6n[h]; B <- setdiff(i6n, A)
  aa <- gm(i6p) - gm(A); bb <- gm(i12n) - gm(B)
  c(r = stats::cor(aa, bb),
    infl_a = stats::sd(aa) / stats::sd(a_shared),
    infl_b = stats::sd(bb) / stats::sd(b_shared))
})))
r_split <- mean(split_tab$r)
# splitting halves the baseline's n, which adds INDEPENDENT noise to both axes and
# so pulls |r| toward zero. The inflation factors say by how much, and undoing it
# is the conservative direction to check.
r_split_dis <- r_split * mean(split_tab$infl_a) * mean(split_tab$infl_b)

stopifnot(nrow(split_tab) == 20L,
          # THE LICENCE FOR THE SURROGATE: shared-baseline VST reproduces the
          # DESeq2 correlation the panel draws
          abs(r_shared_vst - gene_r) < 0.01)

# =============================================================================
# the panel
# =============================================================================
# Strip labels stay short: a facet strip clips at the panel edge (Fig. 1E's trap).
FAC <- c(sprintf("mitochondrial  (%d pathways)", nrow(r143)),
         sprintf("whole transcriptome  (%s genes)",
                 format(nrow(g), big.mark = ",")))
dd <- rbind(
  data.frame(facet = FAC[1], x = r143$p_m6, y = r143$p_tn),
  data.frame(facet = FAC[2], x = g$a,       y = g$b))
dd$facet <- factor(dd$facet, levels = FAC)

# Windowed for display, as Fig. 1G windows its facets; every number is computed on
# the complete set and the counts left out are in the legend block.
WIN <- list(c(-0.55, 0.55), c(-2.2, 2.2))
inw <- (dd$facet == FAC[1] & abs(dd$x) <= WIN[[1]][2] & abs(dd$y) <= WIN[[1]][2]) |
       (dd$facet == FAC[2] & abs(dd$x) <= WIN[[2]][2] & abs(dd$y) <= WIN[[2]][2])
n_out <- c(sum(dd$facet == FAC[1] & !inw), sum(dd$facet == FAC[2] & !inw))
dd <- dd[inw, ]

# Two lines on the gene facet: what is observed, and what is left once the shared
# baseline is broken. Both pass through the centroid of the complete data.
fit_line <- function(x, y, r_use = NULL) {
  s <- if (is.null(r_use)) stats::cor(x, y) * stats::sd(y) / stats::sd(x) else
    r_use * stats::sd(y) / stats::sd(x)
  c(slope = s, intercept = mean(y) - s * mean(x))
}
L <- rbind(
  data.frame(facet = FAC[1], t(fit_line(r143$p_m6, r143$p_tn)), kind = "observed"),
  data.frame(facet = FAC[2], t(fit_line(g$a, g$b)),             kind = "observed"),
  data.frame(facet = FAC[2], t(fit_line(g$a, g$b, r_split)),
             kind = "baseline not shared"))
L$facet <- factor(L$facet, levels = FAC)
L$kind  <- factor(L$kind, levels = c("observed", "baseline not shared"))

ann <- data.frame(
  facet = factor(FAC, levels = FAC),
  lab = c(sprintf('italic(r)~"%.2f"', mito_r_prio),
          sprintf('italic(r)~"%+.2f"', gene_r)))

p <- ggplot2::ggplot(dd, ggplot2::aes(x, y)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey85") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey85") +
  ggplot2::geom_point(colour = "grey35", size = 0.3, alpha = 0.13, stroke = 0) +
  ggplot2::geom_abline(data = L, ggplot2::aes(slope = slope, intercept = intercept,
                                              linetype = kind),
                       linewidth = 0.35, colour = "grey10") +
  ggplot2::geom_text(data = ann, ggplot2::aes(label = lab), parse = TRUE,
                     x = -Inf, y = Inf, hjust = -0.45, vjust = 1.6,
                     size = 1.9, colour = "grey20", inherit.aes = FALSE) +
  ggplot2::scale_linetype_manual(values = c("observed" = "solid",
                                            "baseline not shared" = "22"),
                                 name = NULL, drop = FALSE) +
  ggplot2::facet_wrap(~ facet, nrow = 1, scales = "free") +
  ggplot2::scale_x_continuous(labels = lab_signed) +
  ggplot2::scale_y_continuous(labels = lab_signed) +
  ggplot2::labs(x = "Myc effect at 6W  (myc_6W)",
                y = "wild-type 6>12W  (6>12W_wt)") +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "plain", size = 6, hjust = 0,
                                       margin = ggplot2::margin(0, 0, 1, 0, "mm")),
    strip.clip = "off",
    legend.position = "bottom",
    legend.key.size = ggplot2::unit(3, "mm"),
    legend.margin = ggplot2::margin(-1, 0, 0, 0, "mm"),
    panel.spacing.x = ggplot2::unit(3, "mm"),
    plot.margin = ggplot2::margin(1.5, 2.5, 0.5, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
LEGEND <- panel_legend(
  slot = "Fig. S2B",
  what = paste0(
    "What the wild-type gland does between six and twelve weeks, against what ",
    "Myc does at six weeks, at two scales: one point per mitochondrial pathway ",
    "on the mitoPPS priority ruler (left) and one point per expressed gene on ",
    "the DESeq2 log2 fold change (right). The solid line is the observed ",
    "relationship; the dashed line is what remains once the two contrasts stop ",
    "sharing a baseline."),
  detail = c(
    sprintf("MITOCHONDRIAL: %d MitoPathways (the synthetic mtDNA-encoded pathway excluded, as in Figs. 1E, 1F and 2F). Correlation on the priority ruler %+.3f, and on the content ruler %+.3f -- both negligible.",
            nrow(r143), mito_r_prio, mito_r_cont),
    sprintf("WHOLE TRANSCRIPTOME: %s genes with baseMean >= 20. The observed correlation is %+.3f, and it is a property of the design rather than of the biology.",
            format(nrow(g), big.mark = ","), gene_r),
    sprintf("WHY: `myc_6W` is 6W_myc minus 6W_wt and `6>12W_wt` is 12W_wt minus 6W_wt, so THE SAME SIX WILD-TYPE ANIMALS ARE SUBTRACTED IN BOTH and a fluctuation in that baseline moves the two contrasts together. If the four group means are independent with equal variance v, the covariance the sharing forces is Var(6W_wt) = v = SE^2/2, which predicts a correlation of %+.3f -- most of the %+.3f observed.",
            r_from_se, gene_r),
    sprintf("THE CONTROL, AND IT IS EXHAUSTIVE: estimate the Myc effect against one half of the six wild-type animals and the temporal effect against the other, so no animal is subtracted twice. There are choose(6,3) = %d such assignments and all %d are used, so there is no seed and nothing to drift. The correlation falls to a mean of %+.3f (median %+.3f, range %+.3f to %+.3f, %.0f%% of assignments negative) -- the same value as the mitochondrial scope.",
            nrow(split_tab), nrow(split_tab), r_split, stats::median(split_tab$r),
            min(split_tab$r), max(split_tab$r), 100 * mean(split_tab$r < 0)),
    sprintf("The control runs on the VST matrix rather than on re-fitted DESeq2 models, because re-fitting is an analysis and not a figure. What licenses the surrogate is asserted in the script: with the baseline SHARED, the VST version reproduces the drawn DESeq2 correlation to %.3f (%+.3f against %+.3f).",
            abs(r_shared_vst - gene_r), r_shared_vst, gene_r),
    sprintf("Splitting halves the baseline's sample size, which adds independent noise to both axes and pulls the correlation toward zero; the standard deviations grow by %.0f%% and %.0f%%, and undoing that inflation gives %+.3f rather than %+.3f. The conclusion does not depend on which of the two is used.",
            100 * (mean(split_tab$infl_a) - 1), 100 * (mean(split_tab$infl_b) - 1),
            r_split_dis, r_split),
    sprintf("THE ARTEFACT PUSHES THE WRONG WAY FOR THE MITOCHONDRIAL PANEL TOO, which makes that negative conservative: the priority contrasts share the same wild-type baseline, so the sharing forces a POSITIVE correlation, and what is observed is %+.3f.",
            mito_r_prio),
    sprintf("Script 40 reports the UNCENTRED cosine of these same two vectors, and the script asserts the recomputation matches it exactly (%+.4f content, %+.4f priority). On the content ruler cosine and correlation disagree (%+.3f against %+.3f) because both content vectors have a large positive mean and the cosine picks up that shared offset. The sentence's word is CORRELATION, so the correlation is what is drawn.",
            geo_all$c_cos_wt_myc, geo_all$p_cos_wt_myc, geo_all$c_cos_wt_myc,
            mito_r_cont)),
  bounds = c(
    "BATCH = TIMEPOINT. The wild-type temporal axis is confounded with extraction batch, so it is DESCRIBED, not claimed. That caveat does not undo this negative, though: a batch effect would have to be orthogonal to the Myc programme to leave a correlation of zero, which is a stronger coincidence than the result it would explain away.",
    "THIS IS NOT THE STATEMENT THAT THE Myc+ GLAND'S OWN TEMPORAL CHANGE IS INDEPENDENT OF THE WILD-TYPE ONE. That regression (`shared`, Myc+time ~ WTtime, script 40) has a slope of 0.749 -- much of the Myc+ gland's drift IS the wild-type drift. The two claims are about different pairs of vectors and must not be run together.",
    "A correlation of zero is not independence in the mechanistic sense. It says the two programmes do not move the same pathways in the same direction; it does not say they cannot interact, and the interaction contrast is where that question is asked.",
    "The split-baseline control is a control on the CORRELATION only. It uses VST group means, not DESeq2's negative-binomial fits, so it must not be quoted as an effect size.",
    "n = 6 per group. Every correlation here is over pathways or genes, not over animals, so its precision is not the n = 6 precision -- but the underlying contrasts are, and the split control's spread (sd 0.15 across assignments) is a direct picture of that."),
  source = c(
    "results/background_vs_myc.rds (scripts/40_background_vs_myc_decomposition.R) -- $ruler for both rulers on both contrasts; $geometry for the cosine identity check; $split_runs is the same split-baseline technique applied to a different statistic",
    "results/interaction_results.rds (scripts/03_deseq_results_qc.R) -- raw (unshrunken) DESeqResults for myc_6W and timepoint_neg",
    "results/gsva_scores.rds (scripts/15_gsva_scoring.R) -- $expr_mat, the VST matrix the exhaustive split-baseline control is computed on"))

save_panel_p(p, "figS2_reallocation_independence", height = 52)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the 20 split assignments, one row each -- the spread is the n = 3 baseline
  split_tab[order(split_tab$r), ] |> print(row.names = FALSE, digits = 3)

  ## the same decomposition for the other three contrast pairs. Pairs sharing a
  ## term with the SAME sign are forced positive, with OPPOSITE signs negative --
  ## and the two involving timepoint_pos have real residuals, because
  ## timepoint_pos IS the interaction plus timepoint_neg by construction.
  gl <- function(n) as.data.frame(ir[[n]])$log2FoldChange[ok]
  se <- function(n) as.data.frame(ir[[n]])$lfcSE[ok]
  for (pr in list(c("myc_6W_raw", "timepoint_neg_raw", "+1"),
                  c("myc_12W_raw", "timepoint_neg_raw", "-1"),
                  c("myc_6W_raw", "timepoint_pos_raw", "-1"),
                  c("myc_12W_raw", "timepoint_pos_raw", "+1"))) {
    x <- gl(pr[1]); y <- gl(pr[2])
    v <- (mean(se(pr[1])^2) + mean(se(pr[2])^2)) / 4
    cat(sprintf("%-12s vs %-18s obs %+.3f  forced %+.3f  residual %+.3f\n",
                pr[1], pr[2], stats::cor(x, y),
                as.numeric(pr[3]) * v / (stats::sd(x) * stats::sd(y)),
                (stats::cov(x, y) - as.numeric(pr[3]) * v) /
                  (stats::sd(x) * stats::sd(y))))
  }

  ## script 40's own geometry, which the panel checks itself against
  as.data.frame(bv$geometry) |> print(row.names = FALSE, digits = 3)
  as.data.frame(bv$geometry_null) |> print(row.names = FALSE, digits = 3)
  as.data.frame(bv$artifact_ledger) |> print(row.names = FALSE, digits = 3)
}
