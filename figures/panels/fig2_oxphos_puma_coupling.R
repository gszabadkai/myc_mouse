# =============================================================================
# fig2_oxphos_puma_coupling.R -- the respiratory axis and the PUMA ratio move
# together in one genotype and apart in the other
# -----------------------------------------------------------------------------
# SLOT: Fig. 2I. (Author's call, 2026-08-04: this gets a panel in the main
# figure, not a supplementary.)
#
#   "... the interaction between MYC and OXPHOS subunit coupling to the
#    PUMA/Bcl-XL ratio is highly significant (p = 0.0052), whereas no such link
#    exists for other redox and metabolic axes."
#
# WHAT AN INTERACTION BETWEEN A GENOTYPE AND A COUPLING LOOKS LIKE: two slopes
# that differ. The panel draws them. Each point is one animal; the horizontal
# axis is that animal's OXPHOS-subunit mitoPPS score (left) or its redox score
# (right), and the vertical axis is its PUMA:Bcl-xL log ratio. A line is fitted
# within each genotype, and the statistic the sentence quotes is the difference
# between the two.
#
#   OXPHOS subunits   wild type -4.80, Myc+ +4.35   -- OPPOSITE SIGNS
#   redox             wild type -7.70, Myc+ -6.72   -- parallel
#
# REDOX IS THE CONTROL AND IT IS NOT A NULL AXIS. Both genotypes couple to redox,
# and strongly; what redox does not do is couple DIFFERENTLY in the two. That is
# the right control for an interaction, and it is a stronger one than an axis
# nothing couples to, because it shows the interaction is not simply a
# consequence of the ratio correlating with mitochondrial scores in general.
#
# THE PANEL IS A RECONSTRUCTION AND SAYS SO. Script 43 fits on its own log
# matrix; this rebuilds the ratio from the VST matrix, which reproduces the
# recorded interaction to within 1.2% (asserted below). The number quoted in the
# text should be script 43's, not the panel's.
#
# READ IT AS A LEAD. The empirical null in the same object puts this interaction
# at the 91.7th percentile of 5,000 permutations (p = 0.083), and adding a
# timepoint term takes p from 0.0052 to 0.088. It is the one axis-by-genotype
# interaction in the corpus that reaches nominal significance, and it is a
# ranking statement, not a confirmed one.
#
# Reads (read-only, no re-run):
#   results/substrate_specificity_tradeoff.rds (script 43) -- $tradeoff, the
#                                            statistic; $tradeoff_perm, its null
#   results/priming_arm_teb.rds                (script 42) -- $axis_scores,
#                                            $purity
#   results/gsva_scores.rds                    (script 15) -- $expr_mat, to
#                                            rebuild the per-animal ratio
# Output: outputs/figures/panels/fig2_oxphos_puma_coupling.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

ss_path <- here::here("results", "substrate_specificity_tradeoff.rds")
pa_path <- here::here("results", "priming_arm_teb.rds")
gs_path <- here::here("results", "gsva_scores.rds")
require_fresher_than(ss_path)
require_fresher_than(pa_path)
ss <- readRDS(ss_path); pa <- readRDS(pa_path); gs <- readRDS(gs_path)

tr  <- as.data.frame(ss$tradeoff)
trp <- as.data.frame(ss$tradeoff_perm)
ax  <- as.data.frame(pa$axis_scores)
pu  <- as.data.frame(pa$purity)
E   <- gs$expr_mat
sm  <- as.data.frame(gs$sample_meta)

OUT <- "Bbc3:Bcl2l1 (PUMA priming)"
stopifnot(identical(as.character(ax$sample), colnames(E)),
          identical(as.character(pu$sample), colnames(E)),
          OUT %in% tr$outcome,
          all(c("Bbc3", "Bcl2l1") %in% rownames(E)))

# =============================================================================
# the per-animal quantities
# =============================================================================
# The outcome is the log2 ratio of the two transcripts, z-scored exactly as
# script 43 z-scores it. The axes are script 42's saved mitoPPS scores, used as
# they are.
d <- data.frame(
  sample = colnames(E),
  group  = as.character(sm$group),
  myc    = factor(as.character(sm$myc_status), levels = c("neg", "pos")),
  y      = as.numeric(scale(as.numeric(E["Bbc3", ] - E["Bcl2l1", ]))),
  epi    = as.numeric(pu$epithelial),
  imm    = as.numeric(pu$immune),
  oxphos_ppd = ax$oxphos_ppd,
  redox_ppd  = ax$redox_ppd,
  stringsAsFactors = FALSE)
stopifnot(nrow(d) == 24L)

# --- the reconstruction, proved against the analysis of record ----------------
# Script 43's model is lm(y ~ myc * a + epi + imm) and the statistic is the
# `mycpos:a` term. Refitting it here must land on the recorded value; it will not
# land exactly, because script 43 fits on its own log matrix and this rebuilds
# the ratio from the VST one, so the check is a tolerance and the legend says
# which number belongs in the text.
refit <- function(a) {
  m <- summary(stats::lm(y ~ myc * a + epi + imm,
                         data.frame(y = d$y, myc = d$myc, a = d[[a]],
                                    epi = d$epi, imm = d$imm)))$coefficients
  c(est = m["mycpos:a", 1], p = m["mycpos:a", 4])
}
chk <- vapply(c("oxphos_ppd", "redox_ppd"), refit, numeric(2))
rec <- tr[tr$outcome == OUT, ]
rec <- rec[match(c("oxphos_ppd", "redox_ppd"), rec$axis), ]
rel <- abs(chk["est", ] - rec$myc_x_axis) / abs(rec$myc_x_axis)
# The CLAIM is checked tightly and the CONTROL qualitatively, because a relative
# error on a near-null coefficient is not a meaningful quantity: redox's 2.02
# against 2.16 is 6.8% of a number whose p is 0.73. What has to hold for the
# control is that it still reads as nothing.
stopifnot(rel[["oxphos_ppd"]] < 0.05,            # observed 1.2%
          chk["p", "redox_ppd"] > 0.5,           # the control, in the rebuild
          rec$p[1] < 0.01, rec$p[2] > 0.5)       # the claim, from the record

# =============================================================================
# the panel
# =============================================================================
FAC <- c("OXPHOS subunits", "redox")
long <- rbind(
  data.frame(d[, c("group", "myc", "y")], a = d$oxphos_ppd, facet = FAC[1]),
  data.frame(d[, c("group", "myc", "y")], a = d$redox_ppd,  facet = FAC[2]))
long$facet <- factor(long$facet, levels = FAC)
long$group <- factor(long$group, levels = names(group_cols))

p <- ggplot2::ggplot(long, ggplot2::aes(a, y)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey85") +
  # one line per genotype: the interaction IS the difference between them
  ggplot2::geom_smooth(ggplot2::aes(colour = myc), method = "lm", formula = y ~ x,
                       se = FALSE, linewidth = 0.45) +
  ggplot2::geom_point(ggplot2::aes(fill = group), shape = 21, size = 1.6,
                      stroke = 0.25, colour = "grey25") +
  ggplot2::scale_colour_manual(values = geno_cols, guide = "none") +
  ggplot2::scale_fill_manual(values = group_cols, labels = group_labels,
                             breaks = names(group_cols), name = NULL) +
  ggplot2::facet_wrap(~ facet, nrow = 1, scales = "free_x") +
  ggplot2::scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
  ggplot2::scale_y_continuous(labels = lab_signed) +
  # "mitoPPS", not "priority score" -- the author's naming, 2026-08-09, so this
  # axis and Fig. 2H+I (alt)'s call the same quantity the same thing. The label
  # stays generic because the facets are OXPHOS and redox.
  ggplot2::labs(x = "mitoPPS score, per animal",
                y = "PUMA:Bcl-xL  (log2 ratio, z)") +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1,
                                               override.aes = list(size = 1.7))) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "plain", size = 6, hjust = 0,
                                       margin = ggplot2::margin(0, 0, 1, 0, "mm")),
    strip.clip = "off",
    legend.position = "bottom",
    legend.key.size = ggplot2::unit(2.6, "mm"),
    legend.margin   = ggplot2::margin(-1.5, 0, 0, 0, "mm"),
    panel.spacing.x = ggplot2::unit(3, "mm"),
    plot.margin     = ggplot2::margin(1.5, 2.5, 0.5, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
slope <- function(a, g) {
  k <- d$myc == g
  stats::coef(stats::lm(d$y[k] ~ d[[a]][k]))[2]
}
perm <- trp[trp$outcome == OUT, ]
perm <- perm[match(c("oxphos_ppd", "redox_ppd"), perm$axis), ]
amb  <- as.data.frame(ss$ambient)

LEGEND <- panel_legend(
  slot = "Fig. 2I",
  what = paste0(
    "One point per animal: its mitoPPS score on the horizontal axis -- ",
    "OXPHOS subunits on the left, redox on the right -- against its ",
    "PUMA:Bcl-xL log2 ratio, standardised. A line is fitted within each ",
    "genotype; the statistic the text quotes is the difference between the two ",
    "slopes."),
  detail = c(
    sprintf("n = 24 animals, 6 per group. The outcome is log2(Bbc3) - log2(Bcl2l1) per animal, z-scored; the axes are the mitoPPS pairwise-ratio scores for OXPHOS subunits and for ROS/glutathione metabolism. The model behind the quoted statistic is `ratio ~ genotype * axis + epithelial + immune`, so it is adjusted for cell composition; the drawn lines are unadjusted within-genotype fits and are the picture of the same thing.",
            NULL),
    sprintf("THE INTERACTION IS THE PANEL: on the OXPHOS axis the two genotypes couple in OPPOSITE directions -- wild type %+.2f, Myc+ %+.2f -- and on the redox axis they do not (%+.2f and %+.2f). Script 43's adjusted interaction terms are %+.2f (p = %.4f) for OXPHOS and %+.2f (p = %.2f) for redox.",
            slope("oxphos_ppd", "neg"), slope("oxphos_ppd", "pos"),
            slope("redox_ppd", "neg"),  slope("redox_ppd", "pos"),
            rec$myc_x_axis[1], rec$p[1], rec$myc_x_axis[2], rec$p[2]),
    sprintf("REDOX IS NOT A NULL AXIS, WHICH IS WHY IT IS THE RIGHT CONTROL. Both genotypes couple to it, and strongly (wild-type slope %+.2f). What it does not do is couple DIFFERENTLY between them. An axis nothing couples to would not have excluded the possibility that the PUMA ratio simply tracks mitochondrial scores in general.",
            slope("redox_ppd", "neg")),
    sprintf("THE EMPIRICAL NULL PUTS IT AT THE %.1fST PERCENTILE, not past it. Script 43 permutes the genotype labels 5,000 times: the observed OXPHOS interaction sits at percentile %.1f (empirical p = %.3f) against a null median of %+.2f, and the redox one at percentile %.1f (p = %.2f). Adding a timepoint term moves the OXPHOS p from %.4f to %.4f.",
            perm$percentile[1], perm$percentile[1], perm$p_emp[1],
            perm$null_median[1], perm$percentile[2], perm$p_emp[2],
            rec$p[1], rec$p_with_tp[1]),
    sprintf("The raw Spearman correlations are not the statistic and should not be quoted as one: pooled across genotypes the PUMA ratio correlates %+.2f with the OXPHOS axis and %+.2f with redox, which is the pooling artefact two opposite slopes produce. Design-adjusted they are %+.2f and %+.2f.",
            rec$rho_raw[1], rec$rho_raw[2], rec$rho_adj[1], rec$rho_adj[2]),
    sprintf("Ambient scale for the axis: the OXPHOS-subunit score correlates with an arbitrary one of the other %d mitoPPS pathways at a median |rho| of %.2f (90th percentile %.2f). A coupling of that size is what the compartment hands you, which is why the claim here is about a DIFFERENCE BETWEEN GENOTYPES rather than about a correlation.",
            amb$n_pathways[1], amb$ambient_median_abs_rho[1], amb$ambient_q90[1]),
    "THE PANEL IS A RECONSTRUCTION. Script 43 fits on its own log-expression matrix; this rebuilds the ratio from the VST matrix, and the refitted interaction reproduces the recorded one to within 1.2% on the OXPHOS axis (asserted in the script). The number for the text is script 43's."),
  bounds = c(
    "A LEAD, NOT A RESULT. This is the one axis-by-genotype interaction in the corpus that reaches nominal significance, but it does not survive its own permutation null (p = 0.083) and it weakens when timepoint is added (p = 0.088). At n = 6 per cell an interaction between a genotype and a slope is the least powered thing the design can be asked for.",
    "\"NO SUCH LINK FOR OTHER REDOX AND METABOLIC AXES\" IS TESTED ON ONE COMPARATOR, NOT SEVERAL. Script 43's trade-off analysis carries exactly two axes -- OXPHOS subunits and redox -- so the sentence should say \"for the redox axis\" unless a wider panel of axes is run. What IS broader is the outcome side: the same two axes were tested against five outcomes, and only this one reaches p < 0.05.",
    "mitoPPS is a RELATIVE score: a high OXPHOS-subunit value means the compartment spends more of its budget there, not that it respires more. The coupling is therefore between a resource-allocation state and a transcript ratio, and neither is a rate measurement.",
    "CORRELATION AT n = 6 PER GENOTYPE. Each drawn line is fitted on six animals within an age-mixed genotype group of twelve; the quoted model pools both ages and adjusts for composition. Neither slope should be read as an effect size.",
    "PUMA:Bcl-xL is a transcript ratio, not priming. The measurement of priming is BH3 profiling, and this panel is one of the reasons to do it."),
  source = c(
    "results/substrate_specificity_tradeoff.rds (scripts/43_substrate_specificity_and_tradeoff.R) -- $tradeoff for the interaction terms, $tradeoff_perm for the permutation null, $ambient for the axis's ambient coupling",
    "results/priming_arm_teb.rds (scripts/42_priming_arm_and_teb_substrate.R) -- $axis_scores (the per-animal mitoPPS axes) and $purity (the epithelial and immune composites the model adjusts for)",
    "results/gsva_scores.rds (scripts/15_gsva_scoring.R) -- $expr_mat, used to rebuild the per-animal PUMA:Bcl-xL ratio"))

save_panel_p(p, "fig2_oxphos_puma_coupling", height = 56)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the reconstruction against the record
  data.frame(axis = c("oxphos_ppd", "redox_ppd"),
             recorded = rec$myc_x_axis, refitted = chk["est", ],
             rel_diff = rel) |> print(row.names = FALSE, digits = 3)

  ## every outcome x axis script 43 tested -- this panel is one row of ten
  tr |> print(row.names = FALSE, digits = 3)
  trp |> print(row.names = FALSE, digits = 3)

  ## the within-genotype slopes the panel draws
  for (a in c("oxphos_ppd", "redox_ppd"))
    for (g in c("neg", "pos"))
      cat(sprintf("%-11s %s  slope %+.2f\n", a, g, slope(a, g)))
}
