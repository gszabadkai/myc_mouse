# =============================================================================
# fig2_oxphos_puma_coupling.R -- the respiratory axis against the PUMA ratio, one
# line per genotype, in both fits: a lead that asserts nothing
# -----------------------------------------------------------------------------
# SLOT: Fig. 2I. Rebuilt 2026-09-21 on the author's ruling 7.
#
#   Ruling 7's sentence: "Respiratory priority related to the PUMA to BCL-XL
#   balance differently in the two genotypes, and the difference did not depend
#   on adjustment for epithelial and immune composition."
#
# THE ROBUSTNESS CLAUSE IS WITHDRAWN (author, 2026-09-21, second round). Against
# the within-timepoint permutation null the UNADJUSTED interaction clears and the
# ADJUSTED one does not, so the fit that clears is the one exposed to
# composition. At n = 24 that is an impasse, not a choice between fits: in
# Figure 1 the coupling is hypothesis-generating and asserts nothing. The legend
# carries both permutation results.
#
# WHY IT WAS REBUILT. The 2026-08-05 version drew UNADJUSTED per-genotype lines
# and printed the ADJUSTED model's interaction p beside them. The two fits
# disagree about the wild-type slope -- unadjusted it is negative and nearly as
# steep as the Myc+ one, adjusted it is flat -- so the panel showed one fit and
# quoted the other. Its legend also said six animals per line (there are twelve)
# and that the null permutes genotype labels (it shuffles the axis within each
# timepoint). None of that may ship.
#
# WHAT IT DRAWS NOW. Both fits, side by side, on ONE scale: each axis is a z-score
# over all 24 animals, so a slope is SD of the ratio per SD of the axis and the
# two halves can be read against each other.
#
#   left   unadjusted -- ordinary least squares within each genotype, which is
#          exactly lm(ratio ~ genotype * axis)
#   right  the PRE-SPECIFIED model, lm(ratio ~ genotype * axis + epithelial +
#          immune), covariate coefficients shared across genotypes. Drawn as
#          partial residuals, so the line through each genotype's twelve points IS
#          that model's slope (script 54 asserts it to 1e-10)
#
# The reader sees opposite slopes on the left become flat against positive on the
# right, and a difference between the two lines in each half. That difference is
# printed once per half, from that half's own fit, with its parametric and its
# permutation p side by side (author, 2026-09-21) -- the impasse, on the page --
# and neither half is preferred over the other.
#
# THE QUANTITY IS THE INTERACTION, NOT EITHER SLOPE. The per-genotype slopes were a
# drawing choice, specified in neither form before they were seen; choosing one
# of them now, knowing that the choice decides the sentence, is the move the
# ruling exists to prevent. Both are drawn and neither is quoted in the text.
#
# Reads (read-only, no re-run):
#   results/two_timeline_verification.rds (script 54) -- $coupling_panel,
#       $coupling_lines, $coupling_fits, $coupling_interaction, $coupling_perm,
#       $coupling_decomp, $coupling_cor, $coupling_timepoint
# Output: outputs/figures/panels/fig2_oxphos_puma_coupling.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

tv_path <- here::here("results", "two_timeline_verification.rds")
require_fresher_than(tv_path)
tv <- readRDS(tv_path)

pd <- as.data.frame(tv$coupling_panel)
ln <- as.data.frame(tv$coupling_lines)
cf <- as.data.frame(tv$coupling_fits)
ix <- as.data.frame(tv$coupling_interaction)
pm <- as.data.frame(tv$coupling_perm)
dc <- as.data.frame(tv$coupling_decomp)
cr <- as.data.frame(tv$coupling_cor)
tq <- as.data.frame(tv$coupling_timepoint)

stopifnot(nrow(pd) == 24L, nrow(ln) == 4L,
          all(table(pd$genotype) == 12L),
          setequal(ln$fit, c("unadjusted", "adjusted")))

# --- the numbers the page and the legend carry, read, never typed -------------
FU <- "U unadjusted, within genotype"
FW <- "W epi + imm, within genotype"
FP <- "P pooled, shared covariates: epi + imm"
sl <- function(axis, fit, g, col = "slope")
  cf[[col]][cf$axis == axis & cf$fit == fit & cf$genotype == g]
ia <- function(axis, cv) ix[ix$axis == axis & ix$covariates == cv, ]
i_u <- ia("ox_ppd", "none"); i_a <- ia("ox_ppd", "epi + imm")
r_u <- ia("redox_ppd", "none"); r_a <- ia("redox_ppd", "epi + imm")

# =============================================================================
# the panel
# =============================================================================
# Built by coupling_two_fits() in _panel_common.R, which Fig. 2H+I (alt) part B
# shares, so the two cannot drift apart. It asserts that the drawn lines ARE the
# reported fits and that their difference IS each half's interaction.
p <- coupling_two_fits(tv, key = TRUE)

# =============================================================================
# the legend text (never drawn)
# =============================================================================
pw  <- function(cv) pm[pm$covariates == cv, ]
dwt <- dc[dc$axis == "ox_ppd" & dc$genotype == "neg", ]
cc  <- function(g, pr_) cr$r[cr$genotype == g & cr$pair == pr_]
tqp <- tq[tq$axis == "ox_ppd" & tq$genotype == "pos", ]
# the impasse the legend states, asserted so a re-run cannot flip it silently
stopifnot(pw("none")$p_emp < 0.05, pw("epi + imm")$p_emp >= 0.05)

LEGEND <- panel_legend(
  slot = "Fig. 2I",
  what = paste0(
    "Each animal's OXPHOS mitoPPS score against its PUMA:Bcl-xL log ratio, both ",
    "standardised over the 24 animals, with a line fitted within each genotype: ",
    "without covariates (left) and in the pre-specified model adjusted for the ",
    "epithelial and immune composites (right). Each half prints the difference ",
    "between its two lines, from that half's own fit, with its parametric and its ",
    "permutation p side by side. The panel is a lead and asserts nothing."),
  detail = c(
    sprintf("n = 24 animals; each line is fitted on the TWELVE animals of one genotype (six per age). Both axes are z-scores over all 24, so a slope is SD of the ratio per SD of the axis and the two halves are on one scale. Lines: wild type in blue, Myc+ in vermilion; points in the four-group palette."),
    sprintf("LEFT, UNADJUSTED: ordinary least squares within each genotype, identical to lm(ratio ~ genotype x axis). Wild type %+.2f (p %.3f), Myc+ %+.2f (p %.3f); difference %+.2f (p %.4f).",
            sl("ox_ppd", FU, "neg"), sl("ox_ppd", FU, "neg", "p"),
            sl("ox_ppd", FU, "pos"), sl("ox_ppd", FU, "pos", "p"),
            i_u$interaction, i_u$p),
    sprintf("RIGHT, ADJUSTED: the pre-specified model lm(ratio ~ genotype x axis + epithelial + immune), covariate coefficients shared across genotypes, %d residual df. Points are partial residuals -- the ratio with the two covariate terms removed, covariates centred -- so the line through each genotype's points is exactly that model's slope. Wild type %+.2f (p %.2f), Myc+ %+.2f (p %.4f); difference %+.2f (p %.4f).",
            i_a$df_resid, sl("ox_ppd", FP, "neg"), sl("ox_ppd", FP, "neg", "p"),
            sl("ox_ppd", FP, "pos"), sl("ox_ppd", FP, "pos", "p"),
            i_a$interaction, i_a$p),
    sprintf("THE ADJUSTMENT MOVES BOTH SLOPES. Unadjusted the two slopes are opposite and nearly equal; adjusted the wild-type slope is flat and the Myc+ slope steepens, which is why neither slope is quoted. Adjusting within each genotype separately (own covariate coefficients, %d residual df per genotype) gives the same picture: wild type %+.2f, Myc+ %+.2f.",
            sl("ox_ppd", FW, "neg", "df_resid"), sl("ox_ppd", FW, "neg"), sl("ox_ppd", FW, "pos")),
    sprintf("THE PERMUTATION NULL shuffles the OXPHOS score WITHIN EACH TIMEPOINT (%d times), which keeps the design and breaks only the animal-to-animal pairing. The observed interaction sits at the %.1fth percentile unadjusted (empirical p %.3f) and the %.1fth adjusted (p %.3f), against null medians of %+.2f and %+.2f.",
            pw("none")$n_perm, pw("none")$percentile, pw("none")$p_emp,
            pw("epi + imm")$percentile, pw("epi + imm")$p_emp,
            pw("none")$null_median, pw("epi + imm")$null_median),
    sprintf("THE CONTROL IS THE REDOX mitoPPS AXIS, and what it shows is that its two slopes do NOT differ: interaction %+.2f (p %.2f) unadjusted, %+.2f (p %.2f) adjusted. It is not drawn. Its per-genotype slopes are not a comparator either way: unadjusted wild type %+.2f (p %.4f) against Myc+ %+.2f (p %.2f); adjusted %+.2f and %+.2f, neither distinguishable from zero.",
            r_u$interaction, r_u$p, r_a$interaction, r_a$p,
            sl("redox_ppd", FU, "neg"), sl("redox_ppd", FU, "neg", "p"),
            sl("redox_ppd", FU, "pos"), sl("redox_ppd", FU, "pos", "p"),
            sl("redox_ppd", FP, "neg"), sl("redox_ppd", FP, "pos"))),
  bounds = c(
    sprintf("THE ADJUSTMENT IS FRAGILE, and the legend has to say how. The two covariates are RNA surrogates for composition built from the same count matrix as the ratio (no measured purity exists); they correlate %+.2f with each other in wild type and %+.2f in Myc+; a within-genotype adjusted fit has 8 residual degrees of freedom; and no positive control was carried through the adjustment.",
            cc("neg", "epi ~ imm"), cc("pos", "epi ~ imm")),
    sprintf("WHAT THE ADJUSTMENT REMOVES from the wild-type slope is carried by the immune composite: of a %+.2f change, %+.2f runs through it and %+.2f through the epithelial one. Within the wild-type animals the ratio correlates %+.2f with the immune composite and the OXPHOS score %+.2f.",
            dwt$difference, dwt$via_imm, dwt$via_epi,
            cc("neg", "ratio ~ imm"), cc("neg", "ox_ppd ~ imm")),
    sprintf("BATCH = TIMEPOINT, and it reaches this panel: each line is fitted across both ages. Adjusted for timepoint alone, the unadjusted Myc+ slope falls to %+.2f (p %.2f), so much of it is a difference between the ages rather than between animals. Recorded as an open item, not pursued.",
            tqp$slope_with_tp, tqp$p_with_tp),
    sprintf("AN IMPASSE AT n = 24, NOT A CHOICE. Against the within-timepoint permutation null the UNADJUSTED interaction clears (%.1fth percentile, empirical p %.3f) and the ADJUSTED one does not (%.1fth, p %.3f). The fit that clears is the one exposed to composition, and the fit that removes composition does not clear, so neither can be preferred on these data: the panel is hypothesis-generating and asserts nothing. Both p-values are printed on the page, side by side, because the parametric p alone is anti-conservative at this n and the permutation p alone would hide the discrepancy. Adding a timepoint-by-genotype term to the adjusted model also moves its p from %.4f to %.3f (script 43).",
            pw("none")$percentile, pw("none")$p_emp,
            pw("epi + imm")$percentile, pw("epi + imm")$p_emp,
            tv$tradeoff_43$p, tv$tradeoff_43$p_with_tp),
    "PUMA:Bcl-xL is a TRANSCRIPT ratio. It is not a measurement of how close a cell sits to the apoptotic threshold; BH3 profiling is that measurement.",
    "mitoPPS is a RELATIVE score: a high OXPHOS value means the compartment spends more of its budget there, not that it respires more."),
  source = c(
    "results/two_timeline_verification.rds (scripts/54_two_timeline_verification.R) -- $coupling_panel and $coupling_lines (drawn), $coupling_fits and $coupling_interaction (the three fits on one scale), $coupling_perm (the within-timepoint null), $coupling_decomp, $coupling_cor, $coupling_timepoint"))

save_panel_p(p, "fig2_oxphos_puma_coupling", height = 66)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## all three fits, three axes, one scale
  cf |> print(row.names = FALSE, digits = 3)
  ix |> print(row.names = FALSE, digits = 3)

  ## which covariate carries the wild-type move
  dc |> print(row.names = FALSE, digits = 3)
}
