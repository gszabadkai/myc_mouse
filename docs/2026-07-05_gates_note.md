---
title: Block A Day-1 gates - results and framing decisions
date: 2026-07-05
tags: [myc_mouse, block-A, gates, decisions, trajectory]
status: decision-note
---

# Block A Day-1 gates: results and framing decisions

Day-1 decision layer (scripts 13, 14) sourced in Positron. Gate 1 = SHRINKS
(genotypes maximally divergent at 6W, largely converged by 12W). Gate 2 =
DECISION_POINT (pro-apoptotic module coordinately shifted). Consequences below.
Companion inputs: results/gate1_divergence_timing.rds,
results/gate2_apoptosis_readout.rds; outputs/gates/gate{1,2}_summary.csv.

## Gate 1 - divergence timing (script 13)

Facts:
- Genotype (Myc+ vs Myc-) effect at 6W: 2777 genes padj<0.1 (1967 at <0.05).
- Genotype effect at 12W: 239 genes padj<0.1. Trend: SHRINKS (-91%).
- Interaction direction, 357 unique-to-Myc+ (p<0.05) set: 55.5% negative
  (Myc effect weakens at 12W), median interaction LFC -0.299, binomial p=0.044.
- Broader context, all genes interaction p<0.05: 61.9% negative.

Reading: the genotype's transcriptional footprint is front-loaded - largest at
the earliest timepoint and attenuating by 12W, consistent with the
predominantly-negative interaction. Myc dose is stable across 6W-12W
(mRNA + blot + literature), so this is NOT a driver-dose decline. Constant
driver, moving substrate.

Effect-size check (script 13 PART 2b) - RESOLVED: convergence is BIOLOGICAL,
not a power artifact.
- Concern: myc_6W_raw is the reference-level contrast (name =
  myc_status_pos_vs_neg); myc_12W_raw is a composite (list) contrast with a
  larger SE and therefore LOWER power, so the 2777 -> 239 COUNT drop is
  confounded with reduced 12W power. LFC POINT estimates are SE-independent, so
  the check compares |LFC| across FIXED gene sets (not re-thresholded at 12W).
- Results (median |LFC| 6W -> 12W; slope of lm(lfc_12W ~ 0 + lfc_6W); Pearson r):
    Set A  6W-divergent (padj<0.1, n=2777): 0.62 -> 0.27, ratio 0.44,
           slope 0.46, r 0.70, 97% of genes attenuated. NOTE selection-on-6W
           biases this set TOWARD apparent attenuation.
    Set B  expressed (baseMean>=median, n=9262, LFC-INDEPENDENT/UNBIASED):
           0.22 -> 0.16, ratio 0.72, slope 0.45, r 0.64, 64% attenuated.
- Decisive point: the two SLOPES are identical (0.46 biased vs 0.45 unbiased).
  Regression-to-the-mean would make the 6W-selected set A attenuate MORE than
  the LFC-independent set B; they agree, so the ~2x shrinkage is NOT a
  selection/power artifact. The genotype effect genuinely attenuates to ~45% of
  its 6W magnitude by 12W, direction preserved (positive r -> attenuation, not
  sign reversal). Verdict class: COLLAPSED_BIOLOGICAL.
- Honest residual caveats: (a) "converge" is shorthand - the effect HALVES, it
  does NOT vanish; quote the SLOPE (~0.45), not set A's median ratio (0.44),
  which is composition-inflated vs set B's noise-floor-padded 0.72. (b)
  Regression dilution biases the slope downward but hits both sets equally, and
  the 97%-attenuated figure is dilution-independent - the conclusion is robust.
- Outputs: outputs/gates/gate1_effect_size_check.csv + .pdf;
  gate1_divergence_timing.rds$effect_size_{check,class,verdict}.

Framing consequence:
- The "developmental change licenses Myc" (permission) verb is NOT supported,
  and neither is a "progressive divergence" trajectory. Data show early-max
  divergence that then attenuates ~2-fold (biologically) across the window.
- Convergence IS now supported at the EFFECT-SIZE level (quantified: genotype
  effect ~halves 6W->12W, direction preserved) - not merely a significant-gene
  count artifact. State it as "~2-fold attenuation, not to zero".
- Title verb: move off both "licenses" and "drives progressive"; an
  attenuation/convergence verb is now defensible. The MECHANISM behind the
  attenuation (H1 vs H3 vs H4) is still OPEN, pending Day 2+. No committed
  mechanistic framing yet.

## Candidate trajectory hypotheses (OPEN - to test Day 2+, none adopted)

The author has no committed hypothesis yet; these are candidates, each tied to a
planned analysis that can discriminate it.

- H1 Front-loaded footprint / substrate shift: Myc's impact is maximal at 6W and
  attenuates as the tissue substrate matures; stable dose, moving substrate.
  Test: AP-retention / AP-abund (script 21), per-pathway interaction direction
  (script 19).
- H2 Power artifact: LARGELY REJECTED (script 13 PART 2b, 2026-07-05). The
  count drop is NOT primarily a 12W-power artifact: genotype |LFC| point
  estimates (SE-independent) attenuate to ~45% on a fixed 6W-divergent set
  (slope 0.46), and the LFC-independent expressed set gives the SAME slope
  (0.45), ruling out a selection/power explanation. Residual power effects may
  still inflate the raw COUNT ratio, but the underlying effect-size collapse is
  real. Not a live trajectory hypothesis; retained only as a QC footnote.
- H3 Developmental catch-up (convergence from the WT side): Myc- controls move
  toward the Myc+ state by 12W (WT development advancing), closing the gap rather
  than Myc+ reverting. Test: WT temporal shift AP1 (mitoPPS + MG_* GSVA on Myc-
  6W->12W, scripts 15/18); compare timepoint_neg vs timepoint_pos contrasts.
- H4 Selection / clonal dynamics: early Myc+ divergence carries a stress /
  apoptotic component (see Gate 2 PRO shift) that is selected away by 12W,
  leaving a converged survivor state. Test: cell-death both branches (16/12);
  AP8 cross-sample CV if the selection arm survives.

These are not mutually exclusive (esp. H1+H3). Day 2+ should aim to weight them,
not just pick one.

## Gate 2 - apoptosis readout (script 14)

Facts (raw / Wald interaction LFC; T3 resolved - see flags):
- Apoptosis-PRO (23 genes present): mean interaction LFC -0.130, median -0.080,
  one-sample t p=0.0227 -> coordinately shifted negative (pro-apoptotic Myc
  drive fades at 12W).
- Apoptosis-ANTI (7 genes): mean -0.076, t p=0.25 -> not shifted.
- Bbc3/PUMA: interaction LFC -0.541, rank 2/23, |LFC| > 2x PRO median -> outlier,
  but the PRO group is significant even with Bbc3's influence.
- Class: DECISION_POINT - a module-wide coordinated shift, not an isolated Bbc3
  event.

Reading: supports mitochondrion-as-decision-point (principle 5). The
pro-apoptotic module's coordinated fade is consistent with EITHER coordinated
transcriptional remodelling OR selection against apoptosis-prone cells; Gate 2
alone cannot separate these (Bbc3 is both a group member and an outlier).

Caveat: small n (23 / 7). The PRO t-test rests on a modest set; treat p=0.023 as
suggestive. Corroborate with the two cell-death branches (Day 2, scripts 16/12)
before headlining.

Framing consequence:
- Apoptosis arm STAYS in the main-figure candidate set; AP4 retained. Cell death
  is NOT demoted to Supp 1 (the READOUT branch did not fire).
- Selection-vs-remodelling remains OPEN; carry it into AP-CD and AP8.

## APs / framing surviving Day 1

- Retained/promoted: AP4 (Felsher/Myc co-variation), AP-CD (both branches), the
  decision-point apoptosis arm, AP5 divergence-timing (now a SHRINKS story).
- Reframed / OPEN: trajectory and title verb (permission and
  progressive-divergence both out); depends on H1-H4.
- Unchanged in plan, now also tasked to explain the 6W-max / 12W-converged
  shape: AP1/AP2 dev composition, AP6 preferential alteration, AP7 MB-fork
  (centrepiece).

## Figure lock (provisional; revisit after Day 3)

- Gate 1 SHRINKS panel (6W vs 12W genotype divergence + interaction direction):
  strong Fig 2 / divergence-timing candidate.
- Gate 2 PRO/ANTI forest: Supp 1 candidate; promote to main only if the
  cell-death branches converge.
- Full figure lock deferred until AP7 / AP2 land (Day 3) - the trajectory
  reframing (H1-H4) drives the narrative spine.

## Flags recorded

- T3 RESOLVED: Gate 2 direction + PRO/ANTI already computed on raw / Wald values
  (interaction_results$interaction_raw, unshrunken MLE; IHW touches only padj).
  No shrinkage recompute anywhere in the gate layer. Do not re-litigate.
- Gate 1 power caveat (composite 12W contrast) - RESOLVED 2026-07-05 by the
  effect-size check (script 13 PART 2b): the SHRINKS is biological ~2-fold
  attenuation (slope ~0.45 on both a biased and an unbiased gene set), not a
  power artifact. See outputs/gates/gate1_effect_size_check.{csv,pdf}.
  "Convergence" wording is now defensible if quantified as ~2-fold attenuation
  (not to zero).
- Gate 2 small-n caveat (23 / 7 genes) - corroborate with cell-death branches.
- 01 gene-set tidy still deferred to Block B; trunk self-containment for the raw
  cell-death inputs still a Block B task.

## Decision: STOP at end of Day 1

Do not proceed to Day 2 (scripts 15 GSVA, 16 cell-death port). The trajectory
reframing needs Day 2+ results (H1-H4). Author to review this note and the
candidate hypotheses before we build the GSVA engine and cell-death port.
