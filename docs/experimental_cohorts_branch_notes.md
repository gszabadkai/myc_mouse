---
date: 2026-09-07
tags: [project/myc_mouse, experimental-cohorts, branch-notes, standing-rules]
status: live -- the working notes for the `experimental-cohorts` branch; update in place
relates-to:
  - data/fatpad_timeline/README.md
  - docs/2026-09-07_fatpad_timeline_oxphos_precheck.md
  - docs/2026-09-07_fatpad_tumour_limb_oxphos_trend.md
  - docs/2026-09-08_orthotopic_vector_series.md
  - docs/handoff.md
---

# `experimental-cohorts` — branch notes

Holds the two mouse datasets that are **not** part of the 6W/12W purified-MEC pipeline:
the orthotopic Bcl-xL / PGC1a tumour series, and the fat-pad progression timeline.

## Disposition

This branch **merges back into `paper-final`**. The orthotopic dataset produces main-panel
work, so this is a staging area, not a parallel line. Build accordingly, from the first
commit:

- **Numbered scripts continue `paper-final`'s sequence.** `49` and `50` are taken
  (`49_fatpad_tumour_limb_oxphos_trend.R`, `50_orthotopic_vector_series_scoring.R`); next is `51`. No separate numbering scheme
  to reconcile at merge.
- **Pinned data and provenance READMEs follow the existing conventions exactly**, since
  they arrive in `paper-final` unchanged: large matrices gitignored with a tracked
  `README.md` beside them (`data/FULL.DAT.csv` is the model), MD5 and source path recorded,
  traps named.
- **Scratch stays under `sandbox/`, which is gitignored**, and is deleted before merge.
  Scratch never merges whatever it concludes; only its dated note does.
- **Keep the branch rebased on `paper-final`** rather than letting it diverge. `paper-final`
  is still moving (the text cut has not happened), so rebase before each significant piece
  of work.

## THE STANDING RULE — scores are not comparable across cohorts, within a species

`ox_rel`, `ox_lvl`, GSVA scores, singscore, mitoPPS and every other z-score or
cohort-ranked construction are **cohort-relative**. Their zero and their unit are set by
the samples that went into the scoring run, and by nothing else.

**Therefore: the 6W/12W MEC cohort, the fat-pad timeline and the orthotopic series each get
their own scoring run, and their values are never pooled, never plotted on one axis, and
never differenced.** Not even when the recipe is byte-identical — an identical recipe over
different samples produces a different scale, which is exactly what makes the values look
comparable when they are not.

What **may** travel between cohorts:

- correlations (`rho(ox_rel, Adipoq)` in one cohort against the same in another),
- orderings and directions (which group is highest; which way a slope points),
- effect sizes expressed in each cohort's **own** SD, with the SD named,
- model *form* and its coefficients' signs (the gate's `bMX > 0`), never their magnitudes.

What may **not**: a score, a group mean of a score, a difference of two scores, a "shared"
axis, or any panel with two cohorts' scores on one scale.

Written down now, before the second dataset arrives, because after the merge all three
scoring scripts sit in one repo and the temptation to put their scores on one axis will be
immediate. `mitoPPS` carries the same prohibition for a second, independent reason already
on record: its pairwise-ratio baseline is composition-dependent.

**Corollary 4 — when the INTERVENTION moves the composition markers, a pooled loading can
look material while the per-group loadings disagree.** In the orthotopic vector series
`adipose` loads +0.62 on `ox_rel` across the 21, which trips the material threshold — but the
per-arm values are +0.04 / −0.39 / +0.71, and over all 67 it falls to **+0.09**. The construct
moves `ox_rel` by +0.915, so anything that also differs by arm correlates with it by
construction. **Check the per-group breakdown before believing a pooled confound**, and
remember that adjusting a group contrast for composites that themselves differ by group is
over-adjustment unless a positive control says otherwise.

**Corollary 3 — on this dataset, adjustment cannot separate confound from signal, and the
negative control is how you find that out.** Script 49 PART G residualised the endpoint on
four adipose markers and re-tested. `Mki67`, which genuinely rises across the limb
(tau +0.313, p 0.045), went to **−0.096 (p 0.715)** after the same adjustment — while the
markers explained only 15% of its variance. On this limb adipose depletion and tumour
progression are the same variable, so covarying out composition removes the biology.
**Always carry a positive-control endpoint through any composition adjustment**, and read
it before the endpoint you care about.

**Corollary 2 — a loading is a property of a SUBSET, not of a ruler.** An admissibility
argument ("the confound runs against my prediction, so a positive result is credible") must
be computed on the **exact subset the test runs on**. Script 49 assumed loadings measured
over all 30 samples and they did not transfer: `rho(ox_sub, Adipoq)` is +0.538 over all 30
and **+0.144** within the 12W-to-tumour limb, and the mitoribosome's loading flips sign
between the two. Compute the loading on the analysis subset before relying on its direction.

**FLAG FOR THE HUMAN ARM'S TRAP LIST.** This generalises past this branch. A human cohort
stratified by PAM50, TP53 status or purity is exactly the situation in which a loading
measured on the full cohort stops holding in the stratum being tested — and the human plan
pre-specifies stratified analyses. Every admissibility or specificity argument there must be
recomputed inside each stratum, and any composition/purity adjustment must carry a positive
control through it (Corollary 3). Tell that session; it does not read this branch.

**Corollary that has already bitten twice.** A cohort-relative *construction* also carries
cohort-specific *behaviour*. The compartment-relative share `ox_rel` cancels a confound
only where that confound loads equally on numerator and denominator; in fat pad it loads
more on the respiratory arm than on the rest of MitoCarta, so the share **amplifies** it
(`docs/2026-09-07_fatpad_timeline_oxphos_precheck.md` section 2). Re-check the ruler in
each cohort; do not inherit its properties along with its recipe.

The sharper case is `ox_nuc_mtrib` (script 49): its two halves load on adipose with
**opposite signs** (+0.14 nuclear OXPHOS, −0.24 mitoribosome), so the difference **adds**
their adipose components and the endpoint becomes the most confounded construction on the
limb (+0.58 to +0.78 across four adipose markers) — worse than either half and worse than
`ox_rel`. A difference of composites cancels a confound only when it loads on both with the
same sign and similar magnitude.

## Dataset status

| dataset | pinned | verdict | contributes |
|---|---|---|---|
| fat-pad timeline (Chandan, 30 samples) | `data/fatpad_timeline/` | **CLOSED 2026-09-07.** FAIL for respiratory-axis work, twice: within-tumour correlation (pre-check) and the 12W-to-tumour ordered trend (script 49, T1 FAIL / T3 VOID), and the adjustment that would have rescued it fails its own negative control | the 6W genotype panel (`PUMA:BCL-XL` +0.77 log2, p = 0.016 -- but **+0.65, p = 0.071** once two epithelium-poor `6WK_POS` samples are dropped; cite with the caveat) and `Bcl2l1` flat across progression while `Myc` rises ~11x (unaffected) |
| orthotopic vector series (EV/BclxL/Pgc1a, n = 21 of 67) | `data/orthotopic_series/` | **C1 HOLDS** (`Bcl2l1` +121.8 CPM, p 0.00029, survives adjustment with the positive control intact) and is **SELECTIVE** -- `Mcl1` flat, `Bbc3` falls (script 51); **C2 holds as an observation but is NORMALISATION-DEPENDENT** -- p 0.0029 on CPM, **0.0385 on median-of-ratios**, quote the latter; **C3 UNINFORMATIVE**, for two measured reasons and neither is power | the causal reading of the human `BCL2L1` result, and the in vivo MYC-ceiling observation. **The directed-edge argument is BLOCKED** pending one fact -- see the specificity note section 3 |

**Open, and now the branch's leading question:** does relative OXPHOS share rise again as
tumours establish? The human arm requires a reversal somewhere. Script 49 could not test it
in fat pad — the endpoint is adipose-determined (76% of its variance), a null was
pre-declared uninterpretable, and adjustment removes the progression along with the
confound. **It goes to the orthotopic series**, where the tumour is not embedded in fat pad
— subject to that dataset's own composition check first, computed on the exact subset the
test will use, and carrying a positive control through any adjustment.

Open, carried from the pre-check and **not** resolved: `Myc` rises 7.5x (medians) between
`6WK_POS` and `12WK_POS` in fat pad, against the MEC dataset's stable-dose premise. Three
readings remain live; see the pre-check note section 4. Do not import either dataset's
`Myc` scale into the other.

---

## Two corollaries from the orthotopic series (script 51, 2026-09-08)

**1. The gate's pooled |rho| threshold has a blind spot in BOTH directions.**
The 09-08 note's section 2 showed the first: **pooled-material with the
per-group loadings disagreeing** -- `adipose` +0.621 pooled against
+0.036 / -0.393 / +0.714 by arm, because the intervention itself moves the
composition markers. `ox_lvl` against proliferation shows the mirror: **pooled
-0.323, below threshold, with EV / BclxL / Pgc1a at -0.643 / -0.679 / -0.964** --
three arms agreeing strongly and the pooled value not clearing the bar.
**Read the per-group breakdown before believing either verdict.** A pooled
loading can be material where nothing is consistent, and immaterial where
everything is.

**2. The same lesson on the ENDPOINT rather than the covariate: a near-zero
pooled correlation across designed arms can be two opposed components
cancelling, not an absence.** In the vector series the `Bcl2l1`-vs-`ox_rel`
pooled rho is -0.065, made of a between-arm +0.500 and a within-arm mean of
-0.381; for `log2(Mcl1/Bcl2l1)` it is +0.064 from between-arm -0.500 and
within-arm +0.524. EV and BclxL sit at the SAME point on the respiratory axis
(`ox_rel` -0.292 against -0.291) with `Bcl2l1` 2.1 log2 apart -- a vertical
displacement at constant `ox_rel`, which is the construct.
**Decompose into between-arm and within-arm before reading a pooled rho as a
null.** And the corollary of the corollary: **dropping an arm does not fix it.**
On EV + Pgc1a alone the same rho is +0.657 with an interval excluding zero, but
with two arms and a designed `ox_rel` gap that correlation IS the group
contrast -- circular, and opposite in sign to the within-arm reading.
