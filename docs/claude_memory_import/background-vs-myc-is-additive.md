---
name: background-vs-myc-is-additive
description: "Myc's mito effect is NOT shaped by the background: it rescales x0.55 without reshaping, the background moves elsewhere, and Issue #6's global WT-convergence is a shared-baseline artifact that collapses to chance"
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-27T16:28:45.710Z
---

**Script 40 / fig03, 2026-07-24** (`docs/2026-07-24_background_vs_myc_interpretation.md`,
`results/background_vs_myc.rds`). The author asked whether Myc's mitochondrial effect is
shaped by a changing background. It is not. The supportable claim is **additive, not causal**:
a reallocation both genotypes undergo over the timeline, plus a Myc programme that decays
proportionally without changing shape.

1. **Rescaled, not reshaped.** `Myc@12W ~ Myc@6W`: slope **0.55** (boot 0.51-0.60), R2 0.80,
   rho 0.84 (content); 0.64 / 0.79 / 0.93 (priority); per-tier slopes 0.53-0.83. R2 at the
   **100th pct** of overlap-preserving label shuffles. Pathway-level form of Issue #4.
2. **The Myc-specific temporal term is a DOWNWARD SHIFT, near-constant across pathways.**
   `Myc+time ~ WTtime`: slope 0.75 (null 0.63, 82nd pct -- **not distinguishable**), intercept
   **-0.17** (null -0.05, **beyond all 500 shuffles**). The intercept is the batch-clean
   quantity (a common offset cancels). **Four precision points, added 2026-07-27 when the
   author asked for this to be explained to a non-expert -- doc section 2.2 carries them all:**
   - **"A UNIFORM offset" is an approximation**, exact only at slope 1; at 0.75 there is also a
     term `(slope-1)*c_tn`. Say **the intercept is the extreme, batch-clean term**.
   - **The mirror error is worse:** *"the Myc change is not due to the WT change"* is WRONG --
     slope 0.75 says much of the movement IS shared; the point is that the shared part is
     unremarkable while an **additional** downward component is not.
   - **NEVER FLIP THE AXES.** Regression is not symmetric: `c_tp ~ c_tn` = **0.749/-0.172**,
     `c_tn ~ c_tp` = **0.559/+0.118** (the geometric mirror would be 1.335). Flipping also
     destroys the quantity of interest. WT is the substrate, hence the predictor.
   - **The null is a GENE-LABEL SHUFFLE** -- it tests against the data's own structure, NOT
     against sampling new animals. 500 draws, so "0th percentile" = **p < 0.002, not p = 0**.
     The bootstrap CI (0.57-0.94) resamples **pathways, not mice**.
2b. **The plain-language version, and the one to quote.** WT drifts **UP** over the window
   (median +0.041, 103/143 pathways rise); Myc+ drifts **DOWN** (median -0.135, 22/143 rise).
   Per pathway `c_tp - c_tn` is negative in **133 of 143**, median **-0.18**. Only 96/143 are
   larger in ABSOLUTE size, so **"further down" is robust and "bigger" is not**. That downward
   drift IS the attenuation, not Myc doing something extra. `figureS2` panel C plots this
   difference ranked (a flat band = the near-uniformity shown, not inferred) with the
   regression as an inset -- but the difference must **never** be plotted against either
   temporal contrast (artifact ledger below).
3. **The background is not moving toward Myc.** Projection of WT-time on Myc(6W) = **0.064**
   vs null 0.24 (**9.8th pct** -- less than an arbitrary matched panel); cos **-0.34/-0.39**
   in the biogenesis+OXPHOS arms (moves AWAY). `cos(Myc+time, Myc@6W) = -0.84` (0th pct):
   Myc+ retraces its own programme backwards.
4. **Issue #6's global WT-convergence is an ARTIFACT.** `myc_6W = pos6 - neg6` and
   `timepoint_neg = neg12 - neg6` share the 6W WT baseline with opposite signs (cov =
   +var(neg6)), manufacturing convergence. **The comparison must be MATCHED**: splitting
   halves the baseline, and the extra noise in `d = sign(m6)` pulls `frac_wt_toward` toward
   0.5 from BOTH sides -- so each of the 20 partitions runs twice off the SAME half baseline,
   shared vs split (identical d, universe, noise). Global artifact **+0.242 +/- 0.138, right
   sign in 19/20**: `frac_wt_toward` **0.720 -> 0.478 = chance**, conv_pct **+61% -> -12%**;
   mammary-luminal +0.234. **SURVIVES:** amino acid 0.74->0.65, lipid 0.77->0.63,
   TCA/nucleotide ~0.63. **Artifact ~0 on the DIVERGENT arms** (MYC-core +0.011, OXPHOS
   subunits -0.010) -> their published 0.27 stands; the 0.27->0.37 published-vs-halfshared
   gap there is baseline NOISE, not a correction. Programme-specificity was the whole result;
   the global number was the artifact. Script 31 unchanged on disk.
4b. **Per-tier geometry** (`$geometry`): "biogenesis+OXPHOS" is NOT homogeneous -- protein
   import (-0.63 content / -0.77 priority) and OXPHOS (-0.65) move AGAINST Myc; central dogma
   (+0.74) and metabolism (+0.48) move WITH it. cos(Myc+time, Myc@6W) is -0.80..-0.98 in every
   substantial tier.

**THE ARTIFACT LEDGER (never cite these).** cor(WTtime, delta-Myc-effect) = -0.185 content /
-0.327 priority -- and the value the marginals FORCE is *exactly* the same (delta = tpos-tneg);
cor(retention, WTtime) = -0.52 (retention = 1 + (tpos-tneg)/m6); cor(WT baseline, Myc effect)
= -0.79/-0.88 (mitoPPS is compositional). Carried in `$artifact_ledger` with verdicts.

**Ruler tension:** same gene sets, the share ruler retains **77-83%** of the genotype gap at
12W, the LFC rulers **54-63%**. Sign agrees, size is normalisation-dependent -> quote the range.

**Why:** the manuscript needed a mechanism for the attenuation, and the one on the books
("the WT background matures toward Myc's state") rested on a statistic that shares a baseline
with the contrast it is compared to. Sharpening it made the result stronger, not weaker.

**How to apply:** whenever a claim couples two DESeq contrasts, check first whether they share
a group mean -- if they do, split the shared group before believing the coupling. And compute
the value the algebra forces before quoting a correlation between derived quantities. Related:
[[issue4-attenuation-scoped]], [[correlation-ceiling]],
[[sample-level-couplings-need-a-null]], [[dont-assert-a-mechanism-without-checking-it-can-work]],
[[batch-equals-timepoint]], [[figure-phase-architecture]].
