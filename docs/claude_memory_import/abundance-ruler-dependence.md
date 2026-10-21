---
name: abundance-ruler-dependence
description: The three abundance rulers give different answers for the OXPHOS subunits (+0.061 / +0.201 / +0.226) - a sentence quoting one must name which; plus the tissue-paralog caveat
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-08-15T23:08:37.390Z
---

**A sentence quoting an OXPHOS abundance number must name its ruler.** The same 87 nuclear OXPHOS
subunits on the same contrast give three different answers, and they differ most on precisely the
arm the manuscript's closing claim rests on:

| ruler | 6W_wt>12W_myc |
|---|---|
| unweighted per-gene mean log2FC (script 40's content ruler, the figures' default) | **+0.061** |
| expression-weighted mean | **+0.201** |
| log2 of summed DESeq2-normalised counts (the compartment ruler) | **+0.226** |

The gap IS the weighting — a mean of log ratios versus a log of the sum — and script 45 PART H
asserts the ordering (weighted lands between the other two) so a membership or reconciliation
difference cannot masquerade as it. **Why it matters:** "the respiratory chain returns to baseline"
is a *per-gene* statement. It is much weaker on the abundance-weighted ruler.

**Two further qualifications on that sentence, both load-bearing:**

1. **"Baseline" means relative to the mitochondrial COMPARTMENT, not the transcriptome.** Against
   2000 expression-matched sets the OXPHOS subunits still *rise* on the diagonal, at the **95.4th
   percentile**. What survives is *"the arm that gains least — the only one that does not clear its
   compartment"* (every other arm ≥ 99.6th). Never write a bare "no change"; and a percentile is
   not an effect size (PROLIF pooled hits the 100th on +0.057, because 731 genes give a tight null).
2. **The five lowest expressers are tissue-restricted paralogs** — Cox6a2/Cox7a1/Cox8b (heart,
   muscle), Cox4i2 (lung), Cox6b2 (testis) — at expression **ranks 83–87 of 87**, falling hard. In
   an enzymatically dissociated MEC prep their signal is plausibly residual non-epithelial tissue,
   so their fall may be **less contamination rather than less respiration**. Dropping all five:
   unweighted mean **+0.061 → +0.105**, gradient +0.39 → +0.30, weighted mean unmoved. They inflate
   the discrepancy without creating it — quote **+0.105** beside +0.061.

**The expression gradient itself is real and specific, and is NOT wild-type-only.** All four
contrasts beat an expression-matched null with the null median at the *opposite* sign (matched
random genes have high expressers moving down; the respiratory chain does the reverse):
`wt_time` +0.330 (100th pct), `cross` +0.388 (100th), **`myc_6W` +0.387 (100th)**, `myc_12W` +0.197
(99.7th). So the maturing gland withdraws disproportionately from the *minor* subunits and holds the
high-abundance core — and **Myc's induction is itself expression-graded**, preferentially amplifying
already-abundant subunits. My first reading ("the gradient comes from the wild-type timeline, not
Myc") was wrong and came from comparing `wt_time` against `myc_12W` only.

**The levels ruler corroborates and slightly enlarges the content claim** — genotype betas agree
with script 32's share model at r = 1.000 but are **larger on 11 of 11** matched panels (median
offset +0.048). Expected direction: the share divides by a denominator that itself rises with Myc.
This puts a number on script 32's own "share is a LOWER BOUND" (see [[mito-content-bounded]]).

Everything above is from `scripts/45_state_readings.R` → `results/state_readings.rds`, drawn in the
three panels described in [[figure-phase-architecture]]. Related: [[intro-oxphos-restraint-tension]],
[[mitocarta-sets-are-membership-loose]], [[dont-assert-a-mechanism-without-checking-it-can-work]].
