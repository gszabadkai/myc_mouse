# What survives the global shift — the author's three pushbacks, answered

**Date:** 2026-07-18. **Context:** after script 36 established that a single broad whole-sample axis
("the global shift") inflates every per-sample pathway *correlation*, the author pushed back on the
framing and on three specific results. This doc records that exchange and the resulting
recalibration. It is a companion to `2026-07-18_five_questions_narrative_guidance.md` and corrects
one claim in that doc's section 4 (the mtDNA/nuDNA ratio).

---

## The framing correction (accepted)

The global shift is a **caution scoped to one class of analysis** — per-sample correlations
*between* pathway scores — **not** a headline finding, and not a caveat on the backbone of the
paper. Most load-bearing results are **contrasts** (fold-changes between groups) or **ratios /
normalised values**, and those do not ride the shift. The earlier summary over-weighted the shift
by letting it colour everything; it should be scoped, not foregrounded.

---

## The three points

### Point 2 — the OXPHOS decline over the timeline. **Author is right. Safe.**

It is a DESeq2 contrast (a fold-change with a standard error), not a per-sample correlation. The
ceiling is a property of cross-sectional *correlations*; a fold-change is immune to it. The decline
is **72-90% Myc-specific** (a genotype x time effect), and the genotype dimension is the one the
shift cannot touch (the axis is genotype-independent, p=0.49). **Keep it as a quantified fact.**

### Point 3 — mitoPPS values. **Author is right, with one border.**

mitoPPS is pairwise-ratio-based on linear counts, so a whole-transcriptome shift cancels by
construction. mitoPPS **reallocation among nuclear pathways** (the Issue #4 reprioritisation) is
background-robust. **The one exception:** any mitoPPS measure that includes the **mtDNA-encoded
arm**. mitoPPS cancels a *global* shift but not variation in the *mitochondrial fraction itself*
(mt% swings 3.4->40% across samples); subtracting the mtDNA arm imports that variance. That
exception is exactly the mitonuclear imbalance -> point 1.

### Point 1 — the mtDNA:nuDNA ratio should not be dismissed. **Principle right; this ratio, tested, fails. Earlier "artefact of arithmetic" wording withdrawn.**

A ratio *can* be the correct correction: if numerator and denominator share the background
multiplicatively, the ratio cancels it. But this particular ratio, tested, is **dominated by the
background, not freed from it.** The imbalance is built (`scripts/33_mtdna_axis_and_coupling_null.R:261`)
as `mean(nuclear-OXPHOS mitoPPS) - mtDNA mitoPPS`, and:

| imbalance vs… | rho | r² |
|---|---|---|
| **mt% (mitochondrial fraction)** | **-0.95** | **0.73** |
| contamination | -0.83 | 0.28 |
| IEG / prep-stress | -0.68 | 0.27 |

Its own contrast test (`imbalance_own_stats`):

| term | beta | p |
|---|---|---|
| timepoint (12W vs 6W) | -0.59 | 0.31 |
| genotype (Myc vs WT) | +0.37 | 0.52 |
| interaction | -0.09 | 0.92 |
| **timepoint, adjusted for the IEG axis** | **-0.01** | **0.985** |

Group means: `6W_neg 0.13 · 6W_pos 0.50 · 12W_neg -0.46 · 12W_pos -0.17`.

So the apparent 6W->12W *decline* in the imbalance is **entirely the prep axis** — adjust for it and
the time slope collapses from -0.59 to -0.01 (p=0.99). The ratio did not cancel the background
because the nuclear and mtDNA arms load in **opposite** directions on the contamination/prep axis,
so the *difference* amplifies that axis instead of cancelling it. The subtraction removes the global
shift but keeps — doubles down on — the mt-fraction axis, which is the one unresolved
technical-vs-biological variable.

**Consequence for the narrative.** The imbalance is not "dismissed for being a ratio." It is tested,
does not reach significance even raw (n=6/group), and what trend it has is prep, not time. The one
axis-independent piece is the **genotype offset** (Myc higher at both timepoints, +0.37) — but it is
ns here. **mtDNA qPCR** measures the physical ratio directly, cancels all of this, and would test
that genotype offset cleanly. Honest positive: the imbalance is a **live hypothesis for qPCR**, not
a refuted one — but it is **not readable from this RNA-seq**, in either the raw or the ratio form.

---

## Inventory — what does and does not depend on the global shift

**Independent (the backbone — contrasts, ranks, and true ratios):**

| result | why it's clean |
|---|---|
| Myc raises mitochondrial **content +21-27%** | genotype contrast on read-*shares*; adjustment-robust |
| Myc raises **absolute nuclear OXPHOS** (d=1.27/1.15) | genotype contrast |
| **OXPHOS declines 6W->12W** (72-90% Myc-specific) | DESeq2 LFC contrast — *point 2* |
| The **attenuation** (66% fade + 34% convergence) | genotype x time contrast (Issue #4/#6 identity) |
| **mitoPPS reallocation** among nuclear pathways | pairwise-ratio, cancels the global shift — *point 3* |
| Myc **amplifies** the pubertal programme (d=2.3-3.0); **BMYO suppression** | genotype main effects |
| Myc **bends off the WT developmental trajectory** | correlation of *gene* LFCs, not a per-sample pathway coupling |
| **fGSEA** enrichments per contrast (incl. Tang RCD modalities) | ranks on the Wald statistic of a contrast |
| Developmental-trajectory validation (4-method floor) | contrasts |

**Dependent on the global shift (the exploratory coupling layer only):**

| result | the exposure |
|---|---|
| "OXPHOS is central" / biogenesis / metabolic-axis **couplings** | per-sample GSVA correlations — scripts 28/35/36 |
| The **mitonuclear imbalance as a *time* claim** | dominated by mt%/IEG axis (-0.59 -> -0.01) — *point 1* |
| **mt%** as a standalone readout | entangled with prep; technical-vs-biological unresolved |
| The **mito->death coupling** | circular + rides the axis (the script-36 negative control) |
| Any single-timepoint cross-sectional per-sample correlation | the ceiling |

**The pattern:** fold-changes, rank-enrichments, and content shares are safe; correlations between
two pathway scores are not. Points 2 and 3 sit in the safe column. Point 1 is the one measure that
*looks* like it should be safe (it is a ratio) but empirically is not — because the thing it fails
to cancel is the very axis we cannot yet resolve.

---

## Open for the author's proposal

1. How much narrative weight the mitonuclear story should carry, given it rests on **mtDNA qPCR**
   rather than the RNA (the RNA gives a genotype offset, ns at n=6, and no time trend once prep is
   removed).
2. Whether the **redox -> MB2 fork** lead (the one coupling that survives the correction) earns a
   place alongside the content/reallocation story.
3. The positive spine the author wants to build — expected to lean on the backbone table above
   (content, nuclear OXPHOS, the OXPHOS timeline, mitoPPS reallocation, the developmental bending),
   with the couplings and the imbalance framed as qPCR/bench-testable hypotheses rather than
   results.
