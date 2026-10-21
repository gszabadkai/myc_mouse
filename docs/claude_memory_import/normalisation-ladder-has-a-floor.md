---
name: normalisation-ladder-has-a-floor
description: "Instrument choice for the figures - every quantifier divides by something, the message needs rungs 0-2, but below rung 2 a large gene set stops being arm-specific (script 46)"
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-08-16T22:56:12.083Z
---

**Instruments were chosen on "which can be best quantified" — which optimises for DEFENSIBILITY,
and defensibility is bought with NORMALISATION. Clarity is bought with the opposite thing: units the
reader already owns.** The two criteria pull in a predictable direction, and the corpus sits on one
side of it. The ladder (`docs/2026-08-17_instrument_choice_for_the_reduced_figure_set.md`):

| rung | instrument | divided by | destroys |
|---|---|---|---|
| 0 | summed normalised counts | library size | nothing beyond depth |
| 1 | % share | the transcriptome | magnitude |
| **2** | **log2FC** | the other group | level |
| 3 | mitoPPS | the rest of the compartment | content |
| 4 | GSVA | the cohort | absolute position |
| 5 | fGSEA NES | its own ranked list **and** its own null | units, cross-contrast sign |

Correcting this costs no rigour, because **since scripts 21/34/36/40/43/45 the defence lives in the
NULL, not in the normalisation.** A raw fold change with a matched-set null is both more defensible
and more legible than an NES.

**BUT THE LADDER HAS A FLOOR, and `scripts/46_axis_ruler_test.R` (2026-08-17) measured it.** The test:
refit Fig. 2I's interaction with the x axis swapped from OXPHOS-subunit **mitoPPS** (a budget share)
to **absolute summed normalised counts** (a magnitude), because the closing hypothesis says "high
respiratory state" and mitoPPS explicitly cannot say that. Decision rule fixed in advance. It
reproduces script 43 PART B **exactly** (max |dbeta| = 0) and script 45's levels to 3.6e-15 first.

| axis | ruler | beta | p | p_emp | ambient | **resid_frac** |
|---|---|---|---|---|---|---|
| `oxphos_ppd` | mitoPPS | **+0.779** | **0.0052** | 0.080 | 0.49 | **0.550** |
| `oxphos_lvl` | levels | +0.535 | 0.156 | 0.154 | 0.81 | **0.071** |
| `redox_ppd` | mitoPPS | +0.207 | 0.72 | 0.486 | 0.26 | **0.998** |
| `redox_lvl` | levels | +0.362 | 0.353 | 0.229 | 0.72 | 0.207 |

**The DIRECTION survives the ruler change; the SPECIFICITY does not.** The two rulers correlate
**r = +0.928** and the interaction keeps its sign at 69% of its size — but only **7.1%** of the
absolute level's variance is separable from the global expression factor (r = **+0.964**), against
**55.0%** for mitoPPS. **On the absolute ruler the axis is "overall expression", not "the
respiratory chain"**, so it could not have carried an arm-specific claim however the p came out.

**RULE: move down one rung, not two, and check separability before the last step.** Rungs 0-1 belong
where the claim genuinely IS about magnitude (Fig. 1E's compartment share), not where it is about
one arm among others. **KEEP mitoPPS for Fig. 2I** — the reason is now measured, not inherited.

**TEXT CONSEQUENCE, not cosmetic:** the supported claim is respiratory **PRIORITY** (allocation),
not respiratory **CAPACITY**. "High respiratory state" is a rate no ruler here measures — hand it to
the PGC1a perturbation explicitly. **By-product worth keeping: `redox_ppd`'s resid_frac is 0.998**,
almost perfectly orthogonal to the global factor — the cleanest axis in the corpus, and a better
reason for it being the control than the one in the legend.

**Proposed census:** log2FC is the spine; **one** fGSEA panel (1D, not three); mitoPPS only where the
ratio IS the claim (1F, 2I); **GSVA to zero as a named ruler** (fine for "is this higher in Myc+", a
liability for "does this track that"); the share promoted to lead beat 1. Figure 1 currently makes a
reader re-learn the y axis **five times** across seven panels.

Related: [[correlation-ceiling]], [[abundance-ruler-dependence]], [[figure-phase-architecture]],
[[fgsea-takes-one-list-not-two]].
