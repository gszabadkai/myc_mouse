---
name: intro-oxphos-restraint-tension
description: "Manuscript-intro \"reduced abundance of OXPHOS\" resolved by AXIS and by FRAME: on genotype Myc RAISES OXPHOS; the temporal decline is largely dose; and the WT gland DOES decline (the old \"72-90% Myc-specific / WT barely declines\" gloss was a frame error)"
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-25T23:10:34.841Z
---

The drafted Introduction says MYC leaves "OXPHOS restrained" / "Reduced abundance of OXPHOS
complexes ... associated with desensitisation to MYC-induced apoptosis." Author clarified
(2026-07-13) this means the **Myc+ 6W->12W TIME-axis decline**, NOT Myc-vs-WT.

**Resolved by AXIS** (still stands): on the GENOTYPE axis Myc **RAISES** nuclear-OXPHOS mRNA
(d ~1.25, content +21-27%); on the TIME axis, within Myc+, nuclear-OXPHOS subunit mRNA DECLINES
6->12W. Orthogonal, both true. **Never write that Myc reduces OXPHOS on the genotype axis.**

**CORRECTED 2026-07-26 — the "72-90% Myc-specific / WT OXPHOS barely declines" gloss is WRONG,
and it is a FRAME error, not an arithmetic one.** Issue #6 (`31_attenuation_mechanism.R`,
`$decomp_conv`) decomposes the **attenuation of the genotype GAP** into WT-convergence + Myc-fade.
For `MITOCARTA_OXPHOS_SUBUNITS`: gap 0.674 -> 0.388, atten 0.286 = **-0.126 WT convergence + 0.412
Myc fade** (conv -44%, fade 144%). **That is right: OXPHOS is the arm where the background does NOT
converge toward Myc.** But "does not converge" is about the GAP. Myc+ sits *above* WT on OXPHOS, so
a WT gland that falls is moving *further below* Myc — anti-convergence and a large WT decline are
the same fact from two sides. And the WT decline is large: `MITOCARTA_OXPHOS_SUBUNITS` WT 6->12W
**-0.2548** vs Myc+ **-0.4047** (**63% shared**), fGSEA NES **WT -2.670 vs Myc+ -2.715** (98%).
Script 29's `$decomp` and script 40's `$ruler` agree to six decimals — one measurement, two code
paths. **OXPHOS is the ONLY mitochondrial arm the normal gland withdraws from** (every other arm's
WT NES is near zero or positive: TCA +0.58, lipid +1.49, amino-acid +1.84, mitoribosome -0.61) and
it has the largest shared fraction of all seven Level-1 tiers (0.535 content / 0.517 priority).

**Which frame to use is not taste.** The attenuation frame asks whether the background closes the
genotype gap (Issue #6's question). The **substrate frame** asks what state the gland is IN — and
that is the question the two-input model asks, because MYC-ER holds MYC fixed and samples the
gland's own state. The Introduction needs the substrate frame. See [[dose-vs-competence]].

**Also superseding the 2026-07-13 wording:** the Myc-specific half of the temporal decline is now
just the **x0.55 dose effect** (MYC protein -50% at constant transcript), so "a largely Myc-specific
de-amplification" is an artifact, not a finding. What survives is (a) the WT/developmental
withdrawal above and (b) the **de-prioritisation on the mitoPPS ratio ruler where dose cancels**.
Riders that still hold: it is a NUCLEAR-subunit transcript claim while mtDNA-encoded OXPHOS rises,
so "abundance of assembled complexes" belongs to the protein blots (which now exist); and every WT
temporal statement is batch-confounded ([[batch-equals-timepoint]]), mitigated only by the blots
being off the RNA batch.

**How to apply:** full clause-by-clause verdict and a drafted Introduction are in
`docs/2026-07-26_introduction_alignment_and_the_question.md`, which supersedes sections 3-4 of
`docs/2026-07-13_BlockA_revision_walkthrough_and_intro_alignment.md:884-1066`. Related:
[[issue4-attenuation-scoped]], [[background-vs-myc-is-additive]], [[mito-content-bounded]].
