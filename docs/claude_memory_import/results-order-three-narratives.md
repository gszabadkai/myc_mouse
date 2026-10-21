---
name: results-order-three-narratives
description: "Order D is the FINAL Results order for the attenuation/substrate/death block, written out in docs/2026-08-02_results_narrative_final_order_D.md; it supersedes the A/B/C draft"
metadata:
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-08-02T00:30:34.240Z
---

The Results are written paragraph by paragraph in the Google Doc `MK_Myc_Paper` tab
**Gyorgy Writing**. Paragraphs 1-2 are done (Figure 1). The next block covers attenuation ->
substrate -> death.

**SETTLED 2026-08-02. The order is D, and it is written out in
`docs/2026-08-02_results_narrative_final_order_D.md`** -- two sections as continuous Results
prose, plus the corrections, the seven-item dose list, both answered questions with citations,
the panel table and the trap list. **Read that doc, not the A/B/C one.**

The author did not pick from A/B/C (`docs/2026-08-01_results_narrative_three_orders.md`, kept as
the record of the choice). They wrote a fourth spine instead. **D = C's skeleton with B's hinge**:
section 1 is the oncogene (content -> reprioritisation -> apoptosis demoted on BOTH arms ->
x0.55 rescaling -> dose), section 2 is the tissue (which half of the prolif/death balance moves ->
WT withdraws respiration -> priming balance still flat -> PUMA reverses -> Foxo3 -> MYC-ER ->
prediction). Two author improvements: apoptosis moves INTO section 1, and the mechanism questions
are explicit at the end.

**Why:** the choice governs what the panels must show. Figures were deferred until it was fixed;
they are now commissioned in `figures/panels/PANELS.md` as Fig. 2A-2D and 3A-3D, and **most are
ports of the existing `figures/fig01`-`fig05` / `figS1`-`figS8`**, not new builds. Fig. 1E is
superseded and is not a slot.

**How to apply -- three corrections to the author's own spine, all verified this session:**
**PUMA is NOT deprioritised in the background** (`Bbc3` WT 6->12W +0.06 padj 0.88; falls only in
Myc+, -0.48 padj 0.016, int p 0.0081) -- so OXPHOS is the Myc-INDEPENDENT half and PUMA the
Myc-DEPENDENT half, which is the two-input model and is better than the draft. **The background is
NOT the mirror of Myc's reprioritisation** -- cor = -0.02 across 143 pathways, 50% sign-discordant;
they are *independent*, and the specific sign flips need the matched null (OXPHOS subunits 0th pct,
assembly factors 50th) to be quotable. **"Favours OXPHOS" is true at tier level only** -- the
biggest leaf promotions are one-carbon/amino-acid (glycine cleavage +0.41), and it is the subunits
not the assembly factors.

Also settled and not to be re-derived: the attenuation IS proportional and that is a POSITIVE
result; the WT metabolic rise is CATABOLIC; "collapses" -> **"reverses"**; 12W take-off has only
weak transcriptomic support. Traps: never lead on the DEG count, and the genome-wide interaction
has ZERO genes at FDR 10%.

Related: [[foxo3-puma-pair-and-the-nonlinear-dedosing]], [[dose-vs-competence]],
[[mitochondrial-tradeoff-model]], [[panels-follow-the-narrative]],
[[hypothesis-generating-not-failed-tests]], [[attenuation-competing-hypotheses]],
[[background-vs-myc-is-additive]].
