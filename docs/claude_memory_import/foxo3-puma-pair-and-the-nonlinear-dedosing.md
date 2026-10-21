---
name: foxo3-puma-pair-and-the-nonlinear-dedosing
description: "Foxo3 sits at the 0.456th percentile of the departure-from-dose scan next to Bbc3 at the 0.513th; PUMA suppresses MPC in the literature so the naive arrow is backwards, and the de-dosing non-linearity is Lorenzin affinity-ordered occupancy"
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-08-02T00:30:59.164Z
---

Two things established 2026-08-02 that will otherwise be re-derived from scratch, and one
literature reconciliation the author asked for.

**1. FOXO3 IS THE ANSWER TO "is it Foxo3 driven?", AND IT WAS ALREADY ON DISK.** In script 44's
departure-from-the-x0.55-line scan over 8774 genes (`collapse_module_ownership.rds$collapse_genes`,
`resid_aligned = sign(LFC6) * (LFC12 - RATE*LFC6)`), `Foxo3` sits at the **0.456th percentile**
(z -2.53, retention -0.995) and `Bbc3` at the **0.513th** (z -2.49, retention -1.095) -- adjacent.
DESeq2 interaction p **0.0074** and **0.0081**, the two lowest of any gene in the death/biogenesis
rosters (both padj ~0.84-1.0 genome-wide, so pre-specification is what licenses the test).
**Direction is the informative part: the WT gland RAISES Foxo3 with age (+0.47, padj 0.0078) and
the Myc+ gland does not (-0.01)**, while `Bbc3` is flat in WT (+0.06, padj 0.88) and falls only in
Myc+ (-0.48, padj 0.016). Control that was run and passed: within the 4199 **Myc-induced** genes
alone both hold the same 0.45th/0.52nd percentile, so neither rides the repressive-arm asymmetry
(Myc's repressed arm retains 0.33 vs 0.51 induced among confidently-moved genes; genome-wide there
is no asymmetry, 0.445 vs 0.471). FOXO3 is AKT-regulated post-translationally -> **the western is
the test, message alone cannot settle it.**

**2. THE PUMA -> MPC LITERATURE POINTS THE OTHER WAY, AND THAT IS THE INTERESTING PART.** PUMA
binds the mitochondrial pyruvate carrier and disrupts MPC1-MPC2 hetero-oligomerisation, a
non-apoptotic BH3-independent function (Kim et al., Cancer Cell 2019). So PUMA down should raise
OXPHOS -- the opposite of the co-decline we see. **Our data say MPC is not a PUMA-controlled node
here**: `Mpc2` +0.69 (padj 0.0084) and `Mpc1` +0.56 at 6W, both tracking Myc, both falling with
age in Myc+ -- *with* PUMA. Resolution: respiration is UPSTREAM of the trigger, not downstream.
Anchors: cells lacking OXPHOS with low membrane potential are LESS apoptosis-susceptible
(Dey & Moraes JBC 2000); FOXO3 directly transactivates BBC3/PUMA, ChIP-confirmed, AKT-controlled
(Fitzwalter Dev Cell 2018; Obexer Mol Cancer Ther 2010). The circuit that reconciles both is
**negative feedback: respiration -> FOXO3 -> PUMA -> MOMP, and PUMA -| MPC brakes the state that
induced it.** PGC-1a raising respiration should restore both -- which is the next section's
experiment.

**3. "THE DE-DOSING IS NOT LINEAR" HAS A PUBLISHED MECHANISM.** Myc target promoters differ in
affinity and high-affinity ones are saturated at physiological Myc (Lorenzin et al., eLife 2016),
so halving the protein cannot halve every target equally. Our retentions order exactly as
predicted: **Myc core 0.74-0.77 > global 0.55 > mito arms 0.50-0.58 >> PUMA -0.09.** Affinity
explains everything above the line and **cannot produce a sign reversal** -- which is precisely
why PUMA is the residual worth an experiment. And the reason the genome-wide interaction finds
nothing: median lfcSE **0.333 vs 0.233** (factor 1.43, the root-two of a difference-of-differences)
on a *smaller* effect (0.208 vs 0.257) -> ~3x higher detectable effect; plus **the null is wrong**
-- under a global rescaling every Myc-responsive gene has a non-zero interaction by construction,
so the specific question is departure from the x0.55 line, not difference from zero.

**4. THE BETTER FRAMING, FOUND 2026-08-09 AND NOW IN THE WRITTEN TEXT.** Do not quote the
8774-gene percentile: among script 44's **26 curated biogenesis/cell-death mechanism genes**
(`collapse_module_ownership.rds$mech_genes`, 7 groups: FOXO, BH3-only pre-specified, BH3-only rest,
ISR/ATF4, IMS payload, cristae/OPA1, cardiolipin) **`Foxo3` has the LOWEST interaction p (0.0074)
and `Bbc3` the second (0.0081), with a clear gap to the third (`Pgs1`, 0.066)**. Rank 1 and 2 of 26
in a meaningful universe is a far stronger and more legible statement than "0.46th percentile of
8,774", and it is what the reworded Results now claims. Two companions worth carrying with it:
**`Bmf` rises in BOTH timelines** (+1.02 WT padj 0.0002, +0.84 Myc+ padj 0.003) -- the one death
gene that is purely developmental, and the natural contrast to `Bbc3`; and **`Bax` also falls
significantly on the Myc timeline** (-0.376, padj 0.0077), so "specific to `Bbc3`" holds *among the
BH3-only sensors* and the sentence should say so (Bax is an effector).

**5. THE TWO GENOTYPE CONTRASTS ARE MIRROR IMAGES, AND THAT IS WHY THE INTERACTION IS RANK 1**
(asked and answered 2026-08-09; do not re-derive). Myc RAISES `Foxo3` at six weeks (**+0.243,
p 0.057, padj 0.21**) and LOWERS it at twelve (**-0.242, p 0.059, padj 0.53**) by almost exactly as
much. **Neither half clears 0.05** -- so do not label either one on a panel; their DIFFERENCE is the
strongest term in the roster (-0.486, p 0.0074). On the per-animal z scale the same pair is +0.92
(p 0.076) and -0.89 (p 0.092).

**6. THE OXPHOS -> PUMA:Bcl-xL COUPLING IS AN INTERACTION AND BOTH HALVES ARE REAL.** Within the
Myc+ animals slope **+4.35, R2 0.354, p 0.041**; within the wild-type animals **-4.80, R2 0.328,
p 0.052** -- the wild-type half is as strong as the Myc half and points the OTHER way. The
p = 0.0052 everyone quotes is script 43's *difference* between them. **A figure that draws only the
Myc+ fit reads as a Myc-specific correlation and is misleading**; the interaction p has to travel
with it.

**6a. CORRECTION (2026-09-02, script 48): those two slopes are COVARIATE-FREE, and the wild-type
half does not survive adjustment.** Both fits above omit the epithelial and immune composites that
every other model in scripts 42/43 carries. Refitted with `epi + imm`: the wild-type slope goes to
**-0.01 (p 0.95)** and the Myc+ slope to **+0.87 (p 0.0019)** on the compartment-relative axis
(mitoPPS: **+0.14, p 0.88** and **+2.97, p 0.0036**). So it is a **GATE, not a crossover** --
respiration does NOTHING without Myc rather than protecting against death. Do not write
"respiration is protective in the wild type"; the composition-adjusted data do not support it. The
interaction itself is unharmed and gets *stronger* on adjustment. See [[myc-oxphos-priming-gate-model]].

**How to apply:** these three are the spine of section 2 of the final narrative
(`docs/2026-08-02_results_narrative_final_order_D.md`). Do not re-run the scan to check them; do
not write `Bbc3` as a per-gene significant interaction. If any of the three new numbers is cited
in the manuscript, give it a home in a numbered script rather than the figure layer.

Related: [[results-order-three-narratives]], [[dose-vs-competence]],
[[oxphos-apoptosis-module-verdict]], [[prespecify-from-the-perturbations]],
[[dont-assert-a-mechanism-without-checking-it-can-work]].
