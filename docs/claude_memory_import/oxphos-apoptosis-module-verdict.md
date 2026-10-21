---
name: oxphos-apoptosis-module-verdict
description: "There is no co-regulated OXPHOS/apoptosis module in the mouse data, and the PGC1a axis does not own the apoptotic genes -- MYC's regulon does; the reframe (dependence, not co-transcription) is stronger than the module claim"
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-27T01:23:34.920Z
---

**The author asked whether the title could name "a PGC-1a/Nrf1 regulated common OXPHOS/apoptosis
module". Answer: no, and four checks say so — three in our own data.** Script 44
(`44_collapse_module_and_ownership.R`, `results/collapse_module_ownership.rds`). Numbers below are a
verified scratchpad dry run; **reconcile against the author's own run before they leave the docs.**

**1. OWNERSHIP IS THE OTHER WAY ROUND.** `CORE_MITO` (ESRRA/NRF1/GABPA targets n MitoCarta, n=300,
and it holds 66% of the OXPHOS subunits so it IS a fair test) overlaps the 25 pro-apoptotic MitoCarta
genes in **5 vs 6.6 expected — DEPLETED, fold 0.76, p 0.83** — and those five (**Aifm1, Aifm2, Cycs,
Endog, Ifi27**) are the releasable IMS **payload**, no BH3-only gene. `MYC_MITO` (n=119) is
**ENRICHED 2.68x, p 0.011**: Bax, Casp8, Casp9, Diablo, Pmaip1, Aifm3, Sphk2. Neither contains
`Bbc3`. **PGC1a loads the gun; MYC pulls it.** Set membership in a curated resource — constrains the
TITLE, does not establish regulation.

**2. The shortlist cannot be cited** — see [[dose-vs-competence]]. Base rate, BH FDR 0.99998. The
axis stands on the westerns.

**3. THE WT GLAND DOES NOT MOVE ITS APOPTOTIC TRANSCRIPTS.** Content: OXPHOS subunits **-0.255 (0th
pct)** vs Apoptosis-PRO **+0.025 (63rd)**, ANTI +0.045 (65th). Priority (mitoPPS): **OPPOSITE SIGNS**
— OXPHOS -0.109, PRO **+0.054**. Of 25 PRO + 7 ANTI genes exactly one is significant on
`timepoint_neg`: `Bnip3` +0.635 (padj 0.021), **rising**. Apoptosis is de-prioritised on the **Myc**
axis at 6W (-0.210, padj 0.022), not the developmental one. **So the author's premise — a common
OXPHOS/apoptosis de-prioritisation in WT — does not hold on either ruler.**

**4. THE PGC1a REGULON DOES NOT WITHDRAW AS A UNIT.** `CORE_MITO` +0.013 (75th pct); its 59
OXPHOS-subunit members **-0.216 (0.05th)**, the other 241 **+0.069 at the 99.6th**. Selective, not
global. (`Ppargc1a` baseMean 29.9, WT +0.397 n.s.; `Nrf1` -0.249 n.s.; only `Esrra` moves, on the
GENOTYPE axis +0.684 padj 1.8e-06. A transcript is a poor proxy for a coactivator.)

**THE REFRAME, and it out-novels the module.** *The maturing gland does not dismantle its death
machinery; it withdraws the respiratory state that machinery depends on.* **What is shared is not
transcription but DEPENDENCE.** That is what the Bcl-xL rescue localises (execution, not induction),
and it makes PGC1a co-induction in MYC-ER a prediction rather than a hope.

**THE COLLAPSE SCAN, author-run 2026-07-27. Rank on the residual from the global rate** —
`sign(LFC6) x (LFC12 - rate x LFC6)`, standardised — because under a global x0.55 rescaling every
Myc-induced gene has a negative interaction. **`Bbc3` lands at the 0.51st percentile of 8774
Myc-responsive genes** (retention -1.10, a sign reversal), **pre-specified from the westerns**,
reproducing script 42's PUMA:Bcl-xL result with a statistic that never touches the Bcl-xL
denominator. **`Bcl2l11` lands at the 91.6th — the pair splits.**

**THE ANSWER: PUMA IS A SOLO. No pre-registered question returns a module** — OXPHOS subunits
+1.55 (padj 0.063, if anything *retained*), apoptosis null on every lens, AP_TEB null, ISR null.
What collapses genome-wide is **luminal/hormone-sensing LINEAGE IDENTITY** (`MG_HS_GRAY` NES
**-2.79, padj 1.7e-17**; LHOR/LHS/GarciaSola luminal sets; HS and AP TEB signatures) and what is
retained is **the MYC programme itself** (every MYC signature NES > +2.5) — independently
reproducing Issue #6's "MYC-core protected". 56 sets collapse, 161 retain at BH<0.05.
**So the title is decided: the two-input one, NOT a mechanism-naming module title.**

**CUTS AGAINST THE TIDY CARRIER STORY:** the PGC1a-axis `AP_LE` regulons are among the *best-retained*
sets in the genome (`TFT_GABPA_GRAY_AP_LE` **+2.08, padj 7.3e-06**; `TFT_ESRRA_GRAY_AP_LE` +1.77,
padj 0.008). Myc's grip on the ESRRA/GABPA programme does not weaken. The respiratory withdrawal is
the **gland's**, not a failure of Myc to drive it.

**UNPLANNED INTERNAL CONTROL, and it confirms the premise:** `Myc` itself is at the **100th
percentile, retention 1.09** (+1.679 -> +1.822). **The transgene message does not attenuate**, so the
x0.55 fade downstream cannot be transcript-level dose — the blot's premise ([[dose-vs-competence]])
confirmed from inside the RNA.

**ONE LEAD SURVIVES, ONE DIES.** **`Foxo3`** is the most collapsed gene in the mechanism panel
(0.87th pct) and its profile matches PUMA's: baseMean 3265, **+0.243 at 6W (p 0.057) -> -0.242,
interaction p 0.0074** vs `Bbc3` +0.258 (p 0.064) -> -0.283, p 0.0081. FOXO3 is PUMA's canonical
p53-independent activator and the axis is in `PUMA-and-its-relationships.md` §4 — **though our
direction is OPPOSITE to the canonical one there** (PGC-1a sequesters FOXO3a; losing it *activates*
PUMA). **Caveats must travel: both 6W effects n.s., both interactions genome-wide null, two genes are
not a module, and `Foxo3` RISES in the WT gland (+0.473). A western, not a Results sentence.**
**The ISR route is now TESTED AND NEGATIVE** — `Chac1` looked second-most-collapsed but has no 6W
effect to lose (-0.300, **p 0.47**). Do not propose the retrograde/ISR route again.

**THE PUMA COLLAPSE HAS NO WILD-TYPE DEFINITION.** WT `Bbc3:Bcl2l1` is **+0.168** (drifts UP); the
whole collapse is in the Myc effect (+0.676 at 6W -> -0.044 at 12W). It is **developmental in cause,
Myc-dependent in measurement** — a loss of INDUCIBILITY, not of level. You cannot write it as "the
gland turns down PUMA priming".

**TRAPS THAT BIT DURING THIS WORK, do not repeat them:**
- **25 of the 87 OXPHOS subunits ARE AP LE-marker genes.** Any lineage claim built on the raw sets is
  the OXPHOS claim restated, and adjusting OXPHOS for a raw AP-LE score is self-adjustment (beta
  collapses -0.459 -> -0.044 for arithmetic reasons; the disjoint version gives -0.151). Script 44
  carries mito-removed columns everywhere. See [[le-is-low-expressing]].
- **A padj filter excludes the pre-specified genes** (`Bbc3` 6W padj 0.22, `Bcl2l11` 0.46). Script 44
  reports on a wide set and ranks on a narrow one (`in_ranking` flag). PUMA's finding was never in
  its level.
- **A 2-gene group median describes neither gene** when they sit at opposite ends. Report min/max.

**How to apply:** doc is sections 0.7a / 0.8 / 0.9 of
`docs/2026-07-26_introduction_alignment_and_the_question.md` (module verdict, the pre-registered
results table, mechanism menu). **TITLE DECIDED: "Mitochondria integrate oncogenic and developmental
signals to time Myc-driven mammary tumourigenesis"** — the mechanism-naming option needed a collapse
module and there is none; the dependence reframe goes in the Discussion. **Deferred and wanted: MYAZ + PGC1a RNA-seq**, the only design
that defines the module by perturbation; noted in `paper/myc_mito.qmd`. Related:
[[mitochondrial-tradeoff-model]], [[death-priming-change-untested]],
[[mitocarta-sets-are-membership-loose]], [[prespecify-from-the-perturbations]].
