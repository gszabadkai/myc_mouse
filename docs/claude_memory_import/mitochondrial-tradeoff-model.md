---
name: mitochondrial-tradeoff-model
description: "Respiration is an ASSET and a LIABILITY to a MYC-driven cell (PGC1a kills MYAZ; PGC1a+Bcl-xL grows better than Bcl-xL alone) -- the paper's frame, its question, and its closing claim"
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-27T01:04:17.577Z
---

**The frame that replaced "permissive gate", 2026-07-26.** The MYAZ implants with BOTH comparator
arms: **PGC1a alone KILLS MYAZ cells; PGC1a + Bcl-xL grows BETTER than Bcl-xL alone.** A sign
reversal of the same intervention's effect on fitness, conditional on apoptotic competence. So
mitochondrial respiratory capacity is simultaneously the strongest **pro-growth** and the strongest
**pro-death** input, and which limb is expressed depends only on whether death can be executed.

**THE QUESTION (use this form):** *Mitochondrial respiratory capacity is both an asset and a
liability to a MYC-driven cell: it supports growth, and it primes death. When, and how, does an
incipient tumour capture the asset without paying the liability?* The older "what licenses a
MYC-driven cell to survive its own oncogene" is **subsumed** — it asked only about the liability
limb and could not explain why respiration comes back for metastasis.

**THREE SOLUTIONS, and the paper has all three:** (1) **development gives it away free** — the
maturing gland withdraws the respiratory arm and apoptotic competence with it (the window,
transient, no mutation); (2) **selection finds it expensively** — SS and passaged PGC1a cells shed
the biogenesis and lose the ASSET with the liability; (3) **buffering solves it properly** — block
execution and respiration becomes pure asset (Bcl-xL rescue; PGC1a+Bcl-xL implants).

**THE CLOSING CLAIM (unifies initiation and metastasis into one mechanism):** *what changes across
progression is not what the tumour NEEDS from its mitochondria, but what it can AFFORD.* Early it
cannot afford respiration and development spares it the cost; late, once death is buffered, it can,
and respiration is re-engaged for dissemination. **Human test, cheap and in hand: do the
mitochondria-high METABRIC biclusters carry higher BCL2/BCL2L1/MCL1?**

**Q3 ANSWERED, and it CONTRADICTS the proliferation hypothesis (this is good for the model).** WT
6->12W set-mean raw LFC: **OXPHOS subunits -0.255** vs **PROLIF_* pooled -0.047** (5.4x), nucleotide
metabolism **+0.0002**, mitoribosome **-0.024**, **OXPHOS assembly factors +0.001** (the internal
control -- same complexes, no movement), TEB-vs-ductal **-0.422**; amino-acid +0.184 and lipid
+0.122 RISE. Matched-random-set null: OXPHOS subunits 0th pct, PROLIF 1st, nucleotide/mitoribosome/
assembly n.s. **Paired null on the DIFFERENCE (the actual claim): OXPHOS-minus-PROLIF 0th pct.**
Markers agree (Mki67 -0.04, Ccnb1 +0.15, Plk1 +0.11, Aurka +0.31; only Top2a/E2f1 fall). **The gland
DE-RESPIRES WITHOUT DE-PROLIFERATING.** Physiological reading is **morphogenetic**: it loses the TEB,
and the library already records `HS_TEB` as "apoptosis-loaded (lumen clearance by apoptosis)" -- the
pubertal gland keeps a death-ready mitochondrion because duct cavitation needs apoptosis. **Had
respiration fallen only because proliferation fell, the two-input model would collapse into "less
proliferation, less death" and MYC-ER's 80%-vs-30% would have no explanation.**

**Q1 -- what else respiration buys.** Growth (the MYAZ limb). In vivo proliferation and respiration
are DISSOCIABLE (above); the cross-sectional coupling is **untestable**, not weak (script 36: OXPHOS
IS the global axis, r=0.968) -- never quote a coupling there. In mitoPPS ratio space (script 43 PART
B, standardised outcomes so rows are comparable -- NOT comparable with script 42's +2.79) **ONLY ONE
OUTCOME HAS A CONTROL THAT BEHAVES:** `myc x OXPHOS-priority` on the **gene-level PUMA ratio**
**+6.09, p 0.005, perm 91.7th**, `redox` control **51.2nd = dead centre**. `Bax:Bcl2l1` (67.7 vs
control 65.4) and the **mitoPPS priming COMPOSITE** (9.3 vs control 11.5) are **NOT readable** --
the control does what the axis does. Proliferation is **negative** (-3.68/-3.23, perm 6.3/9.2) with
controls shifted the same way but less far (30/28) = a RANKING observation. **The defensible claim
is a SIGN CONTRAST:** on the one outcome with a working control, respiratory priority makes MYC
better at raising the death trigger, and on no outcome does it make MYC better at raising
proliferation. **Distrust the `with_tp` column** -- adding tp*myc makes the proliferation and
composite terms STRONGER (p 0.107->0.031, 0.060->0.011) while the PUMA term weakens (0.005->0.088);
coefficients that grow on adding a collinear term are collinearity, and only PUMA moves honestly.
And the in-vivo proliferation negative is what the model PREDICTS -- the asset limb is only visible
when death is blocked, and the high-respiratory cells are the ones dying and missing from the
library (survivor bias). Third limb, hypothesis only: `redox -> MB2 fork -0.63` = a ROS/mutagenesis
arm.

**AMBIENT CORRECTION, and it is the error CLAUDE.md warns about.** I quoted ~0.42 as the mitoPPS
ambient for OXPHOS couplings. 0.42 is script 38's **global median across pathways**; the ambient
**for oxphos_ppd specifically is 0.493** (q90 0.838), vs redox_ppd 0.260. **Never borrow one axis's
ambient for another.** No conclusion changes (the permutation null is axis-specific by
construction), but quote 0.49.

**Q2 -- the second hit is a STATE, not a lesion.** The field's answer to "MYC alone is not enough"
is always a lesion that disables apoptosis (Eu-Myc + p53/ARF loss or BCL-2; beta-cell MycER +
Bcl-xL -- **check both citations**, the second is the direct ancestor of this paper's Bcl-xL
rescue). Here normal development supplies a **transient, mutation-free phenocopy**. That is why it
is a *window* and not a switch, and it predicts non-stochastic timing.

**The gland does NOT buffer:** Bcl2 (padj 0.45), Bcl2l1 (0.79), Mcl1 (0.44), Bcl2l2, Xiap, Birc2/3/5
all n.s. on `timepoint_neg`. De-prioritisation is the tissue's route; buffering is the tumour's.

**MODULE VERDICT AND ONE RE-SOURCING (2026-07-27).** Q1's *trigger-supply* route above must cite the
author's **PGC1a westerns (PUMA + BIM)**, NOT the library's `NRF1 -> Bbc3` note — that row is base
rate and fails its own gate. The axis stands; the citation was wrong. And there is **no co-regulated
OXPHOS/apoptosis module** to name in the title: the PGC1a axis is *depleted* of pro-apoptotic genes
and owns the MOMP **payload**, MYC's regulon owns the **effectors**, and the WT gland moves no
apoptotic transcript at all. The reframe that replaces it is stronger — *the gland withdraws the
state the death machinery depends on, so what is shared is not transcription but DEPENDENCE*. Full
verdict, the LE-carrier hypothesis and the mechanism menu: [[oxphos-apoptosis-module-verdict]],
[[le-is-low-expressing]], sections 0.7a/0.8/0.9 of the doc.

**How to apply:** the doc is `docs/2026-07-26_introduction_alignment_and_the_question.md` (question,
sections 0.5/0.5a/0.6/0.7 answering the three questions, clause verdicts, a full drafted
Introduction at axis level -- no PGC1a/PUMA named). All `[43]` numbers are RECONCILED against
`results/substrate_specificity_tradeoff.rds` (author-run 2026-07-26; built-in positive control
returned -0.2548228 against Issue #4's -0.2548). `batch = timepoint` bounds every WT temporal
statement. Related: [[dose-vs-competence]], [[intro-oxphos-restraint-tension]],
[[death-priming-change-untested]], [[correlation-ceiling]].
