---
name: mitonuclear-imbalance-untested
description: "The mitonuclear imbalance has no supported contrast (all p>0.24) and its death coupling fails a null - the mito->death bridge of the manuscript's arc is not established"
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-20T09:38:49.988Z
---

**The manuscript's mito->death bridge is NOT established.** Found 2026-07-16/17 (scripts 32/33,
commits `50ad21c`..`2c5b79e` on `paper-figures`). Two independent failures:

**1. The mitonuclear imbalance has NO supported contrast.** Script 24 reports `mtnuc_index` as
FOUR GROUP MEANS (0.130 / 0.499 / -0.456 / -0.173, *"peaks at 6W_pos, inverts by 12W"*) and the
corpus quoted that downstream as a demonstrated effect. **A contrast was never fitted.** When
fitted (script 33 PART E1): time **p=0.313**, genotype **p=0.523**, interaction **p=0.916**,
Myc-at-6W **p=0.249**. Nothing significant on any term. Adjusting for the IEG axis takes the time
beta -0.586 -> -0.011.

**2. The death couplings FAIL an empirical null** (script 33 PART C). See
[[sample-level-couplings-need-a-null]]. bio_comp r=0.83 vs **null median 0.79**.

=> **"At 6W the mitonuclear-imbalanced substrate is coupled to pro-apoptotic priming, so Myc's
induction lands as killing; by 12W the imbalance resolves, the coupling collapses, killing fades"
is retracted as a mechanism.** With it: *"the OXPHOS reprioritisation and the apoptosis
desensitisation are two faces of one developmental mitochondrial maturation -- the causal spine of
the Introduction's closing sentence."*

**CEILING: unsupported is NOT refuted.** n=6/group with mitoPPS variance. The pattern may be real
and underpowered. But it cannot be reported as a finding, and "peaks at 6W_pos, inverts by 12W"
reads as a demonstrated effect when it is four averages with p=0.31 between them.

**What SURVIVES (do not over-correct):**
- Myc raises **mitochondrial content +21-27%** (script 32, p=0.0006 adjusted) — new and strong.
- Myc raises absolute nuclear OXPHOS (d=1.27 script 29 / d=1.15 script 32, two methods).
- **Issue #4's attenuation is PROTECTED** (script 33 PART E2): 12W is NOT noisier (median SD
  0.2935 vs 0.3515), Myc-contrast SEs ~identical (0.251 vs 0.259), so the DE collapse 1967->135
  tracks |LFC| falling 0.268->0.204 at constant SE = real magnitude loss, not power/batch. The
  IEG axis is genotype-independent so it cancels within each timepoint's contrast.
- Issues #1/#2/#3/#5/#6 (genotype/interaction).
- **The death PHENOTYPE** (Myc kills at 6W not 12W) — external IHC/inducible, never depended on
  this RNA-seq.

**How to apply:** the paper **keeps both halves and loses the bridge**. Myc builds mitochondria
and the program attenuates; the gland becomes death-resistant; whether the mito state is *why* is
an OPEN QUESTION. Do not restore the bridge without new data. **mtDNA qPCR** (copy number) +
**TOMM20/VDAC/CS blot** close four open items between them, including tension A. Written up in
`docs/2026-07-17_evidence_audit_and_narrative.md` (read first) and marked in place throughout the
walkthrough (§0 banner, §1B, §2, §3/4, §5). See [[mito-content-bounded]],
[[dont-assert-a-mechanism-without-checking-it-can-work]]. Do not re-litigate.

**UPDATE 2026-07-19 (tested directly).** A reviewer/author asked whether the mtDNA<->priming
correlation is citable now that mitoPPS gives a clean-null space. TESTED: mtDNA-encoded-OXPHOS <->
priming in mitoPPS space is **-0.71 Spearman (95th pctile, clears the null)** -- looks strong -- BUT
`cor(mtDNA-arm, mt%) = 0.952` (the mtDNA-arm mitoPPS score IS mt%), and adjusting for mt% collapses
it **-0.71 -> -0.20**. The imbalance<->priming is identical (`cor(imbalance, mt%) = -0.94`). So the
mtDNA/imbalance<->priming coupling is **NOT citable -- it is the mt% confound**, exactly the
mtDNA-arm border (mitoPPS cancels a global shift but NOT the mtDNA fraction). CONTRAST: **nuclear**
OXPHOS<->priming is mt%-ROBUST (+0.45 -> +0.49 adjusted; `cor(nuclear OXPHOS, mt%) = -0.71`) -- the
one weak-but-real transcriptomic lead (mid-pack, script 38) -> BH3 profiling. mtDNA questions ->
**mtDNA qPCR**. Full reasoning: `docs/2026-07-19_oxphos_axis_biology_and_mtdna_priming.md`. mt% is
three-way-confounded and death-agnostic -- see [[mtpct-is-three-way-confounded]].

**UPDATE 2026-07-20 (external validation, NOT a rescue).** Lesner/Simon (bioRxiv 2026.07.13.738248,
`docs/MYC_mtDNA_paper.pdf`) show the MYC nuclear-up/mtDNA-down discordance in TCGA (flash-frozen, no
digest) + mouse HCC, with a non-mito MYC control. This **externally VALIDATES our imbalance
PHENOMENON and kills the "our dissociation artifact" worry** -- but it REFRAMES the death arm as
protective MITOPHAGY / turnover (MYC->NRF1->DRP1/FUNDC1), apoptosis only via forced fusion + BH3
mimetics, so it does NOT rescue our mt%<->priming coupling (consistent with our null). Honest nuance:
their discordance is a MYC CONTRAST; our nuclear-up half reproduces as a clean GENOTYPE effect but our
mtDNA-down half is genotype-INDEPENDENT (p~0.49). Validates the QUESTION, pinpoints our underpowered
half. See [[lesner-mtdna-external-validation]].
