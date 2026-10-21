---
name: death-priming-change-untested
description: "The death arm has NO independent transcriptomic support: pro_comp IS a MitoCarta set (circular), 'Myc raises priming' is the mito-content effect, and the priming interaction was fitted+saved+never read"
metadata:
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-20T09:38:55.951Z
---

**Found 2026-07-17 (`scripts/34_death_priming_reassessment.R`, `paper-figures`).** The author
challenged script 33's null: it is computed **entirely within 6W** (`33:278-294`, `idx = i6`), so
it tests a CROSS-SECTIONAL claim and **cannot see** the manuscript's LONGITUDINAL claim (priming
collapses by 12W). Correct challenge. Script 34 built the nulls that can. Outcome: **the death
arm has no independent transcriptomic support.**

**1. THE CIRCULARITY (the deepest problem, and it was flagged in 2026-07-08 and never applied).**
`pro_comp`/`priming` are built from `MITOCARTA_APOPTOSIS_PRO`/`_ANTI` (25/9 genes,
`08_mitoPPS_analysis.R:285-290`) — **mitochondrial** gene sets. The imbalance is a mitoPPS ratio.
"The mito state couples to death priming" has been **correlating mito genes with mito genes**.

**2. "Myc raises apoptotic priming" IS THE MITO-CONTENT EFFECT.** Script 32: Myc raises nuclear
MitoCarta ~19% **wholesale**, so any MitoCarta subset rises whatever its label. Conditional test
(vs expression-matched OTHER MitoCarta genes): APOPTOSIS_PRO **+0.196 vs background +0.364,
z=-2.42, p=0.012** — it rises **BELOW average for a mito gene**. On the **404 non-mito** pro-death
genes the effect **does not reproduce** (-0.106, p=0.078, wrong sign). Independently reproduces
script 23's mitoPPS de-prioritisation (-0.21, padj 0.022) — two methods agree. The composite rises
because its genes are MITOCHONDRIAL, not because they are PRO-APOPTOTIC. Same trap as
[[mitocarta-sets-are-membership-loose]].

**3. THE INTERACTION WAS FITTED, SAVED, AND NEVER READ.** `death_timing_substrate.rds
$h1$state_tests` has carried **PRO int_p=0.177 / priming int_p=0.263** (ANTI 0.953, p53 0.799)
since script 23. The narrative quoted **group means** instead ("~3x more at 6W" = +0.48 vs +0.15).
Worse than the imbalance case, where the contrast was never fitted at all. See
[[mitonuclear-imbalance-untested]].

**4. Gate 2's p=0.023 (the sole positive) does not survive.** `14:91` = one-sample t vs mu=0 on 23
genes. `mu=0` is the WRONG null under a global attenuation (every Myc-induced module has negative
mean interaction LFC) -> matched null **z=-1.27, p=0.21**. Measured inter-gene cor **0.08** ->
effective n **8.4, not 23**. Its own question asked correlation-aware = CAMERA: **p=0.298,
FDR=0.926**. Script 17 already ran it.

**5. THE ATTENUATION IS GLOBAL — the mechanism answer.** vs expression x dispersion-matched nulls,
**nothing** is in excess: death z=0.97, non-mito death z=0.71, **OXPHOS z=0.90, MYC-core z=0.29**
(all BH p>0.14). The comparators failing too is the **control that validates the test**, not the
test failing. => **priming fades because the MYC PROGRAMME fades** — canonical Evan/Lowe
BH3-only biology. No death-specific gate, no mitochondrial gate.

**6. The Tang NES -1.38 is a SHARED developmental decline.** `temporal_neg` = **-1.31** (95% as
much). By Issue #6's identity a shared decline **cancels in the genotype gap**. Never cite
`temporal_pos` alone. (Also: interaction fGSEA CD_Apoptosis padj=0.914, all 15 modalities >0.12.)

**7. THE ONE NUANCE THAT IS NOT FLAT — do not over-correct.** The *collapse* (r 0.75->0.07) sits
at the **96th percentile** of Δrho vs all 884 sets, **perm p=0.066**; Fisher z p=0.055. **Marginal,
not null** — better than the 6W coupling (p=0.49) ever was. BUT: **Pearson-only** (Spearman
p=0.31), fails multiplicity, and **evaporates on every non-mito death axis** (p=0.19-0.92). The
only axis that gets this far is the **circular** one. Ambient context: the median set's coupling
already falls **0.71 -> 0.35** (74% -> 11% of sets above |r|=0.5) while **PC1 is UNCHANGED (43% ->
46%)** — so "12W has less structure" is NOT the explanation, but everything decouples anyway.

**ATTRIBUTION:** r=0.75/p=0.005 is `23_death_timing_substrate.R:326-327` (**Pearson**), NOT script
25. Script 25 Part B (`25:140-148`) uses `stats::cor` = Spearman, **no p-value**, and gets 0.671.
Fixed in 33's header.

**DESIGN LIMIT (state it, don't test it):** the design is **cross-sectional** — different animals
at 6W and 12W. **No within-animal delta exists**, so "does the collapse correlate with mito
reprioritisation" is **not estimable per-sample**; only two 4-point group trajectories, no
inference. Both interactions are null anyway (priming 0.263, imbalance 0.916).

**8. THE DIRECTION IS OPPOSITE TO THE NARRATIVE.** All **4/4 non-mito/general death axes** have a
**NEGATIVE** Myc genotype effect; **2 significant after BH** — `APOP_BH3_REACTOME` **-0.367
(p=0.031)**, `CDC_PRODEATH_APOPTOSIS` **-0.149 (p=0.015)**. On every death axis that is not
mito-defined, **Myc LOWERS the death programme at 6W**. Read WITH the survivor bias, not around
it: we sequence the cells that did NOT die, so a suppressed survivor pool is what OIA would
leave. The walkthrough calls survivor bias "a conservative floor"; this suggests it can **INVERT
the sign**, not just shrink it.

**9. The binomial's null was never 50%, and fixing it flips branch 1.** Genome-wide **59.6%** of
expressed genes have dLFC>0 — **93.9% among Myc-INDUCED genes**. So pro-death at 58.6% looks
significant vs 50% (binom p=0.010) but sits **1.1 points BELOW** the genome-wide baseline. Script
16's 52.7%/p=0.27 got the right answer through the wrong null.

**10. MY OWN MATCHED NULL IS ANTI-CONSERVATIVE — do not quote its p.** `matched_null()` draws
genes **independently** from bins, so the null of the mean is too narrow for co-regulated sets.
**Script 21 documented exactly this for AP6.2** ("raw z inflated by inter-gene correlation → use
EFFECT-SIZE ranking, not absolute z"). I criticised Gate 2 for assuming gene independence while
quoting p=0.012 from a null that assumes gene independence — the same error one level up.
**Non-significant matched-null results are SAFE** (a correction only widens the null), so #5's
"nothing is in excess" and #4's Gate 2 verdict are firm. **Significant ones are not safe on the p
alone** — #2 survives on **effect size** (half the background rise) + the independent mitoPPS
convergence, not its p.

**How to apply:** the RNA-seq claims **no death result**. The death PHENOTYPE is external
(IHC/inducible) and untouched. The 6W-vs-12W timing rests entirely on it. **BH3 profiling** is the
experiment — it measures the apoptotic gateway directly and is independent of the content claim.
Distinguish the paper's **two mitochondria**: (i) biosynthetic organelle (content +21-27%, SOLID —
[[mito-content-bounded]]) vs (ii) apoptotic gateway (BH3/MOMP/PUMA). The bridge assumed (i)->(ii);
(ii) never needed (i). `Bbc3`/PUMA is the one gene moving right on the interaction (z=-2.6, padj
0.84). CEILING: n=6/group — **unsupported is not refuted**. See
[[sample-level-couplings-need-a-null]], [[dont-assert-a-mechanism-without-checking-it-can-work]].
Author had NOT yet sourced 34 when this was written.

**UPDATE 2026-07-19 (`scripts/38_mitopps_priming_and_pgc1a_biogenesis.R`, author-run).** The
circularity (#1) was GSVA/absolute-space-specific. The author's point: mito sub-pathways are
SEPARATELY regulated, and a mitoPPS RATIO cancels both the global shift AND the shared mito-content
mode -- so OXPHOS<->priming in mitoPPS space is NEITHER circular NOR ceiling-bound. Tested there
(Apoptosis-PRO/-ANTI + OXPHOS are nuclear-encoded mitoPPS pathways, script 08 PART 2c, no mt- genes
so no mtDNA-arm border): **mitoPPS space IS readable** (ambient 0.42 vs GSVA 0.80; OXPHOS<->global
0.46 vs 0.97). OXPHOS<->priming(PRO-ANTI) = raw Spearman **+0.45** (null 0.42) / design-adj **+0.70**
(null 0.48) -- but only the **52nd percentile** of OXPHOS's couplings (MID-PACK), genotype
delta-interaction marginal (p=0.08). A **WEAK but CLEAN lead for BH3 profiling**, not a result --
priming still != death. AND: **MYC- and PGC1a/ESRRA-biogenesis are NOT separable** (r=0.93, both
Myc-induced beta +0.57/+0.75; the global-adjusted -0.03 is the [[correlation-ceiling]] OVER-REMOVAL
trap, resid_frac 0.01/0.06), so the hypothesised MYC-vs-PGC1a death contrast is UNTESTABLE in bulk
-> needs an ESRRA/PGC1a PERTURBATION. See [[batch-equals-timepoint]].

**UPDATE 2026-07-20 (the perturbation was DONE externally).** Lesner/Simon (bioRxiv
2026.07.13.738248, `docs/MYC_mtDNA_paper.pdf`) ran exactly the perturbation script 38 said was needed:
**sgNRF1/sgDRP1/sgFUNDC1** on the MYC mito axis. Their reading: MYC->NRF1 drives mtDNA depletion as
protective MITOPHAGY; apoptosis is engaged by forced FUSION + BH3 mimetics -- consistent with our
finding that the death arm has no transcriptomic support and that priming != death. This ENDORSES
**BH3 profiling** as our experiment and reframes the death mechanism away from priming. See
[[lesner-mtdna-external-validation]], [[mtpct-is-three-way-confounded]].
