---
name: myc-death-timing-question
description: "The 6W-vs-12W Myc cell-death phenotype (IHC, inducible model) and the deferred integrated transcriptomic analysis to address it"
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
---

**Phenotype (solid, INDEPENDENT of this RNA-seq).** In a parallel Myc-ER
tamoxifen-inducible model (NO RNA-seq there), acute Myc induction at 6W causes
MUCH more cell death than induction at 12W, measured by IHC. Constant acute Myc
pulse, different tissue age => the death is SUBSTRATE-GATED by developmental
stage, not by Myc dose. This is the core biological question the myc_mouse arm
should ultimately explain.

**Why our (chronic MMTV-Myc) RNA-seq can address it.** The inducible model
isolates the substrate effect; our data characterises the 6W vs 12W substrate.
Bridge: the tissue the acute pulse hits is ~normal (Myc off pre-tamoxifen), so
the **WT / Myc-negative 6W->12W trajectory (the `timepoint_neg` contrast)** is the
cleanest readout of the death-permissive state. Anchor the future analysis there,
then ask whether chronic Myc+ modulates that priming.

**IN BUILD as script 23 `23_death_timing_substrate.R` (2026-07-06b).** Pulled
forward BEFORE the figures (result shapes the final interpretation); author chose
the BROAD scope = all four hypotheses below. Figure scripts renumbered 23/24/25 ->
24/25/26. Build one 6W-death SUBSTRATE model that integrates cell-death branch 1
(script 16, raw binomial) + branch 2 (script 12, Tang RCD fGSEA) + Gate 2
(script 14) + script-17 apoptosis/intersection trajectories, anchored on the WT
`timepoint_neg` contrast (chronic-Myc+ layer secondary/survivor-biased). Testable
hypotheses in the EXISTING data, highest value first:
1. **BH3-only : anti-apoptotic BCL2-family balance + p53/ARF readiness** in the WT
   6W->12W substrate (the canonical Myc-apoptosis axis; Bbc3/Puma already a Gate-2
   outlier). Engaged at 6W, rebalanced anti-apoptotic by 12W?
2. **Biogenesis-death DECOUPLING over time** - links to the biogenesis thread
   (MYC_SPECIFIC_MITO_AND_CDC_PRODEATH front-loads then decouples in script 17) and
   the ER/PGC1a-vs-Myc coupling thread in [[block-a-day1-gates]]. The uncoupling
   may BE the loss of killing.
3. **Proliferation-apoptosis coupling** - 6W pubertal/proliferative state as
   permissive context (oncogene-induced apoptosis).
4. **Selection (H4)** - sensitive cells culled by 6W, resistant survivors by 12W.

**Honest ceiling (state as a manuscript limitation).** Bulk RNA, one timepoint per
age, n=6, and SURVIVOR BIAS (you sequence the cells that did NOT die) mean we can
characterise the death-permissive STATE (association) but cannot prove causation
or fully resolve mechanism; single-cell or functional assays (e.g. BH3 profiling,
caspase-by-state) would be needed for that. The best current transcriptional lead
is pro-death CICD/CDC priming front-loaded at 6W (script 17, CDC_PRODEATH_CICD,
four-method convergent but not FDR-significant at n=6/group).

**RESULT (script 23, commit 2140b37, 2026-07-07).** The model resolved better than
the framing predicted - and CORRECTED one hypothesis:
- **"WT substrate de-primes with age" is REFUTED.** WT baseline BH3:BCL2 priming is
  flat/slightly-rising (6W_neg -0.24 -> 12W_neg -0.08). The substrate-gating does
  NOT act on the WT baseline. Do not anchor the story on WT de-priming.
- **The signal is Myc's COUPLING, and it is powered.** Myc raises pro-apoptotic
  priming (PRO-ANTI genotype effect +0.32, p=0.038; myc_6W PRO module +0.20,
  p=0.041) ~3x more at 6W than 12W (group-mean gap +0.48 vs +0.15); attenuates
  (Gate 2 PRO interaction -0.13, t p=0.023; branch-2 apoptosis temporal_pos NES
  -1.38). Survivor bias is CONSERVATIVE here (dead cells unsequenced -> true
  engagement higher). p53/ARF NULL (p53 composite flat) - rheostat is BCL2-family.
- **Mechanism = MITONUCLEAR IMBALANCE (H2, the standout).** nuclear-minus-mtDNA
  OXPHOS mitoPPS peaks at 6W_pos (+0.50, the death-permissive group); its coupling
  to PRO priming r=0.75 (p=0.005) at 6W -> 0.07 (p=0.82) at 12W. biogenesis~death
  0.77->0.56. Two-lens nuance: Myc raises PRO expression absolutely (+0.20) but
  RELATIVELY de-prioritises Apoptosis-PRO within mito (mitoPPS diff -0.21, padj
  0.022) - biogenesis crowds it out = the imbalance. Ties to script 22 (6W selective/
  imbalanced -> 12W uniform mito) and [[block-a-day1-gates]] reprioritisation.
- **H3 proliferation~death coupling INVARIANT (0.68->0.65)** -> the timing is a
  mito/biogenesis phenomenon, not generic proliferation. H4 CV: Myc+ narrows
  0.017->0.007 (weak, on-narrative). Branch 1 raw binomial null (52.7%, p=0.27).
- **Reframed model:** Myc death-coupling is substrate-gated via mitochondrial state.
  At 6W Myc forces a mitonuclear-imbalanced biogenesis state tightly coupled to
  death priming (oncogene-induced apoptosis); by 12W biogenesis is uniform/balanced
  and DECOUPLED from death, Myc's death drive faded. The 12W DECOUPLING IS the loss
  of killing (this settles the [[block-a-day1-gates]] ER/PGC1a-vs-Myc thread: Myc
  biogenesis STARTS coupled and BECOMES uncoupled - not constitutively uncoupled).
- **Review impact:** narrative D promoted from contingent to a secondary thread with
  a MAIN-figure hook (the imbalance-at-6W_pos panel bridges Fig 1 mito thesis ->
  death); open decision #2 recommended hybrid (imbalance panel main, full death
  analysis Supp). See docs/2026-07-06_block_A_review.md Tier-1 item 6.

**RESULT 2 (scripts 24+25, commits ef6b04a + b50e00a, 2026-07-07).** The review
discussion (docs/BlockA_review_discussion.md) added the DEVELOPMENTAL "why". Key
corrections/additions to the model above:
- **The developmental SUBSTRATE (not chronic-Myc adaptation) sets the timing** -
  matures on TWO coupled axes. (1) MITO axis (robust): the mtDNA-encoded vs
  nuclear-OXPHOS mitoPPS imbalance peaks at 6W_pos (+0.50/+0.56 by defn), inverts
  by 12W; its resolution is ~88% mtDNA-RISE-driven in BOTH genotypes (WT 89%, Myc
  87%) = genotype-shared developmental mito maturation. (2) CELL-STATE axis: WT
  lineage shift (luminal/progenitor down, basal up) is DIRECTIONAL but UNDERPOWERED
  at n=6 (all p>0.18) - do NOT claim a significant composition shift.
- **Imbalance is a WITHIN-mito quantity = mitoPPS ALONE.** fGSEA is a SEPARATE
  importance-vs-transcriptome axis, never merged into the imbalance (author
  correction). mtDNA/nuclear OXPHOS ARE separated in the mitoPPS set (13 mt-* in a
  synthetic pathway). The WT substrate matures from nuclear-led ASSEMBLY (6W:
  complexes/import/chaperones prioritised, mtDNA lagging) to mtDNA-completed
  RUNNING metabolism (12W). Chaperones FALL with age -> the UPR^mt/ISR frame is
  REFUTED; it is assembly->running maturation, not a proteostasis-stress response.
- **CORRECTION to which axis tracks death (script 25 Part B, Spearman vs pro_comp,
  within tp):** at 6W bio_comp 0.83 > MASC/stem 0.78 > mitonuclear_imbalance 0.67;
  ALL decouple by 12W. The mito imbalance is NOT the sole death-coupled axis -
  STEM/progenitor (MASC) content is comparably coupled and decouples with age
  (vindicates the "composition matters" intuition). Coupled death-permissive
  STATE, not separable causes (n=12/tp, survivor-biased).
- **Biogenesis lanes (script 24): Myc drives MYC and ER/PGC1a biogenesis
  near-identically** (Cat-7 not separable by trajectory) - Myc CO-OPTS the
  ER/PGC1a/PGC1a-axis-TF machinery (ESRRA/NRF1/GABPA co-move with MYC; ESR1/ER
  antagonised, re-emerges at 12W). Balanced CORE/MITO_NU biogenesis-death
  intersections are MOST death-coupled at 6W (r 0.77/0.83) - where "ER/PGC1a kills"
  lands.
- **Model closeout:** young 6W substrate (stem-rich, mitonuclear-imbalanced,
  biogenesis-primed) is death-permissive; matured 12W substrate (mtDNA caught up,
  differentiated, converged) is decoupled. Myc-ER anchor: acute pulse hits the
  young imbalanced substrate -> mismatch -> death; by 12W matured -> less death.
  Ceiling unchanged: n=6, survivor bias, correlated axes -> association not
  causation.

**p19ARF/p53 axis re-tested at gene level (2026-07-13).** The canonical Myc escape (ARF loss
/ p53 suppression -> survival) does NOT explain the 12W death drop: Cdkn2a is MORE Myc-induced
at 12W (not lost), Mdm2 flat, p53 target panel not coordinately reduced (p21 less suppressed;
Zmat3/Eda2r/Phlda3/Trp53inp1/Sesn2 up). Only Bbc3/Puma (int z -2.6, ns after FDR) + Bax decline
= BCL2-family EFFECTORS (the BH3 rheostat downstream of p53), not a p53-activity failure.
Corroborates "p53/ARF null." CANONICALISED as script 23 PART 2b (arf_p53_genes + data-driven
verdict -> death_timing_substrate.rds$h1$arf_p53; figure h1c_arf_p53_axis.pdf); logic validated
read-only (verdict=REJECTED, p53-target mean int +0.192 t p=0.18). Author must SOURCE script 23 to
regenerate the rds/pdf. Lives on paper-figures; not yet ported to analysis-exploratory/block-a-reviewed.
Detail in docs/2026-07-13 walkthrough Sec 1B H1. Do not re-litigate.

**p53-INDEPENDENT PUMA regulators tested (2026-07-13; script 23 PART 5b, per
docs/library_reference/PUMA-and-its-relationships.md).** POSITIVE: a FOXO3/PUMA/HTRA2 module
co-moves -- Foxo3 mirrors Bbc3 (dMyc -0.49, int z -2.7 ~ Bbc3 -2.6), Htra2 most 6W-induced
(+1.49 -> dMyc -0.55): a p53-independent candidate for the PUMA loss. NEGATIVE/bounding: PGC1a
BRAKE rejected (ESRRA/PGC1a co-decline, don't rise); FOXO output PUMA-selective (Bim/Noxa flat);
and the FOXO3 target-ACTIVITY signature does NOT track Bbc3 at 6W (r -0.14) nor decline (int
+0.07) -> the Foxo3-gene co-movement is NOT a confirmed FOXO3-activity cascade -> reads as
DE-AMPLIFICATION of the death-primed compartment (compositional spine), not a TF switch.
Directional, n=6, no FDR. -> death_timing_substrate.rds$puma_regulators. Do not re-litigate.

**BMYO/LASP death-dropout confound BOUNDED (2026-07-15; script 26 PART E ->
dev_program_myc_integration.rds$selection_bound).** Author asked whether the Issue #1 BMYO
REPRESSION (and a negative LASP interaction) could be compositional dropout of dying high-MYC cells
rather than a per-cell change. Never tested head-on before in Block A (survivor-bias/composition
machinery only ever aimed at the proliferative/TEB compartment + death-priming program). Two lenses
on on-disk data: (1) per-sample coupling of each MEC lineage composite to pro-death priming
(pro_comp), (2) marker-overlap enrichment vs the 512 pro-death genes. VERDICT: **BMYO = PER-CELL,
not dropout** -- BMYO is OFF the death axis (rho ~0.15/0.13 at 6W/12W; script-19 BASAL_comp even
lower -0.04) and its program is not pro-death-enriched (OR 0.92, p 0.81); plus MMTV drives MYC in
LUMINAL not basal cells, and Myc CLOSES basal chromatin per-cell (CHUNG ATAC). LASP is heterogeneous:
STATE interaction net POSITIVE, but the 4-set luminal-PROGENITOR arm (LAPRO/LP_OPEN/ALVSEC/ALVPROG)
has an all-negative interaction (Myc induces it +0.31@6W then FADES to +0.24), is the MOST
death-coupled lineage (rho ~0.87@6W) and pro-death-enriched (OR 1.60, p 0.006) -> death-dropout NOT
excludable for that arm = the deconvolution-plan target (basal + luminal-progenitor fractions).
Bounds not settles (n=6, survivor bias, set-overlap != death). Script 26 SOURCED by author
(2026-07-15; numbers reproduce read-only validation exactly, PDFs generated). Lives on paper-figures.
Detail in docs/2026-07-13 walkthrough Issue #1 + Sec 5. Do not re-litigate.

**THREE-ROUTE lineage adjudication (2026-07-15; synthesis of the PART E run + the exact identity
int = mycpos_shift - wt_shift).** Every Myc-driven lineage shift decomposes into SELECTION (dropout of
a death-primed subpop) vs FADING (per-cell program change) vs CONVERGENCE (WT background maturing).
Each MEC lineage lands on a DIFFERENT dominant route: BMYO = FADING (per-cell de-basalisation; WT flat
so no convergence, off death axis so no selection; ATAC closes basal chromatin per-cell); LASP state =
CONVERGENCE (int +0.135 is ~100% WT luminal decline, Myc+ flat); LASP luminal-progenitor arm =
SELECTION candidate (the ONLY place dropout survives); LHS = per-cell GAIN (mycpos +0.19, ATAC opens
ML) + convergence. Net proof: no single mechanism; selection is REAL but LOCAL (progenitor arm only) --
extends Issue #6's convergence+fade decomposition by adding the third axis. Written into docs/2026-07-13
walkthrough Issue #1 (three-route table) AND the intro alignment Sec 3 (claim #4 STRENGTHENED) + Sec 4
(rewrite gained an optional cell-intrinsic-reprogramming clause). All committed+pushed on paper-figures
(through 3e8d0bf). Do not re-litigate.

See [[block-a-day1-gates]] for the gate/validation context.
