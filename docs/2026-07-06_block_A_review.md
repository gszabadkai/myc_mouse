---
title: Block A review — myc_mouse finalisation
date: 2026-07-06
status: draft for author review (narrative decisions pending)
relates_to:
  - "docs/2026-07-04_block_A_build_spec.md (the build sequence)"
  - "docs/2026-07-05_gates_note.md (Day-1 gates)"
  - "docs/myc_mouse_finalisation_plan.md (authoritative plan)"
  - "memory: block-a-day1-gates, myc-death-timing-question"
purpose: >
  Synthesises the whole Block A exploratory sweep (scripts 13-21) against the
  original hypotheses, tiers the claims by evidential strength, lays out the
  honest statistical boundaries, and proposes narrative alternatives + a figure
  split for the author to decide. Newer supersedes older.
---

# Block A review — myc_mouse

## Where we are

The Block A analysis spine is complete: gates + the full GSVA/fork/fGSEA/
permutation-null menu, all run in Positron (Option A) and committed. This note
is the review checkpoint: what survived, what the paper can claim, and the open
narrative/figure decisions. The figure scripts (Block B) and script 22 are not
yet built. Cell-death main-vs-supp placement and the deferred death-timing
analysis are OPEN (see below).

## Analysis inventory (scripts 13-21, with the load-bearing result)

| # | Script | Commit | Key result |
|---|--------|--------|-----------|
| 13 | gate1_divergence_timing | c6c7519 | Genotype effect SHRINKS 6W->12W; ~2-fold attenuation, biological (effect-size, LFC-independent) not power artifact |
| 14 | gate2_apoptosis_readout | 9939a37 | Apoptosis-PRO module coordinately shifts (t p=0.023); cell death stays main-figure candidate (DECISION_POINT) |
| 15 | gsva_scoring | 9939a37 | GSVA engine: 884 library sets x 24 samples (VST, cohort-relative). Also exports expr_mat + pathways |
| 16 | cell_death_binomial_raw | 9939a37 | Branch-1 raw port (directional binomial on raw dLFC) |
| 17 | gsva_overview | 8826696 / 3c3058c / e5c8201 | Trajectory viz + 4-method validation. Cross-check GSVA-vs-count rho=0.66. Between-category permutation + per-sample mechanism (both informative negatives) |
| 18 | ap7_mb_fork_projection | a37ae0d | CENTREPIECE. Myc drives the METABRIC biogenesis fork; genotype main effect p=0.005; MB2_UF-specific Wilcoxon p=0.008 |
| 19 | dev_composition | 6f9c128 | Myc=biogenesis (Felsher~bio 0.90, geno p=3.7e-6); suppresses basal (p=0.016); anti-parallel to lineage development; WT mito reallocation |
| 20 | fgsea_percategory | 2f278c9 | MitoCarta NES 2.18 (padj 1e-55); attenuation-via-convergence (timepoint_pos mito down, timepoint_neg mito up) |
| 21 | ap6_permutation_null | e370d1e | HARD "why mito": mito/OXPHOS most preferentially altered vs expression-matched null; proliferation does NOT survive |
| 22 | reframe_mtdna_abund_retention | afcddfe | mtDNA-up/nuclear-down WT reallocation; abund lens-split time-dependent (6W selective SD 0.14 -> 12W uniform 0.10, fGSEA NES flat 2.2); retention: mito/biogenesis convergence per-category |
| 23 | death_timing_substrate | 2140b37 | Myc death-coupling substrate-gated: priming p=0.038, ~3x at 6W, Gate2 p=0.023; mitonuclear imbalance peaks 6W_pos, death coupling 0.75->0.07; prolif coupling invariant; p53 null |

## Hypotheses -> verdict

- **Development "licenses/permits" Myc (permission verb):** KILLED (Gate 1 -
  attenuation, not permission).
- **Progressive divergence trajectory:** KILLED (convergence; ~2-fold attenuation
  to non-zero, direction preserved).
- **H2 attenuation = power/selection artifact:** REJECTED (effect-size slope 0.46,
  97% attenuate; LFC-independent set same slope).
- **Attenuation mechanism (H1 front-loaded / H3 WT catch-up / H4 selection):**
  RESOLVED as **convergence** = WT developmental catch-up on a shared biogenesis
  program + modest Myc fade. NOT development antagonising Myc (script 19 mechanism
  test refuted it: bio~dev POSITIVELY correlated; adjusting for dev does not remove
  the Myc+ biogenesis decline). H4 selection remains UNTESTED (needs CV / AP8).
- **"Why mitochondria?":** ANSWERED (AP6.2 - preferential mito targeting, expression
  confound defeated).
- **MB-fork cross-species validation:** CONFIRMED (AP7 - Myc -> MB2_UF/LP).

## The claim spine (tiered by evidence)

### Tier 1 - powered, significant, robust (the backbone)
1. **Myc preferentially alters the mitochondrial/OXPHOS compartment** - not generic
   proliferation. AP6.2 effect size (obs-null, myc_6W): OXPHOS 0.201, MYC_TARGETS
   0.190, MitoCarta 0.164 >> generic Metabolism 0.095 >> Proliferation 0.032,
   Mammary 0.023. Proliferation's fGSEA enrichment does NOT survive the magnitude-
   matched null -> the focus is specifically metabolic/mitochondrial. In-vivo "why
   mitochondria."
2. **The Myc footprint IS a mito-biogenesis program.** Felsher(Myc activity) ~
   biogenesis composite rho=0.90; MitoCarta fGSEA padj 1e-55; genotype p=3.7e-6.
3. **Myc drives the METABRIC biogenesis fork toward MB2_UF/LP.** AP7 cross-species:
   genotype main effect p=0.005; MB2_UF-vs-MB1_UF Wilcoxon p=0.008. Attenuates
   6W->12W (0.62->0.34) like Gate 1.
4. **Constant driver, ~2-fold attenuation via convergence.** Myc dose constant (MYC
   signatures stable; mRNA/blot/literature). Myc+ biogenesis fades while WT rises
   (fGSEA: timepoint_pos mito -1.5 vs timepoint_neg +0.8; AP1 mito reallocation).
5. **Myc holds lineage differentiation against development.** Suppresses basal
   (p=0.016); its lineage effect is anti-parallel to the WT developmental direction
   (AP3 rho=-0.28) - decoupling the metabolic arm (rides it) from the lineage arm
   (opposes it).
6. **Myc's death-coupling is substrate-gated through mitochondrial state (script
   23).** Myc raises pro-apoptotic priming (genotype p=0.038; myc_6W PRO module
   p=0.041) ~3x more at 6W than 12W, and it attenuates (Gate 2 p=0.023). The
   mechanism: a MITONUCLEAR IMBALANCE (nuclear - mtDNA OXPHOS mitoPPS) peaks at
   6W_pos (+0.50) and its coupling to death priming is strong at 6W (r=0.75) and
   GONE at 12W (r=0.07); biogenesis-death coupling decays 0.77->0.56; proliferation-
   death coupling is INVARIANT (0.68->0.65 - so the timing is mito, not generic
   proliferation). The 12W decoupling IS the loss of killing. Powered anchors =
   the priming genotype effect + Gate 2; the couplings are striking within-timepoint
   directional (n=12, survivor-biased, uncorrected). Bridges the external IHC 6W-
   death phenotype to the mito thesis. p53/ARF null; the rheostat is BCL2-family.

### Tier 2 - directional / method-independent, NOT per-set significant
- Program-trajectory texture, cross-validated GSVA-vs-count rho=0.66.
- Apoptosis pro-death front-loading at 6W (Gate 2 + CDC_PRODEATH_CICD, four-method
  convergent: int_p 0.036 / camera 0.037 / roast 0.020 / mean_int_stat -1.06).
- Biogenesis-death decoupling over time (intersections category; now mechanised by
  script 23 - the mitonuclear-imbalance coupling, see Tier-1 item 6).

### Tier 3 - NOT claimed / refuted
- Per-set and interaction significance (n=6/group floor, confirmed across lm /
  CAMERA / ROAST / permutation). fGSEA-on-interaction reports 461 "sig" but is
  anti-conservative under inter-gene correlation - the correlation-aware verdict
  (weak) stands.
- The H1-vs-H3 per-set mechanism split as individual claims.
- "Development inhibits Myc" (refuted, script 19).
- Permission / divergence framing (killed, Gate 1).

## Statistical boundaries (state these in the paper)

- **The interaction is underpowered at n=6/group** - four methods agree no per-set
  interaction survives. Powered claims are LEVELS / main effects (genotype, timepoint
  within genotype), which is where Tier 1 lives.
- **Library sets are strongly co-regulated** (inter-gene correlation 0.25-0.36,
  measured) - so rank-enrichment (fGSEA) and independence-assuming nulls are
  anti-conservative. Lean on effect sizes and correlation-aware tests.
- **GSVA adds no power** - it reframes onto programs; the trajectory decomposition is
  directional, corroborated by the count-level cross-check, not independently
  significant.
- **AP6.2 matched null defeats the expression confound but not correlation** - use
  the effect-size RANKING (mito >> comparators), not the absolute z.

## Narrative alternatives (to decide)

Not mutually exclusive; the question is the ORGANISING lead. All are supported by
Tier 1.

- **A. Preferential mitochondrial amplification (title-aligned "integrate").** Myc's
  oncogenic signal is read out PREFERENTIALLY through the mitochondrial compartment;
  the myc_mouse arm supplies the in-vivo preferential-mito (AP6.2) + cross-species
  fork (AP7) evidence. Lead: AP6.2 + AP7. Strongest, most defensible, matches the
  manuscript title verb. RECOMMENDED lead.
- **B. Constant driver, converging substrate (trajectory frame).** Organise around
  the attenuation-via-convergence: constant Myc, WT developmentally catches up, the
  footprint attenuates ~2-fold. Mito is the substrate. Lead: Gate 1 + convergence +
  AP1. Elegant but the "attenuation" risks reading as "Myc weakens" unless carefully
  bounded (it is not - MYC signatures stable).
- **C. Oncogene-development decoupling (lineage frame).** Myc accelerates the
  metabolic/biogenesis arm of development while holding lineage differentiation
  against it; MB2_UF/LP as the dedifferentiation endpoint. Lead: 19 + AP7. Novel and
  mechanistic, but the lineage results are partly descriptive (basal p=0.016 is the
  one powered anchor).
- **D. Cell-death / mito-decision-point frame.** The mito state gates apoptotic
  sensitivity (6W-permissive). UPGRADED by script 23 (2026-07-07): no longer
  contingent - it now has a powered anchor (Myc pro-death priming p=0.038; Gate 2
  p=0.023) AND a mechanism (mitonuclear imbalance peaks at 6W_pos, its death
  coupling r=0.75->0.07). D connects to A: the SAME preferential-mito amplification
  (AP6.2) that defines the Myc footprint produces, at 6W, the imbalanced state that
  gates death. "Mitochondria integrate oncogenic and metabolic programs to shape
  progression" - this IS the integration, mechanised. Still bounded (causal death
  phenotype is external/IHC; this is substrate characterisation, survivor-biased).
  Now a strong secondary thread / main-figure hook, not just discussion.

Leading recommendation: **A as the lead**, B and C woven in as the mechanism, and
**D promoted from held-for-discussion to a secondary thread with a main-figure
hook** (the mitonuclear-imbalance-at-6W_pos panel bridges Fig 1's mito thesis to
the death phenotype). Title verb "integrate" holds; drop "licenses/permits."

## Proposed figure split (draft)

- **Fig 1 - why mito + cross-species:** AP6.2 preferential alteration (21) +
  MitoCarta fGSEA (20) + AP7 MB-fork projection (18).
- **Fig 2 - the biogenesis program and its attenuation:** Felsher~bio + Category-7
  discrimination (07) + attenuation-via-convergence (20 timepoint contrasts, AP1
  mito reallocation 19).
- **Supplementary:** cell death (Gate 2 + branches 12/16 + apoptosis front-loading)
  + the death-timing substrate model (BH3:BCL2/p53 rheostat, decoupling, H1-H4,
  script 23); lineage composition (basal suppression, luminal crossover, 19); GSVA
  trajectory validation + the honest boundary (17); mtDNA / abund / retention
  (script 22).

## Open decisions for the author

1. **Framing** - confirm A as lead ("integrate / preferentially amplify"), drop
   permission?
2. **Cell death - main figure or supplementary?** (REFRAMED by script 23,
   2026-07-07.) The death story now has a mechanistic spine (mitonuclear imbalance
   peaks 6W_pos; death coupling 0.75->0.07) + powered anchors (Myc priming p=0.038;
   Gate 2 p=0.023), so it is no longer merely Tier-2. RECOMMENDED = hybrid (c+):
   the mitonuclear-imbalance-at-6W_pos panel goes MAIN (Fig 1 or 2, as the mito->
   death bridge for the "integrate" title), the full death analysis (BH3:BCL2
   rheostat, branches 12/16, Gate 2, H1-H4) goes to a strong Supp. The CAUSAL death
   phenotype stays external (IHC); we frame the transcriptomics as substrate. Author
   to confirm which main figure hosts the imbalance panel.
3. **Script 22** (mtDNA / abund / retention -> Supp 3) - RESOLVED 2026-07-06b: run
   now to finish Block A (being built).
4. **Deferred death-timing analysis** - RESOLVED 2026-07-06b: run BEFORE the
   figures (its result shapes the final interpretation). Built as the new
   **script 23** (`23_death_timing_substrate.R`), covering all four hypotheses
   (H1 BH3:BCL2/p53 rheostat, H2 biogenesis-death decoupling, H3
   proliferation-apoptosis coupling, H4 selection/CV), anchored on the WT
   `timepoint_neg` substrate. See memory `myc-death-timing-question` for the full
   design + the honest bulk-RNA ceiling.
5. **Block B** - move to the figure scripts (renumbered **24/25/26** after the
   death-timing insert) once 22 + 23 land and the narrative is locked.

## The deferred death-timing thread (do not lose)

The core biological question - why Myc kills at 6W not 12W (IHC, Myc-ER inducible
model, INDEPENDENT of this RNA-seq) - is a SEPARATE, under-served thread. The bridge:
the inducible model isolates the substrate; our chronic RNA-seq characterises the
6W-vs-12W substrate (anchor on the WT/Myc- `timepoint_neg` contrast - the tissue the
acute pulse hits). Design + hypotheses (BH3/BCL2 balance first) + the survivor-bias
ceiling are in memory `myc-death-timing-question`. Not started; slotted after the
current flow by author decision.

## Remaining build items

- **Script 22** - reframe mtDNA / abund / retention (Supp 3 inputs). Off
  mitopps_scores + fgsea_percategory. IN BUILD (2026-07-06b).
- **Script 23** - death-timing substrate model (H1-H4, WT-substrate-anchored;
  integrates cell-death branches + Gate 2 + apoptosis/intersection trajectories).
  IN BUILD (2026-07-06b), before the figures.
- **Scripts 24/25/26** - Figure 1, Figure 2, Supplementaries (Block B; after the
  narrative is locked here). Renumbered from 23/24/25 by the death-timing insert.
- **AP8 / CV** - the selection (H4) lens now lives inside script 23; the optional
  mito-fork CV tail stays a script-26 option only if the selection arm is pursued.
