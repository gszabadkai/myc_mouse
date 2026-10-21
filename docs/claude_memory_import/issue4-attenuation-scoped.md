---
name: issue4-attenuation-scoped
description: "Block A revision Issues #4/#5/#6 (6W->12W attenuation): DONE (scripts 29/30/31, committed). #4 NES paradox (magnitude not rank), absolute nuclear OXPHOS up. #5 no MYC-independent axis absorbs it (dev/TF co-drive). #6 CITABLE MECHANISM: attenuation = ~66% Myc-fade + ~34% WT-convergence (biosynthetic-arm-specific), OXPHOS/MYC-core fade-only, MYC-core protected; E2F/prolif co-decline hints composition (deconvolution deferred)"
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
---

Block A revision, branch `BlockA-revision-step-by-step`. Issues #1/#2/#3/#4 DONE +
committed. **Issue #4 DONE 2026-07-11** (`scripts/29_attenuation_decomposition.R`, sourced
clean; results/attenuation_decomposition.rds + outputs/attenuation_decomposition/6 pdfs).
Full outcome in `docs/2026-07-08_BlockA_revision_plan.md` Issue #4 section.

**RESULTS (the LFC data vindicated the correction):** (A) OXPHOS FALLS in BOTH genotypes
(OXPHOS_SUBUNITS WT temporal -0.25 / Myc+ -0.39); biosynthetic RISES in WT (amino-acid
+0.18) -- the aggregate artifact confirmed. Magnitude attenuation carried MORE by mito/
metabolic arms (attn ratio ~0.53) than MYC-target core (0.74). NES flat/high at both ages.
(B) Absolute nuclear OXPHOS mRNA: Myc RAISES it (genotype d~1.27 p0.004, robust across VST/
normcnt/logCPM, agreement Spearman >0.997); WT/Myc+ decline with age directional/ns
(n=6/tp). So mitoPPS "OXPHOS down" is NOT an absolute-mRNA drop -- it is reprioritisation.
(C) Mitonuclear shift REAL in absolute mRNA (nuclear peaks 6W_pos, mtDNA rises to 12W --
both lenses agree; confirms script-24). (D) within-WT (n=12) OXPHOS tracks developmental
luminal axis (partial 0.71) + ER/PGC1a TF activity (0.70-0.83), R^2 0.86 -- hypothesis-
generating, Issue #3 non-MYC OXPHOS axis stays OPEN.

Issue #4 = what IS the 6W->12W attenuation (the "NES paradox"). Resolved conceptually:
fGSEA NES is rank-based / magnitude-blind so it is FLAT across age (mito stays top-enriched
at both timepoints); what SHRINKS is MAGNITUDE (DE count 2777->239, |LFC| slope ~0.46,
Issue #3 matched-null excess). Does NOT shrink: (i) Myc transgene level, (ii) program
rank (NES).

**CORRECTION (2026-07-11, author caught it) -- overrides the old '+0.82 WT mito rise'
convergence story in [[block-a-day1-gates]] (line ~107 'timepoint_neg mito UP +0.8').** That
+0.82 was the MEDIAN over 63 heterogeneous MitoCarta fGSEA sets, carried by the BIOSYNTHETIC
arm (amino-acid metab +1.86, lipid +1.51, metabolism(n=414) +1.66). The actual OXPHOS
pathways FALL in WT with age: MITOCARTA_OXPHOS timepoint_neg -1.94, OXPHOS_SUBUNITS -2.47,
Complex I -2.03 -- and fall FURTHER in Myc+ (OXPHOS timepoint_pos -2.61). So there is **NO
WT developmental OXPHOS rise**; OXPHOS declines in BOTH genotypes 6W->12W (Myc+ faster). Any
gap-shrink at OXPHOS is Myc+ falling toward WT, not WT rising. => **whole analysis is now
PATHWAY-RESOLVED (OXPHOS / nucleotide / TCA -- the Issue #3 detached set -- each on its own);
NEVER report an all-MitoCarta aggregate as "mito."** The tidy "Issue #4 closes the Issue #3
loop via a WT OXPHOS rise" is RETRACTED: the non-MYC OXPHOS axis stays GENUINELY OPEN (what
Part D probes, not something Issue #4 pre-answers).

**mtDNA contained (verified read-only):** fGSEA `MITOCARTA_OXPHOS` (n=134) == `MITOCARTA_
OXPHOS_NU`, leading edge holds ZERO mt-* genes -> fGSEA pathway NES already nuclear-clean.
mtDNA bias lives only in the GSVA `MITOCARTA_ALL` composite (Issue #3) and mitoPPS (which
isolates the 13 mt genes into a synthetic pathway -- handled). Part B excludes mt-* from
nuclear OXPHOS and reports mtDNA separately.

**TPM decision:** SKIPPED as statistically REDUNDANT, not for lack of data. Gene length is
a per-gene constant that cancels under z-scoring in a within-gene between-group comparison
(TPM==log-CPM in shape). Lengths ARE obtainable (biomaRt cached; counts are RSEM-style per
`FULL.DAT.csv` gene_id+transcript_id(s)) if a cross-gene abundance question ever arises.
Part B normalization set = VST + DESeq2 normalized-counts + log-CPM.

Script 29 (reframe): Part A decomposition + NES-paradox figure PATHWAY-RESOLVED; **Part B =
ABSOLUTE nuclear-OXPHOS mRNA across 4 groups, normalization-robust** (the key un-done check
-- mitoPPS "OXPHOS down" is only a reprioritisation ratio, never verified in absolute mRNA);
Part C mitoPPS-vs-absolute reconciliation; Part D first-pass within-WT regulator regression
(OXPHOS trajectory = a DECLINE ~ developmental axis + ER/PGC1a TFs ESRRA/NRF1/GABPA + mtDNA;
hypothesis-generating, first pass INSIDE script 29). See [[block-a-day1-gates]] and
[[mec-consensus-nomenclature]].

**Issue #5 DONE 2026-07-11 (script 30 `30_attenuation_moderation.R`, sourced clean, committed).**
Author asked: can MYC-independent factors (developmental luminal + ER/PGC1a TF activity) FULLY
EXPLAIN the OXPHOS attenuation? **ANSWER = NO -- no axis absorbs it.** The attenuation IS the
genotype x time INTERACTION (moderation); the right test = does conditioning on a candidate axis
SHRINK b_int? Reframe on script 29's `$defs` (no re-run); baseline b_int reproduces Issue #4
absolute_stats EXACTLY (oxphos_abs -0.279 p0.607 directional; oxphos_gsva ~0).
- **RESULTS:** dev+tf held -> b_int essentially UNCHANGED (-0.279->-0.322, delta -0.043);
  dev_luminal alone makes it LARGER (-0.721, delta -0.442 = script-17 co-drive signature, removing
  a co-driven program EXPOSES more Myc effect). EVERY delta_b_int bootstrap CI straddles 0 (n=24)
  -> no RELIABLE movement either way. absorption_frac CIs uninterpretably wide (+-10..18) b/c
  baseline is directional -> "fully explains" nowhere distinguishable. **Level-association != 
  moderation:** all axis LRT sig (5e-4..1e-10) = axes DO add OXPHOS-LEVEL variance but do NOT
  absorb the INTERACTION (exactly why Issue #4 Part D level-assoc was the wrong instrument).
  Within-WT (n=12) OXPHOS slope only shrinks ~29%/19% under dev+TF -> PERSISTS.
- **CEILING:** endogenous mediators (Myc drives dev/TF -> biased), n=6/group, collinear, baseline
  int ns -> BOUNDS not settles; identifying test is experimental (inducible Myc / TF perturbation /
  single-cell / ATAC-ChIP, listed in $settle_it). Full outcome in the revision-doc Issue #5 section.
  Artifacts: script 30 + results/attenuation_moderation.rds + outputs/attenuation_moderation/ (4 pdf).
  Metrics primary = delta_b_int + nested LRT (NOT the unstable fraction).

**Issue #6 DONE 2026-07-12 (script 31 `31_attenuation_mechanism.R`, sourced clean, committed).**
Author: compression re-describes, need a CITABLE MECHANISM. Restatement "MYC signal not reduced,
acts on a different background" -> 3 mechanisms (moving background / Myc-fade / breadth
compression). Breadth compression real but MODEST (SD ratio 0.83) + concentrated in strong
responders (DE-count -91% is a tail artifact overstating -20% magnitude) -> not the mechanism.
- **THE IDENTITY (cor 1.0000):** attenuation per gene = interaction = `myc_12W-myc_6W ==
  timepoint_pos-timepoint_neg`. Aligned to d=sign(myc_6W): **WT-convergence** (d*timepoint_neg,
  WT matures toward Myc) + **Myc-fade** (d*timepoint_pos, oncogenic retreat). atten = wt_conv -
  myc_fade.
- **PART A (POWERED, robust across divergent + effect|LFC6|>0.5 universes):** global ~66% fade /
  34% convergence. PROGRAM-SPECIFIC: biosynthetic arm CONVERGES (amino-acid 37% / lipid 41% /
  nucleotide 31% / mammary-luminal 36% -- WT ramps onto Myc's biosynthetic axis); OXPHOS +
  MYC-target core do NOT (conv NEGATIVE -11..-63%, WT diverges, gap closes by fade only);
  MYC-core also LEAST attenuated (0.17-0.19) = PROTECTED; proliferation = pure fade (conv ~0).
  = powered gene-level resolution of H1(fade, dominant) vs H3(WT catch-up, biosynthetic-only).
- **PART B (BOUNDED, Issue #5 ceiling holds):** broadened TF panel (biogenesis/E2F/ESR1/
  lipogenic/data-driven mito 45-set) -- NO axis reliably absorbs (all delta_b_int boot CIs
  straddle 0); ESR1 doesn't even associate (lrt_p 0.83). tf_e2f numerical standout (b_int
  -0.279->-0.028, ~90% point est, NOT sig CI[-0.64,2.17]) = hint not result. Compositional
  DIAGNOSTIC: OXPHOS fade Myc+ -0.91; co-declining most = tf_mito_panel -0.60 / lipogenic -0.46 /
  E2F -0.36; ESR1 RISES +0.42. Proliferation/mito-TF fades WITH OXPHOS -> SUPPORTS (not proves)
  compositional dilution of the shrinking proliferative/TEB compartment.
- **PART C DEFERRED** (settle_it): reference deconvolution (Bach 2017 GSE106273 + MuSiC/Bisque,
  NOT on disk) / single-cell / E2F perturbation = the definitive test. CEILING: Part A
  descriptive of transcriptional change (composition confounds the fade term, not per-cell
  causal); Part B endogenous + n=6/group. Full outcome in revision-doc Issue #6 section.
  Artifacts: script 31 + results/attenuation_mechanism.rds + outputs/attenuation_mechanism/ (6 pdf).
- **FADE REFRAME (author Q, 2026-07-12, narration only -- no new analysis):** for the fade-led
  programs (both arms decline) the Myc+ slope = SHARED-developmental (=WT slope) + MYC-SPECIFIC
  (=the attenuation). A shared decline CANCELS in the genotype gap, so the WT decline gives
  CONTEXT but cannot explain the gap-shrink; the attenuation IS the Myc-specific excess. OXPHOS
  ~72-90% Myc-specific; proliferation ~100% (WT flat = purest de-amplification); MYC-core highest
  shared frac ~38% (= endogenous pubertal program winding down in WT, ties Issue #2 -> flip side
  of MYC-core "buffered"). Reframe = guardrail (don't read fade as Myc EXHAUSTION) + bridge (the
  excess = DE-AMPLIFICATION of the Myc-built proliferative/high-OXPHOS compartment -> composition).
  Folded into the figure-plan doc (Supp Fig B fade-narration note). Convergence programs are the
  opposite-sign case (WT rises to meet Myc), where this reframe does not apply.
