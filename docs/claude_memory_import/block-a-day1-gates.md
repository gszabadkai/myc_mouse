---
name: block-a-day1-gates
description: "Block A Day-1 gate results (SHRINKS + DECISION_POINT), T3 resolution, and the open trajectory reframing"
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
---

Block A Day-1 decision layer ran 2026-07-05 (scripts 13, 14). Full record:
`docs/2026-07-05_gates_note.md`.

- **Gate 1 = SHRINKS, and it is BIOLOGICAL.** Genotype (Myc+/-) effect: 2777
  genes padj<0.1 at 6W vs only 239 at 12W (-91%). Interaction predominantly
  negative (Myc effect weakens at 12W). This kills BOTH the "permission/licenses"
  verb and a "progressive divergence" trajectory. **Power caveat RESOLVED
  2026-07-05** (script 13 PART 2b effect-size check): genotype |LFC| point
  estimates (SE-independent) attenuate to ~45% of 6W by 12W on a fixed
  6W-divergent set (slope 0.46, 97% of genes attenuate), and the LFC-INDEPENDENT
  expressed set gives the SAME slope (0.45) - so it is NOT a selection/power
  artifact. The genotype effect genuinely ~halves (direction preserved, not sign
  reversal); "convergence" = ~2-fold attenuation, not to zero. H2 (power
  artifact) largely rejected. The MECHANISM (H1 front-loaded / H3 WT catch-up /
  H4 selection) is still OPEN pending Day 2+. Outputs:
  outputs/gates/gate1_effect_size_check.{csv,pdf}.
- **Gate 2 = DECISION_POINT.** Apoptosis-PRO module coordinately shifted negative
  (one-sample t p=0.023); ANTI not (p=0.25); Bbc3 an outlier but PRO significant
  even so. Apoptosis arm STAYS in main-figure candidates; cell death NOT demoted
  to Supp 1. Small-n caveat (23/7) - corroborate with cell-death branches.

**T3 RESOLVED - do not re-litigate:** Gate 2 direction + PRO/ANTI were computed
on raw/Wald values (interaction_results$interaction_raw, unshrunken MLE; IHW
touches only padj). No shrinkage recompute anywhere in the gate layer. See also
CLAUDE.md "Shrunken vs raw LFCs".

**Day 2 DONE (2026-07-06, commit 9939a37).** Script 15 GSVA engine +
script 16 raw cell-death port written and sourced clean in Positron; GSVA added
to 00_setup_packages.R. Script 15 scores all 902 GSVA-tagged library sets
(method gsva/both) per-sample from VST (log, kcdf=Gaussian - OPPOSITE of mitoPPS
linear); output results/gsva_scores.rds (list: scores/set_meta/sample_meta).
Script 16 = faithful port of archive 09 reading combined_df_annotated_raw.rds
(the raw df uses _raw-suffixed LFC cols + has no `gene` col; port aliases +
drops accordingly), outputs suffixed _raw under outputs/cell_death_raw/.

**Day 3 pre-step WRITTEN (2026-07-06), awaiting Positron run.** The author-requested
GSVA visualisation became the new **script 17 `17_gsva_overview.R`** (build spec
17-24 renumbered +1 -> 18-25; doc-only, no files moved). It reads gsva_scores.rds
and per set fits `lm(score ~ timepoint*myc_status)` (vectorised OLS): beta_int =
d12-d6 = the Gate-1 attenuation (RANK metric); beta_time = Myc-/WT slope; myc_slope
= Myc+ slope. This is the **H1-vs-H3 discriminator**: interaction ranks the movers,
but the 4-point two-line profile shows WHICH genotype trajectory moved (H1 = Myc+
descends / front-loaded; H3 = Myc- ascends / WT catch-up). Outputs per category:
`<cat>_landscape.pdf` (top~35 by interaction p, row-z-scored, WT-slope+interaction
row annots), plus `<cat>_profiles.pdf` + `<cat>_dumbbell.pdf` (d6->d12) for the 4
trajectory categories (03_mammary_development, 02_myc_signatures,
07_biogenesis_discrimination, 08_apoptosis). Writes results/gsva_overview.rds
(coef_table + cond_means). Author sources 17 next; report before Day-3 proper.

**GSVA trajectory validation (script 17 PART 6, commit 3c3058c) - READ THIS
BEFORE writing any GSVA trajectory claim.** Added CAMERA (inter.gene.cor=NA) +
ROAST + category-union CAMERA + a cross-check of GSVA beta_int vs mean DESeq2
interaction Wald per set. Outcome:
- **Nothing survives FDR under ANY method** - per-set camera/roast/lm and
  per-category CAMERA all 0 significant (min FDR ~0.93). Individual-set AND
  category-level interaction significance is NOT supported.
- **Inter-gene correlations 0.25-0.36** (camera estimate) confirm the library
  sets are strongly non-independent -> per-set "consistency across many sets" was
  overlap-inflated, as suspected.
- **Cross-check STRONG: Spearman rho=0.663, p<2.2e-16 overall**, sign-consistent
  per category (Apoptosis 0.83, intersections 0.81, Mammary_dev 0.80, ... Biogenesis_
  discrimination lowest 0.29). Count-level interaction reproduces BOTH the
  biogenesis/myc attenuation (negative) AND the mammary-lineage crossover
  (positive). So the trajectory DIRECTIONS are real/method-independent (not a GSVA
  artifact) - but note same 24 samples, so this kills "method artifact", not
  "sampling noise".
- **Verdict:** global attenuation is solid (Gate-1 effect-size, was significant);
  the per-program H1/H3 decomposition is DIRECTIONAL-ONLY (corroborated in sign,
  not individually significant). The H1-vs-H3 biogenesis split (based on tiny non-
  sig beta_time vs myc_slope signs; biogenesis has the LOWEST cross-check rho) is
  the weakest part - report as illustrative, not a claim. Best single directional
  lead = Apoptosis CDC_PRODEATH_CICD (all 4 methods agree, mean_int_stat -1.06,
  int_p 0.036 / camera 0.037 / roast 0.020, not FDR-sig). For figure text: say
  "directionally consistent / method-independent, not significant at n=6/group",
  drop per-set p-values from the narrative.

**AP6.2 permutation null (script 21, commit e370d1e) - the HARD 'why mitochondria'
result.** Expression x dispersion-matched null on mean |raw LFC| defeats the
housekeeping/expression confound behind AP6.1's fGSEA. Effect size (obs-null,
n-indep, myc_6W): OXPHOS 0.201, MYC_TARGETS 0.190, MitoCarta 0.164 >> generic
Metabolism 0.095 >> Proliferation 0.032, Mammary 0.023. Mito/OXPHOS MOST
preferentially altered (~2x metab, ~5x prolif) AFTER removing the expression
advantage = real preferential targeting. KEY: proliferation's fGSEA enrichment
(NES 2.08) does NOT survive the magnitude-matched null (0.032, NS at 12W) -> it
was the expression confound; the preferential alteration is specifically
metabolic/mitochondrial, not generic proliferation ('preferential not just
absolute'). 12W: same ranking, attenuated. CAVEATS: all p_emp hit floor (Myc broad
+ null only removes expression) -> claim is RELATIVE mito>>comparators; raw z
inflated by inter-gene correlation -> use EFFECT-SIZE ranking not absolute z.

**Per-category fGSEA (script 20, commit 2f278c9) - AP6.1.** One fGSEA per
(ranking x category), BH within category; 5 unshrunken Wald rankings.
- Myc genotype (myc_6W/12W): strong coherent enrichment - MitoCarta NES 2.18 (top
  padj 1e-55), biogenesis 2.65, Myc-sig 2.40, prolif 2.08, TF 2.24, metab 1.62 UP;
  mammary_development -1.39 DOWN. Strong at BOTH timepoints (fGSEA reads rank/
  direction; magnitude attenuates per Gate 1 - the CLAUDE.md 'fGSEA unaffected').
- ATTENUATION via CONVERGENCE (powered corroboration): timepoint_pos (Myc+ 6W->12W)
  mito/biogenesis DOWN (-1.5) vs timepoint_neg (WT) mito UP (+0.8) -> gap closes.
  **CORRECTION (Issue #4, 2026-07-11): the "WT mito UP +0.8" here is a MEDIAN over 63
  heterogeneous MitoCarta sets, carried by the BIOSYNTHETIC arm; the actual OXPHOS
  pathways FALL in WT with age (MITOCARTA_OXPHOS timepoint_neg -1.94, OXPHOS_SUBUNITS
  -2.47) and fall further in Myc+ (-2.61). There is NO WT developmental OXPHOS rise. Do
  NOT cite this aggregate as "WT mito rises"; report pathways individually. See
  [[issue4-attenuation-scoped]].**
- AP6 preview: MitoCarta 2.18 > metabolism 1.62, >> Hallmark -1.18, ~ proliferation
  2.08 => preferential vs GENERIC comparators (prolif co-equal); hard test = perm
  null script 21.
- **CORRECTION to the 'interaction floor' framing:** fGSEA-on-interaction reports
  461 sig pathways, BUT fGSEA assumes gene independence and the sets are co-regulated
  (corr 0.25-0.36, measured in 17) -> ANTI-CONSERVATIVE / inflated. Correlation-aware
  CAMERA-on-interaction (17) found ~nothing; that verdict stands. Direction (biogenesis
  interaction NES negative = attenuation) is consistent; the SIG COUNT is not
  trustworthy. Strong Myc-genotype effects robust regardless (huge + Gate-1 backed).
  Precise statement: interaction floor holds for correlation-aware/per-set tests;
  fGSEA over-calls it.

**Dev composition (script 19, commit 6f9c128) - AP1/2/3/4.**
- AP2 POWERED: Myc drives biogenesis composite genotype p=3.7e-6; SUPPRESSES basal
  p=0.016. Luminal (HR/LP) pooled null = the 6W-suppress/12W-elevate crossover (17);
  diff_axis (HR-LP) null (HR/LP move together -> wrong instrument for the LP claim).
- AP4: Felsher(Myc activity) ~ biogenesis composite rho=0.90 => the Myc footprint IS
  a biogenesis program, not lineage (felsher~diff_axis 0.08).
- AP3 descriptive: Myc lineage effect ANTI-parallel to WT developmental direction
  (rho=-0.28 p=1e-4; basal -0.66, HR -0.35) => Myc holds lineage differentiation
  AGAINST development while accelerating biogenesis WITH it (17).
- AP1: WT 6W->12W mitoPPS reallocation - mtDNA-encoded OXPHOS UP (+0.52), nuclear
  OXPHOS subunits DOWN, FAO/carnitine UP = the moving mito substrate (ties to
  [[myc-death-timing-question]]).
- The AP7 LP link lands via biogenesis drive + basal suppression, NOT a clean
  LP-composite rise - do not overclaim a literal LP-marker shift.
- Lineage tagging coarse (109/179 OTHER = cluster/pseudotime sets); powered results
  don't depend on it. Diagnostic plots in outputs/dev_composition/ (7). Polished
  narrative figures deferred to the figure scripts (renumbered 23/24/25) post-review.

**AP7 MB-fork projection (script 18, commit a37ae0d) - the CENTREPIECE, and the
first POWERED/significant result in Block A.** Projects mouse onto Menegollo/
Bentham 2024 MCbiclust forks (docs/Menegollo_et_al_2024.pdf; sets already mouse-
mapped + GSVA-scored in 15). Fork topology: MB1&MB2 share lower fork MB12_LF and
diverge into MB1_UF (mature-luminal, non-Myc) vs MB2_UF (LP, ER-neg, Myc/miRNA,
biogenesis-high); MB3 independent. Score = GSVA(UF=group1) - GSVA(LF=group2).
- Treatment VALIDATED: UF/LF anticorrelated (MB1 -0.77, MB2 -0.84); MB12_LF r=0.98.
  MB3 not a switch in mouse (r=0.18) -> control only.
- **STRONG (powered main effect): Myc drives the biogenesis-high upper fork** -
  MB2_score genotype p=0.0045, MB1_score p=0.0045, robust across HI_CV/BICLUSTER/
  METAB. Attenuates 6W->12W (0.62->0.34, mirrors Gate 1). Direct cross-species
  validation of the mito-biogenesis thesis; genotype MAIN EFFECTS are powered
  (unlike the interaction floor).
- **MB2-specificity SUPPORTED but MODEST:** MB2_over_MB1_UF (=MB2_UF-MB1_UF)
  Wilcoxon p=0.008, lm p=0.026, small effect (0.08 vs 0.8 shared drive); holds
  HI_CV/BICLUSTER NOT METAB => the Myc-specific MB2_UF identity is in REGULATORY
  not metabolic genes (metab genes can't tell MB1_UF from MB2_UF - both just
  biogenesis-high). Myc leans toward LP/ER-neg over mature-luminal (ties to 17's
  luminal-suppression). Falsifiable test PASSES; don't oversell as "Myc->MB2 only"
  - shared biogenesis drive dominates.
- Manuscript: Tier-1 = Myc drives the biogenesis upper-fork switch (validated,
  significant, robust, attenuates); Tier 1-2 = modest MB2_UF/LP preference. Named
  fig mb2_over_mb1_boxplot.pdf. See [[myc-death-timing-question]] (MB2_UF LP/ER-neg
  aggressive state may connect to the 6W death-permissive state).

**Trajectory MECHANISM resolved (script 17 PART 6b/6c, commit e5c8201) - two
informative NEGATIVES; do not re-run these hoping for significance.**
- **Between-category permutation contrast NON-significant** (biogenesis vs mammary
  p=0.236; oncogenic-metabolic vs mammary p=0.184). Cause = the categories are NOT
  a gene partition: Biogenesis_discrimination 469 genes -> only 43 EXCLUSIVE of
  mammary (~90% overlap), and Mammary_development union ~9024 genes (~whole
  transcriptome). So the test collapsed to tiny-set-vs-genome. Lesson: biogenesis
  genes ARE developmental genes; the crossover is a set-composite phenomenon, not a
  clean gene-level dichotomy. A cleaner biogenesis-vs-SPECIFIC-luminal-signature
  contrast is possible but won't change the n=6 significance floor - low priority.
- **Per-sample mechanism test REFUTES 'development inhibits/absorbs Myc'** (the
  window-closing hypothesis): bio~dev composites POSITIVELY correlated (+0.49 neg,
  +0.31 pos - NOT anti-correlated); adjusting for dev does NOT shrink the Myc+
  biogenesis 6W->12W decline (-0.075 -> -0.143, grows). => Development CO-DRIVES
  biogenesis (shared maturation axis). The attenuation = WT developmental
  CONVERGENCE onto a shared biogenesis program + modest fade of Myc's front-loaded
  head start. Development is NOT a culprit/antagonist. "Constant driver, moving
  substrate" holds but substrate CONVERGES WITH Myc, not opposes it.
- **Significance floor now confirmed across FOUR methods** (lm / CAMERA / ROAST /
  permutation): NO interaction test significant at n=6/group. Powered claims =
  levels only (Gate 1). NEW defensible Tier-1 claim = the clean NEGATIVE: the
  attenuation is not developmental antagonism (biogenesis is a shared developmental
  program). Manuscript: assert constant drive + global attenuation + this negative;
  keep the program decomposition descriptive; DROP the H1/H3 split and the
  'development closes the oncogenic window' story (data refute it).

**Interpretive thread (author, 2026-07-06) - biogenesis/death coupling:** the
`Biogenesis_apoptosis_intersections` category (26 sets) was added to script 17's
trajectory panels for a specific reason - it discriminates ER/PGC1a-driven
biogenesis (which appears to COUPLE to cell death in the experimental setup) from
Myc-driven biogenesis (which probably does NOT). Read alongside the `Apoptosis`
category: if intersection programs are Myc+ elevated but their d6->d12 trajectory
does not track the apoptosis-PRO shift (Gate 2), that supports "Myc biogenesis is
uncoupled from death" - biogenesis without the PGC1a/ER death penalty. This shapes
Gate 2 interpretation and the Day-3 dev-composition read. Ties to the companion
Menegollo/Bentham Cancer Res 2024 mito-decision-point theme.

**Next: Day 3 proper** (renumbered): script 18 AP7 MB-fork projection (centrepiece,
resolve fork-set species first line), 19 dev composition AP1/2/3/4, 20 per-category
fGSEA (BH within run, fresh Hallmark via msigdbr). Build spec:
`docs/2026-07-04_block_A_build_spec.md`. The MECHANISM of the Gate-1 attenuation
(H1/H3/H4) is still OPEN; 17 + 18/19 discriminate it.

**Block A analysis complete through script 21 (commit 418ee96); review doc
`docs/2026-07-06_block_A_review.md`.** Scripts 22 + 23 DONE (2026-07-07), before
figures. Figure scripts renumbered 23/24/25 -> 24/25/26.
- **Script 22 (commit afcddfe) mtDNA/abund/retention reframe.** AP-mtDNA: WT mtDNA-
  OXPHOS mitoPPS up (0.89->1.42) / nuclear down = reallocation (cf 19 AP1), mtDNA-
  vs-nuclear per-complex anticorrelated r~-0.5. AP-abund (IMPORTANT): fGSEA vs
  mitoPPS lenses split TIME-DEPENDENTLY - fGSEA absolute NES FLAT (2.18->2.23) but
  mitoPPS reprioritisation dispersion NARROWS (SD 0.14->0.10, median push
  0.052->0.015) = Myc selectively rearranges the mito compartment at 6W, uniform
  biogenesis by 12W. NOT the Gate-1 magnitude attenuation (fGSEA flat) - a distinct
  axis (allocation selectivity). AP-retention: mito/biogenesis retention_abs
  strongly negative (biogenesis -2.71) = attenuation-via-convergence per-category
  on powered temporal contrasts.
- **Script 23 (commit 2140b37) death-timing substrate model.** See
  [[myc-death-timing-question]] for the full result. Headline: Myc death-coupling
  is substrate-gated - priming p=0.038 (~3x at 6W), Gate2 p=0.023; MITONUCLEAR
  IMBALANCE peaks 6W_pos, death coupling r=0.75->0.07; the 12W decoupling IS the
  loss of killing. Promotes narrative D to a main-figure hook. Corrected: "WT
  de-primes with age" REFUTED (baseline flat).
- Review OPEN DECISIONS: #3 (script 22) + #4 (death-timing) DONE; #2 cell-death now
  RECOMMENDED hybrid (imbalance panel main + full death Supp). Still to lock: #1
  framing (A lead), #2 which main figure hosts the imbalance panel, #5 Block B
  figures (now 26/27/28).

**Scripts 24 + 25 DONE (2026-07-07, commits ef6b04a + b50e00a) - resolve the
review discussion (docs/BlockA_review_discussion.md); figure scripts renumbered
24/25/26 -> 26/27/28.** See [[myc-death-timing-question]] RESULT 2 for full detail.
- **Script 24 biogenesis discrimination.** IMBALANCE = mitoPPS ALONE (fGSEA is a
  SEPARATE importance axis - author correction; do NOT define imbalance via a
  fGSEA-sign x mitoPPS-sign quadrant, that is a bulk-vs-mito DISSOCIATION). Within-
  mito imbalance map at Level2-group resolution, mtDNA split out of the "OXPHOS
  subunits" node so the mitonuclear split is not masked. mtnuc imbalance peaks
  6W_pos (+0.56), inverts by 12W. WT substrate = maturation from nuclear-led
  ASSEMBLY (6W) to mtDNA-completed RUNNING metabolism (12W); Chaperones FALL ->
  UPR^mt/ISR REFUTED. Cat-7: Myc drives MYC and ER/PGC1a biogenesis identically
  (co-opts ER/PGC1a-axis TFs ESRRA/NRF1/GABPA; ESR1 antagonised). Cat-9: balanced
  CORE/MITO_NU death-coupling highest at 6W.
- **Script 25 developmental substrate of death (augments 23; 23 intact).** Two
  coupled substrate axes: MITO maturation (imbalance resolution ~88% mtDNA-driven,
  genotype-shared = developmental) + CELL-STATE (WT lineage shift underpowered n=6
  p>0.18, BUT per-sample stem/MASC content death-coupled 0.78 at 6W, decouples by
  12W). Part B CORRECTION: at 6W bio_comp 0.83 > MASC 0.78 > imbalance 0.67 couple
  to death - mito imbalance is NOT the sole death axis, stem content matters too
  (coupled state, not separable). Culling: sd_pro/sd_bio narrow 6W->12W, Myc+ more.

**BLOCK A ANALYSIS COMPLETE (scripts 13-25; HEAD 5b8fdf1 on analysis-exploratory,
2026-07-07).** CANONICAL reference now = `docs/2026-07-07_block_A_synthesis.md`
(full hypothesis->outcome->results/.rds/outputs map, the 4 method lenses +
attenuation-vs-reprioritisation-vs-imbalance-vs-dissociation distinctions, powered-
vs-directional split, the A-lead+D-bridge model). `docs/2026-07-06_block_A_review.md`
is SUPERSEDED for the inventory (its tiered claim-spine still valid). NEXT =
Block B figures (scripts 26/27/28) on a NEW `paper-figures` branch off the reviewed
Block A commit - do NOT start until the author has read the synthesis and locked
the open figure decisions (framing A-lead+D-bridge; which main figure hosts the
mitonuclear-imbalance panel; cell-death main-vs-Supp = recommended hybrid). Author
is reviewing the synthesis; paused 2026-07-07.

**BRANCH STATE UPDATE 2026-07-12 (git consolidation done):** Block A revision (Issues #1-6,
scripts 26-31) is COMPLETE + synthesized. `analysis-exploratory` was fast-forwarded to hold the
reviewed Block A (scripts 00-31); reviewed commit tagged **`block-a-reviewed` @ 4b82a82**.
**CURRENT working branch = `paper-figures`** (off the tag) for Block B publication figures.
`BlockA-revision-step-by-step` is subsumed (safety ref only). `new-analysis` still = 00-12 trunk.
Wherever notes below say the current branch is `BlockA-revision-step-by-step`, it is now
`paper-figures`. NEXT = build 2 main + supplementary figures (re-render Block A outputs, do NOT
re-run 26-31); figure-content partition AWAITS the author's narrative. Plan:
`docs/2026-07-12_BlockA_revision_synthesis_and_figure_plan.md`. See CLAUDE.md "Branches".

**PIVOT 2026-07-08: figures DEFERRED; step-by-step Block A REVISION opened.** On
reading the synthesis the author found substantive problems and chose a depth-first,
one-issue-at-a-time revision on branch `BlockA-revision-step-by-step` (off `5b8fdf1`).
Living log: `docs/2026-07-08_BlockA_revision_plan.md`. Do NOT start figures / the
`paper-figures` branch; the framing decisions above are PARKED (decide-later).

**Issue #1 DONE (script 26 `26_dev_program_myc_integration.R`) - Myc x developmental
program; supersedes script 19's superficial answer.** Reframe on already-scored data
(no re-run): all 179 MG_* dev sets resolved individually across the 4-contrast
trajectory (WT_6W->12W / Myc@6W / Myc@12W / Myc+_6W->12W + interaction), dual lens
GSVA (per-sample state) + fGSEA (vs transcriptome; WT-axis crosswalk rho 0.92).
- **State = Gray/Kessenbrock/Khaled 2025 consensus MEC nomenclature** (see
  [[mec-consensus-nomenclature]]) via the author-curated `data/dev_mec_annotation.csv`
  (state final BMYO/LASP/LHS/other; state_fine = analysis subcategory; lineage kept
  even when demoted). Counts 26 BMYO / 29 LASP / 29 LHS / 95 other (14 subgroups).
- **Result (3 independent lenses converge).** **WT (Myc-) 6W->12W substrate: luminal
  (LASP/LHS) DECLINE with BMYO FLAT - a relative rebalancing by luminal loss, NOT a
  basal expansion** (state_stats: BMYO wt_d 0.07 p=0.90; LASP d-0.60, LHS d-0.49, ~79%
  sets down, per-set t p~1e-5, sample-level composite ns at n=6/tp; MATCHES literature
  - no pubertal->adult BMYO expansion; do NOT claim one). Myc bends OFF that substrate:
  convergence rho(Myc,WT) -0.08 at 6W (orthogonal) -> -0.44 at 12W (oppositional). By
  state: **BMYO suppression POWERED** (Myc effect -0.15/-0.21, d~0.8-1.2, genotype main
  p=0.020 - matches script 19 basal p=0.016); **LHS flip DIRECTIONAL** (-0.15/+0.17,
  interaction p=0.089 trend, credible via cross-modality not the p); LASP near-null
  (Myc main p=0.81). The pooled main-effect null for LHS/LASP is the flip CANCELLING
  over time - must read timepoint-resolved (the point of dropping (d6+d12)/2). LHS flip
  corroborated by GRAY directional nets (Myc -> LESS-DIFFERENTIATED/low-expression +
  TEB/proliferative end-bud, strongest in BMYO) AND CHUNG ATAC OPEN-CLOSED (Myc closes
  basal chromatin, opens LP early, opens ML/LHS by 12W). At 6W = broad off-axis
  identity suppression; by 12W = specific anti-BMYO / pro-LHS reprogramming. NOTE: GRAY
  HEVSLE = High-Expression vs Low-Expression subgroups (Gray 2023) = differentiation-
  state axis, NOT estrogen (the library catalog mislabels it "high/low estrogen").
- Honest ceiling: main effects powered (24 samples); interaction directional at n=6/tp;
  composite p's indicative (sets correlated r0.25-0.36 - not independent replication).
  Association not causation. Output = figure substrate for Block B.

**Issue #2 DONE (script 27 `27_myc_endogenous_amplification.R`) - endogenous Myc &
the pubertal TEB/proliferative phenotype.** Reframe on already-scored GSVA (no re-run);
6 program composites (myc=MYC_signatures 17-set, felsher, hallmark_v2 [V1 NOT scored],
prolif=Proliferation 14, myc_in_teb=TFT_MYC_GRAY_*_TEB, teb_ductal=MG_TEB_VS_DUCTAL
UP-DN). Decomposition: endogenous=WT 6W->12W temporal; transgene=Myc genotype main.
- **Endogenous Myc GATES the pubertal program in WT** (Myc~prolif r0.93, ~myc_in_teb
  0.88, ~teb_ductal 0.87 within WT, tightest at 6W). All 6 programs 6W>12W in WT
  (directional, ns n=6). Supports the literature claim endogenous Myc drives the
  pubertal proliferative/TEB phenotype.
- **Transgene AMPLIFIES it (POWERED):** Myc-target activity d2.3-3.0 p<1e-5, prolif
  d1.0 p0.02, myc_in_teb d1.6 p6e-4; 6W_pos most TEB/proliferative. TEB phenotype
  amplified SPECIFICALLY at 6W (6W_pos +0.41 vs 6W_neg +0.05, ~4x boost vs 12W -> pooled
  TEB geno p0.14 weak b/c time-concentrated; TEBs = pubertal structure, ties Issue #1
  6W substrate).
- **SATURATION nuance:** within Myc+ coupling loosens (r0.87->0.55) + slope flattens
  (1.65->0.83) = endogenous Myc rate-limiting in WT; transgene exceeds range ->
  amplifies LEVELS to saturating high, not a proportional dose-response. Interaction ns
  (additive boost). Ceiling: transgene powered, endogenous directional, correlated
  sets -> indicative, co-variation+amplification not causation.
- Committed on BlockA-revision-step-by-step. Artifacts: script 27 +
  results/myc_endogenous_amplification.rds + outputs/myc_endogenous_amplification/.
- Committed with `data/dev_mec_annotation.csv` + revision doc.

**Issue #3 DONE (script 28 `28_myc_mito_centrality.R`) - is the MYC effect REALLY
"primarily mitochondrial"? REFINES/SOFTENS the theme-B/lead-A synthesis claim + the
AP6.2 "mito MOST preferentially altered" phrasing above.** Reframe on already-scored
GSVA + saved composites (no re-run). Splits "primarily" into Q1 preferential / Q2
central.
- **Q1 PREFERENTIAL = overclaim, softened.** AP6.2's mito-#1 rank is a SIZE artifact of
  the permutation-null z (z inflates with set size: MitoCarta n=1041 vs OXPHOS n=185).
  Size-fair relative-excess (obs-null)/null: MYC-targets 0.65 ~ OXPHOS-core 0.64 >
  MitoCarta-wholesale 0.51 (mito #1 by z, #3 by excess); MYC-targets OVERTAKE at 12W
  (0.69 vs 0.44 vs 0.29). Powered genotype Cohen's d: MYC identity biggest (hallmark_v2
  3.0, felsher 2.7) > mito arms 2.6 ~ amino-acid/glutamine/TCA/glycolysis/one-carbon/
  nucleotide (2.0-2.35, CO-EQUAL with mito) > ... > prolif 1.0 > redox 0.08. Mito is a
  BROAD top tier, MYC-targets ahead, central metabolism alongside - NOT uniquely biggest.
  -> say "the OXPHOS/biogenesis core is top-tier, co-equal with the MYC-target core," NOT
  "the single most preferential." (The AP6.2 effect-size numbers above already hinted this
  - OXPHOS 0.201 ~ MYC_TARGETS 0.190 > MitoCarta 0.164; Issue #3 formalises it.)
- **Q2 CENTRAL = real, robust, SPECIFIC to the OXPHOS + biosynthetic core.** Method: raw
  coupling of each mito/metabolic axis to phenotypic OUTCOMES (MB2 fork, mito death
  priming, TEB/dedifferentiation, proliferation), then PARTIAL correlation deconfounding
  against the mito-stripped Felsher core (rho_pcore) AND the full 17-set MYC composite
  (rho_pmyc, strongest proxy), cross-checked within-WT (rho_wt). 10 ROBUST central
  couplings (survive BOTH proxies + present within WT), ALL bioenergetic-biosynthetic:
  mito_oxphos -> prolif (0.86/|fullMYC 0.59/WT 0.89), teb (0.67/0.44/0.86), priming
  (0.57/0.32/0.34); nucleotide -> prolif (0.85/0.54/0.92), priming (0.63/0.49), teb
  (0.63/0.32/0.91); oxphos_met (independent set) + tca -> prolif. Mito BIOGENESIS/
  translation + MitoCarta-WHOLESALE = BYSTANDERS (collapse under both proxies). MB2 fork
  collapses for all mito arms (fork = MYC-dose-driven) but is robustly tracked by a
  DISTINCT signature: cholesterol/mevalonate (+0.40) + redox (-0.50, within-WT; redox not
  a Myc target, d0.08). Dissociation: biogenesis tracks MYC-signature dose; OXPHOS +
  nucleotide + TCA co-vary with phenotype ORTHOGONAL to MYC-signature + within WT.
- **INTERPRETATION CEILING (critical - do NOT overclaim).** Surviving the partial =
  phenotype signal ORTHOGONAL to the MYC-SIGNATURE axis, NOT MYC-INDEPENDENT or causal.
  Proxy under-capture / non-linear MYC route / common cause (cell-state composition) all
  survive too; DIRECTION unidentified; co-variation says NOTHING about NECESSITY (a
  MYC-readout arm may still be required). Author's proposed "something MYC-independent
  regulates OXPHOS -> regulates the phenotype" is NOT licensed - stop at "orthogonal to
  MYC-signature activity, detectable within WT." Definitive test is genetic. n=6/group,
  correlated sets; priming mito-defined (mito<->priming partly circular; nucleotide<->
  priming is the non-circular corroborator). Paper framing + full numbers in the revision
  log Issue #3 outcome.
- Committed on BlockA-revision-step-by-step. Artifacts: script 28 +
  results/myc_mito_centrality.rds + outputs/myc_mito_centrality/ (5 pdfs). Harmless
  select.Rd Rd-parse warnings on source (namespace shadowing) - not a runtime error.
- Downstream flag: mito-arm decomposition + central-metabolic-axes are NEW framings; when
  the synthesis/model-lead wording is revised (parked), replace "primarily mitochondrial"
  with the dissociable-arms claim. NEXT = author sends Issue #4 (figures still DEFERRED).
