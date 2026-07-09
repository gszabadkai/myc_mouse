# Block A revision — step-by-step log

Branch: `BlockA-revision-step-by-step` (off the reviewed Block A commit `5b8fdf1`).

On reading `docs/2026-07-07_block_A_synthesis.md` the author found substantive
problems and chose a **depth-first, one-issue-at-a-time** revision: reframe the
question, add the exploratory analysis it needs, verify, commit, then move on. This
doc is the running log. Figures (the earlier 26/27/28 plan) are deferred until the
revision settles the science.

Scripts 13-25 and the synthesis remain on `analysis-exploratory` as the completed
Block A record; this branch revises forward from there.

Conventions unchanged: Option A (Claude writes numbered scripts, author sources them
in Positron; every script ends with an `if (FALSE)` sandbox); commit per verified
run; honest ceiling (bulk, n=6/group, survivor bias) stated on every new claim.

---

## Issue #1 — How does Myc integrate with / modify the background developmental program?

**Status:** analysis script written; awaiting Positron run. (2026-07-08)

### The problem (what was wrong)

`19_dev_composition.R` -> `results/dev_composition.rds`, `outputs/dev_composition/`
answered the Myc-vs-development question too superficially. Four defects:

1. **No integration question.** It only asked "does Myc push toward progenitor?"
   (AP2 genotype main effect). It never asked where the WT developmental trajectory
   is, how Myc+ is *already diverged at 6W*, or where 6W_pos -> 12W_pos lands.
2. **Consensus composites over-collapse** ~179 signatures into 5 lineage `_comp`
   colMeans via a coarse name heuristic, after a consensus that dropped "key sets" —
   averaging away cross-source disagreement and developmental-stage information.
3. **`diff_axis = HR_comp - LP_comp`** collapses the MASC->BASAL->LP->ALV->HR
   hierarchy + stage into one scalar.
4. **`(d6+d12)/2`** averages out the timepoint; the Myc effect must be resolved on the
   four-contrast trajectory (Myc-(6->12W), Myc+(6->12W), Myc@6, Myc@12 + interaction).

### The reframe (decisions locked)

Resolve **all 179 `MG_*` sets individually**, labelled by **source + MEC state +
modality**, across the **four-contrast trajectory**, using a **GSVA + fGSEA dual
lens**. Read the WT axis, Myc's 6W divergence, and the 12W_pos landing point.
**New script** (`26_dev_program_myc_integration.R`); scripts 19 and 25 left intact.

**State = consensus MEC nomenclature** (Gray/Kessenbrock/Khaled, Dev Cell 2025;
`docs/MEC_types.pdf`): three discontinuous resting-gland types **BMYO / LASP / LHS**
(not a hierarchy). Script keeps the 6-way `state_fine` for QC/overrides and collapses
to `state` (BMYO=MASC+BASAL+BASAL_PROG; LASP=LP+ALV; LHS=HR). See
[[mec-consensus-nomenclature]]. Split on 179 sets: 34 BMYO / 35 LASP / 34 LHS / 76 other.

### Grounding (already-scored; no upstream re-run)

- 179 `MG_*` GSVA-scored per-sample (`gsva_scores.rds$scores`).
- `gsva_overview.rds$coef_table`: `beta_time`/`myc_slope`/`d6`/`d12`/`beta_int` per set.
- `fgsea_percategory.rds$fgsea` (cat `03_mammary_development`): NES + padj for 167 MG_
  sets x 5 rankings (`timepoint_neg`,`timepoint_pos`,`myc_6W`,`myc_12W`,`interaction`).
- Contrast crosswalk: WT 6->12W = beta_time/timepoint_neg; Myc+ 6->12W =
  myc_slope/timepoint_pos; Myc@6 = d6/myc_6W; Myc@12 = d12/myc_12W; interaction.

### Deliverable

`scripts/26_dev_program_myc_integration.R` (Parts 0-D): set annotation table;
four-contrast dual-lens matrix; WT-axis + Myc-divergence geometry
(accelerate/divert/reverse per set); differentiation-continuum profiles (replaces the
HR-LP scalar); per-source small-multiple trajectories + master heatmap.
Writes `results/dev_program_myc_integration.rds`, `outputs/dev_program_myc_integration/*.pdf`.

### Downstream dependency (later issue)

Script 25's composition-death coupling consumes `dev_composition.rds` composites
(HR/LP/MASC). If Issue #1's richer view supersedes them, re-point 25 in a later issue.
Do not touch 19/25 now.

### Rev 2 (2026-07-08) — annotation-driven

After the consensus run, the author hand-curated the full annotation ->
`data/dev_mec_annotation.csv` (canonical, committed): `state` = final MEC type
(BMYO/LASP/LHS) or `other` (via `new_state` demote/promote); `state_fine` = analysis
subcategory (mains: source/lineage; `other`: SIG / MATRIX / IMMUNE / PAL2017 / HENRY /
SCHEELE_pub / CHUNG_C6 / FETAL + GRAY directional families); `lineage` = consensus
lineage kept even when demoted (for the nets). Script 26 rewritten to consume it
(algorithmic tagger dropped): convergence rho on the 3 MAIN states only (`other`
plotted as squares for observation); CHUNG ATAC OPEN-CLOSED nets added alongside the
GRAY UP-DN nets (lineage-tagged); `other` subdivided by `state_fine` (subgroup_profile
+ one trajectory PDF per subgroup). Final counts: 26 BMYO / 29 LASP / 29 LHS / 95 other
(14 subgroups). `data/dev_state_overrides.csv` retired (superseded).

### Outcome (2026-07-08, sourced clean in Positron; committed)

Three independent lenses converge. **The WT (Myc-) 6W->12W substrate loses luminal
(LASP/LHS) programs while BMYO stays FLAT -- a relative rebalancing by luminal
DECLINE, NOT a basal expansion** (BMYO wt_d=0.07, p=0.90; LASP d=-0.60, LHS d=-0.49,
~79% of sets down, per-set t p~1e-5 but sample-level composite ns at n=6/tp; consistent
with the literature -- no reported pubertal->adult BMYO expansion). Myc bends OFF that
substrate: convergence rho(Myc,WT) = -0.08 at 6W (orthogonal) -> -0.44 at 12W
(oppositional). Per state (state_stats): **BMYO suppression POWERED** (Myc effect
-0.15/-0.21, d~0.8-1.2, genotype main p=0.020); **LHS flip DIRECTIONAL** (-0.15/+0.17,
interaction p=0.089 -- trend, not the p-value but cross-modality agreement carries it);
LASP near-null. Corroborated by GRAY directional nets (Myc -> LESS-DIFFERENTIATED /
low-expression state + TEB/proliferative end-bud, strongest in BMYO) and CHUNG ATAC
OPEN-CLOSED nets (Myc closes basal chromatin, opens LP early, opens ML/LHS by 12W).
NOMENCLATURE NOTE: GRAY `HEVSLE` = **High-Expression vs Low-Expression** subgroups
within each lineage (Gray et al. 2023) = a more/less DIFFERENTIATION-state axis, NOT
high/low estrogen (the library catalog `mammary_dev_sets_catalog.md` mislabels it as
estrogen -- read it as expression/differentiation). Reading: at 6W a broad off-axis
identity suppression; by 12W a specific **anti-BMYO / pro-LHS reprogramming**. Ceiling:
main effects powered (24 samples); interaction directional at n=6/tp; composite p's
indicative (sets correlated); association not causation. Figure substrate for Block B.

Artifacts: `scripts/26_dev_program_myc_integration.R`,
`data/dev_mec_annotation.csv` (curated MEC annotation),
`results/dev_program_myc_integration.rds`, `outputs/dev_program_myc_integration/*`.

---

## Issue #2 — Endogenous Myc & the pubertal TEB/proliferative phenotype

**Status:** DONE (script 27, sourced clean, committed). 2026-07-08.

**Question.** It is reported that endogenous Myc is instrumental in the WT for the
pubertal proliferative phenotype (esp. TEBs). Quantify that endogenous role and show
the Myc+ transgene AMPLIFIES the same program. Reframe on already-scored GSVA (no
re-run). Six per-sample program composites: `myc` (MYC_signatures 17-set), `felsher`,
`hallmark_v2`, `prolif` (Proliferation 14), `myc_in_teb` (TFT_MYC_GRAY_*_TEB = MYC
targets in the Gray TEB context), `teb_ductal` (MG_TEB_VS_DUCTAL UP-DN). Decomposition:
ENDOGENOUS = WT 6W->12W temporal; TRANSGENE = Myc genotype main effect.

### Outcome

- **Endogenous Myc gates the pubertal program (WT).** Within WT, Myc activity tightly
  tracks proliferation (r=0.93), Myc-in-TEB (0.88), TEB-ductal (0.87), tightest at 6W
  (0.77-0.92). All six programs 6W>12W in WT (directional, ns at n=6/tp). Consistent
  with the literature.
- **Transgene amplifies it (POWERED).** Myc-target activity d=2.3-3.0 (p<1e-5),
  proliferation d=1.0 (p=0.02), Myc-in-TEB d=1.6 (p=6e-4); 6W_pos = the most TEB/
  proliferative state. The pure TEB phenotype is amplified SPECIFICALLY at 6W (6W_pos
  +0.41 vs 6W_neg +0.05; Myc boost ~4x larger at 6W than 12W -> pooled TEB geno p=0.14
  weak because time-concentrated), fitting TEBs as a pubertal structure (ties to the
  Issue #1 6W substrate).
- **Saturation nuance.** Within Myc+ the coupling loosens (r 0.87->0.55) and the
  Myc->phenotype slope flattens (1.65->0.83): endogenous Myc is rate-limiting in WT;
  the transgene exceeds that range -> amplification of LEVELS to a saturating high,
  not a proportional continuation of the endogenous dose-response. Interaction ns ->
  the transgene boost is ~additive over time (except TEB, 6W-concentrated).
- Ceiling: transgene powered (main effect); endogenous directional (n=6/tp);
  composite p's/correlations use correlated sets (indicative); endogenous-establishes-
  TEB is a literature claim -> WT co-variation + amplification shown, not causation.
  Hallmark_V1 not scored (V2 + Felsher + 15 targets used).

Artifacts: `scripts/27_myc_endogenous_amplification.R`,
`results/myc_endogenous_amplification.rds`, `outputs/myc_endogenous_amplification/*`.

---

## Issue #3 — Is the MYC effect REALLY "primarily mitochondrial"?

**Status:** DONE (script 28, sourced clean, committed). 2026-07-10.

### The problem (what was wrong)

The synthesis leads with "Myc preferentially amplifies mitochondrial biogenesis"
(theme B / model lead A), resting on AP6.1 (`20`, MitoCarta NES 2.18) and AP6.2 (`21`,
MitoCarta the top permutation-null z). The critique: mito NES is high, but the MYC
signature is ALSO high and spans many arms (ribosome biogenesis, proliferation,
metabolism) -- does mito have independent substance, or is "MYC effect = mitochondrial"
just restating "MYC effect = MYC targets"?

### The reframe (decisions locked)

The one word "primarily" does two jobs; split them.
- **Q1 PREFERENTIAL?** AP6.2's mito-#1 rank is size-confounded: z inflates with set
  size (MitoCarta n=1041 z=19.3 vs OXPHOS n=185 z=10.9 at ~equal effect). Report the
  size-fair metrics -- effect MAGNITUDE and housekeeping-corrected RELATIVE EXCESS over
  the matched null -- and add a powered per-sample Cohen's-d ranking across all arms
  incl. central metabolism.
- **Q2 CENTRAL?** Not "is mito enriched" but "is mito a CORE, phenotype-coupled arm or
  a bystander riding MYC dose?" Couple each mito / metabolic axis to the phenotypic
  OUTCOMES already measured -- MB2 tumorigenic fork (AP7/`18`), mito death priming
  (`23`), TEB/dedifferentiation (Issue #1/#2), proliferation -- then test whether the
  coupling SURVIVES removing the generic MYC axis (PARTIAL correlation vs the Felsher
  core = the 67-gene MYC phenotype core with mito stripped: 61 genes, 8% mito, 0%
  OXPHOS). Central = coupled AND survives deconfounding. Metabolic axes added per the
  cancer-metabolism remit. Full reframe on already-scored GSVA + saved composites.

### Deliverable

`scripts/28_myc_mito_centrality.R` (Parts 1-D): per-sample panel (3 mito arms:
all/OXPHOS-core/biogenesis; 3 MYC identity; proliferation; 12 central metabolic axes;
3 outcomes). Part A1 = reframe the AP6.2 null (magnitude + relative excess, no
recompute); A2 = powered genotype Cohen's d ranking. Part B = coupling (axis x outcome,
overall + by genotype + by timepoint) + partial correlation | Felsher core / prolif +
metabolic-centrality summary. Writes `results/myc_mito_centrality.rds`,
`outputs/myc_mito_centrality/*.pdf`.

### Outcome (2026-07-10, sourced clean in Positron; committed)

**Q1 -- "primarily mitochondrial" does NOT hold as literally-#1 on either lens.**
Size-fair matched-null excess (6W): MYC-targets 0.65 ~ OXPHOS-core 0.64 >
MitoCarta-wholesale 0.51 (mito #1 by the size-confounded z but #3 by excess) > E2F >
Metabolism > Proliferation; MYC-targets OVERTAKE at 12W (0.69 vs OXPHOS 0.44 vs mito
0.29). Powered genotype Cohen's d: hallmark_v2 3.0 > felsher 2.7 > mito_biogenesis 2.6 ~
mito_all 2.57 > amino-acid 2.35 ~ glutamine 2.34 ~ myc_sig 2.32 > mito_oxphos 2.11 ~ tca
2.09 ~ one-carbon 2.05 ~ glycolysis 2.01 ~ nucleotide 1.99 ~ oxphos_met 1.94 > ... >
prolif 1.01 > redox 0.08. Mito is a BROAD top tier with MYC-targets ahead and central
metabolism co-equal -- NOT uniquely the biggest arm. -> soften to "the OXPHOS/biogenesis
core is top-tier, co-equal with the MYC-target core," not "the single most preferential."

**Q2 -- centrality is REAL, ROBUST, and specific to the OXPHOS + biosynthetic core.**
Raw couplings are all high (everything rides MYC dose); the test is what survives
partialling out MYC dose. Deconfounded against BOTH the Felsher core AND the full 17-set
MYC composite (the strongest proxy), and cross-checked within WT (no transgene varying):
**10 ROBUST central couplings, all in the bioenergetic-biosynthetic core** --
`mito_oxphos` -> prolif (raw 0.86 / |fullMYC 0.59 / WT 0.89), teb-dediff (0.67/0.44/0.86),
priming (0.57/0.32/0.34); `nucleotide` -> prolif (0.85/0.54/0.92), priming
(0.63/0.49/0.53), teb (0.63/0.32/0.91); `oxphos_met` (independent set source) + `tca` ->
prolif (0.37, 0.32). The mito BIOGENESIS/translation arm and MitoCarta-WHOLESALE are
BYSTANDERS -- coupling collapses under both proxies (biog->teb 0.52->-0.03; ->prolif
0.74->0.11). The MB2 tumorigenic fork collapses for ALL mito arms (fork = MYC-dose-driven
for mito) but is robustly tracked by a DISTINCT metabolic signature: `cholesterol`/
mevalonate (+0.40) and `redox` (-0.50, even within WT; redox barely moved by Myc, d=0.08
-> a coupled-but-not-a-target axis). Dissociation: biogenesis tracks MYC-signature dose;
OXPHOS + nucleotide + TCA co-vary with the phenotype along an axis ORTHOGONAL to
MYC-signature activity and detectable within WT -- a sharper, better-supported version of
"mitochondria integrate (not merely read out) oncogenic + metabolic programs."

**Interpretation ceiling (critical -- baked into script notes + the Q2 figure).**
Surviving the partial means the phenotype signal is ORTHOGONAL to the MYC-signature axis,
NOT that it is MYC-INDEPENDENT or causal. Proxy under-capture, a non-linear MYC route, and
a common cause (e.g. cell-state composition) all survive too; DIRECTION is unidentified.
rho_pmyc + within-WT tighten "beyond MYC dose" but cannot establish independence or
mediation. Co-variation says NOTHING about NECESSITY (a MYC-readout arm may still be
required). Q1 powered (24 samples); Q2 per-sample association at n=6/group, correlated
sets; priming mito-defined (mito<->priming partly circular; nucleotide<->priming is the
non-circular corroborator). Definitive test is genetic. Association layer -> Block B.

Paper-ready framing: *Myc raises all mitochondrial arms, but their coupling to the
tumorigenic phenotype is dissociable -- ribosome/mtDNA biogenesis tracks MYC-signature
dose, whereas OXPHOS (with nucleotide synthesis and the TCA cycle) co-varies with
proliferation, dedifferentiation and death priming along an axis orthogonal to
MYC-signature activity and detectable within wild-type tissue. This identifies the
respiratory/biosynthetic core, not mitochondrial biogenesis per se, as the arm most
tightly linked to the phenotype. Whether this reflects a MYC-independent input or is
required for tumorigenesis is not resolved by these associations.*

Artifacts: `scripts/28_myc_mito_centrality.R`, `results/myc_mito_centrality.rds`,
`outputs/myc_mito_centrality/*` (q1_preferentiality_reframed, q1_genotype_effect_ranking,
q2_coupling_heatmap, q2_partial_correlation, metabolic_centrality).

### Downstream flag (later issue)

The mito-arm decomposition (OXPHOS core vs biogenesis/translation) and the "central
metabolic axes" (nucleotide/TCA/cholesterol/redox) are new framings not in the synthesis
model. When the synthesis / model-lead wording is revised (parked framing decision),
replace "primarily mitochondrial" (theme B / lead A) with the dissociable-arms claim
above. Do not touch the synthesis doc now.
