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

---

## Issue #4 — What IS the 6W->12W attenuation? (the "NES paradox")

**Status:** DONE (script 29, sourced clean, committed). 2026-07-11.

### The problem (what was fuzzy, and an error caught mid-scoping)

Gate 1 established the genotype effect SHRINKS ~2-fold 6W->12W and that it is biological
(effect-size, not power). But "attenuation" stayed fuzzy because two rulers get read as
one, and there was a live confusion: fGSEA NES of the Myc-vs-WT contrast is FLAT across
age yet the divergence shrinks (the "NES paradox"). During scoping the author caught a real
error in the working explanation: an aggregate "WT mito rises +0.82 with age -> converges
onto Myc" claim. That +0.82 is a MEDIAN over 63 heterogeneous MitoCarta fGSEA sets and is an
ARTIFACT of averaging opposite-direction arms. It is RETRACTED (see below).

### The resolution (two rulers, pathway-resolved)

- **NES is RANK-based / magnitude-blind** -- it asks "are these genes at the TOP of this
  contrast?" (yes at both ages). It measures identity/priority, NOT size. Using NES flatness
  to argue "no attenuation" is a category error; using an all-MitoCarta NES median hides
  opposite-direction pathways.
- **What SHRINKS is MAGNITUDE** (raw |LFC|). DE count 2777->239 (Gate 1); per-pathway here.
- **mtDNA is contained** (verified): the fGSEA MITOCARTA_OXPHOS set is already nuclear (zero
  mt-* genes); mtDNA bias lives only in the GSVA MITOCARTA_ALL composite and mitoPPS.
- **TPM skipped for a STATISTICAL reason** (not a data gap): within-gene between-group
  z-scores make gene length a per-gene constant that cancels, so TPM == log-CPM in shape.
  Part B uses VST + normalized-counts + log-CPM.

### Deliverable

`scripts/29_attenuation_decomposition.R` (reframe on fitted DESeq2 contrasts, GSVA, mitoPPS,
fGSEA; VST/log-CPM are deterministic transforms, no re-fit). Part A decomposition + two-ruler
NES-paradox figure, PATHWAY-RESOLVED. Part B ABSOLUTE nuclear-OXPHOS levels across 4 groups on
3 normalizations (the un-done check). Part C mitoPPS-vs-absolute reconciliation. Part D
first-pass within-WT regulator regression. Artifacts: `results/attenuation_decomposition.rds`,
`outputs/attenuation_decomposition/*` (A_nes_paradox_two_rulers, A_convergence_temporal_by_pathway,
B_nuclear_oxphos_heatmap_vst, B_absolute_composite_trajectory, C_mitopps_vs_absolute,
D_regulator_partials_within_wt).

### Outcome (2026-07-11, sourced clean in Positron; committed)

**A -- the correction is vindicated in the LFC data, not just NES.** Signed temporal LFC:
OXPHOS FALLS in BOTH genotypes (OXPHOS_SUBUNITS WT -0.25 / Myc+ -0.39; OXPHOS -0.13 / -0.28),
while the biosynthetic arm RISES in WT (amino-acid +0.18, lipid +0.12) and is flat/down in
Myc+. There is NO WT developmental OXPHOS rise; the aggregate "+0.82" was averaging these
opposite arms. Genotype NES is flat/high at both ages (2.2-3.5, higher at 12W) confirming the
rank ruler is magnitude-blind. **The magnitude attenuation is carried MORE by the mito/
metabolic arms than the MYC-target core** (attn ratio 12W/6W on the 6W-divergent set ~0.53
for OXPHOS/TCA/nucleotide vs 0.74-0.77 for MYC-targets; DE collapse OXPHOS 63->2, ribosome
54->3 vs MYC-targets 40->17) -- a matter of DEGREE, quantified on |LFC|, not exclusivity.

**B -- absolute nuclear-OXPHOS mRNA (the crux answer).** Myc RAISES nuclear OXPHOS in absolute
mRNA: genotype Cohen's d ~1.27 (VST/normcnt, p~0.004; logCPM d 0.84 p 0.048), ROBUST across all
3 normalizations (per-gene group-mean Spearman 0.997-0.99998). Nuclear OXPHOS DECLINES with age
in both genotypes but only DIRECTIONALLY (WT beta -0.45 p 0.21; Myc+ -0.73 p 0.11 -- ns at
n=6/tp). So **mitoPPS "OXPHOS down" is NOT an absolute-mRNA drop at the genotype level** -- Myc
elevates nuclear OXPHOS mRNA; the mitoPPS signal is within-compartment reprioritisation. This is
the defensible "mitoPPS = ratio, absolute = level" argument.

**C -- reconciliation: the mitonuclear shift is REAL in absolute mRNA.** Both lenses agree in
shape: nuclear OXPHOS peaks at 6W_pos (mitoPPS z +1.21 / absolute z +1.31) and drops by 12W;
mtDNA-encoded OXPHOS is lowest at 6W_pos and rises to 12W in both lenses (mitoPPS z -1.22->+0.37
/ absolute -1.14->+1.23). So script-24's "nuclear-led assembly at 6W -> mtDNA-completed running
metabolism at 12W" is confirmed in absolute levels. Honest divergence: at 12W_pos nuclear OXPHOS
is mitoPPS-deprioritised (0.955, <1) while absolute sits near baseline (+0.05) = reprioritisation
without a level drop.

**D -- first-pass regulators (HYPOTHESIS-GENERATING, within WT, n=12).** The WT OXPHOS trajectory
tracks BOTH the developmental luminal axis (partial r 0.71) and ER/PGC1a TF activity
(ESRRA/NRF1/GABPA; partial r 0.70-0.83), model R^2 0.86-0.89; mtDNA reprioritisation negative/ns.
Candidate MYC-independent inputs = developmental + ER/PGC1a program. The Issue #3 non-MYC OXPHOS
axis stays OPEN (not pre-answered); TF causality would need ATAC/ChIP not in these data.

**Retraction logged:** the "Issue #4 closes the Issue #3 loop via a WT developmental OXPHOS rise"
gloss is WRONG (OXPHOS falls in WT) and is retracted; it did not affect the committed Issue #3
result (partial-correlation centrality does not depend on it). Memory
[[issue4-attenuation-scoped]] and [[block-a-day1-gates]] corrected.

Paper-ready framing: *The 6W->12W attenuation is a shrinking of effect MAGNITUDE, not of program
identity: Myc engages the same top programs at both ages (flat rank enrichment) but the
wild-type<->Myc gap roughly halves, carried disproportionately by the mitochondrial/metabolic
arms (OXPHOS, TCA, nucleotide) rather than the MYC-target core. In absolute mRNA Myc raises
nuclear OXPHOS subunits (robust across normalizations); the mitoPPS "OXPHOS reduction" is a
within-compartment reprioritisation (nuclear-led at 6W, mtDNA-completed by 12W -- confirmed in
absolute levels), not a fall in nuclear OXPHOS transcript. Within wild-type tissue the OXPHOS
trajectory co-varies with the developmental and ER/PGC1a-biogenesis programs, nominating a
partly MYC-independent developmental input as a hypothesis for what re-prioritises the
mitochondrial compartment.*

---

## Issue #5 — Can MYC-independent factors EXPLAIN the OXPHOS attenuation?

**Status:** DONE (script 30, sourced clean, committed). 2026-07-11.

### The problem (why Issue #4 Part D was the wrong instrument)

Issue #4 Part D found the WT OXPHOS trajectory co-varies within WT with the developmental
luminal axis (partial 0.71) and ER/PGC1a TF activity (0.70-0.83). The author then asked the
sharper question: can we "blame" MYC-independent factors (developmental, other TFs) to FULLY
EXPLAIN the reduction? Part D tested LEVEL association (does the axis track OXPHOS?), NOT
whether the axis ABSORBS the attenuation. The attenuation IS a MODERATION phenomenon -- the
genotype (Myc) effect depends on timepoint = the genotype x time INTERACTION (b_int).
"Explaining the reduction" = does conditioning on a candidate MYC-independent axis SHRINK
b_int toward 0? Issue #5 runs that test and states the ceiling.

### Deliverable

`scripts/30_attenuation_moderation.R` (reframe on script 29's saved `$defs` set-lists; same
composite()/comp_gsva() construction so the baseline b_int reproduces absolute_stats; no
DESeq/GSVA re-run). Part A baseline b_int per outcome. Part B covariate absorption per axis
(refit `outcome ~ tp*myc + axis + tp:axis`; report delta_b_int as the STABLE primary metric,
absorption_frac = 1 - b_int_adj/b_int_base with the caveat it is unstable at a small
baseline, int_survives_p, delta_R2, nested LRT). Part C stratified case bootstrap (R=2000,
`boot`; `mediation` not installed) -> percentile CI on delta_b_int and absorption_frac. Part
D within-WT (n=12) OXPHOS slope before/after holding dev + TF (the script-17 analogue). Part
E "what would settle it" (inducible Myc, TF perturbation, single-cell, ATAC/ChIP). Outcomes:
`oxphos_abs` (VST nuclear-OXPHOS, powered/mtDNA-clean) + `oxphos_gsva` (secondary). Axes:
dev_luminal + tf_biogenesis (primary), prolif + mtdna_reprior (secondary), dev+tf combined.
Artifacts: `results/attenuation_moderation.rds`, `outputs/attenuation_moderation/*`
(B_bint_before_after, C_absorption_fraction_bootCI, C_delta_bint_bootCI,
D_withinwt_slope_absorption).

### Outcome (2026-07-11, sourced clean in Positron; committed)

**The answer is NO: no MYC-independent axis absorbs the attenuation.** Baseline b_int on the
powered outcome oxphos_abs = -0.279 (p 0.607, directional -- reproduces Issue #4 absolute_stats
exactly); on oxphos_gsva ~0 (0.003), so its absorption fraction is undefined and read via
delta only.

- **Part B.** Conditioning on the developmental + TF axes together leaves the attenuation
  essentially UNTOUCHED (b_int -0.279 -> -0.322, delta -0.043). dev_luminal ALONE makes it
  LARGER (-0.279 -> -0.721, delta -0.442) -- the same co-drive signature as script 17
  (removing a program Myc co-drives EXPOSES more Myc effect, it does not absorb it).
  tf_biogenesis / prolif / mtdna flip or shrink the point estimate but none toward a clean 0.
- **Part C (the honest ceiling).** EVERY delta_b_int bootstrap CI straddles 0 (e.g. dev+tf
  [-1.39, 1.31]; dev_luminal [-1.17, 0.40]) -- so at n=24 there is no RELIABLE movement of the
  attenuation in either direction: no evidence of absorption, and even the "enlarges it"
  point estimate is not significant. The absorption FRACTION CIs are uninterpretably wide
  (spanning +-10 to +-18) because the baseline interaction is directional -- "fully explains"
  (frac~1) is nowhere distinguishable.
- **Level-association is NOT moderation.** Every axis LRT is significant (lrt_p 5e-4 to
  1e-10): the axes DO add explanatory variance for OXPHOS LEVEL -- they just do not absorb the
  genotype x time INTERACTION. This is exactly why Part D of Issue #4 (level association) was
  the wrong instrument for the "explains the attenuation" question.
- **Part D.** Within WT (n=12) the OXPHOS 6W->12W slope only shrinks ~29% (oxphos_abs
  -0.454 -> -0.323) / ~19% (gsva) when dev + TF are held -- it PERSISTS, the OXPHOS analogue of
  script 17's biogenesis non-absorption.

**CEILING (critical).** Candidate axes are ENDOGENOUS (Myc drives dev/TF too) -> absorption
estimates are BIASED; n=6/group; collinear axes; the baseline interaction is itself ns. This
test BOUNDS the claim -- it cannot prove a MYC-independent cause, and a large absorption
fraction must NOT be read as "development explains it" (report the CI + endogeneity caveat
together). The identifying evidence is experimental (Part E: inducible Myc off/on, TF
perturbation, single-cell/deconvolution, ATAC/ChIP).

Paper-ready framing: *We tested directly whether conditioning on MYC-independent developmental
or ER/PGC1a-biogenesis programs absorbs the genotype x time interaction that constitutes the
attenuation. It does not: holding both axes leaves the interaction essentially unchanged, and
conditioning on the developmental axis alone enlarges rather than shrinks it -- the programs
co-vary with mitochondrial output rather than antagonising Myc (consistent with the earlier
finding that development co-drives, not absorbs, the biogenesis trajectory). These axes remain
associated with OXPHOS level but do not account for the attenuation, and at this sample size
the analysis bounds rather than resolves a MYC-independent contribution; distinguishing a
fading transgene effect from a fixed developmental ceiling requires inducible and perturbation
experiments.* This CLOSES the "can we blame MYC-independent factors" question at the level
these bulk data can answer.

---

## Issue #6 — A citable MECHANISM for the OXPHOS attenuation

**Status:** DONE (script 31, sourced clean, committed). 2026-07-12.

### The problem (compression re-describes; we need a mechanism)

Issues #4/#5 fixed WHAT the attenuation is (magnitude compression, not rank; Myc raises
absolute nuclear OXPHOS) and that no dev/TF axis ABSORBS the interaction. But "compression"
only renames the attenuation -- the author wanted a CITABLE mechanism. The restatement ("the
MYC signal is not reduced, it acts on a different background") splits into three non-exclusive
mechanisms: moving background (WT converges), Myc-fade (front-loaded), and global response
compression (breadth reduction). A read-only check showed the breadth reduction is real but
modest (SD ratio 0.83) and concentrated in the strong responders (the DE-count collapse -91%
is a tail/threshold artifact overstating a -20% magnitude compression) -- so compression alone
is not the citable mechanism.

### The resolution (an exact identity separates two mechanisms)

The attenuation per gene IS the interaction: `myc_12W - myc_6W == timepoint_pos -
timepoint_neg` (verified, cor = 1.0000). Aligned to Myc's induction direction d = sign(myc_6W),
the gap-shrink splits cleanly into **WT-convergence** (`d*timepoint_neg`, WT matures toward
Myc -- the moving background) and **Myc-fade** (`d*timepoint_pos`, oncogenic retreat).

### Deliverable

`scripts/31_attenuation_mechanism.R` (reframe on fitted contrasts + GSVA; no re-run). Part A =
convergence/fade decomposition (POWERED, contrast-level, all genes) over two universes
(6W-divergent padj<0.1; and effect |LFC6|>0.5, selection-independent) x roster programs + pooled
proliferation/mammary-luminal. Part B = broadened-TF absorption (biogenesis/E2F/ESR1/lipogenic/
data-driven mito panel) extending Issue #5, plus a compositional diagnostic (each TF axis's
Myc+ 6W->12W slope vs the OXPHOS fade). Part C (deconvolution) DEFERRED. Artifacts:
`results/attenuation_mechanism.rds`, `outputs/attenuation_mechanism/*` (A1 convergence_fade_by_
program, A2 conv_fade_gene_map, A3 conv_fade_shares, B1 tf_bint_before_after, B2 tf_delta_bint_
bootCI, B3 tf_temporal_vs_oxphos_fade).

### Outcome (2026-07-12, sourced clean in Positron; committed)

**A -- the mechanism (powered, robust across both universes).** Globally the attenuation is
**~66% Myc-fade, ~34% WT-convergence** (42.8% convergence on the effect universe). It is
**PROGRAM-SPECIFIC:** WT-convergence is real for the biosynthetic arm (amino-acid conv 37%,
lipid 41%, nucleotide 31%; TCA 20%) and for pooled mammary-luminal (36%) -- the wild-type gland
matures onto the same biosynthetic axis Myc drives. But it is ABSENT/NEGATIVE for OXPHOS (conv
-11 to -39%; WT diverges, the gap closes purely by fade, fade% >100) and for the MYC-target
core (conv -60 to -63%), which is ALSO the LEAST attenuated (atten 0.17-0.19 vs 0.28-0.50) --
i.e. the MYC identity core is PROTECTED. Pooled proliferation attenuates by PURE fade (conv
~0%, fade 100%). Every sign reproduces on the selection-independent effect universe -> not a
significance-selection artifact. This is the powered, gene-level resolution of the old H1
(front-loaded fade) vs H3 (WT catch-up) question (Gate 1 / script 17 left it directional-only):
H1-dominant (66% fade) with a real, biosynthetic-arm-specific H3 minority.

**B -- broadened-TF absorption (BOUNDED) + a compositional hint.** No broadened TF axis reliably
absorbs the attenuation -- every delta_b_int bootstrap CI straddles 0 (Issue #5 ceiling holds);
ESR1 does not even associate (lrt_p 0.83). The numerical standout is tf_e2f (proliferation/cell-
cycle): it moves b_int on oxphos_abs from -0.279 to -0.028 (point-estimate ~90% absorption) --
but NOT significant (CI [-0.64, 2.17]), a hint not a result. Compositional diagnostic: the
OXPHOS fade is Myc+ slope -0.91; the axes co-declining most are tf_mito_panel (-0.60), lipogenic
(-0.46), E2F (-0.36), while ESR1 RISES (+0.42). So proliferation/mito-TF activity fades
alongside OXPHOS -- SUPPORTING (not proving) that the dominant Myc-fade is partly compositional
dilution of the shrinking proliferative/TEB compartment.

**C -- DEFERRED (settle_it).** Reference deconvolution of the proliferative/TEB FRACTION (needs
an external mouse-mammary sc atlas e.g. Bach 2017 GSE106273 + MuSiC/Bisque; not on disk) is the
next step; single-cell/snRNA-seq or TF perturbation is the definitive test.

**CEILING.** Part A is POWERED and the identity is exact -> the convergence/fade split is solid
but DESCRIPTIVE of the transcriptional change (composition still confounds the fade term; not a
per-cell causal claim). Part B is BOUNDED exactly as Issue #5 (endogenous TFs, n=6/group,
directional baseline). The E2F/proliferation co-decline is a hypothesis-generating pointer to
composition, to be settled by deconvolution/single-cell.

Paper-ready framing: *The genotype gap does not close uniformly. Decomposing the exact
genotype x time interaction into a wild-type-convergence term and a Myc-fade term shows the
attenuation is ~two-thirds attenuation of the oncogenic program on the aging substrate and
~one-third wild-type developmental convergence -- the latter confined to the biosynthetic/
metabolic arm (amino-acid, lipid, nucleotide), where the maturing wild-type gland ramps onto
the same program Myc drives. OXPHOS and the core MYC-target program show no such convergence
(wild-type diverges), and the MYC identity core is the least attenuated of all -- selectively
buffered. The dominant Myc-fade co-occurs with a collapse of proliferative and mitochondrial-TF
activity, consistent with (but not proof of) dilution of the shrinking proliferative
compartment; resolving per-cell fade from compositional dilution requires single-cell or
deconvolution data.* This gives a CITABLE two-mechanism account of the attenuation and scopes
the one remaining ambiguity to a defined next experiment.
