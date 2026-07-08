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

Three independent lenses converge. WT gland matures toward BMYO / loses luminal;
Myc bends OFF that axis: convergence rho(Myc,WT) = -0.08 at 6W (orthogonal) ->
-0.44 at 12W (oppositional). Per state: BMYO -0.27 -> -0.71 (Myc suppresses basal,
deepening), LASP +0.35 -> +0.21 (spared/parallel), LHS -0.14 -> -0.44. Profile
Myc-effect: BMYO -0.15/-0.21, LASP -0.05/+0.09, LHS -0.15/+0.17 = the **LHS flip**
(6W suppress -> 12W elevate, against the WT LHS decline). Corroborated by the GRAY
directional nets (Myc -> low-estrogen + TEB/proliferative, strongest in BMYO) and
the CHUNG ATAC OPEN-CLOSED nets (Myc closes basal chromatin, opens LP early, opens
ML/LHS by 12W). Reading: at 6W a broad off-axis identity suppression; by 12W a
specific **anti-BMYO / pro-LHS reprogramming**. Ceiling: per-set GSVA contrasts
powered (24 samples); n=6 interaction directional; association not causation.
Output = the developmental figure substrate for Block B.

Artifacts: `scripts/26_dev_program_myc_integration.R`,
`data/dev_mec_annotation.csv` (curated MEC annotation),
`results/dev_program_myc_integration.rds`, `outputs/dev_program_myc_integration/*`.
