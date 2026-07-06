---
title: Block A build spec — myc_mouse finalisation
date: 2026-07-04
tags: [project/myc_mouse, planning, block-a, build-spec]
status: draft for author review
relates_to:
  - "docs/myc_mouse_finalisation_plan.md (the authoritative plan)"
  - "docs/branch_manifest.md (script inventory)"
  - "CLAUDE.md (workflow, conventions, data traps)"
purpose: >
  Turns the finalisation plan's Day 1-5 sequence and action-point menu into
  concrete, ordered, Option-A work: for each step, which script is written or
  modified, what it reads, what it writes, and which action point / gate it
  serves. Gates-first dependency logic.
---

# Block A build spec — myc_mouse finalisation

## Context

We are finalising the myc_mouse arm of the manuscript "Mitochondria integrate
oncogenic and metabolic transcriptional programs to shape breast cancer
progression." The authoritative plan (`docs/myc_mouse_finalisation_plan.md`)
reframes the science (two-trajectory, permission-vs-selection, stable Myc dose)
and lays out a Day 1-5 draft build across an action-point menu (AP0-AP8 +
AP-CD/mtDNA/abund/retention). This spec sequences that menu into numbered
scripts under the gates-first dependency logic, on the `analysis-exploratory`
(Block A) branch. Block A is broad-exploratory; Block B (later) re-codes only
the narrative-selected subset for submission.

**Option A throughout:** Claude Code *writes/edits* the numbered scripts and
runs infra (git, snapshots, docs). The author *sources them in Positron*
interactively. Every numbered script ends with an `if (FALSE) { ... }` sandbox
block. Claude never auto-runs analysis scripts.

## Decisions locked (2026-07-04 planning session)

1. **fGSEA scope** — per-category (one run per `by_category/*.gmt`), **BH within
   each run**. The existing pooled runs in script 10 stay, labelled
   exploratory-omnibus only; publication claims cite per-category padj.
2. **AP6 comparators** — **both** panels: generic MSigDB Hallmark (broad) +
   mammary-specific `04_metabolism` / `05_proliferation` library GMTs (sharper,
   tissue-matched).
3. **Hallmark source** — pulled **fresh via `msigdbr(species="Mus musculus",
   category="H")`** inside the new fGSEA script; self-contained, decoupled from
   the legacy `gene_sets_list.rds`. (Not in the library — Hallmark is generic
   MSigDB "H".)
4. **01 gene-set tidy** — **deferred to Block B** (or dropped entirely; decide
   later). Rerouting `01`'s legacy `mitocarta_pathways.csv` / `.gmx` loads would
   rename `MC_*`/`MYC_*`/`MSigDB_*` sets that scripts 10/11 hard-reference. Block
   A leaves `01` untouched; the new per-category fGSEA consumes the snapshot
   GMTs directly.

## Code-read check results (Day 1, no execution)

**Check 1 — Gate 2 inputs (tension T3): RESOLVED. Raw / Wald-based; no recompute.**
- `03_deseq_results_qc.R:77-79` builds `interaction_raw` =
  `results(dds_int, name="timepoint12W.myc_statuspos", filterFun=ihw)` — the
  **unshrunken MLE fit**. `interaction_shrunk` (ashr) is a separate object
  (`03:112-115`). IHW reweights only `padj`, never `log2FoldChange`/`stat`.
- `11_interaction_gene_characterisation.R:126` builds `interaction_raw_df` from
  `interaction_results$interaction_raw`. **Part 2 direction** (binomial, median
  LFC) and **Part 5B PRO/ANTI** (`get_interaction_stats(..., interaction_raw_df)`
  at `11:908-911`; one-sample t-tests `11:1102-1108`) both flow from that raw
  object. The g:Profiler ranking uses `stat` (Wald) from the same fit.
- **Conclusion:** Gate 2 is trustworthy read straight from
  `interaction_gene_characterisation.rds`. The plan's T3 caution is moot for
  this data — recorded so it is not re-litigated.

**Check 2 — raw cell-death input + raw heatmaps: confirmed, one wrinkle.**
- `results/combined_df_annotated_raw.rds` **exists** (1.07 MB, Feb 7), beside
  shrunken `combined_df_annotated.rds`. Archived branch-1
  (`archive_main_pipeline/09_cell_death_pathway_summary.R:30`) reads the
  **shrunken** file — the port's substantive change is repointing that read to
  the `_raw` file. Metric `myc_6W_log2FC - myc_12W_log2FC` (`:71`), keys
  `mgi_symbol`/`effect`.
- `outputs/heatmaps_int/*_raw.pdf` + `*_shrunk.pdf` confirmed present per
  pathway (e.g. `MC_Apoptosis_Pro_raw.pdf`). Use `_raw` for figures.
- **Wrinkle (flag):** `combined_df_annotated_raw.rds` and the `_raw` heatmaps
  are products of the **archived main pipeline** (`archive_main_pipeline/02,03,
  04,07`), not the new-analysis trunk (00-12). For Block A they are stable dated
  inputs — reading them is fine and matches the plan. But the trunk is **not
  self-contained** for these until the builder is ported (Block B task). The
  branch-1 port must verify columns at source in Positron (`head()` the raw df).

**Contradiction with plan assumptions:** the plan's Day-2 phrasing "port main 09
forward onto the interaction model" is really "port branch 1 to read an existing
archived-pipeline raw artifact," not "recompute from the trunk." No hard data
contradictions.

## Script inventory (Block A additions, numbered 13+)

New scripts land at 13+ so they never collide with the 00-12 trunk. Reframe
scripts read existing `.rds`; build scripts add new computation.

| # | Script | Type | Serves | Day |
|---|--------|------|--------|-----|
| 13 | `13_gate1_divergence_timing.R` | reframe | AP5 / Gate 1 | 1 |
| 14 | `14_gate2_apoptosis_readout.R` | reframe | Gate 2 | 1 |
| 15 | `15_gsva_scoring.R` | build | AP1/2/3/4/7 infra | 2 |
| 16 | `16_cell_death_binomial_raw.R` | port | AP-CD branch 1 | 2 |
| 17 | `17_gsva_overview.R` | build | GSVA read / H1-vs-H3 discrimination | 3 |
| 18 | `18_ap7_mb_fork_projection.R` | build | AP7 (centrepiece) | 3 |
| 19 | `19_dev_composition.R` | build | AP1, AP2, AP3, AP4 | 3 |
| 20 | `20_fgsea_percategory.R` | build | AP6.1, AP-retention, AP-abund | 3 |
| 21 | `21_ap6_permutation_null.R` | build | AP6.2 | 4 |
| 22 | `22_reframe_mtdna_abund_retention.R` | reframe | AP-mtDNA, AP-abund, AP-retention | 4-5 |
| 23 | `23_death_timing_substrate.R` | build | deferred death-timing (H1-H4) -> Supp | 5 |
| 24 | `24_figure1.R` | build | Fig 1 draft | 5 |
| 25 | `25_figure2.R` | build | Fig 2 draft | 5 |
| 26 | `26_supplementaries.R` | build | Supp 1-3 draft | 5 |

**Renumber note (2026-07-06):** the GSVA-overview / trajectory-visualisation
script was inserted as the new **17** (author-requested read of `gsva_scores.rds`
before the AP scripts consume it — the H1-vs-H3 discriminator). Everything that
was 17-24 shifted **+1** to 18-25. No script files moved (none existed yet). The
Day 3/4/5 sections below and all cross-references use the new numbering.

**Renumber note (2026-07-06b):** the deferred death-timing analysis was pulled
forward (author decision) to run *before* the figures, because its result shapes
the final interpretation/narrative. It is inserted as the new **23**
(`23_death_timing_substrate.R`, all four hypotheses H1-H4, anchored on the WT
`timepoint_neg` substrate — see memory `myc-death-timing-question`). The figure
scripts shifted **+1**: `23_figure1`->24, `24_figure2`->25,
`25_supplementaries`->26. Script **22** is unchanged. No script files moved.

`12_cell_death_pathway_analysis.R` (branch 2, Tang sets on Wald z) is **done**
(`cell_death_fgsea.rds`); run as-is, no edit. `08` mitoPPS is done
(`mitopps_scores.rds`). `AP8` (cross-sample CV) is optional/last, only if the
selection arm survives Gate 2 — added as a tail block in `24` if reached.

## Ordered build — Day by day

### Day 1 — Decision layer + gates (existing data, lowest risk)

**Infra (Claude runs, after asking):**
- Confirm on `analysis-exploratory` (already are). Library snapshot already in
  `data/genesets_from_library/` incl. `by_category/` + `provenance_table.csv` +
  README — **already present**, no snapshot step needed. Verify pin note in
  README references `v1.0`.

**13 — `13_gate1_divergence_timing.R` (AP5 / Gate 1).**
- Reads: `interaction_results.rds` (`myc_6W_raw` for the 6W genotype contrast;
  `interaction_raw` for direction), `interaction_gene_characterisation.rds`
  (script-11 Part 2 direction summary + the p<0.05 directional set).
- Computes/collates: (i) genotype difference at 6W alone (how many genes
  already divergent at earliest sample); (ii) interaction direction across
  6W->12W (grows / flat / shrinks — negative = Myc effect weakens at 12W); (iii)
  a **placeholder hook** for the stage-annotated dev-set projection, which
  depends on GSVA (AP2) and is filled Day 3 from `15`/`18` outputs.
- Writes: `results/gate1_divergence_timing.rds`, a short
  `outputs/gates/gate1_summary.csv`, one direction/timing plot.
- Gate logic: if divergence *grows* across the window, "developmental change
  licenses Myc" weakens -> title moves off "licenses"; if already-large-and-flat
  at 6W, the soft "already diverged" claim holds. Feeds the figure lock.

**14 — `14_gate2_apoptosis_readout.R` (Gate 2).**
- Reads: `interaction_gene_characterisation.rds` (Part 5B objects —
  Apoptosis-PRO 25 genes, Apoptosis-ANTI 9 genes, the one-sample t-tests; on
  **raw** interaction LFC per Check 1).
- Computes/collates: confirm both t-tests, PRO/ANTI group shift, Bbc3 outlier
  status; interpret selection (module-wide coordinated shift) vs readout
  (isolated/absent — current padj ~0.91 tilts readout/metabolic).
- Writes: `results/gate2_apoptosis_readout.rds`,
  `outputs/gates/gate2_summary.csv`.
- Gate logic: decides whether the mitochondrion is a decision point (principle
  5) or a readout; if null, AP4 + apoptosis framing shrink, lean on the
  metabolic-readout interpretation; cell death held at Supp 1.

**Dated decision note (Claude, on approval -> a `docs/2026-07-0X_gates_note.md`):**
what the gates returned, which framing/APs survive, figure lock, title-verb call.

*End Day 1: permission-vs-divergence fork resolved; selection arm confirmed or
demoted; figure plan locked.*

### Day 2 — Cell death (both branches, raw) + GSVA engine

**16 — `16_cell_death_binomial_raw.R` (AP-CD branch 1, ported).**
- Source: `archive_main_pipeline/09_cell_death_pathway_summary.R`, ported
  forward. **Single substantive change:** read
  `results/combined_df_annotated_raw.rds` (not the shrunken file). Keep metric
  `myc_6W_log2FC - myc_12W_log2FC`, binomial directional test on
  `cell_death_genes_consolidated`.
- Reads: `combined_df_annotated_raw.rds`,
  `data/cell_death_genes_consolidated.csv`.
- Writes (regenerated on raw, suffixed to avoid clobber): e.g.
  `results/cell_death_hypothesis_results_raw.csv`,
  `results/cell_death_genes_full_raw.csv`, figures under
  `outputs/cell_death_raw/`.
- **Positron verification first line:** `head()` the raw df, confirm it carries
  `mgi_symbol`, `myc_6W_log2FC`, `myc_12W_log2FC`, `effect` before trusting the
  join.
- Then run **12 as-is** (branch 2) and assemble a convergence read (where the
  two branches agree/diverge) into `results/cell_death_convergence.rds`.

**15 — `15_gsva_scoring.R` (GSVA infrastructure — the core new build).**
- Input: **log-scale** VST from `dds_int_run.rds` (`DESeq2::vst`); `kcdf=
  "Gaussian"`; **all samples in one run** (cohort-relative). This is the
  *opposite* of mitoPPS (linear) — guard the trap in a header comment.
- Sets: route by method tag from `provenance_table.csv` — GSVA-tagged sets only.
  Primary consumers: `03_mammary_development` (`MG_*` dev/state sets, AP1/2/3),
  `02_myc_signatures` (Felsher, AP4), `07_biogenesis_discrimination` (Fig 2),
  and the MB-fork sets for AP7 (see 17). Ortholog-map via the cached
  `results/ortholog_table.rds` + `dplyr::inner_join` (not `deframe()`), except
  the documented exemptions.
- Writes: `results/gsva_scores.rds` (set x sample matrix + set metadata).
- Validate in Positron (score distributions, no all-NA sets, sample order).

*End Day 2: two-branch cell-death convergence on raw; validated GSVA matrix.*

### Day 3 — GSVA overview read, then GSVA-dependent results + per-category fGSEA

**17 — `17_gsva_overview.R` (GSVA read; H1-vs-H3 discriminator) [runs first].**
- Reads: `results/gsva_scores.rds` (902 GSVA-tagged sets x 24 samples).
- Computes: per set, `lm(score ~ timepoint * myc_status)` on the 24 samples ->
  coefficient table (`beta_time` = WT/Myc- slope; `beta_myc` = d6; `beta_int` =
  d12-d6 = the Gate-1 attenuation; derived `d12`, Myc+ slope), interaction p
  BH-adjusted **within `category_primary`**; 4-column condition-mean matrix
  (`6W_neg,6W_pos,12W_neg,12W_pos`).
- Writes: `results/gsva_overview.rds`; `outputs/gsva_overview/` -
  `<category>_landscape.pdf` (top ~35 sets by interaction p, row-z-scored, WT-slope
  + interaction row annotations, ComplexHeatmap) per category, plus
  `<category>_profiles.pdf` (4-point two-line profile) and `<category>_dumbbell.pdf`
  (d6->d12 gap collapse) for the four trajectory categories
  (03_mammary_development, 02_myc_signatures, 07_biogenesis_discrimination,
  08_apoptosis).
- Rationale: interaction = the attenuation metric (rank on it); the two-line
  profile shows *which* genotype trajectory moved -> discriminates H1 (Myc+
  descends) from H3 (Myc- ascends). Feeds the Gate-1 projection hook (13) and
  sanity-frames dev-composition (19).

**18 — `18_ap7_mb_fork_projection.R` (AP7, CENTREPIECE).**
- Reads: `results/gsva_scores.rds` (or projects fresh if fork sets excluded from
  15), `data/genesets_from_library/metabric_sets.rds` and/or the mouse fork GMTs.
- **Build-time check (flag):** confirm the MB1/MB2/MB12 fork signatures'
  **species** — `metabric_sets.rds` is likely human (METABRIC); AP7 needs
  mouse-mapped fork sets. Resolve at first line: either the library exported
  mouse fork GMTs (check `by_category/07_*` and provenance), or ortholog-map the
  human MB sets. Record which.
- Computes: GSVA-project mouse tumours onto the fork signatures; test whether
  Myc+ progression shifts samples toward **MB2_UF** specifically (falsifiable).
- Writes: `results/ap7_mb_fork.rds`, projection plot(s).

**19 — `19_dev_composition.R` (AP1, AP2, AP3, AP4).**
- Reads: `gsva_scores.rds`, `gsva_overview.rds`, `coldata.rds`, `interaction_results.rds`.
- AP1: WT developmental shift — `MG_*` GSVA + mitoPPS on Myc- 6W vs 12W (mitoPPS
  exists in `mitopps_scores.rds`).
- AP2: lineage-composition genotype effect — `MG_*` GSVA on all samples,
  Myc+ vs Myc- controlling timepoint (the central permission test).
- AP3: interaction on dev-set GSVA + vector correlation (Myc-vs-WT difference
  vector vs WT developmental-shift vector). Corroboration only.
- AP4: Felsher membership lookup + one GSVA co-variation corr (dose constant).
- Writes: `results/dev_composition.rds`; feeds the Gate 1 projection hook (13).

**20 — `20_fgsea_percategory.R` (AP6 part 1 + AP-retention/AP-abund inputs).**
- Reads: `interaction_results.rds` (rankings from `stat` on the **unshrunken**
  `*_raw` contrasts: `interaction_raw`, `timepoint_pos_raw`, `timepoint_neg_raw`,
  `myc_6W_raw`, `myc_12W_raw`); `by_category/*.gmt`; fresh Hallmark via
  `msigdbr`; ortholog table for symbol mapping.
- Engine: loop rankings x categories; **one fGSEA call per (ranking, category
  GMT)**, `minSize`/`maxSize` per library method tags (GSVA large-set exemption
  N/A here — this is fGSEA), **BH within each run**. Directional NES retained.
- Category -> AP routing (plan section 5): `01_mitocarta`->AP6+AP-mtDNA;
  `04_metabolism`/`05_proliferation` + fresh Hallmark->AP6 comparators;
  `06_tf_targets`->AP6 TF/dev lanes; `02_myc_signatures`->AP-retention;
  `03_mammary_development`->dev enrichment; `08_apoptosis`->Gate 2/AP-CD context;
  `07`,`09`->Fig 2 discrimination/intersection.
- Writes: `results/fgsea_percategory.rds` (tidy: ranking, category, pathway,
  NES, pval, padj_within_category, leadingEdge).

*End Day 3: validation (AP7), composition (AP1/2), decomposed enrichment (AP6.1).*

### Day 4 — Finish justification + Figure 1

**21 — `21_ap6_permutation_null.R` (AP6 part 2).**
- Reads: `interaction_results.rds` (mean |LFC| per gene), MitoCarta membership
  (from `mitopps_scores.rds` gene_to_pathway or `01_mitocarta` GMT), comparator
  sets.
- Computes: expression x dispersion-matched permutation null on mean |LFC| — bin
  all genes by baseMean x dispersion, sample matched random sets, locate
  MitoCarta (and each comparator) in the null. Direction-agnostic magnitude
  test; defeats the housekeeping confound. Optional ROAST-mixed/CAMERA is a
  cut-line (section 9 item 3) — include a guarded stub only.
- Writes: `results/ap6_permutation_null.rds`.

**24 — `24_figure1.R` (draft).** MB-fork projection (AP7) + preferential-
alteration panel (AP6) + WT developmental mito shift (AP1). Reads the Day 2-4
`.rds`. Writes `outputs/figures/figure1_draft.pdf`.

*End Day 4: focus justification complete; Fig 1 in draft.*

### Day 5 — Figure 2 + supplementaries + draft text

**22 — `22_reframe_mtdna_abund_retention.R` (reframe, feeds Supp 3).**
- AP-mtDNA: off `mitopps_scores.rds` synthetic mtDNA pathway — mtDNA-encoded vs
  nuclear OXPHOS split, per-complex correlations (per-complex detail is cut-line
  5).
- AP-abund: fGSEA (rank-relative) vs mitoPPS (within-compartment) contrast;
  **time-dependent** answer (6W selective reprioritisation -> 12W uniform
  biogenesis). Do not collapse to one word.
- AP-retention: retention index NES_pos-NES_neg (absolute, the meaningful one)
  and NES_pos/NES_neg on temporal contrasts; carry the Felsher guardrail
  ("developmental decline + stable Myc offset", not loss of Myc activity).
- Reads: `fgsea_percategory.rds`, `mitopps_scores.rds`, `fgsea_results.rds`,
  `fgsea_xs_results.rds`. Writes `results/reframe_supp3.rds`.

**23 — `23_death_timing_substrate.R` (build; deferred death-timing, pulled
forward before the figures).** One integrated 6W-death SUBSTRATE model spanning
cell-death branch 1 (16), branch 2 (12), Gate 2 (14), and the script-17
apoptosis/intersection trajectories. Anchored on the WT `timepoint_neg` substrate
(the tissue the acute inducible-Myc pulse hits — Myc off pre-tamoxifen); the
chronic-Myc+ layer (`timepoint_pos`, `myc_6W`) is secondary and survivor-biased.
Covers all four hypotheses:
- **H1** BH3-only : anti-apoptotic BCL2 rheostat + p53/ARF readiness on the WT
  substrate (priming index falls 6W->12W = less-primed older tissue). Lead.
- **H2** biogenesis-death decoupling over time (ER/PGC1a-couples vs Myc-uncouples).
- **H3** proliferation-apoptosis coupling (oncogene-induced apoptosis at 6W).
- **H4** selection / survivor culling via cross-sample CV narrowing (weakest arm).
- Reads: `interaction_results.rds`, `mitopps_scores.rds`,
  `gate2_apoptosis_readout.rds`, `cell_death_fgsea.rds`,
  `cell_death_*_raw.csv`, `gsva_scores.rds`, `gsva_overview.rds`,
  `ortholog_table.rds`, fresh Hallmark P53_PATHWAY (msigdbr). Writes
  `results/death_timing_substrate.rds`, `outputs/death_timing/`. Honest ceiling:
  bulk RNA + survivor bias + n=6 -> substrate association, not causation. Feeds a
  Supp panel and the discussion; see memory `myc-death-timing-question`.

**25 — `25_figure2.R` (draft).** Lineage-composition GSVA (AP2) + Category 7
discrimination (MYC_SPECIFIC/CORE/DEVELOPMENTAL from `07_*` GMT) + divergence
timing (AP5). Writes `outputs/figures/figure2_draft.pdf`.

**26 — `26_supplementaries.R` (draft).** Supp 1 (apoptosis PRO/ANTI + two-branch
cell-death convergence + death-timing substrate model, from 14/16/12/23) - Supp 2
(permutation null + comparator panel + QC/provenance, from 21/20) - Supp 3 (mtDNA
split + abund panel, from 22). Optional **AP8** tail block (cross-sample CV of
mito-fork score) only if the selection arm survived Gate 2 / H4 (script 23).
Writes `outputs/figures/supp{1,2,3}_draft.pdf`.

**Draft text (author):** AP0 intro framing, results prose with gaps flagged,
figure legends.

*End Day 5: Figs 1-2 + 3 supplementaries in draft; section text drafted.*

## Cut-lines (plan section 9, drop first if behind)

1. AP3/AP4 fine detail - 2. AP8 (script 25 tail) - 3. ROAST/CAMERA (stub in 21)
- 4. comparator breadth -> proliferation only - 5. Supp 3 per-complex detail -
6. figure polish. **Not a cut-line:** the branch-1 raw re-run (16) — it feeds
Gate 2 and runs regardless.

## Verification (how each piece is checked end-to-end, in Positron)

- **Gates (13/14):** run script; inspect `gate{1,2}_summary.csv`; confirm
  direction sign convention (negative interaction = Myc effect weakens at 12W)
  matches script 11 comments. No recompute needed (Check 1).
- **Branch-1 port (16):** first line `head(combined_df_annotated_raw)` — confirm
  `mgi_symbol`/`myc_6W_log2FC`/`myc_12W_log2FC`/`effect` present; then confirm
  the binomial output differs from the shrunken run only as expected.
- **GSVA (15):** score-distribution sanity (Gaussian, no all-NA sets); confirm
  log-scale input (VST), all-samples-one-run; spot-check a known dev set.
- **GSVA overview (17):** confirm `beta_myc == d6` (lm vs group-means agree); the
  four trajectory categories' two-line profiles separate H1 (Myc+ down) from H3
  (Myc- up); landscapes show <=~35 rows after the interaction-p cut.
- **AP7 (18):** confirm fork-set species resolved to mouse; check MB2_UF shift
  direction is testable/falsifiable.
- **Per-category fGSEA (20):** confirm BH is within-run (not pooled); spot-check
  one category's padj against a manual `p.adjust`.
- **Permutation null (21):** confirm bins are expression x dispersion; null is
  direction-agnostic (|LFC|).
- Each script's `if (FALSE)` sandbox holds the line-by-line inspection calls.

## Flags recorded (no Block A action unless noted)

- **T3 resolved** (Check 1) — record in a project memory / gates note so it is
  not re-litigated: Gate 2 direction + PRO/ANTI already run on raw/Wald values.
- **Trunk self-containment** — `combined_df_annotated_raw.rds` and
  `outputs/heatmaps_int/*_raw.pdf` are archived-pipeline outputs; porting their
  builder into the trunk is a Block B task.
- **01 tidy** — deferred (or dropped); `gene_sets_list.rds` stays alive for
  10/11 in Block A.
- **MB-fork species** — resolve mouse mapping at 18's first line.
- Housekeeping (plan section 10): `NES_paradox_explanation.md` "liver" mislabel;
  mitoPPS outlier call; stray `output/` (singular) dir; stale `branch_manifest.md`
  — tidy when convenient, not gating.

## First execution steps when the author approves the build

1. Record the T3-resolution + trunk-self-containment flags (gates note / memory).
2. Write scripts **13** and **14**; author sources them in Positron; Day-1
   gates note. Then proceed Day 2+.
