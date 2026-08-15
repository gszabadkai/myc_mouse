# CLAUDE.md — MMTV-Myc Mouse Timecourse Analysis

Rules and context for Claude Code sessions in this repo. Read this first, then the
plan: `docs/myc_mouse_finalisation_plan.md`.

## What this repo is

Bulk RNA-seq of **purified mammary epithelial cells (MECs)** from MMTV-Myc transgenic
mice — **not whole tissue**. Purification is by **enzymatic dissociation only (no FACS)**.
This matters and was undocumented until 2026-07-16: any stromal/endothelial/adipocyte
signal is **contamination**, not tissue composition (~3–12% residual, consistent with no
sort), and the warm enzymatic digest is itself the stressor that induces the immediate-early
response (van den Brink 2017) — see the mt-axis note below. No dissociation-batch,
viability or RIN metadata exists on disk. **BATCH = TIMEPOINT** (undocumented until
2026-07-19): the 6W and 12W cohorts were extracted as two separate batches, so batch is
perfectly confounded with timepoint. Genotype is balanced within each batch, so the
`~timepoint*myc_status` design absorbs the whole batch effect into the timepoint main
effect — **genotype/interaction contrasts stay clean; the pure between-timepoint
(developmental) reading is confounded with batch and not separable** (no QC metadata).
Four groups x six replicates
(`6W_neg`, `6W_pos`, `12W_neg`, `12W_pos`), DESeq2 genotype x timepoint interaction
model. Investigates Myc-driven transcriptional and mitochondrial changes across the
two timepoints. Feeds the manuscript "Mitochondria integrate oncogenic and metabolic
transcriptional programs to shape breast cancer progression". Analytical companion
paper: Menegollo, Bentham et al., Cancer Res 2024 (CAN-23-3172).

## Current phase

Two-pipeline consolidation is done: `new-analysis` is the trunk (scripts 00-12);
main's removed scripts are archived byte-identical in `scripts/archive_main_pipeline/`.
The finalisation (`docs/myc_mouse_finalisation_plan.md`) runs in two phases:
- **Block A** — broad exploratory (scripts 13-25) + the step-by-step revision (Issues
  #1-6, scripts 26-31). **DONE and reviewed** (2026-07-12): synthesised in
  `docs/2026-07-12_BlockA_revision_synthesis_and_figure_plan.md`; the reviewed commit is
  tagged `block-a-reviewed` (4b82a82) and now lives on `analysis-exploratory` (00-31).
- **Block B** — publication figures, on `paper-figures`. **FROZEN 2026-08-14** at tag
  `block-b-full`: 25 panels, three circulation PDFs, and `paper/analysis_record.qmd`, the
  record of the whole analysis before the cut.
- **Block C** — CURRENT: the **halved paper**, on `paper-final` (off `block-b-full`).
  Collaborators judged the three written Results sections too detailed; they are being cut
  to roughly half, keeping the minimum chain to the closing conclusion, and the figures
  reorganised to match. Before cutting, three **state readings** were added
  (`scripts/45_state_readings.R`, 2026-08-16) because the corpus is built entirely out of
  contrasts and the reduction needed the states:
  - **The diagonal** (`6W_wt>12W_myc`) — never computed; `docs/myc_mouse_finalisation_plan.md
    :379-382` had deferred it as "probably adds little", which is superseded. It needs **no new
    model**: both designs are saturated over the same four groups, so gene-wise
    `cross == myc_12W + 6>12W_wt == myc_6W + 6>12W_myc` to 1.7e-08. **The compartment is up
    almost everywhere (96.5% of 143 MitoPathways, median +0.247) and the respiratory chain is
    the one arm that gains least** — the only drawn arm whose diagonal NES misses BH<0.05
    (2.73/3.28 at the two ages, 1.27 on the diagonal), and demoted on mitoPPS. **"Returns to
    baseline" is relative to the COMPARTMENT, not the transcriptome**: vs 2000
    expression-matched sets it still rises, 95.4th pct.
  - **The absolute levels** — summed DESeq2-normalised counts per tier, four states, per
    animal. Reproduces the content claim on a ruler with no denominator (betas agree with
    script 32 at r=1.000 but are **larger on 11 of 11 panels**, median offset +0.048 — the
    share divides by a denominator that itself rises with Myc, so script 32's "lower bound"
    now has a number). 13 mtDNA genes carry ~21x the transcript of the 152 nuclear OXPHOS
    genes with a **flat genotype bracket (p=0.75)**.
  - **The two abundance rulers disagree on the load-bearing arm.** Same 87 OXPHOS subunits:
    **+0.061 unweighted, +0.201 expression-weighted, +0.226 as a summed-count ratio.** The gap
    IS the weighting. The gradient beats a matched null at the 100th pct (null median has the
    OPPOSITE sign) — but the five lowest expressers are **tissue-restricted paralogs**
    (Cox6a2/Cox7a1/Cox8b/Cox4i2/Cox6b2, ranks 83-87), so their fall may be less contamination
    rather than less respiration: dropping them gives **+0.105**, which belongs beside +0.061.
    The gradient is **not wild-type-only** — as strong in myc_6W (+0.39) as in 6>12W_wt (+0.33).
  - **Rule that follows: a sentence quoting an OXPHOS abundance number must name its ruler.**
  - **Script 32** (`32_mito_content_proxies.R`) — answers "is the biogenesis claim about
    mitochondrial CONTENT?" in absolute compartment shares. **Myc raises mitochondrial
    content ~20–27%**; survives adjustment for prep-stress + contamination. Not a figure
    script. **Figure scripts are now 39+** (34/35/36/37 are analysis scripts, 38 planned — see below).
  - **Script 33** (`33_mtdna_axis_and_coupling_null.R`) — **read this before using ANY
    mt-transcript number.** (1) The mt-transcript share **is not contamination** (mixing
    would need a contaminant mt-share >100%; and nuclear MEC genes track contamination
    *negatively* while mt% tracks it *positively* — mixing cannot give opposite signs). It
    tracks a dominant axis marked by the dissociation IEG signature (rho +0.87 at 6W).
    **Whether mt% is MITOCHONDRIAL is NOT settled** — and do NOT resurrect the withdrawn
    "direction test" (that nuclear mito genes anticorrelate with mt% is *equally* what a
    real mitonuclear imbalance predicts, so it is evidence for neither; circular, withdrawn
    2026-07-17). **PART A does not refute the imbalance.** (2) **The per-sample death couplings fail a null**: script
    25 Part B's couplings to pro_comp (bio_comp 0.83, MASC 0.78, imbalance 0.67) sit at the
    54th–90th percentile of couplings to all 884 library sets (perm p 0.10–0.50). At n=12
    everything correlates with everything at ~0.6–0.8. **Any sample-level composite
    coupling needs an empirical null — the AP6 logic of script 21, which the per-sample
    analyses never got.** The structural reason: within 6W, **PC1 = 43% of variance and
    every measure loads on it** (IEG −0.71, proliferation +0.62, contamination −0.62, MYC
    +0.59) — these are not independent signals. The axis is **mixed** (contamination cannot
    be biological; proliferation can), so it is neither cleanly technical nor cleanly
    biological. It is *not* the TEB/proliferative programme (IEG ~ proliferation −0.76 at
    6W; IEGs *lower* at 6W — both TEB predictions fail).
  - **Script 34** (`34_death_priming_reassessment.R`) — **read this before writing ANY death
    sentence.** Answers the author's challenge that script 33's null is computed *within 6W*
    and so cannot see the *longitudinal* claim. Adds the nulls that can. Outcome: **the death
    arm has no independent transcriptomic support.**
    - **`pro_comp` IS a MitoCarta set** (`MITOCARTA_APOPTOSIS_PRO`, 25 genes, `08:285-290`), and
      so is the imbalance — so "the mito state couples to death priming" correlated **mito genes
      with mito genes**. Every death coupling in the corpus is exposed to this.
    - **"Myc raises apoptotic priming" (geno p=0.037/0.038) is the mito-CONTENT effect.** vs an
      expression-matched **MitoCarta** background, APOPTOSIS_PRO rises **+0.196 vs +0.364
      (z=−2.42, p=0.012)** — i.e. *below average* for a mito gene. Independently reproduces
      script 23's mitoPPS de-prioritisation (−0.21, padj 0.022). On the **404 non-mito**
      pro-death genes the effect does not reproduce (**−0.106, p=0.078**, wrong sign).
    - **The priming interaction was FITTED, SAVED, and NEVER READ**: `death_timing_substrate.rds
      $h1$state_tests` has carried **PRO int_p=0.177 / priming int_p=0.263** all along. The
      narrative quoted group means ("~3x more at 6W" = +0.48 vs +0.15) instead.
    - **Gate 2's p=0.023 is the sole positive and it does not survive**: `mu=0` is the wrong null
      under a global attenuation (matched null z=−1.27, p=0.21); measured inter-gene cor 0.08 →
      effective n **8.4, not 23**; CAMERA (its own question, correlation-aware) p=0.298, FDR=0.926.
    - **The attenuation is GLOBAL**: nothing — death, OXPHOS, or MYC-core — attenuates more than
      expression-matched genes (all BH p>0.14). So **priming fades because the Myc programme
      fades** (canonical Evan/Lowe BH3-only biology); no death-specific or mitochondrial gate.
    - **The Tang NES −1.38 is a SHARED developmental decline** (`temporal_neg` = **−1.31**, 95% as
      much) → cancels in the genotype gap by Issue #6's identity. Never cite `temporal_pos` alone.
    - **The one nuance that is NOT flat:** the *collapse* (r 0.75→0.07) sits at the **96th
      percentile** of Δrho vs all 884 sets, **perm p=0.066**, Fisher z p=0.055 — **marginal, not
      null**, and better than the 6W coupling ever was. But it is **Pearson-only** (Spearman
      p=0.31), fails multiplicity, and **evaporates on every non-mito death axis** (p=0.19–0.92).
      The only axis that gets this far is the circular one. Ambient context: the median set's
      coupling already falls **0.71 → 0.35** (74% → 11% of sets above |r|=0.5) while **PC1 is
      unchanged (43% → 46%)**.
  - **Script 35** (`35_ambient_corrected_couplings.R`) — **THE CORRELATION CEILING. Read before
    quoting any rho.** Every composite here correlates with almost every other, and nobody had
    measured how much. Median |rho| of an axis to an **arbitrary** one of the 884 library sets:
    **mito_oxphos 0.82 (6W) / 0.78 (12W); nucleotide 0.82/0.71; tca 0.80/0.78** — but **redox
    0.41/0.23**. So **a published coupling of 0.8 is what the window hands you.**
    - **The ceiling is mostly MYC DOSE** (the 6-vs-6 genotype split within each timepoint), with
      the IEG/prep axis adding at 6W. **`redox` is the control**: not a Myc target (d=0.08), so
      half the ceiling — **the axes Myc does NOT drive are the only readable ones.**
    - **"6W is stressed, 12W is clean" is WRONG** — script 34's 0.71→0.35 was the mitonuclear
      **imbalance**, a mitoPPS **ratio** (a ratio cancels the shared component). Per axis the
      ceiling barely moves. **Never borrow one axis's ambient for another.**
    - **Q4 (Issue #3 Q2) — half holds.** A raw ambient is the wrong null (Q2's claim is the
      **partial**), so 35 builds a **partial ceiling**. **"OXPHOS is CENTRAL" FAILS**: partial
      →prolif **+0.59 vs a partial ceiling of 0.45** (68th pct, p=0.32); →TEB +0.44 vs 0.45
      (p=0.51). OXPHOS retains residual co-variation with *everything*. **"biogenesis is a
      BYSTANDER" HOLDS, harder than published**: +0.11 vs its ceiling 0.22 (27th pct); →TEB −0.03
      (7th pct) — **below** an arbitrary programme. **The only near-specific couplings are the
      non-Myc ones**: cholesterol→MB2 fork **+0.40 vs a ceiling of 0.13** (96th pct, p=0.043),
      redox→fork −0.50 vs 0.22. **1 of 96 tests beats BH<0.05** — nothing is significant at n=24.
    - **Q5:** Issue #2's "endogenous Myc gates the WT programme" (r=**0.93**) sits at the **92nd
      percentile of its own ceiling (0.75), p=0.08 — marginal**, and is n=6/cell disaggregated.
      Issue #2's **powered half (transgene amplification, a genotype main effect) is untouched.**
    - **The five-question synthesis is `docs/2026-07-18_narrative_synthesis_five_questions.md`** —
      verdicts, numbers, named PDFs, and the gene-set appendix. Read it before writing any figure.
  - **Script 36** (`36_linear_pathway_coupling.R`) — **THE COUPLING METHOD OF RECORD. Read before
    quoting any coupling.** Implements the author's two method notes (`docs/GSVA_global_background_
    and_pathway_coupling.md`; `docs/Gene_set_quantification_for_pathway_correlation.md`). The
    ceiling is **one global per-sample common-mode axis**: `cor(mito_oxphos, mean of all 884
    scores)=0.968`; regress that mean out and the ambient drops 0.80→0.13 (~chance). NOT low n
    (chance floor 0.13), NOT curation (random pairs 0.54), NOT the genotype split (0.83 within
    `6W_pos`). **Methodological in FORM (GSVA doesn't centre within-sample), biological+technical in
    ORIGIN** — the residual (design-removed) axis tracks epithelial +0.79 / contamination −0.70
    (composition/prep) AND proliferation +0.65 / MYC +0.44 (biological). **"Ceiling is mostly Myc
    dose" (script 35's framing) is REFINED**: the RAW ceiling is part four-group contrast (PC1~group
    r=0.48), part a within-group common-mode axis — script 36 decomposes them (PART C).
    - **Method:** primary quantifier is now the **linear mean gene-wise z-score** (its correlation
      IS average cross-gene covariance); **singscore** = cohort-independent robustness; **GSVA** =
      continuity. Correction = **design→PC1→residual, factor estimated AFTER the design** (doc 1 §7)
      + **leave-pair-out** (§4). **NOT mean-subtraction** (equal-loading + sum-to-zero artifact — the
      script-35-era instability). Gene-level vs score-level design-adjustment agree to ~0.02.
    - **THE INTERPRETABILITY TRAP:** when an axis IS the global factor (mito_oxphos r²=0.94 with the
      global mean), removing it leaves ~4% noise — its global-adjusted −0.62 is **noise-on-noise,
      NOT antagonism**. Read a global-adjusted coupling ONLY where `resid_frac ≥ ~0.10`.
    - **RESULT (re-ranks script 35):** "OXPHOS is central" is not just at the ceiling, it is
      **UNTESTABLE** (no separable variance) — the cleanest statement of its failure. **`redox → MB2
      fork` is the ONE robust testable lead** (global-adj −0.63, same sign −0.62…−0.73 across all 3
      quantifiers, axis r²=0.03, Myc-independent). **`cholesterol/mevalonate → fork` is DEMOTED** to
      +0.14 ("separable but weak") — its script-35 96th-pct headline was a raw-ambient artifact.
      biogenesis-bystander unchanged. **Ranking, not proof** (n=24; nothing beats BH<0.05). §0.5/§0.6
      of the synthesis doc carry the write-up + hypothesis ledger + the sentence-level fix pattern.
    - **INFERENCE WEIGHT: ranking + global-factor-as-phenotype, NOT confirmatory CIs** (at n=24 a
      residual r~0.29 arises under the null). **SCOPE unchanged:** re-ranks EXPLORATORY couplings
      only; the content result (raw count shares) and attenuation (DESeq2 LFCs) never touch GSVA/
      these couplings. Author runs 36 (Option A); needs `singscore` (added to `00_setup_packages.R`).
  - **Script 37** (`37_pathway_loading_and_technical_resolution.R`) — **BEHIND the global axis,
    Parts 1+2 of the author's low-dimensional-structure doc** (`docs/Low_dimensional_pathway_
    structure_and_global_axis_analysis.md`). Treats the axis as a phenotype to dissect, not noise to
    delete. Author-run 2026-07-19 (`results/pathway_loading.rds`).
    - **Low-dimensional, quantified:** linear-z-score PC1 = **76%** of variance, effective
      dimensionality **1.7** of a possible 23 (n=24), 79% one-signed loadings. `cor(mito_oxphos,
      global mean)=0.968` (reproduces 36). The 884 pathways move along one dominant axis.
    - **THE POSITIVE REFRAME (answers "OXPHOS in proliferation"):** as a composite, **OXPHOS is the
      single strongest-loading programme — it IS the dominant axis** (r=0.968). The loading LEVEL is
      the ROBUST, defensible statement (near-tautological → a *description* of the axis, not a
      coupling). `redox` is a weak loader (42nd pct), the one axis Myc does not drive.
    - **The loading-CHANGE is a contrast but only an EXPLORATORY lead:** `S ~ G*genotype*timepoint`
      makes loading-change a genotype contrast (readable where the residual coupling was untestable).
      OXPHOS hints at a genotype increase (rel +71%) but **near-collinear at n=6/group and
      mito_biogenesis shows none — not a result, a lead.**
    - **Both-lenses technical correction:** design-only axis (primary) vs contamination+IEG-removed
      (2nd lens) correlate r=0.57. **NAMED over-correction risk:** the RNA-derived anchors absorb
      **r²=0.67 of design-residual proliferation variance** (IEG↔prolif −0.75), so the 2nd lens
      removes some biology; design-only stays primary, prep-driven loadings flagged
      (`r2_removed_by_tech`). No external QC metadata exists → `G_tech` is an RNA surrogate, not
      measured technical.
    - **SCOPE unchanged:** loading layer only; contrasts (content shares, DESeq2 LFCs, Wald fGSEA)
      untouched. Ranking + loading-as-phenotype, NOT confirmatory CIs (n=24).
  - **Script 38** (`38_mitopps_priming_and_pgc1a_biogenesis.R`) — **the death-priming + PGC1a
    questions, moved into spaces where they are answerable.** Author-run 2026-07-19
    (`results/mitopps_priming_pgc1a.rds`).
    - **PART A — mitoPPS space IS readable** (validates the author's point 3): the pairwise ratio
      cancels the global shift, so the mitoPPS ambient is **0.42** (OXPHOS↔global-mean 0.46, PC1 44%)
      vs GSVA's **0.80**. The honest null for a mitoPPS coupling is ~0.42, not ~0.8.
    - **PART B — OXPHOS ↔ death-PRIMING in mitoPPS ratio space is a WEAK but CLEAN lead.** Apoptosis-
      PRO/-ANTI + OXPHOS are nuclear-encoded ratio pathways (script 08 PART 2c; no mt- genes → no
      mtDNA-arm border). OXPHOS↔priming (PRO−ANTI) = raw Spearman **+0.45** (null 0.42), design-adj
      **+0.70** (design-adj null 0.48) — but only the **52nd percentile** of OXPHOS's couplings
      (mid-pack); genotype δ-interaction marginal (p=0.08). **NEITHER circular NOR ceiling-bound**
      (unlike script 34's GSVA-space negative control) — a ranking lead for **BH3 profiling**, not a
      result. priming ≠ death (molecular substrate; phenotype = external IHC/BH3).
    - **PART C — MYC- and PGC1a/ESRRA-biogenesis are NOT separable.** PGC1a = ESRRA/NRF1/GABPA lanes
      (coactivator proxy). vs MYC: raw r=**0.93**, both strongly Myc-induced (β +0.57/+0.75). The
      global-adjusted **−0.03 is the script-36 OVER-REMOVAL trap** (resid_frac 0.01/0.06 — both ARE
      the global factor), NOT independence. **The hypothesised MYC-vs-PGC1a death contrast (PART D)
      is UNTESTABLE here** (all biogenesis axes → priming identically, all resid_frac<0.10) — it
      needs an **ESRRA/PGC1a perturbation** vs the Myc substrate, not n=24 bulk. That is the answer.
  - **The axis is genotype-INDEPENDENT (p=0.49) but time-associated (p=0.0003).** So every
    genotype/interaction result (Issues #1/#2/#3/#5/#6, the content claim) is safe; the
    mt-related TIME story (script 24's mitonuclear narrative, 22/29's mtDNA rise, 25's
    death couplings) is exposed and unresolved. Technical vs biological is **not
    adjudicable** without dissociation-batch/viability metadata.

## Workflow — "Option A" (do not deviate)

- Claude Code **writes and edits** the numbered pipeline scripts. It does **not run
  them.** The author sources them in Positron interactively.
- Infrastructure tasks (git, snapshotting files, provenance READMEs, editing this
  file, planning docs) Claude Code may execute directly. The line: numbered analysis
  scripts are hand-run by the author; plumbing is not.
- Every numbered script ends with an `if (FALSE) { ... }` sandbox block — skipped by
  `source()`/`Rscript`, run line-by-line in Positron for inspection.
- Plan first: in a planning session produce the build spec and wait for approval
  before writing scripts. Commit per verified phase — git is the safety net.
- When in doubt, ask. The cost of clarifying is small; the cost of destroying
  ambiguous state is large.

## Branches, worktrees, git discipline

Branch model (updated 2026-07-12 after the Block A revision consolidation):
- `new-analysis` — stable trunk (scripts 00-12). Unchanged.
- `analysis-exploratory` — **reviewed Block A, scripts 00-31** (fast-forwarded to include
  the revision Issues #1-6). Tag `block-a-reviewed` anchors the reviewed commit (4b82a82).
- `paper-figures` — **Block B, FROZEN 2026-08-14** at tag `block-b-full` (0cb07c2), pushed.
  Scripts 32-44 and the 25-panel `figures/panels/` layer plus `paper/analysis_record.qmd`.
  Do not add to it: it is the named "full version" the reduction is measured against.
- `paper-final` — **Block C, CURRENT working branch** (created off `block-b-full`). The
  halved paper. First work landed 2026-08-16: `scripts/45_state_readings.R` and the three
  state-reading panels (see below). Not yet pushed.
- `BlockA-revision-step-by-step` — the ad-hoc revision branch, now subsumed by
  `analysis-exploratory`; retained as a safety ref, delete once the author is satisfied.
- `main` — earlier pipeline, historical reference only.

Worktrees:
- This working directory — the active `new-analysis`-derived branch.
- `../myc_mouse_main` — clean checkout of `main`, sharing this repo's `.git`. (There
  is no `../myc_mouse_new` worktree.)

Git rules:
1. Do not check out `main` in this working folder. Use the `../myc_mouse_main`
   worktree or `git show main:<path>`. Switching among the `new-analysis`-derived
   branches (`new-analysis`, `analysis-exploratory`, `paper-figures`, `paper-final`)
   here is fine.
2. Do not rename, move, or delete scripts without consulting `docs/branch_manifest.md`
   and confirming with the author. Prefer `git mv`.
3. Read-only git ops (`git show`, `git diff`, `git log`, `git ls-tree`) are always
   fine. Stop-and-check before any destructive action; never force-push a shared
   branch.

## Project structure

- `scripts/` — numbered R pipeline (00-45), `archive_main_pipeline/`, and
  `R_CODING_INSTRUCTIONS.md`. **45 is the newest** (`45_state_readings.R`, the diagonal
  contrast, the absolute levels, and the expression-weighting decomposition).
- `data/` — inputs, plus the library snapshot in `data/genesets_from_library/`.
- `docs/` — `myc_mouse_finalisation_plan.md`, `branch_manifest.md`, `gene_set_history.md`.
- `results/` — intermediate `.rds` (gitignored, generated at runtime).
- `outputs/` — figures/tables (gitignored, generated at runtime).
- `functions/` — shared utilities (`generate_heatmap.R`, `reconcile_gene_symbols.R`).
- `figures/` — the exploratory figure layer (`fig01`-`fig05`, `figS1`-`figS8`, four assembled
  `figure*.R` composites) plus **`figures/panels/`, the publication layer**: one script per
  manuscript panel, `_panel_common.R` (declared scheme + helpers), `rebuild_panels.R`,
  `panels_to_pdf.R`, and **`PANELS.md`, the manifest of record** for which slot each panel holds.
- `paper/` — Quarto. **`analysis_record.qmd` is the record of the analysis as it stands**
  (2026-08-14, extended 2026-08-16): the written Results verbatim, every panel with its legend block, the
  retracted corpus, the open questions. It re-runs no analysis and transcribes no number — it
  sources each panel with `myc.fig.capture = TRUE` and renders `legend_md()`, so it cannot
  drift from the panels. `archive/myc_mito.qmd` is the superseded July walkthrough.
- `external/mitotyping/` — Monzel et al. 2025 mitoPPS reference code (reference only,
  not code to edit).

## Gene sets — consume the library snapshot, do not rebuild

- Finished sets live in `data/genesets_from_library/` — a snapshot of
  `mammary_geneset_library` **v1.0** (nine category GMTs mouse+human, master GMT,
  provenance CSV, `metabric_sets.rds`, `gray_chea_mito_tf_shortlist.csv`), pinned to
  the `v1.0` tag. Provenance in `data/genesets_from_library/README.md`.
- **Do not rebuild gene sets.** Load the GMTs. Each set carries a method tag (`fgsea`,
  `gsva`, or `both`) — route it to the method its tag names. The GSVA large-set
  size-filter exemption applies as recorded in the library decisions.
- The AP6 comparator panel already exists in-repo as `gene_sets_list.rds` (89 sets
  incl. 50 MSigDB Hallmark) — use it, do not re-fetch Hallmark.
- The legacy loads in `01_load_data.R` (`mitocarta_pathways.csv`,
  `myc_signature_genesets.gmx`) predate the library; routing them through the snapshot
  is a pending Step 1/Step 2 tidy (see the plan's open decisions).

## Pipeline structure — Step 1 / Step 2 split

- **Step 1** = data load + DESeq2 (basic + interaction) + QC (basic and interaction).
- **Step 2** = pathway analyses (fGSEA, mitoPPS, GSVA, cell death, characterisation),
  built on Step 1 outputs + the library GMTs.
- Keep both the interaction and the group-based analyses. The split lets Step 2 be
  swapped without re-running DESeq; do not fold set loading or pathway logic back into
  the Step 1 scripts.

## Key methods and analytical conventions

- DESeq2 design: `~ timepoint * myc_status`, five interaction + two group contrasts
  (script 03 is the central hub; changes there propagate downstream).
- **Myc dose is stable across 6W-12W** (mRNA + blot + literature): no pubertal ramp,
  no adult decline. Constant driver, moving substrate. Do not write code or comments
  assuming a dose ramp.
- **fGSEA** — Wald-statistic ranked lists (`stat` = LFC/lfcSE, from the unshrunken MLE
  fit); reads as absolute transcriptional enrichment.
- **mitoPPS** (Monzel 2025, `external/mitotyping/`) — pairwise ratio-based; reads as
  relative resource reallocation within the mitochondrial compartment. Uses
  **linear-scale DESeq2 normalised counts** (not VST). Gene-to-pathway mapping from
  MitoCarta 3.0 Sheet 4 via `splitstackshape::cSplit`, matching the Monzel source. mt-*
  genes (13 mtDNA-encoded) sit in a separate synthetic pathway, never merged into
  nuclear OXPHOS.
- **GSVA** — **log-scale** input (log-CPM or VST, `kcdf="Gaussian"`; or raw counts,
  `kcdf="Poisson"`). The **opposite** of mitoPPS; do not mix them up. Score all samples
  in one run (cohort-relative). GSVA adds no statistical power; it reframes onto
  programs.
- Ortholog mapping — use the cached biomaRt table (`results/ortholog_table.rds`) with
  `dplyr::inner_join`, not `deframe()` + named-vector lookup. Exceptions: `GS_metabolic`
  (own pre-mapped sheet) and `cell_death_genes_consolidated` (MyGene.info + HomoloGene)
  — do not remap either.

## Shrunken vs raw LFCs — important

- fGSEA is **unaffected** (ranks on the Wald stat from the unshrunken fit).
- ΔLFC-difference methods and averaged-LFC visuals must run on **raw (unshrunken)**
  LFCs. This covers the main-09-style cell-death binomial test (metric
  `myc_6W_log2FC - myc_12W_log2FC`) — point it at the existing
  `results/combined_df_annotated_raw.rds` — and the pathway-average heatmaps
  (`outputs/heatmaps_int/` already holds `_raw` and `_shrunk` variants; use `_raw`).
- Before trusting `interaction_gene_characterisation.rds` (script 11) for the PRO/ANTI
  and direction analyses, confirm those parts used raw / Wald-based values — a
  code-read check, no execution.
- Use shrinkage (apeglm/ashr) only for MA-plot QC and ranking-agnostic visualisation.

## Cell death — both branches, neither canonical

**Read `scripts/34_death_priming_reassessment.R` first — both branches are null on the
interaction, and the death composites are mito-defined (circular).** Two standing traps:
`pro_comp`/`priming` are built from `MITOCARTA_APOPTOSIS_PRO`/`_ANTI`, so any coupling to a
mito axis is mito-vs-mito; and branch 1's `binom.test(..., p = 0.5)` is the **wrong null**
under a global attenuation (every Myc-induced gene has ΔLFC>0 by construction — use the
genome-wide supporting fraction, script 34 PART G). Also: **`results/death_timing_substrate.rds`
on disk is STALE** — script 23 as committed builds an 11-row `convergence` table (`23:489-530`);
the saved object has 9, so Parts 2b/5b have never been run. Re-source 23 before reading it.

Both go into the final analysis; report convergence and divergence, do not pick one.
- Branch 1 (main `09`, `archive_main_pipeline/09_cell_death_pathway_summary.R`, ported
  forward) — binomial directional test on `cell_death_genes_consolidated`, on **raw**
  ΔLFC. Directional-hypothesis lens.
- Branch 2 (new `12`) — fGSEA across the Tang 2024 15 RCD modalities (`data/cell-death/`)
  on the Wald z, per contrast including the interaction. Broad-enrichment lens.

## Data traps (a glob of `data/` will otherwise bite)

- **Count matrix:** `01_load_data.R` loads `FULL.DAT.csv` (sample-code headers).
  `FULL.DAT_copy.csv` is the same numeric data with group-label headers and different
  column order — not a byte-duplicate.
- **MitoCarta:** `08_mitoPPS_analysis.R` loads `Mouse.MitoCarta3.0.xls` (dotted,
  new-only) from Sheet 4. `Mouse_MitoCarta3_0.xls` (underscore) is a distinct older
  file, different MD5. Load by exact name; never glob `Mouse*MitoCarta*`.

## R coding rules

Canonical file: `scripts/R_CODING_INSTRUCTIONS.md`. Key points:
- Never `print(n=X)` after `head()` (it may coerce a tibble to data.frame and read `n`
  as `na.print`). Use `head(x) %>% print()` or `x %>% print(n=X)` separately.
- Always `dplyr::count()`, not bare `count()`.
- ASCII-only strings; handle latin1/cp1252 CSV encoding at the read stage.
- Dependencies are centralised in `00_setup_packages.R` (loaded by all scripts via
  `source()`); scripts should not add redundant `library()` calls.
- Use `here::here()` for paths.

## Planning and reference docs

- `docs/myc_mouse_finalisation_plan.md` — the plan. Read as authoritative.
- `docs/branch_manifest.md` — two-pipeline / script inventory (dated 2026-04-22; its
  "consolidation ongoing" framing is superseded by the plan and this file's current
  phase).
- `docs/library_reference/` — read-only reference docs from `mammary_geneset_library`
  v1.0: the set-design rationale (roster, fGSEA/GSVA split, Category 7 algebra, TF
  lanes, dev-set catalog, dropped sets, provenance). Consult when a gene set's meaning
  or method matters for a decision. It is **reference for interpreting the GMTs, not a
  build spec** — do not use it to reconstruct sets. See its `README.md` for an index.
- Decisions are recorded as dated markdown notes; newer supersedes older.
