# CLAUDE.md — MMTV-Myc Mouse Timecourse Analysis

Rules and context for Claude Code sessions in this repo. Read this first, then the
plan: `docs/myc_mouse_finalisation_plan.md`.

## What this repo is

Bulk RNA-seq of MMTV-Myc transgenic mice. Four groups x six replicates
(`6W_neg`, `6W_pos`, `12W_neg`, `12W_pos`), DESeq2 genotype x timepoint interaction
model. Investigates Myc-driven transcriptional and mitochondrial changes across the
two timepoints. Feeds the manuscript "Mitochondria integrate oncogenic and metabolic
transcriptional programs to shape breast cancer progression". Analytical companion
paper: Menegollo, Bentham et al., Cancer Res 2024 (CAN-23-3172).

## Current phase

Two-pipeline consolidation is done: `new-analysis` is the trunk (scripts 00-12);
main's removed scripts are archived byte-identical in `scripts/archive_main_pipeline/`.
Active work is the finalisation in `docs/myc_mouse_finalisation_plan.md`, built in two
phases:
- **Block A** — broad exploratory: gates + full action-point menu across all
  sets/methods/contrasts, plus the raw re-runs. Review against the hypotheses, then
  decide what the paper says.
- **Block B** — focused: re-code only the narrative-selected analyses and figures for
  the writeup. The later submission-quality pass, not the current draft.

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

Branch model:
- `new-analysis` — stable trunk (scripts 00-12).
- `analysis-exploratory` — Block A active development (off `new-analysis`). Current
  working branch; broad runs and raw re-runs land here.
- `paper-figures` — Block B, created later off the reviewed Block A commit. Do not
  create until after the Block A review.
- `main` — earlier pipeline, historical reference only.

Worktrees:
- This working directory — the active `new-analysis`-derived branch.
- `../myc_mouse_main` — clean checkout of `main`, sharing this repo's `.git`. (There
  is no `../myc_mouse_new` worktree.)

Git rules:
1. Do not check out `main` in this working folder. Use the `../myc_mouse_main`
   worktree or `git show main:<path>`. Switching among the `new-analysis`-derived
   branches (`new-analysis`, `analysis-exploratory`, `paper-figures`) here is fine.
2. Do not rename, move, or delete scripts without consulting `docs/branch_manifest.md`
   and confirming with the author. Prefer `git mv`.
3. Read-only git ops (`git show`, `git diff`, `git log`, `git ls-tree`) are always
   fine. Stop-and-check before any destructive action; never force-push a shared
   branch.

## Project structure

- `scripts/` — numbered R pipeline (00-12), `archive_main_pipeline/`, and
  `R_CODING_INSTRUCTIONS.md`.
- `data/` — inputs, plus the library snapshot in `data/genesets_from_library/`.
- `docs/` — `myc_mouse_finalisation_plan.md`, `branch_manifest.md`, `gene_set_history.md`.
- `results/` — intermediate `.rds` (gitignored, generated at runtime).
- `outputs/` — figures/tables (gitignored, generated at runtime).
- `functions/` — shared utilities (`generate_heatmap.R`).
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
