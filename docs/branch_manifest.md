# Branch manifest — `main` vs `new-analysis`

Generated: 2026-04-22
Purpose: Read-only structural reference for the ongoing pipeline
consolidation. The working directory holds a mix of files from both
branches; this document is the authoritative map of what belongs where.
	
important: this file is a previous state, here only for reference - HISTORICAL. Superseded by
the 2026-07-12 consolidation: `analysis-exploratory` now holds the reviewed Block A (scripts
00-31, tag `block-a-reviewed` @ 4b82a82) and `paper-figures` (Block B) is branched off it. See
CLAUDE.md "Branches" for the current model.


> **No files were moved, renamed, or deleted to produce this manifest.**
> All claims are derived from `git show`, `git ls-tree`, `git diff`,
> `diff`, `md5`, and reading file headers.

---

## 1. Worktree layout

| Path | Branch | HEAD | Notes |
|---|---|---|---|
| `/Users/gs/G/data/MK_myc_2022/myc_mouse` | `new-analysis` | `653bb16` | Current working directory; scripts 00-12 + `archive_main_pipeline/` |
| `/Users/gs/G/data/MK_myc_2022/myc_mouse_main` | `main` | `7394bb5` | Clean checkout of main |

**Discrepancy with CLAUDE.md**: CLAUDE.md states a third worktree
`../myc_mouse_new` exists as a clean checkout of `new-analysis`. It does
**not** exist on disk. The `new-analysis` branch is checked out in the
primary working directory (`myc_mouse`). Either create the worktree
(`git worktree add ../myc_mouse_new new-analysis`) or update CLAUDE.md.

---

## 2. Scripts comparison

Main has 18 files in `scripts/` (10 numbered 00-09 with **two** `04_*`
scripts, plus 7 unnumbered legacy files). New-analysis has 13 numbered
scripts 00-12 plus an `archive_main_pipeline/` subdirectory that preserves
main's removed scripts verbatim (verified byte-identical via `diff`).

`archived` in the table below means "present in `scripts/archive_main_pipeline/`
on new-analysis". `renamed` means the script moved forward into the new
numbering (and is therefore **not** in the archive).

| Filename | main | new-analysis | Status | One-liner |
|---|---|---|---|---|
| `00_setup_packages.R` | yes | yes (modified) | **differs** | Package install/load. New adds `here`, `tidyr`, `forcats`, `msigdbr`, `readxl`, `splitstackshape`. |
| `01_load_data.R` | yes | yes (modified) | **differs** | Load count matrix + coldata. New uses `here::here()` and caches biomaRt orthologs. |
| `02_deseq_interaction_model.R` | yes | archived | removed → archive | DESeq2 interaction model with IHW-filtered contrasts. Superseded by new `03_deseq_results_qc.R`. |
| `02_qc.R` | no | yes | **new** (**renamed from** `06_diagnostics_QC.R`) | QC: total counts, library complexity, normalisation diagnostics. |
| `03_lfc_classification.R` | yes | archived | removed → archive | LFC classification / gene ranking. |
| `03_deseq_results_qc.R` | no | yes | **new** | Runs DESeq2; extracts 5 interaction + 2 group contrasts; ashr/apeglm MA-plot QC. **Central pipeline hub.** |
| `04_group_comparison.R` | yes | archived | removed → archive | Group-based (non-interaction) DESeq2. |
| `04_fgsea_timepoint_compare.r` | yes | archived | removed → archive | Early fGSEA draft. **Note lowercase `.r` extension.** |
| `04_fgsea_pathway_analysis.R` | no | yes | **new** | Temporal fGSEA (12W-vs-6W in Myc+/Myc-) on Wald statistics. |
| `05_heatmap_utils.R` | yes | archived | removed → archive | Heatmap helpers for LFC variants. |
| `05_fgsea_visualisation.R` | no | yes | **new** | Dotplots / bar plots / enrichment curves for fGSEA. |
| `06_diagnostics_QC.R` | yes | renamed | **renamed → `02_qc.R`** | QC (PCA, library metrics). Renamed; **not** in archive. |
| `06_fgsea_cross_sectional.R` | no | yes | **new** | Cross-sectional fGSEA (Myc+ vs Myc- at each timepoint). |
| `07_heatmap_oxphos_annotated.R` | yes | archived | removed → archive | OXPHOS subunit heatmap. |
| `07_fgsea_xs_visualisation.R` | no | yes | **new** | Visualisation of cross-sectional fGSEA. |
| `08_pathway_summary.R` | yes | archived | removed → archive | Broad MitoCarta pathway summary. |
| `08_mitoPPS_analysis.R` | no | yes | **new** | Pairwise-ratio mitoPPS scoring (Monzel et al. 2025). |
| `09_cell_death_pathway_summary.R` | yes | archived | removed → archive | **Main's cell-death analysis** (binomial test on curated gene list). |
| `09_mitoPPS_vs_fgsea_comparison.R` | no | yes | **new** | Compares fGSEA NES vs mitoPPS per pathway across contrasts. |
| `10_interaction_fgsea_mitopps.R` | no | yes | **new** | Interaction-contrast fGSEA + mitoPPS overlay. |
| `11_interaction_gene_characterisation.R` | no | yes | **new** | Gene-level dissection of interaction genes (volcano, gprofiler). |
| `12_cell_death_pathway_analysis.R` | no | yes | **new** | **New cell-death analysis** (fGSEA across 15 Tang-2024 RCD modalities). |
| `Heatmaps.R` | yes | archived | legacy → archive | Legacy heatmap generation. |
| `k693_rerun_analysis.R` | yes | archived | legacy → archive | k693 reanalysis. |
| `k693_rerun_analysis_GS.R` | yes | archived | legacy → archive | k693 reanalysis variant. |
| `Myc_timecourse_analysis_GS.R` | yes | archived | legacy → archive | Original timecourse analysis. |
| `Myc_timecourse_analysis_GS_sandbox.R` | yes | archived | legacy → archive | Sandbox version. |
| `NES_paradox_analysis.R` | yes | archived | legacy → archive | NES paradox exploration. |
| `rebuild_all.R` | yes | archived | legacy → archive | Pipeline rebuild orchestrator. |

**Archive semantics**: `scripts/archive_main_pipeline/` holds the 15
scripts from main that were *removed outright* in new-analysis. Scripts
that were *updated in place* (`00`, `01`) or *renamed forward* (`06` →
`02_qc.R`) are **not** archived — their history is preserved only via git.

---

## 3. `data/` and `functions/` comparison

All shared files are MD5-identical between the two branch worktrees.

### `data/`

| File | main | new-analysis | Status | Notes |
|---|---|---|---|---|
| `cell_death_genes_consolidated.csv` | yes | yes | identical | Curated cell-death list used by main's 09. |
| `cell_death_genes_consolidated.rds` | no | yes | new-only | Binary companion. |
| `cell-death/` (15 CSVs) | no | yes | **new-only** | Tang et al. 2024 RCD modalities: Alkaliptosis, Apoptosis, Autophagy_dependent_cell_death, Cuproptosis, Disulfidptosis, Entotic_cell_death, Ferroptosis, Immunogenic_cell_death, Lysosome_dependent_cell_death, MPT_driven_necrosis, Necroptosis, NETotic_cell_death, Oxeiptosis, Parthanatos, Pyroptosis. Used by new 12. |
| `coldata.csv` | yes | yes | identical | Sample metadata. |
| `felsher_integrative_signature.csv` | yes | yes | identical | Felsher signature gene set. |
| `FULL.DAT.csv` | yes | yes | identical | Count matrix; columns are sample codes (`MYBS102_1i` etc.). |
| `FULL.DAT_copy.csv` | yes | yes | identical | Count matrix; columns are group labels (`6 Wk MYC neg` etc.). **Same numeric data, different column naming/ordering** — not a byte-duplicate. |
| `FULL.DAT.COL.DATA.txt` | yes | yes | identical | Column metadata for full dataset. |
| `Human_MitoCarta3_0.xls` | yes | yes | identical | MitoCarta 3.0 human reference. |
| `Mouse_MitoCarta3_0.xls` | yes | yes | identical | MitoCarta 3.0 mouse reference. |
| `Mouse.MitoCarta3.0.xls` | no | yes | **new-only** | **Distinct file** — MD5 differs from the underscore variant. Name-clash risk. |
| `mitocarta_pathways.csv` | yes | yes | identical | MitoCarta pathway annotations. |
| `mitocarta_pathways.xlsx` | yes | yes | identical | Excel version. |
| `myc_signature_genesets.gmx` | yes | yes | identical | Myc signature sets (GMX). |
| `OXPHOS_subunits.csv` | yes | yes | identical | OXPHOS subunit list. |
| `signature_go_mouse_cells.xlsx` | yes | yes | identical | GO signature file. |
| `README.md` | no | yes | new-only | Data directory documentation. |

### `functions/`

| File | main | new-analysis | Status | Notes |
|---|---|---|---|---|
| `generate_heatmap.R` | yes | yes | identical | Only shared utility. |

### `external/`

| Path | main | new-analysis | Notes |
|---|---|---|---|
| `external/mitotyping/` | no | yes | Monzel et al. 2025 mitoPPS reference implementation (cloned from `annamonzel/mitotyping`). Used by new `08_mitoPPS_analysis.R`. |

### Top-level

| File | main | new-analysis | Status |
|---|---|---|---|
| `README.md` | yes | yes | identical |
| `CLAUDE.md` | no | yes | new-only |
| `MYC_mouse_analysis_summary_before_revision_20251002.md` | yes | yes | identical |
| `NES_paradox_explanation.md` | yes | yes | identical |
| `myc_mouse.Rproj` | yes | yes | identical |
| `myc_mouse.code-workspace` | yes | no | main-only (VS Code artefact) |

---

## 4. Cell-death analysis — side-by-side

The central scientific comparison for the consolidation.

| Aspect | main `09_cell_death_pathway_summary.R` | new-analysis `12_cell_death_pathway_analysis.R` |
|---|---|---|
| DE input | `results/combined_df_annotated.rds` (pre-computed DE + annotation) | `results/interaction_results.rds` (from new `03`) |
| Ortholog mapping | None (mouse symbols already) | `results/ortholog_table.rds` (biomaRt human→mouse) |
| Gene sets | `data/cell_death_genes_consolidated.csv` — **custom curated** list (categories: apoptosis, CICD; source tags: GO / KEGG / Reactome / Hallmark) | `data/cell-death/*.csv` — **Tang et al. 2024** 15 regulated cell death modalities |
| Statistical method | **Binomial test** on counts of supporting vs opposing genes (threshold \|ΔLFC\| > 0.2) | **fGSEA** (`minSize=5`, `maxSize=2000`, `nPermSimple=10000`) |
| Ranking metric | Δ LFC: `myc_6W_log2FC − myc_12W_log2FC` | Wald z: `log2FC / lfcSE` |
| Contrasts | 2 implicit: Myc effect on pro-death (6W vs 12W), Myc effect on pro-survival (6W vs 12W) | 5 explicit: Myc@6W, Myc@12W, 12-vs-6 in Myc-, 12-vs-6 in Myc+, interaction |
| RDS output | None | `results/cell_death_fgsea.rds` (gene sets, fGSEA per contrast, NES matrices, LE genes, verdicts) |
| Figure output | 2 PDFs (hypothesis_all_categories, cell_death_hypothesis_scatter) | 6 PDFs (heatmap NES, dotplot temporal / cross-sectional / interaction, paired-temporal dotplot, verdict heatmap) |
| Table output | 5 CSVs (hypothesis_results, genes_full, supporting_genes, opposing_genes + console category summary) | 3 CSVs (fgsea_all_results, leading_edge_annotated, interpretation_summary) |
| Interpretation | Directional hypothesis test: is delta significant and in expected direction? | Multi-layered: NES per contrast + LE extraction + pro-/anti-death keyword tagging + "more resistant" / "more susceptible" verdict per modality |

**Summary**: main tests a **specific directional hypothesis** on a
curated list with a simple binomial test. New-analysis runs **unbiased
enrichment** across a published 15-modality panel, with leading-edge
characterisation and per-pathway verdicts. The two approaches ask
overlapping but distinct questions.

---

## 5. New-analysis dependency graph

Arrows show which RDS / CSV outputs flow between scripts. Script 02
produces QC figures only (no downstream consumers); scripts 05, 07, 09,
11 are terminal visualisation / comparison branches.

```
01_load_data.R
├─ count_matrix.rds, coldata.rds, dds_int.rds
├─ gene_sets_list.rds           (89 sets: 22 MitoCarta + 17 MYC + 2 apoptosis + 50 MSigDB Hallmark)
└─ ortholog_table.rds
      └──> 02, 03, 04, 06, 08, 10, 11, 12

02_qc.R (QC figures only — terminal)

03_deseq_results_qc.R
├─ interaction_results.rds       (myc_6W, myc_12W, timepoint_neg, timepoint_pos, interaction)
├─ group_results.rds             (pos_12W_vs_6W, neg_12W_vs_6W)
├─ dds_int_run.rds, dds_group_run.rds
└─ extended_qc_summary.rds, lfc_comparison_timepoint.rds, results_summary.csv
      └──> 04, 06, 08, 09, 10, 11, 12

04_fgsea_pathway_analysis.R
└─ fgsea_results.rds     ──> 05, 09

05_fgsea_visualisation.R         (figures only — terminal)

06_fgsea_cross_sectional.R
└─ fgsea_xs_results.rds  ──> 07, 09

07_fgsea_xs_visualisation.R      (figures only — terminal)

08_mitoPPS_analysis.R
└─ mitopps_scores.rds    ──> 09, 10, 11

09_mitoPPS_vs_fgsea_comparison.R (comparison tables/figures — terminal)

10_interaction_fgsea_mitopps.R
└─ interaction_fgsea_mitopps.rds ──> 11

11_interaction_gene_characterisation.R (characterisation — terminal)

12_cell_death_pathway_analysis.R
└─ cell_death_fgsea.rds          (terminal; reads 01 + 03 + external cell-death CSVs only)
```

### Dependency summary

| Script | Depends on | Feeds |
|---|---|---|
| 01 | external data | 02, 03, 04, 06, 08, 10, 11, 12 |
| 02 | 01 | — (QC only) |
| 03 | 01 | 04, 06, 08, 09, 10, 11, 12 |
| 04 | 01, 03 | 05, 09 |
| 05 | 01, 04 | — (viz only) |
| 06 | 01, 03 | 07, 09 |
| 07 | 04, 06 | — (viz only) |
| 08 | 01, 03 | 09, 10, 11 |
| 09 | 03, 04, 06, 08 | — (comparison only) |
| 10 | 01, 03, 08 | 11 |
| 11 | 01, 03, 08, 10 | — (characterisation only) |
| 12 | 01, 03 + external cell-death CSVs | — (terminal) |

**Hubs**: `03` (feeds 7 downstream) and `01` (feeds 8 downstream). Any
change to the DESeq2 model or contrast extraction in `03` propagates to
almost every downstream script.

**Pipeline is acyclic**; no dangling reads were found.

---

## 6. Ambiguities and surprises

1. **Missing worktree**. CLAUDE.md lists `../myc_mouse_new` as a clean
   new-analysis checkout; it does not exist. The primary working
   directory is itself on `new-analysis`.
2. **Two `04_*` scripts on main**. Main contains both
   `04_group_comparison.R` and `04_fgsea_timepoint_compare.r` (note
   lowercase `.r`). Both are archived on new-analysis.
3. **Rename vs. removal**. Main's `06_diagnostics_QC.R` became
   new-analysis `02_qc.R`. It is **not** in `archive_main_pipeline/`
   because the archive only holds scripts that were *removed*, not those
   updated or renamed forward. Anyone reviewing history of QC logic
   should use `git log --follow` on the new path.
4. **Two similarly-named MitoCarta files** in new-analysis:
   `Mouse_MitoCarta3_0.xls` (shared with main, MD5 `cf5fcded…`) and
   `Mouse.MitoCarta3.0.xls` (new-only, MD5 `805385fe…`). **Different
   content** despite near-identical names — worth clarifying or
   consolidating before the Quarto writeup.
5. **`FULL.DAT_copy.csv` is not a byte-copy** of `FULL.DAT.csv`. Both are
   54,838 lines with the same gene IDs, but `FULL.DAT.csv` has
   sample-code column headers (`MYBS102_1i`, `MYCF52_5f`, …) while
   `FULL.DAT_copy.csv` has group-label headers (`6 Wk MYC neg`, …).
   Column *order* also differs. The numeric payload is the same. The
   naming "_copy" is misleading.
6. **Script 12 is quasi-independent**. It reads only `01` + `03` outputs
   plus the external cell-death CSVs — it skips the whole fGSEA /
   mitoPPS stack (04-11). This makes it safe to re-run in isolation but
   also means it does not benefit from the new pipeline's curated
   gene-set list or mitoPPS scoring. Worth deciding whether cell death
   should stay parallel or be integrated with the mitoPPS axis for the
   final writeup.
7. **Archive is faithful but partial**. `diff` against `git show
   main:scripts/X` confirms the archived copies are byte-identical to
   main. However, because updated (00, 01) and renamed (06) scripts are
   omitted, `archive_main_pipeline/` alone is **not sufficient to
   reconstruct main** — the `main` branch itself remains the source of
   truth.
