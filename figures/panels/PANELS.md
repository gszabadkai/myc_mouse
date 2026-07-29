# Panel manifest

One script per panel under `figures/panels/`, each built to support a specific sentence of the
Results section. **Figure numbers are recorded here, not in filenames** — the Google Doc's
numbering is still moving, and this table is a one-line edit where a rename is not.

Narrative source: Google Doc `MK_Myc_Paper`, tab **Gyorgy Writing**.

## Convention

- No explanatory text on the page. `theme_panel()` blanks title, subtitle and caption.
- Every script carries a `panel_legend()` block; `rebuild_panels.R` collects them into
  `outputs/figures/panels/legends.md`, which is the raw material for the figure legends.
- The composite object is always named `p`.
- PDFs go to `outputs/figures/panels/`; `options(myc.fig.nosave = TRUE)` gives a dry run.

## Built

| slot | script | supports | inputs |
|---|---|---|---|
| **Fig. S1A** | `figS1A_geneset_library.R` | "a large custom library of ~900 genesets, covering MYC identity … metabolic pathways"; "fGSEA ranking and GSVA scoring" | `data/genesets_from_library/provenance_table.csv` |
| **Fig. S1B** | `figS1B_design_contrasts.R` | "timeline: 6W vs 12W of the WT or Myc+ genotypes, or cross sectional: WT vs Myc+ at 6W or 12W" | none (schematic) |

## Next, blocked on a re-run

Both read objects that predate the 2026-07-24 gene-symbol reconciliation. **Author action,
Option A — re-source in Positron, in order: `scripts/17_gsva_overview.R` →
`scripts/26_dev_program_myc_integration.R` → `scripts/27_myc_endogenous_amplification.R`.**
Script 17 is the slow one (`limma::mroast`, `nrot = 9999`, 902 sets). The panels call
`require_fresher_than()` and stop rather than draw a stale number.

| slot | script | supports | inputs |
|---|---|---|---|
| **Fig. 1B** | `fig1B_cell_state_composition.R` | "no major changes in the overall composition of major mammary cell states (BMYO, LHS, LASP), but a clear TEB-ductal shift away from the pubertal proliferative state at 12W" | `dev_program_myc_integration.rds`, `gsva_scores.rds` |
| **Fig. 1C** | `fig1C_myc_teb_proliferation.R` | "an endogenous Myc program contributed to the pubertal TEB state in WT … the transgene amplified the TEB and proliferation effects and suppressed the BMYO lineage in favor of LHS" | `myc_endogenous_amplification.rds`, `dev_program_myc_integration.rds` |

## Not built here

| slot | why |
|---|---|
| **Fig. 1A** | Whole-mount / histology of the hyperplastic ductal expansion — the author's bench image, assembled outside R. |
| **Fig. 1D** | The library ranked by normalised effect size (paragraph 2). Next increment; `fgsea_percategory.rds` and `pathway_loading.rds` are on disk and reconciled. |
| sample PCA | Deferred 2026-07-29 — belongs to the global-axis section (alignment along the mitochondrial axis), not to the opening overview. |
| DE counts per contrast | Deferred 2026-07-29 — belongs to the attenuation section. Numbers already checked: 6W genotype contrast = **2777** genes at FDR 10%, **1967** at FDR 5%; the `interaction` contrast has **0** at either threshold. |

## Assembly

Deferred until the narrative fixes the numbering. Then a `figures/figure*_overview.R` composes
the signed-off panels with patchwork `tag_levels = "A"` and this table records the mapping. The
four already-assembled manuscript figures (`figures/figure1_myc_mitochondrion.R`,
`figure2_developmental_window.R`, `figureS1_compartment_detail.R`, `figureS2_controls.R`) are
untouched by this layer and still rebuild via `figures/rebuild_manuscript_figures.R`.
