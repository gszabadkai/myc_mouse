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
| **Fig. 1B** | `fig1B_cell_state_composition.R` | "no major changes in the overall composition of major mammary cell states (BMYO, LHS, LASP), but a clear TEB-ductal shift away from the pubertal proliferative state at 12W" | `dev_program_myc_integration.rds`, `myc_endogenous_amplification.rds`, `gsva_scores.rds` |
| **Fig. 1C** | `fig1C_myc_teb_proliferation.R` | "an endogenous Myc program contributed to the pubertal TEB state in WT … the transgene amplified the TEB and proliferation effects and suppressed the BMYO lineage in favor of LHS" | `myc_endogenous_amplification.rds`, `dev_program_myc_integration.rds` |
| **Fig. S1A** | `figS1A_geneset_library.R` | "a large custom library of ~900 genesets, covering MYC identity … metabolic pathways"; "fGSEA ranking and GSVA scoring" | `data/genesets_from_library/provenance_table.csv` |
| **Fig. S1B** | `figS1B_design_contrasts.R` | "timeline: 6W vs 12W of the WT or Myc+ genotypes, or cross sectional: WT vs Myc+ at 6W or 12W" | none (schematic) |

### The 17 → 26 → 27 re-run, done 2026-07-29 15:33

`dev_program_myc_integration.rds` (Jul 15) and `myc_endogenous_amplification.rds` (Jul 8) predated
the 2026-07-24 gene-symbol reconciliation that rebuilt `gsva_scores.rds`, and script 26 also reads
`gsva_overview.rds` (Jul 6) at `26:68`. All three re-sourced with `set.seed(1)`; **every number in
`state_stats` and `prog_stats` is unchanged to three significant figures**, which is what the
reconciliation's scope predicted — its large recovery was ATP synthase / OXPHOS, and these panels
read mammary-development, MYC and proliferation sets. The panels call `require_fresher_than()` and
stop rather than draw a stale number, so the guard stays useful for the next such rebuild.

**One reproducibility note.** `limma::mroast` (`17:379`) is rotation-based and script 17 does *not*
seed it — its only `set.seed` is inside `perm_trajectory_contrast` (`17:495`). The roast columns of
`gsva_overview$coef_table` therefore differ slightly between runs. They feed nothing in these panels
(which use the OLS coefficients and the GSVA scores), but a `set.seed()` before line 379 would be a
one-line durable fix if the roast p-values are ever quoted.

### Where the panels and the written sentences disagree

Both 1B and 1C draw what the objects say, which is not in every respect what paragraph 1 currently
says. Carried in each script's `LEGEND` block and repeated here because it needs a decision:

- **"no major changes in the overall composition"** — in within-group SD, the wild-type 6→12W shift
  is BMYO **+0.10**, LASP **−0.70**, LHS **−0.60**, TEB-ductal **−0.91**. BMYO is flat; the two
  luminal states fall by two thirds of what the TEB axis does. Script 26's own PART C2 header says
  the same ("the data do NOT support a BMYO expansion (flat); the WT change is luminal (LASP/LHS)
  DECLINE"). All four are non-significant at n=6 (p 0.24–0.86).
- **"a clear TEB-ductal shift"** — not supported by the sample-level test (**p = 0.24** in wild type).
  What licenses "clear" is the matched-null ruler: TEB-vs-ductal (HS) **−0.422 at the 0th percentile**
  (`docs/2026-07-26_introduction_alignment_and_the_question.md` §0.7). In Myc+ the timeline shift
  *is* significant (**p = 0.001**).
- **"the transgene amplified the TEB … effects"** — TEB-ductal has **no** genotype main effect
  (+0.63 SD, **p = 0.14**). Proliferation does (+1.01 SD, p = 0.019) and MYC-in-TEB does (+1.61 SD,
  p = 6.2e-04), but the latter is MYC *targets* scored in the TEB context, not the TEB phenotype.
- **"suppressed the BMYO lineage in favor of differentiation to LHS"** — the BMYO half **holds** and
  is the powered state result (−0.98 SD, p = 0.023). The LHS half does not: no genotype main effect
  (p = 0.86), but a time-dependent flip (interaction p = 0.097; below wild type at 6W, above at 12W).

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
