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

## The sample palette (fixed 2026-07-30, project-wide)

Set once in `figures/theme_myc.R`, Okabe-Ito, and used **wherever individual sample values are
shown**:

| group | key on disk | colour |
|---|---|---|
| 6W wt | `6W_neg` | `#0072B2` dark blue |
| 12W wt | `12W_neg` | `#56B4E9` sky blue |
| 6W myc | `6W_pos` | `#D55E00` vermilion |
| 12W myc | `12W_pos` | `#E69F00` orange |

Hue is genotype, lightness is age — and lightness runs the **opposite** way to the palette this
replaced (6W is now the saturated end). `geno_cols` is the 6W pair, so a two-level genotype key and
a four-level group key cannot disagree. The four already-assembled manuscript figures consume
`group_cols` by name, so they change colour on their next rebuild; all four were dry-rendered
against the new palette and still build.

Contrast vocabulary, also fixed 2026-07-30 and declared once in `_panel_common.R`:

- **genotype** — `myc_6W`, `myc_12W`
- **development** — `6>12W_wt`, `6>12W_myc`

The legend blocks map these to the `results/interaction_results.rds` slot names
(`6>12W_wt` = `timepoint_neg`, `6>12W_myc` = `timepoint_pos`).

## Built

| slot | script | supports | inputs |
|---|---|---|---|
| **Fig. 1B** | `fig1B_cell_state_composition.R` | "no major changes in the overall composition of major mammary cell states (BMYO, LHS, LASP), but a clear TEB-ductal shift away from the pubertal proliferative state at 12W" | `dev_program_myc_integration.rds`, `myc_endogenous_amplification.rds`, `gsva_scores.rds` |
| **Fig. 1C** | `fig1C_myc_teb_proliferation.R` | "an endogenous Myc program contributed to the pubertal TEB state in WT … the transgene amplified the TEB and proliferation effects and suppressed the BMYO lineage in favor of LHS" + the de-differentiation reading | `myc_endogenous_amplification.rds`, `dev_program_myc_integration.rds`, `gsva_scores.rds` |
| **Fig. S1A** | `figS1A_geneset_library.R` | "a large custom library of ~900 genesets, covering MYC identity … metabolic pathways" | `provenance_table.csv`, `pathway_loading.rds` (assertion), `gsva_scores.rds` (count) |
| **Fig. S1B** | `figS1B_design_contrasts.R` | "timeline: 6W vs 12W of the WT or Myc+ genotypes, or cross sectional: WT vs Myc+ at 6W or 12W" | none (schematic) |

### Revision of 2026-07-30 (author's review of the first four panels)

**Fig. S1B** — contrast labels are now the author's vocabulary; the linetype key says "genotype"
and "development" and **no longer names the batch** (that is a Methods statement and moved to the
legend block). `WT` / `Myc+` pulled in against the scheme by tightening the x limits, with `HALF_W`
reduced so the nodes keep their size on the page. Node ink follows the palette's lightness, not the
genotype: white at 6W, near-black at 12W.

**Fig. S1A** — the fill was the fGSEA/GSVA routing tag; it is now the **mitochondrial definition of
each set** in the three classes of `outputs/pathway_loading/E_library_coverage.pdf`. The routing tag
is minor and belongs in `paper/myc_mito.qmd` (`method_cols` is kept in `_panel_common.R` for that).
The classification is script 37's rule (`37:445-462`) re-derived over the whole 986-set library and
**asserted identical to `pathway_loading.rds$mito_classification` on the 885 sets they share**. The
point the fill carries: the library is 40 % mitochondrial by label but **13 % once the
construction lanes are set aside** — 260 of the 393 mito-labelled sets are Gray `_MITO` /
`_LE_MITO` TF lanes, which are MitoCarta subsets *by build*.

**Fig. 1B** — was 24 points against a group axis, and the trajectories were unreadable. Now a
**square heatmap laid out as Fig. S1B is**: one 2×2 grid per programme, age across, genotype down,
filled by the group-mean z-score in within-group SD. A trajectory is read left-to-right and the Myc
effect top-to-bottom, the same motion in both panels. Fill is PRGn, deliberately not blue-red, so
the fill cannot be mistaken for genotype. The per-animal spread the tiles average over is kept in
the script's sandbox.

**Fig. 1C** — two structural changes and one addition.
1. **The pooled genotype column is gone**, replaced by `myc_6W` and `myc_12W`. It was hiding the
   result: TEB-ductal is **+1.02 SD at 6 weeks and +0.23 SD at 12**, and the pooled effect (+0.63 SD,
   p = 0.14) is their average and reads as a null.
2. **Full width (183 mm), one shared x scale.** All four columns share the SD scale because the
   comparison the sentence makes is between a genotype effect (up to +3.0 SD) and a development
   effect (down to −1.7 SD); at 120 mm a −0.6 SD point sat 3 mm off zero and read as nothing.
3. **Block 3 is new — the Gray HE/LE lineage-identity axes**, six composites (AP, BA, HS × HE, LE)
   built here as the mean over `TFT_<TF>_GRAY_<lineage>_<HE|LE>`, with the `_MITO` promotions
   excluded. `HE` = high-expressing = differentiated end, `LE` = low-**expressing** =
   lineage-suppressed end (the provenance table's gloss is wrong; see
   `docs/library_reference/Gray_et_al_developmental_TFS_selection.md` §2a).

The Fig. 1C additions were computed **in the figure script**, not in a new numbered analysis
script, because the figure layer already fits exactly these contrasts (`fig01_mito_content.R:82-84`,
`figS1_mitocarta_survey_share.R`, `figS2_oxphos_complex_share.R`, `figS8_myc_network_levels.R`) with
script 27's idiom (`27:103-106`). The shared machinery — `composite_of()`, `wsd_of()`,
`contrast_table()` — now lives in `_panel_common.R` so 1B and 1C cannot drift apart, and every
number that also exists in `results/` is asserted against it (`prog_stats$wt_d`, `$mycpos_shift`,
`$int_beta`, `$within_sd`; `state_stats$d6`, `$d12`).

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

## Where the panels and the written sentences disagree

Carried in each script's `LEGEND` block and repeated here because it needs a decision.

- **"the transgene amplified the TEB … effects" — RESOLVED 2026-07-30, in the author's favour.**
  The pooled genotype main effect on TEB-ductal is +0.63 SD, p = 0.14, which is why the first
  version of Fig. 1C could not support the sentence. Split by age it is **+1.02 SD at 6 weeks**
  (p = 0.10) and **+0.23 SD at 12** (p = 0.70). The Myc effect on the TEB axis exists where the
  phenotype is and has gone by 12 weeks; the pooled test was averaging it away. Still n = 6 per
  cell, so p = 0.10 — a 1 SD effect, not a significant one.
- **"no major changes in the overall composition"** — in within-group SD, the wild-type 6→12W shift
  is BMYO **+0.10**, LASP **−0.70**, LHS **−0.60**, TEB-ductal **−0.91**. BMYO is flat; the two
  luminal states fall by two thirds of what the TEB axis does. Script 26's own PART C2 header says
  the same ("the data do NOT support a BMYO expansion (flat); the WT change is luminal (LASP/LHS)
  DECLINE"). All four are non-significant at n = 6 (p 0.24–0.86).
- **"a clear TEB-ductal shift"** — not supported by the sample-level test (**p = 0.24** in wild
  type). What licenses "clear" is the matched-null ruler: TEB-vs-ductal (HS) **−0.422 at the 0th
  percentile** (`docs/2026-07-26_introduction_alignment_and_the_question.md` §0.7). In Myc+ the
  timeline shift *is* significant (**p = 0.001**).
- **"suppressed the BMYO lineage in favor of differentiation to LHS"** — the BMYO half **holds** and
  deepens with age (−0.80 SD at 6W, −1.16 SD at 12W, p = 0.066). The LHS half is a **reversal, not a
  main effect**: −0.63 SD at 6W, +0.79 SD at 12W (pooled p = 0.86, interaction p = 0.097). The
  sentence is true at 12 weeks and false at 6.
- **The de-differentiation reading is the stronger version of that sentence.** Myc lowers every
  differentiated (HE) composite and raises every lineage-suppressed (LE) composite, in all three
  lineages, at both ages, most strongly at 6 weeks (AP −0.82 / +1.51, BA −1.00 / +1.55,
  HS −0.85 / +1.31 SD). **Bound:** the LE composites sit almost on the global common-mode axis
  (r 0.85–0.90 with the mean of all 885 scores) and the HE composites do not (−0.37, −0.57, +0.07);
  the LE lanes are also larger (165–788 genes vs 26–159) and carry more mitochondrial genes. What
  a common-mode axis cannot produce is **opposite signs within a lineage**, and that opposition is
  what the reading rests on.
- **This panel set is not the attenuation result.** On the MYC programme the genotype effect is the
  same at both ages (+2.36 vs +2.28 SD), so in cohort-relative GSVA space there is nothing to
  attenuate. The attenuation result is a DESeq2 effect-size result (×0.55, scripts 29–31 and 40).
  What attenuates *here* is the developmental arm — TEB and every lineage axis — not the MYC arm.
  Worth saying explicitly in the text, because the two rulers will otherwise look contradictory.

### Two text numbers to make exact

- **"~900 genesets"** is four defensible counts: **986** in the library snapshot, **885** actually
  scored per sample by GSVA (902 were routed; 17 fell to the size filter), **817** fGSEA-eligible,
  **885** in the coupling and loading analyses. 986 for "the library", 885 for "scored".
- **"~3K DEGs"** at 6W is **2777** at FDR 10 % and **1967** at FDR 5 %; the `interaction` contrast
  has **0** at either threshold.

## Not built here

| slot | why |
|---|---|
| **Fig. 1A** | Whole-mount / histology of the hyperplastic ductal expansion — the author's bench image, assembled outside R. |
| **Fig. 1D** | The library ranked by normalised effect size (paragraph 2). Next increment; `fgsea_percategory.rds` and `pathway_loading.rds` are on disk and reconciled. |
| sample PCA | Deferred 2026-07-29 — belongs to the global-axis section (alignment along the mitochondrial axis), not to the opening overview. |
| DE counts per contrast | Deferred 2026-07-29 — belongs to the attenuation section. Numbers already checked: 6W genotype contrast = **2777** genes at FDR 10 %, **1967** at FDR 5 %; the `interaction` contrast has **0** at either threshold. |

## Assembly

Deferred until the narrative fixes the numbering. Then a `figures/figure*_overview.R` composes
the signed-off panels with patchwork `tag_levels = "A"` and this table records the mapping. Note
that Fig. 1C is now **double-column (183 mm)** while 1B and the S1 panels are single-column, which
constrains the layout: 1C wants its own full-width row.

The four already-assembled manuscript figures (`figures/figure1_myc_mitochondrion.R`,
`figure2_developmental_window.R`, `figureS1_compartment_detail.R`, `figureS2_controls.R`) are
untouched by this layer and still rebuild via `figures/rebuild_manuscript_figures.R` — but they now
pick up the new sample palette, so they should be re-rendered and re-checked at the same time.
