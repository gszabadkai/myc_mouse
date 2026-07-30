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

## The diverging scale (fixed 2026-07-30, project-wide)

Every diverging quantity in the manuscript uses one ramp, `ms_diverging` in
`_panel_common.R`: **`#4A3525` deep espresso brown → `#FAFAFA` stark white → `#2A8A6D` crisp mint
green**. Two rules travel with it:

- **Zero is written `0`**, never `0.0` or `+0.0`. `lab_signed()` enforces it and formats every other
  tick with an explicit sign.
- **White is pinned to zero, and the two arms are scaled independently** where the data are
  lopsided. `heat_fill()` takes either one number (symmetric) or the observed range. On the
  +3.0 / −1.8 range of Fig. 1B a symmetric ramp left every negative effect in the first third of the
  brown, where BMYO at −1.16 was indistinguishable from LASP at −0.20 — the fill had stopped
  carrying magnitude on that side. The cost is that **equal ink does not mean equal magnitude across
  the sign change**, so the asymmetric form is only for panels where a quantitative axis carries the
  magnitude anyway, and the colour bar is shown so the asymmetry is inspectable.

## Built

| slot | script | supports | inputs |
|---|---|---|---|
| **Fig. 1B** | `fig1B_myc_teb_proliferation.R` | "an endogenous Myc program contributed to the pubertal TEB state in WT … the transgene amplified the TEB and proliferation effects and suppressed the BMYO lineage in favor of LHS" + the de-differentiation reading | `myc_endogenous_amplification.rds`, `dev_program_myc_integration.rds`, `gsva_scores.rds` |
| **Fig. S1A** | `figS1A_geneset_library.R` | "a large custom library of ~900 genesets, covering MYC identity … metabolic pathways" | `provenance_table.csv`, `pathway_loading.rds` (assertion), `gsva_scores.rds` (count) |
| **Fig. S1B** | `figS1B_design_contrasts.R` | "timeline: 6W vs 12W of the WT or Myc+ genotypes, or cross sectional: WT vs Myc+ at 6W or 12W" | none (schematic) |
| **Fig. S1C** | `figS1C_mb_fork_specificity.R` | specificity control for Fig. 1B's `Human BRCA-MYC` row: is the resemblance to the MYC arm of the human switch, or to biogenesis generally? | `gsva_scores.rds`, `ap7_mb_fork.rds` (assertion) |

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

**One panel where there were two.** `fig1B_cell_state_composition.R` (a 2×2-per-programme heatmap of
group-mean z-scores) was **deleted**: it plotted the same GSVA z-scores as the contrast panel, one
lens further back. `fig1C_myc_teb_proliferation.R` was `git mv`d to
**`fig1B_myc_teb_proliferation.R`** and now holds the 1B slot. The deleted panel's view is preserved
in the surviving script's sandbox block, because a contrast plot cannot show a *level* and that is
occasionally worth checking.

*Two orphaned PDFs sit in `outputs/figures/panels/` — `fig1B_cell_state_composition.pdf` and
`fig1C_myc_teb_proliferation.pdf`, both from 12:10 on 2026-07-30. Their scripts no longer exist, so
`rebuild_panels.R` will never refresh them and `fig1C_...pdf` in particular is a near-twin of the
current 1B. Delete both.*

**Fig. 1B (the surviving panel)** — the Cleveland dot form is kept; six changes.

0. **Nomenclature.** The gene and every programme scored on these mouse samples is **Myc**, not
   MYC. The single exception is the human anchor, which keeps the human spelling.
0b. **The top block is now `Myc-tumourigenesis`** and gains the **METABRIC MB2 fork** as
   `Human BRCA-MYC`: the MCbiclust switch that stratifies human breast cancer, whose upper fork is
   biogenesis-high, proliferative, luminal-progenitor, ER-negative and MYC-activated. Not yet in the
   narrative. **Which score:** `MB2_UF` is `METABRIC_MB2_HI_CV_GROUP1`, but UF and LF are
   anticorrelated poles of one switch, so the project's fixed convention (`scripts/18:96-101`,
   "fixed by the paper + author") is the **contrast UF − LF**. The author's call (2026-07-30) is
   the **raw upper fork**, because this row is a *resemblance* statement — how far the mouse tissue
   sits toward a human MYC-driven state — not a claim about where a switch is thrown. It is also the
   stronger number (+1.63 p = 0.039 / +1.59 p = 0.005 against the fork contrast's +1.46 p = 0.053 /
   +1.09 p = 0.045); the fork contrast is quoted in the legend block and asserted against
   `ap7_mb_fork.rds$fork_df$MB2_score`. **Fig. S1C carries the specificity test that has to go with
   it.**
1. **The pooled genotype column is gone**, replaced by `myc_6W` and `myc_12W`. It was hiding the
   result: TEB-ductal is **+1.02 SD at 6 weeks and +0.23 SD at 12**, and the pooled effect (+0.63 SD,
   p = 0.14) is their average and reads as a null.
2. **Effect size is encoded twice** — dot position *and* dot fill on the manuscript diverging scale.
   Significance is a small asterisk, not a second fill colour: at n = 6 per cell the effect size is
   the interesting quantity and the p-value is the footnote. This is what lets the panel go back to
   **single column (89 mm)**, where 4.7 SD spans ~15 mm per column and a −0.8 SD dot sits 2.5 mm off
   zero — the fill is what makes its magnitude read. All four columns still share one x scale,
   because the comparison the sentence makes is between a genotype effect (up to +3.0 SD) and a
   development effect (down to −1.8 SD).
3. **Two blocks, not three.** Developmental state and lineage identity were the same question asked
   twice, so they are merged under **lineage identity**, each MEC state followed by the LE composite
   of its own lineage — BMYO with BA, LASP with AP, LHS with HS — and the paired row labelled `-LE`.
4. **Only the LE arm is drawn.** HE and LE are near-mirror images (they anticorrelate at −0.71 /
   −0.76 / −0.22 within a lineage), so the HE rows cost six lines to say one thing. `HE` =
   high-expressing = differentiated end, `LE` = low-**expressing** = lineage-suppressed end (the
   provenance table's gloss is wrong; see
   `docs/library_reference/Gray_et_al_developmental_TFS_selection.md` §2a). **The HE numbers must
   still be reported in the text** — see the bounds below, because the opposition of the two signs is
   what licenses the de-differentiation reading and the drawn arm is the one exposed to the
   common-mode axis.

The Fig. 1B additions were computed **in the figure script**, not in a new numbered analysis
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
- **The de-differentiation reading is the stronger version of that sentence.** Myc raises every
  lineage-suppressed (LE) composite and lowers every differentiated (HE) composite, in all three
  lineages, at both ages, most strongly at 6 weeks (AP +1.51 / −0.82, BA +1.55 / −1.00,
  HS +1.31 / −0.85 SD). **Bound, and it applies to the arm that is on the panel:** the LE composites
  sit almost on the global common-mode axis (r 0.85–0.90 with the mean of all 885 scores) and the
  undrawn HE composites do not (−0.37, −0.57, +0.07); the LE lanes are also larger (165–788 genes vs
  26–159) and carry more mitochondrial genes. What a common-mode axis cannot produce is **opposite
  signs within a lineage**, so the reading rests on the opposition — which means the HE numbers have
  to appear in the Results text even though they are not drawn.
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

**Fig. S1C (new)** — the specificity control the `Human BRCA-MYC` row needs. MB1 and MB2 share a
lower fork and diverge into two upper forks that are *both* biogenesis-high and proliferative; the
only feature distinguishing them is MYC, so MB2_UF against MB1_UF is what separates "Myc drives the
Myc fork" from "Myc drives biogenesis generically".

- **Myc raises both forks, by almost the same amount** — MB2_UF +1.63 / +1.59 SD, MB1_UF +1.73 /
  +1.34 SD — and the two scores **correlate 0.982** (230 shared genes of 626 and 419).
- **Only the difference is drawn** (author, 2026-07-30), as a small boxplot with the animals on it,
  because the difference *is* the test. The two component scores are in the legend block and in the
  script's sandbox. It is in **raw GSVA units**: the difference's within-group SD is 0.082 against
  0.335 and 0.297, about a quarter, so putting it on a standardised axis would make the residual look
  as large as the thing it is a residual of.
- Groups are ordered **age-major** so the two boxes a reader compares are adjacent, with a bracket
  per age carrying an asterisk only where the genotype contrast clears p < 0.05.
- **The cancellation is the point:** MB2_UF and MB1_UF each correlate 0.94 / 0.95 with the global
  common-mode axis; their difference falls to **0.49**. The difference is the only one of the three
  rows that is not largely the common mode.
- **The timing contradicts the prediction the test was built for.** Script 18 predicted MB2
  resemblance would peak at 6 weeks (biogenesis being front-loaded) and reported a *pooled* genotype
  main effect, p = 0.026. Split by age the specificity is **absent at 6 weeks (+0.39 SD, p = 0.59)
  and present at 12 (+1.62 SD, p = 0.003)**, while the resemblance to MB2_UF itself is flat with age.
  The interaction is not significant (p = 0.15) at n = 6 per cell, so this is a discrepancy to
  report, not a reversal to claim. It also runs opposite to the TEB result, which is a 6-week effect.

## Not built here

| slot | why |
|---|---|
| **Fig. 1A** | Whole-mount / histology of the hyperplastic ductal expansion — the author's bench image, assembled outside R. |
| **Fig. 1C onward** | Vacated by the 1C→1B renumber. Paragraph 2's panel (the library ranked by normalised effect size) takes it. Next increment; `fgsea_percategory.rds` and `pathway_loading.rds` are on disk and reconciled. |
| MEC state *levels* | The four groups as a 2×2 grid per programme, group-mean z. Built, reviewed, dropped as redundant with the contrasts; kept in `fig1B_myc_teb_proliferation.R`'s sandbox. |
| Gray HE arm | Computed and asserted, reported in the legend block, not drawn — see point 4 above. |
| sample PCA | Deferred 2026-07-29 — belongs to the global-axis section (alignment along the mitochondrial axis), not to the opening overview. |
| DE counts per contrast | Deferred 2026-07-29 — belongs to the attenuation section. Numbers already checked: 6W genotype contrast = **2777** genes at FDR 10 %, **1967** at FDR 5 %; the `interaction` contrast has **0** at either threshold. |

## Assembly

Deferred until the narrative fixes the numbering. Then a `figures/figure*_overview.R` composes
the signed-off panels with patchwork `tag_levels = "A"` and this table records the mapping. Every
panel is now **single column (89 mm)**, so the layout is unconstrained.

The four already-assembled manuscript figures (`figures/figure1_myc_mitochondrion.R`,
`figure2_developmental_window.R`, `figureS1_compartment_detail.R`, `figureS2_controls.R`) are
untouched by this layer and still rebuild via `figures/rebuild_manuscript_figures.R` — but they now
pick up the new sample palette, so they should be re-rendered and re-checked at the same time.
