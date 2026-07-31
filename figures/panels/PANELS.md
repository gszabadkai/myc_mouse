# Panel manifest

One script per panel under `figures/panels/`, each built to support a specific sentence of the
Results section. **Figure numbers are recorded here, not in filenames** — the Google Doc's
numbering is still moving, and this table is a one-line edit where a rename is not.

Narrative source: Google Doc `MK_Myc_Paper`, tab **Gyorgy Writing**.

## Convention

- **Filenames do not carry the slot letter** (2026-07-31). A panel script is
  `fig1_<what it is>.R` or `figS1_<what it is>.R`; which letter it holds is the Built table
  below and nothing else. The letters were in the filenames until the S1A/S1B swap made every
  renumber a rename, a `git mv`, a stranded PDF and four string edits. `rebuild_panels.R`
  derives the PDF slug from the filename, so the `save_panel_p()` slug must match it; and
  because name order is no longer figure order, `legends.md` is now sorted by the `slot` each
  legend block declares (`slug_slot_order()` in `_panel_common.R`).
- No explanatory text on the page. `theme_panel()` blanks title, subtitle and caption.
- Every script carries a `panel_legend()` block; `rebuild_panels.R` collects them into
  `outputs/figures/panels/legends.md`, which is the raw material for the figure legends.
- The composite object is always named `p`.
- PDFs go to `outputs/figures/panels/`; `options(myc.fig.nosave = TRUE)` gives a dry run.

## The sample palette (fixed 2026-07-30, project-wide)

Set once in `figures/theme_myc.R`, Okabe-Ito, and used **wherever individual sample values are
shown**:

| drawn label | key on disk | colour |
|---|---|---|
| `6W_wt` | `6W_neg` | `#0072B2` dark blue |
| `12W_wt` | `12W_neg` | `#56B4E9` sky blue |
| `6W_myc` | `6W_pos` | `#D55E00` vermilion |
| `12W_myc` | `12W_pos` | `#E69F00` orange |

The **drawn labels are the author's own names** (fixed 2026-07-31): `group_labels` in
`theme_myc.R` prints `6W_wt` / `12W_wt` / `6W_myc` / `12W_myc`, not a prettified `6W WT`. The keys
stay the on-disk factor levels (`neg`/`pos`) because that is what the data carry.

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
| **Fig. 1B** | `fig1_myc_teb_proliferation.R` | p1: "an endogenous Myc program contributed to the pubertal TEB state in WT … Transgenic Myc activation amplified the TEB and proliferation state, suppressed the BMYO lineage, and strongly promoted a lineage suppressed (LE) de-differentiation pattern" | `myc_endogenous_amplification.rds`, `dev_program_myc_integration.rds`, `gsva_scores.rds` |
| **Fig. 1C** | `fig1_pathway_pca.R` | p2: "the dominant principal component axis of these genesets across the whole dataset …" — that there IS one, and where the four groups sit on it | `gsva_scores.rds`, `pathway_loading.rds` (assertions) |
| **Fig. 1D** | `fig1_nes_ranking.R` | p2: "Mitochondrial biogenesis and OXPHOS complex genesets were overall top ranked based on normalised changes in effect size, along with core (non-mitochondrial) Myc target genesets, followed by …" | `fgsea_percategory.rds`, `ap6_permutation_null.rds` (legend), `provenance_table.csv` |
| **Fig. S1A** | `figS1_design_contrasts.R` | p1: "both longitudinal (6W versus 12W for WT or Myc+ genotypes) and cross-sectional (WT versus Myc+ at 6W or 12W) comparisons" | none (schematic) |
| **Fig. S1B** | `figS1_geneset_library.R` | p1: "fGSEA ranking and GSVA scoring based on custom-curated genesets"; p2: "we extended the custom library to a total of 986 genesets" | `provenance_table.csv`, `pathway_loading.rds` (assertion), `gsva_scores.rds` (count) |
| **Fig. S1C** | `figS1_axis_loadings.R` | p2: "… aligned almost entirely with variability in mitochondria related terms" — what the axis is made of, with its bound | `gsva_scores.rds`, `pathway_loading.rds` (`mito_classification`, `mito_enrichment`) |
| **Fig. S1D** | `figS1_mb_fork_specificity.R` | not cited yet — specificity control for Fig. 1B's `Human BRCA-MYC` row: is the resemblance to the MYC arm of the human switch, or to biogenesis generally? | `gsva_scores.rds`, `ap7_mb_fork.rds` (assertion) |

### Revision of 2026-07-31 (paragraph 1 rewritten, paragraph 2 written)

**The numbering moved and the letters left the filenames.** Paragraph 1 as rewritten cites the
design schematic before the library, so **S1A and S1B swapped**; the MB-fork panel is cited
nowhere yet and **rolled down to S1D** to make room for the enrichment ranking paragraph 2 does
cite. Rather than rename for the letters a second time, the letters came out of the filenames —
see the Convention section.

**Fig. 1C and Fig. S1C are one PCA, drawn twice.** `pathway_axis()` in `_panel_common.R` rebuilds
script 37's linear z-score score matrix and its pathway-level PCA — `results/pathway_loading.rds`
saves the per-set loadings and the variance percentages but not the score matrix or the sample
scores — and then proves the rebuild before returning: 885 sets, PC1/PC2/PC3 = 76.538 / 8.027 /
4.254 % to 1e-6, `cor(mito_oxphos, global mean)` = script 37's own saved gate (0.9720) to 1e-6,
and concordance ≥ 0.9 with the saved design-residual ranking (observed **0.968**). Runs in 0.5 s
and has no RNG in it.

**Pathway-level PCA only** (author). Script 37 draws the gene-level PCA beside it, because the
77 % is a property of the *composites* (gene-level PC1 is 44 %, and is a composition axis the
compositing correctly drops). That is a methods point; it is in Fig. 1C's legend block and is not
drawn, and the supplementary version of it is explicitly not wanted yet.

**The one thing that makes Fig. 1C readable: PC1 is a genotype axis, not a time axis.**
`PC1 ~ timepoint * genotype` gives genotype **p = 0.0031**, timepoint p = 0.70, interaction
p = 0.43. Since `batch = timepoint`, an age-separating axis would have been uninterpretable. Group
means: 6W Myc+ **+15.2**, 12W Myc+ +4.4, 6W WT −8.4, 12W WT −11.2 — so the genotype gap is 23.6 at
6 weeks and 15.6 at 12, the same *direction* as the attenuation result and not that result (the
interaction is not significant, and attenuation is a DESeq2 effect-size result).

**Fig. 1D's ruler is the fGSEA NES**, which is what "normalised changes in effect size" means in
the manuscript's own methods wording ("enrichment … normalised to the whole transcriptome"). The
source is script 20's per-category run, `fgsea_percategory.rds`. It is a **different ruler** from
1C/1D on purpose: those are per-sample covariation, this is the genotype contrast on the
unshrunken Wald statistic. Both genotype contrasts are drawn because the claim is about a ranking
and the ranking is age-stable (Spearman **0.93**, median |NES| 2.18 vs 2.17) — which is also the
observation the attenuation section opens with.

**One encoding is new and is declared, not inlined:** `programme_group()` in `_panel_common.R`,
with `mito_class3()` (script 37's classification rule, `37:445-462`) lifted out of
`figS1_geneset_library.R` so the three panels that need it share one copy.

### Author's review, same day — four changes

**1D and S1C swapped.** The enrichment ranking is the more relevant main-figure panel: it measures
what Myc *did* (the genotype contrast on the unshrunken Wald statistic), where 1C and the loadings
measure what covaries with what. So `fig1_nes_ranking.R` holds **Fig. 1D** and
`figS1_axis_loadings.R` holds **Fig. S1C**; the files were `git mv`d to match, which is the last
time a renumber will touch a filename now that the letters are gone.

**THE EXPORT DEVICE CHANGED, and it affects every figure in the project.** The author reported
uneven character spacing. It is `cairo_pdf`: it writes text as individually positioned glyphs, and
at the 6 pt these panels use the inter-word advance rounds away — `"Hallmark comparator"` renders
as `"Hallmarkcomparator"` and `"nucleotide"` as `"nudeotide"` in Preview and every other Core
Graphics viewer, at every font family tried. The PDF's *text* is correct (`pdftotext` reads back the
spaces), so it is a rendering artefact, but it is the one the reader sees. `save_panel()` in
`theme_myc.R` now uses the base `pdf()` device with Helvetica, which writes real strings with font
metrics and is clean at the same size. **The trade:** Helvetica is one of the 14 standard PDF fonts
and is referenced rather than embedded (no Ghostscript on this machine, so `embedFonts()` is not
available). `options(myc.fig.cairo = TRUE)` switches back if a submission demands full embedding.
All four assembled manuscript figures dry-render unchanged through the new device.

**Fig. 1D, five row names are the author's** and one of them shifts work onto the colour key: the
by-construction block is now labelled **"curated mitochondrial"** — which is what those sets *are* —
and the key that colours it is what says they are MitoCarta subsets *by build*. Also
`Myc targets (curated)` → **Myc signatures**, `mitochondrial, other` →
**mitochondrial metabolism & dynamics** (the 37 remaining MitoCarta pathways: amino-acid, lipid,
carbohydrate and one-carbon metabolism, fission and surveillance, chaperones and proteases,
transport, calcium and signalling), `TF target lanes` → **TF target sets**, `Hallmark comparator` →
**MSigDB Hallmarks**. (The key's placement was settled in the second review below: one line, under
the plot.)

**Fig. S1C's labels are placed, not repelled.** Five of the seven anchors sit in the left third of
the ranking and a label is ~200 rank units wide, so any left-hand column puts one label across the
next label's leader — which is what the first version did. The wedge *above* the curve on the right
is provably free (the curve falls monotonically, so a straight line from a high-rank anchor to a
point up and to the right can never re-cross it), and both the anchors and the label slots are
monotone in y, so the fan cannot cross itself. A `stopifnot` enforces rank order = label order.
Content changes: `Myc lane x MitoCarta` → **Myc regulon (MitoCarta subset)**; MB1 → **MB2**
(`METABRIC_MB2_HI_CV_GROUP1`, the same score Fig. 1B draws as Human BRCA-MYC); `nucleotide` dropped;
and the biogenesis anchor moved from `MITOCARTA_MITOCHONDRIAL_RIBOSOME` to
`MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA` so the label **mitochondrial biogenesis** names the set it
points at rather than one member of the programme.

### Author's second review, same day — four more

**Group labels are the author's names, project-wide.** `group_labels` now prints `6W_wt`,
`12W_wt`, `6W_myc`, `12W_myc`. One edit in `theme_myc.R`, so it lands on every panel and on the four
assembled figures at their next rebuild.

**The three mitochondrial classes are renamed, also project-wide** (`mito_class_labels`):
**MitoCarta** / **curated (with mitochondria)** / **non-mitochondrial**. Shorter, and the shortening
is what lets Fig. 1D's key fit **one line under the plot** — it costs one line of height instead of
three and cannot collide with a point. Fig. S1B and Fig. S1C pick the names up automatically.

**Fig. 1C's break is off the ticks and both segments are cropped to their data.** The first version
put the break mark on the `+20` tick, which read as a fault rather than a break. Now the main
segment runs to +22.6 with its break in the empty run past `+20`, and the outlier segment is 3 units
wide with its tick at `+44` and its break to the left of it. The panel is 54 mm tall, not 42: with
`coord_fixed` the height is what makes the plot fill 89 mm of width, and shortening the *axis* (52
drawn units against 69.5 unbroken) is what the break bought.

**Fig. S1C's labels are under the curve, close to their points, in three lanes** (y = 0.80, 0.60,
0.44). The lanes are the whole trick and the arithmetic is in the script: a leader dropping to a
deeper lane crosses a shallower label unless that label *ends before the next anchor begins*, so
alternating lanes reduces the constraint to "each label is narrower than the gap to the next
anchor" — which is why three labels are wrapped over two lines and one is right-aligned. The
wrapping is load-bearing, not decoration. `stopifnot` checks every label sits below its own point.

**Fig. 1C has a broken x axis.** One 6-week Myc+ animal sits at PC1 = +44.3 while the next highest
is +19.1, so an unbroken axis spends a third of the width on nothing. Built as two panels rather
than with a package (`ggbreak` is not installed, and doing it by hand keeps `coord_fixed(1)` alive
in both segments, which is the part that matters — within a segment a PC1 unit and a PC2 unit are
the same length, so the collapse onto PC1 is still the picture and not the panel shape). Segment
widths are proportional to their spans and the break is marked with slashes. **What it costs:** a
distance read *across* the gap is not to scale. Two traps found while building it, both recorded in
the script: `clip = "off"` (needed for the break marks) makes each segment draw *all* 24 animals
over its neighbour unless the data are filtered per segment; and patchwork does not merge two
separately built guides, so only the main segment may carry the key. x ticks are by 10.

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
`_LE_MITO` TF lanes, which are MitoCarta subsets *by build*. The key sits **inside the panel, bottom
right**, in the wedge the sorted bars leave empty — the smallest four categories occupy under a
tenth of the width there — which took the panel from 64 mm to 52 mm high.

**One panel where there were two.** `fig1B_cell_state_composition.R` (a 2×2-per-programme heatmap of
group-mean z-scores) was **deleted**: it plotted the same GSVA z-scores as the contrast panel, one
lens further back. `fig1C_myc_teb_proliferation.R` was `git mv`d to
**`fig1B_myc_teb_proliferation.R`** and now holds the 1B slot. The deleted panel's view is preserved
in the surviving script's sandbox block, because a contrast plot cannot show a *level* and that is
occasionally worth checking.

*The two orphaned PDFs left by the rename — `fig1B_cell_state_composition.pdf` and
`fig1C_myc_teb_proliferation.pdf` — were deleted on 2026-07-30. `rebuild_panels.R` only ever writes
the slugs its scripts name, so a rename or deletion always strands the old PDF; check
`outputs/figures/panels/` against the Built table after either.*

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
   `ap7_mb_fork.rds$fork_df$MB2_score`. **Fig. S1D carries the specificity test that has to go with
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
- **"This effect diminished in both the WT and Myc+ 6W → 12W trajectory, suggesting an endogenous
  suppression of the Myc program"** (paragraph 1, closing) — true of the *trajectory* and not of the
  genotype effect, which is the point immediately above. Paragraph 2's attenuation sentence inherits
  the same ambiguity, and Fig. 1C now shows the third version of it: the genotype gap along PC1 goes
  23.6 → 15.6 with an interaction p of 0.43. Three rulers, one direction, one significant result
  (the DESeq2 one). The text should name which ruler it is on.

### Added 2026-07-31, for paragraph 2

- **"aligned almost entirely with variability in mitochondria related terms"** — supported at the
  top of the ranking and bounded three ways on the face of Fig. S1C. Mito sets are 44 % of the
  library and **100 % of the top 50**, median |loading| 0.90 against 0.61. But **94 of the top 100
  are build-tautological** `_MITO` lanes (MitoCarta subsets by construction); non-mitochondrial
  growth programmes reach **0.98** (the METABRIC biclusters, nucleotide metabolism, pentose
  phosphate, the Myc and E2F regulons); and **redox, a mitochondrial gene set, loads 0.67** in the
  bottom half — and it is the one axis Myc does not drive. Suggested wording: **mitochondria-led**.
  What the axis tracks is a coordinated Myc anabolic–proliferative–mitochondrial growth state.
- **"revealing the leading role of mitochondrial remodeling"** — loading is **covariation, not
  primacy**, and script 37 says so in as many words (`pathway_loading.rds$enrichment_verdict`). At
  n = 24 nothing on these panels can order cause; the primacy claim belongs to the perturbations.
- **Fig. 1D's tier order is right at the top and inverted in the middle.** OXPHOS **+2.42**, Myc
  signatures **+2.38** and mitochondrial biogenesis **+2.37** are a three-way tie at the top, exactly
  as written. But the sentence continues "followed by biosynthetic metabolic pathways, E2F and
  overall proliferation signaling", and in the data **E2F / cell cycle (+2.08) sits above
  biosynthetic metabolism (+1.84)**. One-word fix: swap them.
- **The curated-mitochondrial row tops Fig. 1D** at **+2.47**, above OXPHOS. It is drawn as its own
  row in its own colour rather than deleted, because deleting it would flatter every row beneath it.
  The rows to read as programmes start below it — worth one clause in the legend, not in the Results.
- **Apoptosis is the one row that changes sign** (+0.74 at 6 weeks, −0.94 at 12; every other row's
  median moves *up*), and the author asked what drives it. **It is not a result.** The row is bimodal
  and its median sits on zero, so five weak sets crossing zero move it: `TANG_NECROPTOSIS`,
  `TANG_LYSOSOME_DEPENDENT_CELL_DEATH`, `APOP_INTRINSIC_REACTOME`, `APOP_MODULATION_WP` and
  `TANG_PYROPTOSIS`, each going from about +0.8 to about −1.3, and **none significant at either
  age**. The members that *are* significant do not move: `APOP_REGULATION_REACTOME` stays +2.08 →
  +2.43, `TANG_CUPROPTOSIS` +1.86 → +1.93, while `APOP_HALLMARK` (−1.95 → −1.37) and `APOP_KEGG`
  stay depleted and if anything rise. Three of the five movers are non-apoptotic death modalities,
  so there is no mitochondrial-death reading in it either — consistent with script 34, which found
  the death arm has no independent transcriptomic support.
- **Three set counts are now in play and each panel uses a different one.** **986** in the library
  (Fig. S1B, and the number paragraph 2 gives), **885** scored per sample by GSVA (Figs. 1C and S1C),
  **866** in the enrichment ranking (Fig. 1D: 816 fGSEA-eligible library sets that clear the
  minSize 10 / maxSize 500 filter, plus 50 fresh Hallmark comparators). 986 is correct where it
  stands; the other two need naming where they are cited.

### One text number still to make exact

- **"~3K DEGs"** at 6W is **2777** at FDR 10 % and **1967** at FDR 5 %; the `interaction` contrast
  has **0** at either threshold.
- *(Resolved 2026-07-31: "~900 genesets" is now written as 986 in paragraph 2, which is the library
  count and correct. See the three-counts entry above.)*

**Fig. S1D** (built 2026-07-30 as S1C; rolled down 2026-07-31 because it is not cited yet) — the
specificity control the `Human BRCA-MYC` row needs. MB1 and MB2 share a
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
| **Fig. 1E onward** | Paragraph 2's closing sentence — the partial (Myc-signal-independent) correlation of the library to proliferation, de-differentiation, apoptotic priming and the human MYC-driven tumour signature. Next increment; `ambient_corrected_couplings.rds`, `linear_pathway_coupling` outputs and `mitopps_priming_pgc1a.rds` are on disk. Read `docs/2026-07-18_narrative_synthesis_five_questions.md` §0.5/§0.6 first: the raw ambient is the wrong null, "OXPHOS is central" is untestable rather than false, and `cholesterol → fork` was demoted. |
| Gene-level vs pathway-level PCA | Author's call 2026-07-31: pathway lens only, to keep the figure from turning into a methods argument. The comparison (gene PC1 = 44 %, a composition axis; the Myc programme rides at gene-PC2) is in Fig. 1C's legend block, and script 37's figures A2/A3/A4 hold the drawn version if a supplementary is ever wanted. |
| Scree / dimensionality | Same reason. The numbers that matter (PC1 76.5 %, effective dimensionality 1.7, 79 % one-signed loadings) are in Fig. 1C's legend block. |
| MEC state *levels* | The four groups as a 2×2 grid per programme, group-mean z. Built, reviewed, dropped as redundant with the contrasts; kept in `fig1_myc_teb_proliferation.R`'s sandbox. |
| Gray HE arm | Computed and asserted, reported in the legend block, not drawn — see point 4 above. |
| sample PCA | Deferred 2026-07-29 — belongs to the global-axis section (alignment along the mitochondrial axis), not to the opening overview. |
| DE counts per contrast | Deferred 2026-07-29 — belongs to the attenuation section. Numbers already checked: 6W genotype contrast = **2777** genes at FDR 10 %, **1967** at FDR 5 %; the `interaction` contrast has **0** at either threshold. |

## Assembly

Deferred until the narrative fixes the numbering. Then a `figures/figure*_overview.R` composes
the signed-off panels with patchwork `tag_levels = "A"` and this table records the mapping. All
but one panel is **single column (89 mm)** wide (Fig. S1D is 50 mm), so the layout is
essentially unconstrained. Heights as built: 1B 84, 1C 54, 1D 68, S1A 52, S1B 52, S1C 56,
S1D 52 mm — so Figure 1's three built panels stack in 199 mm and Supplementary 1's four in 212,
both about one page column.

The four already-assembled manuscript figures (`figures/figure1_myc_mitochondrion.R`,
`figure2_developmental_window.R`, `figureS1_compartment_detail.R`, `figureS2_controls.R`) are
untouched by this layer and still rebuild via `figures/rebuild_manuscript_figures.R` — but they now
pick up the new sample palette, so they should be re-rendered and re-checked at the same time.
