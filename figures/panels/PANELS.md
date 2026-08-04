# Panel manifest

One script per panel under `figures/panels/`, each built to support a specific sentence of the
Results section. **Figure numbers are recorded here, not in filenames** — the Google Doc's
numbering is still moving, and this table is a one-line edit where a rename is not.

**The slot map is the section "The slot map (2026-08-04, read off the written Results)" below.**
It supersedes every earlier commission in this file, including the Figure 2 / Figure 3 table of
`docs/2026-08-02_results_narrative_final_order_D.md`. Everything above it is the built work and
the record of how it got its form.

Narrative source: Google Doc `MK_Myc_Paper`, tab **Gyorgy Writing**. Read in full on
**2026-08-04**, down to the `READ UP TO HERE ONLY` marker — four section headings, the last of
which (`Restoring OXPHOS by PGC-1a re-establishes the apoptotic trigger`) is a stub for the cell
work. The slot map below is read off *that* text, not off the pre-writing commission.

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
| **Fig. 1B** | `fig1_myc_teb_proliferation.R` | p1: "an endogenous MYC program contributed to the pubertal TEB state in WT … Transgenic MYC activation amplified the TEB and proliferation state, suppressed the BMYO lineage, and strongly promoted a lineage suppressed (LE) de-differentiation pattern" | `myc_endogenous_amplification.rds`, `dev_program_myc_integration.rds`, `gsva_scores.rds` |
| **Fig. 1C** | `fig1_pathway_pca.R` | p2: "the dominant principal component axis of these genesets across the whole dataset …" — that there IS one, and where the four groups sit on it | `gsva_scores.rds`, `pathway_loading.rds` (assertions) |
| **Fig. 1D** | `fig1_nes_ranking.R` | p2: "Mitochondrial biogenesis and OXPHOS complex genesets were overall top ranked based on normalised changes in effect size, along with core (non-mitochondrial) MYC target genesets, followed by …" | `fgsea_percategory.rds`, `ap6_permutation_null.rds` (legend), `provenance_table.csv` |
| **Fig. 1E** | `fig1_mito_content.R` | p3: "it drove a quantitative expansion by systemically upregulating 94% of 143 nuclear-encoded mitochondrial pathways, increasing the mitochondrial transcriptomic fraction by 21–27% (padj < 0.001)" — **both** clauses, two rulers, two parts | `mito_content_proxies.rds` (`$shares`, `$share_stats`), `background_vs_myc.rds` (`$ruler`, `$ruler_summary`) |
| **Fig. S1A** | `figS1_design_contrasts.R` | p1: "both longitudinal (6W versus 12W for WT or Myc+ genotypes) and cross-sectional (WT versus Myc+ at 6W or 12W) comparisons" | none (schematic) |
| **Fig. S1B** | `figS1_geneset_library.R` | p1: "fGSEA ranking and GSVA scoring based on custom-curated genesets"; p2: "we extended the custom library to a total of 986 genesets" | `provenance_table.csv`, `pathway_loading.rds` (assertion), `gsva_scores.rds` (count) |
| **Fig. S1C** | `figS1_axis_loadings.R` | p2: "… aligned almost entirely with variability in mitochondria related terms" — what the axis is made of, with its bound | `gsva_scores.rds`, `pathway_loading.rds` (`mito_classification`, `mito_enrichment`) |
| *(no slot)* | `figS1_mb_fork_specificity.R` | **displaced 2026-08-04.** It held S1D, and the written S1D is the MYC western blot. Still cited nowhere: it is the specificity control for Fig. 1B's `Human BRCA-MYC` row. It keeps building and keeps its legend block; it is not renumbered a third time until a sentence asks for it. | `gsva_scores.rds`, `ap7_mb_fork.rds` (assertion) |

### Fig. 1E, built 2026-08-04

The first panel of the written third paragraph, and the first of fourteen (see the slot map
below). Two parts because the sentence carries two claims and they are on different rulers.

**TOP — the absolute claim.** Five arms as a percentage of the *nuclear* transcriptome (the
denominator drops only the 13 mt-* genes), every animal drawn, four groups per arm.
`fig01_mito_content.R:118-170` is the idiom, ported to `theme_panel(base_size = 6)`.

**The bracket is the POOLED genotype effect, and that is the load-bearing choice.** The number
the sentence quotes is the genotype main effect from script 32's own model,
`lm(log2 share ~ timepoint + myc_status)` — reproduced here and asserted equal to
`$share_stats$geno_beta`/`geno_p` **to 1e-8**. The groups are therefore drawn **genotype-major**
(both WT boxes, then both Myc+ boxes, which is `names(group_cols)` order), so a single bracket
spanning the two halves *is* that contrast. The additive model is the right one: every arm's
genotype x timepoint interaction is far from significant (p 0.35–0.81), i.e. the content effect
does not differ between the ages. Note that `figS1_mb_fork_specificity.R` uses the **opposite**,
age-major order — deliberately, because its test is the genotype gap *within* an age. The order
follows the test, not the panel.

**BOTTOM — the systemic claim.** One point per MitoPathway over a kernel density, coloured by
sign, with zero and the median (+0.429) marked: **95.1 % of the 143 above zero, 7 below.** Values
are set-average **raw** log2 fold changes, as CLAUDE.md requires of an averaged-LFC visual. The
synthetic mtDNA-encoded pathway carries `is_mtdna` and is excluded exactly as script 40 excludes
it, so the panel and the regressions describe the same 143.

**Deliberately neutral.** The tier colour key belongs to Fig. 1F, which ranks these *same 143
pathways* on the priority ruler. The pair is the two-mechanism sentence: one-sided here,
two-sided there. Spending the tier key twice in one figure would obscure that.

**One new encoding, declared not inlined:** `sig_cols` in `_panel_common.R` (red `#E41A1C` for
p < 0.05, `grey45` otherwise) — the bracket-label ink `fig01_mito_content.R:106` had inline. It
applies to bracket labels only; significance is never a fill or a point colour in this project.

**Two tibble traps, both fixed and worth remembering.** Scripts 32 and 40 save **tibbles**, and a
tibble's `[` returns a tibble rather than a scalar, so `stats[a, "n_genes"]` silently poisons
every `sprintf` downstream (it fails with "unsupported type", which does not name the cause).
Every such object is `as.data.frame()`d on read. And `rownames(x) <- …` on a tibble is deprecated
and will eventually stop working.

Sizes: 89 x 78 mm, top:bottom heights 2.5:1, guides collected to one row at the foot.

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

**And each label and its leader take the colour of the set's own class** (third review, same day),
so the key doubles as a way to tell the labels apart. All seven sit below their own point except
the redox control, which sits above it: by rank 585 the curve has fallen past all three lanes, so
the empty space near that point is above the line. The `stopifnot` names that one exception rather
than being relaxed.

**Fig. 1C has a broken x axis.** One 6-week Myc+ animal sits at PC1 = +44.3 while the next highest
is +19.1, so an unbroken axis spends a third of the width on nothing. Built as two panels rather
than with a package (`ggbreak` is not installed, and doing it by hand keeps `coord_fixed(1)` alive
in both segments, which is the part that matters — within a segment a PC1 unit and a PC2 unit are
the same length, so the collapse onto PC1 is still the picture and not the panel shape). Segment
widths are proportional to their spans, and the break is **one back-slash per side, sitting exactly
on the end of each axis line** (third review): one closes the main segment, one opens the outlier
segment. **What it costs:** a
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
   `ap7_mb_fork.rds$fork_df$MB2_score`. **`figS1_mb_fork_specificity.R` carries the specificity
   test that has to go with it** — it held S1D until 2026-08-04 and now holds no slot, because the
   written S1D is the MYC western blot and this row is still uncited.
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

### Added 2026-08-04, from the written paragraphs 3 onward

Six checks against the text as written. None of them is a result changing; all six are wording.

- **"94% of 143" mixes two numbers.** 94.4 % is the fraction above zero over **all 144** rows of
  script 40's ruler (the synthetic mtDNA-encoded pathway included); **95.1 %** is the fraction
  over the **143** that every regression, null and panel actually uses. Fig. 1E draws the 143.
  Pick one and say which.
- **"21–27%" is narrower than what Fig. 1E draws.** Those two numbers are the median per-gene
  rise across the nuclear compartment (**+21.4 %**, 977 genes) and the mass markers (**+26.5 %**).
  The five *set* arms on the panel run **+20 % to +33 %** (nuclear MitoCarta +20.0, mass markers
  +26.5, nuclear OXPHOS +32.2, biogenesis +32.6, mtDNA +5.9 ns). Either widen the sentence to
  "20–33 % depending on the arm" or name the two quantities the range comes from.
- **"(padj < 0.001)" is `BIOGENESIS_FULL`, p = 0.00049**, and it is a **pooled** genotype main
  effect, not a per-age one. Split by age at n = 6 per cell, the four nuclear arms clear p < 0.05
  at **twelve** weeks (0.0015–0.016) and only biogenesis clears it at **six** (0.034; the others
  sit at 0.08–0.11). The effect *size* is the same at both ages — every interaction p is 0.35–0.81
  — so the pooled test is the right one, but the text should not imply a six-week significance
  that the six-week data alone do not carry.
- **"OXPHOS subunits (all complexes)" is not exact on the priority ruler.** Complex IV sits at
  **+0.017**, flat; CV (+0.165) and CIII (+0.160) carry the tier. And the **assembly factors** of
  those same complexes are at **−0.028** — the subunit-versus-assembly split is the interesting
  half and it reappears in the wild-type withdrawal. Write *subunits, led by complexes V and III*.
- **"MYC significantly promoted … biosynthetic pathways including glycine cleavage, pyruvate, and
  serine metabolism"** is right, and it is worth knowing those are the *largest* single
  promotions (+0.41 / +0.33 / +0.21) — larger than OXPHOS, which leads as a **tier** and not as a
  leaf. The sentence as written is fine; do not let a later draft turn it into "OXPHOS leads".
- **`Bbc3` is written "(p < 0.01)"**; the interaction p is **0.0081** and its genome-wide BH is
  **0.84**. What licenses the test is that `Bbc3` was **pre-specified from the cell experiments**,
  and the text does not yet say so. One clause, once.
- **Fig. 1G's "overall R2, OXPHOS R2, p" placeholders.** The numbers that exist: over the 143
  mitochondrial pathways, slope **0.552**, intercept 0, **R² = 0.801** (bootstrap 0.514–0.605),
  and that R² sits at the **100th percentile** of expression-matched shuffled sets. On the
  priority ruler, slope 0.644, R² 0.789. Within the OXPHOS tier, **R² 0.761 content / 0.945
  priority**. For "the whole transcriptome" the number with a home in a numbered script is script
  44's fitted global rate **0.487** over 8,774 genes
  (`collapse_module_ownership.rds$defs$global_rate_fitted`). Which pair the sentence quotes is
  open — settle it when 1G is built.

### One text number still to make exact

- **"~3K DEGs"** at 6W is **2777** at FDR 10 % and **1967** at FDR 5 %; the `interaction` contrast
  has **0** at either threshold.
- *(Resolved 2026-07-31: "~900 genesets" is now written as 986 in paragraph 2, which is the library
  count and correct. See the three-counts entry above.)*

**`figS1_mb_fork_specificity.R`** (built 2026-07-30 as S1C; rolled to S1D on 2026-07-31; displaced
out of the figure entirely on 2026-08-04, all three times because it is not cited) — the
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

## The slot map (2026-08-04, read off the written Results)

**This supersedes the "Planned — Figures 2 and 3" commission of
`docs/2026-08-02_results_narrative_final_order_D.md`.** That document is still the source for the
*argument* — the corrections, the seven-item dose list, the two answered questions, the trap list
— but its panel table is now wrong. What the author actually wrote:

- the oncogene and its de-dosing are **one** figure: Figure 1 grows from B/C/D to **B–H**;
- everything from the MYC-ER model onward is **Figure 2** (`2A–2I`);
- **Figure 3 does not exist** — "the whole narrative has to fit into two figures with 2
  supplementaries";
- **Fig. 1E is free** (the superseded partial-correlation paragraph is gone from the narrative)
  and is reused for mitochondrial content;
- Supplementary 1 gains **S1D–S1F**, Supplementary 2 is **S2A–S2D**.

`bench` = the author's own data, assembled outside R — not built here.

### Figure 1 — the oncogene, and its de-dosing

| slot | the sentence it supports | script | inputs |
|---|---|---|---|
| **1A** | the hyperplastic ductal expansion; **and now also** "the ratio of the Ki-67 and CC3 positive cells increased significantly at 12W … mainly due to a reduction of the number of cells undergoing apoptosis" | bench | — |
| **1B** | TEB→ductal transition; MYC amplifies TEB/proliferation and the LE de-differentiation | BUILT | — |
| **1C** | the dominant principal component axis | BUILT | — |
| **1D** | biogenesis and OXPHOS top-ranked by NES | BUILT | — |
| **1E** | "a quantitative expansion … 94% of 143 … 21–27% (padj < 0.001)" | BUILT 2026-08-04 | — |
| **1F** | "a resource reallocation altered intra-compartmental priorities, as quantified by MitoPPS … promoted protein import/homeostasis, translation, OXPHOS subunits … demoted dynamics and surveillance, including fission, mitophagy, and apoptosis alongside calcium signalling"; **and, from the next section**, "a similar effect on the mitopathway priority scores was observed (see Fig. 1F)" — so the 12W fade must read in the same panel | `fig1_reallocation_ranked.R` | `mitopps_scores.rds`, or `background_vs_myc.rds$ruler` (`p_m6`/`p_m12`) |
| **1G** | "the entire structure is present at half amplitude and unchanged in shape both in the whole and mitochondrial transcriptome (overall R2, OXPHOS R2)" | `fig1_rescaled_not_reshaped.R` | `background_vs_myc.rds$regressions`, `$regression_null`, `$ruler`; `collapse_module_ownership.rds$defs$global_rate_fitted` |
| **1H** | "this ranking was not just maintained but enhanced: the normalized enrichment score for the mitochondrial OXPHOS set and MYC integrative signatures increased substantially" | `fig1_nes_preserved_sharpened.R` | `fgsea_percategory.rds` |

### Supplementary 1

| slot | supports | script |
|---|---|---|
| **S1A** | the design | BUILT |
| **S1B** | the 986-set library | BUILT |
| **S1C** | what the dominant axis is made of | BUILT |
| **S1D** | "this attenuation was driven by a ~50% reduction in MYC protein levels at 12W compared to 6W" | bench (western) |
| **S1E** | "this decline occurred despite stable transcript levels; the expression gap remained constant between 6W_myc and 12W_myc for MYC …" | `figS1_myc_transcript_stable.R` |
| **S1F** | "… and the proximal MYC/MAX/MXD network" | `figS1_myc_network.R` (port `figures/figS8_myc_network_levels.R`) |

### Figure 2 — the window, the tissue, and the trigger

| slot | the sentence it supports | script | inputs |
|---|---|---|---|
| **2A** | the R26-LSL-CAG-MYCER<sup>T2</sup> x MMTV-Cre model | bench | — |
| **2B/2C** | 48 h tamoxifen: proliferation and apoptosis by IHC, 6W vs 12W | bench | — |
| **2D** | 13-day follow-up: MEC expansion and E-cadherin ductal expansion only at 12W | bench | — |
| **2E** | "the TEB signature was lost as expected in the puberty-adult transition, [while] the overall proliferation signalling remained relatively stable" | `fig2_wt_teb_proliferation.R` | `substrate_specificity_tradeoff.rds$wt_null`, `priming_arm_teb.rds$teb_signatures` |
| **2F** | "MECs withdraw from the respiratory chain; OXPHOS subunit LFC **and** MitoPPS drop across all complexes while biogenesis pathways remain relatively stable … the maturing gland upregulates amino-acid and lipid catabolism" | `fig2_wt_mito_contraction.R` | `substrate_specificity_tradeoff.rds$wt_null` + `$comparator_priority`, `background_vs_myc.rds$ruler` |
| **2G** | "Most pro-anti-apoptotic ratios established by MYC at 6W are reduced according to the global rescaling … the PUMA/BCL-XL ratio exhibited a significant reversal" | `fig2_priming_ratios.R` | `priming_arm_teb.rds$priming`, `$pair_null` |
| **2H** | "closely tracking *Bbc3* was its known p53-independent activator, *Foxo3*" | `fig2_departure_from_dose.R` | `collapse_module_ownership.rds$collapse_genes` |
| **2I** | "the interaction between MYC and OXPHOS subunit coupling to the PUMA/Bcl-XL ratio is highly significant (p = 0.0052), whereas no such link exists for other redox and metabolic axes" — **author's call 2026-08-04: this gets a panel, and it is 2I, not a supplementary** | `fig2_oxphos_puma_coupling.R` | `substrate_specificity_tradeoff.rds$tradeoff`, `priming_arm_teb.rds$axis_scores` + `$coincidence_fit` |

### Supplementary 2

| slot | supports | script |
|---|---|---|
| **S2A** | "these alterations were not a consequence of decreased MYC expression at 12W" | bench |
| **S2B** | "the WT temporal program operated independently of the MYC-driven reallocation, showing a negligible correlation across the mitochondrial and the whole transcriptome" | `figS2_reallocation_independence.R` |
| **S2C** | "the overall apoptotic priming remains stable in the WT timeline" | `figS2_priming_balance.R` |
| **S2D** | "no changes in the transcriptome of other known PUMA inducers were observed" — the roster already exists as `priming_arm_teb.rds$exclusions$puma_inputs` (12 genes: `E2f1`, `Trp73`, `Atf4`, `Ddit3` …) | `figS2_puma_inducers.R` |

**Fourteen new panel scripts, built one at a time in citation order:**
`1E → 1F → 1G → 1H → S1E → S1F → 2E → 2F → S2B → S2C → 2G → 2H → S2D → 2I`.
Most are ports of `figures/fig01`–`fig05` / `figS8` into this directory's idiom
(`theme_panel()`, a `panel_legend()` block, `save_panel_p()`, the declared palette and contrast
vocabulary, no prose on the page); **1H, 2H, S2B, S2C and S2D have no port and are new builds.**

**Deliberately not panels.** The DE counts (1967 → 135 at FDR 5 %, 2777 → 239 at 10 %, median
|LFC| of the survivors *rising* 0.678 → 0.869) and the interaction standard-error point (0.333 vs
0.233). Both are threshold arithmetic; drawn, they would flatter a −93 % that is really a −45 %.
The author has them as a trailing italic paragraph in the Doc, which is the right place.

**Three numbers in the narrative are new to this repo** and were derived read-only on 2026-08-02
(no script re-run): `cor(Myc@6W priority, WT-time priority)` = **−0.02** across 143 pathways, so
the background reprioritisation is *independent* of Myc's, not its mirror; median `lfcSE` **0.333**
(interaction) vs **0.233** (genotype); and `Foxo3` at the **0.456th percentile** of the
departure-from-dose scan, adjacent to `Bbc3` at the 0.513th, holding that position within the 4199
Myc-induced genes alone. The first is now cited in the text (S2B) and the third is a panel (2H),
so **both need a home in a numbered script**, not only in the figure layer.

**One freshness trap.** `results/mitopps_scores.rds` (12:35) is one minute *older* than
`results/gsva_scores.rds` (12:36) although both came out of the same post-reconciliation re-run,
so `require_fresher_than()` would wrongly stop **Fig. 1F**. Gate that panel on
`interaction_results.rds` instead, and say why in the script.

## Not built here

| slot | why |
|---|---|
| **Fig. 1A** | Whole-mount / histology of the hyperplastic ductal expansion — the author's bench image, assembled outside R. **It now carries a second quantity** (written 2026-08-04): the Ki-67 : CC3 ratio, up at 12W "mainly due to a reduction of the number of cells undergoing apoptosis". That ratio is the phenotype the whole of Figure 2 explains, so it is worth its own sub-panel rather than a sentence. |
| **Fig. S1D, S2A, Fig. 2A–2D** | The author's bench data: the MYC western blot; MYC-ER expression at 12W; the MYCER<sup>T2</sup> model schematic; and the tamoxifen IHC and 13-day expansion series. |
| The old Fig. 1E paragraph | **SUPERSEDED as narrative, author 2026-08-02.** The partial-correlation paragraph (proliferation, de-differentiation, apoptotic priming, the human MYC-driven signature) is superseded by scripts 35/36: the raw ambient is the wrong null, "OXPHOS is central" is untestable rather than false, and `cholesterol → fork` was demoted. It stays in the Doc **only as a source of reusable statements**. The *slot* 1E was freed by that and is now the mitochondrial content panel — do not confuse the two. |
| Gene-level vs pathway-level PCA | Author's call 2026-07-31: pathway lens only, to keep the figure from turning into a methods argument. The comparison (gene PC1 = 44 %, a composition axis; the Myc programme rides at gene-PC2) is in Fig. 1C's legend block, and script 37's figures A2/A3/A4 hold the drawn version if a supplementary is ever wanted. |
| Scree / dimensionality | Same reason. The numbers that matter (PC1 76.5 %, effective dimensionality 1.7, 79 % one-signed loadings) are in Fig. 1C's legend block. |
| MEC state *levels* | The four groups as a 2×2 grid per programme, group-mean z. Built, reviewed, dropped as redundant with the contrasts; kept in `fig1_myc_teb_proliferation.R`'s sandbox. |
| Gray HE arm | Computed and asserted, reported in the legend block, not drawn — see point 4 above. |
| sample PCA | Deferred 2026-07-29 — belongs to the global-axis section (alignment along the mitochondrial axis), not to the opening overview. |
| DE counts per contrast | Deferred 2026-07-29 — belongs to the attenuation section. Numbers already checked: 6W genotype contrast = **2777** genes at FDR 10 %, **1967** at FDR 5 %; the `interaction` contrast has **0** at either threshold. |

## Assembly

**The numbering is now fixed** (2026-08-04, the slot map above), so assembly is no longer blocked
on the narrative — it is blocked on the panels. Once they are signed off, a
`figures/figure*_overview.R` composes them with patchwork `tag_levels = "A"` and this table
records the mapping.

Every panel is **single column (89 mm)** wide except the displaced MB-fork one (50 mm). Heights as
built: 1B 84, 1C 54, 1D 68, **1E 78**, S1A 52, S1B 52, S1C 56, (MB fork 52) mm. Figure 1 is now
eight slots and will be a full page: A (bench) plus B–H, of which four are built and stack in
284 mm, so **1F, 1G and 1H together have about 180 mm of column to live in** — worth knowing
before designing them, because the two-panel composite form used for 1E is expensive in height and
1F in particular wants width, not height.

The four already-assembled manuscript figures (`figures/figure1_myc_mitochondrion.R`,
`figure2_developmental_window.R`, `figureS1_compartment_detail.R`, `figureS2_controls.R`) are
untouched by this layer and still rebuild via `figures/rebuild_manuscript_figures.R` — but they now
pick up the new sample palette, so they should be re-rendered and re-checked at the same time.
