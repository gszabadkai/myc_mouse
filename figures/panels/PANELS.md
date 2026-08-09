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
| **Fig. 1F** | `fig1_reallocation_ranked.R` | p3: "a resource reallocation altered intra-compartmental priorities, as quantified by MitoPPS … promoted protein import/homeostasis, translation, OXPHOS subunits … demoted dynamics and surveillance, including fission, mitophagy, and apoptosis alongside calcium signalling"; **and** p4: "a similar effect on the mitopathway priority scores was observed (see Fig. 1F)" | `background_vs_myc.rds` (`$ruler`, `$ruler_summary`), `mitopps_scores.rds` (assertion only) |
| **Fig. 1G** | `fig1_rescaled_not_reshaped.R` | p4: "the entire structure is present at half amplitude and unchanged in shape both in the whole and mitochondrial transcriptome (overall R2, OXPHOS R2, p)" | `collapse_module_ownership.rds` (`$collapse_genes`, `$defs`), `background_vs_myc.rds` (`$ruler`, `$regressions`, `$regression_boot`, `$regression_null`, `$null_draws`) |
| **Fig. 1H** | `fig1_nes_preserved_sharpened.R` | p4: "this ranking was not just maintained but enhanced: the normalized enrichment score for the mitochondrial OXPHOS set and MYC integrative signatures increased substantially" | `fgsea_percategory.rds` |
| **Fig. 2E** | `fig2_wt_teb_proliferation.R` | s2p1: "The 6\>12W_wt comparison showed that while the TEB signature was lost as expected in the puberty-adult transition, the overall proliferation signalling remained relatively stable" | `fgsea_percategory.rds` (`$fgsea`, the `timepoint_neg` ranking), `substrate_specificity_tradeoff.rds` (`$comparator`, `$wt_null`, `$defs`), `priming_arm_teb.rds` (`$teb_signatures`) |
| **Fig. S2B** | `figS2_reallocation_independence.R` | s2p1: "the WT temporal program operated independently of the MYC-driven reallocation, showing a negligible correlation across the mitochondrial and the whole transcriptome" | `background_vs_myc.rds` (`$ruler`, `$geometry`), `interaction_results.rds`, `gsva_scores.rds` (`$expr_mat`, the split-baseline control) |
| **Fig. 2F** | `fig2_wt_mito_contraction.R` | s2p1: "MECs withdraw from the respiratory chain; OXPHOS subunit LFC **and** MitoPPS drop across all complexes while biogenesis pathways remain relatively stable … the maturing gland upregulates amino-acid and lipid catabolism" | `background_vs_myc.rds$ruler` (`c_tn`, `p_tn`), `substrate_specificity_tradeoff.rds` (`$comparator`, `$comparator_priority`, `$wt_null`, `$defs`) |
| **Fig. 2I** | `fig2_oxphos_puma_coupling.R` | s2p2: "the interaction between MYC and OXPHOS subunit coupling to the PUMA/Bcl-XL ratio is highly significant (p = 0.0052), whereas no such link exists for other redox and metabolic axes" | `substrate_specificity_tradeoff.rds` (`$tradeoff`, `$tradeoff_perm`, `$ambient`), `priming_arm_teb.rds` (`$axis_scores`, `$purity`), `gsva_scores.rds` |
| **Fig. S2D** | `figS2_puma_inducers.R` | s2p2: "no changes in the transcriptome of other known PUMA inducers were observed" | `priming_arm_teb.rds$exclusions$puma_inputs` |
| **Fig. 2H** | `fig2_departure_from_dose.R` | s2p2: "closely tracking *Bbc3* was its known p53-independent activator, *Foxo3*" | `collapse_module_ownership.rds` (`$collapse_genes`, `$defs`), `interaction_results.rds`, `combined_df_annotated_raw.rds` |
| **Fig. 2G** | `fig2_priming_ratios.R` | s2p2: "most apoptotic priming ratios … remained stable, since both pro- and anti-apoptotic proteins diminished in accordance with the global rescaling … the PUMA/BCL-XL ratio showed a striking reversal" | `priming_arm_teb.rds` (`$priming`, `$pair_null`), `collapse_module_ownership.rds` (`$defs$global_rate_fitted`, `$wt_genes`) |
| **Fig. S2C** | `figS2_priming_balance.R` | s2p1: "the overall apoptotic priming remains stable in the WT timeline" | `collapse_module_ownership.rds$wt_genes`, `substrate_specificity_tradeoff.rds$buffer`, `background_vs_myc.rds$ruler`, `interaction_results.rds` + `combined_df_annotated_raw.rds` (power control) |
| **Fig. S1A** | `figS1_design_contrasts.R` | p1: "both longitudinal (6W versus 12W for WT or Myc+ genotypes) and cross-sectional (WT versus Myc+ at 6W or 12W) comparisons" | none (schematic) |
| **Fig. S1B** | `figS1_geneset_library.R` | p1: "fGSEA ranking and GSVA scoring based on custom-curated genesets"; p2: "we extended the custom library to a total of 986 genesets" | `provenance_table.csv`, `pathway_loading.rds` (assertion), `gsva_scores.rds` (count) |
| **Fig. S1C** | `figS1_axis_loadings.R` | p2: "… aligned almost entirely with variability in mitochondria related terms" — what the axis is made of, with its bound | `gsva_scores.rds`, `pathway_loading.rds` (`mito_classification`, `mito_enrichment`) |
| **Fig. S1E** | `figS1_myc_transcript_stable.R` | p5: "this decline occurred despite stable transcript levels; the expression gap remained constant between 6W_myc and 12W_myc for MYC" | `dds_int_run.rds`, `interaction_results.rds`, `combined_df_annotated.rds`, `collapse_module_ownership.rds` (legend) |
| **Fig. S1F** | `figS1_myc_network.R` | p5: "… and the proximal MYC/MAX/MXD network" — the negative control on the alternative to dose | `interaction_results.rds`, `combined_df_annotated.rds` |
| *(no slot)* | `figS1_mb_fork_specificity.R` | **displaced 2026-08-04.** It held S1D, and the written S1D is the MYC western blot. Still cited nowhere: it is the specificity control for Fig. 1B's `Human BRCA-MYC` row. It keeps building and keeps its legend block, and its `panel_legend(slot = )` now reads **`not currently cited`** rather than a letter, so anything that lays panels out by slot skips it. It gets a letter back the day a sentence asks for one. | `gsva_scores.rds`, `ap7_mb_fork.rds` (assertion) |

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

### Fig. 1F, built 2026-08-04

Cited **twice** — paragraph 3 for the reallocation itself, paragraph 4 for "a similar effect on
the mitopathway priority scores was observed (see Fig. 1F)" — so the 6W→12W fade has to read in
the same panel. Two parts again, and the second is the pair to Fig. 1E's.

**TOP — the thirteen the sentence names**, ranked, each with its 12W effect joined to its 6W one
by a connector. Colour is the declared contrast vocabulary (`contrast_cols`: `myc_6W`,
`myc_12W`); **shape is significance** — filled at padj < 0.05 at six weeks, open otherwise —
because the significance is heterogeneous in a way that matters to the sentence (below).

**Why named rows and not all 144 ranked.** The sentence names twelve things. A 144-point ranked
chart with twelve leader labels is unreadable at 89 mm, and unlabelled it cannot support the
sentence at all. `figures/fig02_reallocation_ranked.R` is the all-144 form and needed 183 mm to
label seven tiers.

**BOTTOM — the compartment the thirteen came out of**, which is what answers the cherry-picking
objection: all 143 non-mtDNA pathways at 6W, the *same construction* as Fig. 1E's lower strip on
the other ruler. **That pairing is the point.** On content, 95.1 % sit above zero; on priority,
56.6 %. Two panels, two strips, and the difference between them *is* "MYC-driven mitochondrial
alterations manifested through two mechanisms".

**One x scale for both parts**, fixed over everything drawn, with the top part's tick labels
dropped so the axis at the foot serves both. A vertical dropped from any row lands where that
pathway sits in the whole compartment; separate scales would have made that reading wrong.

**n is printed in every row label.** Several of the largest effects rest on four genes (glycine
cleavage 4, serine 4), and MitoCarta sets are membership-loose — a big effect on a four-gene set
is one gene, not a module. That belongs on the panel, not in the legend.

**Provenance asserted, freshness guard deliberately NOT called.** `$ruler$p_m6` is asserted
identical (to 0, not to a tolerance) to script 08's `mitopps_pairwise` diff and padj. But
`results/mitopps_scores.rds` is timestamped one minute *before* `results/gsva_scores.rds` although
both came out of the same 2026-07-24 re-run, so `require_fresher_than()` would stop this panel for
an artefact of the ordering inside that session. The identity check is the stronger guard.

**`tier_cols` was re-keyed in `_panel_common.R`** while building this: its keys were abbreviated
strings that matched nothing on disk, so any lookup would have returned `NA` silently. They are
now the seven values `$ruler$tier` carries, with a separate `tier_labels` for display. **1F does
not use it** — four of those hues are the sample palette, and on a page where panel E spends blue
on wild type and orange on Myc+, a tier key reusing them invites the reader to see genotype in a
tier. Thirteen named rows need no tier colour.

Sizes: 89 x 72 mm, top:bottom heights 3:1.

### Fig. 1G, built 2026-08-04 — and it fills in the `R2=***` placeholders

**Two facets, and that is all** (author's review, same day). The first version also carried the
numbers on the right-hand facet and a null distribution underneath. Both came off: *"too much
statistics detail — slope and R-squared is OK for the whole genome; on the mitochondrial plot it
is enough to label the OXPHOS points with the different colour, there is no need to write out the
numbers, can go in the text. Also, the bottom panel is not necessary, it just shows the
significance level, which can again go in the text."* Everything removed is still **computed and
asserted** in the script and reported in the legend block, and the drawn null is kept in the
sandbox — so the text has its numbers and the page does not spend a strip on one p-value.

Three more cuts in the same review: the **OXPHOS fit line** came off (a second line at slope
0.525 beside one at 0.552 is two lines saying they are the same, which the points say better),
the **overall mitochondrial slope and R2 went on** to the right-hand facet so each facet carries
its own line's numbers, and the **OXPHOS key moved inside the frame** — bottom right of the right
facet, which is the wedge the diagonal cloud leaves empty, since a large six-week effect with a
small twelve-week one is exactly what does not happen. Under the plot it cost a whole row of
height for one dot.

The `R2` is a **real superscript**, via plotmath (`parse = TRUE`), which keeps the source ASCII
per the coding rules and renders correctly through the base `pdf()` device where a literal
character would not. **The number must be quoted inside the plotmath string** (`R^2~"0.80"`): left
unquoted it is parsed as a number and printed as `0.8`, silently dropping the significant figure
the reader is being given.

**The distinction the sentence turns on, now text-only:** "half amplitude" is the **slope**;
"unchanged in shape" is the **R2**. They are independent, and only one of them is beyond its null:

| | observed | null median | null max | draws >= observed | p |
|---|---|---|---|---|---|
| slope | 0.552 | 0.419 | — | **64 of 500** | 0.13 |
| **R2** | **0.801** | 0.298 | **0.719** | **0 of 500** | **< 0.002** |

A shuffled set of the same size and expression can halve an effect. What it cannot do is halve it
*coherently*. **So the p that belongs in the parenthesis is the R2 one, and the sentence should
not attach a p to the slope.**

**The numbers for the placeholders:**

| in the sentence | number | what it is |
|---|---|---|
| whole transcriptome, slope | **0.487** | through the origin, over the 2,648 genes Myc moves at 6W (padj < 0.1, \|LFC\| ≥ 0.2, baseMean ≥ 20) — script 44's `global_rate_fitted`, asserted |
| whole transcriptome, "overall R2" | **0.51** | over the same 2,648. Over **all 8,774** reported genes the slope barely moves (0.450) but R2 falls to **0.381** — regression dilution from six thousand genes Myc does not move |
| mitochondrial, slope | **0.552** | intercept 0.0001, so a pure rescaling and not a shift; bootstrap 0.514–0.605 |
| mitochondrial, "overall R2" | **0.801** | |
| **"OXPHOS R2"** | **0.761** | the 19 OXPHOS pathways, content ruler, slope 0.525. **Not 0.945** — that is the *priority* ruler, and mixing rulers inside one parenthesis would be wrong. The priority figures (slope 0.644, R2 0.789; OXPHOS 0.945) belong with the "see Fig. 1F" clause |
| "p <" | **0.002** | empirical, 0 of 500 expression-matched shuffles reach R2 0.801 |

**"Whole transcriptome" needs one word.** The panel's left facet is labelled *Myc-responsive genes
(2,648)*, because that is the set the rate of record is fitted on. If the sentence keeps "whole
transcriptome" it should quote 0.450 / 0.381; if it quotes 0.487 / 0.51 it should say "the genes
Myc moves". Either is fine; they must not be crossed.

**R2 is the squared Pearson correlation throughout**, so the two facets use one definition. A
through-origin `lm` reports an *uncentred* R2 instead (0.565 rather than 0.513 on the left), which
is not comparable to the pathway facet and must never be the quoted number.

**Both facets are windowed for display** — ±3.0 log2 on the left (22 of 2,648 genes outside) and
+1.05 on the right (2 of 143: Glycine metabolism +1.14 and the four-gene Glycine cleavage system
+1.95). Every fitted line and every number is computed on the complete set; the counts are in the
legend block, as Fig. 1E states its own capped point.

**One trap found here, worth remembering:** `geom_segment(x = , y = )` with no `data =` inherits
the layer data and draws the marker once per row — with a 500-row null that turned the "observed"
label into a black smear. Use `annotate()` for single marks.

Size: 89 x 46 mm. OXPHOS highlight is `ms_diverging[["pos"]]` mint — declared, and not a sample
colour.

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

### Fig. 1H, built 2026-08-04 — and it does not support the sentence as written

866 gene sets, NES at twelve weeks against NES at six, identity line, and the two programmes the
sentence names picked out by `programme_group()` — the same encoding Fig. 1D groups its rows by,
so "OXPHOS" and "Myc signatures" mean the same thing on both panels.

**"Maintained" is strongly true and is the panel's real result.** Spearman **0.933** over 866
sets while the effect size halves (Fig. 1G). Better still: the 92 sets that cross sign between the
ages — the two off-diagonal clouds — are **84 % mammary-development and TF-target sets**. Set them
aside and Spearman over the remaining 774 rises to **0.967**. *What holds its order is the
metabolic, mitochondrial and Myc core; what reshuffles is lineage identity* — which is section 2's
subject, arriving early.

**"Enhanced" needs the sentence changed.** Two independent reasons, and the panel shows the first:

1. **The rise is not specific.** Every point is above the line, not just the named ones —
   **81 % of all 866 sets rise (median +0.312)**, and **93 % of the 403 already enriched at six
   weeks (median +0.319)**. As families the two named rise by *the median amount*: OXPHOS
   **+0.308, the 49th percentile**; Myc signatures **+0.386, the 60th**. The individual sets the
   Order D doc quoted are the best members of their families, not typical of them —
   `MITOCARTA_OXPHOS` +0.52 (76th pct), `MYC_felsher_integrative_signature` +0.55 (78th),
   `HALLMARK_MYC_TARGETS_V1` +0.75 (87th).
2. **A global NES rise is what a weaker ranking mechanically produces.** fGSEA normalises each
   enrichment score against a permutation null built from *that same ranked list*. The twelve-week
   Wald list is much flatter than the six-week one — **SD 1.221 vs 1.744, IQR 1.503 vs 2.176,
   9.5 % of genes past |stat| = 2 vs 22.5 %** — so random sets reach smaller enrichment scores,
   the normaliser shrinks, and the same relative enrichment scores a higher NES. **No claim of
   increased Myc activity can rest on this panel.**

**What is true, and is a better sentence:** NES is enrichment *relative to the rest of the
transcriptome*, so these programmes **hold their position at the top of the ranking while the
amplitude halves** — the residual Myc signal is more concentrated on them than it was. That is the
retention asymmetry of scripts 29–31 (Myc core retains 0.74–0.77 against the mitochondrial arms'
0.50–0.58) read on a different instrument. It is a statement about **position, not magnitude**,
and it sits naturally beside Fig. 1G rather than against it.

**One encoding note:** `ms_diverging` is a *ramp* everywhere else in the manuscript; this is the
one panel that spends its two poles as categorical colours (mint = OXPHOS, matching Fig. 1G;
espresso = Myc signatures). Neither is a sample colour. Flagged in the legend block, not left to
be noticed.

Also visible and unplaced: the upper-left cloud is the **release of Myc's suppression of luminal
hormone-sensing identity** (`MG_LHS_CONSENSUS` −3.20 → +2.05, `MG_LHOR_SAEKI` −3.45 → +2.00) —
the largest single departure from the dose line in the transcriptome, and still a beat the Results
has not placed.

Size: 89 x 62 mm; key inside, bottom right.

### Fig. S1E, built 2026-08-04

The control the western blot needs. If the message fell too, the attenuation would need no
further explanation and the dose argument would collapse into "the transgene was silenced". It
does not fall: the genotype gap is **+1.68 → +1.82 log2** and the Myc+ gland does not decline
across the window (**−0.10, padj 0.76**). Measured as retention, `Myc` is the **most retained gene
in the transcriptome** — 1.085, the 99.99th percentile of 8,774.

**Age-major group order, the opposite of Fig. 1E, and for the same reason: the order follows the
test.** 1E draws one pooled genotype main effect, so both wild-type boxes sit together and one
bracket spans it. Here the tests are the two *within-age* gaps, so each age's pair must be
adjacent. Neither order is a default.

**The y axis is log2 of normalised counts, labelled in counts** — on that scale "the gap remained
constant" is a constant vertical distance and can be read off the panel. On a linear count axis it
could not be.

**The brackets carry DESeq2, not a t-test on the drawn points.** Effects and adjusted p-values come
from the raw unshrunken interaction model, and the script asserts the empirical log2 difference
between drawn group means reproduces each DESeq2 fold change to **within 0.1** — so the picture and
the statistic are the same quantity. Neither this panel nor `dds_int_run.rds` is gated by
`require_fresher_than()`: the Step 1 objects predate the gene-symbol reconciliation because that
was a set-membership fix which never touched counts or per-gene results.

**New shared geometry in `_panel_common.R`:** `bracket_frame()` and `bracket_layers()`, the bar +
two end ticks + label idiom `fig01_mito_content.R:69-103` established. Additive on whatever scale
is drawn, so it is correct on a log axis. **The caller owns the x positions** — they follow the
drawn order, which follows the test, and there is no default right for both panels. *Fig. 1E still
has its own inline multiplicative copy from before the helper existed; it should adopt this the
next time it is touched.*

Size: 55 x 55 mm.

### Fig. S1F, built 2026-08-04

The second half of S1E's sentence, and a **negative control on the most obvious alternative to
dose**: MYC cannot bind an E-box without MAX and competes for the same sites with the MXD/MNT
family, so a fall in MAX or a rise in the repressors would attenuate its output at constant MYC.
Neither happens.

**It draws the CONTRAST, not the level** — the sentence's claim is about the gap, so the gap is on
the axis: each gene's Myc genotype effect at six weeks and at twelve, joined, with the roles as
facet strips. The levels form needs ten facets and 183 mm to say the same thing less directly
(`figS8_myc_network_levels.R` panel B).

**Nothing clears significance at either age** — closest are `Mxi1` (+0.33, padj 0.096) and `Mlx`
(+0.39, padj 0.156) — and no interaction is below **padj 0.83**. So there is no significance
encoding on the panel: a key whose filled entry never appears is worse than saying it once in the
legend.

**The error bar is on the six-week point only, and it is drawn OVER the connector in the six-week
colour.** Drawn underneath and in grey it read as spanning both points, which makes a twelve-week
point falling outside it look like a significant difference between the ages. One bar is enough
because **the two ages are equally precise: standard errors agree to within 0.024 log2 across all
ten genes, median ratio 1.01** — asserted in the script. That equality is worth having on its own:
*twelve weeks is not the noisier cohort*, which is the objection every attenuation result in the
paper has to survive.

**Two genes cannot carry a negative** and the bars show it rather than leaving it to be assumed:
`Mxd3` (baseMean 14, SE 0.52) and `Mlxipl` (baseMean 45, SE 0.33) are wide enough that a real
half-log2 effect would be missed.

The key is **under** the plot here, unlike 1G and 1H — this panel has no empty corner, the MLX rows
run to the right edge.

Size: 89 x 55 mm.

### Fig. 2E, built 2026-08-04 — the first panel of section 2

The substrate panel. Everything after it depends on the wild-type gland having changed
*developmentally* between six and twelve weeks, because the MYC-ER animals never saw Myc. It is
also the negative limb that keeps Fig. 2F honest: had proliferation collapsed here, the
respiratory withdrawal would be its shadow rather than a finding.

**The ruler is the fGSEA NES on `timepoint_neg`** (author's call, 2026-08-04) — the same
instrument as Figs. 1D and 1H, so "enrichment" means one thing across the figure set. **One
ranked list, so the normaliser is shared**, which is what makes this comparison safe in exactly
the way Fig. 1H's *cross-contrast* comparison was not. No statement about the size of an NES
crosses a contrast anywhere on this panel.

**Two rows, and nothing else** — author's review the same day: *"this looks quite messy.
Certainly there is no need to separately label the TEBvsDuctal(HS) and the median bars. TEB-down
can be omitted. I'm not sure that the TEBxMitoCarta adds to the interpretation."* The first
version drew four rows, a median rule per row and a label on the strongest set. What is left is
the sentence and only the sentence:

| row | n | median NES | padj < 0.05 |
|---|---|---|---|
| TEB programmes | 38 | **−1.58** | **29** — and every one of the 38 is depleted |
| proliferation | 14 | **−0.89** | **0** |

**The median rules came off because the ±1 guides do their job better.** fGSEA divides the
enrichment score by the mean of the same-signed permutation null, so |NES| ≈ 1 is where a set
that has not moved lands. Against that reference the proliferation row needs no summary statistic
— it *sits on the band* — and the TEB row's real statement is not its median but that **all 38 of
its members are below −1.1**, which the cloud shows directly.

**Two families came off the page and stayed in the legend block, computed and asserted.** The
three `_VS_DUCTAL_*_DN` sets (the direction control: they rise where the TEB sets fall — basal
+1.45 against its UP set's −2.34, hormone-sensing +1.06 against −2.76; the alveolar pair shows
neither) and the twelve `*_TEB_MITO` construction lanes (median +0.49, none significant; 17–59
genes each, MitoCarta subsets by build). **Omitting them from the page is not the same as folding
them into the TEB row**, which is what the script's structure enforces: the `_DN` sets have their
sign inverted by construction and would cancel part of the effect, and the `_MITO` lanes are
Fig. 1D's by-construction class. The four families are identified in one place and two are drawn.

**Three rulers agree on the TEB arm and two are independent of composition.** Gene-level
set-average raw log2FC **−0.4220** over 148 genes, the **0th percentile** of 2000
expression-matched random sets; purity-adjusted per-sample composite (`score ~ tp * genotype +
epithelial + immune`) **−0.4200, p 0.0053**. They agree to **0.002**, so residual stromal or
immune contamination does not explain the loss. Asserted in the script, reported in the legend
block, not drawn.

**`programme_group()` is deliberately not the grouping here** — it would scatter the TEB lanes
across five of its rows, and the sentence is about the TEB *context*, which cuts across them. The
proliferation roster is not re-typed either: it is `substrate_specificity_tradeoff.rds
$defs$prolif_sets`, asserted `setequal`, so the panel and the arm-level statistic quoted in the
legend describe one set of genes.

Size: 89 x 30 mm — the smallest panel in the set, which is what dropping two rows, the median
rules and the label bought. Shape is significance (Fig. 1F's encoding); one ink; key inside, top
right — every TEB programme is depleted, so the top row ends well left of zero and that corner is
empty.

### Fig. 2F, built 2026-08-04 — the two rulers against each other

**The form is `figures/fig04_substrate_specificity.R` panel C** (author's review, same day: *"I
think 2F worked better on fig04 panel C, showing the PPS and logFC plotted against each other and
colour coded up and down, according to the continuous colour scale"*). The first version was
eleven named pathways ranked in two side-by-side columns; this one puts **all 143** on the plane
and lets the named ones sit in it.

- **x = content**, set-average raw log2FC — what the compartment *has*.
- **y = priority**, mitoPPS pairwise-ratio, content-blind — what it *spends its budget on*.
- **fill = the manuscript diverging ramp on x**, espresso down / white at zero / mint up. It
  repeats the x axis rather than adding a variable — Fig. 1B's "encoded twice" idiom — so **no
  colour bar is drawn**: the x axis is the key, and the asymmetric arm scaling is inspectable
  there.

**The plane says three things two ranked lists could not.** The cloud runs along a **diagonal**
(Spearman 0.82, Pearson 0.79 over the 143), so the two rulers agree and the lower-left quadrant is
the conjunction the claim needs — a content fall alone could be normalisation, a content-blind
ratio fall alone could be a reshuffle inside a growing compartment. The cloud's **centre sits up
and right of the origin** (median +0.041 content, 72 % above zero), so the respiratory arm falls
*against* a rising compartment rather than with it — the clause that licenses "withdraw", now
visible instead of asserted. And the **OXPHOS assembly factors sit at the origin, a few
millimetres from their own subunits**: that distance is the entire specificity claim.

| pathway | n | content | priority | compartment rank |
|---|---|---|---|---|
| CIV subunits | 17 | **−0.411** | **−0.180** | 0th / 1st |
| CI subunits | 36 | −0.253 | −0.155 | 3rd / 2nd |
| CIII subunits | 9 | −0.251 | −0.126 | 3rd / 4th |
| CV subunits | 19 | −0.216 | −0.069 | 4th / 15th |
| mitoribosome | 83 | −0.024 | −0.041 | 24th / 30th |
| **OXPHOS assembly factors** | 66 | **+0.001** | **−0.002** | 29th / 47th |
| central dogma | 230 | +0.036 | −0.011 | 47th / 42nd |
| **CII subunits** | 4 | **+0.088** | **+0.026** | 69th / 66th |
| lipid metabolism | 111 | +0.122 | +0.063 | 76th / 79th |
| amino acid metabolism | 79 | +0.184 | +0.082 | 88th / 83rd |
| fatty acid oxidation | 39 | +0.227 | +0.144 | 93rd / 94th |

**Significance — the author asked, and here is everything that exists.** Nothing is drawn, and the
reason is on the record:

1. **Per pathway, priority ruler.** Script 08 *does* test this contrast (`Temporal_Myc-` in
   `mitopps_scores.rds$mitopps_pairwise`, and its `diff` is asserted **identical** to the ruler's
   `p_tn`). **0 of 143 clear BH < 0.05**; the smallest adjusted p is 0.12. What it does show is
   **coherence**: 23 pathways clear an unadjusted 0.05 against 7.2 expected, **binomial
   p = 8e−7** — the compartment moves together and no single pathway carries it.
2. **Per pathway, content ruler.** None exists in the saved objects. A one-sample t over
   member-gene fold changes could be computed but ignores inter-gene correlation and would be
   anti-conservative for exactly the coherently-regulated sets this panel is about. Deliberately
   not done.
3. **Arm level, expression-matched null (script 43, 2000 draws).** This is where it is decisive:
   pooled OXPHOS subunits **percentile 0.0** (p < 0.0005) against their own assembly factors at
   **50.2**; mitoribosome 26.4; amino-acid **100.0** and lipid **99.9** upward.
4. **THE NUMBER TO QUOTE IN THE TEXT, and it was not in the panel brief:** script 43's **paired
   null** redraws *both* sets of a contrast together and tests the **difference** — which is what
   "the chain falls *while* biogenesis does not" actually claims. **OXPHOS subunits minus
   mitoribosome = −0.230, percentile 0.05, p = 0.0005**; minus nucleotide metabolism −0.255,
   p = 0.0085; minus the pooled proliferation set −0.208, p < 0.0005.

A glyph would have to mark 8 tested pathways among 135 untested ones, which reads as 135
negatives — so the panel carries the pattern and the legend carries the tests.

**"Across all complexes" has one exception, and it needs a word in the text.** **CII subunits
+0.088 / +0.026** — the only respiratory complex with no mtDNA-encoded subunit, outside the proton
circuit, and also a TCA enzyme. It is the only one that does not fall. But **n = 4**, and MitoCarta
sets are membership-loose, so it is a direction to note and not a mechanism to claim. Suggested:
*"across the four complexes that carry mtDNA-encoded subunits"*, or simply *"CI, CIII, CIV and
CV"*.

**The mtDNA-encoded subunits are excluded** by the same `is_mtdna` rule Figs. 1E and 1F use — and
they are the largest movement in the compartment (**+0.633 content, +0.525 priority**). That is not
tidiness: the mtDNA-encoded read fraction is confounded three ways (real content, the proliferation
denominator, dissociation leak) and it is the one quantity that is **time-associated rather than
genotype-associated**, so a temporal contrast is exactly where it cannot be read. A nuclear-down /
mtDNA-up mitonuclear discordance is what it looks like; on this axis it is not adjudicable. Stated
in the legend block so the omission is visible.

**Windowed for display, and what falls outside is the membership-loose caveat made visible.** Two
pathways sit below the drawn priority range (Vitamin D metabolism −0.403, 4 genes; Selenoproteins
−0.318, 5) and are named in the legend; every number is computed on the complete 143. **Every
pathway at the periphery of this cloud is a set of three to six genes** — GABA metabolism 6,
catechol 3, molybdenum cofactor 5, the carnitine pair 5 and 6 — because a set mean's noise scales
as 1/√n. The core of the cloud is where the large sets are.

**Every label is placed and none is repelled — and the placement is arithmetic.** Two reasons, one
of them the author's (review of 2026-08-04): **the leader must leave the label on the label's own
line**, and ggrepel cannot promise that — it draws from wherever the box edge happens to fall, so
a repelled label and a placed one do not look like the same object. And seven of the eleven sit in
the knot at the origin or on the crowded upper diagonal, where repel had nowhere local to put them
and stacked them on each other (it did, twice). **ggrepel is no longer a dependency of this panel.**

Three lanes, each in a region the script proves is empty:

| lane | region | holds |
|---|---|---|
| upper left | nothing gains priority while losing content | the knot's four |
| bottom band | nothing sits below −0.135 right of −0.10 | the four complexes |
| right band | past +0.285 content, nothing sits near zero | the catabolic three |

**The emptiness is re-derived from the data** with `in_box()` inside a `stopifnot`, over boxes that
cover the text extent rather than the anchor — not trusted from a render. **And the no-crossing
guarantee is a direct test, not a rule of thumb**: an orientation test over all 55 pairs of leaders.
That change was forced by the author's second note — *"swap OXPHOS assembly with central dogma to
avoid crossing lines"* — which is correct and which the old heuristic (lane order = the points' own
vertical order) would have forbidden: OXPHOS assembly sits *above* central dogma in y, but its
point lies further **left**, so ordering by y crosses the two leaders. Ordering by what actually
matters needed the real test. Same principle as Fig. S1C's three lanes.

**Scripts 40, 43 and 08 all describe these numbers and the script asserts they agree**: script 43's
arms to script 40's ruler at 1e−6 via `$defs$arms`, and script 08's mitoPPS diff to the ruler's
`p_tn` **exactly** (max difference 0). `require_fresher_than()` is deliberately *not* called on
`mitopps_scores.rds` — the same one-minute timestamp artefact Fig. 1F documents.

Size: 89 x 70 mm.

### Fig. S2B, built 2026-08-04 — and the whole-transcriptome half needed a control

Two scatters: what the wild-type gland does between six and twelve weeks against what Myc does at
six weeks. Left, 143 MitoPathways on the mitoPPS **priority** ruler (the ruler the word
*reallocation* refers to) — **r = −0.022**, and on the content ruler −0.038. Right, 15,191 genes on
the DESeq2 log2 fold change — **r = +0.228**.

**That +0.23 is a property of the design, not of the biology, and the panel had to say so.**
`myc_6W` = 6W_myc − 6W_wt and `6>12W_wt` = 12W_wt − 6W_wt: **the same six wild-type animals are
subtracted in both**, so a fluctuation in that baseline pushes the two contrasts the same way.
Two independent estimates of how much:

| | value |
|---|---|
| observed | **+0.228** |
| predicted from the standard errors (Cov = Var(6W_wt) = SE²/2) | **+0.198** |
| **measured, split baseline** — Myc effect against one half of the six wild-type animals, temporal effect against the other | **−0.025** (median −0.021, range −0.44 to +0.16, 60 % negative) |

**The split-baseline control is exhaustive, not sampled.** There are exactly `choose(6,3) = 20`
assignments and all 20 are used, so the panel has **no RNG in it** — nothing to seed and nothing to
drift. It runs on the VST matrix rather than on re-fitted DESeq2 models, because re-fitting is an
analysis and not a figure; **what licenses the surrogate is asserted** — with the baseline *shared*,
the VST version reproduces the drawn DESeq2 correlation to **0.002** (+0.230 against +0.228).
Splitting halves the baseline's n, which adds independent noise and pulls |r| toward zero; undoing
that inflation gives −0.032 rather than −0.025, so the conclusion does not depend on which is used.

**So the sentence is right on both scopes** — −0.02 mitochondrial, −0.03 whole transcriptome — but
the whole-transcriptome number only becomes negligible once the shared baseline is broken, and the
raw +0.23 is what a reader would otherwise compute. The dashed line on the right-hand panel is the
corrected relationship; the solid line is the observed one. Script 40's `$artifact_ledger` records
this class of artefact as **"STRUCTURAL — do not cite"**; this panel is the constructive version of
that warning, and script 40's `$split_runs` is the same technique applied to a different statistic
(Issue #6's wild-type convergence).

**The artefact pushes the wrong way for the mitochondrial panel too**, which makes that negative
conservative: the priority contrasts share the same wild-type baseline, so the sharing forces a
*positive* correlation, and what is observed is −0.02.

**Cosine is not correlation here, and script 40 reports the cosine.** On the content ruler the
uncentred cosine is **+0.210** and the correlation is **−0.038** — both content vectors have a large
positive mean and the cosine picks up that shared offset. The script asserts its recomputation
matches `$geometry` exactly on both rulers, so the panel is provably drawing script 40's vectors;
but the sentence's word is *correlation*, so the correlation is what is drawn. (Script 40's own null
agrees with the reading: the cosine sits at the **41st percentile** of its permutation null, i.e. at
chance, and the projection at the **9.8th**.)

**One bound that must not be dropped:** this is *not* the statement that the Myc+ gland's own
temporal change is independent of the wild-type one. That regression (`shared`, Myc+time ~ WTtime)
has slope **0.749** — much of the Myc+ gland's drift *is* the wild-type drift. Different pairs of
vectors; they must not be run together in the text.

**Two point layers, because the facets differ in density by two orders of magnitude** (author's
review, 2026-08-05: the mitochondrial points were almost invisible). 143 points at full opacity on
the left; 15,191 at 13 % on the right, where anything darker fills solid and the two lines cannot
be read over it. Ink is therefore not comparable between the panels — position is — and the legend
block says so. Both `r` annotations sit in a strip of headroom added by the y-scale expansion,
because at these densities no corner is reliably empty and a label landing on a datum is worse
than a slightly taller panel.

**A better way to check a panel, found here:** `qlmanage -t -s 2400 -o <dir> <file>.pdf` renders
through Core Graphics, which is exactly what the author sees in Preview. `sips` draws small
filled circles as **crosses** — a rasteriser artefact that had me about to change a point shape
that was never wrong. Use `qlmanage` when a glyph looks off.

Size: 89 x 52 mm; 109 of 15,191 genes fall outside the drawn window (none of the 143 pathways).

### Fig. S2C, built 2026-08-05 — a negative, with the control a negative needs

The port of `figures/fig04_substrate_specificity.R` panel D. Every transcript of the mitochondrial
death apparatus across the wild-type window: the **25 pro-** and **7 anti-apoptotic** MitoCarta
genes, plus **five more anti-apoptotic genes MitoCarta does not contain** — the caspase inhibitors
XIAP, cIAP1/2 (`Birc2`/`Birc3`) and survivin (`Birc5`), and the `Bcl2a1b` paralog — as their own
row, because "does the gland move its apoptotic transcripts" and "does it BUFFER against death"
are different questions.

**That row is labelled `IAPs & Bcl2a1`, not `brake (non-MitoCarta)`** (author's review, 2026-08-05:
*"it is not clear what brake (non-MitoCarta) means"*). The old name said how the genes were
*excluded* from the sets above rather than what they *are*. **The script now asserts the row's
membership**, so if `$buffer` ever gains a gene that is neither an IAP nor a Bcl2a1 paralog the
panel stops rather than mislabelling the row. Nothing in it moves either (padj 0.32–0.88).

**One of 37 moves, and it moves the wrong way for a loss of priming:** `Bnip3` **+0.63, padj
0.021** — a pro-apoptotic gene going *up*. Nothing else clears 0.05 on either arm or among the
brakes.

**"Overall priming" is a balance, and both arms move together, so the balance does not move.**
On script 40's ruler — the same instrument as Fig. 2F — across the window the pro arm shifts
**+0.025** and the anti arm **+0.045** on content (difference **−0.020**), and **+0.054** against
**+0.046** on the content-blind priority ruler (difference **+0.008**). For contrast, **Myc moves
the two arms apart**: +0.156 pro against −0.024 anti at six weeks, a difference of **+0.180** —
an order of magnitude larger, and in the opposite geometry. *The balance is something Myc changes
and development does not.*

**The power control, because a negative at n = 6 needs one.** On the same transcripts, in the same
libraries, at the same n, the **genotype** contrast reaches padj < 0.05 for **6 of 35** mapped
genes, against **1** across the window. The measurement can see movement in these genes; there is
none to see here. Adjusted p-values for the genotype contrast are Ensembl-keyed, so the symbols
are mapped through the annotation table and **the mapping is checked against baseMean rather than
trusted**.

**Seven genes are named, in italics** (author's review, 2026-08-05: the major players, not only the
mover). In mouse symbols, with the proteins the text calls them by: `Bax` (BAX) −0.145,
`Bak1` (BAK) +0.023, `Bcl2l11` (BIM) +0.255, `Bbc3` (PUMA) +0.062, `Bnip3` (BNIP3) **+0.635**,
`Bcl2` (BCL-2) +0.461, `Bcl2l1` (BCL-XL) −0.107. **Only `Bnip3` clears padj 0.05**; the rest sit at
0.37–0.99. `Bbc3` and `Bcl2l1` are the numerator and denominator of the ratio Figs. 2G and 2H turn
on, and both are flat across the window.

**Seven labels do not fit beside their own points, so they go into four lanes** — a reserved lane
above the pro row; the gap between the pro and anti rows for the one pro gene whose x-neighbour is
too close to share the top lane (`Bak1`, 0.039 from `Bbc3`); inline to the right for `Bcl2`, which
*is* the rightmost point of its row; and out to the left for `Bcl2l1`. **Every position is checked
against the data**: a half-width in x units is derived from the drawn text size, and `stopifnot`
requires that no two top-lane labels touch, that the rightmost stays inside the panel, that the
left lane is clear of data, and that the inline label really does belong to its row's extreme
point. Nothing here is nudged until it looks right.

**Deterministic vertical offsets, not jitter.** Within a row the genes are ordered by fold change
and the offsets cycle through a fixed ladder, so neighbours in x are separated in y by
construction. Better than random jitter at avoiding overplot, reproducible without a seed, and —
the reason it matters here — it makes every point's position *known*, so the three labels can be
placed rather than repelled. Two go into a lane above the top row that the y limit reserves; Bcl-xL
goes out to the left, into the half of the anti row that is empty. Both regions are **asserted**
empty, and every leader leaves the text on the text's own line (Fig. 2F's convention).

**One trap this panel does not fall into, and the legend says so explicitly.**
`MITOCARTA_APOPTOSIS_PRO`/`_ANTI` are MitoCarta sets, so a coupling between a priming composite and
a mitochondrial axis is mito-vs-mito and circular (script 34). **That warning does not transfer
here**: this is a temporal contrast on the transcripts themselves, not a coupling, and nothing on
the panel is correlated with a mitochondrial score. It would be easy to import the caveat wrongly.

**And the bound that must travel with it:** a negative on *transcripts* is not a negative on
*priming*. Priming is a property of the protein complement and of how close the mitochondrion sits
to the threshold; transcript levels are its substrate, not its measurement. BH3 profiling is the
measurement, and this panel is a reason to do it rather than a substitute for it.

Size: 89 x 38 mm.

### Fig. 2G, built 2026-08-05 — Fig. 1G's plot, one level up

**The sentence was corrected first** (author, 2026-08-05), and the correction is what the panel is
built to: *"While most apoptotic priming ratios established by MYC at 6W **remained stable, since
both pro- and anti-apoptotic proteins diminished** in accordance with the global rescaling of the
MYC effect, the PUMA/BCL-XL ratio showed a striking reversal, deviating significantly from the
expected pattern."* The earlier draft said the *ratios* were reduced; what declines together is
the **members**, which is why the ratios are left where Myc put them. The data say the corrected
version.

**The panel is Fig. 1G's axes at the level of ratios.** x = the Myc effect on a pro-apoptotic:
Bcl-xL log2 ratio at six weeks, y = the same at twelve, and the line is the **expected pattern** —
the Myc effect rescaled by the **global rate 0.487** (script 44's `global_rate_fitted`, the number
Fig. 1G quotes). A ratio that merely follows the rescaling lands on the line.

| ratio | 6W | 12W | retention | |
|---|---|---|---|---|
| `Bax` | +0.901 | +0.495 | **0.55** | on the line |
| `Bid` | +0.946 | +0.396 | 0.42 | on the line |
| `Bak1` | +0.688 | +0.185 | 0.27 | a little under it |
| **`Bbc3`** | **+0.674** | **−0.061** | **−0.09** | **crosses zero** |

**One shared denominator, which is what makes the comparison internal** — script 42's own design
note. Every ratio is against Bcl-xL, so *"it is just the global attenuation"* is refuted from
inside the panel rather than against an outside null: BAX priming retains the global rate while
PUMA priming reverses **against the same denominator**.

**The members are in the legend, not on the panel, and they are exact.** `Bax` retains **0.52** of
its six-week effect and `Bcl2l1` **0.57** — both at the global 0.49 — so their ratio retains 0.55.
`Bbc3` retains **−1.10** against that same denominator, which is the whole of the difference.
Drawing the members too would double the panel to make a point the ratios already carry.

**Shape is "significantly induced at 6W", and it is load-bearing rather than decorative** (author's
review, 2026-08-05: *"established" not very clear — "induced" or "significantly induced" would be
better*; the criterion is the unadjusted p of the six-week genotype coefficient). Retention is a
*quotient*, so it is meaningless where the six-week effect is not distinguishable from zero. Three
of the seven pairs are in that state (p6 = 0.40–0.71) and their positions carry no information
about loss — `Pmaip1`'s retention is +2.08 and `Bmf`'s −0.83 purely because their denominators are
near zero. **No residual is drawn for them**, and the four that are drawn are exactly the four the
sentence's scope ("established by MYC at 6W") names.

**"Significantly" needs its scope, and this is the item for the text.** The PUMA:Bcl-xL interaction
is **p = 0.037 raw, 0.017 purity-adjusted** — but its **BH across the nine pairs is 0.33**, and
against script 42's matched-pair null (random pro-like/anti-like pairs matched on expression,
conditioned on a six-week effect at least as large) the **empirical p is 0.082**, the 92nd
percentile of 510 matched pairs. It is a nominal result with a **pre-specified** licence — PUMA was
named in advance from the PGC1a cell experiments — not a multiplicity-surviving one. The honest
counterweight is in the legend: on the empirical null the most extreme pair is not PUMA but `Bmf`
(p = 0.025), and Myc never established that ratio (p6 = 0.40), which is precisely why the
*conditional* null is the one to read.

**Two numbers must not be interchanged.** The gene-level `Bbc3` interaction (**p = 0.0081**, Fig.
2H's subject) is a different statistic from the *ratio* interaction quoted here (**p = 0.037**).

**The line is imposed, not fitted**, and a bare line through a scatter reads as a regression — so it
is **named on the page in the manuscript's own words**: `global rescaling (slope 0.49)` (author's
review: `x 0.49` was unclear).

**That label lies ALONG the line, and it has to.** A horizontal caption cannot sit parallel to a
sloped line — over the width of the text the line climbs further than the text is tall, so one end
always drifts away from it (author: it *"slipped too high, not close enough to the line"*).
Rotating fixes it, and the angle is the line's angle **on the page**, which depends on the panel's
aspect: `tan(angle) = slope × (plot height / y range) / (plot width / x range)` ≈ 22°. The two
plot-area constants are measured for this panel size and the script says so — the clearance
assertion catches a label that has drifted off the line, but not one drawn at the wrong angle, so
they must be re-checked if the panel is resized. Clearance itself is a **perpendicular** distance
to the line, so it does not depend on the rotation being exactly right.

**The key sits in the bottom-right wedge** (author's review: at top left its text crossed the zero
line). The expected line climbs to the right, so beneath it in that corner is the one region no
point and neither reference line can reach — which also meant moving `Bbc3`'s label from below its
point to its left. The two Mcl-1 pairs are quoted in the legend rather
than drawn: a second denominator needs a second encoding, and the comparison only works with the
denominator held fixed. (`Bax:Mcl1` +0.487 → +0.408, retention 0.84, established; `Bbc3:Mcl1`
+0.260 → −0.149, −0.57, **not** established. PUMA reverses against both, but only the Bcl-xL one
was there to begin with.)

**One R trap, and it cost a run:** `wg[[col]][wg$gene == g]` where `$gene` contains NA returns
extra `NA` elements rather than dropping them — a logical index with NA does not subset, it
propagates. Use `which()`.

Size: 89 x 58 mm.

### Fig. 2H, built 2026-08-05 — and the panel carries the half that failed

Script 44 turns Fig. 1G's global rescaling into a **per-gene residual**: for each of 8,774 genes,
how far its twelve-week Myc effect departs from what the rate predicts, over its own standard
error, signed so that negative means *collapsed further than the dose explains*. The panel is that
distribution, with a rug of every gene in the bottom percentile and three genes marked.

| | z | genes below | percentile |
|---|---|---|---|
| `Foxo3` | **−2.53** | 39 of 8,774 | **0.46** |
| `Bbc3` | **−2.49** | 44 | **0.51** |
| `Bcl2l11` | +1.01 | 8,036 | 91.6 |

**"Closely tracking" has a sharp form and the script asserts it: only FOUR genes of 8,774 lie
between them.**

**Why `Bcl2l11` is on the panel.** The scan was **pre-registered with two genes** — script 44's
`pre_specified_genes` are `Bbc3` and `Bcl2l11`, named in advance from the cell experiments — and
**they split**: `Bbc3` collapses at percentile 0.51, `Bcl2l11` sits in the ordinary middle at 91.6
with an interaction p of 0.86. Drawing only the half that worked would be the wrong panel, so both
are marked.

**No key, and therefore no shape encoding either** (author's review, 2026-08-05: the panel is
obvious without one). The first version distinguished pre-specified from found-in-the-scan by
filled and open marks; with the key gone that difference would be unexplained ink, which is worse
than no difference. All three marks are now identical, and *which* gene was pre-specified — a fact
about the analysis, not about a gene's position — is in the legend block where it belongs.

**The axis reads `difference from global Myc scaling (z)`** (author's wording; kept in the
project's `Myc` casing so it matches Figs. 2E/2G/S1E).

**`Foxo3` was NOT pre-specified — it was found here**, and that is the panel's main caveat rather
than a footnote. Forty-four genes sit below the 0.5th percentile with it and they are mostly
unrelated (`Adh7`, `Agpat2`, `Aldh3a2`, `Ccnjl`, `Cntnap2`, `Inhba`, `Wfdc2`, `Vtcn1` …). **Position
alone is therefore not evidence.** What makes `Foxo3` a lead rather than one of forty-four names is
the **prior** — FOXO3 is PUMA's canonical p53-independent activator, with direct ChIP evidence —
and the text should introduce it that way round.

**The same departure is reached by two different routes, which the text should keep.** `Foxo3`:
Myc raises it at six weeks (+0.243) and lowers it at twelve (−0.242), and the **wild-type gland
raises it with age (+0.473) while the Myc+ gland does not (−0.013)**. `Bbc3`: +0.258 → −0.283, with
the wild-type gland flat (+0.062) and the Myc+ gland falling (−0.479). One is a failure to follow
the normal developmental rise; the other is an active loss.

**None of the interactions survives correction** — `Bbc3` p = 0.0081 (BH 0.84), `Foxo3` p = 0.0074
(BH 1.00), `Bcl2l11` p = 0.86 — which is the expected outcome for an interaction at n = 6 per cell
and is why the licence is pre-specification rather than the p-value. **A z of −2.5 among 8,774
genes is a rank statement, not a test**, and the legend says so.

**And the arrow is not the obvious one.** PUMA restrains the mitochondrial pyruvate carrier (Kim,
*Cancer Cell* 2019), so "PUMA falls, therefore respiration falls" is backwards; respiration sits
upstream (Dey & Moraes), and FOXO3 → BBC3 closes it into a negative-feedback circuit rather than a
linear chain. Nothing on this panel establishes a direction.

**One nicety worth keeping:** neither gene is in the 2,648-gene set the rate was fitted on (their
six-week adjusted p-values are 0.22 and 0.21), so neither contributes to its own expectation.

Size: 89 x 44 mm; windowed at z = +4, with the four genes beyond it named in the legend (they are
genes Myc *retains* better than the dose predicts — `Myc` itself is there at z +5.03, the 99.99th
percentile, which is Fig. S1E's result read on this scan).

### Fig. S2D and Fig. 2I, built 2026-08-05 — the last two, run here rather than in Positron

The author was away and authorised a one-off deviation from Option A: these two were written, **run
and their PDFs generated in this session**. Everything else about them is unchanged — the same
assertions, the same legend blocks, the same `rebuild_panels.R`.

**Fig. S2D — the alternative to Foxo3, tested.** Fig. 2H's obvious objection is that some *other*
upstream input to PUMA moved and Foxo3 is a bystander. The roster is script 42's own
`exclusions$puma_inputs` — the p53 family (`Trp73`), the E2F arm (`E2f1`), the integrated stress
response (`Atf4`, `Ddit3`, `Trib3`, `Chac1`, `Eif2ak3`, `Nupr1`, `Sesn2`) and the FOXO family
(`Foxo1`, `Foxo3`, `Foxo4`) — drawn as the Myc effect at each age, joined. **Nothing reaches
padj < 0.05 at either age and no interaction survives** (asserted; smallest adjusted p 0.19 at six
weeks, 0.70 at twelve, all interaction padj = 1).

- **`Foxo3` is in the roster and is drawn, in bold**, because it is the one gene here that is *not*
  "other": +0.243 → −0.242, interaction −0.486, the same sign change as `Bbc3`. Leaving it out
  would make the panel look like a clean negative when it is a negative **with one exception**, and
  the exception is the preceding sentence. **The text should name it as the exception** rather than
  let the reader assume the roster excludes it.
- **The other two FOXO paralogues do not do it** — `Foxo1` −0.184 → −0.300, `Foxo4` +0.071 →
  −0.072 — so the Foxo3 result is not a family-wide effect.
- **One gene moves on a different axis and is named in the legend so the panel is not read as
  flatter than it is:** `Trp73` rises **+2.14 across the wild-type window at padj 0.0039** — on a
  mean expression of **12 counts**, at the floor of what this design can measure, and on the
  batch-confounded contrast. Not a Myc effect, not on this panel's axis.
- **Bound that matters:** a transcript-level negative is not a pathway-level one. The ISR and the
  p53 axis act substantially through protein stability, phosphorylation and localisation — ATF4 in
  particular is translationally controlled and can be fully active with an unchanged message. What
  this excludes is a *transcriptional* re-routing of PUMA's inputs, which is the alternative the
  text raises.

Size: 89 x 52 mm; key under the plot (Fig. S1F's solution — twelve full-width rows leave no empty
corner, and inside it landed on the bottom row).

**Fig. 2I — what an interaction between a genotype and a coupling looks like: two slopes that
differ.** One point per animal, its mitoPPS priority score against its PUMA:Bcl-xL log ratio, a
line fitted within each genotype.

| axis | wild type | Myc+ | script 43's interaction |
|---|---|---|---|
| **OXPHOS subunits** | **−4.80** | **+4.35** | **+6.09, p = 0.0052** |
| redox | −7.70 | −6.72 | +2.16, p = 0.72 |

**On the OXPHOS axis the two genotypes couple in opposite directions; on redox they are parallel.**

- **Redox is the right control precisely because it is not a null axis.** Both genotypes couple to
  it, and strongly (wild-type slope −7.70). What it does not do is couple *differently*. An axis
  nothing couples to would not have excluded the possibility that the PUMA ratio simply tracks
  mitochondrial scores in general.
- **The panel is a reconstruction and says so.** Script 43 fits on its own log matrix; this rebuilds
  the ratio from the VST matrix and **reproduces the recorded interaction to 1.2 %** (asserted).
  The claim is checked tightly and the control qualitatively — a relative error on a near-null
  coefficient is not a meaningful quantity, so redox is required only to still read as nothing
  (refit p 0.74). **The number for the text is script 43's.**
- **Read it as a lead.** The permutation null in the same object puts it at the **91.7th percentile
  of 5,000 draws (empirical p = 0.083)**, and adding a timepoint term moves p from 0.0052 to 0.088.
  It is the one axis-by-genotype interaction in the corpus that reaches nominal significance.
- **A wording item:** *"no such link exists for other redox and metabolic axes"* — script 43's
  trade-off analysis carries **exactly two axes**, OXPHOS subunits and redox. The sentence should
  say *"for the redox axis"* unless a wider panel is run. What *is* broader is the outcome side:
  two axes against five outcomes, and only this one reaches p < 0.05.
- **Do not quote the raw correlations.** Pooled across genotypes the PUMA ratio correlates +0.04
  with the OXPHOS axis and −0.48 with redox — the pooling artefact two opposite slopes produce.

Size: 89 x 56 mm; four-group sample palette on the points, genotype on the lines, key under the
plot.

### ROUND 2 (2026-08-09) — four alternatives for Figure 2, on one grammar

**Why.** The section is written and every commissioned panel is built, but the panels read as
supplementary material: each answers its own sentence carefully and none carries the argument.
The author's requirement — *"the figures have to talk by themselves, people do not read the
text"* — and the message they must deliver:

> A complex interplay between the WT and Myc-driven timelines leads to **OXPHOS subunit repression,
> more than the overall MYC effect decline**; that decline is coupled to the **selective loss of the
> PUMA/Bcl-xL ratio**, presumably through **Foxo3**; the lost trigger is why the adult gland cannot
> die.

The reworded Doc section puts **both timelines side by side in every claim**, and that is exactly
what the built panels do not do. **Every current panel is kept.** These four are alternatives; the
author picks per slot.

**THE GRAMMAR: the two timelines.** `x = 6>12W_wt`, `y = 6>12W_myc`, identity line = *development
alone*. Declared once as **`two_timeline_base()`** in `_panel_common.R` so the four cannot drift
apart. Quadrants gain a fixed meaning the reader learns once, and the device is **exact, not
approximate**: script 40's ruler satisfies `c_tp = c_tn + c_int` to floating point (1.2e-16 on
content, bit-identical on priority, asserted), so **a vertical drop from the diagonal IS the
interaction**. `coord_equal()` is not cosmetic — both axes are the same quantity, so equal scaling
is honest *and* it pins the diagonal at exactly 45°, which is what lets its label lie along it
without the aspect arithmetic Fig. 2G needed.

**Talking within the publication rule.** `theme_panel()` still blanks title/subtitle/caption. What
carries meaning is **naming drawn elements** — `development alone` along the diagonal,
`below the line = lost under Myc as well` in a corner — the device Fig. 2G already uses for its
`global rescaling (slope 0.49)` line. Naming a drawn element is what a key does.

| slot | script | what it draws |
|---|---|---|
| **Fig. 2F (alt)** | `fig2_two_timelines_mito_alt.R` | 143 MitoPathways on the plane; arrows from the diagonal down to eleven named arms = the Myc-specific loss |
| **Fig. 2G (alt)** | `fig2_death_two_timelines_alt.R` | the death machinery, per gene, same plane |
| **Fig. 2H+I (alt)** | `fig2_puma_chain_alt.R` | **A** the chain per animal (heatmap) + **B** the OXPHOS→PUMA coupling |
| **Fig. S2D (alt)** | `figS2_p53_arm_alt.R` | the p53 arm, same plane, with `Bbc3`/`Foxo3` as the movers |

**Fig. 2F (alt).** Development **raises** the compartment (median +0.041, 72 % above zero) while
**stripping OXPHOS**; the Myc+ gland lowers everything (93 % below the diagonal); OXPHOS ends
deepest because the two terms **add** — subunits −0.255 in development, −0.150 more from Myc,
−0.405 under Myc. The **assembly-factor control is on the panel**: +0.001 in development, and only
the global Myc term (−0.175) after. One ruler deliberately — drawing both rulers *and* both
timelines at 89 mm is what made the first version unreadable; the priority ruler is in the legend
and the sandbox redraws the panel on it with one substitution.

**Fig. 2G (alt).** `Bbc3` (+0.06 → **−0.48**, padj 0.016) is the only BH3-only sensor that is flat
in development and lost under Myc — asserted. `Bcl2l1` sits at the origin, so the ratio's reversal
is **all numerator**. `Bmf` (+1.02 → +0.84, significant on both) sits **on** the diagonal: the one
death gene the window itself moves, and the panel's own negative control. **`Bax` also falls
significantly under Myc (−0.376, padj 0.008)** and is visible — it is an effector, so "specific to
`Bbc3`" holds *among the BH3-only sensors*, and the sentence should say so.

**Fig. 2H+I (alt).** Part A's premise had to be corrected against the data mid-build: it is **not**
a block that falls together. Group means of the standardised rows:

| | OXPHOS | `Foxo3` | `Bbc3` | `Bcl2l1` | PUMA:Bcl-xL |
|---|---|---|---|---|---|
| 6W_wt | +0.07 | −0.90 | −0.12 | +0.74 | −0.55 |
| 12W_wt | −0.78 | **+0.89** | +0.18 | +0.39 | −0.16 |
| 6W_myc | +0.95 | +0.01 | +0.84 | −0.71 | +0.94 |
| 12W_myc | −0.25 | −0.00 | **−0.90** | −0.42 | −0.23 |

Read **down the two timelines**: in the wild-type gland OXPHOS falls, **`Foxo3` rises** and PUMA is
untouched; under Myc OXPHOS falls from a higher start, **`Foxo3` fails to rise** and **PUMA
collapses**. That is the text's own sentence, and `Bcl2l1` follows neither. Three build decisions
are on the record: the fill is **clipped at 2 SD** (one animal at |z| 2.7 otherwise leaves every
tile pale); each group gets a **bordered mean column with its value printed**, because a mean of
six z-scores is necessarily paler than its animals and would be the weakest mark on the panel
exactly where the claim lives; and the number is **rotated 90°**, since a 28-column tile is 2.5 mm
wide and `+0.95` is not. Part B is the current 2I, reconstructed and asserted to 1.2 %.

**Fig. S2D (alt).** The author's point — *the p53-dependent arm is not moving* — needs
`exclusions$p53_axis`, not the `puma_inputs` roster the current S2D draws (which contains one
p53-family gene). The arm forms a **horizontal band at y ≈ 0**, not a cluster at the origin: several
drift on the wild-type axis (`Eda2r` −0.99 at padj 0.019, the one p53-arm gene to clear 0.05 on
either timeline, on 101 counts — named in the legend so the panel is not read as flatter than it
is) while the whole arm stays within ±0.26 under Myc. **`Bbc3` leaves the band; `Foxo3` leaves the
diagonal** — two different ways of moving, and the grammar shows both.

**Circulation.** `panels_to_pdf.R` now writes a **third** document, `Figure2_alternatives.pdf`
(4 pages), matched on the `(alt)` slot suffix, and the two committed documents **exclude** that
suffix — so they are unchanged at 12 and 8 pages, verified.

**Text items this round turned up.** `Foxo3` is **rank 1 and `Bbc3` rank 2 of the 26
biogenesis/cell-death mechanism genes by interaction p** (0.0074, 0.0081; gap to third `Pgs1` at
0.066) — the framing the reworded text already uses, and far stronger than the 8,774-gene
percentile the current 2H draws. `Bmf` rises in **both** timelines and belongs in the text as the
contrast to `Bbc3`. `Bax` falls significantly under Myc, so "specific to `Bbc3`" needs "among the
BH3-only sensors". `Trp73` moves on the wild-type timeline (+2.14, padj 0.004) on **12 counts**.

### ROUND 2 REVIEW (2026-08-09) — the author's changes, panel by panel

The four alternatives were reviewed together and each got a list. Recorded here as they land.

**Fig. 2F (alt) — done.** Five changes, all made.

1. **No drop arrows.** The vertical distance from the diagonal is still the interaction and still
   exact; twelve arrowheads over a 143-point cloud cost more ink than they returned.
2. **One continuous fill, on the wild-type log2FC** (`heat_fill`, the manuscript ramp, symmetric
   ±0.558, white pinned to zero), replacing the two-colour respiratory/catabolic rosters. The
   rosters survive **only as the twelve labels**. The fill is **redundant with the x axis by
   construction** — that is the point, and the legend block says so: the plane sorts by what
   development did without a roster deciding it for the reader.
   *Consequence that forced a geometry change:* zero on this ramp is `#FAFAFA`, so a solid point at
   the compartment median (+0.041) would be invisible. Points are now **outlined circles**
   (`shape = 21`, grey45 stroke), labelled arms larger with a darker stroke.
3. **`development alone` moved to the top of the line.** Its glyph run is a 14 mm segment lying
   along the diagonal, so the placement is checked with a real **point-to-segment distance in panel
   fractions**, not a bounding box. The binding constraint is catechol metabolism, 1.9 mm below the
   run; at this font size **no placement on the upper diagonal clears both** that point by ≥2.6 mm
   *and* the corner — pushed further up, the last letter clips (it did, and `qlmanage` caught it).
   Anchor `diag_at = 0.79`, assertion threshold 0.025 of the panel, with the reason in the script.
4. **`below the line = lost under MYC`** (was "… under Myc as well").
5. **Both axes carry `(set average log2FC)`.**
6. **The key is horizontal and inside, top left** (second pass): laid the same way as the axis it
   repeats, so the eye reads the two as one scale. The quadrant it occupies is the empty one —
   x below zero with y above it would be a pathway the normal gland strips and Myc restores, and
   there is none — asserted rather than assumed. Moving it off the right also let the square data
   region grow to the full column width. One trap: `theme_classic` sets `axis.text` to `rel(0.8)`
   of `base_size`, so a legend text size copied from `base_size` prints **larger than the axis it
   repeats**; the key's numbers are set to 4.8 explicitly.

**Fig. 2G (alt) — done.** Four changes, all made.

1. **`below the line = lost under MYC` moved to the bottom right** (`quadrant_at = c(0.99, 0.02)`),
   away from the labels.
2. **`development alone` moved to the bottom of the line** (`diag_at = 0.15`). `Bmf` at +1.02 sets
   the limits, so the seventeen other transcripts occupy the middle third and **the outer thirds
   are genuinely empty** — asserted (`!any(fx < 0.36)`), which is what makes both moves safe.
3. **Continuous fill on the Myc timeline** (the vertical axis, which is the axis the sentence is
   about), `heat_fill` symmetric ±0.925.
   *This forced the significance encoding to move from the SHAPE to the RING* — an open shape has
   no fill, and only four of eighteen transcripts are significant, so the old `21`/`1` pair would
   have left the ramp on four points. Now every point is `shape = 21` and the ring is dark for
   `padj < 0.05` under Myc, pale for n.s.
4. **Both axes read `(log2FC)`** (the x axis said "raw log2FC"; "raw (unshrunken)" stays in the
   legend block, where the shrinkage rule belongs).

Second pass, two more:

5. **The bar is vertical**, laid the same way as the y axis it repeats — the mirror of 2F (alt),
   whose fill is the x axis and whose bar is horizontal. The **wide** key (the rings, with its long
   `padj` label) takes the top row where there is room; the **narrow** bar hangs below it down the
   empty left edge.
6. **Significant gene names print in the declared significance ink** (`sig_cols[["sig"]]`,
   `#E41A1C`) — `Bmf`, `Bax`, `Bbc3`, `Htra2`, exactly the four with `padj < 0.05` on the Myc
   timeline. Passed as a per-row constant in the data's own order, not as an aesthetic, so the
   panel does not grow a third key.
   *This widened `sig_cols`' declared scope*, which said it applied **only** to bracket labels. The
   rule behind it is unchanged and still binds — where the geometry carries the effect size,
   significance stays a footnote and never becomes a mark's fill. A per-gene scatter is the case
   where the geometry *is* the effect size on both axes, so the ring and the red name are
   annotation rather than the reading. Recorded at the declaration.

**Not changed, and why:** the seven transcripts crowding the origin (`Cycs`, `Xiap`, `Casp3`,
`Apaf1`, `Diablo`, `Bak1`, `Bcl2l1`) force one long `Apaf1` leader. Four repel seeds were rendered
and compared; the leader is structural, not a seed accident, so the approved layout (`seed = 6`)
stands.

**Fig. 2H+I (alt) — rebuilt.** The author's verdict on the per-animal heatmap: *"the heatmap does
not work, the individual data are too much variable"*. Right about the data, and the display was
the problem — at n = 6 the reader had to average six noisy tiles by eye in four places at once.

**Part A left — the chain as distributions, in Fig. 1E's idiom** (box, every animal, brackets).
Three facets in **the order the chain runs** — OXPHOS mitoPPS, `Foxo3`, PUMA:Bcl-xL — not the order
they were listed in, so the panel reads the way the mechanism does. Each facet carries **the two
within-genotype temporal contrasts**, which is the whole argument in six brackets:

| | wild type 6>12W | Myc+ 6>12W |
|---|---|---|
| OXPHOS mitoPPS | −0.11, p 0.108 | −0.15, **p 0.027** |
| `Foxo3` | +0.46, **p 0.0044** | −0.00, p 0.982 |
| PUMA:Bcl-xL | +0.17, p 0.517 | −0.50, **p 0.019** |

The respiratory arm falls on **both** timelines; `Foxo3` rises **only** in the wild-type gland; the
priming ratio collapses **only** under Myc. Significant brackets in the declared significance ink,
as in 2G (alt).

**One instrument across the panel, and it is checked.** Neither the mitoPPS axis nor the ratio has
a DESeq2 test, so all six brackets are OLS on the drawn values (Fig. 1E's `fit_simple` idiom) —
drawing a Wald p on the gene facet and a t on the other two would put two instruments on one panel.
Where a record exists it agrees: `Foxo3` +0.473 padj 0.0078 / −0.013 padj 0.977, `Bbc3` +0.062 padj
0.881 / −0.479 padj 0.016, `Bcl2l1` n.s. on both. Same signs, same side of 0.05 — **asserted**.

**Part A right — the two members of the ratio as group medians** (2 × 4, the manuscript ramp, value
printed in each cell). `Bbc3` tracks the ratio (+0.7 → −0.8 under Myc, flat in the wild type);
`Bcl2l1` follows neither timeline. The reversal is all numerator.
*Two layout facts worth keeping:* written-out column labels ("12W_myc" rotated under a 6 mm column)
cost **10 mm of height**, and patchwork then pads the boxplots' blank x axis to match — that gap
between the two halves was the first version's. The columns instead carry a **strip of the same
sample colours**, keyed by the one legend under part A; `colour` is free because the tiles use
`fill`, so it needs no second scale. And `coord_fixed()` is required, or the two rows stretch to
whatever height the panel beside them needs and a 6 × 20 mm tile stops reading as a heat map.

**Part B — one line, two numbers.** The WT fit is gone, the Myc+ fit is dashed and carries its R²
(0.35), the y-axis title sits closer to the axis, and part B's duplicate group key is dropped (part
A names the four groups directly above, in the same palette).

**The interaction p is printed beside the R²** (author, 2026-08-09), and **the ink says whose
number it is**: the R² in the Myc+ colour, because it belongs to the dashed fit; `interaction
p = 0.0052` in neutral ink, because it belongs to both genotypes and is what the text quotes.
This was needed — the wild-type animals trend the **opposite way just as strongly** (R² 0.33,
p 0.052), so with one line drawn and one number, a viewer reads a Myc-specific correlation where
the result is a genotype-by-slope interaction.

**Asked and answered: is `Foxo3` 12W_wt vs 12W_myc significant? No — and the reason is the
result.** On the drawn values −0.89 z, **p = 0.092**; DESeq2 −0.242 log2, p = 0.059, **padj = 0.53**.
The 6W genotype contrast is its **mirror**: +0.92 z (p = 0.076), DESeq2 +0.243 (p = 0.057).
Myc raises `Foxo3` at six weeks and lowers it at twelve by almost exactly the same amount; neither
half clears 0.05 alone, and that is precisely why the **interaction** is the strongest term
(−0.486, **p = 0.0074**, rank 1 of the 26 mechanism genes). Labelling either genotype contrast
would print a failed test next to the panel's actual claim. Not drawn.

**Added: one interaction p per facet** (author, 2026-08-09) — the batch-clean quantity, and the
difference the two lower brackets already draw. Drawn as a **second-tier bracket from 1.5 to 3.5**,
i.e. from the midpoint of one lower bracket to the midpoint of the other, because that is what a
difference of differences is:

| | wild type 6>12W | Myc+ 6>12W | difference |
|---|---|---|---|
| OXPHOS mitoPPS | p 0.108 | **p 0.027** | p 0.615 |
| `Foxo3` | **p 0.0044** | p 0.982 | **p 0.013** |
| PUMA:Bcl-xL | p 0.517 | **p 0.019** | **p 0.042** |

Same instrument as the six below it, so the three numbers in a facet are comparable; DESeq2's
record for the two transcripts agrees (`Foxo3` −0.486, p 0.0074; `Bbc3` −0.541, p 0.0081 — ranks 1
and 2 of the 26 mechanism genes, neither clearing genome-wide BH, which is why pre-specification is
what licenses them). The `Foxo3` agreement is asserted.

**Naming.** `OXPHOS priority` → **`OXPHOS mitoPPS`** wherever it is drawn, and
`fig2_oxphos_puma_coupling.R`'s axis follows ("mitoPPS score, per animal", generic because its
facets are OXPHOS and redox). `priority` survives only as the word for the **ruler** in prose,
where it is the counterpart of `content` — so Figs. 1F and 2F, which draw the pathway-level ruler
rather than this per-animal score, are untouched.

**Fig. S2D (alt) — done.** Two changes.

1. **The two movers carry their own padj on the page**, each for the timeline it actually moves on:
   `Foxo3` **padj 0.0078** (wild-type window — its displacement is along x) and `Bbc3` **padj
   0.016** (Myc+ window — its displacement is down y). Which timeline is asserted, so the legend
   block cannot drift from the drawing. Built as a **plotmath expression on one line**
   (`italic("Bbc3")~"padj 0.016"`) so the gene stays italic and the number does not — a single
   `fontface` cannot do both. `atop()` was tried first and set the two halves so far apart that the
   two labels interleaved and you could not tell which number belonged to which gene.
2. **The class label carries the threshold**: `the two that move` → **`moves (padj < 0.05)`**.
   Not `padj < 0.001` — the two values are **0.0078 and 0.016**, so 0.001 would be a claim neither
   supports.
3. **Both axes read `(log2FC)`** (the x axis said "raw log2FC").

Note on the red convention: 2G (alt) marks significant genes with `sig_cols[["sig"]]` because it
has to pick four out of eighteen. Here the movers are already **their own class with their own
colour and the threshold in the key**, so the class colour *is* the significance mark and a second
one would be redundant. Same rule, different instrument.

**One shared bug fixed on the way:** `two_timeline_base()` declared `quadrant_hjust` and then
hard-wired `hjust = 1`, so a caller asking for a centred note silently got a right-aligned one.
Both existing callers now pass 1 explicitly, which is what they already render — **nothing moved**.
Had it been fixed silently, 2F's note would have centred straight onto the glycine-cleavage point.

**A note on `MYC` vs `Myc`.** The quadrant text is the author's wording verbatim. The axes on the
same panel say `Myc+`, so the panel currently spells the oncogene two ways; say the word and it
becomes `Myc` everywhere.

### The circulation document — `panels_to_pdf.R`, added 2026-08-04

**Two documents as of 2026-08-05**, one per figure, because that is how the Results is read:

| file | pages | page size |
|---|---|---|
| `outputs/figures/panels/Figure1_and_S1_panels.pdf` | 12 | 105 x 111 mm |
| `outputs/figures/panels/Figure2_and_S2_panels.pdf` | **8** | 105 x 97 mm |

Every built panel, **one to a page, in slot order, with the slot as the page title**, for
collaborators who want to see the panel each sentence cites without assembling a figure. One
`source()` writes both:

```
source(here::here("figures", "panels", "panels_to_pdf.R"))
```

**Which document a panel lands in is decided by the slot its own `panel_legend()` declares**, not
by a list kept here — so a panel that changes figure moves documents by editing its own script.
Anything whose slot matches neither figure is reported and skipped (currently
`figS1_mb_fork_specificity`, which is still cited nowhere).

**It re-renders rather than concatenating the existing PDFs**, because this machine has no
Ghostscript, qpdf, pdftk, pdfjam or LaTeX and neither the `qpdf` nor the `pdftools` R package is
installed — nothing available can place one PDF page inside another. Re-rendering is also the
better answer: it stays vector, it cannot go stale against the scripts, and the page titles come
from the `panel_legend()` blocks rather than from a list kept by hand.

**Every panel is drawn at its designed physical size**, centred on a page just big enough for the
largest. Scaling a ggplot to fill a page does not magnify it — the type stays at 6 pt while the
plot area stretches, which is a different picture from the one that was signed off. At true size
each page looks exactly like the manuscript figure will, and the reader magnifies with the
viewer's zoom.

**Sizes are read from each panel's own `save_panel_p()` call** through a new `myc.fig.capture`
hook in `_panel_common.R`, which records the designed size and writes nothing. So this file keeps
no second copy of a number that lives in the panel script. The registry sits in the global
environment because each panel is sourced into a child of it.

**One R trap, and it cost a debugging round:** `on.exit()` at the *top level* of a sourced script
attaches to the frame evaluating that single expression, so it fires immediately. Here it restored
`myc.fig.capture` before a single panel had run, and every panel wrote its own PDF instead of
reporting its size. There is no `on.exit()` anywhere in that file now; the option and the device
are restored explicitly at the foot.

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
- **"MYC *significantly* promoted … OXPHOS subunits (all complexes)" — the word `significantly`
  does not reach that far, and this is the one item on the list that needs a real change.**
  Checked against script 08's own BH-adjusted p across all 144 pathways at 6W (33 clear 0.05):

  | in the sentence | priority at 6W | padj | verdict |
  |---|---|---|---|
  | translation | +0.117 | **0.042** | supported (translation factors +0.143, padj 0.0043) |
  | pyruvate | +0.328 | **0.042** | supported |
  | serine | +0.212 | **0.016** | supported |
  | protein import / homeostasis | +0.117 | 0.053 | **marginal at tier level** — but SAM 0.016, chaperones 0.039, protein homeostasis 0.044 all clear it |
  | glycine cleavage | +0.405 | 0.118 | **not significant**, and it is the *largest* promotion on the list |
  | **OXPHOS subunits** | +0.112 | **0.274** | **not significant** |
  | — CI / CII / CIII / CIV / CV | +0.055 / +0.021 / +0.160 / +0.017 / +0.165 | 0.37 / 0.54 / 0.16 / 0.86 / 0.15 | **not one complex is significant** |
  | dynamics & surveillance | −0.156 | **0.028** | supported |
  | fission | −0.202 | **0.016** | supported |
  | apoptosis | −0.247 | **0.043** | supported |
  | mitophagy | −0.182 | 0.084 | **not significant** |
  | calcium signalling | −0.075 | **0.026** | supported (calcium homeostasis; the uniporter is −0.306 at padj 0.060) |

  So: no part of the OXPHOS arm is significant on this ruler, and the **one OXPHOS row that is,
  is `CIII assembly factors` at −0.141 (padj 0.022) — demoted.** Suggested repair: move OXPHOS out
  of the "significantly promoted" clause and into a separate one — *"the OXPHOS subunits moved in
  the same direction (+0.112) while their assembly factors did not (−0.028), a split that returns
  in the wild-type timeline"* — which is a more interesting sentence and is the one section 2
  needs. Fig. 1F draws all of this with filled/open points, so the panel and the text must agree.
- **"MYC significantly promoted … biosynthetic pathways including glycine cleavage, pyruvate, and
  serine metabolism"** is right, and it is worth knowing those are the *largest* single
  promotions (+0.41 / +0.33 / +0.21) — larger than OXPHOS, which leads as a **tier** and not as a
  leaf. The sentence as written is fine; do not let a later draft turn it into "OXPHOS leads".
- **"the expression gap remained constant *between 6W_myc and 12W_myc*"** — a gap is between
  genotypes, not between ages. Both readings happen to be true and Fig. S1E draws both (the
  genotype gap +1.68 → +1.82; no decline within Myc+, −0.10 at padj 0.76), but the sentence should
  say which it means.
- **`Bbc3` is written "(p < 0.01)"**; the interaction p is **0.0081** and its genome-wide BH is
  **0.84**. What licenses the test is that `Bbc3` was **pre-specified from the cell experiments**,
  and the text does not yet say so. One clause, once.
- **Fig. 1G's "overall R2, OXPHOS R2, p" placeholders — RESOLVED 2026-08-04**, see the Fig. 1G
  entry above for the table. Short form: whole transcriptome slope **0.487** / R2 **0.51** over
  the 2,648 genes Myc moves (0.450 / **0.381** over all 8,774 — say which); mitochondrial slope
  **0.552** / R2 **0.801**; **OXPHOS R2 0.761** on the content ruler (*not* 0.945, which is the
  priority ruler and belongs with the "see Fig. 1F" clause); **p < 0.002**, and that p attaches to
  the **R2 only** — the slope is reached by 64 of 500 shuffles (p = 0.13).

### Added 2026-08-04, from Fig. 2E

- **"the overall proliferation signalling remained relatively stable" — right, and it needs one
  qualifier.** On the enrichment ruler it is as clean as it gets: **not one of the 14 curated
  proliferation sets is significant** on `6>12W_wt`, and the family median (−0.89) sits inside
  the ±1 band an unmoved set occupies — against the TEB programmes' −1.58 with 29 of 38
  significant. But on the gene-level ruler `PROLIF_*` pooled (731 genes) moves **−0.047 log2 at
  the 1.3rd percentile of its own expression-matched null, p = 0.013**. It is *small, not
  immobile* — script 43's own phrase. Suggested: keep "relatively stable" and make it
  quantitative — **the TEB arm moves 8.9 times as far**.
- **The generic proliferation regulons agree with the curated sets and the lineage-context lanes
  of the same transcription factors do not.** `TFT_E2F1_DOROTHEA_ABC` −1.18 and
  `TFT_E2F1_CHIPATLAS` +0.86 are both n.s., but the Gray `_LE` lanes fall hard (FOXM1_HS_LE
  −2.10, TFDP1_HS_LE −2.08, MYBL2_HS_LE −2.04, E2F1_HS_LE −1.72, all padj < 1e−3). Those are
  lineage-identity lanes, not cell-cycle sets. **What the normal gland withdraws from is lineage
  and morphogenesis; the cell cycle itself does not move** — the same reading Fig. 1H's
  off-diagonal clouds gave, arriving where the narrative wants it. Worth a clause.
- **The direction is controlled, in the legend block rather than on the panel** (author's call;
  see the Fig. 2E entry): the ductal `_DN` sets rise where the TEB `_UP` sets fall, in two of the
  three lineages. If the sentence wants a stronger verb than "lost", *"the gland moved from the
  end-bud programme toward the ductal one"* is what the data say — and that phrasing is worth
  more here than the control being drawn, because the text is where it now lives.

### Added 2026-08-04, from Fig. 2F

- **"across all complexes" has one exception.** CI, CIII, CIV and CV subunits all fall on both
  rulers; **CII subunits rise (+0.088 content, +0.026 priority)**. Complex II is the only
  respiratory complex with no mtDNA-encoded subunit and the only one outside the proton circuit —
  but the set is **four genes**, so this is a direction to note, not a mechanism. One-word fix:
  *"across the four complexes that carry mtDNA-encoded subunits"*, or name them.
- **"OXPHOS subunit LFC and MitoPPS drop" is the strongest form of the claim and it is worth
  saying why both are quoted:** a content fall alone could be normalisation, a content-blind
  ratio fall alone could be a reshuffle inside a growing compartment. Across the 143 pathways the
  two rulers agree at **Spearman 0.82**.
- **The clause the sentence is missing is the background.** The median MitoPathway **gains**
  content over this window (+0.041, 72 % above zero), so the respiratory arm is falling *against*
  a rising compartment, not with it. That is what licenses "withdraw", and it is one clause.
- **The assembly-factor control deserves a half-sentence in the text**, not just the panel: the
  assembly factors of the same complexes sit at **+0.001 / −0.002**, percentile 50.2 of the
  matched null. *The gland withdraws the stoichiometry of the chain, not the machinery that
  builds it.*
- **The statistic to put in the sentence is the PAIRED null, not either arm's own.** "The
  respiratory chain falls *while* biogenesis does not" is a comparative claim, and script 43
  tests it by redrawing both sets together: **OXPHOS subunits minus mitoribosome = −0.230,
  p = 0.0005** (minus nucleotide metabolism −0.255, p = 0.0085; minus the pooled proliferation
  set −0.208, p < 0.0005). Per pathway nothing survives BH on either ruler — 0 of 143 — so the
  text should not attach a p-value to any single MitoPathway.
- **The mtDNA arm is not on the panel and the text should not quote it either.** It is the largest
  movement in the compartment (+0.633 / +0.525) and it is the one quantity that is
  time-associated rather than genotype-associated, so the wild-type temporal contrast is exactly
  where it cannot be read.

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
| **2E** | "the TEB signature was lost as expected in the puberty-adult transition, [while] the overall proliferation signalling remained relatively stable" | **BUILT 2026-08-04** `fig2_wt_teb_proliferation.R` | `fgsea_percategory.rds` (drawn), `substrate_specificity_tradeoff.rds$wt_null` + `$defs`, `priming_arm_teb.rds$teb_signatures` |
| **2F** | "MECs withdraw from the respiratory chain; OXPHOS subunit LFC **and** MitoPPS drop across all complexes while biogenesis pathways remain relatively stable … the maturing gland upregulates amino-acid and lipid catabolism" | **BUILT 2026-08-04** `fig2_wt_mito_contraction.R` | `background_vs_myc.rds$ruler` (drawn), `substrate_specificity_tradeoff.rds$wt_null` + `$comparator_priority` + `$defs` |
| **2G** | **sentence corrected 2026-08-05:** "While most apoptotic priming ratios established by MYC at 6W remained stable, since both pro- and anti-apoptotic proteins diminished in accordance with the global rescaling of the MYC effect, the PUMA/BCL-XL ratio showed a striking reversal, deviating significantly from the expected pattern" | **BUILT 2026-08-05** `fig2_priming_ratios.R` | `priming_arm_teb.rds$priming` + `$pair_null`, `collapse_module_ownership.rds$defs$global_rate_fitted` + `$wt_genes` |
| **2H** | "closely tracking *Bbc3* was its known p53-independent activator, *Foxo3*" | **BUILT 2026-08-05** `fig2_departure_from_dose.R` | `collapse_module_ownership.rds$collapse_genes` + `$defs`, `interaction_results.rds`, `combined_df_annotated_raw.rds` |
| **2I** | "the interaction between MYC and OXPHOS subunit coupling to the PUMA/Bcl-XL ratio is highly significant (p = 0.0052), whereas no such link exists for other redox and metabolic axes" — **author's call 2026-08-04: this gets a panel, and it is 2I, not a supplementary** | **BUILT 2026-08-05** `fig2_oxphos_puma_coupling.R` | `substrate_specificity_tradeoff.rds` (`$tradeoff`, `$tradeoff_perm`, `$ambient`), `priming_arm_teb.rds` (`$axis_scores`, `$purity`), `gsva_scores.rds$expr_mat` |

### Supplementary 2

| slot | supports | script |
|---|---|---|
| **S2A** | "these alterations were not a consequence of decreased MYC expression at 12W" | bench |
| **S2B** | "the WT temporal program operated independently of the MYC-driven reallocation, showing a negligible correlation across the mitochondrial and the whole transcriptome" | **BUILT 2026-08-04** `figS2_reallocation_independence.R` |
| **S2C** | "the overall apoptotic priming remains stable in the WT timeline" | **BUILT 2026-08-05** `figS2_priming_balance.R` |
| **S2D** | "no changes in the transcriptome of other known PUMA inducers were observed" — the roster already exists as `priming_arm_teb.rds$exclusions$puma_inputs` (12 genes: `E2f1`, `Trp73`, `Atf4`, `Ddit3` …) | **BUILT 2026-08-05** `figS2_puma_inducers.R` |

**Fourteen new panel scripts, built one at a time in citation order:**
**ALL FOURTEEN BUILT** (2026-08-05): `1E → 1F → 1G → 1H → S1E → S1F → 2E → 2F → S2B → S2C → 2G → 2H → S2D → 2I`. Only the author's bench panels remain (1A, S1D, S2A, 2A–2D).
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
