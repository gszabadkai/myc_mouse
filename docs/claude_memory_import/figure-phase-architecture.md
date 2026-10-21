---
name: figure-phase-architecture
description: "Block B publication-figure phase: the figures/panels/ layer, paper/analysis_record.qmd (the frozen record + its capture/!expr mechanism), and the git structure for the reduction (tag block-b-full -> branch paper-final, retire via retired/)"
metadata:
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-08-14T08:58:55.818Z
---

**Started 2026-07-22** (branch `paper-figures`, commit `5940c11`). Moving from exploratory
drafts (scripts 32-39 each write draft panels) to the ACTUAL paper: **2 main + 2 supplementary
figures, ~6-8 panels**, built one panel at a time, plus a reflective Quarto document.

**Two locked decisions (author):**
1. **Figures first, thin Quarto stub now** (not a parallel Quarto build) -- fill the doc as panels land.
2. **Build panel-by-panel, DEFER the full 2+2 partition** -- start with the strongest result and let
   the spine emerge. (Two candidate spines were offered: "Content->Progression" title-aligned, or
   "Program->Methods" foregrounding the ceiling; author picked neither yet -- deferred.)

**Architecture (NEW layer, separate from analysis scripts):**
- **`figures/theme_myc.R`** -- reusable backbone: `theme_myc()` (theme_classic base), genotype palette
  `geno_cols` (neg=#4575B4 WT blue / pos=#D73027 Myc+ red), `group_labels`, `fig_w` (Nature 89/120/183mm),
  `save_panel()` = ggsave via `cairo_pdf` (Type-42 fonts, mm). Does NOT source 00_setup (light layer).
- **`figures/fig01_mito_content.R`** -- reads `results/mito_content_proxies.rds` (+ script-33 rds for the
  adjusted stat), writes `outputs/figures/fig01_mito_content.pdf`. Save GUARDED by
  `options(myc.fig.nosave=TRUE)` so the Quarto doc displays `p` without re-exporting.
- **`paper/myc_mito.qmd` + `paper/_quarto.yml`** -- Quarto stub. knitr engine, `execute-dir: project`,
  `freeze: auto`, self-contained HTML. Section spine = exploratory arc (content -> global axis/ceiling ->
  attenuation -> death null -> Lesner -> leads); read-only setup chunk; one worked section embeds fig01.
  **`here::i_am("paper/myc_mito.qmd")`** pins root (here mis-anchors to paper/ otherwise). Renders clean
  (quarto 1.9.38). Build artifacts gitignored via `paper/.gitignore`. Render: `quarto render paper/myc_mito.qmd`.
- **THE DRIVE `Icon\r` TRAP (2026-07-27, fixed).** The repo is on a **Google Drive for Desktop mount**
  (`/Users/gs/G/...`), which stamps a zero-byte `Icon\r` (resource-forked) into **every folder it syncs**.
  Quarto walks `paper/myc_mito_files/` with readdir and dies on the stub -- `NotADirectory (os error 20)`
  -- **after every chunk has run**. **A pre-render clean is NOT enough:** Quarto creates that directory
  mid-render and Drive stamps it within seconds, so the killer file is written *inside* the render window
  (one render succeeded, the next failed, identical inputs -- a race, not a state). Fix = `paper/
  clean_icon_files.sh` wired as `pre-render: ./clean_icon_files.sh watch` / `post-render: ... stop`: a
  token-locked background sweeper that runs for the render's duration (no PIDs -> nothing stale to kill).
  Verified by two consecutive renders. **Generalise: any tool here that walks a directory tree can hit
  this**, and hooks take args + run with cwd = project dir (`QUARTO_PROJECT_ROOT` is set for post-render
  only). Supersedes the old "validate with pandoc, quarto render is broken" workaround.
- **RESTART / SESSION RECOVERY (2026-07-28, `docs/2026-07-28_restart_runbook.md`).** **Nothing needs
  re-running after a reboot** -- `results/*.rds` are files, not session state. `scripts/99_session_
  bootstrap.R` (read-only) restores the session: audits the **11 objects the whole figures/ layer
  reads**, loads them into a list `S` (not globals -> cannot shadow a figure script), prints git +
  Icon-stub state. `figures/rebuild_manuscript_figures.R` re-sources the 4 manuscript figures (own env
  each; under `myc.fig.nosave` it also DRAWS to `pdf(NULL)` -- nosave skips `save_panel()`, which is
  where patchwork assembly/fonts actually render, so a dry run without this only tests the data paths).
  `clean_icon_files.sh` gained **`repo` mode** (sweeps the whole repo incl. `.git`) + `-size 0` in the
  match, which is what makes pointing it at `.git` safe. **Two durable traps found:** (1) `quarto
  render` re-executes all 18 figure scripts but writes **no PDFs** (qmd setup sets `myc.fig.nosave=
  TRUE`) -- refreshing `outputs/figures/` is a separate action; (2) **`combined_df_annotated{,_raw}.rds`
  have NO writer in `scripts/`** -- they come from `archive_main_pipeline/02_deseq_interaction_model.R`
  (bare relative paths, needs cwd=root), and `results/` is gitignored, so they exist in ONE place.
  Also: a 9.7 MB `.RData` from 2024-10-01 sits in the root with `RestoreWorkspace: Default`.
- **`00_setup_packages.R`** gained `ggbeeswarm` (quasirandom, n=6/group) + `patchwork` (composition).

**THE FIVE ARMS (fig01 + fig01b, settled 2026-07-23):** mass markers (12) / nuclear MitoCarta (1027)
/ nuclear OXPHOS (155) / **`BIOGENESIS_FULL` (314)** / mtDNA-encoded (13). The biogenesis arm =
`union(MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA` 230, `MITOCARTA_PROTEIN_IMPORT_SORTING_AND_HOMEOSTASIS`
88`)`, **replacing the mitoribosome facet** (83/83 nested inside it). WHY: the claim is OXPHOS vs
GENERAL biogenesis and only **7/314 (2%) are in OXPHOS_NU** -> genuinely disjoint, which a
mitoribosome-vs-OXPHOS comparison never was; and it reproduces the project's own `mito_biogenesis`
axis (`28:382`, used by 35/36) so figure and coupling analyses name ONE thing. **REJECTED: the
`Biogenesis_discrimination` panel** (`CORE_MITO`=PGC1a regulon 300, `MYC_MITO` 119, the `bio_comp` of
19/23/24) -- those are REGULONS and MYC_MITO is Myc targets by construction => circular on a genotype
share axis (cf. script 38 PART C). Caveat noted in-script: 5 of the 12 mass markers (Tomm20/22/40,
Timm23/44) sit inside BIOGENESIS_FULL. Script 32's `panel_roster` gained 3 rows (2 gmt + the derived
union; components retained as sub-arms) -> **author must re-run 32**; rows are purely additive
(shares computed per panel, no loops/perms) so every existing number is unchanged, and script 33
filters by panel name so needs no re-run. fig01b needs NO re-run (builds membership from the GMT).
LFC result: biogenesis is the STRONGEST arm at 6W (94.7% up, 155 sig) and attenuates like OXPHOS
(87.8% vs 87.3% up at 12W, but 14 sig vs 1 -> counts track power/size, not biology).
**Author ran 32 + fig01 2026-07-23; pre-existing arms bit-identical (additive guarantee HELD).
Biogenesis is the ONLY arm significant on genotype at BOTH timepoints: 6W p=0.035, 12W p=0.0015
(vs mass 0.085/0.0028, nuclear 0.11/0.018, OXPHOS 0.12/0.015, mtDNA 0.56/0.48). Share space and LFC
space now RANK THE SAME ARM FIRST -> qmd tension 1 rewritten: part of the "6W marginal" problem was
ARM CHOICE (a sub-arm), not only n=6. Sole significant time bracket = mtDNA within Myc+ p=0.0079,
batch-confounded -> described, not claimed.**

**Panel 1 = the load-bearing result** ([[mito-content-bounded]]), revised 2026-07-23 on author review
(commit `c7cf91f`, NOT pushed): Myc raises the nuclear mito arms vs mtDNA flat. **5 facets** (see
above), free_y, quasirandom
points n=6. **x-axis dropped**; 4 groups carry a bottom colour legend via `group_cols` (hue=genotype,
lightness=age; RdBu endpoints). **p-value brackets** for 4 contrasts (genotype@6W / @12W, 6W-vs-12W in
WT / in Myc+), recomputed by script 32's OWN log2-share simple-effect method so they reconcile;
**significant p<0.05 in red**. Key finding surfaced: per-timepoint genotype is **marginal at 6W (p
0.05-0.12) but firm at 12W (p 0.003-0.018)** -- disaggregating costs power vs the pooled +27%/p=0.003
headline (author kept the per-TP split; batch-confound of the time brackets handled in text, not greyed).
Brackets compressed to top strip -> data keeps ~2/3 of axis. One 12W WT mtDNA sample (67%) drawn
off-scale (cap 40). **Y-axis = "fraction of the nuclear transcriptome (%)"** -- NOT "non-mitochondrial":
`share_nomt` drops ONLY the 13 mtDNA genes (`32:253-256`), leaving all nuclear MitoCarta in the
denominator. **Panel 1 not yet signed off / not pushed.**

**Panel 1b = LFC companion** (`figures/fig01b_mito_content_lfc.R`, commit `e628a2f`, NOT pushed): per-gene
RAW log2FC for the SAME four comparisons as fig01's brackets (from `interaction_results.rds` = the 4 raw
DESeqResults `myc_6W_raw`/`myc_12W_raw`/`timepoint_neg_raw`/`timepoint_pos_raw`; LFC agreement w/ combined
cor=1). Rebuts fig01's compositional critique: nuclear arm genes rise absolutely + coordinately (Myc@6W
83-98% up), mtDNA NOT Myc-induced (0-15% up), rises only over time. **Now a patchwork COMPOSITE** (commit `9090cde`): TOP
track = per arm x contrast 100% stacked direction bar (up `#D6604D` / down `#4393C3` -- MUTED RdBu mid
tones, author request -- by % of set genes) with the padj<0.05 count inside each segment; count label
BLANKED when nsig==0 or its segment is <6% wide (a 0 is noise; a count on a sliver spills outside the
fill). BOTTOM track = the LFC boxplots (text removed, arm labels on top
track), columns aligned. x-axis dropped -> bottom legend (Myc@6W/Myc@12W/WT-time/Myc+-time). Needs
`patchwork`. **NESTING kept as-is (author 2026-07-23):** Nuclear MitoCarta is an UMBRELLA fully containing
OXPHOS_NU (155/155) + BIOGENESIS_FULL (314/314) + mass (12/12); mtDNA disjoint. But the two CLAIM facets
(OXPHOS vs biogenesis) are effectively disjoint of each other, 7 shared genes. Facets double-count; a future
"Other nuclear" (877) partition is deferred (would also touch fig01). Documented in-script. Arm
memberships = GMT
(`mammary_mito_myc_metab_v1_mouse.gmt`) + mass roster, symbol-matched to mgi_symbol (reproduces script 32
counts). Needs `DESeq2` to coerce the DESeqResults. **TWO CROSS-PANEL TENSIONS (noted in qmd sec-content +
sec-attenuation):** (1) 6W LFC-strong (369 sig) vs fig01 share-marginal (p~0.09) = per-gene induction big
but n=6 aggregate noisy; (2) 12W counts undersell the coordinate shift (OXPHOS 87% up but 1 sig) = the
%-up line (power-independent) fixes this, previews attenuation. Integrated: 6W strong-per-gene/noisy-
aggregate; 12W attenuated-per-gene/consolidated-share. Competitive NES deferred. Fig01/fig01b = the natural
main+supp pair (share vs LFC). Parked idea: 2x2 contrast-square diagram (4 contrasts as sides).

**SUPPLEMENTARY SURVEY (figS1/figS1b, built 2026-07-24).** Answers the author's question: in the
WT 6->12W timeline the TOTAL nuclear MitoCarta + biogenesis are ~0 LFC but OXPHOS is mostly DOWN
-- so what rises to hold the total flat? Surveys **all main MitoCarta groups** the same two ways as
fig01/fig01b. **Roster = the 7 canonical MitoPathways top-level categories with METABOLISM split
into its 9 depth-2 children + mtDNA-encoded ref = 16 groups.** Hierarchy verified against
`data/Mouse.MitoCarta3.0.xls` Sheet 4 `MitoPathway Hierarchy` (7 depth-1: OXPHOS, Metabolism,
Central dogma, Import/sorting/homeostasis, Dynamics & surveillance, Signaling, Small-molecule
transport; Metabolism's 9 children: amino-acid, carbohydrate, lipid, nucleotide, vitamin,
metals&cofactors, detox, sulfur, electron-carriers). All GMT sets already mt-* free ->
`MITOCARTA_OXPHOS` == `MITOCARTA_OXPHOS_NU` (155, NO mt-* genes; confirmed). **OXPHOS is shown
as its two halves: Nuclear OXPHOS (OXPHOS_NU 155) + mtDNA OXPHOS (OXPHOS_MT 13, IDENTICAL genes
to MITOCARTA_MTDNA_ENCODED -- the 13 mtDNA genes ARE the 13 mtDNA OXPHOS subunits)** -> so the
nuclear category carries no mtDNA bias and there is NO separate "mtDNA reference" facet (author
correction 2026-07-24). `share_nomt` throughout. The two figures now have **7 top-level + 9
metabolism = 16 groups** (was 6+9+ref). mtDNA OXPHOS moves OPPOSITE to nuclear OXPHOS over the WT
timeline (nuclear-down/mtDNA-up = the Lesner mitonuclear discordance, on the confounded time axis).
- **`figures/figS1_mitocarta_survey_share.R`** -- share_nomt boxplots. **2 patchwork blocks with a
  DIVIDER** (bold block subtitle): Top-level (7, OXPHOS pair leads) / Metabolism level-2 (9), each
  facet_wrap ncol=3, one collected legend, y-label on top block; needs patchwork. **p<0.05-ONLY
  brackets** (author request): fig01's 4 comparisons via script 32's log2-share simple-effect
  method (pval_for), drawn only where p<0.05, genotype=red / time=grey; reconciles fig01 (Nuclear
  OXPHOS geno 0.015, mtDNA OXPHOS time 0.0079). Rendered for REAL (32 re-run done 2026-07-24).
- **`figures/figS1b_mitocarta_survey_lfc.R`** -- reads `interaction_results.rds` + GMT directly, **no
  re-run**, verified in scratchpad. The fig01b bar-atop-box composite does NOT scale to 16 groups
  (aligns only for 1 facet row) -> redesigned as **A: median-LFC heatmap** (16 groups x 4 contrasts,
  diverging fill, sig-count text) over **B: the LFC boxplot grid**. Also writes the **divergence
  table** `outputs/figures/figS1b_diverse_genes.csv` -- per group the top-15 genes by LFC SPREAD
  (max-min) across the 4 contrasts (author's "most diverse", = context-dependent movers).
- **RESULT (the answer):** WT timeline OXPHOS 38% up (6 sig DOWN) while substrate-metabolism groups
  RISE (sulfur 80/carb 77/amino 76/lipid 68% up; metals, detox, transport, e-carriers >55%). Flat
  total = REALLOCATION not stasis. Myc@6W induces every group (red), fades by 12W (attenuation).
  Time columns batch-confounded -> described not claimed; genotype cols clean; % up is the
  power-independent read.
- **Script 32 gained 12 additive gmt roster rows** (PART 2b) for the share panel; author re-runs 32
  (pre-existing numbers bit-identical). fig01/fig01b untouched. Takes a `figS1` slot -- the `fig02`
  reserved for the reallocation main panel is unaffected.

**COMPLEX BREAKDOWN (figS2/figS2b, built 2026-07-24).** Granular supplementary drilling into nuclear
OXPHOS: the 5 respiratory complexes (MITOCARTA_COMPLEX_I..V; CI 57/CII 8/CIII 15/CIV 47/CV 24; all
nuclear-encoded, mt-* removed). figS2 = share boxplots (1 row, p<0.05 brackets); figS2b = median-LFC
heatmap (5x4) over LFC boxplot grid + divergence table (`figS2b_diverse_genes.csv`).
- **RESULT = UNIFORMITY not selectivity:** Myc raises all 5 complexes coordinately at 6W (share
  genotype-sig for CI/CIII/CIV/CV; CII the small underpowered exception, the one complex w/ no mtDNA
  subunit), all attenuate in step by 12W. No complex singled out -> respiratory chain driven as one
  block (coordinate biogenesis, consistent w/ sec-content).
- **CRITICAL MAPPING FIX (reusable):** MitoCarta complex sets carry OLD ATP-synthase names
  (Atp5a1/b/c1/...); the count matrix/combined_df use CURRENT names (Atp5f1a.., Atp5mc*, Atp5pd/pf/po).
  A current-SYMBOL match silently drops **17 of Complex V's 24 genes** (CV showed 7/24!). Fix: map
  MitoCarta Symbol -> its OWN EnsemblGeneID (**Sheet 2** of `data/Mouse.MitoCarta3.0.xls`, split on
  "|") -> intersect the DE/count universe. Recovers CV 24/24. **This also means script 32's ens_of
  (symbol match) undercounts OXPHOS_NU 134/155 etc. -- a KNOWN pre-existing loss in fig01/figS1; NOT
  fixed there (would change signed-off numbers).** figS2/figS2b are therefore **SELF-CONTAINED**
  (compute shares from count_matrix directly; read mito_content_proxies only for sample metadata) ->
  no script 32 change, NO re-run, both verified in scratchpad. Needs readxl.

**REALLOCATION -- MAIN fig02 + SUPP figS3 (built 2026-07-24).**
**MAIN = `figures/fig02_reallocation.R`** -- the "reprioritisation EXISTS" panel, an ORDERED DOT
CHART (author corrected TWICE: NOT per-sample scatter, NOT tier-fan bars -- a ranked Cleveland dot
chart). Reads `mp$mitopps_pairwise` (per pathway x contrast diff/padj; Myc_effect_6W / _12W) +
`pathway_tier1_map`. All ~144 MitoPathways RANKED by their Myc mitoPPS effect (Myc+ - WT, from
`diff`), most demoted bottom -> most promoted top; y = pathway (labelled, 4.2pt, tall ~250mm),
x = effect, one dot/pathway coloured by tier-1, **filled=6W / open=12W joined by a thin line**
(the 6W->12W change). Structured by function: demoted end = dynamics/apoptosis/signalling + mtDNA
OXPHOS; promoted end = biosynthetic metabolism + import (SAM/chaperones); **OXPHOS visibly SPLIT
(subunits high, assembly low) -> ~neutral as a tier**. Metabolism greyed. Tier-MEAN story
(import +0.10 / dynamics -0.11) is cleaner than individual extremes (small metab pathways noisy)
-> subtitle kept descriptive ("structured by function"), not "biogenesis at top". No re-run;
verified in scratch.
**MAIN = `figures/fig02_reallocation_ranked.R`** (renamed from fig02b; the tall one became figS4).
(author: the all-labelled tall version was too detailed; wanted a
horizontal, NO per-pathway labels, only asterisks + colour + labels on the top hierarchy). Same
ranking HORIZONTAL (x = rank demoted->promoted, y = 6W Myc effect S-curve), dots tier-coloured,
`*` on the 33 padj<0.05 pathways, and ggrepel leader-labels on ONLY the **7 Level-1 top-of-
hierarchy pathways** (`pathway_levels$Level=="Pathway_Level1"`; they are themselves scored
composites -> anchor where each top-level category sits: Dynamics most demoted ... Import most
promoted, Metabolism ~0). Compact ~100mm. Ranked DESCENDING (promoted left -> demoted right); 12W shown as small points
joined to 6W by thin connectors (the FADE); anchor labels on a SINGLE y line per strip
(direction="x") -- promoted top row / demoted bottom row, off the data. Both kept: **fig02
(`fig02_reallocation_ranked.R`) = clean MAIN**, **figS4 (`figS4_reallocation_full.R`, was the
tall fig02_reallocation.R) = the fully-labelled tall companion (supplementary)**; figS3 =
two-rulers detail. needs ggrepel.
**SYNTHESIS = `figures/figS5_reallocation_scatter.R`** (author idea 2026-07-24): the two rulers
as ONE scatter -- x = content (DESeq set-avg log2FC), y = priority (mitoPPS effect), per pathway,
6W, tier-coloured, ringed = padj<0.05, priority-extreme pathways labelled, + **MARGINAL DENSITIES**
(patchwork). **NO trend line** (author: OLS fit misleading since mitoPPS is ~content-minus-mean,
so a slope is partly built-in). The message the marginals carry (author's point): CONTENT is
one-sided-positive (~94% up, top hump right of 0) while PRIORITY is WIDE both ways (~56% up, right
density straddles 0) -> the wide priority spread on a near-uniform content rise IS the reallocation.
FINDINGS
worth the paper: (1) r=0.76, cloud in content-UP half -> reprioritisation is a SECOND-ORDER skew
on a broad rise, not the main axis; (2) reallocation is DIRECTIONAL -- the content-down/priority-up
quadrant is EMPTY (0 pathways; 55 are content-up/priority-down) -> demoted pathways aren't shut
off, they rise slower than biogenesis and lose relative share; (3) residuals name candidates:
ABOVE = chaperones/immune/anaplerotic metab (creatine/pyruvate/serine); BELOW = OXPHOS assembly,
cAMP-PKA, dynamics/apoptosis. CAVEAT: mitoPPS is a ratio so part of r is built in -> the SIGNAL is
the residual, not the trend. Reads mitopps_scores + interaction_results via reconciler, no re-run.
qmd sec-reallocation gained the synthesis paragraph + figS5.
**SUPP = `figures/figS3_reallocation_two_rulers.R`** (renamed from fig02_..., author designated it
supplementary) -- the pathway-level detail below.

The two-rulers heatmap across **all ~148 MitoPathways**
(tier-blocked, facet_grid(tier~.)): LEFT "Content (absolute)" = DESeq set-average raw log2FC
(via the reconciler on `mp$gene_to_pathway` genes + interaction_results; * = one-sample t vs 0,
BH<0.05); RIGHT "Priority (mitoPPS)" = `mp$mitopps_pairwise$diff` (* = padj<0.05). Two genotype
contrasts (6W/12W); native units + separate legends so the discordance = a COLOUR FLIP. Reads
`mitopps_scores.rds` + `interaction_results.rds`, no re-run (mitoPPS regenerated by the author's
post-reconciliation re-run). Verified in scratch.
- **RESULT:** content positive almost everywhere (broad Myc induction, attenuating 12W), but
  priority reshuffled. **OXPHOS = the exception: subunits PROMOTED (CIII-sub mitoPPS +0.20,
  CV-sub +0.17), assembly factors DEMOTED (CII −0.26, CIII −0.14 padj 0.022, respirasome −0.11)
  + mtDNA arm −0.34 -- all with POSITIVE content.** Subunit-over-assembly split only visible at
  all-148 granularity. Dynamics/apoptosis + signalling demoted; import chaperones/SAM up on both.
- **Framing:** clean genotype axis, ratio cancels common-mode (ambient 0.42, script 38), but
  DESCRIPTIVE/EXPLORATORY (n=6; few tiles clear padj<0.05 -> read the pattern). raw LFCs
  (CLAUDE.md). Tall ~285mm portrait; a compact main-figure summary can be extracted later.
- **`raw_pairwise` (absolute) is UNUSABLE cross-pathway** (raw normalized-count scale, mtDNA-OXPHOS
  −19560 dominates) -> the absolute ruler MUST be DESeq set-LFC. qmd sec-reallocation embeds it.

**TEMPORAL reallocation figS6+figS7 (built 2026-07-24).** mitoPPS change over the TIMELINE
(within-genotype 6->12W): WT (`Temporal_Myc-`) and Myc+ (`Temporal_Myc+`) from `mp$mitopps_pairwise`.
Both mitoPPS-ONLY (read only mitopps_scores; no content/reconciler; no re-run).
- **`figS6_reallocation_temporal_ranked.R`** (author: rank by WT, dropline to Myc+): compact
  horizontal, pathways ranked descending by the WT change (filled dot), open Myc+ dot + thin
  dropline = the Myc-specific bend; 7 Level-1 labels in top/bottom strips (fig02 treatment). NO
  asterisks (temporal padj batch-confounded). Point: the WT BACKGROUND reprioritises over the
  timeline on its own (developmental), distinct from the Myc genotype effect.
- **`figS7_reallocation_temporal_divergence.R`**: x=WT-time, y=Myc+-time mitoPPS, coord_equal,
  **y=x dashed = developmental/shared**; off-diagonal = Myc bends it. 7 Level-1 anchors labelled.
  cor(WT,Myc+)=**0.75**, 70% same direction -> mostly shared development; the biogenesis/OXPHOS arms
  sit BELOW y=x (Myc fades them = the attenuation in mitoPPS space).
- **PRIMING negative (verified):** the apoptosis/mitophagy/autophagy pathways have tiny WT-vs-Myc+
  interaction (median |inter| 0.015 vs 0.047 all-pathway; Apoptosis-PRO ~0) -> they track the
  developmental background, **no Myc-specific temporal priming signal in the mitoPPS lens**
  (consistent w/ [[death-priming-change-untested]]).
- **BATCH-CONFOUND framing:** time axis = batch, so individual WT/Myc+ temporal contrasts are
  DESCRIPTIVE; the clean quantity is the divergence (dropline / y=x deviation) where batch cancels.
  Both SUPPLEMENTARY, exploratory (n=6). qmd sec-reallocation-time (callout-warning on the confound).

**BACKGROUND-vs-MYC fig03 + script 40 (built 2026-07-24, `sec-attenuation` now WRITTEN).** Answers
"is Myc's mito effect shaped by a background change?" -- see [[background-vs-myc-is-additive]] for
the result and the artifact ledger. **`scripts/40_background_vs_myc_decomposition.R`** (author-run,
writes `results/background_vs_myc.rds`): PART A three-ruler pathway table (identity asserted on both
rulers); B the two regressions + TWO nulls (primary = within-decile gene-LABEL SHUFFLE preserving
set size, expression AND pathway overlap -- MitoPathways nest, so overlap alone creates cross-pathway
correlation; secondary = decile-matched resampling); C the **sample-split control**; D geometry +
null; E ruler reconciliation (4 quantifiers); F the 7x4 state table for the figure; G the artifact
ledger. **`figures/fig03_background_vs_myc.R`** = 4 panels (A state-space arrows per tier / B rescale
scatter / C shared-vector scatter / D convergence lollipop), patchwork, 183x175mm.
- **Traps hit and fixed while building** (all reusable): (i) `tibble()` MASKS an outer object of the
  same name as a column being created -- bind `x <- ruler[[v]]` *before* the `tibble(ruler = ...)`
  call; (ii) `dplyr::summarise()` sees earlier results, so `sd(x)` after `x = mean(x)` returns NA --
  order matters; (iii) `conv_pct` is a ratio of noisy quantities -> per-split sd reaches 4 digits;
  aggregate as a **ratio of means**, and read `frac_wt_toward` as primary; (iv) the synthetic
  **mtDNA-encoded pathway is a 4-6.5 SD outlier** on every temporal contrast and alone moves the
  shared-vector slope 0.75->0.84 -> flagged `is_mtdna`, excluded from every fit/null/geometry;
  (v) the same 13 genes inside the OXPHOS **tier** made its content share an mt% readout (sd 23 on a
  mean of 10) -> `setdiff(tier_genes, mt_ens)`; (vi) ggplot draws the arrow key TWICE unless the
  `colour` and `linetype` scales share name+labels+order (no `override.aes` on one of them).

**SUBSTRATE fig04 + DEATH fig05, and the NUMBERING CONVENTION (built 2026-07-27).** The argument order
was decided as **substrate first**, which promoted script 43 PART A to the opening figure. **Build
numbers stay stable; the manuscript order lives in a manifest table at the head of `paper/myc_mito.qmd`
and in the section order** -- renaming 14 scripts to encode an argument order would break every commit,
dated doc and `outputs/` reference for no analytical gain (author's decision). Mapping: **fig04 = MS
Fig 1**, fig01/01b = Fig 2, fig02 = Fig 3, fig03/figS8 = Fig 4, **fig05 = MS Fig 5**.
- **`figures/fig04_substrate_specificity.R`** (43 `$comparator`/`$wt_null`/`$paired_null`/
  `$comparator_priority`/`$buffer` + 44 `$wt_genes`): A ranked arms w/ matched-null percentile
  (OXPHOS subunits -0.255 0th pct, PROLIF pooled -0.047 **1.3rd**, assembly +0.001 50th) / B the
  PAIRED null on the difference / C content-vs-priority / D the 32 apoptosis + 5 brake transcripts,
  1 mover and it RISES. **Design point: the percentile must be ENCODED** -- PROLIF pooled is also
  beyond its own null, so "proliferation doesn't move" would overstate it; the claim is comparative.
  Guards: field-existence + assert `$defs$oxphos_wt_reference` (-0.2548) so a stale rds fails at
  render.
- **`figures/fig05_death_arm.R`** (42 `$machinery`/`$priming`/`$pair_null` + 44 `$collapse_genes`):
  A the 6W roster (Htra2 +1.49) / B retention vs the x0.55 line, **open = 6W effect n.s. so its
  retention is a ratio of noise** / C the pro:anti ratios vs the CONDITIONAL matched-pair null
  (its own median 0.551 reconstructs x0.55 = the positive control) / D the collapse distribution
  over 8774 genes w/ Bbc3 0.51 pct, Bcl2l11 91.59, Myc 99.99. **Honesty encoding that matters:
  `Bmf:Bcl2l1` is FURTHER out than Bbc3 (97.5th, p 0.025) but has no 6W effect (p6 0.40)** -> the
  fill encodes the precondition rather than the caption carrying it.
- `sec-substrate` (new, before `sec-content`) + `sec-death` (rewritten from the nine-step order) are
  now written; qmd title updated to the decided one. Remaining stubs: `sec-axis`, `sec-lesner`.

**THE PAPER FIGURES, 2026-07-27 — and the reason fig01-fig05 are NOT them.** The author ran fig04/fig05,
agreed they carry the conclusions, then rejected them for the paper: *"the argumentation became too
complex ... readers of a leading journal need more clear cut explanations ... the story has to be more
simple and use less jargon."* Root cause: the panels were designed to answer *"does this survive its
null?"*, so they lead with hedges. **Re-scoped to 2 main figures (<=6 panels) + 2 supplementary**, ordered
by ¶4 of the drafted Introduction, which is **MYC-FIRST** (this SUPERSEDES the substrate-first figure
order in the death-narrative §9 — the nine steps stay the argument, the figures partition by actor).
- **`figures/figure1_myc_mitochondrion.R`** (32 + 08 + 42): A content share 4 arms w/ genotype brackets
  only / B all 144 MitoPathways ranked by mitoPPS effect / C content-vs-priority, the **content-down/
  priority-up quadrant is EMPTY** (0 pathways; 55 the other way) / D the roster at 6W (HTRA2 +1.49, BAX
  +0.48, BCL-xL -0.42).
- **`figures/figure2_developmental_window.R`** (43 + 44 + 42): A the withdrawal + null pct / B paired
  null / C both rulers / D machinery intact / E PUMA reverses vs the x0.55 line, w/ the rest of the
  roster as a rug / **F THE THREE PROPOSITIONS** (see below).
- **PANEL F is the new idea and it came from the author's question** *"why is it important that the
  gland withdraws the respiratory arm — does it contribute to the lost sensitivity? can it be shown
  after?"*. Answer: structurally load-bearing (without it the mouse is a pure dose story and MYC-ER is
  unexplainable), mechanistically **a hypothesis this dataset cannot promote**. So the panel presents
  **H1 respiration gates PUMA priming** (91.9th pct, redox control 50.9th = dead centre) / **H2 the TEB
  transition moves both** (90.8th; the honest rival, not separable at n=24) / **H3 respiration costs
  proliferation UNLESS death is blocked** (6.3rd/9.2nd pct) — each annotated with the experiment that
  adjudicates it, and **H3's sign flip drawn**: mouse (death intact) respiration costs growth; MYAZ +
  Bcl-xL (death blocked) PGC-1a buys it. **H3 had been mis-filed as a control; it is the hand-off.**
  Plotted on the PERMUTATION PERCENTILE because the two source tables (42 `perm_axis`, 43
  `tradeoff_perm`) are on different coefficient scales. The script ASSERTS the two independent
  implementations of H1 still agree (91.94 vs 91.72) before drawing.
- **qmd restructured**: `sec-substrate` + `sec-death` MERGED into `# Development opens the window
  {#sec-window}`; `sec-content` renamed "MYC builds a death-competent mitochondrion" and carries
  figure1; bounds moved into `collapse="true"` callouts so the main read is the claim. Old fig01-fig05
  chunks stay as the exploratory layer. **¶4 of the Introduction rewritten** to name PUMA/PGC-1a/BAX/
  HTRA2/BCL-xL, with PGC-1a **only in cell clauses** (not separable from MYC here, r 0.93).
- **`figures/figureS1_compartment_detail.R`** (Supp 1, controls for Fig 1): A per-gene DIRECTION bars
  across the 4 contrasts = the answer to "a share can rise because everything else falls" / B the
  16-group median-LFC heatmap, WT column = what rises while OXPHOS falls (14 of 16 rise) / C the five
  complexes. **Guard that earned itself immediately:** two of panel A's arms (`MASS_MARKERS_NOCHAP`,
  `BIOGENESIS_FULL`) are NOT GMT sets — they are script-32 constructs (mass roster minus chaperones;
  the central-dogma + import UNION) — so the tabulator takes a named list of SYMBOL VECTORS, not set
  names. Separate assert: Complex V must resolve to >=20 genes or the reconciler has regressed.
- **`figures/figureS2_controls.R`** (Supp 2, controls for Fig 2): A driver-vs-output (drivers widen
  +0.17/+0.32/+0.25, output narrows -0.36 => attenuation is DOWNSTREAM of dose; asserted in-script)
  / B the x0.55 rescale + its 100th-pct null / C shared slope (82nd, not beyond) vs Myc-specific
  intercept (0th) / D the genome-wide collapse. The WT-convergence control stays in fig03.
- Both supplementaries embedded in the qmd (`sec-mito-survey`, `sec-attenuation`); manifest updated.

**THE APPROACH TURNS AROUND, 2026-07-29 (commit `ff15bb6`).** The four assembled figures are "a good
starting point but still missing points from the narrative", so the author now writes the Results
**in the Google Doc `MK_Myc_Paper`, tab "Gyorgy Writing"**, and panels are built **one at a time to
support the sentences as they are written** -- see [[panels-follow-the-narrative]]. New layer
**`figures/panels/`** (own PDFs in `outputs/figures/panels/`), on the same `paper-figures` branch;
`theme_myc.R` and the four assembled figures are UNTOUCHED and still rebuild.
- **THE NATURE RULE, made mechanical:** no explanatory information on the page. `theme_panel()` in
  `figures/panels/_panel_common.R` BLANKS title/subtitle/caption, so a panel cannot acquire prose by
  accident (the Gen-2 scripts bake `title = "A   MYC raises..."` into the plot -- that idiom is now
  forbidden for new panels). Each script carries a **`panel_legend(slot, what, detail, bounds, source)`**
  block; `rebuild_panels.R` collects them into `outputs/figures/panels/legends.md` = the raw material
  for the figure legends. `_panel_common.R` also declares the three copy-pasted palettes ONCE and gives
  method routing a purple ramp (never spend the genotype blue/red on a non-genotype key).
- **Figure numbers live in `figures/panels/PANELS.md`, not in filenames** -- the Doc's Figure 1 is the
  GLOBAL OVERVIEW, while the repo's `figure1_myc_mitochondrion.R` is the mito-content figure, which sits
  later in the Doc's numbering. A slot change is a table edit, not a rename.
- **Paragraph 1 maps to:** 1A = the author's whole-mount micrograph (NOT built here); **S1A** =
  the 986-set library by category x quantifier (BUILT); **S1B** = the 2x2 design + five contrasts,
  batch=timepoint encoded solid-vs-dashed (BUILT); **1B** = GSVA cell-state composition (BMYO/LASP/LHS
  flat) + the TEB->ductal shift; **1C** = endogenous vs transgene on MYC/prolif/TEB + BMYO->LHS.
- **1B/1C ARE BLOCKED on a re-run:** `dev_program_myc_integration.rds` (Jul 15) and
  `myc_endogenous_amplification.rds` (Jul 8) predate the Jul-24 reconciliation that rebuilt
  `gsva_scores.rds`; script 26 also reads `gsva_overview.rds` (Jul 6) at `26:68`. Chain = **17 -> 26 ->
  27**, author-run; 17 is slow (`mroast` nrot=9999 x 902 sets). `require_fresher_than()` stops a panel
  rather than draw a stale number.
- **Two text fixes surfaced:** "~900 genesets" is four defensible numbers (**986** library / **902**
  GSVA-scored / **817** fGSEA-eligible / **884** in the coupling analyses); "~3K DEGs" at 6W is **2777
  at FDR 10%**, 1967 at FDR 5%, and the `interaction` contrast has **0**.
- **The honesty problem inside 1B:** "a clear TEB-ductal shift" fails the sample-level test
  (`teb_ductal` WT 6->12W = **-0.320, p=0.243**, n=6). What licenses "clear" is the matched-null
  effect size in `docs/2026-07-26...` §0.7 (**-0.422, 0th percentile**). Panel shows the honest n=6
  distribution; the legend block carries the statistic that supports the word.

**STATE 2026-08-09 (the `figures/panels/` layer, 25 scripts).** Round 1 delivered **all fourteen
commissioned panels** (1E-1H, S1E, S1F, 2E-2I, S2B-S2D), each written against one written sentence,
run by the author, signed off one at a time, and recorded in `figures/panels/PANELS.md` — which is
now the live slot map and also carries the running list of sentences the data contradicts. Round 2
added **four alternatives** on the two-timeline plane ([[two-timeline-grammar]],
[[main-figure-panels-must-carry-the-argument]]): `fig2_two_timelines_mito_alt.R`,
`fig2_death_two_timelines_alt.R`, `fig2_puma_chain_alt.R`, `figS2_p53_arm_alt.R`, commit `37c3e61`.
**`panels_to_pdf.R` now writes THREE circulation documents** — `Figure1_and_S1_panels.pdf` (12 pp),
`Figure2_and_S2_panels.pdf` (8 pp), `Figure2_alternatives.pdf` (4 pp) — the first two matched with
`not = "\\(alt\\)"` so adding alternatives cannot disturb what has already circulated. **Verify a
panel PDF with `qlmanage -t -s 2000 -o <dir>`, NOT `sips`** — `sips` renders small filled circles as
crosses, which invites "fixing" a shape that was never wrong. Still outstanding: the author's own
bench panels (1A, S1D, S2A, 2A-2D) and figure assembly.

**How to apply:** every new panel = a `figures/figNN_*.R` sourcing `theme_myc.R`, reading a `results/*.rds`,
guarded save. **patchwork trap hit twice now:** `&  guides(...)` applies to EVERY panel in the band and
resurrects guides that individual panels set to `guide = "none"` — put `guides()` on the one plot that
should emit the legend, and use `&` only for `theme()`. Also: never build a gene-pair label by
`sub("Bcl2l1", "BCL-xL", pair)` — it turns `Bcl2l11` (BIM) into "BCL-xL1"; build from the pro/anti columns. Verify via scratchpad-redirected render (`options(myc.fig.nosave=TRUE)` + ggsave into the
scratchpad, then `sips -s format png` and READ the image -- clipped titles/subtitles and collided repel
labels only show up visually), never touch `outputs/`/`results/` (Option A: author runs canonical +
pushes). Scripts 34/35 were never run -> their draft panel dirs are absent; re-source before any
death-null / ceiling-correction panel.

**BLOCK B FROZEN, 2026-08-14 (commit `0cb07c2`, tag `block-b-full`, branch `paper-final` created).**
Collaborators said the paper does not need this detail, so the three written Results sections are
being cut to **half**. Before the cut the full state was recorded, because most of what goes is not
wrong -- it is the evidence, bounds and negative controls that LICENSE the surviving sentences.
- **The narrative source moved**: Google Doc **`11.8.26_MK_Myc_Paper`, tab `CURRENT VERSION`**
  (not `MK_Myc_Paper` / "Gyorgy Writing"). Working title *"Mitochondrial content sets an apoptotic
  threshold for Breast cancer transformation"*. The doc export starts `# CURRENT VERSION` and the
  three bioinformatics sections sit at char offsets ~7.3k / 12.1k / 14.8k of the tab.
- **`paper/analysis_record.qmd`** = the record. Ten chapters: how to read it; the measurement + the
  FOUR standing constraints (no sort / batch=timepoint / n=6 / the 0.8 correlation ceiling); the
  three sections **verbatim** with all 25 panels; **what did not survive, grouped by WHY** (circular
  by construction, at the ceiling or untestable, the wrong null, shared-baseline artifact,
  confounded beyond repair, fitted-and-never-read, doesn't transfer); a 26-row sentence-vs-data
  table; open questions each paired with its adjudicating experiment; a MARKED recommendation for a
  halved paper; generated provenance. Renders ~5 min cold, 5.2 MB self-contained HTML.
- **THE MECHANISM, reusable.** Source each panel with `options(myc.fig.capture = TRUE)`:
  `save_panel_p()` then records the DESIGNED size into `.myc_panel_sizes` in `globalenv()` and
  writes nothing, while `p` and `LEGEND` are still built -- **one pass gives figure + legend +
  size**. Chunk options read it via **`#| fig-width: !expr panel_w("slug")`** (knitr evaluates
  `!expr` BEFORE the chunk body, so cache the env). Legends render through the project's own
  `legend_md()` inside `::: {.callout-note collapse="true"}`. **Nothing is transcribed**, so the
  document cannot drift from the panels. Document ASSERTS its own completeness (25 slugs once each,
  every `(alt)` contiguous with its original, sizes == registry) -> a render that completes passed.
- **Three render traps.** (1) **PNG at 300 dpi, not SVG** -- `svglite` and `ragg` are NOT installed,
  the only vector route is cairo, and cairo is what mangles inter-word spacing at 6 pt; rasterising
  puts it below a pixel (~146 KB/panel). (2) **No bare `@`-token in a legend string** -- pandoc read
  `figS1_design_contrasts`'s "paper/myc_mito.qmd @sec-overview" as a cross-reference. (3)
  **`system2()` pastes args into a shell string** -> `--format=%h %ad` splits and git exits 128;
  `shQuote()` it.
- **GIT STRUCTURE for the reduction (decided, executed).** The full version is a NAME, not a memory:
  commit on `paper-figures` -> annotated tag **`block-b-full`** -> branch **`paper-final`** off the
  tag (created, NOT switched to; author pushes first). **Retire panels with
  `git mv figures/panels/<x>.R figures/panels/retired/`** -- the glob in `rebuild_panels.R:44` and
  `panels_to_pdf.R:68` is `list.files(pattern="^fig.*\\.R$")` with **no `recursive=TRUE`**, so a
  subdirectory drops out of both runners with NO code change and `git log --follow` still works.
  Then delete the stranded PDFs in `outputs/figures/panels/` (gitignored housekeeping).
- **`paper/myc_mito.qmd` -> `paper/archive/`** (recorded in `docs/branch_manifest.md` **section 7**,
  which is now the log for post-manifest renames). It walks the superseded `figures/fig01`-`fig05`
  layer; kept for its global-axis / ceiling prose, with a banner saying its own manifest is stale.
  Its build artifacts were DELETED not moved -- a stale `_freeze/` under the old path would silently
  short-circuit a later render. `paper/.gitignore` `/*_files/` -> `*_files/` (any depth).

## BLOCK C OPENED, 2026-08-16 -- `paper-final` is now the working branch

`paper-figures` is **frozen** at `block-b-full` (0cb07c2, pushed) and must not be added to: it is
the named "full version" the reduction is measured against. Work happens on **`paper-final`**,
which was created off that tag and has now had its first eight commits. Not yet pushed.

- **The panel count is no longer 25.** Three state-reading panels landed (28 total), so the
  record document's completeness assertion is **glob-driven and deliberately has no literal count**
  -- the one hard-coded 25 was a comment and is now a note saying why. A render printing
  "28 panels, 4 alternatives, all paired and contiguous" is the pass signal.
- **The retirement mechanism above is still unused and still correct.** Nothing has been retired
  yet; the halving of the text had not happened when Block C opened.
- **Unslotted panels stay out of the circulation PDFs automatically.** `panels_to_pdf.R` groups by
  Figure 1/2 slot, so a panel declaring `slot = "not currently cited"` is listed as "in no
  document" and the three PDFs stay at 12 / 8 / 4 pages. That is the right default for a panel
  built before its sentence exists.
