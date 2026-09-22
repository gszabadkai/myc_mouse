---
date: 2026-09-22
tags: [project/myc_mouse, handoff, main, figure-1, figure-2, foxo3, verification]
status: live handoff -- overwrite in place at the end of each session
relates-to:
  - docs/2026-09-21_two_timeline_verification.md       (READ FIRST -- script 54, the rulings of record, the four rates; section 9 is the 2026-09-22 Foxo3 record)
  - figures/panels/PANELS.md                            (the manifest of record; from "The Figure-1 verification pass" on)
  - scripts/54_two_timeline_verification.R              (writes results/two_timeline_verification.rds)
  - scripts/47_biogenesis_axis_and_the_developmental_oxphos_decline.R (PART I: the FOXO3 target programme)
  - docs/2026-09-10_orthotopic_escape_series.md         (scripts 52 and 53)
  - docs/experimental_cohorts_branch_notes.md           (the retired branch's notes; its cohort rule still holds)
  - docs/2026-09-09_orthotopic_identity_correction.md   (voids the conclusions of scripts 50 and 51)
  - docs/2026-09-02_myc_oxphos_priming_gate_model.md    (the gate model)
  - docs/2026-08-27_human_validation_plan.md            (the parallel arm; untracked)
  - docs/2026-08-02_results_narrative_final_order_D.md  (the narrative the cut is measured against)
---

# Handoff -- restart from here

**2026-09-22: work continues on `main` from `b7253f7`; `experimental-cohorts` is retired and identical to `main`.**

Read this file first. Sections 1 to 3 are current as of **2026-09-22**. Section 3a is the
2026-09-21 evening session. Sections 3b to 3g and 4 to 9 are carried over, and where today's
work supersedes something in them, it is marked in place.

**The author has the read-only report in 3.3, on whether the FOXO3 chain has a timing
problem, and will decide how the paragraph reads.** Nothing has been built from it. Start
there, and wait.

---

## 1. State

| | |
|---|---|
| checkout | `main`, tracking `origin/main`. **It is the only branch**, local and remote |
| HEAD | this handoff's commit, on top of `1a43f04` ("PANELS.md: the second Fig. 2H row no longer lists combined_df_annotated_raw.rds") |
| pushed | everything through `1a43f04`. **This handoff's commit is local only:** push when the author says so |
| branches | every other branch was deleted on 2026-09-22 (3.1). Five are tagged `archive/*`: `analysis-exploratory`, `BlockA-revision-step-by-step`, `new-analysis`, `paper-figures`, `paper-final`. `experimental-cohorts` was not tagged; its tip `20104a3` is an ancestor of `main`. The other tags are `block-a-reviewed`, `block-b-full` and `orthotopic-pre-identity-correction`. All eight tags are on the remote |
| scripts | 49 to 54, all run; nothing new this session. 52 and 53 are the orthotopic escape series (53's PDF is `outputs/orthotopic/53_escape_dose.pdf`); 54 is the Figure-1 verification. **The next number is 55** |
| results | **no object was written this session.** All were fresh by the rule as of 2026-09-21: Task D's re-sources (3b; the note's section 7), and 33 and 42 at 16:29 UTC. Checked today: script 47's object (2026-08-18) postdates `334ebab`. Script 23's object (2026-07-17) postdates its last commit, but predates `gsva_scores.rds` and `mitopps_scores.rds` (07-24) and `gsva_overview.rds` (07-29) |
| panels | still 28 manuscript slots. **Fig. 2H's legend gained a bound** (3.2). **`legends.md` was regenerated** at 2026-09-22 11:44 by a real `rebuild_panels.R` run, 28 of 28 ok (3.4) |
| repo | `/Users/gs/code/myc_mouse`, off Drive since 2026-09-21 (`50a5638`). `.git` holds no Icon stubs; the six working-tree stubs under `data/` are gitignored |
| scratch | none in the repo. This session's computations ran as read-only R in the session; 3.3 carries the code |

**Untracked, and staying that way:**
- the six `.txt` session exports at the repo root;
- `docs/library_reference/gray_chea_mito_tf_shortlist.csv`;
- `docs/library_reference/2026-08-22_consensus_myc_double_hit_thread.md`;
- `docs/2026-08-27_human_validation_plan.md`.

**THE TEXT CUT HAS STILL NOT HAPPENED**, as far as the repo records. Nothing here blocks it
or depends on it.

---

## 2. Do this first

1. **Read 3.3 and wait for the author's ruling on the FOXO3 paragraph.**
   - The manuscript's chain is: respiratory priority falls, Foxo3 follows, PUMA follows.
   - Bbc3's fall is entirely an interaction. The FOXO3 target programme is lower under MYC at
     both ages, with no interaction detected.
   - So either the link is written as transcript-level, or the chain has a gap at the activity
     step. 3.3 says what the numbers allow.
   - Build nothing until the author decides.
2. **Still waiting from 2026-09-21: the seven-gene table** (3a.3). The author's suggestions for
   what to build from it have not come, and nothing has been built.
3. **Open, and the author's call: Fig. 2G (alt)'s restored claim item** (3a.2). Since today's
   rebuild it is in `legends.md`, so it can be reviewed there.
4. **Offered, not done: the isoform caveat in Fig. S2D (alt)'s legend.** The legend carries
   Cdkn2a's low-count caveat, but not the fact that a gene-level count cannot separate p19Arf
   from p16Ink4a (note section 9.3). This is the author's call.
5. **Standing from 2026-09-10: the orthotopic conclusions of scripts 50 and 51 are void.**
   Read `docs/2026-09-09_orthotopic_identity_correction.md` before touching that dataset;
   3e and section 5 have the detail.

---

## 3. What happened this session (2026-09-22)

Three commits, all pushed: `dd58b8d`, `d3d7f3b` and `1a43f04`. There were also branch and tag
operations, and two read-only reports: the `analysis-exploratory` resolution (3.1) and the
FOXO3 timing report (3.3).

### 3.1 The move recorded, and every other branch retired

- **`dd58b8d`:** one line each in `CLAUDE.md`, this file and
  `docs/experimental_cohorts_branch_notes.md`: work continues on `main` from `b7253f7`, and
  `experimental-cohorts` is retired. The trees were checked and are identical.
- **`analysis-exploratory` held nothing new.**
  - Its two commits (`bedb15e`, `46bcba7`) are script-23-only copies of `07b6516` and
    `fe6d004`, which have been on `main` since 13 July.
  - Script 23 is byte-identical on both (blob `2994a74`), so a cherry-pick would have been
    empty.
  - The author accepted the correction. The tip was tagged `archive/analysis-exploratory` and
    the branch deleted.
- **What the resolution report established,** in case PARTs 2b or 5b are asked about again:
  - PART 5b's Foxo3 numbers are the same estimate as script 54's check 2: the same object,
    contrasts and gene ID. 35 of 35 roster genes map identically by script 23's route and by
    the reconciler.
  - PART 2b's gene table is duplicated, for the shared genes, by script 42's
    `exclusions$p53_axis`. Its extra genes (Perp, Aen, Ei24) do not move.
  - Both saved tables reproduce from `interaction_results.rds` with zero deviation.
  - Nothing downstream reads `$h1$arf_p53` or `$puma_regulators`.
- **Then the author deleted the rest:** `experimental-cohorts`, `new-analysis`,
  `paper-figures`, `paper-final` and `BlockA-revision-step-by-step`, from the remote.
  - Every tip was checked first, and all are ancestors of `main`.
  - The local `experimental-cohorts` was deleted with `-d`.

### 3.2 Foxo3's target programme does not follow its transcript -- `d3d7f3b`

**The finding.** The record is the verification note, section 9.1.
- **The source is script 47 PART I** (`biogenesis_axis_developmental.rds$foxo_lanes`,
  `TFT_FOXO3_CHUNG`, 38 genes):
  - NES -1.76 at 6W and -1.80 at 12W, adjusted p 0.003 for each within the TF category;
  - interaction NES +1.09, p 0.30;
  - FOXO1 null as the control.
- **A second method agrees:** script 23 PART 5b's programme score gives genotype -0.351
  (p 0.006) and interaction +0.072 (p 0.76).
- **PART 5b's per-animal correlations are not cited.**
- **Against the transcript:** Foxo3 is +0.473 in wild type and -0.013 in Myc+, interaction
  -0.486. At 6W, Foxo3's message is higher under MYC while its output is lower.

**What changed.**
- **Fig. 2H** (`fig2_departure_from_dose.R`):
  - a new bound reads the programme's NES from script 47's object, with asserts;
  - the prior now reads "a p53-independent activator of PUMA", without "canonical".
  - It was rendered, and the drawn panel is unchanged.
  - The author asked for "lower under MYC at both ages without a genotype difference". Since
    "lower under MYC" is itself a genotype difference, it is written as "no difference between
    the ages is detected".
- **One manuscript sentence (the author's),** note 9.2: "Foxo3, which switches PUMA on without
  p53" becomes "Foxo3, a p53-independent activator of PUMA". The clause asserted activity in
  this tissue, and the programme shows output going the other way. The manuscript is outside
  the repo, so the change is recorded, not applied.
- **Two ARF caveats,** note 9.3 and CLAUDE.md:
  - Cdkn2a's baseMean is 5.6 counts;
  - a gene-level count cannot distinguish p19Arf from p16Ink4a.
- **Do not cite script 23 PART 2b's module test** (note 9.4, CLAUDE.md). Its target set
  contains Bbc3 and Bax, and it tests against mu = 0.
- **CLAUDE.md's stale claim is corrected.** It said script 23's object had 9 convergence rows
  and that 2b/5b had never run. The object has 11 rows and both outputs.
- **PANELS.md:** the Fig. 2H entry and its slot row are updated. The second 2H row followed in
  `1a43f04`.

### 3.3 The read-only report: does the FOXO3 chain have a timing problem?

**The author's framing.** Bbc3's fall is entirely an interaction: flat in wild type, -0.541
under MYC, appearing between 6 and 12 weeks. A programme depressed under MYC at both ages
cannot explain a change within the window. So either the link is transcript-level, which is
weak, or the chain has a gap. Three questions were asked. Nothing was built or saved.

**1. The programme's interaction is not detected, and the interval is wide.** On 5b's score
(38 genes, per-animal z units, 20 df):

| | estimate | 95% CI | p |
|---|---|---|---|
| MYC effect, 6W | -0.387 | -0.74 to -0.04 | 0.032 |
| MYC effect, 12W | -0.315 | -0.66 to +0.04 | 0.075 |
| interaction | **+0.072** (SE 0.237) | **-0.42 to +0.57** | 0.76 |

- **The interval spans a lot.** Its lower end is a further drop about the size of the whole 6W
  effect. Its upper end is the suppression erased and reversed. This is "not detected", not
  "absent".
- **The point estimate has the wrong sign for the chain:** slightly less suppression at 12W.
- **On a log2 scale** (the same genes, centred per gene, not scaled):
  - the programme's interaction is +0.079 (-0.15 to +0.31);
  - Bbc3 in the same model is -0.493 (-0.83 to -0.15).
  - So the interval excludes a fall the size of Bbc3's, but not a smaller one. A 38-gene
    average can also hide a subset that fell.
- **The 12W MYC effect alone misses 0.05 on this score.** "Lower at both ages" rests on script
  47's NES and on the pooled genotype effect (p 0.006).
- **Script 47's NES has no interval.** fGSEA returns an NES and a p. An interval would need a
  bootstrap over animals, which is new analysis.

**2. The FOXO3 targets** (`interaction_results.rds`, raw MLE):

| gene | baseMean | in the Chung set | 6W MYC: LFC (SE), p | 12W MYC: LFC (SE), p | interaction: LFC (SE), p |
|---|---|---|---|---|---|
| **Pdk4** | 2085 | no | **+2.001 (0.465), 1.7e-5** | +0.232 (0.465), 0.62 | **-1.769 (0.657), 0.0071** |
| Sod2 | 1423 | yes | +0.303 (0.176), 0.086 | +0.230 (0.177), 0.19 | -0.073 (0.250), 0.77 |
| Cat | 897 | yes | +0.177 (0.138), 0.20 | +0.009 (0.139), 0.95 | -0.168 (0.196), 0.39 |
| Gadd45a | 547 | yes | +0.038 (0.181), 0.84 | -0.163 (0.183), 0.37 | -0.201 (0.258), 0.44 |
| Txnip | 3982 | yes | -0.201 (0.276), 0.47 | -0.183 (0.276), 0.51 | +0.018 (0.390), 0.96 |

- **None follows the programme with any support.**
  - Only Txnip has its shape, and it is nowhere near significant. It is also a known
    MYC-repressed gene.
  - Sod2 and Cat point the other way.
- **Pdk4 does not fall under MYC.** It rises about four-fold at 6W and is back near baseline
  by 12W. It carries an interaction, with the transcript's shape, not the programme's.
- **Three things weaken reading Pdk4 as FOXO3 output:**
  - **Most of its interaction is the global fade.** On Fig. 2H's scan it sits at z -1.44, the
    6.0th percentile, against about the 0.5th for Foxo3 and Bbc3. The 0.487 rate alone
    predicts a 12W effect of about +0.97.
  - **Pdk2, which is not a FOXO target, does the same:** interaction -0.85, p 0.008.
  - **Pdk4 is also an ERRa/PGC-1a and PPAR target,** and Esrra fades similarly (+0.68 at 6W,
    +0.37 at 12W).
- **It is not adipocyte contamination:** Spearman -0.11 against script 33's seven adipocyte
  markers over the 24 animals. Only adipocytes were checked.

**3. The pyruvate arm is off the diagonal, below it: MYC-specific by the declared magnitude.**
- **Where it is computed:** script 40, `background_vs_myc.rds$ruler`, row "Pyruvate
  metabolism".
  - It has 14 expressed genes: the PDH complex, Pdk1 to Pdk4, the phosphatases, and Mpc1/2.
  - Fig. 1F draws its 6W and 12W values on the priority ruler.
- **Priority (mitoPPS):**
  - The MYC effect is +0.328 at 6W (padj 0.042) and +0.006 at 12W. The "12W deprioritisation"
    is the loss of the 6W promotion, not a demotion below the compartment.
  - On the plane: wild type +0.166, Myc+ -0.156, interaction **-0.321** (95% CI -0.55 to
    -0.09, raw p 0.009, from a linear model on the per-animal scores).
  - That is beyond script 54's 0.20. It is rank 1 of 143 pathways, and does not survive
    correction over the 143.
  - Script 54's eight arms all sit within 0.048 of the diagonal.
- **Content:** the interaction is -0.365. The median pathway's is -0.183 (the global fade), so
  pyruvate ranks 7th.
- **Which way it points.** At the level of the arm it cuts against FOXO3, but only if the score
  is read as flux, and it is not flux: it sums the enzyme and its inhibitors.
  - **The loss is carried by the brakes.** The interactions run Pdk4 -1.77 and Pdk2 -0.85,
    then Pdha1 -0.54, Pdhb -0.48, Mpc2 -0.40 and Dlat -0.31.
  - **The fraction of the 6W effect kept at 12W** is 0.12 to 0.14 for the kinases, 0.33 to
    0.63 for the PDH subunits, and 0.78 for Mpc1.
  - Gene by gene, then, the brake retreats faster than the enzyme. That is the direction
    reduced FOXO3-type output predicts, but it is weak and not FOXO3-specific.

**What it means for the paragraph.**
- **The link is transcript-level.** Foxo3, Bbc3 and Pdk4 share "up at 6W, lost by 12W".
  Foxo3's and Bbc3's go beyond the global fade; Pdk4's mostly does not.
- **At the activity step,** the programme shows no detected change within the window, and its
  point estimate has the wrong sign. The interval leaves room for a moderate fall.
- **The chain has a gap at activity** that these data neither fill nor rule out.

**To regenerate** (read-only, about 5 s):

```r
suppressPackageStartupMessages(library(DESeq2))
source(here::here("functions", "reconcile_gene_symbols.R"))
h   <- function(...) here::here(...)
ir  <- readRDS(h("results", "interaction_results.rds"))
gs  <- readRDS(h("results", "gsva_scores.rds"))
gmt <- fgsea::gmtPathways(h("data", "genesets_from_library", "mammary_mito_myc_metab_v1_mouse.gmt"))
sm  <- as.data.frame(gs$sample_meta); sm <- sm[match(colnames(gs$expr_mat), sm$sample), ]
tp  <- factor(sm$timepoint, c("6W", "12W")); myc <- factor(sm$myc_status, c("neg", "pos"))
# 1. the programme score (script 23 PART 5b's construction): the interaction and its CI
g   <- intersect(gmt[["TFT_FOXO3_CHUNG"]], rownames(gs$expr_mat))
sc  <- colMeans(t(scale(t(gs$expr_mat[g, ]))))
fit <- lm(sc ~ tp * myc); print(round(cbind(est = coef(fit), confint(fit)), 3))
# 2. five FOXO3 targets: 6W MYC effect, 12W MYC effect, interaction -- LFC, SE, raw p
D <- lapply(ir[c("myc_6W_raw", "myc_12W_raw", "interaction_raw")], as.data.frame)
G <- c("Pdk4", "Sod2", "Cat", "Gadd45a", "Txnip")
e <- vapply(G, function(x) recon_to_ensembl(x, rownames(D[[1]])), character(1))
print(signif(do.call(cbind, lapply(D, function(d) d[e, c("log2FoldChange", "lfcSE", "pvalue")])), 3))
# 3. the pyruvate arm on both rulers (script 40), and its per-animal interaction
r <- as.data.frame(readRDS(h("results", "background_vs_myc.rds"))$ruler)
print(r[r$pathway == "Pyruvate metabolism", ], digits = 3)
ps <- as.data.frame(readRDS(h("results", "mitopps_scores.rds"))$mitopps_scores)
print(round(confint(lm(ps[["Pyruvate metabolism"]] ~ factor(ps$timepoint, c("6W", "12W")) *
                         factor(ps$myc_status, c("neg", "pos")))), 3))
```

### 3.4 Legends rebuilt, and the `combined_df_annotated_raw` audit -- `1a43f04`

- **`rebuild_panels.R`, a real run:** 28 of 28 panels ok in 13 s, and `legends.md` rewritten
  with 28 blocks, ASCII-clean.
  - Ten blocks changed since 16 August: 2E, 2G, 2G (alt), 2H, 2H+I (alt), 2I, S1E, S2C, S2D and
    S2D (alt). They match the verification pass and today's 2H bound.
  - The runner's glob is `^fig`, so the `plane_`, `biogax_`, `explainer_` and `fatpad_` panels
    were not re-rendered.
  - One ggplot2 warning ("Vectorized input to `element_text()`") came from a panel other than
    2H. It is not traced.
- **`grep combined_df_annotated_raw`:** 49 hits, and one was stale: PANELS.md's second Fig. 2H
  row, fixed in `1a43f04`.
  - **Fig. 2H+I (alt) and S2C** still read the file, as a symbol-to-Ensembl dictionary only.
    Every value they draw comes from `interaction_results.rds`.
  - **By design:** script 16 (required by CLAUDE.md), script 39, `99_session_bootstrap.R`, and
    four older exploratory figures outside `figures/panels/`.
  - **The rest** are records and history.

---

## 3a. 2026-09-21, evening -- Fig. 2G, the legend guard, and the seven-gene table

Two commits, both pushed: `fea4d60` and `98b61e0`. Then a read-only report, 3a.3, which the
author has not yet built on.

### 3a.1 Fig. 2G rebuilt on ruling 4 -- `fea4d60`, signed off

`figures/panels/fig2_priming_ratios.R` now reads `results/two_timeline_verification.rds` and
nothing else. It was built to the note's sections 1.2, 1.3, 1.5 and 3. It is 89 x 64 mm,
up from 58.

**What it draws.**
- **The line:** through the origin, fitted to the eight ratios other than PUMA:Bcl-xL. Slope
  0.429 (SE 0.137), residual SD 0.221 on 7 df. It is printed "slope 0.43" in Fig. 1G's form.
  No rate is imported, and no hard-coded 0.55 enters.
- **The band:** the line's 95% prediction interval (`$ratio_band`), drawn in full.
- **The points:** all nine ratios (ruling 3). Shape is the denominator: circles Bcl-xL,
  triangles Mcl-1. PUMA:Bcl-xL is the one open point, keyed **"not used to fit the line"**.
  - The key first said "not in the fit". The author read that as a verdict ("does not fit"),
    against a point that lies inside the band.
  - The key now names the method, and the legend's first sentence says PUMA:Bcl-xL "is read
    against a line it did not help to draw".

**What the legend says.**
- **The result, without a verdict:**
  - PUMA:Bcl-xL is at +0.674 at 6W and -0.061 at 12W, against +0.289 predicted by the line.
  - Its residual is -0.350, or -1.58 residual SDs, the largest of the nine.
  - It lies inside the prediction interval, -0.276 to +0.855.
  - No "off the line" criterion was declared, so none is applied, and the panel does not
    support "significantly".
- **What else is reported, not drawn:** the seven-ratio sensitivity fit (slope 0.456,
  residual -0.368, still inside), check 3 (IQR 0.444 >= 0.30), and the nine six-week
  p-values.
- **Bax:Bcl-xL's retention** is given as a computed 0.5495, not compared with any rate.
- **The ratio's own interaction is kept apart:** -0.735, raw p 0.037, Benjamini-Hochberg 0.33
  across the nine. Its null is no attenuation, and eight of the nine ratios attenuate, so it
  is not this panel's comparison.
- **The members are not on this panel** (added at sign-off). The manuscript's clause about
  the pro- and anti-apoptotic members cites the four-gene plane. The legend says so in its
  first sentence and in a bound, and it cites the plane by filename (`plane_four_genes.R`),
  which is asserted. Renaming that file when it gets a slot will stop 2G, on purpose.
- **"Priming" is cleared** from the block (N3), and an assert enforces it.

**Asserted, so a re-run of script 54 that changes one stops the panel:**
- the drawn line is the refit of `$ratios` (to 1e-12);
- the band is that line's closed-form prediction interval (to 1e-10);
- PUMA:Bcl-xL is rank 1 of 9, below the line, and inside the interval, on both fits;
- check 3's verdict;
- retention = d12/d6 for all nine.

**Gone from the 2026-08-05 panel, and the second is reversible on the author's word:**
- the dashed residual drops: Bak1:Bcl-xL sits 0.014 from PUMA:Bcl-xL on x, so a drop would
  pass about 0.5 mm from Bak1's point;
- the filled/open "significantly induced at 6W" encoding: it scoped a retention, which
  ruling 4 retires. The p-values are in the legend instead.

**Also this commit:**
- the note's lower interval edge is corrected from -0.277 to **-0.276** (-0.27645, which had
  been rounded twice);
- script 42's matched-pair null (p 0.082 in the old legend) is not carried over, because it
  is not in script 54's object.

### 3a.2 `panel_legend()` guards every panel -- `98b61e0`

**The guard.**
- **The failure.** An item built from a zero-length value evaluates to `character(0)`, and
  `c()` drops it before `panel_legend()` sees the vector. This is how the PUMA:Bcl-xL item
  went missing from 2G's first draft, with no error.
- **Why it reads the call.** The value alone cannot show that an item is missing.
  `panel_legend()` therefore evaluates each argument of `detail`, `bounds` and `source`
  separately, and errors on NULL, zero-length, NA or empty text, naming the item. A NULL
  passes only from an explicit `if (...)`.
- **Tests.** Twelve deliberate cases, including a reproduction of the 2G bug. All 39 panels
  were then sourced under the guard.
- **Not caught:** a zero-length value inside one item, as in `paste0("x ", sprintf(...))`.
- 2G's local count check is removed, because the guard covers it.

**It found another, in Fig. 2G (alt), and it was the panel's claim.**
- `fig2_death_two_timelines_alt.R`'s item "THE CLAIM, AND THE SCRIPT ASSERTS IT" named Bik
  among the other BH3-only sensors. Script 42's roster has never carried Bik.
- So `say("Bik")` returned `character(0)`, and so did the whole item. It has been missing
  from every rendered legend since the panel was built (`37c3e61`, 2026-08-09). The
  16 August `legends.md` has the block without it.
- **The fix:** `OTHER_SENSORS` (Bcl2l11, Bid, Pmaip1) is declared and asserted equal to the
  roster's BH3-only sensors other than Bbc3 and Bmf.
- **The restored item is unreviewed text, and the author's call.** It says Bbc3 is -0.479
  under Myc (padj 0.016), "the ONLY BH3-only sensor that is flat in development and
  significantly down under Myc". That is a claim on the temporal arm, which the panel's own
  first bound says is described, not claimed. It is restored as written, minus Bik.

### 3a.3 The seven-gene table: Bbc3 against six others on the 6W MYC effect and the interaction

**The author's question:** from `interaction_results.rds` and script 42's saved object, report
the six-week MYC effect and the interaction (log2FC, SE and raw p for each) for Bbc3,
Bcl2l11, Bax, Bcl2l1, Bak1, Bmf and Htra2. The purpose is to see whether "induced at 6W and
lost by 12W" separates Bbc3 from the others on those two axes, before deciding whether it is
a panel. A read-only report; nothing was built or saved to the repo.

**Sources.**
- `results/interaction_results.rds` (script 03): the contrasts `myc_6W_raw`, `interaction_raw`
  and `myc_12W_raw`. Raw (unshrunken) MLE log2FC, Wald SE and raw Wald p; n = 6 per group.
- `results/priming_arm_teb.rds` (script 42, re-sourced 16:29): `$machinery`, for the arm
  labels and as a second symbol route.
- Symbols are mapped to Ensembl IDs with `recon_to_ensembl()`
  (`functions/reconcile_gene_symbols.R`), one ID per gene.

**Controls.**
- Script 42's saved `lfc_myc_6W`, `lfc_myc_12W` and `lfc_interaction` equal
  `interaction_results.rds` at those IDs exactly: the maximum deviation is 0.
- Script 42's `baseMean` is `round(baseMean)` (`42:160`) and matches after rounding.
- The identity 12W = 6W + interaction holds to 1.1e-16.

| gene | arm (script 42) | baseMean | 6W MYC effect: LFC (SE), raw p | interaction: LFC (SE), raw p | 12W MYC effect (= 6W + interaction): LFC, raw p |
|---|---|---|---|---|---|
| **Bbc3** | BH3-only trigger (PUMA) | 152 | +0.258 (0.139), 0.064 | **-0.541 (0.204), 0.0081** | -0.283, 0.059 |
| Bcl2l11 | BH3-only trigger (BIM) | 2305 | -0.263 (0.154), 0.088 | -0.038 (0.218), 0.86 | -0.301, 0.051 |
| Bax | effector (pro) | 481 | +0.484 (0.103), 2.8e-6 | -0.231 (0.150), 0.12 | +0.253, 0.019 |
| Bcl2l1 | brake (Bcl-xL) | 870 | -0.418 (0.150), 0.0054 | +0.179 (0.213), 0.40 | -0.239, 0.12 |
| Bak1 | effector (pro) | 472 | +0.258 (0.200), 0.20 | -0.302 (0.284), 0.29 | -0.044, 0.83 |
| Bmf | BH3-only trigger | 303 | -0.243 (0.215), 0.26 | -0.188 (0.304), 0.54 | -0.431, 0.045 |
| Htra2 | execution (IMS protease) | 367 | +1.489 (0.226), 4.0e-11 | -0.546 (0.322), 0.089 | +0.943, 3.9e-5 |

The 12W column is derived; its p is the 12W genotype contrast's own.

**Derived values.**

| gene | Bbc3 | Bcl2l11 | Bax | Bcl2l1 | Bak1 | Bmf | Htra2 |
|---|---|---|---|---|---|---|---|
| interaction z | -2.65 | -0.17 | -1.54 | +0.84 | -1.06 | -0.62 | -1.70 |
| 12W/6W | -1.10 | +1.14 | +0.52 | +0.57 | -0.17 | +1.77 | +0.63 |

**The reading: "induced at 6W and lost by 12W" does not separate Bbc3 at the gene level.**
- **Bbc3 is not significantly induced at 6W:** +0.258, p 0.064.
- **Bak1 has the identical 6W effect (+0.258)** and has also fallen to about zero by 12W
  (-0.044). On this pattern Bak1 looks like a noisier Bbc3.
- **The two genes MYC clearly induces at 6W keep about half of it at 12W:** Bax (p 2.8e-6)
  keeps 0.52 and Htra2 (p 4e-11) keeps 0.63. They weaken; they are not lost.
- **Bcl2l1, Bcl2l11 and Bmf are lower in Myc+ at 6W,** so the pattern does not apply to them.

**What does separate Bbc3.**
- Its interaction relative to its SE: -0.541, z -2.65, raw p 0.0081. It is the only one of
  the seven below 0.05; the genome-wide IHW value is 0.84.
- Htra2's interaction is the same size (-0.546), but its SE is larger (z -1.70, p 0.089), and
  it comes from a 6W effect six times bigger.
- On the two axes, Bbc3's interaction is about twice its 6W effect. That takes it past zero
  at 12W (-0.283, p 0.059), which no other gene that starts positive does. But this reading
  rests on a 6W effect that is not distinguishable from zero (z 1.85). Script 42's correction
  #1, as script 44 restates it: you cannot lose an effect you never had.

**The "induced at 6W" half belongs to the PUMA:Bcl-xL ratio, and mostly to its denominator.**
- **At 6W**, script 42's `$priming` (per-animal OLS on log2(count + 1)) gives PUMA:Bcl-xL at
  +0.674, p 0.0089. On the DESeq2 members, Bbc3 minus Bcl2l1 is +0.258 - (-0.418) = +0.676:
  **38% Bbc3 and 62% Bcl-xL**, which is lower under MYC at 6W (p 0.0054).
- **The ratio's interaction** is -0.735 (p 0.037). The members give -0.541 - (+0.179) =
  -0.720: **75% Bbc3 and 25% Bcl-xL**. Bcl-xL's MYC repression weakens at about the rate Bax
  does: -0.418 to -0.239, retention 0.57.
- So 2G (alt)'s legend line, that the ratio's reversal is "entirely its numerator", is
  about three-quarters true on these numbers.

**Licence.**
- Bbc3 and Bcl2l11 are script 44's two `pre_specified_genes`, per Fig. 2H's legend. They
  split: Bcl2l11 has no interaction (-0.038, p 0.86).
- Bcl2l1 is the fixed denominator of the pre-specified pair.
- Bax, Bak1, Bmf and Htra2 are the comparison set chosen for this question. So "the only one
  below 0.05 of the seven" is a descriptive ranking, not a test.

**What any build has to respect.** These are constraints, not a recommendation; the author's
suggestions come next.
- A gene-level panel must not mark Bbc3 as induced at 6W (p 0.064).
- "Lost by 12W" applies to Bak1 as well.
- The Bbc3 interaction is nominal (raw p 0.0081, IHW 0.84), and its licence is
  pre-specification.
- Any class such as "induced" or "lost" defined now would be defined after the numbers were
  seen. It must be labelled descriptive, as the note's rules require.
- The ratio's 6W induction is mostly the denominator.
- If both readings are drawn, the ratio (2G) and the members (the four-gene plane) must not
  contradict each other.

**To regenerate** (read-only, about 5 s):

```r
suppressPackageStartupMessages(library(DESeq2))
source(here::here("functions", "reconcile_gene_symbols.R"))
ir <- readRDS(here::here("results", "interaction_results.rds"))
pa <- readRDS(here::here("results", "priming_arm_teb.rds"))
G  <- c("Bbc3", "Bcl2l11", "Bax", "Bcl2l1", "Bak1", "Bmf", "Htra2")
D  <- lapply(ir[c("myc_6W_raw", "myc_12W_raw", "interaction_raw")], as.data.frame)
U  <- rownames(D$interaction_raw)
ens <- vapply(G, function(g) recon_to_ensembl(g, U), character(1))
stopifnot(!anyNA(ens))
m  <- as.data.frame(pa$machinery); m <- m[match(G, m$gene), ]
v  <- function(k, col) D[[k]][ens, col]
stopifnot(max(abs(m$lfc_myc_6W      - v("myc_6W_raw",      "log2FoldChange"))) < 1e-9,
          max(abs(m$lfc_interaction - v("interaction_raw", "log2FoldChange"))) < 1e-9)
data.frame(gene = G, arm = m$arm, baseMean = round(v("myc_6W_raw", "baseMean")),
           lfc6 = v("myc_6W_raw", "log2FoldChange"), se6 = v("myc_6W_raw", "lfcSE"),
           p6   = v("myc_6W_raw", "pvalue"),
           int  = v("interaction_raw", "log2FoldChange"), seI = v("interaction_raw", "lfcSE"),
           pI   = v("interaction_raw", "pvalue"),
           lfc12 = v("myc_12W_raw", "log2FoldChange"), p12 = v("myc_12W_raw", "pvalue"))
```

---

## 3b. Earlier the same day (2026-09-21) -- script 54 and three rounds of rulings

The record is `docs/2026-09-21_two_timeline_verification.md`: every outcome against rules
fixed before retrieval, the rulings of record in section 3, and the four rates in section 8.
Commits run from `932d592` to `b21c6b0`. In short:

**The three checks.**
- **Check 1 FAILED.** Bax's interaction is -0.231, beyond the declared 0.20 and with Bbc3's
  sign. "Only the PUMA to BCL-XL balance behaved differently" is **withdrawn**.
- **Check 2 passes on the amended rule.** Foxo3 "rose in the normal gland and did not rise
  under MYC": +0.473 in wild type, -0.013 in Myc+, interaction -0.486. It is exploratory:
  its position is shown and its p is labelled, and the manuscript quotes no p for it.
- **Check 3 passes.** The IQR is 0.444, so 2G stays a scatter, with all nine ratios.

**The coupling (Fig. 2I): an impasse at n = 24.**
- The unadjusted interaction clears its within-timepoint permutation null (p 0.026); the
  adjusted one does not (p 0.085).
- Withdrawn: "only where MYC was present", the robustness clause, and the redox sentence.
- In Figure 1 the coupling asserts nothing. Both p-values are printed on the panel, side by
  side.

**The arms.**
- "MYC neither causes nor prevents the respiratory withdrawal" is **withdrawn**. The
  replacement: most of the withdrawal occurs on both timelines (content ruler, 63% shared).
- All seven mitochondrial arms lie beyond their matched nulls on the content ruler. On
  mitoPPS each is within 0.048 of the diagonal.

**The figure layer.**
- The adjusted p on disk is IHW, not BH; five legends are relabelled.
- "Priming" is cleared from the legend blocks (N3).
- The third round: exploratory p-values (Foxo3, Bax) are **shown and labelled** in legends,
  never quoted in the text.

**Task D and the rate audit.**
- The Figure-1 paragraph's numbers are asserted in scripts 32, 43, 44 and 33, all re-sourced,
  and every assert passed.
- There are **four rates, and none substitutes for another:** 0.450 (all reported genes),
  0.487 (the genes MYC moves), 0.552 (the mitochondrial compartment) and 0.43 (the eight
  ratios).
- The hard-coded 0.55 is none of them. `figures/fig05_death_arm.R` and
  `figures/figure2_developmental_window.R` are marked SUPERSEDED.

**Recorded and not pursued** (the note, section 5): the unadjusted Myc+ coupling slope is
largely timepoint.

## 3c. 2026-09-10 -- scripts 52 and 53, the orthotopic escape series

This is the first handoff written after them. The record is
`docs/2026-09-10_orthotopic_escape_series.md`.
- **Script 52** re-scored all 49 vector and CRISPR samples in one run, as a dose-of-escape
  series.
- **It is a selection experiment:** every arm is a survivor population, so nothing there says
  PGC1a regulates anything.
- **Outcomes:**
  - R4, endogenous `Bcl2l1` null, **closes C1**.
  - R2 and R3: `Bbc3` falls in the construct-free arm too, and only 1 of 5 p53 targets moves
    with it. The escape-route explanation is not supported, and `Bbc3` tracks retained PGC1a
    dose.
  - R1 holds on the letter of its rule, not on its substance.
- **Script 53** draws the series as a dose figure: three panels off script 52's object, with
  x the RETAINED PGC1a dose, which is an outcome of selection and not an assigned dose.
- `docs/experimental_cohorts_branch_notes.md` was last updated with script 52's run
  (`faf1e8a`), so it records 52 but not 53. For the detail, read the note.

---

## 3d. History (2026-09-08 to 09-10) -- script 51, specificity and normalisation

> **PARTIALLY VOID as of 2026-09-09 (sample identity).** The arm this section calls `Pgc1a` is
> PGC1a + Bcl-xL, so C1 and everything resting on `Bcl2l1` in that arm fall. **C3, the
> composition gate with its positive control, R3 (demoted) and R5 stand.** Preserved as the
> record of what was concluded and why; see section 2 item 0.

**Not a fourth claim.** C1 holds, C2 holds as an observation, C3 is uninformative,
and none of those verdicts was reopened. Script 51 rebuilt the 21-sample cohort
from source and **proved it against 50 before computing anything new** --
median-CPM to 0.784%, 50's five C1 point estimates to 0.0238%, set sizes asserted
at 87 / 880 / 13. Full note:
`docs/2026-09-08_orthotopic_specificity_checks.md`. Figure:
`outputs/orthotopic/51_arm_geometry.pdf`.

**Gate 0, and the limitation that matters more than the length offset.** No
transcript-level quantification exists anywhere -- and the delivery folder is a
**results folder, not a pipeline tree**, so these files were never delivered
rather than merely lost. The consequence is **Bcl-xS**: gene-level `Bcl2l1` pools
two functionally opposite isoforms, C1 rests entirely on `Bcl2l1`, and nothing
here separates *"PGC1a raises the guardian"* from *"PGC1a shifts splicing toward
the short pro-apoptotic isoform"*. **A limitation on C1's INTERPRETATION, not its
measurement**, and the one a referee finds first.

**C1 got stronger and more specific.**

- `Mcl1` does not move: **-0.020 [-0.154, +0.127]** against `Bcl2l1` +1.539.
- `Bbc3` **falls**: -0.727 [-1.299, -0.376]. The gap widens on both limbs.
- Among the anti-apoptotic members **only `Bcl2l1` rises** -- `Mcl1` flat,
  `Bcl2l2` falls. **Not a general anti-apoptotic response.**
- **But `Bik` rises** +0.765, so the gap sentence is scoped to `Bbc3` and not to
  the sensitiser arm. That is in the note rather than omitted.
- **PGC1a is a broad biogenesis intervention with a 2:1 OXPHOS tilt** (numerator
  +1.925, denominator +1.011), not an OXPHOS-selective one. Describe it that way.

**C2 is normalisation-dependent, and the mechanism is measured.** The
80th-percentile `Myc` contrast against EV goes **0.0029 (CPM) to 0.0385
(median-of-ratios)**. The Pgc1a arm carries **14.7% of its library in MitoCarta
genes against 8.0% in EV** -- PGC1a nearly doubles the mitochondrial share, which
deflates every other gene's CPM in that arm, `Myc` included. **C2 stands as an
observation; quote 0.0385, or quote both.** The 09-08 note's 0.0029 must not
travel alone.

**C3 now has two independent reasons for being uninformative, and neither is
power.**

- **The arms cancel.** EV and BclxL sit at the *same* `ox_rel` (-0.292 against
  -0.291) with `Bcl2l1` 2.1 log2 apart -- a vertical displacement at constant
  respiration, which is the construct. Between-arm +0.500 against within-arm
  -0.381 cancel to the pooled -0.065.
- **Dropping the arm does not fix it.** C1's spec excluded BclxL from pooling and
  C3's spec pooled it, so 7 of 21 `Bcl2l1` values are set by a construct. On
  EV + Pgc1a alone the rho is **+0.657 excluding zero** -- but with two arms and a
  designed `ox_rel` gap that correlation **IS the group contrast**, circular and
  opposite in sign to the within-arm reading. Reported as a demonstration, not a
  result.

**R5 is uninformative and says so.** `Ppargc1a` spans only **1.51-fold** within
the Pgc1a arm, so the flatness of `Bcl2l1` cannot discriminate a threshold from a
dose-response. The rule was written in advance to say that rather than read the
flatness.

**One correction to the brief, made rather than applied silently.** It asked for
50's PART F C3 fit "on median-of-ratios rather than CPM"; 50 **already ran on
median-of-ratios**. Both are fitted; the endpoint carrying the sign flip is
normalisation-stable (-1.442 against -1.431).

**The 09-08 vector-series note gained a cross-reference banner with no number
changed**, and the branch notes gained two corollaries -- see section 5.

---

## 3e. History -- the orthotopic vector series (script 50)

`scripts/50_orthotopic_vector_series_scoring.R`, author-run 2026-09-08, reproducing the dry
run exactly; the note is `docs/2026-09-08_orthotopic_vector_series.md` and all 33 quoted
numbers reconcile against `results/orthotopic_vector_series.rds`.

**Gate 1, answered before anything was written.** **No transcript-level quantification exists
anywhere on disk** -- no `quant.sf`, no salmon dirs, no `tx2gene` -- so `tximport` is
unavailable and the transcript-length offset is lost; that qualifies every number. **MYAZ is
identifiable** (the parental implant line) but **`FaMY1` and `BoMY` are not determinable from
files on disk**, and three file-level facts put all three in a separate series: a different
animal-ID namespace (`BFGE`/bare against `BRNS###-#x`), n = 6 against 7, and the originating
analysis pairing them only with each other -- **there is no MYAZ-vs-EV contrast on disk**.
**Scoring set fixed at the 21 vector-series samples** before any number was computed. p21
series out of scope.

**C1 HOLDS -- the one clean causal result.** `Bcl2l1` +121.8 CPM [98.1, 151.2] in Pgc1a
against EV, one-sided **p = 0.00029**, which is 1/3432, the FLOOR for 7 v 7: the arms separate
perfectly. `Ppargc1a` 0.04 -> 302 CPM and both OXPHOS measures move with it, so the
manipulation worked where it claims to. **Survives composition adjustment** at +120.0
(p 0.0031) **with the positive control at 2.9e-9**. PGC1a alone raises the guardian in vivo
with no Bcl-xL construct present -- the human `BCL2L1` coefficient gets its causal reading.

**C2 HOLDS as an observation, and both statistics were needed.** The **max** contrast against
EV is the **weakest** cell (p 0.096) because EV carries one animal at 3170 CPM; the **80th
percentile**, which that animal cannot move, is the strongest (**p 0.0029**). So the ceiling
is not an artefact of one PGC1a tumour. Zero PGC1a tumours above 2000 `Myc` CPM against 3 EV
and 4 BclxL; `Bcl2l1` converges 1.25x against 1.58x and 1.83x. Exact null, all 3432 splits.
**Written as an observation, not a test.**

**C3 IS UNINFORMATIVE, and that is the item needing a decision.** `guardian_ratio ~ ox_rel`
= **+0.064 [-0.347, +0.493]**; every interval spans zero and the **sign flips** across the
three readings (pooled +, group-adjusted +, composition-adjusted -0.144). Reference
directions, never pooled: gland **+0.351**, human tumours **-0.31 to -0.46** -- the estimate
is compatible with both and with zero. **The reason is design, not power:** the vector series
was built to MOVE `ox_rel` between arms (+0.915), which makes it excellent for a group
contrast and poor for a within-state correlation.

**The gate, read more carefully than its own threshold.** `adipose` and `endothelial` tripped
|rho| >= 0.5 against `ox_rel` -- but per arm the loadings disagree (`adipose` +0.04 / -0.39 /
+0.71) and over all 67 `adipose` falls to **+0.09**. The pooled loading is substantially the
construct effect, so PART F is partly adjusting for the intervention itself. That is why the
positive control mattered, and it is now **Corollary 4** in the branch notes.

**Not licensed by this dataset:** the MYC x OXPHOS interaction (no MYC-low arm; it lives in
the iMMEC rtTA-MYC +/-dox x +/-PGC1a design), the abstract sentence above, anything about
apoptotic competence as a cellular state (N3 held throughout -- the word is never applied to
a transcript), and any pooling with the timeline or the human cohorts.

---

## 3f. History -- the fat-pad timeline, closed (script 49)

**New branch `experimental-cohorts`**, off `paper-final`, for the two datasets outside the
6W/12W purified-MEC pipeline: the orthotopic Bcl-xL / PGC1a series (not started) and the
fat-pad progression timeline (now closed). It **merges back into `paper-final`** -- numbered
scripts continue that sequence, data and READMEs follow the existing conventions, scratch
lives in a gitignored `sandbox/` and never merges.

**The fat-pad timeline is CLOSED. It cannot carry a respiratory claim, and the reason is
structural rather than a matter of power.** Whole mammary fat pad, adipose-majority
throughout (`Adipoq` 59,342 -> 26,415 on group medians), 30 samples, no 12W negative control.
Three tests, all pre-registered, all failed in an informative way:

1. **Pre-check** (scratch, `docs/2026-09-07_fatpad_timeline_oxphos_precheck.md`): can `ox_rel`
   support a within-tumour correlation? **FAIL** -- +0.318 to +0.561 against the adipose
   markers at n = 15. And the compartment-relative share made it **worse** than the absolute
   level, because adipose loads more on the respiratory arm than on the rest of MitoCarta.
2. **Script 49, the ordered trend** on the 12W-to-tumour limb (n = 20), the estimand the
   pre-check never asked: **T1 FAIL** (`ox_nuc_mtrib` tau -0.147, one-sided p 0.800) and
   **T3 VOID** -- the confound runs WITH the drift. The endpoint's two halves load on adipose
   with **opposite signs** (+0.14 nuclear OXPHOS, -0.24 mitoribosome), so their difference
   ADDS the tissue signal: +0.58 to +0.78 across four markers, **76% of its variance**
   adipose-explained against **7%** for `ox_rel`. The construction built to isolate the
   biology is the one the tissue gradient most contaminates.
3. **Script 49 PART G, the adjustment**: the pre-registered **negative control fails**.
   `Mki67` genuinely rises across the limb (tau +0.313, p 0.045) and goes to **-0.096
   (p 0.715)** after covarying out four adipose markers, while those markers explain only 15%
   of its variance. On this limb adipose depletion and tumour progression are the same
   variable, so removing the confound removes the biology. There is no corrected version.

**What the dataset still contributes**, neither needing a respiratory score: the **6W
genotype panel** (`PUMA:BCL-XL` +0.77 log2, Wilcoxon p 0.016) and **`Bcl2l1` flat across the
whole progression** while `Myc` rises ~11x, which forecloses "BCL-XL accumulates as tumours
grow". **The 6W panel now carries a caveat** -- `Adipoq` is matched between the 6W groups but
`Krt8`/`Epcam` are not; two of five `6WK_POS` samples have almost no epithelium (102-fold and
23-fold below their own group median) and carry the group's two highest ratios. Dropping them:
**+0.646, p 0.071** at n = 5 vs 3. Direction survives, significance does not. Do not cite the
p = 0.016 without it.

**Also flagged and NOT resolved:** `Myc` rises 7.5x (medians) from `6WK_POS` to `12WK_POS` in
fat pad, against the MEC dataset's stable-dose premise. Three readings stay live; the premise
is about purified MECs and must not be quietly generalised to tissue.

**The panel for collaborators:** `figures/panels/fatpad_confound.R` ->
`outputs/figures/panels/fatpad_confound.pdf` (183 x 66 mm). Three parts: the score tracks the
fat; its two halves load with opposite signs so the difference amplifies; adjustment destroys
a real trend. Discussion panel, `fatpad_` prefix, outside both runners -- **the 28-panel count
is unchanged** and `PANELS.md` is untouched.

**Two lessons written into the branch notes, and both are flagged for the human arm:**

- **A loading is a property of a SUBSET, not of a ruler.** Compute the admissibility argument
  on the exact subset the test runs on. `rho(ox_sub, Adipoq)` is +0.538 over all 30 samples
  and **+0.144** within the limb, and the mitoribosome's loading flips sign between them.
- **Any composition or purity adjustment must carry a POSITIVE CONTROL, read first.** Without
  `Mki67`, PART G would have read as a clean null.

A human cohort stratified by PAM50 / TP53 / purity is exactly where both bite. That session
does not read this branch -- **tell it**.

---

## 3g. History -- the gate model (script 48)

**The ask:** a model of MYC, OXPHOS subunits and PUMA/BIM that can be tested in human
tumours, verified in the mouse first.

**Finding: the model was already in the repo as a fitted regression, not as a model.**
`scripts/42_priming_arm_and_teb_substrate.R:487` (PART H, "the coincidence model") and
`scripts/43_substrate_specificity_and_tradeoff.R:255` (PART B) both fit
`Bbc3:Bcl2l1 ~ myc * oxphos + epi + imm`. What was missing was the form that survives
leaving the mouse: the mouse `myc` is a **binary genotype** and the axis is **mitoPPS**,
and a human cohort has neither.

### The model

`M` = MYC dose, `X` = respiratory state of the mitochondrial compartment, `T` = a
BH3-sensor : guardian log-ratio, all z-scored:

```
T  = b0 + bM*M + bX*X + bMX*(M x X) + covariates        (1)
M* = -bX / bMX                                          (2)
```

`bMX > 0` is the whole claim. `bX = 0` at `M = 0` collapses (1) to a one-parameter gate
-- the manuscript's title thesis written as an equation. `M*` is the MYC level at which
respiration stops being neutral and starts priming, comparable across species because
both inputs are z-scores. Two further layers (`P(BUFFER) ~ M x X`, then conditional
chemosensitivity) are **human-only by construction**: this dataset has no tumours.

### What the mouse verification says

1. **It is a pure gate, not a crossover.** With the genotype as `M`, dropping the OXPHOS
   main effect costs nothing (R2 0.807 both ways, anova p 0.973).
2. **The ruler question is settled.** `ox_rel` = mean z(87 nuclear OXPHOS subunits) -
   mean z(the other 999 MitoCarta genes) **beats mitoPPS and beats the absolute level**:
   +0.857 (p 0.0010, perm p 0.047) vs +0.779 (0.0052) vs **+0.595 (0.082, fails)**. It is
   also the only ruler separable from MYC (r 0.33 against 0.61). Adding proliferation
   makes the gate stronger; `redox_rel` (+0.09, p 0.78) stays the control that behaves;
   the mitoribosome rivals it marginally and **loses the head-to-head** (+1.086 p 0.040
   against -0.354 p 0.59).
3. **The guardian half carries it, not the sensor.** Fitted alone against
   expression-matched nulls: `Bcl2l1` **-0.670 (2.1st pct, p 0.0026)**, `Mcl1` **+0.806
   (98.9th, p 0.018)**, `Bbc3` +0.746 (94.8th, **p 0.059**). `Bbc3:Mcl1` -- same numerator,
   other denominator -- is **null (+0.183, p 0.63)**. `Mcl1:Bcl2l1` is the most robust
   endpoint anywhere in this corpus (+0.994, p 0.00014; **within-timepoint +0.351,
   p 0.00058**). BIM is dead centre (37.6th pct); NOXA is a clean negative (7.7th).
4. **The transfer is weak, and that is the headline.** Rescaled by each estimator's own
   genotype gap: `Myc` transcript recovers **74%** of the genotype effect, MSigDB MYC
   signature **15%**, DoRothEA **16%**.

### Reconciled and settled

Script 48 was author-run twice and all 47 numbers in its doc reconcile against
`results/gate_model_verification.rds`. The one thing the author's run overturned:
**the collapse to a pure gate holds for the trigger and NOT for the guardian balance** --
`Mcl1:Bcl2l1` has `bX` -0.541 (p 0.0040) and crosses over at `M*` = +0.545. Section 5
carries the scoped correction.

### Three corrections to the human plan

`docs/2026-08-27_human_validation_plan.md` sections 7.1, 7.2 and 7.3 are superseded:

| plan section | change |
|---|---|
| 7.1 MYC estimator | **dose over activity** -- MYC mRNA / 8q24 GISTIC primary; the signature demoted to the concordance check the plan already requires |
| 7.2 mitochondrial axis | **`ox_rel` primary**, absolute level as a reported sensitivity, mitoPPS for shape only |
| 7.3 endpoints | **`MCL1 - BCL2L1` co-primary with `BBC3 - BCL2L1`**, plus mandatory decomposition of both halves whenever `PRIME` is positive |

Power follows: n ~ 170 on a dose estimator for a third of the mouse effect, ~ 1193 on a
signature. TCGA (~1000) and METABRIC (~1900) are adequate on the former only.

**Artifact of the same content, for circulation:**
https://claude.ai/code/artifact/c0bd6218-c635-4e6c-bc7b-0e06c9679cd7

---

## 4. What was already done -- reading list in dependency order

**Read first on this branch now (added 2026-09-21):**
- `docs/2026-09-21_two_timeline_verification.md`: script 54, the rulings of record (section
  3) and the four rates (section 8). Read it before quoting any Figure-1 number.
- `figures/panels/PANELS.md`, from "The Figure-1 verification pass" to "Still open": what
  changed, panel by panel, including this session's 2G rebuild and the `panel_legend()` guard.
- `docs/2026-09-10_orthotopic_escape_series.md`: scripts 52 and 53.

The list below is carried over unchanged. Its PANELS.md line numbers date from 2026-09-10 and
have moved since, so search for the headings instead.

1. `docs/2026-08-02_results_narrative_final_order_D.md` -- the settled narrative (Order D):
   section 1 oncogene, section 2 tissue. The cut is measured against this.
2. `docs/2026-07-26_introduction_alignment_and_the_question.md` -- the trade-off frame,
   the question in its final form, the three solutions, the closing claim.
3. `docs/2026-07-25_death_narrative_dose_vs_competence.md` -- dose against competence;
   why MMTV-Myc and MYC-ER answer different questions.
4. `docs/2026-08-17_developmental_oxphos_decline_and_the_biogenesis_axis.md` (37 kB) --
   the PGC1a question, closed in the negative; what licenses the rescue as an
   *instrument*, not a mimic. Sections 6 and 7 are the paste-ready sentences and their risks.
5. `docs/2026-08-17_instrument_choice_for_the_reduced_figure_set.md` -- which ruler each
   panel should draw, and why Fig. 2I keeps mitoPPS.
6. `figures/panels/PANELS.md` -- the manifest of record. Slot map at line 1792; "Where the
   panels and the written sentences disagree" at 1564; the no-panel sentence list at 1947.
7. `docs/2026-09-02_myc_oxphos_priming_gate_model.md` -- the gate model and the three
   corrections it forces on the human plan.
8. **Read before touching the orthotopic or fat-pad data** (formerly `experimental-cohorts`):
   `docs/experimental_cohorts_branch_notes.md` (rules, **four** corollaries, dataset status),
   then `docs/2026-09-08_orthotopic_vector_series.md` (this session; sections 5 and 6 are the
   manuscript consequences), `docs/2026-09-07_orthotopic_analysis_plan.md` (its spec), and
   `docs/2026-09-07_fatpad_tumour_limb_oxphos_trend.md` (sections 3, 7a, 7b and the
   Discussion paragraph) and `docs/2026-09-07_fatpad_timeline_oxphos_precheck.md`.

---

## 5. Corrections that must not be re-introduced

**Added 2026-09-22: Foxo3, ARF and the pyruvate arm** (the note, section 9, and 3.3).
- **The Foxo3 transcript is not FOXO3 output.** The target programme is lower under MYC at
  both ages, with no interaction detected. Do not write that FOXO3 "switches PUMA on" in this
  tissue; the manuscript now says "a p53-independent activator of PUMA".
- **"No interaction" is not established for the programme.** The interaction is +0.072 with a
  95% CI of -0.42 to +0.57: not detected, not absent.
- **ARF cannot be read from Cdkn2a** (5.6 counts, and gene level cannot separate p19Arf from
  p16Ink4a).
- **Do not cite:** script 23 PART 2b's module test, or PART 5b's per-animal correlations.
- **The pyruvate arm's priority score is not pyruvate flux.** It sums PDH and its kinases, and
  its 12W loss is carried mostly by Pdk4 and Pdk2.
- **Pdk4's interaction (-1.77, p 0.007) is mostly the global fade** of a four-fold 6W effect:
  it sits at the 6.0th percentile of the departure scan. It is not evidence of FOXO3 output on
  its own.

**Added 2026-09-21: script 54, the rulings of record, and this session.**
- **Withdrawn, do not bring back:**
  - "Only the PUMA to BCL-XL balance behaved differently" (Check 1 failed on Bax).
  - "Respiratory priority tracked this ratio only where MYC was present".
  - The coupling's robustness clause, and the redox sentence.
  - "MYC neither causes nor prevents the respiratory withdrawal".
- **The coupling is an impasse at n = 24, not a result.** The unadjusted fit clears its
  permutation null and the adjusted one does not.
- **There are four rates, none substitutes for another,** and the hard-coded 0.55 is none of
  them. A ratio's retention compared with a gene- or pathway-level rate compares two
  estimators; ruling 4 retired that comparison.
- **PUMA:Bcl-xL does not deviate "significantly" from the eight-ratio line.** It is the
  largest departure of the nine, it lies inside the 95% prediction interval, and no
  off-the-line criterion was declared.
- **The adjusted p on disk is IHW, not BH.**
- **At gene level, Bbc3 is not significantly induced at 6W** (+0.258, p 0.064). "Induced at
  6W" belongs to the PUMA:Bcl-xL ratio, 62% of which is Bcl-xL falling under MYC (3a.3).
- **2G (alt)'s "entirely its numerator" is about three-quarters** on the DESeq2 members
  (3a.3).
- **A key or label must name the method, never read as a verdict.** "Not in the fit" was read
  as "does not fit".
- **An empty legend item vanishes silently.** `panel_legend()` now stops on one (3a.2). Do not
  build a panel's legend in a way that bypasses the guard, for example with a pre-built
  vector when the items could be written in `c()`.

**Added 2026-09-08 from script 51, and both are in
`docs/experimental_cohorts_branch_notes.md`:**

- **The composition gate's pooled |rho| threshold has a blind spot in BOTH
  directions.** Pooled-material with the per-group loadings disagreeing
  (`adipose` +0.621 pooled against +0.036 / -0.393 / +0.714 by arm), and the
  mirror -- `ox_lvl` against proliferation, **pooled -0.323 below threshold with
  the three arms at -0.643 / -0.679 / -0.964**. **Read the per-group breakdown
  before believing either verdict.**
- **The same lesson on the endpoint: a near-zero pooled correlation across
  designed arms can be two opposed components cancelling, not an absence.**
  Decompose into between-arm and within-arm before reading a pooled rho as a
  null -- and note that **dropping an arm does not fix it**, it makes the
  remaining correlation the group contrast.


Ten statements have been made in this project, tested, and found false. They keep
coming back.

1. **PRC does not clear its expression-matched null** (7.1st pct). `padj 1.9e-4` at
   baseMean 2102 is precision, not effect size.
2. **Never quote the within-age OXPHOS <-> proliferation r = 0.89 raw.** OXPHOS <->
   `CORE_MITO` is 0.98 in the same samples; 0.89 *is* the ceiling. Use gene-level
   percentiles against all expressed genes.
3. **"Myc blocks the FOXO3 rise" is false.** `Foxo3` rises in both genotypes, interaction
   padj = 1, and `Bbc3` is flat in the wild-type timeline. Neither genotype half clears
   0.05, so neither may be labelled on a panel -- only their difference.
   **SUPERSEDED IN PART, 2026-09-21** (the note, section 1.4, Check 2).
   - "Rises in both genotypes" is wrong on the raw DESeq2 contrasts: wild type +0.473
     (p 0.00022), Myc+ **-0.013** (p 0.92), interaction -0.486 (raw p 0.0074, IHW 1.000).
   - The sentence of record is "rose in the normal gland and did not rise under MYC". Do not
     write "fell under MYC", and keep "blocks" out as a causal verb.
   - What still stands: the IHW value, `Bbc3` flat in wild type, and the rule that the
     genotype pair (6W +0.243, p 0.057; 12W -0.242, p 0.059) may not be labelled.
   - Foxo3 is exploratory: its p is shown in legends, labelled as exploratory, and quoted
     nowhere in the text.
4. **Screen risers as well as fallers.** `Esrrb` was missing from the first ledger. The
   within-animal coupling was the discriminating test both times.
5. **NEW (2026-09-02): "respiration is protective without MYC" is false ON THE TRIGGER, and
   the qualifier is load-bearing.** For `Bbc3:Bcl2l1` the -4.80 wild-type slope on record has
   no covariates; adjusted for epithelial and immune content it is **-0.01 (p 0.95)** while
   Myc+ is **+0.87 (p 0.0019)** -- a **gate**, not a crossover, and the interaction gets
   stronger on adjustment. **But the guardian balance `Mcl1:Bcl2l1` behaves the other way:**
   its wild-type slope is genuinely negative and survives adjustment (`bX` **-0.541,
   p 0.0040**; dropping `X` costs R2 0.828 -> 0.724), so it is a real crossover at
   `M*` = **+0.545**. Never state the collapse without naming the endpoint.
6. **NEW (2026-09-07): a loading is a property of a SUBSET, not of a ruler.** An
   admissibility or specificity argument must be computed on the exact subset the test runs
   on. In fat pad, `rho(ox_sub, Adipoq)` is +0.538 over all 30 samples and **+0.144** within
   the 20-sample limb, and the mitoribosome's loading flips sign between them. A
   pre-registered argument built on the full-cohort numbers did not survive.
7. **NEW (2026-09-07): a composition or purity adjustment must carry a positive control,
   read FIRST.** `Mki67` rises across the fat-pad limb (tau +0.313, p 0.045) and goes to
   -0.096 after adjusting on four adipose markers that explain only 15% of it. Without the
   control that would have read as a clean null.
8. **NEW (2026-09-07): matching one composition marker is not matching composition.** The 6W
   fat-pad panel was cleared on `Adipoq` being matched; `Krt8` was not, and two of five
   `6WK_POS` samples have essentially no epithelium.
9. **NEW (2026-09-08): an interval that spans zero is not a direction.** Script 50's first
   draft called `rho +0.064 [-0.347, +0.493]` "gland-like, the abstract sentence must be
   withdrawn" **on the sign alone**. Caught before the author's run; the reading rule is now
   three-way and fixed in code, with a sign-stability check beside it. Any verdict that turns
   a sign into a direction needs the interval and the stability of that sign across readings.
10. **NEW (2026-09-08): when the intervention moves the composition markers, a pooled loading
   can look material while the per-group loadings disagree.** `adipose` loads +0.62 on
   `ox_rel` across the orthotopic 21 but +0.04 / -0.39 / +0.71 per arm and **+0.09 over all
   67**. Check the per-group breakdown before believing a pooled confound.

---

## 6. The first task when the cut text arrives

**Re-map the slots off the cut text and report before changing anything:** which panels
lose their sentence, which sentences have no panel. Then retire / keep / rebuild per slot,
with sign-off between each.

Retiring a panel is three steps, all of them required:

1. `git mv figures/panels/<x>.R figures/panels/retired/`
2. delete the stranded PDFs in `outputs/`
3. mark the row retired in `paper/analysis_record.qmd`, or its completeness check fails

---

## 7. Open items

**New 2026-09-22:**

- **The FOXO3 paragraph: transcript-level link, or a stated gap at activity?** This is waiting
  on the author (3.3).
- **Fig. S2D (alt)'s legend lacks the p19Arf/p16Ink4a isoform caveat** (section 2, item 4). It
  was offered, not done.
- **One ggplot2 warning in the rebuild** ("Vectorized input to `element_text()`"), from a panel
  other than 2H. It is not traced (3.4).
- **Older documents still say "canonical p53-independent activator":** the 2026-07-26 and
  2026-08-17 notes, narrative v3, `figS2_puma_inducers.R`'s header comment, and
  `paper/analysis_record.qmd`. They are left as written; newer supersedes older (note 9.2).

**Carried from 2026-09-21, still open:**

- **Whether the seven-gene reading (3a.3) becomes a panel, and which reading it draws.** This
  is waiting on the author's suggestions.
- **Fig. 2G (alt)'s restored claim item: keep, reword, or cut.** This is the author's call
  (3a.2). It is now in `legends.md`.
- **Fig. 2G (alt)'s "entirely its numerator"** is about three-quarters on these numbers (3a.3).
  Rewording it is the author's call; it has not been changed.
- **2G against 2G (alt) for the slot.** This was already open (below), and both changed on
  2026-09-21.
- ~~**`legends.md` needs a real `rebuild_panels.R` run.**~~ **Done 2026-09-22** (3.4).
- **The note's section 7 still calls the 33/42 re-source "pending".** Both ran on 2026-09-21 at
  16:29 UTC, and that line of the note has not been edited.
- **The superseded 0.55 drawings.** `figures/fig05_death_arm.R` and
  `figures/figure2_developmental_window.R` are marked SUPERSEDED. Retiring them, or pointing
  them at `defs$global_rate_fitted`, is the author's call.
- **`paper/analysis_record.qmd` still describes the old 2G and 2I,** by instruction. It
  sources panels by slug, so a render shows the new panels and their legend blocks.
- **Recorded, not pursued** (the note, section 5): the unadjusted Myc+ coupling slope is
  largely timepoint.

The lists below are carried over from 2026-09-10 and **were not re-verified this session.**
Script 52 has since closed C1 (3c), so read the orthotopic items against
`docs/2026-09-10_orthotopic_escape_series.md` first.

**The manuscript** (these were the `paper-final` items; that branch is archived, and the work is on `main`):

- **The text cut.** Blocked only on the text arriving. Everything else is done.
- Per-slot choice between each original panel and its alternative: **2F, 2G, 2H+2I, S2D**.
- `MYC` against `Myc` on the 2F/2G quadrant notes -- unresolved.
- A fifth arrow in `figS1_design_contrasts.R:60-93` if the diagonal panel is adopted.
- Whether the gate model earns a panel. It has none, deliberately.

**The orthotopic and fat-pad datasets** (these were the `experimental-cohorts` items; that branch is retired, and the work is on `main`):

- **THE C3 / ABSTRACT DECISION IS OPEN AND IS YOURS.** The orthotopic series cannot test
  *"mouse tumours resemble human tumours rather than the gland they arose from"*. Another
  source, or the sentence comes out. Nothing else in the paper depends on the answer.
- **R3 IS BLOCKED ON ONE FACT AND IT IS YOURS** -- was the standing Bcl-xL/OXPHOS finding
  **transcriptional or functional**? See section 2 item 0. The asymmetry holds and the
  directed-edge argument is available on the transcriptome; whether it is a contradiction or
  an interesting cross-measurement discordance depends entirely on that answer. **Do not
  write the argument until it is settled.**
- **C2's number changed and the old one must not travel alone.** Quote **p = 0.0385**
  (median-of-ratios), or quote both. The 09-08 note's 0.0029 is CPM and the Pgc1a arm's
  library is 14.7% MitoCarta against EV's 8.0%.
- **Does the orthotopic result earn a panel?** C1 is the cleanest causal statement in the
  branch and currently has none -- and script 51 made it **more** quotable, not less: `Mcl1`
  flat, `Bbc3` falling, only `Bcl2l1` rising among the guardians. C2 is explicitly a
  figure-and-observation and would need one if it is used at all. Neither is built; the
  manuscript panel count is still 28.
- **`Bcl2l1` cannot be resolved into Bcl-xL and Bcl-xS** and no file on disk can fix it. If
  the isoform question matters to the paper, it needs new quantification from whoever ran the
  pipeline -- the delivery folder never contained transcript-level output.
- The orthotopic series is **done for the questions asked of it**. If the reversal question is
  still wanted, it needs a design where `ox_rel` varies **within** a common state -- not one
  built to move it between arms. Script 51 section 5 shows *why* from two directions.
- ~~The branch is **pushed** (`origin/experimental-cohorts`, created 2026-09-08) and still
  unmerged. Rebase on `paper-final` before the next substantial piece of work.~~ **Superseded
  2026-09-22:** everything is on `main`, and both branches are deleted (3.1).
- `sandbox/precheck_fatpad_oxphos_confound.R` is scratch on disk, gitignored, and is deleted
  before merge. Only its dated note survives.
- **CLAUDE.md's branch section is stale** -- it describes `new-analysis`/`analysis-exploratory`/
  `paper-figures` and predates `paper-final` and `experimental-cohorts`. Flagged three times
  now, deliberately not fixed.
  Its workflow, coding, gene-set, mitoPPS/GSVA and data-trap sections all still hold.
- The `Myc` 7.5x fat-pad discrepancy is open and needs per-epithelial-cell measurement to
  settle. Do not import either dataset's `Myc` scale into the other.

**Everywhere (updated 2026-09-22):** `main` is the only branch, pushed through `1a43f04`,
and only this handoff's commit is local.

---

## 8. The parallel human session

The human arm runs in a separate Claude Code session and has so far written **only two
docs** into this tree and **no scripts** (`scripts/` stops at 48; nothing new in
`results/`). That is why running the two sessions in parallel has worked.
*(As of 2026-09-10. Updated 2026-09-21: the human arm now lives in two separate repositories
beside this one: `../myc_human_exploratory` and `../myc_human_validation`, each with its own
GitHub remote. They are not worktrees of this repo, and nothing from them is in this tree. On
`experimental-cohorts`, `scripts/` runs to 54, all of it mouse work. `git worktree list` also
no longer shows the `../myc_mouse_main` worktree that CLAUDE.md describes. It was not
recreated after the move off Drive; use `git show main:<path>` instead.)*

Before that session writes scripts or `.rds` files, give it its own worktree:

```
git worktree add ../myc_mouse_human -b human-validation main
```

The pattern is already established by `../myc_mouse_main`. Two sessions sharing one
working tree collide in three places: the git index (`index.lock`, and one session's
commit sweeping up the other's half-finished edits), shared `results/` and `outputs/`
(a silent read of a half-written `.rds`), and same-file edits (last write wins, no merge).
Sequential is only necessary if the human arm needs to *change* mouse scripts or results.

**Tell that session two things.** First, three of its plan's measurement choices are decided
by mouse data (section 3a) and two of the three were wrong. Second, the two traps from
section 3: a loading belongs to the SUBSET it was measured on, and any purity or composition
adjustment needs a positive control read before the endpoint of interest. A cohort stratified
by PAM50 / TP53 / purity is exactly where both bite, and that session does not read
this repo.

Note also that `docs/2026-08-27_human_validation_plan.md` is untracked by rule, so it is
**not on the remote** -- if that session expects to find its own spec on `origin`, it is not
there.

---

## 9. Working rules, verbatim

- **Option A:** Claude Code writes `scripts/`, `figures/`, `paper/`; the author runs the
  numbered scripts in Positron. Running a panel to verify is fine and expected.
- **No explanatory prose on a panel.** `theme_panel()` blanks title/subtitle/caption;
  every script carries a `panel_legend()` block. Panels talk by NAMING DRAWN ELEMENTS.
- **Use the declared scheme;** anything new is declared in `_panel_common.R`, never inlined.
- **Verify with** `qlmanage -t -s 2000 -o <dir> <file>.pdf`. Never `sips`.
- **Run `./paper/clean_icon_files.sh repo` before the first git command** (Drive mount).
  *(2026-09-21: the repo moved off Drive, to `/Users/gs/code/myc_mouse`. The sweep matters
  again only if it goes back onto a Drive mount; check `find .git -name $'Icon\r'` if in
  doubt.)*
  **Never `git add -A`** -- the six `.txt` exports and the untracked docs stay untracked.
  **Commit per signed-off change; do not push without a word.**
- **Check every claim against the data before drawing it,** and say when a sentence in the
  cut does not survive -- add it to the `PANELS.md` list rather than softening it silently.

---

## 10. Follow-up threads, after the manuscript

- **What reshapes mitochondrial priorities in high-mitochondria human tumours** (MYC-high
  against ER-high): `docs/2026-09-22_mitochondrial_priorities_followup.md`. Untestable at
  n = 21, not refuted. A follow-up needs about 116 both-forks tumours.
- **Esr1**, from the FOXO3 dilution test (2026-09-22): the only gene in the FOXO3 target set
  with a raw interaction p below 0.05 (+0.813, raw p 0.036, IHW 1.00).
  - MYC represses it at 6W (-0.813) and not at 12W; most of that is the global fade (7.7th
    percentile of script 44's departure scan).
  - Flagged, not pursued: `docs/2026-09-21_two_timeline_verification.md` section 10.3.
  - In TCGA the ESR1-high and MYC-high high-mitochondria groups overlap (21 tumours), and ESR1
    is inversely associated with the MYC score within them (follow-up note, section 3.1).
