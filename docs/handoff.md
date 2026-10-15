---
date: 2026-10-15
tags: [project/myc_mouse, handoff, main, nat-metab, figure-layer, figure-1, extended-data-2, extended-data-3, script-55, script-56, prism-format, type-5-6pt, box-plot-axes]
status: live handoff -- overwrite in place at the end of each session
relates-to:
  - figures/natmetab/FIGURE_RULES.md                    (READ FIRST for figure work: every rule of the new panel layer, with its reason)
  - figures/natmetab/_style.R                           (the style; figures/natmetab/_style_selftest.R must pass after any edit)
  - Google Doc 11.8.26_MK_Myc_Paper, tabs "Nat metab version" and "fig plans" (Drive id 1M-kEC1d_VJmbTG-SyYYfykZ3F1_aEiDtu0xtwFaPtJo)
  - docs/2026-09-21_two_timeline_verification.md       (script 54, the rulings of record, the four rates; sections 9-10 are the 2026-09-22 Foxo3 record)
  - docs/2026-09-22_mitochondrial_priorities_followup.md (parked human follow-up; section 10 below)
  - figures/panels/PANELS.md                            (the OLD layer's manifest: the August long-form version, kept as the record)
  - scripts/54_two_timeline_verification.R              (writes results/two_timeline_verification.rds)
  - scripts/47_biogenesis_axis_and_the_developmental_oxphos_decline.R (PART I: the FOXO3 target programme)
  - docs/2026-09-10_orthotopic_escape_series.md         (scripts 52 and 53; Fig. 3k-m will redraw script 53's three panels)
  - docs/experimental_cohorts_branch_notes.md           (the retired branch's notes; its cohort rule still holds)
  - docs/2026-09-09_orthotopic_identity_correction.md   (voids the conclusions of scripts 50 and 51)
  - docs/2026-09-02_myc_oxphos_priming_gate_model.md    (the gate model)
  - docs/2026-08-27_human_validation_plan.md            (the parallel arm; untracked)
  - docs/2026-08-02_results_narrative_final_order_D.md  (the narrative the cut is measured against)
---

# Handoff -- restart from here

**2026-10-15, END OF DAY: FIG. 1 IS RESIZED AND SET IN 5/6-PT TYPE (2a.34). Everything is pushed,
this handoff included.**
- **Fig. 1b goes to Prism:** the author redraws it from the R panel's numbers;
  `fig1b_myc_transcript.R` stays as their source (2a.33, 2a.34).
- **Fig. 1d, 100 x 60:** 5/6 pt with a reduced *P* exponent. The y title reads "Fraction of
  nuclear-encoded transcriptome (%)", not "total", because the denominator leaves out the mtDNA
  transcripts. **Box-plot axes**, the author's new rule: from the tick below the lowest animal to
  the tick above the *P* label.
- **Fig. 1e, 78 x 55**, and **Fig. 1f, 100 x 56**, both 5/6 pt. 1f now sits alone in its row.
- **Fig. 1g 66 + 1h 52 + 1i 56, all 56 tall: one 180-mm row.** g's key reads "z-score", has no
  frame, and moved left.
- **"z-score", hyphenated, everywhere**, ED 3a and 3d included.

FIRST: the author re-renders ED 3a and ED 3d, whose PDFs still say "z score"; every Fig. 1 PDF is
current. Then ask what comes next (section 2, item 2): the ED figures are still 6/7 pt.

Read this file first. Sections 1, 2 and 2a are current as of **2026-10-15**. Section 2a covers
Extended Data Fig. 2 (2026-10-12, morning, 2a.1-2a.6), the first Fig. 1g-h proposals
(2026-10-12, afternoon, 2a.7-2a.13), the arguments and proposals of 2026-10-13 (2a.14-2a.22), and the
247-mm page with two rounds of proposals and scripts 55-56 (2026-10-14, 2a.23-2a.31), and the
evening's changes and the promotion (2a.32), Fig. 1b in the Prism format (2026-10-15, 2a.33), and
the afternoon's resizing of Fig. 1d-i in 5/6-pt type (2a.34). All are pushed. Section 2b is the
2026-10-10/11 session (Fig. 1d-f and ED Fig. 2, old letters), 2c the 2026-10-10 session (Fig. 1b
and ED Fig. 1) and 2d the 2026-10-09 session, which started the Nat Metab layer. Section 3 is the
2026-09-22 session and 3a the 2026-09-21 evening session; 3b to 3g and 4 to 9 are carried over.
Section 10 holds the threads parked for after the manuscript.

---

## 1. State

| | |
|---|---|
| checkout | `main`, tracking `origin/main`. It is the only branch, local and remote |
| HEAD | this handoff's commit (2026-10-15, evening), on top of `8bb8242` (the g-h-i row, 5/6 pt, z-score), `7f2d178` (1d's box-plot axes, 1e), `dfa8986` (1d in 5/6 pt) and `c740aae` (the morning's handoff) |
| pushed | **everything, this handoff included**, at the author's word (2026-10-15: "commit and push and handoff") |
| branches | unchanged since 2026-09-22: only `main`. Tags: five `archive/*`, plus `block-a-reviewed`, `block-b-full` and `orthotopic-pre-identity-correction` |
| scripts | **55** changed in the evening (the 99% line, `p_emp_resp`, `outside_99`) and re-run by the author (23:36). **56** run by the author (matches the dry run exactly). **The next number is 57** |
| results | `results/window_and_chain_statistics.rds` (script 56) and `results/departure_from_myc_fade.rds` (script 55, with the 99% columns), both current |
| new figure layer | `figures/natmetab/`: `_style.R` (on 2026-10-15 also `lw_pt()`, `natural_size()`, `label_mm()`, `box_axis()`, the `exponent` argument of `fmt_p()` and `brackets()`, and `txt`/`txt_title`/`line`/`tips` arguments on `theme_nm()` and `brackets()`), `_style_selftest.R` (33 checks, all pass), `FIGURE_RULES.md` (the built table and the page are current), **twenty panels of record**: Fig. 1b, d-i; ED 1a-c; ED 2a-f; ED 3a-d. 20 of 20 pass the dry run. `page_mock_247.R` draws the three pages (Fig. 1's bench row, with 1b, is a placeholder there). `retired/` holds the three replaced panels; `proposals/` holds only superseded ones |
| old figure layer | `figures/panels/` is **untouched and kept as the record**; not used for the Letter |
| outputs | `outputs/` is gitignored. The four stale PDFs (old 1f, narrow 2c, old 2f, old 2g) were moved to `outputs/natmetab/_retired/`. **Every Fig. 1 PDF is current** (the author rendered 1b at 10:09 and 1d-i between 15:42 and 16:27 on 2026-10-15). **ED 3a and ED 3d are stale**: they still say "z score" (2026-10-14). The other ED PDFs are current |
| text | the Doc (tab "Nat metab version") was read on 2026-10-12 (copy in that session's scratchpad); **nothing was written to it**. Its Fig. 1g-i paragraph needs the fixes in 2a.25, and its ED citations renumber (2a.31) |
| repo | `/Users/gs/code/myc_mouse`, off Drive |
| human arm | unchanged, read-only. Fig. 4g-j's data live in `../myc_human_exploratory` |

**Untracked, and staying that way:**
- the six `.txt` session exports at the repo root;
- `docs/library_reference/gray_chea_mito_tf_shortlist.csv`;
- `docs/library_reference/2026-08-22_consensus_myc_double_hit_thread.md`;
- `docs/2026-08-27_human_validation_plan.md`.

---

## 2. Do this first

1. **Re-render ED 3a and ED 3d** (the author, Option A), for their "z-score" labels:
   `source(here::here("figures", "natmetab", "edfig3a_per_animal.R"))` and
   `edfig3d_coupling_both_fits.R` likewise. Every other panel is current. Sourcing a panel's script
   overwrites its PDF in `outputs/natmetab/<Fig>/`; Illustrator updates a linked placement and
   needs an embedded one placed again.
2. **Ask what comes next.** Fig. 1 is done in R (2a.34): 1b is the author's to redraw in Prism, and
   the bench row a-c is theirs. Open:
   - **Fig. 1f sits alone in its row** (100 of 180 mm): does it widen, or share the row?
   - **The ED figures are still 6/7 pt**: do they take 5/6 too, and the box-plot axis rule (ED 3a)?
3. **Legends** (the scripts print every number): 2a.30's list, plus 2a.32's and 2a.34's additions.
   - **Fig. 1b's legend carries the WT window**, *P* = 0.48, which is not drawn; the script prints
     all four *P* values, whether 1b is drawn in R or Prism.
   - **Fig. 1d's legend names the denominator**: the 13 mtDNA-encoded transcripts are left out of
     it, which is why the title says "nuclear-encoded".
   - The box sentence, for R and Prism alike: "boxes show the median and interquartile range,
     whiskers the minimum and maximum, and dots individual animals".
4. **The text** (the author's Doc; Claude has written nothing to it): the fixes in 2a.25, and the
   renumbering in 2a.31. Every legend number is printed by its script (2a.30).
5. **After Fig. 1:** Fig. 3k-m (script 53's three panels, 2d.1; the void orthotopic conclusions of
   scripts 50 and 51 matter again, 3e). Fig. 4g-j is open (2d.1).
6. **How each panel is built** (two bullets added 2026-10-15):
   - one script per panel in `figures/natmetab/`, sourcing `_style.R` and saved with
     `save_panel()`, which refuses a panel that breaks a rule;
   - size it for its page: **183 x 247 mm for the submission** (ED 180 x 247), rows kept modular for
     the 170-mm page at acceptance; **5 pt where font size is the limit** (FIGURE_RULES.md);
   - preview the PDF with Quick Look (`qlmanage -t -s 2000`) and look at it: the checks cannot see
     clipped or colliding labels, and poppler's `pdftoppm` draws the italic *P* upright;
   - for a choice between designs, show rendered proposals and a page mock-up, not a description;
   - a number the text quotes is asserted in the script that computes it, and a panel computes no
     test: a statistic no object holds goes into a numbered script first (as 56 did today);
   - after any edit to `_style.R`, run `Rscript figures/natmetab/_style_selftest.R`, then every
     panel under `options(natmetab.dry_run = TRUE)` (20 of 20 pass). For a change meant to leave
     the other panels alone, render each one before and after and compare them (2026-10-15:
     byte-identical once the PDF dates are stripped);
   - a panel among Prism graphs takes the Prism format (FIGURE_RULES.md; Fig. 1b is the model);
   - **Fig. 1's type is 5/6 pt** (2026-10-15): 5 pt for tick and row labels, keys, *P* values and
     tile numbers; 6 pt for axis and strip titles, group names and a colour bar's title. An
     exponent or *R*2 superscript is reduced (`exponent = "reduced"`; it prints at 4 pt);
   - **a box plot's y axis** runs from the tick closest below the lowest whisker to the first tick
     above the highest mark, a *P* label included (`box_axis()`, on a graph of fixed height);
   - **"z-score"**, hyphenated with the soft hyphen, everywhere;
   - **a size the author reads in Illustrator:** check the PDF's own size first (its MediaBox).
     Twice on 2026-10-15 the reading came from a different file (67.5 for a 66-mm page, 104 for
     100). Base `pdf()` writes page sizes in whole points, rounded down, so a 56-mm page reads 55.74.
7. **At the end of the figure sessions:** turn `FIGURE_RULES.md` and `_style.R` into the GS-lab
   figure skill.
8. **Carried items, possibly moot now: ask before acting.**
   - the seven-gene table (3a.3);
   - Fig. 2G (alt)'s claim item (3a.2);
   - the isoform caveat in Fig. S2D (alt);
   - the void orthotopic conclusions of scripts 50 and 51 (3e), which matter again for Fig. 3k-m.
9. **Do not reopen without the author:** the FOXO3 activity step (closed, 3.5), or the parked human
   threads (section 10).

---

## 2a. What happened on 2026-10-12 to 2026-10-15

- **2026-10-12, morning, 2a.1-2a.6:** Extended Data Fig. 2. Two commits, both pushed: `a927ce8`, and
  that morning's handoff, `ab62730`.
- **2026-10-12, afternoon, 2a.7-2a.13:** the first Fig. 1g-h proposals. Pushed: `6cb04f1` and
  `9fd2b16`.
- **2026-10-13, 2a.14-2a.22:** the arguments behind Fig. 1g-h, and three more proposals. Pushed:
  `8bde9c4` and `3114ec2`.
- **2026-10-14, 2a.23-2a.31:** the 247-mm page, two rounds of proposals, scripts 55 and 56. Pushed:
  `4fb3aa9` and `ccad756`.
- **2026-10-14, evening, 2a.32:** the author's notes after the first handoff, the promotion, and the
  min-max whiskers. Pushed (2026-10-15): `125bb7f`, `969ee5c`, `813993d` and `e37da20`.
- **2026-10-15, 2a.33:** Fig. 1b in the Prism format, two rounds of the author's notes. Pushed:
  `7bc98c3` and that morning's handoff, `c740aae`.
- **2026-10-15, afternoon, 2a.34:** Fig. 1b to Prism; Fig. 1d-i resized and set in 5/6-pt type;
  box-plot axes; "z-score". Pushed: `dfa8986`, `7f2d178`, `8bb8242` and this handoff.

### 2a.1 Extended Data Fig. 2 renumbered (the author's order)

The text cites the independence panel first, so it moved to the front of the window panels.
Scripts renamed with `git mv`; cross-references, the rules file and the PDFs updated.

| old | new | script |
|---|---|---|
| 2g | **2d** | `edfig2d_window_independent.R` |
| 2d | **2e** | `edfig2e_myc_reallocation_ages.R` |
| 2e | **2f** | `edfig2f_window_matched_controls.R` |
| 2f | **2g** | `edfig2g_window_two_rulers.R` |

Citations in the approved paragraph (2b.2) change accordingly:
- "(Fig. 1e, ED Fig. 2d)" becomes **2e**;
- "(ED Fig. 2g)" becomes **2d**;
- "(Fig. 1f, ED Fig. 2e,f)" becomes **2f,g**;
- "(Fig. 1d,f; ED Fig. 2e)" becomes **2f**.

For citation order, the "independently of MYC" sentence (2d) must come before the
MYC-at-both-ages sentence (2e).

### 2a.2 ED Fig. 2a-c built (from the old S1A, S1B, and 1C + S1C)

| slot | script | PDF (mm) | content |
|---|---|---|---|
| ED 2a | `edfig2a_design.R` | 46 x 40 | the 2 x 2 design, n = 6 per group (asserted); arrows for the four contrasts and the interaction, in each contrast's own colour, so the panel is also the colour key; text black |
| ED 2b | `edfig2b_geneset_library.R` | 58 x 40 | the 986-set library by category, filled by mitochondrial class in slate; script 37's classification re-derived over the whole library and asserted identical on the shared sets |
| ED 2c | `edfig2c_pathway_axis.R` | 72 x 40 | the pathway-level PCA (equal units, no broken axis; ringed group means) beside the ranked absolute PC1 loadings with OXPHOS, E2F targets and glutathione named; script 37's PCA rebuilt and asserted (885 sets, PC1-3 %, OXPHOS gate, concordance >= 0.9) |

- **Ported, not sourced:** `pathway_axis()` and `mito_class3()` from the old
  `_panel_common.R` are copied into the scripts that need them; the new layer does not source the
  old one.
- **ED 2c names three sets, not four.** The human breast-cancer MYC bicluster (MB2, 0.98, as high
  as OXPHOS) did not fit at 30 mm; the script prints it for the legend.
- **ED 2c's loading curve has no class key:** it uses ED 2b's, so the two sit side by side.

### 2a.3 The ED 2a-c sentence (approved with the plan, 2026-10-12)

The Doc's version: "...we explored the timeline of the MMTV-MYC-induced transcriptome of MECs
purified from the tumours, which was dominated by a mitochondria-driven trend (ED Fig. 2a,b,c)".
Three problems:
1. **"Purified from the tumours":** the cells came from 6- and 12-week glands.
2. **"Timeline":** PC1 is a genotype axis, not an age axis (2a.4).
3. **"Mitochondria-driven":** a loading is co-variation, not cause; the axis is common-mode; the
   top is mitochondrial by construction.

The replacement:
> ...we profiled MECs purified from 6- and 12-week glands. Across 885 gene sets, a single axis
> carried 77% of the variance; it separated MYC+ from WT glands at both ages and was led by
> mitochondrial gene sets (ED Fig. 2a-c).

**Asserted in `edfig2c_pathway_axis.R`:** 885 sets, PC1 rounds to 77%, MYC+ above WT and Welch
P < 0.05 at each age, and no non-mitochondrial set in the top 50.

### 2a.4 Numbers for the ED 2a-c legends (printed by the scripts)

- **Design:** each age was extracted as one batch, so the 12W vs 6W changes coincide with batch;
  the genotype effects and the interaction are within batch.
- **Library:** 986 sets in 9 categories; set size median 46 genes (5-1253):
  - MitoCarta 133;
  - curated (TF target sets intersected with MitoCarta) 260;
  - not mitochondrial 593.
  
  885 are scored per sample by GSVA, and 817 are fGSEA-eligible.
- **PCA:** PC1 76.5%, PC2 8.0%, PC3 4.3%.
  - `PC1 ~ age + genotype`: genotype P = 0.00067, age P = 0.19; interaction P = 0.43.
  - Within age (Welch): 6W P = 0.022, 12W P = 0.012.
  - Group means: WT 6W -8.4, MYC+ 6W +15.2, WT 12W -11.2, MYC+ 12W +4.4.
- **The bounds:**
  - r(PC1, mean of all sets) = 0.981; 79% of loadings one-signed; median |loading| 0.88;
  - top 50: 47 curated (by construction), 3 MitoCarta;
  - median |loading| by class: MitoCarta 0.90, curated 0.98, not mitochondrial 0.61;
  - MB2 0.98 (rank 156); glutathione 0.67 (rank 585); OXPHOS 0.98 (rank 114); E2F 0.87 (rank 450).

### 2a.5 The Extended Data page is 180 mm wide

**Nature's research figure guide gives 180 x 170 mm for Extended Data.** The 183 mm used before
is the main-figure width. So ED 2d-g were resized: 2d 114 x 60, 2e 60 x 60, 2f 92 x 64 and
2g 84 x 64.

ED 2 now runs in three rows:

| row | panels | width | height |
|---|---|---|---|
| 1 | a + b + c | 176 mm | 40 mm |
| 2 | d + e | 174 mm | 60 mm |
| 3 | f + g | 176 mm | 64 mm |

That is 164 mm of height before the gaps between panels.

### 2a.6 Style changes, all in `FIGURE_RULES.md`

- **`save_panel()` refuses an `EDFig` panel wider than 180 mm** (`NM_ED_MAX_W`).
- **The tick-end check applies to drawn axes only,** so a schematic (ED 2a, canvas in
  millimetres) passes.
- **Mitochondrial classes are slate:** `mito_class_cols` / `mito_class_labels`. The old
  red/orange/blue is out because blue and orange mean genotype.
- **Self-test:** two new cases; 20 of 20 pass.

### 2a.7 Fig. 1g-h: the request, and the paragraph against the data

**The author's request.** g and h "are interrelated and the basis of the paragraph describing the
relationship in the transcriptome between pro-, antiapoptotic proteins and OXPHOS ... a fundamental
part of the paper, where the figure has to be very clear and meaning should be the strongest
possible."
- The closest earlier versions are Fig. 2G (alt) (`outputs/figures/panels/fig2_death_two_timelines_alt.pdf`)
  and Fig. 2H+I (alt). The latter is page 3 of `outputs/figures/panels/Figure2_alternatives.pdf`
  (16 August); its current render is `fig2_puma_chain_alt.pdf`.

**The Doc's paragraph, checked sentence by sentence** (tab "Nat metab version", read 2026-10-12).

| Doc sentence | what the data say | change |
|---|---|---|
| "at 6 weeks MYC shifted the balance toward pro-apoptotic factors" | All nine BCL-2 pro:anti ratios rise under MYC at 6W: median +0.49 log2, five at *P* < 0.01 | Holds |
| "the pro-apoptotic pathway was progressively suppressed" | There are only two ages. And not the whole pathway: *Bmf* rises about 2-fold in both genotypes; *Bcl2l11* and *Pmaip1* rise slightly | Name the genes; drop "progressively" |
| "both *Bax* and *Bbc3* declined compared to the normal gland" | *Bbc3*: −0.48 in MYC+ (IHW *P* 0.016), +0.06 in WT; interaction *P* 0.0081 (pre-specified). *Bax*: −0.38 in MYC+ but also −0.15 in WT; interaction *P* 0.12 (exploratory) | *Bbc3* holds; *Bax* "fell in both, more in MYC+" |
| "the steepest drop was in *Bbc3*" | True within the BCL-2 family (*Bax* −0.38, *Bid* −0.36). False across script 42's wider roster: *Htra2* −0.67, *Birc5* −0.93 | "Of the BCL-2 family" |
| "shifting distinctly the PUMA/BCL-XL ratio" | It goes from +0.67 to −0.06, losing its whole shift, while the other eight keep 0.43 of theirs. It is the largest departure of the nine, but inside the 95% band (2d.6, item 4) | Drop "distinctly" |
| *Foxo3* "most distinct ... among all evaluated genes" | First of script 44's 26 apoptosis-related genes by interaction \|z\| (2.68; *Bbc3* 2.65). Genome-wide IHW *P* is 1.0. Exploratory | "Of the 26 we examined"; no *P* in the text |
| "the resulting drop in PUMA ... indicating" | OXPHOS falls in both genotypes (Fig. 1f), *Bbc3* only in MYC+. Nothing here tests cause | Cite 1f; drop "resulting"; "suggesting" |

**Three facts behind the design.**
1. **The Doc's own citations already match the design.** It cites g for the 6-week balance (the
   nine ratios) and h for *Bax*, *Bbc3*, *Bcl2l1* and *Foxo3* (gene level). 2G (alt) draws no
   6-week effect, so it cannot carry g's sentence; its content went into h instead.
2. **The 6W PUMA:BCL-XL shift is 62% BCL-XL falling** (3a.3). "MYC shifted the balance" is
   accurate, and h shows BCL-XL lower in MYC+ at 6W.
3. **The coupling's MYC+ slope is largely the 6W-12W group difference.** It goes from +0.56 to
   +0.26 (*P* 0.38) with age in the model (verification note, section 5). This is why Claude
   suggested Extended Data for it. The author keeps it (2a.8).

The 2d.6 flags are addressed by the proposals: item 3 (*Foxo3*) is drawn in 1h, item 4
("distinctly") is dropped from the draft, and item 5 (the coupling impasse) prints both fits and
both *P* pairs.

### 2a.8 The author's rulings and questions (afternoon)

- **"I need to see your proposals, difficult to make decisions just by reading them."**
  - Proposals are now rendered at their printed size and shown on page mock-ups.
  - This is a rule in `FIGURE_RULES.md` and a memory.
- **Keep the coupling**, "that is an important part, even if only a suggestion".
- **Keep the draft paragraph for now** (2a.9), to be edited once the graphs have been seen.
- **Character size.** The author asked whether Claude was working at the minimum size. For the
  initial submission, characters "throughout the figure can be smaller" if that gets more
  convincing content in.
  - Claude's answer: the panels are at the house sizes (6 pt labels, 7 pt titles), and Layout B
    shows the whole figure at 5/6 pt.
  - The mock-ups show that panel area, not type, is the limit.
  - Below 5 pt was not done. It reverses the 2026-10-09 rule ("nothing under 5 pt, even for the
    initial submission"), and Nature requires 5 pt at acceptance.
  - **Open: the author's word.**

### 2a.9 The draft paragraph (accepted for now), and a sentence for the coupling

> At 6 weeks, MYC raised all nine pro- to anti-apoptotic transcript ratios of the BCL-2 family
> (Fig. 1g). By 12 weeks, MYC shifted the other eight by less than half as much, but no longer
> shifted the PUMA:BCL-XL ratio at all (Fig. 1g). The change was in PUMA. Of the BCL-2 family,
> *Bbc3* fell furthest between 6 and 12 weeks in the MYC+ gland and did not change in the WT
> gland, whereas *Bcl2l1* changed in neither (Fig. 1h). *Bax* fell in both, more in the MYC+
> gland. *Foxo3*, a p53-independent activator of PUMA, rose across the window in the WT gland but
> not in the MYC+ gland, separating the genotypes more sharply than any of the 26
> apoptosis-related genes we examined (Fig. 1h). OXPHOS was deprioritised in both genotypes
> (Fig. 1f), yet PUMA fell only in the presence of MYC. This suggests that the developmental
> withdrawal of respiration permits, rather than causes, the loss of PUMA. Together, these
> findings suggested that MYC-dependent apoptosis via PUMA relies on high mitochondrial
> respiration, a state withdrawn between puberty and adulthood.

**A sentence for the coupling.** Not yet offered to the author.

> Across individual glands, the relation between OXPHOS priority and the PUMA:BCL-XL ratio
> differed between the genotypes: positive in MYC+ glands, negative or absent in WT glands.

- WT is negative when unadjusted (−0.61, *P* 0.053) and absent when adjusted (+0.05).
- "Only where MYC was present" stays withdrawn (3b).
- The legend carries both fits, both *P* pairs and section 5's age caveat.

Every number in the draft is asserted in a proposal script (2a.10).

### 2a.10 What was built: four proposal panels and three page layouts

The scripts are in `figures/natmetab/proposals/`; the PDFs are in `outputs/natmetab/proposals/`.

| proposal | script | PDF (mm) | content |
|---|---|---|---|
| 1g | `p1g_bcl2_balance.R` | `Fig1g_bcl2_balance.pdf` (56 x 56) | the ruling-4 2G redrawn from script 54: the nine ratios, MYC+ vs WT at 6W against 12W; the line through the other eight (slope 0.43) and its 95% band in full. Grey = the eight that fit the line; black = *Bbc3:Bcl2l1*. It was open in the old panel, but open means "not significant" in Fig. 1e |
| 1h | `p1h_puma_members.R` | `Fig1h_puma_members.pdf` (72 x 62) | *Bbc3* and *Bcl2l1* on top, *Bax* and *Foxo3* below, every animal over a box, genotype-major. Lower brackets: 12W vs 6W in each genotype (DESeq2, IHW-adjusted *P*). Upper bracket: the difference between the two (interaction, raw *P*). The 6W *Bbc3* difference is not bracketed (*P* 0.064, 3a.3) |
| 1h, alternative | `p1h_alt_bcl2_plane.R` | `Fig1h_alt_bcl2_plane.pdf` (60 x 60) | 2G (alt)'s plane: the BCL-2 family plus *Foxo3*, 12W vs 6W in WT against MYC+. Dashed = the same change in both; dotted = 0.20 either side (the declared threshold); black = the four genes the text names |
| coupling | `p1i_coupling.R` | `Fig1i_coupling.pdf` (62 x 56) | script 54's two fits, unadjusted and adjusted for epithelial and immune content, one line per genotype; each half prints the slope difference, *P* and permutation *P* |

**The page mock-ups.** `page_mock.R` builds `Fig1_gh_proposals.pdf`.
- **How it works:** it sources every panel script as a dry run with `natmetab.capture`, and places
  each panel at its printed size on an A4 sheet.
- **What it shows:** the 183 x 170 mm page is outlined, and the bench row a-c is a grey placeholder
  50-52 mm tall.
- **Re-run it** after any change: `Rscript figures/natmetab/proposals/page_mock.R`.

| page | layout | letters after c | verdict |
|---|---|---|---|
| 1 | **A**: 1e to Extended Data. Row 2: d (100) + the window panel (80). Row 3: ratios / genes / coupling at 50 / 70 / 57 | d, e window, f ratios, g genes, h coupling | reads; the window panel's two x axes collide at 80 mm (fixable) |
| 2 | **A2**: 1d to Extended Data. Row 2: 1e (80) + the window panel at its own 100 | d reallocation, e window, f, g, h | **the only layout where every panel keeps its own design at house size** |
| 3 | **B**: every panel kept, 5/6 pt. Row 2: d / e / f at 66 / 55 / 56. Row 3: 52 / 70 / 55 | d-i, coupling = i | d-f too narrow: e and f's labels collide and their data are squeezed; they would need redesigning, not shrinking |
| 4 | the four proposals at 100%, house size | | |

**The page arithmetic.** The bench row (about 50 mm), two panel rows (about 55 mm each) and the
gaps fill the 170 mm. Six panels therefore need two rows of three, so either one panel leaves
(A, A2) or all six narrow (B).

**Asserted in the proposal scripts**, so a re-run that changes any of these stops the panel:
- **1g:**
  - the drawn line and band are script 54's;
  - all nine 6W shifts are above 0, and the slope is below 0.5;
  - PUMA:BCL-XL's 12W effect is within 0.1 of 0;
  - PUMA:BCL-XL is the largest departure, below the line and inside the band.
- **1h:**
  - the values are script 54's `$four_genes`, and the drawn group means are within 0.1 of DESeq2;
  - *Bbc3*: IHW *P* below 0.05 in MYC+ and above 0.5 in WT; interaction raw *P* below 0.01;
  - *Bcl2l1*: IHW *P* above 0.5 in both;
  - *Bax*: IHW *P* below 0.05 in MYC+, falling less than *Bbc3*;
  - *Foxo3*: IHW *P* below 0.05 in WT and above 0.5 in MYC+;
  - *Bbc3* falls furthest of the BCL-2 family in MYC+ (script 42);
  - *Foxo3* has the largest |interaction z| of script 44's 26.
- **Alternative 1h:** the values are script 42's `$machinery`.
- **Coupling:**
  - the drawn lines are the reported fits, and their differences are the interactions;
  - the impasse holds: the unadjusted permutation *P* is below 0.05 and the adjusted one is not.

### 2a.11 Before a proposal becomes a panel of record

- **Move it** into `figures/natmetab/` under its final letter, which depends on the layout. Drop
  `dir =`, and give it the house legend printout and a SANDBOX block.
- **Layout A only:** the window panel at 80 mm needs a coarser x step.
- **1h:** the group labels are at 5 pt, because at 6 pt "MYC+ 6W" and "MYC+ 12W" touch at 72 mm.
  Keep 5 pt or widen the panel.
- **Alternative 1h:** the *Pmaip1* label sits beside *Bcl2l1*'s square.
- **1g:** the author could restore the 6-week significance marks (five of the nine have unadjusted
  *P* < 0.01). Ruling 4 removed them, and that is reversible on the author's word. If restored,
  the adjustment comes from script 54, not the panel.
- **The panel that leaves** under A (1e) or A2 (1d) needs an Extended Data home, and its sentence
  needs the new citation.

### 2a.12 Style changes (in `6cb04f1`)

- **`_style.R`:**
  - `NM_TXT` and `NM_TXT_TITLE` now read `natmetab.txt` and `natmetab.txt_title`. The defaults are
    6 and 7, the 5-pt floor is unchanged, and the title is asserted at 7 pt or less.
  - Under `natmetab.capture`, `save_panel()` keeps the drawable panel in `.nm_captured`.
  - The self-test passes 20 of 20, and all 14 panels pass a dry run.
- **`FIGURE_RULES.md`:**
  - the soft hyphen: `"­"` prints a true hyphen from base `pdf()`, so "BCL-XL" prints
    correctly;
  - the two switches;
  - "show options, don't describe them";
  - the proposals;
  - the two pending items (the layout, and text below 5 pt).

### 2a.13 Numbers for the legends (printed by the proposal scripts)

- **1g:**
  - Line over the other eight: slope 0.429 (SE 0.137).
  - PUMA:BCL-XL: +0.674 to −0.061; residual −0.350 (−1.58 SD); 95% prediction interval −0.276
    to +0.855.
  - The six-week *P* values (per-animal OLS, script 42):

    | ratio | *P* |
    |---|---|
    | *Bax:Bcl2l1* | 0.00012 |
    | *Bid:Bcl2l1* | 0.0036 |
    | *Bak1:Bcl2l1* | 0.0052 |
    | *Bax:Mcl1* | 0.0059 |
    | *Bbc3:Bcl2l1* | 0.0089 |
    | *Bbc3:Mcl1* | 0.12 |
    | *Bmf:Bcl2l1* | 0.40 |
    | *Bcl2l11:Bcl2l1* | 0.46 |
    | *Pmaip1:Bcl2l1* | 0.71 |

- **1h:** DESeq2 Wald test on raw (unshrunken) MLE log2 fold changes. The 12W-vs-6W changes carry
  IHW-adjusted *P*; the interaction carries raw *P*.

  | gene | WT | MYC+ | interaction | licence |
  |---|---|---|---|---|
  | *Bbc3* | +0.062 (*P* 0.88) | −0.479 (*P* 0.016) | −0.541 (SE 0.204), *P* 0.0081 | pre-specified |
  | *Bcl2l1* | −0.107 (0.79) | +0.073 (0.83) | +0.179, *P* 0.40 | the pair's denominator |
  | *Bax* | −0.145 (0.49) | −0.376 (0.0077) | −0.231, *P* 0.12 | exploratory |
  | *Foxo3* | +0.473 (0.0078) | −0.013 (0.98) | −0.486, *P* 0.0074 | exploratory |

  - The BCL-2 family across the MYC+ window, in order: *Bbc3* −0.48, *Bax* −0.38, *Bid* −0.36,
    *Bak1* −0.28, *Mcl1* +0.03, *Bcl2l1* +0.07, *Pmaip1* +0.21, *Bcl2l11* +0.22, *Bcl2* +0.34,
    *Bmf* +0.84.
  - Outside the family (script 42's roster): *Htra2* −0.67, *Birc5* −0.93.
  - *Foxo3*'s |interaction z| is 2.68, first of the 26 (*Bbc3* 2.65).
- **Coupling:**

  | fit | WT slope | MYC+ slope | difference | *P* | permutation *P* |
  |---|---|---|---|---|---|
  | unadjusted | −0.61 | +0.56 | +1.17 | 0.0046 | 0.026 |
  | adjusted | +0.05 | +0.83 | +0.78 | 0.0052 | 0.085 |

### 2026-10-13: the arguments behind Fig. 1g-h

### 2a.14 The day, in order

1. **Htra2:** why it is not the cause, though it falls further than *Bbc3* along the MYC+
   timeline (2a.15).
2. **The balance, made precise:** a pro-death tilt at 6W, its common shrink with MYC's fade, and
   an added effect (2a.16).
3. **A new proposal:** the apoptotic genes ranked by their departure from the fade (2a.17).
4. **A new proposal:** the three steps as rows, set beside the 1g plane (2a.18).
5. ***Bbc3* against *Bmf*:** and *Bmf* added to 1h as the developmental control (2a.19).
6. ***Apaf1*** (2a.20).
7. **Where Fig. 1g-h stands**, with a suggested paragraph (2a.21). The numbers checked today are
   in 2a.22.

### 2a.15 Htra2 is the internal control (script 42, `single_gene_null`)

**The question** (author): *Htra2* falls further than *Bbc3* along the MYC+ timeline (−0.67
against −0.48), so why is it not the cause of the 12W shift?

**The argument was already in the repo.** It is in script 42 (`single_gene_null` and its closing
note) and in `docs/2026-07-25_death_narrative_dose_vs_competence.md`, section 3, point 4.
- **The reasoning:** MYC's whole programme fades, so a gene MYC induced strongly at 6W must fall a
  lot anyway.
- **The test:** each gene's loss is sign(LFC6) x (LFC6 − LFC12), compared with genes matched on
  both baseMean and the size of the 6W effect.

| gene | MYC effect 6W → 12W | loss | matched genes lose (median) | percentile | empirical *P* |
|---|---|---|---|---|---|
| *Htra2* | +1.49 → +0.94 | 0.546 | 0.598 | **36.8** | 0.63 |
| *Bbc3* | +0.26 → −0.28 | 0.541 | 0.116 | **94.6** | 0.054 |
| *Bax* | +0.48 → +0.25 | 0.231 | 0.217 | 51.4 | 0.49 |
| *Bcl2l1* | −0.42 → −0.24 | 0.179 | 0.243 | 34.6 | 0.65 |

- **The two losses are the same size but mean opposite things.** *Htra2* loses slightly less than
  genes it is matched to; *Bbc3* loses about four times as much, enough to change sign.
- **Script 44's departure scan agrees.** *Htra2* is at the 87.5th percentile: it keeps 0.63 of its
  6W effect, against the overall 0.487. *Bbc3* is at the 3.7th and *Foxo3* at the 0.9th.
- ***Htra2* is still raised at 12W:** +0.94 log2 above WT (IHW *P* 0.012). The adult MYC+ gland
  still carries more of it than the WT gland does.
- **Consequence:** "steepest drop" is the wrong measure, and the right one is the departure from
  the fade. On the raw MYC+ timeline *Htra2* (−0.67) and *Birc5* (−0.93) fall further than
  *Bbc3*, so "the steepest drop was in *Bbc3*" must go (2a.7).

### 2a.16 The balance, made precise: tilt, common shrink, added effect

**The author's question:** "If both pro and anti-apoptotic shrinks by the same factor the balance
would not change, I assume?"

**The answer.** What fades is MYC's *effect* on each gene (MYC+ minus WT), not the levels. When
every effect halves, the MYC+ gland moves halfway back toward WT on every gene, and the ratio moves
halfway back too.

| | 6W | 12W |
|---|---|---|
| *Bax*, MYC+ vs WT | +0.48 | +0.25 (about half) |
| *Bcl2l1*, MYC+ vs WT | −0.42 | −0.24 (about half) |
| *Bax:Bcl2l1* | **+0.90** | **+0.50** (the tilt halves) |

A common fade preserves the composition of the balance, not its size: every ratio keeps the same
fraction of its tilt. This is "rescaled, not reshaped", which the text already says of the whole
MYC programme (ED 1b,c).

**The three-step argument** (agreed, 2026-10-13):
1. **MYC builds a pro-death tilt at 6W.** All nine ratios rise: median +0.49 log2, five at
   *P* < 0.01.
2. **The tilt shrinks with MYC's overall fade.** The other eight ratios keep 0.43 of their tilt (the
   1g line). This is the dose part, consistent with MYC protein halving (Fig. 1b).
3. **On top of the shrink, PUMA:BCL-XL loses all of its tilt:** +0.67 to −0.06, against +0.29
   expected.
   - At the gene level, *Bbc3* departs from the fade (+0.26 to −0.28), while *Bcl2l1* fades
     normally (−0.42 to −0.24).
   - Fig. 1c fits this: at constant MYC dose (MYCER^T2) the adult gland still loses the apoptotic
     response, so something beyond dose is expected, and the PUMA departure is the transcriptomic
     candidate.

**Two expected-shrink values, never mixed:**
- 0.43 for ratios: script 54's line through the other eight.
- 0.487 for genes: script 44's rate, fitted over the genes MYC moves.

These are two of the "four rates" (3b; the verification note, section 8).

### 2a.17 New proposal: the fade ranking (`p_fade_ranking.R`)

**The author's idea.** "rank the proapoptotic and antiapoptotic proteins according to their
relative change compared to the timeline ... I guess on this plot Bbc3 should sit at the end."

**The measure.**
- It is script 44's departure, (LFC12 − 0.487 x LFC6) / sqrt(SE12² + 0.487² SE6²), re-signed by
  function. Positive means toward survival: a pro-apoptotic transcript lower, or an
  anti-apoptotic one higher, at 12W than the fade predicts.
- For script 44's reporting set it equals sign(LFC6) x `z_resid` (asserted).

**The band:** the central 95% of the same quantity over the 8,774 MYC-responsive genes
(|z| < 1.73).

**The rosters, declared before plotting:**
- **BCL-2 family:** the anti-apoptotic members, the effectors and the BH3-only proteins, at
  baseMean 20 or more. That leaves 16 genes; *Hrk*, *Bcl2a1d* and *Bcl2l10* fall below the floor.
- **Downstream:** script 42's eight execution and IAP genes.

| group | toward survival → toward death (z) |
|---|---|
| BCL-2 family | ***Bbc3* +2.49** (outside the band), *Bmf* +1.31, *Bcl2l11* +1.01, *Bak1* +0.76, *Bid* +0.49, *Bcl2l2* +0.11, *Bok* −0.10, *Bax* −0.15, *Bcl2a1b* −0.18, *Bcl2l1* −0.21, *Pmaip1* −0.23, *Bik* −0.55, *Bad* −0.55, *Bcl2* −0.81, *Bcl2a1a* −0.95, *Mcl1* −1.16 |
| downstream | ***Apaf1* +2.14** (outside), *Diablo* +0.88, *Casp9* +0.83, *Casp3* +0.04, *Birc5* −0.38, *Xiap* −0.55, *Cycs* −0.84, *Htra2* −0.86 |

**Caveats:**
- **Grey genes:** for genes MYC did not move at 6W (*Bcl2l2*, *Bok*, *Pmaip1*, *Bad*, *Bcl2a1a*,
  *Mcl1*, *Diablo*, *Casp9*, *Casp3*), the departure is computed in the proposal script. **It must
  move to a numbered script before the panel is final.**
- ***Bcl2a1a/b*** are A1 paralogs expressed mostly in immune cells, so their signal may come from
  residual immune cells.
- ***Bbc3* is nominal.** By |z| it is more extreme than 99.3% of the 8,774 genes, but 64 are more
  extreme still; its licence is pre-specification. The band marks "unusual", not "significant".

**Outputs:** `Fig1x_fade_ranking_bcl2.pdf` (62 x 60), `Fig1x_fade_ranking_all.pdf` (66 x 86) and
`Fig1_fade_ranking_sheet.pdf`. The letter is open ("x").

### 2a.18 New proposal: the three steps as rows (`p1g_alt_ratio_rows.R`), beside the plane

**What each row draws:**
- the 6W tilt (large orange dot);
- 0.43 x the tilt, the common shrink (black tick);
- the line's 95% prediction interval (grey bar, the same as the plane's band);
- the 12W value (small amber dot).

*Bbc3:Bcl2l1* is in bold because it was not used to fit the line.
`sheet_1g_compare.R` writes `Fig1g_compare_sheet.pdf`, the plane and the rows side by side at 100%.

| ratio | 6W | expected 12W | 12W | departure |
|---|---|---|---|---|
| *Bid:Bcl2l1* | +0.95 | +0.41 | +0.40 | −0.01 |
| *Bax:Bcl2l1* | +0.90 | +0.39 | +0.50 | +0.11 |
| *Bak1:Bcl2l1* | +0.69 | +0.30 | +0.18 | −0.11 |
| ***Bbc3:Bcl2l1*** | +0.67 | +0.29 | −0.06 | **−0.35** (inside its bar) |
| *Bax:Mcl1* | +0.49 | +0.21 | +0.41 | +0.20 |
| *Bbc3:Mcl1* | +0.26 | +0.11 | −0.15 | −0.26 |
| *Bmf:Bcl2l1* | +0.24 | +0.10 | −0.20 | −0.31 |
| *Pmaip1:Bcl2l1* | +0.19 | +0.08 | +0.39 | +0.31 |
| *Bcl2l11:Bcl2l1* | +0.17 | +0.07 | −0.08 | −0.15 |

**Claude's read:**
- **The rows are easier to read than the plane,** because each step has its own mark.
- **At the ratio level PUMA:BCL-XL is the largest departure, not a lone one.** *Bmf:Bcl2l1* falls
  almost as far, and *Pmaip1* overshoots.
- **Both PUMA ratios fall short,** against both denominators.
- **PUMA stands out cleanly only at the gene level,** so the rows pair with the ranking.

### 2a.19 *Bbc3* against *Bmf*; *Bmf* added to 1h as the developmental control

**The author's question:** "what makes Bbc3 a better candidate for the added effect compared to
Bmf?"

| | *Bbc3* (PUMA) | *Bmf* |
|---|---|---|
| 12W vs 6W, WT | +0.06 (IHW *P* 0.88) | **+1.02** (*P* 0.0002) |
| 12W vs 6W, MYC+ | **−0.48** (*P* 0.016) | **+0.84** (*P* 0.003) |
| interaction | −0.54 (raw *P* 0.0081) | −0.19 (*P* 0.54) |
| MYC+ vs WT, 6W → 12W | +0.26 → −0.28 (changes sign) | −0.24 → −0.43 (below WT at both ages) |
| its ratio's 6W tilt | +0.67 (*P* 0.0089) | +0.24 (*P* 0.40) |
| departure from the fade, \|z\| percentile of 8,774 | 2.49 (99.3rd) | 1.31 (88th) |

**What separates them:**
1. **MYC-specific against developmental.**
   - *Bmf* roughly doubles across the window in both genotypes and sits on the diagonal. 2G (alt)
     already used it as its developmental control.
   - In the MYC+ gland *Bmf* rises, which is toward death. Its ratio falls short only because it
     rose a little less than in WT.
2. **A tilt to lose.** PUMA:BCL-XL starts with a real tilt; Bmf:BCL-XL barely has one, because MYC
   slightly lowers *Bmf* at 6W. *Bmf:Bcl2l1* is also one of the eight ratios the line is fitted
   to, so its shortfall is partly built in. PUMA:BCL-XL's is a prediction.
3. **At the gene level only *Bbc3* is outside the band.**
4. **Both PUMA ratios fall short:** −0.35 against BCL-XL and −0.26 against MCL-1. *Bmf* is paired
   only with BCL-XL.
5. **The prior.** *Bbc3* was pre-specified from the PGC-1α experiments, and the bench figures follow
   PUMA. Context, not evidence: BMF is best known as the BH3-only sensor of detachment (anoikis),
   which fits a developmental rise.

**What cuts the other way:**
- At 12W, *Bmf* is lower in MYC+ than in WT (−0.43, raw *P* 0.045), as *Bbc3* is (−0.28,
  *P* 0.059). The difference is the trajectory.
- *Bbc3*'s own 6W induction is not significant (*P* 0.064), and 62% of PUMA:BCL-XL's tilt is BCL-XL
  falling.
- *Bbc3*'s interaction is nominal (IHW *P* 0.84 genome-wide).

**Done** (author's request): *Bmf* is in `p1h_puma_members.R` as the developmental control.
- The panel has five genes at 96 x 64 mm. Row 1: *Bbc3*, *Bcl2l1*, *Bax*. Row 2: *Bmf*, *Foxo3*.
- **Asserted:** *Bmf* is up in both genotypes (IHW *P* < 0.05 each), and on the diagonal by script
  54's rule (|interaction| < 0.20, raw *P* > 0.05).
- **Layout:** its brackets are staggered, because a *P* with an exponent is too wide to sit beside
  its neighbour (`FIGURE_RULES.md`).

### 2a.20 *Apaf1*

**No argument specific to *Apaf1* is recorded anywhere in the repo.** The author recalled one. The
search covered the docs, the scripts, the old panels and the session exports. The nearest reasoning
is general:
- **Script 17's module scan (July):** "death is primed but execution is held down". MYC lowers the
  caspase/apoptosome, cytochrome-c and SMAC/XIAP modules at 6W.
- **Fig. S2C:** the WT window does not move the death machinery.

The gene-specific internal-control argument is *Htra2*'s (2a.15).

**The case now** (checked against the current objects):
1. **It barely changes in the MYC+ gland:** −0.11 across the window (IHW *P* 0.60); interaction
   −0.13 (raw *P* 0.39).
   - Its departure from the fade (z +2.14) comes from *failing to fade*. The fade predicted its 6W
     repression (−0.24) would halve to −0.12; it went to −0.37 instead.
   - *Bbc3*, by contrast, actually falls (−0.48, *P* 0.016).
2. **It was already lower under MYC at 6W** (−0.24, IHW *P* 0.092), when MYC+ glands show the most
   apoptosis. Low APAF1 did not stop death then.
3. **Its module moves the other way.**

   | module (script 17) | MYC+ vs WT, 6W → 12W | interaction |
   |---|---|---|
   | caspase/apoptosome (Reactome, 6 genes) | −0.55 → −0.29 | +0.26 (*P* 0.28) |
   | cytochrome-c release | −0.41 → −0.15 | +0.26 |
   | SMAC/XIAP | −0.45 → −0.20 | +0.25 |

   All three move toward death, not survival, so *Apaf1* does not represent its module.
4. **It acts after the decision point.** APAF1 works after cytochrome c release. The release is set
   by the BCL-2 balance, which is where the bench experiments act: BCL-xL rescue, and PGC-1α
   inducing PUMA.

**Literature, from memory; verify before citing.** In Myc-driven mouse lymphoma, losing Apaf-1 or
caspase-9 did not accelerate tumour formation (Scott et al., *J Cell Biol* 2004). Losing PUMA did
(Hemann et al., *PNAS* 2004; Garrison et al., *Mol Cell Biol* 2008; Michalak et al., *Cell Death
Differ* 2009). These are lymphoma studies, not mammary.

**Caveat.** On the author's own fade logic *Apaf1* counts, because ending up lower than the fade
predicts is a departure. Claude's suggestion is to keep it visible in the ranking, with a legend
note that it departs by failing to fade, not by falling.

### 2a.21 Where Fig. 1g-h stands

**The argument is agreed** (2a.16). The author will come back to finalise, and "will continue with
ideas what to add".

**The candidates:**
- **1g:** the plane (`p1g_bcl2_balance.R`) or the three-step rows (`p1g_alt_ratio_rows.R`).
- **The added effect at gene level:** the fade ranking (`p_fade_ranking.R`), the BCL-2 family alone
  or with the downstream genes.
- **1h:** five genes in every gland (`p1h_puma_members.R`) or the BCL-2 plane
  (`p1h_alt_bcl2_plane.R`).
- **The coupling** (`p1i_coupling.R`, kept by the author).

**Claude's suggested set:** the rows, then the ranking, then the five-gene panel, then the coupling.
That puts the tilt and the shrink first, then which gene departs, then the animals behind it, then
the link to OXPHOS.

**The page.**
- This set exceeds even yesterday's layouts (2a.10), so the layout has to be redone for the chosen
  set.
- **`page_mock.R` still places the four-gene 1h at 70 x 54 mm.** Resize it before re-running, or
  the five-gene panel is squeezed.

**A suggested paragraph along the three steps.** It is not yet seen by the author, and it keeps
2a.9's wording where it can. Fig. 1x is the ranking, with its letter open.

> At 6 weeks, MYC raised all nine pro- to anti-apoptotic transcript ratios of the BCL-2 family
> (Fig. 1g). By 12 weeks, as the MYC programme faded to about half its amplitude (ED Fig. 1b), MYC
> shifted eight of these ratios by less than half as much, but no longer shifted the PUMA:BCL-XL
> ratio at all (Fig. 1g). The difference lay in PUMA. Of the BCL-2 family, only *Bbc3* departed from
> MYC's overall fade beyond the range of 95% of MYC-responsive genes (Fig. 1x). It fell between 6
> and 12 weeks in the MYC+ gland and not in the WT gland, whereas *Bcl2l1* changed in neither and
> *Bmf* rose in both (Fig. 1h). *Foxo3*, a p53-independent activator of PUMA, rose across the window
> in the WT gland but not in the MYC+ gland, separating the genotypes more sharply than any of the
> 26 apoptosis-related genes we examined (Fig. 1h). OXPHOS was deprioritised in both genotypes
> (Fig. 1f), yet PUMA fell only in the presence of MYC, suggesting that the developmental
> withdrawal of respiration permits, rather than causes, the loss of PUMA. [The coupling sentence,
> 2a.9.] Together, these findings suggested that MYC-dependent apoptosis via PUMA relies on high
> mitochondrial respiration, a state withdrawn between puberty and adulthood.

- **If the ranking shows the downstream genes,** an optional clause: "among the genes acting
  downstream of the BCL-2 family, only *Apaf1* also departed, by keeping its repression rather than
  losing it".
- **Every number in the paragraph is asserted** in the proposal scripts (2a.10, 2a.17-2a.19).

### 2a.22 Checked today against the saved objects (read-only; no object written)

- **Script 42, `$single_gene_null` and `$machinery`:** the *Htra2* control (2a.15); *Apaf1* and
  *Bmf* on all contrasts.
- **Script 44, `$collapse_genes` and `$mech_genes`:**
  - the fade rate 0.4872;
  - the band: the 95th percentile of |departure z| is 1.73;
  - percentiles: *Bbc3* 0.51 (aligned) and 99.3 (|z|), *Bmf* 88.2 (|z|).
- **Script 17, `gsva_overview.rds$coef_table`:** the three execution modules (2a.20).
- **Script 54, `$ratios` and `$ratio_line`:** the rows (2a.18).
- **Script 03, `interaction_results.rds`:**
  - *Apaf1* interaction raw *P* 0.39;
  - *Bmf*: 6W −0.243 (SE 0.215), 12W −0.431 (raw *P* 0.045), interaction −0.188 (SE 0.304,
    *P* 0.54).

### 2026-10-14: the 247-mm page, two rounds of proposals, scripts 55 and 56

### 2a.23 The day, in order

1. **The author's plan** (typed in a pasted block):
   - the page is **183 x 247 mm for the initial submission** (ED 180 x 247), and "always use the
     minimum font size where font size becomes limiting" (5 pt; nothing below it);
   - **1g:** a heatmap of z for the BCL-2 family, as 2H+I (alt)'s top panel;
   - **ED 3b:** the ranking, version B, coloured by a graded 6W MYC effect, axis "Difference from
     the overall MYC effect";
   - **1h:** the simplified ratio rows, "aligned with the common shrink", ordered by the 6W-12W
     difference;
   - **1i:** the adjusted coupling, the WT line dashed and thinner, the MYC+ line's R2 on the panel;
   - **the FOXO3 argument goes to Extended Data:** a dot plot of the S2D (alt) roster plus *Bbc3*,
     *P* < 0.05 labelled, same axis name. ("Foxo1/Foxo1a" in the plan was read as *Foxo3*, FOXO3a:
     *Foxo1* is a control in that roster);
   - **ED 2c:** wider, with PC1 broken for the outlier as in the old Fig. 1C.
2. **The author's answers:** "1. the expected" (1h draws the expected 12W against the observed);
   "3. yes, label them all" (the four regulators at *P* < 0.05); "4. keep them, will discuss later"
   (the ED 3a and 3d additions); "5. write it, I will run it" (script 55); "place the pa values in
   the legend" (1i's *P* values).
3. **Script 55** written; a sign bug on *Diablo* found and fixed; the author re-ran it (2a.24).
4. **Round 1** built and laid out on the 247-mm page (2a.27).
5. **The author's six notes** on the mock-ups; **round 2** built (2a.28).
6. **Script 56** written for the two statistics round 2 needed, dry-run, and the two proposals
   switched to read it (2a.29). Committed `4fb3aa9` and pushed, at the author's request.

### 2a.24 Script 55: the departure from MYC's fade, for every expressed gene

`scripts/55_departure_from_myc_fade.R` writes `results/departure_from_myc_fade.rds`: `$genes`
(15,191 expressed genes), `$band`, `$rosters`, `$defs` (with `$text`, the numbers the text quotes,
asserted where computed).
- **The statistic is script 44's, extended to the genes MYC did not move at 6W:** departure =
  LFC12 - 0.4872 x LFC6, SE from both contrasts, z, nominal *P*. Unsigned; the BCL-2 family and the
  downstream genes also get a function sign (+ = toward survival).
- **Bands** (95th percentile of |z|): **1.73** over the 8,774 responsive genes (script 44's set,
  rebuilt and asserted), **1.71** over the unmoved ones.
- **The interaction is not the departure.** Interaction = fade share + departure:
  *Bbc3* -0.541 = -0.132 + **-0.408** (mostly departure, z -2.49); *Htra2* -0.546 = **-0.764** +
  0.217 (all fade, z +0.86); *Bax* -0.231 = -0.248 + 0.017 (it only fades).
- **Rosters:** BCL-2 family 16 read (roles: effectors, activators BID/BIM/PUMA, sensitisers,
  anti-apoptotic); downstream 8; PUMA regulators 20 read (*Trp73* and *Cdkn2a* below the floor).
- **Findings:** *Bbc3* is the only BCL-2 member outside its band (+2.49 toward survival); *Apaf1*
  the only downstream one (+2.14). Among the PUMA regulators, **four** pass raw interaction
  *P* < 0.05, not two: *Foxo3* (z -2.53, *P* 0.0074) and *Bbc3* (-2.49, 0.0081) below the fade,
  the p53 targets *Zmat3* (+2.74, 0.046) and *Eda2r* (+2.74, 0.014) above it, through a WT fall with
  the MYC+ gland flat. None survives BH across the 20.
- **The Diablo fix:** the first version signed *Diablo* (SMAC, an IAP antagonist) anti-apoptotic.
  It now uses the `^brake` rule and asserts the anti-apoptotic downstream set is {*Xiap*, *Birc5*};
  *Diablo*'s function z is +0.876. The author's re-run carries it.

### 2a.25 What the data say against the draft text (flagged; nothing written to the Doc)

- **"progressively"** spans two ages; say "between 6 and 12 weeks".
- ***Bax*** falls under MYC only by the common fade (departure +0.02); *Bbc3* falls beyond it. Do not
  pair them as one loss.
- **"the steepest drop was in *Bbc3*"** is false on raw values (*Htra2* and *Birc5* drop further).
  Say "fell furthest beyond MYC's overall fade".
- **At 6W only two BCL-2 genes are significant on their own:** *Bax* up (IHW 1.4e-4) and *Bcl2l1*
  down (IHW 0.04). "All nine ratios rise" is the ratio statement (median +0.49), not nine genes.
- **New from round 2 (the merged 1f, 2a.28):**
  - "both glands **demoted OXPHOS in priority**" and "FAO and the carnitine shuttle **gained
    priority**": no window change of these pathways passes script 08's BH test (OXPHOS adjusted *P*
    0.46 WT, 0.17 MYC+), so every priority point is open. The claims rest on ranks among 143:
    OXPHOS 14th and 13th lowest; carnitine shuttle 2nd and 3rd highest; FAO 8th and 25th. Say
    "ranked among".
  - "OXPHOS priority **falls on both timelines**" is not licensed: in ED 3a the WT window falls as
    far but *P* = 0.086 (script 56 asserts it).
  - The proliferation signature's WT point is filled: it fell only -0.05, but beyond 98% of
    matched sets (731 genes give a narrow null). "Fell further" is the bracket's job (*P* < 0.0005).

### 2a.26 The coupling's two *P* values (the author's question 2)

- The parametric *P* (0.0052) tests the slope difference against 0.
- The within-age permutation shuffles OXPHOS within each age, which keeps the age-level
  co-movement, so its null median is **+0.40**, not 0. The observed +0.78 sits at the 91.5th
  percentile (*P* 0.085): **suggestive only**.
- With age and the covariates in the model, the MYC+ slope is +0.62 (*P* 0.043).
- The author: both *P* values go to 1i's legend.

### 2a.27 Round 1: the proposals on the 247-mm page

| slot | script | what |
|---|---|---|
| 1g | `p1g_heatmap.R` | 16 BCL-2 genes by role; group means of VST z; espresso-white-mint, +/-2; 5-pt numbers. MYC+ 6W: effectors and activators up, BCL-XL down; *Bbc3* +0.8 at 6W, -0.9 at 12W |
| 1h | `p1h_ratio_expected.R` | each ratio's expected 12W (0.4291 x 6W, the line through the other eight; SE 0.137) against observed, rows by gap. PUMA:BCL-XL: expected +0.289, observed -0.061, inside the 95% prediction interval [-0.276, +0.855], the largest gap |
| 1i | `p1i_coupling_adjusted.R` | adjusted fit only; WT dashed 0.35; MYC+ R2 0.73 (squared r of the drawn points; 0.35 unadjusted); *P* values in the legend |
| ED 2c | `ped2c_pathway_axis_wide.R` | PCA with PC1 broken before the one 6W MYC+ animal at +44, loadings beside |
| ED 3b | `ped3b_fade_ranking.R` | reads script 55; fill = 6W LFC; *Bbc3* alone outside in the BCL-2 family |
| ED 3c | `ped3c_puma_regulators.R` | unsigned z, shapes by role, the four labelled with raw interaction *P* |
| ED 3a, 3d | `p1h_puma_members.R`, `p1i_coupling.R` | the five genes per animal; both coupling fits |

### 2a.28 Round 2: the author's six notes, and what was done

1. **1h:** "The X-axis title is misleading": now "MYC+ vs WT (normalised to overall MYC activity)",
   key "6W" and "12W". **The legend must say the factor is 0.43, the common shrink of the other
   eight ratios** (the same per-animal estimator; script 54's ruling 4), not the gene-level fade
   0.49. With 0.49 the row order is identical and PUMA:BCL-XL's expected is +0.33.
2. **ED 2b-c colours** ("distinguish the mitochondrial gene sets"): `mito_class_cols` is now
   reddish purple for the two mitochondrial classes (`#8E2C74` MitoCarta, `#CC79A7` curated with
   mitochondria) against grey (`#D9D9D9`). `options(natmetab.mito_class_cols = ...)` tries others.
3. **ED 2c lower, more labels:** 46 mm (was 58); the PCA data area is fixed at 30 mm so PC1 and PC2
   have equal units; the loadings name OXPHOS (rank 114, 0.98), MYC (Felsher) (304, 0.94; the
   Felsher integrative signature), MYC V1 (358, 0.92), Cell cycle (Reactome, 409), E2F targets (450)
   and Glutathione (585). **"MYC V1" is MSigDB's MYC_UP.V1_UP: the library has no Hallmark MYC
   targets V1, only V2** (rank 334, 0.93). The top 100 sets are 94 TF regulons x MitoCarta, 5
   MitoCarta, 1 other (printed).
4. **ED 2d and 2e the same size:** every scatter's data area is a fixed 44-mm square on one shared
   range (`fixed_size()`; `stack_parts()` with grid units). Panels of record, edited in place.
5. **1f and ED 2f merged** (`p1f_window_combined.R`, 100 x 56): 1f's rows and two rulers in one
   panel (no facets, so a bracket can cross the groups); **filled = *P* < 0.05** -- for expression,
   beyond the central 95% of 2,000 expression-matched sets (script 54; the two rows it lacks from
   script 56); for priority, script 08's BH, Fig. 1e's grammar, so every priority point is open;
   **no median ticks**; **a bracket for OXPHOS minus proliferation, WT, *P* < 0.0005** (script 43).
   ED 2f leaves; the old 2g becomes 2f.
6. **ED 3a per animal** (`ped3a_per_animal.R`, 104 x 60): OXPHOS (mitoPPS), *Bbc3:Bcl2l1*, *Bax*,
   *Bmf*, z over the 24 animals, Fig. 2H+I (alt)'s form; lower brackets 12W vs 6W per genotype, upper
   the interaction; raw *P* throughout (the old version's lower brackets were IHW-adjusted).
7. **ED 3b without *Apaf1*:** the omission is asserted (outside the band, interaction *P* > 0.05)
   and printed with its numbers for the legend (2a.30).

### 2a.29 Script 56: the two statistics the round-2 proposals needed

`scripts/56_window_and_chain_statistics.R` writes `results/window_and_chain_statistics.rds`:
`$window_null`, `$window_null_control`, `$values` (the 24 animals' four quantities), `$tests`,
`$welch_08`, `$defs`, `$notes`. **Dry run only (save sent to the scratchpad); the author has not run
it.** Seeded (56), so the author's run reproduces these numbers.
- **PART A, Fig. 1f's two extra rows**, by script 54's recipe (20 bins, 2,000 draws), proved on three
  of script 54's arms (OXPHOS assembly, TCA cycle, nucleotide metabolism: within **1.45** percentile
  points; tolerance 5) and against script 40's ruler:
  - carnitine shuttle (5 genes): WT +0.200, 85.5th percentile; MYC+ +0.104, 76.4th: within;
  - fatty acid oxidation (39): WT +0.227, **100th** (beyond every matched set); MYC+ -0.020, 75.6th.
- **PART B, ED 3a's tests** (12W vs 6W in WT | in MYC+ | interaction; raw *P*):

  | quantity | instrument | WT | MYC+ | interaction |
  |---|---|---|---|---|
  | OXPHOS mitoPPS | lm(value ~ age x genotype) | -0.109, *P* 0.086 | -0.153, *P* 0.020 | -0.044, *P* 0.62 |
  | PUMA:BCL-XL | the same, script 42's model | +0.180, *P* 0.45 | -0.555, *P* 0.027 | -0.735, *P* 0.037 |
  | *Bax* | DESeq2 Wald | -0.145, *P* 0.17 | -0.376, *P* 3.3e-4 | -0.231, *P* 0.12 |
  | *Bmf* | DESeq2 Wald | +1.023, *P* 1.8e-6 | +0.835, *P* 1.1e-4 | -0.188, *P* 0.54 |

  The ratio's fit reproduces script 42's saved values to 1e-8. The OXPHOS score is script 08's
  (window differences to 1e-6), whose Welch test gives *P* 0.12 and 0.040 (BH 0.46, 0.17): Fig. 1f's
  legend quotes Welch, ED 3a draws the model, and each legend names its test.
- **The mock-up** `outputs/natmetab/proposals/Fig1_ED2_ED3_247.pdf` was drawn with these numbers:
  first inside the two proposals (same recipe, another seed: no point changed state), then verified
  against script 56's dry-run object in the scratchpad.

### 2a.30 Legends: what each must say (every script prints its numbers)

- **1f (with ED 2f):** filled = *P* < 0.05 against 2,000 expression-matched sets (two-sided
  empirical; script 54, two rows script 56) for expression, BH-adjusted Welch (script 08) for
  priority; priority claims are ranks among 143; the bracket is the WT paired comparison (script 43),
  with OXPHOS minus mitochondrial ribosome (*P* 0.0005) and minus nucleotide metabolism (*P* 0.0085)
  in words; each window coincides with an extraction batch; the OXPHOS protein decline is Fig. 1b.
- **1h:** the 6W points are scaled by 0.43, the common shrink of the other eight ratios; PUMA:BCL-XL
  is inside the 95% prediction interval (no row is "significant").
- **1i:** slope difference +0.78, *P* 0.0052 parametric and 0.085 against the within-age
  permutation; MYC+ R2 is of the drawn (adjusted) points.
- **ED 2c:** the class colours are ED 2b's key; MYC V1 is MYC_UP.V1_UP; glutathione is the redox
  control; the breast-cancer MYC bicluster loads as high as OXPHOS (0.98).
- **ED 3a:** the instruments per quantity (2a.29), raw *P*; within-genotype brackets coincide with
  batch, the interaction does not.
- **ED 3b:** "*Apaf1* (+2.14 on the panel's axis, outside the band) is not shown: MYC had already lowered it at 6W
  (-0.24), its change across the window does not differ by genotype (interaction *P* 0.39), and it
  acts after the BCL-2 decision" -- or the author's wording; the omission must be visible.

### 2a.31 Where things stand

**Open for the author** (round 2 is not yet reviewed):
- the merged 1f (and so ED 2's renumbering);
- ED 3a's new form, and ED 3d's place ("keep them, will discuss later");
- the purple class colours; ED 2c's six labels, and whether the wide ED 2c replaces the panel of
  record (`edfig2c_pathway_axis.R`);
- the final panel sizes;
- the text fixes (2a.25).

**Renumbering in the Doc when the set is fixed:** the new ED 3 pushes the bench ED Figs 3 and 4 down
one; ED 2f leaves (into 1f) and ED 2g becomes 2f; the FOXO3 regulators are ED 3c.

**The layouts** (`page_mock_247.R`, all checked by eye):
- **Fig. 1** (247 of 247 mm, nothing to spare): bench row a-c 56 (placeholder); d 100 + e 80 (66);
  f 100 + g 80 (56); h 100 + i 80 (60); 3-mm gaps. *Superseded on 2026-10-15 by 2a.34.*
- **ED 2** (about 215 mm): a 46 + b 58 (40); c 176 x 46; d 114 + e 60 (56); f (old g) 84 x 64.
- **ED 3** (about 165 mm, room to spare): a 104 + d 73 (60); b 86 x 100 + c 86 x 92.

### 2a.32 The evening: the author's notes after the first handoff, and the promotion

1. **ED 3b's fill scale** "should be narrowed, the differences are not visible": +/-0.5 (was
   +/-1.5); *Bid* (+0.61), *Birc5* (+1.16) and *Htra2* (+1.49) are drawn at the end colour, and the
   legend names them.
2. **ED 3d, "adjusted to epithelial cells? ... why adjust to that?"** -- answered: the adjustment is
   for sample purity, not for the MECs. Script 42's two RNA composites (epithelial: *Epcam*,
   *Krt8*, *Krt18*, *Krt5*, *Krt14*, *Cdh1*, *Krt17*; immune: *Ptprc*, *Cd52*, *Cd3e*, *Lyz2*,
   *Cd74*, *H2-Aa*, *Itgam*, *Ms4a1*) estimate each unsorted prep's non-epithelial share; they
   correlate -0.84 (WT) and -0.92 (MYC+), one purity axis. In WT the ratio tracks immune content
   (r +0.90), so the unadjusted WT slope (-0.61) is mostly contamination; adjusted it is +0.05, and
   MYC+ goes +0.56 -> +0.83 (script 54's decomposition: the immune score carries most of it). The
   label is now "Adjusted for immune content" (author); 1i's y axis matches ("*Bbc3:Bcl2l1* (z
   score), adjusted for immune content"). **The legend must name both scores**, and the caveat that
   the keratin score also moves with epithelial state. The WT line is dashed and thinner, as in 1i.
3. **The band in ED 3b-c** became two dashed pairs, labelled along their length: **P = 0.05**
   (|z| 1.73) and **P = 0.01** (|z| 2.38), the 95th and 99th percentiles of |z| over the 8,774
   MYC-responsive genes (author: "more important to show that the others are not significant").
   - **The P is empirical** (script 55's new `p_emp_resp`: the fraction of responsive genes
     departing at least as far, two-sided), not the nominal z-test (*Bbc3* 0.013, *Foxo3* 0.011
     nominal): departures among MYC-responsive genes are narrower than a standard normal.
   - Every drawn gene is beyond both pairs or inside both (smallest P inside: 0.10); asserted. Only
     the omitted *Apaf1* (P 0.018) would fall between them.
   - **ED 3c now prints the same empirical P** (*Foxo3* 0.0067, *Bbc3* 0.0074, *Zmat3* 0.0043,
     *Eda2r* 0.0041). The interaction P it printed before (0.0074, 0.0081, 0.046, 0.014) goes to the
     legend: beside a P = 0.01 line, *Zmat3* would have shown "P = 0.046" beyond it. Same four genes.
4. **Saved sizes = the page's slots** (1e 80 x 66, 1g 80 x 56, 1h 100 x 60, 1i 80 x 60, ED 3b
   86 x 100, ED 3c 86 x 92), so each PDF is placed at 100% exactly as reviewed.
5. **The promotion** (author: "promote them to final panels and commit"): see FIGURE_RULES.md's
   built table. New names: `fig1f_window_reallocation.R` (the merged panel; the old one is
   retired under the same name), `fig1g_bcl2_heatmap.R`, `fig1h_ratio_expected.R`,
   `fig1i_coupling.R`, `edfig2c_pathway_axis.R` (the wide one), `edfig2f_window_two_rulers.R` (the
   old 2g), `edfig3a_per_animal.R`, `edfig3b_fade_ranking.R`, `edfig3c_puma_regulators.R`,
   `edfig3d_coupling_both_fits.R`. All moved with `git mv`.
6. **Prism (the author's question):** Fig. 1d's top is a box-and-whisker plot (median and IQR,
   Tukey whiskers at 1.5 x IQR) with every animal drawn. Prism shows all points only with
   min-to-max whiskers, so R switched to min-max (author: "yes change the whiskers to min max to
   all panels with box and whiskers"): Fig. 1b, 1d and ED 3a use `geom_boxplot(coef = Inf)`, and
   `save_panel()` now refuses a box whose whiskers do not run min to max (self-test case added).
   Legend sentence for R and Prism alike: "boxes show the median and interquartile range, whiskers
   the minimum and maximum, and dots individual animals".

### 2a.33 2026-10-15: Fig. 1b in the Prism format

The author's four notes, all done in `fig1b_myc_transcript.R` (FIGURE_RULES.md, new section "The
Prism format"): 1b sits among Prism graphs, and "Prism is less flexible, so will need to adapt".
1. **No x axis and no group labels** (Illustrator sets them).
2. **Graph 20 x 20 mm**: the data area, fixed with `fixed_size()`. The y axis ends at 12,500 (a minor
   tick), with the top animal at 12,180; labels every 5,000. The brackets sit above the axis end.
3. **Axis title 6 pt, tick labels and P values 5 pt.** The title is on two lines (`atop()`, as in
   1i and ED 3b-c), because one line is 30 mm against the 20-mm axis. **Brackets have no tips.**
   The two genotype P values are 12 mm wide at 5 pt, so the brackets take four levels, about 13 mm
   above the axis end.
4. **Every line 1 pt**, read from the PDF: axis, ticks, box outlines, whiskers, medians (ggplot
   doubles the median unless `median.linewidth` is set) and brackets. At 1 pt the two WT boxes,
   whose interquartile ranges are under 1 mm on this axis, print as solid grey bars.

**Second round, the same day** (author: "lines are slightly thicker than prism - try 0.9 point?",
and less space around the panel). Page 37.5 x 34 -> **30 x 28.5 mm**:
1. **Lines 0.9 pt** (0.90 in the PDF).
2. **Brackets over the animals they compare**: the 12W genotype bracket moved down into the graph,
   over the 12W groups; the WT window (P = 0.48) left out, its P for the legend (the script prints
   it, marked "not drawn"); the MYC+ window moved down to just above the 6W bracket.
3. **Tick labels in thousands** (0, 5, 10), title "*Myc* mRNA (x1,000, / normalised counts)". The
   comma in 1,000 follows the rules (the author wrote "x1000"), and "normalised" is the author's
   spelling, as on Fig. 1h. The break after "x1,000," keeps both lines inside the 20-mm axis; one
   line would be 26 mm.
4. **The title is closer to the axis, with its lines closer together**: 2.4 mm baseline to
   baseline, 0.6 mm from the tick labels. `atop()` cannot set that, so the title is a narrow
   column of its own, stacked beside the graph (`stack_parts`, all widths in mm). The x is
   Helvetica's own glyph (byte D7), so Illustrator can edit it.

**`_style.R` gained, with defaults unchanged** (the other 19 panels re-rendered pixel-identical;
20 of 20 pass the dry run; self-test 26 of 26):
- `lw_pt()`, points to linewidth;
- `theme_nm(txt, txt_title, line)` and `brackets(txt, line, tips)`; `tips = FALSE` draws the bar
  alone, with no zero-length tip segments;
- `natural_size()`, the page a fixed-size panel needs. `save_panel()` now refuses a fixed-size
  panel saved smaller than that.

**The author re-renders it:** `source(here::here("figures", "natmetab", "fig1b_myc_transcript.R"))`.
Committed and pushed as `7bc98c3` (2026-10-15).

### 2a.34 2026-10-15, afternoon: Fig. 1 resized, in 5/6-pt type

**Fig. 1b goes to Prism** (author: "easier to get the prism lookalike if I just do it there").
- **The data:** the author was given the 24 values as a Prism column table. They are DESeq2's
  normalised *Myc* counts (`counts(dds, normalized = TRUE)` from `results/dds_int_run.rds`), the
  panel's own, one row per animal and unpaired.
- **The *P* values** are the four IHW-adjusted Wald values the script prints: 6 wk 2.7 x 10^-17,
  12 wk 2.1 x 10^-20, MYC+ window 0.76, WT window 0.48 (legend only). Prism must not test: its
  t-test or ANOVA would not give these.
- **Boxes:** Prism's quartiles differ from R's at n = 6, so its boxes are wider; the medians and
  min-max whiskers are the same.
- `fig1b_myc_transcript.R` stays as the source of the numbers.

**Fig. 1d, 100 x 60** (four rounds):
1. **Type:** tick labels, the key and *P* values 5 pt; strip and axis titles 6 pt. The exponent is
   reduced (`brackets(exponent = "reduced")`: plotmath's scriptstyle, 70%, which base `pdf()`
   rounds to 4 pt).
2. **The y title.** The author asked for "Fraction of total transcriptome (%)".
   - The values are script 32's `share_nomt`, whose denominator leaves out the 13 mtDNA-encoded
     transcripts, 7-16% of each library. So "total" would overstate every value: Nuclear MitoCarta
     is 2.2-4.6% as drawn and 1.8-4.5% of the total.
   - The title reads **"Fraction of nuclear-encoded transcriptome (%)"**; the author did not object.
   - Script 32's `share_all` would make "total" true. It would also move the text's 20-33% to
     21-33%, and the *P* values from 0.00049-0.0074 to 0.0054-0.036, because the mtDNA transcripts
     vary a lot between animals.
3. **Box-plot axes.** The author: "for box and whiskers plots it is acceptable that the y axis does
   not start at 0 ... start the y-axis at the closest tick to the whisker at the bottom, and end it
   above the highest one, in this case including the p value".
   - The axes now run 0.05-0.2, 2-5, 0.3-1.1 and 0.4-1.2 (from 0 to 0.25, 6, 1.2 and 1.5 before).
   - The ends are worked in printed mm on upper graphs fixed at 28 mm (`box_axis()`, `label_mm()`);
     the lower graph takes the rest, 11.3 mm.
   - The graphs do not clip, so the lowest OXPHOS animals (0.31, on an axis from 0.3) draw whole.
4. **Size.** The author asked for 60 mm, then 56, then 95 wide, then "revert back to 100x60". The
   67.5- and 104-mm readings had come from measuring a different file.

**Fig. 1e:** 78 x 55 (was 80 x 66), in 5/6 pt. The author chose it from two previews, against the
6/7-pt version.

**Fig. 1f:** 5/6 pt; its size is unchanged (100 x 56), and it now sits alone in its row.

**Fig. 1g, h and i: one row of 180 mm** (author), all 56 tall: 66 + 52 + 56, with two 3-mm gaps.
- **g** (was 80 wide): the key's title is "z-score" alone, with no frame. The key moved left, to
  the top of the role names' column beside Effectors and Activators, 1 mm past the wider of the
  two. That is a negative `legend.box.spacing`, worked from the names' printed widths. Its 18-mm
  bar ends above Sensitisers.
- **h** (was 100 x 60): narrowed most, by the x axis (about 37 mm). Its title takes two lines,
  because one line would be 46 mm.
- **i** (was 80 x 60): the square data area is as large as the height allows (about 41 mm). The
  *R*2 superscript is reduced like 1d's exponents; that was Claude's choice, and the author has not
  objected.

**"z-score" everywhere** (author): Fig. 1g's key, Fig. 1i's two axis titles, and the axis titles of
ED 3a and ED 3d, written with the soft hyphen. The italic *z* of ED 3b-c is a test statistic, so it
stays.

**`_style.R`, defaults unchanged.** The twelve untouched panels are byte-identical; 20 of 20 pass;
the self-test passes 33 of 33.
- `fmt_p()` and `brackets()` take `exponent`: "full" (displaystyle) or "reduced" (scriptstyle).
  `check_panel()` refuses a superscript in neither style.
- `label_mm()` gives a plotmath label's printed height (5-pt *P* labels: 1.30 mm, or 1.90 with a
  reduced exponent).
- `box_axis()` applies the box-plot axis rule in printed mm, reworking the mm on each axis it finds
  until that axis keeps its own margins.
- `check_panel()` exempts box plots from the from-0 check.
- `natural_size()` reads a null nested in a sum as 0, as grid's layout does. Patchwork centres a
  collected key that way, so a fixed stack with a shared key now measures.

**Two facts about base `pdf()`**, both in FIGURE_RULES.md and read from the PDFs:
- it rounds every font size to a whole point: 5.5 prints at 6, and a 70% exponent of 5 pt at 4;
- it writes the page size in whole points, rounded down: a 56-mm page reads 55.74 mm in Illustrator.

**The page** (`page_mock_247.R`): Fig. 1 takes 237 of 247 mm. Bench row a-c 56; d 100 x 60 + e
78 x 55 (row 60); f 100 x 56; g + h + i (56).

Commits, all pushed: `dfa8986`, `7f2d178`, `8bb8242`, and this handoff.

---

## 2b. What happened on 2026-10-10 evening to 2026-10-11

*ED Fig. 2 letters in this section are the ones in use until 2026-10-12. Renumbered since (2a.1):
2d is now 2e, 2e is now 2f, 2f is now 2g, and 2g is now 2d.*

Three commits, all pushed: `a22a149` (Fig. 1d), `08c0dda` (the 21-27% asserts removed) and
`f5dd876` (Fig. 1e-f and ED 2d-g), plus that day's handoff, `babd188`.

### 2b.1 What was built

| slot | script | PDF (mm) | content |
|---|---|---|---|
| Fig. 1d | `fig1d_mito_content.R` | 100 x 66 | four module shares, every animal, over the MYC effect on 143 MitoPathways at 6W and 12W; one PDF, one key; "20-33%" asserted |
| Fig. 1e | `fig1e_myc_reallocation.R` | 83 x 66 | mitoPPS MYC effect at 6W and 12W (12W smaller): 10 rows chosen by rule (significant at 6W, same sign at 12W) plus OXPHOS subunits |
| Fig. 1f | `fig1f_window_reallocation.R` | 100 x 56 | 12W vs 6W in WT and MYC+ (MYC+ smaller), expression and priority side by side; TEB and proliferation above six mitochondrial arms |
| ED 2d | `edfig2d_myc_reallocation_ages.R` | 62 x 62 | MYC effect at 12W against 6W over 143 pathways: slope 0.64, R2 0.79, r 0.89; Fig. 1e's rows black, unnamed |
| ED 2e | `edfig2e_window_matched_controls.R` | 95 x 66 | 9 programmes against 2,000 expression-matched sets in both windows (MYC+ smaller); below, OXPHOS minus proliferation, ribosome and nucleotide metabolism, paired, WT: P < 0.0005, 0.0005, 0.0085 |
| ED 2f | `edfig2f_window_two_rulers.R` | 88 x 66 | expression against priority over 143 pathways: r 0.79 (WT), 0.86 (MYC+) |
| ED 2g | `edfig2g_window_independent.R` | 120 x 62 | WT window against MYC's 6W effect, r -0.02; MYC+ window against WT window, r 0.69 |

ED Fig. 2 starts at d on the author's numbering; a-c are design, the mitochondrial axis and total
abundance (2d.1's slot map). Each script prints its legend's numbers.

### 2b.2 The paragraph of record (author, 2026-10-11)

> Within the nuclear-encoded mitochondrial transcriptome, scored as pairwise ratios of pathways
> independent of its overall size (mitoPPS19), MYC promoted biogenesis and biosynthetic pathways
> and demoted mitochondrial dynamics, surveillance and signalling, in the same pattern at both ages
> (Fig. 1e, ED Fig. 2d). The pubertal-adult transition, by contrast, reshaped the compartment
> independently of MYC (ED Fig. 2g). Between 6 and 12 weeks, OXPHOS subunit expression fell in both
> the WT and the MYC+ gland, further than all 2,000 expression-matched gene sets, and OXPHOS
> subunits lost more priority than nine in ten mitochondrial pathways, while the carnitine shuttle
> and fatty acid oxidation, which feed electrons to the respiratory chain, ranked among the largest
> gains (Fig. 1f, ED Fig. 2e,f). In the WT gland, OXPHOS subunits fell further than the
> proliferative signature (mean log2 fold change -0.25 versus -0.05), although the terminal end bud
> programme regressed and the mitochondrial share of the transcriptome did not fall (Fig. 1d,f; ED
> Fig. 2e); they also declined at the protein level (Fig. 1b). The adult gland thus withdraws from
> respiration far more than from proliferation.

**Asserted where computed:**
- **Fig. 1f:**
  - 2,000 draws (`NSET` in scripts 43 and 54) and OXPHOS at the 0th percentile in both windows;
  - OXPHOS losing more priority than at least 90% of the 143 pathways, in both windows;
  - the carnitine shuttle and fatty acid oxidation in the top fifth of gains, in both windows;
  - -0.25 and -0.05 to two decimals;
  - the paired OXPHOS-minus-proliferation test, and TEB below all matched sets.
- **Fig. 1d:** the WT nuclear MitoCarta share does not fall across the window (+17.6%, Welch on
  log2 share, P = 0.078).
- **ED 2g:** |r| < 0.1 between the window and MYC's effect, which "independently of MYC" rests on.

### 2b.3 The author's rulings

- **Fig. 1d:** no histogram; both parts in one PDF with one key; the text says 20-33%; 66 mm tall.
- **"Fell further than all 2,000 matched sets"**, with no multiplier.
- **The OXPHOS western is in Fig. 1b**, and the paragraph cites it for the protein-level decline.
- **Development contrasts:** "12W vs 6W, WT" in `#0072B2` and "12W vs 6W, MYC+" in `#D55E00`
  (`geno_cols` with `dev_labels`).
- **Size marks the leading contrast.** Fig. 1e draws 12W smaller than 6W; Fig. 1f and ED 2e draw
  MYC+ smaller than WT, the gland's own programme. ED 1a and ED 1c stay equal, because they compare
  6W and 12W directly.
- **Fig. 1f keeps filled dots** and draws no open ones: there is no test of record on priority
  (2b.4).

### 2b.4 Findings that changed the text

- **Priority across the window is not significant per pathway.** Script 08's Welch test passes BH
  for 0 of 144 pathways in WT (10 in MYC+):

| pathway | WT change | P | adjusted P | rank of 143 (WT, MYC+) |
|---|---|---|---|---|
| OXPHOS subunits | -0.109 | 0.12 | 0.46 | 14th, 13th lowest |
| carnitine shuttle | +0.289 | 0.011 | 0.16 | 2nd, 3rd highest |
| fatty acid oxidation | +0.144 | 0.030 | 0.25 | 8th, 25th highest |

  No priority null exists, so the paragraph states ranks.
- **Expression is the tested half.** OXPHOS subunits fall below all 2,000 matched sets in both
  windows (scripts 43 and 54), and below proliferation in the paired test (WT).
- **-0.26 was wrong:** the WT OXPHOS value is -0.2548, so -0.25.
- **Proliferation's own fall is beyond the matched 95%** (1.3-1.7th percentile). So the closing
  sentence says "far more than from proliferation", not "while sustaining proliferation".
- **The multiplier was dropped** because it depends on the weighting: about 5.4x unweighted and
  2.8x expression-weighted. The fall is not the five tissue-restricted COX paralogs (-0.226 without
  them), and 69 of the 87 subunits fall.
- **Retired from the old sentences** (2d.6, item 2):
  - "completely inverted": MYC's OXPHOS effect stays positive at both ages;
  - glycine cleavage: it gains no priority across the window;
  - "constant": now "did not fall".

### 2b.5 Page budget (183 x 170 mm)

- **ED Fig. 2:** d + g fill one row (182 x 62) and e + f the next (183 x 66), leaving 42 mm for
  a-c.
- **Fig. 1:** three rows:
  - d (100) beside e (83), 66 mm;
  - f (100) beside g/h (83), 56 mm;
  - a 56-mm top row (b's height).
  
  That is 178 mm, 8 mm over.

### 2b.6 Style changes and lessons, all in `FIGURE_RULES.md`

- **Axis breaks are rounded to 10 decimals** in `.nm_axis()`. An unrounded 0.4000000000000001
  fell outside an axis ending at 0.4, and ggplot dropped that tick (Fig. 1f). A self-test case
  covers it.
- **`stack_parts()` takes `ncol` and `widths`** for parts side by side (Fig. 1f, ED 2g).
- **Patchwork laid out the shared key with `theme_get()` at 11 pt.** The helper now gives it the
  panel theme and no outer margin, and `theme_nm()` keeps keys tight. All 11 panels were
  re-rendered after the change.
- **Base `pdf()` draws the ASCII hyphen as a minus.** So labels write "1C", not "1-C", while
  "OXPHOS - proliferation" prints a true minus (ED 2e).
- **Preview with Quick Look, not `pdftoppm`,** which draws Helvetica-Oblique upright.
- **Labels on scatters:** give the `ggrepel` layer every point, the unnamed ones with an empty
  label. ED 2d's eleven names did not fit at 62 mm, so its points are black and the legend names
  them.

---

## 2c. What happened on 2026-10-10

Four commits, all pushed: `b01da95`, `880e466`, and that day's handoff, `931e0a9` and `7505b39`.

### 2c.1 The author rearranged Fig. 1b and Extended Data Fig. 1

| slot | content | from the old layer |
|---|---|---|
| Fig. 1b | *Myc* mRNA per animal, beside the MYC western | S1E (the style demo) |
| ED 1a | the MYC/MAX/MXD network | S1F |
| ED 1b | MYC effect at 12W against 6W, MYC-responsive genes only; no MitoCarta facet | 1G, gene half |
| ED 1c | NES of every gene set at 6W and 12W, by programme | 1D |

- **This settled the three Fig. 1b questions** of 2d.5:
  - the blot lane order: the demo is used as is, age-major;
  - "target gene sets": ED 1c, the library-wide NES;
  - the programme panel: genes only, in ED 1b.
- **The author has updated the text:** the amplitude sentence now cites ED 1b,c, and the MYCER^T2
  model, previously cited as "ED Fig. 1b", has another letter.

### 2c.2 What was built

| slot | script | PDF (mm) |
|---|---|---|
| Fig. 1b | `fig1b_myc_transcript.R` | `Fig1/Fig1b_myc_transcript.pdf` (48 x 56) |
| ED 1a | `edfig1a_myc_network.R` | `EDFig1/EDFig1a_myc_network.pdf` (80 x 64) |
| ED 1b | `edfig1b_myc_effect_rescaled.R` | `EDFig1/EDFig1b_myc_effect_rescaled.pdf` (58 x 58) |
| ED 1c | `edfig1c_nes_by_programme.R` | `EDFig1/EDFig1c_nes_by_programme.pdf` (112 x 80) |

All four pass every check, and the author signed them off.
- **Fig. 1b:** DESeq2 IHW-adjusted P passed in. The script asserts that the drawn group means give
  DESeq2's log2 fold changes to within 0.1.
- **ED 1a:**
  - 6W and 12W are offset vertically, each with a 95% interval (log2FC +/- 1.96 SE). The old
    panel showed +/-1 SE at 6W only.
  - It shows the MYC+ vs WT gap at each age. Claude offered to draw the 6W-to-12W change within
    MYC+ instead, to match the text's word "fall"; the author did not take it up. The smallest
    interaction P is 0.83.
- **ED 1b:**
  - slope 0.487, asserted to be script 44's rate;
  - R2 0.513 (squared Pearson);
  - 22 genes outside +/-3 log2 are not drawn but enter every number, and the legend says so.
- **ED 1c:**
  - points coloured by age, in ED 1a's colours;
  - no median and no gridlines;
  - Spearman 0.933 across 866 sets between the ages.

### 2c.3 Decisions and findings, all recorded in `FIGURE_RULES.md`

- **Genotype-effect labels and colours:** "MYC+ vs WT, 6W" in `#D55E00` and "MYC+ vs WT, 12W" in
  `#E69F00`, as `contrast_cols` and `contrast_labels` in `_style.R`.
- **Intervals are 95% and unadjusted.** In ED 1a three intervals exclude 0 although no adjusted P
  passes, and the legend must say so:

| gene | age | raw P | IHW-adjusted P |
|---|---|---|---|
| *Mxi1* | 6W | 0.020 | 0.096 |
| *Mlx* | 6W | 0.039 | 0.16 |
| *Mlxip* | 12W | 0.001 | 0.061 |

- **No median or mean on a row that has sets on both sides of zero.**
  - In ED 1c the Apoptosis median goes +0.74 to -0.94, and its mean +0.17 to -0.35.
  - The means of three more rows change sign: TF target sets, MSigDB Hallmarks and Mammary
    development.
  - The shape barely moves. Apoptosis has 8 negative sets of 19 at 6W and 13 at 12W.
  - The author chose larger points and no median.
- **A correction to the old Fig. 1D legend.** Five Apoptosis sets change sign, and the old legend
  said none of them is significant at either age. In fact TANG_NECROPTOSIS reaches a
  within-category BH P of 0.033 at 12W (NES -1.53); the other four are not significant.
  - Suggested legend sentence: "At 12 weeks, 5 of the 19 apoptosis sets change sign (NES +0.6 to
    +0.9 at 6 weeks; -0.9 to -1.5 at 12 weeks). Four are not significant at either age;
    TANG_NECROPTOSIS reaches a within-category adjusted P of 0.033 at 12 weeks."
- **A display window that hides points** is allowed only if every number uses all points and the
  legend gives the count.

### 2c.4 Style changes, and a lesson

- **`_style.R`:** gained `contrast_cols` and `contrast_labels`, and a 2.5 mm right margin. Labels
  at the end of an axis overhang it by half their width, and "+3" and "+4" were being clipped. The
  self-test still passes 14 of 14.
- **The checks cannot see layout.** Clipped end labels, and ED 1c's "+4" running into the next
  facet's "-4" (now 6 mm apart), were caught only by rendering the PDF and looking at it. Render
  every panel before calling it done.

---

## 2d. What happened on 2026-10-09

Two commits, both pushed: `2cdd00e`, the figure layer (2d.4), and that day's handoff, `df9390d`.

### 2d.1 The paper is now a Nature Metabolism Letter

- **Source of record:** the Google Doc `11.8.26_MK_Myc_Paper` (Drive id
  `1M-kEC1d_VJmbTG-SyYYfykZ3F1_aEiDtu0xtwFaPtJo`), tab **"Nat metab version"**, edited 2026-10-09.
  - Title: *Mitochondrial respiration constrains MYC-driven mammary tumourigenesis*. Four figures.
  - The **"fig plans"** tab lists Figs. 2-4 panel by panel, mostly bench panels. **3k-m are marked
    "script 52, analysed, no panel". Its Fig. 1 table is empty.**
  - "CURRENT VERSION", the tab PANELS.md names, is the previous long form.
  - Reading it: the Doc is ~260k characters, so the Drive connector spills it to a JSON file; take
    `.fileContent` with jq. This session read the current tabs in full and only searched the three
    old drafts (OLD, OLD 2, Version 1 with figures).
- **Letter limits** (nature.com, 2026-10-09):
  - 2-4 display items; up to 10 Extended Data figures, one page each;
  - up to 40 references and 2,500 words; no headings;
  - a referenced introductory paragraph of up to 200 words instead of an abstract. The draft's
    abstract runs ~330 words; that is not a figure matter and was flagged once.
- **Figure specifications:**
  - 89 or 183 mm wide, at most 170 mm tall;
  - 5-7 pt Helvetica or Arial, panel letters 8 pt bold lowercase, no coloured text;
  - editable, embedded fonts in the final files.
- **The repo panels the text cites**, mapped from its sentences. The empty Fig. 1 table in "fig
  plans" can be filled from this once agreed:

| new slot | what the text cites it for | built in the OLD layer as |
|---|---|---|
| 1b, RNA-seq part | programme at half amplitude, same shape; Myc transcript steady (the western is bench) | `fig1_rescaled_not_reshaped` (old 1G), `figS1_myc_transcript_stable` (old S1E) |
| ED 1a | the MYC/MAX/MXD transcripts | `figS1_myc_network` (old S1F) |
| 1d | mitochondrial content +21-27% | `fig1_mito_content` (old 1E) |
| 1e | mitoPPS: OXPHOS at 6W, demoted on both timelines; catabolic arms gain | `fig1_reallocation_ranked` (old 1F, MYC effect only); the timelines are in `plane_arms_mitopps` or `fig2_wt_mito_contraction` (old 2F) |
| 1f | OXPHOS falls >5x proliferation against matched sets; Bax, Bbc3, Bcl2l1, Foxo3 | `plane_arms_content` and `plane_four_genes` (unslotted; built 2026-09-21 for this paragraph) |
| 1g | apoptosis ratios; PUMA:Bcl-xL | `fig2_priming_ratios` (old 2G, the ruling-4 rebuild) |
| 1h | PUMA drops only under MYC | probably `fig2_oxphos_puma_coupling` (old 2I); unconfirmed |
| ED 2a-e | design; the mitochondrial axis; total abundance; the adult programme | old S1A, old 1C, `fig1_compartment_levels`, old S2B |
| 3k-m | the escape series | script 53's panels A-C, an analysis figure; to be redrawn in the new layer |
| 4g-j | TCGA, SCAN-B, the chemotherapy cohorts, METABRIC | `../myc_human_exploratory`; no figure layer |

The old panels are **content references only**. Every Nat Metab panel is a new script in the new
layer.

### 2d.2 The author's rulings on the figure layer

All of these are recorded, with reasons, in `figures/natmetab/FIGURE_RULES.md`.
- **Panel by panel.** R makes single panels at their final printed size. The author assembles each
  figure in Illustrator with blots, images and Prism graphs, placing panels at 100%. Panel letters
  go on in Illustrator.
- **Start clean.** A new folder; no `panel_legend()`; no rebuild runner. `figures/panels/` stays
  untouched as the record, so section 6's retirement procedure does not apply. Stay on `main`. PDFs
  go in the repo under `outputs/natmetab/<Fig>/`.
- **Type:** Helvetica throughout (R, Prism, Illustrator), 6-7 pt. **Nothing under 5 pt, even for
  the initial submission.** Black text only.
- **Axes:** the y axis starts at 0 whenever an x axis is drawn, so magnitudes get no log axes. Both
  ends of an axis sit on a drawn tick.
- **P values:** exact, from the analysis of record, never computed by the plot, with an italic
  capital *P*. DESeq2 values are IHW-adjusted, and the legends say so.
- **Labels:** "WT / 6W" and "MYC+ / 12W" (MYC+, as in the text). Contrast labels are pending.
- **Goal:** at the end of the figure sessions, turn `FIGURE_RULES.md` and `_style.R` into a GS-lab
  figure skill.

### 2d.3 What was checked, and what was ruled out

- **The July font problem was spacing, not resolution** (commit a58fa1c). Re-checked on R 4.6.1:
  `cairo_pdf` embeds the font but collapses word spaces at 5-7 pt in Apple's viewers, while base
  `pdf()` is clean.
  - So panels use base `pdf()` with Helvetica, not embedded; Illustrator embeds the font on save.
  - Ghostscript, which is not installed, is not needed.
- **Dead end:** writing "Arial" from base `pdf()` by renaming Helvetica's metric files. Apple's
  viewers draw it in Times.
- **ggpubr is used only for `geom_bracket()`.** Its `stat_compare_means()` tests the plotted points:
  Wilcoxon by default, which cannot go below P = 0.0022 at 6 against 6, so it would contradict
  DESeq2. **ggstatsplot** is not used, for the same reason. **data-to-viz** is the checklist for
  chart forms.
- **P's case:** Nature's written guides are silent. One published Nat Metab paper prints an italic
  capital *P* 48 times and lowercase never: house style.
- **Skills:** 38 personal skills in `~/.claude/skills/` (installed 2026-07-22). They are optional,
  and their `cairo_pdf` advice is overridden here. No third-party Nature-figure skill improves on
  this layer: `nature-figure` exports with `cairo_pdf`.
- **Numbers read for 1b** (read-only):
  - *Myc*, MYC+ against WT: +1.68 at 6W and +1.82 at 12W (IHW P 2.7e-17 and 2.1e-20).
  - *Myc*, 6W to 12W: -0.10 in MYC+ (P 0.76), -0.25 in WT (P 0.48).
  - The programme across 2,648 genes: slope 0.49, R2 0.51.
  - Across the 143 MitoPathways: slope 0.552, R2 0.801. The R2 is at the 100th percentile of 500
    matched shuffles and the slope at the 87th. Tier slopes run 0.52-0.83.
  - Library-wide NES, 6W against 12W: 866 sets, Spearman 0.933.

### 2d.4 What was written: `2cdd00e`

`figures/natmetab/` holds three files.
- **`_style.R`:** the theme, palettes, P formatting, brackets, tick-ending axes, and `save_panel()`
  with its checks.
- **`_style_selftest.R`:** 14 checks. A correct panel passes, and nine broken panels are each
  refused.
- **`FIGURE_RULES.md`:** the rules, with reasons.

The self-test caught two bugs before the commit: an argument passed to the wrong function, and an
axis check that accepted an ending on an undrawn minor tick.

### 2d.5 Fig. 1b: three questions for the author (settled 2026-10-10, 2c.1)

1. **The lane order of the MYC western**, so that the mRNA panel can match it.
2. **"Without selective repression of specific target gene sets":** the LFC scatter covers
   MitoPathways only. Three options:
   - (a) the text says "mitochondrial pathways";
   - (b) add the library-wide NES ranking (866 sets, Spearman 0.93) as an Extended Data panel;
   - (c) draw a new scatter of every library set's mean MYC effect, 12W against 6W, with the MYC
     target sets marked. This needs a numbered script (55), run first.

   Claude leans to (c) if MYC targets in general are meant.
3. **The programme panel:** both halves (genes and pathways, about 89 x 45 mm) in 1b, or one half in
   1b and the other in Extended Data.

### 2d.6 Five sentences that do not match the data, for the walkthrough

1. **ED 2a,b** says the MECs were "purified from the tumours". They came from 6- and 12-week glands. *Resolved 2026-10-12 (2a.3).*
2. **1e** (the mitoPPS ruler, `background_vs_myc.rds$ruler`). *Resolved 2026-10-11 (2b.4).*
   - OXPHOS subunits' MYC effect is +0.112 at 6W (padj 0.27, not significant) and still +0.069 at
     12W, so it is not inverted. What falls is the change across the window: -0.109 in WT and
     -0.153 in MYC+.
   - Glycine cleavage does not gain priority across the window (+0.045 WT, -0.156 MYC+). It rises
     only on the content ruler, in WT (+0.240).
3. **1f**, "Foxo3 ... most distinct among all evaluated genes" (*2026-10-11: the new 1f does not draw
   the four genes, so this sentence has no panel yet*): the panel draws four genes. By the
   2026-09-21 ruling, Foxo3 is first of script 44's 26 by P and |z|, and fifth by size.
4. **1f,g**, "shifting distinctly the PUMA/BCL-XL ratio": the ruling-4 2G shows it as the largest of
   nine departures, inside the 95% prediction interval, and reports no verdict.
5. **1h**, if it is the coupling panel: the impasse ruling stands (permutation P 0.026 unadjusted,
   0.085 adjusted), so the panel asserts nothing.

*Items 3-5 are addressed by the 2026-10-12 proposals (2a.7):*
- *Foxo3* is drawn in 1h;
- "distinctly" is dropped from the draft paragraph;
- the coupling prints both fits and both P pairs.

---

## 3. What happened on 2026-09-22

Seven commits, all pushed: `dd58b8d`, `d3d7f3b`, `1a43f04`, `b881522` (this handoff's first
version), `8e4dd96`, `531d531` and `8114213`. This handoff's update is local.

There were also branch and tag operations, and these read-only analyses:
- the `analysis-exploratory` resolution (3.1);
- the FOXO3 timing report (3.3) and the dilution test that closed it (3.5);
- three TCGA checks on the human arm's data (3.6 to 3.8).

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
- **Ruled on the same day** (3.5): the activity step was tested again (the subset) and closed
  as a measured negative.

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

### 3.5 Dilution rejected, and the FOXO3 activity step closed -- `8e4dd96`

**The question** (the author's): does the programme score's undetected interaction hide a
responsive subset, targets falling under MYC within the window and averaged away?

**The rules, declared before retrieval:**
- at least 5 of the 38 genes with an interaction below -0.20;
- the subset mean below -0.20, **and** the remainder within 0.10 of zero;
- 1,000 seeded, baseMean-matched random sets.

**Both tests fail.**
- **k = 9 against a null median of 10** (32nd percentile). There is no hidden falling subset.
- **The remainder mean is +0.161** (94th percentile).

**The set skews the other way.**
- 13 genes sit above +0.20, against a null median of 5, and the set's mean interaction is at
  the 95th percentile.
- Across the 38, interaction and the 6W MYC effect correlate at **-0.89**. That is the 0.487
  rescaling.
- The set is rich in genes MYC represses at 6W (Esr1, Bcl6, Ar, Plau, Tnfsf10, Cdkn1a, Id1),
  whose repression relaxes by 12W. **The programme's interaction is the mirror of the global
  fade, not a FOXO3 signal.**

**Recorded as flaws or limits:**
- **75% of random matched sets pass the rule as written.** The null carried the test.
- **The six "drifting lower" genes** were defined after the numbers were seen. They are
  descriptive only and were not tested.

**The author's ruling.**
- The activity step has been tested three ways: the score, the canonical targets and the
  subset. None supports it.
- **Figure 1's Foxo3-to-PUMA link is transcript-level, and it stays there.** This is a measured
  negative, not an untested gap.

**Where it is recorded.**
- The verification note, **section 10**, set out like Check 1's failure. It has a pointer at the
  head of section 1, and 9.1's "not ruled on" sentence now points to it.
- **Esr1 is flagged for a future session:** the only raw interaction p below 0.05 in the set
  (+0.813), mostly the fade (7.7th percentile of script 44's scan). See note 10.3 and section 10.

### 3.6 TCGA: the ESR1-high and MYC-high high-mitochondria groups (read-only)

This is cross-sectional, association only. It is recorded in the follow-up note, section 3.1.
- **The data:** Menegollo's published fork calls on 849 TCGA tumours.
  - MB2 upper is "MYC-high" (M_a median +0.56, 38% ER+).
  - MB1 upper is "ESR1-high" (99% ER+).
- **The two are not mutually exclusive:** 21 tumours are in both, against 9.4 expected (odds
  ratio 3.1, 95% CI 1.7 to 5.5). The 21 are all ER+ and mostly LumB, with intermediate MYC
  scores.
- **ESR1 against the MYC score within the 158 high-mitochondria tumours:** rho -0.69, and -0.22
  with the three groups held fixed.

### 3.7 TCGA: do the two allocate the compartment differently? Parked -- `531d531`

**Circularity was checked before any scoring.**
- The fork-defining gene blocks were rebuilt from the upstream Menegollo scripts. MB2's is
  confirmed (rho 0.975), and MB1's is the best-supported candidate (0.88).
- The blocks are 7 to 13% MitoCarta. MB2's high side over-represents translation and
  respiratory genes, at 1.4 to 3.2 times the denominator rate.

**The primary ER+ contrast** (21 both-forks against 65 MB1-upper-only):
- import +0.027, translation +0.027, OXPHOS subunits -0.007;
- no interval excludes zero, and none lies below it.
- **The verdict: UNTESTABLE, not refuted.**

**Power:**
- The design detects 0.70 SD, and the observed effects are about 0.3 SD.
- **A follow-up needs about 116 both-forks tumours.**
- Script 07's null is shifted for this contrast (the median pathway sits at the 86th
  percentile).

**The secondary contrast is confounded by subtype.** 98 of 142 pathways differ, which is what a
basal-versus-luminal difference looks like.

**Parked** as `docs/2026-09-22_mitochondrial_priorities_followup.md`, with section 10 of this
handoff created to point at it.

### 3.8 Subtype as the subject, and what separating it costs -- `8114213`

**The author's framing** is in the follow-up note, section 6.
- For a future paper, subtype is the hypothesis: developmental state sets mitochondrial
  priorities, with the oncogene on top.
- **The prediction is declared in advance:** the subtype signal large, the MYC signal small.
  Today's result must not be cited as support.

**The cost check** (TCGA and SCAN-B, on M_a): the covariance breaks on the **luminal side only**.
- MYC-high luminal tumours are plentiful: 140 and 623.
- MYC-low basal tumours barely exist: 8 and 11.
- **SCAN-B is enough for the main effects.** The detectable MYC slope with subtype held fixed is
  0.060 SD per SD, but about 90% of that estimate comes from luminal tumours.
- **SCAN-B is not enough for the basal side.** A within-basal contrast detects only d 0.86, and
  d 0.5 would need a cohort of about 9,400.
- The scarcity may be biological.

**Everything through `8114213` was pushed at the author's word.**

### 3.9 Fig. 2G's decomposition, and a scope word for the two alts -- `ea0de33`

**The session opened on a false premise, and it is recorded so it cannot return.** The request
was to rebuild Fig. 2G on ruling 4, on the understanding that the committed panel was still the
16 August one. **It was not.** The rebuild landed on 2026-09-21 as `fea4d60` (3a.1), and the
committed file was checked against the whole spec before anything was written: it reads only
script 54's object, fits the line inside the ratio set, draws the band in full and all nine
ratios, carries no imported rate, and its block is clear of "priming" under an assert. The
author had carried the item over from the interrupted session without checking `fea4d60`.
- **The key's wording stands as the author set it:** "not used to fit the line", not "not in
  the fit", which had read as a verdict.

**What was actually outstanding, and is now done.**
- **2G gains one bound:** what the ratio's interaction is made of. On the members' own
  interactions it is about three quarters Bbc3 (-0.541) and one quarter Bcl2l1 (+0.179). The
  denominator enters with the opposite sign, so both members push the ratio down, and their
  difference (-0.720) is the members' version of the ratio's own -0.735.
  - Read from `$four_genes`, with asserts on the opposite signs (so the shares partition the
    difference) and on the two estimators agreeing in sign and within 0.05.
- **2G (alt) and 2H+I (alt) take a scope word.** "Entirely its numerator" now says it is about
  the TIMELINES, where Bcl-xL does not move (p 0.53 both), and each points at 2G's block for
  the genotype INTERACTION, which is a different contrast.
  - **The fraction is deliberately not in either panel** (the author's ruling). The two claims
    are not in conflict: one is about the temporal arms, the other about the interaction.
- All three panels render. 2H+I (alt)'s change is in its header comment only.

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

**Added 2026-09-22, afternoon: the dilution test and the human checks** (the verification
note's section 10; the follow-up note; 3.5 to 3.8).
- **The FOXO3 activity step is a measured negative, not an untested gap.** It was tested three
  ways, and none supports it. There is no hidden falling subset in the target set.
- **The programme's flat interaction is the mirror of the global fade,** not a FOXO3 signal.
- **A subset rule without the null built in discriminates nothing.** The dilution rule as
  written passes 75% of random matched sets.
- **The six "drifting lower" genes are descriptive.** Never cite them as a finding.
- **Human allocation (MYC-high against ER-high): untestable at n = 21, not refuted.** Never read
  those near-zero estimates as evidence of similarity.
- **The ESR1-high and MYC-high high-mitochondria groups overlap** (21 tumours). They are not
  two exclusive routes, and the whole comparison is association only.
- **In the human cohorts, a "MYC effect adjusted for subtype" is about 90% a luminal
  estimate.** MYC-low basal tumours barely exist (8 in TCGA, 11 in SCAN-B).

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

- ~~**The FOXO3 paragraph: transcript-level link, or a stated gap at activity?**~~ **Resolved
  2026-09-22:** the author ruled that the link is transcript-level and the activity step is a
  measured negative (3.5).
- **Parked, after the manuscript:** the human mitochondrial-priorities thread and Esr1
  (section 10).
- **Fig. S2D (alt)'s legend lacks the p19Arf/p16Ink4a isoform caveat** (section 2, item 3). It
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
- ~~**Fig. 2G (alt)'s "entirely its numerator"** is about three-quarters on these numbers
  (3a.3).~~ **Closed 2026-09-22** (3.9): the phrase is scoped to the timelines in 2G (alt) and
  2H+I (alt), where Bcl-xL does not move (p 0.53 both), and each now says the genotype
  interaction is a different contrast. The fraction stays out of both, by the author's ruling:
  the two claims describe different contrasts and both hold. The decomposition lives in 2G's
  own block.
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

**Everywhere (updated 2026-09-22, end of day):** `main` is the only branch, pushed through
`8114213`, and only this handoff's commit is local.

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

**For the Nat Metab panels, `figures/natmetab/FIGURE_RULES.md` replaces the panel bullets
below** (`theme_panel()`, `panel_legend()`, `_panel_common.R` belong to the old layer).

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
  n = 21, not refuted. A follow-up needs about 116 both-forks tumours. Its section 6 has the
  subtype framing and the cost of separating subtype from MYC.
- **Esr1**, from the FOXO3 dilution test (2026-09-22): the only gene in the FOXO3 target set
  with a raw interaction p below 0.05 (+0.813, raw p 0.036, IHW 1.00).
  - MYC represses it at 6W (-0.813) and not at 12W; most of that is the global fade (7.7th
    percentile of script 44's departure scan).
  - Flagged, not pursued: `docs/2026-09-21_two_timeline_verification.md` section 10.3.
  - In TCGA the ESR1-high and MYC-high high-mitochondria groups overlap (21 tumours), and ESR1
    is inversely associated with the MYC score within them (follow-up note, section 3.1).
