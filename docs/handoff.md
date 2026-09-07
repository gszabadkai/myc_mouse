---
date: 2026-09-02
tags: [project/myc_mouse, handoff, block-c, paper-final]
status: live handoff -- overwrite in place at the end of each session
relates-to:
  - docs/2026-09-02_myc_oxphos_priming_gate_model.md   (this session's deliverable)
  - docs/2026-08-27_human_validation_plan.md            (the parallel arm; three of its sections are superseded)
  - figures/panels/PANELS.md                            (the manifest of record)
  - docs/2026-08-02_results_narrative_final_order_D.md   (the narrative the cut is measured against)
---

# Handoff -- restart from here

Read this file first, then the two docs named under "What happened this session" and
"What was already done". Everything below is current as of 2026-09-02.

---

## 1. State

| | |
|---|---|
| branch | `paper-final` (Block C, the halved paper), created off tag `block-b-full` |
| HEAD | `e348dd8` -- "The gate collapses for the trigger and not for the guardian -- I said it collapses" |
| commits ahead of `block-b-full` | **19** |
| pushed | **yes, 2026-09-03** -- `paper-final` is on `origin` at `e348dd8`, tracking set, 0 ahead / 0 behind. First push of this branch; no PR opened, since Block C is not finished. |
| panels | 28 `figures/panels/fig*.R` + 6 `biogax_*.R` discussion panels (deliberately outside both runners and the 28-count) |
| circulation PDFs | 12 / 8 / 4 pages |
| scripts | 00-48. **48 RUN by the author 2026-09-02**, positive control passing, reproducing the scratchpad run bit-for-bit. `results/gate_model_verification.rds` exists. It has since gained a two-endpoint `model_form` (below) and needs one ~2-minute re-run to carry it. |
| working tree | clean apart from the nine untracked-by-rule items |

**Untracked and staying that way:** the six `.txt` session exports at the repo root,
`docs/library_reference/gray_chea_mito_tf_shortlist.csv`,
`docs/library_reference/2026-08-22_consensus_myc_double_hit_thread.md`,
`docs/2026-08-27_human_validation_plan.md` (the parallel arm's spec).

**New this session, committed and pushed:**

- `scripts/48_gate_model_mouse_verification.R` -- author-run 2026-09-02; the committed
  version is **one re-source ahead** of the saved object (two-endpoint `model_form`).
- `docs/2026-09-02_myc_oxphos_priming_gate_model.md`
- `docs/handoff.md` (this file)

**Not on the remote, by the untracked rule:** `docs/2026-08-27_human_validation_plan.md`.
If the parallel human session expects to find its own spec on `origin`, it is not there.

**THE TEXT CUT HAS STILL NOT HAPPENED.** It has been the pending main task for five
sessions. Block C exists to halve three written Results sections and reorganise the
figures to match; the cut text has not yet arrived.

---

## 2. Do this first

1. **Re-source script 48 once** (~2 min). It has not changed in any way that could move a
   number -- `model_form` now loops over both co-primary endpoints instead of `PUMA` alone,
   and every `PUMA` row is identical to your run. The re-run only makes the saved object
   carry the `BUFFER` rows that section 3's crossover finding rests on. Verified in the
   scratchpad; the positive control still prints `+6.089 / 0.00523`.
2. **Then the cut**, if the text is ready. See section 6.

---

## 3. What happened this session -- the gate model

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

### The reconciliation (2026-09-02, second pass)

Script 48 was re-run in the scratchpad with the save path redirected, so `results/` is
untouched. The positive control passes and **every previously-reported number reproduces
bit-for-bit**. Five numbers in the doc were not carried by the saved object, which its own
provenance section requires; the script now stores all five (`separability$r_myc_dor`,
`axis_cor`, `set_sizes$n_genes_tested`, `simple_slopes$r2_wt`/`r2_myc`, `gene_gate$bX_myc`
-- additive fields only, no existing computation touched). Three confirmed as written
(16,819 genes tested; `mitorib_rel` <-> `ox_rel` = 0.755; the covariate-free wild-type
`ox_ppd` fit R2 0.325, p 0.053). **Two were wrong and are corrected in the doc:**

- the two MYC signatures do **not** both sit at 0.61 with `ox_rel` -- MSigDB 0.61,
  DoRothEA **0.64** (the transcript is 0.33). This strengthens section 4.1 rather than
  weakening it: the signature contains more of the axis, not less.
- section 7's `Mcl1` simple slopes were **-0.93 / -0.26**, from a probe rather than the
  script. The fitted values are **WT -1.00, Myc+ -0.20** (`gene_gate`, epi + imm adjusted).
  The claim -- both negative, so the ratio rises because BCL-XL falls faster -- is unchanged.

Script 44's percentiles cited in section 3.3 were checked against
`results/collapse_module_ownership.rds` and are exact (`Bbc3` 0.51, `Bcl2l11` 91.6).

**And the author's run then contradicted a third statement, the one that mattered most.**
Section 3.2 claimed the model collapses to a pure gate. That is true of `Bbc3:Bcl2l1` and
**false of `Mcl1:Bcl2l1`** -- the endpoint this same document promotes to co-primary. The
guardian balance has a real, covariate-surviving negative wild-type slope, so it crosses
over at `M*` = +0.545 instead of switching on at zero. `model_form` now fits both endpoints
so the object cannot hide it again, section 3.2 carries both tables, and the human
pre-specification (4.3) now says to retain the `X` main effect for `BUFFER`. The crossover
reaches p < 0.05 only on the genotype estimator -- transcript p 0.20, signature p 0.18, same
sign and consistent `M*` -- so read it as power, not as contradiction.

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
7. **New:** `docs/2026-09-02_myc_oxphos_priming_gate_model.md`.

---

## 5. Corrections that must not be re-introduced

Five statements have been made in this project, tested, and found false. They keep
coming back.

1. **PRC does not clear its expression-matched null** (7.1st pct). `padj 1.9e-4` at
   baseMean 2102 is precision, not effect size.
2. **Never quote the within-age OXPHOS <-> proliferation r = 0.89 raw.** OXPHOS <->
   `CORE_MITO` is 0.98 in the same samples; 0.89 *is* the ceiling. Use gene-level
   percentiles against all expressed genes.
3. **"Myc blocks the FOXO3 rise" is false.** `Foxo3` rises in both genotypes, interaction
   padj = 1, and `Bbc3` is flat in the wild-type timeline. Neither genotype half clears
   0.05, so neither may be labelled on a panel -- only their difference.
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

- Per-slot choice between each original panel and its alternative: **2F, 2G, 2H+2I, S2D**.
- `MYC` against `Myc` on the 2F/2G quadrant notes -- unresolved.
- A fifth arrow in `figS1_design_contrasts.R:60-93` if the diagonal panel is adopted.
- Script 48 run and reconciled. One cheap re-run outstanding, for the two-endpoint `model_form`.
- Whether the gate model earns a panel. It currently has none, deliberately -- it is a
  response-to-reviewers and human-arm asset, and the 28-panel count is unchanged.
- Nothing unpushed. `paper-final` is on `origin`; further pushes still need a word.
- The artifact watch on the gate-model page dropped (connection lost, 2026-09-07) and was
  not restarted. The page itself is unaffected and still current at the URL in section 3.

---

## 8. The parallel human session

The human arm runs in a separate Claude Code session and has so far written **only two
docs** into this tree and **no scripts** (`scripts/` stops at 48; nothing new in
`results/`). That is why running the two sessions in parallel has worked.

Before that session writes scripts or `.rds` files, give it its own worktree:

```
git worktree add ../myc_mouse_human -b human-validation paper-final
```

The pattern is already established by `../myc_mouse_main`. Two sessions sharing one
working tree collide in three places: the git index (`index.lock`, and one session's
commit sweeping up the other's half-finished edits), shared `results/` and `outputs/`
(a silent read of a half-written `.rds`), and same-file edits (last write wins, no merge).
Sequential is only necessary if the human arm needs to *change* mouse scripts or results.

**Tell that session about section 3.** Three of its plan's measurement choices are now
decided by mouse data, and two of the three were wrong.

---

## 9. Working rules, verbatim

- **Option A:** Claude Code writes `scripts/`, `figures/`, `paper/`; the author runs the
  numbered scripts in Positron. Running a panel to verify is fine and expected.
- **No explanatory prose on a panel.** `theme_panel()` blanks title/subtitle/caption;
  every script carries a `panel_legend()` block. Panels talk by NAMING DRAWN ELEMENTS.
- **Use the declared scheme;** anything new is declared in `_panel_common.R`, never inlined.
- **Verify with** `qlmanage -t -s 2000 -o <dir> <file>.pdf`. Never `sips`.
- **Run `./paper/clean_icon_files.sh repo` before the first git command** (Drive mount).
  **Never `git add -A`** -- the six `.txt` exports and the untracked docs stay untracked.
  **Commit per signed-off change; do not push without a word.**
- **Check every claim against the data before drawing it,** and say when a sentence in the
  cut does not survive -- add it to the `PANELS.md` list rather than softening it silently.
