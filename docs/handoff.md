---
date: 2026-09-07
tags: [project/myc_mouse, handoff, block-c, paper-final]
status: live handoff -- overwrite in place at the end of each session
relates-to:
  - docs/experimental_cohorts_branch_notes.md          (the CURRENT branch's rules and status)
  - docs/2026-09-07_fatpad_tumour_limb_oxphos_trend.md (the fat-pad verdict, dataset now closed)
  - docs/2026-09-07_fatpad_timeline_oxphos_precheck.md
  - docs/2026-09-02_myc_oxphos_priming_gate_model.md   (the gate model)
  - docs/2026-08-27_human_validation_plan.md            (the parallel arm; three of its sections are superseded)
  - figures/panels/PANELS.md                            (the manifest of record)
  - docs/2026-08-02_results_narrative_final_order_D.md   (the narrative the cut is measured against)
---

# Handoff -- restart from here

Read this file first, then the docs named under "What happened this session" and "What was
already done". Everything below is current as of 2026-09-07.

**Two branches are live.** `paper-final` holds the manuscript and is where the text cut will
happen; `experimental-cohorts` is the current checkout and stages the two datasets outside
the MEC pipeline. Section 1 covers both.

---

## 1. State

**Current checkout: `experimental-cohorts`** (created off `paper-final`, merges back into it).

| | |
|---|---|
| HEAD | `3bf29cc` -- "A panel for collaborators: why the fat-pad series cannot carry a respiratory claim" |
| commits ahead of `paper-final` | **6** |
| pushed | **no.** `experimental-cohorts` is local only. |
| scripts | `49_fatpad_tumour_limb_oxphos_trend.R` -- author-run 2026-09-07, all 33 quoted numbers reconciled. Next number is **50**. |
| data | `data/fatpad_timeline/` pinned (CSV gitignored, README tracked) |
| panels | 28 manuscript panels UNCHANGED. `fatpad_confound.R` is a discussion panel outside both runners, like the six `biogax_*.R`. |
| scratch | `sandbox/` exists and is gitignored; it never merges |

**`paper-final`** (the manuscript branch, not currently checked out):

| | |
|---|---|
| HEAD | `b691133` -- "Script 48 re-run: 47 numbers checked back against the object, none moved" |
| ahead of `block-b-full` | 20, **all pushed**; `origin/paper-final` is 0 ahead / 0 behind |
| scripts | 00-48, all run and reconciled. Nothing outstanding on the mouse side. |

**Untracked and staying that way:** the six `.txt` session exports at the repo root,
`docs/library_reference/gray_chea_mito_tf_shortlist.csv`,
`docs/library_reference/2026-08-22_consensus_myc_double_hit_thread.md`,
`docs/2026-08-27_human_validation_plan.md` (the parallel arm's spec -- note it is therefore
**not on the remote**).

**THE TEXT CUT HAS STILL NOT HAPPENED.** It has been the pending main task for six
sessions. Block C exists to halve three written Results sections and reorganise the figures
to match; the cut text has not yet arrived. Nothing on `experimental-cohorts` blocks it or
depends on it.

---

## 2. Do this first

1. **Decide which branch.** If the cut text is ready, `git checkout paper-final` and go to
   section 6 -- that is the main line and it is blocked only on the text. If the next move is
   the **orthotopic Bcl-xL / PGC1a series**, stay on `experimental-cohorts` and read
   `docs/experimental_cohorts_branch_notes.md` first; it carries the standing rules and the
   two hard-won corollaries that dataset must obey before anything is fitted.
2. **Rebase before starting anything substantial** on this branch:
   `git rebase paper-final`. `paper-final` is still moving and the branch must not diverge.
3. **Nothing is outstanding on either dataset already done.** Scripts 48 and 49 are run,
   their objects are in step with their scripts, and every number in their docs has been
   checked back against the objects.

---

## 3. What happened this session -- the fat-pad timeline, closed

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

## 3a. The session before -- the gate model

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
8. **New, and read before touching `experimental-cohorts`:**
   `docs/experimental_cohorts_branch_notes.md` (rules, corollaries, dataset status), then
   `docs/2026-09-07_fatpad_tumour_limb_oxphos_trend.md` (sections 3, 7a, 7b and the
   Discussion paragraph) and `docs/2026-09-07_fatpad_timeline_oxphos_precheck.md`.

---

## 5. Corrections that must not be re-introduced

Eight statements have been made in this project, tested, and found false. They keep
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

**On `paper-final`:**

- **The text cut.** Blocked only on the text arriving. Everything else is done.
- Per-slot choice between each original panel and its alternative: **2F, 2G, 2H+2I, S2D**.
- `MYC` against `Myc` on the 2F/2G quadrant notes -- unresolved.
- A fifth arrow in `figS1_design_contrasts.R:60-93` if the diagonal panel is adopted.
- Whether the gate model earns a panel. It has none, deliberately.

**On `experimental-cohorts`:**

- **The orthotopic Bcl-xL / PGC1a series -- not started.** It is expected to produce
  main-panel work, and it is where the reversal question goes. Before anything is fitted it
  needs its own composition check, computed on the exact subset the test will use, and any
  adjustment must carry a positive control. Pin the data with a README first, following
  `data/fatpad_timeline/README.md`.
- The branch is **unpushed** and unmerged. Rebase on `paper-final` before the next
  substantial piece of work.
- `sandbox/precheck_fatpad_oxphos_confound.R` is scratch on disk, gitignored, and is deleted
  before merge. Only its dated note survives.
- The `Myc` 7.5x fat-pad discrepancy is open and needs per-epithelial-cell measurement to
  settle. Do not import either dataset's `Myc` scale into the other.

**Everywhere:** no push authorised beyond `paper-final`, which is current.

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

**Tell that session two things.** First, three of its plan's measurement choices are decided
by mouse data (section 3a) and two of the three were wrong. Second, the two traps from
section 3: a loading belongs to the SUBSET it was measured on, and any purity or composition
adjustment needs a positive control read before the endpoint of interest. A cohort stratified
by PAM50 / TP53 / purity is exactly where both bite, and that session does not read
`experimental-cohorts`.

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
  **Never `git add -A`** -- the six `.txt` exports and the untracked docs stay untracked.
  **Commit per signed-off change; do not push without a word.**
- **Check every claim against the data before drawing it,** and say when a sentence in the
  cut does not survive -- add it to the `PANELS.md` list rather than softening it silently.
