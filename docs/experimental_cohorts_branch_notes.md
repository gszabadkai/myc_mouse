---
date: 2026-09-07
tags: [project/myc_mouse, experimental-cohorts, branch-notes, standing-rules]
status: live -- the working notes for the `experimental-cohorts` branch; update in place
relates-to:
  - data/fatpad_timeline/README.md
  - docs/2026-09-07_fatpad_timeline_oxphos_precheck.md
  - docs/handoff.md
---

# `experimental-cohorts` — branch notes

Holds the two mouse datasets that are **not** part of the 6W/12W purified-MEC pipeline:
the orthotopic Bcl-xL / PGC1a tumour series, and the fat-pad progression timeline.

## Disposition

This branch **merges back into `paper-final`**. The orthotopic dataset produces main-panel
work, so this is a staging area, not a parallel line. Build accordingly, from the first
commit:

- **Numbered scripts continue `paper-final`'s sequence.** The next one is `49`. No separate
  numbering scheme to reconcile at merge.
- **Pinned data and provenance READMEs follow the existing conventions exactly**, since
  they arrive in `paper-final` unchanged: large matrices gitignored with a tracked
  `README.md` beside them (`data/FULL.DAT.csv` is the model), MD5 and source path recorded,
  traps named.
- **Scratch stays under `sandbox/`, which is gitignored**, and is deleted before merge.
  Scratch never merges whatever it concludes; only its dated note does.
- **Keep the branch rebased on `paper-final`** rather than letting it diverge. `paper-final`
  is still moving (the text cut has not happened), so rebase before each significant piece
  of work.

## THE STANDING RULE — scores are not comparable across cohorts, within a species

`ox_rel`, `ox_lvl`, GSVA scores, singscore, mitoPPS and every other z-score or
cohort-ranked construction are **cohort-relative**. Their zero and their unit are set by
the samples that went into the scoring run, and by nothing else.

**Therefore: the 6W/12W MEC cohort, the fat-pad timeline and the orthotopic series each get
their own scoring run, and their values are never pooled, never plotted on one axis, and
never differenced.** Not even when the recipe is byte-identical — an identical recipe over
different samples produces a different scale, which is exactly what makes the values look
comparable when they are not.

What **may** travel between cohorts:

- correlations (`rho(ox_rel, Adipoq)` in one cohort against the same in another),
- orderings and directions (which group is highest; which way a slope points),
- effect sizes expressed in each cohort's **own** SD, with the SD named,
- model *form* and its coefficients' signs (the gate's `bMX > 0`), never their magnitudes.

What may **not**: a score, a group mean of a score, a difference of two scores, a "shared"
axis, or any panel with two cohorts' scores on one scale.

Written down now, before the second dataset arrives, because after the merge all three
scoring scripts sit in one repo and the temptation to put their scores on one axis will be
immediate. `mitoPPS` carries the same prohibition for a second, independent reason already
on record: its pairwise-ratio baseline is composition-dependent.

**Corollary that has already bitten once.** A cohort-relative *construction* also carries
cohort-specific *behaviour*. The compartment-relative share `ox_rel` cancels a confound
only where that confound loads equally on numerator and denominator; in fat pad it loads
more on the respiratory arm than on the rest of MitoCarta, so the share **amplifies** it
(`docs/2026-09-07_fatpad_timeline_oxphos_precheck.md` section 2). Re-check the ruler in
each cohort; do not inherit its properties along with its recipe.

## Dataset status

| dataset | pinned | verdict | contributes |
|---|---|---|---|
| fat-pad timeline (Chandan, 30 samples) | `data/fatpad_timeline/` | **FAIL** for respiratory-axis work, 2026-09-07 | the 6W genotype panel (`PUMA:BCL-XL` +0.77 log2, Wilcoxon p = 0.016, `Adipoq` matched) and `Bcl2l1` flat across progression while `Myc` rises ~11x |
| orthotopic Bcl-xL / PGC1a series | not yet | not started | expected main-panel work |

Open, carried from the pre-check and **not** resolved: `Myc` rises 7.5x (medians) between
`6WK_POS` and `12WK_POS` in fat pad, against the MEC dataset's stable-dose premise. Three
readings remain live; see the pre-check note section 4. Do not import either dataset's
`Myc` scale into the other.
