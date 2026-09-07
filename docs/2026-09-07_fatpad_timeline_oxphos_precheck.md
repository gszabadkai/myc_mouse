---
date: 2026-09-07
tags: [project/myc_mouse, experimental-cohorts, fat-pad-timeline, precheck, oxphos, confound]
status: settled -- FAIL on the pre-registered rule; the dataset's contribution is re-scoped, not cancelled
relates-to:
  - data/fatpad_timeline/README.md
  - docs/experimental_cohorts_branch_notes.md
  - scripts/48_gate_model_mouse_verification.R   (the ox_rel recipe this reuses)
---

# Can the fat-pad timeline carry a respiratory axis? No.

A pre-check, run before anything was built for this dataset. Scratch code under
`sandbox/`, nothing written to `results/`, no numbered script. The scratch does not merge;
this note is what survives.

---

## 1. The verdict

**FAIL**, on the rule fixed before the data were looked at:

> PASS iff `|rho(ox_rel, Adipoq)| < 0.3` in tumours (n = 15), with `Fabp4`, `Plin1` and
> `Cidec` agreeing in magnitude.

| ruler, tumours only (n = 15) | `Adipoq` | `Fabp4` | `Plin1` | `Cidec` |
|---|---|---|---|---|
| **`ox_rel`** | **+0.318** | **+0.486** | **+0.357** | **+0.561** |
| `ox_lvl` | +0.200 | +0.346 | +0.254 | +0.371 |
| `ox_mt` (13 mtDNA genes) | +0.675 | +0.821 | +0.679 | +0.732 |

Spearman, on `log2(DESeq-normalised counts + 1)`. `ox_rel` misses on the primary marker
(+0.318 against a 0.3 threshold) and misses clearly on the other three. Over all 30
samples it is worse still (`Adipoq` +0.640). At n = 15 no adjustment recovers a ruler
confounded this hard, so **the timeline dataset gets no respiratory-axis script.**

The same table read against the other composition axes says what the ruler is actually
measuring:

| ruler, tumours / all 30 | `Krt8` | `Epcam` | `Ptprc` | `Mki67` |
|---|---|---|---|---|
| `ox_rel` | +0.054 / −0.412 | −0.286 / −0.586 | −0.264 / +0.008 | **−0.461 / −0.722** |
| `ox_lvl` | +0.375 / −0.281 | +0.064 / −0.436 | −0.446 / −0.067 | −0.318 / −0.604 |

Respiratory score up, proliferation and epithelium down, adipose up. That is a tissue-
composition axis wearing a mitochondrial label.

## 2. The result that was not expected: the compartment-relative share made it worse

The pre-check was set up expecting `ox_rel` to clear the confound where `ox_lvl` did not —
that would have been direct evidence for the construction, which the MEC cohort cannot
give (it has no composition gradient to test against). **The opposite happened, in both
subsets.** The mechanism:

| subset | `rho(ox_sub, Adipoq)` | `rho(rest_mito, Adipoq)` | `rho(ox_rel, Adipoq)` |
|---|---|---|---|
| tumours | +0.200 | +0.157 | **+0.318** |
| all 30 | +0.538 | +0.326 | **+0.640** |

A within-compartment share cancels a confound **only when the confound loads equally on
numerator and denominator.** Here adipose loads *more* on the nuclear respiratory subunits
(+0.538) than on the rest of MitoCarta (+0.326), so the subtraction removes shared
non-adipose variance and leaves the adipose component proportionally **larger**. The share
amplifies this confound instead of cancelling it.

**This does not retract `ox_rel` for the MEC cohort or the human arm.** It is a statement
about a confound whose loading is unequal across the compartment, and it names the
condition under which the bridge ruler is safe — one that has to be checked per cohort
rather than assumed. Adipose is the worst case for it: fat is the one tissue whose
mitochondrial content is both large and respiratory-weighted.

## 3. The nu/mt separation, documented in this dataset

Carried over from the MEC work until now. Measured here:

| subset | n | `rho(ox_lvl, ox_mt)` | `rho(ox_rel, ox_mt)` | `rho(ox_mt, Adipoq)` |
|---|---|---|---|---|
| tumours | 15 | +0.464 | +0.632 | **+0.675** |
| all 30 | 30 | +0.763 | +0.762 | **+0.731** |

All 13 mtDNA-encoded genes are detected; they sit in the `ox_rel` **denominator** by
construction, and the nuclear numerator (`MITOCARTA_OXPHOS_SUBUNITS`, 87 of 89 symbols
mapped) contains none of them — both asserted in code, not assumed.

**The mtDNA proxy is the adipose reader.** The scoping figure of rho +0.729 against
`Adipoq` reproduces here as +0.731 over all 30 samples. Disregarding it was correct, and
now there is a number for why: it is the most adipose-loaded of the three rulers on every
marker. Note also that in *this* dataset the nuclear and mtDNA arms run **together**
(+0.464 to +0.763), not opposite as in the human arm's panel — because here both are
partly reading the same tissue fraction. That is a property of fat pad, not a
contradiction of the MEC/human finding.

## 4. The flag: `Myc` dose across 6W-12W

The MEC dataset's premise — **transgene dose is stable across 6W-12W** (mRNA + blot +
literature; a constant driver against a moving substrate) — is load-bearing for the
interaction argument. These data look different:

| gene | stat | 6WK_NEG | 6WK_POS | 12WK_POS | SMALL_TUM | LARGE_TUM |
|---|---|---|---|---|---|---|
| `Myc` | median | 759.0 | 988.5 | 7,406.2 | 9,062.7 | 10,980.9 |
| `Myc` | mean | 1,113.0 | 2,444.0 | 8,457.3 | 8,856.7 | 12,288.5 |
| `Krt8` | median | 7,494.7 | 8,908.6 | 11,495.7 | 13,087.8 | 14,822.1 |
| `Adipoq` | median | 59,342.3 | 59,189.7 | 31,595.1 | 29,515.0 | 26,414.8 |

Kruskal on `Myc` across the five groups, p = 0.0054. `6WK_POS -> 12WK_POS` is **7.5x** on
medians (989 -> 7,406) and **3.5x** on means (2,444 -> 8,457). `Krt8` rises only 1.3x over
the same interval, so epithelial expansion does not account for it.

**Reported, not resolved.** Three readings are live and this dataset cannot separate them:

1. **A real disagreement.** Transgene expression genuinely rises 6W->12W in this cohort,
   and the MEC premise does not generalise beyond purified MECs.
2. **Composition, in the direction that is not obvious.** `Adipoq` halves over the same
   interval (59,342 -> 31,595). A per-cell-constant transgene in an expanding epithelial
   compartment rises in whole-tissue units simply because the adipose denominator shrinks —
   but `Krt8`'s 1.3x sets a ceiling on how much of a 7.5x that can explain, so this cannot
   be the whole story either.
3. **Different animals, different assay.** Separate cohort, separate library prep,
   upstream normalisation this repo did not perform, and no shared samples with the MEC
   dataset. The two `Myc` numbers are not on a common scale and were never intended to be.

What would settle it: the MEC cohort's own `Myc` counts are on record and stable; the
decisive measurement is transgene expression **per epithelial cell** in this tissue —
sorted, or single-cell, or an IHC/blot on matched fat pad. Nothing in the pinned matrix
can do it.

**Consequence for the manuscript: none yet, and it must not be quietly imported.** The
stability premise is about purified MECs and is stated that way. A sentence claiming
stability *in tissue* would now need this dataset to agree, and it does not.

## 5. What the dataset is cleared to contribute

**Cleared: the two readings that need no respiratory score at all.** Both are single- or
two-gene readings on the linear counts, so the FAIL verdict does not touch them.

1. **The 6W genotype panel — an independent replication of the gate configuration in whole
   tissue.** `Bbc3` 90.0 -> 117.5, `Bcl2l1` 1,185.5 -> 1,001.2 (medians), with `Adipoq`
   matched at 59,342 / 59,190 — i.e. the adipose fraction does not move between the two
   groups, so this contrast is the one place in the dataset where composition is not a
   competing explanation. The PUMA:BCL-XL log2 ratio rises **+0.77** (median of the
   per-sample ratios; +0.62 as the mean), **Wilcoxon p = 0.016** at n = 5 vs 5.
   This is the mouse gate's sensor:guardian configuration reproducing in a second cohort,
   a different tissue preparation, and a different laboratory pipeline.

   **CAVEAT ADDED 2026-09-07, from script 49's per-sample record.** `Adipoq` is matched
   between the two 6W groups, and that was the argument that composition is not a competing
   explanation. **It is not sufficient: `Krt8` and `Epcam` are not matched.** Two of the five
   `6WK_POS` samples carry almost no epithelial signal -- `6WK_POS_R1` at `Krt8` 6.67 and
   `6WK_POS_R2` at 8.85, i.e. **102-fold and 23-fold below that group's own median** -- and
   `6WK_POS` is the only group in the dataset with a wide `Krt8` spread (6.85 log2, against
   0.46-1.41 in every other group). At 6 weeks the ductal tree occupies only part of the fat
   pad, so a sample can miss it; that is the likely origin. **Those two samples carry the
   group's two highest PUMA:BCL-XL values**, so they push the effect in its claimed
   direction. Recomputed without them:

   | subset | n | delta median | delta mean | Wilcoxon p |
   |---|---|---|---|---|
   | all 5 vs 5 (as claimed) | 5 v 5 | +0.772 | +0.624 | **0.016** |
   | `Krt8`-low dropped | 5 v 3 | **+0.646** | +0.482 | **0.071** |

   **The direction survives and the magnitude falls only ~16%, but the significance does
   not.** The honest statement is now: the configuration reproduces in the same direction in
   whole tissue, at n = 5 vs 3 once epithelium-poor samples are excluded, with p = 0.071 --
   supporting, not independently significant. Do not cite the p = 0.016 without the
   epithelial caveat. Recorded in code at `results/fatpad_tumour_limb_trend.rds$six_week_qc`.

2. **`Bcl2l1` is flat across the entire progression while `Myc` rises.** Medians 1,185.5 /
   1,001.2 / 948.4 / 986.1 / 1,010.9 against a `Myc` rise of ~11x from `6WK_NEG` to
   `LARGE_TUMOUR`; Kruskal p = 0.065 for `Bcl2l1`, p = 0.0054 for `Myc`. **This forecloses
   "BCL-XL accumulates as tumours grow"** — the guardian is demand-driven, set at the 6W
   genotype step and held, not progression-driven. A negative worth having: it is the
   obvious alternative reading of the human arm's guardian result, and it is now excluded
   in mouse tissue.

**Not cleared:** any respiratory score (`ox_rel`, `ox_lvl`, `ox_mt`, mitoPPS or GSVA on a
mito set), any statement about OXPHOS across this progression, and any use of these data
as the tissue-level counterpart to the MEC respiratory decline. If a respiratory question
must be asked of this material, it needs deconvolution against an adipose reference or a
sorted/single-cell replacement — not an adjustment at n = 15.

**Also not available, for a separate reason:** no `12WK_NEG` group exists, so the
genotype contrast this dataset can support is 6W only.

## 6. Method notes

- `ox_rel` / `ox_lvl` to script 48's recipe, read from the code and not from prose:
  `mito_all` = union of every `MITOCARTA_*` set in the library GMT (1,083 genes here);
  `ox_sub` = `MITOCARTA_OXPHOS_SUBUNITS` (87 mapped); `rest_mito` = the other 996, which
  is where the 13 mtDNA-encoded genes land. `comp_e(e)` = colMeans of row z-scores of
  `log2(counts + 1)`.
- Set membership through `functions/reconcile_gene_symbols.R` (`recon_to_ensembl`), as
  CLAUDE.md requires for any new set-based analysis.
- 15,963 of 54,838 genes at mean normalised count >= 10, zero-variance rows dropped.
- mitoPPS was **not** run: it needs the linear counts and a separate input object, and the
  decision rule does not depend on it.
- **Every number in the session brief was a group MEDIAN.** Means differ materially at
  n = 5-9 (`Myc` at `6WK_POS`: 989 median against 2,444 mean, a 2.5x gap). Both are given
  above wherever a claim depends on the choice. Name the statistic in any sentence built
  on these data.
