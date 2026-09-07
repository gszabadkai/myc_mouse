---
date: 2026-09-07
tags: [project/myc_mouse, experimental-cohorts, fat-pad-timeline, oxphos, ordered-trend, confound]
status: settled -- T1 FAIL and T3 VOID; the reversal question stays OPEN and this cohort cannot close it
relates-to:
  - scripts/49_fatpad_tumour_limb_oxphos_trend.R  (writes results/fatpad_tumour_limb_trend.rds)
  - docs/2026-09-07_fatpad_timeline_oxphos_precheck.md  (the different, earlier estimand)
  - docs/experimental_cohorts_branch_notes.md
  - data/fatpad_timeline/README.md
---

# Does relative OXPHOS share rise again as tumours establish? Not measurably — and the test could not be made admissible

A between-group **ordered trend** across `12WK_POS` < `SMALL_TUMOUR` < `LARGE_TUMOUR`
(n = 20), pre-specified in full before the fit. This is a **different estimand** from the
pre-check's within-tumour correlation, which stays closed and is not reopened.

**Headline: the primary endpoint does not rise (T1 FAIL), and the hard admissibility
requirement fails too (T3 VOID) — so the null must not be read as evidence against the
reversal hypothesis.** The reason T3 fails is the substantive finding and is set out in
section 3.

---

## 1. The admissibility argument, and where it broke

Fixed in advance, in the session brief and in the script header:

> The pre-check measured adipose loading at **+0.538** on the nuclear OXPHOS subunits and
> **+0.326** on the rest of MitoCarta, so adipose pushes `ox_rel` **up**; adipose **falls**
> across the series. Composition alone therefore predicts the score **falling**. The
> prediction is that it **rises** after 12W. A positive result is credible because it runs
> against the confound; **a null is uninterpretable.**

The argument is sound in form. Two things about it did not survive contact with the data,
and the script was built to find out rather than assume (PART D):

1. **Those loadings were measured over all 30 samples. They do not hold within the limb.**
   Restricted to the 20 samples actually tested, `rho(ox_sub, Adipoq)` is **+0.144**, not
   +0.538.
2. **The argument was never established for the primary endpoint.** It was written for
   `ox_rel`. `ox_nuc_mtrib` is a different construction and, as section 3 shows, behaves
   the opposite way.

**The transferable lesson: a loading is a property of a subset, not of a ruler.** An
admissibility argument has to be computed on the exact subset the test runs on. Carried
into the branch notes.

## 2. T1 and T2 — the trend

Jonckheere–Terpstra, one-sided for an increase, 20,000 label permutations; Kendall's
tau-b with a 5,000-sample percentile bootstrap as the effect estimate.

| endpoint | tau [95% boot] | J-T one-sided p | Kruskal p |
|---|---|---|---|
| **`ox_nuc_mtrib`** (primary) | **−0.147** [−0.480, +0.241] | **0.800** | 0.716 |
| `ox_rel` (concordance) | −0.121 [−0.442, +0.304] | 0.758 | 0.082 |
| `ox_lvl` | −0.070 [−0.418, +0.337] | 0.663 | 0.145 |
| `mtrib_lvl` | +0.083 [−0.287, +0.445] | 0.336 | 0.374 |
| `mito_total` | −0.019 [−0.351, +0.354] | 0.554 | 0.282 |
| `ox_mt` | −0.313 [−0.579, 0.000] | 0.959 | 0.191 |
| `Mki67` | +0.313 [−0.058, +0.625] | **0.045** | 0.181 |
| `Adipoq` | −0.249 [−0.552, +0.104] | 0.919 | 0.352 |

- **T1: FAIL.** `ox_nuc_mtrib` does not rise. Its point estimate is *negative* and its
  interval spans zero; nothing here is distinguishable from flat.
- **T2: concordant.** `ox_rel` agrees in direction (also weakly negative), so the two
  independently-constructed rulers do not disagree.
- The only endpoint that moves at all on this limb is **`Mki67`** (tau +0.313, p = 0.045) —
  proliferation rises from 12W through small to large tumours, as it must. It is one
  nominal p among 21 endpoints and carries no weight on its own; it is reported because it
  shows the limb is not simply inert.

## 3. T3 — the hard requirement, and why it fails

Measured within the limb (n = 20), Spearman:

| marker | on `ox_sub` | on `mtrib` | on **`ox_nuc_mtrib`** | on `ox_rel` |
|---|---|---|---|---|
| `Adipoq` | +0.144 | −0.241 | **+0.576** | +0.220 |
| `Fabp4` | +0.182 | −0.227 | **+0.659** | +0.259 |
| `Plin1` | +0.080 | −0.328 | **+0.723** | +0.159 |
| `Cidec` | +0.113 | −0.307 | **+0.776** | +0.236 |

**The two halves of the primary endpoint load on adipose with opposite signs.** Nuclear
OXPHOS is weakly positive (+0.08 to +0.18); the mitoribosome is negative (−0.23 to −0.33).
Subtracting one from the other therefore **adds** their adipose components instead of
cancelling them, and `ox_nuc_mtrib` ends up the **most adipose-loaded construction on the
limb** — +0.58 to +0.78, well above `ox_rel` (+0.16 to +0.26) and above either of its own
parts. This is the pre-check's amplification mechanism again, in a sharper form: a
difference of two composites cancels a confound only if the confound loads on both with
the same sign and similar magnitude, and here it does neither.

The consequence for admissibility, computed as `tau(marker) x loading(marker on endpoint)`:

| marker | tau on the limb | loading on primary | predicted push |
|---|---|---|---|
| `Adipoq` | −0.249 | +0.576 | **−0.144** |
| `Fabp4` | −0.134 | +0.659 | −0.088 |
| `Plin1` | −0.083 | +0.723 | −0.060 |
| `Cidec` | −0.019 | +0.776 | −0.015 |

Adipose falls across the limb and loads positively on the primary, so **composition
predicts the primary falling — which is the direction observed** (tau −0.147, against a
predicted push of −0.144 from `Adipoq` alone). The confound runs **with** the observed
trend, not against it. By the rule fixed in advance, **T1 is void whatever its p-value.**

**Not an outlier artefact.** `LARGE_TUMOUR_R8` sits at `Adipoq` 10.7 (log2) where every
other limb sample is 13.9–15.6, and it carries the lowest primary score (−1.94). Dropping
each sample in turn:

| quantity | full | leave-one-out range |
|---|---|---|
| tau, primary | −0.147 | **−0.244 to −0.057** (never positive) |
| `rho(primary, Adipoq)` | +0.576 | **+0.505 to +0.688** |

The confound is structural and the absence of a rise is not one animal.

## 4. T4 — nothing moves on this limb

| endpoint | tau | one-sided p |
|---|---|---|
| `mito_total` (whole nuclear MitoCarta) | −0.019 | 0.554 |
| `ox_lvl` | −0.070 | 0.663 |
| `mtrib_lvl` | +0.083 | 0.336 |

Neither the compartment as a whole nor its respiratory arm nor its translation arm changes
detectably from 12W through small to large tumours. So the T4 question — reprioritisation
(share moves, total flat) against a content claim (both move) — **does not arise**: nothing
moves. The per-sample figure shows the same thing directly: every panel falls from 6W to
12W and is then flat across the three limb groups.

## 5. T5 — the guardian, observed

`Bcl2l1` tau **+0.147** (p 0.223) against the primary's −0.147; `Bbc3` +0.160 (0.204);
`log2(Bbc3/Bcl2l1)` +0.147 (0.224); `Mcl1` +0.045; `Myc` +0.147 (0.219). All inside the
noise at n = 20, and no threshold was set for this. **The hoped-for reading — respiratory
share rising into tumours with the guardian tracking it, the human configuration appearing
as a mouse trajectory — is not available**, because the antecedent does not hold: the share
does not rise. `Bcl2l1`'s flatness here is consistent with the pre-check's finding that it
is flat across the whole progression while `Myc` rises ~11x, which remains this dataset's
solid negative.

## 6. Group summaries — both statistics

Medians above, means below; the pre-check established they diverge materially at these
group sizes, so any sentence must name which it rests on.

| endpoint | stat | `12WK_POS` | `SMALL` | `LARGE` |
|---|---|---|---|---|
| `ox_nuc_mtrib` | median | −0.110 | −0.186 | −0.221 |
| | mean | −0.048 | −0.051 | −0.375 |
| `ox_rel` | median | −0.005 | −0.496 | −0.320 |
| | mean | −0.024 | −0.332 | −0.331 |
| `ox_lvl` | median | −0.098 | −0.757 | −0.454 |
| | mean | −0.144 | −0.499 | −0.468 |
| `mtrib_lvl` | median | +0.005 | −0.778 | −0.339 |
| | mean | −0.096 | −0.449 | −0.093 |
| `mito_total` | median | −0.107 | −0.320 | −0.154 |
| | mean | −0.124 | −0.191 | −0.157 |
| `ox_mt` | median | −0.146 | −0.200 | −0.676 |
| | mean | +0.149 | −0.187 | −0.674 |
| `Adipoq` (log2) | median | 14.95 | 14.84 | 14.69 |
| `Mki67` (log2) | median | 10.01 | 10.70 | 10.72 |
| `Bcl2l1` (log2) | median | 9.891 | 9.947 | 9.983 |

Per-sample values are in `results/fatpad_tumour_limb_trend.rds$per_sample`; the figure is
`outputs/fatpad_limb/49_per_sample_trends.pdf` (eight endpoints, all 30 samples, limb
samples filled and 6W samples open so the excluded groups are visible but marked).

## 7. What the 12W-to-tumour limb now supports

**Nothing about the reversal, in either direction.** The pre-registered logic was explicit
that a positive result would be credible because it ran against the confound and that a
null would be uninterpretable; the result is a null, and T3 shows the confound in fact runs
*with* the weak downward drift that was observed. So this limb neither supports nor
refutes the claim that relative OXPHOS share rises again as tumours establish. **The
question raised by the human arm — MYC-high tumours are respiration-high, so the
developmental decline must reverse somewhere — remains open, and this cohort cannot close
it.**

What would: material without a majority-adipose background, or a deconvolution against an
adipose reference — sorted tumour epithelium, single-cell or spatial data from the same
model, or the orthotopic series, where the tumour is not embedded in fat pad. The
orthotopic dataset is the next thing this branch takes on and is the natural place for the
question to go, subject to its own composition check first.

**What is unchanged:** the MEC cohort's 6W-to-12W reprioritisation, which is clean and does
not depend on anything here; and this dataset's two score-free contributions (the 6W
genotype panel, `PUMA:BCL-XL` +0.77 log2 at Wilcoxon p = 0.016 with `Adipoq` matched; and
`Bcl2l1` flat across the progression while `Myc` rises ~11x).

## 8. Method notes

- **Sets are script 08's authoritative splits**, rebuilt here from `Mouse.MitoCarta3.0.xls`
  Sheet 4 with `splitstackshape::cSplit` and every `mt-*` gene stripped from every pathway,
  exactly as 08 does. They are *not* read back from `mitopps_scores.rds`, whose gene lists
  were already filtered to the MEC matrix — a cohort-dependent filter that must not travel.
  Mapped in this cohort: OXPHOS subunits **87**, mitochondrial ribosome **83**, nuclear
  MitoCarta **975**, mtDNA-encoded **13**. `ox_rel` keeps script 48's GMT recipe (87 and
  996) so T2 compares two independently-constructed rulers.
- Membership through `functions/reconcile_gene_symbols.R`, per CLAUDE.md.
- 15,963 of 54,838 genes at mean normalised count >= 10, zero-variance rows dropped.
  z-composites on `log2(counts + 1)`. **mitoPPS was not run** — it needs the linear counts
  and its own input object, and no decision rule depends on it.
- **Jonckheere–Terpstra is implemented in the script**, not taken from a package, so it adds
  no dependency `00_setup_packages.R` does not already load. Verified against
  `PMCMRplus::jonckheereTest` on five endpoints: the standardised statistic agrees to four
  decimal places (`ox_nuc_mtrib` z −0.8066, `ox_rel` −0.6664, `Mki67` +1.7185, `Adipoq`
  −1.3678, `ox_mt` −1.7185); the permutation p and the package's normal approximation differ
  only in the third decimal. `PMCMRplus` is not a dependency of the script.
- No FDR across endpoints, per repo standard: estimates with intervals are reported and the
  pattern carries it.
- **Cross-cohort comparison is forbidden** and nothing here does it. These scores are
  cohort-relative; only the directions and orderings above may travel.
