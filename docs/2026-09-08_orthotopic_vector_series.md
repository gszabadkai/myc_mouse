---
date: 2026-09-08
tags: [project/myc_mouse, experimental-cohorts, orthotopic, bcl2l1, pgc1a, myc-ceiling, reversal]
status: settled -- C1 holds and survives adjustment; C2 is an observation and holds; C3 is UNINFORMATIVE and the abstract sentence is left untested
relates-to:
  - scripts/50_orthotopic_vector_series_scoring.R  (writes results/orthotopic_vector_series.rds)
  - docs/2026-09-07_orthotopic_analysis_plan.md    (the spec; Gate 1 answered below)
  - data/orthotopic_series/README.md
  - docs/2026-09-02_myc_oxphos_priming_gate_model.md  (the gland reference direction)
  - docs/experimental_cohorts_branch_notes.md
---

# The orthotopic vector series: one clean causal result, one observation, and one question this dataset cannot answer

> **QUALIFIED AND CORRECTED 2026-09-08 by
> `docs/2026-09-08_orthotopic_specificity_checks.md` (script 51).** No number in
> this note changes and no verdict is reopened -- script 51 rebuilt this cohort
> from source and reproduces the section 3 C1 table to 0.0238%. Four things it
> adds, and the third and fourth must travel with any quotation from here:
>
> 1. **C1 is MORE SPECIFIC than this note could show.** `Mcl1` does not move
>    (-0.020 [-0.154, +0.127]) and `Bbc3` FALLS (-0.727 [-1.299, -0.376]) while
>    `Bcl2l1` rises. Among the anti-apoptotic members only `Bcl2l1` rises;
>    `Bcl2l2` falls. C1 is a selective event, not a general anti-apoptotic
>    response. (`Bik` rises, so the gap sentence is scoped to `Bbc3`.)
> 2. **PGC1a is a BROAD BIOGENESIS intervention with a 2:1 OXPHOS tilt** --
>    numerator +1.925, denominator +1.011 -- not an OXPHOS-selective one.
> 3. **C2's section 4 `p = 0.0029` IS NORMALISATION-DEPENDENT.** On
>    median-of-ratios the same contrast is **p = 0.0385**. The mechanism is
>    measured: the Pgc1a arm carries **14.7% of its library in MitoCarta genes
>    against 8.0% in EV**, which deflates every other gene's CPM in that arm.
>    **Quote 0.0385, or quote both.**
> 4. **C3's CONSTRUCTION was inconsistent with C1's.** C1's spec said BclxL is
>    reported and never pooled; section 5 pooled it, so 7 of the 21 `Bcl2l1`
>    values are set by an expression construct rather than by the tumour.
>    **C3 stays UNINFORMATIVE** -- script 51 gives two independent reasons,
>    neither of which is power -- but section 5's readings should be read with
>    that correction in view.

Author-run 2026-09-08; every number below is read from
`results/orthotopic_vector_series.rds`. Scoring set **fixed before any number was computed**
at the 21 vector-series samples (EV / BclxL / Pgc1a, 7 each). The p21 series is out of scope;
MYAZ / FaMY / BoMY are excluded as a separate series.

| | verdict |
|---|---|
| **GATE** | `adipose` and `endothelial` reach the material threshold against `ox_rel` — but see §2, the pooled loading is substantially the construct effect |
| **C1** | **HOLDS.** `Bcl2l1` +121.8 CPM [98.1, 151.2], one-sided p = 0.00029, and survives composition adjustment with the positive control intact |
| **C2** | **HOLDS as an observation.** 80th-percentile contrast −903 CPM (p = 0.0029); zero PGC1α tumours above 2000 `Myc` CPM against 3 EV and 4 BclxL |
| **C3** | **UNINFORMATIVE.** rho +0.064 [−0.347, +0.493], and the sign flips across the three readings. The abstract sentence is **neither supported nor withdrawn** |

---

## 1. Gate 1, answered

- **Quantification level.** No transcript-level quantification exists anywhere on disk — no
  `quant.sf`, no salmon output directories, no `tx2gene`, no transcript matrix, under
  `/Users/gs/G/data/MK_myc_2022` or in this repo. **`tximport` is unavailable.** Counts go
  through `DESeqDataSetFromMatrix` on rounded salmon estimated counts, so the
  **transcript-length offset is lost**. That applies to every number in this note. It is a
  real limitation and not a formality: gene-level counts without the offset misstate genes
  whose isoform usage differs between arms, and nothing here can detect that.
- **MYAZ / FaMY1 / BoMY.** `MYAZ` is identifiable — the MMTV-Myc tumour-derived implant line
  that is the parental background for the vector series. **`FaMY1` and `BoMY` are not
  determinable from files on disk**: no sample sheet, no metadata, no coldata exists for this
  cohort, and the only text mention anywhere is the analysis plan. Three file-level facts put
  all three in a separate series: a different animal-ID namespace (bare `410-6h`,
  `BFGE516-1h` against `BRNS###-#x` for all 49 vector and p21 samples), n = 6 against 7, and
  the originating analysis pairing them only with each other — **there is no MYAZ-vs-EV
  contrast on disk**. Excluded, not pooled.
- **Cohort.** Vector series alone, n = 21. Wider set as a named sensitivity in §2 only.
- **Composition.** §2.
- **Load control.** The median-CPM table for `Bcl2l1` / `Ppargc1a` / `Myc` across all ten
  groups reproduces the plan's §2 to **0.78%**, so the load and the group parsing are the
  ones the preliminary pass used. Group tokens are parsed longest-first and asserted; `EV`
  and `Pgc` matched before `NTEV`/`KOEV`/`NTPgc`/`KOPgc` would collapse four groups into two.
  `ox_sub` maps **87** nuclear subunits here, not the plan's 88 — reported, not forced.

## 2. The gate, and a nuance that changes how it reads

Spearman inside the 21, `ox_rel`:

| compartment | pooled (n = 21) | EV | BclxL | Pgc1a | all 67 (sensitivity) |
|---|---|---|---|---|---|
| **endothelial** | **+0.692** | 0.000 | −0.357 | +0.786 | +0.425 |
| **adipose** | **+0.621** | +0.036 | −0.393 | +0.714 | **+0.093** |
| immune | +0.457 | +0.679 | +0.107 | +0.071 | +0.081 |
| proliferation | +0.126 | −0.250 | +0.393 | +0.179 | −0.215 |
| epithelial | −0.068 | +0.536 | −0.643 | −0.321 | +0.080 |
| stromal | −0.155 | +0.536 | +0.214 | +0.643 | −0.158 |

Two compartments clear the pre-registered |rho| ≥ 0.5, so PART F was load-bearing. **But the
pooled loading is substantially the construct effect, not a cellularity gradient**, and the
table says so itself:

- the **per-group** correlations do not agree — `adipose` runs +0.036 / −0.393 / +0.714 across
  EV / BclxL / Pgc1a, and `endothelial` 0.000 / −0.357 / +0.786. A genuine cellularity
  confound would load consistently within arms;
- **`adipose` collapses to +0.093 over all 67**, where the construct contrast is diluted by
  seven other groups;
- the constructs move `ox_rel` enormously by design (Pgc1a vs EV **+0.915**, §3), so anything
  that also differs by arm correlates with `ox_rel` across the pooled 21 by construction.

Per-group rhos at n = 7 are individually near-uninterpretable (|rho| = 0.71 is p ≈ 0.07), so
this is a statement about *pattern*, not about any one cell. The consequence is that PART F's
adjustment is partly an adjustment **for the intervention itself** — which is exactly the
over-adjustment risk the positive control exists to detect. It did its job: see §3.

## 3. C1 — the licensing relationship, by intervention

`Bcl2l1` in Pgc1a against EV, one-sided by pre-specification, n = 7 v 7:

| measure | median Pgc1a | median EV | difference [95% boot] | one-sided p |
|---|---|---|---|---|
| **`Bcl2l1` (CPM)** | 192.8 | 71.0 | **+121.8 [98.1, 151.2]** | **0.00029** |
| `Ppargc1a` (CPM) — manipulation check | 301.6 | 0.037 | +301.6 [268.6, 356.4] | 0.00029 |
| `ox_lvl` (OXPHOS level) | +1.216 | −0.708 | +1.925 [1.663, 2.296] | 0.00029 |
| `ox_rel` (respiratory share) | +0.607 | −0.308 | +0.915 [0.761, 1.007] | 0.00029 |
| `Bcl2l1` — BclxL arm (construct present) | 326.4 | 71.0 | +255.4 [181.9, 292.7] | 0.00029 |

p = 0.00029 is **1/3432, the floor for 7 v 7 one-sided** — the two groups separate perfectly
on all four measures. The manipulation worked at the level it claims to: `Ppargc1a` goes from
undetectable to 302 CPM, and both the absolute OXPHOS level and the compartment-relative
share move with it.

**And it survives the gate.** Adjusted for `adipose` + `endothelial`:

| measure | β Pgc1a vs EV | p |
|---|---|---|
| `Bcl2l1` (CPM) | **+120.0** | **0.0031** |
| **`Ppargc1a` (CPM) — POSITIVE CONTROL** | **+275.2** | **2.9e-9** |
| `ox_lvl` | +2.07 | 2.5e-7 |

The control holds emphatically, so the adjustment is not eating the group effect wholesale
(the fat-pad failure mode), and `Bcl2l1` retains 98% of its unadjusted effect. **This is the
one clean causal result in the dataset.**

**What it licenses.** PGC1α overexpression alone, with no Bcl-xL construct present, raises
`Bcl2l1` in vivo by 2.7×. That is a positive guardian–respiration association **produced by
manipulating respiration**, not observed by correlating it — the argument that the human
`BCL2L1` coefficient reflects a causal link rather than a cross-sectional association. Write
it as an intervention result at n = 7 v 7 in one model, not as an effect size to transfer.

## 4. C2 — the MYC ceiling, as an observation

Exact null: all `choose(14,7)` = **3432** relabellings of the two arms contrasted.

| statistic | contrast | observed | one-sided p |
|---|---|---|---|
| max | Pgc1a vs EV | −1207 CPM | 0.096 |
| **80th percentile** | **Pgc1a vs EV** | **−903 CPM** | **0.0029** |
| max | Pgc1a vs BclxL | −613 CPM | 0.035 |
| 80th percentile | Pgc1a vs BclxL | −802 CPM | 0.012 |

| group | `Myc` max | n over 2000 CPM | `Bcl2l1` min–max | fold spread |
|---|---|---|---|---|
| EV | 3170 | 3 | 61.8–97.7 | 1.58 |
| BclxL | 2576 | 4 | 200.0–365.5 | 1.83 |
| **Pgc1a** | **1963** | **0** | 181.2–227.3 | **1.25** |

**Reporting both statistics was necessary, and it settles the obvious objection in the
opposite direction to the one expected.** The *max* contrast against EV is the **weakest**
cell in the table (p = 0.096) because EV carries one animal at 3170 CPM that makes its own
maximum unstable; the 80th percentile, which is not hostage to that animal, is the strongest
(p = 0.0029). So the ceiling is **not** an artefact of a single PGC1α tumour — and the
`Bcl2l1` convergence (1.25× against 1.58× and 1.83×) is a floor as well as a ceiling.

**This is a figure and an observation, not a test.** n = 7 per arm, four contrasts, no
multiplicity correction, and the groups differ in the intervention that defines them. It is
consistent with M and O not both being high in an established tumour; it does not demonstrate
it.

## 5. C3 — which side of the reversal: this dataset does not say

`Mcl1:Bcl2l1` and both members against `ox_rel`:

| endpoint | reading | rho | 95% boot | leave-one-out | sign stable? |
|---|---|---|---|---|---|
| **`guardian_ratio`** | pooled (n = 21) | **+0.064** | **[−0.347, +0.493]** | −0.054 to +0.144 | **no** |
| `guardian_ratio` | group-adjusted | +0.300 | [−0.255, +0.737] | +0.223 to +0.505 | yes |
| `Mcl1` | pooled | +0.279 | [−0.155, +0.640] | +0.187 to +0.362 | yes |
| `Mcl1` | group-adjusted | +0.138 | [−0.353, +0.581] | +0.041 to +0.317 | yes |
| `Bcl2l1` | pooled | −0.065 | [−0.524, +0.387] | −0.161 to +0.053 | no |
| `Bcl2l1` | group-adjusted | −0.194 | [−0.616, +0.295] | −0.335 to −0.113 | yes |

And composition-adjusted (PART F): `guardian_ratio` β = **−1.44 (p 0.111)**, `Mcl1` −0.214
(0.135), `Bcl2l1` +1.228 (0.145).

**Every interval spans zero, and the sign of the primary endpoint flips across the three
readings** — pooled **+**, group-adjusted **+**, composition-adjusted **−**. The
leave-one-out range on the pooled reading also crosses zero.

**So this dataset places these tumours on neither side of the reversal.** Reference
directions, quoted for direction only and never pooled with these cohort-relative scores: the
gland gives **+0.351** on the same endpoint and ruler (within-timepoint,
`docs/2026-09-02` §3.4), human tumours **−0.31 to −0.46**. The orthotopic estimate is
compatible with both and with zero.

**Consequence for the manuscript, and it is the one that needs a decision.** The abstract
sentence

> *"Mouse tumours resemble human tumours rather than the gland they arose from, placing this
> reversal at the transition into malignancy rather than between species."*

is **not tested by this dataset**. It is not refuted either. If it is currently written as
supported, it needs a different source or it comes out — **it must not be cited to the
orthotopic series.** The honest form, if the sentence is kept at all, is that mouse tumours
were examined and did not resolve the question at n = 21.

Why it fails is worth stating, because it is not simple lack of power. The vector series was
built to *move* `ox_rel` between arms (+0.915, §3), which makes it excellent for C1 — a group
contrast — and poor for C3, which needs `ox_rel` to vary *within* a common state. The pooled
reading mixes the construct contrast into the association; the group-adjusted reading removes
it and is left with n = 7 per arm. **A design optimised for the intervention question is the
wrong design for the correlation question**, and no analysis choice fixes that here.

## 6. What this licenses, and what it does not

**Licenses:**
- the causal reading of the human `BCL2L1` result — PGC1α alone raises the guardian in vivo
  (§3), adjustment-robust with a positive control;
- the in vivo selection observation — no PGC1α tumour reaches the `Myc` levels EV and BclxL
  tumours reach, on the upper tail rather than the mean (§4), as an observation.

**Does not license:**
- **anything about the MYC × OXPHOS interaction.** Every arm is MYAZ-derived and MYC-high;
  there is no MYC-low arm, so it is not identifiable here at any n. It lives in the iMMEC
  rtTA-MYC ±dox × ±PGC1α design on the death readout. No group term was allowed to stand in
  for it;
- **the abstract's "mouse tumours resemble human tumours" sentence** (§5);
- **anything about apoptotic competence as a cellular state.** N3: these are transcript
  associations, and the word is not applied to a transcript here;
- **any pooling** with the 6W/12W timeline or the human cohorts. Species = cohort; these
  scores are cohort-relative and only directions travel.

## 7. Method notes

- Scoring set fixed at n = 21 before any number was computed; `ox_rel` is relative, so the
  cohort choice sets every value. The all-67 column in §2 is a sensitivity computed in **its
  own scoring run** and is not numerically comparable with the rest of the table.
- `ox_rel` = nuclear MitoCarta OXPHOS subunits minus the rest of nuclear MitoCarta, built from
  `Mouse.MitoCarta3.0.xls` **Sheet 4 by exact filename** (never a glob) with every `mt-*` gene
  stripped from every pathway as script 08 does. **87** subunits map here (the plan states 88),
  880 in the denominator, all **13** mtDNA-encoded genes present and in the denominator —
  asserted in code, not assumed. Membership through `functions/reconcile_gene_symbols.R`.
- 15,721 of 78,334 genes at mean normalised count ≥ 10. Two input objects that never mix:
  linear normalised counts and `log2(+1)`. **mitoPPS and GSVA were not run** — no
  pre-specified claim needs either, and a cohort-relative score that nothing rests on is how
  such a number ends up quoted by accident.
- Marker panels require **≥ 3 detected genes**; a two-gene "composite" is a gene. All six
  compartments cleared it after `endothelial` was widened to `Pecam1`/`Cdh5`/`Cldn5`/`Emcn`/`Tek`.
- No per-gene FDR; estimates with bootstrap intervals throughout (5,000 resamples), and the
  C2 null is exact rather than sampled.
- **A correction made during the build:** the first draft's C3 verdict called rho +0.064 with
  an interval from −0.347 to +0.493 "gland-like, sentence must be withdrawn" on the **sign
  alone**. The reading rule is now three-way and fixed in code — an interval spanning zero
  places the tumours on neither side — with the sign-stability check stored alongside.

## 8. Carried to the branch notes

The two fat-pad corollaries were applied here and both earned their place: the loading table
was computed **inside the 21 samples the tests run on** (§2), and the composition adjustment
carried a **positive control read first** (§3). The control held, which is what makes §3
readable — and §2 adds a third pattern worth keeping: **when the intervention itself moves the
composition markers, a pooled loading can look "material" while the per-group loadings
disagree.** Check the per-group breakdown before believing a pooled confound.
