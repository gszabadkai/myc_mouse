---
date: 2026-09-08
tags: [project/myc_mouse, experimental-cohorts, orthotopic, bcl2l1, mcl1, bbc3, pgc1a, specificity, normalisation]
status: VOID (sample identity, 2026-09-09) -- the `Pgc1a` arm is PGC1a + Bcl-xL. Body preserved as the record of what was concluded; see docs/2026-09-09_orthotopic_identity_correction.md
relates-to:
  - scripts/51_orthotopic_specificity_checks.R  (writes results/orthotopic_specificity_checks.rds)
  - docs/2026-09-08_orthotopic_vector_series.md  (the result this qualifies and corrects)
  - docs/2026-09-08_script51_specificity_prompt.md  (the brief)
  - scripts/50_orthotopic_vector_series_scoring.R
  - docs/experimental_cohorts_branch_notes.md
---

# Specificity, normalisation, and one correction: what script 50's three claims look like under a closer look

> **VOID — SAMPLE IDENTITY CORRECTION, 2026-09-09.** The arm this note calls `Pgc1a` is
> **PGC1a + Bcl-xL**, not PGC1a alone, so its contrast against `EV` differs by two
> manipulations and every conclusion resting on `Bcl2l1` in that arm is void. The code is a
> correct execution of a wrong premise and is deliberately unaltered. Read
> `docs/2026-09-09_orthotopic_identity_correction.md` before this document. The body below is
> preserved as written, as the record of what was concluded and why.


Author-run 2026-09-08. Every number below is read from
`results/orthotopic_specificity_checks.rds`, which **reproduces the pre-run dry
run in all 16 non-date objects**. Same fixed cohort as 50 -- the 21 vector-series
samples, EV / BclxL / Pgc1a, 7 each -- rebuilt from source and proven against
50 before anything new was computed.

**THIS IS NOT A FOURTH CLAIM.** C1 holds, C2 holds as an observation, C3 is
uninformative, and none of those verdicts is reopened. This pass qualifies C1's
**wording**, states C2's **normalisation dependence**, and corrects C3's
**construction**. Nothing here becomes a manuscript claim.

**N3.** Transcript associations throughout. The word "primed" is not applied to
a transcript here.

| rule | outcome |
|---|---|
| **R1** `Mcl1` specificity | **SELECTIVE.** `Bcl2l1` +1.539 log2, `Mcl1` **-0.020 [-0.154, +0.127]**. The guardian-dissociation sentence is available |
| **R2** `Bbc3` gap | **GAP WIDENS.** `Bbc3` **-0.727 [-1.299, -0.376]** -- it falls while `Bcl2l1` rises |
| **R3** asymmetry | **HOLDS.** BclxL vs EV: `ox_lvl` +0.093, `ox_rel` -0.023, `ox_mt` -0.246, all three intervals covering zero |
| **R4** C2 under median-of-ratios | **WEAKENS.** 80th-pct vs EV p 0.0029 -> **0.0385**. Survives at 0.05; the mechanism is measured in section 4 |
| **R5** threshold or dose | **UNINFORMATIVE**, and said so rather than read. `Ppargc1a` spans **1.51-fold** within the arm |
| **PART A** cancelling triangle | **CONFIRMED**, and quantified |
| **PART F** C3 construction | **CORRECTION RECORDED.** C3 pooled the arm C1 excluded |

---

## 0. Gate 0, and the limitation that matters more than the length offset

`/Users/gs/G/data/MK_myc_2022` was searched exhaustively, no depth limit, for
`*length_scaled*`, `*transcript_counts*`, `*gene_lengths*`, `quant.sf` and
`tx2gene*`. **None exist.** The only salmon artefact anywhere is
`salmon.merged.gene_counts.tsv` in **three byte-identical copies**, MD5
`fba3f8fb41de4d4a2c2b76091f933c69`: the delivery folder
(`MK_paper_2026/OneDrive_1_15-03-2026/Bulkseq/`), the working copy
(`orth_tumour_data/`), and this repo's snapshot.

**This is stronger than "not found". The delivery folder is a RESULTS folder,
not a pipeline tree** -- the gene matrix plus thirteen downstream MSigDB CSVs
and two empty directories. No `star_salmon/`, no per-sample directories, no
`.rds`. There is no truncated pipeline output to return to; **these files were
never delivered**, and re-obtaining them means going back to whoever ran the
pipeline. Section 1 of the 09-08 note stands unchanged.

### The Bcl-xS problem, which the search establishes rather than merely raises

Gene-level `Bcl2l1` **pools Bcl-xL and Bcl-xS**. They arise from alternative 5'
splice site usage at the same locus and are functionally **opposite** -- long
anti-apoptotic, short pro-apoptotic. Transcript-level quantification would
separate them and it does not exist.

**C1 rests entirely on `Bcl2l1`.** So *"PGC1a raises the guardian"* cannot be
separated from *"PGC1a shifts splicing toward the short pro-apoptotic isoform"*
in these data, and no analysis choice here can separate them. It also qualifies
the BclxL arm: the construct is Bcl-xL cDNA, but the readout cannot distinguish
construct from endogenous Bcl-xS.

**This is a limitation on C1's INTERPRETATION, not on its measurement.** The
+121.8 CPM is what it is. What it means depends on an isoform ratio nothing here
observes. **For a paper resting a causal claim on `Bcl2l1` this is a sharper
limitation than the lost length offset, and it is the one a referee finds
first.** Sections 2 and 6 both carry it.

---

## 1. C1 is more specific than script 50 could show (R1, R2)

Pgc1a vs EV, n = 7 v 7, all on the log2 normalised scale so magnitudes are
comparable:

| measure | median Pgc1a | median EV | difference [95% boot] | p |
|---|---|---|---|---|
| **`Bcl2l1`** -- anchor | 12.44 | 10.90 | **+1.539 [+1.122, +1.799]** | 0.00029 (1-s) |
| **`Mcl1`** -- R1 | 12.94 | 12.96 | **-0.020 [-0.154, +0.127]** | 0.69 (1-s) |
| **`Bbc3`** -- R2 | 8.24 | 8.96 | **-0.727 [-1.299, -0.376]** | 1.00 (1-s) |
| `log2(Mcl1/Bcl2l1)` -- two-sided | 0.54 | 2.05 | **-1.506 [-1.777, -1.175]** | 0.00058 |
| `Bcl2l1` -- BclxL arm | 13.11 | 10.90 | +2.211 [+1.586, +2.470] | 0.00029 |
| `Mcl1` -- BclxL arm | 12.87 | 12.96 | -0.086 [-0.234, +0.033] | 0.95 |

**R1: SELECTIVE.** `Mcl1` does not move -- the interval is `[-0.154, +0.127]`,
tight and centred on zero, against a `Bcl2l1` move of +1.539. PGC1a moves the
guardian and leaves its partner alone. **The guardian-dissociation sentence is
available.** `Mcl1` is equally flat in the BclxL arm, which is the internal
consistency check.

**R2: the gap widens on both limbs.** `Bbc3` does not merely stay flat, it
**falls** by 0.727 log2 with an interval excluding zero. So the intervention
raises the guardian *and* lowers the sensitiser. The `log2(Mcl1/Bcl2l1)` group
contrast that 50 never ran comes out at **-1.506**, driven entirely by
`Bcl2l1` since `Mcl1` is flat.

### The family panel, and the honest complication in it

All thirteen, Pgc1a vs EV, two-sided, **no per-gene FDR** -- this is a panel and
not a set of tests. `Bcl2a1b` and `Bcl2a1d` are reported separately as the
non-one-to-one mouse paralogues they are and are never summed. Only the six with
intervals excluding zero are listed; the other seven cover zero.

| gene | side | difference [95% boot] |
|---|---|---|
| **`Bcl2l1`** | anti | **+1.539 [+1.122, +1.799]** |
| **`Bik`** | pro | **+0.765 [+0.304, +1.215]** |
| `Bad` | pro | +0.183 [+0.039, +0.473] |
| `Bcl2l11` | pro | -0.212 [-0.409, -0.004] |
| `Bcl2l2` | anti | -0.268 [-0.401, -0.100] |
| **`Bbc3`** | pro | **-0.727 [-1.299, -0.376]** |

**Two readings, and both go in.**

**The specificity is real and it is sharper than R1 alone says.** Among the
**anti-apoptotic** members, only `Bcl2l1` rises: `Mcl1` is flat and **`Bcl2l2`
(Bcl-w) falls**. This is therefore *not* a general anti-apoptotic response --
the other guardians go down or nowhere while `Bcl2l1` goes up by nearly three
fold. `Bcl2l1` moves twice as far as the next-largest member in either
direction.

**But the sensitisers do not move as a block, and that must be said.** `Bbc3`
falls, `Bcl2l11` falls slightly -- and **`Bik` rises by +0.765**. R2 was
pre-specified on `Bbc3` and its verdict stands as scored, but "the
guardian-sensitiser gap widens" is a statement about `Bbc3`, not about the
sensitiser arm as a whole. **Do not write it as the latter.**

---

## 2. What PGC1a actually did to the compartment (PART C decomposition)

The 09-08 note left the OXPHOS-selectivity of the intervention as an inference
from the gap between `ox_lvl` +1.925 and `ox_rel` +0.915. Decomposed:

| component | Pgc1a vs EV | BclxL vs EV |
|---|---|---|
| numerator: nuclear OXPHOS subunits | **+1.925 [+1.663, +2.296]** | +0.093 [-0.298, +0.581] |
| denominator: rest of nuclear MitoCarta | **+1.011 [+0.801, +1.369]** | +0.218 [-0.224, +0.565] |
| `ox_rel` = numerator - denominator | **+0.915 [+0.761, +1.007]** | -0.023 [-0.165, +0.096] |

**PGC1a raised the whole mitochondrial compartment, with roughly a 2:1 OXPHOS
tilt.** The denominator moves +1.011, more than half as far as the numerator.
So this is a **broad biogenesis intervention with an OXPHOS bias**, not an
OXPHOS-selective one, and C1's manipulation should be described that way.
`ox_rel` is positive because the tilt exists, not because only OXPHOS moved.

---

## 3. R3 -- Bcl-xL does not move the respiratory transcriptome (blocker deleted 2026-09-09)

BclxL vs EV, two-sided:

| ruler | difference [95% boot] |
|---|---|
| `ox_lvl` | +0.093 [-0.298, +0.581] |
| `ox_rel` | -0.023 [-0.165, +0.096] |
| `ox_mt` (mtDNA-encoded) | -0.246 [-0.547, +0.303] |

**All three cover zero.** Overexpressing Bcl-xL does not move respiration
transcriptionally, on any of the three rulers.

> **CORRECTION, 2026-09-09 — the blocker that stood here has been deleted.** It
> asserted a "standing experimental finding that Bcl-xL overexpressing tumours
> have higher OXPHOS" and forbade the argument until that was reconciled. **No
> such finding exists.** It entered in the brief that commissioned this script,
> `docs/2026-09-08_script51_specificity_prompt.md:149-152`, and was never in the
> data or in the manuscript. The manuscript has Bcl-xL as the **permissive**
> partner throughout — it rescues the Ndi1 growth defect and makes raised
> respiration advantageous — and nowhere claims it raises respiration. R3 is
> therefore consistent with the manuscript's own model and there is nothing to
> reconcile. See `docs/2026-09-09_orthotopic_identity_correction.md` section 1.
>
> **R3 is unblocked but DEMOTED.** The sentence above it — *"the edge runs one
> way, and the directed-edge argument is available"* — is withdrawn: with C1 void
> there is no forward edge for R3 to be the reverse of. R3 stands as a
> **standalone negative control** (Bcl-xL does not move the respiratory
> transcriptome), **transcript-level only**, while the manuscript's respiration
> claims are functional.

---

## 4. R4 -- C2 weakens under median-of-ratios, and the mechanism is measured

Same exact permutation null, all `choose(14,7)` = 3432 relabellings, both scales:

| statistic | contrast | CPM p | median-of-ratios p |
|---|---|---|---|
| max | Pgc1a vs EV | 0.0962 | 0.0962 |
| **80th pct** | **Pgc1a vs EV** | **0.0029** | **0.0385** |
| max | Pgc1a vs BclxL | 0.0350 | 0.0350 |
| 80th pct | Pgc1a vs BclxL | 0.0120 | 0.0288 |

**Every cell keeps its sign and stays under 0.05 where it was under 0.05, but
the strongest cell drops by thirteen-fold.** The max contrasts are unchanged;
both 80th-percentile contrasts weaken.

**And the mechanism is not a matter of argument, it is measured:**

| arm | size factor | library (M) | sf / libM | **MitoCarta share of library** |
|---|---|---|---|---|
| EV | 0.998 | 26.8 | 0.0371 | **7.96%** |
| BclxL | 1.034 | 27.4 | 0.0368 | **7.75%** |
| **Pgc1a** | 0.985 | 28.4 | **0.0350** | **14.73%** |

**PGC1a nearly doubles the mitochondrial share of the transcriptome.** With
14.7% of the library sitting in MitoCarta genes, **every other gene's CPM in
that arm is deflated relative to EV** -- and `Myc` is one of those genes. That
is exactly the compositional route by which a CPM-based ceiling can be partly
manufactured by the intervention that defines the arm.

**So C2 stands as an observation, with the dependence stated.** It is not
abolished -- the 80th-percentile contrast against EV is still p = 0.0385 on the
normalisation that does not have this problem, and the against-BclxL contrast
survives too. But **the number to quote is 0.0385, not 0.0029**, and the
sentence must say which normalisation it is on. The 09-08 note's `p = 0.0029`
should not be quoted without this beside it.

---

## 5. C3: two independent reasons it is uninformative, and neither is power

### 5.1 PART A -- the arms cancel, and that is now a number

Arm means:

| arm | `ox_rel` | `ox_lvl` | `Bcl2l1` | `Mcl1` | `log2(Mcl1/Bcl2l1)` |
|---|---|---|---|---|---|
| EV | **-0.292** | -0.694 | 11.0 | 13.0 | +1.98 |
| BclxL | **-0.291** | -0.556 | 13.0 | 12.8 | -0.22 |
| Pgc1a | **+0.583** | +1.250 | 12.5 | 13.0 | +0.49 |

**EV and BclxL sit at the same point on the respiratory axis** -- `ox_rel`
-0.292 against -0.291, a difference of 0.001 -- **while their `Bcl2l1` differs
by 2.1 log2 units.** That is the construct, and it is a vertical displacement at
constant `ox_rel`. Pgc1a then sits high on `ox_rel` with intermediate `Bcl2l1`.
**The three arm-means are NON-MONOTONE: the cancelling-triangle explanation is
CONFIRMED.**

Decomposed:

| endpoint | pooled (n = 21) | between arms (n = 3) | within EV | within BclxL | within Pgc1a | within mean |
|---|---|---|---|---|---|---|
| `Bcl2l1` | **-0.065** | +0.500 | -0.143 | -0.357 | -0.643 | **-0.381** |
| `Mcl1` | +0.279 | -0.500 | +0.714 | +0.107 | +0.250 | +0.357 |
| `guardian_ratio` | **+0.064** | -0.500 | +0.750 | 0.000 | +0.821 | **+0.524** |

**The between-arm and within-arm components have opposite signs and cancel to
approximately zero.** The pooled near-zero is not an absence of association --
it is two components pointing in opposite directions.

**And this identifies which reading C3 actually wanted.** C3 needs `ox_rel` to
vary *within* a common state, which is the within-arm component -- and that is
script 50's group-adjusted reading, **+0.300 [-0.255, +0.737]**, whose interval
spans zero at n = 7 per arm. **So the verdict does not change.** What changes is
that the near-zero pooled value is now explained rather than merely reported.
The between-arm rho is on n = 3 and is descriptive only; no interval is invented
for it.

Figure: `outputs/orthotopic/51_arm_geometry.pdf`.

### 5.2 PART F -- the construction correction

**C1's spec said BclxL is reported and NEVER pooled. C3's spec pooled it.**
Seven of the 21 `Bcl2l1` values are therefore set by an **expression construct**
rather than by the tumour. Same endpoint, opposite handling, caught only because
C1 and C3 were pre-specified as separate claims.

Recomputed on EV + Pgc1a alone -- **scores subset from the 21-sample cohort and
never re-scored**, because `ox_rel` is cohort-relative:

| endpoint | reading | n | rho [95% boot] | excludes zero |
|---|---|---|---|---|
| `guardian_ratio` | pooled, all 21 | 21 | +0.064 [-0.360, +0.489] | no |
| `guardian_ratio` | **EV + Pgc1a only** | 14 | **-0.560 [-0.689, -0.101]** | **yes** |
| `Mcl1` | pooled, all 21 | 21 | +0.279 [-0.146, +0.644] | no |
| `Mcl1` | EV + Pgc1a only | 14 | +0.099 [-0.407, +0.592] | no |
| `Bcl2l1` | pooled, all 21 | 21 | -0.065 [-0.511, +0.378] | no |
| `Bcl2l1` | **EV + Pgc1a only** | 14 | **+0.657 [+0.250, +0.823]** | **yes** |

> **DO NOT READ THIS AS C3 RESOLVING.** With the third arm removed there are two
> arms separated by a large deliberate `ox_rel` gap, so the "correlation"
> collapses into the **EV-to-Pgc1a group contrast, which IS C1**. The rho of
> +0.657 is C1 wearing a correlation's clothes, and its sign is **opposite to
> the within-arm reading** (-0.381) for the same gene, which is the tell.
> **It is circular and it is reported as a demonstration of the construction
> problem, not as a result.**

**C3 REMAINS UNINFORMATIVE. The abstract sentence stays withdrawn.** Sections
5.1 and 5.2 now give two independent reasons why, and neither of them is lack of
power: pooled across three arms the components cancel; pooled across two arms
the reading is circular. **A design built to move `ox_rel` between arms cannot
answer a question about `ox_rel` varying within one.**

### 5.3 The normalisation does not rescue it either

The composition-adjusted C3 fit under both normalisations. **A correction to the
brief, stated rather than applied silently: script 50's PART F already ran on
median-of-ratios** -- its endpoints are `log2(counts(dds, normalized = TRUE) + 1)`
-- so the brief's direction was inverted. Both are fitted here.

| endpoint | median-of-ratios (script 50) | CPM (the mirror) |
|---|---|---|
| `guardian_ratio` | -1.442 (p 0.111) | -1.431 (p 0.112) |
| `Mcl1` | -0.214 (p 0.135) | -0.384 (p 0.0185) |
| `Bcl2l1` | +1.228 (p 0.145) | +1.047 (p 0.203) |

**The endpoint carrying the sign flip is normalisation-stable** (-1.442 against
-1.431). Only `Mcl1` moves, and it is not the primary endpoint. So the
normalisation choice is not what makes C3 uninformative.

---

## 6. R5 -- the dose question is not answerable in this arm, and that is the answer

Within the Pgc1a arm, n = 7:

| measure | min | max | fold spread | CV |
|---|---|---|---|---|
| `Ppargc1a` (CPM) | 242 | 364 | **1.51** | 0.150 |
| `Bcl2l1` (CPM) | 181 | 227 | 1.25 | 0.092 |

`rho(Ppargc1a, Bcl2l1)` = **+0.286 [-0.765, +0.882]**, n = 7.

**UNINFORMATIVE, and it is written that way rather than read as a threshold.**
The transgene spans only 1.51-fold within the arm, so flatness of `Bcl2l1`
cannot discriminate a viability threshold from a dose-responsive induction. The
interval on the correlation runs from -0.77 to +0.88. **Neither the threshold
reading nor the dose-response reading is available**, and R5 was written in
advance to say so rather than to read the flatness.

**The external dose point, hedged exactly as specified.** `NTPgc1a` carries
`Ppargc1a` at **129 CPM** against `NTEV` at 0.16, with `Bcl2l1` **79.3 against
81.6**. **Raw CPM observation only** -- never scored, never tested, never in a
claim, and a **different genetic background**. It does not touch the p21 series'
out-of-scope status. It is recorded because it is in the deposited matrix and a
reader will find it: at roughly 43% of the vector series' transgene level,
`Bcl2l1` is unmoved. That is *suggestive* of a threshold and it is **not
evidence**, because the background differs and nothing was scored.

---

## 7. What this licenses, and what it does not

**Licenses, and these are wording changes to C1 rather than new claims:**

- **the specificity sentence.** PGC1a raises `Bcl2l1` and does not move `Mcl1`
  (-0.020 [-0.154, +0.127]); among the anti-apoptotic members only `Bcl2l1`
  rises, `Bcl2l2` falls. C1 is a **selective** event, not a general
  anti-apoptotic response;
- **the gap sentence, scoped to `Bbc3`.** The guardian rises and `Bbc3` falls
  (-0.727 [-1.299, -0.376]), so the gap widens on both limbs **by
  intervention**. Scoped to `Bbc3`: `Bik` rises;
- **describing the manipulation correctly.** PGC1a is a **broad mitochondrial
  biogenesis intervention with roughly a 2:1 OXPHOS tilt** (numerator +1.925,
  denominator +1.011), not an OXPHOS-selective one;
- **C2 with its normalisation named.** The 80th-percentile contrast against EV
  is **p = 0.0385 on median-of-ratios**, and that is the number to quote.

**Does not license:**

- **the directed-edge / asymmetry argument**, until section 3's open question is
  answered -- whether the standing Bcl-xL/OXPHOS finding was transcriptional or
  functional;
- **any reading of C3 in either direction**, including the +0.657 in section
  5.2, which is circular;
- **quoting C2's `p = 0.0029` without the median-of-ratios value beside it**;
- **any threshold or dose-response reading** of the Pgc1a arm (section 6);
- **"PGC1a raises the guardian" as an isoform-resolved statement.** Section 0:
  `Bcl2l1` pools Bcl-xL and Bcl-xS and nothing here separates them;
- **anything about the `MYC x OXPHOS` interaction.** Not identifiable here at
  any n; it lives in the iMMEC rtTA-MYC +/-dox x +/-PGC1a design on the death
  readout;
- **any pooling** with the 6W/12W timeline or the human cohorts. Species =
  cohort; only directions travel.

---

## 8. Method notes

- **PART 0 is the gate and it passed.** The cohort was rebuilt from source --
  load, longest-first group parsing, MitoCarta Sheet 4 by exact filename with
  every `mt-*` stripped, both input objects, both rulers -- and then proven:
  the median-CPM table reproduces the plan's section 2 to **0.784%**, and script
  50's five C1 point estimates reproduce to **0.0238%** with all five at the
  7 v 7 one-sided floor p = 1/3432. Set sizes asserted at **87 / 880 / 13** with
  15,721 genes carried. Nothing new was computed until all of that passed.
- Bootstrap intervals are **reported beside** 50's published ones but **not
  asserted**: they depend on RNG position in the stream and 51 does not consume
  randomness in 50's order. The point estimates and the Wilcoxon p are
  deterministic and are what the control tests.
- 5,000 bootstrap resamples throughout; the C2 nulls are **exact**, all 3432
  relabellings. **No per-gene FDR** across the thirteen family members.
- The composition gate reproduced 50's material axes, `adipose` and
  `endothelial`, and section 5.3's adjusted fits carry them.
- **mitoPPS and GSVA were not run.** No pre-specified reading needs either.
- Two input objects that never mix: linear normalised counts, and
  `log2(+1)` of them.

---

## 9. Carried to the branch notes

**The gate's pooled |rho| threshold has a blind spot in BOTH directions**, and
this dataset now shows both halves.

- Section 2 of the 09-08 note showed the first: **pooled-material with the
  per-group loadings disagreeing** -- `adipose` +0.621 pooled against
  +0.036 / -0.393 / +0.714 by arm, because the intervention itself moves the
  composition markers.
- `ox_lvl` against proliferation shows the mirror: **pooled -0.323, below
  threshold, with EV / BclxL / Pgc1a at -0.643 / -0.679 / -0.964** -- three arms
  agreeing strongly and the pooled value not clearing the bar.

**Read the per-group breakdown before believing either verdict.** A pooled
loading can be material where nothing is consistent, and immaterial where
everything is.

**And a second corollary from PART A, which is the same lesson on the endpoint
rather than the covariate:** a near-zero pooled correlation across designed arms
can be **two opposed components cancelling**, not an absence. Decompose into
between-arm and within-arm before reading a pooled rho as a null.

---

## 10. Where the numbers live

| | |
|---|---|
| script | `scripts/51_orthotopic_specificity_checks.R` |
| object | `results/orthotopic_specificity_checks.rds` -- `$load_control` `$arm_geometry` `$c1_panel` `$asymmetry` `$norm_diagnostic` `$c2_mor` `$c3_partf_mor` `$dose` `$c3_ev_pgc1a` `$verdict` `$gate0` `$notes` |
| figure | `outputs/orthotopic/51_arm_geometry.pdf` |
| the gate | `$load_control` -- 0.784% on the median-CPM table, 0.0238% on 50's C1 |
| what it qualifies and corrects | `docs/2026-09-08_orthotopic_vector_series.md` |
| the brief | `docs/2026-09-08_script51_specificity_prompt.md` |
