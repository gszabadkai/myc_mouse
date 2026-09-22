---
date: 2026-09-10
tags: [project/myc_mouse, experimental-cohorts, orthotopic, escape-series, pgc1a, bbc3, p53, selection]
status: settled -- R4 CLOSES C1; R2 and R3 overturn the correction document's section 2 explanation and replace it with something better; R1 holds on the letter of its rule but not on its substance
relates-to:
  - scripts/52_orthotopic_escape_series.R  (writes results/orthotopic_escape_series.rds)
  - docs/2026-09-09_orthotopic_identity_correction.md  (the spec this tests)
  - data/orthotopic_series/README.md
  - docs/experimental_cohorts_branch_notes.md
---

# The orthotopic cohort as a dose-of-escape series: C1 closes, and `Bbc3` tracks PGC1a dose rather than the escape route

Author-run 2026-09-10. Every number is read from `results/orthotopic_escape_series.rds`, which
reproduces the pre-run dry run in **23 of 23 non-date objects**. All 49 vector + CRISPR samples
in one scoring run, seven arms of seven.

**This is a selection experiment.** Every arm is a survivor population: PGC1a kills MYAZ cells,
so a PGC1a-expressing line is selected before it is ever injected. Nothing below is a statement
that PGC1a regulates anything; it is a statement about what tumours that survived PGC1a
expression look like.

| rule | outcome |
|---|---|
| **GATE** | no composition axis is material; but the control passes **trivially** — see §1 |
| **R1** | holds on the letter of the rule (`ox_lvl` spans zero) and **not on its substance** — §2 |
| **R2** | `Bbc3` **falls in the construct-free arm too**: the §2 p53 account is insufficient — §3 |
| **R3** | only 1 of 5 p53 targets moves with `Bbc3`: **the escape-route explanation is not supported** — §3 |
| **R4** | endogenous `Bcl2l1` is **null**: this **closes C1** — §4 |

---

## 1. The gate is clean, and the positive control passes for a reason worth stating

No compartment reaches |rho| >= 0.5 against `ox_rel` in the 49: `endothelial` +0.472,
`adipose` +0.329, `immune` +0.219, `proliferation` −0.148, `stromal` −0.136, `epithelial`
−0.011. `material_axes` is empty.

**Which means the positive control passes trivially, and that must not be over-read.** With no
material axis, `adjusted_for` is `"nothing"` and the "adjustment" is a no-op — the control
confirms the *unadjusted* contrast (`Ppargc1a` +301 CPM, p 6.2e-10 in the vector series; +137,
p 3.5e-04 in NT). It is the right outcome and it is not the test it would have been had
something been material. The per-arm breakdown is reported regardless, per the standing lesson,
and it is uninformative in the honest way: at n = 7 per arm the per-arm rhos swing from −0.857
to +0.75 for `adipose` alone.

## 2. R1 — the rule holds as written; the data are less comfortable than that

`NTPgc1a` vs `NTEV`, the construct-free PGC1a contrast:

| ruler | difference [95% boot] | Wilcoxon p |
|---|---|---|
| `ox_lvl` (absolute) | +0.257 [**−0.003**, +0.773] | 0.128 |
| **`ox_rel`** (share) | **+0.176 [+0.039, +0.348]** | **0.011** |
| **mitoPPS OXPHOS** | **+0.053 [+0.030, +0.116]** | **0.004** |
| `ox_mt` | −0.338 [−0.728, +0.229] | 0.259 |
| `ox_den` (the denominator) | −0.023 [−0.138, +0.405] | 0.710 |

R1 was specified on `ox_lvl`, and on `ox_lvl` the interval spans zero — by the letter, **§3's
dose-response reading holds** and that is what the verdict table prints.

**But two of the five rulers exclude zero, and they are the share rulers.** With `ox_den` flat,
the pattern in `NTPgc1a` is a **tilt**, not mass expansion: the respiratory arm rises *relative
to* the rest of the compartment without the compartment growing. Contrast `Pgc1a_BclxL`, where
everything moves together (`ox_lvl` −0.359 -> +2.06, `ox_rel` −0.241 -> +1.06, `ox_mt` −0.824 ->
+1.88) — that is mass expansion.

And `Myc` does not respond in `NTPgc1a`: median 2286 against `NTEV`'s 2208, with C2 giving
p = 0.50 and p = 0.37 (§5). So the honest statement is **a small but real respiratory *share*
increase with no MYC cap** — which is closer to R1's failure condition than the verdict line
conveys.

**Recommendation: §3's table stands, with its middle row reworded.** `NTPgc1a` is not "no
respiratory rise, therefore no pressure"; it is "a tilt without mass expansion, and no cap".
Whether a tilt of +0.18 constitutes selective pressure is a biological judgement this dataset
cannot make. The 0.5 threshold that produced the word "SMALL" in the verdict table was chosen
in code by Claude and is not from the spec.

## 3. R2 and R3 — the section 2 explanation fails, and what replaces it is stronger

**R2. `Bbc3` falls in the construct-free arm too**: −0.345 [−0.890, −0.010], p 0.026 in
`NTPgc1a` vs `NTEV`. By the pre-declared rule that makes the §2 p53 account **insufficient** —
the fall is not confined to the arm with a downstream buffer.

**R3. And the p53 programme is not moving with it.** In `Pgc1a_BclxL` vs `BclxL`:

| gene | difference [95% boot] | direction |
|---|---|---|
| **`Bbc3`** | **−0.731 [−1.120, −0.450]** | **down** |
| `Cdkn2a` | −0.210 [−0.595, −0.100] | down |
| `Bax` | −0.495 [−1.248, −0.135] | down |
| `Trp53` | −0.151 [−0.210, +0.001] | spans zero |
| `Zmat3` | −0.293 [−2.158, +0.281] | spans zero |
| `Cdkn1a` | +0.191 [−0.335, +0.930] | spans zero |
| `Mdm2` | +0.024 [−0.369, +0.374] | spans zero |
| `Ccng1` | +0.432 [−0.943, +0.932] | spans zero |

**One of the five named targets moves down.** `Bbc3` has the largest fall in the panel by a
factor of three and is not accompanied by the programme it is supposed to belong to. The
escape-route explanation is **not supported**.

**What the data show instead is a dose-response, and it is the better finding.** `Bbc3` across
the three within-series contrasts, ordered by the PGC1a dose the arm retained:

| arm | `Ppargc1a` CPM | `Bbc3` vs its own control |
|---|---|---|
| `Pgc1a_BclxL` | 301.6 | **−0.731 [−1.120, −0.450]** |
| `NTPgc1a` | 128.6 | **−0.345 [−0.890, −0.010]** |
| `KOPgc1a` | 0.17 (transgene lost) | +0.044 [−0.231, +1.180] |

**Monotone in dose, across three independent backgrounds, each against its own matched control**
— a polyclonal pool and two separate CRISPR clones. `Bbc3` tracks retained PGC1a dose, and it
falls to nothing where the transgene is lost. That is what R2's alternative branch predicted,
and it is a far more robust statement than the single-arm escape-route reading it replaces.

**Say it as selection.** The correct form is *tumours that retained more PGC1a show less `Bbc3`*
— not that PGC1a represses `Bbc3`. Nothing here can distinguish a regulatory effect from
selection against high-`Bbc3` cells at high PGC1a dose, and the survivor framing makes the
second at least as likely.

**One curiosity to park.** `Cdkn2a` goes **up** in `KOPgc1a` vs `KOEV` (+0.587, p 0.004) — in
the p19^ARF-knockout background. Mouse `Cdkn2a` pools p16^INK4a and p19^ARF from alternative
first exons and a CRISPR indel does not remove transcript, so gene-level counts cannot read this
and the western separates what these counts cannot.

## 4. R4 — endogenous `Bcl2l1` is null, and that closes C1

`NTPgc1a` vs `NTEV` is the only place in the cohort where `Bcl2l1` means endogenous BCL-XL:

**`Bcl2l1` −0.113 [−0.283, +0.174], p 0.535 — spans zero.** `Mcl1` is null in the same contrast
(−0.100 [−0.229, +0.047]).

C1 claimed PGC1a raises the guardian. In the one arm where the claim is testable without a
construct in the way, **it does not.** This **closes C1** rather than leaving it merely void on
a technicality: the identity error made the original contrast two-factor, and the construct-free
test that replaces it is null.

For completeness and behind its flag: `Bcl2l1` in `Pgc1a_BclxL` vs `BclxL` is −0.667
[−0.890, −0.195]. **Uninterpretable** — construct and endogenous transcript are not separable at
gene level, and BCL-XL protein is equal between those arms by western, so the transcript
difference is not differential construct expression either.

The rest of the construct-free family, none of which was pre-specified and all of which is
reported without per-gene FDR: `Bcl2` +0.344 [+0.052, +0.552] up, `Bcl2a1a` +0.810
[+0.159, +1.250] up, `Bid` −0.182 [−0.234, −0.018] down, and eight spanning zero including
`Foxo3` (−0.081) and `Bcl2l11` (+0.165). `Bik` (0.40 CPM) and `Hrk` (0.03 CPM) are below the
1 CPM floor and were reported, not fitted.

## 5. C2 on the escape series — the cap appears only where the dose was kept

Exact null, all 3432 relabellings, both statistics, both normalisations:

| contrast | max (CPM) | p | 80th pct (CPM) | p | 80th pct p, median-of-ratios |
|---|---|---|---|---|---|
| **`Pgc1a_BclxL` vs `BclxL`** | **−613** | **0.035** | **−802** | **0.012** | **0.029** |
| `NTPgc1a` vs `NTEV` | −317 | 0.500 | −24 | 0.366 | 0.467 |
| `KOPgc1a` vs `KOEV` | +2036 | 0.769 | −84 | 0.298 | 0.427 |

**The MYC ceiling is present only in the arm that kept the full PGC1a dose under a downstream
buffer**, and it survives median-of-ratios (p 0.029), which matters because that arm carries
**14.7% of its library in MitoCarta genes against 7.8-8.6% everywhere else** — it nearly doubles
its mitochondrial share, deflating every other gene's CPM without any gene changing. Both
normalisations agree, so the cap is not CPM deflation.

`KOPgc1a` carries a `Myc` maximum of **4611**, far above every other arm's, which is what makes
its max contrast positive and uninterpretable. Park it as an outlier, not a result.

**Written as an observation, not a test.** n = 7 per arm, three contrasts, two statistics, two
scales, no multiplicity correction.

## 6. What this licenses, and what it does not

**Licenses:**
- **C1 is closed**, on the construct-free contrast, not merely voided by the identity error (§4).
- **`Bbc3` falls in dose-order with retained PGC1a across three backgrounds** (§3) — as a
  selection statement, not a regulatory one.
- **The MYC cap tracks the arm that kept the dose** (§5), normalisation-robust, as an observation.
- **The four-arm escape reading of §3 survives**, with its `NTPgc1a` row reworded from "no
  respiratory rise" to "a tilt without mass expansion, and no cap" (§2).

**Does not license:**
- **The §2 p53 / escape-route explanation of `Bbc3`.** One of five targets moves. It needs a
  different account, and §3's dose-response is the candidate.
- Anything resting on `Bcl2l1` in the 14 construct-carrying samples, or on Bcl-xL vs Bcl-xS
  anywhere — no transcript-level quantification exists, so `tximport` is unavailable and the
  length offset is lost for every number here.
- **C3, the guardian ratio and the guardian-sensitiser gap** — closed, and not computed here in
  any form.
- **The MYC x OXPHOS interaction** — every arm is MYAZ-derived and MYC-high, no MYC-low arm, not
  identifiable at any n. It lives in the iMMEC rtTA-MYC +/-dox x +/-PGC1a design.
- Any comparison with the 21-sample scoring behind scripts 50 and 51, or with the timeline or
  human cohorts. Species = cohort.

## 7. Method notes

- **Scoring set fixed at 49 before any number**: 7 arms x 7, plus 18 excluded = 67, **asserted in
  PART 0** rather than described. The spec had carried "35" through five documents before anyone
  multiplied 7 by 7; the count now fails the run rather than a reader.
- **The ruler is cohort-dependent and that is asserted too.** 86 nuclear OXPHOS subunits over the
  49 against 87 over the 21, because `Cox6b2` sits on the expression floor (mean normalised count
  10.16 -> 8.17). It is a tissue-restricted paralog CLAUDE.md already flags among the five lowest
  expressers, so nothing respiratory rests on it — but it is one more reason a value from this
  run is never quotable beside one from the 21-sample run.
- 15,788 of 78,334 genes at mean normalised count >= 10. Two input objects that never mix: linear
  normalised counts for mitoPPS, `log2(+1)` for every z-composite.
- **mitoPPS** computed here for the first time on this cohort, on 141 pathways after the >=3-gene
  filter, following Monzel's algorithm as script 08 implements it. The vectorised implementation
  is **asserted equal to the definition** on a brute-force subset rather than assumed.
- MitoCarta by exact filename, Sheet 4, every `mt-*` gene stripped from every pathway as script 08
  does. Membership through `functions/reconcile_gene_symbols.R`.
- Bootstrap intervals at 5,000 resamples; the `Myc` null exact. No per-gene FDR anywhere.
- **The one threshold not from the spec:** R1's "small" was encoded as `|diff| < 0.5` or an
  interval spanning zero. That choice was Claude's, it is what produced the word "SMALL" in the
  verdict table, and §2 sets out why the substance is less comfortable than the label.
- N3 throughout: transcript associations, and the word is not applied to a transcript.

## 8. Carried to the branch notes

**Corollary 7 — a rule keyed to one ruler can pass while the cohort says otherwise.** R1 was
specified on `ox_lvl` and passes there; `ox_rel` and mitoPPS both exclude zero in the same
contrast. Pre-registration fixes the decision rule, which is its purpose, but the pre-registered
ruler is itself a choice and the other rulers still have to be read. **Report every ruler in the
contrast the rule is about, not only the one the rule names.**
