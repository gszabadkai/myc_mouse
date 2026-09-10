# CORRECTION v2 — orthotopic series identity, provenance and interpretation

**2026-09-09.** Supersedes v1 of the same date. v1 was written before the
derivation, the western blots, the growth data and the surviving PGC1α-only arms
were known. Two of its central claims are withdrawn below.

Still standing from v1: the arm labelled `Pgc1a` is **PGC1α + Bcl-xL**, and every
conclusion in scripts 50 and 51 that rests on `Bcl2l1` in that arm is void.

---

## 1. Withdrawn from v1

**"The empty cell is the result."** v1 stated that PGC1α-only tumours do not form
unless the transgene is silenced. **False.** Three PGC1α-only arms exist and all
three formed tumours. What varies is how much PGC1α they kept. The argument moves
from *existence* to *dose* — weaker as a headline, considerably more informative
as data. See §3.

**"Standing finding that Bcl-xL tumours have higher OXPHOS."** I introduced this
and it is not supported anywhere. The manuscript has Bcl-xL as the *permissive*
partner throughout — it rescues the Ndi1 growth defect, recovers the Ndi1
metastatic defect, and makes raised respiration advantageous. It is never claimed
to raise respiration. **The blocker in `docs/2026-09-08_orthotopic_specificity_checks.md`
§3 must be deleted**; there is nothing to reconcile. R3 (`ox_lvl` +0.093,
`ox_rel` −0.023, `ox_mt` −0.246, all covering zero) is consistent with the
manuscript's own model.

R3 is also demoted. With C1 void there is no forward edge, so it is a standalone
negative control — Bcl-xL does not move the respiratory transcriptome — not a
directed-edge argument. And it is transcript-level only; the manuscript's
respiration claims are functional.

**Also withdrawn, from chat rather than v1:** the proposal to anchor the human
arm on CDKN2A. iMMEC carries p53DN + E1A + SV40-LT, and PGC1α still induces PUMA
and MYC plus respiration still kills — so the p53–ARF route is one path to the
endpoint, not the mechanism. It is also known biology, and Menegollo already
reported CDKN2A loss with the OXPHOS-high biclusters, so running it in TCGA
confirms the companion paper rather than advancing past it.

---

## 2. Resolved

**Clonality — closed.** `BclxL` and `Bcl-xL PGC1a` are the **same pool**: MYAZ
transduced with Bcl-xL, then split, then EV or PGC1α, polyclonally selected.
BCL-XL protein is equal on the blot. So **`Bcl-xL PGC1a` vs `Bcl-xL EV` is a
valid PGC1α contrast** for every gene except `Bcl2l1` itself.

Densitometry is still worth having, since a blot does not resolve the 1.6× the
transcript shows — but the derivation fact does the work clonality testing was
going to do.

**`Bcl2l1` 326 → 193 — deprioritised.** Equal protein means the transcript
difference is not differential construct expression. Gene-level counts cannot
separate construct from endogenous in either arm, so the number is not
interpretable and does not need to be.

**`Bbc3` −0.727 — explained, and it stops being an anomaly.** MYAZ is p53
proficient; the PGC1α survivors have **p19^ARF and p21 down** by western.
Unrestrained MDM2, reduced p53 output, and `Bbc3` — a canonical p53 target —
comes down with the rest of the programme. The fall is a readout of *which escape
route this arm took*, not a PGC1α regulatory effect.

That also resolves the transcript–protein inversion. PUMA protein is up (Supp
Fig 5A) while the transcript is down, and sequestration is excluded because
BCL-XL protein is equal between arms. Reduced transcriptional induction with
retained post-transcriptional accumulation fits both observations without either
giving way.

**Prediction, testable now in the existing counts.** If the p53 explanation
holds, `Cdkn2a`, `Cdkn1a`, `Mdm2`, `Ccng1` and `Zmat3` should move down
**together** with `Bbc3` in `Bcl-xL PGC1a` vs `Bcl-xL EV`. Caveat: mouse
`Cdkn2a` encodes p16^INK4a and p19^ARF from alternative first exons and
gene-level counts pool them, so the WB separates what the RNA-seq cannot.

---

## 3. The design as it actually is — a dose-of-escape series

Not a 2 × 2. Four PGC1α arms, each with a matched control, differing in how much
of the death pathway was broken and therefore how much PGC1α survived:

| arm | escape route | `Ppargc1a` CPM | OXPHOS | growth | `Myc` median vs its control |
|---|---|---|---|---|---|
| MYAZ + PGC1α *(no RNA-seq)* | transgene silenced | — | none | = control | — |
| `KOPgc1a` (p19^ARF KO) | transgene lost entirely | **0.2** | ? | = NT; KO advantage cancelled | 2176 vs 2368 |
| `NTPgc1a` | partial dose retained | **128.6** | slight ↑ (WB) | = NT-EV | **2286 vs 2208** |
| `Pgc1a` = **Bcl-xL PGC1α** | buffer downstream, full dose kept | **301.6** | large ↑ | **> Bcl-xL EV** | **1440 vs 2033** |

Read down the table: **the amount of PGC1α a tumour keeps is set by how much of
the death pathway it broke.** Only the arm with a downstream buffer retains the
full dose and the full respiratory increase — and it is the only arm that grows
faster and the only arm where `Myc` comes down.

**What they were escaping from is now named.** The LbNOX result (2026-09-09)
identifies the lethal variable as **mitochondrial matrix NADH oxidation on a MYC
background**: mtLbNOX kills MYAZ + Tre-MYC while making no ATP, cytLbNOX kills
much less, and Ndi1 behaves as mtLbNOX does. So the table above resolves into
**two escape classes** — three arms lowered the oxidative driver, one blocked the
effector downstream and kept the driver. Both are escapes from the same gate. See
`2026-09-09_reprioritisation_narrative_v3.md` §2 and §5.

**This is a better C2 than the two-arm contrast ever was.** The MYC reduction
appears only where OXPHOS actually rose, across four arms with three separate
controls. Medians at n = 6–7 in differing backgrounds, so an observation rather
than a test — but a dose-response observation, and one that survives the
objection that killed the original C2.

`NTPgc1a` not being capped becomes a supporting datum rather than a problem: no
respiratory pressure, nothing to select against. **Pre-declared reading rule:**
this holds only if `ox_lvl` in `NTPgc1a` is small. If it moves substantially with
no MYC effect, the dose-response reading fails.

The `KOPgc1a` cell — PGC1α cancels the ARF-loss growth advantage while its own
transcript is null — is curious and should be parked.

---

## 4. `NT-EV` vs `NT-PGC1a` is now the highest-value pull in the cohort

It is the **only construct-free PGC1α contrast** in the dataset, and therefore the
only place `Bcl2l1` means endogenous BCL-XL.

From the load table: `Bcl2l1` **79.3 vs 81.6** — flat. It needs a fitted contrast
with an interval, but medians that close at n = 7 will not become a 2.9× rise.
Whatever remained of C1's direction, the construct-free test looks negative.

The pull, in order:

1. **`ox_lvl`, `ox_rel`, `ox_mt`, mitoPPS** in `NTPgc1a` vs `NTEV`. Decides §3's
   reading rule, and says whether this arm is mass expansion like the double or
   something closer to the timeline's reprioritisation.
2. **`Bbc3`.** Pre-declare: rising or flat here means the double's fall requires
   sustained high dose under a buffer and is an escape-route signature. Falling
   here too means it tracks PGC1α at any dose and the p53 explanation is
   insufficient.
3. **The twelve, construct-free.**
4. **`Foxo3`**, plus the p53 target panel as a background control across both
   series.
5. **`Myc` max and 80th percentile**, same exact permutation null, so C2 is
   tested where the model's prediction is now explicit.

**Scoring rule.** `ox_rel`, GSVA and mitoPPS are cohort-relative. Re-score **all
35 vector + CRISPR samples in one run** as the primary for anything involving NT,
and leave the 21-sample scoring untouched as the record behind scripts 50 and 51.
Never quote a value from one run beside a value from the other.

**Contrast rule.** Clean contrasts are **within series** — `NTPgc1a` vs `NTEV`,
`Bcl-xL PGC1a` vs `Bcl-xL EV`. NT is a CRISPR clone and the vector series a
polyclonal pool, so cross-series contrasts carry a derivation difference and are
descriptive only. Two backgrounds, one manipulation: close to a replication.

---

## 5. What the orthotopic can and cannot say

**Can:**
- PGC1α raises the mitochondrial transcriptional programme in vivo,
  dose-dependently across arms.
- Tumour growth increases when PGC1α is combined with a downstream apoptotic
  block.
- Tumours carrying MYC and high respiration exist only where something in the
  death pathway was broken, and the transcriptome records which.
- `Myc` falls only in the arm where respiration actually rose.

**Cannot:**
- Anything resting on `Bcl2l1` in the 14 construct-carrying samples.
- The guardian ratio, the guardian–sensitiser gap, or C3 in any form.
- A clean PGC1α → PUMA regulatory statement. These are survivors, and the
  transcript reads the escape route.
- The human arm's **V4** ("move OXPHOS, read the twelve"). `NTPgc1a` is the
  nearest available answer and it is a small-dose, single-clone arm.

**The interpretive rule this episode establishes, for the branch notes:** every
orthotopic arm is a survivor population. PGC1α kills MYAZ cells, so a
PGC1α-expressing line is selected before it is ever injected. Matched EV handling
removes the protocol asymmetry but not this one. **Read the orthotopic
transcriptome as a record of what survived, not as a response to a manipulation.**

Second corollary: a construct arm is not identified by its label. The manuscript
named these arms correctly — *EV EV, Bcl-xL EV, Bcl-xL Pgc1a*, Figure 5 — the
whole time. The check that would have caught this is reconciling sample labels
against the manuscript's own figure legends before pre-specifying anything.

---

## 6. Repo actions — `myc_mouse`, `experimental-cohorts`

Unchanged from v1: **do not edit or delete scripts 50 and 51.** Tag the current
commit. Do not merge to `paper-final` until 52 is signed off.

Updated:

1. **`data/orthotopic_series/README.md`** — rename `Pgc1a` → `Pgc1a_BclxL`;
   record the derivation (single Bcl-xL pool, split, polyclonal selection);
   record equal BCL-XL protein; record that `NTEV`/`NTPgc1a` are a CRISPR-clone
   background carrying a surviving PGC1α-only arm; record the WB findings per arm.
2. **`docs/2026-09-09_orthotopic_identity_correction.md`** — this document.
3. **Banner** the two 09-08 notes with a pointer. **Delete** the fabricated
   blocker in the specificity note §3 and replace it with a one-line correction
   naming its source.
4. **`scripts/52_orthotopic_escape_series.R`** — renamed from v1's plan, because
   the job changed:
   - All 35 vector + CRISPR samples, scored in one run. Labels corrected and
     asserted.
   - Primary contrasts within series.
   - The p53 target panel as the pre-declared test of §2's explanation.
   - The construct-free twelve in `NTPgc1a` vs `NTEV`.
   - The four-arm `Myc` and respiratory table as one output object.
   - Reading rules in the header before any number, as for 51.
5. **Scope change, recorded deliberately:** the p21/ARF series was out of scope in
   the 09-07 plan. `NTEV`/`NTPgc1a` carry the only construct-free PGC1α contrast
   and `KOEV`/`KOPgc1a` complete the escape series. Both return, with the reason
   stated.

**`myc_human_exploratory`:** V4 reopened. Remove any citation of orthotopic C1 as
causal support for the human `BCL2L1` coefficient. No CDKN2A analysis.

**`myc_human_validation`:** frozen, no action.

---

## 7. Open items

- Whether the in vitro PUMA qPCR and WB were on the **drug-selected population
  that went into the fat pad** or on an acute post-transduction population. This
  localises the transcript reversal to the dish or to the animal, and it is the
  single most informative unanswered question.
- Banked pre-injection cells: a three-gene qPCR (`Ppargc1a`, `Bbc3`, `Bcl2l1`)
  would settle the above directly.
- Densitometry on the BCL-XL blot.
- Which comparator Supp Fig 5A's PUMA western used — EV EV or Bcl-xL EV.
- Whether PUMA protein and `Bbc3` were measured in the **same animals**. A
  within-animal inversion is a result; two group-level observations are a puzzle.
- `Cdkn2a` in the mouse timeline — expected flat, and a stated negative protects
  the novelty claim by showing the developmental gate operates without it.
- Bcl-xL / Bcl-xS: unresolvable, no transcript-level quantification, 14 of 21
  samples affected.
