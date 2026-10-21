---
name: pgc1a-axis-did-not-change
description: "The developmental 6->12W OXPHOS decline is NOT a PGC1a-axis activity change - the regulon splits in two directions; and the TF-target layer reports cell-state CONTEXT, not factor identity, so no TF may be named from it"
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-08-17T23:50:29.636Z
---

**Closed 2026-08-17** (author's question during the Block C cut; `docs/2026-08-17_developmental_
oxphos_decline_and_the_biogenesis_axis.md`, `scripts/47_biogenesis_axis_and_the_developmental_
oxphos_decline.R`). It had been raised repeatedly across Block A and never resolved.

**THE DECISIVE NUMBER WAS ALREADY ON DISK** since July, in script 44's `core_decomp`, and had never
been read as the answer to this question. Inside `CORE_MITO` (the ERRa u NRF1 u GABP regulon,
300 genes):

| part | n | WT 6->12W content | matched-null pct |
|---|---|---|---|
| n OXPHOS subunits | 59 | **-0.216** | **0.05** |
| the rest | 241 | **+0.069** | **99.55** |
| whole regulon | 300 | +0.013 | 74.90 |

**An activity change acts on a REGULON, not on a functional subset of it.** Script 47 PART B
generalises it: **all four regulons that actually contain the respiratory chain split the same way**
(pct_OXPHOS / pct_rest -- `CORE_MITO` 0.0/99.2, `ESRRA_MITO` 0.4/97.6, `NRF1_MITO` 1.0/84.8,
`GABPA_MITO` 4.1/79.4). The other three hold 2-9 OXPHOS genes and cannot show the split.
**The one regulon that DOES move as a unit is the ER one** -- `DEVELOPMENTAL_MITO` at the 96.0th
percentile, i.e. UP. That is the shape an axis change makes; the PGC1a axis does not make it.

Three more lines agree:
the *build* arms do not follow the *run* arm down (**OXPHOS assembly factors +0.0007, 50.25th pct
-- dead centre of their own null**; mitoribosome, central dogma, mtDNA maintenance flat-to-up;
catabolic arms rise, FAO +1.89); the factors are flat (`Esrra` +0.23, `Nrf1` -0.25, `Gabpa` +0.19,
`Tfam` +0.10, all ns); and the **generic** (non-Gray) regulons go UP (`GABPA_CHIPATLAS` +1.20,
`ESRRA_CHIPATLAS` +1.19; Category 7 `CORE_MITO` +0.85 p 0.90).

**`Ppargc1a` is at baseMean 30 -- effectively absent from purified MEC.** The only PGC-1 family
member expressed at a real level is **`Pprc1` (PRC), baseMean 2102**, and it is the only one that
moves: -0.471 padj 1.9e-4 (WT), -0.585 padj 1.1e-6 (Myc+), interaction ns. **PRC appears nowhere
else in the corpus.** BUT it does **NOT** clear its expression-matched null (**7.1st percentile**,
inside the pre-set [5,95] bar) -- the tiny padj is **precision, not effect size**. A lead for an
antibody, not a mechanism. **Nothing in the 57-gene TF/coactivator roster clears the bar on the
DOWN side** (`Mybl2` 5.6th, `Nrip1` 9.2nd). What clears it is all UP: **`Foxo3` 95.7th
(padj 0.0078), `Sirt1` 97.9th (0.019)**, `Srebf1` 98.1st, `Elf5` 98.3rd, `Esrrb` 99.6th.

## The trap that kept it open: the TF layer cannot rank TFs

**Never name a transcription factor from a `_GRAY_` or `_MITO` lane in this project.**
- `_MITO` lanes are `TF programme n MitoCarta` -- an OXPHOS-heavy overlap falls **by construction**.
- `_GRAY_<context>` lanes report **cell state**: adjusted R2 **0.472 for context vs 0.267 for TF
  identity** over 352 lanes; within `AP_LE` the 87 mito lanes span an IQR of 0.50 NES against a
  between-context spread of 3.34. Inside `AP_LE`, **`CHCHD3` (a MICOS structural protein) -2.56
  and `GTF3A` (general transcription factor IIIA) -2.53 both outrank `ESRRA` -2.33.**
- The layer **ambient is -1.461** (median of 421 lanes; the whole layer is negative-going). The
  PGC1a-axis lane median is **-1.144, i.e. ABOVE ambient** (script 24). Always read a TF lane
  against ambient, never against zero.
- Script 47's two-sided reading rule passes **0 of 218** mito lanes. The only readable sub-layer is
  the 69 **context-free** lanes, where **`TFT_FOXO3_CHUNG` +1.43 is the single highest**.

## What the decline IS (hypothesis, not result)

A **cell-state switch, TEB / lineage-suppressed AP -> mature duct**: `MG_TEB_VS_DUCTAL_HS_GRAY_UP`
**-0.422 (pct 0.0)**, a bigger fall than OXPHOS's own -0.255; AP LE markers -0.294, still **-0.238
with every mito gene removed**. **NOT IDENTIFIED**: carrier and cargo correlate at **r = 0.93**, so
script 44's -0.459 -> -0.151 adjustment is a **bound, not mediation**. Sorted AP/HS/BA or
single-cell is the design that settles it.
## FOXO3: the strongest readable SIGNAL, and NOT the driver (script 47 PART I)

It earns the first half: its target set contains **0 OXPHOS subunits** (no content artifact) and
carries **no Gray context** (no state artifact), so it is the top of the 69 context-free lanes
(+1.43). `Foxo3` +0.473 padj 0.0078 (95.7th pct), `Sirt1` +0.560 (97.9th), `Bnip3` +0.635;
`FOXO1` null throughout as the control.

**It loses the second half on REACH -- 0 of the 89 declining subunits, against CORE_MITO's 59.**
A factor that touches none of the affected genes cannot be their proximal cause. FOXO3 is the
readable MARKER of the adult state, corroborating the state hypothesis rather than competing with
it; its non-mito composite is only at the **11.5th percentile** against OXPHOS and anti-correlates
with the carrier just as hard (-0.53), so it is not a separable axis either.

**Three things that are NOT true and will otherwise be said:**
1. "Myc blocks the developmental FOXO3 rise" -- **NO.** The programme rises in BOTH genotypes
   (+1.43 WT, +1.71 Myc+) and the interaction is ns. What Myc does is **suppress it by a constant
   offset at both ages** (-1.76 / -1.80, padj 0.003) -- the clean genotype result worth keeping.
2. "The adult gland arms PUMA through FOXO3" -- **NO.** `Bbc3` is **flat in the wild-type
   timeline** (+0.062, padj 0.88). PUMA falls in the *Myc+* timeline (-0.479, padj 0.016).
3. "It is a quiescence/arrest programme" -- only the **atrophy/turnover** arm moves (`Fbxo32`
   +1.16 padj 5e-10, `Bnip3`, `Sirt1`); the arrest arm (`Cdkn1a/1b`, `Gadd45a`, `Ccng2`) is flat.

Also: **`Foxo3` interaction padj = 1**, and NES is scale-free so equal NES at both ages is NOT
equal effect size ([[fgsea-takes-one-list-not-two]]).

## The tension to keep straight -- and the ceiling that goes with it

Between the ages OXPHOS falls **-0.255 against proliferation's -0.047**, five times more, so
"the gland de-respires without de-proliferating" is a **between-age** statement and needs that
scope word. Within an age the two covary at r = +0.887 -- **but never quote that raw r**: in the
same 12 samples OXPHOS <-> `CORE_MITO` is **0.980** and <-> mitoribosome **0.967**. A set-level
0.89 is the ceiling ([[correlation-ceiling]]), not a finding. **Only the gene-level percentile
against all 16,922 expressed genes is readable:** `Tfdp1` 99.8th, `Mterf4` 99.5th, `Gata3` 95.5th,
`Mybl2` 94.2nd at the top; **`Sirt1` 3.8th, `Bnip3` 5.7th at the bottom**; and the biogenesis
factors in the middle (`Gabpa` 79th, `Tfam` 65th, `Nrf1` 45th, `Pprc1` 44th, `Esrra` 43rd).

## What licenses the PGC1a/NRF1 rescue anyway (the author's follow-up)

Frame it as an **instrument, not a mimic** -- and the justification is measured:
1. The mouse **cannot de-confound** respiration from the state carrying it (r = 0.93); the
   intervention can. The limitation IS the licence.
2. **Reach:** the axis covers **59 of the 89 declining subunits (66%)**.
3. **Two factors are a control:** ERRa n NRF1 = 22 genes, **50% OXPHOS subunits, 2.07x enriched
   over their union, p = 0.0052** -- convergence localises the effect to the chain.
4. **Not circular:** the death effectors are **NOT** in the PGC1a regulon (5/25, **0.76x**, no
   PUMA/BAX/caspase); they are in MYC's (**2.68x, p = 0.011**). The one exception, `Cycs`, is
   itself a respiratory carrier.
5. **Calibration:** the deficit is **x1.19** (WT) to **x1.32** (Myc+) at subunit level. Since
   PGC1a is absent in MEC, calibrate on OUTPUT, not transgene.
6. **Pre-specify the overshoot:** assembly factors / mitoribosome / central dogma never fell but
   are 29-35% inside the regulon, so they will rise. Then the claim stays "raising respiratory
   capacity restores death competence", not "reverting the developmental change".

## ERRb: the riser that had to be tested (author's challenge, 2026-08-18)

`Esrrb` is the **top riser of the 57-gene roster** (+1.578, padj 0.040, **99.6th
expression-matched pct**, up in both genotypes, Myc-suppressed at 6W by -1.04) and was **missing
from the first ledger** -- because the ledger was organised around what could LOWER the chain and
ERRb goes UP, so it needed a repression/competition model nobody had considered. **Coherent, and
reach does NOT disqualify it** (ERRa/b/g share the ERRE, so it would touch all 49 of ERRa's 89).
**Demoted on three tests:** within-animal coupling to OXPHOS **+0.36 (WT) vs -0.29 (Myc+)** --
opposite signs; sibling `Esrrg` also rises (+0.49, 87.8th) and couples **+0.75**; per-animal spread
inside one group is **12-fold** (58/85/101/312/332/696) against `Esrra`'s 1.2-fold. No ESRRB target
set exists in the library. Tracks no state composite (AP_LE +0.13, TEB +0.25, basal +0.09).
**Named test that would settle it: knock ERRb down in the 12W-like state; respiratory subunits
should rise while the rest of the ERRa regulon stays put.**

**RULE: screen RISERS as well as fallers, then discriminate with the within-animal coupling.** That
one test demoted both ERRb and PRC.

**Six discussion panels** (`figures/panels/biogax_*.R`, prefix keeps them outside both runners'
globs -- no slot, no panel count, absent from `analysis_record.qmd`) are collected by
`paper/biogenesis_axis_discussion.qmd` for the collaborator conversation.

Related: [[oxphos-apoptosis-module-verdict]], [[mitocarta-sets-are-membership-loose]],
[[le-is-low-expressing]], [[figure-b-loaders-are-curated-not-top]], [[batch-equals-timepoint]],
[[normalisation-ladder-has-a-floor]], [[mito-content-bounded]].
