# The developmental OXPHOS decline, the PGC1a question, and what licenses the rescue

**2026-08-17.** Written during the Block C cut, in answer to a question raised repeatedly
across Block A and never closed: *is the wild-type 6W->12W OXPHOS decline a PGC1a-axis change?*
The Results describe that decline, the decline carries the anti-apoptotic window, and the paper
then models it experimentally with **PGC1a and NRF1**. A reviewer will ask whether the
developmental change is a PGC1a change at all. It is not, and saying so is what makes the
experiment defensible rather than vulnerable.

Companion: `scripts/47_biogenesis_axis_and_the_developmental_oxphos_decline.R` ->
`results/biogenesis_axis_developmental.rds`. Numbers below carry their source. Values marked
**[47]** are produced by that script; everything else is read from an object already on disk.

---

## 0. Scope, said once

**BATCH = TIMEPOINT** (CLAUDE.md). Every between-age number here is confounded with cohort.
Two things survive that and one does not:

- **survives** -- the **arm-selectivity**. A shared batch effect cannot move one arm of one
  regulon to the 0.05th percentile of a matched null while moving the rest of the same regulon
  to the 99.55th. The comparison is internal.
- **survives** -- the **within-timepoint couplings** (section 4), which remove the timepoint
  means entirely and therefore contain no cohort contrast.
- **does not survive** -- the **magnitude** of the developmental decline.

n = 24 (n = 12 per genotype). This is ranking plus a set of negatives, not confirmatory inference
-- the standing rule since scripts 35/36/37.

---

## 1. Why the question stayed open

Not power. **Both instruments that look like they should answer it are broken for this contrast,
and nobody had measured how badly.**

**(i) The `_MITO` TF lanes are made of OXPHOS genes.** Category 6 promotes, for any significant
TF, the set `TF programme n MitoCarta` (`docs/library_reference/Gray_et_al_developmental_TFS_selection.md`,
lane 2c: 218 sets). A lane whose overlap is OXPHOS-heavy falls because OXPHOS falls. Reading
`TFT_ESRRA_GRAY_AP_LE_MITO` = -2.33 as "ERRa activity fell" is circular.

**(ii) The `_GRAY_<context>` lanes report a cell state, not a factor.** Over the 352 Gray lanes
in this contrast, **adjusted R2 = 0.472 for context and 0.267 for TF identity**. Within `AP_LE`
the 87 mito lanes span an IQR of **0.50 NES**; the spread between context medians is **3.34**.
The demonstration: inside `AP_LE` the most negative "TFs" are `TCF15` -2.74, `SNAPC5` -2.71,
`THAP3` -2.66, `HMGN3` -2.59, `ARNT` -2.59, **`CHCHD3` -2.56 (a MICOS structural protein, not a
transcription factor)**, `GTF3A` -2.53 (general transcription factor IIIA) -- all *below* `ESRRA`
at -2.33. A layer in which a MICOS subunit outranks ERRa is not ranking transcription factors.

**The ambient nobody applied.** Median NES over all 421 TF lanes in this contrast is **-1.461**;
the whole layer is negative-going. Script 24 puts the PGC1a-axis lane median at **-1.144**
(`biogenesis_discrimination.rds$tf_fgsea_summary`, WT_6W->12W) -- **above** ambient. Against its
own layer the biogenesis axis is not preferentially down at all.

**The reading rule that follows** (fixed in `47` before the numbers were read): a `_MITO` lane is
never TF activity; a Gray lane is readable only if it departs from its own context median by more
than that context's IQR, in the same direction as its departure from ambient; a non-Gray lane
carries no context confound and is readable if it departs from ambient by more than the layer
IQR (0.941). The rule passes **0 of 218** mito lanes, 34 of 134 Gray lanes and 19 of 69
context-free lanes **[47]**.

---

## 2. The answer: four lines, and they agree

### 2.1 The regulon splits in two directions -- the decisive test

`results/collapse_module_ownership.rds$core_decomp` (script 44, already on disk):

| part of the PGC1a axis (`CORE_MITO`) | n | WT 6->12W content | matched-null pct |
|---|---|---|---|
| n OXPHOS subunits | 59 | **-0.216** | **0.05** |
| the rest | 241 | **+0.069** | **99.55** |
| whole regulon | 300 | +0.013 | 74.90 |

**A change in a factor's activity acts on its regulon; it cannot act on a functional subset of
it.** The regulon does not move. Its respiratory arm falls at the 0.05th percentile while its
other 241 members rise at the 99.6th. Script 47 PART B generalises this to seven regulons
(`CORE_MITO`, `ESRRA_MITO`, `NRF1_MITO`, `GABPA_MITO`, `E2F1_MITO`, `MYC_MITO`,
`DEVELOPMENTAL_MITO`); if they all split the same way, the split belongs to the OXPHOS genes and
to no factor **[47]**.

`paper/analysis_record.qmd:1906` already carries the sentence -- *"The PGC1a regulon is flat
while only its OXPHOS arm falls"* -- it has simply never been used to answer this question.

### 2.2 The arm profile is "run", not "build"

fGSEA, `timepoint_neg` (`results/fgsea_percategory.rds`), with the content ruler from script 44
where available:

| falls | NES | | does not fall | NES |
|---|---|---|---|---|
| OXPHOS subunits | **-2.670** (padj 1.3e-11) | | OXPHOS **assembly factors** | -0.736 (p 0.94) |
| CI subunits | -2.359 | | CI assembly factors | -0.513 |
| CIV subunits | -2.168 | | CIV assembly factors | **+0.781** |
| CV subunits | -1.979 | | mitoribosome | -0.606 (p 1.00) |
| CIII | -1.586 | | mitochondrial central dogma | **+0.945** |
| chaperones | -1.665 | | mtDNA maintenance | **+1.096** |
| protein import | -1.541 | | mt-tRNA synthetases | +0.833 |

And every catabolic arm **rises**: fatty-acid oxidation **+1.885**, amino acid **+1.835**,
autophagy +1.651, pyruvate +1.632, carbohydrate +1.581, BCAA +1.537, dynamics/surveillance +1.435.

On the content ruler the split is sharper still: OXPHOS subunits **-0.255 at the 0.00 percentile**,
OXPHOS assembly factors **+0.0007 at the 50.25 percentile** -- dead centre of its own null
(`collapse_module_ownership.rds$wt_content`).

**A biogenesis-axis withdrawal predicts subunits AND assembly AND mitoribosome AND import moving
together.** Only the structural subunits move. The corpus already separates the two programmes --
OXPHOS ~188 genes that *run* respiration against biogenesis ~249 genes that *build* the organelle,
17 shared (`docs/2026-07-19_oxphos_axis_biology_and_mtdna_priming.md`, section 3) -- and it is the
"run" half that falls.

*(Reported, and nothing rests on it: the 13 mtDNA-encoded subunits go the other way in the
wild-type timeline, `mt-Co2` +1.68 padj 0.039, `mt-Atp6` +1.76 padj 0.036. That is the axis
CLAUDE.md flags as time-associated and three-way confounded.)*

### 2.3 The factors themselves do not move -- and PGC1a is barely expressed

`results/interaction_results.rds`, wild-type timeline, raw LFC:

| gene | baseMean | WT 6->12W | padj |
|---|---|---|---|
| `Esrra` | 333 | +0.227 | 0.29 |
| `Nrf1` | 489 | -0.249 | 0.35 |
| `Gabpa` | 806 | +0.194 | 0.13 |
| `Tfam` | 277 | +0.097 | 0.81 |
| `Ppargc1b` | 150 | +0.204 | 0.69 |
| **`Ppargc1a`** | **30** | +0.397 | 0.30 |
| **`Pprc1` (PRC)** | **2102** | **-0.471** | **1.9e-4** |

Two facts the paper needs. **`Ppargc1a` sits at baseMean 30 -- effectively absent from purified
MEC.** And the only PGC-1 family member expressed at a real level is **PRC (`Pprc1`), 70x higher,
and it is the only one that moves** -- -0.471 (padj 1.9e-4) in the wild type, -0.585 (padj 1.1e-6)
in the Myc+ gland, interaction padj 1. A genotype-independent developmental decline of the
growth-coupled coactivator. **`Pprc1` appears nowhere else in the corpus.**

### 2.4 The generic regulons of the same factors go UP

Free of the Gray context artifact: `TFT_GABPA_CHIPATLAS` **+1.20**, `TFT_ESRRA_CHIPATLAS`
**+1.19**, `TFT_NRF1_CHIPATLAS` +0.98, `TFT_NRF1_DOROTHEA_ABC` +0.95. And in Category 7 --
the sets built to discriminate exactly this -- `CORE_MITO` **+0.85 (p 0.90)**, `ESRRA_MITO` +0.80,
`GABPA_MITO` +0.79, `NRF1_MITO` -0.98, `E2F1_MITO` +0.93 (the proliferation-confounder control),
`DEVELOPMENTAL_MITO` +1.27. **Not one of them is down.**

---

## 3. The hypothesis ledger

**H0 -- the PGC1a/ERRa/NRF1/GABP axis declines. REJECTED** on all four lines above.

**H1 (primary) -- the respiratory chain is cargo of a cell-state switch: TEB / lineage-suppressed
AP -> mature duct.** The TF changes that underlie it are developmental, not mitochondrial: down
`HMGA2 ARNTL2 ELK3 GRHL3 FOXQ1 IRF6 FOXM1 TFDP1 CENPA` (all `_HS_TEB` / `_HS_LE`, -2.1 to -2.4),
up `MKX OSR1 OSR2 HEY2 ATOH8` (all `_HS_DUCT` -- **the only positive context in the whole layer**,
median +1.27 against -1.46 ambient). Content ruler, script 44 `le_content`:
`MG_TEB_VS_DUCTAL_HS_GRAY_UP` **-0.422 at the 0.0 percentile** -- a larger fall than OXPHOS's own
-0.255 -- and `MG_HEVSLE_AP_GRAY_DN` -0.294 (pct 0.0), still **-0.238 (pct 0.00) with every
MitoCarta gene removed**, so the carrier state falls independently of its mitochondrial content.

**H2 (runner-up) -- a coactivator-supply change: PRC, not PGC1a.** The one version of "the
biogenesis axis changed" the data leave standing, because it changes coactivator *supply* on the
NRF1/GABP respiratory promoters without changing factor abundance or regulon-wide output.
**Already half-falsified:** within timepoint, `Pprc1`'s per-sample coupling to the OXPHOS
composite is ambient. If it acts, it acts as a switch, not a rheostat **[47]**.

**H3 -- the reciprocal arm: FOXO3/SIRT1 quiescence.** `TFT_FOXO3_CHUNG` **+1.43** (p 0.059,
padj 0.10) is **the single highest of all 69 context-free TF lanes**, with `TFT_FOXO1_CHUNG`
-0.63 as the specificity control. `Foxo3` **+0.473 (padj 0.0078)**, `Sirt1` +0.560 (padj 0.019),
`Bnip3` +0.635 (padj 0.021); autophagy +1.651; dynamics and surveillance +1.435. Within
timepoint, `Sirt1` and `Bnip3` are among the strongest *inverse* correlates of the OXPHOS
composite in the whole transcriptome **[47]**.

> **The caveat that must travel with H3, every time.** The `Foxo3` **interaction padj = 1**. The
> supported statement is *"the wild type rises (padj 0.0078) and the Myc+ gland does not (-0.013,
> ns)"*. It is **not** *"Myc blocks the FOXO3 rise"*. Two halves that differ in significance do
> not differ significantly -- the same trap already recorded for the Foxo3/PUMA pair.

---

## 4. The tension H1 has to carry honestly

Within a timepoint the OXPHOS composite tracks the cell-cycle programme at **r = +0.887 in the
wild type and +0.842 in the Myc+ gland** (87-gene composite, timepoint means removed, so no
cohort contrast remains). Yet *between* the ages OXPHOS falls **-0.255 against proliferation's
-0.047** -- five times more.

Both are true, and together they say something neither says alone: the decline is **not
proportional growth withdrawal** and **not an axis-wide biogenesis change**. It is a state switch,
and a state switch is not a dial. This protects the record's load-bearing negative
(*"the gland de-respires without de-proliferating"*, `analysis_record.qmd:1901`), which stays true
on the between-age ruler and now gains its within-age companion.

**The limit that stops H1 being a result.** Carrier and cargo correlate at **r = 0.93**
(`collapse_module_ownership.rds$le_within`, 12 wild-type samples, mito genes removed from the
carrier). Script 44's adjustment takes the OXPHOS time-beta from -0.459 to -0.151 -- but at that
collinearity the adjusted beta is **a bound on how much of the decline could be carried, not an
estimate of how much is**. Bulk cannot separate them. **Sorted AP/HS/BA or single-cell is the
design that settles it**, and that is the honest answer to a reviewer who asks.

---

## 5. What licenses reverting the decline with PGC1a and NRF1

If the axis did not cause the decline, the experimental part needs its justification in numbers.
Five, and the first is the strongest because it turns the mouse's limitation into the
experiment's licence.

### H1. The mouse cannot de-confound; the intervention can

In the gland, respiratory-chain expression is inseparable from the state that carries it
(r = 0.93 to the carrier; r = 0.887 to the cell-cycle programme within an age). **No observational
contrast in this animal can distinguish "respiration fell" from "the state that respires became
rarer"** -- that is precisely why the question stayed open for a whole analysis block. An
intervention that moves respiration *while cell state is held fixed* is the only instrument that
breaks the confound. PGC1a is not chosen because it reproduces the developmental cause. It is
chosen because it does something the gland cannot do.

### H2. Reach, not authorship

PGC1a did not lower these genes, but it reaches them.

| regulon | n | covers the declining subunits |
|---|---|---|
| `CORE_MITO` (ERRa u NRF1 u GABP) | 300 | **59 of 89 (66.3%)** |
| `ESRRA_MITO` | 177 | 49 of 89 (55.1%) |
| `NRF1_MITO` | 81 | 19 of 89 (21.3%) |
| `GABPA_MITO` | 81 | 10 of 89 (11.2%) |

Reverting a deficit requires **reach over the affected genes**, not responsibility for the deficit.
This is the ordinary logic of a rescue and it should be stated as such.

### H3. Two factors are a control, not a redundancy

`ESRRA_MITO` and `NRF1_MITO` overlap in only **22 genes** -- and that intersection is where the
respiratory subunits concentrate:

| part | n | OXPHOS-subunit fraction |
|---|---|---|
| ESRRA not NRF1 | 155 | 24.5% |
| NRF1 not ESRRA | 59 | 13.6% |
| **shared** | **22** | **50.0%** |

Against the union of the two regulons the shared part is **2.07x enriched for respiratory
subunits (11 observed vs 5.31 expected, hypergeometric p = 0.0052)**. Two interventions with
largely non-overlapping reach that produce the same phenotype **localise the effect to their
intersection**, and their intersection is the respiratory chain. This converts "PGC1a raises many
things" from the reviewer's objection into the design's specificity control.

### H4. PGC1a cannot restore death by transcribing the death machinery

`collapse_module_ownership.rds$ownership`, hypergeometric against the MitoCarta universe:

| regulon | pro-apoptotic overlap | expected | fold | p_enrich | genes |
|---|---|---|---|---|---|
| `CORE_MITO` | 5 of 25 | 6.58 | **0.76** | 0.83 | Aifm1, Aifm2, **Cycs**, Endog, Ifi27 |
| `MYC_MITO` | 7 of 25 | 2.61 | **2.68** | **0.011** | Aifm3, **Bax**, Casp8, Casp9, Diablo, **Pmaip1**, Sphk2 |

**The PGC1a regulon is if anything depleted for the death effectors -- no PUMA, no BAX, no
caspase.** So if PGC1a overexpression restores PUMA-dependent death, it is not doing it by
directly transcribing the machinery; it has to run through the respiratory state, which is the
claim under test. Had the effectors been inside the regulon, the experiment would have been
circular. **They are inside MYC's regulon instead** -- which is the paper's own result
(`analysis_record.qmd`, "no co-regulated OXPHOS/apoptosis module").

The single exception is worth a sentence rather than a footnote: the one death gene PGC1a
reaches is **cytochrome c**, which is simultaneously a respiratory carrier -- the paper's
asset-and-liability duality inside one gene.

### H5. The mouse sets the rescue's target effect size

The deficit is a defined size at OXPHOS-subunit level (content ruler, script 44):

| context | deficit (log2) | remaining | lost | restoration needed |
|---|---|---|---|---|
| wild type 6W->12W | -0.255 | 0.838 | **-16.2%** | **x1.19** |
| Myc+ 6W->12W | -0.405 | 0.755 | **-24.5%** | **x1.32** |

Because `Ppargc1a` is at the floor of expression in MEC, **the transgene level is not a
calibration -- the output is**. Report where the intervention lands against x1.19-x1.32 at
subunit level. That is what separates *"we reverted the developmental change"* from *"we pushed
respiration supraphysiologically"*, and either is publishable as long as it is the one claimed.

### H6. Pre-specify the overshoot, because it will happen

The developmental change was arm-selective. **The rescue will not be** -- the arms that never fell
are still inside the regulon:

| arm that did NOT fall | its WT content (pct of null) | fraction inside `CORE_MITO` |
|---|---|---|
| OXPHOS assembly factors | +0.0007 (50.25) | 31.3% |
| mitoribosome | flat | 34.9% |
| mitochondrial central dogma | rises | 29.1% |

So predict, in writing, that PGC1a overexpression raises them. If death competence returns
alongside a broader mitochondrial expansion, the supported sentence is **"raising respiratory
capacity restores death competence"**, not "reverting the developmental change restores it".
Naming the overshoot before the experiment is what keeps the claim intact after it.

---

## 6. The sentences for the paper

**On the developmental change (Results, section 2):**

> The maturing gland withdraws from the respiratory chain specifically: the structural subunits of
> complexes I, III, IV and V fall, while the assembly factors that build the same complexes, the
> mitochondrial ribosome and the mitochondrial central dogma do not, and the catabolic arms rise.
> The withdrawal is not a change in mitochondrial biogenesis programme: the ERRa/NRF1/GABP regulon
> as a whole is unchanged across the two ages, and only its respiratory subset moves.

**On the intervention (Results or Discussion, wherever PGC1a first appears):**

> Because respiratory gene expression in the gland is inseparable from the epithelial state that
> carries it, the developmental data cannot establish whether respiration itself is the operative
> variable. We therefore raised it directly. PGC1a and NRF1 were chosen not to reproduce the
> developmental mechanism -- which is not an axis-activity change -- but because their regulons
> reach two thirds of the subunits that declined while excluding the apoptotic effectors, so a
> restored death phenotype cannot be direct transcription of the death machinery.

**On what the two factors buy (Methods or legend):**

> The ERRa and NRF1 target sets overlap in 22 genes, of which half are respiratory-chain subunits
> (2.1-fold enrichment over their union, p = 0.005); convergence of the two interventions therefore
> localises the effect to the respiratory arm rather than to either regulon as a whole.

---

## 7. Sentence risks this creates (for the `PANELS.md` list)

1. Any sentence implying the developmental OXPHOS loss is a **decline in PGC1a-driven biogenesis**
   -- reword. The regulon sits at the 74.9th percentile of its null while its OXPHOS arm sits at
   the 0.05th.
2. *"The gland de-respires without de-proliferating"* -- true **between the ages** (-0.255 against
   -0.047); within an age the two covary at **r = 0.89**. Keep the sentence, add the scope word.
3. Any sentence attributing the decline to a **named mitochondrial TF read off a `_GRAY_` or
   `_MITO` lane** -- those lanes report cell-state context and gene content. A MICOS subunit and a
   general transcription factor outrank ERRa in the same context.
4. Any sentence pairing the **FOXO3** rise with Myc **blocking** it -- interaction padj = 1.

---

## 8. What this hands to the perturbation

- **Pre-specified primary readout:** the 22-gene ERRa n NRF1 intersection (50% respiratory
  subunits), against the 155- and 59-gene exclusive parts as the specificity contrasts.
- **Pre-specified calibration:** x1.19 (wild-type deficit) to x1.32 (Myc+ deficit) at
  OXPHOS-subunit level.
- **Pre-specified overshoot check:** assembly factors, mitoribosome, central dogma -- flat in the
  gland, inside the regulon, expected to rise.
- **Pre-specified negative:** the death effectors are not in the PGC1a regulon; `Bbc3`, `Bax`,
  `Casp9`, `Pmaip1` should not be direct PGC1a targets in the perturbation transcriptome. If they
  are, the rescue is circular and the interpretation changes.
- **The question bulk cannot answer:** whether the state or the respiration carries the
  developmental decline (r = 0.93). Sorted AP/HS/BA or single-cell.
- **The open runner-up:** PRC (`Pprc1`), the only PGC-1 family coactivator expressed in MEC, down
  in both genotypes. Never examined in this project; the obvious comparator arm if the PGC1a
  result raises the question of which coactivator the gland actually uses.

---

## 9. Relevant sets and analyses -- the inventory

**Gene sets.** Category 7 `07_biogenesis_discrimination_mouse.gmt`, built for exactly this
question: `CORE_MITO` (ERRa u NRF1 u GABP targets n MitoCarta), `ESRRA_MITO`, `NRF1_MITO`,
`GABPA_MITO`, `E2F1_MITO` (proliferation-confounder control), `DEVELOPMENTAL_MITO` (ER),
`MYC_MITO`, `MYC_SPECIFIC_MITO`. Category 1 MitoCarta arms for the run/build split. Category 6 TF
lanes -- **read only under the rule in section 1**. `MG_HEVSLE_*` and `MG_TEB_VS_DUCTAL_*` for the
carrier. `gray_chea_mito_tf_shortlist.csv` for the TF-to-mito overlap provenance.

**Analyses already on disk.** Script 20 (`fgsea_percategory.rds`) -- the NES layer, the only object
carrying all five rankings. Script 24 (`biogenesis_discrimination.rds`) -- Category 7 lanes per
contrast, and the TF lane medians. Script 26 -- the developmental trajectory. Script 32 -- the
content claim. Script 43 -- substrate specificity. **Script 44 (`collapse_module_ownership.rds`) --
the one that already contains the answer**: `core_decomp`, `wt_content`, `le_content`, `le_within`,
`ownership`. Script 45 -- the state readings. Script 46 -- the ruler test behind "priority, not
capacity".

**New in script 47.** The arm map with matched nulls across all MitoCarta arms; the within-regulon
split generalised to seven regulons; the TF-layer decomposition and reading rule; the factor roster
against an expression-matched null; the within-timepoint couplings; the carrier's identifiability
limit; and the rescue pre-specification of section 5.
