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
The demonstration: inside `AP_LE`, ranked by NES over 87 lanes, `ESRRA` comes **21st** and `GABPA`
**25th**. Above them: `TCF15` -2.74, `SNAPC5` -2.71, `THAP3` -2.66, `HMGN3` -2.59, `ARNT` -2.59,
**`CHCHD3` -2.56 (rank 6 -- a MICOS structural protein, not a transcription factor)**, `POU5F1B`,
`NR2C2`, **`GTF3A` -2.53 (rank 9 -- general transcription factor IIIA)**, `ZNF232`, `ZKSCAN8`,
`ZNF706`. **A layer in which a MICOS subunit and a general transcription factor both outrank ERRa
is not ranking transcription factors.**

**The ambient nobody applied.** Median NES over all 421 TF lanes in this contrast is **-1.461**;
the whole layer is negative-going. Script 24 puts the PGC1a-axis lane median at **-1.144**
(`biogenesis_discrimination.rds$tf_fgsea_summary`, WT_6W->12W) -- **above** ambient. Against its
own layer the biogenesis axis is not preferentially down at all.

**The reading rule that follows** (fixed in `47` before the numbers were read): a `_MITO` lane is
never TF activity; a Gray lane is readable only if it departs from its own context median by more
than that context's IQR, in the same direction as its departure from ambient; a non-Gray lane
carries no context confound and is readable if it departs from ambient by more than the layer IQR
(0.941). It passes **0 of 218 mito lanes**, 34 of 134 Gray lanes (25%) and 19 of 69 context-free
lanes (28%) -- it discriminates rather than waving everything through.

**One caveat the rule surfaces, and it should be stated rather than buried.** `TFT_ESRRA_GRAY_AP_LE`
(-2.32) *does* clear its own context median by 0.54, so ERRa's alveolar-progenitor programme is a
genuine within-context outlier. It is still not evidence of an ERRa activity change: that programme
is precisely where the ERRa-mito overlap was detected in the first place (47 mito genes,
p = 1.1e-22, the promotion that created `ESRRA_MITO`), so it fails the *content* filter even while
passing the *context* one. Both filters have to pass, and the generic ERRa regulon -- which fails
neither -- goes **up** (+1.19).

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
other 241 members rise at the 99.6th.

**Script 47 PART B generalises it to seven regulons, and the four that contain the respiratory
chain all split the same way:**

| regulon | n OXPHOS | n rest | **pct OXPHOS** | **pct rest** | pct whole | verdict |
|---|---|---|---|---|---|---|
| `CORE_MITO` | 59 | 241 | **0.0** | **99.2** | 73.4 | arm-selective |
| `ESRRA_MITO` | 49 | 128 | **0.4** | **97.6** | 54.5 | arm-selective |
| `NRF1_MITO` | 19 | 62 | **1.0** | **84.8** | 40.8 | arm-selective |
| `GABPA_MITO` | 10 | 71 | **4.1** | **79.4** | 59.2 | arm-selective |
| `MYC_MITO` | 9 | 109 | 24.3 | 91.0 | 87.6 | inconclusive |
| `DEVELOPMENTAL_MITO` | 9 | 84 | 30.8 | 97.6 | **96.0** | **moves as a unit (UP)** |
| `E2F1_MITO` | 2 | 42 | -- | 80.6 | 74.2 | not testable |

**Read it as four of four, not four of seven.** The three that fail the split are exactly the
three with almost no respiratory content -- `E2F1_MITO` holds 2 OXPHOS subunits (below the
3-gene floor, so the test cannot run at all), `MYC_MITO` and `DEVELOPMENTAL_MITO` hold 9 each.
A regulon cannot show an arm-selective split in an arm it barely contains.

**One regulon does move as a unit, and it is the ER one:** `DEVELOPMENTAL_MITO` sits at the
**96.0th percentile** of its null as a whole -- the ER-reachable mitochondrial genes *rise* in the
maturing gland. That is the shape an axis change makes, and it is the shape the PGC1a axis does
not make.

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

**On the content ruler with a matched null, by arm class** (script 47 PART A, all 95 arms):

| class | n arms | median content | **median null pct** |
|---|---|---|---|
| **chain_structural** | 13 | **-0.159** | **1.75** |
| chain_assembly | 3 | +0.0024 | **53.2** |
| organelle_build | 17 | +0.036 | 67.8 |
| metabolic | 25 | +0.089 | 75.6 |
| other | 35 | +0.046 | 78.0 |
| *chain_mtDNA (confounded, see below)* | 2 | +0.633 | 100 |

**The structural chain sits at the 1.75th percentile as a class; the assembly factors that build
those same complexes sit at the 53rd.** Individually: `CIV_SUBUNITS` **-0.411 (pct 0.5)** is the
most negative arm in the compartment, then OXPHOS subunits -0.255 (pct 0.0), CI -0.253 (0.05),
CIII -0.251 (2.9), CV -0.216 (1.75). Against `OXPHOS_ASSEMBLY_FACTORS` +0.0007 (**pct 50.9**),
`CI_ASSEMBLY` +0.0024 (53.2), `CIV_ASSEMBLY` +0.018 (55.2).

**Two details worth keeping.** `COMPLEX_II` **does not fall** (+0.021, pct 57.5) -- succinate
dehydrogenase is the one respiratory complex that neither pumps protons nor contains an
mtDNA-encoded subunit, and it is also a TCA enzyme. And `ELECTRON_CARRIERS` **rises** (+0.155,
pct 94.6). The withdrawal is from the *proton-pumping, mtDNA-containing* complexes specifically.

**Honest qualification.** The build arms do not follow the chain down, but two of them are low
rather than central: chaperones **pct 9.5** and protein import **pct 8.65**, with the mitoribosome
at **pct 26.0** and the central dogma at **86**. So the accurate sentence is *"the assembly
factors do not move at all and the build arms do not follow"*, not *"nothing else moves"*.

**A biogenesis-axis withdrawal predicts subunits AND assembly AND mitoribosome AND import moving
together.** They do not. The corpus already separates the two programmes -- OXPHOS ~188 genes that
*run* respiration against biogenesis ~249 genes that *build* the organelle, 17 shared
(`docs/2026-07-19_oxphos_axis_biology_and_mtdna_priming.md`, section 3) -- and it is the "run" half
that falls, with the assembly factors as the sharp internal control.

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
in the Myc+ gland, interaction padj 1. **`Pprc1` appears nowhere else in the corpus.**

> **But it does not clear its own null, and that matters.** Against expression-matched genes PRC's
> wild-type fall sits at the **7.1st percentile** -- inside the [5, 95] bar fixed before the test.
> Its tiny padj is **precision, not effect size**: a gene at baseMean 2102 is measured well enough
> for a modest fold change to clear FDR. So PRC is a *lead worth an antibody*, not a result, and
> the same test disposes of the other candidate fallers: `Mybl2` 5.6th percentile, `Nrip1` 9.2nd,
> `Ppard` 9.2nd. **Nothing in the roster clears the bar on the down side.**
>
> What does clear it, all on the up side: **`Foxo3` 95.7th** (padj 0.0078), **`Sirt1` 97.9th**
> (0.019), `Srebf1` 98.1st, `Elf5` 98.3rd, `Esrrb` 99.6th. The developmental signal in this roster
> is a *gain* of the quiescence arm, not a loss of the biogenesis arm.

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

**H2 (demoted after 47) -- a coactivator-supply change: PRC, not PGC1a.** This was the one version
of "the biogenesis axis changed" that looked survivable: PRC changes coactivator *supply* on the
NRF1/GABP respiratory promoters without changing factor abundance or regulon-wide output.
**It does not survive its own null.** PRC's fall is at the **7.1st expression-matched percentile**
(inside the pre-set bar) and its within-timepoint coupling to the OXPHOS composite is ambient.
Keep it as **a lead for an antibody, not a mechanism** -- and keep it mainly because PRC, not
PGC1a, is what this tissue actually expresses.

**H3 (promoted after 47) -- the reciprocal arm: FOXO3/SIRT1 quiescence.** This is now the
strongest positive signal in the analysis, and it is the only one that clears every bar set for it.
`TFT_FOXO3_CHUNG` **+1.43** is **the single highest of all 69 context-free TF lanes** and is
readable under the lane rule, with `TFT_FOXO1_CHUNG` -0.63 as the specificity control.
`Foxo3` **+0.473, padj 0.0078, 95.7th expression-matched percentile**; `Sirt1` +0.560, padj 0.019,
**97.9th**; `Bnip3` +0.635, padj 0.021, 93.0th. Autophagy +1.651; dynamics and surveillance +1.435.
And within timepoint, against all 16,922 expressed genes, **`Sirt1` (3.8th percentile) and `Bnip3`
(5.7th) are among the strongest inverse correlates of the OXPHOS composite in the transcriptome**
**[47]**.

> **The caveat that must travel with H3, every time.** The `Foxo3` **interaction padj = 1**. The
> supported statement is *"the wild type rises (padj 0.0078) and the Myc+ gland does not (-0.013,
> ns)"*. It is **not** *"Myc blocks the FOXO3 rise"*. Two halves that differ in significance do
> not differ significantly -- the same trap already recorded for the Foxo3/PUMA pair.

---

## 3A. The FOXO3 story, in full -- and whether it is the strongest hypothesis

**Short answer: FOXO3 is the strongest readable SIGNAL in this contrast, and it is not the
strongest hypothesis for the OXPHOS decline. The distinction is the whole of the section.**

### Why it is the strongest readable signal

Every other candidate in the TF layer dies on one of two filters. FOXO3 passes both.

- **Content.** `TFT_FOXO3_CHUNG` is 39 genes and **contains no OXPHOS subunit at all** (0 of 89).
  Its six MitoCarta members are `Bcl2l11`, `Bnip3l`, `Cat`, `Pink1`, `Prdx3`, `Sod2` -- antioxidant,
  mitophagy and apoptosis, not the respiratory chain. So its rise cannot be mitochondrial content.
- **Context.** It is a Chung lane, not a Gray lane, so it carries no cell-state context and cannot
  be the artifact that destroys the mammary-context layer. It is **the single highest of all 69
  context-free lanes** at **+1.43** (p 0.059, padj 0.10).
- **Gene level agrees, against an expression-matched null:** `Foxo3` +0.473 (padj 0.0078,
  **95.7th pct**), `Sirt1` +0.560 (0.019, **97.9th**), `Bnip3` +0.635 (0.021, 93.0th).
  `FOXO1` is null in every contrast (-0.63 in the wild-type timeline) -- the specificity control.

### The lane across all four contrasts, which reframes it

| ranking | NES | padj |
|---|---|---|
| `myc_6W` | **-1.76** | **0.0030** |
| `myc_12W` | **-1.80** | **0.0031** |
| `6>12W_wt` | **+1.43** | 0.103 |
| `6>12W_myc` | **+1.71** | **0.0061** |
| interaction | +1.09 | 0.37 |

Two readings, and the second is the one worth keeping:

1. **The programme rises with age in BOTH genotypes** (+1.43 wild type, +1.71 Myc+), and the
   interaction is not significant. So **"Myc blocks the developmental FOXO3 rise" is false.**
2. **Myc suppresses the programme at both ages, by about the same amount** (-1.76, -1.80), with
   FOXO1 null throughout. That is a clean, stable genotype effect on the PUMA regulator's
   programme -- and it is a genotype contrast, which this design measures well.

> **Do not read the two genotype rows as "the FOXO3 programme does not attenuate".** NES is
> scale-free and cannot see an amplitude change -- that is the four-list result recorded under
> Fig. 1H. Equal NES at both ages is not equal effect size.

### Which arm of the programme moves -- named before the values were read

FOXO3 has an arrest arm, an atrophy/turnover arm, an antioxidant arm and an apoptotic arm. They do
not have to move together, and here they do not:

| arm | moves? | evidence |
|---|---|---|
| **atrophy / turnover** | **yes** | `Fbxo32` **+1.156 (padj 5.4e-10)**, `Bnip3` +0.635 (0.021), `Sirt1` +0.560 (0.019), `Pink1` +0.370 |
| antioxidant | marginal | `Cat` +0.384 (0.068); `Sod2`, `Prdx3`, `Txnip`, `Sesn1` flat |
| cell-cycle arrest | **no** | `Cdkn1a` -0.213, `Cdkn1b` +0.028, `Gadd45a` -0.034, `Ccng2` -0.097, all ns |
| **apoptotic** | **no** | `Bcl2l11` +0.255 (ns), **`Bbc3` +0.062 (padj 0.881)** |

Leading edge, 15 of 38: `Fbxo32, Fasl, Sirt1, Cat, Pink1, Bnip3l, Rbl2, Cited2, Bcl2l11, Pik3ca,
Abcb1a, Klf4, Ar, Hbp1, Dusp5` -- broad enough that this is not one gene carrying a small set.

**The consequence for the death story: `Bbc3` is FLAT in the wild-type timeline.** The attractive
sentence -- *the adult gland raises FOXO3, which arms PUMA* -- **is not supported**. PUMA falls in
the *Myc+* timeline (-0.479, padj 0.016), which is a different claim, and the interaction is ns
(padj 0.844).

### Why it is not the driver: reach

The same test that disqualified the PGC1a axis as the *cause* disqualifies FOXO3, and it has to be
applied consistently or it is not a test:

| regulon | reaches, of the 89 declining subunits |
|---|---|
| `CORE_MITO` (ERRa/NRF1/GABP) | **59** |
| `ESRRA_MITO` | 49 |
| **`TFT_FOXO3_CHUNG`** | **0** |

**A factor that touches none of the affected genes cannot be their proximal transcriptional
cause**, however clean its own signal is. Whatever lowered the respiratory subunits, FOXO3 did not
do it directly.

### And it is not separable from the state

Within timepoint, on the non-mitochondrial part of the programme (the 6 MitoCarta members stripped
first, or the coupling would be partly self-correlation):

| | vs OXPHOS composite | percentile among all genes | vs the AP_LE carrier |
|---|---|---|---|
| wild type | -0.691 | **11.5th** | -0.533 |
| Myc+ | -0.521 | 15.2nd | -0.697 |

A mid-range percentile, and an anti-correlation with the carrier state as strong as the one with
OXPHOS. **FOXO3 is not an independent axis; it is another reading of the adult state.**

### The verdict

**FOXO3 is the cleanest TF-level READOUT of the state in which respiration is lower.** It
corroborates H1 rather than competing with it: the state hypothesis predicts exactly this -- a
quiescence/turnover programme up, respiration down, no direct regulatory link between them.

What FOXO3 genuinely adds to the paper is **not** a mechanism for the OXPHOS decline. It is:

1. **A stable, Myc-suppressed programme on the PUMA regulator** (padj 0.003 at both ages, FOXO1
   null). That is a genotype result, which this design measures well, and it sits directly beside
   the death story.
2. **The one place a TF statement is possible at all in this dataset**, which is worth saying in
   the methods rather than leaving a reader to wonder why no factor is named.
3. **A named perturbation target** for the cell system: if the adult state's turnover arm is what
   accompanies low respiration, FOXO3 is how you test that -- not by asking whether it lowers
   OXPHOS transcripts, which it cannot, but by asking whether it changes death competence at fixed
   respiratory state.

---

## 4. The tension H1 has to carry honestly

Within a timepoint the OXPHOS composite tracks the cell-cycle programme at **r = +0.887** in the
wild type and +0.842 in the Myc+ gland (87-gene composite, timepoint means removed, so no cohort
contrast remains). Yet *between* the ages OXPHOS falls **-0.255 against proliferation's -0.047** --
five times more.

> **The set-level r must not be quoted on its own, and this is the corpus's oldest trap.** In the
> same 12 samples the OXPHOS composite correlates with `CORE_MITO` at **0.980** and with the
> **mitoribosome at 0.967** -- both *higher* than with proliferation. Everything correlates with
> everything; a set-level 0.89 is what the window hands you. **Only the gene-level percentile
> against all 16,922 expressed genes is readable**, and there the ordering is informative:
> `Tfdp1` **99.8th**, `Mterf4` 99.5th, `Gata3` 95.5th, `Mybl2` **94.2nd**, `Foxm1` 88.1st at the
> top; `Sirt1` **3.8th** and `Bnip3` **5.7th** at the bottom; and the biogenesis factors in the
> middle (`Gabpa` 79th, `Tfam` 65th, `Nrf1` 45th, `Pprc1` 44th, `Esrra` 43rd). **The growth
> machinery and the quiescence arm are the two poles; the biogenesis axis is not at either.**

Taken together the two rulers say something neither says alone: the decline is **not proportional
growth withdrawal** and **not an axis-wide biogenesis change**. It is a state switch, and a state
switch is not a dial. This protects the record's load-bearing negative (*"the gland de-respires
without de-proliferating"*, `analysis_record.qmd:1901`), which stays true on the between-age ruler
and now has its within-age companion -- stated as the percentile, not the raw r.

**The limit that stops H1 being a result.** Carrier and cargo correlate at **r = 0.930** with every
mitochondrial gene stripped from the carrier (12 wild-type samples), giving a **variance inflation
factor of 7.4**. Script 44's adjustment takes the OXPHOS time-beta from -0.459 to -0.151 -- but at
that collinearity the adjusted beta is **a bound on how much of the decline could be carried, not
an estimate of how much is**. The specificity control behaves as it should (`AP_HE` non-mito
**-0.129**, i.e. no association), and the TEB composite is a looser carrier (**+0.608**), which is
the more honest anchor of the two. Bulk cannot separate carrier from cargo. **Sorted AP/HS/BA or
single-cell is the design that settles it**, and that is the answer to a reviewer who asks.

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
| `ESRRA_MITO` | 4 | 3.88 | 1.03 | 0.56 | Aifm2, **Cycs**, Endog, Ifi27 |
| `NRF1_MITO` | 1 | 1.78 | **0.56** | 0.85 | Aifm1 |
| `MYC_MITO` | 7 of 25 | 2.61 | **2.68** | **0.011** | Aifm3, **Bax**, Casp8, Casp9, Diablo, **Pmaip1**, Sphk2 |

**The PGC1a regulon sits at or below the base rate and contains no BCL-2 family member, no BAX and
no caspase** -- the four or five genes it does hold are the release/effector-adjacent set
(Aifm1/2, Endog, Ifi27) plus cytochrome c. So if PGC1a overexpression restores PUMA-dependent
death, it is not doing it by directly transcribing the machinery; it has to run through the
respiratory state, which is the claim under test. Had the effectors been inside the regulon, the
experiment would have been circular. **They are inside MYC's regulon instead** -- which is the
paper's own result (`analysis_record.qmd`, "no co-regulated OXPHOS/apoptosis module").

*(State the fold honestly: `CORE_MITO` at 0.76 and `NRF1_MITO` at 0.56 are below expectation but
neither is significantly depleted -- with 25 pro-apoptotic genes in the MitoCarta universe the
test has no power to prove depletion. The claim that carries is the **specific absence of PUMA,
BAX and the caspases**, which is a membership fact and needs no p-value.)*

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
   -0.047, and the five-fold gap is the point). Add the scope word. Within an age the two do
   covary, but **do not quote the raw r (0.89)**: in the same samples OXPHOS correlates with
   `CORE_MITO` at 0.98 and the mitoribosome at 0.97, so a set-level 0.89 is the ceiling, not a
   finding. If the within-age point is needed, quote the **gene-level percentile** -- `Tfdp1`
   99.8th, `Mybl2` 94.2nd at the top, `Sirt1` 3.8th, `Bnip3` 5.7th at the bottom.
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
- **The open lead, at its true weight:** PRC (`Pprc1`) is the only PGC-1 family coactivator
  expressed in MEC (baseMean 2102 against PGC1a's 30) and falls in both genotypes (padj 1.9e-4) --
  but it does **not** clear its expression-matched null (7.1st percentile) and its per-sample
  coupling is ambient. It is a reason to check a protein, not a mechanism to test. Its real weight
  in the paper is the *expression* fact: if a reviewer asks which coactivator this tissue uses,
  the answer is not PGC1a.

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
