---
date: 2026-11-10
tags: [project/myc_mouse, nat-metab, reviewer-prep, myc-dose, attenuation, figure-1b, extended-data-1]
status: >
  Reviewer preparation, assembled from the existing record at the author's request
  (collaborator comments, 2026-11-10). NO NEW ANALYSIS: every number is read from a
  committed object, a dated note or the author's Prism file, and each is sourced below.
  Three numbers are new READ-OUTS of existing files rather than new results, and are
  flagged as such. Two items need the author to settle them (section 6).
purpose: >
  The question a reviewer is expected to ask -- "why does MYC protein fall between 6
  and 12 weeks?" -- and what the corpus can and cannot answer. Written so the reply
  can be drafted without re-deriving anything.
relates-to:
  - docs/2026-07-25_attenuation_competing_hypotheses.md   (the four hypotheses and the 2026-07-25 scans; its banner is now partly superseded, see section 5)
  - outputs/myc_stability_panel/README.md                 (script 41: the turnover machinery, generated from the run)
  - docs/2026-07-25_death_narrative_dose_vs_competence.md (the dose-vs-competence frame after the blots)
  - paper/analysis_record.qmd                             ("The mechanism: dose", the long-form reasoning)
  - docs/2026-09-21_two_timeline_verification.md          (section 8: the four rescaling rates)
  - figures/natmetab/edfig1a_myc_network.R               (ED Fig. 1a, the network negative control)
  - results/lesner_mechanism_test.rds                     (script 39: autophagy and mitophagy)
  - results/interaction_results.rds                       (DESeq2, raw LFC + IHW-adjusted P)
---

# Why does MYC protein fall 6W -> 12W? What the record answers, and what it does not

## 0. The distinction the whole reply turns on

Two different declines are easy to conflate, and the manuscript explains only one of them.

| | what falls | is it explained? |
|---|---|---|
| **The programme** | the amplitude of MYC's transcriptional output, by about half | **Yes.** It tracks the protein, and MYC-ER shows the biology does not depend on it |
| **The protein** | MYC protein itself, at constant message | **No.** We measured it and excluded three routes; no mechanism is established |

The *Myc* **transcript does not fall at all** -- this is the control the dose argument needs, and
it holds. A reviewer asking "why does MYC protein fall" is asking about row 2, where the corpus
has negatives and a hypothesis, not an answer. Sections 1-4 are what we can say; section 5 is the
proposal; section 7 is how to reply.

## 1. The three anchored observations

**1.1 The message is flat and the genotype gap widens** (DESeq2, `interaction_results.rds`, raw
LFC with IHW-adjusted *P*, n = 6 per group):

| contrast | *Myc* |
|---|---|
| MYC+ vs WT, 6W | **+1.68** [+1.33, +2.03], padj 2.7e-17 |
| MYC+ vs WT, 12W | **+1.82** [+1.47, +2.17], padj 2.1e-20 |
| 12W vs 6W within MYC+ | −0.10 [−0.45, +0.25], padj **0.76** |
| 12W vs 6W within WT | −0.25 [−0.60, +0.11], padj 0.48 |

The author's own Prism sheet `1b Myc mRNA` agrees independently: MYC+ 8,140 -> 7,580 normalised
counts, Tukey adjP **0.91**; WT 2,540 -> 2,140, adjP 0.97. Measured as retention across the
transcriptome, `Myc` is the **most retained gene of 8,774** (1.085, 99.99th percentile).
**So the transgene is not silenced.** If the message had fallen, the attenuation would need no
further explanation and the dose argument would collapse into "the transgene switched off".

**1.2 The programme rescales by about half.** ×**0.487** over the 2,648 genes MYC moves at 6W
(R² 0.51) -- the rate of record, and the 0.49 that ED Fig. 1h's legend already quotes. Shape is
preserved: per-tier slopes 0.53-0.83, coherence at the 100th percentile of gene-label shuffles.

> **Quote the right rate.** Four exist, all reading "about a half", none substituting for another
> (`2026-09-21_two_timeline_verification.md` section 8): **0.450** whole transcriptome, **0.487**
> the genes MYC moves, **0.552** the 143 MitoPathways on the content ruler, **0.43** the eight
> pro:anti ratios. The hard-coded `GLOBAL_RATE = 0.55` in scripts 42 and 44 is **none of these** --
> it is a typed, rounded copy of the mitochondrial slope, and the old "0.55-fold transcriptomic
> rescaling" sentence is wrong.

**1.3 The protein falls in BOTH genotypes, and harder in the wild type.** From the author's
`1b Myc WB Quant` sheet (`NM_all_figs_corrected.prism`, 2026-11-06; one-way ANOVA, Tukey):

| group | n | mean (AU) | s.d. |
|---|---|---|---|
| WT 6W | 8 | 146 | 35.2 |
| MYC+ 6W | 8 | 197 | 21.9 |
| WT 12W | 6 | **62.4** | 14.3 |
| MYC+ 12W | 8 | **133** | 15.3 |

- **MYC+ falls 197 -> 133 = −32%** (adjP 6.3e-05).
- **WT falls 146 -> 62.4 = −57%** (adjP 4.1e-06).
- The genotype ratio **widens** at protein level, 1.35 -> 2.13, mirroring the widening message gap.
- **12W MYC+ is indistinguishable from 6W WT** (adjP **0.70**): by twelve weeks the transgenic
  gland carries about the MYC protein the normal pubertal gland had.

**This is a read-out of the author's own file, not a new analysis, and it is the single most
useful fact for the reply.** It is currently not stated in the manuscript, which frames the
decline as something happening to the transgenic gland.

## 2. The MYC/MAX/MXD network: excluded

MYC cannot bind an E-box without MAX and competes for the same sites with MXD/MNT, so a fall in
MAX or a rise in the repressors would attenuate output at constant MYC. **Neither happens.**

- **ED Fig. 1a** (ten genes, both ages): nothing clears significance at either age, and no
  interaction is below **padj 0.83**. `Max` is +0.22 (padj 0.21) at 6W and +0.06 (0.87) at 12W.
- **The defensive sub-point is worth as much as the negative.** The standard errors of the two
  ages agree to within **0.024 log2 across all ten genes**, median ratio 1.01 -- so twelve weeks
  is not simply the noisier cohort. This is asserted in the panel script, not assumed.
- **Every E-box binder is flat** (2026-07-25 scan: Usf1/2, Bhlhe40/41, Tfe3/Tfeb/Mitf, the whole
  MAX/MXD/MNT/MGA/MLX network). Closest movements are Srebf1 (+0.61, padj 0.050) and Usf2
  (+0.21, p 0.051), neither Myc-specific.
- **A 111-gene repressor scan** (nuclear receptors, NCoR/SMRT, Sin3, CoREST, Polycomb, repressive
  E2Fs, DREAM, circadian, mammary TFs): six move in the WT background at FDR < 0.05, **zero have a
  significant interaction**.

### 2.1 "Notwithstanding some factor significantly change" -- the three that did, and why each was demoted

| lead | what it does | why it does not carry the explanation |
|---|---|---|
| **`Id4`** +1.53, padj 0.0013 (WT) | an HLH sequestrator; the most attractive hit | `Id1/2/3` are flat (−0.40, +0.23, −0.28, all ns); `Id4` rises inside a **myoepithelial maturation signature** (`Oxtr` +1.72 padj 2.7e-4; `Myh11` +0.85 padj 0.024; `Trp63` +1.17) -- composition, not a MYC antagonist. Mechanistically wrong anyway: **ID proteins sequester E-proteins (class A bHLH); MYC:MAX is bHLH-LZ and is not an ID target** |
| **`Mlxipl`** −1.03, padj 0.037 (MYC+ time) | ChREBP, an E-box competitor | baseMean **45**; WT falls too (−0.39); interaction padj 1.00; mammary is not a ChREBP tissue |
| **`Mlxip`** −0.27 padj 0.061 at 12W | MondoA, the one shaped lead | **the sign is backwards** -- fewer E-box competitors means *less* competition, so MYC should hold up better, not worse. Interaction padj 0.83, and `Txnip`, MondoA's most sensitive output, is flat |
| *(two of the six repressor movers point the wrong way)* | | `Esrrb` +1.58 is an **activator** going up; `E2f7` −0.83 is a **repressive** E2F going down = *less* repression |

**One route survives in principle and cannot be tested by RNA-seq:** MLX is shared between
MondoA/ChREBP and MXD1/MXD4/MNT, so losing MondoA could free MLX into **MXD:MLX repressive
complexes** at E-boxes -- with no repressor transcript changing. That is a protein-stoichiometry
hypothesis; only co-immunoprecipitation or occupancy can read it.

## 3. Proteolysis / turnover machinery: excluded at transcript level, and pointing the wrong way

Script 41 (`outputs/myc_stability_panel/README.md`), 31 genes: phospho-degron kinases and
phosphatases, E3 ligases, deubiquitinases, acetyl/deacetylases, SUMO enzymes, O-GlcNAc.

Five down at FDR < 0.05 (`Trim6` −1.01, `Pias4` −0.40, `Fbxw7` −0.35, `Usp36` −0.33,
`Huwe1` −0.24) and two up (`Ep300` +0.25, `Sirt1` +0.76). **None of it supports the hypothesis**,
for four reasons:

1. **Nothing is Myc-specific.** Interaction padj > 0.05 for **every** gene in the panel. Several
   hits move in WT as well (`Pias4` −0.27 p 0.047; `Huwe1` −0.13 p 0.049; `Sirt1` +0.56 p < 0.001),
   i.e. they sit on the batch-confounded time axis, where a contrast is described, not claimed.
2. **The panel does not move more than the genome.** Against expression-matched background genes:
   on the MYC+ time axis the panel median is −0.087 vs a matched null of −0.071 (**42nd percentile,
   p 0.83**); on the interaction, −0.111 vs −0.077 (**24th percentile, p 0.49**). This is the global
   attenuation passing through the panel, not a targeted change.
3. **The directions predict a MORE stable MYC.** `Fbxw7` (−0.35) and `Huwe1` (−0.24), the two
   best-established degradative E3 ligases for MYC, are both **down**. The destabilising and
   stabilising arms also move together (E3 ligases median −0.25; deubiquitinases −0.24), so the
   ubiquitin balance does not shift.
4. **The canonical T58/S62 phospho-degron axis shows no transcriptional change at all** --
   `Mapk1`, `Mapk3`, `Gsk3b`, `Pin1`, `Ppp2ca`, `Ppp2cb`: median **+0.04**, smallest p 0.182.

**Read-only extension (2026-11-10, not in a numbered script -- re-derive before citing).** The
26S proteasome subunits and the calpain arm were read on the same contrasts and behave the same
way: `Psma3` −0.41 (padj 0.024) on the MYC+ axis but −0.48 (padj 0.0083) in WT, i.e. shared;
`Psmb1/5/6/7`, `Psmc1`, `Psmd1/4` flat; `Capn1`, `Capn2`, `Cast` flat or shared. **No interaction
below padj 1.00 anywhere in either set.**

## 4. Autophagy: excluded, and it moves the other way

Script 39 (`results/lesner_mechanism_test.rds`), fGSEA NES with BH-adjusted *P*:

| set | MYC+ vs WT 6W | 12W vs 6W in MYC+ | 12W vs 6W in WT | interaction |
|---|---|---|---|---|
| AUTOPHAGY | **+1.62** (0.021) | +0.95 (**0.62**) | +1.66 (0.024) | −1.58 (0.031) |
| MITOPHAGY | **+1.56** (0.021) | −0.75 (**0.87**) | +1.34 (0.20) | −1.46 (0.067) |

Three readings, all against the hypothesis:

- **MYC raises autophagy; it is not a route that removes MYC.** The genotype effect is positive at
  both ages.
- **It does not rise over the window in MYC+** (padj 0.62 and 0.87) -- it **fades with the rest of
  the programme** (interaction −1.58). Where it does rise, it rises in **WT** (+1.66, padj 0.024),
  i.e. developmentally and not Myc-specifically.
- **The genotype rise is not autophagy-specific.** Against an expression-matched MitoCarta
  background it is *below* average for a mitochondrial gene: MITOPHAGY observed +0.24 vs background
  +0.39 (**z −1.37, p 0.92**); DYNAMICS +0.27 vs +0.39 (**z −2.69, p 1.00**). It is the ~20%
  wholesale mitochondrial content effect. The matched null is anti-conservative, so the
  non-specific verdict is safe.

**Read-only extension (2026-11-10 -- re-derive before citing).** Gene level on the MYC+ time axis:
`Sqstm1`, `Becn1`, `Ulk1`, `Tfeb`, `Atg5`, `Atg7`, `Atg12` all flat; `Atg16l1` **−0.32**
(padj 0.016), i.e. **down**; `Map1lc3b` +0.28 (padj 0.029) but WT +0.41 (padj 0.00053) -- shared.
Lysosomal `Lamp1`, `Lamp2`, `Ctsb`, `Ctsd`, `Ctsl` flat. No interaction below padj 1.00.

## 5. The limitation all three share, and the proposal that follows

**Every scan above is transcript-level, and MYC protein stability is controlled
post-translationally.** A constant `Fbxw7` mRNA with altered T58 phosphorylation would look
exactly like this result. Script 41's own README states it: *"No transcriptional evidence" is not
"not happening."* What the three scans establish is that **no large transcriptional reprogramming
of the degradation machinery, the E-box network or the autophagic apparatus accompanies the
decline** -- which is a real and quotable exclusion, and is also precisely what a
post-translational mechanism predicts.

The proposal, in order of what the data support:

### 5.1 The decline is developmental, not oncogenic (strongest, and already in the author's data)

Endogenous MYC in the wild-type gland falls **−57%** while the transgenic protein falls **−32%**
(section 1.3). Whatever lowers MYC therefore acts on the protein **irrespective of which promoter
made the message**, and the variable is the pubertal-to-adult transition of the gland, not the
transgene. This reframes the reviewer's question from "why did your model weaken" to "what does
the maturing mammary gland do to MYC protein" -- a question about normal biology, to which the
honest answer is that it is known to happen and we did not set out to map it.

### 5.2 The mechanism class is post-translational control of MYC half-life

The canonical axis -- mitogenic input -> ERK -> S62 (stabilising), GSK3β -> T58 -> FBXW7
(destabilising) -- is the obvious candidate, because the window is exactly when terminal end buds
regress and pubertal mitogenic drive winds down in the wild type (`Aldh1a3` −2.22 padj 5.5e-5;
TEB-versus-ductal signature −0.45, p 0.0028; `Oxtr` +1.72 padj 2.7e-4 marking the differentiation
shift). **It fits every negative above by construction:** it changes no transcript, so sections
2-4 are blind to it. It also sits comfortably beside section 3's finding that the degradative E3s
are *down* -- if the cell is not raising its degradative capacity and MYC protein still falls, the
control is on the substrate's phospho-state (or on a route the panel does not carry), not on the
amount of ligase.

**This is a hypothesis, not a result.** It is the right thing to offer a reviewer as the likely
explanation, and the wrong thing to assert in the manuscript.

### 5.3 Compositional dilution (conditional -- see section 6)

The 12W gland holds a different cell mix: fewer proliferative/TEB cells, more differentiated
myoepithelial and luminal cells (the `Oxtr`/`Myh11`/`Trp63` signature). If MYC is quantified in
**whole lysate**, a shrinking MYC-high compartment lowers the measured protein with no per-cell
change. **Standing tension:** dilution should shrink the *Myc* message gap too, and the gap
*widens* (+1.68 -> +1.82). **This route is largely removed if the quantification is per-nucleus**
-- which section 6 asks the author to settle.

### 5.4 The fallback, and the strongest card in the reply

**The paper's conclusion does not rest on the mechanism.** The argument runs: dose halves -> the
programme rescales -> **but MYC-ER delivers equal or higher MYC at 12W and still loses the
apoptotic response while proliferation is maintained** -- in animals that had never seen MYC
before induction, so the change cannot be MYC-driven selection. The permissive window is therefore
**substrate competence, not dose**, and why the dose changed becomes an interesting side question
rather than a load-bearing one. A reviewer who presses on the protein fall is pressing on a step
the design was built to make non-essential, and the reply should say so plainly.

> **Banner correction to `2026-07-25_attenuation_competing_hypotheses.md`.** That note's banner
> reads "H3 wins outright ... a ~50% fall in effective MYC and the x0.55 rescaling are the same
> number". Both halves need softening: the protein fall in MYC+ is **−32%** on the current
> quantification (section 1.3) and the transcriptome rate of record is **×0.487** (section 1.2).
> The qualitative claim -- a roughly multiplicative reduction in effective MYC produces a roughly
> proportional rescaling with preserved shape -- is unaffected. The numerical coincidence is not
> as tight as the note's rhetoric, and should not be quoted as such.

## 6. Two things the author should settle before a reviewer finds them

1. **"~50%" against the figure.** The manuscript says MYC protein "fell by ~50%". The current
   Prism quantification gives **−32% in MYC+** (and −57% in WT). Either the text takes the
   measured number, or the discrepancy is explained (different quantification, different cohort).
   As it stands the claim and the panel behind it do not agree.
2. **What Fig. 1b's left graph actually measures.** The axis reads
   **"MYC (mean nuclear intensity, AU)"** while the Prism sheet is named **"1b Myc WB Quant"**,
   with n = 8/8/6/8 (the RNA-seq is n = 6). This decides section 5.3: **per-nucleus IHC largely
   removes compositional dilution** and pushes the explanation to 5.2, whereas **whole-lysate
   western keeps dilution live**. The reply changes depending on the answer, so it should be
   pinned down first.

## 7. The reply, in the order it should be made

1. The *Myc* **transcript does not fall** -- the transgene is intact (1.1), and the standard
   errors show twelve weeks is not the noisier cohort (2).
2. **MYC protein falls in both genotypes, and more steeply in the wild type** (1.3). This is a
   feature of the maturing gland, not a failure of the model.
3. We **excluded** three routes at transcript level, two of which point the wrong way: the E-box
   network is unmoved (2), the degradative machinery is *down* and no more changed than background
   (3), and autophagy **falls** with the programme rather than rising (4).
4. We **state the limitation**: these are transcript measurements and the regulation is
   post-translational, so the likely explanation is reduced stabilising input to MYC protein as
   pubertal mitogenic signalling winds down (5.2). We offer it as a hypothesis.
5. We **note that the conclusion is independent of it**: MYC-ER holds dose constant and reproduces
   the asymmetry, so the permissive window is competence (5.4).
6. The experiments that would settle it: **phospho-T58/S62 immunoblot at both ages, a
   cycloheximide chase, and MYC CUT&RUN or ChIP** to separate "less protein" from "same protein,
   less occupancy". If offered as a revision experiment, the cycloheximide chase is the cheapest
   decisive one.

## 8. Supporting literature, with DOIs for Paperpile

- **Lorenzin, F. et al.** Different promoter affinities account for specificity in MYC-dependent
  gene regulation. *eLife* **5**, e15161 (2016). doi:10.7554/eLife.15161 -- why halving MYC cannot
  halve every target equally. Already cited in the manuscript, and the reason the retentions order
  as they do: MYC core 0.74-0.77 > global ~0.49 > mitochondrial arms 0.50-0.58 >> PUMA −0.09.
  Affinity-ordered occupancy explains everything above the line and **cannot produce a sign
  reversal**, which is why PUMA is the residual worth an experiment.

## 9. Provenance

- DESeq2 values: `results/interaction_results.rds`, raw (unshrunken) log2 fold change with
  IHW-adjusted *P*, `~ timepoint * myc_status`, n = 6 per group. Batch = timepoint, so the
  within-genotype time contrasts are **batch-confounded and descriptive**; the interaction is the
  batch-clean term, and it is null for every gene cited here.
- Protein and mRNA group values: `NM_all_figs_corrected.prism` (2026-11-06 17:36) on the local
  Drive mount `/Users/gs/G/data/MK_myc_2022/paper figures/NM_Figs_to_submit/`, sheets
  `1b Myc WB Quant` and `1b Myc mRNA`, read with
  `myc_human_exploratory/docs/legends/tools/prism_inventory.py`. Percentages in 1.3 are computed
  from the group means printed there.
- Matched-null percentiles, the per-arm tables and the T58/S62 summary: script 41's generated
  `outputs/myc_stability_panel/README.md` (2026-07-25), which is written from its own run.
- fGSEA NES and the mitophagy specificity null: `results/lesner_mechanism_test.rds` (script 39).
- The 2026-07-25 E-box, 111-gene repressor, identity-marker and MLX scans were **read-only
  scratchpad explorations**. They carry no provenance beyond
  `2026-07-25_attenuation_competing_hypotheses.md` section 4 and **must be re-derived inside a
  numbered script before any of their numbers is cited in the manuscript**. The same applies to
  the two read-only extensions marked in sections 3 and 4 of this note.
- Blot-derived statements about MYC-ER (equal or higher dose, loss of the apoptotic response with
  proliferation maintained) are the author's bench data, reported here as given.
