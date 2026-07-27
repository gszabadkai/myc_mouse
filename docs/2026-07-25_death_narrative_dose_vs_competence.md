# Dose versus competence: a single line through the death data

**Status:** synthesis note, written 2026-07-25 after the MYC and apoptotic-protein blots.
Supersedes the framing (not the content) of `docs/2026-07-25_attenuation_competing_hypotheses.md`,
which was written before the blots. **Reconciled against `results/priming_arm_teb.rds`
(`scripts/42_priming_arm_and_teb_substrate.R`, author-run 2026-07-25)** — every transcriptome
number below is that script's output. Four claims in the pre-run draft did not survive
reconciliation and are corrected in place: the purity-adjusted interaction p is **0.017**, not
0.006; the BAX absorption figure was computed on a different axis and is restated; the
AP_TEB regulon finding is a **context** statement, not a TF-specificity one; and the absorption
test is dropped because its control failed. All four are flagged where they occur.

---

## 0. The epistemic contract

Everything below is written to one rule, and the manuscript should be too:

> **The in-vivo transcriptome generates the hypothesis. The cell perturbations prove it.**

No correlation in these 24 mice is asked to carry a causal claim. That is not a concession — it
is the correct division of labour, and it is what makes the section defensible. The mouse data has
n=6 per group, a batch-confounded time axis, and a correlation ceiling that has defeated every
cross-sectional coupling this project has attempted. The cell data has both directions of a
perturbation (PGC1a gain sensitises; selected and passaged cells lose biogenesis and resist) in
three systems. Causality belongs to the second, convergent observation to the first.

Read this way, the in-vivo material has a real and defensible job: it shows that **the two things
the cell work links causally — the OXPHOS/PGC1a arm and the PUMA trigger — both go down in the
gland over the window in which death is lost**, and it shows that on more than one ruler. That is
an indication sufficient to motivate the experiment, which is exactly what it is used for.

---

## 1. The two experiments answer two different questions

The confusion dissolves once the designs are put side by side.

| | MMTV-Myc, 6W vs 12W | MYC-ER, induced at 6W vs at 12W |
|---|---|---|
| MYC protein | falls ~50%, transcript constant | equal or higher at 12W |
| prior MYC exposure | 12 weeks of it | **none** — acute induction |
| what varies | the dose | the tissue |
| what it measures | what *less* MYC does | what the *same* MYC does to a different substrate |
| the answer it gives | **dose** | **competence** |

**MMTV-Myc is a dose experiment.** A ~50% fall in MYC protein at constant message, and a x0.55
rescaling of the transcriptional programme. Those two numbers are the same number. And a uniform
multiplicative reduction in effective MYC is exactly what predicts the pattern script 40 measured:
the 12W profile is the 6W profile at 55% amplitude, same pathways, same rank order, per-tier
slopes 0.53-0.83, coherence at the 100th percentile of gene-label shuffles. It also retrospectively
explains every negative in the corpus — script 34 found no programme attenuating more than
expression-matched background; the 2026-07-25 scans found no repressor induced and no E-box
competitor moving; Issue #5 found no Myc-independent axis absorbing the attenuation. There was
nothing to find, because **output per unit MYC is unchanged**. (This reads the programme as
linear in MYC protein; the coincidence of ~50% and x0.55 is the evidence for that, not an
assumption made in advance.)

The mito-protein blots and the BAX/BIM/PUMA protein falls in Myc+ are the same x0.55 propagating
from message to protein. They need no separate explanation.

**MYC-ER is a competence experiment**, and one design feature does the decisive work: those 12W
animals had **never been exposed to MYC before induction**. Whatever changed in them cannot be
MYC-driven selection. It has to be a property of the maturing gland.

So propositions 1 and 2 in the original question are not rivals. Both are true, of different
observables. The MMTV attenuation is dose. The loss of death sensitivity at constant dose is
substrate.

---

## 2. The coincidence model, and how the sub-hypotheses stand

The right frame is a **two-input (coincidence) model**: MYC-driven death requires **MYC *and* a
competent mitochondrial substrate**, and the substrate input declines developmentally, so the
same MYC no longer clears the threshold.

| hypothesis | verdict |
|---|---|
| (i) MYC engages PGC1a more at 6W than 12W | genuinely **dose**-flavoured — cannot by itself explain MYC-ER, where dose is equal or higher |
| (ii) the PGC1a/biogenesis input is lower at 12W, so the same MYC does not clear the threshold | **the coincidence model. MYC-ER is its test, not its refutation** — it holds MYC fixed and varies exactly the input this hypothesis names |
| (iii) developmental / compositional | **compatible and parallel** — a second candidate for the same second input, non-exclusive with (ii) |

The coincidence frame is better than "the programme decays" or "the background opposes" on three
counts, and all three are worth keeping in the paper:

- **It explains the asymmetry with no extra assumption.** Proliferation needs one input; death
  needs two. A fall in either input therefore costs death superlinearly and proliferation
  linearly — which is what MYC-ER measures (80% versus 30%).
- **A threshold explains the *shape* of the PUMA result.** PUMA priming does not merely attenuate,
  it collapses *through zero* (retention -0.09) while everything else retains ~0.55. A product of
  two declining inputs crossing a threshold gives that; a uniform dose reduction does not.
- **It is the manuscript's title thesis stated literally** — mitochondria *integrating* an
  oncogenic and a metabolic input.

One correction to the framing of (iii): the direction in our data is **not** "an increase in
luminals". It is loss of the pubertal/TEB and luminal-progenitor compartment and maturation of the
differentiated one — `Aldh1a3` -2.22 (padj 5.5e-5), `Nrg1` -1.82 (0.014), `Areg` -1.30 down;
`Oxtr` +1.72 (padj 2.7e-4), `Myh11` +0.85 (0.024), `Trp63` +1.17 up.

---

## 3. MYC builds a death-ready mitochondrion — and loses the trigger

At 6W, on the clean genotype axis, MYC does four things at once:

| what | genes | Myc effect at 6W (padj) |
|---|---|---|
| loads the executioner | **Htra2** | **+1.49 (2.6e-08)** — the strongest MYC-induced apoptotic gene in the dataset |
| loads the effector | Bax | +0.48 (1.4e-4) |
| removes the brake | Bcl2l1 (Bcl-xL) | **-0.42 (0.04)**; Xiap -0.32 (0.078) |
| builds the organelle | Tomm22 / Tomm40 / Atp5f1a / Hspd1 / Esrra | +0.94 (4.5e-5) / +0.90 (2.4e-5) / +0.87 (5.9e-4) / +1.12 (8.7e-4) / +0.68 (1.8e-6) |

HTRA2 is worth pausing on: it is the intermembrane-space protease that
`docs/library_reference/PUMA-and-its-relationships.md` places as the *physical executioner*
downstream of the PGC1a/PUMA cascade, released through the BAX/BAK pore. MYC induces it more
strongly than anything else in the apoptotic machinery here.

At 12W MYC still does all of it, at ~0.55-0.8 of the amplitude — the dose effect. **Except the
trigger.**

| priming ratio (per mouse, log2) | Myc effect 6W | Myc effect 12W | retention | conditional matched-pair null |
|---|---|---|---|---|
| Bax : Bcl-xL | +0.90 (p 1.2e-4) | +0.50 | **0.550** | 50th pct (p_emp 0.50) — *exactly* the global rate |
| Bid : Bcl-xL | +0.95 (p 0.0036) | +0.40 | 0.419 | 70th pct (p_emp 0.30) |
| Bak1 : Bcl-xL | +0.69 (p 0.0052) | +0.18 | 0.268 | 82nd pct (p_emp 0.18) |
| **PUMA : Bcl-xL** | +0.67 (p 0.0089) | **-0.06** | **-0.090** | **92nd pct (p_emp 0.082)**, interaction p 0.037 |
| PUMA : Bcl-xL, purity-adjusted | +0.61 (p 0.0018) | | | **interaction p 0.017** |

Four features make this more than a lone p-value.

1. **It contains its own control.** The panel shares a denominator. BAX priming fades at the
   global rate; PUMA priming collapses. "Isn't this just the x0.55 attenuation?" is answered from
   inside the panel, not by appeal to an external null.
2. **The null reproduces the global rate exactly.** Random expression-matched gene pairs retain a
   median of **0.55** — script 40's rescaling slope, rebuilt from arbitrary genes. The instrument
   is calibrated, and this was the script's built-in positive control.
3. **Adjusting for composition strengthens it.** Epithelial and immune covariates move the
   interaction from p 0.037 to **p 0.017**, while BAX stays non-significant (0.146 -> 0.185). So
   it is not contamination — and PUMA's raw correlation with the immune composite (+0.34) was the
   obvious worry.
4. **A second, ratio-free framing agrees.** `Bbc3`'s own aligned attenuation (script 34's
   `sign(lfc6)*(lfc6-lfc12)`) sits at the **94.6th percentile** of genes matched on both baseMean
   and |lfc6| (p_emp 0.054), while `Bax` sits at the 51st — and `Htra2`, the most strongly
   Myc-induced apoptotic gene here, at the **37th**, i.e. attenuating *less* than matched genes.
   So this is not a property of Myc-induced apoptotic genes in general.

**The specificity is graded, not binary, and the paper should say so.** Two other ratios have more
extreme retentions than PUMA — `Bmf:Bcl-xL` (-0.83) and `Bcl2l11:Bcl-xL` (-0.45) — but both have
**non-significant 6W effects** (p 0.41 and 0.46), and a retention is a ratio of two noisy
quantities: you cannot lose an effect you never had. Restricted to the pairs with a real 6W effect
— Bax, Bid, Bak1, PUMA, Bax:Mcl1 — the retentions form a gradient (0.55, 0.42, 0.27, **-0.09**,
0.84) and **PUMA is the only one that reverses sign** and the only one whose conditional p_emp
falls below 0.1. The honest claim is that the BH3-only arms lose priming to differing degrees and
PUMA loses it completely, not that PUMA is uniquely affected.

**And PUMA was pre-specified.** The PGC1a experiments induce OXPHOS and PUMA "but not other
BH3-only proteins". The one BH3-only the cell work singles out is the one whose in-vivo priming
collapses. This was not a fishing expedition through the BCL-2 family.

**One sentence for the paper: MYC builds a death-competent mitochondrion at both ages; what it
loses by twelve weeks is the trigger.**

That maps onto the cell data with unusual precision. Bcl-xL rescue shows the machinery is still
built — OXPHOS and BIM/BAX go up, more so than without it — and the cells still do not die: the
block is at execution, not at transcription. PGC1a supplies PUMA, and they die.

---

## 3.5 The three convergent strands — what generates the hypothesis

Section 3 gives one of them. Here they are together, because it is their convergence, not any one
of them, that motivates the cell experiments.

### Strand 1 — OXPHOS goes down on *both* rulers, and the priority ruler cannot be dose

Per Level-1 tier, medians:

| tier | content Myc@6W | content WT 6-12W | **priority WT 6-12W** | **priority Myc+ 6-12W** |
|---|---|---|---|---|
| **OXPHOS** | +0.456 | -0.142 | **-0.069** | **-0.134** |
| Protein import | +0.518 | -0.056 | -0.055 | -0.085 |
| Central dogma | +0.469 | +0.051 | -0.011 | -0.050 |
| Metabolism | +0.397 | +0.084 | +0.037 | +0.012 |
| Dynamics | +0.236 | +0.030 | +0.003 | +0.022 |

**OXPHOS is the most de-prioritised Level-1 tier on both temporal axes** — figS6's reading,
tabulated. And within OXPHOS it is the **structural subunits** that carry it (priority -0.109 WT /
-0.153 Myc+; content -0.255 / -0.405) while the **assembly factors do not move** (-0.002 /
-0.029). The subunits are the PGC1a / NRF1 / ERRalpha output. So what loses priority is precisely
the arm PGC1a builds.

**The simplest form of this, added 2026-07-26.** On the temporal fGSEA NES for
`MITOCARTA_OXPHOS_SUBUNITS` (`results/attenuation_decomposition.rds$decomp`, script 29): **WT
-2.670 versus Myc+ -2.715.** The normal gland withdraws from the respiratory subunits essentially
as strongly as the Myc+ gland does, and **no other arm in that table behaves this way** (MYC-target
core WT -1.40, mitoribosome -0.61, nucleotide -0.66, TCA **+0.58**, lipid **+1.49**, amino-acid
**+1.83**). Use this pair, not the tier medians, when the point is that the *substrate* falls.

**And read it in the right frame.** Issue #6 (`31_attenuation_mechanism.R$decomp_conv`) reports
OXPHOS as the arm where the wild-type background does **not converge** toward Myc (conv -44%, fade
144%), and that is correct — but it is a statement about the **genotype gap**. Myc+ sits *above* WT
on OXPHOS, so a wild-type gland that falls is moving *further below* it: anti-convergence and a
large wild-type decline are the same fact seen from two sides. The **gap frame** is right for
Issue #6's question (why the Myc effect halves); the **substrate frame** is the one this section
needs, because MYC-ER holds MYC fixed and samples the gland's own state, not the gap. The gloss
"WT OXPHOS barely declines" was a frame error and does not survive — see
`docs/2026-07-26_introduction_alignment_and_the_question.md` section 1.1.

**Why this strand carries weight the content ruler cannot.** mitoPPS is a *pairwise ratio* within
the mitochondrial compartment, so a uniform x0.55 scaling of the whole programme **cancels** in it.
The content ruler is dose-dominated — that is section 1's whole point. The priority ruler is not,
and OXPHOS still falls on it. **The de-prioritisation is not the dose effect.**

### Strand 2 — the trigger goes down

Section 3's table: PUMA:Bcl-xL priming retains **-0.09** against a global rate of 0.55 (92nd
percentile of expression-matched pairs; interaction p 0.017 purity-adjusted), while Bax:Bcl-xL
retains exactly 0.55.

### Strand 3 — the two are coupled, in the one space where a coupling is readable

Script 38 PART B, already run and on disk: **OXPHOS-ratio couples to priming-ratio (PRO-ANTI) at
design-adjusted Spearman +0.70 against a design-adjusted null of 0.48** (raw +0.45 against a raw
null of 0.42), beating its null on both. `redox` is the control axis and **fails** (+0.12, does
not beat its null). Crucially this coupling is **neither circular nor ceiling-bound** — a ratio
cancels the shared mito-content mode, which is exactly what made the GSVA-space version
(script 34's negative control) uninterpretable.

### The bounds on all three

- **The apoptosis *composites* show nothing temporally.** Apoptosis-PRO priority +0.054, ANTI
  +0.046 — they move together, so the balance is flat. The death finding is **gene-specific to
  PUMA**. That is why script 34's composite analysis found nothing, and it is what the cell work
  predicts: PGC1a induces PUMA "but not other BH3-only proteins", so a 25-gene PRO composite is
  the wrong instrument.
- Strand 3's coupling is **mid-pack** — 52nd percentile of OXPHOS's couplings to all 141 mito
  pathways — and its genotype delta is p 0.080. A ranking lead.
- Strand 1's temporal axes are **batch-confounded**. Their mitigation is that the protein blots
  move the same way, off the RNA batch entirely.

**Together:** the two things the cell work links causally both go down, on rulers that fail in
different ways, and they co-vary where co-variation is measurable. That is an indication. It is
not proof, and section 0 says what proves it.

## 4. Why PUMA, and why developmentally

Three routes, in decreasing order of what our data can support.

**4.1 The PUMA promoter is reachable in the TEB context and, in this resource, nowhere else.**
*(Corrected after the run — the pre-run draft over-read this as TF specificity.)* The Gray/ChEA
shortlist (`data/genesets_from_library/gray_chea_mito_tf_shortlist.csv`) puts `Bbc3` in the target
set of **71 transcription factors — and every one of them is `AP_TEB`.** It appears in no other
context. Seventy-one TFs including TP53, ESR1, JUN, FOS, SOX2 and POU5F1 is not TF specificity;
it is what a **broadly bound, accessible promoter in one cell context** looks like. That is still
the claim the model needs — a context in which PUMA is reachable at all — but it is a statement
about chromatin, not about any TF's target list.

The specific claim that does survive is sharper and is in the biogenesis-TF table: **MYC and E2F1
reach `Bbc3` only in `AP_TEB`. In every other context they reach `Bid`, `Bok`, `Bak1` or `Pmaip1`
instead** (MYC BA_LE -> Bid, AP_LE -> Bok, HS_LE -> Bok, HS_TEB -> Pmaip1; E2F1 BA_LE ->
Bak1/Bid, AP_LE -> Bak1, HS_TEB -> Pmaip1). So **which BH3-only MYC can touch is
context-dependent**, and PUMA is the one it can touch in the terminal end bud. And
`docs/library_reference/Gray_et_al_developmental_TFS_selection.md` already records, months ago and
for an unrelated purpose:

> "Roster TFs lean pro-apoptotic: MYC -> Bbc3/Bok/Bid; E2F1 -> Bak1/Bid/Pmaip1; ESR1/ESRRA ->
> Aifm2; GABPA -> Bak1; **NRF1 -> Bbc3**. I.e. the biogenesis TFs that build mito mass also touch
> pro-apoptotic effectors — a candidate co-regulation of apoptotic sensitivity."

and

> "the `HS_TEB` context is apoptosis-loaded (TEB lumen clearance by apoptosis during pubertal
> elongation)."

That is the PGC1a-to-PUMA hypothesis, in an independent resource, written before the question was
asked. NRF1 is a core PGC1a-axis biogenesis TF and PUMA is among its targets.

**4.2 The TEB compartment demonstrably regresses in our own WT data.**
`MG_TEB_VS_DUCTAL_HS_GRAY_UP` falls **-0.42 (p 0.0053)** purity-adjusted between 6W and 12W, and
every TEB set that moves is an **HS_TEB** one: ARNTL2 -0.36 (0.0090), HMGA2 -0.34 (0.023), GRHL3
-0.28 (0.011), MYC_HS_TEB -0.21 (0.020), ELK3 -0.28 (0.073), E2F7 -0.27 (0.074), E2F1_HS_TEB
-0.19 (0.070). The markers agree (section 2). The HS_TEB context — the one the library flags as
apoptosis-loaded — is the only one that declines, which is the coherence the claim needs.

So the proposal is: **the promoter context in which MYC and the biogenesis TFs reach PUMA is a
property of the pubertal gland, and it goes away with the terminal end buds.** At 12W the same
MYC, in a matured gland, builds the same organelle and cannot pull the trigger. That is
developmental by construction, which is what MYC-ER requires.

**4.3 PGC1a itself is not readable here.** `Ppargc1a` has baseMean ~30 in this tissue — too low to
interpret. The ESRRA/NRF1/GABPA lanes are the activity proxy, as in script 38 PART C. `Esrra` is
strongly MYC-induced at 6W (+0.68, padj 1.8e-6) and not significantly at 12W (+0.37, p 0.12).

---

## 5. Three negatives that must travel with the story

Stated as prominently as the positives, because each is a place the paper could otherwise
overclaim.

1. **The cross-sectional OXPHOS-PUMA bundle fails.** Per-mouse, the nuclear-OXPHOS composite
   correlates with PUMA at **+0.00 at 6W (50th percentile of the ambient distribution) and -0.40
   at 12W (25th percentile)** — at or *below* what an arbitrary gene gives. BAX does better
   (+0.66 / +0.45, ~78-80th pct) but the ambient median is already 0.52/0.41. **The OXPHOS-PUMA
   co-regulation is carried by the cell perturbations, never by in-vivo correlation.** The mouse
   contributes the *contrast*, not the coupling. (Six near-identical mice vary in mitochondrial
   content for compositional and technical reasons that need not carry PUMA, so this is a weak
   test — but it is the test, and it failed, so the claim is not available.)
2. **PUMA is not a TEB transcript.** `Bbc3` is absent from every TEB-vs-ductal signature. The TEB
   claim in 4.1 is about **TF-binding context**, not about where PUMA is expressed. Writing it the
   other way would be wrong.
3. **The compartment shift is marker-level, not compartment-level.** The MEC consensus composites
   (LASP, LHS, BMYO, LP) all move in the expected directions but none is significant. Only single
   markers clear FDR. So: maturation is documented; a wholesale change in cell proportions is not.
4. **The mediation test failed its own control and is not usable.** *(Corrected after the run.)*
   Script-30-style absorption of the PUMA interaction by each candidate axis gives -27%
   (OXPHOS level), -25% (OXPHOS priority), +13% (TEB) — and **+41% for `redox`, the negative
   control axis**, which "absorbs" more than any real one. Adding any axis with the right noise
   structure moves the interaction, so absorption is measuring nothing here. **Do not cite it in
   either direction** — including not as evidence against the biogenesis route.
5. **In vivo, BAX tracks mitochondrial mass and PUMA does not** — the inverse of the cell result.
   Per-mouse, the OXPHOS composite correlates with BAX (+0.66 at 6W, +0.46 at 12W, ~75-80th
   percentile of ambient) and not with PUMA (+0.00, -0.40, at or below ambient). The cell
   experiments say PGC1a induces PUMA specifically. State this tension rather than let a reader
   find it; the most likely reconciliation is that mitochondrial **mass** and PGC1a **activity**
   are not the same variable — and note that `Ppargc1a` itself is baseMean ~30 here and
   unreadable, so nothing in these data tests the activity version at all. That is where the cell
   perturbations live, and they are perturbations.
6. **Nothing fully separates from timepoint at n=24.** The best coincidence term (§6.5) drops from
   p 0.005 to p 0.088 once `tp*myc` is in the model; the TEB axis drops from p 0.018 to p 0.27.
   The candidate second inputs can be **ranked** and not established.

And the standing bound on all of section 4: **batch = timepoint.** Every WT temporal statement is
described, not claimed. Its one mitigation is real, though: the mitochondrial and PUMA protein
blots move the same way, off the RNA batch entirely.

---

## 6. What is excluded — and what that narrows

| candidate | verdict |
|---|---|
| anti-apoptotic buffering rises | **No.** Bcl2l1, Mcl1, Bcl2 all flat across the timeline. The brake does not go up; the trigger stops being pulled |
| p53 / ARF pathway attenuates | **No.** Trp53, Mdm2, Cdkn1a, Zmat3, Phlda3, Ccng1 all flat (Cdkn2a unmeasurable). And the iMMECs are p53-null and die anyway |
| an induced repressor / background antagonism | **No.** Every E-box binder flat; 111-gene repressor scan with zero significant interactions (2026-07-25 note) |
| PUMA's known transcriptional inputs | **No.** E2f1, Atf4, Ddit3, Foxo1 all flat. Foxo3 does rise in WT (+0.47, padj 0.0078) but PUMA does not follow |
| MYC-driven selection | **Not for MYC-ER** — those animals never saw MYC. Selection belongs to the tumour-progression half of the story |

These exclusions are worth as much as the positive: the mechanism is narrowed to **the
inducibility of the pro-apoptotic arm**, not the buffer, not p53, not a repressor.

---

## 6.5 The coincidence model, tested — and its control works

This is the one place the in-vivo data speaks directly to the two-input model rather than around
it, and it is the strongest thing in the analysis. The test is `lm(priming ratio ~ Myc x axis +
epithelial + immune)`: does MYC raise priming **only** where the substrate is competent?

| axis for PUMA:Bcl-xL | Myc x axis | p | with `tp*myc` in the model | permutation p_emp |
|---|---|---|---|---|
| **OXPHOS priority (mitoPPS ratio)** | **+2.79** | **0.005** | **+2.06 (p 0.088)** | **0.081 (92nd pct)** |
| TEB context | +1.15 | 0.018 | +0.53 (p 0.27) | 0.092 (91st pct) |
| OXPHOS level (log-expression) | +0.36 | 0.14 | +0.18 (p 0.44) | 0.153 |
| `redox` priority — **negative control** | +0.99 | **0.72** | +0.95 (p 0.66) | **0.491 (51st pct)** |

Four things to read from it. The **mitoPPS OXPHOS-priority axis is the strongest** — three times
the TEB term and eight times the level term. It is the **only one that partially survives the
timepoint it is confounded with**, where the TEB axis collapses. The **within-timepoint
permutation** (shuffling the axis inside each timepoint, so only the mouse-to-mouse link breaks)
puts it at the 92nd percentile. And the **control behaves**: `redox` is null on every reading and
sits at the 51st percentile of its own permutation — which is what makes the positives readable at
all. On BAX the same axis gives +1.21 (p 0.069), weaker and not surviving `tp*myc`.

**And the size objection is answered in the space where it should be asked.** The obvious challenge
to the model is that the developmental decline looks too small to matter: on the log-expression
level ruler, the WT 6->12W fall in OXPHOS is **-0.25** against a Myc effect of **+0.96** — 26%. But
on the mitoPPS priority ruler, where the dose effect cancels, the WT fall is **-0.109** against a
Myc effect of **+0.112** — **97%**. In the space where the second input should be measured, **the
developmental decline is as large as the entire effect of the oncogene.** The `redox` control axis
gives 70% by the same calculation, so this is a ranking statement rather than a clean one — but it
removes the objection that the substrate change is too small to be the second input.

The verdict: **the mitoPPS OXPHOS-priority axis is the best-supported candidate second input**,
marginal on every individual test (p 0.005 -> 0.088 with timepoint, permutation 0.081), with a
control that works. A ranking lead, in a project where every cross-sectional coupling has landed
in the same place — and, per section 0, exactly the weight the in-vivo data should carry.

## 7. Alternatives the paper must address

**7.1 Threshold uncoupling — the reviewer's first move.** Apoptosis requires a higher MYC
threshold than proliferation (the Murphy/Evan dose-threshold work; check the citation before
using). A 50% fall in MYC protein would therefore remove death preferentially, and the MMTV
asymmetry follows with no substrate change at all. Our data is consistent with it: proliferation
retains ~0.6 (Mki67 +1.04 -> +0.50, Plk1 +1.46 -> +0.58, Aurka +1.13 -> +0.67) while PUMA priming
retains ~0. **The answer is MYC-ER**: equal or higher dose, still 80% less death. The paper should
raise this objection itself and answer it, rather than wait for it.

**7.2 Replication-stress coupling.** MYC induces the DDR at 6W (Parp1 padj 5.6e-5, Cdc25a 0.0017,
Brca2 0.024, Atm 0.048) and several members collapse well below the global rate (Cdc25a 0.20,
Brca2 0.19, Wee1 0.15). A quiescent 12W gland offers less replication stress for MYC to convert
into death. This is a **co-mechanism, not a rival** — and it is compatible with everything above.

**7.3 Survivor bias, which cuts in our favour.** At 6W the Myc+ gland is actively dying, and the
cells that died are not in the RNA library. Bulk RNA therefore systematically **under-reports** the
6W pro-apoptotic induction, making the PUMA collapse a **lower bound**. This is the in-vivo
analogue of why Bcl-xL was needed in vitro to reveal the full BIM/BAX induction — the rescue
experiment removes the same bias experimentally. It probably also explains an awkward earlier
result (script 34: APOPTOSIS_PRO rising *less* than an average mitochondrial gene, z = -2.42).

**7.4 The experiment that would settle the coincidence model, and it is cheap.** The model makes
one prediction that no amount of further analysis can substitute for: **co-induce PGC1a (or
otherwise raise mitochondrial biogenesis) in MYC-ER at 12W, and death should be restored.** That
is the in-vivo form of the iMMEC PGC1a+MYC experiment, it uses a system already in hand, and it
converts the whole section from inference to test. The converse is equally available — knock the
biogenesis arm down at 6W and death should fall at constant MYC. Alongside it, **BH3 profiling at
both ages** is the phenotypic readout of the priming claim, since priming is what we measure and
death is what we infer. Note also that script 38 PART D concluded this question "needs an
ESRRA/PGC1a perturbation vs the Myc substrate, not n=24 bulk" — **that perturbation now exists in
the cell arm**, and the in-vivo version is the natural completion.

**7.5 The death is invisible to bulk RNA, and that must be said.** There is no efferocytosis or
phagocyte-clearance signature at 6W — Gas6 is actually *down* (-0.80, padj 0.0086). The death
phenotype comes from IHC and from the MYC-ER counts, and the transcriptome speaks only to the
molecular substrate for it. Priming is not death.

---

## 8. The contradiction the paper has to confront

`docs/library_reference/PUMA-and-its-relationships.md` records the canonical model:

> "When PGC-1a is abundant, it complexes with FOXO3a, guiding it to turn on antioxidant survival
> genes while suppressing its ability to trigger apoptosis. When PGC-1a is **lost** or
> downregulated ... unchecked FOXO3a freely binds to the PUMA promoter, causing massive PUMA
> upregulation."

That is the opposite of the experimental result: PGC1a **overexpression** induces PUMA and
sensitises to MYC-driven death, in iMMECs and in MYAZ. **The manuscript is inverting a
literature** — PGC1a as sensitiser rather than protector — and must say so explicitly and early,
with the muscle-atrophy and neurodegeneration contexts named as where the protective model comes
from. Reviewers will otherwise supply the contradiction themselves.

Our data does at least say the FOXO3-uncoupling model is not what is operating here: `Foxo3` rises
in the WT background (+0.47, padj 0.0078) and PUMA does not follow.

---

## 9. Proposed narrative order — SUBSTRATE FIRST

*(Reordered 2026-07-27 on the author's decision: start with what the gland does to itself, then what
Myc does on top of it. The previous Myc-first order is kept at the foot of this section, because the
reason it was rejected is worth remembering.)*

A reordering of what exists, not new work. The current write-up is organised by *analysis*; the
paper should be organised by *argument* — and, following section 0, the order should still be
**observation -> hypothesis -> causal test -> consequence**, so that the in-vivo correlation is never
standing where a causal claim belongs.

1. **WHAT THE NORMAL GLAND DOES TO ITSELF.** Between 6W and 12W the wild-type mammary epithelium
   withdraws from respiration and from nothing else that would explain it: OXPHOS subunits
   **-0.255**, at the **0th percentile** against expression-matched random sets, while the
   proliferative programme barely moves (`PROLIF_*` pooled **-0.047**, 1.3rd) and the assembly
   factors of the very same complexes sit at the **50th**. The paired null on the *difference* is at
   the 0th percentile — the claim is comparative and is tested comparatively. **The gland
   de-respires without de-proliferating.** Its physiological reading is morphogenetic: the TEB
   regresses (-0.422), and the TEB is the apoptosis-loaded compartment, because duct cavitation
   needs apoptosis. What may carry the respiratory arm is the lineage-suppressed (LE) progenitor
   state, which declines in the AP and BA compartments (section 0.7a of the alignment doc).
2. **AND IT DOES NOT DISMANTLE THE DEATH MACHINERY.** Every pro- and anti-apoptotic transcript sits
   where it was at 6W — of 25 PRO + 7 ANTI genes exactly one moves, and it moves *up*. Nor does the
   gland buffer: `Bcl2` padj 0.45, `Bcl2l1` 0.79, `Mcl1` 0.44. **What is withdrawn is the state the
   machinery depends on, not the machinery.** This is what makes the next step a *competence*
   argument rather than an abundance one.
3. **MYC ARRIVES AND BUILDS A DEATH-READY ORGANELLE ON THAT SUBSTRATE** — content +21-27%, all five
   complexes, import and central dogma, at both ages, plus HTRA2 (+1.49, padj 2.6e-08), BAX up,
   BCL-xL down. Asset and liability, constructed together, which is the trade-off stated as
   molecules.
4. **WHAT MYC LOSES BY 12W IS THE TRIGGER, NOT THE MACHINERY.** `PUMA:Bcl-xL` priming retains
   **-0.09** against a global rate of 0.55 (92nd percentile of expression-matched pairs; purity-
   adjusted interaction p 0.017), while `BAX:Bcl-xL` retains exactly 0.55. Independently, on a
   genome-wide statistic that never touches the Bcl-xL denominator, `Bbc3` sits at the **0.51st
   percentile of 8774 Myc-responsive genes** — pre-specified from the PGC1a westerns, not a scan hit
   (`[44]`). Note what this is: **the collapse has no wild-type definition** (the wild-type
   `Bbc3:Bcl2l1` ratio drifts *up*, +0.168). It is a loss of *inducibility* — developmental in
   cause, Myc-dependent in measurement.
5. **DOSE VERSUS COMPETENCE, and why both mouse experiments are needed.** MMTV-Myc 6->12W is a
   **dose** experiment: MYC protein -50% at constant transcript gives a uniform x0.55 rescaling of
   the whole programme (scripts 40, 41, figS8), which retrospectively explains every negative in the
   corpus. MYC-ER is a **competence** experiment: same or higher MYC, **80% less death against only
   30% less proliferation**, in animals that had never seen MYC — so the change cannot be selection
   and must be developmental. Both inputs fall in MMTV-Myc; only the substrate falls in MYC-ER. The
   asymmetry is what a two-input requirement predicts and a one-input model cannot produce.
6. **HYPOTHESIS, stated as generated and not proven.** The mitochondrial OXPHOS/PGC1a state is the
   **second input**: MYC supplies the oncogenic signal, the respiratory substrate supplies apoptotic
   competence, and death requires both.
7. **CAUSAL TEST, in cells, both directions.** PGC1a raises OXPHOS and PUMA and does not kill on its
   own; PGC1a plus MYC kills more than MYC alone; the same in MYAZ, where **PGC1a alone kills and
   PGC1a + Bcl-xL grows better than Bcl-xL alone** — a sign reversal of one intervention's fitness
   effect, conditional on apoptotic competence. Conversely, long-term MYC selects SS cells that lose
   OXPHOS and resist, and passaged PGC1a cells lose the biogenesis and with it the sensitivity, with
   PGC1a still overexpressed — so it is the *output*, not the coactivator, that was selected
   against. Bcl-xL rescue localises the step to **execution**: the programme is still induced and
   the cells still do not die.
8. **THE WINDOW.** Two hits on the death arm (less MYC, and no trigger even at full MYC), one
   partial hit on proliferation. Expansion outruns death and tumourigenesis takes off.
   **What the cells achieve by selection, the gland achieves by development** — a transient,
   mutation-free phenocopy of the second hit every other MYC model needs as a lesion.
9. **CLOSING.** What changes across progression is not what the tumour *needs* from its
   mitochondria, but what it can *afford*.

Points 8 and 9 are the sentences the paper is for, and 7 is where the in-vitro selection experiments
find their place: SS cells and passaged PGC1a cells are the *tumour-progression* model, not the
6W-to-12W substrate model. Conflating the two would be the easiest mistake to make here, because
MYC-ER rules selection out for the developmental window — the gland gets there without it.

**Figure consequence.** Substrate-first promotes **script 43 PART A to an early main figure** — the
de-respires-without-de-proliferating comparator with its matched and paired nulls now opens the
argument rather than arriving late in support of it. The mito-content panels (fig01/fig01b) move to
step 3, and the dose material (fig03/figS8) to step 5.

**Why not the other two orders.** *Myc-first* (the previous version of this section, which opened on
"MYC drives mitochondrial biogenesis") makes the substrate arrive as a late qualifier, and the reader
has already formed the impression that the paper is about what Myc does. *Death-phenotype-first* —
leading with the killing and then showing it correlates with OXPHOS de-prioritisation — says the same
thing but puts the correlation before the perturbation and invites the reader to weigh it as
evidence. Neither does what the order above does.

---

## 10. Provenance

Transcriptome numbers are from `results/priming_arm_teb.rds`
(`scripts/42_priming_arm_and_teb_substrate.R`, author-run 2026-07-25). Underlying sources:
`results/interaction_results.rds` (raw/unshrunken LFC, genome-wide BH padj),
`results/dds_int_run.rds` (median-of-ratios normalised counts),
`results/mitopps_scores.rds` (per-sample pairwise-ratio scores),
`data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt` and
`data/genesets_from_library/gray_chea_mito_tf_shortlist.csv`. Blot values (MYC -50%; BAX/BIM/PUMA
and the mitochondrial panel) are the author's, reported here as given.

Strand 1's tier tables come from `results/background_vs_myc.rds$ruler` (script 40) and strand 3
from `results/mitopps_priming_pgc1a.rds$priming_couplings` (script 38 PART B, run 2026-07-19) —
both already on disk; script 42 PARTS G and H re-tabulate them alongside the new tests.

Related: `docs/2026-07-24_background_vs_myc_interpretation.md` (script 40, the x0.55),
`docs/2026-07-25_attenuation_competing_hypotheses.md` (how the question was bounded before the
blots), `outputs/myc_stability_panel/README.md` (script 41),
`scripts/34_death_priming_reassessment.R` (why composite death scores are circular and why this
analysis works gene-by-gene instead).
