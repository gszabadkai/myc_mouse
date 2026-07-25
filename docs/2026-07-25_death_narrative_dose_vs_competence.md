# Dose versus competence: a single line through the death data

**Status:** synthesis note, written 2026-07-25 after the MYC and apoptotic-protein blots.
Supersedes the framing (not the content) of `docs/2026-07-25_attenuation_competing_hypotheses.md`,
which was written before the blots. All transcriptome numbers below are **provisional**: they come
from read-only exploration and are re-derived by `scripts/42_priming_arm_and_teb_substrate.R`
(author-run, writes `results/priming_arm_teb.rds`). Nothing here may be cited until that script
has been run and this document reconciled against it.

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

## 2. This adjudicates the three sub-hypotheses

| hypothesis | verdict |
|---|---|
| (i) MYC engages PGC1a more at 6W than 12W | a **dose** mechanism — cannot explain MYC-ER, where dose is equal or higher |
| (ii) PGC1a is always there; MYC pushes it over a threshold at 6W, and less MYC at 12W falls under it | a **dose** mechanism — **excluded for MYC-ER by construction** |
| (iii) developmental / compositional | **the only one that survives the MYC-ER design** |

(i) and (ii) may still be true *in MMTV-Myc* — indeed (ii) is close to a restatement of what the
dose result says. They simply cannot be the explanation for the MYC-ER phenotype, and it is the
MYC-ER phenotype that forces a substrate story into the paper.

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

| priming ratio (per mouse, log2) | Myc effect 6W | Myc effect 12W | retention | vs expression-matched pairs |
|---|---|---|---|---|
| Bax : Bcl-xL | +0.90 (p 1.2e-4) | +0.50 (p 0.034) | **0.55** | 46th pct — *exactly* the global rate |
| **PUMA : Bcl-xL** | +0.67 (p 0.0089) | **-0.06 (p 0.77)** | **-0.09** | **92.4th pct**, interaction p 0.037 |
| PUMA : Bcl-xL, purity-adjusted | +0.73 (p 0.0007) | | | **interaction p 0.0059** |

Three features make this more than a lone p-value.

1. **It contains its own control.** The panel shares a denominator. BAX priming fades at the
   global rate; PUMA priming collapses. "Isn't this just the x0.55 attenuation?" is answered from
   inside the panel, not by appeal to an external null.
2. **The null reproduces the global rate independently.** Random expression-matched gene pairs
   retain a median of **0.53** — script 40's x0.55, rebuilt from arbitrary genes. The instrument
   is calibrated.
3. **Adjusting for composition strengthens it.** Epithelial and immune covariates move the
   interaction from p 0.037 to **p 0.006-0.017**, while BAX stays non-significant (0.146 ->
   0.185). So it is not contamination — and PUMA's raw correlation with the immune composite
   (+0.34) was the obvious worry.

Supporting: `Bbc3`'s own aligned attenuation (script 34's `sign(lfc6)*(lfc6-lfc12)`) sits at the
**94.6th percentile** of genes matched on both baseMean and |lfc6| (p_emp 0.054) — a second,
ratio-free framing landing in the same place.

**And PUMA was pre-specified.** The PGC1a experiments induce OXPHOS and PUMA "but not other
BH3-only proteins". The one BH3-only the cell work singles out is the one whose in-vivo priming
collapses. This was not a fishing expedition through the BCL-2 family.

**One sentence for the paper: MYC builds a death-competent mitochondrion at both ages; what it
loses by twelve weeks is the trigger.**

That maps onto the cell data with unusual precision. Bcl-xL rescue shows the machinery is still
built — OXPHOS and BIM/BAX go up, more so than without it — and the cells still do not die: the
block is at execution, not at transcription. PGC1a supplies PUMA, and they die.

---

## 4. Why PUMA, and why developmentally

Three routes, in decreasing order of what our data can support.

**4.1 PUMA sits in the TEB-context regulon of MYC and of the biogenesis TFs.** The Gray/ChEA
shortlist (`docs/library_reference/gray_chea_mito_tf_shortlist.csv`) lists `Bbc3` as a target of
**MYC, MYCN, E2F1, KDM5B and MITF — and only in the `AP_TEB` context.** And
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
`MG_TEB_VS_DUCTAL_HS_GRAY_UP` falls **-0.45 (p 0.0028)** purity-adjusted between 6W and 12W; the
BA lane -0.36 (p 0.053); multiple HS_TEB TF lanes decline coherently (ARNTL2 -0.39 p 0.0054,
GRHL3 -0.32 p 0.0038, HMGA2 -0.34 p 0.022). The markers agree (section 2). The HS_TEB lane — the
one the library flags as apoptosis-loaded — is the one that falls hardest.

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

**7.4 The death is invisible to bulk RNA, and that must be said.** There is no efferocytosis or
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

## 9. Proposed narrative order

A reordering of what exists, not new work. The current write-up is organised by *analysis*; the
paper should be organised by *argument*.

1. **MYC drives mitochondrial biogenesis in mammary epithelium** — content +21-27%, all five
   complexes, import and central dogma, at both ages.
2. **Mitochondrial biogenesis sets apoptotic priming** — established by **perturbation, in cells**:
   PGC1a raises OXPHOS and PUMA and does not kill; PGC1a plus MYC kills more than MYC alone; the
   same in MYAZ. And the converse: long-term MYC selects SS cells that lose OXPHOS and resist, and
   passaged PGC1a cells lose the biogenesis and with it the sensitivity. Bcl-xL rescue localises
   the step — the mitochondrion is the execution node, not the transcriptional one.
3. **In vivo, MYC builds the same death-ready mitochondrion** — HTRA2, BAX up, BCL-xL down,
   biogenesis up — at both ages.
4. **What is lost by 12W is the trigger** — PUMA priming collapses while BAX priming fades only at
   the global rate.
5. **Why the programme fades at all: dose** — MYC protein -50% at constant transcript, uniform
   x0.55 rescaling (scripts 40, 41, figS8).
6. **Dose is not the whole story** — MYC-ER: same dose, 80% less death, 30% less proliferation,
   and no prior MYC exposure, so the change is developmental.
7. **The developmental substrate** — TEB regression, and the TEB is the apoptosis-loaded
   compartment.
8. **The window** — two hits on the death arm (less MYC, and no trigger even at full MYC), one
   partial hit on proliferation. Expansion outruns death and tumourigenesis takes off.
   **What the cells achieve by selection, the gland achieves by development.**

Point 8 is the sentence the paper is for, and it is also where the in-vitro selection experiments
(SS cells, passaged PGC1a cells) find their place: they are the *tumour-progression* model, not
the 6W-to-12W substrate model. Conflating the two would be the easiest mistake to make here,
because MYC-ER rules selection out for the developmental window.

---

## 10. Provenance

Transcriptome numbers are **provisional** until `scripts/42_priming_arm_and_teb_substrate.R` is
run; they come from read-only exploration on 2026-07-25 and are re-derived there with their nulls.
Sources: `results/interaction_results.rds` (raw/unshrunken LFC, genome-wide BH padj),
`results/dds_int_run.rds` (median-of-ratios normalised counts),
`data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt`,
`docs/library_reference/gray_chea_mito_tf_shortlist.csv`. Blot values (MYC -50%; BAX/BIM/PUMA and
the mitochondrial panel) are the author's, reported here as given.

Related: `docs/2026-07-24_background_vs_myc_interpretation.md` (script 40, the x0.55),
`docs/2026-07-25_attenuation_competing_hypotheses.md` (how the question was bounded before the
blots), `outputs/myc_stability_panel/README.md` (script 41),
`scripts/34_death_priming_reassessment.R` (why composite death scores are circular and why this
analysis works gene-by-gene instead).
