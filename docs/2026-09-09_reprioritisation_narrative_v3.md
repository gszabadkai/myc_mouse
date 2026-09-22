# Narrative v3 — "Developmental reprioritisation of the mitochondrial transcriptome gates MYC-driven mammary tumourigenesis"

**Supersedes `2026-09-09_reprioritisation_narrative_v2.md`.** Movement 2 is
rewritten around the LbNOX result, which names the lethal variable and adds a
functional demonstration of the MYC interaction. Movements 5 and 7 change as a
consequence.

⬚ marks numbers I have from conversation rather than from a document.
N3 throughout.

---

## 0. The spine

Development reprioritises the mitochondrial transcriptome, withdrawing the
respiratory chain from a compartment that is not itself shrinking. **The lethal
variable is mitochondrial NADH oxidation on a MYC background** — not
mitochondrial mass, not ATP, not membrane potential — so the developmental
withdrawal opens a window in which MYC-driven transformation is survivable.
Tumours that later combine MYC with high respiration exist only because they
lowered the oxidative driver or blocked the effector, and their transcriptomes
record which. In human breast tumours the respiratory axis, not MYC, orders the
apoptotic machinery: the post-gate state of the same relationship.

**The novelty claim:** respiratory *priority* within the compartment, rather than
compartment size, gates oncogene-induced death — and the proximate variable is
the matrix NADH/NAD⁺ ratio.

---

## 1. Movement 1 — The gland reprioritises its mitochondria, and MYC does not cause it

*(Unchanged from v2.)*

Signed temporal LFC falls in both genotypes (`OXPHOS_SUBUNITS` WT −0.25, Myc+
−0.39) while the biosynthetic arm rises in WT (amino acid +0.18, lipid +0.12).
Against expression-matched controls the WT window gives proliferation −0.047
(1.3rd percentile) against OXPHOS subunits −0.255 (**0th percentile**) — both
real, OXPHOS fivefold larger. That pair opens the paper because it pre-empts
"rediscovered puberty."

*Reprioritisation* is load-bearing and exact: at 12W_pos nuclear OXPHOS is
mitoPPS-deprioritised (0.955) while its absolute level sits at +0.05. In absolute
mRNA MYC *raises* nuclear OXPHOS (d ≈ 1.27, p ≈ 0.004). Never write that MYC
reduces OXPHOS.

Nuclear OXPHOS peaks at 6W_pos (mitoPPS z +1.21) and falls by 12W while
mtDNA-encoded subunits run the other way (−1.22 → +0.37). **Complex I strongest**
— which now matters more than it did, since complex I is the NADH dehydrogenase
arm and Movement 2 names NADH oxidation as the lethal variable. The developmental
axis and the manipulation series converge on the same node.

Mass separates from function: OXPHOS (~188 genes) and biogenesis (~249) overlap in
17; after the shared axis is removed OXPHOS retains coupling to outcomes and
biogenesis is a MYC-dose bystander.

**Limits.** Absolute age decline directional only (WT β −0.45 p 0.21; Myc+ −0.73
p 0.11). Depth confound disclosed (6W 18.8–29.0M vs 12W 9.9–16.1M; cancelled by
construction for shares and size factors). The within-WT regression on the
developmental luminal axis (partial r 0.71) and ER/PGC1α TF activity (0.70–0.83,
R² 0.86–0.89) is hypothesis-generating, and conditioning on those axes does not
absorb the MYC-effect attenuation (−0.279 → −0.322). `Cdkn2a` flat ⬚ — a stated
negative that protects the novelty claim.

---

## 2. Movement 2 — The lethal variable is mitochondrial NADH oxidation, and it requires MYC

**Rewritten. This is the novelty pillar and it now names a mechanism rather than
excluding one.**

| manipulation | compartment mass | ATP-coupled | NADH/NAD⁺ | ΔΨm | kills at MYAZ baseline | kills with extra MYC |
|---|---|---|---|---|---|---|
| PGC1α | **expands** | yes | ↓ ⬚ | ↑ ⬚ | **yes** | yes |
| Ndi1 | unchanged | yes, downstream | **↓ (measured)** | **↓** | no ⬚ | **yes** |
| **mtLbNOX** | unchanged | **no** | ↓ ⬚ | ↓ ⬚ | **no** | **yes** |
| cytLbNOX | unchanged | no | ↓ cytosolic ⬚ | — | no | **much less** ⬚ |

**Three exclusions, each from a manipulation that kills.**

*Not mass* — Ndi1 and mtLbNOX expand nothing and kill. *Not ATP* — LbNOX is a
water-forming oxidase that makes none, and kills. *Not membrane potential* — Ndi1
lowers ΔΨm, mtLbNOX diverts electrons from the pumping chain and lowers it too,
PGC1α raises it, and all three kill.

**One inclusion, and one compartment.** All four manipulations oxidise NADH. Only
the matrix-targeted ones kill efficiently: **mtLbNOX kills, cytLbNOX much less.**
So the variable is the **mitochondrial matrix NADH/NAD⁺ ratio**, not whole-cell
redox.

That compartment result also eliminates the leading transcriptional candidate.
CtBP1/2 are direct NADH sensors that repress a pro-apoptotic set including PUMA
and NOXA, but they respond to free *cytosolic* NADH — and cytLbNOX is the arm
that barely kills. **CtBP is excluded by the data**, which is worth stating,
because it is the mechanism a reader who knows the redox literature will reach for
first. SIRT3 is the compartment-matched candidate that remains; the relay from a
matrix NAD⁺ signal to nuclear apoptotic transcription is unspecified and should be
presented as such.

**And here is the interaction, demonstrated functionally.** Ndi1 and mtLbNOX do
not kill MYAZ at baseline — MYAZ already carries MMTV-Myc — but both kill
MYAZ + Tre-MYC. **The respiratory manipulation is lethal only above a MYC
threshold.** That is the `MYC × OXPHOS` interaction the gate model requires,
shown by intervention in the form the human cross-section cannot show, and
reproduced in a second background: Ndi1 potentiates MYC-induced death in
iMMEC + Tre-MYC as PGC1α does.

PGC1α is the exception — it kills at MYAZ baseline. Either it is a stronger
manipulation of the same variable or it carries additional death inputs (ROS,
lipid, compartment expansion), and the current data cannot separate those. **Say
so.** It is also the reason PGC1α alone cannot establish the mechanism and the
mimics can.

⬚ **Needs pinning before this is written:** whether NADH/NAD⁺ was actually
measured in the mtLbNOX and PGC1α lines or only in Ndi1 (the manuscript documents
it for Ndi1); a number for "much less" in cytLbNOX; and whether mitochondrial
content was measured directly in the Ndi1 and LbNOX lines (mtDNA copy number,
TOM20, citrate synthase), since the "mass unchanged" column is currently inference
from what the constructs are.

**The dissociation this creates, and keep it explicit.** Death sensitivity is set
by matrix NADH oxidation and is ATP-independent. The *growth* advantage requires
ATP-coupled respiration — cyto- and mitoLbNOX do not increase in vivo growth,
Ndi1 + Bcl-xL does. **Two variables, two phenotypes**, which is why Bcl-xL plus
Ndi1 grows faster and Bcl-xL plus LbNOX would not be expected to.

---

## 3. Movement 3 — PGC1α induces PUMA independently of p53

*(Unchanged from v2.)*

iMMEC carries p53DN, E1A and SV40-LT; MYAZ is p53 proficient. PGC1α induces PUMA
in both, higher basal in MYAZ ⬚. The canonical route is one path, not the
mechanism — which narrows the candidate inducers to the p53-independent set.

**FOXO3, prior first.** FOXO3 is PUMA's canonical p53-independent activator with
direct ChIP evidence; PGC1α is its coactivator; the induction is p53-independent
in our own hands. *Therefore* we asked whether `Foxo3` behaves as the model
predicts across the window, and it does — rising in the maturing WT gland
(+0.473, padj 0.0078), falling along the MYC timeline, interaction −0.486,
p 0.0074. Written this way the scan rank is unnecessary and the ~44-genes
objection does not apply.

State that `Foxo3` was not pre-specified and that the FOXO3-repression arm did not
transfer to human tumours.

**Note the layer gap explicitly:** sirtuins act on FOXO3 by deacetylation, while
the timeline measures transcript. The NAD⁺ arm of Movement 2 and the `Foxo3`
transcript behaviour are **two leads that converge on one protein**, not one
chain. Do not write them as a chain.

**What would make it a pillar:** FOXO3 knockdown in iMMEC or MYAZ, asking whether
it blocks PGC1α-induced PUMA.

---

## 4. Movement 4 — The apoptotic transcriptome reorganises with MYC across the window

*(Unchanged from v2.)*

The wild-type negative with its power control: pro +0.025 vs anti +0.045
(difference −0.020) across the window, against MYC moving the arms apart at 6W
(+0.156 vs −0.024, **+0.180**) — and the genotype contrast reaching padj < 0.05
for **6 of 35** genes against **1** across the window.

`Bbc3` crosses zero: +0.674 → −0.061 against a ×0.55 rescaling line, retention
−0.09, the only ratio in the panel to cross; constant in WT, −0.48 (padj 0.016)
under MYC, interaction p 0.0081, 0.51st percentile of a scan across 8,774 genes.

Specificity with exceptions: `Bcl2l1` at the origin; `Bcl2l11` in the ordinary
middle (91.6th percentile, p 0.86) as the pre-registered pair split; `Bmf` rising
in both timelines as internal control; `Bax` also falling (−0.376, padj 0.008), so
"specific" holds among the BH3-only sensors.

Cross-sectionally: `Bcl2l1` withdrawn (`bMX` −0.670, 2.1st percentile), `Mcl1`
retained (+0.806, 98.9th), `Mcl1:Bcl2l1` +0.994 (p 0.00014; within-timepoint
+0.351, p 0.00058). M\* ≈ +0.007 — respiration does nothing without MYC, which is
now the transcriptional echo of Movement 2's functional threshold.

**Limits, once and plainly.** Nothing survives genome-wide multiplicity (ratio BH
0.33, `Bbc3` BH 0.84); the licence is pre-specification. The coupling is p 0.0052
but empirical p 0.083 on permutation, 0.088 with a timepoint term.

---

## 5. Movement 5 — Escape is either lowering the driver or blocking the effector

**Now with a named thing to escape from.**

PGC1α kills MYAZ cells, so every PGC1α-expressing line is a survivor population
before it is ever injected. Four arms, each with a matched control:

| arm | escape route | `Ppargc1a` CPM | OXPHOS | growth | `Myc` vs control |
|---|---|---|---|---|---|
| MYAZ + PGC1α *(no RNA-seq)* | **driver silenced** | — | none | = control | — |
| p19^ARF KO + PGC1α | **driver lost entirely** | 0.2 | ? | = NT; KO advantage cancelled | 2176 vs 2368 |
| NT + PGC1α | **driver reduced to a partial dose** | 128.6 | slight ↑ | = NT-EV | 2286 vs 2208 |
| **Bcl-xL + PGC1α** | **effector blocked; full dose kept** | **301.6** | **large ↑** | **> Bcl-xL EV** | **1440 vs 2033** |

**Two escape classes, both observed.** Three arms lowered the oxidative driver.
One blocked the effector downstream and kept the driver — and it is the only arm
that grows faster and the only arm where `Myc` comes down.

The buffered arm shows a *second* effector block on top of Bcl-xL: p19^ARF and
p21 down by western, with `Bbc3` −0.727 in the same tumours. A canonical p53
target falling with a dismantled p53 axis, not a PGC1α regulatory effect.
Prediction testable in existing counts: `Cdkn2a`, `Cdkn1a`, `Mdm2`, `Ccng1`,
`Zmat3` down together. The iMMEC survivors take a different route (PUMA slightly
down by WB) in a p53-null background — same argument, different route.

**Which unifies the paper.** The gate is matrix NADH oxidation on a MYC
background. It is passed by lowering the driver or blocking the effector.
**Development does the former by programme; tumours do either by selection.**
Which effector gets broken is arm-specific and incidental — that is why CDKN2A is
an instance rather than the story.

**What the orthotopic cannot say:** anything resting on `Bcl2l1` in the 14
construct-carrying samples; the guardian ratio or gap in any form; a clean
PGC1α → PUMA regulatory statement.

---

## 6. Movement 6 — In human tumours the respiratory axis orders the death machinery

*(Unchanged from v2.)*

| | TCGA | SCAN-B |
|---|---|---|
| OXPHOS split | 0.453 | 0.489 |
| ...MYC removed | **0.485** | **0.525** |
| MYC split | 0.187 | 0.137 |
| ...OXPHOS removed | **−0.043** | **−0.058** |

Spread 1.6× wider on OXPHOS (SD ratio 1.640 [1.473, 1.784]; 1.558 [1.461,
1.650]); all five contrasts exclude their null in both cohorts; per-gene
replication 0.87–0.92. The only rung whose falsifier was written before the answer
was seen. Control: MYC still tracks the mitoribosome at z 8.66 / 6.54 while the
death machinery sits at 0.67 / 0.13.

Configuration: `BAD` +0.50, `BCL2L1` +0.39, `BIK` +0.28, `BBC3` +0.27, `BID`
+0.21 up; `MCL1` −0.27, `BMF` −0.21, `BCL2L11` −0.17, `BCL2` −0.12, `PMAIP1`
−0.11 down. Plus the apoptosis-specific negative: cytosolic effectors run against
the mitochondrial axis (−0.106 / −0.094, **0th percentile of 30 pathways**).

**Forbidden:** balance or ratio language (grid 92–95% additive; 0 of 35 ratios
gain on OXPHOS in both cohorts); "pro-apoptotics rise as a block"; any
`MYC × OXPHOS` interaction phrasing on the human side. Never set a mouse `bMX`
beside a human partial Spearman — different estimands, orthogonal by
construction. Audit item across the whole draft.

---

## 7. Movement 7 — The interaction is real, demonstrated in vitro, and absent from human tumours for the reason the model predicts

**Now an affirmative argument, not only a defensive one.**

The human `MYC × OXPHOS` interaction broke under three of four falsifiers and the
sibling pre-registered study found the functional version null. **Nothing here
overturns that.**

What changes is that the interaction no longer rests on the human cross-section at
all. **Movement 2 demonstrates it directly:** Ndi1 and mtLbNOX kill only above a
MYC threshold, in two backgrounds. That is the interaction, measured functionally,
in the pre-transformation state where the model locates it.

The human null is then exactly what the model predicts. A gate is a property of
the window before transformation; every human tumour has passed it; a
cross-section of survivors cannot show a gate. Movement 5's four-arm table is that
shown by selection.

Position this in the main text so a reviewer meets the argument before the null.

---

## 8. Figures

| | content | source |
|---|---|---|
| **1** | Developmental reprioritisation: OXPHOS vs proliferation against matched nulls; mitoPPS four groups with nuclear/mtDNA crossover; mitoPPS-vs-absolute at 12W_pos; within-WT regression | timeline |
| **2** | **The lethal variable**: the four-manipulation grid (mass / ATP / NADH / ΔΨm / death); mtLbNOX vs cytLbNOX; the MYC-threshold requirement in MYAZ and iMMEC; Bcl-xL + Ndi1 growth rescue | in vitro ⬚ |
| **3** | MYC-conditional apoptotic reorganisation: WT negative with power control; two-timeline plane; `bMX` forest with expression-matched percentiles | timeline |
| **4** | p53-independent PUMA: iMMEC vs MYAZ; FOXO3 across the window | in vitro ⬚ + timeline |
| **5** | The escape series: four-arm driver/OXPHOS/growth/`Myc` table; p53 target panel; construct-free NT contrast | orthotopic + WB |
| **6** | The human axis: conditioning asymmetry with mitoribosome control; the twelve with intervals; regulon-up / cytosolic-down against 30 nulls | TCGA + SCAN-B |

Extended Data: per-complex nuclear/mtDNA split; absolute-vs-mitoPPS across three
normalisations; depth disclosure; ratio grid additivity and zero-gain; MOM/MOM
null; the four falsifiers; orthotopic composition gate with positive control; C2
normalisation dependence; the Bcl-xS limitation.

---

### 8b. Alternative allocation — four figures, Letter format

**Not a replacement for §8. An alternative to evaluate against it.** The case for
Letter is not only that the molecular mechanism is incomplete — it is that the
four-figure constraint cuts precisely the material with known weaknesses, and
forces the single claim that seven movements avoid making.

Pushed to Extended Data under this scheme: the mouse apoptotic transcriptome arm
(BH 0.33 and 0.84, empirical p 0.083, pre-specification-licensed), the human
cohorts (post-gate, no interaction, ratio language forbidden, estimand audit
pending), and FOXO3 (not pre-specified, did not transfer). Each is a place a
reviewer would push; none is load-bearing for the title.

| | content | why essential |
|---|---|---|
| **1** | Developmental reprioritisation. OXPHOS vs proliferation against matched nulls (0th vs 1.3rd percentile); mitoPPS across four groups with the nuclear/mtDNA crossover; **mitoPPS 0.955 against absolute +0.05 at 12W_pos** | The title's first clause. The third panel *is* the title. |
| **2** | The lethal variable. Four-manipulation grid (mass / ATP / ΔΨm / NADH / death); mtLbNOX vs cytLbNOX; the MYC threshold in MYAZ and iMMEC | Three exclusions, one inclusion. Carries the novelty and the interaction. |
| **3** | The gate. MYC-ER transformation at 12W but not 6W, with the respiratory state of the two windows | The verb. Without it the title says "correlates with". |
| **4** | Escape and relevance. Four-arm driver/OXPHOS/growth/`Myc` table; Bcl-xL + Ndi1 rescue; **one human panel** — the conditioning asymmetry with its mitoribosome control | Tumours confirm the gate by escaping it; human relevance in one bar chart. |

Reads as: development does this → this is why it matters → this is the
consequence → this is what tumours do about it. Nothing in the chain is nominal.

Extended Data lands at roughly ten items — per-complex nuclear/mtDNA;
absolute-vs-mitoPPS across three normalisations plus depth; the within-WT
regression; the apoptotic transcriptome arm; the `bMX` forest; FOXO3 and
p53-independence; the human twelve plus the cytosolic negative; ratio-grid
additivity and the falsifiers; orthotopic composition gate and C2 normalisation;
full respirometry and NADH/NAD⁺. At the ceiling, so anything added later displaces
something.

**The binding constraint is words, not figures.** ~2,500 for a Nature-family
Letter, carrying the reprioritisation-versus-content distinction, three exclusions
and one inclusion, the MYC threshold, the gate, two escape classes and human
relevance. The reprioritisation subtlety is the expensive one: "priority within
the compartment, not compartment size" needs room or a reviewer reads it as
content, and that misreading collapses the novelty claim. Budget that paragraph
first.

**Risk to watch:** human data in Extended Data may read as a mouse paper to a
metabolism editor. The Fig 4 panel plus an abstract clause should hold it; if an
editor pushes back, that is the signal to go long rather than to argue. Expanding
from a Letter with full ED to Nat Comms discards nothing.

---

## 9. Computation queue

1. **p53 target panel across all orthotopic arms** — tests Movement 5's mechanism.
2. **`NTPgc1a` vs `NTEV`** — respiratory rulers, mitoPPS, construct-free twelve,
   `Foxo3`. The only construct-free PGC1α contrast.
3. **Four-arm `Myc` and OXPHOS table**, all 49 vector + CRISPR samples in one scoring run.
4. **`Cdkn2a` in the timeline** — the protective negative.
5. **`Sirt3`, `Ctbp1`, `Ctbp2`, `Nampt`, `Nmnat1`, `Nmnat3`** across timeline and
   orthotopic. Won't test activity, but a coordinated shift in NAD⁺ synthesis
   capacity across the window would be worth seeing given Movement 2.
6. **Whether the in vitro PUMA measurement was on the injected population.**

Not queued, deliberately: CDKN2A in human; any recovery of C3; any mouse–human
comparison before the estimand audit.
