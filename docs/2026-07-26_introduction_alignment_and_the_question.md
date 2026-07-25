# The Introduction, the question, and the trade-off the paper is actually about

**Status:** alignment note, written 2026-07-26, after the MYC/apoptotic blots, scripts 40-42, and
the MYAZ implant result. Supersedes section 3-4 of
`docs/2026-07-13_BlockA_revision_walkthrough_and_intro_alignment.md` (the clause table and the
proposed rewrite at `:884-1066`). Companion to
`docs/2026-07-25_death_narrative_dose_vs_competence.md`, whose section 9 is the agreed argument
order; this note is the layer above it — what the paper is asking, and how the Introduction has to
be re-pointed to ask it.

**Scope decisions taken with the author:** all three stages stay in the paper (initiation, primary
growth, metastasis — the metastasis arm has data not represented in this repo); the Introduction
**names the axis, not the genes** (no PGC1a, no PUMA); the metastasis unification is the paper's
**closing claim**; this note supplies an evaluation plus a complete drafted Introduction, not an
abstract.

**Numbers flagged `[43]` are provisional** — they come from a read-only reconciliation on
2026-07-26 and from a scratchpad dry run of `scripts/43_substrate_specificity_and_tradeoff.R` at
reduced draw counts. They become citable when the author runs 43 at its committed settings.

---

## 0. The question

The current Introduction contains three candidate questions and commits to the weakest.

- **Paragraph 1** asks *how are metabolic priorities remodelled across the stages of progression?*
  Descriptive, and it stakes the novelty on longitudinal-versus-cross-sectional.
- **Paragraph 2** asks *how does MYC engage mitochondrial programmes in vivo?* Also descriptive,
  and largely answered in humans by Menegollo et al. 2024.
- **Paragraph 3** asks the real one, and then abandons it: *MYC drives biogenesis and apoptotic
  sensitivity through the same organelle; how is that conflict resolved so transformation can
  proceed?*

Paragraph 3 is the right question, and the MYAZ implant result sharpens it into the form the paper
should use. PGC1a in MYAZ cells **kills** them; PGC1a **plus Bcl-xL grows better than Bcl-xL
alone**. That is a **sign reversal of the same intervention's effect on fitness, conditional on
apoptotic competence** — so respiration is not simply a liability to be shed. It is simultaneously
the strongest pro-growth and the strongest pro-death input, and which limb is expressed depends
entirely on whether death can be executed. The model is therefore not a permissive gate but a
**coupled trade-off**:

> **Mitochondrial respiratory capacity is both an asset and a liability to a MYC-driven cell: it
> supports growth, and it primes death. When, and how, does an incipient tumour capture the asset
> without paying the liability?**

The earlier form of this question — *what licenses a MYC-driven cell to survive its own oncogene?*
— is subsumed: it asked only about the liability limb, and it could not explain why respiration
would ever come back.

**Three solutions to one trade-off, and the paper has all three.**

| solution | how it works | cost | where in the paper |
|---|---|---|---|
| **Development gives it away free** | the maturing gland withdraws the respiratory arm, and apoptotic competence with it | none to the cell — but transient | MMTV-Myc 6->12W; MYC-ER at fixed dose |
| **Selection finds it, expensively** | long-term MYC selects OXPHOS-low cells; passaged PGC1a cells lose the biogenesis | gives up the **asset** along with the liability | SS cells, passaged PGC1a cells |
| **Buffering solves it properly** | block execution and respiration becomes pure asset | requires a lesion | Bcl-xL rescue in iMMECs; PGC1a + Bcl-xL in MYAZ |

And the closing claim, which turns the two halves of the paper into one mechanism:

> **What changes across progression is not what the tumour needs from its mitochondria, but what it
> can afford.** Early it cannot afford respiration, and development spares it the cost. Late, once
> death is buffered, it can — and respiration is re-engaged for dissemination.

That is the manuscript title (mitochondria *integrating* an oncogenic and a metabolic programme)
stated as a mechanism rather than as an association, and it is the only framing under which the
initiation and metastasis arms are one paper rather than two.

**What each paragraph then has to do.** ¶1 establishes that metabolic state is treated as an output
read at one point, and that an endpoint cannot say *when* a commitment was made or *what else* had
to be true. ¶2 shows, in human data, that the axis is branched and MYC-linked — so the driver alone
does not determine the state. ¶3 names the trade-off, names the blind spot that has hidden it, and
asks the question. ¶4 answers it in the order of section 9, closing on affordability.

---

## 0.5 Why respiration is coupled to death, and what else it buys

*(Answers the author's first question: is it only priming that travels with the respiratory arm?)*

**Why coupled at all.** Two routes are available and **the paper does not have to choose between
them.** The *trigger-supply* route is transcriptional: the biogenesis transcription factors that
build mitochondrial mass also reach pro-apoptotic effectors — the library's own developmental-TF
note records `NRF1 -> Bbc3` among them, in a resource built for an unrelated purpose. The
*execution-threshold* route is structural: more respiratory apparatus is more of the machinery that
performs MOMP, so the same trigger clears a lower bar. The Bcl-xL rescue localises the block to
**execution** — the programme is still induced and the cells still do not die — while the
transcript work supplies the trigger, and both observations are consistent with either route
operating. Say so, rather than over-committing.

**What else respiration buys: growth.** This is the limb the earlier framing missed, and the MYAZ
implants demonstrate it directly, with both comparator arms present. On a death-competent
background PGC1a is lethal; on a death-blocked background it is a growth advantage. **The same
intervention has opposite effects on fitness depending only on whether the cell can die.** No
correlational analysis can match that, and it is why the trade-off, not the gate, is the right
frame.

**Is proliferation supported by respiration in vivo? The two are dissociable — and that is what
the model needs.**

- In the normal gland, the respiratory arm withdraws while the proliferative programme barely
  moves (section 0.7) — de-respiration without de-proliferation.
- MYC-ER loses ~80% of death but only ~30% of proliferation at equal dose. Death needs two inputs,
  proliferation needs roughly one.
- Cross-sectionally the coupling is **not testable**, and script 36 says why: in GSVA space OXPHOS
  *is* the global axis (r=0.968 with the global mean), so there is no separable variance to test.
  This is an honest negative, not a weak positive — do not quote a coupling here.
- `[43]` In mitoPPS ratio space, where the ambient is ~0.4 rather than ~0.8, the asymmetry appears
  on a single axis with a working control: `myc x OXPHOS-priority` is **+6.09 (p 0.005)** on the
  PUMA priming ratio with the `redox` control null (+2.16, p 0.72, 50th permutation percentile),
  and **negative** on both proliferation composites (−3.68, p 0.11; −3.23, p 0.17). A **sign**
  asymmetry is stronger than a magnitude one: high respiratory priority makes MYC better at raising
  the death trigger and no better at raising proliferation.

**And the in-vivo negative on proliferation is what the trade-off predicts.** The asset limb is
only visible when death is blocked. In a gland where death is *not* blocked, the cells with the
highest respiratory priority are the ones dying, and they are not in the library — the survivor
bias already documented in the death narrative (section 7.3). So the growth advantage should be
invisible or inverted in vivo and visible in the Bcl-xL arms, which is exactly the pattern. State
this as an interpretation, not a result.

**A third limb, hypothesis only.** The project's one robust Myc-independent coupling is
**redox -> MB2 tumorigenic fork (global-adjusted −0.63**, same sign across three quantifiers,
`docs/2026-07-18_narrative_synthesis_five_questions.md`): lower antioxidant capacity tracks the
tumorigenic fork. If respiration also supplies ROS, then mutagenesis is a third consequence of the
same asset, with its own cost. Worth one sentence in the Discussion; nothing more.

---

## 0.6 MYC alone, and a second hit that is a state rather than a lesion

*(Answers the author's second question.)*

MYC is rarely studied as a solitary driver because on its own it usually fails: the apoptotic
response it triggers has to be removed first. The field's answer to *what removes it* is always a
**lesion** — p53/ARF loss or BCL-2 in the Eu-Myc lymphoma models, and, in the closest single
precedent, Bcl-xL co-expression in beta-cell MycER, where suppressing apoptosis is what exposes
MYC's oncogenic properties. **Both citations must be checked before use**; the beta-cell precedent
in particular is the direct ancestor of the Bcl-xL rescue in this paper and should be cited as
such, not rediscovered.

The claim here is different in kind: **normal mammary development supplies a transient,
mutation-free phenocopy of that second hit.** The gland withdraws the respiratory arm as
morphogenesis ends, apoptotic competence goes with it, and the same MYC that killed at six weeks
no longer can. Nothing was mutated; nothing was selected — the MYC-ER animals had never seen MYC.

Three consequences worth stating explicitly, because they are what makes the claim useful rather
than merely tidy:

1. **It explains why it is a window and not a switch.** A lesion is permanent; a state is not. The
   permissive period has a beginning and an end set by the tissue's own programme, which is why
   transformation in this model is stage-dependent at all.
2. **It predicts timing.** A stochastic second lesion predicts stochastic latency; a developmental
   state predicts that the permissive period is the same in every animal. That is testable.
3. **It reframes the cell arm.** SS cells and passaged PGC1a cells are the *tumour-progression*
   model — cells finding by selection what the gland handed over for free. Conflating the two would
   be the easiest error here, and MYC-ER is what rules selection out for the developmental window.

---

## 0.7 What the adult gland is actually doing — and it is not de-proliferating

*(Answers the author's third question. The proposed explanation — that OXPHOS falls because the
adult gland has lost the TEB proliferative axis — is **contradicted** by the data, and the paper is
better for it.)*

Wild-type 6->12W, set-mean raw log2 fold change, same code path as script 40 `[43]`:

| arm, wild-type 6->12W | WT | Myc+ | WT as % of Myc+ | vs matched random sets |
|---|---|---|---|---|
| **OXPHOS subunits** | **-0.255** | -0.405 | 63% | **0th percentile** |
| TEB vs ductal (HS) | -0.422 | -0.220 | 191% | 0th percentile |
| OXPHOS (all) | -0.145 | -0.306 | 47% | 0.5th percentile |
| **PROLIF_\* pooled (731 genes)** | **-0.047** | -0.149 | 32% | 1st percentile |
| mitoribosome | -0.024 | -0.252 | 10% | 29th (n.s.) |
| nucleotide metabolism | +0.0002 | -0.183 | 0% | 45th (n.s.) |
| **OXPHOS assembly factors** | **+0.001** | -0.175 | 0% | 49th (n.s.) |
| amino-acid metabolism | +0.184 | -0.057 | — | 100th |
| lipid metabolism | +0.122 | -0.028 | — | 100th |

Four things to read from it.

1. **The gland de-respires without de-proliferating.** The respiratory withdrawal is **5.4x** the
   proliferative one. The proliferative programme does fall (1st percentile against matched genes —
   it is not literally flat) but by a fifth as much.
2. **The two mitochondrial arms most tied to growth do not move at all** — nucleotide metabolism
   (+0.0002) and the mitoribosome (-0.024), both non-significant against matched random sets. If
   this were growth withdrawal, those are the arms that should have led it.
3. **The internal control is inside the same complex.** OXPHOS *assembly factors* do not move
   (+0.001, 49th percentile) while OXPHOS *subunits* fall hardest. So it is not "mitochondria"
   generically, and not a mass effect — it is the structural respiratory chain specifically, which
   is the PGC1a/NRF1/ERRalpha output.
4. **The comparison beats its own null.** `[43]` The paired null — both sets redrawn together, so
   the null is on the *difference*, which is what the claim actually is — puts OXPHOS-minus-
   proliferation at the **0th percentile** (p_emp 0.000), OXPHOS-minus-nucleotide at the 1st, and
   OXPHOS-minus-mitoribosome at the 0th. The dissociation is not an eyeballed contrast.

Individual markers agree: Mki67 -0.04, Ccnb1 +0.15, Plk1 +0.11, Aurka +0.31 in the wild-type
gland; only Top2a (-0.41) and E2f1 (-0.35) fall. `[43]` Within the twelve wild-type mice the
respiratory score and the proliferation score are strongly collinear (r 0.83) and the timepoint
term nonetheless survives adjustment for proliferation (-0.50 -> -0.42) — descriptive only, since
`batch = timepoint`.

**So what is the gland doing?** Losing the **terminal end bud**, not losing proliferation. The TEB
signature is the largest mover in the table (-0.42), and
`docs/library_reference/Gray_et_al_developmental_TFS_selection.md` already records that the
`HS_TEB` context is *"apoptosis-loaded (TEB lumen clearance by apoptosis during pubertal
elongation)"*. The physiological reading is therefore **morphogenetic, not proliferative**: the
pubertal gland keeps a death-ready mitochondrion because duct cavitation requires apoptosis, and
stands the machinery down when morphogenesis is complete. The permissive window is the interval
after the gland stops needing to kill its own cells.

**Why this matters for the model rather than just being tidy.** Had respiration fallen only because
proliferation fell, the two-input model would collapse into "less proliferation, less death", and
the MYC-ER asymmetry — ~80% less death against ~30% less proliferation — would have no explanation.
The dissociation is what licenses the asymmetry.

**And the gland does not solve the trade-off by buffering.** `[43]` No anti-apoptotic gene moves in
the wild-type timeline: Bcl2 (padj 0.45), Bcl2l1 (0.79), Mcl1 (0.44), Bcl2l2 (0.88), Xiap (0.79),
Birc2/3/5 all non-significant. The gland de-prioritises; tumours buffer. Three solutions, and the
normal tissue uses only the first.

**The standing bound.** Every wild-type temporal statement above is exposed to `batch = timepoint`
and is described, not claimed. Its one mitigation is that the mitochondrial protein blots move the
same way, off the RNA batch entirely.

---

## 1. What changed since 2026-07-13 — and one correction that matters more than the rest

### 1.1 The correction is a frame error, not an arithmetic one — and OXPHOS is the arm the normal gland withdraws from

The 2026-07-13 alignment resolved the Introduction's "reduced abundance of OXPHOS complexes" onto
the **Myc-specific** 6->12W decline, on the strength of an Issue #6 decomposition reported as
*"OXPHOS subunits -0.39 = -0.11 shared (28%) + -0.28 Myc-specific (72%) ... **WT OXPHOS barely
declines** ... the biosynthetic arm does converge; OXPHOS is precisely the arm where it does not."*
That was carried into memory and into the proposed rewrite.

**Issue #6's arithmetic is right. The gloss attached to it is wrong, and the Introduction was built
on the gloss.** Read read-only on 2026-07-26 from `results/attenuation_mechanism.rds$decomp_conv`
(script 31, Issue #6), what that decomposition actually decomposes is the **attenuation of the
genotype gap** into a wild-type convergence term and a Myc fade term, on the divergent-gene
universe. For `MITOCARTA_OXPHOS_SUBUNITS` (n=45 divergent): gap 0.674 -> 0.388, attenuation 0.286 =
**-0.126 wild-type convergence + 0.412 Myc fade**, i.e. conv -44%, fade 144%. That is correct and
it stands: **OXPHOS is the arm where the wild-type background does *not* converge toward Myc.**

But "does not converge" is a statement about the **gap**, and it does not license "WT OXPHOS barely
declines". The Myc+ gland sits *above* wild type on OXPHOS; a wild-type gland that falls is moving
*further below* the Myc+ level — away from it. Anti-convergence and a large wild-type decline are
the same observation seen from two sides. And the decline is large:

| `MITOCARTA_OXPHOS_SUBUNITS` | WT 6->12W | Myc+ 6->12W | WT as % of Myc+ |
|---|---|---|---|
| set-mean raw LFC, `attenuation_decomposition.rds$decomp` (script 29, Issue #4) | **-0.254823** | **-0.4047** | **63%** |
| set-mean raw LFC, `background_vs_myc.rds$ruler` (script 40, independent code path) | **-0.254823** | **-0.40468** | **63%** |
| fGSEA NES, script 29 | **-2.670** | **-2.715** | **98%** |

Two independent code paths agree to six decimal places, so this is one measurement, not two.

**Which frame the Introduction needs is not a matter of taste.** The attenuation frame asks whether
the background closes the genotype gap — the right question for Issue #6, which was about why the
Myc effect halves. The substrate frame asks what state the gland is *in* — and that is the question
the two-input model asks, because MYC-ER holds MYC fixed and samples the gland's own state, not the
genotype gap. In the substrate frame the wild-type gland's withdrawal from the respiratory arm is
large, and it is the largest of any mitochondrial programme.

(Nothing else moves: `background-vs-myc-is-additive` records that Issue #6's convergence term is
inflated by a shared-baseline artifact in the *biosynthetic* arms; the divergent arms, OXPHOS among
them, carry an artifact of ~0. The anti-convergence reading is unaffected.)

**And OXPHOS is not merely developmental, it is the *only* mitochondrial arm that is.** On the
temporal NES from script 29's `decomp`, no other arm behaves this way:

| arm | WT 6->12W NES | Myc+ 6->12W NES |
|---|---|---|
| **OXPHOS subunits** | **-2.670** | **-2.715** |
| OXPHOS (all) | -2.096 | -2.711 |
| Complex I | -2.022 | -2.492 |
| MYC-target core (Hallmark V2) | -1.398 | -2.257 |
| nucleotide metabolism | -0.660 | -1.396 |
| mitochondrial ribosome | -0.606 | -2.306 |
| TCA cycle | +0.579 | -1.755 |
| lipid metabolism | +1.492 | -0.683 |
| amino-acid metabolism | +1.835 | +0.708 |

Script 40's tier table says the same thing on both rulers: of the seven Level-1 tiers, **OXPHOS has
the largest shared fraction** (WT/Myc+ = 0.535 on content, 0.517 on priority), while every other
tier's wild-type background is near zero or opposite-signed.

**Consequence for the Introduction.** The clause has to be re-pointed, and it now points somewhere
better than before: from *"the Myc+ gland de-amplifies its own elevated OXPHOS"* (which the blot
has since reduced to a dose artifact) to **"the maturing gland itself withdraws from the
respiratory arm, and it is the only mitochondrial arm it withdraws from — the MYC-expressing gland
merely does the same thing slightly harder."** That is the substrate fall the coincidence model
needs, in animals whose gland is normal, which is exactly the population MYC-ER samples.

### 1.2 The other supersessions

| 2026-07-13 claim | status now |
|---|---|
| "a largely Myc-specific de-amplification" as the centrepiece of the OXPHOS clause | **Superseded twice** — it rests on the gloss corrected in 1.1 (the wild-type gland declines, and by 63% as much), and the Myc-specific residual is itself the x0.55 dose effect (MYC protein -50% at constant transcript) |
| the mitonuclear-imbalance -> priming coupling as the mechanism | **Retracted 2026-07-17** (fails an empirical null; the imbalance has no supported contrast) |
| the optional tail "~2/3 Myc fade + ~1/3 WT convergence" | **Dead** — Issue #6's global WT-convergence is a shared-baseline artifact (`background-vs-myc-is-additive`); on the matched control it collapses to chance |
| "reduced abundance of assembled OXPHOS complexes [COMPANION REF]" | **Now answerable** — the author's mitochondrial protein blots show the decline at protein level, off the RNA batch entirely |
| "the death spine gives it a mechanism" | **Replaced** by the two-input model plus the cell perturbations (`2026-07-25_death_narrative`, sections 2, 3, 6.5) |

---

## 2. Clause-by-clause verdict on the current Introduction

Same SUPPORTED / NUANCE / REWORD idiom as the 2026-07-13 table, so the two can be read together.

### Paragraph 1 (clinical funnel; the temporal axis)

| clause | verdict |
|---|---|
| clinical framing, metabolic adaptation as a target, refs 1-7 | **SUPPORTED** — keep verbatim |
| "the metabolic phenotype of a given tumour remains difficult to predict" | **SUPPORTED**, and it is the best sentence in the paragraph — it is the seed of the second-input idea. Promote it |
| "a further, underexplored layer ... along the temporal axis of oncogenesis" | **NUANCE.** Keep (the metastasis arm earns it), but it is currently doing the work of the *gap*, and as a gap it is weak: it invites the reader to notice we have two timepoints |
| "providing cross-sectional snapshots of what is fundamentally a dynamic process" | **REWORD.** The defensible gap is not that snapshots are static but that **an endpoint cannot say when a metabolic commitment was made, or what else had to be true of the tissue for an oncogene to make it.** Same observation, and it points at the paper's actual claim |

### Paragraph 2 (METABRIC / Menegollo 2024)

| clause | verdict |
|---|---|
| the biclustering result, its cohorts, and its stratification | **SUPPORTED** — the strongest justification in the Introduction; it is the human anchor and it is the authors' own |
| "MYC signalling emerged as one of the principal oncogenic drivers associated with distinct mitochondrial bicluster types" | **SUPPORTED**, and under-used. That MYC associates with **distinct** bicluster types is the point: the axis is **branched**, and the driver alone does not determine which branch |
| "whose engagement of these programmes in vivo we set out to resolve" | **REWORD.** "Engagement" promises description. Replace with the two things a cross-sectional cohort cannot resolve — *when* the branch is taken, and *what else* is required for it to be taken one way |

### Paragraph 3 (MYC's two demands; the barrier)

| clause | verdict |
|---|---|
| MYC deregulation, prognosis, the growth/metabolic programme, refs 19-25 | **SUPPORTED** |
| "the same deregulation also sensitises cells to apoptosis, a barrier that must be overcome" | **SUPPORTED — this is the spine.** Promote it |
| "Mitochondria sit at the centre of both processes ... MYC therefore places competing demands on the same organelle" | **SUPPORTED**, and it is the title thesis. Sharpen "competing demands" to the concrete form: the organelle must be **built for growth and disarmed for survival** |
| "how mammary epithelial cells reconcile them ... remains poorly understood" | **NUANCE.** "Reconcile" presumes the cell solves it. Our answer is that the cell does not — the tissue changes underneath and the conflict dissolves. Keep the sentence but do not pre-commit the agent |
| "Much of what is known comes from ... tumours and cell lines that have already adapted ... or models in which it is pre-empted by loss of p53 or amplification of Bcl-2" | **SUPPORTED and badly under-weighted.** This is the best gap statement in the Introduction and it is currently a throwaway. It is a *systematic blind spot*, not a coverage gap: where the apoptotic barrier is already removed, the mitochondrion can only look like a supplier, never like the gate |
| "we use deregulated MYC-driven mouse models to temporarily resolve the mitochondrial remodelling that supports tumour initiation, tumour development and subsequent metastasis" | **REWORD.** ("temporarily" is a typo for "temporally".) It pivots back to description and drops the question just posed. Replace with the question |

### Paragraph 4 (findings)

| clause | verdict |
|---|---|
| "MYC is reported to stimulate overall mitochondrial biogenesis ... its effect is not uniform" | **SUPPORTED** — content rises +21-27% (script 32) *and* the compartment is reprioritised |
| "the dynamics of mitochondrial pathway reprioritisation is central to how transformation proceeds" | **NUANCE.** "Central" reads causal. It is now defensible in a specific and better form: the reprioritisation is what **removes the second input to death**. Say that, not "central" |
| "MYC plays a part and interacts with the transcriptional trajectory of normal mammary development" | **REWORD.** Our result is closer to **additive superposition** than interaction — MYC rescales without reshaping (`background-vs-myc-is-additive`). "Interacts" invites a reviewer to look for an interaction we explicitly did not find. Use "superimposed on" |
| "leading to a robust increase in mitochondrial biogenesis" | **SUPPORTED** with a measured referent (script 32: mass markers +27%, nuclear MitoCarta +19%, import +33%; a share-of-transcriptome **lower bound**) |
| "key mitochondrial pathways are selective[ly] reprioritised, holding OXPHOS restrained" | **REWORD — the most important fix.** Agent ambiguity: as written it reads as *MYC represses OXPHOS*, which contradicts our own genotype result (MYC **raises** nuclear OXPHOS, d ~1.25). The restraint is (i) relative within the compartment and (ii) **substantially developmental** — the normal gland withdraws from OXPHOS almost as steeply as the Myc+ gland (section 1.1) |
| "Reduced abundance of OXPHOS complexes ..." | **REWORD.** Re-point from the Myc-specific decline (now dose) to the **de-prioritisation that survives on the mitoPPS ratio ruler, where a uniform dose scaling cancels** — OXPHOS the most de-prioritised Level-1 tier on both temporal axes, carried by structural subunits (-0.109 WT / -0.153 Myc+) while assembly factors do not move (-0.002 / -0.029). At transcript level this is a *nuclear-subunit* claim; "abundance of complexes" belongs to the protein blots, which now exist |
| "... is associated with a desensitisation to MYC-induced apoptosis" | **REWORD.** The in-vivo association **fails**: per-mouse, the nuclear-OXPHOS composite correlates with the PUMA arm at +0.00 at 6W (50th percentile of ambient) and -0.40 at 12W (25th) — at or below what an arbitrary gene gives. The link is proven by the cell perturbations. State the in-vivo material as **co-decline**, and let causality come from the perturbations, per the epistemic contract |
| "opening a permissive window for tumourigenesis at the adult developmental stage" | **SUPPORTED** (external IHC + the MYC-ER counts; never depended on this RNA-seq) — and it is the paper's payoff sentence |
| "while restrained OXPHOS favours initial transformation, higher OXPHOS activity is subsequently required to support metastasis" | **SUPPORTED** (metastasis arm, not represented in this repo) — and it is what makes the stage framing earn its place. Keep as the closing inversion |
| the final "Together ..." sentence | **SUPPORTED** — keep essentially verbatim; it already says the right thing |

**Missing entirely, and it is the best thing the paper has:** *development, not adaptation.* MYC-ER
animals were never exposed to MYC before induction, so the substrate change cannot be selection.
Set against paragraph 3's own list of the field's answers (p53 loss, BCL-2 amplification, adapted
lines), that is the claim that separates this paper — and we have **both** halves, selection in the
cell arm and development in the mouse arm.

---

## 3. Emphasis map

With all three stages retained, the stage framing keeps real weight; what changes is which *gap* it
states, and how much of the Introduction the trade-off gets.

| element | now | proposed | why |
|---|---|---|---|
| stage / temporal framing (¶1) | ~40% | ~28% | keeps its place because the metastasis arm delivers on it; its gap changes from "studies are cross-sectional" to "an endpoint cannot say when, or what else was required" |
| METABRIC human anchor (¶2) | ~20% | ~20% | unchanged in size, re-pointed: the axis is **branched** and MYC-linked, so the driver alone does not fix the state |
| the trade-off + the barrier + the blind spot (¶3) | ~25% | **~35%** | the spine; ends on the question |
| findings preview (¶4) | ~15% | ~25% | rebuilt in section 9's order, and it has to carry three things the current draft does not: the specificity of the respiratory withdrawal, the sign reversal in cells, and development-not-adaptation |

The single biggest structural change is that ¶3's closing gap sentence — *in a system that cannot
die, only the asset limb of the trade-off is observable* — now does double duty. It is the field
gap, and it is the explanation of why the trade-off has gone unseen. In the committed draft that
sentence only said the models had adapted; it did not say what that cost the field.

---

## 4. Drafted Introduction

Axis-level naming throughout, per the author's decision — no PGC1a, no PUMA, no gene names. Existing
reference callouts preserved where the sentence survives. Square brackets mark optional clauses.

> Breast cancer remains the primary cause of cancer-related mortality in women, with metastasis
> associated with the worst prognosis.^1-3^ It is now established that bioenergetic and metabolic
> adaptations are fundamental for tumourigenesis and progression, making tumour-specific metabolic
> features an attractive therapeutic target.^4-7^ A plethora of studies have sought to identify what
> drives this highly heterogeneous metabolic adaptation, implicating the intrinsic genetic makeup of
> the cancer cell, the tumour microenvironment, and the systemic metabolic state of the host,^8-15^
> yet the metabolic phenotype of a given tumour remains difficult to predict. Part of the reason may
> be that metabolic state is treated as an output of the tumour genotype, read at a single point:
> most studies compare established tumours with normal tissue, tumour-derived cells of high versus
> low metastatic potential, or tumours under pharmacological challenge,^16-18^ and each of these is a
> snapshot of an endpoint. An endpoint cannot say when a metabolic commitment was made, nor what
> else had to be true of the tissue for an oncogene to make it. Metabolic priorities are in fact
> remodelled along the temporal axis of oncogenesis, reprioritised as disease advances from initial
> transformation through primary tumour growth to metastatic dissemination, and how and why they are
> remodelled remains poorly understood, limiting the development of stage-appropriate metabolic
> therapies.
>
> Among these metabolic pathways, mitochondrial content and function stand out as a particularly
> heterogeneous axis, with direct consequences for tumour bioenergetics and survival. Our previous
> unbiased biclustering of the breast cancer transcriptome, seeded on the mitochondrial proteome,
> across more than 2,000 human breast tumours (METABRIC, with independent confirmation in TCGA and
> Oslo2) identified the largest coordinately regulated transcriptional switches in the disease,
> stratifying tumours by cell-of-origin, proliferative state, histology, and clinical outcome
> (Menegollo et al. 2024). Among these switches, MYC signalling emerged as one of the principal
> oncogenic drivers associated with distinct mitochondrial bicluster types. Two things follow.
> Mitochondrial gene programmes constitute a principal, oncogene-linked axis of breast cancer
> transcriptional state, rather than a passive readout of proliferative rate. And that axis is
> branched: tumours sharing a driver do not share a mitochondrial state, so the driver alone does not
> determine it. What a cross-sectional cohort cannot resolve is when, over the course of the disease,
> the branch is taken, and what — besides MYC — has to be true for it to be taken one way rather than
> the other.
>
> The oncogene c-MYC is deregulated, through gene amplification or constitutive upstream signalling,
> in about half of all cancers, and its deregulation carries a worse prognosis in breast cancer.^19^
> MYC drives a transcriptional programme that promotes cell growth, proliferation, and
> pro-tumorigenic metabolic change, including mitochondrial biogenesis.^14,20-25^ Yet the same
> deregulation also sensitises cells to apoptosis, a barrier that must be overcome for transformation
> to proceed.^26-30^ Mitochondria sit at the centre of both processes: they generate ATP through
> oxidative phosphorylation (OXPHOS) and act as the gatekeepers of intrinsic apoptosis through
> BAX/BAK-mediated cytochrome c release.^31^ The organelle that a MYC-driven cell must build up to
> support growth is therefore the same organelle that executes its death, and mitochondrial
> respiratory capacity becomes simultaneously an asset and a liability. How an incipient tumour
> resolves that trade-off — and when — remains poorly understood, and the systems in which
> MYC-driven metabolic change has mostly been characterised cannot resolve it. Established tumours
> and tumour-derived lines have already adapted to this apoptotic stress, and the models used in
> their place pre-empt it by loss of tumour suppressors such as p53 or amplification of
> BCL-2.^32,33^ This is a systematic blind spot rather than a gap in coverage: in a system that
> cannot die, only the asset limb of the trade-off is observable, and the mitochondrion can appear
> only as a supplier of building blocks, never as the gate. It also settles the answer in advance,
> as something the transformed cell acquires for itself. Here we combine deregulated and acutely
> inducible MYC-driven mouse models with defined manipulation of the mitochondrial biogenic
> programme in mammary epithelial cells to ask when, and how, a MYC-driven cell captures the asset
> without paying the liability — and whether the solution is one the cell finds, or one the tissue
> hands it.
>
> Although MYC is reported to stimulate overall mitochondrial biogenesis, we find that its effect on
> the mitochondrial transcriptome is not uniform, and that the resulting reprioritisation is what
> decides whether transformation proceeds. Superimposed on the transcriptional trajectory of normal
> mammary development, MYC raises mitochondrial content and builds an organelle that is not only
> bioenergetically competent but death-competent, loading the pro-apoptotic machinery while
> releasing the mitochondrial brake. Against this, the respiratory arm is the one mitochondrial
> programme from which the maturing gland itself withdraws — and it does so as terminal end bud
> morphogenesis ends, not as part of a proliferative exit: the structural respiratory subunits fall
> in the normal gland several-fold more steeply than the proliferative programme does, while the
> assembly factors of the same complexes, the mitochondrial ribosome and nucleotide metabolism do
> not move at all. The withdrawal persists on a ratio-based measure in which a change in oncogene
> dose cancels. Apoptotic competence is lost over the same window while the remainder of the MYC
> programme is merely scaled down. In mammary epithelial cells the trade-off is then demonstrated
> directly, and it reverses sign: raising the mitochondrial biogenic programme kills MYC-driven
> cells, whereas raising it in cells whose mitochondrial death execution is blocked accelerates
> their growth — the same intervention is lethal or advantageous according only to whether the cell
> can die. [This assigns the mitochondrial biogenic programme a sensitising rather than the
> cytoprotective role usually ascribed to it.] Conversely, prolonged MYC expression selects cells
> that have shed the respiratory programme and with it the sensitivity. In vivo, acute activation of
> MYC in the adult gland at undiminished dose fails to kill, in animals never previously exposed to
> it, so the change is a property of the tissue and not an adaptation of the transformed cell:
> **what tumour cells achieve by selection, the normal gland achieves by development.** Together
> these define a permissive window for tumourigenesis at the adult developmental stage. The
> requirement then inverts with progression, whereas higher OXPHOS activity is subsequently required
> to support metastasis — and on this account it is not the mitochondrial requirement that changes
> across progression but the capacity to afford it, respiration being re-engaged once the death it
> primes can no longer be executed. Our in vitro and in vivo findings therefore define the
> mitochondrial transcriptome as an active and stage-dependent axis of MYC-driven breast cancer,
> reprioritised to permit transformation early and re-engaged to enable progression.

**Per-change rationale, against the section 2 table.**

- ¶1: gap swapped from "snapshots of a dynamic process" to "an endpoint cannot say *when* or *what
  else*". Same observation, but it now points at the second-input claim instead of at our timepoint
  count. "temporarily" -> "temporally" is resolved by rewriting the sentence out of ¶3.
- ¶2: the hand-off becomes the two questions a cohort cannot answer. "Tumours sharing a driver do
  not share a mitochondrial state" is the load-bearing sentence — it is what makes a *second* input
  necessary rather than merely interesting. Check it against the prior paper before submission
  (section 6).
- ¶3: "competing demands" replaced by the **trade-off** — asset and liability in one organelle,
  which is what the MYAZ sign reversal licenses and what makes the metastasis half follow. The
  p53/BCL-2 sentence is promoted from a throwaway to the blind spot, and the reason is now stated
  as *"in a system that cannot die, only the asset limb is observable"* — which is both the sharpest
  form of the gap and an explanation of why the trade-off has not been seen. "Reconcile" is dropped
  so the agent is not pre-committed to the cell; the paragraph ends on the question.
- ¶4: "interacts with" -> "superimposed on" (additive, not interaction). "Holding OXPHOS restrained"
  replaced by a sentence with an explicit agent — *the maturing gland* withdraws — which removes the
  contradiction with our own genotype result, and now carries the **specificity** (assembly factors,
  mitoribosome and nucleotide metabolism do not move; the withdrawal is not a proliferative exit),
  which is section 0.7's result and the thing that licenses the two-input asymmetry. "Associated
  with a desensitisation" becomes "apoptotic competence is lost over the same window", i.e.
  co-decline, with causality carried by the cell sentence that follows — the epistemic contract
  enforced at sentence level. The cell sentence now states the **sign reversal** rather than a
  one-directional sensitisation, because that is the stronger and more surprising fact.
  "In animals never previously exposed to it" is four words that do the entire work of excluding
  selection. The closing sentence adds the affordability clause — the paper's closing claim — while
  keeping the author's own final sentence intact.
- **The protein corroboration is deliberately not in ¶4.** The mitochondrial blots confirm that the
  decline is real off the RNA batch, which is the standing mitigation for `batch = timepoint` — but
  they are a *level* measurement and cannot corroborate a *ratio-based* de-prioritisation, and the
  wild-type limb of the claim needs blots from the wild-type gland specifically. It belongs in the
  Results, stated for the genotypes actually blotted.
- The bracketed clause discharges the canonical-cytoprotection inversion early, at axis level,
  without naming the coactivator. If it is cut, see section 5.

---

## 5. What stays out of the Introduction, and where it goes instead

1. **The dose result — MYC protein -50% at constant transcript, equivalent to the x0.55 rescaling of
   the whole programme.** It is deliberately absent above: an Introduction should not carry a
   control. But it must land **early in the Results**, before any 6->12W comparison, because
   unstated it invites the reader to read the entire attenuation as a substrate story when it is a
   dose story. It probably also needs one clause in the Abstract.
2. **The epistemic contract** (*the in-vivo transcriptome generates the hypothesis, the cell
   perturbations prove it*) belongs at the head of the death section, not the Introduction.
3. **The PGC1a inversion.** If the bracketed clause in ¶4 is cut, the inversion has to be discharged
   explicitly in the Results or early Discussion, with the muscle-atrophy and neurodegeneration
   contexts named — see `2026-07-25_death_narrative` section 8. Do not leave it for a reviewer.
4. **Threshold uncoupling** (apoptosis requires a higher MYC threshold than proliferation) is the
   reviewer's first move against the MMTV arm. It is answered by the inducible model and belongs in
   the Results/Discussion, raised by us rather than by them.
5. **`batch = timepoint`.** ¶4's sentence about the normal gland's withdrawal is a wild-type temporal
   statement, and the wild-type temporal axis is batch-confounded in this design. The Introduction
   may state it plainly; the Results must carry the caveat and the mitigation — the mitochondrial
   and apoptotic protein blots move the same way, off the RNA batch entirely.
6. **Priming is not death.** The transcriptome measures the molecular substrate; the death phenotype
   comes from IHC and the inducible-model counts. ¶4 says "apoptotic competence", not "death", for
   this reason. Keep that discipline in the Results.
7. **The three-solutions table** (section 0) is a Discussion structure, not an Introduction one. The
   Introduction states the trade-off and the question; the Discussion is where development,
   selection and buffering are laid out as three routes to the same end, with the cell arm supplying
   two of them.
8. **The ROS/redox limb** (section 0.5) is a single Discussion sentence at most. It rests on one
   exploratory coupling and would over-extend the model if promoted.

---

## 6. Open items

1. **For the author, on ¶2.** Do the METABRIC mitochondrial biclusters contain both
   mitochondria-high and mitochondria-low MYC-associated classes? If they do, ¶2's claim that
   "tumours sharing a driver do not share a mitochondrial state" can be stated directly from the
   prior paper and becomes the strongest possible hand-off into ¶3. If they do not, soften to the
   fork result we can cite in the mouse (AP7 cross-species: genotype main effect p=0.005;
   MB2_UF-specific Wilcoxon p=0.008). I cannot settle this from the repo.
1b. **The human test of the closing claim, in the authors' own cohort.** If respiration is
   re-engaged once death can no longer be executed, then the mitochondria-high biclusters should
   carry **higher anti-apoptotic buffering** — `BCL2`, `BCL2L1`, `MCL1` — than the mitochondria-low
   ones, and the association should be independent of proliferative index. METABRIC is already in
   hand and the gene sets are in `data/genesets_from_library/metabric_sets.rds`. This is the single
   cheapest test of the paper's closing claim and it needs no new material.
1c. **The two citations in section 0.6 must be checked** before the blind-spot passage is submitted:
   the Eu-Myc/p53-ARF/BCL-2 requirement, and the beta-cell MycER + Bcl-xL precedent. The second is
   the direct ancestor of the Bcl-xL rescue here and should be cited as such.
2. **Resolved, no longer open:** the OXPHOS shared-versus-Myc-specific question (section 1.1). The
   two decompositions were never in conflict — Issue #6 measured convergence of the genotype gap,
   Issue #4 and script 40 measure the wild-type gland's own decline — and the Introduction needs the
   second. No re-run is required; both objects are current.
3. **Propagation.** Section 3.5 (strand 1) of `2026-07-25_death_narrative_dose_vs_competence.md`
   currently makes the developmental case on the tier-median LFCs. The NES pair for
   `MITOCARTA_OXPHOS_SUBUNITS` — **WT -2.670 versus Myc+ -2.715** — is a stronger and simpler
   statement of the same thing and should be added there, with the frame distinction of 1.1 stated
   alongside it so the two readings are never confused again. The memory
   `intro-oxphos-restraint-tension` must be corrected: its rule (a) (never say MYC reduces OXPHOS on
   the genotype axis) stands, but "WT OXPHOS barely declines" does not.

---

## 7. Provenance

Section 1.1 read read-only on 2026-07-26 from three objects: `results/attenuation_mechanism.rds
$decomp_conv` (`scripts/31_attenuation_mechanism.R`, Issue #6 — the gap-convergence decomposition),
`results/attenuation_decomposition.rds$decomp` (`scripts/29_attenuation_decomposition.R`, Issue #4 —
the per-arm temporal LFCs and NES), and `results/background_vs_myc.rds$ruler`
(`scripts/40_background_vs_myc_decomposition.R`, author-run 2026-07-25). All three resolve set
membership through `functions/reconcile_gene_symbols.R`. Tier medians and the priority ruler are
`…$ruler` grouped by Level-1 tier, excluding the flagged mtDNA pathway. Priming, coupling and TEB numbers are
`results/priming_arm_teb.rds` (script 42) as reported in
`docs/2026-07-25_death_narrative_dose_vs_competence.md`. Content figures (+21-27%) are script 32.
The METABRIC fork statistics are `docs/2026-07-06_block_A_review.md` (AP7). The `redox -> MB2 fork`
coupling is `docs/2026-07-18_narrative_synthesis_five_questions.md` (scripts 35/36). The
`NRF1 -> Bbc3` and `HS_TEB` apoptosis-loading statements are
`docs/library_reference/Gray_et_al_developmental_TFS_selection.md`.

**Everything marked `[43]`** — the wild-type comparator table and its matched and paired nulls, the
within-wild-type fit, the buffer padj values, the mitoPPS ambient and the trade-off asymmetry — is
**provisional**. It comes from a read-only reconciliation and from a scratchpad dry run of
`scripts/43_substrate_specificity_and_tradeoff.R` executed at reduced draw counts (200 set draws,
300 permutations against the committed 2000 and 5000), with output redirected away from `results/`.
The committed script writes `results/substrate_specificity_tradeoff.rds`; when the author runs it,
every `[43]` number must be reconciled against that object and the flags removed. Percentiles and
empirical p-values in particular will shift with the larger draw counts.

Blot and cell-culture values — MYC -50%, the mitochondrial and pro-apoptotic panels, the iMMEC and
MYAZ experiments including the PGC1a / Bcl-xL implant arms — are the author's, reported as given.
The prior clause table this note supersedes is
`docs/2026-07-13_BlockA_revision_walkthrough_and_intro_alignment.md:884-1066`.
