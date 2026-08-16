# Instrument choice for the reduced figure set

**2026-08-17, branch `paper-final` (Block C).** Which quantifier should carry each beat of the
final message, judged on *clarity and impact* rather than on *what can be best quantified*.

Every number here was read from the panels' own computed `panel_legend()` blocks
(`fig2_priming_ratios`, `fig2_departure_from_dose`, `fig2_oxphos_puma_coupling` re-run
2026-08-17) or from `results/fgsea_percategory.rds`. Nothing is transcribed from memory.

---

## 0. The criterion changed, and the change has a direction

The instrument for each panel was chosen on **"which one can be best quantified."** That
criterion optimises for *defensibility*, and defensibility is bought with **normalisation** —
NES normalises against the ranked list, GSVA against the cohort, mitoPPS against the compartment,
a share against the transcriptome, a fold change against the other group.

The new criterion — clarity and impact — is bought with the opposite thing: **units the reader
already owns.**

So the two criteria do not merely differ, they pull in a predictable direction, and the whole
corpus sits on one side of it:

| rung | instrument | divided by | what the division destroys |
|---|---|---|---|
| **0** | summed DESeq2-normalised counts | library size only | nothing beyond depth |
| **1** | % share of transcriptome | the transcriptome | magnitude |
| **2** | log2 fold change | the other group | level |
| **3** | mitoPPS priority | the rest of the compartment | content |
| **4** | GSVA score | the cohort | absolute position |
| **5** | fGSEA NES | its own ranked list *and* its own permutation null | units, and sign comparability across contrasts |

**The message needs rungs 0–2. The corpus drew rungs 3–5.** That is the single correction this
document proposes, and it is not a loss of rigour — because since scripts 21, 34, 36, 40, 43 and
45, **the defence now lives in the null, not in the normalisation.** A raw fold change with a
matched-random-set null is more defensible *and* more legible than an NES. Both properties moved
the same way; only the panels did not.

> **The ladder has a floor, and script 46 found it** (§4). For a *large* gene set, the absolute
> level is not arm-specific: abundance summed over 87 genes tracks the global expression axis at
> **r = 0.964**, leaving 7% separable variance. Rungs 0–1 are legible but they stop distinguishing
> *this arm* from *everything*. So the optimum is usually **rung 2** — which is where most of the
> paper already sits — and rungs 0–1 belong where the claim genuinely *is* about magnitude
> (Fig. 1E's compartment share) rather than about one arm among others. **Move down one rung, not
> two, and check separability before the last step.**

---

## 1. What each instrument is actually for

| instrument | the one question it answers | the price | what a reviewer does with it |
|---|---|---|---|
| **fGSEA NES** | *"Is this programme involved — judged against the whole transcriptome, with no cut-off I chose?"* | no units; sign flips between contrasts for arithmetic reasons; p anti-conservative for correlated sets; BH within category, so rows are not comparable | **Skims it.** It is insurance against "did you cherry-pick", not persuasion. Every genomics paper has one. |
| **GSVA score** | *"Where does each animal sit on this programme?"* | cohort-relative (adding samples changes it); does not centre within sample, which is *why* the ambient coupling ceiling is 0.80 (script 36) | Reaches for the objection the corpus has already confirmed: **everything correlates with everything.** |
| **log2 fold change** | *"How big, in units I know?"* | a ratio — no level, no magnitude | **Believes it.** Fold change is the field's native unit. |
| **DESeq2 normalised counts** | *"How much of this is there?"* | vulnerable to "that's composition, not biology" | *"Good — so this is real abundance, not a ratio artefact."* Per-animal points also kill the n=6 objection on sight. |
| **% share of transcriptome** | *"How much of the cell is this compartment?"* | a lower bound under global amplification | Accepts it — a percentage is a unit everybody owns. |
| **mitoPPS priority** | *"Is the *balance inside* the compartment shifting?"* | content-blind **by construction**; must be taught to the reader | Has to learn it. Costs a legend paragraph **per use**. |

Two consequences fall out immediately.

**fGSEA earns exactly one panel.** Its unique contribution — "we scanned everything and did not
choose" — is a claim about *method*, and a claim about method needs to be made once. Three NES
panels (1D, 1H, 2E) make it three times and persuade no more than one.

**mitoPPS earns exactly one panel.** Its unique contribution is *reallocation*, and that is
Fig. 1F's whole subject. Everywhere else in the set it is paying the teaching cost without
collecting the benefit — and in one place (Fig. 2I) it is actively working against the sentence
it supports (§4).

---

## 2. Beat by beat

### Beat 1 — "Myc drives an oncogenic programme aligned to a mitochondrial axis"

**Claim type:** a *scan* result plus a *geometry*. Two different things, currently three panels.

| | now | proposed |
|---|---|---|
| unbiasedness | **1D** fGSEA NES, whole library, matched null | keep — this is the one place NES is irreplaceable |
| the axis | **1C** pathway PCA (76%, r = 0.968) | keep — "aligned" is *literally* a geometric word, and PCA is the geometry |
| the magnitude | **1E** share + set-average logFC | **promote.** This is the beat's strongest panel |
| the ranking, again | **1H** NES at 12W vs 6W | **demote to supplementary** (see Beat 2) |

**1E is the panel that should lead.** "Myc raises the mitochondrial share of the transcriptome by
21–27%, per animal, p < 0.001" converts an *alignment* into a *physical quantity*. Alignment is
an abstraction; more mitochondrion is a fact. A reader remembers the second.

> **Wording.** `docs` §37 and the loading audit are explicit that **loading is covariation, not
> primacy**, and that the axis is *mito-led, not mito-specific* (non-mito growth programmes load
> 0.95–0.98). Script 36 goes further: because OXPHOS **is** the global factor (r² = 0.94), "OXPHOS
> is central" is **untestable**, not false. So write **"aligned with"**, never **"primarily acts
> on"** — the latter asserts targeting, which is exactly the thing that cannot be shown.

### Beat 2 — "parallel to the WT developmental programme, partly interacting with it"

**Claim type:** a *geometry* — two vectors and the angle between them. This is the beat where the
instrument choice is most obviously wrong at present, and where the best available idea is already
built but unslotted.

**Draw it on the two-timeline plane** (`two_timeline_base()`): x = `6>12W_wt`, y = `6>12W_myc`,
both raw log2FC, the diagonal = *development alone*. Then:

- **"parallel"** is the distance along the diagonal;
- **"partly interacting"** is the **drop below it** — and that drop **is** the interaction, exactly,
  because `c_tp = c_tn + c_int` to floating point.

One plane, one ruler, both halves of the sentence, no second panel. Nothing else in the set gets
two claims for one axis.

Two numbers now available that make this beat much stronger than when it was written:

- **Spearman(`myc_6W`, `6>12W_wt`) = −0.068 over 866 gene sets.** "Parallel" is now *measured*, not
  asserted. It belongs in the text.
- **The diagonal** (`6W_wt > 12W_myc`, script 45): 96.5% of 143 MitoPathways rise, median +0.247,
  and the respiratory chain is the one arm that gains least. That is the closing conclusion reached
  in a single contrast.

> **Consequence: 2E (fGSEA NES of the wild-type window) can go.** Its claim — TEB lost,
> proliferation stable — sits on the same plane as a labelled point, on the same ruler as
> everything else in the figure.

### Beat 3 — "the apoptotic arm is the primary driver"

**This is the beat most exposed to review, and the fix makes it stronger, not weaker.**

The sentence cannot stand as a *set-level* claim. Script 34 established, and nothing since has
overturned:

- `pro_comp` and the priming composite **are MitoCarta sets**, so every mito↔death coupling in the
  corpus correlated mitochondrial genes with mitochondrial genes;
- the priming **interaction was fitted, saved and never read** — `int_p` 0.177 / 0.263;
- **the attenuation is global**: nothing — death, OXPHOS or MYC-core — attenuates more than
  expression-matched genes (all BH p > 0.14). So there is no death-specific gate to be primary.

**What survives is sharper than what is claimed.** The strongest evidential move in the whole
paper is not a pathway result at all:

> A single global rate was fitted to the transcriptome, and then the question asked was **which
> genes leave it**. Two genes went in pre-specified from the cell work. `Bbc3` came out at the
> **0.51st percentile of 8,774**; `Bcl2l11` at the **91.6th**. `Foxo3` — PUMA's canonical
> p53-independent activator — sits at the **0.46th**, with **only four genes between them.**

That is a **genome-wide test with a named prior**. It is worth more than any enrichment score in
the corpus, and it is a kind of evidence reviewers trust. Sell it as such.

**Replace "the apoptotic arm" with "one apoptotic node, named in advance."** Weaker-sounding;
much stronger evidence; survives review.

**Instrument: per-gene log2FC residual against a global rate (Fig. 2H).** No set method can touch
this claim, and none should be asked to.

> **Promote 2H.** It is currently fourth in its section, behind two panels that set up substrate.
> It is the panel the section should be built around.

### Beats 4 and 5 — PUMA, and why the dose does not explain it

Already correct in instrument. Both are log2FC. Keep both, and keep **2G's rhetorical structure
explicit**, because it is doing something unusual and valuable:

> **Every ratio in 2G shares the same denominator (Bcl-xL).** So "it is just the global
> attenuation" is refuted **from inside the panel**, with no external null required. Bax retains
> 0.52 and Bcl-xL 0.57 against a global 0.49, so their ratio retains 0.55 — the rescaling passes
> straight through. `Bbc3` retains **−1.10** against the same denominator; the ratio retains
> **−0.09**. That single contrast is the argument.

Bounds that must travel: interaction p = 0.037 raw / **0.017 purity-adjusted**, **BH across the
nine pairs 0.33**, and against the matched-pair conditional null **empirical p = 0.082** (92nd
percentile of 510 pairs). And **draw `Bcl2l11`** — the pre-registration had two genes and they
split; drawing only the one that worked would be the wrong panel.

### Beat 6 — the OXPHOS × genotype interaction on PUMA:Bcl-xL

**Claim type:** a per-animal interaction. This *requires* a per-sample score, so it is the one
place in the paper where a sample-level quantifier is unavoidable. See §4 — it needs one change
and it is the highest-value change in this document.

### Beat 7 — the hypothesis

Not a data panel. A **schematic**, and per the standing preference it should carry the *competing*
reading and the experiment that adjudicates it — respiration upstream of PUMA (Dey & Moraes;
FOXO3→BBC3 ChIP) versus PUMA upstream of the pyruvate carrier (Kim, *Cancer Cell* 2019). Nothing
in the transcriptome establishes the direction, and a schematic that pretends otherwise is the
easiest thing in the paper for a reviewer to attack.

---

## 3. One grammar per figure

Counted by *family*, not by panel:

| figure | rulers now | panels |
|---|---|---|
| **Figure 1** | per-sample composite (1B, 1C) · NES (1D, 1H) · share (1E) · logFC (1E, 1F, 1G) · mitoPPS (1F) — **five** | seven, two of which carry two rulers each |
| **Figure 2** | NES (2E) · mitoPPS (2F, 2I) · logFC (2F, 2G, 2H) · per-animal (2I) — **four** | five |

A reader re-learns the y axis five times in Figure 1. That is invisible when panels are judged one
at a time, which is exactly how they were built, and it is the largest clarity cost in the set.

**Target: two rulers per figure, three at the outside.**

- **Figure 1** → *share + log2FC*, with **one** NES panel and **one** mitoPPS panel as the declared
  exceptions, each earning its keep for a reason stated in its legend.
- **Figure 2** → *log2FC* throughout, with the per-animal scatter as the single declared exception
  at the end. Moving 2E onto the two-timeline plane (Beat 2) achieves this on its own.

---

## 4. The highest-value single change: Fig. 2I's x axis

**Now:** the x axis is the animal's **OXPHOS-subunit mitoPPS score** — a *budget share*.

**The problem is not aesthetic, it is that the axis contradicts the sentence it supports.** The
closing hypothesis is *"high respiratory state is required for MYC–PUMA mediated death."* mitoPPS
**cannot** say that, and the panel's own legend already has to warn:

> *"mitoPPS is a RELATIVE score: a high OXPHOS-subunit value means the compartment spends more of
> its budget there, **not that it respires more**."*

**Proposed:** re-run the same interaction with x = **summed DESeq2-normalised counts of the 87
nuclear OXPHOS subunits, per animal** — already computed in `scripts/45_state_readings.R` PART G.
The axis becomes *"respiratory-chain transcript, normalised counts"*: rung 0, a quantity the reader
owns, and the closest transcriptomic proxy to the words in the hypothesis.

**This is a test, not a fix, and it has a real chance of failing.** A ratio cancels the global
common-mode axis; an absolute level does not, and script 36 showed that the OXPHOS composite **is**
that axis (r = 0.968 with the global mean). The absolute version could be dominated by common-mode
variance and give a different answer.

**Decision rule, fixed in advance:**

| outcome | what to do |
|---|---|
| the interaction **reproduces** on absolute counts | draw absolute counts. Clearer axis **and** independent corroboration on a second quantifier |
| it **does not reproduce** | keep mitoPPS, and **report the ruler-dependence in the bounds** |

---

### The test was run. **It rejected this proposal.**

`scripts/46_axis_ruler_test.R`, written 2026-08-17. It reproduces script 43 PART B's model exactly
(max |Δβ| = 0, max |Δp| = 0 over 10 fits) and script 45's OXPHOS levels to 3.55e-15 **before**
changing anything, then refits with only the axis swapped, over the same five outcomes with the same
5,000 within-timepoint permutation null.

| axis | ruler | β | p | p<sub>emp</sub> | ambient | **resid_frac** | verdict |
|---|---|---|---|---|---|---|---|
| `oxphos_ppd` | mitoPPS | **+0.779** | **0.0052** | 0.080 | 0.49 | **0.550** | reference |
| `oxphos_lvl` | levels | +0.535 | 0.156 | 0.154 | 0.81 | **0.071** | **partial** |
| `redox_ppd` | mitoPPS | +0.207 | 0.72 | 0.486 | 0.26 | 0.998 | fails |
| `redox_lvl` | levels | +0.362 | 0.353 | 0.229 | 0.72 | 0.207 | fails |

**The direction survives the ruler change; the specificity does not.** The two OXPHOS rulers
correlate **r = +0.928**, so they are largely the same variable, and the interaction keeps its sign
at about 69% of its size (+0.535 against +0.779). But **only 7.1% of the absolute level's variance is
separable from the global expression factor** (r = **+0.964**), against **55.0%** for mitoPPS, and
the ambient nearly doubles (0.49 → 0.81).

**So on the absolute ruler the axis is "overall expression", not "the respiratory chain."** It could
not have carried an arm-specific claim however the p had come out. **Keep mitoPPS** — and the reason
is now *measured* rather than inherited from script 43.

**Consequence for the text, and it is not cosmetic.** The transcriptomic claim is about respiratory
**priority** (allocation), not respiratory **capacity**. The hypothesis sentence's "high respiratory
state" is a statement about a *rate*, which no ruler here measures — that is the perturbation
experiment's job, and the sentence should hand it over explicitly.

**A by-product worth keeping:** `redox_ppd`'s resid_frac is **0.998** — it is almost perfectly
orthogonal to the global factor, the cleanest axis in the corpus. That is *why* it is the right
control, and it is a better reason than the one currently given.

**The smaller change still stands:** the y axis is the z-scored PUMA:Bcl-xL ratio, and drawing the
raw log2 ratio costs nothing and buys a unit the reader owns.

**One smaller change to the same panel:** the y axis is currently the z-scored PUMA:Bcl-xL ratio.
Drawing the **raw log2 ratio** costs nothing and buys a unit the reader owns.

**And keep redox as the control** — but see §5, because the sentence currently describes it
backwards and throws away its own best argument.

---

## 5. Sentences that do not survive as written

Listed with the replacement, in the order they appear in the message. These go on the `PANELS.md`
sentence list.

**1. *"primarily acts on a mitochondrial transcriptome aligned axis"***
→ **"aligned with"**, not "acts on". Loading is covariation, not primacy; the axis is mito-*led*,
not mito-*specific*; and because OXPHOS is the global factor, centrality is **untestable** rather
than demonstrated.

**2. *"unequivocally identified the apoptotic arm as the primary driver"***
→ Neither word holds. "Unequivocally" — the death composites are MitoCarta-defined and the priming
interaction is null (`int_p` 0.177/0.263). "Primary driver" — **nothing attenuates more than
expression-matched genes** (all BH p > 0.14), so there is no arm to be primary.
**Replacement:** *"one apoptotic node, pre-specified from the cell work, leaves the global dose
line: PUMA:Bcl-xL retains −0.09 where the transcriptome retains 0.49."*

**3. *"apoptotic sensitivity declined ... due to loss of PUMA (Bbc3) in 12W_myc mice"***
→ The claim is the **departure from the dose line**, not a loss at one age. `Bbc3`'s six-week
genotype effect has **padj 0.22** — the half the sentence presupposes is not itself established —
and the sharp statistic is the *ratio*'s reversal, not the gene's level.

**4. *"transcriptional repression of Foxo3"***
→ **Not repression — failure to rise.** The wild-type gland raises `Foxo3` with age (**+0.473**);
the Myc+ gland does not (**−0.013**). The genotype halves are mirror images at **+0.243 (p 0.057)**
and **−0.242 (p 0.059)** — *neither is significant*, which is why the interaction is rank 1 of the
26 mechanism genes, and why **neither half may be labelled on a panel.**
Also: **`Foxo3` was not pre-specified.** It was found in the scan with 44 genes below it. Its
licence is the **prior** — FOXO3 is PUMA's canonical p53-independent activator with direct ChIP
evidence — and the text must introduce it that way round.

**5. *"the interaction ... is highly significant (p = 0.0052)"***
→ **"Highly" cannot stand.** p = 0.0052 is the composition-adjusted term; against its own
5,000-permutation null it sits at the **91.7th percentile, empirical p = 0.083**, and adding a
timepoint term moves it to **p = 0.088**.
**Replacement:** *"the one axis-by-genotype interaction in the corpus that reaches nominal
significance (p = 0.0052 adjusted for cell composition; empirical p = 0.083)"* — a ranking
statement, which is what n = 6 per cell can support.

**6. *"whereas no such link exists for other redox and metabolic axes"***
→ **Two errors, and the second one discards the sentence's own strongest argument.**
(a) **One comparator was tested, not several** — script 43 carries exactly two axes. Say "for the
redox axis."
(b) **Redox is not a null axis. It links, strongly, in both genotypes** (wild-type slope −7.70,
Myc+ −6.72; pooled ρ −0.48). What it does *not* do is link **differently** (interaction p = 0.72).
**Replacement:** *"the redox axis couples to the same ratio in both genotypes and to a similar
degree; what is specific to the respiratory chain is not the coupling but its **reversal** between
genotypes."*
That is a far better control than an axis nothing couples to, because it excludes the obvious
alternative — that the PUMA ratio simply tracks mitochondrial scores in general. **And the pooled
OXPHOS correlation is +0.04**: there is no main-effect coupling at all. The entire signal is the
difference of slopes, so the panel must never be drawn with a single fitted line.

**7. *"reduction of respiratory capacity (OXPHOS subunits)"***
→ **Capacity is not measured**; transcript is. And the rulers disagree on the load-bearing arm:
the same 87 subunits give **+0.061 / +0.201 / +0.226** on the diagonal depending on weighting.
**Script 46 now makes this concrete rather than pedantic:** the only ruler on which the Fig. 2I
coupling is arm-specific is mitoPPS, which is an *allocation* score. So the supported claim is a
change in respiratory **priority**, and "capacity" — a rate — is what the PGC1a perturbation is for.
Name the ruler, or say "respiratory-chain transcript".

---

## 6. The proposed census

| instrument | panels now | proposed | rationale |
|---|---|---|---|
| **log2 fold change** | 1E · 1F · 1G · 2F · 2G · 2H | **the spine — plus 2E's slot** | the field's native unit; carries a per-gene identity; the only ruler that can show a null geometrically |
| **% share** | 1E | **1E, promoted to lead Beat 1** | the one physical, unit-bearing mitochondrial claim in the corpus |
| **normalised counts** | none slotted | **optionally one supplementary** — *not* Fig. 2I: §4's test failed | rung 0. Legible, but for a large set it is not arm-specific (resid_frac 0.071), so it belongs only where the claim is about magnitude |
| **fGSEA NES** | 1D · 1H · 2E | **1D only** | the unbiasedness claim is about method and needs making once |
| **mitoPPS** | 1F · 2F · 2I | **1F only** (2I pending §4) | reallocation is 1F's subject; elsewhere it pays the teaching cost without the benefit |
| **GSVA** | 1B · (1C via scores) | **zero as a named ruler** | superseded as the coupling instrument by the linear z-score (script 36 is the method of record). Fine for *"is this programme higher in Myc+"*; a liability for *"does this track that"* |

**What the demotions cost, stated honestly.** 1H carries a real and non-obvious result — the
programme holds its rank while the amplitude halves (Spearman 0.933). Moving it to supplementary
means Fig. 1G's legend must carry that sentence, and 1G is the right home: it already draws the
rescaling. 2E's TEB/proliferation result is not lost either — it becomes two labelled points on
the two-timeline plane.

---

## 7. What to run

1. ~~The 2I ruler test~~ — **written and run: `scripts/46_axis_ruler_test.R`** (author to source it
   in Positron; it writes `results/axis_ruler_test.rds`). Verified in a write-free harness,
   20 s, every assertion passing. **Outcome: keep mitoPPS** (§4), with the reason now measured, and
   the "respiratory capacity" wording corrected to "respiratory priority" throughout.
2. **Nothing else.** Every remaining recommendation is a slot decision or a wording change, and
   needs no new computation.

## 8. The one-line version

**Move every panel one rung down the normalisation ladder — but only one, and check separability
before the last step — keep exactly one fGSEA panel and one mitoPPS panel, put the two-timeline
plane under Beat 2, and build the death section around the genome-wide dose-line scan rather than
around the apoptotic set.** The defence has already moved from the normalisation to the null; the
panels have not followed it.
