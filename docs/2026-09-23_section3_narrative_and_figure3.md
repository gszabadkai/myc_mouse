# Section 3 — narrative text and Figure 3

**2026-09-23, revised.** Assembles Section 3 of the Nature Metabolism Letter from
four sources: the figure reconciliation document of 12 September; §4 of the
reprioritisation narrative v5; the existing manuscript prose on the orthotopic
series; and **the orthotopic RNA-seq analysis itself, scripts 50 to 52**.

**The first version of this note was built without script 52**, which was run and
signed off at `faf1e8a`. Two things it said are wrong, and both are corrected
here. See §2.3 and §4.1.

Companion to `2026-09-22_section2_narrative_and_figure2.md` and to
`2026-09-12_figure_reconciliation.md`.

⬚ marks a number or statement held from conversation rather than from a document
or a panel.

---

## 1. What Section 3 claims

One sentence: **tumours that combine deregulated MYC with a high respiratory
state exist only where the death pathway has been broken, there are two classes
of break, and the transcriptome records which one was taken.**

The argument runs in four steps.

1. Cells selected to survive continuous MYC induction, with MYC protein
   unchanged, split by whether they carry a buffer. Without Bcl-xL they lower
   the chain and PUMA together. With Bcl-xL they keep both raised.
2. The same split appears with the biogenesis driver. Passaged iMMECs carrying
   MYC and Pgc1a lower the chain and PUMA, whilst MYAZ cells surviving Pgc1a
   instead lower p19^ARF, p21 and MYC itself.
3. Tumours do the same in vivo. PGC-1α alone is silenced and growth matches
   control. With Bcl-xL the transgene is retained, the chain rises, PUMA rises
   and the tumours outgrow both controls.
4. Across three sequenced PGC-1α arms, each against its own matched control, the
   retained dose scales with the respiratory share and `Bbc3` falls
   monotonically with it, going flat where the transgene was lost.

**The step that makes it a section rather than a list.** Every arm in the
orthotopic series is a surviving population, so the transcriptome there reads
the route taken rather than a response to the manipulation. That is a limitation
of the design read as a feature, and it is the honest form of the claim.

---

## 2. Figure 3 panel map

The in vitro block keeps the reconciliation's letters, with one move: the MYAZ
qPCR sat at 3J, after the in vivo block, and belongs beside the MYAZ survivor
western it explains. **The in vivo block is rebuilt around what script 52
actually produced**, which is different from and better than what the
reconciliation anticipated.

| new | source | shows | status |
|---|---|---|---|
| **3a** | Fig 3C | ETC subunits and c-MYC, 48 h against SS, EV and Bcl-xL | exists |
| **3b** | Supp 3D *(promoted)* | PUMA and BIM, same comparison | exists |
| **3c** | Supp 3C *(promoted)* | cell death, 48 h against SS | exists |
| **3d** | Supp 3E *(promoted)* | **MYC unchanged in SS**, the control that licenses the block | exists |
| **3e** | thesis Fig 4.3 C, D | iMMEC MYC + PGC-1α after passaging: ETC down, PUMA down, death no longer elevated | exists ⬚ |
| **3f** | thesis Fig 4.2 I, J | MYAZ + PGC-1α surviving fraction: p19^ARF down, p21 down | exists ⬚ |
| **3g** | MYAZ qPCR *(was 3J)* | `Myc` and `p19arf` down under PGC-1α and NRF1, **unchanged under Ndi1** | new figure ⬚ |
| **3h** | Fig 4A, 4B, 4C | PGC-1α alone in vivo: transgene silenced, no ETC change, growth = control | exists |
| **3i** | Fig 4D, 4E | Bcl-xL + PGC-1α: transgene retained, ETC up, growth up | exists |
| **3j** | Supp 4A *(promoted)* | PUMA protein **up** in the Bcl-xL + PGC-1α tumours | exists |
| **3k** | script 52 | **the retained-dose series**: `Ppargc1a` and the respiratory rulers across the three PGC-1α arms, **ordered by route**, each against its own control | **analysed, no panel** |
| **3l** | script 52 | **the `Bbc3` dose-response**, monotone across three backgrounds, with `Cdkn2a` beside it and four canonical p53 targets as non-movers | **analysed, no panel** |
| **3m** | script 52 | **the MYC ceiling**, present only in the buffered arm | **analysed, no panel** |

Thirteen panels, as before, but three of them are now analysed results with
intervals rather than things that might be built.

**Two panels are cited twice in the draft and carry two claims each**, which
changes how they must be drawn.

**3k** carries the dose-against-respiratory-share relationship *and* the
head-to-head route comparison. For the second to be readable, **order the three
arms by route rather than alphabetically or by construct**: effector blocked
(301.6 CPM), pathway intact (128.6), checkpoint broken (0.17). Label each arm by
what was done to the death pathway with the construct as a subscript, so a reader
sees that the deliberately broken checkpoint sits at the bottom. Drawn as three
unordered bars with construct names, the comparison is present in the data and
invisible on the page.

**3l** carries the `Bbc3` dose-response *and* its specificity, with `Cdkn2a`
beside it and the four canonical p53 targets as non-movers.

Every other panel is cited once and can be drawn to its own logic.

### 2.1 The in vitro block is the strongest content in the dataset

Panels 3a to 3g are four escape routes in two cell lines, with nothing re-cut and
MYC protein constant where it was measured. They are engineered rather than
inferred, which is what the orthotopic arm on its own cannot be.

| route | panel | what the cell gave up |
|---|---|---|
| lower the respiratory driver | 3a EV arm, 3e | the chain, and PUMA with it |
| block the effector downstream | 3a Bcl-xL arm | nothing, the buffer absorbs it |
| break the upstream checkpoint | 3f, 3g | p19^ARF and p21 |
| lower the oncogene | 3g | MYC itself |

**3d is not optional.** Without the demonstration that MYC protein is unchanged
in the survivors, a reader can attribute the whole survivor state to MYC
silencing, and the block collapses.

**3g carries a control that is easy to overlook.** `Myc` and `p19arf` fall under
PGC-1α and NRF1 but not under Ndi1. The checkpoint break therefore tracks the
comprehensive biogenesis drivers rather than the NADH-oxidation mimics, which is
the same split between intervention classes that Figure 2 relies on.

**What the MYC control licenses, and what it does not.** MYC is driven from a
synthetic promoter that holds expression high, so 3d rules out the obvious
reading, that survivors lowered their chain because they lowered the driver.
The second exclusion is the design rather than the panel: both arms went through
the same fourteen-day protocol and only the unbuffered one lowered the chain, so
it is not a general effect of prolonged culture. **Neither exclusion supports
"selected" over "adapted".** A fourteen-day survivor population is compatible
with differential survival of pre-existing variants and with a state change
inside the same cells, and nothing here separates them. The draft therefore
states the two exclusions and stops. Constant MYC protein is also not constant
MYC activity, which depends on MAX availability, MXD competition and T58/S62
phosphorylation, none of it measured; that one is for the Discussion if a
reviewer raises it.

**3e makes the in vivo block a replication rather than a separate observation.**
Passaged iMMECs carrying MYC and PGC-1α lower the chain, which is the in vitro
analogue of the transgene silencing in 3h.

### 2.2 The design of the sequenced series

**Seven arms, seven animals each, 49 of the 67 sequenced samples.** The other 18
are MYAZ, FaMY1 and BoMY, a separate series in a different animal-ID namespace
with no MYAZ-against-EV contrast on disk.

| arm | PGC-1α | background | matched control |
|---|---|---|---|
| EV | no | parental | — |
| BclxL | no | Bcl-xL | EV |
| **Pgc1a_BclxL** | yes | Bcl-xL | BclxL |
| NTEV | no | non-targeting CRISPR | — |
| **NTPgc1a** | yes | non-targeting CRISPR | NTEV |
| KOEV | no | p19^ARF KO | — |
| **KOPgc1a** | yes | p19^ARF KO | KOEV |

**Three PGC-1α arms, each against its own control, is the design feature the
whole in vivo block rests on.** `NTPgc1a` against `NTEV` is the only contrast in
the dataset where PGC-1α is present without a Bcl-xL construct, and it is
therefore the only place `Bcl2l1` means endogenous BCL-XL.

### 2.3 What script 52 found, and what it overturns

Every number below is from the signed-off run, which reproduced its dry run in 23
of 23 non-date objects with 38 numbers checked back against the object.

**The dose-response, and it is the panel of the section.** `Bbc3` against its own
control in each background:

| arm | retained `Ppargc1a` CPM | `Bbc3` against its control |
|---|---|---|
| Pgc1a_BclxL | 301.6 | **−0.731 [−1.120, −0.450]** |
| NTPgc1a | 128.6 | **−0.345 [−0.890, −0.010]** |
| KOPgc1a | 0.17 | +0.044 [−0.231, +1.180] |

Monotone across three independent backgrounds, and flat where the transgene was
lost. **Write it as selection**: tumours that retained more PGC-1α show less
`Bbc3`, and nothing here separates regulation from selection against
high-`Bbc3` cells.

**The construct-free arm tilts the compartment without expanding it.**
`NTPgc1a` against `NTEV`: `ox_rel` +0.176 [+0.039, +0.348] and mitoPPS
+0.053 [+0.030, +0.116], both excluding zero, **with the denominator flat**. This
is the in vivo tumour version of Figure 1's priority-against-content distinction,
and it arrives in the one arm with no construct in the way.

**The MYC ceiling survives only where the dose was kept.** 80th-percentile
contrast −802, p 0.012 on CPM and 0.029 on median-of-ratios, against nothing in
NT or KO. That arm carries 14.7% of its library in MitoCarta against 7.8 to 8.6%
elsewhere, so CPM deflation was the obvious objection and both normalisations
agree.

**Two overturned readings, and both corrections improve the section.**

*The p53 programme, and why this is a split result rather than a negative.* The
p53 target panel was run. `Cdkn2a` moves (−0.210) whilst the four canonical p53
targets do not, and `Bbc3` falls three times harder (−0.731). **Read as a test of
a coordinated p53-axis break it fails.** Read as a test of which part of the
checkpoint arm is engaged it succeeds, and it says something the in vitro block
cannot: the sensor moves and the output does not.

**That maps onto the in vitro result and makes the section's premise carry.**
MYAZ survivors in culture lowered both p19^ARF and p21 (3f, 3g). In vivo only the
`Cdkn2a` half moves. **The checkpoint route is therefore taken in both systems
but used differently**, which is what "escape routes" should look like if they
are routes rather than a programme. Section 2 closes by predicting that
MYC-plus-respiration tumours arise only where the pathway was interrupted, and
dropping this leaves the in vivo half showing two interruptions where the in
vitro half shows four.

**So 3l carries both, in one panel.** `Bbc3` and `Cdkn2a` against retained dose,
with the four non-movers plotted beside them. A reader sees the effect and its
specificity in the same object, rather than meeting the specificity later in
Extended Data.

**One thing to verify before drawing it.** `Cdkn2a` is a designed variable in
this series, since one background is a p19^ARF knockout. The −0.210 must be the
buffered arm against its own control and not a pooled estimate, or the knockout
arm contaminates it. ⬚

*`Myc` and respiration.* Narrative v5 read the ceiling as "`Myc` falls only in
the arm where respiration actually rose". Respiration rose, in relative terms,
in the construct-free arm too, and `Myc` was not capped there. **The ceiling
tracks retained dose and the buffer, not the respiratory rise.**

**Two things not to over-read, from the note.** The composition gate's positive
control passes trivially, because no axis was material and the adjustment is
therefore a no-op. And `KOPgc1a`'s `Myc` maximum of 4611 is an outlier that makes
that arm's max contrast positive and uninterpretable.

### 2.4 Two escape solutions, one constraint — the section's spine

**The objection this answers, because a reviewer will raise it.** The paper argues
that raising matrix NADH oxidation on a MYC background raises PUMA and kills. The
orthotopic series shows the opposite sign: more retained PGC-1α, less `Bbc3`.

**The frame belongs before the results, not after them.** Both systems report on
survivors, so what either measures is the state a cell had to reach in order to
exist, not its response to the manipulation. Stated first, as sentence 2 of the
draft now does, the `Bbc3` fall arrives as the expected reading. Stated
afterwards, it arrives as damage control. Same sentence, opposite effect.

**The two systems solved the same constraint differently, and that is the
argument.** Section 2 closes by predicting that MYC and a high respiratory state
coexist only where the route to death has been suppressed. There are two ways to
satisfy that, and each system took one.

| system | what happened | how the constraint was met |
|---|---|---|
| iMMEC survivors, no buffer | ETC subunits **and** PUMA fall together (3b, 3c) | the combination was never reached, respiration was given up |
| iMMEC survivors, Bcl-xL | both stay raised (3b, 3c) | the effector was already blocked, so nothing had to give |
| passaged iMMEC + Pgc1a | chain and PUMA fall together (3e) | as the unbuffered survivors |
| orthotopic tumours | respiration **kept**, `Bbc3` falls with dose (3k, 3l) | the combination was reached and the trigger suppressed |

**Two different solutions is stronger evidence for the constraint than either
alone.** A single escape route is compatible with a generic culture artefact. Two
routes that are mechanically opposite, converging on the same requirement, are
not. The draft makes the pivot explicit in one sentence: *these tumours had
therefore kept the combination that the cultured survivors gave up.*

### The p19^ARF knockout is the paper's only head-to-head route comparison

**What it shows.** Retained `Ppargc1a` was 301.6 CPM with the effector blocked,
128.6 with an intact pathway and 0.17 with the checkpoint knocked out. The
knockout arm kept **less than the intact control**. Breaking the checkpoint in
advance did not allow the respiratory state to be held; blocking the effector
did. Nowhere else in the paper are two escape routes set against each other on
the same endpoint.

**It resolves an apparent tension in the in vitro data.** MYAZ survivors lower
p19^ARF and p21 spontaneously (3f, 3g), so the checkpoint route is clearly taken.
Yet engineering the break in advance bought nothing. The reading that fits both:
**p19^ARF loss is a correlate of the adaptation rather than the licensing step.**
Panel 3g supports that directly, since `p19arf` falls together with `Myc` under
PGC-1α and NRF1 but not under Ndi1. The closing sentence of the draft now says
so, that the checkpoint falls alongside rather than permitting.

**Four caveats, and the third is the serious one.** n = 7 per arm. Retention is
an outcome of selection, not a designed variable. ⬚ **An ARF-null line may simply
grow fast enough that the transgene becomes dispensable, which removes the need
rather than failing to license** — a growth-rate comparison between `KOEV` and
`NTEV` would separate these and the data to do it exist. And the knockout itself
needs validation on the record. ⬚

**It feeds Section 4 directly.** If two routes differ in whether they permit the
respiratory state to be kept, they should differ in what the resulting tumour can
do, which is the prediction Section 4 tests.

**The flat arm is what makes this testable rather than circular.** Where the
transgene was lost (0.17 CPM) no pressure was carried, nothing needed
suppressing, and `Bbc3` does not move (+0.044). Had `Bbc3` fallen at every dose
the selection reading would collapse into a generic PGC-1α effect.

**What would falsify it.** `Bbc3` falling in the p19^ARF knockout arm despite the
lost transgene. Or acute PGC-1α lowering PUMA in MYAZ rather than raising it,
which Figure 2c,d rules out in the same cell line. **That acute comparison is the
falsification anchor and it is deliberately not in the 486 words**, because with
the frame up front the section no longer needs to argue against itself. Keep it
for the response to reviewers.

**Two limits to keep on the record.** Nothing here separates suppression as
regulation from selection against high-`Bbc3` cells, and the draft says
"suppressed" rather than "downregulated" for that reason. And the human β₂ result,
where OXPHOS runs positive with `BBC3` at fixed MYC, must **not** be set beside
these coefficients as though the signs were comparable. Partial Spearman
correlations across established human tumours and a mouse dose-of-escape contrast
are different estimands, and the standing rule forbids quoting them as comparable.
The compatible statement is the general one: an acute response and a selected
endpoint need not share a sign on the same axis.

### 2.5 What the orthotopic transcriptome cannot say

For Methods and the Extended Data legend, stated once:

- Every arm is a survivor population, so a "PGC-1α effect" measured there is an
  escape-route signature rather than a response to a manipulation.
- **No transcript-level quantification exists.** The directory was searched
  exhaustively for `quant.sf`, `tx2gene`, length-scaled and transcript counts;
  the only salmon artefact is `salmon.merged.gene_counts.tsv` in three
  byte-identical copies. Gene-level `Bcl2l1` therefore pools Bcl-xL with Bcl-xS,
  which are functionally opposite products of the same locus.
- `Bcl2l1` is construct-contaminated in every arm carrying the Bcl-xL vector,
  which is why `NTPgc1a` against `NTEV` is the only readable guardian contrast.
- The retained dose in 3k is an outcome of selection. **Label the x-axis retained
  dose, not assigned dose.**
- The intervention moves the composition markers, so pooled loadings mislead:
  `adipose` is +0.621 pooled against +0.036, −0.393 and +0.714 by arm. Read the
  per-group breakdown before believing either verdict.

### 2.6 Panel budget

**First out: the MYAZ against iMMEC baseline.** It was 3M in the reconciliation
and it is now displaced by three analysed in vivo panels. This is a real loss,
because complex I is the one variable appearing in all three figures: part of the
subunit set the developmental window lowers in Figure 1, lowered in this tumour
line, and the arm Ndi1 replaces in Figure 2 (see §4.3 on the superlative that
used to sit here). **Recommendation: Extended Data, cited from the main text**, unless a
panel has to be cut for space, in which case reconsider 3g.

**The p53 target panel stays in the main figure**, folded into 3l, per §2.3. The
per-gene detail goes to Extended Data.

### 2.7 Supplementary for Figure 3

Supp 3H in the draft's numbering (p53 protein, qPCR of p53 targets, proliferation
with pifithrin), currently an unbuilt placeholder · the p53 target panel as a
negative · the MYAZ against iMMEC baseline · the orthotopic composition gate with
its positive control and the no-op caveat · the C2 normalisation dependence · the
Bcl-xS limitation · the `BclxL` against `EV` asymmetry, which is the control
showing the buffer alone does not move respiration (`ox_lvl` +0.093, `ox_rel`
−0.023, `ox_mt` −0.246, all three intervals covering zero).

**The Supp 3 numbering conflict recorded in the Section 2 note applies here
too**, and it matters more for this figure because three of its panels are
promotions out of Supp 3.

---

## 3. Draft text, 508 words

Metrics against `STYLE.md`: 23 sentences, mean 22.1 words, median 24, range
13–32, 39% under 20 words. This is the longest of the four sections, which is
defensible because it carries two systems and three sequenced backgrounds.
**Decision taken 2026-09-23: length stands for the initial submission**, now 508
after the survivor frame moved to the front, the p19^ARF knockout result was
added (§2.4) and PUMA protein was restored (§4.4). §3.1 is held in reserve
for revision or for a length request from the editor; the first cut there brings
it to 460. No em-dashes, no "priming", British spelling, no
serial comma, nothing under the twelve-word floor.

> If raising mitochondrial NADH oxidation kills MYC-expressing cells,
> tumours that combine MYC with a high respiratory state should exist only
> where the death pathway has been broken. We therefore examined survivors
> in two systems, where the measured state is what a cell had to reach in
> order to exist rather than its response to the manipulation. In culture,
> iMMECs were kept under continuous MYC induction until they no longer died
> (Fig. 3a). Survivors carrying an empty vector had lowered their electron
> transport subunits and their PUMA together, whilst survivors carrying Bcl-
> xL kept both raised (Fig. 3b,c). MYC protein was unchanged in both arms,
> on the same protocol, so the lower respiratory state is neither a loss of
> the driver nor an effect of culture (Fig. 3d). Without a buffer the cells
> gave up respiration to survive MYC, and with one they kept it. iMMECs
> carrying MYC and Pgc1a reached the same state after passaging, lowering
> the chain and PUMA together (Fig. 3e). MYAZ cells surviving Pgc1a took two
> further routes, lowering p19^ARF and p21, and lowering MYC itself (Fig.
> 3f,g).
>
> To ask what tumours do, we implanted MYAZ cells expressing Pgc1a into the
> mammary fat pad on backgrounds with the death pathway intact, blocked at
> the effector or broken upstream. Three were sequenced against their own
> controls: Bcl-xL overexpression, a non-targeting guide and a p19^ARF
> knockout. On the parental background the transgene was carried as
> transcript but silenced as protein, and growth matched controls (Fig. 3h).
> With Bcl-xL co-expressed the transgene was retained, the subunits and PUMA
> protein rose and the tumours outgrew both controls (Fig. 3i,j). These
> tumours had therefore kept the combination that the cultured survivors
> gave up. Retained Ppargc1a ranged from near zero to three hundredfold, and
> the relative respiratory share rose with it whilst the compartment did not
> expand (Fig. 3k). Bbc3 fell monotonically with that dose, and went flat
> where the transgene had been lost (Fig. 3l). The tumours that kept most
> respiration damped the message most whilst the protein stayed up, the same
> constraint met by the opposite route in culture. Cdkn2a fell in the same
> tumours whilst four canonical p53 targets did not, so the checkpoint route
> was taken only in part (Fig. 3l). The p19^ARF knockout arm nonetheless
> kept least of the three, so breaking the checkpoint in advance did not
> allow the respiratory state to be held, whereas blocking the effector did
> (Fig. 3k). Myc was capped only in the arm carrying both the highest
> retained dose and the Bcl-xL buffer, and not in the construct-free arm
> (Fig. 3m).
>
> The parental MYAZ line is itself such a survivor, with complex I lowered
> and the matrix NADH to NAD⁺ ratio high (Extended Data Fig. 7). Complex I
> oxidises matrix NADH, falls in the developmental window, and is the arm
> Ndi1 replaces. These data indicate that the gate is passed by giving up
> respiration or by blocking the effector, whilst the checkpoint falls
> alongside rather than permitting it. Which route a tumour takes should
> therefore decide what it can do next.

The closing sentence hands off to Section 4 (consequence), in the same shape as
Sections 1 and 2.

### 3.1 Trim order, held in reserve

**Not to be applied at initial submission**, per the decision above. This is the
order to work in if a revision or an editorial length request forces a cut.

**Cut the `Myc` ceiling sentence first**, which takes the section to 483 words.
It is the weakest of the three sequenced results, at n = 6 to 7 with an outlier
already excluded from one arm, and nothing else depends on it.

**Second, merge the two in vitro MYAZ routes into the dose-response sentence.**
Costs a little, because "lowering MYC itself" is a route the reader otherwise
never meets.

**Do not cut the MYC-unchanged control or the two design sentences that open the
in vivo paragraph.** The first licenses the whole in vitro block. The second pair
states the strategy before the results, which is what lets a reader follow panels
3h to 3m without assembling the design from the findings.

**Two sentences changed in this revision.** The first version read "MYC protein
was unchanged in the survivors, so the respiratory change was selected and not
imposed", which asserted selection over adaptation and was inconsistent with §2.3,
where the same distinction is left open for the in vivo arms. It now states the
two exclusions the design supports and stops, and it has moved after the result
it licenses.

The earlier draft also read "That arm carried a second break upstream, with
p19^ARF and p21 lowered by western and Bbc3 falling in the same tumours." The
p53 target panel rejects the coordinated break that sentence asserts. It has been
replaced by two sentences: the `Bbc3` dose-response, and a narrower checkpoint
sentence saying that `Cdkn2a` fell whilst the canonical p53 targets did not. The
closing verdict now names three routes rather than two, so the in vivo half
matches the prediction Section 2 ends on.

---

## 4. What is not yet supported

### 4.1 The orthotopic identity correction, and how C1 actually closed

**The arm labelled `Pgc1a` in the vector series was PGC-1α + Bcl-xL
double-expressing cells.** Correction Phase 1 was committed at `02d3acb` and
tagged `orthotopic-pre-identity-correction`; Phase 2 reanalysed all 49 vector and
CRISPR samples as a dose-of-escape series in script 52, committed at `faf1e8a`.

**C1 was closed on evidence, not voided on a technicality, and the distinction
matters.** C1 held that PGC-1α alone raises `Bcl2l1` in vivo (+121.8 CPM,
one-sided p = 0.00029). Script 52 tested it where it is testable, in `NTPgc1a`
against `NTEV`, the only contrast with no Bcl-xL construct in the way. `Bcl2l1`
is **−0.113 [−0.283, +0.174]**, null, with `Mcl1` null beside it. PGC-1α does not
raise the guardian. That is a reportable negative rather than a missing result,
and it should be written as one.

**Consequence for the draft.** Any sentence giving the human `BCL2L1` coefficient
a causal reading from this dataset comes out. The licensing model, in which
PGC-1α raises its own buffer, is not supported in vivo.

**What survives, and now with numbers.** The `Bbc3` dose-response across three
backgrounds; the construct-free respiratory tilt without mass expansion; the MYC
ceiling in the buffered arm at both normalisations. All three are §2.3.

**The open question, load-bearing for three panels.** Whether the labelling error
affects only the RNA-seq samples or also the Prism panels in current Figure 4.
Panels 3h, 3i and 3j are those Prism panels. If the arm drawn as "Pgc1a" in
Fig 4A–C is the double-expressing line, 3h is not the transgene-silencing panel
the section needs. **Confirm before drafting the legend.**

### 4.2 Panels held from conversation rather than from a document

**3e and 3f come from thesis figures** (4.3 C,D and 4.2 I,J). The reconciliation
recorded that the adapted-PGC-1α iMMEC panels are described in the current Fig 3
legend but could not be located in the layout. ⬚ Confirm they exist and are
publication-grade.

**3g is a supplied new figure** ⬚ rather than one located in the existing set.

**The MYAZ against iMMEC baseline is flagged "blot? confirm"** ⬚ and now sits in
Extended Data, so the stakes are lower than in the previous version of this note.

### 4.3 "Complex I strongest" has never been computed

**Provenance, traced.** The claim originates as an eyeball reading of
`MC_OXPHOS _ OXPHOS subunits_raw.pdf` in the April bioinformatics notes,
recorded as "Strongest effects are on Complex I - follow up???" with an unticked
checkbox. It entered the finalisation plan and its v2 as a **task** ("Complex I
shows the strongest effect ... a plotting/interpretation task, with per-complex
correlations in a supplementary"), then hardened into an assertion in narrative
v3 Movement 1 and into "the strongest component of the developmental
reprioritisation" in v5 §4b. **No per-complex statistic exists.** What is
computed sits one tier up: OXPHOS is the most de-prioritised Level-1 tier on both
temporal axes, carried by structural subunits (−0.109 WT, −0.153 Myc+) whilst
assembly factors do not move (−0.002, −0.029).

**The draft no longer uses it.** The sentence now says complex I oxidises matrix
NADH, falls with the subunit set the window lowers, and is the arm Ndi1 replaces.
All three are supported, the loop across the three systems survives intact, and
the superlative was carrying rhetorical weight rather than argumentative weight.

**It is cheap to compute and worth computing**, because if it holds the stronger
sentence is available and it strengthens Figure 1 more than Figure 3. The
per-complex sets nest exactly inside `OXPHOS subunits` in MitoCarta Sheet 4,
under the `OXPHOS > Complex N > CN subunits / CN assembly factors` hierarchy, so
this is a set definition and a ranking, not a new analysis.

**Two things to control for when it is run.** Complex I is by far the largest
complex, so a set-level "strongest" may be reading set size rather than biology;
report against an expression-matched or size-matched null. And batch equals
timepoint in this dataset, so the **magnitude** of any developmental decline does
not survive, whilst a per-complex **ranking** does, on the same arm-selectivity
logic that licenses the rest of Figure 1.

**This belongs to Section 1 more than to Section 3.** Narrative v5 places it in
Figure 1 and Figure 3 only borrows it.

### 4.4 The transcript–protein inversion, resolved in an earlier session

**Correction to the previous draft of this note.** I speculated that the
comparators differ, protein against EV and transcript against `BclxL`. They do
not. `BclxL` and `Bcl-xL PGC1a` are the **same pool**: MYAZ transduced with
Bcl-xL, then split, then EV or PGC-1α, polyclonally selected. Both readouts run
against that matched control, and the contrast is valid for every gene except
`Bcl2l1` itself.

**Stabilisation by BCL-XL is excluded, not supported.** BCL-XL protein is equal
on the blot between the two arms, so buffering capacity is equal and no
differential sequestration is available. That route was closed cleanly in an
earlier session rather than left open. ⬚ Densitometry would turn "equal" into a
number, since a blot does not resolve the 1.6× the transcript shows.

**What stands: reduced transcriptional induction with retained
post-transcriptional accumulation.** Put the observations in time order.

| | `Bbc3` transcript | PUMA protein |
|---|---|---|
| in vitro, acute PGC-1α | **up** (qPCR) | up |
| in vivo, established tumours | **down** (−0.731) | up |

The transcript reverses between acute and chronic; the protein does not, and
that is the shape selection produces. PGC-1α acutely raises both. Bcl-xL makes
the residual survivable. Among Bcl-xL-protected cells, those that additionally
damped `Bbc3` transcription are the ones that grew out, whilst the
post-transcriptional component is not something transcriptional damping removes.
**The transcript is the selection signature and the protein is the persisting
effect of the manipulation.** It also explains the clean 7-against-7 separation
better than noise does, since selection produces uniform outcomes.

**The candidate for the post-transcriptional half, and it is speculative.** PUMA
is short-lived and mitochondrially associated, so its half-life depends on
available outer-membrane binding surface. PGC-1α nearly doubles the
mitochondrial share of the transcriptome, 7.96% to 14.73%. More membrane, more
retained PUMA, longer effective half-life, more protein on less message. It has
to be independent of BCL-XL given the protein is equal. A cycloheximide chase
would settle it.

**One earlier explanation has since failed its own test.** The first account
attributed the `Bbc3` fall to a damped p53 programme, since MYAZ is
p53-proficient and the survivors have p19^ARF and p21 down by western. It
predicted that `Cdkn1a`, `Mdm2`, `Ccng1` and `Zmat3` would fall with `Bbc3`.
Script 52 ran that panel and they do not (§2.3). **The selection reading
survives; the p53-programme mechanism does not.**

**What the draft now says.** PUMA protein is restored to the Bcl-xL sentence,
and the interpretive sentence reads "damped the message most whilst the protein
stayed up", which states the inversion without spending words explaining it.
The explanation belongs in the legend or the Methods.

**⬚ One reference to confirm.** The PUMA protein panel is Supp 4A in the
reconciliation and Supp Fig 5A in the earlier session. Given the Supp numbering
conflict already recorded for Section 2, settle which panel before the legend is
written.

### 4.5 One route that is asserted but not measured

**BCL-XL protein in the adapted survivor lines was not measured.** If the
survivors raise their own BCL-XL, that is a fourth route and the one most
directly comparable to the human data, where `BCL2L1` amplification co-occurs
with `MYC` amplification. It is a western, not an experiment. **It also now
carries more weight**, because script 52 has closed the transcript-level version
of the same question with a null.

---

## 5. Used, moved, omitted

**Used.** The iMMEC selected-survivor block with its MYC-unchanged control; the
passaged PGC-1α iMMEC panels; the MYAZ PGC-1α surviving fraction; the MYAZ qPCR;
the orthotopic PGC-1α-alone and Bcl-xL + PGC-1α arms with PUMA; the script 52
retained-dose series, `Bbc3` dose-response and MYC ceiling.

**Moved to Extended Data.** The p53 target panel, as a negative. The MYAZ against
iMMEC baseline. The composition gate, the C2 normalisation dependence, the
Bcl-xS limitation and the `BclxL`-against-`EV` asymmetry control.

**Moved to Section 4.** Everything on metastasis: spontaneous lung metastasis,
the IV colonisation series, Twist1 with MMP2 and MMP9, circulating tumour cells,
the LbNOX in vivo growth result and the spontaneous metastatic MMTV-Myc lines.

**Closed, and stays closed.** C3, the `Mcl1:Bcl2l1` reversal against `ox_rel`,
uninformative for two measured reasons and neither is power. The abstract
sentence about mouse tumours resembling human tumours must not be cited to this
dataset.

**Open, and it belongs to Section 4.** The PGC-1α-high and PGC-1α-low cooperation
model with the 1:1 mixing experiment. Narrative v5 proposed reducing it to one
discussion sentence; that decision has not been taken.

**Not claimable.** That PGC-1α alone raises `Bcl2l1` in vivo (§4.1). That the
retained dose was assigned. That `Bbc3` falls as part of a p53 programme. Any
juxtaposition of a mouse `bMX` coefficient with a human partial Spearman. No
"primed" or "priming", anywhere.

---

## 6. What has to be built for Figure 3

1. **Panels 3k, 3l and 3m from script 52.** The analysis is done and signed off,
   so these are drawing jobs against a saved object, not computations.
2. **Confirm the Prism in vivo panels against the identity correction** (§4.1).
   Nothing else in the in vivo block can be trusted until this is settled.
3. **Locate 3e, 3f and 3g** and confirm they are publication-grade.
4. **Locate or abandon the MYAZ against iMMEC blot** for Extended Data.
5. **A western for BCL-XL protein in the adapted survivor lines** (§4.5), which
   is the one new experiment this section would benefit from and the only one.
6. **The per-complex ranking across the developmental window** (§4.3). Cheap, it
   belongs to Figure 1, and it decides whether the stronger complex I sentence is
   available.
