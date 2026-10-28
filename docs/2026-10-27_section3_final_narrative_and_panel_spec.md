---
date: 2026-10-27
version: v2 (amends v1 of the same date)
status: >
  Narrative settled. Panel inventory settled apart from one conditional panel,
  listed in section 7. Four numbers still open, listed in section 6. Supersedes
  2026-09-23_section3_narrative_and_figure3.md for the narrative, the panel
  inventory and the Extended Data contents. The September note remains
  authoritative for the script 52 analysis itself, the design of the sequenced
  series, and the limitations recorded in its section 2.5.
purpose: >
  Context for the script 57 panel build. Gives the prose the panels have to
  match, what each panel is for, what was dropped, and what must not be claimed.
changes-from-v1: >
  Four p53 non-movers with Cdkn2a reinstated to Extended Data; MYAZ against
  iMMEC complex baseline reinstated to Extended Data, now with quantified blots;
  composition adjustment moved to Methods as prose rather than a panel; Myc
  ceiling held conditional on the KOPgc1a Myc result. Each reinstated panel now
  carries a main-text sentence, because legends in this manuscript describe
  panels only.
relates-to:
  - 2026-09-23_section3_narrative_and_figure3.md (superseded in part, see above)
  - 2026-10-23_MYAZ_proteomics_mass_redox_and_bcl2_family.md (section 6 flag
    against the p19ARF measurement in Fig. 3f)
  - script 52, committed at faf1e8a
---

# Section 3: final narrative and panel spec

## 0. Standing rule that governs this build

**Legends describe the panels on the figure. No new result and no
interpretation goes in a legend.** Every panel therefore needs a main-text
citation and whatever sentence the result requires. This is why the Extended
Data contents below are shorter than the available space would allow: the
constraint is main-text words, not Extended Data room.

---

## 1. What the section claims, and what changed

**One sentence.** Deregulated Myc and a high respiratory state coexist only where
the mitochondrial apoptotic effector has been blocked, and cultured cells and
orthotopic tumours give the same answer when buffer status is matched.

**What changed from the September note.** That note built the section around
*two escape solutions, one constraint*: cultured cells gave up respiration,
tumours kept it and suppressed the trigger. That framing existed to defend the
`Bbc3` fall against the acute PUMA rise in Fig. 2, and it rested on a comparison
between the **buffered** in vivo tumours and the **unbuffered** in vitro
survivors. Buffer status is the one variable matched across both systems, so
dropping it manufactured a contrast that is not in the data.

Aligned by buffer status the two systems agree at every level:

| apoptotic status | in vitro | in vivo |
|---|---|---|
| effector available | iMMEC+EV survivors lower OXPHOS subunits and PUMA together (3b); passaged iMMEC+Pgc-1α do the same (3d,e) | NT+Pgc-1α silences the protein, OXPHOS subunits unchanged, growth = control (3g) |
| effector blocked (Bcl-xL) | survivors keep both raised, no death (3a,b) | transgene retained, OXPHOS up, PUMA protein up, tumours outgrow both controls (3h) |
| checkpoint removed | MYAZ surviving fraction lowers p19ARF, Cdkn1a and Myc (3f) | p19ARF KO arm retains least of the three (3i) |

The section is therefore a replication across systems with one discriminating
negative, not two different solutions. It hands off directly to Fig. 4, which
uses only the Bcl-xL route.

**The route-outcome ladder**, which is the section's verdict shown as a tumour
phenotype rather than inferred from a transcript:

| route | retained dose | growth against its own control |
|---|---|---|
| effector blocked (Bcl-xL) | highest | up |
| pathway intact (NT) | partial | unchanged |
| checkpoint removed (p19ARF KO) | lost | down |

The third row closes the obvious objection. If the transgene were simply
superfluous in a fast-growing ARF-null background, losing it should leave the
ARF-null growth advantage intact, so `KOPgc1a` should grow like `KOEV`. It does
not: `KOPgc1a` grows at the `NTEV` rate, significantly below `KOEV`.

---

## 2. Final prose

Panel calls in this text are the ones the build must produce. Paragraph 1 is the
author's existing text with three corrections applied and marked. Paragraphs 2
and 3 are drafted against `STYLE.md`: 18 sentences, 383 words, mean 21.3,
median 22, range 15–33, 39% under 20 words, no sentence under the 12-word floor,
no em-dashes, British spelling, no serial comma. The section totals roughly 511
words with paragraph 1, against 508 for the September version, which was
accepted for initial submission. There is no slack left.

> We tested this hypothesis by evaluating survival in cellular and orthotopic
> mouse tumour models with varying Myc and OXPHOS levels to uncover the required
> survival state. *In vitro*, we kept iMMECs under sustained Myc induction until
> resistance to cell death was acquired (Fig. 3a). Survivors carrying an empty
> vector had lowered their OXPHOS subunits and their PUMA levels, whilst
> survivors carrying Bcl-xL kept both raised (Fig. 3b). Myc levels remained
> constant across both groups, thus the reduced respiratory capacity is neither
> a loss of the driver nor an effect of prolonged culture (Fig. 3c). iMMECs
> carrying Myc and Pgc-1α (*see Fig. 2e*) reached the same state after several
> passages, lowering OXPHOS, PUMA, Bax together with apoptotic sensitivity
> (Fig. 3d,e). On the other hand, the Myc-adapted MYAZ cells surviving Pgc-1α
> induced apoptosis took two further routes, lowering *Cdkn2a*/p19ARF, *Cdkn1a*
> and *Myc* transcription itself (Fig. 3f).

> To ask whether tumours take the same routes *in vivo*, Pgc-1α expressing MYAZ
> cells were implanted orthotopically on three backgrounds. These differed in
> where the apoptotic pathway was interrupted: (i) a non-targeting guide
> control, leaving the pathway intact, (ii) Bcl-xL overexpression, blocking the
> mitochondrial effector, and (iii) p19ARF knockout, removing the upstream
> checkpoint. On the non-targeting control background the *Ppargc1a* transcript
> was retained whilst Pgc-1α protein was no longer detectable. OXPHOS subunit
> protein was accordingly unchanged and the tumours grew at the rate of the
> control cohort (Fig. 3g). Conversely, co-expression of Bcl-xL sustained the
> transgene, raised OXPHOS subunits and PUMA protein, and allowed the tumours to
> outgrow both control groups (Fig. 3h). At matched Bcl-xL status the two
> systems therefore behaved alike, lowering the respiratory programme where
> apoptosis could still proceed and retaining it where the mitochondrial
> effector was blocked. Sequencing the three arms against their own controls
> showed that retained *Ppargc1a* was accompanied by a corresponding rise in the
> relative OXPHOS share (Fig. 3i). The highest dose was carried by the Bcl-xL
> arm, whilst the p19ARF knockout retained less than the intact background
> (Fig. 3i). Loss of the transgene there was not neutral, since these tumours
> also gave up the growth advantage that the knockout confers on its own
> (Extended Data Fig. 5). In the same tumours *Bbc3* fell with the retained dose
> and did not move where the transgene had been lost, whilst PUMA protein
> remained elevated (Extended Data Fig. 5). Four canonical p53 targets did not
> move in these tumours, whilst *Cdkn2a* fell alongside *Bbc3* (Extended Data
> Fig. 5). The fall in *Bbc3* is therefore unlikely to reflect a change in the
> p53 transcriptional output.

> Notably, the parental MYAZ line had itself lowered complex I relative to
> immortalised mammary epithelial cells, whilst the other respiratory complexes
> were unchanged (Extended Data Fig. 5). Complex I oxidises matrix NADH, falls
> with the subunit set deprioritised in the developing gland, and is the arm
> that Ndi1 replaces. Altogether, the gate that the developing gland passes by
> deprioritising OXPHOS (*see Fig. 1f*) is passed in established tumours only
> where the mitochondrial effector has been blocked. Removal of the upstream
> p19ARF checkpoint in advance did not substitute for that blockade. The
> surviving fraction in culture nonetheless lowered p19ARF alongside *Myc*
> itself (*see Fig. 3f*), suggesting that its loss accompanies the adaptation
> rather than licensing it. We therefore asked whether progression and
> metastasis of high-OXPHOS tumours depend on the same blockade at the
> mitochondrial checkpoint.

**The three corrections in paragraph 1.** The panel call for the passaged
Pgc-1α iMMECs was `3e` and becomes `3d,e`, since the built figure splits the
western from the death assay. `Cdkna2` and `Cdkna1` are typos for `Cdkn2a` and
`Cdkn1a`. The Myc-unchanged sentence had a stray full stop before its citation,
and claimed "selective pressure", which the design does not support: constant
Myc excludes driver loss, and both arms running the same protocol excludes
culture drift, but neither separates selection of pre-existing variants from a
state change within the same cells.

---

## 3. Panel inventory

`a` to `h` exist in `NM_Fig3.pdf`. One main panel and three Extended Data panels
are to be built, with a fourth held conditional.

| panel | content | status |
|---|---|---|
| 3a | iMMEC+Tre-Myc, % apoptotic cells, GFP against Myc at 48 h and continuous, in EV and Bcl-xL | built |
| 3b | western: Myc, Bcl-xL, ETC subunits, PUMA, BIM, same comparison | built |
| 3c | Myc expression relative to GFP, showing it unchanged in survivors | built |
| 3d | western: passaged iMMEC+Tre-Myc, EV against Pgc-1α, ±dox, with Bax 6A7 | built |
| 3e | passaged iMMEC apoptosis, death no longer elevated by Pgc-1α | built |
| 3f | MYAZ qPCR, surviving fraction: p19ARF, Cdkn1a, Myc | built |
| 3g | NT orthotopic: growth, western, **Ppargc1a bar (placeholder data, rebuild)** | **rebuild one sub-panel** |
| 3h | vector series orthotopic: growth, western, PUMA protein | built |
| **3i** | **route comparison: retained *Ppargc1a* and relative OXPHOS share, three PGC-1α arms each against its own control, plus BclxL against EV as a zero arm** | **build** |
| **ED a** | ***Bbc3*** **against retained dose, three backgrounds, with *Cdkn2a* and four canonical p53 targets as a specificity strip** | **build** |
| **ED b** | **four-arm growth: NTEV, NTPgc1a, KOEV, KOPgc1a** | **build** |
| **ED c** | **MYAZ against iMMEC respiratory complex baseline, from quantified blots** | **build** |
| *ED d* | *Myc across arms* | **conditional, see §7** |

**3g is one of the three sequenced arms, not a separate experiment.** It is
`NTEV` against `NTPgc1a`, so it shares data with 3i. The three backgrounds named
in the design sentence map one to one onto the three sequenced PGC-1α arms.

**3g and 3h sit on different control backgrounds.** 3g's control is `NTEV` from
the CRISPR series; 3h's are `EV+EV` and `Bcl-xL+EV` from the vector series. Each
contrast is against its own matched control, which is correct, but the legend
must say so, since that is description rather than interpretation.

**ED c is not a job against the sequencing object.** It comes from quantified
blots. Confirm where the quantification lives before building it.

### To Methods as prose, with no panel

| item | what Methods should say |
|---|---|
| composition adjustment | Scores recomputed with composition-marker adjustment; estimates for `ox_rel` and `Bbc3` unchanged. Record the marker loadings **per arm**, not pooled: `adipose` runs +0.036, −0.393 and +0.714 across the three arms, so the pooled value of +0.621 misleads. |
| Bcl-xS limitation | No transcript-level quantification exists in this dataset, so gene-level `Bcl2l1` pools Bcl-xL with Bcl-xS, which are functionally opposite products of the same locus. |
| survivor caveat | Every arm is a survivor population, so a "PGC-1α effect" measured there is an escape-route signature rather than a response to a manipulation. |

---

## 4. Drawing requirements that come from the prose

1. **3i must be ordered by route, not by construct name and not alphabetically**:
   effector blocked, pathway intact, checkpoint removed. Label each arm by what
   was done to the death pathway, with the construct as a subscript. Drawn as
   unordered bars with construct names, the head-to-head route comparison is
   present in the data and invisible on the page.
2. **The x-axis of 3i reads "retained", never "assigned".** Retention is an
   outcome of selection, not a designed variable.
3. **One unit across 3g and 3i.** Script 52 quotes CPM; the existing Prism axis
   says DESeq2-normalised. Pick one and state which.
4. **3g needs per-animal values, not summary statistics**, since Prism computes
   its own statistics from the pasted replicates.
5. **The legend must cite one statistical test consistently.** A Prism unpaired
   t-test on exported counts will not reproduce the DESeq2 Wald P for the same
   contrast. Have both numbers before the legend is written.
6. **ED a: the four p53 targets are `Cdkn1a`, `Mdm2`, `Ccng1` and `Zmat3`**, the
   panel script 52 already ran. Label all four. Do not drop `Cdkn1a` to avoid
   the collision described in §5.
7. **ED c must be labelled descriptively**, as a difference between two lines
   rather than an effect of selection. See §5.
8. **Extended Data Fig. 5 is a placeholder number.** The manuscript has ED 1–4
   plus two unnumbered figures in the human section; settle the numbering first.

---

## 5. What must not be claimed from these panels

- **No causal verb for the retained-dose relation.** Every arm is a survivor
  population at endpoint, so 3i reports an association across escape routes.
  "Was accompanied by", not "drove".
- **ED c reports a difference between lines, not an effect of selection.** MYAZ
  and iMMECs differ in origin (tumour against normal), in Myc status, and in
  immortalisation method (E1A, p53DN, SV40-LT). The complex I difference is
  consistent with the section's argument and cannot be attributed to it. The
  September draft's phrasing, "the parental MYAZ line is itself such a
  survivor", claims a cause this comparison does not support and is not used.
- **`Cdkn2a` must be the `Pgc1a_BclxL` arm against its own control**, not a
  pooled estimate, since one background in the series is a p19ARF knockout.
- **No "primed" or "priming"** for any transcript, cell line or animal (N3).
- **Never set these mouse coefficients beside the human partial Spearman
  correlations** in Fig. 4. Different estimands; the standing rule forbids it.
  The compatible statement is the general one, that an acute response and a
  selected endpoint need not share a sign on the same axis.
- **PGC-1α does not raise `Bcl2l1` in vivo.** Tested where it is testable, in
  `NTPgc1a` against `NTEV`, the only contrast with no Bcl-xL construct in the
  way: −0.113 [−0.283, +0.174], with `Mcl1` null beside it. This is a reportable
  negative. Any sentence giving the human `BCL2L1` coefficient a causal reading
  from this dataset comes out.

**One collision that will be visible and is deliberately not addressed in the
prose.** `Cdkn1a` is one of the four p53 non-movers in ED a and also appears in
3f, falling in the in vitro surviving fraction. Down in culture, flat in
tumours. That is what the narrative predicts, since those arms differ in buffer
status and the summary paragraph already states that the checkpoint route was
taken in culture, but it is not spelled out. Expect the question; do not spend
words pre-empting it, and do not drop `Cdkn1a` from the panel.

---

## 6. Open items that could still change a panel

1. **Within the `Pgc1a_BclxL` arm only, does `Bbc3` track `Ppargc1a` across its
   animals?** This decides whether ED a is a dose-response or a three-point
   comparison between three genotypes, since the three arms differ in background
   as well as in dose. Report it whichever way it comes out and describe the
   panel accordingly. If it is flat, the Extended Data sentence must say the
   three backgrounds differ rather than implying a dose-response within tumours.
2. **`NTPgc1a` `ox_rel`.** The September note records +0.176 [+0.039, +0.348],
   excluding zero, at transcript level, in the same arm whose western shows ETC
   subunits unchanged. The prose handles this by saying "OXPHOS subunit protein"
   in the 3g sentence and by making 3i a dose-scaled relation, so a modest
   transcriptomic tilt with no visible protein change sits where it should.
   Confirm the number before the legend is written.
3. **Whether the four growth cohorts were implanted in the same experiment.** The
   Methods record mixed genotypes per cage to avoid cage effects, so this is
   probably fine, but ED b depends on it. If the tumour volumes exist only in
   Prism and not in the repo, say so.
4. **The p19ARF measurement behind 3f.** The October proteomics note found
   CDKN2A protein unchanged in all three MYAZ construct arms, and `Cdkn2a`
   encodes p16INK4a and p19ARF from alternative reading frames sharing exons 2
   and 3. Confirm the qPCR primer set is exon-1β specific. The panel axis says
   p19ARF, so it probably is, but the checkpoint sentence in the summary
   paragraph depends on it.

---

## 7. The one conditional panel

**ED d, `Myc` across arms, is held pending one number: `Myc` in `KOPgc1a`
against `KOEV`.**

*Hypothesis, not a conclusion.* If ARF-null cells grow faster because they
tolerate a higher Myc level, then adding PGC-1α to that background creates
exactly the Myc-plus-high-respiration combination the paper argues is lethal,
and the only escape is to lose the transgene *and* come down on Myc, the second
half of which costs the growth advantage. That would also tie 3f together, where
p19ARF and `Myc` fall as a pair in the surviving fraction.

**If the number is positive**, `Myc` comes down under PGC-1α pressure in three
places: the buffered arm (the September ceiling, −802 at the 80th percentile),
the KO arm, and 3f in culture. That is a secondary route replicated across both
systems, worth one panel and one main-text sentence.

**If it is null**, the ceiling stays what the September note called it, the
weakest of the three sequenced results at n = 6 to 7 with an outlier already
removed, and ED d is not built.

**Either way it must be an arm median or mean.** Script 52 flags a `Myc` outlier
of 4611 in `KOPgc1a` that makes the 80th-percentile and maximum contrasts
uninterpretable.

**And either way there is a cost to naming it.** The section's premise is that
these tumours kept Myc and respiration together, and 3c has just established Myc
unchanged in the in vitro survivors. "Myc is capped" softens that, so the
sentence would have to say ceiling rather than abolition, precisely.
