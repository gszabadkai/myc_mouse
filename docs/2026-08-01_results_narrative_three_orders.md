# Three narrative orders for the attenuation / substrate / death sections

**Status:** drafting note, 2026-08-01. Written for the author to read away from the session.
The three orders are written out **as continuous Results prose**, in the paper's voice, so the
*shape* of each can be judged rather than its beat list. They cover the eleven content items the
author listed on 2026-08-01, and all three end at the hand-off to the cell and MYC-ER work
(PGC-1a -> PUMA -> death on a Myc background; PGC-1a with death blocked -> growth).

**Nothing here is a figure plan.** Where a panel would obviously be needed I have written
`[panel]` and moved on. The figure partition comes after the order is chosen.

**Recommendation: Order C**, with one sentence borrowed from B (see the comparison at the foot).
But read all three first -- that is what this file is for.

**Provenance.** Every number is from an object already on disk and already reconciled against its
author-run script: `interaction_results.rds` (raw, unshrunken), `background_vs_myc.rds` (script 40),
`attenuation_decomposition.rds` (script 29 / Issue #4), `priming_arm_teb.rds` (script 42),
`substrate_specificity_tradeoff.rds` (script 43), `collapse_module_ownership.rds` (script 44),
`mitopps_priming_pgc1a.rds` (script 38), `ap7_mb_fork.rds` (script 18), `gsva_overview.rds`.
Blot, IHC and MYC-ER values are the author's, reported as given. No script was re-run.

---

# ORDER A -- subtract the dose, then see what is left

*Elimination-first. Establish the null explanation, exhaust it, and let everything that survives
subtraction be the finding.*

## Section 1 -- "The Myc transcriptional programme is scaled down, not rewritten, in the adult gland"

The Myc+ and wild-type transcriptomes were already widely separated at six weeks -- 2,777 genes
differentially expressed at 10% FDR and 1,967 at 5% -- and by twelve weeks that separation had
largely closed gene by gene, to 239 and 135. The collapse in gene counts overstates what changed.
A differential-expression count is a threshold statistic, and under a uniform reduction in effect
size it falls far faster than the effect size itself: the genes that still clear the threshold at
twelve weeks are simply the largest ones, and their median absolute fold change is accordingly
*higher* than at six weeks, not lower. Measured directly, the change is a rescaling. Regressing the
Myc effect at twelve weeks on the Myc effect at six weeks across the 143 nuclear-encoded
mitochondrial pathways returns a slope of 0.55, an intercept indistinguishable from zero and
R2 = 0.80: the adult profile is the pubertal profile at just over half amplitude, not a different
profile. Consistent with this, the rank of every Myc-driven programme was preserved or sharpened
rather than lost -- the normalised enrichment score for the OXPHOS complexes rose from 2.81 to 3.34
and for the core Myc signature from 2.96 to 3.50 -- and not one gene in the genome reached
significance on the genotype-by-timepoint interaction. The attenuation is a change of amplitude,
spread evenly across the programme, and it is a property of the distribution rather than of any set
of genes. **[panel]**

The reason is not transcriptional. The MMTV-Myc transgene message itself does not attenuate at all:
`Myc` sits at the 100th percentile of the genome for retention of its own genotype effect, which is
marginally larger at twelve weeks (+1.82) than at six (+1.68). What falls is the protein. MYC is
reduced by approximately 50% in the adult gland by immunoblot, and a 50% reduction in effective MYC
is exactly what a uniform x0.55 rescaling of the downstream programme reports. These are one
measurement made twice -- once at protein level and once across eight thousand transcripts.
**[panel]** Everything else we looked for was absent: no Myc-antagonising repressor was induced, no
E-box competitor moved, and no Myc-independent transcriptional axis absorbed the change. There was
nothing to find, because output per unit MYC is unchanged.

The same is true inside the mitochondrial compartment, which carries the largest Myc effect in the
dataset and is the most tightly co-regulated network in the library. At six weeks Myc raises 94% of
the 143 mitochondrial pathways (median +0.43 log2) and, within the compartment, reprioritises them:
on a pairwise ratio measure that is blind to how many mitochondria a cell contains, protein import
(+0.106), the respiratory chain (+0.052) and mitochondrial central dogma (+0.036) are promoted,
while dynamics and surveillance (-0.182), signalling (-0.100) and small-molecule transport (-0.092)
are demoted. The apoptotic machinery sits inside the demoted arm (Apoptosis -0.247, padj 0.043),
with its pro- and anti-apoptotic branches moving together so that the balance itself is unchanged.
By twelve weeks that entire structure -- promotion and demotion alike -- is reproduced at reduced
amplitude, the priority profile regressing on itself with a slope of 0.64 and R2 up to 0.95 for the
respiratory arm. Myc's mitochondrial programme is rescaled, not reshaped. **[panel]**

If nothing in the Myc effect changes shape across this window, then any change in the biology of the
gland has to come from what Myc is acting on. That is what the wild-type timeline shows.

## Section 2 -- "The maturing gland withdraws the respiratory chain, and Myc loses the apoptotic trigger"

Between six and twelve weeks the wild-type mammary epithelium is not doing a smaller version of
anything. It withdraws one arm of the mitochondrial compartment specifically. The structural
subunits of the respiratory complexes fall by -0.255, at the 0th percentile against
expression-matched random gene sets, and they fall complex by complex (complex I -0.253, complex
III -0.251, complex IV -0.411, complex V -0.216), together with the import channels that deliver
them (TOM -0.169, TIM23 -0.202, MIA40 -0.273). The assembly factors of those same complexes do not
move at all (+0.001, 50th percentile), so this is not a loss of mitochondrial mass and not a general
contraction of the organelle. Nor is it a withdrawal from oxidative metabolism: over the same window
the gland *raises* the machinery that supplies electrons to the chain -- fatty-acid oxidation +0.227,
the carnitine shuttle +0.220, branched-chain amino-acid dehydrogenase +0.266, pyruvate metabolism
+0.237, the glycine cleavage system +0.240 -- while lowering the chain itself. The maturing gland
keeps and expands everything that feeds the respiratory chain and stands the chain down.
**[panel]** The wild-type temporal axis is batch-confounded in this design and we describe rather
than claim it; the mitigation is that the mitochondrial protein blots move the same way, off the RNA
batch entirely.

The obvious reading -- that the adult gland simply proliferates less -- does not survive. The pooled
proliferative programme falls by -0.047, a fifth of the respiratory fall, and a paired null computed
on the *difference* between the two puts the dissociation at the 0th percentile. Individual markers
agree (Mki67 -0.04, Ccnb1 +0.15, Plk1 +0.11, Aurka +0.31; only Top2a and E2f1 fall). The gland
de-respires without de-proliferating. What it loses instead is the terminal end bud, the compartment
that clears its own lumen by apoptosis during pubertal elongation: the TEB-versus-ductal signature is
the largest single mover in the whole comparator (-0.42). The physiological reading is morphogenetic,
not proliferative -- the pubertal gland keeps a death-ready mitochondrion because duct cavitation
requires apoptosis, and stands it down when morphogenesis is complete.

Against that background Myc goes on building a death-competent mitochondrion at both ages. It raises
mitochondrial content by 21-27%, induces the intermembrane-space protease Htra2 more strongly than
any other apoptotic gene in the dataset (+1.49, padj 2.6e-08) and Bax (+0.48), and lowers the
mitochondrial brake Bcl-xL (-0.42) -- and it does all of this at twelve weeks as well, at the x0.55
amplitude. What it loses is the trigger. Myc's induction of Puma relative to Bcl-xL retains -0.09 of
its pubertal value: it does not fade, it reverses. In the same animals, sharing the same denominator,
Bax:Bcl-xL retains exactly 0.55 -- the global rate -- so the panel contains its own control, and the
objection that this is merely the x0.55 attenuation is answered from inside it. Adjusting for
epithelial and immune composition strengthens rather than weakens the interaction (p 0.017), and an
independent genome-wide statistic that never touches the Bcl-xL denominator places `Bbc3` at the
0.51st percentile of 8,774 Myc-responsive genes. The gland does not achieve this by buffering: no
anti-apoptotic transcript moves across the window (Bcl2 padj 0.45, Bcl2l1 0.79, Mcl1 0.44). What is
withdrawn is the state the death machinery depends on, not the machinery. **[panel]**

Is the loss of killing simply the loss of MYC? The transcriptome says not -- everything else retains
about 0.55 while Puma priming reverses, and Bim, the other BH3-only protein induced by PGC-1a in
cells, sits at the opposite extreme (91.6th percentile). The decisive test is the inducible model.
Acute MYC activation in the adult gland, at equal or higher dose and in animals that had never
previously been exposed to MYC, produces approximately 80% less death against only 30% less
proliferation. Selection is excluded by design, and what changed is the tissue. This also resolves
the apparent paradox of the timeline. Transformation accelerates at twelve weeks, when the Myc
transcriptional effect is at half strength, because far less of the expansion is being cancelled --
and because the part of the programme that does not attenuate is the tumour-relevant part: the core
Myc target sets retain 0.74-0.77 of their pubertal effect against 0.50-0.58 for the mitochondrial
arms, the human MYC-associated mitochondrial bicluster signatures retain 0.91-0.98 while
proliferative programmes retain 0.63-0.87, and what loses its Myc-inducibility genome-wide is
luminal lineage identity rather than the Myc programme itself. None of these differences reaches
significance at n = 24 and we report them as a consistent direction, not as a result.

The two things that decline together over this window are the respiratory arm and the Puma trigger.
In the one space where a coupling between them is measurable -- a within-compartment ratio, which
cancels the shared change in mitochondrial content -- they co-vary (design-adjusted Spearman +0.70
against a null of 0.48, with the redox arm as a negative control that fails), and respiratory
priority moderates Myc's effect on Puma priming at the 92nd percentile of a within-timepoint
permutation null while that control axis sits at the 51st. Nothing here is significant at n = 24,
and the two candidate second inputs -- respiratory state and the terminal end bud -- are not
separable in a two-timepoint cross-sectional design. This is where the mouse hands over. If
respiratory capacity is the second input that Myc-driven death requires, then supplying it should
restore the trigger and removing it should remove the death at constant MYC. Only a perturbation can
test that.

---

# ORDER B -- the paradox first

*Question-driven. Open on the puzzle, let the dose result earn its place as the answer to half of
it, and use the other half to motivate the substrate.*

## Section 1 -- "Tumourigenesis accelerates while the Myc transcriptional effect halves"

Transformation in this model accelerates over exactly the window we sampled: the hyperplastic
expansion of the ductal epithelium is progressive, and by twelve weeks the Myc+ gland is
substantially further along it than at six. The transcriptome does the opposite. Between six and
twelve weeks the genotype effect contracts on every summary we computed. The number of
differentially expressed genes falls from 2,777 to 239 at 10% FDR, and from 1,967 to 135 at 5%.
Taken alone that would suggest the oncogenic programme had largely switched off. It has not: a
differential-expression count is a threshold statistic and exaggerates a uniform reduction in effect
size, and the genes still clearing the threshold at twelve weeks are simply the largest ones.
Measured directly, the Myc effect at twelve weeks regresses on the Myc effect at six weeks with a
slope of 0.55, an intercept of zero and R2 = 0.80 across the 143 nuclear-encoded mitochondrial
pathways. The ranking of the programme is not lost but sharpened -- normalised enrichment for the
OXPHOS complexes rises from 2.81 to 3.34 and for the core Myc signature from 2.96 to 3.50 -- and no
gene in the genome reaches significance on the genotype-by-timepoint interaction. The programme did
not change. Its amplitude did. **[panel]**

The cause is not transcriptional. The transgene message does not attenuate: `Myc` sits at the 100th
percentile of the genome for retention of its own effect. What falls is the protein, by
approximately 50% on immunoblot, and a 50% loss of effective MYC is precisely what a uniform x0.55
rescaling of eight thousand downstream transcripts reports. **[panel]** No repressor was induced, no
E-box competitor moved, and no Myc-independent axis absorbed the change; there was nothing to find,
because output per unit MYC is unchanged.

That leaves a genuine paradox, because a halved oncogenic input predicts less of everything,
including less transformation. Two things reconcile it, and only one of them is about Myc. The
first is that the attenuation, while close to uniform, is not indifferent to what it attenuates:
the core Myc target programme retains 0.74-0.77 of its pubertal effect against 0.50-0.58 for the
mitochondrial arms; the human MYC-associated mitochondrial bicluster signatures retain 0.91-0.98
while proliferative programmes retain 0.63-0.87; and genome-wide, what loses its Myc-inducibility
across the window is luminal lineage identity, while the Myc programme itself is among the
best-retained. None of this is significant at n = 24 and we report it as a direction rather than a
result, but its direction is consistent: the part of the programme that survives the dose fall is
the tumour-associated part. The second, and the subject of the rest of this section, is that far
less of the expansion is being cancelled. Myc-induced apoptosis, prominent in the pubertal gland, is
largely gone by twelve weeks. So the question is not why the Myc effect gets smaller -- that is dose
-- but why the same oncogene stops killing.

## Section 2 -- "The adult gland withdraws the respiratory chain, and with it Myc's apoptotic trigger"

Myc's largest and most coherent effect in this dataset is on mitochondria, and it is not simply an
increase. At six weeks Myc raises 94% of the 143 mitochondrial pathways (median +0.43 log2) and
mitochondrial content by 21-27%, but within the compartment it reprioritises: on a pairwise ratio
measure blind to how many mitochondria a cell contains, protein import (+0.106), the respiratory
chain (+0.052) and mitochondrial central dogma (+0.036) are promoted while dynamics and surveillance
(-0.182), signalling (-0.100) and small-molecule transport (-0.092) are demoted. The apoptotic
machinery sits in the demoted arm (-0.247, padj 0.043), with both branches moving together so the
balance is unchanged. By twelve weeks that whole structure is reproduced at reduced amplitude
(priority slope 0.64, R2 up to 0.95 for the respiratory arm). Nothing that Myc does to mitochondria
changes shape between the two ages, so if the biology changes, the change is in the substrate.
**[panel]**

And the substrate changes. Between six and twelve weeks the wild-type epithelium withdraws one arm
of the mitochondrial compartment specifically: the structural subunits of the respiratory complexes
fall by -0.255, at the 0th percentile against expression-matched random sets, complex by complex
(CI -0.253, CIII -0.251, CIV -0.411, CV -0.216), together with the import channels that deliver them
(TOM -0.169, TIM23 -0.202, MIA40 -0.273). The assembly factors of those same complexes do not move
(+0.001, 50th percentile), so this is not a loss of mitochondrial mass. Nor is it a withdrawal from
oxidative metabolism: the gland simultaneously raises everything that supplies electrons to the
chain -- fatty-acid oxidation +0.227, the carnitine shuttle +0.220, branched-chain amino-acid
dehydrogenase +0.266, pyruvate +0.237, the glycine cleavage system +0.240. It keeps the fuel and
stands down the chain. **[panel]** These wild-type temporal statements are batch-confounded and are
described rather than claimed; their mitigation is that the mitochondrial protein blots move the
same way, off the RNA batch entirely.

Nor is this a proliferative exit. The pooled proliferative programme falls by -0.047, a fifth as
much, and a paired null on the difference between the two places the dissociation at the 0th
percentile; individual markers agree (Mki67 -0.04, Ccnb1 +0.15, Plk1 +0.11, Aurka +0.31). The gland
de-respires without de-proliferating. What it loses is the terminal end bud (-0.42, the largest
mover in the comparator) -- the compartment that clears its own lumen by apoptosis during pubertal
elongation. The withdrawal is morphogenetic, and it ends when morphogenesis does.

Against that background Myc continues to build a death-competent mitochondrion at both ages, raising
Htra2 (+1.49, padj 2.6e-08) and Bax (+0.48) and lowering Bcl-xL (-0.42) at twelve weeks as well as
at six, at the x0.55 amplitude. What it loses is the trigger. Myc's induction of Puma relative to
Bcl-xL retains -0.09 of its pubertal value -- it reverses -- while Bax:Bcl-xL, in the same animals
and over the same denominator, retains exactly 0.55, the global rate. The panel therefore contains
its own control. Composition adjustment strengthens the interaction (p 0.017), and an independent
genome-wide statistic that never touches the Bcl-xL denominator puts `Bbc3` at the 0.51st percentile
of 8,774 Myc-responsive genes, while Bim sits at the 91.6th. The gland does not buffer: no
anti-apoptotic transcript moves across the window (Bcl2 padj 0.45, Bcl2l1 0.79, Mcl1 0.44). What is
withdrawn is the state the death machinery depends on, not the machinery. **[panel]**

That the killing is lost for a reason other than the fall in MYC is settled by the inducible model
rather than by the transcriptome. Acute MYC activation in the adult gland, at equal or higher dose
and in animals never previously exposed to MYC, gives approximately 80% less death against only 30%
less proliferation. Selection is excluded by design; what changed is the tissue. The respiratory arm
and the Puma trigger decline together over this window, and in the one space where a coupling
between them can be measured -- a within-compartment ratio, which cancels the shared change in
content -- they co-vary (design-adjusted Spearman +0.70 against a null of 0.48; the redox control
axis fails), with respiratory priority moderating Myc's effect on Puma priming at the 92nd
percentile of a within-timepoint permutation null against a control at the 51st. Nothing is
significant at n = 24, and respiratory state and the terminal end bud are not separable in this
design. If respiratory capacity is the second input that Myc-driven death requires, then supplying
it should restore the trigger and removing it should abolish the death at constant MYC. That is the
prediction the cell and inducible-model experiments test.

---

# ORDER C -- two actors  *(recommended)*

*Agent-partitioned. Section 1 is the oncogene; section 2 is the tissue. The paradox sits on the
hinge between them.*

## Section 1 -- "Myc reprioritises the mitochondrial compartment, and the whole programme is scaled down in the adult gland"

The MMTV-Myc gland expands progressively across this window and is substantially more hyperplastic
at twelve weeks than at six. As we show below, the Myc transcriptional effect at twelve weeks is
about half what it is at six. We first establish what that effect is, and then what happens to it.

Myc's largest and most coherent effect in these mammary epithelial cells is on mitochondria, and it
is not simply an increase in mitochondrial biogenesis. At six weeks Myc raises 94% of the 143
nuclear-encoded mitochondrial pathways (median +0.43 log2) and mitochondrial content by 21-27%, but
within the compartment the arms are ordered. On a pairwise ratio measure that is blind to how many
mitochondria a cell contains, protein import (+0.106), the respiratory chain (+0.052) and
mitochondrial central dogma (+0.036) are promoted, while dynamics and surveillance (-0.182),
signalling (-0.100) and small-molecule transport (-0.092) are demoted. The apoptotic machinery sits
inside the demoted arm (Apoptosis -0.247, padj 0.043; pro-apoptotic -0.210, anti-apoptotic -0.273),
and because both branches move together the apoptotic balance itself is unchanged -- Myc de-emphasises
the death machinery relative to the rest of the compartment without shifting it towards survival.
**[panel]**

That compartment-level reprioritisation is what changes least across the window, and what changes is
its amplitude. Between six and twelve weeks the number of differentially expressed genes falls from
2,777 to 239 at 10% FDR and from 1,967 to 135 at 5%. The gene counts overstate it: a count is a
threshold statistic and exaggerates a uniform reduction in effect size, so that the genes still
clearing the threshold at twelve weeks are simply the largest. Measured directly, the Myc effect at
twelve weeks regresses on the Myc effect at six weeks with a slope of 0.55, an intercept
indistinguishable from zero and R2 = 0.80. The ranking of the programme is preserved or sharpened
rather than lost -- normalised enrichment for the OXPHOS complexes rises from 2.81 to 3.34 and for
the core Myc signature from 2.96 to 3.50 -- and no gene in the genome reaches significance on the
genotype-by-timepoint interaction. The programme did not change; its amplitude did. **[panel]**

The cause is not transcriptional. The transgene message does not attenuate at all: `Myc` sits at the
100th percentile of the genome for retention of its own genotype effect, marginally larger at twelve
weeks (+1.82) than at six (+1.68). What falls is the protein, by approximately 50% on immunoblot,
and a 50% reduction in effective MYC is exactly what a uniform x0.55 rescaling of the downstream
programme reports; these are one measurement made twice. **[panel]** No Myc-antagonising repressor
was induced, no E-box competitor moved and no Myc-independent axis absorbed the change -- there was
nothing to find, because output per unit MYC is unchanged.

The mitochondrial compartment is no exception, and that is worth stating positively rather than as
an absence. The promotion-and-demotion structure described above is reproduced at twelve weeks in
full, at reduced amplitude: the priority profile regresses on itself with a slope of 0.64, and with
R2 up to 0.95 for the respiratory arm. Myc's mitochondrial programme is rescaled, not reshaped, and
the same holds arm by arm, with one departure -- the core Myc target sets retain 0.74-0.77 of their
pubertal effect against 0.50-0.58 for the mitochondrial arms, so the parts of the programme most
directly bound by Myc fade least. **[panel]**

This leaves the timeline paradoxical. Transformation accelerates over exactly the interval in which
the oncogenic input halves, and a lower dose predicts less of everything, including less
transformation. Two observations bear on it. The first is weak but consistent in direction: what
survives the dose fall is disproportionately the tumour-relevant part of the programme -- besides
the core Myc targets, the human MYC-associated mitochondrial bicluster signatures retain 0.91-0.98
of their pubertal genotype difference while proliferative programmes retain 0.63-0.87, and
genome-wide, what loses its Myc-inducibility is luminal lineage identity rather than the Myc
programme. None of this reaches significance at n = 24. The second is decisive and is the subject of
the next section: far less of the expansion is being cancelled, because Myc-induced apoptosis, which
is prominent in the pubertal gland, is largely gone by twelve weeks. Since nothing Myc does changes
shape across the window, the explanation has to lie in what Myc is acting on.

## Section 2 -- "The maturing gland withdraws the respiratory arm, and Myc's apoptotic trigger goes with it"

Between six and twelve weeks the wild-type mammary epithelium is not doing a smaller version of
anything -- it is doing something different, and to one arm of the compartment in particular. The
structural subunits of the respiratory complexes fall by -0.255, at the 0th percentile against
expression-matched random gene sets, complex by complex (complex I -0.253, complex III -0.251,
complex IV -0.411, complex V -0.216), together with the import channels that deliver them (TOM
-0.169, TIM23 -0.202, MIA40 -0.273). The assembly factors of those same complexes do not move at all
(+0.001, 50th percentile), so this is not a loss of mitochondrial mass. Nor is it a withdrawal from
oxidative metabolism: over the same window the gland raises the machinery that supplies electrons to
the chain -- fatty-acid oxidation +0.227, the carnitine shuttle +0.220, branched-chain amino-acid
dehydrogenase +0.266, pyruvate metabolism +0.237, the glycine cleavage system +0.240. The maturing
gland keeps and expands the fuel supply and stands down the chain itself. **[panel]** The wild-type
temporal axis is batch-confounded in this design and we describe rather than claim it; the mitigation
is that the mitochondrial protein blots move the same way, off the RNA batch entirely.

It is also not a proliferative exit, which matters because that is the reading the model cannot
afford. The pooled proliferative programme falls by -0.047, a fifth of the respiratory fall, and a
paired null computed on the difference between the two puts the dissociation at the 0th percentile;
individual markers agree (Mki67 -0.04, Ccnb1 +0.15, Plk1 +0.11, Aurka +0.31, with only Top2a and
E2f1 falling). The gland de-respires without de-proliferating. What it loses is the terminal end bud
-- the single largest mover in the comparator (-0.42) -- which is the compartment that clears its own
lumen by apoptosis during pubertal elongation. The withdrawal is morphogenetic rather than
proliferative, and it ends when morphogenesis ends.

Against that changing background Myc goes on building a death-competent mitochondrion at both ages.
It induces the intermembrane-space protease Htra2 more strongly than any other apoptotic gene in the
dataset (+1.49, padj 2.6e-08) and Bax (+0.48), and lowers the mitochondrial brake Bcl-xL (-0.42) --
at twelve weeks as well as at six, at the x0.55 amplitude. What it loses is the trigger. Myc's
induction of Puma relative to Bcl-xL retains -0.09 of its pubertal value: it does not attenuate, it
reverses. In the same animals, over the same denominator, Bax:Bcl-xL retains exactly 0.55 -- the
global rate -- so the objection that this is merely the x0.55 dose effect is answered from inside the
panel. Adjusting for epithelial and immune composition strengthens the interaction rather than
weakening it (p 0.017), and an independent genome-wide statistic that never touches the Bcl-xL
denominator places `Bbc3` at the 0.51st percentile of 8,774 Myc-responsive genes, while Bim -- the
other BH3-only protein induced by PGC-1a in cells -- sits at the 91.6th. Nor does the gland solve the
problem by buffering: no anti-apoptotic transcript moves across the window (Bcl2 padj 0.45, Bcl2l1
0.79, Mcl1 0.44). **What is withdrawn is the state the death machinery depends on, not the
machinery.** **[panel]**

Whether the loss of killing is simply the loss of MYC is not something the transcriptome can settle,
and the inducible model does. Acute MYC activation in the adult gland, at equal or higher dose and in
animals that had never previously been exposed to MYC, produces approximately 80% less death against
only 30% less proliferation. Selection is excluded by the design, so the change is a property of the
tissue and not an adaptation of the transformed cell. The asymmetry is itself informative:
proliferation requires one input and death requires two, so a fall in the second costs death
superlinearly and proliferation linearly.

The two things that decline together over this window are therefore the respiratory arm and the Puma
trigger. In the one space where a coupling between them is measurable -- a within-compartment ratio,
which cancels the shared change in mitochondrial content -- they co-vary (design-adjusted Spearman
+0.70 against a null of 0.48, with the redox arm as a negative control that fails), and respiratory
priority moderates Myc's effect on Puma priming at the 92nd percentile of a within-timepoint
permutation null while the control axis sits at the 51st. Nothing here is significant at n = 24, and
respiratory state and the terminal end bud are not separable in a two-timepoint cross-sectional
design; we treat the mouse as generating the hypothesis rather than testing it. If respiratory
capacity is the second input that Myc-driven death requires, then supplying it should restore the
trigger and removing it should abolish the death at constant MYC. Those are the predictions the cell
and inducible-model experiments were designed to test.

---

# Comparison

| | **A -- subtract the dose** | **B -- paradox first** | **C -- two actors** |
|---|---|---|---|
| §1 opens on | the DEG collapse and the x0.55 rescaling | hyperplasia at 12W | Myc's mitochondrial reprioritisation |
| the paradox sits | at the end of §2, as a resolution | at the top of §1, as the question | on the hinge between §1 and §2 |
| §1 / §2 division | control / finding | puzzle / answer | **oncogene / tissue** |
| headings already in the Doc | need re-titling | need re-titling | **fit almost verbatim** |
| protagonist | shifts from dose to tissue | shifts from phenotype to tissue | Myc, then the gland -- one handover |
| strongest asset | reviewer-proof; every later number is read against an established baseline | the best opening sentence in the material | opens on your own novel result, not on a control |
| main cost | a page of control before any biology; the best hook is spent last | opens on histology/IHC and leans on the weakest evidence early | dose arrives second, so one paragraph could be misread as covering both ages |

**Recommendation: C, with one sentence of B.** C is the only order whose two sections each carry
exactly one claim, and it is the only one that does not require re-titling what is already written.
The signpost sentence at the head of §1 -- *the gland is substantially more hyperplastic at twelve
weeks, when the Myc effect is at half strength* -- borrows B's hook without making the weak
retained-signature evidence load-bearing, because the paradox is not developed until the hinge,
where the death arm is already in view.

**If you prefer B**, one thing has to be fixed first: paragraph 1 as written says the model "showed
hyperplastic expansion of the developing ductal epithelial structures **at these time points**",
i.e. at both ages. B's opening sentence needs the histology stated as progression (or quantified),
or the paradox has nothing to stand on.

---

# Answers to the four questions in the list

**1. "Is the attenuation all proportional in the mitochondria, the most coupled gene network?"**
**Yes -- and this is a positive result, not a null.** Content slope 0.552 (R2 0.80, rho 0.84);
priority slope 0.644 (R2 0.79, rho 0.93); per-tier priority R2 up to 0.945 for OXPHOS. Per-tier
content slopes run 0.53-0.83. The Myc mitochondrial programme is *rescaled without being reshaped*.
There are exactly two departures: the core Myc targets fade least (0.74-0.77 against 0.50-0.58), and
at gene level `Bbc3` alone reverses. Saying the attenuation is proportional **strengthens** the
section, because it is what makes the two non-proportional things -- the wild-type withdrawal and the
Puma reversal -- stand out as findings rather than as noise.

**2. "Tumourigenesis takes off at 12W -- any evidence in the transcriptome, or did it already take
off at 6W?"** The primary evidence is histological and IHC, and the paper should say so. The
transcriptome adds one weak but internally consistent line: **what fails to attenuate is the
tumour-relevant part of the programme.** Core Myc targets retain 0.74-0.77 against 0.50-0.58 for the
mitochondrial arms; the METABRIC MYC-associated mitochondrial bicluster signatures retain 0.91-0.98
of their pubertal genotype difference while proliferative programmes retain 0.63-0.87 and
stem/alveolar sets 0.12-0.66; the MB2-over-MB1 fork genotype gap widens from 0.032 to 0.133; and the
genome-wide collapse scan finds the Myc programme retained while luminal identity collapses
(`MG_HS_GRAY` NES -2.79, padj 1.7e-17). **None of it is significant** -- the fork interaction is
p 0.15, and the bicluster retentions are a ratio of group-mean GSVA differences with no null behind
them. One sentence, framed as a direction. And no: it has not already happened at 6W in the
transcriptome. What is different at 6W is how much of the expansion death is cancelling.

**3. "The WT timeline: metabolism up -- is it anabolic?"** **No -- the opposite, and the correct
answer is better for the argument.** What rises in the wild-type gland is fuel catabolism and
substrate oxidation: fatty-acid oxidation +0.227, carnitine shuttle +0.220, branched-chain
amino-acid dehydrogenase +0.266, pyruvate +0.237, glycine cleavage +0.240. What falls is the
respiratory chain itself, subunit class by subunit class. **The maturing gland raises everything
that feeds electrons into the chain and withdraws the chain.** That is a sharper statement than
"metabolism up", and it reinforces the specificity claim -- this is not a general mitochondrial
contraction, and it is not a shift out of oxidative metabolism.

**4. "Collapses -- find a better synonym."** The measurement is a **sign reversal**, not a decline
to zero (`Puma:Bcl-xL` retention -0.09 against a global rate of 0.55). Best single word:
**reverses**. Also accurate: *inverts*; *is lost outright*; *does not merely attenuate -- it
disappears and turns over*. Avoid *collapses* only because it is silent about direction, and the
sign change is the whole point -- it is why the result beats its null.

---

# Traps to keep out of the text

1. **Do not lead on the DEG count.** 1,967 -> 135 is -93%, but the effect size falls only ~45%. The
   median absolute fold change among the survivors actually *rises* (0.678 -> 0.869), for exactly the
   threshold reason. Lead with the regression; give the counts as illustration with the caveat
   stated. All three narratives above do this.
2. **The genome-wide interaction contrast has zero genes at 5% or 10% FDR.** "The Myc effect
   attenuates" is a set-level, distributional statement and must never be written as a per-gene
   significant result -- including for `Bbc3`, whose interaction is genome-wide null.
3. **"The ranking showed no change" understates it.** fGSEA NES *rises* at 12W for every arm. Write
   *preserved or sharpened*.
4. **"The only arm the gland withdraws from" is not exact.** The protein-import channels also fall
   (TOM -0.169, TIM23 -0.202, MIA40 -0.273; tier median -0.056). Use *by far the steepest, with the
   import machinery a distant second* -- which is what the narratives do by naming both.
5. **The Myc+ gland's own OXPHOS decline is mostly the dose effect** and is not evidence for the
   substrate model. The two legs that are evidence are the **wild-type** withdrawal (no Myc to lose)
   and the **Puma reversal** (the one Myc-dependent measurement departing from the x0.55 rate). The
   summary sentence in the list -- "Myc acting on a different background" -- should be built from
   those two, not from the Myc+ OXPHOS trace.
6. **`batch = timepoint`.** Every wild-type temporal statement is described, not claimed; the
   mitigation is the protein blots, off the RNA batch. Stated once per narrative above.
7. **The epistemic contract.** Every causal verb sits in a cell-experiment sentence. The mouse
   supplies the contrast and the co-decline; the perturbations supply causality.
8. **Priming is not death.** The transcriptome measures the molecular substrate; the phenotype is
   IHC and the MYC-ER counts. There is also no efferocytosis signature at 6W (`Gas6` is down), so
   the death is invisible to bulk RNA by construction.
9. **Housekeeping in the Doc.** The partial-correlation paragraph ("only select mitochondrial
   pathways, such as OXPHOS, nucleotide metabolism and the Krebs cycle...") currently sits orphaned
   under the second heading. It is paragraph 2's closing sentence -- the queued **Fig. 1E** -- and
   belongs back under the first heading, otherwise whichever order is chosen opens on a topic it
   does not develop. Separately, its "central coupling" claim is in tension with scripts 35/36:
   OXPHOS partial -> proliferation is +0.59 against a partial ceiling of 0.45 (68th percentile,
   p 0.32), and script 36 finds the centrality claim **untestable** rather than false, because
   OXPHOS *is* the global axis (r = 0.968 with the per-sample global mean). That is the Fig. 1E
   conversation, not this one.

---

# The hand-off, common to all three

All three narratives end on a prediction rather than a conclusion: *if respiratory state is the
second input, supplying it should restore the trigger.* The next section then tests it in both
directions -- PGC-1a raises OXPHOS and PUMA and kills on a Myc background, and once execution is
blocked the same intervention becomes a growth advantage. That sign reversal is the strongest single
fact in the paper, and it belongs to the cells, not to the mouse. The canonical-literature inversion
(PGC-1a as cytoprotective, from the muscle-atrophy and neurodegeneration contexts) should be
discharged there or early in the Discussion, not left for a reviewer to supply.
