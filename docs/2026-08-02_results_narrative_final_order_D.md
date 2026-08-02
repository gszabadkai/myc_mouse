# The final narrative -- Order D, with its panels

**Status:** 2026-08-02. This supersedes `docs/2026-08-01_results_narrative_three_orders.md`,
which stays as the record of how the choice was made.

The author read the three orders and wrote a fourth, tighter spine: four beats plus two open
mechanistic questions. **Order D is Order C's skeleton with B's hinge**, plus two improvements
that are the author's own -- the apoptosis observation moves *into* section 1 as part of the
reprioritisation list, so "competent but not more primed" is established before the window
argument needs it; and the mechanism questions are made explicit at the end rather than buried.

**What is new in this document, beyond the prose.** Every assertion in the author's spine was
checked against the objects on disk. Three of them needed changing, and two of the three changes
make the argument stronger. Six things are missing from the spine and are added. Both open
questions are answered as far as the data and the literature reach. The panel list is at the foot.

**Provenance.** `interaction_results.rds` (raw, unshrunken), `background_vs_myc.rds` (40),
`attenuation_decomposition.rds` (29), `priming_arm_teb.rds` (42),
`substrate_specificity_tradeoff.rds` (43), `collapse_module_ownership.rds` (44),
`mito_content_proxies.rds` (32), `mitopps_priming_pgc1a.rds` (38). Blot, IHC and MYC-ER values are
the author's, reported as given. No script was re-run. Where a panel is needed I have written
`[2A]`, `[3C]` and so on, keyed to the table at the foot.

---

# Section 1 -- "Myc reprioritises the mitochondrial compartment, and the whole programme is scaled down by dose in the adult gland"

Myc does two separable things to the mitochondrial compartment, and they need two rulers to see.
The first is quantitative. Across the 143 nuclear-encoded mitochondrial pathways the transgene
raises set-average expression in 94% of them, with a median effect of +0.43 log2, and the absolute
mitochondrial share of the transcriptome rises with it -- by 21 to 27% depending on the proxy, an
effect that survives adjustment for dissociation stress and for residual non-epithelial
contamination (adjusted p = 0.0006) `[2B]`. The gland is not merely re-labelling its
mitochondria; it is building more of them.

The second is a change in the internal priorities of the compartment, and to see it the content
increase has to be removed. On the pairwise ratio ruler, which normalises each pathway against the
rest of the compartment and is therefore blind to how much mitochondrion there is, the Myc effect
resolves into a clear reallocation `[2A]`. Promoted: protein import, sorting and homeostasis
(+0.106 at tier level, +0.117 for the 87-gene set, padj 0.052); mitochondrial translation and the
mitochondrial central dogma (mitoribosome +0.121; central dogma +0.092, padj 0.044); the oxidative
phosphorylation subunits (+0.112, led by complex V at +0.165 and complex III at +0.160); and, at
the level of individual pathways, one-carbon and amino-acid metabolism most strongly of all --
glycine cleavage +0.41, creatine metabolism +0.37, pyruvate metabolism +0.33, S-adenosylmethionine
+0.30, serine +0.21. Demoted: mitochondrial dynamics and surveillance as a whole (-0.182), and
within it fission (-0.202, padj 0.016), autophagy (-0.191), mitophagy (-0.182) and apoptosis
(-0.247, padj 0.043); mitochondrial signalling (-0.100); small-molecule transport (-0.092); the
calcium uniporter (-0.306); and the thirteen mitochondrially encoded subunits (-0.341).

Two features of that reallocation are worth stating because they recur. First, it is the OXPHOS
*subunits* that are promoted and not the *assembly factors* of the same complexes, which sit at
-0.028; the same subunit-versus-assembly split reappears in the wild-type timeline below, so one
sub-compartment carries both halves of this account. Second, the promotion of respiration comes
with the promotion of everything needed to build a respiring organelle -- import channels,
mitochondrial ribosomes, the mitochondrial central dogma -- which is what a biogenesis programme
looks like, rather than a selective induction of the chain.

Apoptosis is in the demoted list, and the way it is demoted matters. Both arms move together and
in the same direction: the pro-apoptotic MitoCarta set falls by 0.210 (padj 0.022) and the
anti-apoptotic set by 0.273 (padj 0.088). The *balance* between them therefore does not move.
Myc is not raising apoptotic priming at six weeks; it is holding the apoptotic machinery at a
constant relative weight while it expands everything else. And it builds an organelle that is
nonetheless fully competent to die: HTRA2 rises by 1.49 log2 (padj 2.6e-8), BAX by 0.48
(padj 1.4e-4) and BCL-xL falls by 0.42 (padj 0.040) -- and cell death is what the six-week gland
in fact shows by immunohistochemistry `[3B]`. Asset and liability are constructed in one move.

At twelve weeks the entire structure is present at half amplitude and unchanged in shape `[2C]`.
Regressing the Myc effect at twelve weeks on the Myc effect at six weeks across the 143 pathways
returns a slope of 0.552 with an intercept indistinguishable from zero and R2 = 0.80 (bootstrap
95% interval 0.514 to 0.605); on the priority ruler the slope is 0.644 with R2 = 0.79, and within
the OXPHOS tier alone R2 reaches 0.945. The coherence of that rescaling sits at the 100th
percentile of expression-matched shuffled sets, so it is not what an arbitrary collection of genes
would do. The ranking is not merely preserved but sharpened: the normalised enrichment score of
the mitochondrial OXPHOS set rises from 2.81 to 3.34 between the two ages, and of the Myc
integrative signature from 2.96 to 3.50. The adult programme is the pubertal programme, played
quieter.

*(The gene-level counts tell the same story and exaggerate it, which is why they are reported and
not led on. Differential expression falls from 1,967 genes to 135 at 5% FDR, and from 2,777 to 239
at 10%. A differential-expression count is a threshold statistic: under a uniform reduction in
effect size it falls far faster than the effect size does, and the survivors at twelve weeks are
simply the largest effects, so their median absolute fold change is higher than at six weeks --
0.869 against 0.678 -- not lower. The rescaling is a 45% reduction reported by a 93% one.)*

The reason is dose. MYC protein in the twelve-week gland is about half of what it is at six weeks
by western blot, and that is the primary evidence; everything else is consistency with it `[2D]`.
The transgene message does not fall -- the genotype gap is +1.68 log2 at six weeks (padj 2.7e-17)
and +1.82 at twelve (padj 2.1e-20), with no within-genotype decline across age (-0.10, padj 0.76);
measured as the ratio of its twelve-week to its six-week effect, *Myc* is the single most retained
gene in the transcriptome, at the 99.99th percentile of 8,774. A halving of protein predicts the
x0.55 rescaling that is observed. The shape is preserved on both rulers and at every tier. The
most retained gene sets on a ranking of departure from that rescaling are the Myc regulon itself
(Muhar +2.83, Hallmark MYC targets V2 +2.81, Felsher +2.81, all padj < 1e-8). And every
alternative has been scanned and is negative: no E-box binder changes, a 111-gene survey of
nuclear receptors, corepressors and the DREAM complex yields no significant interaction, no
Myc-independent axis absorbs the attenuation, and the post-translational MYC stability machinery
is not transcriptionally down -- what does move in it predicts a *more* stable protein, not a less
stable one.

Dose does not act linearly, and it is not expected to. Myc target promoters differ in affinity and
the high-affinity ones are already saturated at physiological Myc, so halving the protein cannot
halve every target equally (Lorenzin et al., eLife 2016). Our retention values order exactly as
that predicts: the Myc core retains 0.74 to 0.77 of its six-week effect while the mitochondrial
arms retain 0.50 to 0.58. Affinity accounts for everything that sits above the global line. It
cannot account for anything that crosses zero -- which is what makes the one gene that does worth
an experiment, and what section 2 turns on.

**The hinge.** If the oncogene's programme only shrinks -- same pathways, same ranking, less of it
-- then nothing in it explains a gland that becomes overtly hyperplastic at twelve weeks, when the
driver's transcriptional output is at half strength. Whatever changed is not in what Myc is doing.
It is in what Myc is doing it to.

---

# Section 2 -- "The maturing gland withdraws respiration, and Myc's apoptotic trigger goes with it"

Histologically the balance between proliferation and cell death shifts across this window. The
transcriptome cannot see the phenotype, but it can say which half of the balance is moving, and
the answer is unambiguous: the death half. Measured against expression-matched random sets, the
wild-type gland's six- to twelve-week change in proliferation signatures is -0.047, at the 1.3rd
percentile -- real, but five times smaller than its change in OXPHOS subunits, which is -0.255 at
the 0th percentile `[3A]`. Proliferation barely tracks the tumourigenic timeline. What changes
underneath it is the substrate for killing.

What the normal gland does to itself in this window is specific, and the specificity is what makes
it a finding rather than a general contraction. It withdraws from the respiratory chain -- OXPHOS
subunits at the 0th percentile of matched sets, complex by complex (CI -0.253, CIII -0.251,
CIV -0.411, CV -0.216) -- while the assembly factors of those same complexes do not move at all,
sitting at the 50th percentile (+0.0007), and the mitochondrial ribosome sits at the 26th. The
protein import channels fall too, though far less steeply (TOM -0.169, TIM23 -0.202, MIA40 -0.273).
And at the same time the gland *raises* fuel catabolism: amino-acid metabolism at the 100th
percentile (+0.184) and lipid metabolism at the 99.9th (+0.122), with fatty-acid oxidation +0.227,
the carnitine shuttle +0.220, branched-chain amino-acid dehydrogenase +0.266, pyruvate +0.237 and
glycine cleavage +0.240. The maturing gland raises everything that feeds electrons into the chain
and withdraws the chain. It is a catabolic shift, not an anabolic one, and it is the only
mitochondrial arm the normal gland withdraws from.

This is a different programme from Myc's, not its mirror. Across the 143 pathways the correlation
between the Myc reallocation and the wild-type temporal reallocation is -0.02, with half the
pathways discordant in sign -- which is what chance gives. The two reprioritisations are
independent. On the arms that matter to this paper they nonetheless point opposite ways: OXPHOS
subunits go from +0.112 under Myc to -0.110 with age, import from +0.117 to -0.056, the
mitochondrial ribosome from +0.121 to -0.041, and apoptosis from -0.247 to +0.049. Two programmes
acting on one compartment, and the tissue is moving the arm Myc most depends on.

*(A caveat that belongs here and is stated once. The two ages were extracted as two batches, so
batch is aligned with timepoint by design. Genotype is balanced within each batch, so every
genotype and interaction statement in this paper is clean; the pure between-age comparison is not
separable from batch. The wild-type temporal changes are therefore described rather than claimed
as developmental, and the mitigation is that the protein measurements move the same way, off the
RNA batch entirely.)*

The apoptotic balance still does not move. In the wild-type timeline the pro- and anti-apoptotic
arms rise together (+0.054 and +0.046), as they fell together under Myc. Neither the oncogene nor
the tissue changes the aggregate priming of the compartment. What changes is one gene.

Of the pro-apoptotic ratios that Myc establishes at six weeks, every one fades at the global rate
except one `[3B]`. BAX over BCL-xL is +0.90 at six weeks (p = 1.2e-4) and +0.50 at twelve, a
retention of 0.549 -- which is the global rescaling rate, arriving independently on a pair of
genes. BAK1 over BCL-xL retains 0.27, BID over BCL-xL 0.42. PUMA over BCL-xL is +0.67 at six weeks
(p = 0.0089) and **-0.06** at twelve: a retention of -0.09, a sign reversal, on the same
denominator as the internal control. Measured as departure from the x0.55 line across 8,774 genes,
*Bbc3* sits at the 0.51st percentile with z = -2.49. It is not a class effect: *Bcl2l11* (BIM)
sits at the 91.6th percentile and moves the other way. And it is not a background effect either --
*Bbc3* does not change in the wild-type gland across this window (+0.06, padj 0.88); it falls only
in the Myc+ gland (-0.48, padj 0.016), which is to say it is an interaction (p = 0.0081).

That is the structure of the result, and it is worth naming: **the respiratory withdrawal is
Myc-independent and the trigger loss is Myc-dependent.** The tissue does the first to itself; the
second happens only where the oncogene is. The phenotype needs both.

The gene that sits beside *Bbc3* on the same ranking is *Foxo3* `[3C]`. It is at the 0.456th
percentile with z = -2.53, immediately adjacent, and its interaction p-value (0.0074) is the
lowest of any gene in the death or biogenesis rosters. Its direction is the informative part: the
wild-type gland *raises* Foxo3 across the window (+0.47, padj 0.0078) while the Myc+ gland does
not (-0.01). FOXO3 is the canonical transcriptional activator of PUMA, binding and transactivating
the *BBC3* locus directly, and it is controlled post-translationally by AKT -- so a message-level
observation nominates it and cannot settle it. That is what the FOXO3 western across both ages,
and in the PGC-1a arm, is for.

None of this follows from the fall in Myc dose. A lower dose predicts less of everything, PUMA
included, in proportion -- and PUMA is the one measurement in the death arm that departs from
proportion. The decisive evidence is external to the transcriptome: in the MYC-ER model, at equal
or higher Myc dose and in animals never previously exposed to the oncogene, the adult gland shows
roughly 80% less death against roughly 30% less proliferation. The change between six and twelve
weeks is in the competence of the tissue to execute, not in the strength of the signal.

*(One statistical point, because it will be asked. The genome-wide interaction contrast has zero
genes at 5% or 10% FDR, and the attenuation is therefore reported as a distributional, set-level
result -- 143 pathways, slope 0.552, 95% interval 0.514 to 0.605 -- and never as a per-gene
significant interaction, *Bbc3* included, whose adjusted p is 0.84. Two things make that expected
rather than contradictory. The interaction is a difference of two differences at n = 6 per group:
its median standard error is 0.333 against 0.233 for a genotype contrast, a factor of 1.43, while
the effect it is asked to detect is smaller, so the detectable effect size is about three times
higher. And under a global rescaling every Myc-responsive gene has a non-zero interaction by
construction, so testing against zero measures the dose fade rather than the biology; the
biologically specific question is departure from the x0.55 line, which is how the PUMA and FOXO3
results above are computed. The PUMA test is licensed by pre-specification: *Bbc3* was named in
advance from the cell experiments, not selected from this scan.)*

**The hand-off, stated as a prediction rather than a conclusion.** If respiratory state is the
second input to the trigger, then supplying it should restore the trigger. The mouse can rank that
hypothesis but not test it: on the ratio ruler the interaction between Myc and OXPHOS priority on
PUMA-over-BCL-xL has p = 0.0052, while the same test on a redox axis that Myc does not drive gives
p = 0.72 `[3D]` -- a specificity contrast at n = 24, not a proof. What settles it is the
perturbation, in both directions: PGC-1a raises respiration and PUMA and kills on a Myc
background; and once execution is blocked, the same intervention becomes a growth advantage.

---

# The corrections

Three assertions in the author's spine needed changing. Two of the changes strengthen the
argument.

**1. PUMA is not deprioritised in the background.** The spine reads "the overall priming balance is
not changing, but OXPHOS and PUMA are deprioritised", with both attributed to the background. The
balance half is right on both contrasts. The OXPHOS half is right and is a wild-type statement. The
PUMA half is not a wild-type statement: *Bbc3* in wild type across the window is +0.06 (padj 0.88),
flat. It falls only in Myc+ (-0.48, padj 0.016), interaction p = 0.0081. **This is better than the
draft**: it makes the two inputs cleanly separable -- one from the tissue, Myc-independent; one
from the oncogene, an interaction -- which is the two-input model the paper argues, and it turns
the hand-off into a real prediction instead of a restatement.

**2. The background is not the mirror of Myc's reprioritisation.** Correlation across the 143
pathways is -0.02 (Spearman -0.03), 50.3% sign-discordant. The two reallocations are
*independent*, not opposite. The sign flips on OXPHOS, import, mitoribosome and apoptosis are real
but selecting those four is selection unless a null is quoted -- and the content ruler has one
(OXPHOS subunits at the 0th percentile of matched sets, assembly factors at the 50th). Write the
withdrawal from the matched-null test; write the priority flip as description.

**3. "The reprioritisation favours OXPHOS" is true at tier level and not at pathway level.** The
largest single promotions are one-carbon and amino-acid metabolism (glycine cleavage +0.41,
creatine +0.37, pyruvate +0.33). OXPHOS leads as a tier, not as a leaf, and it is the subunits
rather than the assembly factors.

# What was missing from the spine

1. **Mitochondrial content rises 21-27%** (script 32, adjusted p = 0.0006). The only solid new
   result in the corpus, and the ruler that makes "reprioritisation" mean what it says -- priority
   is a content-blind ratio, so without the content panel a reader takes reprioritisation for a
   content claim. Both rulers or neither.
2. **Foxo3 answers the author's own question and the answer was already on disk** -- 0.456th
   percentile, adjacent to *Bbc3* at the 0.513th, interaction p = 0.0074, raised by age in wild
   type and not in Myc+.
3. **Myc's repressive arm fades faster than its inductive arm** -- retention 0.33 against 0.51
   among confidently moved genes, while genome-wide there is no such asymmetry (0.445 against
   0.471). A property of the genes Myc actually moves.
4. **The largest single departure from the dose line in the whole transcriptome is the release of
   luminal suppression.** `MG_HS_GRAY` NES -2.79, padj 1.7e-17 over 299 genes: the Myc effect goes
   from -0.466 to +0.090. Note the direction -- this is not "luminal identity collapses", it is
   "Myc stops suppressing it". Partly the repressive-arm asymmetry above, but far beyond it in
   size. It fits section 2's account of what changes by twelve weeks, but it is a
   de-differentiation statement and Fig. 1B already owns that topic. **Flagged for the author's
   decision, not placed.**
5. **Bcl2l11/BIM is the negative control inside the BH3-only family** (91.6th percentile,
   retention 1.14). The pair splits; PUMA is not a class effect.
6. **The batch caveat**, once, as written into section 2.

# The dose evidence, as a list

The author asked for it explicitly.

1. The transgene message does not fall: gap +1.68 -> +1.82 log2, within-genotype -0.10 (padj 0.76),
   retention at the 99.99th percentile of 8,774 genes.
2. MYC protein falls about 50% by western blot. **This is the proof.**
3. -50% protein is quantitatively the fitted rescaling: slope 0.552, 95% interval 0.514-0.605,
   intercept zero.
4. The shape is preserved: R2 0.80 (content) and 0.79 (priority), rho 0.84 and 0.93, per-tier
   priority R2 up to 0.945; the R2 sits at the 100th percentile of matched nulls.
5. The ranking is preserved or sharpened, not degraded: NES rises at twelve weeks for every arm,
   and the most-retained sets on the departure ranking are the Myc regulon itself.
6. Every alternative mechanism has been scanned and is negative: no E-box binder moves; a 111-gene
   nuclear-receptor / corepressor / DREAM scan gives no significant interaction; no Myc-independent
   axis absorbs it; the post-translational stability machinery is not transcriptionally down and
   what moves in it predicts a more stable protein.
7. The Myc-specific temporal term is a uniform offset (0th percentile of nulls), not a different
   set of pathways.

# The two questions

## Why do PUMA and OXPHOS go together, when PUMA suppresses MPC?

The author's recollection of the literature is correct and the naive prediction does fail. PUMA
binds the mitochondrial pyruvate carrier, disrupts MPC1-MPC2 hetero-oligomerisation and blocks
pyruvate uptake -- a non-apoptotic, BH3-independent function (Kim et al., *Cancer Cell* 2019,
"Wild-type p53 promotes cancer metabolic switch by inducing PUMA-dependent suppression of oxidative
phosphorylation"). If PUMA sat upstream of respiration, PUMA falling should raise OXPHOS. It does
not.

Our data say MPC is not behaving as a PUMA-controlled node in this tissue: *Mpc2* is +0.69 at six
weeks (padj 0.0084) and *Mpc1* +0.56, both tracking the Myc programme, and both fall with age in
Myc+ (-0.47 and -0.10) -- with PUMA, not against it.

So the arrow runs the other way, and respiration is upstream of the trigger. Two anchors. Cells
without oxidative phosphorylation and with a low membrane potential are *less* susceptible to
apoptosis, and BCL-xL protection is unaffected by that (Dey & Moraes, *JBC* 2000) -- a de-respiring
mitochondrion is intrinsically harder to kill. And FOXO3 directly transactivates *BBC3*/PUMA, with
ChIP confirmation, under AKT control (Fitzwalter et al., *Dev Cell* 2018; Obexer et al.,
*Mol Cancer Ther* 2010).

The circuit that reconciles both readings is negative feedback: respiration -> FOXO3 -> PUMA ->
(a) MOMP and (b) an MPC brake on the very state that induced it. Under that circuit the co-decline
we observe is what is expected -- a de-respiring gland induces less PUMA, so there is less brake
and less death -- and PGC-1a, by raising respiration, should restore both. That is the experiment
the next section runs.

Three candidate mechanisms, each with an experiment that decides it:

| | mechanism | evidence now | what decides it |
|---|---|---|---|
| H1 | respiratory state gates the trigger | Myc x OXPHOS-priority on PUMA:BCL-xL, p = 0.0052; redox control p = 0.72 | the PGC-1a arm |
| H2 | FOXO3 activity | 0.456th and 0.513th percentile, adjacent; interaction p 0.0074 / 0.0081; wild type raises Foxo3 and Myc+ does not | the FOXO3 western -- FOXO3 is AKT-regulated post-translationally, so message alone cannot settle it |
| H3 | context loss: PUMA is a Myc/E2F1 target only in the TEB lane, and the TEB regresses | *Bbc3* appears in the AP_TEB TF lanes only; TEB-versus-ductal at the 0th percentile; the Fig. 1B TEB effect goes +1.02 SD at 6W to +0.23 at 12W | single-cell, or a TEB-staged perturbation; currently unresolved |

## Why does it not come out of the interaction, when it *is* an interaction?

Two reasons, and the author's instinct is right on both.

**Power.** A difference of two differences at n = 6 per group: median standard error 0.333 against
0.233 for a genotype contrast (a factor of 1.43, the root-two the algebra predicts), while the
median absolute effect it is asked to detect is *smaller* (0.208 against 0.257). The detectable
effect is about three times higher. Zero genes pass at 5% or 10% FDR.

**And the null is wrong.** Under a global rescaling every Myc-responsive gene has a non-zero
interaction by construction, so a test against zero measures the dose fade rather than the biology.
The specific question is departure from the x0.55 line -- and asked that way the signal is present:
*Bbc3* at z = -2.49 and *Foxo3* at z = -2.53, both in the bottom 0.5% of 8,774 genes, and in the
bottom 0.5% of the 4,199 Myc-induced genes taken alone, so neither is riding the repressive-arm
asymmetry. This is the same logic as the wrong-null correction already applied to the binomial
death test.

**And the non-linearity has a published mechanism.** Myc target promoters differ in affinity and
high-affinity promoters are saturated at physiological Myc (Lorenzin et al., *eLife* 2016), so
halving the protein cannot halve every target equally. Our retentions order as that predicts: Myc
core 0.74-0.77 > global 0.55 > mitochondrial arms 0.50-0.58 >> PUMA -0.09. Affinity explains
everything above the line and cannot produce a sign reversal, which is exactly why PUMA is the
residual worth an experiment.

---

# The panels

Numbering assumes paragraphs 1-2 keep Figure 1. Section 1 takes **Figure 2** (the oncogene),
section 2 takes **Figure 3** (the tissue and the trigger); the cell and MYC-ER work becomes
Figure 4.

**Most of this already exists.** `figures/fig01`-`fig05` and `figS1`-`figS8` were built for the
Quarto write-up and cover nearly every panel below. The work is porting them into the
`figures/panels/` idiom -- `theme_panel()`, a `panel_legend()` block, `save_panel_p()`, the
declared palette and contrast vocabulary, no prose on the page -- not building from scratch.

| slot | what it shows | port from | inputs |
|---|---|---|---|
| **2A** | the reallocation, ranked, with the twelve-week point joined by a connector so the fade reads in the same panel | `fig02_reallocation_ranked.R` | `background_vs_myc.rds$ruler` |
| **2B** | mitochondrial content rises -- the ruler that makes 2A mean what it says | `fig01_mito_content.R` (group boxplot idiom, lines 118-170) | `mito_content_proxies.rds` |
| **2C** | rescaled, not reshaped: Myc@12W against Myc@6W over 143 pathways, slope 0.552, R2 0.80, with the null percentile | `fig03_background_vs_myc.R` panel B | `background_vs_myc.rds$regressions`, `$regression_null` |
| **2D** | the dose: transgene message flat or rising against target output narrowing, beside the blot quantification | `figS8_myc_network_levels.R` | `interaction_results.rds`, western |
| **3A** | what the normal gland does: wild-type 6->12W per arm against expression-matched nulls -- the withdrawal, the assembly-factor control, the catabolic rise and the proliferation negative, in one panel | `fig04_substrate_specificity.R` panel A | `substrate_specificity_tradeoff.rds$wt_null` |
| **3B** | Myc builds the organelle and the means of its own execution; then everything fades at x0.55 except the trigger | `fig05_death_arm.R` panels A+B | `priming_arm_teb.rds$priming`, `collapse_module_ownership.rds` |
| **3C** | **new** -- the departure-from-dose distribution over 8,774 genes with *Foxo3* (0.46th) and *Bbc3* (0.51st) marked, *Bax* (59th) and *Bcl2l11* (92nd) as controls | none | `collapse_module_ownership.rds$collapse_genes` |
| **3D** | the coupling that names the experiment: Myc x OXPHOS-priority on PUMA:BCL-xL, p = 0.0052, redox control p = 0.72 | `fig04`'s trade-off panel | `substrate_specificity_tradeoff.rds$tradeoff` |

**Supplementary, all ports:** the two-ruler scatter (`figS5`) -- 94% content-up against a
two-sided priority spread; the full 144-pathway reallocation with adjusted p (`figS3` / `figS4`);
the four-gene contrast table (*Foxo3*, *Bbc3*, *Bax*, *Bcl2l1* across five contrasts); and the
luminal-release panel if that beat is taken.

**Not a panel, text only:** the DE counts and the interaction standard-error point. Both are
threshold arithmetic, and drawn they would flatter a 93% that is really a 45%.

---

# Traps that carry forward

1. Never lead on the DE count: -93% against a ~45% fall in effect size.
2. The genome-wide interaction has **zero** genes at 5% or 10% FDR. Attenuation is a set-level
   result. *Bbc3* is licensed by pre-specification, not by its adjusted p.
3. "Ranking showed no change" understates it -- NES *rises*. Write *preserved or sharpened*.
4. "The only arm the gland withdraws from" is not exact: import channels fall too (TOM -0.169,
   TIM23 -0.202, MIA40 -0.273). Write *by far the steepest, with import a distant second*.
5. The Myc+ gland's own OXPHOS decline is mostly dose and is not evidence for the substrate model.
   The two legs that are evidence are the **wild-type** withdrawal and the **PUMA** reversal.
6. `batch = timepoint`: describe, do not claim, every wild-type temporal statement. Once.
7. Every causal verb sits in a cell-experiment sentence. The mouse supplies the contrast; the
   perturbations supply causality.
8. Priming is not death. The transcriptome measures the molecular substrate; the phenotype is IHC
   and the MYC-ER counts.
9. The `pro_comp` and `priming` *composites* are MitoCarta-defined, so any coupling of a
   mitochondrial axis to them is mitochondria against mitochondria. The gene-level PUMA:BCL-xL
   ratio used here is not exposed to that, and neither is the departure-from-dose scan.
