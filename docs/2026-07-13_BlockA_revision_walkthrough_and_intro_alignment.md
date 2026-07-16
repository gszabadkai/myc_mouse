# Block A revision, point by point -- what we did, what we found, and how it aligns with the introduction

Written 2026-07-13, branch `paper-figures` (@ 5c00c5e). Audience: the author, building the
manuscript narrative. This is the TEACHING layer for Block A: the step-by-step revision (Issues
#1-6, scripts 26-31; Section 1) AND the cell-death / death-timing spine (scripts 12/14/16/23-25;
Section 1B) that the Introduction's closing sentence rests on. It restates each piece in plain terms
with the effect sizes, defines every gene set precisely (source + what genes it covers; full glossary
in Section 6), and then confronts the whole body of findings against the broad-overview paragraph
drafted for the Introduction, ending with a proposed reconciled rewrite.

- The dense chronological record is `docs/2026-07-08_BlockA_revision_plan.md` (the running log; all
  numbers here trace to its verified outcome sections and to scripts 26-31).
- The 5-bullet compression + figure inventory is `docs/2026-07-12_BlockA_revision_synthesis_and_figure_plan.md`.
- This doc does NOT supersede those; it is the readable walkthrough + intro-alignment they lack.

---

## 0. How to read this

**The experiment.** Bulk RNA-seq, MMTV-Myc transgenic mouse mammary gland, 4 groups x 6 replicates:
`6W_neg`, `6W_pos`, `12W_neg`, `12W_pos` (age x genotype). DESeq2 `~ timepoint * myc_status`.
"6W" = pubertal, TEB-rich; "12W" = adult. "neg" = wild-type (WT); "pos" = Myc transgene.

**The premise that frames everything.** Myc transgene dose is STABLE 6W->12W (mRNA + blot +
literature): a **constant driver acting on a moving substrate**. So every age-dependent change is a
change in the SUBSTRATE or in how the constant driver reads it -- not a change in Myc dose.

**What the revision is.** Issues #1-6 REFRAME already-fitted data (DESeq2 contrasts, GSVA
per-sample scores, fGSEA, mitoPPS). No script re-runs DESeq/GSVA/fGSEA; they re-read saved results.

**The honest ceiling, stated once, applies throughout.** Bulk tissue, n=6/group. Powered claims =
genotype MAIN effects (24 samples). Directional/indicative claims = within-timepoint or within-WT
contrasts (n=6-12), or composites of correlated gene sets. Co-variation is NOT causation, and bulk
cannot separate a per-cell change from a change in cell-type composition. Each issue names where it
sits on this ladder.

**Two rulers, kept distinct throughout (this trips people up):**
- **fGSEA NES** is RANK-based / magnitude-blind -- "are these genes at the TOP of this contrast?"
  It reads IDENTITY / priority, not size.
- **|LFC| magnitude** (and per-sample Cohen's d) reads SIZE / abundance change.
- **mitoPPS** is a within-mitochondrion RATIO (relative reallocation), on linear-scale normalised
  counts. **GSVA** is a cohort-relative per-sample score on log-scale input. mitoPPS "down" and an
  absolute-level "down" are DIFFERENT statements (this is the crux of Issue #4).

---

## 1. The six issues, point by point

### Issue #1 (script 26) -- How does Myc integrate with the normal developmental program?

**The question / what was wrong.** The earlier answer (`19_dev_composition.R`) only asked "does Myc
push toward progenitor?", collapsed ~179 developmental signatures into 5 lineage averages, reduced
differentiation to one HR-minus-LP scalar, and averaged the Myc effect over time -- discarding the
timepoint resolution that the whole question needs.

**What we did.** Reframe: resolve all **179 developmental sets individually**, labelled by source +
mammary-epithelial state + modality, across the FOUR-contrast trajectory, with a **dual lens**
(GSVA per-sample state + fGSEA importance). Contrasts: WT 6->12W (`timepoint_neg`), Myc+ 6->12W
(`timepoint_pos`), Myc@6W (`myc_6W`), Myc@12W (`myc_12W`), and the interaction. State uses the
consensus MEC nomenclature (Gray/Kessenbrock/Khaled, Dev Cell 2025): three discontinuous resting
types -- **BMYO** (basal-myoepithelial), **LASP** (luminal adaptive secretory precursor), **LHS**
(luminal hormone-sensing) -- from the author's hand-curated `data/dev_mec_annotation.csv`
(26 BMYO / 29 LASP / 29 LHS / 95 other).

**Gene sets used.** The 179 `MG_*` developmental sets (source `GRAY`, category
`Mammary_development`; e.g. `MG_TEB_VS_DUCTAL_*_GRAY_UP/DN`, 150 genes each = terminal-end-bud vs
mature-duct contrasts per lineage), plus the GRAY directional families (UP/DN, high- vs
low-EXPRESSION -- see the HEVSLE note in Section 6) and CHUNG ATAC OPEN/CLOSED chromatin nets. All
are library sets in `mammary_mito_myc_metab_v1_mouse.gmt`.

**Result.**
- The **WT substrate matures by LOSING luminal programs while BMYO stays FLAT** -- a rebalancing by
  luminal DECLINE, not basal expansion (BMYO wt_d = +0.07, p = 0.90; LASP d = -0.60; LHS d = -0.49;
  ~79% of luminal sets down; per-set p ~1e-5 but sample-level composite ns at n=6/tp -- matches the
  literature: no reported pubertal->adult BMYO expansion).
- **Myc bends OFF that axis, increasingly with age:** convergence rho(Myc, WT) = -0.08 at 6W
  (orthogonal) -> -0.44 at 12W (oppositional).
- Per state: **BMYO suppression is POWERED** (Myc effect -0.15/-0.21, d ~0.8-1.2, genotype main
  p = 0.020); **LHS flip is DIRECTIONAL** (-0.15 at 6W -> +0.17 at 12W, interaction p = 0.089 --
  carried by cross-modality agreement, not the p-value). Corroborated by the GRAY
  directional nets (Myc -> less-differentiated / TEB-proliferative end) and CHUNG ATAC (Myc closes
  basal chromatin, opens luminal by 12W).
- **LASP is heterogeneous (corrects the earlier "LASP near-null").** As a whole state the LASP
  interaction is net POSITIVE (GSVA +0.14, Myc@12W +0.09). But its 4-set **luminal-progenitor /
  secretory-precursor arm** (`MG_LAPRO_SAEKI`, `MG_LP_OPEN_CHUNG_ATAC`,
  `MG_LUMINAL_ALVSEC/ALVPROG_GARCIASOLA`) carries a NEGATIVE interaction -- all four sets, mean
  -0.07 GSVA and negative fGSEA importance at 12W -- i.e. Myc INDUCES the progenitor program at 6W
  (Myc@6W +0.31) and that induction FADES by 12W (+0.24). The alveolar-differentiation clusters
  (`MG_ALV_C*`) and the ATAC CLOSED pole net positive and dominate the state average, hiding the arm.
  This progenitor arm is the thread the selection bound (below) picks up.

**Interpretation + ceiling.** At 6W Myc imposes a broad off-axis identity suppression; by 12W a
specific **anti-BMYO / pro-LHS** reprogramming. Main effects powered (24 samples); the interaction
is directional at n=6/tp; composites use correlated sets. Association, not causation.

**Selection / death-dropout bounding (script 26 PART E -> `dev_program_myc_integration.rds$selection_bound`;
figures `e_selection_bound.pdf`, `e_lasp_arm_split.pdf`).** The BMYO suppression is asserted as a
PER-CELL state change, but the competing read is **compositional dropout**: high-MYC BMYO cells DIE
(oncogene-induced apoptosis), so the survivor bulk carries fewer BMYO transcripts -- an apparent
per-cell repression that is really selection. Bulk cannot separate fraction x per-cell, but two
already-computed lenses BOUND it (single-cell / FACS-sorted-basal SETTLES). Reading rule: a lineage
that death culls must sit ON the death-permissive axis (large |coupling|) and carry a death-gene-rich
program; |coupling| ~ 0 means the lineage is OFF that axis and dropout cannot manufacture its
repression. Results:
- **BMYO reads as PER-CELL, not dropout.** Its per-sample composite is essentially UNCOUPLED from the
  pro-death priming score (`pro_comp` from the death spine, script 25): Spearman rho ~ 0.15 (6W) /
  0.13 (12W) MEC-native, corroborating the even lower script-19 `BASAL_comp` coupling (-0.03 / -0.04)
  -- BMYO is OFF the death-permissive axis. And the BMYO program is NOT pro-death-gene-enriched
  (marker overlap vs the 512 pro-death genes of `cell_death_genes_consolidated`: OR 0.92, p = 0.81;
  panel pro-death fraction 0.034 vs base rate 0.036). Two independent facts already cut against a
  BMYO-death route: the MMTV-LTR drives the transgene in LUMINAL epithelium (so the high-MYC cells are
  luminal, not basal), and Myc CLOSES basal chromatin within cells (CHUNG ATAC) -- a per-cell
  remodeling signature, and the WT basal composite RISES with age rather than falling. Convergent:
  the BMYO suppression is a per-cell identity change.
- **The LASP luminal-progenitor arm is the residual candidate (dropout NOT excludable).** That 4-set
  arm is the MOST death-coupled of any lineage (rho ~ 0.87 at 6W, 0.62 at 12W) AND is pro-death-gene-
  enriched (OR 1.60, p = 0.006) -- and it is exactly the arm whose Myc induction FADES over the window
  (the negative interaction above). So for the luminal-progenitor compartment (unlike BMYO) death-
  driven dropout of Myc-induced cells cannot be excluded, and this is precisely the "shrinking
  luminal-progenitor" compartment the deconvolution plan already targets.
- **Ceiling.** These are association-level bounds at n=6/tp on survivor-biased bulk; gene-set overlap
  is not cell death and within-timepoint coupling pools genotype. The definitive separation of
  fraction x per-cell needs single-cell / FACS-sorted compartments / deconvolution (Section 5).

**Hypothesis-proofing: three routes to an apparent lineage shift, adjudicated per lineage.** Any
Myc-driven change in a lineage's bulk score over 6W->12W can arise three ways, and the machinery now
separates them. Using the exact interaction identity `interaction = Myc+ trajectory - WT trajectory`
(state_stats: `int = mycpos_shift - wt_shift`) plus the PART E death lens:
- **CONVERGENCE** -- the WT background matures and the genotype gap opens/closes with Myc doing nothing
  per-cell (the `wt_shift` term; the same WT-convergence axis as Issue #6).
- **FADING (per-cell)** -- the Myc+ tissue's own per-cell program moves (the `mycpos_shift` term):
  genuine transcriptional reprogramming of surviving cells.
- **SELECTION (compositional dropout)** -- part of the Myc+ trajectory is not per-cell but loss of a
  death-primed subpopulation; flagged only where the lineage is death-coupled AND pro-death-enriched.

Each of the three MEC lineages proves out to a DIFFERENT dominant route -- that non-uniformity is the
result (see `e_selection_bound.pdf`, `e_lasp_arm_split.pdf`, `state_stats.pdf`):

| Lineage / arm | interaction = mycpos - wt | Dominant route | Proof |
|---|---|---|---|
| **BMYO** | -0.064 = -0.051 - (+0.013) | **FADING (per-cell repression)** | WT is FLAT (wt_d 0.07 -> no convergence available); OFF the death axis (rho 0.15/0.13) and not pro-death-enriched (OR 0.92 -> selection excluded); powered genotype effect (p 0.020) + Myc CLOSES basal chromatin per-cell (ATAC: WT +0.29 vs Myc -0.33 -> -0.52). A genuine per-cell de-basalisation. |
| **LASP (state)** | +0.135 = -0.007 - (-0.142) | **CONVERGENCE (WT luminal decline)** | ~100% the WT term: WT loses LASP with age (wt_d -0.68, 79% sets down, per-set t p 4e-5) while Myc+ stays flat (mycpos -0.007). The apparent "Myc gain" is a mirage of the WT falling away -- not Myc per-cell. |
| &nbsp;&nbsp;-> **LASP luminal-progenitor arm** | -0.071 (all 4 sets negative) | **SELECTION candidate** | The one compartment where dropout is NOT excludable: most death-coupled of any arm (rho 0.87/0.62), pro-death-enriched (OR 1.60, p 0.006); Myc INDUCES it at 6W (+0.31) and the induction FADES by 12W (+0.24) -- a fade that is dropout-compatible. = the deconvolution target. |
| **LHS** | +0.322 = +0.190 - (-0.132) | **FADING (per-cell GAIN) + convergence** | ~59% a genuine per-cell Myc induction (mycpos +0.19; Myc OPENS ML/luminal chromatin +0.31 at 12W) on top of ~41% WT decline (wt -0.13). Death coupling (0.24 -> 0.73) is secondary and 12W-only -- Myc is actively driving hormone-sensing identity, not culling. |

**The overall proof.** No single mechanism carries the lineage reprogramming; the three MEC identities
are each dominated by a different route -- BMYO by a per-cell fade (de-basalisation), the LASP state by
WT-decline convergence, LHS by a per-cell gain. SELECTION survives adjudication in exactly ONE place:
the LASP luminal-progenitor arm -- consistent with the death spine, where the death-permissive
compartment is stem/progenitor, not basal or mature-luminal. This extends the Issue #6 mito
decomposition (which resolved only convergence + fade, no separable selection term) by adding the third
axis and showing selection is real but LOCAL. The author's original worry -- that the BMYO repression
is dying high-MYC BMYO cells -- is answered: it is not (BMYO is off the death axis, and MMTV drives MYC
in luminal not basal cells); the dropout worry belongs to the luminal-progenitor arm instead.

---

### Issue #2 (script 27) -- Endogenous Myc, the pubertal TEB program, and amplification

**The question.** The literature says endogenous Myc drives the WT pubertal proliferative phenotype
(especially terminal end buds, TEBs). Quantify that endogenous role and test whether the transgene
AMPLIFIES the same program rather than creating a new one.

**What we did.** Reframe on already-scored GSVA. Six per-sample program composites (below).
Decomposition: ENDOGENOUS role = the WT 6->12W temporal change; TRANSGENE effect = the genotype
MAIN effect (powered); plus the interaction.

**Gene sets used (the six composites; see Section 6 for members).** `myc` = the 17-set
**MYC_signatures** composite; `felsher` = **MYC_felsher_integrative_signature** (61 genes, the
mito-stripped MYC phenotype core); `hallmark_v2` = **MYC_HALLMARK_MYC_TARGETS_V2** (54); `prolif` =
the 14-set **Proliferation** composite; `myc_in_teb` = MYC targets in the Gray TEB context
(`TFT_MYC_GRAY_*_TEB`); `teb_ductal` = **MG_TEB_VS_DUCTAL** UP minus DN (the pure TEB-vs-duct axis).

**Result.**
- **Endogenous Myc gates the pubertal program (WT).** Within WT, Myc activity tracks proliferation
  r = 0.93, Myc-in-TEB 0.88, TEB-ductal 0.87 (tightest at 6W, 0.77-0.92); all six programs 6W > 12W
  in WT (directional, ns at n=6/tp) -- consistent with the literature.
- **The transgene amplifies it (POWERED).** Myc-target activity d = 2.3-3.0 (p < 1e-5),
  proliferation d = 1.0 (p = 0.02), Myc-in-TEB d = 1.6 (p = 6e-4). 6W_pos is the most
  TEB/proliferative state; the pure TEB phenotype is boosted SPECIFICALLY at 6W (6W_pos +0.41 vs
  6W_neg +0.05; ~4x larger at 6W than 12W, so the pooled TEB genotype p = 0.14 is weak only because
  it is time-concentrated) -- fitting TEBs as a pubertal structure (ties to the Issue #1 substrate).
- **Saturation nuance.** Within Myc+, the Myc->phenotype coupling loosens (r 0.87 -> 0.55) and the
  slope flattens (1.65 -> 0.83): endogenous Myc is rate-limiting in WT; the transgene pushes past
  that range -- amplification of LEVELS to a saturating high, not a proportional continuation.

**Interpretation + ceiling.** Frame Myc as an oncogene AMPLIFYING a normal pubertal
proliferative/TEB program, not creating a de novo one. Transgene powered; endogenous directional
(n=6/tp); "endogenous ESTABLISHES the TEB" stays a literature claim (we show WT co-variation +
amplification, not causation).

---

### Issue #3 (script 28) -- Is the MYC effect REALLY "primarily mitochondrial"?

**The question / what was wrong.** The pre-revision synthesis led with "Myc preferentially amplifies
mitochondrial biogenesis," resting on a MitoCarta enrichment that ranked #1 by a permutation-null z.
But z INFLATES with set size (MitoCarta n=1140, z=19.3 vs OXPHOS n=155, z=10.9 at ~equal effect), and
the MYC signature is also high and spans many arms. Does mito have independent substance, or is
"MYC effect = mitochondrial" just restating "MYC effect = MYC targets"? The one word "primarily"
does two jobs; the script splits them.

**What we did.** Reframe on scored GSVA + saved composites, no recompute.
- **Q1 (PREFERENTIAL?)** Report the SIZE-FAIR metrics -- effect magnitude and housekeeping-corrected
  relative EXCESS over the matched null -- plus a powered per-sample genotype Cohen's d ranking
  across all arms including central metabolism.
- **Q2 (CENTRAL?)** Couple each mito/metabolic axis to phenotypic OUTCOMES already measured (the MB2
  tumorigenic fork, mito death priming, TEB/dedifferentiation, proliferation), then test whether the
  coupling SURVIVES partialling out the generic MYC axis (partial correlation vs the Felsher core AND
  vs the full 17-set MYC composite), cross-checked within WT. Central = coupled AND survives
  deconfounding.

**Gene sets used.** Three MYC-identity arms (`myc_sig` 17-set, `felsher` 61-gene core -- 8% mito,
0% OXPHOS, `hallmark_v2` 54). Three mito arms built by regex over MitoCarta pathway sets:
`mito_all` = **MITOCARTA_ALL** (1140); `mito_oxphos` = the OXPHOS/complex/subunit/assembly/cristae
core; `mito_biogenesis` = the ribosome / mt-tRNA / mtDNA / import / translation arm. `prolif` (14).
Twelve central-metabolic axes (glycolysis, TCA, oxphos_met, PPP, nucleotide, one-carbon, glutamine,
FAO, FA-synthesis, amino-acid, cholesterol, redox), each a composite of `GS_metabolic` +
KEGG/Reactome/WP/Hallmark sets. Outcomes: `mb2_fork` (AP7), `priming` = MITOCARTA_APOPTOSIS_PRO
(25) minus ANTI (9), `teb_dediff` = MG_TEB UP minus DN.

**Result.**
- **Q1: "primarily mitochondrial" does NOT hold as literally #1.** Size-fair matched-null excess at
  6W: MYC-targets 0.65 ~ OXPHOS-core 0.64 > MitoCarta-wholesale 0.51 > E2F > Metabolism >
  Proliferation (mito is #1 by the size-confounded z but #3 by excess); MYC-targets OVERTAKE at 12W
  (0.69 vs OXPHOS 0.44 vs mito 0.29). Powered genotype Cohen's d: hallmark_v2 3.0 > felsher 2.7 >
  mito_biogenesis 2.6 ~ mito_all 2.57 > amino-acid 2.35 ~ glutamine 2.34 ~ myc_sig 2.32 >
  mito_oxphos 2.11 ~ TCA 2.09 ~ one-carbon 2.05 ~ glycolysis 2.01 ~ nucleotide 1.99 ~ oxphos_met
  1.94 > ... > prolif 1.01 > redox 0.08. Mito is a BROAD top tier, MYC-targets ahead, central
  metabolism co-equal -- not uniquely biggest.
- **Q2: centrality is REAL, ROBUST, and specific to an OXPHOS + biosynthetic core.** Raw couplings
  are all high (everything rides MYC dose); what survives partialling out MYC dose is the test.
  **10 robust central couplings, all in the bioenergetic-biosynthetic core:** `mito_oxphos` ->
  proliferation (raw 0.86 / |full-MYC 0.59 / within-WT 0.89), -> TEB-dediff (0.67/0.44/0.86),
  -> priming (0.57/0.32/0.34); `nucleotide` -> proliferation (0.85/0.54/0.92), -> priming
  (0.63/0.49/0.53), -> TEB (0.63/0.32/0.91); `oxphos_met` and `TCA` -> proliferation (0.37, 0.32).
  The mito **biogenesis/translation** arm and MitoCarta-WHOLESALE are BYSTANDERS -- coupling
  collapses under both proxies (biogenesis -> TEB 0.52 -> -0.03; -> prolif 0.74 -> 0.11). The MB2
  fork collapses for ALL mito arms (MYC-dose-driven) but is tracked by a distinct metabolic
  signature: cholesterol/mevalonate (+0.40) and redox (-0.50, even within WT; redox is barely moved
  by Myc, d = 0.08 -- coupled but not a Myc target).

**Interpretation + ceiling.** Replace "primarily mitochondrial" with **dissociable arms**: mito
biogenesis/translation tracks MYC-signature DOSE (bystander); OXPHOS + nucleotide + TCA co-vary with
the phenotype along an axis ORTHOGONAL to MYC-signature activity and detectable within WT. Surviving
the partial means orthogonal-to-MYC-dose, NOT MYC-independent or causal (proxy under-capture, a
non-linear MYC route, or a common cause such as cell-state composition all survive too; direction
unidentified). Q1 powered; Q2 per-sample association at n=6/group. Definitive test is genetic.

---

### Issue #4 (script 29) -- What IS the 6W->12W attenuation? (the "NES paradox")

**The question / what was fuzzy.** Gate 1 established the genotype effect SHRINKS ~2-fold 6W->12W
and that this is biological (effect-size, not power). But two things were confused: fGSEA NES of the
Myc-vs-WT contrast is FLAT across age even though the divergence shrinks (the "NES paradox"); and a
working gloss claimed "WT mito rises +0.82 with age and converges onto Myc." That +0.82 was a MEDIAN
over 63 heterogeneous MitoCarta sets -- an artifact of averaging opposite-direction arms. It is
RETRACTED here.

**What we did.** Reframe on fitted DESeq2 contrasts, GSVA, mitoPPS, fGSEA (VST / log-CPM are
deterministic transforms of the fitted data -- no re-fit). Part A: pathway-resolved decomposition +
the two-ruler NES-paradox figure. Part B: ABSOLUTE nuclear-OXPHOS levels across the 4 groups on 3
normalisations (the previously undone check). Part C: mitoPPS-vs-absolute reconciliation. Part D:
first-pass within-WT regulator regression.

**Gene sets used.** `MITOCARTA_OXPHOS` (155, nuclear -- zero mt-* genes, verified) and
`MITOCARTA_OXPHOS_SUBUNITS` (89; nuclear subset = subunits minus the 13 mtDNA-encoded genes) for the
absolute-level check; the mito/metabolic MitoCarta pathway sets for the pathway-resolved LFCs;
`luminal_sets` (LASP/LHS dev sets) and `tf_sets` = ESRRA_MITO / GABPA_MITO / NRF1_MITO (the ER/PGC1a
mito-biogenesis TF-activity sets) for Part D.

**Result.**
- **A -- the correction is vindicated in the LFC data.** Signed temporal LFC: OXPHOS FALLS in BOTH
  genotypes (OXPHOS_SUBUNITS WT -0.25 / Myc+ -0.39; OXPHOS -0.13 / -0.28), while the biosynthetic arm
  RISES in WT (amino-acid +0.18, lipid +0.12) and is flat/down in Myc+. There is NO WT developmental
  OXPHOS RISE; the +0.82 was averaging opposite arms. Genotype NES is flat/high at both ages
  (2.2-3.5) -- confirming the rank ruler is magnitude-blind. The magnitude attenuation is carried
  MORE by the mito/metabolic arms than the MYC-target core (attn ratio 12W/6W ~0.53 for
  OXPHOS/TCA/nucleotide vs 0.74-0.77 for MYC-targets; DE collapse OXPHOS 63 -> 2, ribosome 54 -> 3 vs
  MYC-targets 40 -> 17) -- a matter of DEGREE, on |LFC|, not exclusivity.
- **B -- the crux answer.** Myc RAISES nuclear-OXPHOS mRNA in ABSOLUTE terms: genotype Cohen's
  d ~1.27 (VST/normcnt p ~0.004; logCPM d 0.84 p 0.048), ROBUST across all 3 normalisations (per-gene
  group-mean Spearman 0.997-0.99998). Nuclear OXPHOS declines with age in both genotypes but only
  DIRECTIONALLY (WT beta -0.45 p 0.21; Myc+ -0.73 p 0.11). So the mitoPPS "OXPHOS down" is NOT an
  absolute-mRNA drop -- it is within-compartment REPRIORITISATION.
- **C -- reconciliation.** Both lenses agree in shape: nuclear OXPHOS peaks at 6W_pos (mitoPPS z
  +1.21 / absolute z +1.31) and drops by 12W; mtDNA-encoded OXPHOS is lowest at 6W_pos and rises to
  12W (mitoPPS -1.22 -> +0.37 / absolute -1.14 -> +1.23). So "nuclear-led assembly at 6W ->
  mtDNA-completed metabolism at 12W" holds in ABSOLUTE levels. Honest divergence: at 12W_pos nuclear
  OXPHOS is mitoPPS-deprioritised (0.955, <1) while its absolute level sits near baseline (+0.05) =
  reprioritisation WITHOUT a level drop.
  > **IMPRECISE + COHORT-ALIGNED (script 32 QC gate; wording corrected 2026-07-16).** The mtDNA
  > absolute z rise -1.14 -> +1.23 is a TEMPORAL claim and carries two caveats -- neither of which is
  > a depth artifact (**depth is not a confound**: shares and DESeq2 size factors both cancel it by
  > construction, and within cohort depth predicts nothing, 6W rho=-0.06 / 12W rho=-0.32, both ns).
  > (i) The mt-* share is **high-variance** (3.4%-40% across samples), so the rise is **imprecise** at
  > n=6/group. (ii) The 6W/12W depth ranges are **disjoint** (18.8-29.0M vs 9.9-16.1M), so batch is
  > perfectly aligned with timepoint and not separable -- inherent to any cross-sectional design, and
  > no artifact is positively indicated. The *genotype* half of Part B/C is unaffected (depth is
  > genotype-balanced p=0.41, litter-controlled). Report the rise; do not lean hard on its magnitude.
  > See the script-32 section above.
- **D -- first-pass regulators (hypothesis-generating, within WT, n=12).** The WT OXPHOS trajectory
  tracks BOTH the developmental luminal axis (partial r 0.71) and ER/PGC1a TF activity
  (ESRRA/NRF1/GABPA; partial r 0.70-0.83), model R^2 0.86-0.89.

**Interpretation + ceiling.** The attenuation is a shrinking of MAGNITUDE, not of program IDENTITY.
In absolute mRNA Myc RAISES nuclear OXPHOS; do not say Myc "reduces OXPHOS" -- say it REALLOCATES the
mitochondrial compartment. Part B is powered (genotype main effect); the age decline and Part D are
directional / within-WT associations. The retraction did not affect the committed Issue #3 result
(partial-correlation centrality does not depend on it).

---

### Does mitochondrial CONTENT per cell change? (script 32 -- author question, 2026-07-16)

**The question.** The manuscript says Myc drives a "robust increase in mitochondrial biogenesis". A
reviewer will ask whether that is mitochondrial CONTENT/mass per cell -- actual organelle biogenesis --
or only a transcriptional PROGRAM. Block A never answered it, and **three objects that sound like they
did, do not**:

- **`reframe_supp3$abund` (script 22, "AP-abund") is not an abundance measure.** It is the SD/IQR of
  mitoPPS Myc-effect diffs across 19 OXPHOS pathways (0.140 -> 0.100) = dispersion of *reallocation*,
  computed on ratio-normalised scores. The finalisation plan glosses it as "uniform biogenesis (an
  abundance-**like** change)"; the "-like" is load-bearing.
- **mitoPPS is abundance-BLIND by construction** (`08:404`: "cancels out both total mitochondrial
  content..."). So `mtnuc_index` -- the mitonuclear-imbalance headline -- is a *priority ratio*: a
  mitochondrion with 2x of everything scores identically.
- **`total_mito_score` (`08:1153`, captioned "reflects mito content") was plotted to a PDF and never
  saved or tested.** It is also `rowSums` over 142 *overlapping* MitoCarta pathways (double-counts
  genes) and is mtDNA-dominated.

And `bio_comp` (`19:91`) -- the composite behind "the Myc footprint IS a biogenesis program, rho=0.90"
-- is a mean of **GSVA** scores: rank-based, measures program, not organelle. The only real abundance
result was Issue #4 Part B, and only for OXPHOS subunits.

**What we did.** Script 32 computes per-sample **compartment shares** (share of transcriptome held by a
mito gene panel) on raw counts, two denominators (whole transcriptome; and one excluding the 13 mt-*
genes, which are ~77% of all MitoCarta counts and swing wildly). Genotype x time model mirrors Issue
#4's `level_stat_one` column-for-column, on `log2(share)`. Pure reframe -- no re-fit, 26-31 untouched.

**Gene sets used.** Library GMT panels (`MITOCARTA_ALL/_NUCLEAR_ENCODED/_MTDNA_ENCODED`,
`_MITOCHONDRIAL_RIBOSOME`, `_PROTEIN_IMPORT_AND_SORTING`, `_MTDNA_REPLICATION`, `_MTDNA_NUCLEOID`,
`_TRANSCRIPTION`, `_FISSION`, `_FUSION`, `_OXPHOS_NU/_MT`) -- **the mtDNA-machinery sets have never
been reported by any script**. Plus `MASS_MARKERS`, a 14-gene hand roster (Vdac1/2/3, Tomm20/22/40,
Tspo, Cs, Immt, Slc25a3, Timm23/44, Hspa9, Hspd1) = the in-silico stand-in for the TOMM20/VDAC/CS/HSP60
blot we do not have. This is a **flagged exception** to the "do not rebuild gene sets" rule: the
library has no mass-marker equivalent.

**Result (genotype axis; all 17 panels leave-one-out stable in sign and significance).**
- **Myc raises the mitochondrial share of the transcriptome.** Mass markers (chaperone-free) **+27%**
  (d=1.34, p=0.003); nuclear MitoCarta **+19%** (d=1.16, p=0.008); mitoribosome **+35%** (p=0.0014);
  protein import **+33%** (p=0.0016); nuclear OXPHOS **+28%** (d=1.15, p=0.009). The mass panel is
  **coherent**: 13/14 member genes up with Myc, 11 at p<0.05. `Tspo` is the lone non-mover
  (-0.10, p=0.30) and is reported individually rather than absorbed.
- **The chaperones are the largest effects and are handled separately.** `Hspd1` (+1.03 log2,
  p=3.3e-5) and `Hspa9` (+0.73, p=7.5e-5) are standard mass markers **and direct MYC transactivation
  targets**, so the chaperone arm has a route to elevation the structural arm does not. The
  **claim-bearing composite excludes them** (+27%); the full panel (+41%) is corroboration only.
- **Myc does NOT scale mtDNA-encoded output with it** (+6%, p=0.83) -- the script-24 mitonuclear
  imbalance, now in **absolute share** rather than a mitoPPS ratio. A materially harder claim to attack.
- **WHERE THE IMBALANCE IS *NOT* -- an informative negative (PART 3b).** The nuclear mito compartment
  rises by a **median +21%** (934 expressed nuclear-MitoCarta genes), and the mtDNA-dedicated core
  (Tfam/Polg/Polg2/Twnk/Ssbp1/Tfb2m/Polrmt/Tefm) rises **with** it, unremarkably: median +13%, **39th
  percentile** of that distribution, **Wilcoxon p=0.12**. So the imbalance is **not a transcript-level
  deficit anywhere in the machinery** -- nothing in the transcriptome is the bottleneck. That
  relocates the cause **downstream of transcript abundance**: mtDNA copy number, transcription rate,
  or turnover. It is what makes **mtDNA qPCR the discriminating experiment** rather than a
  nice-to-have.
  > **RETRACTED (2026-07-16, same day, before any figure depended on it).** The first version of this
  > section read the set-level `MITOCARTA_TRANSCRIPTION` (+43%, p=8.5e-5) and `MITOCARTA_MTDNA_NUCLEOID`
  > (+22%, d=2.32) results as **"commissioned but unbuilt"** -- Myc builds the mtDNA machinery and the
  > output does not follow. That was a **set-composition artifact**. Gene by gene, both effects are
  > carried by loosely-assigned, strongly-Myc-induced members -- `Mrpl12` **+116%** (a *mitoribosomal*
  > protein), `Atad3a` **+113%**, `Top1mt` +76%, `Poldip2` +58% -- while the actual core barely moves:
  > `Tfam` +11% (ns), `Polg` +3% (ns), `Twnk` +1% (ns), `Polg2` -9% (ns); only `Polrmt` +18% is
  > nominal (p=0.036). The **opposite** story ("Myc *spares* the mtDNA machinery") is equally
  > unsupported by the Wilcoxon above and is **not** substituted for it. **General lesson for this
  > corpus: MitoCarta process sets are membership-loose -- any set-level claim about a MECHANISM must
  > be resolved gene by gene before it is believed.** Script 32 PART 3b now enforces this.
- **Dynamics:** fusion +16% (p=0.022), fission +5% (ns) -- a mild fusion bias, reported because a
  reviewer will ask, not because it carries weight.
- **No interaction anywhere** (all int_p 0.40-0.95): Myc's content effect is constant across the window.

**Reconciliation (PART 5) -- three independent checks.**
1. **vs Issue #4** (different method entirely -- per-gene VST z-composite vs raw-count share):
   nuclear OXPHOS d=1.27 vs **1.15**; mtDNA d=0.03 vs **0.09**. Two lenses, same answer on genotype.
2. **vs script 24's mitoPPS**: the mtDNA component's group-shape agrees perfectly (Spearman 1.0). The
   nuclear component diverges (0.4) -- which is *exactly* Issue #4 Part C's "reprioritisation without a
   level drop" at 12W_pos, reappearing independently.
3. **The chaperone reconciliation.** Script 24's "Chaperones fall -> UPR^mt REFUTED" is a *mitoPPS*
   result (within-budget prioritisation). In absolute share the chaperones are Myc-elevated (+0.88,
   p=4e-5) and age-flat (WT temporal +0.05, p=0.87). **Not a contradiction** -- they deprioritise
   *inside* a reallocating compartment while their absolute share holds. Script 24's result stands.

Script 32 also **rebuilds `total_mito_score` correctly** (`MITOCARTA_ALL` counted once): it comes out
null (p=0.67), reproducing script 08's null -- and explains it. The mt-* genes are **76.9% of all
MitoCarta counts**, so the "total" was a noisy mt-* readout wearing a MitoCarta label. Counted mt-free,
the same compartment gives d=1.16, p=0.008.

**The QC gate -- what limits the TIME axis, and what does not (PART 4).**

> **CORRECTED 2026-07-16, on the author's challenge ("aren't the DESeq2 normalised values accounted
> for depth?"). They are -- and so are these shares. The first version of this section was WRONG.**
> It claimed the mt share "tracks library depth (rho -0.47)" and concluded TIME was **uninterpretable**.
> Both retracted:
> - **Depth cannot be a confound here.** A share is a **proportion** -- sequencing deeper multiplies
>   numerator and denominator alike -- so it is depth-invariant *by construction*, for exactly the same
>   reason DESeq2's median-of-ratios size factors cancel depth. There is no mechanism.
> - **Empirically depth does not act.** Within cohort it predicts nothing: 6W rho=-0.06 (p=0.87), 12W
>   rho=-0.32 (p=0.32). The pooled rho=-0.47 is **entirely** the gap between two clouds that differ in
>   both depth and mt% -- a pooling artifact, cited as if it were a signal.
> - **12W is not the degraded cohort.** Top-100 share of non-mt counts is **31% at 12W vs 40% at 6W**
>   -- 12W libraries are *more* even. The lower detected-gene count at 12W (23,255 vs 24,599) is what
>   lower depth gives you. The one metric bearing on RNA quality argues *against* the worry that was
>   raised.

**What actually limits the temporal claims -- two things, both duller than a confound:**
1. **Cohort-alignment (inherent, not mt-specific).** The depth ranges are **disjoint** (6W 18.8-29.0M
   from MYCF62-65/MYBS10x; 12W 9.9-16.1M from MYCF52-56), so the two ages were near-certainly
   different prep/sequencing batches, and batch is perfectly aligned with timepoint. No common support
   => not separable. **But this is true of any cross-sectional design** (different animals per age);
   the whole project already lives with it, and nothing positively indicates an artifact.
2. **Imprecision (mt-specific -- the real caveat).** The mt-* share spans **3.4%-40%** across samples,
   SD 13 pp within 12W_neg alone, one sample at 40.1%. Temporal mt claims are **imprecise** at n=6.

So Issue #4's "mtDNA absolute z -1.14 -> +1.23" is **imprecise and cohort-aligned, not
uninterpretable**. The genotype axis is unaffected either way: depth is genotype-**balanced** (p=0.41)
and both genotypes occur inside single litters (MYCF62, MYCF52, MYCF56). There is no RIN or batch
column on disk. The gate **diagnoses**; there is nothing mechanical to correct.

*(The scripts 29-vs-32 divergence on the WT temporal beta, -0.454 vs +0.094, was previously cited as
corroborating a confound. It does not: the two use different denominators and different composites, so
they are not expected to agree on a temporal slope. Dropped as evidence.)*

**Interpretation + the three ceilings.** *Myc increases mitochondrial content per cell -- modestly
(~25-27% by transcript proxy), coherently across structural, import, ribosomal and OXPHOS arms, and
without scaling mtDNA output.* This converts the vulnerable "robust increase in mitochondrial
biogenesis" from an enrichment-score claim into a quantified, bounded one. The ceilings:
1. **Bulk polyA with no spike-ins and no cell counts cannot measure per-cell content in absolute
   terms.** DESeq2's median-of-ratios normalisation removes exactly that scale. Every number here is a
   **share of transcriptome**.
2. **Myc globally amplifies total RNA per cell**, so a *constant* share would already imply more
   mitochondria per cell => **+27% is a LOWER BOUND** on the Myc-driven increase, and bulk cannot
   recover the true one.
3. **Transcript share is not protein and not organelle volume.** A TOMM20/VDAC/CS blot, mtDNA qPCR, or
   EM **settles** it; n=6/group. Note this is the *same* orthogonal experiment the protein-level OXPHOS
   question (tension A) already needs -- one blot answers both.

---

### Issue #5 (script 30) -- Can MYC-independent factors EXPLAIN the attenuation?

**The question / why Part D was the wrong instrument.** Issue #4 Part D showed OXPHOS co-varies
within WT with developmental and ER/PGC1a-TF axes (LEVEL association). The author then asked the
sharper question: can we "blame" those MYC-independent factors to FULLY EXPLAIN the reduction? The
attenuation is a MODERATION phenomenon -- the genotype effect DEPENDS on timepoint, i.e. the
genotype x time INTERACTION (b_int). "Explaining the reduction" = does conditioning on a candidate
axis SHRINK b_int toward 0? Level-association is not moderation, so Part D could not answer it.

**What we did.** Reframe reusing script 29's saved set-lists (so the baseline b_int reproduces the
Issue #4 result exactly). Part A: baseline b_int per outcome. Part B: covariate absorption per axis
-- refit `outcome ~ tp*myc + axis + tp:axis`, reporting delta_b_int (the STABLE primary metric),
absorption fraction (with the caveat it is unstable at a small baseline), the interaction-survives
p, delta_R2, and the nested LRT. Part C: stratified case bootstrap (R=2000) for CIs. Part D: within-
WT (n=12) OXPHOS slope before/after holding dev + TF. Part E: what would settle it.

**Gene sets / axes used.** Outcomes: `oxphos_abs` (VST nuclear-OXPHOS composite -- powered,
mtDNA-clean) and `oxphos_gsva` (secondary). Axes: `dev_luminal` (LASP/LHS luminal composite),
`tf_biogenesis` (ESRRA/GABPA/NRF1 _MITO), `prolif` (14-set), `mtdna_reprior` (mtDNA-encoded OXPHOS
mitoPPS), and `dev_plus_tf` (both).

**Result -- the answer is NO.** Baseline b_int on oxphos_abs = -0.279 (p 0.607, directional --
reproduces Issue #4 exactly); on oxphos_gsva ~0 (0.003), so its absorption fraction is undefined.
- **Part B.** Conditioning on dev + TF together leaves the attenuation essentially UNTOUCHED
  (-0.279 -> -0.322, delta -0.043). dev_luminal ALONE makes it LARGER (-0.279 -> -0.721, delta
  -0.442) -- a co-drive signature (removing a program Myc co-drives EXPOSES more Myc effect, it does
  not absorb it).
- **Part C (the honest ceiling).** EVERY delta_b_int bootstrap CI straddles 0 (dev+tf [-1.39, 1.31];
  dev_luminal [-1.17, 0.40]) -- at n=24 there is no reliable movement in either direction. The
  absorption-FRACTION CIs are uninterpretably wide because the baseline interaction is directional;
  "fully explains" is nowhere distinguishable.
- **Level-association is NOT moderation.** Every axis LRT is significant (5e-4 to 1e-10): the axes
  DO add variance for OXPHOS LEVEL -- they just do not absorb the genotype x time INTERACTION.
- **Part D.** Within WT the OXPHOS 6->12W slope only shrinks ~29% (abs -0.454 -> -0.323) / ~19%
  (gsva) when dev + TF are held -- it PERSISTS.

**Interpretation + ceiling.** No MYC-independent axis absorbs the attenuation; the dev/TF programs
CO-VARY with mitochondrial output rather than antagonising Myc. Candidate axes are ENDOGENOUS (Myc
drives dev/TF too), so absorption estimates are biased and this test BOUNDS rather than resolves a
MYC-independent contribution. The identifying evidence is experimental (inducible Myc, TF
perturbation, single-cell, ATAC/ChIP).

---

### Issue #6 (script 31) -- A CITABLE MECHANISM for the attenuation

**The question.** Issues #4/#5 fixed WHAT the attenuation is (magnitude, not rank; Myc raises
absolute OXPHOS; no dev/TF axis absorbs it) -- but "compression" only renames it. The author wanted
a citable MECHANISM. The restatement ("the MYC signal is not reduced; it acts on a different
background") splits into candidate mechanisms; a read-only check showed global response compression
is real but modest (SD ratio 0.83) and cannot be the citable mechanism.

**What we did.** Use an EXACT identity in the interaction model (verified, cor = 1.0000): the
per-gene attenuation IS the interaction,
`myc_12W - myc_6W == timepoint_pos - timepoint_neg`. Aligning each gene to Myc's induction direction
`d = sign(myc_6W)`, the gap-shrink splits into two mechanisms:
- **WT-convergence** = `d * timepoint_neg` -- wild-type tissue matures TOWARD Myc's state (the moving
  background);
- **Myc-fade** = `d * timepoint_pos` -- the oncogenic program retreats on the aging substrate;
- `attenuation = WT-convergence - Myc-fade`.

Part A: the decomposition (POWERED, contrast-level, all genes) over two universes (6W-divergent
padj < 0.1; and a selection-independent effect universe |LFC6W| > 0.5) x a program roster. Part B: a
broadened-TF absorption panel (extending Issue #5) + a compositional diagnostic (each TF axis's Myc+
6->12W slope vs the OXPHOS fade). Part C (reference deconvolution): DEFERRED.

**Gene sets used.** Program roster (MitoCarta pathway sets + two pooled non-mito representatives):
`MITOCARTA_OXPHOS_SUBUNITS` (89) and `MITOCARTA_OXPHOS` (155) = OXPHOS core; `MITOCARTA_TCA_CYCLE`
(20); `MITOCARTA_NUCLEOTIDE_METABOLISM` (36); `MITOCARTA_AMINO_ACID_METABOLISM` (90) +
`MITOCARTA_LIPID_METABOLISM` (126) = biosynthetic; `MITOCARTA_MITOCHONDRIAL_RIBOSOME` (83) =
biogenesis; `MYC_HALLMARK_MYC_TARGETS_V2` (54) + `MYC_felsher_integrative_signature` (61) =
MYC-target core; pooled `PROLIFERATION` (all `PROLIF_*`) and pooled `MAMMARY_LUMINAL` (LASP/LHS).
Broadened TF axes: `tf_biogenesis` (ESRRA/GABPA/NRF1 _MITO); `tf_e2f` (`TFT_E2F[0-9]_GRAY_*_MITO` +
`PROLIF_E2F_HALLMARK`, 200); `tf_esr1` (`TFT_ESR1_GRAY_*` + `ESR1_AND_CORE_MITO`, 59); `tf_lipogenic`
(`TFT_(SREBF1|SREBF2|MLX)_GRAY_*_MITO`); `tf_mito_panel` (data-driven top-20 mito-TFs from
`gray_chea_mito_tf_shortlist.csv`, excluding MYC/E2F/etc).

**Result.**
- **A -- the mechanism (powered, robust across both universes).** Globally the attenuation is
  **~66% Myc-fade, ~34% WT-convergence** (42.8% convergence on the effect universe). It is
  **PROGRAM-SPECIFIC:** WT-convergence is real for the biosynthetic arm (amino-acid 37%, lipid 41%,
  nucleotide 31%; TCA 20%) and pooled mammary-luminal (36%) -- the WT gland matures onto the same
  biosynthetic axis Myc drives. But convergence is ABSENT/NEGATIVE for OXPHOS (-11 to -39%; WT
  DIVERGES, the gap closes purely by fade) and for the MYC-target core (-60 to -63%), which is ALSO
  the LEAST attenuated (atten 0.17-0.19 vs 0.28-0.50) -- the MYC identity core is PROTECTED. Pooled
  proliferation attenuates by PURE fade (conv ~0%, fade 100%). Every sign reproduces on the
  selection-independent universe. This is the powered, gene-level resolution of the old H1
  (front-loaded fade) vs H3 (WT catch-up) question: H1-dominant with a real, biosynthetic-specific
  H3 minority.
- **The fade has two faces (the shared-developmental reframe).** For fade-led programs the Myc+
  temporal slope = a SHARED developmental decline (also seen in WT) + a MYC-SPECIFIC excess. A shared
  decline CANCELS in the genotype gap, so only the MYC-specific excess constitutes the attenuation
  (OXPHOS ~72-90% Myc-specific). This is a guardrail: do not over-read the fade as generic
  exhaustion -- the part that matters is Myc-specific.
- **B -- broadened-TF absorption (BOUNDED) + a compositional hint.** No broadened TF axis reliably
  absorbs the attenuation (every delta_b_int bootstrap CI straddles 0; ESR1 does not even associate,
  lrt_p 0.83). The standout is `tf_e2f`: it moves b_int on oxphos_abs from -0.279 to -0.028
  (point-estimate ~90% absorption) -- but NOT significant (CI [-0.64, 2.17]), a hint not a result.
  Compositional diagnostic: the OXPHOS fade is a Myc+ slope of -0.91; the axes co-declining most are
  `tf_mito_panel` (-0.60), lipogenic (-0.46), E2F (-0.36), while ESR1 RISES (+0.42). Proliferation /
  mito-TF activity fades ALONGSIDE OXPHOS.
- **C -- DEFERRED.** Reference deconvolution of the proliferative/TEB fraction (needs an external
  mouse-mammary sc atlas, e.g. Bach 2017 GSE106273, + MuSiC/Bisque; not on disk) is the next step;
  single-cell or TF perturbation is definitive. See `docs/deconvolution_subproject_plan.md`.

**Interpretation + ceiling.** A citable two-mechanism account: the genotype gap closes by ~2/3
Myc-fade (oncogenic retreat) + ~1/3 WT-convergence (the WT gland maturing onto Myc's biosynthetic
axis), with OXPHOS and the MYC core NOT converging and the MYC core selectively buffered. Part A is
POWERED and the identity is exact -- but it is DESCRIPTIVE of the transcriptional change: the fade
term still confounds per-cell weakening with compositional dilution of the shrinking
proliferative/TEB compartment (the E2F/prolif co-decline is a pointer to composition, not proof).
Resolving that requires single-cell / deconvolution.

---

## 1B. The cell-death / death-timing spine (Block A exploratory: scripts 12, 14, 16, 23-25)

The six revision issues (Section 1) are the mode-of-action half of Block A. The OTHER half --
produced earlier in Block A exploratory, and the direct support for the Introduction's closing
sentence -- is the cell-death / death-timing spine. It is the substrate the revision built on
(Issue #3 uses the death-`priming` axis as a phenotypic outcome), so it belongs in this walkthrough.

### The external phenotype that anchors it (independent of this RNA-seq)

In a PARALLEL Myc-ER tamoxifen-inducible model (no RNA-seq there), an acute Myc pulse at 6W causes
MUCH more cell death than the same pulse at 12W (IHC). Constant acute pulse, different tissue age =>
the death is **substrate-gated by developmental stage**, not by Myc dose. Our chronic MMTV-Myc
RNA-seq cannot see the inducible phenotype but CHARACTERISES the 6W-vs-12W substrate the pulse hits.
Because Myc is off pre-tamoxifen there, the cleanest readout of the death-permissive state is the WT
`timepoint_neg` contrast; the chronic-Myc+ layer is secondary and SURVIVOR-BIASED (we sequence the
cells that did NOT die -- which makes any Myc death-engagement we DO see a conservative floor).

### Two cell-death lenses + a decision gate (scripts 12, 14, 16 -- both branches, neither canonical)

- **Branch 1 -- directional binomial on RAW dLFC (script 16).** Tests the a priori hypothesis that
  pro-death genes are induced more by Myc at 6W and pro-survival more at 12W, on
  `cell_death_genes_consolidated` using the metric `myc_6W_log2FC - myc_12W_log2FC` (RAW/unshrunken,
  per the shrinkage rule). Result: **NULL** (52.7% in the predicted direction, p = 0.27) -- no
  genome-wide directional pro/anti asymmetry. A clean negative that keeps the story honest.
- **Branch 2 -- Tang 2024 15-RCD fGSEA (script 12).** fGSEA of 15 regulated-cell-death modalities
  (Tang et al. 2024) on the Wald z of five contrasts. Result: **DIRECTIONAL** -- the apoptosis
  modality has temporal_pos NES -1.38 (padj 0.016): the apoptotic program DECLINES in Myc+ with age.
- **Gate 2 -- decision-point vs readout (script 14).** Asks whether the pro-apoptotic module shifts
  COORDINATELY on the interaction (selection / active decision) or only as isolated genes (passive
  readout). Result: **CONFIRMED coordinated shift** (Apoptosis-PRO module t p = 0.023) -> the
  DECISION-POINT framing survives; cell death stays a main-figure candidate.

### The death-timing SUBSTRATE model (script 23) -- four hypotheses, WT-anchored

- **"WT substrate de-primes with age" is REFUTED.** WT baseline BH3:BCL2 priming is flat/slightly
  rising (6W_neg -0.24 -> 12W_neg -0.08). The gating is NOT in the WT baseline -- do not anchor the
  story there.
- **H1 (BH3/BCL2 rheostat) CONFIRMED as Myc COUPLING, and powered.** Myc raises pro-apoptotic priming
  (PRO-ANTI genotype effect +0.32, p = 0.038; myc_6W PRO module +0.20, p = 0.041), ~3x more at 6W
  than 12W (group-mean gap +0.48 vs +0.15). p53/ARF is NULL -- the rheostat is the BCL2 family, not
  p53. (Survivor bias makes this a conservative floor.)
  - **p19ARF / p53 axis, explicitly re-tested at gene level (2026-07-13; canonicalised as script 23
    PART 2b -> `death_timing_substrate.rds$h1$arf_p53`, figure `h1c_arf_p53_axis.pdf`).** The
    canonical Myc-antiapoptosis escape -- p19ARF induction driving p53,
    then its loss/relief permitting survival -- does NOT explain the 12W death drop. ARF is not lost:
    `Cdkn2a` is MORE Myc-induced at 12W (myc6 +0.14 -> myc12 +0.53; interaction independent-filtered,
    so no q, but direction is anti-hypothesis; mRNA cannot separate p19ARF from p16INK4a). p53 is not
    disengaged: `Mdm2` is flat (no added braking) and the p53 TARGET panel is not coordinately reduced
    at 12W -- `Cdkn1a`/p21 is less suppressed (+0.33) and `Zmat3`/`Eda2r`/`Phlda3`/`Trp53inp1`/`Sesn2`
    rise; p53 output looks intact-to-higher. The only hit is `Bbc3`/Puma (myc6 +0.26 -> myc12 -0.28,
    interaction z -2.6, padj 0.84) with `Bax` trailing (-0.23) -- both pro-apoptotic BCL2-family
    EFFECTORS, i.e. the BH3 rheostat downstream of p53, not a p53-activity failure. Nothing survives
    FDR (all int_q ~1). => Corroborates "p53/ARF null" and localises the fade to the BCL2-family
    effector arm coupled to the mitonuclear-imbalance spine, not to an upstream ARF/p53 switch.
  - **p53-INDEPENDENT PUMA regulators, gene + signature (2026-07-13; canonicalised as script 23
    PART 5b -> `death_timing_substrate.rds$puma_regulators`, figures `h1d_puma_regulators.pdf`,
    `h1e_puma_tf_signature_coupling.pdf`).** PUMA (Bbc3) has p53-independent drivers
    (`docs/library_reference/PUMA-and-its-relationships.md`): FOXO3a (gated by a PGC-1a brake), FOXO
    BH3 co-targets Bim/Noxa, and the executioner HTRA2. *Positive:* at GENE level a FOXO3 / PUMA /
    HTRA2 module is coordinately Myc-induced at 6W and withdrawn by 12W -- `Foxo3` mirrors `Bbc3`
    almost exactly (Foxo3 dMyc -0.49, interaction z -2.7; Bbc3 z -2.6) and `Htra2` is the most
    6W-induced gene in the panel (+1.49 -> dMyc -0.55). A coherent p53-independent candidate for the
    PUMA loss. *Negative / bounding:* (i) the textbook PGC-1a BRAKE is REJECTED -- ESRRA/PGC1a do NOT
    rise as a brake, they CO-DECLINE (Esrra dMyc -0.32); (ii) the FOXO output is PUMA-SELECTIVE, not
    the full program (Bim/Bcl2l11 dMyc -0.04, Noxa +0.19, antioxidant Sod2/Cat flat); (iii) crucially,
    the FOXO3 target-ACTIVITY signature (TFT_FOXO3_CHUNG, 38 genes) does NOT track Bbc3 at 6W (r -0.14,
    ns) nor decline (trajectory interaction +0.07, p 0.76) -- so the `Foxo3`-gene co-movement is NOT
    confirmed as a FOXO3-activity cascade. => Net: the p53-independent FOXO3/HTRA2 module co-moves with
    PUMA (strengthening "a BH3/effector arm is de-amplified"), but the balance of evidence (activity
    flat, PGC1a co-declining, program PUMA-selective) reads as DE-AMPLIFICATION of the 6W death-primed
    compartment rather than a specific TF-activity switch -- consistent with the compositional spine.
    Directional, n=6, nothing survives FDR.
- **H2 (mitonuclear imbalance) = the MECHANISM (the standout).** The nuclear-minus-mtDNA OXPHOS
  mitoPPS imbalance peaks at 6W_pos (+0.50, the death-permissive group); its coupling to PRO priming
  is r = 0.75 (p = 0.005) at 6W -> 0.07 (p = 0.82) at 12W. Two-lens nuance mirroring Issue #4: Myc
  raises PRO expression ABSOLUTELY (+0.20) but RELATIVELY de-prioritises Apoptosis-PRO within the mito
  compartment (mitoPPS diff -0.21, padj 0.022) -- biogenesis crowds it out.
- **H3 (proliferation-death) REFUTED as the cause.** prolif~death coupling is INVARIANT (0.68 ->
  0.65): the timing is a mito/biogenesis phenomenon, not generic proliferation.
- **H4 (selection/culling) DIRECTIONAL, weakest.** Cross-sample dispersion of priming narrows
  6W->12W, more in Myc+ (sd_pro ratio 0.42 vs 0.60) -- consistent with, but not proof of, culling.

### The developmental "why" (scripts 24-25) -- what matures in the substrate

- **The developmental SUBSTRATE (not chronic-Myc adaptation) sets the timing, on two COUPLED axes.**
  (1) MITO axis (ROBUST): the mitonuclear imbalance resolves 6W->12W, ~88% by mtDNA RISING (WT 89% /
  Myc 87%) -- genotype-shared developmental maturation (nuclear-led assembly at 6W -> mtDNA-completed
  running metabolism at 12W; the same mitonuclear handover as Issue #4). (2) CELL-STATE axis:
  progenitor/luminal down, basal up, but DIRECTIONAL and underpowered (all p > 0.18) -- do not claim a
  significant composition shift.
- **Which background axis tracks death (script 25 Part B).** At 6W, per-sample coupling to pro-death
  priming: biogenesis 0.83 > stem/MASC content 0.78 > mitonuclear imbalance 0.67 -- ALL decouple by
  12W. So the mito imbalance is NOT the sole death-coupled axis; stem/progenitor content is comparably
  coupled (vindicating "composition matters"). A coupled death-permissive STATE, not separable causes.
- **Where "ER/PGC1a biogenesis" lands (script 24).** Myc co-opts the ER/PGC1a machinery
  (ESRRA/NRF1/GABPA co-move with MYC; ESR1/ER antagonised, re-emerges at 12W); balanced CORE/MITO_NU
  biogenesis-death intersections are most death-coupled at 6W (r 0.77 / 0.83).

### The spine, in one line (and how it meets the OXPHOS thread)

At 6W the young substrate is death-PERMISSIVE: stem-rich, mitonuclear-IMBALANCED (nuclear OXPHOS
running ahead of mtDNA), biogenesis-primed, and all three tightly coupled to pro-apoptotic priming --
so Myc's induction lands as oncogene-induced apoptosis (killing). By 12W the substrate has MATURED
(mtDNA caught up, imbalance resolved, cells differentiated, priming/biogenesis converged), the death
coupling COLLAPSES (r 0.75 -> 0.07), and Myc's killing fades. **The 12W decoupling IS the loss of
killing** -- the desensitisation that opens the permissive window. Crucially this is the SAME
mitonuclear maturation as the OXPHOS thread: the 6W nuclear-led-assembly imbalance (Issue #4) is the
death-permissive state, and its resolution -- together with the Myc-specific OXPHOS de-amplification
over the same window (tension A) -- is what desensitises the adult gland. The OXPHOS reprioritisation
and the apoptosis desensitisation are two faces of one developmental mitochondrial maturation.

**Ceiling.** Bulk RNA, one timepoint per age, n=6/group, SURVIVOR BIAS. This characterises the
death-permissive STATE (association) and cannot prove causation or fully resolve mechanism; BH3
profiling / single-cell / caspase-by-state would be needed. Correlated axes (mito, stem, biogenesis)
are not separable at this n. The branch-1 null and the p53/ARF null are firm; the couplings and H4
culling are directional.

---

## 2. The through-line (Issues #1-6 + the death spine as one arc)

A constitutively expressed Myc AMPLIFIES the normal pubertal proliferative/biosynthetic program of
the mammary gland (strongest at 6W, the TEB-rich window) rather than creating a new one, bending the
gland OFF its maturing developmental axis (#1-2). Its transcriptional footprint is BROAD, not
narrowly mitochondrial: the OXPHOS/biosynthetic core is engaged CO-EQUALLY with the MYC-target and
central-metabolic programs, and the arm most tightly COUPLED to the phenotype -- along an axis
orthogonal to generic MYC dose and visible within WT -- is a specific respiratory/biosynthetic core
(OXPHOS + nucleotide + TCA), NOT mitochondrial biogenesis, which merely tracks MYC dose (#3). Within
the mitochondrial compartment Myc REALLOCATES resources: it RAISES absolute nuclear-OXPHOS transcript
while reprioritising the compartment (nuclear-led assembly at 6W -> mtDNA-completed metabolism at
12W), so the "OXPHOS reduction" is relative, not an absolute transcript drop (#4). The 6W->12W
attenuation is a loss of MAGNITUDE (not rank), which no MYC-independent developmental/TF axis absorbs
(they co-drive) (#4-5); decomposed exactly, it is ~2/3 oncogenic fade + ~1/3 wild-type convergence
onto the shared biosynthetic axis, sparing the buffered MYC core, and the fade co-occurs with a
collapse of proliferative/mito-TF activity -- a compositional hint (#6). This mitochondrial
maturation is also the death spine (scripts 12/14/16/23-25): at 6W the mitonuclear-imbalanced,
biogenesis-primed, stem-rich substrate is tightly coupled to pro-apoptotic priming, so Myc's
induction lands as oncogene-induced apoptosis (killing); by 12W the imbalance resolves (mtDNA catches
up), the death coupling collapses (r 0.75 -> 0.07), Myc's killing fades, and the desensitised adult
gland is the permissive window for tumourigenesis. The OXPHOS reprioritisation and the apoptosis
desensitisation are two faces of the SAME developmental mitochondrial maturation -- which is the
causal spine of the Introduction's closing sentence.

---

## 3. Confrontation and alignment with the Introduction paragraph

The drafted paragraph (verbatim):
> Although MYC is reported to stimulate overall mitochondrial biogenesis, we find that its effect on
> the mitochondrial transcriptome is not uniform, and the dynamics of mitochondrial pathway
> reprioritisation is central to how transformation proceeds. MYC plays a part and interacts with the
> transcriptional trajectory of normal mammary development, leading to a robust increase in
> mitochondrial biogenesis. However, key mitochondrial pathways are selective reprioritised, holding
> OXPHOS restrained. Reduced abundance of OXPHOS complexes is associated with a desensitisation to
> MYC-induced apoptosis, opening a permissive window for tumourigenesis at the adult developmental stage.

Clause-by-clause verdict (SUPPORTED / NUANCE / REWORD):

| # | Claim in the paragraph | Issue(s) | Verdict |
|---|---|---|---|
| 1 | "MYC is reported to stimulate overall mitochondrial biogenesis" (prior/literature framing) | #3, #4 | SUPPORTED -- accurate as background; our data also confirm biogenesis rises. |
| 2 | "its effect on the mitochondrial transcriptome is not uniform" | #3, #4 | SUPPORTED -- dissociable arms (#3); OXPHOS falls while biosynthetic rises in WT (#4). |
| 3 | "the dynamics of mitochondrial pathway reprioritisation is central to how transformation proceeds" | #3 (Q2), #4 | NUANCE -- centrality is COUPLING-level (partial-corr, orthogonal to MYC dose), not demonstrated causality/necessity. Keep, but flag as association. |
| 4 | "MYC ... interacts with the transcriptional trajectory of normal mammary development" | #1, #2 | SUPPORTED -- Myc bends off the maturing WT axis; amplifies the endogenous pubertal program. "plays a part and interacts" is well-calibrated. STRENGTHENED (Issue #1 PART E): the interaction is adjudicated per lineage into three routes -- BMYO by a genuine per-cell fade (de-basalisation), the LASP state by WT-decline convergence, LHS by a per-cell gain -- and is a real cell-intrinsic reprogramming, NOT a selection artifact (dropout survives only in the death-coupled luminal-progenitor arm). |
| 5 | "leading to a robust increase in mitochondrial biogenesis" | #3, #4 | NUANCE -- true and POWERED (biogenesis d ~2.6; absolute nuclear OXPHOS d ~1.27), but biogenesis is the MYC-DOSE-tracking BYSTANDER; the phenotype-coupled arm is OXPHOS/nucleotide/TCA. Precision, not contradiction. |
| 6 | "key mitochondrial pathways are selective[ly] reprioritised" | #4 (mitoPPS), #3 | SUPPORTED -- the reprioritisation (nuclear-led -> mtDNA-completed) is a central finding. |
| 7 | "holding OXPHOS restrained" | #4, #6 | NUANCE (author's intent = the TIME axis) -- over 6->12W the Myc+ gland pulls its elevated OXPHOS back down (a de-amplification). Valid on the time axis; NOT a genotype-level suppression (vs WT, Myc RAISES OXPHOS, d ~1.27). Keep, but scope to "over the adult transition." |
| 8 | "Reduced abundance of OXPHOS complexes ..." (author confirmed: the Myc+ 6->12W temporal decline) | #4, #6, death-timing | NUANCE -- the temporal decline is real (subunits -0.39, core -0.34) and its shared/specific split is POWERED, but it is 72-90% MYC-SPECIFIC, NOT a reflection of WT development; and it is a NUCLEAR-subunit decline while mtDNA-encoded OXPHOS RISES (a mitonuclear rebalancing). Window: SUPPORTED. |

### The three tensions to resolve

**(A) "holding OXPHOS restrained" / "reduced abundance of OXPHOS complexes" -- resolved by axis
(author clarified: the TIME axis, Myc+ 6W->12W).** The apparent contradiction with Issue #4 was an
axis ambiguity, now settled. There are two ORTHOGONAL axes and both hold at once: on the GENOTYPE
axis Myc RAISES nuclear-OXPHOS transcript (d ~1.27, powered); on the TIME axis, within Myc+, nuclear-
OXPHOS subunit transcript DECLINES 6W->12W (subunits slope -0.39, core -0.34). The author's "reduced
abundance" is the TIME-axis decline -- so there is NO self-contradiction. Three points must be baked
into the wording, however:

1. **It is largely MYC-SPECIFIC, not a reflection of normal development.** Decomposing the Myc+ 6->12W
   slope into a shared-developmental (= WT) part + a Myc-specific excess (Issue #6 fade-narration):
   OXPHOS subunits -0.39 = -0.11 shared (28%) + **-0.28 Myc-specific (72%)**; OXPHOS core -0.34 =
   -0.03 shared (10%) + **-0.30 Myc-specific (90%)**. WT OXPHOS barely declines over this window, so
   the drop is 72-90% a Myc-specific DE-AMPLIFICATION of the elevated OXPHOS program -- NOT the gland
   "reflecting the WT developmental trajectory." (Contrast the biosynthetic arm -- amino-acid, lipid
   -- where the WT DOES converge; OXPHOS is precisely the arm where it does not.) This de-amplification
   is the compositional-dilution pointer (the Myc-built high-OXPHOS/proliferative compartment shrinking
   with maturation). The shared/specific SPLIT is powered (gene-level, all genes); the sample-level
   composite slope itself is directional/ns at n=6/tp (Myc+ p ~0.11).
2. **It is a NUCLEAR-subunit decline; mtDNA-encoded OXPHOS RISES 6->12W** (absolute z -1.14 -> +1.23).
   At the transcript level this is a mitonuclear REBALANCING (nuclear-led assembly at 6W ->
   mtDNA-completed at 12W), so "reduced abundance of assembled OXPHOS complexes" does not follow
   cleanly from the transcriptome -- nuclear subunit SUPPLY falls while mtDNA subunits rise.
3. **Whether assembled OXPHOS complex abundance actually falls is a PROTEIN/functional question**
   (blot, respirometry, the Menegollo et al. 2024 companion). If the intro's "reduced abundance of
   OXPHOS complexes" is meant at that level, cite that data; if it is meant transcriptionally, scope it
   to the Myc-specific nuclear-subunit de-amplification over the adult transition, per point 1-2.

**(B) "central to how transformation proceeds" vs the Issue #3 ceiling.** Issue #3 establishes
centrality as COUPLING that survives deconfounding against MYC dose and is present within WT -- a
strong association, explicitly NOT a demonstration of necessity or causality (bulk, n=6/group; a
MYC-readout arm could still be required; direction unidentified). "Central to how transformation
proceeds" reads as mechanistic/causal. Either soften to a coupling statement ("tightly and
dissociably coupled to the tumorigenic phenotype") or keep the word "central" but frame it as the
arm most tightly linked, with causality flagged as the open genetic test.

**(C) "robust increase in mitochondrial biogenesis" as the headline vs dissociable arms (#3).**
Not a contradiction -- biogenesis DOES rise robustly -- but leading with biogenesis buries the
revision's main refinement: biogenesis is the arm that merely TRACKS MYC dose (a bystander), whereas
the phenotype-coupled core is respiratory/biosynthetic (OXPHOS + nucleotide + TCA). For precision,
name the respiratory/biosynthetic core as the phenotype-linked arm, and keep the biogenesis increase
as the broad-footprint context.

**STRENGTHENED, and the word "biogenesis" now has a measured referent (script 32).** Until 2026-07-16
this claim rested entirely on **enrichment scores** -- `bio_comp` is a GSVA composite, mitoPPS is
abundance-blind by construction, and the one object captioned "reflects mito content" was never tested.
A reviewer asking "is that mitochondrial content per cell?" would have found nothing. Script 32 answers
it in **absolute compartment share**: Myc raises mass markers +27% (chaperone-free, d=1.34, p=0.003),
nuclear MitoCarta +19%, mitoribosome +35%, import +33% -- 13/14 mass-panel genes coherent, all panels
LOO-stable, genotype axis QC-CLEAN, and independently concordant with Issue #4 (nuclear OXPHOS d=1.15
vs 1.27 by a different method). So **"biogenesis" can be stated as a bounded content claim rather than
an enrichment claim** -- with two riders: it is a *share of transcriptome*, hence a **LOWER BOUND**
(Myc amplifies total RNA per cell, so a flat share would already mean more mitochondria); and it is
transcript, not protein or organelle volume. The refinement in the paragraph above is unaffected --
biogenesis rising and biogenesis being a MYC-dose bystander are statements about different things
(magnitude vs phenotype-coupling), and script 32 sharpens the first without touching the second.

**And it sharpens the imbalance by an informative NEGATIVE.** The nuclear mito compartment rises by a
median +21% and the mtDNA-dedicated core (Tfam/Polg/Twnk/Polrmt/Tfb2m...) rises **with** it,
unremarkably (39th percentile, Wilcoxon p=0.12) -- while mtDNA-encoded output stays flat (p=0.83). So
the mitonuclear imbalance is **not a transcript-level deficit anywhere in the machinery**: nothing in
the transcriptome is the bottleneck, and the cause must sit **downstream of transcript abundance**
(mtDNA copy number, transcription rate, or turnover). Do not write the narrative as though the
transcriptome explains the imbalance -- it explains the *nuclear* half and demonstrates that the mtDNA
half is not transcriptionally accounted for. **mtDNA qPCR is the discriminating experiment.**
(A "commissioned but unbuilt" reading -- Myc builds the mtDNA machinery, output does not follow -- was
drafted and **retracted the same day** as a set-composition artifact; see the script-32 section.)

### The apoptosis / permissive-window clause -- SUPPORTED (the death spine, Section 1B)

The final clause -- "associated with a desensitisation to MYC-induced apoptosis, opening a permissive
window ... at the adult developmental stage" -- is the best-supported half of the sentence, and it is
what makes the OXPHOS thread matter. The death-timing spine gives it a mechanism, not just an
association: at 6W the mitonuclear-imbalanced substrate is tightly coupled to pro-apoptotic priming
(r 0.75, p 0.005) and Myc engages that priming (~3x stronger at 6W; PRO-ANTI p 0.038); by 12W the
imbalance resolves, the coupling collapses (r -> 0.07), the apoptotic program declines (branch-2 NES
-1.38), and Myc's killing fades -- the desensitisation. The permissive window IS that 12W decoupled
state. So the intro's causal chain (mitochondrial state -> apoptosis desensitisation -> permissive
window) is directly backed. Two guardrails for the wording: (i) the desensitisation is carried by the
mitonuclear-imbalance RESOLUTION and biogenesis-death DECOUPLING, of which the OXPHOS de-amplification
(tension A) is one face -- not by "reduced OXPHOS complexes" acting alone; (ii) it is an ASSOCIATION
in survivor-biased bulk (the model characterises the death-permissive STATE), so phrase the link as
"associated with," exactly as the draft already does. Do NOT overstate to causation.

---

## 4. Proposed reconciled rewrite (a suggestion -- accept, trim, or edit)

> Although MYC is reported to stimulate overall mitochondrial biogenesis, we find that its effect on
> the mitochondrial transcriptome is not uniform. A constitutively expressed Myc amplifies the normal
> pubertal proliferative and biosynthetic program of the mammary gland rather than creating a new one,
> reshaping epithelial lineage identity cell-intrinsically -- repressing the basal-myoepithelial and
> driving the hormone-sensing luminal program -- rather than by selective loss of Myc-high cells; and
> within the mitochondrial compartment it does not act as a uniform biogenic switch: it raises
> mitochondrial biogenesis and absolute nuclear-encoded OXPHOS transcript, yet REALLOCATES resources
> across mitochondrial pathways -- deprioritising OXPHOS RELATIVE to the rest of the compartment (a
> mitonuclear shift from nuclear-led assembly at six weeks to mtDNA-completed metabolism at twelve
> weeks) even as OXPHOS transcript abundance is maintained. The respiratory and biosynthetic core
> (OXPHOS together with nucleotide and TCA-cycle metabolism), rather than mitochondrial biogenesis per
> se, is the arm most tightly coupled to the proliferative, dedifferentiated and apoptosis-primed
> phenotype, along an axis separable from generic MYC-target activity. As the gland matures into its
> adult stage, the Myc+ epithelium de-amplifies this elevated OXPHOS program (a largely Myc-specific
> decline in nuclear-encoded OXPHOS-subunit transcript over six-to-twelve weeks, distinct from the
> minimal wild-type developmental change) [and, at the protein level, a reduced abundance of assembled
> OXPHOS complexes -- COMPANION REF / DIRECT DATA, see tension A]; this lower-OXPHOS adult state is
> associated with a desensitisation to MYC-induced apoptosis, opening a permissive window for
> tumourigenesis.

Per-change rationale:
- "not uniform" kept; split the run-on so each claim stands (readability).
- "amplifies ... rather than creating a new one" -- adds the #1-2 result (the developmental
  interaction the paragraph gestures at).
- "reshaping epithelial lineage identity cell-intrinsically ... rather than by selective loss of
  Myc-high cells" -- adds the Issue #1 PART E adjudication: the lineage reprogramming decomposes into
  per-cell fade (BMYO de-basalisation; off the death axis, rho ~0.14, not pro-death-enriched),
  WT-convergence (the LASP state), and a per-cell LHS gain -- selection/dropout survives only in the
  death-coupled luminal-progenitor arm. Answers the "is the basal repression just dying cells?"
  reviewer worry inside the sentence. Drop the clause if the intro should stay strictly mitochondrial;
  it is defensible but optional. See Issue #1 PART E + Section 5.
- "raises ... absolute nuclear-encoded OXPHOS transcript, yet REALLOCATES ... deprioritising OXPHOS
  RELATIVE ... even as ... abundance is maintained" -- fixes tension A at the transcript level (do
  not imply an absolute mRNA drop).
- "respiratory and biosynthetic core ... rather than mitochondrial biogenesis per se ... separable
  from generic MYC-target activity" -- fixes tensions B and C (names the phenotype-coupled core;
  frames centrality as coupling).
- "de-amplifies this elevated OXPHOS program (a largely Myc-specific decline ... distinct from the
  minimal wild-type change)" -- this is the author's intended TIME-axis reading of "reduced abundance,"
  made precise: it is 72-90% Myc-specific (Issue #6), NOT a reflection of WT development, and it is a
  nuclear-subunit decline (mtDNA rises). See tension A.
- "[at the protein level, a reduced abundance of assembled OXPHOS complexes -- COMPANION REF]" -- keep
  ONLY if there is protein/respirometry data; the transcriptome shows a mitonuclear rebalancing, not an
  unambiguous drop in assembled complexes. If there is no protein source, drop this bracket and rely on
  the transcript-level de-amplification clause alone.
- "desensitisation to MYC-induced apoptosis ... permissive window ... adult developmental stage" kept
  -- now SUPPORTED with a mechanism by the death spine (Section 1B): the 6W mitonuclear-imbalanced
  substrate is death-coupled and Myc engages apoptosis there; by 12W the imbalance resolves, the
  coupling collapses (r 0.75 -> 0.07) and killing fades. Optionally make it explicit: "-- as the
  mitonuclear imbalance that couples the young gland to Myc-induced apoptosis resolves, the adult
  gland is desensitised." Keep "associated with" (survivor-biased bulk = association).
- Optional mechanistic tail for the attenuation: "-- a window that widens as the oncogenic program
  fades on the aging substrate (~2/3) while the wild-type gland converges onto the same biosynthetic
  axis (~1/3)." (#6.)

---

## 5. What remains open (carried-forward ceilings; non-blocking for the writeup)

- **Per-cell fade vs compositional dilution (#6).** The dominant Myc-fade term confounds a per-cell
  weakening with dilution of the shrinking Myc-built proliferative/TEB compartment. Bulk bounds it;
  the E2F/proliferation co-decline points to composition. Settled by deconvolution / single-cell --
  parked plan in `docs/deconvolution_subproject_plan.md`.
- **Lineage repression: per-cell vs death-dropout (#1, script 26 PART E).** BOUNDED, not fully
  settled. BMYO suppression is per-cell (BMYO off the death axis, rho ~0.14; program not pro-death-
  enriched, OR 0.92) -- dropout does not explain it. The residual candidate is the LASP luminal-
  progenitor arm, which IS death-coupled (rho ~0.87) and pro-death-enriched (OR 1.60) and whose Myc
  induction fades -- there, death-driven dropout of Myc-induced cells cannot be excluded. This is the
  same shrinking-luminal-progenitor compartment as the #6 fade; `docs/deconvolution_subproject_plan.md`
  already lists basal + luminal-progenitor as fractions to estimate, so both settle together (single-
  cell / FACS-sorted-basal / deconvolution).
- **Protein-level OXPHOS (tension A).** Whether OXPHOS COMPLEX abundance is reduced is a
  protein/functional question; source it to the blot / Menegollo companion, not this RNA-seq.
- **Mitochondrial content is BOUNDED, not measured (script 32).** Myc raises the mito share of the
  transcriptome (+27% mass markers, chaperone-free, d=1.34, p=0.003; coherent across import/ribosome/
  OXPHOS arms; genotype axis QC-clean and LOO-stable), but bulk polyA with no spike-ins and no cell
  counts cannot give per-cell content in absolute terms -- median-of-ratios normalisation removes that
  scale. And because Myc globally amplifies total RNA per cell, **+27% is a LOWER BOUND**. A
  TOMM20/VDAC/CS blot, mtDNA qPCR, or EM settles it -- the SAME orthogonal experiment tension A
  already needs, so one blot answers both.
- **The mitonuclear imbalance is NOT transcriptionally explained (script 32 PART 3b) -- the sharpest
  open mechanism.** The nuclear mito compartment rises a median +21%, the mtDNA-dedicated core rises
  with it unremarkably (39th percentile, Wilcoxon p=0.12), and mtDNA output stays flat (p=0.83). So no
  transcript-level bottleneck exists anywhere in the machinery: the cause is **downstream** -- mtDNA
  copy number, transcription rate, or turnover, which bulk cannot separate. **mtDNA qPCR is a
  one-experiment test** and would turn the imbalance from a description into a mechanism. (A
  "commissioned but unbuilt" reading was retracted the same day as a set-composition artifact --
  Mrpl12/Atad3a/Poldip2/Top1mt carry the set effects, not Tfam/Polg/Twnk.)
- **Do the genes that follow the mtDNA pattern form a module? (author question, 2026-07-16 -- OPEN.)**
  First pass says the mtDNA machinery does **not** (above). A transcriptome-wide scan is worth doing.
  Its real trap is **not** depth (shares cancel depth by construction; see the corrected QC gate) but
  **cohort structure**: the mt share's biggest swings are between two batch-aligned cohorts, so a scan
  correlating genes against the mt composite across all 24 samples will mostly rank genes by *cohort
  membership*. Design it **within timepoint** (where the mt share still spans 3.4-14.9% at 6W), and
  report the cohort-driven module as the named alternative rather than discovering it at the end. Not
  yet built.
- **The TIME axis is cohort-aligned and the mt arm is imprecise -- but NOT depth-confounded (script 32
  QC gate; corrected 2026-07-16).** **Depth is not a confound**: a share is a proportion and DESeq2's
  size factors cancel depth, both by construction, and within cohort depth predicts nothing (6W
  rho=-0.06 p=0.87; 12W rho=-0.32 p=0.32) -- the pooled rho=-0.47 is a pooling artifact of two clouds.
  12W is also not degraded (top100 31% vs 40%). What is real: (i) the 6W/12W depth ranges are
  **disjoint** (18.8-29.0M vs 9.9-16.1M) so batch is perfectly aligned with timepoint and not
  separable -- **inherent to any cross-sectional design**, not an mt-specific defect, no artifact
  indicated; (ii) the mt-* share is **high-variance** (3.4%-40%), so temporal mt claims incl. Issue
  #4's "mtDNA absolute z -1.14 -> +1.23" are **imprecise** at n=6. Genotype contrasts are clean
  (depth-balanced p=0.41, litter-controlled). Worth one methods sentence on the batch alignment; do
  not overstate it as a confound.
- **Causality (#3, #5).** The phenotype-coupled respiratory core is orthogonal-to-MYC-dose, not
  proven MYC-independent or necessary; the identifying tests are inducible Myc, TF perturbation, and
  ATAC/ChIP.
- **The death spine is survivor-biased association (Section 1B).** We sequence the cells that did NOT
  die, at one timepoint per age, n=6 -- so the death-permissive STATE is characterised but not the
  causal act of killing; the mito / stem / biogenesis death-coupled axes are not separable at this n.
  BH3 profiling, caspase-by-state, or single-cell would resolve it. The external Myc-ER inducible
  phenotype (kills more at 6W) is the anchor these bulk data explain but do not themselves reproduce.

---

## 6. Gene-set glossary (provenance and composition)

Sources: `data/genesets_from_library/provenance_table.csv` (per-set source/category/method/size),
the master GMT `mammary_mito_myc_metab_v1_mouse.gmt` (986 sets; N = gene count per set),
`gray_chea_mito_tf_shortlist.csv`, and `docs/library_reference/` catalogs. `N` = mouse gene count.
Composites are built IN the scripts (line-cited) by averaging per-sample GSVA scores over member
sets; they are NOT new gene sets.

### Library sets referenced

| Set / family | Source | Category | N | What it covers |
|---|---|---|---|---|
| `MITOCARTA_ALL` | MitoCarta3.0 Sheet 4 | All | 1140 | Whole mitochondrial compartment (all MitoCarta genes). |
| `MITOCARTA_OXPHOS` | MitoCarta3.0 Sheet 4 | OXPHOS_nu | 155 | Nuclear-encoded OXPHOS (verified 0 mt-* genes). |
| `MITOCARTA_OXPHOS_SUBUNITS` | MitoCarta3.0 | MitoCarta pathway | 89 | Structural OXPHOS subunits (incl. 13 mtDNA-encoded; nuclear subset = 89 - 13). |
| `MITOCARTA_OXPHOS_ASSEMBLY_FACTORS` | MitoCarta3.0 | MitoCarta pathway | 67 | OXPHOS complex assembly factors. |
| `MITOCARTA_TCA_CYCLE` | MitoCarta3.0 | MitoCarta pathway | 20 | TCA-cycle enzymes. |
| `MITOCARTA_NUCLEOTIDE_METABOLISM` | MitoCarta3.0 | MitoCarta pathway | 36 | Mitochondrial nucleotide metabolism. |
| `MITOCARTA_AMINO_ACID_METABOLISM` | MitoCarta3.0 | MitoCarta pathway | 90 | Mitochondrial amino-acid metabolism (biosynthetic arm). |
| `MITOCARTA_LIPID_METABOLISM` | MitoCarta3.0 | MitoCarta pathway | 126 | Mitochondrial lipid metabolism (biosynthetic arm). |
| `MITOCARTA_MITOCHONDRIAL_RIBOSOME` | MitoCarta3.0 | MitoCarta pathway | 83 | Mitoribosome (biogenesis/translation arm). |
| `MITOCARTA_APOPTOSIS_PRO` / `_ANTI` | MitoCarta3.0 | MitoCarta pathway | 25 / 9 | Pro- / anti-apoptotic mito genes; the death-`priming` axis = PRO - ANTI (Issue #3 outcome; death spine H1/H2). |
| Tang 2024 15-RCD modalities | Tang et al. 2024 (Comput Struct Biotechnol J), human -> mouse orthologs | cell death (`data/cell-death/*.csv`) | 15 sets | Regulated-cell-death modalities (Apoptosis, Ferroptosis, Necroptosis, Pyroptosis, Autophagy-dep, Cuproptosis, Disulfidptosis, etc.); branch-2 fGSEA (script 12). |
| `cell_death_genes_consolidated` | Consolidated RCD list, mapped via MyGene.info + HomoloGene (not remapped) | cell death (`data/cell-death/`) | -- | Pro-/anti-death annotated gene list for the branch-1 directional binomial (script 16). |
| `MYC_HALLMARK_MYC_TARGETS_V2` | Felsher MYC compendium GMX (snapshot) | Felsher compendium | 54 | MSigDB Hallmark MYC targets V2. |
| `MYC_felsher_integrative_signature` | Felsher integrative signature (snapshot) | Felsher integrative | 61 | Mito-STRIPPED MYC phenotype core (8% mito, 0% OXPHOS) -- the deconfounding proxy in #3. |
| `MYC_signatures` (category) | Felsher compendium (e.g. `MYC_MUHAR_MYC_SIGNATURE`, 88) | MYC_signatures | 17 sets | The `myc` composite: 17 MYC-target/signature sets. |
| `Proliferation` (category) | MSigDB (KEGG/Reactome/WP/Hallmark: `PROLIF_*`) | Proliferation | 14 sets | The `prolif` composite: cell-cycle / E2F / G2M / DNA-replication sets. |
| `PROLIF_E2F_HALLMARK` | MSigDB | proliferation: E2F targets | 200 | Hallmark E2F targets (member of `tf_e2f`). |
| `MG_*` developmental sets | GRAY (Gray et al.) | Mammary_development | 179 sets | Per-lineage developmental signatures (e.g. `MG_TEB_VS_DUCTAL_*_GRAY_UP/DN`, 150 each = TEB vs mature duct); GRAY UP/DN and CHUNG ATAC OPEN/CLOSED nets. HEVSLE = High-/Low-EXPRESSION (differentiation state), NOT estrogen. |
| `ESRRA_MITO` / `NRF1_MITO` / `GABPA_MITO` | INTERSECTION (TF targets AND MitoCarta) | component: PGC1a axis | 177 / 81 / ~ | ER/PGC1a mito-biogenesis TF activity ON the mito compartment; the `tf_biogenesis` axis. |
| `ESR1_AND_CORE_MITO` | INTERSECTION | discrimination: ER-and-core | 59 | ESR1 targets AND mito core (member of `tf_esr1`). |
| `TFT_<TF>_GRAY_*_MITO` | GRAY developmental TF targets AND MitoCarta | TF-in-mito | varies | Gray developmental-TF target sets restricted to mito; feed `tf_e2f`, `tf_esr1`, `tf_lipogenic`, `tf_mito_panel`. |
| `gray_chea_mito_tf_shortlist.csv` | Gray + ChEA (library) | mito-TF shortlist | -- | Ranked mito-TF list; source of the data-driven `tf_mito_panel` (top-20, excl. MYC/E2F/etc). |
| `GS_METAB_*` (e.g. `GS_METAB_KREBS` 27, `GS_METAB_NUCLEOTIDE` 62) | GS_metabolic (pre-mapped) | metabolic sub-pathways | varies | Curated central-metabolism sets; feed the 12 metabolic axes below. |

### Script-defined composites (per-sample GSVA averages; not new sets)

- **Issue #2 (script 27:88-89).** `myc` = MYC_signatures (17); `felsher` = MYC_felsher_integrative
  (61); `hallmark_v2` = MYC_HALLMARK_MYC_TARGETS_V2 (54); `prolif` = Proliferation (14); `myc_in_teb`
  = `TFT_MYC_GRAY_*_TEB`; `teb_ductal` = MG_TEB_VS_DUCTAL UP - DN.
- **Issue #3 (script 28:90-121).** `mito_all` = MITOCARTA_ALL; `mito_oxphos` = regex over
  OXPHOS|COMPLEX|_SUBUNITS|ASSEMBLY_FACTORS|ELECTRON_CARRIERS|CRISTAE; `mito_biogenesis` = regex over
  RIBOSOME|CENTRAL_DOGMA|MT_TRNA|MT_RRNA|MTRNA|MTDNA|IMPORT|TRANSLATION. Twelve metabolic axes
  (`metab_axes`): glycolysis, tca, oxphos_met, ppp, nucleotide, one_carbon, glutamine, fao, fa_synth,
  amino_acid, cholesterol, redox -- each a composite of GS_METAB_* + KEGG/Reactome/WP/Hallmark sets.
  Outcomes: `mb2_fork` (AP7 MB2/MB1), `priming` (APOPTOSIS_PRO - ANTI), `teb_dediff` (MG_TEB UP - DN).
- **Issue #4/#5 (scripts 29:266-437, 30:104-133).** `oxphos_abs` = VST composite over nuclear OXPHOS
  subunit Ensembl IDs (subunits minus 13 mt); `dev_luminal` = GSVA composite over LASP/LHS luminal
  sets; `tf_biogenesis` = ESRRA/GABPA/NRF1_MITO; `mtdna_reprior` = mtDNA-encoded OXPHOS mitoPPS.
- **Issue #6 (script 31:113-199).** Roster as listed in Issue #6 above; broadened TF axes `tf_e2f`,
  `tf_esr1`, `tf_lipogenic`, `tf_mito_panel` as defined above.
- **Death spine (scripts 23-25).** `priming` = MITOCARTA_APOPTOSIS_PRO minus ANTI (pro-apoptotic
  balance); `mitonuclear_imbalance` = nuclear minus mtDNA-encoded OXPHOS mitoPPS (the death-permissive
  substrate feature, = the Issue #4 mitonuclear axis); biogenesis/stem (MASC) composites for the
  per-sample death-coupling in script 25 Part B. Branch 1 (script 16) uses
  `cell_death_genes_consolidated`; branch 2 (script 12) uses the Tang 15-RCD sets.
- **Mito content (script 32:120-190).** NOT GSVA -- these are **compartment shares**: `100 *
  colSums(raw counts over the panel) / colSums(raw counts)`, reported on two denominators (whole
  transcriptome; and one excluding the 13 mt-* genes, which are 76.9% of all MitoCarta counts). Panels
  are library GMT sets except **`MASS_MARKERS`**, the one hand roster in the corpus: `Vdac1/2/3,
  Tomm20/22/40, Tspo` (OMM_structural) + `Cs, Immt, Slc25a3, Timm23/44` (matrix_IMM) + `Hspa9, Hspd1`
  (chaperone). A **flagged exception** to CLAUDE.md's "do not rebuild gene sets" -- the library has no
  mass-marker equivalent; it is the in-silico stand-in for a TOMM20/VDAC/CS/HSP60 blot. The
  claim-bearing composite is `MASS_MARKERS_NOCHAP` (chaperones excluded because Hspa9/Hspd1 are
  standard mass markers **and** direct MYC targets, so they have a confounded route to elevation).
  Shares are valid BETWEEN samples for a fixed set (gene length cancels), never BETWEEN panels.

---

*This walkthrough is the teaching/writing layer; the archival record stays in
`docs/2026-07-08_BlockA_revision_plan.md`. Numbers trace to that log and to scripts 26-31; no
analysis was re-run to produce this doc.*
