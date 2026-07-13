# Block A revision, point by point -- what we did, what we found, and how it aligns with the introduction

Written 2026-07-13, branch `paper-figures` (@ 5c00c5e). Audience: the author, building the
manuscript narrative. This is the TEACHING layer for the step-by-step Block A revision (Issues
#1-6, scripts 26-31). It restates each issue in plain terms with the effect sizes, defines every
gene set precisely (source + what genes it covers; full glossary in Section 6), and then confronts
the whole body of findings against the broad-overview paragraph drafted for the Introduction,
ending with a proposed reconciled rewrite.

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
  carried by cross-modality agreement, not the p-value); LASP near-null. Corroborated by the GRAY
  directional nets (Myc -> less-differentiated / TEB-proliferative end) and CHUNG ATAC (Myc closes
  basal chromatin, opens luminal by 12W).

**Interpretation + ceiling.** At 6W Myc imposes a broad off-axis identity suppression; by 12W a
specific **anti-BMYO / pro-LHS** reprogramming. Main effects powered (24 samples); the interaction
is directional at n=6/tp; composites use correlated sets. Association, not causation.

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
- **D -- first-pass regulators (hypothesis-generating, within WT, n=12).** The WT OXPHOS trajectory
  tracks BOTH the developmental luminal axis (partial r 0.71) and ER/PGC1a TF activity
  (ESRRA/NRF1/GABPA; partial r 0.70-0.83), model R^2 0.86-0.89.

**Interpretation + ceiling.** The attenuation is a shrinking of MAGNITUDE, not of program IDENTITY.
In absolute mRNA Myc RAISES nuclear OXPHOS; do not say Myc "reduces OXPHOS" -- say it REALLOCATES the
mitochondrial compartment. Part B is powered (genotype main effect); the age decline and Part D are
directional / within-WT associations. The retraction did not affect the committed Issue #3 result
(partial-correlation centrality does not depend on it).

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

## 2. The through-line (Issues #1-6 as one arc)

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
collapse of proliferative/mito-TF activity -- a compositional hint that opens a permissive adult
window for tumourigenesis (#6).

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
| 4 | "MYC ... interacts with the transcriptional trajectory of normal mammary development" | #1, #2 | SUPPORTED -- Myc bends off the maturing WT axis; amplifies the endogenous pubertal program. "plays a part and interacts" is well-calibrated. |
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

---

## 4. Proposed reconciled rewrite (a suggestion -- accept, trim, or edit)

> Although MYC is reported to stimulate overall mitochondrial biogenesis, we find that its effect on
> the mitochondrial transcriptome is not uniform. A constitutively expressed Myc amplifies the normal
> pubertal proliferative and biosynthetic program of the mammary gland rather than creating a new one,
> and within the mitochondrial compartment it does not act as a uniform biogenic switch: it raises
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
- "permissive window ... adult developmental stage" kept (SUPPORTED by #6). Optionally add the
  mechanism: "-- a window that widens as the oncogenic program fades on the aging substrate while the
  wild-type gland converges onto the same biosynthetic axis." (~2/3 fade, ~1/3 convergence; #6.)

---

## 5. What remains open (carried-forward ceilings; non-blocking for the writeup)

- **Per-cell fade vs compositional dilution (#6).** The dominant Myc-fade term confounds a per-cell
  weakening with dilution of the shrinking Myc-built proliferative/TEB compartment. Bulk bounds it;
  the E2F/proliferation co-decline points to composition. Settled by deconvolution / single-cell --
  parked plan in `docs/deconvolution_subproject_plan.md`.
- **Protein-level OXPHOS (tension A).** Whether OXPHOS COMPLEX abundance is reduced is a
  protein/functional question; source it to the blot / Menegollo companion, not this RNA-seq.
- **Causality (#3, #5).** The phenotype-coupled respiratory core is orthogonal-to-MYC-dose, not
  proven MYC-independent or necessary; the identifying tests are inducible Myc, TF perturbation, and
  ATAC/ChIP.

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
| `MITOCARTA_APOPTOSIS_PRO` / `_ANTI` | MitoCarta3.0 | MitoCarta pathway | 25 / 9 | Pro- / anti-apoptotic mito genes; `priming` = PRO - ANTI. |
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

---

*This walkthrough is the teaching/writing layer; the archival record stays in
`docs/2026-07-08_BlockA_revision_plan.md`. Numbers trace to that log and to scripts 26-31; no
analysis was re-run to produce this doc.*
