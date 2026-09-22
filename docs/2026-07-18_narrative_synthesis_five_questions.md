# Narrative synthesis — the five questions, answered

Written 2026-07-18 on `paper-figures` (scripts 32–35, commits `50ad21c`..`022e217`). Read after
`docs/2026-07-17_evidence_audit_and_narrative.md`. Every number is reproducible from a committed
script; nothing is quoted from conversation only.

---

## 0. The one thing that explains why results changed

**Almost every composite in this dataset correlates with almost every other, and nobody had
measured how much.**

The median |rho| of a mitochondrial/metabolic axis to an **arbitrary** one of the 884 library
programmes (script 35 PART A):

| axis | all (n=24) | 6W (n=12) | 12W (n=12) |
|---|---|---|---|
| mito_oxphos | 0.80 | 0.82 | **0.78** |
| nucleotide | 0.77 | 0.82 | 0.71 |
| tca | 0.77 | 0.80 | 0.78 |
| mito_biogenesis | 0.74 | 0.78 | 0.71 |
| **redox** | **0.22** | **0.41** | **0.23** |

**A published coupling of 0.8 is what the window hands you.** It is not a finding.

**Why the ceiling is there.** Within a timepoint the **6-vs-6 genotype split** alone makes every
Myc-responsive programme co-vary with every other. The IEG/dissociation axis (script 33) adds to
it at 6W. `redox` is the control that proves it: **not a Myc target (genotype d=0.08)**, and its
ceiling is half the others'. Readable axes are the ones Myc does not drive.

**Two corrections to earlier framings, both mine, both caught by testing:**
1. *"6W is stressed, 12W is clean"* — **wrong**. Script 34 measured the ambient collapsing
   0.71 → 0.35, but that was the mitonuclear **imbalance**, a mitoPPS **ratio** — and a ratio
   cancels the shared component. Per axis, the ceiling barely moves (0.82 → 0.78). **Never borrow
   one axis's ambient for another.**
2. A **raw** ambient is the wrong null for Issue #3 Q2 anyway, because Q2's claim is about the
   **partial** correlation. The matched null is a **partial ceiling** (script 35 PART B).

**What the ceiling cannot touch: genotype CONTRASTS.** The IEG axis is genotype-**independent**
(p=0.49), so it cancels inside every Myc-vs-WT comparison. That is the dividing line running
through all five answers: **contrasts survived; cross-sectional couplings did not.**

---

## 0.5 Where the ceiling comes from — one global per-sample axis

The author's two method notes (`docs/GSVA_global_background_and_pathway_coupling.md`;
`docs/Gene_set_quantification_for_pathway_correlation.md`) named the ceiling as a **common-mode
latent factor** and prescribed the fix. **Script 36 implements it and reproduces the origin probes
as saved fields** (`linear_pathway_coupling.rds$ceiling_probes`; numbers below are the linear
z-score, GSVA in parentheses):

| candidate | probe | result |
|---|---|---|
| low n | chance floor (permute the axis's sample labels) | **0.13** vs observed **0.80** → NO |
| set curation / overlap | median \|rho\| of **random set pairs** | **0.54** (0.53) → not curation |
| the Myc genotype split | ambient **within `6W_pos` alone** (contrast gone) | **0.83** → NO, survives it |
| **one global per-sample axis** | `cor(mito_oxphos, mean of all 884 scores)` | **0.968** (0.955) |
| " | ambient after regressing the global mean out of every set | **0.13** ≈ chance |

**So the ceiling is a single global per-sample offset that every score rides.** GSVA scores each
sample against the cohort but does **not centre within a sample**, so any broad shift in a mouse
lifts every set score together. **Methodological in FORM, biological + technical in ORIGIN.**

**Biological or methodological? Both — and the decomposition is measurable** (script 36 PART C).
The **raw** PC1 explains **76%** of variance with **79%** of loadings one-signed (a common mode,
not a contrast) and correlates with the group label at **r=0.48** — so a large part of the raw
ceiling **is** the four-group design. The **residual** PC1 (design removed first, doc 1 §7) still
explains **70%**, and its strongest correlates are **epithelial purity +0.79**, **contamination
−0.70**, **proliferation +0.65**, **MYC activity +0.44**, **IEG/prep −0.31**. That is
**composition + prep (technical) AND proliferation/MYC (biological), together** — neither cleanly.
For enzymatically dissociated, non-FACS MECs, an epithelial-purity/stromal-contamination axis is
exactly the expected prep variation. **It is not noise, and a bigger n or better prep would not
remove it** — it is the thing a raw coupling mostly measures. Report it as a phenotype **and**
remove it for specific coupling. Figures: `linear_pathway_coupling/A_ceiling_origin.pdf`,
`B_global_factor_phenotype.pdf`.

---

## 0.6 Which instrument, and do you have to use the floor

**Do you have to use the floor to generate a hypothesis? No.** Prior + plausibility + a suggestive
number is enough when the bench does the testing. **The problem was never r=0.75 — it was p=0.005**
(which asserts a test was passed) and the claim level around it. *But a ceiling-level correlation
did not nominate the hypothesis; the prior did.* Label such leads **prior-driven** — legitimate and
honest. The floor is a **ranking tool, not a filter**, and it changes the ranking.

**The instrument itself was wrong for coupling** (doc 2). GSVA's nonlinear competitive score has no
clean covariance reading. The coupling instrument of record is now the **linear mean gene-wise
z-score** (its correlation *is* average cross-gene covariance); **singscore** is the
cohort-independent robustness arm; **GSVA** is kept for continuity. Gene-level vs score-level design
adjustment agree to ~0.02 (script 36 `gene_vs_score_design`), so the correction is robust. The
correction is **design → PC1 → residual, with the factor estimated AFTER removing the design**
(doc 1 §7) and **leave-pair-out** for each tested pair (§4) — **not** mean-subtraction, which
assumes equal loadings and imposes a sum-to-zero constraint (that was the script-35-era
instability).

**The interpretability trap this exposes (script 36 PART D).** When an axis *is* the global factor
— `mito_oxphos` carries **r²=0.94** with the global mean — removing that factor leaves ~**4%** of
its variance, which is **noise**. Its global-adjusted "coupling" to proliferation prints **−0.62**
and is sign-consistent across all three quantifiers, **but this is noise-on-noise at n=24, not
antagonism.** The honest statement of "is OXPHOS central?" is therefore sharper than Q4's:
**OXPHOS has essentially no variance independent of the broad state, so no separable central
coupling is testable at all** — which *confirms* "OXPHOS-central fails," and supersedes script 35's
raw-ambient framing.

**What the principled correction actually nominates** (`linear_pathway_coupling.rds$robustness`,
global-adjusted, three quantifiers):

| lead | global-adj (z) | across methods | axis r² w/ global mean | verdict |
|---|---|---|---|---|
| **redox → MB2 fork** | **−0.63** | −0.62 … −0.73, same sign | **0.03** (90% survives) | **robust, testable — the readable lead** |
| cholesterol/mevalonate → MB2 fork | +0.14 | +0.14 … +0.27 | 0.54 | separable but **weak** — *demoted* from script 35's 96th-pct headline |
| mito_oxphos / mito_biogenesis → prolif, teb | (noise) | — | 0.78–0.94 | **untestable — axis IS the global factor** |
| mito_oxphos → priming (neg. control) | (noise) | sign flips | 0.94 | untestable + circular |

**Lower redox/glutathione capacity tracks the tumorigenic MB2 fork, independently of Myc dose** —
the one coupling with real independent variance and cross-quantifier agreement. It is a **ranking
lead for the bench, not a result** (n=24; nothing beats BH<0.05).

**The replication trap.** A ceiling-level lead **will "replicate" in any similar bulk dataset**,
because the ceiling is a property of the DESIGN + SCORING, not the biology. **Do not validate with
more bulk — validate orthogonally** (BH3 profiling / mtDNA qPCR / redox assay / the blot).

**The sentence-level fix pattern** (keeps the narrative honest without deleting it):
> was: *"the imbalance couples to pro-apoptotic priming (r=0.75, p=0.005) — the causal spine…"*
> becomes: *"the imbalance and pro-apoptotic priming co-vary at 6W (r=0.75). Couplings of this size
> arise between most programme pairs in this design, so we treat this as hypothesis-generating and
> test it directly by BH3 profiling (Fig X)."*

**Delete the p, keep the r, name the bar, point at the experiment.**

### Hypothesis ledger
- **Data-nominated, testable:** `redox → MB2 fork` (robust, Myc-independent — the strongest lead);
  the `imbalance → priming` **collapse** at 6W→12W (96th pct of Δrho, perm p=0.066 — Pearson-only,
  circular axis, so a *lead* not a finding).
- **Data-nominated, weak:** `cholesterol/mevalonate → MB2 fork` (testable but +0.14 once the factor
  is removed) — still worth a statin/mevalonate bench look given the companion paper's fork.
- **Prior-nominated, data-consistent:** `Bbc3`/PUMA (right BH3 biology; 1-of-23 is chance, but a
  fine BH3-profiling hypothesis).
- **Neither:** "OXPHOS is central" (untestable — axis is the global factor).

---

## Q1 — Do the attenuation mechanisms hold? **YES, with one refinement.**

The attenuation is the most robust mechanistic result in the paper and nothing this week touched
it. It is a **gene-level genotype × time contrast**, so the ceiling and the prep axis are both
irrelevant to it.

- **It is real, not a power artifact** (script 33 PART E2): 12W is **not** noisier (median
  within-group SD **0.2935** vs **0.3515**), and Myc-contrast SEs are ~identical (**0.251** at 6W
  vs **0.259** at 12W). So the DE collapse **1967 → 135** tracks median |LFC| falling
  **0.268 → 0.204** *at constant SE* — genuine magnitude loss. Noise inflates |LFC| rather than
  shrinking it, so it could not manufacture this in any case.
- **The two-mechanism account holds** (Issue #6, exact identity): **~66% Myc-fade + ~34%
  WT-convergence**, with convergence confined to the **biosynthetic arm**; OXPHOS and the
  MYC-target core close by **fade alone**.
- **The fade is mostly Myc-specific, not developmental** (OXPHOS 72–90%, proliferation ~100%) —
  so "Myc+ rides the shared developmental decline and falls **further**; the further part IS the
  attenuation."

**REFINE — "the MYC identity core is selectively BUFFERED".** It attenuates less than the global
average (`atten` **0.168** vs **0.225**) but **not less than expression-matched genes**
(script 34 `atten_excess`: z=**0.29**, p=0.78). Same for OXPHOS (z=0.90, p=0.37). **The buffering
tracks its expression level, not a protective mechanism.** Nothing attenuates more than matched
background — **the attenuation is global**, which is Issue #4's own finding.

**Figures:** `attenuation_decomposition/A_nes_paradox_two_rulers.pdf` (magnitude not rank — the
NES/|LFC| paradox); `attenuation_mechanism/A1_convergence_fade_by_program.pdf` (**headline**: the
66/34 split); `A2_conv_fade_gene_map.pdf`; `A3_conv_fade_shares.pdf` (robustness);
`B1_tf_bint_before_after.pdf` + `B2_tf_delta_bint_bootCI.pdf` (no TF axis absorbs it);
`B3_tf_temporal_vs_oxphos_fade.pdf` (the compositional pointer).

---

## Q2 — Is there mito reprioritisation? **GENOTYPE yes — TIME no.**

**GENOTYPE — holds, and is the strongest new thing in the paper.**
Myc raises mitochondrial **content +21–27%** (script 32; **+21.1%, p=0.0006** after adjusting for
prep-stress + contamination; mass panel coherent **13/14 genes up, 11 at p<0.05**; all 17 panels
leave-one-out stable), raises **absolute nuclear OXPHOS** (d=**1.27** script 29 / d=**1.15**
script 32 — two independent methods), *while* mitoPPS shows a **within-compartment
REALLOCATION**. So: **Myc REALLOCATES the mitochondrial compartment; it does not "reduce OXPHOS".**
All genotype contrasts ⇒ the ceiling and the IEG axis cannot bias them.
*Rider:* every number is a **share of transcriptome**, and Myc amplifies total RNA/cell, so
**+21–27% is a LOWER BOUND**.

**TIME — unsupported.** The mitonuclear imbalance has **no supported contrast** (time p=**0.313**,
genotype p=**0.523**, interaction p=**0.916**, Myc@6W p=**0.249**; adjusting for the IEG axis takes
the time beta **−0.586 → −0.011**). The quoted pattern (0.130 / 0.499 / −0.456 / −0.173, *"peaks
at 6W_pos, inverts by 12W"*) is **four group means that never had a contrast fitted**. And the mt
axis **is** time-associated (mt% p=0.023; IEG p=0.0003), so the scripts 22/24/29 mtDNA-rise /
"nuclear-led assembly → mtDNA-completed running metabolism" story is entangled with the prep axis
and **not adjudicable** on these data. *Unsupported is not refuted* (n=6/group). **mtDNA qPCR
settles it** at the copy-number level the transcriptome cannot see.

**Figures:** `mito_content_proxies/A_compartment_shares.pdf` (**headline**),
`B_myc_vs_age_axes.pdf`, `C_qc_gate.pdf`; `attenuation_decomposition/C_mitopps_vs_absolute.pdf`
(**the reconciliation panel** — reprioritisation vs absolute), `B_absolute_composite_trajectory.pdf`
(3-normalisation robustness), `B_nuclear_oxphos_heatmap_vst.pdf`.
**Do NOT use as claims:** `dev_substrate_death/A1_mito_maturation_imbalance.pdf`,
`biogenesis_discrimination/A2_mitonuclear_imbalance.pdf` — four group means, no contrast.

---

## Q3 — Is there evidence for Bbc3 attenuation? **NO — it is exactly the chance expectation.**

`Bbc3`/PUMA on the interaction: **LFC −0.541, z=−2.65, p=0.008, padj=0.844** — the only BH3 gene
with a nominal interaction. **With 23 genes tested, chance predicts 1.15 nominal hits. Bbc3 *is*
the expected one.**

The panel does not corroborate it — and it should, if a BH3 rheostat were shifting:

| gene | role | baseMean | interaction LFC | z | p |
|---|---|---|---|---|---|
| **Bbc3** (PUMA) | BH3-only | 152 | −0.541 | **−2.65** | 0.008 |
| Bax | effector | 481 | −0.231 | −1.54 | 0.12 |
| Bid | BH3-only | 148 | −0.461 | −1.20 | 0.23 |
| Bak1 | effector | 472 | −0.302 | −1.06 | 0.29 |
| **Bcl2l11** (BIM) | BH3-only | **2310** | −0.038 | **−0.17** | 0.86 |
| **Pmaip1** (NOXA) | BH3-only | 197 | **+0.188** | +0.37 | 0.71 |

**BIM — the highest-expressed BH3 gene in the panel — is flat, and NOXA moves the wrong way.**
Bbc3 remains the **right hypothesis** (BH3-only, direct MOMP trigger, canonical Myc target). It is
a hypothesis, not evidence.

### What else could give "less priming at 12W"

1. **The Myc programme fades** — priming attenuates with everything else; no death-specific gate
   (script 34: the death sets are **not** in excess of an expression-matched null, z=0.97/0.71).
   This is the **parsimonious** account and it is canonical Evan/Lowe Myc-apoptosis biology.
2. **Compositional dilution** — the proliferative/TEB compartment Myc built at 6W shrinks
   (Issue #6's B6 pointer). Fewer high-Myc cells ⇒ less priming, with no per-cell change.
3. **Survivor selection** — the primed cells died; we sequence the survivors (H4: sd_pro ratio
   **0.42** vs **0.60**, Myc+ narrowing more). Script 34 found **all 4 non-mito death axes have a
   NEGATIVE Myc effect** (APOP_BH3_REACTOME **−0.367**, p=0.031; CDC_PRODEATH_APOPTOSIS
   **−0.149**, p=0.015) — a suppressed-looking survivor pool is exactly what OIA leaves behind.
   Survivor bias may **invert** the sign, not merely shrink it.
4. **A protein-level BCL2 rheostat** — transcript ≠ protein ≠ MOMP threshold.
5. **NOT p53/ARF** — refuted (p=0.457 / 0.799).

**None of these are separable in bulk transcriptome.** **BH3 profiling** separates (4) from the
rest and measures the gateway directly; **single-cell/deconvolution** separates (2) from (3).

**Figures (as descriptive substrate, NOT as mechanism):** `death_timing/h1_bh3_bcl2_rheostat.pdf`,
`h1_bcl2_family_wt_forest.pdf`; `death_priming_reassessment/D_attenuation_excess.pdf` (the
mechanism discriminator — nothing in excess).

---

## Q4 — Does "OXPHOS not biogenesis is central" hold? **HALF of it.**

This is the answer that changed most, and the change is not in the direction I expected.

**(1) The raw coupling is at the ceiling in every window.** mito_oxphos → proliferation is
**+0.86** (n=24) against a raw ceiling of **0.80**; **+0.83 vs 0.82** at 6W; **+0.77 vs 0.78** at
12W. A mito axis correlates with ~everything, at both timepoints.

**(2) The 'OXPHOS is CENTRAL' half FAILS its matched null.** Q2's claim was never about the raw
coupling — it was that the coupling **survives removing MYC dose**. So the null must be partial
too (script 35 PART B): partial rho(axis, *every* library set | MYC), with the outcome located in
it.

| coupling (all, n=24) | partial \| MYC | partial ceiling | percentile | perm p |
|---|---|---|---|---|
| mito_oxphos → prolif | **+0.59** | **0.45** | 68 | 0.32 |
| mito_oxphos → teb | +0.44 | 0.45 | 49 | 0.51 |
| nucleotide → prolif | +0.54 | 0.35 | 82 | 0.18 |
| **mito_biogenesis → prolif** | **+0.11** | 0.22 | **27** | 0.73 |
| **mito_biogenesis → teb** | **−0.03** | 0.22 | **7** | 0.93 |
| **cholesterol → MB2 fork** | **+0.40** | **0.13** | **96** | **0.043** |
| redox → MB2 fork | −0.50 | 0.22 | 87 | 0.13 |
| *mito_oxphos → priming* (neg. control) | +0.32 | 0.45 | — | 0.64 |

OXPHOS retains **residual co-variation with everything** after MYC is removed (ceiling 0.45), and
proliferation is not special within that. **"Central to phenotype" implies PREFERENTIAL coupling
to phenotype; that is not what the data show.**

**(3) The 'biogenesis is a BYSTANDER' half HOLDS — harder than published.** Its partial couplings
sit **BELOW its own ceiling** (27th and 7th percentile). **After MYC dose is removed, biogenesis
couples to nothing.** This is the well-supported half and it should lead.

**(4) The only couplings approaching specificity are the NON-Myc ones — and script 28 already
found them while the narrative led with the mito arms.** `cholesterol/mevalonate → MB2 fork`
(partial **+0.40** against a ceiling of only **0.13**, 96th percentile, p=0.043) and
`redox → MB2 fork` (−0.50 vs 0.22). Their ceilings are low **precisely because they are not Myc
targets** — which is why they are readable and the mito arms are not.

**(5) Multiplicity:** **1 of 96** tests beats BH<0.05; cholesterol → fork is p=0.043 raw but
**BH=0.35**. **Nothing here is a significant finding** at n=24 with correlated sets.

**(6) The principled linear correction (script 36) CONFIRMS and SHARPENS this — and re-ranks the
leads.** Under design→PC1→residual on a **linear** score with **leave-pair-out** (see §0.6): the
mito axes carry **r²≈0.94** with the global mean, so they have essentially no variance independent
of the broad state — **"OXPHOS is central" is not merely at the ceiling, it is UNTESTABLE**, which
is the cleanest possible statement of its failure. **`redox → MB2 fork` is promoted to the one
robust, testable lead** (global-adjusted **−0.63**, same sign −0.62…−0.73 across all three
quantifiers, axis r²=0.03), while **`cholesterol/mevalonate → fork` is DEMOTED** to +0.14
("separable but weak") once the factor is removed — its script-35 96th-percentile headline was a
raw-ambient artifact. The biogenesis-bystander call is unchanged.

**Restate as:** *"Mitochondrial biogenesis is a MYC-dose bystander — once MYC dose is removed it
tracks nothing. The OXPHOS/nucleotide arm retains residual phenotype co-variation beyond MYC dose,
but not preferentially."* **Drop the death-priming third entirely** — `priming` is defined at
`28:125` as `MITOCARTA_APOPTOSIS_PRO − MITOCARTA_APOPTOSIS_ANTI`, a **MitoCarta** set, so
"mito → priming" is **mito genes vs mito genes** by construction (flagged in
`2026-07-08_BlockA_revision_plan.md:236`, never acted on).
*The OXPHOS-vs-biogenesis ordering (+0.59 vs +0.11) is real and descriptive — but it was never
tested AS A DIFFERENCE, and two percentiles are not a contrast.*

**Figures:** `myc_mito_centrality/q2_partial_correlation.pdf` (**headline** — but caption the
bystander half, not the central half), `q1_preferentiality_reframed.pdf` (the size-fair Q1
reframe), `q1_genotype_effect_ranking.pdf`, `metabolic_centrality.pdf`;
`ambient_corrected_couplings/B_excess_over_ambient.pdf` (**new** — the ceiling correction).
**Use `q2_coupling_heatmap.pdf` with care:** it shows raw couplings, which are all at the ceiling.

---

## Q5 — Does the Myc × developmental-trajectory effect hold? **YES for the powered half.**

**Holds — genotype main effects (24 samples), which the ceiling cannot touch:**
- The transgene **AMPLIFIES** the pubertal programme: Myc-target activity **d 2.3–3.0, p<1e-5**;
  proliferation d1.0 p0.02; myc_in_teb d1.6 p6e-4. The TEB phenotype is amplified **specifically
  at 6W** (**+0.41** vs **+0.05**, ~4×) — TEBs are a pubertal structure.
- **BMYO suppression is POWERED** (genotype main **p=0.020**, matching script 19's basal p=0.016).
- The WT substrate matures by **luminal (LASP/LHS) decline with BMYO flat** — a relative
  rebalancing by luminal loss. **Do not claim a basal expansion** (literature agrees).

**Holds — gene-level, not per-sample:** Myc bends **OFF** the WT developmental axis — convergence
rho(Myc, WT) **−0.08 at 6W (orthogonal) → −0.44 at 12W (oppositional)**. This is a correlation of
**LFC vectors**, not a sample-level coupling, so the ceiling does not apply.

**Directional only:** the LHS flip (interaction p=0.089) — credible via cross-modality (GRAY
directional nets, CHUNG ATAC), not via the p.

**CAVEAT to add — Issue #2's "endogenous Myc GATES the pubertal programme in WT".** This is a
**within-WT per-sample correlation pooled across both timepoints**: `myc_sig ~ proliferation`
r=**+0.93** — but the ambient in that window is **0.75**, putting it at the **92nd percentile,
perm p=0.08 — marginal**. Disaggregated it is n=6/cell (+1.00 at 6W, ambient 0.83, p=0.27;
+0.94 at 12W, ambient 0.71, p=0.03). **Report it as co-variation with its ceiling stated, not as
r=0.93.** The powered half — the transgene amplification — is untouched.

**Figures:** `myc_endogenous_amplification/endogenous_vs_transgene.pdf` (**headline**),
`program_panels.pdf`, `myc_axis_scatter.pdf`; `dev_program_myc_integration/state_stats.pdf`
(**headline** — BMYO powered, LHS directional), `convergence_vector_myc6.pdf` (the off-axis bend),
`heatmap_gsva_effect.pdf`, `e_selection_bound.pdf` (the dropout bound).
**Use `program_correlation.pdf` with care:** it is the r=0.93 family — caption against ambient.

---

## The ledger, in one screen

| claim | status | why |
|---|---|---|
| Myc raises mito **content +21–27%** | **STANDS** | genotype contrast, adjustment-robust, 2 methods |
| Myc raises **absolute nuclear OXPHOS**, reallocates the compartment | **STANDS** | d=1.27 / 1.15 |
| The **attenuation** (66% fade + 34% convergence) | **STANDS** | genotype×time contrast; SE-protected |
| Myc **amplifies** the pubertal programme; **BMYO suppression** | **STANDS** | genotype main, d 2.3–3.0 |
| Myc bends **off** the WT developmental axis | **STANDS** | LFC-vector correlation |
| **Biogenesis is a MYC-dose bystander** | **STANDS (strengthened)** | partial below its own ceiling |
| "MYC core selectively **buffered**" | **REFINE** | expression effect (z=0.29 vs matched) |
| "**OXPHOS is central** to phenotype" | **DROP as a coupling** | script 36: axis IS the global factor (r²=0.94) → no separable coupling is testable |
| **redox → MB2 fork** (Myc-independent) | **LEAD** | script 36: robust global-adjusted −0.63 across 3 quantifiers, axis r²=0.03 |
| "endogenous Myc **gates** the WT programme" (r=0.93) | **SOFTEN** | 92nd pct of its ceiling, p=0.08 |
| The **mitonuclear imbalance** (time) | **UNSUPPORTED** | no contrast fitted; all p>0.24 |
| The **mito→death bridge** | **NOT ESTABLISHED** | circular + fails every null |
| **Bbc3** attenuation | **NOT EVIDENCE** | = the chance expectation (1 of 23) |
| The **death phenotype** (kills at 6W not 12W) | **STANDS** | external IHC/inducible — never this RNA-seq |

---

## Appendix — gene sets used, by question

All from the library snapshot `data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt`
(986 sets; **884** GSVA-scored = the null universe) unless noted. **Do not rebuild these.**

| set family | source | size | used for | note |
|---|---|---|---|---|
| `MITOCARTA_*` (19 OXPHOS-arm, 21 biogenesis-arm sets) | MitoCarta 3.0 | 120 sets | Q2, Q4 | the mito axes |
| **`MITOCARTA_APOPTOSIS_PRO` / `_ANTI`** | MitoCarta 3.0 | **25 / 9** | Q3, Q4 | **MITO-DEFINED — the circularity. `pro_comp`/`priming` are built from these** (`08:285-290`) |
| `MITOCARTA_NUCLEAR_ENCODED` / `_MTDNA_ENCODED` | MitoCarta 3.0 | — | Q2 | content shares; the mito background for script 34 F3 |
| `GS_METAB_*` / `METAB_*` (12 axes) | Reactome/KEGG/Hallmark/GS_metabolic | 93 sets | Q4 | **NOT MitoCarta** — glycolysis, tca, oxphos_met, ppp, nucleotide, one_carbon, glutamine, fao, fa_synth, amino_acid, **cholesterol**, **redox** |
| `MYC_*`, `felsher` (67-gene, mito-stripped → 61), `hallmark_v2` | MSigDB + Felsher | 17 sets | Q4, Q5 | the MYC-dose proxies for the partial |
| `PROLIF_*` | library | 14 sets | Q4, Q5 | the `prolif` outcome |
| `MG_*` developmental + `data/dev_mec_annotation.csv` | Gray/Kessenbrock/Khaled 2025 | 179 sets (26 BMYO / 29 LASP / 29 LHS / 95 other) | Q5 | consensus MEC nomenclature |
| `MG_TEB_VS_DUCTAL_*_UP/_DN` | Gray | — | Q4, Q5 | `teb_dediff` = UP − DN |
| MB1/MB2/MB3 fork sets (`MB*_UF`/`_LF`) | Menegollo/Bentham 2024 | — | Q4 | `mb2_fork` = GSVA(UF) − GSVA(LF) |
| Tang 2024 15-RCD modalities | `data/cell-death/*.csv` | 15 sets | Q3 | branch-2 fGSEA |
| `cell_death_genes_consolidated` | curated | 512 pro-death (**41 mito / 471 non-mito**), 587 pro-survival | Q3 | branch-1 binomial; the `is_mitochondrial` split breaks the circularity |
| `MASS_MARKERS` / `_NOCHAP` | **hand roster (script 32)** | 14 / 7 | Q2 | **the corpus's only hand-built set** — a flagged exception; Hspa9/Hspd1 are direct MYC targets so the claim-bearing composite excludes them |
| IEG/dissociation panel | van den Brink 2017 | 13 | §0 | Fos/Fosb/Jun/Junb/Jund/Egr1/Ier2/Ier3/Atf3/Hspa1a/Hspa1b/Socs3/Zfp36 — **inferred covariate, not measured** |

**Method routing (CLAUDE.md):** fGSEA on Wald-stat ranks (unshrunken); **mitoPPS** on
linear-scale normalised counts (ratio-based, content-**blind**); **GSVA** on log-scale VST
(cohort-relative). ΔLFC methods use **raw** LFCs.

---

## What would settle the open items

| open item | experiment | also settles |
|---|---|---|
| Is there a mitonuclear imbalance at all? | **mtDNA qPCR** (copy number) | whether mt% means anything; the Q2 time story |
| Is mitochondrial content really up? | **TOMM20 / VDAC / CS blot** (or EM) | tension A; the lower-bound rider |
| Is the 6W substrate death-permissive? | **BH3 profiling** | the death mechanism, independent of survivor bias |
| Is the priming loss compositional? | **single-cell / deconvolution** | Issue #6's B6 pointer |
| Does redox capacity gate the tumorigenic fork? | **redox/glutathione assay** (GSH:GSSG, ROS) on the fork | the one Myc-independent coupling lead (script 36) |
| Is the IEG axis technical? | a recorded **dissociation-batch / viability** variable | whether the time axis is usable |

**Two bench experiments — mtDNA qPCR and the blot — close four items between them.**

*Numbers: `results/{mito_content_proxies, mtdna_axis_and_coupling_null, death_priming_reassessment,
ambient_corrected_couplings}.rds`. Verdicts print on sourcing; sandboxes walk every claim.
Scripts 00–31 unmodified; no re-fits. Figure scripts remain 36+ and gated on this narrative.*
