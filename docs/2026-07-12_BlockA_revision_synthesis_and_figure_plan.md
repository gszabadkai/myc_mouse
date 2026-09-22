# Block A revision synthesis (Issues #1-6) + a 2-supplementary-figure draft plan

Written 2026-07-12, branch `BlockA-revision-step-by-step` (HEAD 4874ca8). Purpose: (1)
synthesize how the step-by-step revision (Issues #1-6) changed the way the project should be
written up, and (2) on that basis specify TWO supplementary figures that tell the revised
story, mapped to EXISTING exploratory PDFs for a review draft. This is a DRAFT-assembly plan;
publication-quality figure scripts come AFTER the author reviews this. Per-issue detail lives in
`docs/2026-07-08_BlockA_revision_plan.md`; the pre-revision synthesis is
`docs/2026-07-07_block_A_synthesis.md` (now partly SUPERSEDED -- see below).

---

## Part 1 -- What the revision changed (synthesis)

### The old line (pre-revision synthesis, 2026-07-07)
"Myc drives a **primarily mitochondrial** program that licenses progression; the effect shows
progressive divergence; A-lead (mitochondria central) + D-bridge (death-timing)." Several
pieces of this did not survive the scrutiny of Issues #1-6.

### The revised line (what Issues #1-6 established)

1. **Constant driver, amplifying a pubertal program (Issues #1-2).** Myc dose is stable
   6W->12W; endogenous Myc GATES the pubertal proliferative/TEB program in wild-type tissue, and
   the transgene AMPLIFIES it (powered; d 1.0-3.0), specifically at 6W (the pubertal, TEB-rich
   window). The wild-type SUBSTRATE matures 6W->12W (luminal LASP/LHS decline, BMYO flat), and
   Myc bends OFF that developmental axis (orthogonal at 6W -> oppositional at 12W): BMYO
   suppression is powered, the LHS flip directional. => Frame as an oncogene amplifying a
   normal developmental program, not creating a de novo one.

2. **"Primarily mitochondrial" is an overclaim -> dissociable arms (Issue #3).** The mito
   "#1 rank" was a set-SIZE artifact of the permutation-null z. Size-fair, the OXPHOS/biogenesis
   core is TOP-TIER but CO-EQUAL with the MYC-target core and central metabolism (amino-acid,
   TCA, nucleotide) -- not uniquely biggest. What is CENTRAL to phenotype (partial-correlation,
   deconfounded against MYC-signature dose, present within WT) is a specific **OXPHOS + nucleotide
   + TCA** core; mito **biogenesis/translation** is a MYC-dose BYSTANDER. => Replace "primarily
   mitochondrial" with "a dissociable bioenergetic/biosynthetic core, separable from
   MYC-dose-driven biogenesis." (The parked lead-A rewrite.)

3. **Myc RAISES absolute OXPHOS; the "mito reduction" is REPRIORITISATION (Issue #4).** In
   absolute mRNA Myc elevates nuclear OXPHOS subunits (d~1.27, normalization-robust). The
   mitoPPS "OXPHOS down" is a within-compartment reprioritisation ratio, not an absolute drop.
   The mitonuclear shift is real: nuclear-led assembly at 6W -> mtDNA-completed running
   metabolism at 12W. => Do not say Myc "reduces OXPHOS"; say it REALLOCATES the mitochondrial
   compartment.

4. **The attenuation is magnitude, not rank, and no MYC-independent axis absorbs it (Issues
   #4-5).** fGSEA NES is flat (rank-based, magnitude-blind); what shrinks is |LFC| magnitude.
   Conditioning on developmental / ER-PGC1a-TF axes does NOT absorb the genotype x time
   interaction (they co-drive, not antagonise) -- level-association is not moderation.

5. **A citable two-mechanism account of the attenuation (Issue #6).** Via the exact identity
   `attenuation = WT-convergence + Myc-fade`: ~66% **Myc-fade** (oncogenic retreat) + ~34%
   **WT-convergence**, the convergence confined to the **biosynthetic arm** (WT matures onto the
   same axis). OXPHOS and the MYC-target core show NO convergence (close by fade alone), and the
   MYC identity core is the LEAST attenuated -- selectively BUFFERED. The fade itself is not
   monolithic: for the fade-led programs it decomposes into a SHARED developmental decline (also
   seen in WT) plus a MYC-SPECIFIC excess, and only the excess constitutes the attenuation (a
   shared decline cancels in the gap; OXPHOS ~72-90% Myc-specific -- see the fade-narration note
   under Supp Fig B). The fade co-occurs with a collapse of proliferation/mito-TF activity, and
   the Myc-specific excess reads as DE-AMPLIFICATION of the Myc-built proliferative/high-OXPHOS
   compartment -> a COMPOSITIONAL hint (dilution of the shrinking
   proliferative/TEB compartment), to be settled by deconvolution/single-cell (deferred;
   `docs/deconvolution_subproject_plan.md`).

### One-paragraph revised abstract-level claim (for the writeup)
*A constitutively expressed Myc amplifies the normal pubertal proliferative/biosynthetic
program of the mammary gland rather than creating a new one. Its transcriptional footprint is
broad -- the mitochondrial OXPHOS/biosynthetic core is engaged co-equally with the MYC-target
and central-metabolic programs -- and within the mitochondrial compartment Myc REALLOCATES
resources (raising absolute nuclear OXPHOS while reprioritising toward a nuclear-led->
mtDNA-completed maturation) rather than simply upregulating it. A dissociable OXPHOS/nucleotide/
TCA core, separable from MYC-dose-driven biogenesis, tracks the oncogenic phenotype and is
detectable even within wild-type tissue. The oncogene's effect attenuates from puberty to
adulthood not by switching programs but by two quantifiable mechanisms -- a two-thirds fade of
the oncogenic program on the maturing (and shrinking) proliferative substrate, and a one-third
convergence of the wild-type gland onto the same biosynthetic axis -- with the MYC identity core
selectively buffered.*

### Status of prior claims
- SUPERSEDED: "primarily mitochondrial", "progressive divergence", "Myc reduces OXPHOS".
- RETAINED (from earlier Block A, still valid): AP7 MB-fork biogenesis switch (script 18,
  powered); death-timing substrate-gating (scripts 23/25); the 4-method significance floor.
- These two supplementary figures carry the REVISED mechanistic story; the death-timing/MB-fork
  results remain their own (earlier-planned) figures.

---

## Part 2 -- Two supplementary figures (draft assembly from existing PDFs)

Design logic: **Supp Fig A = what Myc does (mode of action + mitochondrial engagement);
Supp Fig B = how the effect attenuates with age.** Each panel below names an EXISTING PDF under
`outputs/` to drop into the review draft, the claim it supports, and whether it is core or
optional. All exploratory PDFs are dense; final scripts will trim/re-render -- for the draft,
use as-is.

### Supp Figure A -- "Constant Myc amplifies a pubertal program; broad, reprioritising mitochondrial engagement"

| Panel | Source PDF (`outputs/`) | Claim it carries | Priority |
|---|---|---|---|
| A1 | `myc_endogenous_amplification/endogenous_vs_transgene.pdf` | Endogenous Myc gates the pubertal program in WT; the transgene amplifies it (powered), concentrated at 6W (Issue #2) | core |
| A2 | `dev_program_myc_integration/state_stats.pdf` | Myc modifies the developmental substrate: BMYO suppression powered, LHS flip directional (Issue #1) | core |
| A3 | `dev_program_myc_integration/convergence_vector_myc6.pdf` | Myc bends OFF the wild-type developmental axis (orthogonal @6W -> oppositional @12W) | optional (pairs with A2) |
| A4 | `myc_mito_centrality/q1_preferentiality_reframed.pdf` | Size-fair: OXPHOS core co-equal with the MYC-target core + central metabolism -- NOT uniquely "primarily mitochondrial" (Issue #3 Q1) | core |
| A5 | `myc_mito_centrality/q2_partial_correlation.pdf` | Dissociable arms: OXPHOS/nucleotide/TCA central to phenotype (orthogonal to MYC dose); biogenesis a bystander (Issue #3 Q2) | core |
| A6 | `attenuation_decomposition/C_mitopps_vs_absolute.pdf` | Myc RAISES absolute nuclear OXPHOS; mitoPPS "down" = reprioritisation; mitonuclear shift 6W->12W (Issue #4) | core |

Caption thesis: *Myc acts as a constant amplifier of the normal pubertal proliferative/
biosynthetic program, engaging the mitochondrial OXPHOS/biosynthetic core co-equally with the
MYC-target and central-metabolic programs, and reallocating -- not simply upregulating -- the
mitochondrial compartment.*

### Supp Figure B -- "The oncogenic effect attenuates by two mechanisms (fade + convergence), not a program switch"

| Panel | Source PDF (`outputs/`) | Claim it carries | Priority |
|---|---|---|---|
| B1 | `attenuation_decomposition/A_nes_paradox_two_rulers.pdf` | The attenuation is MAGNITUDE not rank: NES flat while |LFC| shrinks (Issue #4) | core |
| B2 | `attenuation_mechanism/A1_convergence_fade_by_program.pdf` | The exact-identity decomposition: ~66% Myc-fade + ~34% WT-convergence, per program (Issue #6) | core (headline) |
| B3 | `attenuation_mechanism/A2_conv_fade_gene_map.pdf` | Biosynthetic arm converges; OXPHOS + MYC-target core fade (WT diverges) | core |
| B4 | `attenuation_mechanism/A3_conv_fade_shares.pdf` | Convergence share by program, robust across the divergent + selection-independent universes | optional (robustness) |
| B5 | `attenuation_mechanism/B1_tf_bint_before_after.pdf` | No broadened TF axis absorbs the attenuation (Issues #5-6); pair with `B2_tf_delta_bint_bootCI.pdf` for the CI | core |
| B6 | `attenuation_mechanism/B3_tf_temporal_vs_oxphos_fade.pdf` | Compositional pointer: E2F/proliferation + mito-TF activity co-decline with the OXPHOS fade -> motivates deconvolution | core |

Caption thesis: *From puberty to adulthood the Myc effect does not switch programs; it
attenuates in magnitude via a two-thirds fade of the oncogenic program on a maturing, shrinking
proliferative substrate and a one-third convergence of wild-type tissue onto the same
biosynthetic axis, with the MYC identity core selectively buffered and no measured
MYC-independent axis accounting for it.*

### Narration for the fade-led panels (B2/B3): shared-developmental vs Myc-specific

A reviewer/reader can reasonably ask a second question about "fade": if these pathways ALSO
decline in wild-type with age, is the Myc+ decline just the tissue's normal developmental
down-trend? The answer, and the guardrail it imposes on the caption, uses only values already
in `results/attenuation_mechanism.rds` (NO new analysis) via the identity

  `Myc+ temporal slope  =  shared-developmental (= WT slope)  +  Myc-specific (= the attenuation)`

A decline SHARED by both genotypes CANCELS in the genotype gap, so it cannot create the
attenuation; the attenuation is, by construction, the Myc-SPECIFIC excess (Myc+ falls beyond the
shared developmental trend). This split is meaningful only where BOTH arms decline (same sign);
the biosynthetic/convergence programs are the opposite-sign case where WT RISES to meet Myc.

Fade-led programs (both arms decline), from the divergent set:

| program | Myc+ slope | shared-developmental (WT) | Myc-specific (= attenuation) | shared % |
|---|---|---|---|---|
| OXPHOS subunits | -0.39 | -0.11 | **-0.28** | 28% |
| OXPHOS (core) | -0.34 | -0.03 | **-0.30** | 10% |
| MYC-target V2 | -0.31 | -0.12 | **-0.19** | 38% |
| MYC Felsher | -0.27 | -0.11 | **-0.17** | 39% |
| proliferation (pooled) | -0.32 | 0.00 | **-0.31** | ~0% |

Three reads this licenses:
1. **The fade is mostly Myc-specific** (OXPHOS 72-90%, proliferation ~100%): the WT decline
   contextualizes but does not account for the gap-shrink. So the caption should say *"Myc+
   rides the shared developmental decline and falls FURTHER; the further part IS the
   attenuation"* -- NOT "Myc's effect exhausts."
2. **Proliferation is ~100% Myc-specific** (WT proliferation is flat at the gene level) -- the
   purest de-amplification; consistent with Myc having built the proliferative compartment at
   6W that then collapses.
3. **The MYC-target core has the HIGHEST shared fraction (~38%)** and the smallest Myc-specific
   excess -- its shared component is the endogenous pubertal Myc/proliferation program winding
   down in wild-type with age (ties directly to Issue #2). This is the flip side of "MYC core
   is buffered": less of its (already small) decline is Myc-specific.

**Why this ADDS something (not just a relabel):** it (a) prevents over-reading "fade" as
Myc-specific exhaustion, and (b) reframes the Myc-specific excess as DE-AMPLIFICATION of the
Myc-built proliferative/high-OXPHOS compartment -- which is precisely the compositional-dilution
hypothesis that panel B6 hints at and the deconvolution sub-project would test. Implementation:
a caption note + optionally a small stacked-bar re-render of B2 (shared vs Myc-specific per
fade-led program) at final-figure time; values are already computed, no re-analysis.

### Notes for the draft -> final transition
- **Overlap to resolve at final:** B5 can also draw on `attenuation_moderation/B_bint_before_after.pdf`
  (Issue #5 dev/TF panel); pick the broadened Issue #6 version (`B1_tf_bint_before_after.pdf`)
  as primary and cite #5 in the caption.
- **Dense panels to trim:** the `dev_program_myc_integration/trajectories_*` family (25 PDFs) is
  exploratory -- do NOT put in a figure; A2/A3 are the distilled panels. `q2_coupling_heatmap.pdf`
  is an alternative to A5 if a heatmap reads better than the partial-correlation dotplot.
- **Absolute-OXPHOS alternatives for A6:** `attenuation_decomposition/B_absolute_composite_trajectory.pdf`
  (3-normalization robustness) or `B_nuclear_oxphos_heatmap_vst.pdf` (gene-level) if the
  reconciliation panel needs support.
- **What these two figures deliberately EXCLUDE:** the MB-fork biogenesis switch (script 18) and
  death-timing (scripts 23/25) -- those belong to their own (earlier-planned, retained) figures,
  not the revision's mechanistic-correction story.

---

## Part 3 -- Next steps (after author review of this plan)
1. Author reviews the synthesis + the panel selection; adjusts the A/B partition and
   core/optional calls.
2. Lock the lead-A rewrite wording ("dissociable arms" replacing "primarily mitochondrial") --
   the parked framing decision this synthesis now forces.
3. THEN write the publication figure scripts (a new `paper-figures`-style branch/scripts) that
   re-render the selected panels to a consistent style -- the deferred Block B work, now with a
   revised, defensible narrative to serve.
4. Deconvolution sub-project (`docs/deconvolution_subproject_plan.md`) remains parked; it would
   only add a panel to Supp Fig B (the compositional test) if activated.
