---
name: correlation-ceiling
description: "Every composite correlates with every other at ~0.8 — the ceiling is ONE global per-sample axis (cor 0.968), part Myc-contrast part within-group prep/composition. Coupling method of record is script 36's linear z-score + design->PC1->residual, NOT GSVA, NOT mean-subtraction."
metadata:
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-20T09:39:06.472Z
---

**The fact that reorganises the corpus** (`scripts/35_ambient_corrected_couplings.R`, 2026-07-18).
Median |rho| of an axis to an **ARBITRARY** one of the 884 scored library programmes:

| axis | all (n=24) | 6W | 12W |
|---|---|---|---|
| mito_oxphos | 0.80 | 0.82 | **0.78** |
| nucleotide | 0.77 | 0.82 | 0.71 |
| tca | 0.77 | 0.80 | 0.78 |
| mito_biogenesis | 0.74 | 0.78 | 0.71 |
| **redox** | **0.22** | **0.41** | **0.23** |

**A published coupling of 0.8 is what the window hands you.** Not a finding.

**WHERE IT COMES FROM (script 36, 2026-07-18 — the resolution, after the author's two method
notes `docs/GSVA_global_background_and_pathway_coupling.md` + `..._Gene_set_quantification...md`).**
It is **ONE global per-sample common-mode axis**: `cor(mito_oxphos, mean of all 884 scores)=0.968`;
regress that mean out and the ambient drops **0.80->0.13** (~chance). **NOT low n** (chance floor
0.13), **NOT curation** (random pairs 0.54), **NOT the genotype split alone** (0.83 within `6W_pos`,
contrast gone). GSVA doesn't centre within-sample, so any broad shift lifts every score.
**Methodological in FORM, biological+technical in ORIGIN.** The script-35 headline "MOSTLY MYC DOSE"
is **REFINED**: the RAW ceiling = (four-group contrast; raw PC1~group r=0.48) + (within-group axis).
The **residual** axis (design removed first) tracks **epithelial +0.79 / contamination -0.70**
(composition/prep) AND **proliferation +0.65 / MYC +0.44** — neither cleanly. `redox` is still the
control (r^2 with global mean 0.03) — the axes Myc does NOT drive are the only readable ones.

**THE METHOD OF RECORD IS NOW SCRIPT 36, NOT 35.** Primary quantifier = **linear mean gene-wise
z-score** (correlation IS avg cross-gene covariance); singscore = robustness; GSVA = continuity.
Correction = **design->PC1->residual, factor estimated AFTER design + leave-pair-out** — **NOT
mean-subtraction** (equal-loading + sum-to-zero artifact = the script-35-era instability). **THE
TRAP:** when an axis IS the global factor (mito_oxphos r^2=0.94), removing it leaves ~4% noise; its
global-adjusted -0.62 is **noise-on-noise, NOT antagonism** — read a global-adjusted coupling ONLY
where `resid_frac >= ~0.10`.

**TWO ERRORS I MADE AND THIS SCRIPT CAUGHT — both are the same shortcut:**
1. **"6W is stressed, 12W is clean" is WRONG.** Script 34's ambient collapsing 0.71 → 0.35 was
   measured on the mitonuclear **IMBALANCE** — a mitoPPS **RATIO**, and a ratio cancels the
   component every programme shares. Per axis the ceiling barely moves (0.82 → 0.78). **NEVER
   BORROW ONE AXIS'S AMBIENT FOR ANOTHER.** I built a whole plan's headline on the borrowed number
   after explicitly telling the author that borrowing was the error class behind the retractions.
2. **A RAW ambient is the wrong null for a PARTIAL claim.** Issue #3 Q2 never claimed "OXPHOS
   correlates with proliferation" (it does; so does everything) — it claimed the coupling
   **survives removing MYC dose**. The matched null is a **partial ceiling**: partial rho(axis,
   every set | MYC). Match the null to the CLAIM, not to the data type.

**Q4 (Issue #3 Q2) — HALF HOLDS:**
- **"OXPHOS is CENTRAL" FAILS.** Partial →prolif **+0.59 vs a partial ceiling of 0.45** (68th pct,
  p=0.32); →TEB **+0.44 vs 0.45** (p=0.51). OXPHOS retains residual co-variation with
  **everything**; proliferation is not special within it. "Central" implies PREFERENTIAL coupling.
- **"biogenesis is a BYSTANDER" HOLDS, harder than published.** Partial →prolif **+0.11 vs its own
  ceiling 0.22** (27th pct); →TEB **−0.03** (7th pct) — **BELOW** an arbitrary programme. Once MYC
  dose is removed it tracks **nothing**. This is the well-supported half; lead with it.
- **SCRIPT 36 RE-RANKS THE LEADS (linear correction supersedes 35's raw-ambient).** "OXPHOS central"
  is not just at the ceiling, it is **UNTESTABLE** — the mito axes ARE the global factor (r^2=0.94),
  so no separable coupling exists (cleanest statement of the failure). **`redox → MB2 fork` is the
  ONE robust testable lead**: global-adjusted **−0.63**, same sign −0.62…−0.73 across ALL 3
  quantifiers, axis r^2=0.03 (Myc-independent — lower redox/glutathione ↔ tumorigenic fork).
  **`cholesterol/mevalonate → fork` is DEMOTED to +0.14** ("separable but weak") — its script-35
  96th-pct/p=0.043 headline was a **raw-ambient artifact**. biogenesis-bystander unchanged.
- **1 of 96 tests beats BH<0.05**. **Nothing is significant** at n=24 — leads are for the BENCH.
- The OXPHOS-vs-biogenesis ordering was **never tested as a DIFFERENCE**; two percentiles are not a
  contrast (and under script 36 OXPHOS isn't even testable).

**Q5:** Issue #2's "endogenous Myc gates the WT programme" **r=0.93** sits at the **92nd percentile
of its own ceiling (0.75), perm p=0.08 — MARGINAL**; disaggregated it is n=6/cell. Issue #2's
**powered half (transgene amplification, genotype main effect, d=2.3–3.0) is untouched.**

**THE DIVIDING LINE THROUGH THE WHOLE CORPUS:** the prep axis is genotype-**independent** (p=0.49)
so it cancels inside every Myc-vs-WT comparison. **Contrasts survived; per-sample cross-sectional
couplings did not.** Content +21–27%, absolute nuclear OXPHOS, the attenuation, the amplification,
BMYO suppression = all contrasts = safe. The mito→death bridge, the imbalance-as-time-claim,
"OXPHOS central", Bbc3 = all couplings/single genes = gone.

**THE CONSTRUCTIVE NEXT STEP (script 37, 2026-07-19).** Instead of only *removing* the axis, DISSECT
it as a phenotype (author's doc `Low_dimensional_pathway_structure_and_global_axis_analysis.md`).
`scripts/37_pathway_loading_and_technical_resolution.R` (`results/pathway_loading.rds`): the matrix is
**low-dimensional, quantified** — PC1=76% var, effective dimensionality **1.7** of 23 (n=24). **THE
POSITIVE REFRAME that answers "OXPHOS in proliferation":** as a composite, **OXPHOS IS the dominant
axis** (loading r=0.968) — the loading LEVEL is a ROBUST *description* of the axis (near-tautological,
not a coupling); `redox` is a weak loader (42nd pct, the non-Myc control). The loading-CHANGE by
genotype (`S ~ G*genotype*timepoint`) is a CONTRAST (readable where the residual coupling was not) but
here only an **exploratory lead** (OXPHOS rel +71% but near-collinear at n=6/group, biogenesis none).
**Both-lenses tech correction** (design-only primary + contamination/IEG-removed 2nd lens, r=0.57):
NAMED over-correction — the RNA anchors absorb **r²=0.67 of proliferation variance** (IEG↔prolif
−0.75), so they remove biology too. See [[batch-equals-timepoint]]. Script 38 (mitoPPS priming-
coupling in ratio space + MYC-vs-PGC1a/ESRRA separability) is next.

**LOADING FIGURE, INTERROGATED (2026-07-20).** `outputs/pathway_loading/B_axis_loadings.pdf` plots 11
CURATED composites, NOT the empirical top loaders (Gray `TFT_*_LE_MITO` lanes ~0.99). Mito-defined =
44% of the 884-set library but 100% of top-100 — striking, BUT non-mito growth programmes load
0.95–0.98 (nucleotide 0.965, PPP 0.972, MYC-target/METABRIC biclusters 0.98) vs redox control 0.26,
and the very top is CONSTRUCTION-inflated (the 218 `_LE_MITO` sets are MitoCarta subsets by build;
`TFT_MYC_..._LE_MITO` 0.993 rank6 vs same TF's `TFT_MYC_..._AP_HE` 0.013 rank881). => axis is
**mito-LED not mito-SPECIFIC**; loading = covariation NOT primacy. See
[[figure-b-loaders-are-curated-not-top]]. mt% (the mtDNA/imbalance axis) is three-way-confounded and
death-agnostic — [[mtpct-is-three-way-confounded]].

**How to apply:** for coupling, use **`scripts/36_linear_pathway_coupling.R`** — linear z-score,
design->PC1->residual, leave-pair-out, `resid_frac` gate (don't interpret a global-adjusted coupling
where the axis IS the global factor). Report the global factor as a phenotype too, don't only remove
it (doc 1 §14). For fGSEA the ceiling analogue is script 21's permutation null. Full write-up:
`docs/2026-07-18_narrative_synthesis_five_questions.md` §0.5/§0.6 (origin, instrument, hypothesis
ledger, sentence-level fix pattern). **Script 36 was RUN by the author 2026-07-18 and reproduced the
scratchpad numbers EXACTLY** (`results/linear_pathway_coupling.rds` + 3 PDFs) — confirmed, not
provisional. See [[death-priming-change-untested]], [[sample-level-couplings-need-a-null]],
[[dont-assert-a-mechanism-without-checking-it-can-work]], [[mito-content-bounded]].
