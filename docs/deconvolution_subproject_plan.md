# Deconvolution sub-project plan (reference; NOT current priority)

Standalone scoping doc for a future cell-composition deconvolution of the MMTV-Myc bulk
timecourse. Written 2026-07-12 as the basis for later research, not for immediate execution.
It exists because Issue #6 (`scripts/31_attenuation_mechanism.R`) left ONE mechanistic
ambiguity that bulk RNA-seq cannot resolve on its own.

---

## 1. Why this sub-project exists (the question deconvolution would answer)

Issue #6 decomposed the 6W->12W attenuation of the Myc effect into an exact identity:
`attenuation = WT-convergence + Myc-fade`, with the attenuation ~66% **Myc-fade** (the oncogenic
program retreats on the aging substrate) and ~34% **WT-convergence** (biosynthetic arm only).
The dominant Myc-fade term is fundamentally ambiguous at bulk resolution:

- **(a) per-cell:** Myc's transcriptional program genuinely weakens in each cell, OR
- **(b) compositional:** per-cell Myc activity is unchanged but the Myc-responsive
  proliferative / terminal-end-bud (TEB) / luminal-progenitor compartment SHRINKS as a
  fraction of the bulk as the pubertal gland matures to a quiescent adult one -> the
  bulk-averaged effect DILUTES.

Bulk expression is `signal = sum_c (fraction_c * per-cell-expression_c)`; fraction and per-cell
expression are confounded in one sum. Issue #6 found an indirect pointer to (b): the
E2F/proliferation and mito-TF activities co-decline with the OXPHOS fade (script 31 Part B3),
and proliferation attenuates by pure fade. Deconvolution is the bulk-feasible way to **bound**
the compositional contribution; single-cell/snRNA-seq is what would **settle** it.

**Precise deliverable of the sub-project:** per-sample estimates of cell-state fractions
(basal/myoepithelial, luminal-progenitor, mature/hormone-sensing luminal, and a
proliferative/TEB state), then (i) test whether the proliferative/LP fraction drops 6W->12W,
and (ii) test whether conditioning the Myc-fade / the genotype x time interaction on that
fraction shrinks it (the compositional-absorption analogue of Issues #5/#6, but on a FRACTION
rather than a program score).

---

## 2. Hard constraint: nothing is on disk

Verified 2026-07-11 (see Issue #6 log): NO single-cell reference, signature matrix, `.h5ad`/
`.loom`/SCE, and NO deconvolution package are present in `data/`, `external/`, or
`00_setup_packages.R`. So this sub-project requires BOTH an outward-facing data fetch (needs
author approval) AND a package install (author, per Option A). That is why it is scoped
separately and not folded into an Issue.

---

## 3. Candidate single-cell references (critical evaluation)

The reference must be **mouse mammary gland**, ideally spanning the pubertal->adult window that
brackets our 6W and 12W timepoints, and annotated to states compatible with the project's
adopted Gray/Kessenbrock/Khaled 2025 consensus nomenclature (BMYO / LASP / LHS; see memory
`mec-consensus-nomenclature`). Verify every accession before download.

| Reference | Scope | Fit for us | Risks / caveats |
|---|---|---|---|
| **Bach et al. 2017** (Nat Commun; GEO ~GSE106273) | Mouse mammary across development (nulliparous, gestation, lactation, involution), 10x | **Primary candidate.** Mouse, developmental, clean basal / LP / mature-luminal + proliferative annotation; widely used as a deconvolution reference | Nulliparous adult-weighted; strain/age not identical to ours; pubertal TEB state under-represented vs Giraddi |
| **Giraddi et al. 2018** (Cell Rep; ~GSE111113) | Fetal -> postnatal/pubertal mouse mammary trajectory | **Strong complement** for the 6W pubertal / TEB state that Bach under-covers (directly relevant to the Issue #2 TEB amplification) | Early-development-weighted; fewer adult quiescent cells |
| **Pal et al. 2021** (EMBO J; mouse+human developmental atlas) | Mouse & human mammary across development/pregnancy | Broad, consensus-aligned; enables a human cross-check | Verify the exact mouse accession; large, needs subsetting to relevant stages |
| **Saeki et al. 2021** / other normal-gland atlases | Mouse mammary across age/estrous | Adult quiescent states, estrous covariate | Coverage of proliferative compartment varies |
| Human breast atlases (Wu 2021, Kumar/Nguyen, Reed 2024) | Human normal/tumor | Only as ortholog cross-validation | Human; not a valid mouse deconvolution reference |

**Recommendation:** a **two-reference strategy** -- Bach 2017 as the primary adult/developmental
reference, Giraddi 2018 to cover the pubertal/TEB compartment -- combined via an
ensemble/robustness method (SCDC) so results are not hostage to one atlas. Cross-validate
fractions against canonical markers regardless of method: *Krt5/Krt14/Acta2* (BMYO),
*Krt8/Krt18* (luminal), *Elf5/Kit/Aldh1a3* (LASP/LP), *Esr1/Pgr/Prlr* (LHS/HS), *Mki67/Top2a*
(proliferative/TEB).

---

## 4. Candidate deconvolution methods (critical evaluation)

All are **reference-based** (bulk deconvolved against a labeled sc reference) unless noted.
Marker-only "pseudo-deconvolution" is listed last as the zero-dependency fallback.

| Method (pkg) | Approach | Strengths | Weaknesses / fit |
|---|---|---|---|
| **MuSiC** (Wang 2019; R/Bioc) | Cross-subject weighted NNLS; up-weights genes consistent across sc subjects | Handles multi-subject references, no manual signature; solid-tissue standard | Needs a multi-subject sc reference; sensitive to reference cell-type balance |
| **BisqueRNA** (Jew 2020; R) | Reference-based regression with assay-platform decomposition | Robust to bulk-vs-sc platform bias (relevant: our bulk is RSEM, reference is 10x); fast | Assumes overlapping subjects help (we have none); marker-mode simpler but coarser |
| **SCDC** (Dong 2021; R) | ENSEMBLE across multiple references, ENSEMBLE weighting | Directly supports the Bach+Giraddi two-reference plan; robustness | More moving parts; heavier |
| **DWLS** (Tsoucas 2019; R) | Dampened weighted least squares | Accurate for rare/small compartments (the shrinking proliferative fraction is exactly this) | Slower; signature construction matters |
| **BayesPrism** (Chu 2022; R) | Bayesian; jointly infers fractions + per-cell-type expression, separates malignant from reference | **Notable for us:** MMTV-Myc tissue is hyperplastic/transformed -- BayesPrism tolerates a bulk state absent from a normal reference and returns per-state expression (could directly probe per-cell fade) | Compute-heavy; needs care with priors |
| **CIBERSORTx** (Newman 2019) | SVR on a signature matrix + batch correction | Strong, batch-robust | Web/Docker, licensing friction; not a clean scripted R dep |
| **granulator** (Bioc) | Wraps/benchmarks several methods | Good for the BENCHMARKING step -- run several, compare | A harness, not a new method |
| **Scaden** (Menden 2020; Python) | Deep-learning, simulation-trained | Robust, reference-flexible | Python dependency; breaks the R/Positron Option-A flow |
| **Marker-score proxy** (in-repo GSVA composites) | ssGSEA/GSVA of state marker sets as fraction PROXIES | Zero external deps; feasible today | NOT true deconvolution -- conflates fraction with per-cell intensity; a first-pass bound only |

**Recommendation / decision tree:**
1. **Zero-cost first pass (in-repo, now-able):** the marker-score proxy using the library's
   BMYO/LASP/LHS + proliferation/TEB sets -- already what Issue #6 Part B did partially. Honest
   bound, no download. Use it to decide whether the full sub-project is warranted.
2. **Primary reference-based run:** **MuSiC** or **BisqueRNA** with **Bach 2017**; add
   **Giraddi 2018** and ensemble via **SCDC**. Because our bulk is RSEM and the reference is
   10x, Bisque's platform decomposition is attractive; MuSiC is the conventional default.
3. **Robustness/benchmark:** run the panel through **granulator**; report fraction concordance
   across methods and against the marker cross-check. If per-cell expression is needed (to
   attack the per-cell-vs-composition question directly), add **BayesPrism**.
4. **Validate** every run against the canonical markers (Section 3) before interpreting.

---

## 5. Analysis design once fractions exist

1. **Descriptive:** per-sample fractions across the 4 groups; test the proliferative/LP/TEB
   fraction 6W->12W (does it drop, and more in Myc+?). Powered contrast if the effect is large.
2. **Compositional absorption:** refit the Issue #6 / Issue #5 moderation with the proliferative
   fraction as the covariate -- `oxphos ~ timepoint*myc_status + frac + timepoint:frac` -- and
   report `delta_b_int` + bootstrap CI (reuse `scripts/30`/`31` machinery). Does the fraction
   absorb the Myc-fade where TF activity did not?
3. **Reconcile with Issue #6:** map the fraction change onto the `myc_fade` term of
   `results/attenuation_mechanism.rds` -- how much of the fade is compositional vs residual
   per-cell.

---

## 6. Ceilings that deconvolution does NOT remove (state honestly)

- **Endogeneity persists:** the proliferative fraction is partly a Myc consequence (Issue #2:
  Myc drives the pubertal proliferative/TEB program), so conditioning on it still over-controls
  -- deconvolution sharpens the instrument (an explicit fraction) but does not make it causal.
- **Reference-bulk mismatch:** strain, age, platform (10x vs RSEM), and the transgenic/
  hyperplastic state of MMTV-Myc tissue differ from any normal atlas -> deconvolution BIAS;
  BayesPrism mitigates the "state absent from reference" problem but not all of it.
- **n = 6/group** caps the fraction x genotype x time interaction power (same floor as Issues
  #5/#6).
- **Fundamental identifiability:** even perfect fractions leave `fraction x per-cell` only
  partially separable without per-cell measurements. The DEFINITIVE test remains
  **single-cell/snRNA-seq of the 4 groups** or FACS-sorted compartments; deconvolution BOUNDS.

---

## 7. Concrete future scripts (when activated)

- `scripts/32_build_deconv_reference.R` -- fetch + QC + label the sc reference(s) (Bach 2017
  [+ Giraddi 2018]); build the SCE/ExpressionSet and a signature matrix; cache under
  `data/deconv_reference/` (gitignored) with a provenance README (accession, date, cell-type
  map to BMYO/LASP/LHS + proliferative). REQUIRES: author-approved download + package install
  (add MuSiC/BisqueRNA/SCDC[/BayesPrism] to `00_setup_packages.R`).
- `scripts/33_deconvolve_bulk.R` -- deconvolve the 24 bulk samples (primary + ensemble +
  benchmark via granulator); marker cross-validation; per-sample fraction table +
  `results/deconv_fractions.rds`.
- `scripts/34_composition_vs_attenuation.R` -- the Section-5 analyses; integrate with
  `results/attenuation_mechanism.rds`; report the compositional share of the Myc-fade with CIs
  and the full endogeneity/identifiability ceiling.

## 8. Open decisions for when this is picked up
- Single vs ensemble reference (recommend ensemble Bach+Giraddi).
- Method: MuSiC vs Bisque as primary; whether to add BayesPrism for per-cell expression.
- Python (Scaden) allowed, or R-only to preserve the Positron/Option-A flow (recommend R-only).
- Whether the in-repo marker-score first pass is sufficient to answer the paper's needs before
  committing to the full external-reference build.
