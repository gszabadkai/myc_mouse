---
name: mec-consensus-nomenclature
description: The Gray/Kessenbrock/Khaled 2025 consensus 3-type mammary epithelial cell (MEC) nomenclature adopted project-wide for developmental state
metadata: 
  node_type: memory
  type: reference
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
---

**Adopted 2026-07-08 for the myc_mouse developmental analysis.** Source:
Gray, Kessenbrock, Khaled et al., "Defining breast epithelial cell types in the
single-cell era", Developmental Cell 60, Sept 2025 (docs/MEC_types.pdf). A
consensus perspective harmonising ~13 prior nomenclatures onto THREE
phenotypically discontinuous resting-gland epithelial types (plus the lactocyte,
lactation-only, not relevant to the resting gland):

- **BMYO** = basal-myoepithelial (CL:0002324). Subsumes old basal / myoepithelial /
  MaSC-enriched basal/stem. The paper argues a discrete adult MaSC population is
  NOT supported; stem potential is a facultative property residing in BMYO.
- **LASP** = luminal adaptive secretory precursor (CL:4033057). Subsumes old LP /
  alveolar (AV) / luminal progenitor / secretory luminal.
- **LHS** = luminal hormone-sensing (CL:4033058). Subsumes old ML / mature luminal /
  HR+ / hormone-responsive. "Sensing" (not "responding") is the defended term.

Key stance: these are DISCONTINUOUS states, NOT a differentiation hierarchy /
trajectory (they fail lineage-trajectory assumptions). Do not frame them as an
ordered continuum (stem->differentiated).

**How applied in script 26 (dev_program_myc_integration):** the fine catalog tag
(MASC/BASAL_PROG/BASAL/LP/ALV/HR/other) is kept as `state_fine` for QC/overrides;
`state` collapses it to BMYO/LASP/LHS/other and drives the analysis. Collapse map:
{MASC,BASAL,BASAL_PROG}->BMYO; {LP,ALV}->LASP; {HR}->LHS. On the 179 MG_* sets the
consensus split is 34 BMYO / 35 LASP / 34 LHS / 76 other (other = signalling sigs,
mNEC contamination, unassigned clusters). Table 1 of the paper gives the per-source
crosswalk (e.g. Gray 2023 BA/AV/HS = BMYO/LASP/LHS; Henry, Pal, etc.).

See [[block-a-day1-gates]] for the Block A revision context (Issue #1).
