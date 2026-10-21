---
name: panels-follow-the-narrative
description: Build a figure panel only when the Results sentence that needs it is written; do not pre-build panels for sections still to come
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-29T14:23:38.911Z
---

**2026-07-29.** Asked which global-overview panels to build for the opening Results paragraph, the
author selected four, then on seeing the plan cut two of them: *"Leave Fig S1C and S1D for later,
they will come up in the narrative later. PCA shows the alignment along mito, DE can be used for the
attenuation question. Not now."*

**Why.** The author knows where each piece of evidence lands in the argument, and a panel built early
is a panel built without its sentence. The sample PCA is not an opening QC ornament here -- it is the
evidence for *alignment along the mitochondrial axis*, which belongs to the global-axis section. The
DE counts are not a "how big is the dataset" statistic -- they are the setup for *attenuation*
(myc_6W 2777 significant genes at FDR 10% versus myc_12W 239). Building either into the overview
would spend the evidence in the wrong place and force it to be re-argued later.

**How to apply.** Build a panel when the sentence it serves exists, not when the data for it exists.
When scoping ahead, list the candidates and let the author assign them to sections -- do not assume a
conventional figure order (overview -> QC -> PCA -> DE) applies. Record deferred candidates with the
section they belong to and any numbers already checked, so the work is not re-derived:
`figures/panels/PANELS.md` has a "Not built here" table for exactly this. Related:
[[figure-phase-architecture]], [[hypothesis-generating-not-failed-tests]].
