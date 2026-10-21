---
name: le-is-low-expressing
description: "Gray's HE/LE = HIGH- and LOW-EXPRESSING (lineage-program suppression), NOT high/low oestrogen -- the library metadata is wrong, and the LE state is the candidate carrier of the PGC1a/OXPHOS programme"
metadata: 
  node_type: memory
  type: reference
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-27T01:03:50.705Z
---

**THE LIBRARY METADATA IS WRONG AND IT MISLED A WHOLE SESSION.** Gray et al. 2023, Cell Reports
42:113293 (`docs/GRAY_2023_paper.pdf`), **p.7 verbatim**: *"the AP4, HS4, and BA2 clusters constituted
a 'low-expressing' (LE) state within their respective lineages compared to their high-expressing (HE)
counterparts."* `LE` = **LOW-EXPRESSING**, i.e. **lineage-program suppression** — the paper's entire
subject (*"systematic suppression of lineage programs ... as a common feature of mammary epithelial
expansion during puberty, the estrous cycle, and gestation"*). **Nothing to do with oestrogen.**

But `data/genesets_from_library/provenance_table.csv` glosses `MG_HEVSLE_*` as **"HE vs LE (AP)"**
tagged **"adult"**, and `docs/library_reference/mammary_dev_sets_catalog.md:108-110` expands that to
**"High- vs low-estrogen"**. Both wrong. The library reference is read-only, so the correction lives
in script 44 and in section 0.7a of `docs/2026-07-26_introduction_alignment_and_the_question.md`.
`AP` = alveolar progenitor = **LASP** in the consensus nomenclature ([[mec-consensus-nomenclature]]).

**WHY IT MATTERS — the author's hypothesis, and it survives testing.** ESRRA's and GABPA's
mitochondrial programmes are detected specifically in **`AP_LE`**, the lineage-suppressed progenitor
state the pubescent gland is rich in. So the live alternative to "PGC1a activity declines" is **"the
cell state that CARRIES the PGC1a-driven OXPHOS programme becomes less abundant"** — which also
explains why `CORE_MITO` is flat while its respiratory arm falls ([[oxphos-apoptosis-module-verdict]]).

**THE TEST (script 44 PART B2), mito-removed to kill the circularity — 25 of the 87 OXPHOS subunits
ARE AP LE-marker genes, so raw numbers are the OXPHOS result restated:**

| lineage | HE-marker (`_UP`) | LE-marker (`_DN`) | verdict |
|---|---|---|---|
| **AP** (the prediction) | -0.049 (21st pct) | **-0.238 (0th pct)** | **holds** |
| BA (control) | +0.056 (99.7th) | -0.105 (1.4th) | holds |
| HS (control) | -0.275 (0th) | -0.141 (0.15th) | **fails** — HE falls MORE |

So: **LE-specific decline in AP and BA, absent in HS**, where the whole lineage falls (`MG_HS_GRAY`
-0.142, TEB-vs-ductal -0.435). Honest statement is *loss of the lineage-suppressed state in the
alveolar-progenitor and basal compartments*, **not** an AP-exclusive carrier.

**BOUNDS.** Gray's own nuance: `AP4/HS4/BA2` are *"broadly distributed among samples of all ages"*, so
this is a shift in LE **ABUNDANCE**, not a state appearing/vanishing; the pubertally restricted LE
clusters are `AP5` (*"almost completely restricted to 3-week-old glands"*) and `HS5` (*"significantly
enriched in pubescent glands"*, tied by the authors to the TEB). The within-WT version
(`ox_beta_time` -0.459 -> -0.151 on a **disjoint** AP-LE score) is **descriptive only** — n=12,
`batch = timepoint`, cannot separate lineage shift from developmental regulation. Deconvolution is
the design that could (`docs/deconvolution_subproject_plan.md`).

**TWO THINGS WORTH CITING.** Gray argues **LE is NOT a proliferation readout** (*"the majority of
IdU+ cells were in fact HE cells"*) — independent convergence on our own de-respires-without-
de-proliferating result, from a different assay. And **MYC is elevated mostly in ER-negative tumours,
which arise from the LASP/AP lineage** (the author's point) — so the state carrying the PGC1a/OXPHOS
programme is the lineage of origin for the MYC-high human tumours. That is the link to paragraph 2
and the METABRIC arm.

**How to apply:** never read a Gray context label off our provenance table — check the paper. And
whenever a `MG_HEVSLE_*` or lineage set is compared with a mitochondrial one, report the
**mito-removed** version, or the comparison is partly self-referential.
