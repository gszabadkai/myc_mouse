---
name: figure-b-loaders-are-curated-not-top
description: "The pathway-loading figure B plots 11 CURATED composites, not the empirical top loaders; the dominant axis is mito-LED but not mito-SPECIFIC (non-mito growth programmes load 0.95-0.98), and the very-top mito monopoly is construction-inflated"
metadata:
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-20T11:45:02.856Z
---

**Found 2026-07-20** (read-only interrogation of `results/pathway_loading.rds`, script 37).
Prompted by the author's question: are the axes in `outputs/pathway_loading/B_axis_loadings.pdf`
the actual top loaders, and are there any NON-mito sets in the high rankings (a possible argument
for mitochondria-PRIMARY regulation of early tumorigenesis by the Myc-mito axis)?

**1. Figure B plots only the 11 CURATED composites**, NOT the empirical top loaders. The true top
per-set loaders on the dominant axis are the Gray `TFT_*_GRAY_*_LE_MITO` TF lanes at ~0.99 -- none
of which appear in figure B. So B is a designed panel, not a ranking.

**2. Mito enrichment is REAL but not the whole story.** Mito-defined sets are **44% of the library
(385/884)** but **100% of the top 100** and **195/200 of the top 200** -- a striking enrichment.
BUT non-mito growth programmes load nearly as high: MYC-target/METABRIC breast-cancer biclusters
**0.98**, pentose-phosphate **0.972**, nucleotide **0.965**, proliferation max **0.955** (median
0.853), MYC_signatures median **0.876** -- vs the redox control **0.26**. => the dominant axis is a
coordinated Myc anabolic-proliferative-mitochondrial GROWTH state: **mito-LED but not
mito-SPECIFIC.**

**3. The very-top mito monopoly is CONSTRUCTION-INFLATED.** The 218 `_LE_MITO` + biogenesis-
discrimination `_MITO` sets are subsets of MitoCarta BY BUILD (`docs/library_reference/
Gray_et_al_developmental_TFS_selection.md` sec 2c: TF program INTERSECT MitoCarta, promoted at BH
FDR<0.01 & >=15 mito genes). Proof: `TFT_MYC_GRAY_BA_LE_MITO` loads **0.993 (rank 6)** but the same
TF's non-mito lane `TFT_MYC_GRAY_AP_HE` loads **0.013 (rank 881)** -- loading tracks mito gene
CONTENT, not the regulator.

**4. The ranking does NOT stratify sub-programme** (oxphos vs mitoribo/biogenesis vs priming):
across the 218 `_MITO` lanes `cor(|loading|, oxphos_frac) = -0.005`, `cor(|loading|, mitoribo_frac)
= 0.148` -- all saturate at the mito ceiling (median |loading| 0.981). It reads "mito as a block,"
not sub-programme resolution.

**LIBRARY COVERAGE (script 37 PART C4 / figure `E_library_coverage.pdf`, 2026-07-20).** The 884
scored sets divide functionally: TF_targets 360 (218 construction _MITO + 142 non-mito),
Mammary_development 203 (all non-mito), MitoCarta 120, Metabolism 93, Apoptosis 36,
Biogenesis_apoptosis_intersections 25, MYC_signatures 17 (all non-mito), Biogenesis_discrimination 16,
Proliferation 14. Total 385 mito / 499 non-mito, but **259 of the 385 "mito" are construction _MITO
lanes -> genuine MitoCarta = 126.** KEY POINT for the loading argument: **MYC's OWN curated signatures
are NON-mitochondrial gene sets** -- the MYC_signatures category is 0/17 mito, and Felsher's
integrative MYC signature is only 5/61 genes mito (**8%**, at the ~7% transcriptome baseline = NO mito
enrichment). So the Myc-mito coupling is EMERGENT (co-variation of the dominant axis), NOT built into
the MYC sets. **PROVENANCE: only 19 of 884 scored sets are off-the-shelf MSigDB Hallmark**
(Mammary_development 9, Metabolism 4, Proliferation 4, Apoptosis 1, MYC_signatures 1; `_AND_`
intersections excluded) **-> ~98% of the library was CONSTRUCTED for this project** (marked by diamonds
in figure E1). Visual = stacked bar (not a pie), dataviz-skill palette validated (CVD dE 17.4).

**How to apply:** the honest, defensible statement is **"OXPHOS/mito is the single strongest-loading
programme of the dominant transcriptomic axis"** (a DESCRIPTION -- the loading LEVEL, near-tautological
since mito_oxphos IS the global factor at r=0.968) -- NOT "mitochondria PRIMARILY regulate early
tumorigenesis." Loading = covariation, NOT primacy; at n=24 the "OXPHOS central" question is
UNTESTABLE (no separable variance). The mito-vs-non-mito loading-enrichment panel (script 37 ADDENDUM
1, specced but not built as of 2026-07-20) makes this figure-worthy. See [[correlation-ceiling]],
[[mito-content-bounded]], [[dont-assert-a-mechanism-without-checking-it-can-work]].
