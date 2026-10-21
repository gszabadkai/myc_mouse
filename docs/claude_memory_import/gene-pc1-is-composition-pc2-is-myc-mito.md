---
name: gene-pc1-is-composition-pc2-is-myc-mito
description: "Gene-level PC1 (44%) is the epithelial/immune COMPOSITION nuisance (genotype-independent, correctly dropped); the Myc-mito program is gene-PC2, which the pathway analysis foregrounds as the 76% global axis -- the manuscript entry point for mito involvement in the Myc program"
metadata:
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-20T10:56:08.365Z
---

**The gene-level entry point for "mitochondria in the Myc program"** (2026-07-20, script 37 PART A2,
`results/pathway_loading.rds`; Spearman n=24, reproduced read-only). The author asked what the pathway
compositing KEEPS vs DROPS: the pathway global axis is genotype-associated, but the largest *gene*-level
axis is a different direction -- what is it? Decomposed the gene-level PCA (top-500 VST) by covariate
correlation, correspondence to the pathway axis, and fGSEA on genome-wide loadings.

**gene-PC1 (44% of gene variance) is the COMPOSITION nuisance -- correctly DROPPED.**
- **genotype-INDEPENDENT (rho 0.02)**; it is epithelial-purity (+0.84) vs immune-infiltration (-0.90)
  -- the residual enzymatic-dissociation contamination. fGSEA: OXPHOS/MYC-target pole vs
  immune/inflammation pole (allograft rejection NES -3.5 padj 1e-32, IFN-g, IL6/JAK/STAT3).
- Its metabolic/OXPHOS pole is a **cell-IDENTITY** effect (epithelial MECs are mito-dense, immune
  infiltrate is not), NOT Myc dose. This is the axis the design/compositing correctly removes.

**gene-PC2 (13%) IS the Myc-mito program -- FOREGROUNDED by the pathway analysis.**
- MYC/E2F/proliferation (prolif 0.76, myc 0.49) with **OXPHOS + mitochondrial translation/mtRNA
  (biogenesis) baked in**: fGSEA E2F 1e-52, MYC 1e-45, mito-translation 1e-35, OXPHOS, G2M. So
  **mitochondria ride INSIDE the Myc program**, not the composition nuisance.
- The pathway global axis tracks **gene-PC2 more than gene-PC1 (cor 0.60 vs 0.44)** and **SHARPENS the
  genotype signal (gene-PC2 genotype 0.26 -> pathway-PC1 0.66)**. The z-score/compositing step promotes
  the modest gene-PC2 Myc-mito program to the dominant (76%) axis and demotes the 44% composition axis.

**THE NARRATIVE FLIP (author agreed):** NOT "PC1 is the mito axis." It is **"gene-PC1 (44%) is what we
correctly discard (composition/immune); the Myc-mito program is gene-PC2, which the pathway analysis
foregrounds."** Characterising PC1 is what PROVES the global axis is not merely contamination.

**CAVEAT (the ceiling's gene-level root):** OXPHOS enriches at BOTH gene-PC1's metabolic pole
(composition) AND gene-PC2 (Myc program) -- this dual loading is why mito correlates with everything.
See [[correlation-ceiling]], [[figure-b-loaders-are-curated-not-top]], [[batch-equals-timepoint]].

**How to apply:** the defensible mito-Myc statement is gene-PC2 -> pathway-PC1 (genotype-clean); do NOT
cite gene-PC1's OXPHOS pole as Myc biology (it is composition, genotype-independent). Figures in script
37: A3 (PC x covariate heatmap), A4 (fGSEA assigning the axes to pathways), A2 (sample PCA),
A_scree (two-panel gene vs pathway scree). Pairs with the content result [[mito-content-bounded]].
