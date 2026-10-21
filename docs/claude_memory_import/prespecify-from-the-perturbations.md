---
name: prespecify-from-the-perturbations
description: "The author wants hypotheses driven towards what the later experiments actually show -- naming genes/sets in advance from the independent cell perturbations is pre-registration, not circularity, and it is the strongest available position at n=24"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-27T01:04:07.801Z
---

The author, 2026-07-27: *"I know that we are working here backwards — trying to use the RNAseq data
for hypothesis generation, but we could also drive the hypothesis towards what later experiments are
actually shown."*

**This is right and should be the default.** The cell/implant perturbations (PGC1a overexpression,
Bcl-xL rescue, MYAZ implants, MYC-ER) are **independent of the mouse RNA-seq**, so using them to fix
which genes and sets get tested is **pre-registration, not circularity**. At n=24, with a correlation
ceiling of ~0.8 in GSVA space and nothing surviving BH<0.05 in any coupling analysis, a
pre-specified test is worth far more than a scan hit — it converts "we found something in a big
ranking" into "we predicted this and it held".

**Why:** every blind scan in this corpus has died on multiplicity or on a null (scripts 34, 35, 36,
25 Part B). The two results that survived — PUMA's priming collapse and the mito-content effect —
were both anchored outside the transcriptome. Pre-specification is how this dataset earns claims.

**How to apply:**
- When designing a scan, **name the genes/sets from the perturbation data first**, in the script
  header, with the reason. Script 44 does this for `Bbc3` and `Bcl2l11`.
- **Report the pre-specified items even when they fail the scan's own filters** — `Bbc3` (6W padj
  0.22) and `Bcl2l11` (0.46) fail a padj gate, so script 44 reports on a wide gene set and ranks on
  a narrow one (`in_ranking` flag). Excluding them would have removed exactly what the scan was for.
- **Report the halves that fail.** The westerns give PUMA **and** BIM; in vivo `Bbc3` sits at the
  0.51st percentile and `Bcl2l11` at the 91.6th. Saying so is what makes the surviving half
  quotable — see [[oxphos-apoptosis-module-verdict]].
- Keep the epistemic contract intact: the in-vivo transcriptome **generates**, the perturbations
  **prove** ([[dose-vs-competence]]). Pre-specification does not upgrade a correlation to a cause.

Related: [[dont-assert-a-mechanism-without-checking-it-can-work]], [[sample-level-couplings-need-a-null]].
