---
name: sample-level-couplings-need-a-null
description: "At n=6-12 with one dominant sample-level axis, every composite correlates with every other at ~0.6-0.8; any per-sample coupling claim needs an empirical null (the AP6/script-21 logic)"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
---

**In this dataset, any two per-sample composites correlate at ~0.6-0.8. A big r is not a finding.**
Established 2026-07-16 by `scripts/33_mtdna_axis_and_coupling_null.R`.

Null = all 884 GSVA-scored library sets; ask where the published coupling sits in the distribution
of that axis's couplings to *arbitrary* sets. Within 6W (n=12), couplings to `pro_comp`:

| coupling | rho | null median | percentile | perm p | verdict |
|---|---|---|---|---|---|
| mitonuclear imbalance | +0.67 | +0.62 | 61st | 0.49 | NOT SPECIFIC |
| bio_comp | +0.83 | **+0.79** | 54th | 0.50 | NOT SPECIFIC |
| MASC_comp | +0.78 | +0.57 | 90th | 0.10 | marginal |
| IEG/prep-stress axis | -0.78 | -0.55 | 0.5th | 0.03 | beats the null |

`bio_comp` r=0.83 looks compelling until you see the null median is **0.79**. 73-81% of ALL library
sets couple at |rho|>0.5. The only axis that beats the null is the **prep-quality signature**.

**The project already knew this.** `scripts/21_ap6_permutation_null.R` built exactly this kind of
null for fGSEA (AP6: "mito genes are highly expressed, so enrichment could be an artifact — test
against a matched null"). Nobody applied the same logic to the per-sample couplings in scripts
23/25. That gap is what let script 25 Part B's death couplings look mechanistic.

**How to apply:**
1. **Never report a per-sample composite coupling without an empirical null.** Reuse the
   `null_for_axis()` pattern in script 33 PART C.
2. Applies to script 25 Part B, script 24 PART C (cat-9 death coupling), script 23's substrate
   couplings — anything of the form "axis X correlates with composite Y across samples".
3. **Specificity failure is independent of the axis's origin** — it holds whether the driver is
   technical or biological. So it does not wait on the technical-vs-biological adjudication.
4. Genotype/interaction contrasts are a *different* design and do NOT have this problem — they are
   powered, litter-controlled, and the dominant axis is genotype-independent (p=0.49).

Caveat: the 884 sets are internally correlated and some overlap `pro_comp`, so the empirical p is
indicative, not exact — but a coupling at the **median** needs no precision to read.

**The structural reason:** within 6W, **PC1 = 43%** of variance and every measure loads on it (IEG
-0.71, proliferation +0.62, contamination -0.62, MYC +0.59, mt% -0.51). These are not independent
signals — so pairwise composite couplings are structurally uninformative at this n. Not a quirk of
the death sets.

See [[mitonuclear-imbalance-untested]], [[mito-content-bounded]],
[[dont-assert-a-mechanism-without-checking-it-can-work]]. Recorded in CLAUDE.md, the walkthrough
(script-33 section) and `docs/2026-07-17_evidence_audit_and_narrative.md`. Committed on
`paper-figures` through 2c5b79e.
