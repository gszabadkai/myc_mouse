---
name: batch-equals-timepoint
description: "The 6W and 12W cohorts were extracted as two separate batches, so batch is perfectly confounded with timepoint. Genotype is balanced within each batch, so the design absorbs it: genotype/interaction clean, developmental (pure-time) reading confounded."
metadata:
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
---

**Load-bearing technical fact, from the author 2026-07-19** (undocumented until then, like the
enzymatic-dissociation detail). The RNA-seq was **two extractions -- 6W and 12W processed
separately** -- so **batch is perfectly confounded with timepoint.**

**Why it is mostly reassuring, not alarming.** Genotype is **balanced within each batch** (each
timepoint has both WT and Myc+). So the DESeq2 design `~timepoint*myc_status` **absorbs the entire
batch effect into the timepoint main-effect term**. Consequences:
- **Every genotype and interaction contrast stays clean** -- the backbone is safe (content +21-27%,
  nuclear OXPHOS, the attenuation, amplification, BMYO suppression, all genotype/interaction).
- **Only the pure between-timepoint (developmental) reading is confounded with batch** and cannot be
  separated from it -- there is **no external batch/RIN/QC metadata on disk** to adjudicate.

This is the *mechanism* behind the caution we already held from script 33: the global/IEG axis is
**genotype-INDEPENDENT (p=0.49) but time-associated (p=0.0003)**. Batch=timepoint is why. Same
verdict, now with a named cause: **genotype clean, time exposed.**

**How to apply.** Any code or narrative that removes the design first (script 36's design->PC1->
residual; script 37's biological axis) is orthogonal to batch by construction -- fine. Do NOT make a
clean "developmental vs Myc-specific" claim from a pure-timepoint contrast without flagging the
batch confound. The genotype/interaction story never needs this caveat. mtDNA/imbalance TIME claims
already carry it (see [[mitonuclear-imbalance-untested]], [[timepoint-cohort-confound]],
[[correlation-ceiling]]). The clean test for anything time-confounded is an orthogonal measurement
(qPCR, blot), not more of this bulk RNA.
