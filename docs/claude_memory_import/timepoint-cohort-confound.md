---
name: timepoint-cohort-confound
description: Timepoint is batch-aligned in the MMTV-Myc bulk data (disjoint depth ranges) and the mt share is high-variance; but depth itself is NOT a confound - shares and DESeq2 both cancel it
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
---

**Sequencing depth differs sharply between the two ages, and it is NOT a confound.** Corrected
2026-07-16 after the author asked "aren't the DESeq2 normalised values accounted for depth?" --
they are, and so are share-based measures. My first version of this memo claimed the opposite.

**Why depth cannot be the confound:**
- A **share is a proportion** (`colSums(panel)/colSums(all)`). Sequencing deeper multiplies
  numerator and denominator alike => **depth-invariant by construction**. DESeq2's
  median-of-ratios size factors cancel depth for the same reason. **There is no mechanism.**
- **Empirically it does not act**: within cohort depth predicts mt% not at all --
  6W rho=-0.056 (p=0.87), 12W rho=-0.315 (p=0.32).
- The scary pooled **rho=-0.469 (n=24) is a POOLING ARTIFACT** -- two clouds at disjoint depths
  that also differ in mt%. Never quote the pooled correlation.
- **12W is NOT the degraded cohort**: top100 share of non-mt counts is **31% at 12W vs 40% at 6W**
  (12W is *more* even). Its lower detected-gene count (23,255 vs 24,599) follows from lower depth.
  The one quality-bearing metric argues *against* the degradation worry.

**What IS true (both duller than a confound):**
1. **Cohort/batch alignment — inherent, not mt-specific.** Depth ranges are **disjoint**: 6W
   18.8-29.0M (animals MYCF62-65 + MYBS102/104) vs 12W 9.9-16.1M (MYCF52-56). `depth ~ timepoint`
   p=4.2e-09. No common support => batch and time are not separable. **But this is true of ANY
   cross-sectional design** (different animals per age) and the whole project already lives with
   it. Worth one methods sentence; not a reason to distrust mt specifically.
2. **mt imprecision — the real, mt-specific caveat.** The mt-* share spans **3.4%-40%** across
   samples (SD 13 pp within 12W_neg; one sample at 40.1%). Temporal mt claims are **IMPRECISE**
   at n=6/group.

=> Temporal mt numbers, incl. script 29's "mtDNA absolute z -1.14 -> +1.23", are **IMPRECISE and
COHORT-ALIGNED, NOT "uninterpretable"**. Report them; don't lean on their magnitude.

**GENOTYPE is clean** either way: depth is genotype-balanced (p=0.41) and both genotypes occur
inside single litters (MYCF62, MYCF52, MYCF56) — litter-controlled.

There is no RIN or batch column on disk (`data/FULL.DAT.COL.DATA.txt` has only
sample/group/cancer/age).

**How to apply:** don't over-warn. Say "batch is aligned with timepoint (cross-sectional design)
and the mt share is high-variance"; do NOT say "depth-confounded" or "time is uninterpretable".
Script 32 PART 4 + `qc$mt_depth_within` / `qc$depth_ranges` hold the evidence. See
[[mito-content-bounded]], [[dont-assert-a-mechanism-without-checking-it-can-work]].
Committed on `paper-figures`. Do not re-litigate.
