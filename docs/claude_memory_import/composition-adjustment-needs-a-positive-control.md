---
name: composition-adjustment-needs-a-positive-control
description: "Two traps found on the fat-pad cohort that generalise to the human arm: a loading is a property of a SUBSET not a ruler, and any composition/purity adjustment must carry a positive control or it removes the biology"
metadata:
  type: feedback
---

**A loading is a property of a SUBSET, not of a ruler.** An admissibility argument ("the
confound runs against my prediction, so a positive result is credible") must be computed on
the **exact subset the test runs on**. Script 49: `rho(ox_sub, Adipoq)` is +0.538 over all 30
fat-pad samples and **+0.144** within the 20-sample limb actually tested, and the
mitoribosome's loading **flips sign** between the two. The pre-registered admissibility
argument was built on the full-cohort numbers and did not survive.

**A difference of two composites cancels a confound only when it loads on both with the same
sign and similar magnitude.** `ox_nuc_mtrib` (nuclear OXPHOS minus mitoribosome) has halves
loading +0.14 and **-0.24** on adipose, so the subtraction ADDS the tissue signal: the
endpoint ends up +0.58 to +0.78 across four adipose markers and **76% of its variance**
adipose-explained, against **7%** for `ox_rel`. The construction that isolates the biological
claim was the one the tissue gradient most contaminated.

**ANY composition/purity adjustment must carry a POSITIVE CONTROL through it, read FIRST.**
Script 49 PART G: `Mki67` genuinely rises across the limb (tau +0.313, p 0.045) and goes to
**-0.096 (p 0.715)** after adjusting on four adipose markers -- while those markers explain
only 15% of its variance. Adipose depletion and tumour progression are the same variable
there, so covarying out composition removed the biology. Without the control this would have
read as a clean null.

**Why:** all three were caught only because the tests were pre-registered with a hard
admissibility requirement and a negative control, and each would otherwise have produced a
confident wrong sentence. **How to apply:** this is flagged for the **human arm's trap list**
-- a cohort stratified by PAM50 / TP53 / purity is exactly where a full-cohort loading stops
holding in the stratum tested, and the human plan pre-specifies stratified analyses plus
purity adjustment. Related: [[correlation-ceiling]], [[sample-level-couplings-need-a-null]],
[[hypothesis-generating-not-failed-tests]], [[myc-oxphos-priming-gate-model]].
