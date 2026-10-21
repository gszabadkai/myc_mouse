---
name: dont-assert-a-mechanism-without-checking-it-can-work
description: Twice in one session I shipped an interpretation without checking its mechanism was even possible; both times the author caught it. Check the mechanism before writing the claim.
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
---

**Before asserting an interpretation, check that the mechanism it names could actually operate.**
Twice in one session (2026-07-16, script 32) I shipped a confident claim whose mechanism was
impossible or untested, and **the author caught both** — the second time after already running the
script.

**Failure 1 — "commissioned but unbuilt".** Read set-level `MITOCARTA_TRANSCRIPTION` +43% /
`MITOCARTA_MTDNA_NUCLEOID` +22% (d=2.32) as "Myc builds the mtDNA machinery". Never opened the sets.
Gene by gene the effects were carried by Mrpl12 +116% (a *mitoribosomal* protein), Atad3a +113%,
Poldip2 +58% — not by Tfam/Polg/Twnk, which barely move. **Tell I read past: a 9-11 gene set with
the largest effect size in the table.** See [[mitocarta-sets-are-membership-loose]].

**Failure 2 — "depth-confounded, TIME is not interpretable".** Cited mt% ~ depth rho=-0.47 as a
confound. But **a share is a proportion — depth cancels by construction**, exactly as DESeq2's
median-of-ratios does. Depth *cannot* drive a share. The author asked plainly: "aren't the DESeq2
normalised values accounted for depth?" Within cohort the correlation vanishes (rho=-0.06/-0.32,
both ns); the pooled value was a pooling artifact of two clouds. I also asserted 12W was the
degraded cohort while the one metric bearing on it (library complexity) said the opposite.
See [[timepoint-cohort-confound]].

**Failure 3 — "Myc raises nuclear mito but not mtDNA = the mitonuclear imbalance in absolute
share".** Built the headline of script 32 by noting nuclear was p=0.009 and mtDNA was p=0.83.
**That is comparing significance to non-significance and calling it a difference** — the commonest
statistical error there is. Tested as the difference it actually is (nuclear minus mtDNA per
sample): **p=0.57 raw, p=0.82 adjusted**. The author had already run the script twice.

**Failure 4 (near-miss, caught only because the author kept pushing) — the death couplings.**
I was one step from accepting "imbalance couples to pro-death priming, r=0.75, p=0.005" as a real
mechanism. It sits at the **61st percentile** of that axis's couplings to all 884 library sets
(null median +0.62, perm p 0.49). At n=12 everything correlates with everything at ~0.6-0.8.
**The project already knew this**: script 21 built a permutation null for exactly this reason on
the fGSEA side. Nobody applied it to the per-sample couplings. See
[[sample-level-couplings-need-a-null]].

**Failure 5 — the circular "direction test".** Argued mt% is not a mitochondrial readout because
nuclear mito genes anticorrelate with it (rho -0.65): "if mt% measured content, nuclear would rise."
**But the mitonuclear imbalance IS the claim that those two arms dissociate** — so the
anticorrelation is exactly what a *real* imbalance predicts. I assumed the phenomenon's absence in
order to prove the measurement broken. The author: *"that is exactly what imbalance would be,
wouldn't it?"* Withdrawn 2026-07-17. **Before using "X should track Y if Z were real" as a test, ask
whether the hypothesis under study already predicts they DON'T.**

**The common root:** I checked whether a pattern was *there*, never whether the claim I built on it
was the claim the test supported. Association got promoted to causal story; two p-values got
promoted to a difference; a big r got promoted to a finding without a null; and a test whose premise
the hypothesis already denied got promoted to a refutation.

**Every one was caught by the author asking a plain question, never by me.** "Shouldn't the nuclear
be there as well?" broke the contamination story open. "Aren't the DESeq values accounted for
depth?" killed Failure 2 in one line. When the author asks a simple question about a result I am
confident in, **that is the signal to test, not to explain.**

**2026-07-17 — THE pct ARTIFACT, TWICE, TWO DAYS APART.** Script 32's `pct_retained = 388%` (a
percentage where there was no effect to take a percentage of) was fixed on 2026-07-17. In script
34 I reintroduced the identical artifact by carrying script 31's `conv_pct = 100*wt_conv/atten`
across without checking that its premise — a substantial `atten` — holds for small death sets. It
printed **conv_pct -310% / fade_pct +410%**. Script 31's programmes all attenuate substantially,
so it never surfaced there. **Copying a formula copies its unstated preconditions.** Any ratio
needs a floor on its denominator. Caught by reading the author's output table, not by writing it.

**2026-07-17 — the failure recurs INSIDE the fix.** Writing script 34 (the script whose whole
point is "Gate 2 assumed gene independence"), I quoted **p=0.012 from a matched null that itself
assumes gene independence**. Script 21's own header had documented the caveat ("raw z inflated by
inter-gene correlation → use EFFECT-SIZE ranking, not absolute z"). Caught only by reading the
author's full output tables afterwards, not by writing the test. **Diagnosing an error class does
not exempt your own instrument from it.** Ask of every null you build: what does it assume, and
does my query violate it? Non-significant results under an anti-conservative null are safe;
significant ones need corroboration. See [[death-priming-change-untested]].

**2026-07-18 — over-nulling EXPLORATION is itself an error (author had to point it out).** I built
a full permutation-null apparatus to gate per-sample couplings, then framed the couplings as needing
to *pass* it. The author: *"do we HAVE TO use the floor? We are not seeking proof at this stage, we
are exploring."* Correct. A hypothesis needs prior + plausibility + a suggestive number when the
bench does the testing — **the problem was never r=0.75, it was p=0.005** (which asserts a test was
passed) and the claim level. **Calibrate rigor to STAGE:** a null is a RANKING tool for exploration,
a GATE only for a claim. But the mirror still holds — a ceiling-level lead "replicates" in any bulk
dataset (property of design+scoring), so validate ORTHOGONALLY, not with more bulk. And there is a
symmetric trap (script 36): **OVER-correcting** — when an axis IS the global factor (r²=0.94),
removing it leaves noise, and the resulting strong "coupling" is not a finding. Under- and
over-nulling are the same failure to match the instrument to the question. See [[correlation-ceiling]].

**How to apply:**
1. Before writing "X causes/explains Y", state the mechanism in one sentence and ask **could it
   physically work?** Ratios cancel scale. Proportions cancel depth. A median is robust to a few
   genes. If the arithmetic forbids the mechanism, stop.
2. **Aggregate -> mechanism requires disaggregation first.** Print the members with per-gene effect
   and expression before interpreting any set/composite result.
3. **A pooled correlation across groups is not evidence** until it survives within group.
4. **Do not over-warn.** Inability to *exclude* an artifact is not evidence *for* one. Calibrate the
   caveat to the positive evidence ("imprecise and cohort-aligned"), not to the worst case
   ("uninterpretable"). Overstated caveats cost credibility the same way overstated findings do.
5. **After retracting, do not reflexively substitute the mirror claim** — check it with equal
   rigour. Here "Myc *spares* the machinery" was also unsupported (p=0.12); the honest result was an
   informative negative.
6. Read-only validation caught neither: both claims *ran fine*. **Validation checks the code, not
   the interpretation.** Reasoning needs a separate pass.
