---
name: myc-oxphos-priming-gate-model
description: "The MYC x OXPHOS -> priming model written as an equation for the human arm, verified in the mouse by script 48 -- and the three measurement corrections it forces on the human validation plan"
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-09-01T23:59:23.885Z
---

**The model, 2026-09-02.** `T = b0 + bM*M + bX*X + bMX*(M x X) + covariates`, with `T` a
BH3-sensor:guardian log-ratio, `M` MYC **dose**, `X` the respiratory state of the mitochondrial
compartment. Predictions: `bMX > 0` (the whole claim), `bX = 0` at `M = 0`, crossover
`M* = -bX/bMX` in SD of M. Two further layers exist only in tumours and are human-only by
construction: `P(BUFFER) ~ M x X` (the escape) and the conditional-chemosensitivity layer.
Doc: `docs/2026-09-02_myc_oxphos_priming_gate_model.md`. Script: `scripts/48_gate_model_mouse_verification.R`
-> `results/gate_model_verification.rds`, positive control reproduces script 43 PART B (+6.089, p 0.00523).

**THE TRIGGER IS A PURE GATE; THE GUARDIAN BALANCE IS A CROSSOVER. Never state the collapse
without naming the endpoint.** For `T = Bbc3:Bcl2l1`, dropping the OXPHOS main effect costs
**nothing** (R2 0.807 both ways, anova p 0.973) -- so `T = b0 + bM*M + bMX*(M x X)`, one parameter
fewer, and it is the manuscript's title thesis as an equation; the wild-type slope is **zero, not
negative**, once epi+imm are in the model (see [[foxo3-puma-pair-and-the-nonlinear-dedosing]] 6a).
For `B = Mcl1:Bcl2l1` the collapse is **REJECTED** (`bX` -0.541 p 0.0040; R2 0.828 -> 0.724):
the wild-type slope is genuinely negative and survives adjustment, so `B` crosses over at
**`M*` = +0.545**, between the genotypes -- the only place equation (2) returns a positive number
inside the observed range. Significant on the genotype estimator only (transcript p 0.20,
signature p 0.18, same sign) -- power, not contradiction. **The human model must retain the X
main effect for BUFFER.** Caught by the author's own run of 48, not by the dry run.

**THE RULER QUESTION IS SETTLED: use `ox_rel`.** `ox_rel = mean z(87 nuclear OXPHOS subunits) -
mean z(the other 999 MitoCarta genes)` -- a within-compartment share like mitoPPS but portable to
any cohort or species. Gate coefficient **+0.857 p 0.0010 (perm p 0.047)** vs mitoPPS +0.779
(0.0052) vs **absolute level +0.595 (0.082, FAILS)** -- the ruler the human plan makes primary is
the weakest. The reason is separability: `ox_rel` correlates with MYC mRNA at **0.33** where the
absolute level is at 0.61, and with proliferation 0.59 vs 0.80. **The gate is not proliferation** --
adding a proliferation composite makes it STRONGER (0.857 -> 0.913, p 0.0005). `redox_rel` (+0.09,
p 0.78) and `fao_rel` (-0.35) are the controls that behave. The mitoribosome rivals it marginally
but **loses the head-to-head** (myc:ox +1.086 p 0.040 vs myc:mitorib -0.354 p 0.59) -- always
report the head-to-head, or "OXPHOS-specific" is unsupported at n=24.

**THE GATE'S STRONGEST COMPONENT IS THE GUARDIAN, NOT THE SENSOR -- and this changes the human
endpoint.** Fitted alone against expression-matched nulls: `Bcl2l1` **-0.670, 2.1st pct, p 0.0026**;
`Bbc3` +0.746, 94.8th, **p 0.059**; `Mcl1` +0.806, **98.9th**, p 0.018. `Bbc3:Mcl1` -- same
numerator, other denominator -- is **null (+0.183, p 0.63)**. So `PUMA:BCL-XL` is a **balance**,
roughly half carried by BCL-XL falling. **`Mcl1:Bcl2l1` is the most robust endpoint in the whole
analysis** (+0.994 p 0.00014; **within-timepoint +0.351 p 0.00058**; perm 99.4th pct) -- as
MYC-driven respiration rises the gland withdraws BCL-XL and retains MCL1, a guardian switch visible
before any tumour. Pre-specify it as co-primary with PRIME, and it carries the drug prediction
(S63845/AMG-176, DepMap MCL1 Chronos, BH3 profiling MS1 vs HRK). **BIM is flat** (37.6th pct, dead
centre) -- fourth independent confirmation the PUMA/BIM pair splits in vivo; keep it pre-specified
but state that a null is uninformative about BIM protein. **NOXA is a clean negative** (7.7th pct).

**THE TRANSFER IS WEAK AND THAT IS THE HEADLINE.** Rescaled by each estimator's own genotype gap,
`Myc` transcript recovers **74%** of the genotype effect (p 0.0055), the MSigDB MYC signature
**15%** (p 0.25), DoRothEA **16%** (p 0.33) -- both signatures separate the genotypes by only ~1 SD
and correlate with `ox_rel` at 0.61 (MSigDB) / 0.64 (DoRothEA) against the transcript's 0.33, so they
partly CONTAIN the axis. **Prefer DOSE estimators
(MYC mRNA, 8q24 GISTIC) over activity signatures** -- the human plan currently says the opposite.
Caveat that must travel: the mouse dose is experimentally set, so this cannot prove mRNA beats a
signature *in a tumour*; what it fixes is the effect size to power on. n needed at 80%: **170** on
a dose estimator for a third of the mouse effect, **1193** on a signature.

**Cross-sectional regime (timepoint means removed -- the only variation a human cohort has):** PUMA
+0.288 **p 0.053**, BUFFER +0.351 **p 0.00058**. Quote these, not the full-cohort numbers, whenever
the claim is cross-sectional.

**Why:** the human arm was about to be built on three measurement choices that the mouse can
adjudicate for free, and two of the three were wrong. **How to apply:** this does NOT overturn the
manuscript's PUMA claim -- that is the departure-from-the-dose-line statistic between 6W and 12W
(script 44, rank 1-2 of 26 with `Foxo3`), a different statistic on different variation. Say which
statistic a sentence is about, every time. Related: [[dose-vs-competence]],
[[mitochondrial-tradeoff-model]], [[oxphos-apoptosis-module-verdict]], [[correlation-ceiling]],
[[mtpct-is-three-way-confounded]], [[hypothesis-generating-not-failed-tests]].
