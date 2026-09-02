---
date: 2026-09-02
tags: [project/myc_mouse, project/human_validation, model, apoptotic-priming, oxphos, myc, bbc3, bcl2l11, mcl1, bcl2l1]
status: proposed model + mouse verification (scratchpad dry run; reconcile against the author's own run of script 48)
supersedes-in-part: docs/2026-08-27_human_validation_plan.md sections 7.1, 7.2, 7.3
relates-to:
  - scripts/48_gate_model_mouse_verification.R  (writes results/gate_model_verification.rds)
  - scripts/42_priming_arm_and_teb_substrate.R  PART H  (the coincidence test, published form)
  - scripts/43_substrate_specificity_and_tradeoff.R PART B (the trade-off asymmetry, positive control)
  - docs/2026-07-25_death_narrative_dose_vs_competence.md
  - docs/2026-07-26_introduction_alignment_and_the_question.md
---

# The gate model: MYC x OXPHOS -> apoptotic priming, and what the mouse says about testing it in human tumours

## 0. What this document is

The human arm needs a model, not a hypothesis list. A model is a set of equations with
named parameters, predicted signs, and a rule that says which observations would kill it.
This document states that model, verifies the part of it the mouse can verify, and reports
three places where the verification **contradicts the measurement choices in the current
human plan**.

The model is not new biology. It is the coincidence / trade-off frame this project has been
working with since 2026-07-25, written down as an equation so that it can be fitted in a
cohort that has no genotype, no timepoints and no perturbation.

**Scope, stated once.** n = 24; batch = timepoint; every model below is a correlation among
24 mice. This is measurement design and ranking, not confirmatory inference. The causal
claims in the manuscript belong to the perturbations. Numbers below come from a verified
scratchpad run of `scripts/48_gate_model_mouse_verification.R`; the script carries a
positive control that reproduces script 43 PART B to five decimal places
(`+6.089`, `p = 0.00523`), and every number should be re-read from the author's own run.

---

## 1. The model

Per sample `i`, with all scores z-standardised:

| symbol | quantity | mouse | human |
|---|---|---|---|
| `M` | MYC **dose** | genotype; `Myc` transcript | `MYC` mRNA, 8q24 GISTIC, RPPA MYC |
| `X` | respiratory state **of the mitochondrial compartment** | `ox_rel` (below) | `ox_rel`, identical recipe |
| `T` | apoptotic trigger, a sensor:guardian log-ratio | `log2(Bbc3) - log2(Bcl2l1)` | `log2(BBC3) - log2(BCL2L1)` |
| `B` | guardian balance | `log2(Mcl1) - log2(Bcl2l1)` | `log2(MCL1) - log2(BCL2L1)` |
| `C` | covariates | epithelial, immune, proliferation | purity, leukocyte fraction, proliferation, PAM50, TP53, plate |

### 1.1 Layer 1 -- the gate (testable in mouse and human)

```
T_i  =  b0  +  bM * M_i  +  bX * X_i  +  bMX * (M_i x X_i)  +  g' C_i  +  e_i     (1)
```

| parameter | prediction | meaning |
|---|---|---|
| `bMX` | **> 0** | the whole claim: respiration primes **in proportion to MYC** |
| `bX`  | **= 0** at `M = 0`, **on `T`** | respiration alone does not prime; it is not an independent input. **Holds for the trigger, not for the guardian balance `B`** -- see section 3.2 |
| `bM`  | > 0 | MYC raises the trigger -- canonical, not the test |

If `bX = 0` holds, (1) collapses to a **one-parameter gate**

```
T_i  =  b0  +  bM * M_i  +  bMX * (M_i x X_i)  +  g' C_i  +  e_i                  (1')
```

which is the manuscript's title thesis -- *mitochondria integrate an oncogenic and a
metabolic input* -- written as an equation. (1') is the sharper thing to test: it is nested
in (1), so the collapse is a likelihood-ratio test with one degree of freedom.

The scale-free summary that transfers across species and platforms is the point at which
`dT/dX` changes sign:

```
M*  =  -bX / bMX                                                                  (2)
```

`M*` is the MYC level, in SD units, at which respiration stops being neutral and starts
priming -- the **asset-becomes-liability point**. It is a number, not a metaphor, and it is
comparable between a mouse gland and a tumour cohort because both `M` and `X` are z-scores.

### 1.2 Layer 2 -- the escape (human only, by construction)

A cohort of established tumours has been filtered by the gate. Two consequences follow, and
neither can be tested in this mouse dataset because it contains no tumours:

```
P(BUFFER_i)  =  logit^-1( c0 + cM * M_i + cX * X_i + cMX * (M_i x X_i) + h' C_i )  (3)
```

with `cMX > 0`: occupancy of the high-`M`, high-`X` corner **requires** anti-apoptotic
buffering. And a **density** prediction that needs no outcome data at all: in the (`M`, `X`)
plane, the unbuffered high-high cell is *under-occupied* relative to the product of its
margins. That is a chi-square on a 2 x 2 x 2 table and it is the cheapest test in the whole
human arm.

`cMX > 0` also implies that the human `bMX` in (1) is **attenuated** relative to the mouse:
the tumours that expressed the gate most strongly are the ones that are not in the cohort.
A weaker human `bMX` is therefore *predicted*, not a failure. Only a **sign reversal** falsifies.

### 1.3 Layer 3 -- the consequence (human only)

Among tumours in the high-`M`, high-`X` region, the **unbuffered** ones retain the trigger
and are chemo-sensitive; the buffered ones are not. This is the `STATE` variable already
pre-specified in `docs/2026-08-27_human_validation_plan.md` section 7.5, and the model does
not change it.

---

## 2. The measurement bridge: `ox_rel`

The published mouse interaction is fitted on **mitoPPS**, which cannot be compared across
cohorts (the pairwise-ratio baseline is composition-dependent). The human plan therefore
made the **absolute OXPHOS level** primary. The mouse says that swap loses the result.

`ox_rel` is the bridge:

```
ox_rel  =  mean z( nuclear MitoCarta OXPHOS subunits )  -  mean z( the rest of MitoCarta )
```

87 subunits against the other 999 MitoCarta genes. It is a **within-compartment share**,
like mitoPPS -- so the compartment-wide common mode cancels -- but it is a difference of two
z-composites, so it needs nothing but MitoCarta membership and an expression matrix and is
computed identically in any cohort of any species.

**Why it works, measured rather than asserted.** An interaction test needs its two inputs to
be separable, and `ox_rel` is the only ruler that is:

| ruler | r with MYC mRNA | r with MSigDB MYC | r with DoRothEA MYC | r with proliferation |
|---|---|---|---|---|
| `ox_lvl` (absolute level) | 0.61 | 0.81 | 0.88 | 0.80 |
| `ox_ppd` (mitoPPS) | 0.45 | 0.75 | 0.81 | 0.75 |
| **`ox_rel`** | **0.33** | **0.61** | **0.64** | **0.59** |

---

## 3. The mouse verification

### 3.1 The gate holds, and it is respiratory-specific

Equation (1), genotype form, endpoint `Bbc3:Bcl2l1`, covariates epithelial + immune.
`perm` is the within-timepoint permutation null of scripts 42/43 -- shuffling `X` inside each
timepoint keeps the design and breaks only the animal-to-animal link. **Its median is not
zero** (0.19-0.48 depending on the axis), and the excess over that median is the claim.

| axis | `bX` | `bMX` | p | + timepoint | + proliferation | perm null | perm pct | perm p |
|---|---|---|---|---|---|---|---|---|
| **`ox_rel`** | **-0.006** | **+0.857** | **0.0010** | +0.878 (0.0009) | +0.913 (0.0005) | 0.484 | 95.3 | **0.047** |
| `ox_ppd` (mitoPPS) | +0.051 | +0.779 | 0.0052 | +0.806 (0.0045) | +0.973 (0.0012) | 0.395 | 91.6 | 0.084 |
| `ox_lvl` (absolute) | +0.203 | +0.595 | 0.082 | +0.669 (0.051) | +0.953 (0.017) | 0.271 | 85.7 | 0.143 |
| `mitorib_rel` | +0.022 | +0.976 | 0.0061 | +0.912 (0.0088) | +1.155 (0.0021) | 0.189 | 99.7 | 0.003 |
| `oxasm_rel` (assembly) | +0.169 | +0.516 | 0.056 | +0.507 (0.063) | +0.580 (0.040) | 0.318 | 76.2 | 0.238 |
| `tca_rel` | +0.311 | +0.608 | 0.071 | +0.504 (0.114) | +0.590 (0.097) | 0.105 | 92.7 | 0.073 |
| `fao_rel` | -0.052 | **-0.346** | 0.292 | -0.299 (0.327) | -0.306 (0.304) | -0.448 | 65.9 | 0.341 |
| `redox_rel` | -0.678 | **+0.090** | 0.780 | +0.135 (0.673) | +0.097 (0.767) | -0.171 | 78.5 | 0.215 |

Three things to read off it.

1. **The ordering is `ox_rel` > `ox_ppd` > `ox_lvl`, consistently, in every column.** The
   ruler the human plan makes primary is the weakest of the three and does not clear p < 0.05
   without help.
2. **`redox_rel` remains the control that behaves** -- +0.090, p 0.78, exactly as in script 38
   and script 43. `fao_rel` is negative. Neither is a respiratory arm.
3. **The gate is not proliferation.** Adding a proliferation composite *strengthens* every
   respiratory row (`ox_rel` 0.857 -> 0.913, p 0.0010 -> 0.0005). This is the negative control
   that matters most, because `X` and proliferation correlate at 0.59.

**The one arm that rivals OXPHOS is the mitoribosome, and it does not survive the
head-to-head.** `mitorib_rel` correlates with `ox_rel` at 0.755. Fitted together:

```
T ~ myc * ox_rel + myc * mitorib_rel + epi + imm
    myc:ox_rel        +1.086   p = 0.040
    myc:mitorib_rel   -0.354   p = 0.595
```

The respiratory term takes the whole interaction; mitochondrial translation contributes
nothing once respiration is in the model. Report the head-to-head with the marginal row --
without it, "OXPHOS-specific" is not supported at n = 24.

### 3.2 The collapse to a pure gate holds for the trigger -- and **not** for the guardian balance

Endpoint `T` = `Bbc3:Bcl2l1`:

| MYC estimator | `bX` | p for `bX` | R2 additive | R2 interaction | R2 gate-only (drop `X`) | **p for dropping `X`** | `M*` |
|---|---|---|---|---|---|---|---|
| genotype | **-0.006** | 0.973 | 0.643 | 0.807 | **0.807** | **0.973** | +0.007 |
| MYC transcript | +0.372 | 0.022 | 0.627 | 0.760 | 0.677 | 0.022 | -1.08 |
| MYC signature | +0.273 | 0.198 | 0.637 | 0.663 | 0.630 | 0.198 | -1.97 |

With the experimentally-set MYC dose, **dropping the OXPHOS main effect costs exactly
nothing** (identical R2, p = 0.973). Equation (1) collapses to (1'). The trigger is a
**product**: no MYC, no effect of respiration.

Endpoint `B` = `Mcl1:Bcl2l1`, the co-primary promoted in section 4.3 -- **and here the
collapse is rejected**:

| MYC estimator | `bX` | p for `bX` | R2 additive | R2 interaction | R2 gate-only (drop `X`) | **p for dropping `X`** | `M*` |
|---|---|---|---|---|---|---|---|
| genotype | **-0.541** | **0.0040** | 0.607 | 0.828 | **0.724** | **0.0040** | **+0.545** |
| MYC transcript | -0.160 | 0.200 | 0.641 | 0.843 | 0.828 | 0.200 | +0.377 |
| MYC signature | -0.282 | 0.184 | 0.531 | 0.663 | 0.628 | 0.184 | +0.915 |

**The guardian balance is a genuine crossover, not a gate.** In the MYC-null gland the
`MCL1:BCL-XL` ratio *falls* as respiration rises (`slope_wt` -0.512, covariate-adjusted;
`Mcl1` itself at -1.00 against a flat `Bcl2l1`), and MYC reverses the sign. `M*` = **+0.545**
lands between the two genotypes rather than at either of them -- the only place in this
analysis where equation (2) returns a number that is both positive and inside the observed
range. Two honest limits on it: the crossover reaches p < 0.05 only on the **genotype**
estimator (transcript p 0.20, signature p 0.18 -- same sign, consistent `M*`, so this reads
as power rather than contradiction), and `M*` is a ratio of two estimated coefficients at
n = 24, so treat it as a location, not as an interval.

**Consequence for the human model.** Fit (1) with the `X` main effect retained for `B` and
report the collapse test rather than assuming it. The one-parameter form (1') is licensed
for `T` alone.

**Correction to the record.** The internal note carried "within the wild-type animals the
OXPHOS -> PUMA:BCL-XL slope is -4.80, R2 0.328, p 0.052" -- i.e. a genuine crossover, with
respiration *protective* in the absence of MYC. That fit has no covariates. It reproduces
here (`ox_ppd`, no covariates: WT slope R2 0.33, p 0.053), but **it does not survive
adjustment for epithelial and immune content**:

| axis | covariates | WT slope | p | Myc+ slope | p |
|---|---|---|---|---|---|
| `ox_rel` | none | -0.566 | 0.050 | +0.657 | 0.024 |
| `ox_rel` | epi + imm | **-0.012** | 0.95 | **+0.870** | 0.0019 |
| `ox_ppd` | none | -2.193 | 0.053 | +2.017 | 0.042 |
| `ox_ppd` | epi + imm | **+0.138** | 0.88 | **+2.966** | 0.0036 |

So the honest statement is a **gate, not a crossover**: in the MYC-null gland the slope is
zero, not negative. `M*` from equation (2) is +0.007 SD on the genotype fit -- i.e. the sign
change sits essentially at the MYC-null state, which is what (1') says. Do not write
"respiration is protective without MYC" in the paper; the composition-adjusted data do not
support it. Write "respiration does nothing without MYC".

### 3.3 Which half of the ratio carries the gate -- the finding that changes the human endpoint

Every endpoint in the menu shares the `Bcl2l1` denominator. Fitted alone, against an
expression-matched null (840 genes per bin, 20 bins) and against all 16,819 expressed genes:

| gene | `bMX` alone | p | pct of all genes | **pct expression-matched** | within-timepoint p |
|---|---|---|---|---|---|
| `Mcl1` | **+0.806** | **0.018** | 97.1 | **98.9** | **0.027** |
| `Bbc3` (PUMA) | +0.746 | 0.059 | 95.9 | 94.8 | 0.44 |
| `Bmf` | +0.585 | 0.099 | 91.2 | 88.5 | 0.028 |
| `Bak1` | +0.493 | 0.123 | 87.0 | 84.2 | 0.45 |
| `Bid` | +0.342 | 0.362 | 78.0 | 76.6 | 0.78 |
| `Bax` | +0.104 | 0.707 | 57.4 | 53.0 | 0.38 |
| `Myc` | +0.047 | 0.713 | 51.6 | 48.3 | 0.82 |
| `Bcl2l11` (BIM) | **-0.072** | 0.856 | 38.7 | **37.6** | 0.37 |
| `Pmaip1` (NOXA) | -0.430 | 0.033 | 10.0 | 7.7 | 0.32 |
| `Bcl2l1` (BCL-XL) | **-0.670** | **0.0026** | 2.7 | **2.1** | **0.0026** |

And the endpoint menu, full model and within-timepoint (the cross-sectional regime):

| endpoint | `bMX` | p | `bMX` within-timepoint | p within |
|---|---|---|---|---|
| **`Mcl1:Bcl2l1`** | **+0.994** | **0.00014** | **+0.351** | **0.00058** |
| `Bbc3:Bcl2l1` (PUMA priming) | +0.857 | 0.0010 | +0.288 | 0.053 |
| `Bak1:Bcl2l1` | +0.807 | 0.011 | +0.294 | 0.068 |
| `Bmf:Bcl2l1` | +0.805 | 0.031 | +0.591 | 0.0044 |
| `Bid:Bcl2l1` | +0.588 | 0.076 | +0.148 | 0.45 |
| `Bax:Bcl2l1` | +0.443 | 0.022 | +0.124 | 0.29 |
| `Bcl2l11:Bcl2l1` (BIM) | +0.437 | 0.229 | +0.103 | 0.36 |
| **`Bbc3:Mcl1`** | **+0.183** | **0.63** | -0.063 | 0.63 |
| `Pmaip1:Bcl2l1` (NOXA) | -0.102 | 0.52 | +0.050 | 0.70 |

**Four conclusions, and the third is uncomfortable.**

1. **BIM is flat.** `Bcl2l11` sits at the 37.6th percentile of its expression-matched null --
   dead centre -- and the `BIM:BCL-XL` ratio is null in every regime. This is the fourth
   independent confirmation that the PUMA/BIM pair from the PGC1a westerns **splits in vivo**
   at the RNA level (script 44's collapse statistic put `Bbc3` at the 0.51st percentile and
   `Bcl2l11` at the 91.6th). BIM stays in the human pre-specification, because the westerns
   pre-specified it and a positive would strengthen the paper -- but its expected value is
   zero, BIM protein is set post-translationally (ERK-driven degradation), and a null result
   is **uninformative about the protein**. Say that in Methods, before fitting.
2. **NOXA is a clean negative in the opposite direction** (7.7th percentile). Free specificity.
3. **The gate's strongest single component is the guardian, not the sensor.** `Bcl2l1` at the
   2.1st percentile is more extreme than `Bbc3` at the 94.8th; `Bcl2l1` and `Mcl1` are the only
   two members of the roster significant in **both** the full model and the within-timepoint
   restriction, and `Bbc3` is significant in neither (p 0.059 and 0.44). `Bbc3:Mcl1` -- the same
   numerator over a different denominator -- is **null (+0.183, p 0.63)**. So the pre-specified `PUMA:BCL-XL` endpoint is carried roughly half by
   PUMA rising and half by BCL-XL falling, and it is the **balance**, not PUMA, that the gate
   moves.
4. **`Mcl1:Bcl2l1` is the most robust endpoint in the whole analysis** (+0.994, p 0.00014;
   within-timepoint +0.351, p 0.00058; permutation 99.4th percentile, p 0.0056). As MYC-driven
   respiration rises, **BCL-XL is withdrawn while MCL1 is retained** -- a switch of guardian
   dependence, visible in normal epithelium, before any tumour exists. It is also the one
   endpoint with a **two-sided** shape: section 3.2 shows the wild-type slope is genuinely
   negative, so `B` crosses over at `M*` = +0.545 rather than switching on at zero.

**This does not overturn the manuscript's PUMA claim.** That claim is about a *different*
statistic: the departure of `Bbc3` from the x0.55 dose line between 6W and 12W (script 44;
rank 1 and 2 of 26 mechanism genes with `Foxo3`), plus the PGC1a westerns. This section is
about a *cross-sectional coupling within the cohort*. Both can be true, and the endpoint
that transfers to a cross-sectional human cohort is not necessarily the one that carries a
longitudinal collapse. Say which statistic a sentence is about, every time.

### 3.4 The cross-sectional regime -- the honest transfer estimate

Every variable centred inside its timepoint, so the cohort contrast (which is also the batch
contrast) is gone entirely. This is the only variation a human cross-section has:

| axis | `Bbc3:Bcl2l1` | p | `Mcl1:Bcl2l1` | p |
|---|---|---|---|---|
| `ox_rel` | +0.288 | **0.053** | +0.351 | **0.00058** |
| `ox_ppd` | +0.244 | 0.101 | +0.348 | 0.0011 |
| `ox_lvl` | +0.144 | 0.403 | +0.366 | 0.0033 |

The PUMA gate is **marginal** in this regime at n = 24 and only on `ox_rel`; the buffering
switch is solid on all three rulers. Quote these numbers, not section 3.1's, whenever the
claim is cross-sectional.

---

## 4. Three corrections the mouse forces on the human plan

### 4.1 The MYC estimator ordering is wrong (plan section 7.1)

The plan makes the **MYC target signature** primary (M-a) and **8q24 amplification** a third
estimator. Each estimator's coefficient is per SD of itself, so they are not comparable as
printed; multiplying by that estimator's own genotype gap puts all of them on the genotype
scale and answers the only question that matters -- *how much of the experimentally-created
effect does this instrument see?*

| estimator | genotype gap (SD) | `bMX` per SD | p | in genotype units | **fraction recovered** |
|---|---|---|---|---|---|
| `Myc` transcript | 1.85 | +0.343 | **0.0055** | +0.636 | **74%** |
| MSigDB MYC targets (MitoCarta-stripped, 291 genes) | 0.93 | +0.138 | 0.249 | +0.129 | **15%** |
| DoRothEA MYC regulon (stripped, 302 genes) | 1.09 | +0.123 | 0.332 | +0.134 | **16%** |

Both signature scores separate the genotypes by only ~1 SD, and both correlate with `ox_rel`
more strongly than the transcript does -- 0.61 (MSigDB) and 0.64 (DoRothEA) against 0.33 --
so the signature partly *contains* the axis it is supposed to interact with, and the
interaction becomes nearly unidentifiable.

**Consequence for the human arm: prefer dose estimators to activity estimators.** `MYC` mRNA
and 8q24 GISTIC amplification are measured as independently of the mitochondrial axis as
anything available (amplification is a DNA-level quasi-instrument and is the cleanest of
all); a target signature is downstream of both inputs. Keep the signature as the concordance
check the plan requires -- but if the signature and the dose estimator disagree, the mouse
says the dose estimator is the one to believe, and the plan currently says the opposite.

**The caveat that must travel with this.** In this mouse the dose is experimentally set and
the `Myc` transcript reads the transgene directly, so the mouse cannot prove that MYC mRNA
beats a signature *in a tumour*, where mRNA is a much noisier proxy for activity. What the
mouse does establish is the **effect size on the signature scale**, and that number is what
the human power calculation has to use.

### 4.2 The primary axis measure is wrong (plan section 7.2)

Replace "Level (primary): mean z-score of nuclear-encoded MitoCarta OXPHOS subunits" with
**`ox_rel`, the compartment-relative respiratory score** (section 2 above). Keep the absolute
level as a reported sensitivity, not as the primary. mitoPPS stays as the shape reading and
still must never be compared numerically across cohorts.

### 4.3 The endpoint set needs a fourth member, promoted (plan section 7.3)

Pre-specify **two co-primaries**, in this order:

- `BUFFER = log2(MCL1) - log2(BCL2L1)` -- the guardian switch. Strongest and most robust in
  the mouse, survives the cross-sectional restriction, and carries a drug-shaped prediction
  (MCL1 dependence: S63845 / AMG-176; DepMap MCL1 Chronos; BH3 profiling with MS1 against HRK).
  **Fit it with the `X` main effect retained** (section 3.2): unlike `PRIME`, this endpoint
  does not collapse to the one-parameter gate, and `M*` is the quantity to report for it.
- `PRIME = log2(BBC3) - log2(BCL2L1)` -- unchanged, the cross-species continuity endpoint.

and add two **decomposition endpoints** that are mandatory whenever `PRIME` is positive,
because without them a positive `PRIME` cannot be attributed:

- `BBC3` and `BCL2L1` fitted **alone**. If the signal is in `BCL2L1`, the sentence is about
  the balance, not about PUMA.
- `log2(BBC3) - log2(MCL1)` -- the same numerator over the other denominator. Null in the
  mouse; a positive in human would mean the human gate is broader than the mouse gate.

Negative-control endpoints unchanged, with the mouse's expected values now attached:
`BCL2L11:BCL2L1` (null, 0.44 / p 0.23), `PMAIP1:BCL2L1` (**negative**, -0.10), `BAX:BCL2L1`
(weak, rides the denominator), `BID:BCL2L1` (weak).

---

## 5. Falsification

The model is dead if **any** of the following:

- `bMX <= 0` in equation (1) for `ox_rel` in the full cohort **and** in the TP53-mutant,
  PIK3CA-wild-type and within-PAM50 strata, on both co-primary endpoints;
- `bMX` is as large for `redox_rel`, `fao_rel` or `tca_rel` as for `ox_rel` (the gate is
  mitochondrial-in-general, not respiratory);
- `bMX` is as large for endpoints the model does not name (`BAX`, `BID`, `NOXA`) as for the
  co-primaries, with `BCL2L1` alone flat -- the gate would then be a general apoptotic shift;
- `cMX <= 0` in equation (3) **and** the unbuffered high-`M`/high-`X` cell is *over*-occupied;
- OXPHOS-high predicts chemo**resistance** regardless of buffering status (this inverts the
  model rather than merely failing to support it).

**Not falsifying:** a human `bMX` smaller than the mouse's. Section 1.2 predicts attenuation.
Only a sign reversal counts.

---

## 6. Power

SE scales as `1/sqrt(n)` at fixed design correlation; the mouse SE at n = 24 fixes the curve.
2.8 SE is the usual 80% / two-sided 5% rule of thumb.

| estimator | `bMX` per SD (mouse) | SE | n for the full mouse effect | n for half | **n for a third** |
|---|---|---|---|---|---|
| `Myc` transcript | +0.343 | 0.109 | 19 | 76 | **170** |
| MYC signature | +0.138 | 0.116 | 133 | 531 | **1193** |
| DoRothEA regulon | +0.123 | 0.124 | 190 | 759 | **1708** |

TCGA-BRCA (n ~ 1000) and METABRIC (n ~ 1900) are adequate **on a dose estimator** even if the
human effect is a third of the mouse's, and adequate **on a signature estimator** only if the
human effect is at least half the mouse's. That asymmetry is a second, independent reason to
promote the dose estimators. The plan's existing caveat -- interaction tests need roughly four
times the n of a main effect -- is consistent with these figures and still applies.

---

## 7. What this does not claim

- **Nothing here is causal.** Twenty-four correlated mice. The gate is a description of
  covariation; the manuscript's causal claims rest on the MYC-ER, PGC1a and Bcl-xL
  perturbations, and this document does not add to them.
- **Layers 2 and 3 are untested and untestable here.** No tumours, no outcome data.
- **`Ppargc1a` is not in the model.** baseMean ~30 in MEC; coactivator activity is
  post-translational. The model uses the respiratory *state*, which is what is measurable, and
  says nothing about which upstream factor sets it -- consistent with
  `docs/2026-08-17_developmental_oxphos_decline_and_the_biogenesis_axis.md`, which closed that
  question in the negative.
- **The gate is not the developmental decline.** The developmental fall is an arm-selective
  loss of respiratory subunits between two cohorts; the gate is a within-cohort coupling.
  Script 48 PART G separates them and the paper should too.
- **`Mcl1`'s rise is relative.** Both `Mcl1` slopes are negative in absolute terms
  (`gene_gate`, epi + imm adjusted: WT **-1.00**, Myc+ **-0.20**); the ratio rises because
  BCL-XL falls faster still (`Bcl2l1` WT +0.01, Myc+ **-0.66**). Write "the guardian balance shifts
  toward MCL1", not "MCL1 is induced".

---

## 8. Provenance

Every number above is read from `results/gate_model_verification.rds`, written by
`scripts/48_gate_model_mouse_verification.R`, which reads only committed objects
(`dds_int_run.rds`, `combined_df_annotated.rds`, `mitopps_scores.rds`) and the library GMT,
takes set membership exclusively through `functions/reconcile_gene_symbols.R`, and stops on a
positive control against script 43 PART B before computing anything new. `set.seed(1)`,
`NPERM = 5000`, `NBIN = 20`, genes filtered at mean normalised count >= 10 (16,819 genes).

Script 48's own verdict table, on a rule fixed before the tests:

| axis | verdict | `bMX` | p | perm p | fraction recovered by the signature estimator |
|---|---|---|---|---|---|
| `ox_rel` | **transfers_weak** | +0.857 | 0.0010 | 0.047 | 0.150 |
| `ox_ppd` | transfers_weak | +0.779 | 0.0052 | 0.084 | 0.086 |
| `ox_lvl` | **fails** | +0.595 | 0.082 | 0.143 | 0.071 |

"transfers_weak" is the honest headline: **the gate is real in the mouse and the ruler
question is settled, but no MYC estimator available to a human cohort recovers more than
~15% of it.** That is why section 6 exists, and why the human arm should lead with the
endpoint that survives the cross-sectional restriction (`MCL1:BCL2L1`) rather than the one
that is marginal in it.
