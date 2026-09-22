# Is Myc's mitochondrial effect shaped by a changing background?

**2026-07-24. Decision record for `sec-attenuation`.** Supersedes the loose "the WT
background matures toward Myc's state" framing of Issue #6's write-up. Numbers here come
from `results/background_vs_myc.rds` (script 40), built on the 2026-07-24
post-reconciliation run — **not** from the older docs, which predate the symbol fix.

Companion artefacts: `scripts/40_background_vs_myc_decomposition.R`,
`figures/fig03_background_vs_myc.R`.

---

## 1. The three rulers, both axes

144 MitoPathways, medians. content = DESeq2 set-average **raw** log2FC; priority = mitoPPS
pairwise diff (content-blind); share = % of the nuclear transcriptome (script 32).

| ruler | Myc@6W | Myc@12W | WT 6->12 | Myc+ 6->12 | interaction |
|---|---|---|---|---|---|
| content (log2FC) | **+0.43** (94% up) | **+0.23** (94% up) | +0.04 (72% up) | -0.13 (16% up) | **-0.18** (92% down) |
| priority (mitoPPS) | +0.025, sd 0.153 | -0.002, sd 0.111 | +0.003, sd 0.103 | -0.023, sd 0.105 | -0.027, sd 0.073 |
| share (nuclear MitoCarta) | +22.0% | +18.0% | +15.8% | +17.1% | int_p 0.79 |

The identity `(Myc@12W - Myc@6W) == (Myc+ temporal - WT temporal)` holds exactly on **both**
the content and the priority ruler (max deviation 1.1e-16 / 0), so the two decompositions
below describe the same quantity.

---

## 2. What the data support

### 2.1 Myc's effect is RESCALED, not reshaped

`Myc@12W ~ Myc@6W`: **slope 0.55** (bootstrap 0.51-0.60), R2 0.80, rho 0.84 on the content
ruler; slope 0.64, R2 0.79, rho 0.93 on priority. Per-tier slopes 0.53-0.83 — uniform. The
R2 sits at the **100th percentile** of within-decile gene-label shuffles that preserve set
size, expression and pathway overlap; the slope at the 89th.

> The 12W Myc profile *is* the 6W profile at 55% of the amplitude. The same pathways, in the
> same rank order, less of it. Nothing about which pathways Myc favours changes.

This is the powered, pathway-level form of Issue #4 ("magnitude not rank shrinks").

### 2.2 The Myc-specific part of the timeline is a uniform OFFSET

`Myc+time ~ WTtime`: slope **0.75**, intercept **-0.17**, R2 0.42.

- The **slope** (the component both genotypes share) is **not distinguishable from the
  null** (null median 0.63, 82nd percentile).
- The **intercept** — the Myc-specific offset — is at the **0th percentile** of the null
  (null median -0.05). It is also the batch-clean quantity: a batch offset common to both
  temporal contrasts cancels in the difference.

> Both genotypes' mitochondrial compartments move over the timeline in a pattern that is
> not specific to either. What Myc adds is a downward shift, near-constant across pathways,
> rather than a different set of pathways moving.

**FOUR PRECISION POINTS, added 2026-07-27 when the author asked for this panel to be explained
to a non-expert. None changes a number; each stops a reading the numbers do not support.**

1. **"A UNIFORM offset" is an approximation, not an identity.** It would be exact only at slope
   1. At 0.75 the Myc-specific change also carries a component proportional to how far the
   wild-type moved: `c_tp - c_tn = (slope - 1) * c_tn + intercept`. The defensible claim is
   that **the intercept is the extreme, batch-clean term**, not that the offset is literally
   constant. (Ranked, the per-pathway difference *is* a visibly flat band, so "near-constant"
   is fair; "uniform" is not.)
2. **The mirror error is worse: "the Myc change is NOT DUE TO the WT change" is too strong.**
   Slope 0.75 says much of the Myc gland's movement **is** shared with the wild-type's. The
   result is that the shared part is *unremarkable* (82nd percentile) while an **additional**
   downward component is not.
3. **WHAT THE NULL IS A NULL OF.** A gene-label shuffle preserving set size, expression and
   pathway overlap — a test against **the data's own internal structure**, not against sampling
   new animals. That is exactly why the slope's 82nd percentile is reported as a negative and
   the intercept's as a result. **500 draws**, so "0th percentile" means **p < 1/500 = 0.002
   and no smaller** — it is not p = 0. The bootstrap CI on the slope (0.57–0.94) resamples
   **pathways, not mice**; at n=6/group the animal-level question is separate and weaker.
4. **DO NOT FLIP THE AXES.** Regression is not symmetric: `c_tp ~ c_tn` gives **0.749 /
   -0.172**, `c_tn ~ c_tp` gives **0.559 / +0.118** (the geometric mirror would be 1.335, so
   the two are unrelated). Flipping also destroys the quantity of interest — the intercept
   would become "where the wild-type sits when the Myc gland does not move". The orientation is
   forced by the question: the wild-type gland is the substrate, hence the predictor.

**Direction, since it is the message and "bigger" gets it wrong.** Across the window the
wild-type gland drifts **UP** (median +0.041; 103 of 143 pathways rise) while the Myc gland
drifts **DOWN** (median -0.135; 22 of 143 rise). Per pathway `c_tp - c_tn` is negative in
**133 of 143** cases, median **-0.18**. Only 96 of 143 have a larger move in absolute size, so
"further down" is robust and "bigger" is not. That downward drift **is** the attenuation, not
Myc doing something extra. `figures/figureS2_controls.R` panel C plots this difference directly
(ranked) with the regression as an inset — and the difference must **never** be plotted against
either temporal contrast, per the artifact ledger in section 6.

### 2.3 The background is NOT moving toward the Myc state

Projection of the WT 6->12 vector onto the Myc(6W) vector = **+0.064** (1.0 would mean the
12W WT has arrived where 6W Myc+ was), against a null median of **0.24** — the **9.8th
percentile**. cos angle +0.21 overall (41st percentile, i.e. chance), but **-0.34 (content)
and -0.39 (priority) restricted to the biogenesis + OXPHOS arms**: in exactly the arms that
carry the claim, the WT background moves *away*.

Per tier (`geometry`), the "biogenesis+OXPHOS" aggregate is **not homogeneous** — it is two
arms moving away and one moving with:

| tier | cos(WT-time, Myc@6W) content | priority | cos(Myc+time, Myc@6W) content |
|---|---|---|---|
| Protein import / homeostasis | **-0.63** | **-0.77** | -0.98 |
| OXPHOS | **-0.65** | -0.41 | -0.95 |
| Mitochondrial central dogma | **+0.74** | +0.30 | -0.96 |
| Metabolism | +0.48 | +0.07 | -0.80 |
| Dynamics & surveillance | +0.24 | -0.37 | -0.91 |
| Signaling (n=7) | +0.40 | +0.74 | -0.33 |
| Small molecule transport (n=5) | +0.86 | +0.36 | +0.02 |

So the background moves *against* Myc in **import and OXPHOS**, and *with* it in central
dogma and metabolism. The last column is the striking one: the Myc+ temporal trajectory is a
near-exact reversal of the Myc effect in every substantial tier (-0.80 to -0.98).

By contrast `cos(Myc+ temporal, Myc@6W) = -0.84` (0th percentile of the null): the Myc+
trajectory over the timeline is close to an exact reversal of the Myc effect itself. The
attenuation is the Myc programme retracing its own steps, not the background catching up.

### 2.4 The WT-convergence term does not survive its control (the new result)

Issue #6 splits the attenuation into WT-convergence and Myc-fade after aligning each gene to
`d = sign(myc_6W)`. **`myc_6W = pos6 - neg6` and `timepoint_neg = neg12 - neg6` share the 6W
WT baseline with opposite signs**, so `cov = +var(neg6) > 0`: noise in that baseline
manufactures apparent convergence. Script 40 PART C splits the six 6W_neg mice — three anchor
the genotype contrast (and define `d`), three anchor the WT temporal contrast — over all 20
three-versus-three partitions.

**The comparison must be matched.** Splitting halves the baseline (3 mice, not 6), which adds
noise to `d`; sign noise pulls `frac_wt_toward` toward 0.5 *from both sides*, so comparing the
published 6-mouse value against a 3-mouse split confounds the artifact with plain
attenuation-to-chance. Each partition is therefore run **twice off the same half baseline A** —
once with A anchoring both contrasts, once with B anchoring the temporal one. Identical `d`,
identical universe, identical noise; the only difference is the sharing.

`frac_wt_toward` = fraction of genes whose WT temporal change points toward Myc (0.5 = chance):

| arm | published (6 mice) | half-shared | **half-split** | paired artifact (± sd) | expected sign, /20 |
|---|---|---|---|---|---|
| ALL genes (global) | 0.630 | 0.720 | **0.478** | **+0.242 ± 0.138** | 19 |
| Mammary luminal (pooled) | 0.606 | 0.700 | **0.466** | +0.234 ± 0.168 | 18 |
| Amino acid metabolism | 0.769 | 0.742 | **0.653** | +0.089 ± 0.180 | 13 |
| Lipid metabolism | 0.727 | 0.767 | **0.625** | +0.142 ± 0.202 | 13 |
| TCA cycle | 0.615 | 0.724 | **0.638** | +0.087 ± 0.547 | 11 |
| Nucleotide metabolism | 0.722 | 0.666 | **0.630** | +0.036 ± 0.385 | 10 |
| Mitoribosome | 0.509 | 0.568 | **0.516** | +0.052 | 10 |
| OXPHOS | 0.397 | 0.509 | **0.451** | +0.058 | 11 |
| OXPHOS subunits | 0.269 | 0.372 | **0.382** | **-0.010** | 10 |
| MYC targets (Hallmark) | 0.273 | 0.371 | **0.360** | **+0.011** | 11 |
| MYC signature (Felsher) | 0.286 | 0.377 | **0.359** | **+0.017** | 10 |
| Proliferation (pooled) | 0.411 | 0.528 | **0.391** | +0.137 | 13 |

**The artifact is large exactly where convergence was claimed and absent where divergence
was.** Globally, holding noise constant, sharing the baseline adds **+0.24** of apparent
convergence, in 19 of 20 partitions. On the MYC-target core and OXPHOS subunits the artifact
is **~0** — so their published divergence (0.27) is not affected by the sharing at all; the
0.27 -> 0.37 shift between the published and half-shared columns is the *noise* of halving
the baseline, not a correction.

Same picture as percentages (ratio of means, never a mean of ratios — `conv_pct` is a ratio of
two noisy quantities and explodes on individual splits): global **+61% -> -12%**; amino acid
45 -> **31**; lipid 47 -> **22**; TCA 37 -> **22**; nucleotide 31 -> **15**; OXPHOS subunits
-20 -> **-36**; MYC core -39 -> **-57**. `conv_pct_shared` (6-mouse, same quantifier)
reproduces `conv_pct_published` (42.2 vs 42.8 globally), so the pipeline is not the difference.

> **The headline "~34% of the attenuation is WT convergence" does not survive.** Globally the
> background does not move toward Myc at all. Convergence is real but **restricted to the
> biosynthetic/metabolic arms** (amino acid, lipid, TCA, nucleotide), and OXPHOS and the
> MYC-target core genuinely **diverge**. Issue #6's programme-specificity was the whole
> result; its global number was the artifact.

Ranges across the 20 partitions are wide for small arms (TCA n~15 spans 0.06-1.00), and the
paired artifact is only well resolved for the large sets (global sd 0.138 on +0.242; TCA sd
0.547 on +0.087). No p-value is offered against 0.5: genes within an arm are strongly
correlated, so a binomial SE would be anti-conservative, and the split-to-split SD measures
split variability, not sampling error.

---

## 3. The artifact ledger — statistics that are NOT evidence

Every naive test of "the background shapes Myc" is algebraically forced, because the
contrasts share terms. For each, the value the algebra forces given the observed marginals:

| statistic | observed | forced | mechanism |
|---|---|---|---|
| cor(WTtime, delta-Myc-effect), content | -0.185 | **-0.185** | delta = tpos - tneg; tneg enters negatively |
| cor(WTtime, delta-Myc-effect), priority | -0.327 | **-0.327** | same |
| cor(retention, WTtime), content | -0.522 | — | retention = 1 + (tpos - tneg)/m6 |
| cor(WT baseline @6W, Myc effect @6W), priority | -0.789 | — | mitoPPS is compositional (centred at 1); the effect subtracts the baseline |
| cor(WT baseline @12W, Myc effect @12W), priority | -0.885 | — | same |

**None of these may be cited.** They are carried in
`results/background_vs_myc.rds$artifact_ledger` so the next reader meets them with the
verdict attached.

---

## 4. The ruler tension: how big IS the attenuation?

Same gene set, four quantifiers of the genotype gap, `pct_retained` at 12W:

| set | share of nuclear transcriptome | set-avg LFC (unweighted) | set-avg LFC (expr-weighted) | mean gene-wise z |
|---|---|---|---|---|
| MITOCARTA_NUCLEAR_ENCODED | **83%** | 56% | 54% | 59% |
| MITOCARTA_OXPHOS_NU | **83%** | 63% | 57% | 65% |
| MITOCARTA_CENTRAL_DOGMA | 77% | 54% | 55% | 57% |
| PROTEIN_IMPORT_SORTING_HOMEOSTASIS | 83% | 56% | 63% | 57% |

Every ruler agrees on the **sign and the direction**; they disagree on the **size** (17-23%
lost by share, 37-46% by the LFC rulers). The difference is the denominator (per-sample
transcriptome sum vs median-of-ratios size factors) and the weighting. The manuscript should
state the **range**, not pick a ruler.

---

## 5. Verdict for the manuscript

The supportable claim is **not** "a changing background shapes Myc's mitochondrial effect".
It is the additive one:

> **Two components.** A reallocation of the mitochondrial compartment that both genotypes
> undergo over the timeline, and a Myc programme that decays proportionally (to ~55% of its
> 6W amplitude) without changing shape. They are close to independent: the background's own
> move is 57-62% as large as the Myc effect, points in a direction the Myc effect does not
> (projection 0.064, below an arbitrary matched panel), and in the biogenesis/OXPHOS arms
> points against it. Where the background does converge on the Myc state — the biosynthetic
> arms — it survives the sample-split control; globally it does not.

## 6. What is NOT separable

- **BATCH = TIMEPOINT.** The batch offset cancels in the interaction (the panel C intercept,
  the whole of 2.1 and 2.2's Myc-specific term) but **not** in the convergence/fade split.
  Calling the shared temporal vector "developmental" rather than "batch" is an
  interpretation, not a measurement. The sample-split removes the shared-baseline artifact;
  it cannot remove the batch confound.
- **Composition.** The Myc-fade term is per-cell weakening OR dilution of a shrinking
  Myc-responsive compartment; script 31 PART C (deconvolution) remains deferred.
- n=6/group, exploratory. This is a decomposition with nulls, not confirmatory inference.

## 7. Downstream consequences

- Script 31's `decomp_conv` **stays on disk unchanged**; script 40 tests it and reports the
  correction. Rewriting 31 is a separate decision.
- Any sentence quoting "34% WT-convergence / 66% Myc-fade" as a global statement must be
  rewritten to the arm-specific form above.
- `docs/2026-07-13_..._walkthrough...md`, `docs/2026-07-12_BlockA_revision_synthesis...md`
  and `docs/2026-07-18_narrative_synthesis_five_questions.md` all carry the global number;
  they are historical records and are not edited, but this note supersedes them.
