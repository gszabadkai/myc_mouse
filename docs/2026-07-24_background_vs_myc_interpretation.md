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
> not specific to either. What Myc adds is a uniform downward offset, not a different set
> of pathways moving.

### 2.3 The background is NOT moving toward the Myc state

Projection of the WT 6->12 vector onto the Myc(6W) vector = **+0.064** (1.0 would mean the
12W WT has arrived where 6W Myc+ was), against a null median of **0.24** — the **9.5th
percentile**. cos angle +0.21 overall (36th percentile, i.e. chance), but **-0.34 (content)
and -0.39 (priority) restricted to the biogenesis + OXPHOS arms**: in exactly the arms that
carry the claim, the WT background moves *away*.

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

`frac_wt_toward` = fraction of genes whose WT temporal change points toward Myc (0.5 = chance):

| arm | published | shared baseline (same quantifier) | **split baseline** | verdict |
|---|---|---|---|---|
| ALL genes (global) | 0.630 | 0.630 | **0.478** | collapses to chance |
| Mammary luminal (pooled) | 0.606 | 0.610 | **0.466** | collapses to chance |
| Amino acid metabolism | 0.769 | 0.718 | **0.653** | converges (survives) |
| Lipid metabolism | 0.727 | 0.738 | **0.625** | converges (survives) |
| TCA cycle | 0.615 | 0.583 | **0.638** | converges (survives) |
| Nucleotide metabolism | 0.722 | 0.647 | **0.630** | converges (survives) |
| Mitoribosome | 0.509 | 0.545 | **0.516** | chance, as published |
| OXPHOS | 0.397 | 0.478 | **0.451** | diverges |
| OXPHOS subunits | 0.269 | 0.400 | **0.382** | diverges (survives) |
| MYC targets (Hallmark) | 0.273 | 0.265 | **0.360** | diverges (survives) |
| MYC signature (Felsher) | 0.286 | 0.283 | **0.359** | diverges (survives) |
| Proliferation (pooled) | 0.411 | 0.366 | **0.391** | diverges |

Same picture as percentages (ratio of means, never a mean of ratios — `conv_pct` is a ratio
of two noisy quantities and explodes on individual splits): global **+42% -> -12%**; amino
acid 39 -> **31**; lipid 38 -> **22**; TCA 32 -> **22**; nucleotide 26 -> **15**; OXPHOS
subunits -25 -> **-36**; MYC core -50 -> **-57**.

`conv_pct_shared` reproduces `conv_pct_published` (42.2 vs 42.8 globally), so the split is
the only thing that changes.

> **The headline "~34% of the attenuation is WT convergence" does not survive.** Globally the
> background does not move toward Myc at all. Convergence is real but **restricted to the
> biosynthetic/metabolic arms** (amino acid, lipid, TCA, nucleotide), and OXPHOS and the
> MYC-target core genuinely **diverge**. Issue #6's programme-specificity was the whole
> result; its global number was the artifact.

Ranges across the 20 splits are wide for small arms (TCA n~15 spans 0.06-1.00). No p-value is
offered against 0.5: genes within an arm are strongly correlated, so a binomial SE would be
anti-conservative, and the split-to-split SD measures split variability, not sampling error.

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
