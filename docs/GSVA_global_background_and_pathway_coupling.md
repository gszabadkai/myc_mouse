# Interpreting and adjusting the global background in GSVA pathway coupling analysis

## Study design

- **Samples:** 24
- **Groups:** 4 groups defined by:
  - timepoint: 6W versus 12W
  - genotype: WT versus MYC+
- **Gene sets:** 884
- **Important prior check:** the gene sets were tested for redundancy, and redundancy is minimal.

Two observations motivate the analysis:

- Random pairs of gene-set scores correlate at approximately **0.678**.
- The `mito_oxphos` score correlates with the per-sample mean of all 884 GSVA scores at **0.979**.

After removing the global per-sample mean, the apparent correlation ceiling falls to approximately **0.287**.

The central problem is therefore not ordinary gene-set overlap. It is the presence of a strong sample-level axis that moves many largely non-redundant gene-set scores together.

---

## 1. What is the background?

The most useful statistical description is a **common-mode latent factor** in the GSVA score matrix.

For sample \(s\) and gene set \(k\):

\[
S_{sk} = \mu_k + X_s\beta_k + \lambda_k G_s + \varepsilon_{sk}
\]

where:

- \(S_{sk}\) is the GSVA score;
- \(\mu_k\) is the average score of gene set \(k\);
- \(X_s\) contains known experimental variables, particularly timepoint, genotype and their interaction;
- \(G_s\) is a global sample-level axis;
- \(\lambda_k\) is the loading of gene set \(k\) on that axis;
- \(\varepsilon_{sk}\) is pathway-specific variation not explained by the design or the global axis.

For two gene sets \(k\) and \(l\), the observed covariance includes:

\[
\lambda_k\lambda_l\operatorname{Var}(G)
\]

even if their pathway-specific residual activities are unrelated.

Thus, the high pairwise correlations are consistent with a dominant, approximately rank-one covariance component shared across many gene sets.

Because redundancy between the 884 gene sets is minimal, this common covariance is unlikely to be explained primarily by duplicated biological annotations or extensive gene overlap. The more likely explanation is that the samples differ along a broad transcriptomic axis that influences many distinct biological programmes simultaneously.

---

## 2. What might generate this global axis?

The global axis is likely to be a mixture of biological and technical sources.

### 2.1 Experimental group structure

Timepoint, MYC status and their interaction may affect a large fraction of the transcriptome. If many pathways shift in the same direction between groups, correlations calculated across all 24 samples will be high even when there is no additional coupling between the pathways.

For example, two pathways can correlate strongly simply because both are high in 6W MYC+ samples and low in 12W WT samples.

This is coordinated differential pathway activity, but it is not necessarily evidence that the pathways are coupled between individual mice.

### 2.2 Broad biological state

The latent factor may represent a genuine whole-sample biological state, such as:

- proliferation or differentiation;
- MYC-associated transcriptional state;
- inflammatory or immediate-early-gene activation;
- stromal, immune or epithelial composition;
- mitochondrial content;
- developmental progression;
- tissue integrity or stress.

The very high correlation between `mito_oxphos` and the global mean indicates that `mito_oxphos` is acting primarily as a readout of this broad state in the present dataset.

That does not mean the `mito_oxphos` score is invalid. It means that most of its between-sample variation is shared with many other pathway scores.

### 2.3 Technical or compositional variation

Residual technical effects can also contribute:

- RNA quality;
- batch;
- sequencing depth or library complexity not fully handled upstream;
- differences in the fraction of highly expressed gene classes;
- sample composition;
- normalization-related compositional effects.

These effects can be real properties of the measured libraries without representing the pathway-specific biology of interest.

### 2.4 Cohort-dependent scoring

GSVA scores are estimated relative to the expression distributions represented by the analysed cohort. Broad sample-to-sample expression differences can therefore propagate into many pathway scores.

However, the background should not be described simply as an additive GSVA offset. A perfectly uniform increase in every transcript would not necessarily raise all competitive enrichment scores equally.

The empirical finding is better described as a **broad relative-expression programme captured by many GSVA scores**.

---

## 3. Is this MYC-driven global RNA amplification?

Possibly in part, but this cannot be concluded from standard normalized bulk RNA-seq alone.

Conventional library-size normalization measures relative transcript abundance. A uniform increase in total RNA per cell is generally not identifiable without an external absolute reference, such as:

- spike-ins added per cell or per fixed number of cells;
- direct measurements of RNA per cell;
- suitable absolute transcript quantification;
- carefully controlled cell-number-normalized assays.

Therefore, the latent factor may be associated with MYC and may reflect widespread MYC-dependent transcriptional changes, but it should initially be called something neutral, such as:

> global GSVA factor, common pathway-score axis, broad transcriptomic state, or latent sample-level pathway activity.

Its biological interpretation should be based on associations with known covariates and gene-level data.

---

## 4. Why subtracting the per-sample mean is informative but not ideal

For each sample, subtracting the mean score across all 884 gene sets removes the dominant common shift:

\[
S'_{sk} = S_{sk} - \overline{S}_{s}
\]

The collapse of the correlation ceiling from approximately 0.678 to 0.287 is strong evidence that a sample-level common component generates much of the raw correlation structure.

However, mean subtraction assumes that every gene set loads equally on the global factor:

\[
\lambda_1 = \lambda_2 = \cdots = \lambda_{884}
\]

This is unlikely to be exactly true.

It also forces the average adjusted pathway score within every sample to zero. This creates a compositional constraint and can introduce negative correlations: when some adjusted scores rise, others must fall relative to the sample mean.

Because gene-set redundancy is minimal, the per-sample mean is more interpretable here than it would be in a highly redundant collection. It is not being dominated by repeated copies of the same biological programme. Nevertheless, it remains an equal-weight summary of 884 variables, so regression on an empirically estimated factor is preferable for formal inference.

Mean-centred scores are therefore useful for:

- visualisation;
- diagnosing the magnitude of the common mode;
- sensitivity analysis;
- an intuitive first-pass estimate of specific pathway behaviour.

They should not be the only basis of formal coupling claims.

---

## 5. Separate three different questions

### 5.1 Total coordination

\[
r_{\mathrm{raw}} =
\operatorname{cor}(S_k,S_l)
\]

This asks whether two gene-set scores move together across the complete experiment.

It includes:

- timepoint effects;
- genotype effects;
- interaction effects;
- the global transcriptomic axis;
- pathway-specific coupling;
- residual technical variation.

### 5.2 Coordination beyond the experimental group means

\[
r_{\mathrm{design}} =
\operatorname{cor}
\left(
S_k \mid \text{timepoint}*\text{genotype},
S_l \mid \text{timepoint}*\text{genotype}
\right)
\]

This asks whether mice with unusually high activity of pathway \(k\), relative to other mice in the same experimental design, also have unusually high activity of pathway \(l\).

This is closer to within-condition biological coupling.

### 5.3 Specific coupling beyond the global sample state

\[
r_{\mathrm{specific}} =
\operatorname{cor}(\varepsilon_k,\varepsilon_l)
\]

where both scores have been adjusted for:

1. timepoint;
2. genotype;
3. timepoint-by-genotype interaction;
4. the global latent factor;
5. justified technical covariates.

This asks whether the pathways remain associated after removing broad coordinated transcriptomic variation.

All three quantities can be biologically meaningful, but they answer different questions.

---

## 6. Recommended analysis workflow

## Step 1: model the known experimental design

Create a sample-by-gene-set matrix and remove the fixed design effects.

```r
# gsva_scores:
#   rows = gene sets
#   columns = samples
#
# meta:
#   one row per sample
#   row names match columns of gsva_scores

Y <- t(gsva_scores)  # samples x gene sets

meta <- meta[colnames(gsva_scores), , drop = FALSE]

X <- model.matrix(
  ~ timepoint * genotype,
  data = meta
)

# Residual GSVA scores after removing the known design
R_design <- qr.resid(qr(X), Y)
```

Include technical covariates only when they are independently justified, for example:

```r
X <- model.matrix(
  ~ batch + RIN + timepoint * genotype,
  data = meta
)
```

Avoid automatically including variables that could be biological mediators of MYC or developmental effects.

---

## Step 2: inspect the global factor before removing it

Estimate the first principal component of the design-residualised GSVA matrix.

```r
pc <- prcomp(
  R_design,
  center = FALSE,
  scale. = FALSE
)

global_factor <- pc$x[, 1]
global_loadings <- pc$rotation[, 1]

summary(pc)
```

Examine:

- variance explained by PC1;
- correlation between PC1 and the mean GSVA score;
- correlation between PC1 and `mito_oxphos`;
- the sign and distribution of pathway loadings;
- association of PC1 with timepoint, genotype and batch;
- association with RNA quality and library metrics;
- association with cell-composition estimates;
- association with selected gene-level signatures.

Example:

```r
global_mean <- rowMeans(Y)

cor(global_factor, global_mean)
cor(global_factor, Y[, "mito_oxphos"])

plot(global_mean, global_factor)
plot(Y[, "mito_oxphos"], global_factor)
```

If most gene sets load in the same direction, PC1 is functioning as a common score-level factor.

If gene sets split into strong positive and negative loading groups, PC1 is instead a broad biological contrast between programmes. It can still be adjusted for, but its interpretation should be explicit.

---

## Step 3: regress each pathway on the design and global factor

```r
X_global <- cbind(
  X,
  global_factor = global_factor
)

R_specific <- qr.resid(
  qr(X_global),
  Y
)
```

The columns of `R_specific` contain pathway scores adjusted for both the experimental design and the common latent factor.

Pairwise coupling can then be calculated as:

```r
r_raw <- cor(
  Y[, "mito_oxphos"],
  Y[, "myc_signature"]
)

r_design <- cor(
  R_design[, "mito_oxphos"],
  R_design[, "myc_signature"]
)

r_specific <- cor(
  R_specific[, "mito_oxphos"],
  R_specific[, "myc_signature"]
)

c(
  raw = r_raw,
  design_adjusted = r_design,
  global_adjusted = r_specific
)
```

---

## Step 4: reduce circularity when testing a selected pair

If the global factor is estimated from all 884 pathways, the two pathways being tested contribute slightly to the factor used to adjust them.

With 884 largely non-redundant sets, this influence is probably small for most individual pairs. Nevertheless, a leave-pair-out procedure is preferable for selected biological hypotheses.

```r
pairwise_global_adjusted_correlation <- function(
    set_a,
    set_b,
    Y,
    X,
    R_design
) {
  background_sets <- setdiff(
    colnames(R_design),
    c(set_a, set_b)
  )

  pc_pair <- prcomp(
    R_design[, background_sets, drop = FALSE],
    center = FALSE,
    scale. = FALSE
  )

  G_pair <- pc_pair$x[, 1]

  X_pair <- cbind(
    X,
    global_factor = G_pair
  )

  E_pair <- qr.resid(
    qr(X_pair),
    Y[, c(set_a, set_b), drop = FALSE]
  )

  c(
    raw = cor(Y[, set_a], Y[, set_b]),
    design_adjusted =
      cor(R_design[, set_a], R_design[, set_b]),
    global_adjusted =
      cor(E_pair[, 1], E_pair[, 2])
  )
}
```

Example:

```r
pairwise_global_adjusted_correlation(
  set_a = "mito_oxphos",
  set_b = "myc_signature",
  Y = Y,
  X = X,
  R_design = R_design
)
```

For a small set of pre-specified mito–MYC comparisons, this is a defensible primary analysis.

---

## Step 5: compare mean subtraction with PC1 adjustment

The per-sample mean and PC1 may be nearly identical in this dataset. This should be demonstrated rather than assumed.

```r
Y_mean_centered <- sweep(
  Y,
  MARGIN = 1,
  STATS = rowMeans(Y),
  FUN = "-"
)

cor(
  rowMeans(Y),
  global_factor
)
```

For each selected pathway pair, compare:

- raw correlation;
- design-adjusted correlation;
- sample-mean-adjusted correlation;
- PC1-adjusted correlation;
- leave-pair-out-PC1-adjusted correlation.

Agreement between the approaches would show that the result is not an artefact of one particular correction.

---

## 7. Should the global factor be estimated before or after removing the design?

The preferred primary sequence is:

\[
\boxed{
\text{known experimental design}
\rightarrow
\text{global latent factor}
\rightarrow
\text{specific coupling}
}
\]

Estimating PC1 from the unadjusted score matrix risks making the global factor equivalent to the four-group experimental contrast.

Residualising timepoint, genotype and their interaction first defines the latent factor as the broad variation that remains between mice beyond the expected group means.

However, an unadjusted PC1 remains useful descriptively. Comparing adjusted and unadjusted PC1 can reveal how much of the global axis is explained by the experimental design.

---

## 8. How to interpret `mito_oxphos`

Given:

\[
\operatorname{cor}
(
\texttt{mito\_oxphos},
\text{mean of 884 scores}
)
= 0.979
\]

the `mito_oxphos` score is currently almost indistinguishable from the global score axis.

This means that its raw pairwise correlations mainly answer:

> Does the second pathway track the broad sample-level transcriptomic state that is also captured by `mito_oxphos`?

They do not initially answer:

> Is OXPHOS specifically coupled to the second pathway independently of the broader state?

The latter requires global-factor-adjusted residuals.

This distinction is particularly important for MYC biology. MYC may simultaneously affect proliferation, ribosome biogenesis, metabolism, stress, differentiation and tissue composition. A raw OXPHOS–MYC correlation can therefore represent participation in a common MYC-associated state rather than a direct or selective OXPHOS–MYC relationship.

---

## 9. How to report the results

For each biologically important pair, report at least:

| Measure | Interpretation |
|---|---|
| Raw correlation | Overall coordination across all samples |
| Design-adjusted correlation | Coupling beyond timepoint/genotype group means |
| Global-adjusted correlation | Specific coupling beyond the common sample state |
| PC1 loading of each pathway | Participation of each pathway in the common programme |
| Bootstrap interval | Uncertainty with \(n=24\) |

Example interpretation:

### Pattern A: global-only association

- raw \(r = 0.82\)
- design-adjusted \(r = 0.65\)
- global-adjusted \(r = 0.06\)

Interpretation:

> The pathways strongly track the same broad sample-level state, but there is little evidence for specific coupling after the common factor is removed.

### Pattern B: global plus specific association

- raw \(r = 0.82\)
- design-adjusted \(r = 0.68\)
- global-adjusted \(r = 0.48\)

Interpretation:

> Part of the association reflects the global state, but substantial pathway-specific coupling remains.

### Pattern C: masked antagonism

- raw \(r = 0.10\)
- design-adjusted \(r = -0.15\)
- global-adjusted \(r = -0.52\)

Interpretation:

> A shared positive global component masked an underlying negative pathway-specific relationship.

---

## 10. Statistical limitations with 24 samples

There are:

\[
\frac{884 \times 883}{2}
= 390{,}286
\]

possible gene-set pairs.

With only 24 samples, it is not realistic to infer a stable unrestricted network of all pairwise couplings.

Even before fitting covariates, the sampling uncertainty of a correlation is substantial. A residual correlation near 0.287 can easily arise under the null for an individual pair.

The adjusted analysis also consumes degrees of freedom:

- intercept;
- timepoint;
- genotype;
- interaction;
- global factor;
- any technical covariates.

Therefore, the most reliable strategy is:

1. define a limited set of biologically motivated pathway pairs;
2. calculate raw, design-adjusted and global-adjusted correlations;
3. quantify uncertainty by bootstrap;
4. treat the full 884-by-884 correlation matrix as exploratory;
5. avoid estimating correlations separately within each \(n=6\) group, except descriptively.

---

## 11. Bootstrap recommendation

Resample mice within the four experimental groups so that the study design is preserved.

For every bootstrap iteration:

1. resample samples within each timepoint–genotype group;
2. refit the design model;
3. recalculate the global factor;
4. residualise the selected pathway pair;
5. recalculate the adjusted correlation.

The global factor must be re-estimated inside every bootstrap iteration. Otherwise, uncertainty introduced by estimating the factor is ignored.

A percentile or bias-corrected bootstrap interval can then be reported.

---

## 12. A better null model for pairwise coupling

Since gene-set redundancy is minimal, overlap-driven correlation is not the primary concern. Nevertheless, random pairs of pathways are not automatically a clean null because distinct pathways can still respond to the same global biological programmes.

For an empirical background distribution, compare an observed pathway pair with random pairs matched on:

- gene-set score variance;
- global-factor loading;
- proportion of variance explained by the experimental design;
- gene-set size;
- mean score, when relevant.

Matching on global-factor loading is especially important. Two pathways with large positive loadings are expected to have high raw correlation even without specific coupling.

For adjusted correlations, the empirical null should be generated after applying the complete design and global-factor correction.

---

## 13. Alternative factor models

PC1 is the simplest and probably most transparent adjustment. It is appropriate when the score matrix is dominated by one common axis.

Check the scree plot:

```r
plot(
  pc$sdev^2 / sum(pc$sdev^2),
  type = "b",
  xlab = "Principal component",
  ylab = "Proportion of variance explained"
)
```

If there is one dominant component, a one-factor model is well justified.

If several components are substantial, possible approaches include:

- adjusting for the first two or three PCs;
- surrogate variable analysis;
- factor analysis;
- probabilistic PCA;
- mixed or latent-factor models.

However, with only 24 samples, removing several components can rapidly remove meaningful biology and destabilise correlations.

The number of factors should therefore be pre-specified using diagnostics and sensitivity analysis, not selected separately for every pathway pair.

A reasonable primary analysis is one global factor, with two-factor adjustment as a sensitivity analysis only if PC2 is clearly substantial and interpretable.

---

## 14. Avoid overcorrection

The global factor may itself be central to the biology. Removing it answers a narrower question but does not make it irrelevant.

For example, if MYC drives a broad coordinated state that includes OXPHOS, proliferation and biosynthesis, then the common factor may be the strongest biological result.

Therefore, preserve two parallel interpretations:

### Global programme biology

Which pathways load most strongly on the common factor?

What sample characteristics predict the factor?

Does it distinguish timepoints or genotypes?

Does it associate with mitochondrial abundance, proliferation, developmental state or cell composition?

### Pathway-specific coupling

After accounting for the global programme, which pathway pairs retain a reproducible association?

The analysis should not replace the raw score matrix with residuals for every purpose. Residuals are specifically for testing selective coupling beyond the global programme.

---

## 15. Recommended primary analysis for this dataset

1. Fit each GSVA score to:

\[
\text{score}
\sim
\text{timepoint} * \text{genotype}
\]

2. Estimate PC1 from the 884 design-residualised scores.

3. Characterise PC1 using:
   - explained variance;
   - pathway loadings;
   - correlation with the per-sample score mean;
   - correlation with `mito_oxphos`;
   - associations with sample and QC variables.

4. For each pre-specified mito–MYC or mito–developmental pair, calculate:
   - raw correlation;
   - design-adjusted correlation;
   - PC1-adjusted correlation;
   - leave-pair-out-PC1-adjusted correlation.

5. Bootstrap the entire adjustment procedure within experimental groups.

6. Use mean-subtracted scores only as a transparent sensitivity analysis.

7. Interpret:
   - raw correlation as shared global coordination;
   - residual correlation as specific coupling beyond the global state.

8. Do not treat all 390,286 pairwise correlations as confirmatory. Use the complete matrix for exploration and hypothesis generation only.

---

## 16. Concise conclusion

The background is best understood as a strong **sample-level common pathway-activity factor**. Because the gene sets have minimal redundancy, it is unlikely to arise mainly from duplicated or overlapping pathway definitions. It more plausibly represents broad biological and technical variation that simultaneously affects many distinct pathways.

The very high correlation of `mito_oxphos` with the mean GSVA score shows that raw `mito_oxphos` variation is dominated by this common factor. Raw pairwise correlations therefore measure participation in a shared whole-sample programme more than selective pathway coupling.

The recommended analysis is:

\[
\boxed{
\text{remove timepoint/genotype effects}
\rightarrow
\text{estimate the common latent factor}
\rightarrow
\text{test correlations between residual pathway activities}
}
\]

The global factor should also be analysed as a biological phenotype in its own right, rather than treated only as unwanted noise.
