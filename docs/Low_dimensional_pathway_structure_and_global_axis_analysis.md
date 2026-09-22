# Understanding low-dimensional pathway structure and the global axis

## Study context

The dataset contains:

- **24 samples** in four groups: 6W WT, 6W MYC+, 12W WT and 12W MYC+
- **884 gene-set scores**
- minimal gene-set redundancy
- a strong common sample-level component in the pathway-score matrix

The main observations are:

- apparently unrelated gene sets correlate strongly;
- `mito_oxphos` correlates almost perfectly with the mean score across all 884 sets;
- removal of one global sample-level axis eliminates most of the pairwise correlation.

This document explains what **low-dimensional** means in this setting and how to analyse the global axis without discarding potentially important MYC biology.

---

# 1. What low-dimensional means

## 1.1 Intuitive definition

Although 884 pathway scores were measured, the samples may not vary independently along 884 separate biological directions.

Instead, many pathways may rise and fall together because they respond to one or a few shared forces, such as:

- MYC activity;
- developmental stage;
- proliferation or differentiation;
- cell-type composition;
- tissue stress;
- mitochondrial abundance;
- RNA quality or library properties.

If most pathway variation is explained by only one, two or three such axes, the dataset is **low-dimensional**. The number of measured variables is large, but the number of underlying independent patterns is small.

## 1.2 Simple statistical representation

For sample \(s\) and pathway \(k\):

\[
S_{sk} = \mu_k + \lambda_k G_s + \varepsilon_{sk}
\]

where:

- \(S_{sk}\) is the pathway score;
- \(\mu_k\) is the mean score of pathway \(k\);
- \(G_s\) is a shared sample-level factor;
- \(\lambda_k\) is the loading of pathway \(k\) on that factor;
- \(\varepsilon_{sk}\) is pathway-specific variation.

If one factor \(G\) explains most of the covariance among the 884 scores, the data are approximately one-dimensional. If two or three factors are needed, the data are still low-dimensional relative to 884 measured pathways.

## 1.3 Geometrical interpretation

Each mouse can be represented as a point in an 884-dimensional space:

\[
(\text{score}_1,\text{score}_2,\ldots,\text{score}_{884})
\]

In a genuinely high-dimensional dataset, samples would vary along many independent directions. In a low-dimensional dataset, the 24 sample points lie close to:

- a line, if one factor dominates;
- a plane, if two factors dominate;
- a small subspace, if a few factors dominate.

Thus, the samples occupy only a small part of the nominal 884-dimensional space.

## 1.4 Biological interpretation

Low-dimensionality does **not** mean that only one biological pathway is active.

It means that many pathways change together and therefore cannot be statistically separated into independent programmes in this dataset. MYC may simultaneously influence ribosome biogenesis, proliferation, nucleotide synthesis, mitochondrial metabolism, translation and stress responses. These pathways are biologically distinct, but if they move together across the 24 mice they form one broad transcriptomic axis.

This is different from gene-set redundancy:

- **gene-set redundancy** asks whether the sets contain similar genes or annotations;
- **low-dimensional sample structure** asks whether distinct pathways respond together across samples.

Minimal redundancy among the 884 sets therefore does not prevent strong low-dimensional structure.

## 1.5 Why the current observations suggest low-dimensionality

The observations are consistent with a dominant common factor:

- random pathway pairs correlate at approximately 0.678;
- `mito_oxphos` correlates with the mean score across all pathways at approximately 0.979;
- removal of one global axis lowers the apparent correlation ceiling to approximately 0.287.

This implies that much of the covariance between two pathways can be approximated as:

\[
\operatorname{Cov}(S_k,S_l)
\approx
\lambda_k\lambda_l\operatorname{Var}(G)
\]

The remaining pathway-specific covariance is much smaller.

In practical terms, many pathway pairs appear correlated because they track the same broad state, not because each pair has a unique direct relationship.

## 1.6 Quantifying low-dimensionality

### Principal-component variance

Calculate the proportion of variance explained by PC1, PC2 and later components.

```r
pc <- prcomp(
  pathway_scores,
  center = TRUE,
  scale. = FALSE
)

variance_explained <- pc$sdev^2 / sum(pc$sdev^2)
variance_explained[1:10]

plot(
  variance_explained,
  type = "b",
  xlab = "Principal component",
  ylab = "Proportion of variance explained"
)
```

A dominant PC1 indicates a strong one-dimensional component.

### Effective dimensionality

An approximate effective dimension can be calculated from the covariance eigenvalues:

\[
D_{\mathrm{eff}}
=
\frac{\left(\sum_i \lambda_i\right)^2}
{\sum_i \lambda_i^2}
\]

```r
eig <- pc$sdev^2

effective_dimension <- sum(eig)^2 / sum(eig^2)
effective_dimension
```

With 24 samples, the centred matrix can have at most 23 non-zero principal components, regardless of the 884 measured pathways. In practice, it may be dominated by far fewer.

---

# 2. Why simply removing the global axis can be misleading

## 2.1 The global axis can contain both biology and technical variation

A more realistic model is:

\[
S_{sk}
=
\mu_k
+
X_s\beta_k
+
\lambda_kG^{\mathrm{bio}}_s
+
\theta_kG^{\mathrm{tech}}_s
+
\varepsilon_{sk}
\]

where:

- \(X_s\) contains timepoint, genotype and their interaction;
- \(G^{\mathrm{bio}}\) is a broad biological programme;
- \(G^{\mathrm{tech}}\) is unwanted technical variation;
- \(\varepsilon_{sk}\) is pathway-specific activity.

For two pathways:

\[
\operatorname{Cov}(S_k,S_l)
=
\lambda_k\lambda_l\operatorname{Var}(G^{\mathrm{bio}})
+
\theta_k\theta_l\operatorname{Var}(G^{\mathrm{tech}})
+
\operatorname{Cov}(\varepsilon_k,\varepsilon_l)
\]

Removing one undifferentiated PC1 may remove all these components together. That can be too aggressive if the broad biological programme is central to the MYC phenotype.

## 2.2 A disappearing residual correlation does not make the raw association false

Suppose:

\[
\mathrm{MYC}
\rightarrow
G
\rightarrow
\begin{cases}
\mathrm{OXPHOS}\\
\mathrm{Pathway\ B}
\end{cases}
\]

MYC induces a broad state \(G\), and both pathways respond to it. Conditioning on \(G\) removes the OXPHOS–pathway B correlation.

This means that the association is **mediated through the shared global state**, rather than representing an additional pathway-specific relationship. It does not mean that the original relationship was artefactual.

The analysis should distinguish:

1. total coordination;
2. coordination mediated through the global state;
3. residual pathway-specific coupling.

---

# 3. Ways to get behind the global axis

## 3.1 Separate technical and biological variation using anchors

The best way to distinguish technical from biological global variation is to use variables external to the pathway-score matrix.

Potential technical anchors include:

- batch and processing date;
- sequencing lane;
- RNA integrity;
- mapping rate;
- duplication rate;
- insert size;
- library complexity;
- total read count;
- percentage of mitochondrial reads;
- negative-control genes;
- technical replicates.

Possible approaches include:

- direct regression on known QC covariates;
- supervised surrogate-variable analysis;
- RUV methods using negative-control genes;
- factor estimation using technical replicates.

The preferred sequence is:

\[
\boxed{
\text{remove anchored technical variation}
\rightarrow
\text{estimate the remaining biological axis}
}
\]

This is preferable to treating the entire PC1 as unwanted background.

### Identifiability limitation

If technical effects are completely confounded with genotype or timepoint, they cannot be separated reliably from biology using the expression data alone. For example, if every MYC+ sample was processed in one batch and every WT sample in another, no statistical correction can unambiguously identify which component is MYC and which is batch.

## 3.2 Decompose the global factor into design-driven and within-group components

After technical adjustment, estimate the remaining global factor \(G\), then fit:

\[
G
\sim
\text{timepoint}
*
\text{genotype}
+
\text{technical covariates}
\]

Decompose:

\[
G = G_{\mathrm{design}} + G_{\mathrm{within}}
\]

where:

- \(G_{\mathrm{design}}\) is the portion predicted by timepoint, genotype and interaction;
- \(G_{\mathrm{within}}\) is the remaining between-mouse variation.

```r
global_model <- lm(
  global_factor ~ timepoint * genotype + batch + RIN,
  data = meta
)

global_design <- fitted(global_model)
global_within <- residuals(global_model)
```

### Interpretation

`global_design` captures broad coordinated response to:

- developmental age;
- MYC status;
- differential MYC effects at 6W and 12W.

`global_within` may reflect:

- heterogeneity in MYC activity;
- variable cell composition;
- mitochondrial abundance;
- inflammation or stress;
- mouse-to-mouse biological variability;
- residual technical effects.

These components should be analysed separately.

## 3.3 Analyse pathway loadings on the global axis

If most pathways participate in one broad programme, the loading on that programme becomes an important phenotype.

For pathway \(k\):

\[
S_k = \alpha_k + \lambda_k G + \varepsilon_k
\]

The loading \(\lambda_k\) describes how strongly pathway \(k\) follows the global state.

Important questions include:

- Is OXPHOS one of the strongest-loading pathways?
- Are mitochondrial pathways enriched among high-loading pathways?
- Are MYC-target pathways disproportionately aligned with the factor?
- Does pathway loading change between 6W and 12W?
- Does MYC strengthen or weaken alignment between OXPHOS and the global programme?

A condition-dependent loading model is:

\[
S_{\mathrm{OXPHOS}}
=
\alpha
+
\lambda G
+
\gamma_1(G\times\mathrm{MYC})
+
\gamma_2(G\times\mathrm{time})
+
\gamma_3(G\times\mathrm{MYC}\times\mathrm{time})
+
\varepsilon
\]

```r
fit_loading <- lm(
  mito_oxphos ~
    global_factor * genotype * timepoint +
    batch + RIN,
  data = dat
)

summary(fit_loading)
anova(fit_loading)
```

This can be more informative than asking only whether OXPHOS retains a residual correlation after its dominant biological signal has been removed.

## 3.4 Test differential coupling across conditions

A pooled correlation assumes that coupling is the same in all four groups. It may instead:

- exist only at 6W;
- exist only in MYC+;
- weaken between 6W and 12W;
- reverse direction;
- be induced or lost by MYC.

Four separate correlations with \(n=6\) per group are unstable. A better strategy is an interaction model using all 24 samples:

\[
B
=
\alpha
+
\beta A
+
\gamma_1\mathrm{MYC}
+
\gamma_2\mathrm{time}
+
\gamma_3(\mathrm{MYC}\times\mathrm{time})
+
\delta_1(A\times\mathrm{MYC})
+
\delta_2(A\times\mathrm{time})
+
\delta_3(A\times\mathrm{MYC}\times\mathrm{time})
+
\varepsilon
\]

The \(\delta\) terms test whether the relationship between pathways A and B changes by genotype or timepoint.

```r
fit_coupling <- lm(
  pathway_B ~
    pathway_A * genotype * timepoint +
    technical_factor,
  data = dat
)

summary(fit_coupling)
anova(fit_coupling)
```

Because the sample size is small, this should be limited to a small number of pre-specified pathway pairs.

## 3.5 Look below the whole-pathway average

An aggregate pathway score can hide biologically specific substructure.

A pathway may contain:

- one subgroup specifically coupled to OXPHOS;
- another subgroup responding only to the global MYC programme;
- a third subgroup moving in the opposite direction.

The whole-set score may show no residual coupling even when a coherent submodule is present.

For selected pathway pairs, examine the cross-set gene-correlation matrix:

\[
C_{gh}
=
\operatorname{cor}(Z_g,Z_h),
\qquad
g\in A,\ h\in B
\]

Useful exploratory approaches include:

- average cross-set gene correlation;
- hierarchical clustering;
- biclustering;
- sparse partial least squares;
- sparse canonical correlation;
- predefined functional submodules;
- leading-edge or core genes;
- leave-one-sample-out stability.

With 24 samples, learned submodules should be considered exploratory unless validated in an independent dataset.

## 3.6 Use a bifactor or sparse-plus-dense model

The pathway-score matrix may contain:

- one dense general factor affecting most pathways;
- several smaller local factors affecting selected pathway families.

A conceptual model is:

\[
S_{sk}
=
\lambda_kG_s
+
\sum_{j=1}^{J}\phi_{kj}L_{sj}
+
\varepsilon_{sk}
\]

where:

- \(G_s\) is the general global factor;
- \(L_{sj}\) are local factors;
- \(\phi_{kj}\) are sparse local loadings.

Potential local factors might include:

- mitochondrial biogenesis;
- immune infiltration;
- extracellular matrix;
- cell cycle;
- apoptosis;
- ribosome biogenesis.

For this dataset, a practical implementation would require:

- reducing the analysis to perhaps 20–50 pathway families;
- strong regularisation;
- one fixed general factor;
- only a few candidate local factors;
- bootstrap stability assessment;
- preferably external validation.

Attempting to infer a complex latent network across all 884 pathways with 24 samples would be unreliable.

## 3.7 Use external biological measurements

The strongest way to interpret the global axis is to relate it to measurements not derived from the same RNA-seq matrix.

Potential anchors include:

- MYC protein;
- validated MYC-target activity;
- total RNA per cell;
- mitochondrial mass;
- mtDNA copy number;
- respiratory activity;
- OXPHOS protein abundance;
- histological cell composition;
- epithelial, immune or stromal fractions;
- proliferation markers;
- independent MYC mouse datasets.

A model could be:

\[
G
\sim
\text{MYC protein}
+
\text{mitochondrial mass}
+
\text{cell composition}
+
\text{technical QC}
\]

This turns PC1 from an anonymous statistical axis into a partially interpretable biological phenotype.

## 3.8 Treat the global axis as a mediator

For selected hypotheses, mediation may be more appropriate than complete residualisation.

For example:

\[
\mathrm{MYC}
\rightarrow
G
\rightarrow
\mathrm{OXPHOS}
\]

or:

\[
\mathrm{MYC}
\rightarrow
G
\rightarrow
\begin{cases}
\mathrm{OXPHOS}\\
\mathrm{Pathway\ B}
\end{cases}
\]

The total MYC effect can then be separated into:

- a component mediated by the global programme;
- a component independent of the global programme.

With 24 samples, the mediation model must remain small and strongly hypothesis-driven.

---

# 4. Recommended reporting framework

For each selected pathway pair, report three levels of association.

## 4.1 Total association

\[
r_{\mathrm{total}}
=
\operatorname{cor}(S_A,S_B)
\]

This captures the full experimental biology and technical variation.

## 4.2 Global-state-mediated association

For a one-factor model, the global covariance component is:

\[
\operatorname{Cov}_{\mathrm{global}}(A,B)
=
\lambda_A\lambda_B\operatorname{Var}(G)
\]

This quantifies how much of their coordination is explained by shared participation in the global programme.

## 4.3 Residual-specific association

\[
r_{\mathrm{specific}}
=
\operatorname{cor}(\varepsilon_A,\varepsilon_B)
\]

This captures coupling not explained by the known design or global programme.

A near-zero residual correlation should be interpreted as:

> The pathways are coordinated primarily through the shared global state, with little evidence for an additional selective relationship.

---

# 5. Recommended workflow for this dataset

## Step 1: characterise the dimensionality

Calculate:

- PC1, PC2 and later variance;
- effective dimensionality;
- correlation between the global score mean and PC1;
- distribution of pathway loadings.

## Step 2: identify technical anchors

Assemble:

- batch and processing variables;
- RNA quality measures;
- mapping and complexity metrics;
- negative-control genes;
- replicate information.

## Step 3: remove anchored technical variation

Perform technical correction before pathway scoring or before estimating the biological factor.

## Step 4: estimate the remaining biological global factor

Use PC1 or a one-factor model after technical adjustment.

## Step 5: decompose the factor

Fit:

\[
G \sim \text{timepoint} * \text{genotype}
\]

Retain both:

- design-predicted global activity;
- within-group global variability.

## Step 6: analyse pathway loadings

Determine which pathways participate most strongly in the global programme. Test whether OXPHOS loading changes with MYC or timepoint.

## Step 7: test selected differential couplings

For a small number of mito–MYC or mito–developmental hypotheses, fit pathway-by-genotype and pathway-by-time interaction models.

## Step 8: inspect gene-level substructure

Determine whether selected pathway pairs contain coherent residual gene modules even when aggregate scores do not correlate.

## Step 9: validate across scoring methods

Compare:

- mean gene-wise z-score;
- singscore;
- GSVA.

## Step 10: bootstrap the complete procedure

Resample mice within the four experimental groups and re-estimate:

- pathway scores;
- technical factors;
- the global factor;
- pathway loadings;
- selected coupling estimates.

## Step 11: seek external validation

Use independent measurements or an external dataset wherever possible.

---

# 6. Interpretation if almost no specific coupling remains

If, after:

- technical adjustment;
- preservation of MYC and time effects;
- alternative scoring;
- condition-dependent analysis;
- gene-level inspection;
- bootstrap assessment;

almost no residual coupling remains, the analysis has not failed.

The defensible conclusion is:

> Most pathway coordination in these samples is mediated through a broad, low-dimensional transcriptomic state. There is little evidence for additional pathway-pair-specific covariance beyond this shared programme.

This would suggest that MYC and development reorganise the transcriptome mainly through a coordinated systems-level response rather than through many independently varying pathway relationships.

That is a substantive biological result.

---

# Final conceptual summary

The dataset is low-dimensional because many of the 884 pathway scores move together along one or a few dominant sample-level axes.

The global axis should not automatically be treated as unwanted noise. It likely combines:

- MYC biology;
- developmental biology;
- cell composition;
- whole-sample state;
- technical effects.

The central question should not be only:

\[
\text{What correlation survives after removing PC1?}
\]

The more informative questions are:

\[
\boxed{
\begin{aligned}
&\text{What generates the global axis?}\\
&\text{Which pathways participate most strongly in it?}\\
&\text{Does MYC alter pathway loading on the axis?}\\
&\text{Does coupling change by timepoint or genotype?}\\
&\text{Are there local gene modules beyond the global programme?}
\end{aligned}
}
\]

The recommended strategy is:

\[
\boxed{
\text{technical anchoring}
\rightarrow
\text{global biological factor}
\rightarrow
\text{loading and mediation analysis}
\rightarrow
\text{selected residual coupling tests}
}
\]

This preserves the global biological signal while still allowing specific pathway relationships to be investigated.
