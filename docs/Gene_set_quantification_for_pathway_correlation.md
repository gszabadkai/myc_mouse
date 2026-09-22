# Choosing a gene-set quantifier for pathway–pathway correlation analysis

## Study context

The dataset contains:

- **24 samples**
- **4 experimental groups**
  - 6W WT
  - 6W MYC+
  - 12W WT
  - 12W MYC+
- **884 gene sets**
- minimal gene-set redundancy
- a strong common sample-level axis in the current GSVA score matrix

The main analytical aim is not simply to detect differential pathway activity. It is to identify **gene sets whose activities are correlated across samples**, while distinguishing specific pathway coupling from a broad transcriptomic background shared by many pathways.

The current GSVA analysis shows:

- high correlations between apparently unrelated gene sets;
- `mito_oxphos` correlating almost perfectly with the mean score across all 884 sets;
- a major reduction in pairwise correlations after removal of the global sample-level score axis.

This raises two separate questions:

1. Is GSVA the most suitable pathway quantifier for estimating pathway–pathway correlations?
2. Regardless of the quantifier, how should the broad common pathway factor be handled?

Changing the quantifier may improve interpretability, but it will not automatically remove genuine global biological variation.

---

## 1. Recommended primary quantifier: mean gene-wise z-score

For the specific purpose of estimating pathway–pathway coupling, the most interpretable primary score is a simple linear summary of standardized gene expression.

For gene \(g\) in sample \(s\):

\[
Z_{sg}
=
\frac{E_{sg}-\overline{E}_{g}}
{\operatorname{SD}(E_g)}
\]

where \(E_{sg}\) is appropriately transformed expression, such as VST, rlog or voom/logCPM.

For gene set \(k\), containing \(m_k\) genes:

\[
A_{sk}
=
\frac{1}{m_k}
\sum_{g\in k} Z_{sg}
\]

This score represents the mean relative expression of the genes belonging to the pathway.

The GSVA package's `zscore` method uses a closely related scaling:

\[
A_{sk}
=
\frac{\sum_{g\in k} Z_{sg}}
{\sqrt{m_k}}
\]

For correlations between pathways, dividing by \(m_k\) or by \(\sqrt{m_k}\) gives the same result because each pathway is only multiplied by a constant.

---

## 2. Why a linear score is preferable for correlation analysis

A linear pathway score has a direct covariance interpretation.

For two non-overlapping gene sets \(k\) and \(l\):

\[
\operatorname{Cov}(A_k,A_l)
=
\frac{1}{m_km_l}
\sum_{g\in k}
\sum_{h\in l}
\operatorname{Cov}(Z_g,Z_h)
\]

Therefore, a correlation between pathway scores reflects the average cross-gene covariance between the two sets.

Because the gene-set collection has already been checked and has minimal redundancy, correlations are unlikely to be dominated by repeated pathway definitions or extensive reuse of the same genes.

The resulting biological question is clear:

> Do genes belonging to pathway A collectively covary with genes belonging to pathway B across mice?

This is more transparent than correlating two nonlinear enrichment statistics.

### Advantages of a mean z-score

- linear and easy to interpret;
- symmetric between samples and pathways;
- directly traceable to gene-level covariance;
- straightforward to residualise for design and technical covariates;
- easy to bootstrap;
- does not require an enrichment walk;
- less dependent on parameter choices than GSVA or ssGSEA;
- simple to calculate separately for total and design-adjusted expression.

---

## 3. Why GSVA is not necessarily optimal for this purpose

GSVA is valuable for:

- sample-level pathway scoring;
- differential pathway activity;
- ranking samples according to enrichment;
- detecting coordinated changes across members of a gene set.

However, it was not specifically designed to provide an unbiased covariance structure between pathways.

GSVA scores are:

- nonlinear;
- based on transformed gene-level distributions;
- cohort-dependent;
- competitive relative to genes outside the set;
- potentially affected by broad changes in within-sample expression ranks.

Consequently, the correlation between two GSVA scores is less directly interpretable than the correlation between two linear pathway summaries.

This does not make GSVA invalid. It means that GSVA is better retained as:

- a differential pathway analysis method;
- an alternative scoring method;
- a sensitivity analysis.

For pathway–pathway coupling, a linear score is a cleaner primary measure.

---

## 4. Two different ways to construct the linear score

## 4.1 Total pathway activity score

Start from transformed expression without removing the experimental design.

This score retains:

- developmental effects;
- MYC effects;
- the timepoint-by-genotype interaction;
- broad biological state;
- cell-composition differences;
- technical residual variation.

The correlation between two such scores measures **total coordination across the complete experiment**.

---

## 4.2 Design-adjusted pathway activity score

First regress each gene on the known design:

\[
E_g
\sim
\text{timepoint}
*
\text{genotype}
+
\text{technical covariates}
\]

Then standardize the gene residuals and calculate the pathway mean.

This correlation asks:

> Among mice after accounting for the expected four-group structure, do the genes in pathway A still covary with the genes in pathway B?

This is closer to within-condition pathway coupling.

Both total and design-adjusted scores should be retained because they answer different biological questions.

---

## 5. Recommended R implementation

Assume:

- `expr` is a genes-by-samples VST, rlog or voom/logCPM matrix;
- `meta` contains one row per sample;
- `gene_sets` is a named list of gene symbols or identifiers.

### 5.1 Align metadata and expression

```r
meta <- meta[colnames(expr), , drop = FALSE]
```

### 5.2 Create the experimental design

```r
X <- model.matrix(
  ~ timepoint * genotype,
  data = meta
)
```

Add independently justified technical covariates where appropriate:

```r
X <- model.matrix(
  ~ batch + RIN + timepoint * genotype,
  data = meta
)
```

Avoid adjusting for variables that may be biological mediators of MYC or development unless the analytical question explicitly requires it.

### 5.3 Residualise each gene

```r
# expr is genes x samples
# qr.resid requires samples x genes

gene_residuals <- qr.resid(
  qr(X),
  t(expr)
)

gene_residuals <- t(gene_residuals)
```

### 5.4 Standardize each gene across samples

```r
gene_z <- t(scale(t(gene_residuals)))
```

Remove genes with zero variance or non-finite values:

```r
keep <- apply(
  gene_z,
  1,
  function(x) all(is.finite(x))
)

gene_z <- gene_z[keep, , drop = FALSE]
```

### 5.5 Calculate the gene-set means

```r
score_gene_set <- function(
    genes,
    gene_z,
    minimum_genes = 5L
) {
  genes <- intersect(
    genes,
    rownames(gene_z)
  )

  if (length(genes) < minimum_genes) {
    return(
      rep(NA_real_, ncol(gene_z))
    )
  }

  colMeans(
    gene_z[genes, , drop = FALSE],
    na.rm = TRUE
  )
}

pathway_scores <- vapply(
  gene_sets,
  score_gene_set,
  FUN.VALUE = numeric(ncol(gene_z)),
  gene_z = gene_z
)

rownames(pathway_scores) <- colnames(gene_z)
```

The resulting matrix is:

- rows: samples;
- columns: gene sets.

---

## 6. Directional and bidirectional signatures

A simple mean works best when the member genes are expected to change in the same direction.

For signatures containing separate upregulated and downregulated components, use:

\[
A_{sk}
=
\operatorname{mean}(Z_{\mathrm{up}})
-
\operatorname{mean}(Z_{\mathrm{down}})
\]

Example:

```r
signed_score <- function(
    up_genes,
    down_genes,
    gene_z
) {
  up_genes <- intersect(
    up_genes,
    rownames(gene_z)
  )

  down_genes <- intersect(
    down_genes,
    rownames(gene_z)
  )

  up_score <- colMeans(
    gene_z[up_genes, , drop = FALSE],
    na.rm = TRUE
  )

  down_score <- colMeans(
    gene_z[down_genes, , drop = FALSE],
    na.rm = TRUE
  )

  up_score - down_score
}
```

This is particularly important for transcription-factor or perturbational signatures containing both activated and repressed targets.

Without signed scoring, opposing components may cancel even when the pathway is strongly active.

---

## 7. Secondary method: singscore

`singscore` is the most useful independent sensitivity analysis.

It calculates pathway scores from gene ranks within each individual sample. A sample can therefore be scored without using the expression distributions of the other samples.

### Advantages

- cohort-independent scoring;
- rank-based robustness;
- less sensitive to a small number of extreme expression values;
- suitable for one-sample-at-a-time scoring;
- useful for checking whether the GSVA background depends on cohort estimation;
- supports directional signatures.

### Limitations

- bounded rank-based scores;
- reduced sensitivity to expression magnitude;
- correlations may be nonlinear;
- large coordinated changes in absolute expression can be compressed;
- pathway scores can still share a biological global factor.

For singscore matrices, examine:

- Pearson correlation;
- Spearman correlation;
- global factor structure;
- correlation with corresponding linear z-scores.

---

## 8. Interpreting comparison across scoring methods

Generate three pathway-score matrices:

1. current GSVA scores;
2. mean gene-wise z-scores;
3. singscore scores.

For each matrix calculate:

- median correlation among randomly selected pathway pairs;
- distribution of all pairwise correlations;
- correlation between `mito_oxphos` and the mean pathway score;
- variance explained by PC1;
- distribution of PC1 loadings;
- correlation of corresponding pathway scores across methods;
- pairwise correlations after design adjustment;
- pairwise correlations after global-factor adjustment.

### Interpretation table

| Observation | Likely interpretation |
|---|---|
| Strong global factor in GSVA, z-score and singscore | Predominantly biological or compositional global state |
| Strong in GSVA but weak in singscore | Important cohort-dependent GSVA component |
| Strong in z-score but weaker in rank-based methods | Driven partly by expression magnitude |
| Strong in both linear and rank-based methods | Robust broad ordering of samples |
| Specific residual correlation appears in all methods | Strong evidence for selective pathway coupling |
| Specific residual correlation appears only in one method | Method-sensitive and not yet robust |

---

## 9. Other possible gene-set quantifiers

## 9.1 ssGSEA

ssGSEA ranks genes within each sample and calculates an enrichment statistic.

### Potential advantages

- sample-level score;
- rank-based;
- less dependent on absolute expression magnitude;
- widely used.

### Limitations for pathway coupling

- nonlinear enrichment statistic;
- competitive scoring;
- still permits broad common-mode effects;
- less directly interpretable as cross-gene covariance;
- unlikely to provide a major conceptual advantage over singscore.

### Recommendation

Use only as an optional sensitivity analysis. It is not the preferred primary method for the coupling problem.

---

## 9.2 PLAGE or pathway eigengene

PLAGE uses singular value decomposition to derive the dominant expression pattern within each gene set.

For each pathway, it is similar to taking the first principal component of its member genes.

### Potential advantages

- captures the dominant coordinated expression programme;
- can downweight noisy pathway genes;
- useful for highly coherent complexes;
- may work well for selected mitochondrial respiratory-chain modules.

### Limitations

- sign is arbitrary;
- pathway weights are learned from the same 24 samples;
- unstable when the sample size is small;
- PC1 may represent only a subset of genes;
- different pathways can have very different internal structures;
- pathway scores become less directly comparable;
- strong cohort dependence.

### Recommendation

Do not use as the primary score for all 884 pathways.

It may be informative for selected, internally coherent pathways such as:

- respiratory-chain complexes;
- mitochondrial ribosome;
- TCA-cycle modules;
- well-defined MYC target modules.

---

## 9.3 PROGENy and footprint-based methods

Footprint methods estimate signalling-pathway activity from genes that respond downstream of pathway perturbation rather than from the genes physically belonging to the pathway.

### Advantages

- often closer to functional signalling activity;
- can be more causal in interpretation;
- avoids some problems of pathway-member expression;
- useful for selected signalling pathways.

### Limitations

- available only for a restricted set of pathways;
- cannot replace scoring of 884 general gene sets;
- not suited to many mitochondrial, structural or developmental pathways;
- depends on the quality and context of the footprint model.

### Recommendation

Use as an independent validation for selected signalling pathways, not as a universal replacement for the gene-set library.

---

## 9.4 Simple unstandardized pathway mean

An alternative is:

\[
A_{sk}
=
\frac{1}{m_k}
\sum_{g\in k} E_{sg}
\]

This gives highly expressed genes more influence and can cause pathways containing abundant transcripts to dominate.

### Recommendation

Do not use as the primary measure.

Gene-wise standardization before aggregation is preferable because every gene contributes on a comparable scale.

---

## 9.5 Weighted pathway scores

A pathway score can be generalized to:

\[
A_{sk}
=
\sum_{g\in k} w_g Z_{sg}
\]

Weights may represent:

- prior biological importance;
- direction of regulation;
- reproducibility;
- gene reliability;
- loadings estimated in an independent dataset.

### Caution

Weights learned from the same 24 samples risk overfitting and can inflate apparent correlations.

### Recommendation

Use weights only when they are:

- biologically pre-defined;
- learned from independent data;
- validated by cross-validation or external replication.

---

## 10. The global factor remains a separate issue

A better pathway quantifier does not remove the need to adjust for the broad sample-level axis.

After calculating linear pathway scores:

\[
S_{sk}
=
\mu_k
+
X_s\beta_k
+
\lambda_k G_s
+
\varepsilon_{sk}
\]

where:

- \(X_s\) contains timepoint, genotype and interaction;
- \(G_s\) is the common pathway-score factor;
- \(\varepsilon_{sk}\) is pathway-specific activity.

The recommended sequence remains:

\[
\boxed{
\text{gene-level design adjustment}
\rightarrow
\text{linear pathway scoring}
\rightarrow
\text{global-factor adjustment}
\rightarrow
\text{specific pathway correlation}
}
\]

---

## 11. Estimate and remove the global factor

Once the pathway-score matrix has been calculated:

```r
pc <- prcomp(
  pathway_scores,
  center = TRUE,
  scale. = FALSE
)

global_factor <- pc$x[, 1]
global_loadings <- pc$rotation[, 1]
```

Remove the factor:

```r
X_global <- model.matrix(
  ~ global_factor
)

pathway_specific <- qr.resid(
  qr(X_global),
  pathway_scores
)
```

If the experimental design was not already removed at gene level, include it here:

```r
X_global <- model.matrix(
  ~ timepoint * genotype + global_factor,
  data = meta
)

pathway_specific <- qr.resid(
  qr(X_global),
  pathway_scores
)
```

For selected confirmatory pairs, estimate PC1 after excluding the two pathways being tested.

---

## 12. Recommended scoring and correlation hierarchy

For each selected pathway pair, calculate the following.

### Level 1: raw correlation

Correlation between pathway scores from the unadjusted expression matrix.

Interpretation:

> total coordination across the full experiment.

### Level 2: design-adjusted correlation

Correlation between scores calculated from gene-level residuals after removing timepoint, genotype and interaction.

Interpretation:

> coupling beyond the four expected group means.

### Level 3: global-factor-adjusted correlation

Correlation between pathway residuals after removal of the common pathway-score factor.

Interpretation:

> specific pathway coupling beyond the broad transcriptomic state.

### Level 4: cross-method robustness

Compare the result using:

- mean gene-wise z-score;
- singscore;
- GSVA.

Interpretation:

> robustness to the definition of pathway activity.

---

## 13. Statistical considerations with 24 samples

With 24 samples, pathway correlations have substantial sampling uncertainty.

After fitting:

- timepoint;
- genotype;
- interaction;
- global factor;
- technical covariates;

the effective residual degrees of freedom become limited.

Therefore:

- avoid constructing a fully confirmatory 884-by-884 network;
- pre-specify biologically important mito–MYC and developmental pairs;
- treat the complete network as exploratory;
- bootstrap the entire scoring and correction procedure;
- avoid pathway correlations within individual groups of six mice except as descriptive results;
- report confidence intervals rather than only point estimates.

---

## 14. Bootstrap recommendation

Resample mice within the four experimental groups.

For each bootstrap iteration:

1. resample mice within each group;
2. refit the gene-level design model;
3. calculate gene residuals;
4. standardize genes;
5. recalculate pathway scores;
6. re-estimate the global factor;
7. calculate the selected pathway correlation.

This captures uncertainty from:

- sample selection;
- gene-level residualisation;
- pathway scoring;
- estimation of the global factor;
- final correlation.

The pathway scores and global factor should not be held fixed across bootstrap iterations.

---

## 15. A practical comparison pipeline

### Matrix A: GSVA

Use the current GSVA matrix.

### Matrix B: linear z-score

Calculate mean gene-wise standardized expression.

Produce:

- total scores;
- design-adjusted scores.

### Matrix C: singscore

Calculate within-sample rank-based pathway scores.

For all three matrices, record:

```text
1. Correlation of mito_oxphos with the global mean
2. Median random-pair correlation
3. PC1 variance explained
4. Distribution of PC1 loadings
5. Raw mito–MYC correlation
6. Design-adjusted mito–MYC correlation
7. Global-factor-adjusted mito–MYC correlation
8. Bootstrap confidence interval
```

This comparison will determine whether the common background is:

- intrinsic to the expression data;
- strengthened by GSVA;
- dependent on expression magnitude;
- robust across scoring frameworks.

---

## 16. Recommended primary method for this dataset

The preferred primary analysis is:

### Step 1

Transform the RNA-seq count matrix using VST or voom/logCPM.

### Step 2

Fit each gene to:

\[
\text{expression}
\sim
\text{timepoint}
*
\text{genotype}
+
\text{justified technical covariates}
\]

### Step 3

Standardize the gene-level residuals across samples.

### Step 4

Calculate the mean standardized residual expression for each gene set.

### Step 5

Estimate PC1 of the resulting pathway-score matrix.

### Step 6

For selected pathway pairs, remove PC1 using a leave-pair-out procedure.

### Step 7

Calculate pathway-specific Pearson and Spearman correlations.

### Step 8

Validate the main results with singscore.

### Step 9

Retain GSVA for:

- differential pathway analysis;
- comparison with earlier results;
- demonstrating method robustness.

---

## 17. Final recommendation

For the aim of finding correlated gene sets, use the following hierarchy:

1. **Mean gene-wise z-score**
   - preferred primary measure;
   - most interpretable for pathway covariance.

2. **singscore**
   - preferred independent sensitivity analysis;
   - useful for testing cohort dependence and rank robustness.

3. **GSVA**
   - retain for differential pathway analysis and comparison;
   - do not rely on it alone for pathway-coupling inference.

4. **PLAGE**
   - use selectively for highly coherent pathways;
   - not recommended for the complete 884-set collection.

5. **ssGSEA**
   - optional;
   - unlikely to add major value beyond singscore and GSVA.

6. **PROGENy or other footprint methods**
   - useful for selected signalling pathways;
   - not a replacement for the complete pathway library.

The most defensible coupling measure is therefore:

\[
\boxed{
\text{correlation between linear pathway scores}
\text{ after design and global-factor adjustment}
}
\]

with singscore and GSVA used to establish robustness.

A strong association that survives:

- design adjustment;
- global-factor adjustment;
- bootstrap analysis;
- alternative pathway quantification;

would provide substantially stronger evidence for selective biological coupling than a raw GSVA correlation alone.
