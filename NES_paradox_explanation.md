# NES Paradox: Higher Pathway Enrichment Despite Lower Gene-Level Effects

## The Observation

In our MYC-driven liver cancer model, we observed an apparent paradox:

- **Gene-level**: Individual mitochondrial genes show **lower** log2 fold changes (LFC) at 12W compared to 6W
- **Pathway-level**: Mitochondrial pathway NES scores are **higher** at 12W compared to 6W

This document explains why this is not a contradiction but reflects the underlying biology and statistics of pathway enrichment analysis.

---

## Key Statistics

| Metric | 6W | 12W |
|--------|-----|-----|
| Mean NES (MitoCarta pathways) | 2.255 | 2.526 |
| Variance of genome-wide Wald statistics | 3.04 | 1.49 |
| % genes with extreme stats (\|stat\| > 2) | 22.5% | 9.5% |
| Mean sequencing depth | ~23M reads | ~13M reads |

---

## The Explanation

### 1. NES is a Rank-Based Metric, Not Magnitude-Based

GSEA/fgsea calculates NES by:
1. Ranking all genes by their Wald statistic
2. Measuring whether pathway genes cluster at the extremes of this ranked list
3. **Normalizing** for gene set size and background distribution

The key insight: **NES reflects relative position, not absolute effect size**.

### 2. The "Quieter Background" Phenomenon

At 6W (early MYC activation):
- Broad transcriptional chaos: 22.5% of genes show extreme statistics
- MitoCarta genes are upregulated, but compete with many other dysregulated genes
- High variance in background makes it harder to detect specific pathway enrichment

At 12W (established MYC phenotype):
- Transcriptional response is more focused: only 9.5% of genes are extreme
- MitoCarta genes **maintain their relative ranking** despite lower absolute effects
- With less competition from other genes, the enrichment signal is stronger

### 3. Analogy

Imagine a race:
- At 6W: Your team runs fast (10 mph), but so does everyone else. Hard to stand out.
- At 12W: Your team runs slower (7 mph), but others are walking. You dominate.

---

## Formal Tests

We conducted four formal tests to validate this explanation:

### Test 1: Permutation Test on Rank Preservation

**Question**: Do MitoCarta genes maintain their relative ranking between timepoints?

**Result**: 
- Observed rank change: +0.006 (essentially unchanged)
- P-value: 0.508 (not unusual compared to random gene sets)

**Interpretation**: MitoCarta genes maintain nearly identical relative positions in the ranked list at both timepoints.

### Test 2: Paired NES Comparison

**Question**: Is the NES increase at 12W statistically significant?

**Result**:
- Mean NES at 6W: 2.255
- Mean NES at 12W: 2.526
- Paired t-test P-value: 3.27 × 10⁻⁷
- 35/39 pathways (89.7%) show higher NES at 12W

**Interpretation**: The NES increase is highly significant and consistent across pathways.

### Test 3: Uniform Scaling Simulation

**Question**: If we uniformly compress all Wald statistics (mimicking reduced effect sizes), does NES change?

**Result**: Correlation between original and scaled NES = 1.000

**Interpretation**: Uniform scaling doesn't change NES because it preserves ranks. This confirms NES is purely rank-based.

### Test 4: Downsampling Analysis (Critical Test)

**Question**: Could lower sequencing depth at 12W explain the higher NES?

**Method**: Downsample 6W reads to match 12W depth, re-run DESeq2 and fgsea.

**Result**:
| Comparison | NES Change |
|------------|------------|
| 6W → 12W (biological) | **+0.270** |
| 6W → downsampled 6W (technical) | **−0.111** |

**Interpretation**: 
- Lower sequencing depth **decreases** NES (more noise → weaker signal)
- The observed 12W NES increase occurs **despite** the depth disadvantage
- This directly contradicts the "lower counts cause higher NES" hypothesis

---

## Conclusion

The higher NES at 12W compared to 6W is a **biological phenomenon**, not a technical artifact:

1. **Early MYC activation (6W)** causes widespread transcriptional dysregulation
2. **Established MYC phenotype (12W)** shows a more focused response
3. Mitochondrial genes maintain their relative importance as background noise decreases
4. This results in **stronger enrichment scores** despite smaller absolute effect sizes

The downsampling analysis provides definitive evidence: if lower depth caused higher NES, downsampling 6W should have increased NES. Instead, it decreased NES by 0.111, while the biological 6W→12W transition increased NES by 0.270.

---

## Files

- **Analysis script**: `scripts/NES_paradox_analysis.R`
- **Test results**: `results/fgsea/NES_paradox_tests.rds`
- **Original scatter plot**: `results/fgsea/NES_scatter_sig.pdf`

---

## References

- Subramanian et al. (2005). Gene set enrichment analysis. PNAS.
- Korotkevich et al. (2021). Fast gene set enrichment analysis. bioRxiv.
