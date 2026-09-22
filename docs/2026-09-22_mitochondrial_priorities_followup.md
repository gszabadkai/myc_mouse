---
date: 2026-09-22
tags: [project/myc_mouse, follow-up, future-paper, human, tcga, mitopps, exploratory, parked]
status: parked follow-up thread -- after the manuscript. Nothing here is for the current paper
relates-to:
  - docs/handoff.md                                  ("Follow-up threads, after the manuscript")
  - docs/2026-09-21_two_timeline_verification.md     (section 10: the FOXO3 record; 10.3, Esr1)
  - /Users/gs/code/myc_human_validation              (frozen at d3ac60e; read-only input)
  - /Users/gs/code/myc_human_exploratory             (data/menegollo_biclusters; read-only input)
---

# Mitochondrial priorities in high-mitochondria human breast tumours: a parked follow-up

**Parked 2026-09-22.** Exploratory, cross-sectional and read-only. No panel was built, no
manuscript text was touched, and nothing was written to either human repository.

---

## 1. The question, for a future paper

**What reshapes mitochondrial priorities during tumourigenesis?**
- The mouse arm of the current manuscript shows a developmental programme that reallocates the
  nuclear-encoded mitochondrial transcriptome substantially, and independently of MYC.
- Whether ER-driven and MYC-driven routes to a high-mitochondria state in human breast tumours
  allocate the compartment differently is the natural extension. It is unanswered.

## 2. What the literature has and has not established

The author's summary, 2026-09-22.
- **MYC-high and ER+ tumours are each independently described as OXPHOS-high.** Both are
  reported as vulnerable to complex I inhibition, and IACS-010759 appears in both literatures.
  Those references were not pulled in this session.
- **Alles et al. 2009** (*PLoS ONE* 4(3):e4710; PMID 19270750; PMC2650420;
  doi:10.1371/journal.pone.0004710; verified on PubMed 2026-09-22). In ER-negative basal
  tumours, MYC transcriptional activity reproduced much of the indirect oestrogen response seen
  in ER+ cells.
- **Nobody appears to have compared the two high-mitochondria populations on HOW the compartment
  is allocated,** only on whether it is large. MitoPPS (Monzel et al. 2025) postdates the
  Menegollo analysis (*Cancer Res* 2024), which is why.

---

## 3. What was run on 2026-09-22, and what it returned

**The headline is a negative.** The question was asked and the design could not answer it. The
numbers are all here, not summarised.

### 3.1 The groups

- **Cohort:** TCGA-BRCA, the 849 tumours with Menegollo's published fork calls (snapshot pinned
  at upstream `8fdbb34`). The other 188 of 1,037 cannot be recovered from the snapshot.
  - **METABRIC's** forks cannot be attached to patients: the sample-identifier vector is
    missing.
  - **SCAN-B** has no fork calls.
- **Both upper forks are the mitochondria-high ones:** OXPHOS GSVA medians +0.37 (MB1 upper)
  and +0.34 (MB2 upper).
  - **MB2 upper, "MYC-high":** n 86; MYC score M_a median +0.56; 38% ER+.
  - **MB1 upper, "ESR1-high":** n 93; 99% ER+.
- **The groups are not mutually exclusive.** 21 tumours are in both, against 9.4 expected
  (odds ratio 3.1, 95% CI 1.7 to 5.5, Fisher p 0.00016).

  | group | n | ER+ | main PAM50 | M_a (median) |
  |---|---|---|---|---|
  | MB2-upper only | 65 | 10/60 | Basal | +0.59 |
  | both forks | 21 | 21/21 | LumB (17) | +0.41 |
  | MB1-upper only | 72 | 65/66 | LumB/LumA | -0.07 |

- **ESR1 against the MYC score, within the 158 high-mitochondria tumours:**
  - Spearman -0.69 with M_a (-0.76 to -0.60), and -0.72 with M_b;
  - -0.22 (p 0.006) with the three groups held fixed.
  - This is an association only.

### 3.2 Circularity, checked before any scoring

**What defines a TCGA fork call** (read from the upstream scripts). It is PC1 over the top
1,000 genes of a METABRIC bicluster correlation vector, by |correlation| and split by sign.
- **MB2's block is confirmed.** A score built from it tracks the published PC1 at Spearman
  0.975, and a PC1 projection at -0.991.
- **MB1's block is the best-supported candidate:** the ICT1 vector, at 0.88 against 0.73 for the
  next candidate.

**The blocks are mostly not mitochondrial.**
- MitoCarta genes are **7.4% of the MB1 block and 13.0% of MB2's.**
- Block genes are about 6.1% (MB1) and 10.4% (MB2) of every pathway's MitoPPS denominator.

**Numerator overlap for the pathways that matter here.** "High side" is the arm of the MB2
block that is high in MB2-upper tumours.

| pathway | n | in MB1 block | in MB2 block, high side | times MB2's denominator rate |
|---|---|---|---|---|
| Protein import and sorting | 47 | 10.6% | 14.9% | 1.4 |
| Translation | 154 | 10.4% | 18.2% | 1.7 |
| OXPHOS subunits | 88 | 6.8% | 18.2% (plus 1.1% on the other side) | 1.9 |
| Mitochondrial ribosome | 82 | 17.1% | 28.0% | 2.7 |
| CIII subunits | 9 | 11.1% | 33.3% | 3.2 |

- **29 of 142 pathways** have at least 20% of their numerator genes in a block. Most are 3- to
  9-gene pathways.
- **The primary contrast is exposed to MB2's block only.** Both of its groups are MB1-upper;
  only the 21 are MB2-upper.
- Translation and respiratory genes sit **1.4 to 3.2 times above the denominator rate on MB2's
  high side.** So **any positive difference in those pathways is partly a restatement of the
  fork call.**

### 3.3 Power, stated before the results

- **The primary contrast (21 against 65)** detects **0.70 SD** at 80% power (alpha 0.05), and
  **1.11 SD** after correction across 142 pathways.
- **In MitoPPS units,** 0.70 SD is 0.056 for import, 0.062 for translation and 0.176 for OXPHOS
  subunits.
- **The secondary contrast (65 against 65)** detects 0.49 SD.
- **The observed primary effects are about 0.3 SD.**

> **THE DESIGN REQUIREMENT FOR ANY FOLLOW-UP: about 116 tumours in the both-forks group to
> detect 0.3 SD.** This assumes 80% power, two-sided alpha 0.05, one pre-declared pathway, and a
> comparison group 3.1 times larger (about 359). Across 142 pathways with Bonferroni correction
> it is 287. With equal groups it is 175 per group.

**What that means for cohort size.** This is arithmetic at TCGA's prevalence (21 of 849, 2.5%),
assuming fork calls could be made elsewhere and the prevalence held.
- 116 such tumours need a cohort of about 4,700.
- SCAN-B (3,207) would give about 79, and METABRIC (1,981) about 49.
- Neither reaches 116 alone. The three cohorts together (about 149) would, but only by
  **meta-analysing per-cohort effects**. MitoPPS is cohort-relative and is never pooled across
  cohorts.

### 3.4 The primary contrast: 21 ER+ both-forks tumours against 65 ER+ MB1-upper-only

Seven of the 72 MB1-upper-only tumours lack a positive ER call and were excluded. Using all 72
gives import +0.026, translation +0.026 and OXPHOS subunits -0.001, the same.

The prediction under test came from the mouse and was stated before looking: MYC promotes
protein import, mitochondrial translation and the respiratory chain. Its three pathways were
fixed before any group was scored.

| pathway | both | MB1-only | difference | 95% CI | d | p | BH | null pct (all genes) | null pct (within compartment) |
|---|---|---|---|---|---|---|---|---|---|
| Protein import and sorting | 1.025 | 0.997 | **+0.027** | -0.011 to +0.066 | **0.34** | 0.16 | 0.34 | 97.7 | 90.4 |
| Translation | 1.082 | 1.055 | **+0.027** | -0.018 to +0.072 | **0.31** | 0.23 | 0.43 | 99.7 | 98.0 |
| OXPHOS subunits | 1.114 | 1.121 | **-0.007** | -0.098 to +0.084 | **-0.03** | 0.88 | 0.93 | 94.7 | 93.5 |

- **No interval excludes zero, and none lies below zero.**
- **The verdict, against the rule declared before looking: UNTESTABLE, not refuted.**
  - **HOLDS** needed all three intervals above zero and above the 95th null percentile.
  - **FAILS** needed an interval below zero, or all three estimates at or below zero.
  - Neither happened.
- **OXPHOS subunits show nothing.** That agrees with the mouse's own record: its OXPHOS-subunit
  promotion on MitoPPS was not significant either (+0.112, adjusted p 0.27).
- **Children of the three pathways, for context only:**
  - protein import, sorting and homeostasis (the tier): +0.105 (+0.002 to +0.209), d 0.86, BH
    0.20;
  - mitochondrial ribosome: +0.028 (-0.046 to +0.102);
  - mt-tRNA synthetases: +0.109 (+0.014 to +0.204), BH 0.17;
  - CIII subunits: +0.192 (+0.026 to +0.358), BH 0.17, with 33% of its genes in the MB2 block.
- **Across all 142 pathways,** 32 intervals exclude zero, against about 7 expected. **Seven
  clear BH 0.05:**
  - Six are lower in the both-forks group, and none has a block gene: ABC transporters -0.236,
    Q-linked reactions -0.168, tetrahydrobiopterin synthesis -0.400, autophagy -0.104, sulfur
    metabolism -0.171 and trafficking -0.148.
  - One is higher: Cytochromes +0.370, a 3-gene pathway with 33% in a block.

### 3.5 Null inflation

- **Script 07's all-gene null is shifted for this contrast.** The median pathway sits at the
  **86th percentile**, and 115 of 141 null-tested pathways sit above 50.
- **Content does not simply differ between the groups.**
  - Nuclear OXPHOS content is higher in the 21 (+0.11 log2).
  - The total MitoCarta share is slightly lower (0.112 against 0.122, p 0.08), because mtDNA
    transcripts are lower (MitoPPS -0.47).
- **A within-compartment null is shifted too:** medians -0.06 to -0.23 for the three pathways
  (expression-matched MitoCarta genes only, 10 bins over 1,025 genes).
- **Translation's 98th percentile must be read against that.** Its interval includes zero, and
  18% of its genes are on MB2's high side.

### 3.6 The secondary three-group contrast (confounded by subtype)

MB2-upper-only (65; basal, mostly ER-negative) against ER+ MB1-upper-only (65).

| pathway | difference | 95% CI | null percentile |
|---|---|---|---|
| Protein import and sorting | +0.132 | +0.098 to +0.166 | 86.7 |
| Translation | +0.060 | +0.027 to +0.094 | 75.1 |
| OXPHOS subunits | -0.012 | -0.094 to +0.071 | 47.8 |

- **98 of 142 pathways differ at BH < 0.05,** and the median pathway sits at the 56.7th null
  percentile. 23 of the 98 are block-flagged, against 29 of 142 overall.
- **The groups differ across most of the compartment,** which is what a basal-versus-luminal
  difference looks like.
- **This contrast cannot separate allocation by MYC from allocation by subtype.**

### 3.7 Descriptive only, not tested

The both-forks group sits between the other two, in the same order as the groups' MYC scores
(+0.59, +0.41, -0.07):

| pathway | MB2-upper only | both forks | MB1-upper only |
|---|---|---|---|
| Protein import and sorting | 1.129 | 1.025 | 0.997 |
| Translation | 1.116 | 1.082 | 1.055 |
| mt-tRNA synthetases | 1.178 | 1.062 | 0.953 |

This ordering was noticed after the numbers were seen. It is not a result.

---

## 4. Routes to adequate power -- options, not decisions

- **SCAN-B** (3,207 tumours; about 79 both-forks tumours at TCGA's prevalence). It has no fork
  calls. Two options:
  - run MCbiclust on it;
  - find another way to separate MYC-driven from ER-driven high-OXPHOS tumours that does not
    depend on the Menegollo biclusters. Those were seeded on mitochondrial genes, and their TCGA
    defining blocks carry 7 to 13% MitoCarta genes. This option would also **dissolve** the
    circularity problem rather than quantify it.
- **METABRIC** (about 49). Its forks exist but cannot be attached to patients in the current
  snapshot: `METABRIC_DATA.RData` or the compile script's output is missing (the snapshot
  README, "The METABRIC gap"). The upstream README links the file on Google Drive. Worth
  checking whether the identifiers are recoverable.
- **Gene level rather than pathway level** (raised by the author). Pathway means may be
  averaging away real differences.
  - The FOXO3 target-set analysis earlier today is the worked example of what a gene-level view
    finds behind a flat mean. There the 38-gene mean (+0.072) concealed a heavy positive tail of
    genes MYC represses at 6W (verification note, section 10.1).
  - There, the concealed structure turned out to be the global fade, not a responsive subset.
    That is the caution the same approach carries here.
  - MitoPPS is pathway-level by construction, so a gene-level version needs its own ruler, such
    as each gene's share of the compartment.
- **Across cohorts:** meta-analyse per-cohort effects; never pool cohort-relative scores.

---

## 5. The framing to preserve

**The question has not been asked at adequate power. That is different from the answer being
no.** A future session must not read today's near-zero estimates as evidence of similarity.
- At 21 against 65, a true allocation difference of 0.3 SD is detected only 22% of the time
  (alpha 0.05). It is missed about three times in four.
- The design that can see it needs about 116.

---

## 6. Link back

This connects to the mouse result that a developmental programme reallocates mitochondrial
priorities independently of MYC. That is the current manuscript's Figure 1 (the author's
reference).
- **In the panel manifest** the result is drawn as Fig. 2F (`fig2_wt_mito_contraction.R`) and
  Fig. S2B (`figS2_reallocation_independence.R`), both cited to "s2p1".
- **The number, checked today:** across 143 MitoPathways, MYC's 6W priority effect and the
  wild-type 6-to-12W priority change correlate at **-0.022** (Spearman -0.033), with 72 of 143
  discordant in sign (`background_vs_myc.rds$ruler`).

---

## 7. Provenance

**Inputs, all read-only.**
- **`myc_human_validation` at `d3ac60e`:**
  - `results/tcga_brca_mito_scores.rds` (MitoPPS universe of 142 pathways; `mito_paths`; content
    scores);
  - `results/tcga_brca_linear.rds`, `results/tcga_brca_vst.rds`,
    `results/tcga_brca_covariates.rds` (`er_call`, `PAM50`) and `results/tcga_brca_myc_scores.rds`
    (`M_a`, declared primary there; `M_b`).
- **`myc_human_exploratory/data/menegollo_biclusters/`:** `TCGA.all.biclusters.RNAseq.Rdata`
  (the fork calls), `TCGA_MB{1,2}_RNASeq_data_nonorm.RData` and `METABRIC_sort_data.RData`.

**Fetched read-only with `gh`:**
- `gszabadkai/Menegollo_Bentham` at `8fdbb34`: `R scripts/TCGA MCbiclust analysis/TCGA_RNASeq_Analysis_MB1.R`
  (blob `8f1f9b7`) and `..._MB2.R` (blob `b375142`);
- `gszabadkai/mammary_geneset_library` at `v1.0` (`cbd8f16`):
  `data/raw/user_curated/METABRIC_allSampleParameters_allBiclustergenelists_corrected_groups2.xlsx`
  (blob `b17a7a9`).
- The workbook's `MB{n}.allgenome.cv` sheets equal the snapshot's positional `ICT1.cv[[1]]` and
  `mito.cv[[2]]` in order, which is what lets the blocks be named.

**Construction choices:**
- **Block:** the top 1,000 by |CV|, restricted to the old RNASeqV2 names in
  `TCGA.MB{n}.CV.RNAseq.df$Genes` (exact, as upstream), then reconciled to current symbols
  (exact, then HGNC alias via `org.Hs.eg.db`).
- **MitoPPS:** the saved universe, rebuilt to 4.4e-16.
- **Intervals:** Welch 95% CIs.
- **All-gene null:** script 07's. 2,000 sets per pathway, drawn without replacement within 20
  ventiles of mean linear expression over all 18,115 genes, and scored as a query against the
  universe with the pathway held out. The mtDNA pathway is not null-tested, by script 07's rule.
- **Within-compartment null:** the same, with 10 bins over the 1,025 MitoCarta pathway genes.
- **Seed and ER status:** seed 20260922; ER+ means `er_call == "Positive"`.

**To regenerate.** The computation ran in the session scratchpad and is reproduced verbatim in
the appendix.
1. Fetch the workbook to a folder and point `SC` at it:
   `gh api repos/gszabadkai/mammary_geneset_library/git/blobs/b17a7a99ecf4fa27e8075374a74080337c5fdf05 --jq .content | tr -d '\n' | base64 -d > "$SC/metabric_biclusters.xlsx"`
2. Run the main script, about 2.5 minutes.
3. Then run the within-compartment chunk.

---

## Appendix: the computation, verbatim

Run on 2026-09-22 from the session scratchpad, read-only. Set `SC` to the folder holding the
fetched workbook; the main script also writes its result object there, and nowhere else.

### A.1 The main script: circularity, groups, power, both contrasts, and the all-gene null

```r
# Read-only exploratory check (2026-09-22): do the high-mitochondria MYC-high and
# ER-high TCGA tumours differ in HOW the mitochondrial compartment is allocated?
# Reads the frozen myc_human_validation objects (d3ac60e), the Menegollo snapshot
# in myc_human_exploratory, and the library workbook fetched to this scratchpad.
# Writes nothing outside the scratchpad.
suppressPackageStartupMessages(library(Matrix))
options(width = 220)
V  <- "/Users/gs/code/myc_human_validation/results/"
X  <- "/Users/gs/code/myc_human_exploratory/data/menegollo_biclusters/"
SC <- "/private/tmp/claude-501/-Users-gs-code-myc-mouse/4b5805d8-f9c0-49f2-b1e5-6f67f78be89d/scratchpad/"
set.seed(20260922)
PRED <- c(import = "Protein import and sorting", translation = "Translation",
          respiratory_chain = "OXPHOS subunits")
NSET <- 2000L; NBIN <- 20L; SKIP_NULL <- "mtDNA-encoded OXPHOS subunits"

# ---- load -----------------------------------------------------------------------
ms <- readRDS(paste0(V, "tcga_brca_mito_scores.rds"))
U  <- ms$mitopps_universe
L  <- readRDS(paste0(V, "tcga_brca_linear.rds"))$mat
cv <- as.data.frame(readRDS(paste0(V, "tcga_brca_covariates.rds"))$covariates)
my <- as.data.frame(readRDS(paste0(V, "tcga_brca_myc_scores.rds"))$estimators)
e  <- new.env(); load(paste0(X, "TCGA.all.biclusters.RNAseq.Rdata"), envir = e)
mb <- e$TCGA.all.biclusters.RNAseq.df2; mb$patient <- substr(mb$X, 1, 12)
stopifnot(identical(colnames(U), colnames(L)))

# ---- rebuild the pathway score matrix and prove it reproduces the saved mitoPPS ----
pres <- lapply(ms$mito_paths, function(g) intersect(g, rownames(L)))
pres <- pres[lengths(pres) >= 3]
stopifnot(setequal(names(pres), rownames(U)))
pres <- pres[rownames(U)]
Su <- t(vapply(pres, function(g) colMeans(L[g, , drop = FALSE]), numeric(ncol(L))))
q_pps <- function(Sq, Su) {                       # script 07's .mitopps_query, verbatim
  N <- ncol(Su); P <- nrow(Su); Bi <- 1 / Su
  A <- (Sq %*% t(Bi)) / N
  Sq * ((1 / A) %*% Bi) / P
}
chk <- vapply(c(1L, 50L, nrow(Su)), function(p)
  max(abs(q_pps(Su[p, , drop = FALSE], Su[-p, , drop = FALSE]) - U[p, ])), numeric(1))
cat(sprintf("identity: held-out query reproduces the saved mitoPPS to %.1e\n", max(chk)))
stopifnot(max(chk) < 1e-8)

# ---- the fork-defining blocks (the TCGA scripts' own construction) -------------------
wb <- paste0(SC, "metabric_biclusters.xlsx")
block_of <- function(b) {
  ag <- as.data.frame(suppressMessages(readxl::read_excel(wb, sheet = paste0(b, ".allgenome.cv"))))
  top <- order(abs(ag[[2]]), decreasing = TRUE)[1:1000]
  list(pos = ag[[1]][top][ag[[2]][top] > 0], neg = ag[[1]][top][ag[[2]][top] < 0], all_cv = ag)
}
B1 <- block_of("MB1"); B2 <- block_of("MB2")
tc <- function(n) { f <- new.env(); load(sprintf("%sTCGA_MB%d_RNASeq_data_nonorm.RData", X, n), envir = f)
  get(sprintf("TCGA.MB%d.CV.RNAseq.df", n), envir = f) }
T1 <- tc(1); T2 <- tc(2)
# provenance: which METABRIC vector does each TCGA fork's own correlation vector track?
prov <- sapply(c(MB1 = 1, MB2 = 2), function(n) { Tn <- if (n == 1) T1 else T2
  sapply(c(MB1 = "MB1", MB2 = "MB2"), function(b) { ag <- (if (b == "MB1") B1 else B2)$all_cv
    k <- intersect(Tn$Genes, ag[[1]]); cor(Tn$CV[match(k, Tn$Genes)], ag[[2]][match(k, ag[[1]])], use = "complete.obs") }) })
cat("\ncorrelation of each TCGA fork's CV (columns) with each METABRIC CV (rows):\n"); print(round(prov, 3))
# genes actually used in TCGA = block symbols present among the old RNASeqV2 names (exact, as upstream)
used <- function(B, Tn) lapply(B[c("pos", "neg")], function(g) intersect(g, Tn$Genes))
U1 <- used(B1, T1); U2 <- used(B2, T2)
# reconcile old symbols to the matrix's current symbols (exact first, then HGNC alias)
cur <- rownames(L)
to_cur <- function(g) {
  ex <- g[g %in% cur]; rest <- setdiff(g, cur)
  al <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, rest, "SYMBOL", "ALIAS", multiVals = "first"))
  unique(c(ex, intersect(stats::na.omit(unname(al)), cur)))
}
R1 <- lapply(U1, to_cur); R2 <- lapply(U2, to_cur)
cat(sprintf("\nblocks: MB1 used %d+%d -> in matrix %d+%d | MB2 used %d+%d -> in matrix %d+%d\n",
    length(U1$pos), length(U1$neg), length(R1$pos), length(R1$neg), length(U2$pos), length(U2$neg), length(R2$pos), length(R2$neg)))

# ---- circularity, per MitoPathway ---------------------------------------------------
allmito <- unique(unlist(pres))
circ <- do.call(rbind, lapply(names(pres), function(p) {
  g <- pres[[p]]; den <- unique(unlist(pres[setdiff(names(pres), p)]))
  f <- function(x, s) length(intersect(x, s)) / length(x)
  data.frame(pathway = p, n = length(g),
    MB1_num = f(g, c(R1$pos, R1$neg)), MB1_num_pos = f(g, R1$pos), MB1_num_neg = f(g, R1$neg), MB1_den = f(den, c(R1$pos, R1$neg)),
    MB2_num = f(g, c(R2$pos, R2$neg)), MB2_num_pos = f(g, R2$pos), MB2_num_neg = f(g, R2$neg), MB2_den = f(den, c(R2$pos, R2$neg)))
}))
circ$flag <- circ$MB1_num >= 0.20 | circ$MB2_num >= 0.20
cat(sprintf("\ncompartment: %d MitoCarta pathway genes in the matrix; in MB1 block %d (%.1f%%), in MB2 block %d (%.1f%%)\n",
    length(allmito), length(intersect(allmito, unlist(R1))), 100 * mean(allmito %in% unlist(R1)),
    length(intersect(allmito, unlist(R2))), 100 * mean(allmito %in% unlist(R2))))
cat(sprintf("the blocks themselves: MB1 %.1f%% MitoCarta, MB2 %.1f%% MitoCarta\n",
    100 * mean(unlist(R1) %in% allmito), 100 * mean(unlist(R2) %in% allmito)))
cat("pathways with >= 20% of numerator genes in a block:", sum(circ$flag), "of", nrow(circ), "\n")

# ---- groups ------------------------------------------------------------------------------
d <- merge(mb[, c("patient", "MB1.fork", "MB2.fork")], cv[, c("patient", "er_call", "PAM50")], by = "patient")
d <- merge(d, my[, c("patient", "M_a")], by = "patient")
d$both    <- d$MB1.fork == "Upper" & d$MB2.fork == "Upper"
d$mb1only <- d$MB1.fork == "Upper" & d$MB2.fork != "Upper"
d$mb2only <- d$MB2.fork == "Upper" & d$MB1.fork != "Upper"
ERp <- !is.na(d$er_call) & d$er_call == "Positive"
gA <- d$patient[d$both & ERp]; gB <- d$patient[d$mb1only & ERp]; gC <- d$patient[d$mb2only]
gB72 <- d$patient[d$mb1only]
iA <- match(gA, colnames(U)); iB <- match(gB, colnames(U)); iC <- match(gC, colnames(U)); iB72 <- match(gB72, colnames(U))
cat(sprintf("\ngroups: both & ER+ %d | MB1-upper-only & ER+ %d (of %d; excluded %d not ER+ called) | MB2-upper-only %d\n",
    length(iA), length(iB), length(gB72), length(gB72) - length(gB), length(iC)))
stopifnot(!anyNA(c(iA, iB, iC, iB72)))

# ---- power ---------------------------------------------------------------------------------
dmin <- function(n1, n2, a) (qnorm(1 - a / 2) + qnorm(0.8)) * sqrt(1 / n1 + 1 / n2)
cat(sprintf("\npower (80%%): primary %d vs %d detects d = %.2f SD at alpha 0.05, %.2f SD at 0.05/%d | secondary %d vs %d: %.2f SD\n",
    length(iA), length(iB), dmin(length(iA), length(iB), 0.05), dmin(length(iA), length(iB), 0.05 / nrow(U)), nrow(U),
    length(iC), length(iB), dmin(length(iC), length(iB), 0.05)))

# ---- per-pathway differences, CIs, matched null ---------------------------------------------
gm <- rowMeans(L)
bin_of <- setNames(cut(rank(gm, ties.method = "first"), breaks = NBIN, labels = FALSE), rownames(L))
by_bin <- split(names(bin_of), bin_of)
draw <- function(g) { b <- bin_of[g]; unlist(lapply(split(b, b), function(k) {
  pool <- by_bin[[as.character(k[[1]])]]; pool[sample.int(length(pool), length(k))] }), use.names = FALSE) }
welch <- function(x, y) { t <- stats::t.test(x, y); c(diff = mean(x) - mean(y), lo = t$conf.int[1], hi = t$conf.int[2], p = t$p.value,
  d = (mean(x) - mean(y)) / sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) / (length(x) + length(y) - 2))) }
gi <- setNames(seq_len(nrow(L)), rownames(L))
res <- do.call(rbind, lapply(rownames(U), function(p) {
  pr <- welch(U[p, iA], U[p, iB]); se <- welch(U[p, iC], U[p, iB]); pr72 <- welch(U[p, iA], U[p, iB72])
  pct_p <- NA_real_; pct_s <- NA_real_
  if (!p %in% SKIP_NULL) {
    g <- pres[[p]]; k <- length(g)
    idx <- unlist(lapply(seq_len(NSET), function(s) gi[draw(g)]), use.names = FALSE)
    G <- sparseMatrix(i = rep(seq_len(NSET), each = k), j = idx, x = 1 / k, dims = c(NSET, nrow(L)))
    Sq <- as.matrix(G %*% L)
    nul <- q_pps(Sq, Su[setdiff(rownames(Su), p), , drop = FALSE])
    nd_p <- rowMeans(nul[, iA, drop = FALSE]) - rowMeans(nul[, iB, drop = FALSE])
    nd_s <- rowMeans(nul[, iC, drop = FALSE]) - rowMeans(nul[, iB, drop = FALSE])
    pct_p <- 100 * mean(nd_p < pr[["diff"]]); pct_s <- 100 * mean(nd_s < se[["diff"]])
  }
  data.frame(pathway = p, n = length(pres[[p]]),
    mA = mean(U[p, iA]), mB = mean(U[p, iB]), mC = mean(U[p, iC]),
    diff = pr[["diff"]], lo = pr[["lo"]], hi = pr[["hi"]], p = pr[["p"]], d = pr[["d"]], null_pct = pct_p,
    diff72 = pr72[["diff"]], lo72 = pr72[["lo"]], hi72 = pr72[["hi"]],
    s_diff = se[["diff"]], s_lo = se[["lo"]], s_hi = se[["hi"]], s_p = se[["p"]], s_null_pct = pct_s)
}))
res$bh <- p.adjust(res$p, "BH"); res$s_bh <- p.adjust(res$s_p, "BH")
res <- merge(res, circ, by = c("pathway", "n"))
saveRDS(list(res = res, circ = circ, prov = prov, groups = list(A = gA, B = gB, C = gC, B72 = gB72), d = d),
        paste0(SC, "mitopps_forks_result.rds"))
cat("\nsaved to scratchpad\n")
```

### A.2 The within-compartment null (run after A.1)

This reproduces the four within-compartment percentiles in 3.4 and 3.5 exactly (checked on
2026-09-22): 90.4, 98.0, 93.5 and 99.9.

```r
# Within-compartment null for the three prediction pathways and the import tier.
# Run after the main script (it uses L, pres, Su, q_pps, U, iA, iB).
mito <- unique(unlist(pres))
gmm  <- rowMeans(L[mito, ])
binm <- setNames(cut(rank(gmm, ties.method = "first"), breaks = 10L, labels = FALSE), mito)
bybm <- split(names(binm), binm)
set.seed(20260922)
for (p in c("Protein import and sorting", "Translation", "OXPHOS subunits",
            "Protein import, sorting and homeostasis")) {
  g <- pres[[p]]; k <- length(g); b <- binm[g]
  idx <- unlist(lapply(seq_len(2000), function(s) match(unlist(lapply(split(b, b), function(kk) {
    pool <- bybm[[as.character(kk[[1]])]]; pool[sample.int(length(pool), length(kk))] }),
    use.names = FALSE), rownames(L))), use.names = FALSE)
  G   <- sparseMatrix(i = rep(seq_len(2000), each = k), j = idx, x = 1 / k, dims = c(2000, nrow(L)))
  nul <- q_pps(as.matrix(G %*% L), Su[setdiff(rownames(Su), p), , drop = FALSE])
  nd  <- rowMeans(nul[, iA]) - rowMeans(nul[, iB])
  obs <- mean(U[p, iA]) - mean(U[p, iB])
  cat(sprintf("%-40s observed %+.3f | null median %+.3f (95%% %+.3f to %+.3f) | percentile %.1f\n",
              p, obs, median(nd), quantile(nd, 0.025), quantile(nd, 0.975), 100 * mean(nd < obs)))
}
```
