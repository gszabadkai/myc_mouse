# scripts/40_background_vs_myc_decomposition.R
# =============================================================================
# Block B -- IS MYC'S MITOCHONDRIAL EFFECT SHAPED BY A CHANGING BACKGROUND?
# =============================================================================
#
# The question the figures so far only hint at. fig01/fig02 read Myc WITHIN a
# timepoint; figS6/figS7 read each genotype ACROSS the timeline. Nothing confronts
# the two, and "the WT background changes and that is why the Myc effect fades" has
# been asserted (Issue #6's WT-convergence term) without the controls it needs.
#
# WHAT THIS SCRIPT IS FOR. Three rulers, two axes:
#   content  = DESeq2 set-average RAW log2FC   (absolute transcriptional amount)
#   priority = mitoPPS pairwise diff           (relative allocation, content-blind)
#   share    = % of the nuclear transcriptome  (script 32's absolute compartment share)
# Read off disk (2026-07-24 post-reconciliation run), the picture is:
#   1. Myc's effect is RESCALED, not reshaped: Myc@12W ~ Myc@6W has slope 0.55
#      (content) / 0.65 (priority), R2 ~0.8, per-tier slopes 0.53-0.83.
#   2. The temporal move is largely SHARED: Myc+time ~ WTtime slope ~1 in the
#      import/central-dogma/OXPHOS arms, plus a constant negative intercept.
#   3. But the background is NOT moving toward the Myc state: projection of the WT
#      6->12 vector on the Myc(6W) vector = +0.06; cos angle -0.32 (content) /
#      -0.51 (priority) in exactly the arms that carry the claim.
#
# WHY CONTROLS ARE MANDATORY HERE. Every naive test of "the background shapes Myc"
# is algebraically rigged, because the contrasts share terms:
#   cor(WTtime, delta-Myc-effect) = -0.33 -- the value FORCED by the observed sds and
#     r(tneg,tpos)=0.75 is -0.33 exactly (delta = tpos - tneg).
#   cor(retention, WTtime) = -0.52 -- built in, since retention = 1 + (tpos-tneg)/m6.
#   cor(WT baseline, Myc effect) = -0.79/-0.88 -- shared-baseline + compositional.
#   ** Issue #6's WT-convergence carries the same exposure: myc_6W = pos6 - neg6 and
#      timepoint_neg = neg12 - neg6 share neg6 with OPPOSITE signs, so
#      cov = +var(neg6) > 0 -- noise in the 6W WT baseline MANUFACTURES apparent
#      convergence. The bias pushes frac_wt_toward UP, so the OXPHOS DIVERGENCE
#      (frac 0.31) is robust against it, but the biosynthetic convergence (0.74) and
#      the headline 34% convergence / 66% fade are the exposed half. **
# PART C is the fix: split the six 6W_neg mice so the genotype contrast and the
# temporal contrast no longer share a baseline.
#
# Input:  results/interaction_results.rds   (raw LFC contrasts + baseMean)
#         results/mitopps_scores.rds        ($mitopps_pairwise, $mitopps_group_means,
#                                            $gene_to_pathway, $pathway_tier1_map,
#                                            $pathway_levels, $mtdna_genes_separated)
#         results/dds_int_run.rds           (size factors, colData)
#         results/count_matrix.rds          (raw counts -> shares)
#         results/attenuation_mechanism.rds ($decomp_conv = the published Issue #6 split)
#         results/attenuation_decomposition.rds ($defs$luminal_sets for the pooled arms)
#         data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt
#         functions/reconcile_gene_symbols.R  (MANDATORY -- vintage-aware membership)
# Output: results/background_vs_myc.rds
#
# CEILING. Descriptive/exploratory at n=6/group. This is a DECOMPOSITION with nulls,
# not confirmatory inference; no CIs are offered as tests. BATCH = TIMEPOINT, so both
# temporal contrasts carry the same batch offset: it CANCELS in the interaction (the
# intercept of PART B's second regression) but it does NOT cancel in the
# convergence/fade split, and "developmental" remains an INTERPRETATION of the shared
# vector, never a measurement. The sample-split of PART C removes the shared-baseline
# ARTIFACT; it cannot remove the batch confound.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

set.seed(1)
NPERM <- 500L          # matched-random-set draws (PART B/D null); 500 is enough to rank
NBOOT <- 2000L         # pathway bootstrap for the regression CIs

group_levels <- c("6W_neg", "6W_pos", "12W_neg", "12W_pos")

# =============================================================================
# PART 1: LOAD
# =============================================================================
ir  <- readRDS(here::here("results", "interaction_results.rds"))
mp  <- readRDS(here::here("results", "mitopps_scores.rds"))
dds <- readRDS(here::here("results", "dds_int_run.rds"))
cts <- readRDS(here::here("results", "count_matrix.rds"))
am  <- readRDS(here::here("results", "attenuation_mechanism.rds"))
ad  <- readRDS(here::here("results", "attenuation_decomposition.rds"))
gmt <- fgsea::gmtPathways(
  here::here("data", "genesets_from_library", "mammary_mito_myc_metab_v1_mouse.gmt"))

sm <- as.data.frame(SummarizedExperiment::colData(dds))
sm$group <- factor(sm$group, levels = group_levels)
samples  <- colnames(dds)
cts      <- cts[, samples, drop = FALSE]
stopifnot(identical(colnames(cts), samples), all(table(sm$group) == 6))

universe_all <- rownames(ir$myc_6W_raw)
ens_of <- function(syms, universe = universe_all) recon_to_ensembl(syms, universe)

# raw (unshrunken) LFC vectors -- the identity below only holds on the MLE fit
V  <- function(k) stats::setNames(ir[[k]]$log2FoldChange, rownames(ir[[k]]))
m6  <- V("myc_6W_raw");      m12 <- V("myc_12W_raw")
tn  <- V("timepoint_neg_raw"); tp <- V("timepoint_pos_raw")
bm  <- stats::setNames(ir$myc_6W_raw$baseMean, rownames(ir$myc_6W_raw))
p6  <- stats::setNames(ir$myc_6W_raw$padj,     rownames(ir$myc_6W_raw))

# =============================================================================
# PART A: THE THREE-RULER PATHWAY TABLE
# =============================================================================
# One row per MitoPathway: the content ruler (set-average raw LFC) and the priority
# ruler (mitoPPS pairwise diff) for the same four contrasts, plus its Level-1 tier.
pwise <- mp$mitopps_pairwise
tier1 <- mp$pathway_tier1_map
gp    <- mp$gene_to_pathway
pcol  <- names(gp)[1]; gcol <- names(gp)[2]

paths <- intersect(unique(pwise$pathway), names(tier1))
path_ens <- lapply(stats::setNames(paths, paths), function(p)
  ens_of(unique(gp[[gcol]][gp[[pcol]] == p])))
message(sprintf("PART A: %d MitoPathways; median genes resolved = %.0f",
                length(paths), stats::median(vapply(path_ens, length, integer(1)))))

set_mean <- function(v) vapply(paths, function(p) {
  x <- v[path_ens[[p]]]
  if (length(x)) mean(x, na.rm = TRUE) else NA_real_
}, numeric(1))

grab <- function(k) {
  s <- pwise[pwise$contrast == k, ]
  list(diff = s$diff[match(paths, s$pathway)], padj = s$padj[match(paths, s$pathway)])
}
g_m6 <- grab("Myc_effect_6W"); g_m12 <- grab("Myc_effect_12W")
g_tn <- grab("Temporal_Myc-"); g_tp  <- grab("Temporal_Myc+")

ruler <- tibble::tibble(
  pathway = paths,
  tier    = unname(tier1[paths]),
  n_genes = vapply(path_ens, length, integer(1)),
  c_m6 = set_mean(m6), c_m12 = set_mean(m12),
  c_tn = set_mean(tn), c_tp  = set_mean(tp),
  p_m6 = g_m6$diff,    p_m12 = g_m12$diff,
  p_tn = g_tn$diff,    p_tp  = g_tp$diff,
  p_m6_padj = g_m6$padj, p_m12_padj = g_m12$padj) |>
  dplyr::mutate(c_int = c_m12 - c_m6, p_int = p_m12 - p_m6,
                is_mtdna = pathway == mp$mtdna_pathway_name) |>
  dplyr::filter(!is.na(c_m6), !is.na(p_m6))

# THE IDENTITY, on both rulers: (Myc@12W - Myc@6W) == (Myc+ temporal - WT temporal)
id_c <- max(abs(ruler$c_int - (ruler$c_tp - ruler$c_tn)))
id_p <- max(abs(ruler$p_int - (ruler$p_tp - ruler$p_tn)))
message(sprintf("PART A identity: content max|dev| = %.2e ; priority max|dev| = %.2e", id_c, id_p))
stopifnot(id_c < 1e-8, id_p < 1e-8)

tier_levels <- c("Protein import, sorting and homeostasis", "Mitochondrial central dogma",
                 "OXPHOS", "Metabolism", "Signaling",
                 "Mitochondrial dynamics and surveillance", "Small molecule transport")
stopifnot(all(ruler$tier %in% tier_levels))
bio_tiers <- tier_levels[1:3]              # the claim-bearing biogenesis + OXPHOS arms

# THE mtDNA PATHWAY IS EXCLUDED FROM EVERY FIT. mitoPPS keeps the 13 mtDNA-encoded
# subunits as their own synthetic pathway, tiered under OXPHOS. It is a 4-6.5 SD
# outlier on EVERY temporal contrast (mt% is three-way confounded: real content x
# proliferation denominator x prep leak) and on its own inflates the shared-vector
# slope from 0.75 to 0.84 and its R2 from 0.42 to 0.48. It stays in `ruler` (flagged
# `is_mtdna`) and is reported separately; it never enters a regression, a null or the
# geometry.
stopifnot(sum(ruler$is_mtdna) == 1L)
ruler_fit <- ruler[!ruler$is_mtdna, ]

ruler_summary <- dplyr::bind_rows(lapply(
  c("c_m6", "c_m12", "c_tn", "c_tp", "c_int", "p_m6", "p_m12", "p_tn", "p_tp", "p_int"),
  function(v) {
    x <- ruler[[v]]                       # bind first: `ruler` is masked inside tibble()
    tibble::tibble(
      ruler  = ifelse(startsWith(v, "c_"), "content", "priority"),
      metric = sub("^[cp]_", "", v),
      median = stats::median(x), mean = mean(x),
      sd     = stats::sd(x),     pct_up = 100 * mean(x > 0),
      q25    = unname(stats::quantile(x, 0.25)),
      q75    = unname(stats::quantile(x, 0.75)))
  }))

# =============================================================================
# PART B: THE TWO REGRESSIONS, WITH AN EMPIRICAL NULL
# =============================================================================
# Regression 1  m12 ~ m6    -- "rescaled, not reshaped": slope = the retained fraction.
# Regression 2  tp  ~ tn    -- "shared vector + Myc offset": slope = how much of the
#                              temporal move is common to both genotypes; the INTERCEPT
#                              is the batch-clean interaction.
fit2 <- function(y, x, label, scope) {
  ok <- !is.na(x) & !is.na(y)
  m  <- stats::lm(y[ok] ~ x[ok])
  tibble::tibble(model = label, scope = scope, n = sum(ok),
                 intercept = unname(stats::coef(m)[1]),
                 slope     = unname(stats::coef(m)[2]),
                 slope_se  = summary(m)$coefficients[2, 2],
                 r2        = summary(m)$r.squared,
                 rho       = stats::cor(x[ok], y[ok], method = "spearman"))
}

regressions <- dplyr::bind_rows(
  fit2(ruler_fit$c_m12, ruler_fit$c_m6, "rescale (Myc@12W ~ Myc@6W)",   "content, all"),
  fit2(ruler_fit$p_m12, ruler_fit$p_m6, "rescale (Myc@12W ~ Myc@6W)",   "priority, all"),
  fit2(ruler_fit$c_tp,  ruler_fit$c_tn, "shared (Myc+time ~ WTtime)",   "content, all"),
  fit2(ruler_fit$p_tp,  ruler_fit$p_tn, "shared (Myc+time ~ WTtime)",   "priority, all"),
  local({ i <- ruler_fit$tier %in% bio_tiers; dplyr::bind_rows(
    fit2(ruler_fit$c_m12[i], ruler_fit$c_m6[i], "rescale (Myc@12W ~ Myc@6W)", "content, biogenesis+OXPHOS"),
    fit2(ruler_fit$p_m12[i], ruler_fit$p_m6[i], "rescale (Myc@12W ~ Myc@6W)", "priority, biogenesis+OXPHOS"),
    fit2(ruler_fit$c_tp[i],  ruler_fit$c_tn[i], "shared (Myc+time ~ WTtime)", "content, biogenesis+OXPHOS"),
    fit2(ruler_fit$p_tp[i],  ruler_fit$p_tn[i], "shared (Myc+time ~ WTtime)", "priority, biogenesis+OXPHOS")) }),
  dplyr::bind_rows(lapply(tier_levels, function(t) {
    i <- ruler_fit$tier == t
    if (sum(i) < 5) return(NULL)
    dplyr::bind_rows(fit2(ruler_fit$c_m12[i], ruler_fit$c_m6[i], "rescale (Myc@12W ~ Myc@6W)", paste0("content, ", t)),
                     fit2(ruler_fit$p_m12[i], ruler_fit$p_m6[i], "rescale (Myc@12W ~ Myc@6W)", paste0("priority, ", t)))
  })))

# --- bootstrap CI over pathways (both rulers) --------------------------------
boot_slope <- function(y, x) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  s <- vapply(seq_len(NBOOT), function(b) {
    i <- sample.int(n, n, replace = TRUE)
    if (stats::sd(x[i]) == 0) return(NA_real_)
    unname(stats::coef(stats::lm(y[i] ~ x[i]))[2])
  }, numeric(1))
  stats::quantile(s, c(0.025, 0.975), na.rm = TRUE)
}
regression_boot <- tibble::tibble(
  model = c("rescale content", "rescale priority", "shared content", "shared priority"),
  lo = NA_real_, hi = NA_real_)
regression_boot[1, c("lo", "hi")] <- as.list(boot_slope(ruler_fit$c_m12, ruler_fit$c_m6))
regression_boot[2, c("lo", "hi")] <- as.list(boot_slope(ruler_fit$p_m12, ruler_fit$p_m6))
regression_boot[3, c("lo", "hi")] <- as.list(boot_slope(ruler_fit$c_tp,  ruler_fit$c_tn))
regression_boot[4, c("lo", "hi")] <- as.list(boot_slope(ruler_fit$p_tp,  ruler_fit$p_tn))

# --- matched-random-set null (CONTENT ruler only) ----------------------------
# The priority ruler cannot be nulled this way: mitoPPS is a pairwise ratio built
# ON the MitoCarta partition, so a random gene set has no mitoPPS. Its honest
# ambient is script 38's mitoPPS ceiling (~0.42), not a re-drawn set.
uni <- names(m6)[!is.na(m6) & !is.na(m12) & !is.na(tn) & !is.na(tp) &
                 !is.na(bm) & bm > 0]
dec_of <- stats::setNames(cut(rank(bm[uni], ties.method = "first"),
                              breaks = 10, labels = FALSE), uni)
pool   <- split(uni, dec_of[uni])
Mlfc   <- cbind(m6 = m6[uni], m12 = m12[uni], tn = tn[uni], tp = tp[uni])

path_uni <- lapply(path_ens[ruler_fit$pathway], function(e) intersect(e, uni))
path_dec <- lapply(path_uni, function(e) table(factor(dec_of[e], levels = names(pool))))

# NULL 1 (PRIMARY) -- within-decile LABEL SHUFFLE. One permutation of gene identities
# inside each expression decile, applied to all 144 pathways at once. Preserves set
# SIZE, the expression profile AND the pathway OVERLAP structure (MitoPathways nest:
# a Level-1 set contains its children, which by itself creates cross-pathway
# correlation); destroys only the gene-to-pathway assignment.
draw_shuffled <- function() {
  pm <- stats::setNames(uni, uni)
  for (k in names(pool)) { g <- pool[[k]]; pm[g] <- sample(g) }
  t(vapply(path_uni, function(e) colMeans(Mlfc[pm[e], , drop = FALSE]), numeric(4)))
}
# NULL 2 (secondary) -- expression-decile-matched RESAMPLING. Independent pseudo-sets,
# so it does NOT preserve overlap; kept as the more conservative comparator on slope.
draw_matched <- function() {
  t(vapply(path_dec, function(tb) {
    g <- unlist(lapply(names(tb)[tb > 0], function(k)
      sample(pool[[k]], tb[[k]], replace = TRUE)), use.names = FALSE)
    colMeans(Mlfc[g, , drop = FALSE])
  }, numeric(4)))
}

null_stats <- function(X) {
  r1 <- stats::lm(X[, "m12"] ~ X[, "m6"]); r2 <- stats::lm(X[, "tp"] ~ X[, "tn"])
  ca <- function(a, b) sum(a * b) / sqrt(sum(a^2) * sum(b^2))
  c(rescale_slope = unname(stats::coef(r1)[2]), rescale_r2 = summary(r1)$r.squared,
    shared_slope  = unname(stats::coef(r2)[2]), shared_r2  = summary(r2)$r.squared,
    shared_int    = unname(stats::coef(r2)[1]),
    cos_wt_myc    = ca(X[, "tn"], X[, "m6"]),
    proj_wt_myc   = sum(X[, "tn"] * X[, "m6"]) / sum(X[, "m6"]^2),
    cos_mycT_myc  = ca(X[, "tp"], X[, "m6"]))
}
null_draws     <- vapply(seq_len(NPERM), function(i) null_stats(draw_shuffled()), numeric(8))
null_draws_rs  <- vapply(seq_len(NPERM), function(i) null_stats(draw_matched()),  numeric(8))

pctile <- function(obs, nulls) 100 * mean(nulls <= obs, na.rm = TRUE)
null_summary <- tibble::tibble(
  statistic     = rownames(null_draws),
  null_median   = apply(null_draws, 1, stats::median, na.rm = TRUE),
  null_lo       = apply(null_draws, 1, stats::quantile, 0.025, na.rm = TRUE),
  null_hi       = apply(null_draws, 1, stats::quantile, 0.975, na.rm = TRUE),
  null_median_rs = apply(null_draws_rs, 1, stats::median, na.rm = TRUE),
  null_lo_rs     = apply(null_draws_rs, 1, stats::quantile, 0.025, na.rm = TRUE),
  null_hi_rs     = apply(null_draws_rs, 1, stats::quantile, 0.975, na.rm = TRUE))

# =============================================================================
# PART C: THE SAMPLE-SPLIT NULL FOR WT-CONVERGENCE  (the key new control)
# =============================================================================
# Issue #6 aligns each gene to d = sign(myc_6W) and splits the attenuation into
#   WT-convergence = mean(d * timepoint_neg)   and   Myc-fade = mean(d * timepoint_pos).
# myc_6W and timepoint_neg BOTH contain the 6W WT mean with opposite signs, so noise
# in that baseline inflates the convergence term. Here the six 6W_neg mice are split:
# three anchor the genotype contrast (they define d), the other three anchor the WT
# temporal contrast. All C(6,3) = 20 splits are enumerated; the complement pairing
# makes each split its own control.
nc  <- DESeq2::counts(dds, normalized = TRUE)[, samples, drop = FALSE]
lg  <- log2(nc + 1)
idx <- split(seq_along(samples), sm$group)
gm_of <- function(cols) rowMeans(lg[, cols, drop = FALSE])

# shared-baseline estimates from the SAME quantifier (the control: this is what the
# split is compared against, so the only thing that changes is the baseline sharing)
sh_m6  <- gm_of(idx[["6W_pos"]])  - gm_of(idx[["6W_neg"]])
sh_m12 <- gm_of(idx[["12W_pos"]]) - gm_of(idx[["12W_neg"]])
sh_tn  <- gm_of(idx[["12W_neg"]]) - gm_of(idx[["6W_neg"]])
sh_tp  <- gm_of(idx[["12W_pos"]]) - gm_of(idx[["6W_pos"]])
qc_agree <- stats::cor(sh_m6[intersect(names(sh_m6), names(m6))],
                       m6[intersect(names(sh_m6), names(m6))], use = "complete.obs")
message(sprintf("PART C: log-mean vs DESeq2 MLE agreement on myc_6W: r = %.3f", qc_agree))
stopifnot(qc_agree > 0.9)

# the arms Issue #6 reports on (script 31 roster), rebuilt here so the comparison is
# gene-for-gene the same partition
arm_roster <- tibble::tribble(
  ~program,                            ~arm,
  "MITOCARTA_OXPHOS_SUBUNITS",         "OXPHOS core",
  "MITOCARTA_OXPHOS",                  "OXPHOS core",
  "MITOCARTA_TCA_CYCLE",               "TCA",
  "MITOCARTA_NUCLEOTIDE_METABOLISM",   "nucleotide",
  "MITOCARTA_AMINO_ACID_METABOLISM",   "biosynthetic",
  "MITOCARTA_LIPID_METABOLISM",        "biosynthetic",
  "MITOCARTA_MITOCHONDRIAL_RIBOSOME",  "biogenesis/translation",
  "MYC_HALLMARK_MYC_TARGETS_V2",       "MYC-target core",
  "MYC_felsher_integrative_signature", "MYC-target core")
arm_ens <- lapply(stats::setNames(arm_roster$program, arm_roster$program),
                  function(p) ens_of(gmt[[p]]))
arm_ens[["PROLIFERATION_pooled"]] <-
  ens_of(unique(unlist(gmt[grep("^PROLIF_", names(gmt), value = TRUE)])))
arm_ens[["MAMMARY_LUMINAL_pooled"]] <-
  ens_of(unique(unlist(gmt[intersect(ad$defs$luminal_sets, names(gmt))])))
arm_roster <- dplyr::bind_rows(arm_roster, tibble::tibble(
  program = c("PROLIFERATION_pooled", "MAMMARY_LUMINAL_pooled"),
  arm     = c("proliferation", "mammary-dev")))

# one decomposition given a (possibly split) set of per-gene contrast vectors
decomp_split <- function(v_m6, v_tn, v_tp, genes) {
  e <- intersect(genes, names(v_m6))
  e <- e[!is.na(v_m6[e]) & !is.na(v_tn[e]) & !is.na(v_tp[e]) & v_m6[e] != 0]
  if (length(e) < 5) return(NULL)
  d  <- sign(v_m6[e])
  wc <- d * v_tn[e]; mf <- d * v_tp[e]
  atten <- mean(wc) - mean(mf)
  tibble::tibble(n = length(e), atten = atten,
                 wt_conv = mean(wc), myc_fade = mean(mf),
                 conv_pct = 100 * mean(wc) / atten,
                 fade_pct = 100 * (-mean(mf)) / atten,
                 frac_wt_toward = mean(wc > 0), frac_myc_retreat = mean(mf < 0))
}

# universe: effect-based and defined from the SAME baseline that defines d, so the
# selection cannot re-import the shared-baseline coupling
run_universe <- function(v_m6) names(v_m6)[!is.na(v_m6) & abs(v_m6) > 0.5]

one_run <- function(v_m6, v_tn, v_tp, tag, split_id) {
  u  <- run_universe(v_m6)
  g0 <- decomp_split(v_m6, v_tn, v_tp, u)
  stopifnot(!is.null(g0))
  rows <- list(dplyr::bind_cols(tibble::tibble(program = "ALL", arm = "global"), g0))
  rows <- c(rows, lapply(seq_len(nrow(arm_roster)), function(i) {
    r <- decomp_split(v_m6, v_tn, v_tp, intersect(arm_ens[[arm_roster$program[i]]], u))
    if (is.null(r)) return(NULL)
    dplyr::bind_cols(tibble::tibble(program = arm_roster$program[i],
                                    arm = arm_roster$arm[i]), r)
  }))
  dplyr::bind_rows(rows) |> dplyr::mutate(mode = tag, split = split_id, .before = 1)
}

shared_run <- one_run(sh_m6, sh_tn, sh_tp, "shared_baseline", 0L)

# THE MATCHED CONTROL. Splitting the baseline halves it (3 mice, not 6), which adds
# noise to d = sign(myc_6W). Sign noise pulls frac_wt_toward toward 0.5 from BOTH
# sides, so a "shared 6 vs split 3" comparison confounds the sharing artifact with
# plain attenuation-to-chance. Each split therefore runs TWICE off the SAME half
# baseline A -- once with A anchoring both contrasts (shared) and once with B
# anchoring the temporal contrast (split). Identical d, identical universe,
# identical noise level; the ONLY difference is whether the baseline is shared.
splits <- utils::combn(idx[["6W_neg"]], 3, simplify = FALSE)   # 20 configurations
split_runs <- dplyr::bind_rows(lapply(seq_along(splits), function(k) {
  a  <- splits[[k]]                      # anchors the GENOTYPE contrast (defines d)
  b  <- setdiff(idx[["6W_neg"]], a)      # anchors the WT TEMPORAL contrast (split only)
  m6a <- gm_of(idx[["6W_pos"]]) - gm_of(a)
  dplyr::bind_rows(
    one_run(m6a, gm_of(idx[["12W_neg"]]) - gm_of(a), sh_tp, "half_shared", k),
    one_run(m6a, gm_of(idx[["12W_neg"]]) - gm_of(b), sh_tp, "half_split",  k))
}))

# NOTE ON AGGREGATION. conv_pct = 100 * wt_conv / atten is a RATIO of two noisy
# quantities: on individual splits `atten` can pass near zero and the percentage
# explodes (per-split sd reaches 4 digits). The split estimate is therefore a RATIO
# OF MEANS (mean wt_conv over mean atten), never a mean of ratios, and the
# distribution-free summary `frac_wt_toward` -- the fraction of genes whose WT
# temporal change points toward Myc, with 0.5 = chance -- is the primary read.
agg <- function(mode_tag, suffix) split_runs |>
  dplyr::filter(mode == mode_tag) |>
  dplyr::group_by(program, arm) |>
  dplyr::summarise(n_splits    = dplyr::n(),
                   conv_pct    = 100 * mean(wt_conv) / mean(atten),
                   frac_toward = mean(frac_wt_toward),
                   frac_sd     = stats::sd(frac_wt_toward),
                   frac_min    = min(frac_wt_toward),
                   frac_max    = max(frac_wt_toward),
                   wt_conv     = mean(wt_conv),
                   wt_conv_sd  = stats::sd(wt_conv),
                   myc_fade    = mean(myc_fade),
                   atten       = mean(atten),
                   n_genes     = mean(n), .groups = "drop") |>
  dplyr::rename_with(~ paste0(.x, suffix), -c(program, arm))

# THE PAIRED DIFFERENCE is the estimate of the artifact: within each split, the same
# half baseline, shared vs split. Its sd across the 20 splits is a real uncertainty
# on the artifact size (unlike a comparison against the 6-mouse shared run, which
# also changes the noise level).
paired <- split_runs |>
  dplyr::select(program, arm, mode, split, frac_wt_toward, wt_conv) |>
  tidyr::pivot_wider(names_from = mode, values_from = c(frac_wt_toward, wt_conv)) |>
  dplyr::group_by(program, arm) |>
  dplyr::summarise(
    frac_artifact      = mean(frac_wt_toward_half_shared - frac_wt_toward_half_split),
    frac_artifact_sd   = stats::sd(frac_wt_toward_half_shared - frac_wt_toward_half_split),
    wt_conv_artifact   = mean(wt_conv_half_shared - wt_conv_half_split),
    frac_artifact_frac_pos = mean(frac_wt_toward_half_shared > frac_wt_toward_half_split),
    .groups = "drop")

split_summary <- agg("half_split", "_split") |>
  dplyr::left_join(agg("half_shared", "_halfshared"), by = c("program", "arm")) |>
  dplyr::left_join(paired, by = c("program", "arm")) |>
  dplyr::left_join(
    shared_run |> dplyr::select(program, arm,
                                conv_pct_shared = conv_pct,
                                frac_toward_shared = frac_wt_toward,
                                wt_conv_shared = wt_conv,
                                atten_shared = atten, n_genes_shared = n),
    by = c("program", "arm")) |>
  # No z or p is offered against the 0.5 chance line: genes within an arm are strongly
  # correlated, so a binomial SE would be anti-conservative, and the split-to-split sd
  # measures split variability, not sampling error. The range across the 20 splits is
  # reported instead and read descriptively.
  dplyr::mutate(conv_pct_bias    = conv_pct_shared - conv_pct_split,
                frac_toward_bias = frac_toward_shared - frac_toward_split) |>
  dplyr::left_join(                                  # the PUBLISHED Issue #6 numbers
    am$decomp_conv |> dplyr::filter(universe == "effect") |>
      dplyr::select(program, conv_pct_published = conv_pct,
                    frac_toward_published = frac_wt_toward),
    by = "program")

# =============================================================================
# PART D: THE GEOMETRY (is the background moving TOWARD the Myc state?)
# =============================================================================
cosang <- function(a, b) sum(a * b, na.rm = TRUE) /
  sqrt(sum(a^2, na.rm = TRUE) * sum(b^2, na.rm = TRUE))
projn  <- function(a, b) sum(a * b, na.rm = TRUE) / sum(b^2, na.rm = TRUE)

geom_one <- function(d, scope) tibble::tibble(
  scope = scope, n = nrow(d),
  c_cos_wt_myc   = cosang(d$c_tn, d$c_m6), c_proj_wt_myc  = projn(d$c_tn, d$c_m6),
  p_cos_wt_myc   = cosang(d$p_tn, d$p_m6), p_proj_wt_myc  = projn(d$p_tn, d$p_m6),
  c_cos_mycT_myc = cosang(d$c_tp, d$c_m6), p_cos_mycT_myc = cosang(d$p_tp, d$p_m6),
  c_sd_ratio_bg  = stats::sd(d$c_tn) / stats::sd(d$c_m6),
  p_sd_ratio_bg  = stats::sd(d$p_tn) / stats::sd(d$p_m6))

geometry <- dplyr::bind_rows(
  geom_one(ruler_fit, "all pathways"),
  geom_one(ruler_fit[ruler_fit$tier %in% bio_tiers, ], "biogenesis + OXPHOS"),
  geom_one(ruler_fit[!ruler_fit$tier %in% bio_tiers, ], "metabolism + other"),
  dplyr::bind_rows(lapply(tier_levels, function(t) {
    d <- ruler_fit[ruler_fit$tier == t, ]; if (nrow(d) < 5) NULL else geom_one(d, t) })))

geometry_null <- tibble::tibble(
  statistic = c("cos_wt_myc", "proj_wt_myc", "cos_mycT_myc"),
  observed  = c(geometry$c_cos_wt_myc[1], geometry$c_proj_wt_myc[1],
                geometry$c_cos_mycT_myc[1]),
  null_median = c(stats::median(null_draws["cos_wt_myc", ]),
                  stats::median(null_draws["proj_wt_myc", ]),
                  stats::median(null_draws["cos_mycT_myc", ])),
  percentile = c(pctile(geometry$c_cos_wt_myc[1],  null_draws["cos_wt_myc", ]),
                 pctile(geometry$c_proj_wt_myc[1], null_draws["proj_wt_myc", ]),
                 pctile(geometry$c_cos_mycT_myc[1], null_draws["cos_mycT_myc", ])),
  percentile_rs = c(pctile(geometry$c_cos_wt_myc[1],  null_draws_rs["cos_wt_myc", ]),
                    pctile(geometry$c_proj_wt_myc[1], null_draws_rs["proj_wt_myc", ]),
                    pctile(geometry$c_cos_mycT_myc[1], null_draws_rs["cos_mycT_myc", ])))

regression_null <- tibble::tibble(
  statistic   = c("rescale_slope", "rescale_r2", "shared_slope", "shared_r2", "shared_int"),
  observed    = c(regressions$slope[regressions$scope == "content, all" &
                                      grepl("^rescale", regressions$model)],
                  regressions$r2[regressions$scope == "content, all" &
                                   grepl("^rescale", regressions$model)],
                  regressions$slope[regressions$scope == "content, all" &
                                      grepl("^shared", regressions$model)],
                  regressions$r2[regressions$scope == "content, all" &
                                   grepl("^shared", regressions$model)],
                  regressions$intercept[regressions$scope == "content, all" &
                                          grepl("^shared", regressions$model)])) |>
  dplyr::mutate(
    null_median = vapply(statistic, function(s) stats::median(null_draws[s, ]), numeric(1)),
    percentile  = vapply(seq_along(statistic),
                         function(i) pctile(observed[i], null_draws[statistic[i], ]), numeric(1)),
    null_median_rs = vapply(statistic, function(s) stats::median(null_draws_rs[s, ]), numeric(1)),
    percentile_rs  = vapply(seq_along(statistic),
                            function(i) pctile(observed[i], null_draws_rs[statistic[i], ]), numeric(1)))

# =============================================================================
# PART E: RULER RECONCILIATION (how big IS the attenuation?)
# =============================================================================
# Same gene sets, four quantifiers of the genotype gap at each age. The share ruler
# and the DESeq LFC ruler disagree on the SIZE of the attenuation (not on its sign):
# they differ in denominator (per-sample transcriptome sum vs median-of-ratios size
# factors) and in weighting (expression-weighted sum vs unweighted per-gene mean).
mt_ens   <- ens_of(mp$mtdna_genes_separated, rownames(cts))
den_nomt <- colSums(cts[setdiff(rownames(cts), mt_ens), , drop = FALSE])
share_of <- function(e) 100 * colSums(cts[intersect(e, rownames(cts)), , drop = FALSE]) / den_nomt

gene_z <- {                                    # script 36's linear quantifier (36:113-125)
  z <- t(scale(t(lg)))
  z[is.finite(rowSums(z)), , drop = FALSE]
}
gap_by <- function(v) {                        # log2 genotype gap at each age from a per-sample vector
  c(g6  = mean(v[idx[["6W_pos"]]])  - mean(v[idx[["6W_neg"]]]),
    g12 = mean(v[idx[["12W_pos"]]]) - mean(v[idx[["12W_neg"]]]))
}

recon_sets <- c("MITOCARTA_NUCLEAR_ENCODED", "MITOCARTA_OXPHOS_NU",
                "MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA",
                "MITOCARTA_PROTEIN_IMPORT_SORTING_AND_HOMEOSTASIS")
ruler_reconciliation <- dplyr::bind_rows(lapply(recon_sets, function(s) {
  e  <- ens_of(gmt[[s]])
  ec <- intersect(e, rownames(cts)); ez <- intersect(e, rownames(gene_z))
  w  <- bm[intersect(e, names(bm))]; w[is.na(w)] <- 0
  q  <- list(
    `share of nuclear transcriptome` = gap_by(log2(share_of(ec))),
    `set-average LFC (unweighted)`   = c(g6 = mean(m6[e], na.rm = TRUE),
                                         g12 = mean(m12[e], na.rm = TRUE)),
    `set-average LFC (expr-weighted)`= c(g6 = stats::weighted.mean(m6[names(w)], w, na.rm = TRUE),
                                         g12 = stats::weighted.mean(m12[names(w)], w, na.rm = TRUE)),
    `mean gene-wise z-score`         = gap_by(colMeans(gene_z[ez, , drop = FALSE])))
  dplyr::bind_rows(lapply(names(q), function(k) tibble::tibble(
    set = s, quantifier = k, n_genes = length(e),
    gap_6W = unname(q[[k]]["g6"]), gap_12W = unname(q[[k]]["g12"]),
    attenuation = unname(q[[k]]["g12"] - q[[k]]["g6"]),
    pct_retained = 100 * unname(q[[k]]["g12"] / q[[k]]["g6"]))))
}))

# =============================================================================
# PART F: THE STATE TABLE (per Level-1 tier x group) -- what fig03 panel A draws
# =============================================================================
# Two absolute STATES, not contrasts: priority = mitoPPS group mean (centred at 1),
# content = share of the nuclear transcriptome. Tier membership comes from the mitoPPS
# partition itself, so the two axes describe the same gene sets.
# The 13 mtDNA-encoded genes are dropped from the NUMERATOR as well as the denominator.
# They sit inside the OXPHOS tier, are 3-40% of the library, and swing wildly between
# samples (mt% is three-way confounded: content x proliferation denominator x prep leak),
# so leaving them in makes the OXPHOS content axis an mt% readout (sd 23 vs 5 on a mean
# of 10). Every nuclear arm in this project is read mt-free; the tiers follow.
tier_genes <- lapply(stats::setNames(tier_levels, tier_levels), function(t) {
  ps <- names(tier1)[tier1 == t]
  setdiff(ens_of(unique(gp[[gcol]][gp[[pcol]] %in% ps]), rownames(cts)), mt_ens)
})
message("PART F tier gene counts: ",
        paste(sprintf("%s=%d", names(tier_genes), lengths(tier_genes)), collapse = ", "))
stopifnot(all(lengths(tier_genes) > 20))

state_samples <- dplyr::bind_rows(lapply(tier_levels, function(t) tibble::tibble(
  tier = t, sample = samples, group = sm$group,
  timepoint = sm$timepoint, myc_status = sm$myc_status,
  content_share = share_of(tier_genes[[t]]))))

l1_paths <- mp$pathway_levels$Pathway[mp$pathway_levels$Level == "Pathway_Level1"]
gmeans   <- mp$mitopps_group_means
prio <- gmeans[gmeans$pathway %in% l1_paths, ]
prio$tier <- unname(tier1[prio$pathway])
stopifnot(all(tier_levels %in% prio$tier))

state_table <- state_samples |>
  dplyr::group_by(tier, group) |>
  dplyr::summarise(content_sd    = stats::sd(content_share),
                   content_share = mean(content_share), .groups = "drop") |>
  dplyr::left_join(prio[, c("tier", "group", "mean_mitopps", "sd_mitopps")],
                   by = c("tier", "group")) |>
  dplyr::rename(priority = mean_mitopps, priority_sd = sd_mitopps) |>
  dplyr::mutate(group = factor(group, levels = group_levels)) |>
  dplyr::arrange(tier, group)

# =============================================================================
# PART G: THE ARTIFACT LEDGER
# =============================================================================
# For every correlation this section could quote, the value the ALGEBRA forces given
# the observed marginals, next to the observed value. A statistic whose forced value
# equals its observed value carries NO evidence and must never be cited.
forced_cor_diff <- function(a, b) {          # cor(a, b - a) when b - a is the difference
  sa <- stats::sd(a); sb <- stats::sd(b); r <- stats::cor(a, b)
  (r * sa * sb - sa^2) / (sa * stats::sd(b - a))
}
artifact_ledger <- tibble::tibble(
  statistic = c("cor(WTtime, delta-Myc-effect) [content]",
                "cor(WTtime, delta-Myc-effect) [priority]",
                "cor(retention, WTtime) [content]",
                "cor(WT baseline @6W, Myc effect @6W) [priority]",
                "cor(WT baseline @12W, Myc effect @12W) [priority]"),
  observed = c(
    stats::cor(ruler$c_tn, ruler$c_int),
    stats::cor(ruler$p_tn, ruler$p_int),
    local({ k <- ruler$c_m6 > 0.15
            stats::cor(ruler$c_m12[k] / ruler$c_m6[k], ruler$c_tn[k]) }),
    local({ w <- gmeans[gmeans$group == "6W_neg", ]
            stats::cor(w$mean_mitopps[match(ruler$pathway, w$pathway)], ruler$p_m6,
                       use = "complete.obs") }),
    local({ w <- gmeans[gmeans$group == "12W_neg", ]
            stats::cor(w$mean_mitopps[match(ruler$pathway, w$pathway)], ruler$p_m12,
                       use = "complete.obs") })),
  forced = c(forced_cor_diff(ruler$c_tn, ruler$c_tp),
             forced_cor_diff(ruler$p_tn, ruler$p_tp),
             NA_real_, NA_real_, NA_real_),
  mechanism = c(
    "delta = tpos - tneg : tneg enters with a negative sign",
    "delta = tpos - tneg : tneg enters with a negative sign",
    "retention = 1 + (tpos - tneg)/m6 : tneg is inside the statistic",
    "mitoPPS is compositional (centred at 1) and the effect subtracts the baseline",
    "mitoPPS is compositional (centred at 1) and the effect subtracts the baseline"),
  verdict = c("STRUCTURAL -- do not cite", "STRUCTURAL -- do not cite",
              "STRUCTURAL -- do not cite", "STRUCTURAL -- do not cite",
              "STRUCTURAL -- do not cite"))

# =============================================================================
# SAVE
# =============================================================================
notes <- paste(
  "Script 40. Is Myc's mitochondrial effect shaped by a changing background?",
  "",
  "PART A. Three rulers, 144 MitoPathways. The identity (Myc@12W - Myc@6W) ==",
  "  (Myc+ temporal - WT temporal) holds exactly on BOTH the content and the priority",
  "  ruler, so the two decompositions below are of the same quantity. The synthetic",
  "  mtDNA-encoded pathway is flagged `is_mtdna` and EXCLUDED from every fit, null and",
  "  geometry: it is a 4-6.5 SD outlier on every temporal contrast and alone moves the",
  "  shared-vector slope 0.75 -> 0.84. Its row stays in `ruler`, reported separately.",
  "PART B. Regression 1 (rescale): the 12W Myc profile IS the 6W profile times a",
  "  constant -- read `regressions` for the slope and `regression_null` for where it",
  "  sits against expression-matched random sets. Regression 2 (shared): the slope is",
  "  the shared temporal component, the INTERCEPT is the batch-clean interaction.",
  "  TWO nulls, both CONTENT-ruler only (mitoPPS has no random-set analogue -- it is",
  "  a ratio built on the MitoCarta partition; its ambient is script 38's ~0.42):",
  "  PRIMARY = within-decile LABEL SHUFFLE, which preserves set size, expression and",
  "  the pathway OVERLAP structure (MitoPathways nest, and that nesting alone creates",
  "  cross-pathway correlation); SECONDARY (`_rs` columns) = decile-matched",
  "  resampling, which does not preserve overlap.",
  "PART C. THE CONTROL THAT MATTERS. Issue #6's WT-convergence shares the 6W WT",
  "  baseline with the genotype contrast, which manufactures apparent convergence.",
  "  `split_summary` recomputes it over all 20 three-versus-three splits of the",
  "  6W_neg mice. `conv_pct_shared` is the same-quantifier control (log-mean, no",
  "  split) and tracks `conv_pct_published` (script 31, effect universe). But a",
  "  6-mouse shared baseline vs a 3-mouse split one ALSO changes the noise in",
  "  d = sign(myc_6W), and sign noise pulls frac_wt_toward toward 0.5 from BOTH",
  "  sides -- so the artifact is estimated from the MATCHED PAIR instead: the same",
  "  half baseline A, once shared between both contrasts (`_halfshared`) and once",
  "  split (`_split`). `frac_artifact` = the paired difference, with its across-split",
  "  sd and the fraction of splits in which it has the expected sign.",
  "  PRIMARY READ = `frac_wt_toward`, the",
  "  fraction of genes whose WT temporal change points toward Myc (0.5 = chance);",
  "  conv_pct is a ratio of noisy quantities and is reported as a ratio of means,",
  "  never a mean of ratios. No p is offered against 0.5 -- genes within an arm are",
  "  correlated -- so the 20-split range is what is read.",
  "PART D. Geometry: cos and projection of the WT 6->12 vector on the Myc(6W) vector,",
  "  whole compartment and per tier, against the same null.",
  "PART E. The attenuation's SIZE is normalisation-dependent (share vs DESeq LFC vs",
  "  expression-weighted LFC vs gene-wise z). Quote the range, not one number.",
  "PART F. `state_table` = the 7 Level-1 tiers x 4 groups in (priority, content)",
  "  space; fig03 panel A draws the WT 6->12, Myc+ 6->12 and genotype arrows on it.",
  "PART G. `artifact_ledger`: statistics whose observed value equals the value the",
  "  algebra forces. They are not evidence. Do not cite them. It is computed on the",
  "  FULL 144-pathway table on purpose -- these are the numbers a reader would get.",
  "",
  "SCOPE. Descriptive/exploratory, n=6/group. BATCH = TIMEPOINT: the batch offset",
  "cancels in the interaction (regression 2's intercept) but NOT in the",
  "convergence/fade split, and 'developmental' stays an interpretation of the shared",
  "vector. The split removes the shared-baseline ARTIFACT, not the batch confound.",
  sep = "\n")

out <- list(
  ruler = ruler, ruler_summary = ruler_summary,
  regressions = regressions, regression_boot = regression_boot,
  regression_null = regression_null, null_summary = null_summary,
  null_draws = null_draws, null_draws_rs = null_draws_rs,
  shared_run = shared_run, split_runs = split_runs, split_summary = split_summary,
  geometry = geometry, geometry_null = geometry_null,
  ruler_reconciliation = ruler_reconciliation,
  state_table = state_table, state_samples = state_samples,
  artifact_ledger = artifact_ledger,
  defs = list(tier_levels = tier_levels, bio_tiers = bio_tiers,
              path_ens = path_ens, tier_genes = tier_genes, arm_ens = arm_ens,
              n_perm = NPERM, n_boot = NBOOT, n_splits = length(splits),
              logmean_vs_mle_cor = qc_agree, identity_dev = c(content = id_c, priority = id_p)),
  analysis_date = Sys.Date(), notes = notes)

saveRDS(out, here::here("results", "background_vs_myc.rds"))
message("wrote results/background_vs_myc.rds")

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  bm40 <- readRDS(here::here("results", "background_vs_myc.rds"))

  ## the three-ruler summary (the table the paper quotes)
  bm40$ruler_summary |> print(n = nrow(bm40$ruler_summary))

  ## 1. rescaled, not reshaped -- and 2. shared vector + Myc offset
  bm40$regressions |> dplyr::filter(scope %in% c("content, all", "priority, all",
                                                 "content, biogenesis+OXPHOS",
                                                 "priority, biogenesis+OXPHOS")) |> print()
  bm40$regression_boot |> print()
  bm40$regression_null |> print()

  ## 3. THE CONTROL: does WT-convergence survive the sample split?
  ##    PRIMARY READ = frac_wt_toward (0.5 = chance). The percentage is a ratio of
  ##    noisy quantities; conv_pct_split is a ratio of means, read as secondary.
  bm40$split_summary |>
    dplyr::select(program, arm, frac_toward_published, frac_toward_shared,
                  frac_toward_halfshared, frac_toward_split,
                  frac_artifact, frac_artifact_sd, frac_artifact_frac_pos) |>
    dplyr::arrange(frac_toward_split) |> print(n = 20)
  bm40$split_summary |>
    dplyr::select(program, arm, conv_pct_published, conv_pct_shared,
                  conv_pct_halfshared, conv_pct_split,
                  wt_conv_halfshared, wt_conv_split, wt_conv_artifact) |>
    print(n = 20)

  ## 4. geometry: is the background moving toward the Myc state?
  bm40$geometry |> print(n = nrow(bm40$geometry))
  bm40$geometry_null |> print()

  ## 5. how big is the attenuation, by ruler?
  bm40$ruler_reconciliation |> print(n = nrow(bm40$ruler_reconciliation))

  ## 6. the state table fig03 panel A draws
  bm40$state_table |> print(n = 28)

  ## 7. the artifact ledger -- observed vs the value the algebra forces
  bm40$artifact_ledger |> print()

  cat(bm40$notes, "\n")
}
