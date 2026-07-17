# scripts/34_death_priming_reassessment.R
# =============================================================================
# Does the 6W death priming actually COLLAPSE by 12W -- and is that collapse
# death-specific, or just the Myc programme fading? (Block B, author challenge
# 2026-07-17)
# =============================================================================
#
# WHY THIS EXISTS. Script 33 PART C showed the per-sample death couplings fail a
# specificity null. The author's challenge: that null is computed ENTIRELY WITHIN
# 6W (`33:278-294`, idx = i6), so it tests a CROSS-SECTIONAL claim ("at 6W, the
# more imbalanced animals are the more death-primed ones"). But the manuscript's
# actual claim is LONGITUDINAL: priming is raised at 6W and COLLAPSES by 12W, and
# that collapse tracks mitochondrial reprioritisation. The 6W null speaks to
# neither. This script builds the nulls that do.
#
# THE CORRECTED NULL DESIGN -- three levels, replacing the single 6W null:
#
#   L1 (PART A) AMBIENT-CORRELATION DIAGNOSTIC. Within 6W, PC1 = 43% of variance
#      and everything loads on it. PC1 WITHIN 12W WAS NEVER COMPUTED. If the
#      ambient correlation is weaker at 12W then EVERY set's coupling collapses
#      6W->12W and the imbalance's collapse means nothing. The 6W null median is
#      +0.62; if the 12W null median is ~0.15, the "collapse" IS the ambient
#      structure. This is the cheapest possible test of the headline claim.
#
#   L2 (PART B) THE DELTA-RHO NULL -- the null the collapse claim always needed.
#      For each of the 884 library sets: d_s = rho_6W(X, s) - rho_12W(X, s).
#      Locate d_obs = rho_6W(X, pro_comp) - rho_12W(X, pro_comp) in that
#      distribution. Corroborated by a Fisher z test of r6 vs r12.
#
#   L3 (PART C) THE INTERACTION-PERCENTILE NULL. Does the death programme
#      attenuate MORE than an arbitrary programme? Script 17's coef_table already
#      holds beta_int for all 884 sets; locate the apoptosis sets in it.
#
# A STRUCTURAL LIMIT, STATED NOT TESTED. The design is CROSS-SECTIONAL: different
# animals at 6W and 12W. There is NO within-animal delta. So "the priming collapse
# correlates with mito reprioritisation" is NOT estimable per-sample at all -- it
# can only ever be two 4-point group trajectories co-moving, which carries no
# inference. What IS estimable is whether the two share a genotype x time
# interaction. Both are null (priming p=0.263, imbalance p=0.916). See
# `cross_sectional_limit` in PART C.
#
# THREE PROBLEMS WITH THE DEATH ARM THAT PREDATE THIS SCRIPT:
#
#   (1) THE INTERACTION WAS FITTED AND SAVED ALL ALONG, AND NOBODY READ IT.
#       `death_timing_substrate.rds$h1$state_tests` carries int_p = 0.177 (PRO) /
#       0.263 (priming). The narrative quoted the GROUP MEANS instead ("~3x more
#       at 6W" = +0.48 vs +0.15). This is worse than the mitonuclear-imbalance
#       case, where the contrast was never fitted at all -- here it was fitted,
#       saved, and stepped over. PART C surfaces it.
#
#   (2) THE COUPLING IS CIRCULAR. `pro_comp` is MITOCARTA_APOPTOSIS_PRO -- 25
#       hand-curated MITOCHONDRIAL apoptosis genes (08_mitoPPS_analysis.R:285-290).
#       The imbalance is a mitoPPS ratio. "The mitochondrial state couples to death
#       priming" is mitochondrial genes correlating with mitochondrial genes. This
#       was flagged in docs/2026-07-08_BlockA_revision_plan.md:236 ("priming
#       mito-defined ... partly circular") and never applied to the death spine.
#       PART F breaks it.
#
#   (3) GATE 2's p=0.023 IS A SELF-CONTAINED TEST WHERE A COMPETITIVE ONE IS
#       NEEDED. `14:91` runs a one-sample t-test on 23 genes' interaction LFCs
#       against mu=0, treating genes as independent draws. Script 17 MEASURED the
#       inter-gene correlation in these sets at 0.25-0.36. And under a GLOBAL
#       attenuation (Issue #4: median |LFC| 0.268 -> 0.204), EVERY Myc-induced
#       module has a negative mean interaction LFC -- so mu=0 is the wrong null.
#       PART E replaces it with a matched null + the correlation-corrected n.
#
# NOTE ON THE ALIGNMENT BIAS (PART D). The Issue #6 identity aligns genes by
# d = sign(m6), which makes gap6 positive BY CONSTRUCTION, so `atten` is biased
# positive by regression-to-the-mean alone. That is precisely why PART D uses a
# MATCHED null: the null genes receive the SAME alignment, so the bias cancels.
# Do not read a positive `atten` as attenuation without its matched null.
#
# ATTRIBUTION NOTE. r=0.75/p=0.005 is 23_death_timing_substrate.R:326-327
# (PEARSON), not script 25. Script 25 Part B (25:140-148) uses stats::cor
# (Spearman, NO p-value) and reports 0.671. Script 33's header credits it to 25.
# Because the published claim is Pearson, PARTs A/B run BOTH methods.
#
# Input:  results/developmental_substrate_death.rds  (master per-sample frame)
#         results/death_timing_substrate.rds         (h1$state_tests, h2$coupling)
#         results/gsva_scores.rds                    (884 sets = the null universe)
#         results/gsva_overview.rds                  (coef_table: beta_int, camera, roast)
#         results/interaction_results.rds            (raw contrasts)
#         results/interaction_gene_characterisation.rds (Gate 2's 23 PRO genes)
#         results/attenuation_mechanism.rds          (Issue #6 identity defs)
#         results/cell_death_fgsea.rds               (Tang 15 RCD, 5 contrasts)
#         results/dds_int_run.rds, results/combined_df_annotated.rds
#         data/cell_death_genes_consolidated.csv     (512 pro-death, is_mitochondrial)
# Output: results/death_priming_reassessment.rds
#         outputs/death_priming_reassessment/*.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

out_dir <- here::here("outputs", "death_priming_reassessment")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

group_levels <- c("6W_neg", "6W_pos", "12W_neg", "12W_pos")
set.seed(1)

# =============================================================================
# PART 1: LOAD + ALIGN
# =============================================================================
ds  <- readRDS(here::here("results", "developmental_substrate_death.rds"))$per_sample
dt  <- readRDS(here::here("results", "death_timing_substrate.rds"))
gs  <- readRDS(here::here("results", "gsva_scores.rds"))
ov  <- readRDS(here::here("results", "gsva_overview.rds"))
ir  <- readRDS(here::here("results", "interaction_results.rds"))
ch  <- readRDS(here::here("results", "interaction_gene_characterisation.rds"))
cdf <- readRDS(here::here("results", "combined_df_annotated.rds"))
cdg <- readr::read_csv(here::here("data", "cell_death_genes_consolidated.csv"),
                       show_col_types = FALSE)

samples <- ds$sample
scores  <- gs$scores[, samples, drop = FALSE]
stopifnot(identical(colnames(scores), samples))

i6  <- ds$timepoint == "6W"
i12 <- ds$timepoint == "12W"
stopifnot(sum(i6) == 12, sum(i12) == 12)

sym2ens <- cdf[!is.na(cdf$mgi_symbol) & !duplicated(cdf$mgi_symbol), c("mgi_symbol", "gene")]
universe_all <- rownames(ir$myc_6W_raw)
ens_of <- function(syms) {
  e <- sym2ens$gene[match(intersect(syms, sym2ens$mgi_symbol), sym2ens$mgi_symbol)]
  intersect(e[!is.na(e)], universe_all)
}

# raw (unshrunken) LFC vectors -- the identity and every dLFC metric need raw
V   <- function(k) stats::setNames(ir[[k]]$log2FoldChange, rownames(ir[[k]]))
m6  <- V("myc_6W_raw"); m12 <- V("myc_12W_raw")
tn  <- V("timepoint_neg_raw"); tp <- V("timepoint_pos_raw")
lint <- V("interaction_raw")

# the axes under test (script 23's published per-sample values)
imbalance <- ds$mitonuclear_imbalance
pro_comp  <- ds$pro_comp

# =============================================================================
# PART A: THE AMBIENT-CORRELATION DIAGNOSTIC (L1)
# =============================================================================
# If the 12W cohort simply has less shared variance, every coupling collapses and
# the imbalance's collapse is not a finding. 33's D2 computed PC1 within 6W only.
vst_mat <- SummarizedExperiment::assay(
  DESeq2::vst(readRDS(here::here("results", "dds_int_run.rds")), blind = TRUE))
vst_mat <- vst_mat[, samples, drop = FALSE]

pc_one <- function(idx, label) {
  v <- vst_mat[order(-matrixStats::rowVars(vst_mat[, idx, drop = FALSE]))[1:2000],
               idx, drop = FALSE]
  s <- summary(stats::prcomp(t(v), scale. = TRUE))$importance[2, 1:3]
  tibble::tibble(window = label, pc1_pct = 100 * s[1], pc2_pct = 100 * s[2],
                 pc3_pct = 100 * s[3])
}
pc1_by_timepoint <- dplyr::bind_rows(pc_one(i6, "within 6W"), pc_one(i12, "within 12W"))

# the null universe's own coupling level, per timepoint. The KEY comparison.
null_level <- function(x, idx, label, method = "pearson") {
  r <- apply(scores[, idx, drop = FALSE], 1,
             function(s) suppressWarnings(stats::cor(s, x[idx], method = method)))
  r <- r[is.finite(r)]
  tibble::tibble(window = label, method = method, n_sets = length(r),
                 null_median_abs = stats::median(abs(r)),
                 null_median = stats::median(r),
                 frac_above_0.5 = mean(abs(r) > 0.5))
}
ambient_by_timepoint <- dplyr::bind_rows(
  null_level(imbalance, i6,  "6W",  "pearson"),
  null_level(imbalance, i12, "12W", "pearson"),
  null_level(imbalance, i6,  "6W",  "spearman"),
  null_level(imbalance, i12, "12W", "spearman"))

ambient_verdict <- sprintf(paste0(
  "AMBIENT-CORRELATION DIAGNOSTIC -- AND IT IS NOT WHAT PC1 WOULD HAVE TOLD YOU. PC1 is %.0f%% ",
  "of variance within 6W vs %.0f%% within 12W -- essentially UNCHANGED, so 'the 12W cohort has ",
  "less structure' is NOT the explanation. But the ambient COUPLING level collapses anyway: the ",
  "imbalance's couplings to ALL %d library sets have median |r| %.2f at 6W vs %.2f at 12W, and ",
  "%.0f%% vs %.0f%% of sets exceed |r|=0.5. So at 6W almost everything correlates with almost ",
  "everything and by 12W most of that is gone. ANY axis's coupling to ANY set collapses over ",
  "this window. That is the baseline the published r 0.75 -> 0.07 must be judged against, and ",
  "it is why PART B compares the SIZE of the collapse to the size of an arbitrary set's ",
  "collapse rather than asking whether a collapse happened."),
  pc1_by_timepoint$pc1_pct[1], pc1_by_timepoint$pc1_pct[2],
  ambient_by_timepoint$n_sets[1],
  ambient_by_timepoint$null_median_abs[1], ambient_by_timepoint$null_median_abs[2],
  100 * ambient_by_timepoint$frac_above_0.5[1], 100 * ambient_by_timepoint$frac_above_0.5[2])

# =============================================================================
# PART B: THE DELTA-RHO NULL (L2) -- the null the collapse claim never had
# =============================================================================
# Observed: d_obs = rho_6W(X, pro_comp) - rho_12W(X, pro_comp).
# Null:     d_s   = rho_6W(X, s)        - rho_12W(X, s)        for all 884 sets.
# Same method for observed and null, or the comparison is meaningless.
delta_rho_null <- function(x, label, method = "pearson") {
  r6  <- apply(scores[, i6,  drop = FALSE], 1,
               function(s) suppressWarnings(stats::cor(s, x[i6],  method = method)))
  r12 <- apply(scores[, i12, drop = FALSE], 1,
               function(s) suppressWarnings(stats::cor(s, x[i12], method = method)))
  d_s <- r6 - r12
  d_s <- d_s[is.finite(d_s)]
  o6  <- suppressWarnings(stats::cor(pro_comp[i6],  x[i6],  method = method))
  o12 <- suppressWarnings(stats::cor(pro_comp[i12], x[i12], method = method))
  d_obs <- o6 - o12
  pp <- mean(abs(d_s) >= abs(d_obs))
  tibble::tibble(
    axis = label, method = method, rho_6W = o6, rho_12W = o12, delta_obs = d_obs,
    null_median_delta = stats::median(d_s), null_q90_delta = stats::quantile(d_s, 0.90),
    n_sets = length(d_s),
    empirical_percentile = 100 * mean(d_s < d_obs),
    perm_p = pp,
    verdict = dplyr::case_when(
      pp > 0.10 ~ "NOT SPECIFIC (any gene set collapses this much)",
      pp > 0.05 ~ "marginal",
      TRUE      ~ "beats the null"))
}

collapse_null <- dplyr::bind_rows(
  delta_rho_null(imbalance,     "mitonuclear_imbalance", "pearson"),
  delta_rho_null(imbalance,     "mitonuclear_imbalance", "spearman"),
  delta_rho_null(ds$bio_comp,   "bio_comp",              "pearson"),
  delta_rho_null(ds$bio_comp,   "bio_comp",              "spearman"),
  delta_rho_null(ds$MASC_comp,  "MASC_comp",             "pearson"),
  delta_rho_null(ds$MASC_comp,  "MASC_comp",             "spearman"))

# Fisher z: is r6 vs r12 a DIFFERENCE at all, ignoring specificity? (n=12 each)
fisher_z <- function(r1, r2, n1 = 12, n2 = 12) {
  z1 <- atanh(r1); z2 <- atanh(r2)
  z  <- (z1 - z2) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  c(z = z, p = 2 * stats::pnorm(-abs(z)))
}
fisher_z_collapse <- collapse_null |>
  dplyr::filter(method == "pearson") |>
  dplyr::rowwise() |>
  dplyr::mutate(fisher_z = fisher_z(rho_6W, rho_12W)[["z"]],
                fisher_p = fisher_z(rho_6W, rho_12W)[["p"]]) |>
  dplyr::ungroup() |>
  dplyr::select(axis, rho_6W, rho_12W, fisher_z, fisher_p) |>
  dplyr::mutate(reads_as = "direct test of r6 vs r12, ignoring specificity; n=12/12")

# The full delta distributions, for the figure
collapse_distributions <- purrr::map_dfr(
  list(mitonuclear_imbalance = imbalance, bio_comp = ds$bio_comp, MASC_comp = ds$MASC_comp),
  function(x) {
    r6  <- apply(scores[, i6,  drop = FALSE], 1,
                 function(s) suppressWarnings(stats::cor(s, x[i6],  method = "pearson")))
    r12 <- apply(scores[, i12, drop = FALSE], 1,
                 function(s) suppressWarnings(stats::cor(s, x[i12], method = "pearson")))
    d <- r6 - r12
    tibble::tibble(delta = d[is.finite(d)])
  }, .id = "axis")

collapse_verdict <- {
  z  <- collapse_null |> dplyr::filter(axis == "mitonuclear_imbalance", method == "pearson")
  zs <- collapse_null |> dplyr::filter(axis == "mitonuclear_imbalance", method == "spearman")
  b  <- collapse_null |> dplyr::filter(axis == "bio_comp", method == "pearson")
  mm <- collapse_null |> dplyr::filter(axis == "MASC_comp", method == "pearson")
  fz <- fisher_z_collapse |> dplyr::filter(axis == "mitonuclear_imbalance")
  sprintf(paste0(
    "THE COLLAPSE NULL (the test the claim never had) -- AND THIS IS THE DEATH ARM'S BEST ",
    "SHOWING, SO READ IT CAREFULLY IN BOTH DIRECTIONS. The published collapse (imbalance~PRO ",
    "r=%.2f at 6W -> %.2f at 12W) gives delta=%.2f against a null median delta of %.2f -- i.e. ",
    "the ambient coupling ALREADY falls a long way, but the imbalance falls about twice as far, ",
    "landing at the %.0fth percentile of its collapses against all %d sets. perm p=%.3f (%s); ",
    "Spearman p=%.3f. As a direct difference, Fisher z=%.2f, p=%.3f at n=12/12. TWO ",
    "INDEPENDENT TESTS CONVERGE AT p~0.06. SO: this is NOT the flat null that the WITHIN-6W ",
    "coupling gave (perm p=0.49, script 33 PART C) -- the COLLAPSE is a better claim than the ",
    "6W coupling ever was, and it should not be dismissed as noise. But it does NOT reach ",
    "significance on either test, and it does not survive multiplicity (3 axes x 2 methods). ",
    "The comparators do not come close: bio_comp p=%.2f, MASC_comp p=%.2f. VERDICT: MARGINAL ",
    "-- suggestive, unproven, and NOT reportable as a finding at n=12/12. CRITICAL CAVEAT: the ",
    "only axis that gets even this far is the CIRCULAR one (pro_comp IS a MitoCarta set, like ",
    "the imbalance). PART F re-runs this on non-mito death axes -- read that before believing ",
    "this number."),
    z$rho_6W, z$rho_12W, z$delta_obs, z$null_median_delta, z$empirical_percentile,
    z$n_sets, z$perm_p, z$verdict, zs$perm_p, fz$fisher_z, fz$fisher_p,
    b$perm_p, mm$perm_p)
}

# =============================================================================
# PART C: DOES PRIMING COLLAPSE AT ALL? (L3 + surfacing the buried contrast)
# =============================================================================
# C1. The contrast that was fitted, saved, and never read.
priming_own_stats <- tibble::as_tibble(dt$h1$state_tests) |>
  dplyr::mutate(int_significant = int_p < 0.05,
                geno_significant = geno_p < 0.05) |>
  dplyr::mutate(source = "death_timing_substrate.rds$h1$state_tests (script 23, already on disk)")

priming_group_means <- dt$h1$group_means

# C2. The interaction-percentile null: does the death programme attenuate MORE
# than an arbitrary programme? coef_table already holds beta_int for all 884 sets.
ct <- tibble::as_tibble(ov$coef_table)
death_cats <- c("Apoptosis", "Biogenesis_apoptosis_intersections")
interaction_percentile <- ct |>
  dplyr::mutate(
    is_death = category_primary %in% death_cats,
    abs_beta_int = abs(beta_int)) |>
  dplyr::group_by(category_primary) |>
  dplyr::summarise(
    n_sets = dplyr::n(),
    median_abs_beta_int = stats::median(abs_beta_int, na.rm = TRUE),
    median_beta_int = stats::median(beta_int, na.rm = TRUE),
    n_camera_fdr05 = sum(camera_fdr < 0.05, na.rm = TRUE),
    n_roast_fdr05  = sum(roast_fdr  < 0.05, na.rm = TRUE),
    .groups = "drop") |>
  dplyr::mutate(
    pct_of_all_sets_below = 100 * vapply(median_abs_beta_int, function(v)
      mean(abs(ct$beta_int) < v, na.rm = TRUE), numeric(1))) |>
  dplyr::arrange(dplyr::desc(median_abs_beta_int))

# per-set percentile for the named apoptosis sets that carry the narrative
named_death_sets <- c("MITOCARTA_APOPTOSIS_PRO", "MITOCARTA_APOPTOSIS_ANTI",
                      "CDC_PRODEATH_APOPTOSIS", "CDC_PRODEATH_APOPTOSIS_MITO",
                      "CDC_PRODEATH_CICD", "APOP_BH3_REACTOME",
                      "APOP_INTRINSIC_REACTOME", "APOP_HALLMARK")
death_set_percentile <- ct |>
  dplyr::filter(set_name %in% named_death_sets) |>
  dplyr::mutate(
    abs_beta_int_percentile = 100 * vapply(abs(beta_int), function(v)
      mean(abs(ct$beta_int) < v, na.rm = TRUE), numeric(1)),
    camera_cor_measured = camera_cor,
    n_eff = camera_ngenes / (1 + (camera_ngenes - 1) * pmax(camera_cor, 0))) |>
  dplyr::select(set_name, camera_ngenes, beta_int, d6, d12, int_p,
                abs_beta_int_percentile, camera_cor_measured, n_eff,
                camera_p, camera_fdr, roast_p, roast_fdr) |>
  dplyr::arrange(dplyr::desc(abs_beta_int_percentile))

# C3. The five-method concordance table, assembled from disk (no new fits).
# NOTE: the 5 contrasts live under $fgsea_results / $fgsea_combined (75 rows =
# 15 modalities x 5 contrasts), NOT at the top level of the rds.
cd_fg  <- readRDS(here::here("results", "cell_death_fgsea.rds"))
fg_all <- tibble::as_tibble(cd_fg$fgsea_combined)
fg_int <- fg_all |> dplyr::filter(contrast == "interaction")
apop_int <- fg_int |> dplyr::filter(pathway == "CD_Apoptosis")

method_concordance <- tibble::tribble(
  ~method, ~unit, ~what_it_tests, ~result, ~verdict,
  "per-sample composite interaction (23:162-167)", "n=24 lm",
  "priming ~ timepoint*genotype",
  sprintf("PRO int_p=%.3f; priming int_p=%.3f; ANTI %.3f; p53 %.3f",
          priming_own_stats$int_p[priming_own_stats$metric == "PRO"],
          priming_own_stats$int_p[priming_own_stats$metric == "priming (PRO-ANTI)"],
          priming_own_stats$int_p[priming_own_stats$metric == "ANTI"],
          priming_own_stats$int_p[priming_own_stats$metric == "p53"]),
  "NULL",
  "fGSEA on the DESeq2 interaction (12:158)", "gene ranks",
  "are death modalities enriched in the interaction?",
  sprintf("Apoptosis NES=%.2f padj=%.3f; min padj over 15 modalities=%.3f",
          apop_int$NES[1], apop_int$padj[1], min(fg_int$padj, na.rm = TRUE)),
  "NULL",
  "branch-1 binomial (16:133)", "gene signs",
  "pro-death dLFC sign asymmetry vs 50%",
  "52.7% supporting, p=0.273 (see PART G -- 50% is the WRONG null)",
  "NULL (and mis-nulled)",
  "CAMERA/ROAST on apoptosis sets (17:367-383)", "correlation-aware",
  "interaction, accounting for inter-gene correlation",
  sprintf("%d/%d apoptosis-category sets at camera_fdr<0.05",
          sum(ct$camera_fdr[ct$category_primary == "Apoptosis"] < 0.05, na.rm = TRUE),
          sum(ct$category_primary == "Apoptosis")),
  "NULL",
  "Gate 2 (14:91)", "one-sample t, 23 genes",
  "is the PRO module's mean interaction LFC != 0?",
  "p=0.023 -- the SOLE positive; see PART E",
  "POSITIVE (but self-contained + assumes gene independence)")

# C4. The structural limit -- stated, not tested.
cross_sectional_limit <- tibble::tibble(
  question = "Does the priming collapse correlate with mitochondrial reprioritisation?",
  estimable = FALSE,
  why = paste(
    "The design is CROSS-SECTIONAL: 6W and 12W are DIFFERENT ANIMALS. There is no",
    "within-animal delta, so a per-sample correlation between two CHANGES does not",
    "exist and cannot be computed. Any such 'correlation' can only be two 4-point",
    "group trajectories co-moving, which carries no inference (n=4, and the two",
    "trajectories are not independent draws)."),
  what_is_estimable = paste(
    "Whether the two share a genotype x time INTERACTION. Both are null:",
    sprintf("priming int_p=%.3f (23$h1$state_tests); imbalance int_p=%.3f (23$h2$imbalance_test).",
            priming_own_stats$int_p[priming_own_stats$metric == "priming (PRO-ANTI)"],
            dt$h2$imbalance_test$int_p[1]),
    "Two null interactions cannot co-vary informatively."))

priming_verdict <- sprintf(paste0(
  "DOES PRIMING COLLAPSE? THE CONTRAST WAS FITTED, SAVED, AND NEVER READ. ",
  "`death_timing_substrate.rds$h1$state_tests` has carried the answer since script 23 was ",
  "written: PRO int_p=%.3f, priming(PRO-ANTI) int_p=%.3f, ANTI %.3f, p53 %.3f -- NOTHING ",
  "significant on the interaction. The narrative instead quoted the GROUP MEANS ('~3x more at ",
  "6W' = +0.48 vs +0.15). That is the retracted-error pattern a THIRD time (cf. the mtDNA ",
  "'difference' and the mitonuclear imbalance) -- and it is the worst instance, because here ",
  "the test existed on disk. WHAT IS POWERED AND SURVIVES: the GENOTYPE main effect -- PRO ",
  "geno_beta=%+.3f (p=%.3f), priming %+.3f (p=%.3f). Myc RAISES apoptotic priming. That is a ",
  "genotype contrast, so the genotype-independent IEG axis (p=0.49) cannot bias it. CEILING: ",
  "at n=6/group, the interaction being null is UNSUPPORTED, not REFUTED."),
  priming_own_stats$int_p[priming_own_stats$metric == "PRO"],
  priming_own_stats$int_p[priming_own_stats$metric == "priming (PRO-ANTI)"],
  priming_own_stats$int_p[priming_own_stats$metric == "ANTI"],
  priming_own_stats$int_p[priming_own_stats$metric == "p53"],
  priming_own_stats$geno_beta[priming_own_stats$metric == "PRO"],
  priming_own_stats$geno_p[priming_own_stats$metric == "PRO"],
  priming_own_stats$geno_beta[priming_own_stats$metric == "priming (PRO-ANTI)"],
  priming_own_stats$geno_p[priming_own_stats$metric == "priming (PRO-ANTI)"])

# =============================================================================
# PART D: EXCESS OVER THE GLOBAL ATTENUATION -- THE MECHANISM DISCRIMINATOR
# =============================================================================
# THE QUESTION. Myc drives apoptotic priming (PART C, powered). The Myc programme
# attenuates 6W->12W (~66% fade + ~34% WT-convergence, Issue #6). If priming is
# simply PART OF that programme, it attenuates WITH it -- no death-specific gate,
# no mitochondrial-state mechanism, and the canonical Evan/Lowe account (Myc
# sensitises to apoptosis via BH3-only proteins) covers the whole phenomenon.
# That hypothesis makes a sharp prediction: THE DEATH ATTENUATION IS NOT IN EXCESS
# OF THE GLOBAL ATTENUATION. This is that test.
#
# Metric, per gene, aligned to Myc's induction direction (the Issue #6 identity):
#   atten_g = sign(m6) * (m6 - m12)
# Positive = the Myc effect shrank. NOTE: alignment makes this biased positive by
# regression-to-the-mean alone -- which is exactly why the null genes get the SAME
# alignment. Reading `obs` without `null_mean` is meaningless.
dds_int <- readRDS(here::here("results", "dds_int_run.rds"))
disp <- DESeq2::dispersions(dds_int); names(disp) <- rownames(dds_int)

gene_tbl <- tibble::tibble(
  ensembl  = rownames(ir$myc_6W_raw),
  baseMean = ir$myc_6W_raw$baseMean,
  lfc6     = m6[rownames(ir$myc_6W_raw)],
  lfc12    = m12[rownames(ir$myc_6W_raw)],
  lfc_int  = lint[rownames(ir$myc_6W_raw)],
  dispersion = disp[rownames(ir$myc_6W_raw)]) |>
  dplyr::filter(baseMean > 0, !is.na(dispersion), dispersion > 0,
                !is.na(lfc6), !is.na(lfc12), lfc6 != 0) |>
  dplyr::mutate(
    atten_g  = sign(lfc6) * (lfc6 - lfc12),      # aligned attenuation
    bm_bin   = dplyr::ntile(log10(baseMean), 10),
    disp_bin = dplyr::ntile(log10(dispersion), 10),
    bin      = (bm_bin - 1L) * 10L + disp_bin)

B <- 5000L

# `background` = the pool the null is drawn from. NULL = the whole transcriptome
# (asks "is this set special among all genes"). Supply a gene set to ask a
# CONDITIONAL question -- e.g. "is APOPTOSIS_PRO special among MitoCarta genes",
# which is the only way to separate a death effect from the mito-content effect.
# Query genes are excluded from their own pool.
#
# *** READ THIS BEFORE QUOTING ANY p FROM THIS FUNCTION. *** It draws null genes
# INDEPENDENTLY from their bins, so the null distribution of the mean is too NARROW
# for real (co-regulated) gene sets -> p_emp is ANTI-CONSERVATIVE, exactly as
# script 21 documented for AP6.2 ("raw z inflated by inter-gene correlation -> use
# EFFECT-SIZE ranking, not absolute z"). It matches expression x dispersion; it does
# NOT match inter-gene correlation. Consequences here:
#   - A NON-significant result is SAFE (a proper correction only widens the null).
#     PART D's "nothing is in excess" and PART E's Gate 2 verdict are therefore firm.
#   - A SIGNIFICANT result is NOT safe on its p alone. Quote the EFFECT SIZE and the
#     DIRECTION, and require independent corroboration. The two significant results
#     below (PART F3) both have it: APOPTOSIS_PRO rising HALF as much as the mito
#     background reproduces script 23's mitoPPS de-prioritisation (-0.21, padj 0.022)
#     by an independent method.
# This is the same trap as Gate 2 (PART E), so do not fall into it while reporting
# Gate 2's failure.
matched_null <- function(ens, label, metric_col, background = NULL,
                         excess_label = "the global attenuation") {
  qi <- which(gene_tbl$ensembl %in% ens & !is.na(gene_tbl[[metric_col]]))
  if (length(qi) < 5) return(NULL)
  pool <- gene_tbl |>
    dplyr::filter(!is.na(.data[[metric_col]]), !ensembl %in% ens)
  if (!is.null(background)) pool <- dplyr::filter(pool, ensembl %in% background)
  if (nrow(pool) < 50) return(NULL)
  bin_vals <- split(pool[[metric_col]], pool$bin)
  obs    <- mean(gene_tbl[[metric_col]][qi])
  bins_q <- as.character(gene_tbl$bin[qi])
  samp   <- vapply(bins_q, function(b) {
    v <- bin_vals[[b]]
    if (is.null(v) || !length(v)) rep(NA_real_, B) else sample(v, B, replace = TRUE)
  }, numeric(B))
  null_means <- rowMeans(samp, na.rm = TRUE)
  nm <- mean(null_means); ns <- stats::sd(null_means)
  pe <- (1 + sum(abs(null_means - nm) >= abs(obs - nm))) / (B + 1)
  zz <- (obs - nm) / ns
  tibble::tibble(
    set = label, n = length(qi), n_pool = nrow(pool), observed = obs,
    null_mean = nm, null_sd = ns, z = zz,
    p_emp_two_sided = pe, excess = obs - nm,
    # DIRECTION-AWARE. The p is two-sided, so a bare "in excess" would be printed
    # for a NEGATIVE z -- i.e. for a set sitting BELOW its background, which is the
    # opposite claim. Say which side.
    verdict = dplyr::case_when(
      pe >= 0.05 ~ paste0("not distinguishable from ", excess_label),
      zz > 0     ~ "ABOVE the matched background",
      TRUE       ~ "BELOW the matched background"))
}

# the death rosters -- including the non-mito split that breaks the circularity
pro_mito_syms <- ch$interaction_by_geneset$apoptosis_pro$gene_symbol
prodeath_all  <- cdg$mouse_symbol[cdg$effect == "pro-death"]
prodeath_nomito <- cdg$mouse_symbol[cdg$effect == "pro-death" & !cdg$is_mitochondrial]
prodeath_mito   <- cdg$mouse_symbol[cdg$effect == "pro-death" &  cdg$is_mitochondrial]
prosurv_all   <- cdg$mouse_symbol[cdg$effect == "pro-survival"]

gmt_path <- here::here("data", "genesets_from_library", "mammary_mito_myc_metab_v1_mouse.gmt")
gmt <- fgsea::gmtPathways(gmt_path)
set_ens <- function(nm) if (!is.null(gmt[[nm]])) ens_of(gmt[[nm]]) else character(0)

mitocarta_bg <- ens_of(gmt[["MITOCARTA_NUCLEAR_ENCODED"]])

# labels carry the MAPPED gene count, not the roster count -- ens_of() drops
# unmapped symbols and a stale label is how a claim starts drifting.
lab <- function(stem, ens) sprintf("%s (n=%d)", stem, length(ens))
death_rosters <- list()
add_roster <- function(stem, ens) {
  if (length(ens) >= 5) death_rosters[[lab(stem, ens)]] <<- ens
}
add_roster("MITOCARTA_APOPTOSIS_PRO -- pro_comp, MITO",  ens_of(pro_mito_syms))
add_roster("pro-death ALL (consolidated)",               ens_of(prodeath_all))
add_roster("pro-death NON-MITO -- non-circular",         ens_of(prodeath_nomito))
add_roster("pro-death MITO-ONLY",                        ens_of(prodeath_mito))
add_roster("pro-survival ALL",                           ens_of(prosurv_all))
add_roster("APOP_BH3_REACTOME",                          set_ens("APOP_BH3_REACTOME"))
add_roster("APOP_INTRINSIC_REACTOME",                    set_ens("APOP_INTRINSIC_REACTOME"))
# comparators: programmes whose attenuation IS established (Issue #4/#6)
add_roster("MITOCARTA_OXPHOS_SUBUNITS (comparator)",     set_ens("MITOCARTA_OXPHOS_SUBUNITS"))
add_roster("MYC_HALLMARK_MYC_TARGETS_V2 (comparator)",   set_ens("MYC_HALLMARK_MYC_TARGETS_V2"))

# multiplicity matters: this is ~9 sets tested against the same null
atten_excess <- dplyr::bind_rows(lapply(names(death_rosters), function(n)
  matched_null(death_rosters[[n]], n, "atten_g"))) |>
  dplyr::mutate(p_bh = stats::p.adjust(p_emp_two_sided, method = "BH"),
                verdict = dplyr::case_when(
                  p_bh >= 0.05 ~ "not distinguishable from the global attenuation",
                  z > 0        ~ "attenuates MORE than matched background (BH)",
                  TRUE         ~ "attenuates LESS than matched background (BH)"))

# D2. The Issue #6 identity applied to the death sets: shared-developmental vs
# Myc-specific. A decline SHARED by both genotypes CANCELS in the genotype gap and
# therefore cannot create an attenuation.
decomp_one <- function(genes, label, universe_genes) {
  e <- intersect(genes, universe_genes)
  e <- e[!is.na(m6[e]) & !is.na(m12[e]) & !is.na(tn[e]) & !is.na(tp[e]) & m6[e] != 0]
  if (length(e) < 5) return(NULL)
  d <- sign(m6[e])
  gap6 <- d * m6[e]; gap12 <- d * m12[e]
  wc <- d * tn[e];   mf <- d * tp[e]
  atten <- mean(gap6) - mean(gap12)
  tibble::tibble(program = label, n = length(e),
                 gap6 = mean(gap6), gap12 = mean(gap12), atten = atten,
                 wt_conv = mean(wc), myc_fade = mean(mf),
                 conv_pct = 100 * mean(wc) / atten, fade_pct = 100 * (-mean(mf)) / atten,
                 frac_wt_toward = mean(wc > 0), frac_myc_retreat = mean(mf < 0))
}
p6adj <- stats::setNames(ir$myc_6W_raw$padj, rownames(ir$myc_6W_raw))
universes <- list(
  divergent = names(m6)[!is.na(p6adj) & p6adj < 0.1],
  effect    = names(m6)[!is.na(m6) & abs(m6) > 0.5],
  expressed = gene_tbl$ensembl)

death_decomp <- dplyr::bind_rows(lapply(names(universes), function(u) {
  rows <- c(list(decomp_one(universes[[u]], "ALL (the global attenuation)", universes[[u]])),
            lapply(names(death_rosters), function(n)
              decomp_one(death_rosters[[n]], n, universes[[u]])))
  dplyr::bind_rows(rows) |> dplyr::mutate(universe = u, .before = 1)
}))

# D3. The Tang NES: is the apoptosis decline Myc-SPECIFIC or shared with WT?
tang_shared_vs_specific <- fg_all |>
  dplyr::filter(pathway == "CD_Apoptosis") |>
  dplyr::select(contrast, pathway, NES, pval, padj) |>
  dplyr::mutate(contrast = factor(contrast, levels = c(
    "myc_effect_6W", "myc_effect_12W", "temporal_neg", "temporal_pos", "interaction"))) |>
  dplyr::arrange(contrast)

tang_verdict <- {
  tneg <- tang_shared_vs_specific$NES[tang_shared_vs_specific$contrast == "temporal_neg"]
  tpos <- tang_shared_vs_specific$NES[tang_shared_vs_specific$contrast == "temporal_pos"]
  sprintf(paste0(
    "THE TANG NES -1.38 IS A SHARED DEVELOPMENTAL DECLINE, NOT A MYC EFFECT. The walkthrough ",
    "cites apoptosis temporal_pos NES=%.2f as death-spine evidence ('the apoptotic programme ",
    "DECLINES in Myc+ with age'). But temporal_neg (WT) = %.2f -- the WT gland declines %.0f%% ",
    "as much. By Issue #6's OWN exact identity, a decline SHARED by both genotypes CANCELS in ",
    "the genotype gap and therefore cannot create an attenuation. The Myc-specific excess is ",
    "%.2f NES units. Do not cite temporal_pos alone."),
    tpos, tneg, 100 * abs(tneg / tpos), abs(tpos) - abs(tneg))
}

atten_verdict <- {
  a <- atten_excess |> dplyr::filter(grepl("MITOCARTA_APOPTOSIS_PRO", set))
  n <- atten_excess |> dplyr::filter(grepl("NON-MITO", set))
  o <- atten_excess |> dplyr::filter(grepl("OXPHOS_SUBUNITS", set))
  y <- atten_excess |> dplyr::filter(grepl("MYC_HALLMARK", set))
  sprintf(paste0(
    "THE MECHANISM DISCRIMINATOR: IS THE DEATH ATTENUATION IN EXCESS OF THE MYC PROGRAMME ",
    "FADING? Aligned per-gene attenuation vs an expression x dispersion-matched null (script ",
    "21's idiom; the null genes get the SAME sign-alignment, so the regression-to-the-mean bias ",
    "cancels). pro_comp's %d mito genes: observed %+.3f vs matched null %+.3f, z=%.2f, ",
    "p_emp=%.3f => %s. The non-circular %d non-mito pro-death genes: observed %+.3f vs null ",
    "%+.3f, z=%.2f, p_emp=%.3f => %s. THE CONTROL THAT VALIDATES THE TEST: OXPHOS subunits -- ",
    "whose attenuation Issue #4 ESTABLISHED -- are ALSO not in excess (z=%.2f, p_emp=%.3f), and ",
    "neither is the MYC-target core (z=%.2f, p=%.3f). That is not the test failing; it is the ",
    "point. The attenuation is GLOBAL (Issue #4: median |LFC| 0.268 -> 0.204 across the ",
    "transcriptome), so NO programme attenuates more than expression-matched genes -- death ",
    "included. CONCLUSION: priming fades because the MYC PROGRAMME fades. That is canonical ",
    "Evan/Lowe Myc-induced apoptosis via BH3-only proteins, and it needs no death-specific gate ",
    "and no mitochondrial-state gate. Neither is supported here."),
    a$n[1], a$observed[1], a$null_mean[1], a$z[1], a$p_emp_two_sided[1], a$verdict[1],
    n$n[1], n$observed[1], n$null_mean[1], n$z[1], n$p_emp_two_sided[1], n$verdict[1],
    o$z[1], o$p_emp_two_sided[1], y$z[1], y$p_emp_two_sided[1])
}

# =============================================================================
# PART E: A COMPETITIVE NULL FOR GATE 2
# =============================================================================
# Gate 2 (14:91) runs t.test(PRO_interaction_LFC, mu = 0) -> p=0.023 and calls
# DECISION_POINT. Two things are wrong with mu=0 as the null:
#   (i)  Under a GLOBAL attenuation, every Myc-induced module has a negative mean
#        interaction LFC. mu=0 is not the no-death-specific-effect hypothesis.
#   (ii) It treats 23 genes as independent draws. Script 17 MEASURED the inter-gene
#        correlation of these sets at 0.25-0.36 -> the effective n is far below 23.
gate2_genes <- ch$interaction_by_geneset$apoptosis_pro
gate2_anti  <- ch$interaction_by_geneset$apoptosis_anti

gate2_selfcontained <- tibble::tibble(
  module = c("Apoptosis-PRO", "Apoptosis-ANTI"),
  n = c(nrow(gate2_genes), nrow(gate2_anti)),
  mean_int_lfc = c(mean(gate2_genes$log2FoldChange), mean(gate2_anti$log2FoldChange)),
  t_p_mu0 = c(stats::t.test(gate2_genes$log2FoldChange, mu = 0)$p.value,
              stats::t.test(gate2_anti$log2FoldChange,  mu = 0)$p.value),
  # script 14 messages the nominal-p count; padj is the honest one. Report both,
  # so this table cannot be read as disagreeing with 14's own console output.
  n_genes_p05 = c(sum(gate2_genes$pvalue < 0.05, na.rm = TRUE),
                  sum(gate2_anti$pvalue  < 0.05, na.rm = TRUE)),
  n_genes_padj05 = c(sum(gate2_genes$padj < 0.05, na.rm = TRUE),
                     sum(gate2_anti$padj  < 0.05, na.rm = TRUE)),
  reads_as = "the ORIGINAL Gate 2 test, reproduced -- self-contained, assumes gene independence")

gate2_competitive <- dplyr::bind_rows(
  matched_null(ens_of(gate2_genes$gene_symbol), "Apoptosis-PRO (Gate 2's 23 genes)",
               "lfc_int", excess_label = "the global attenuation"),
  matched_null(ens_of(gate2_anti$gene_symbol),  "Apoptosis-ANTI (Gate 2's 7 genes)",
               "lfc_int", excess_label = "the global attenuation"),
  matched_null(ens_of(prodeath_nomito), lab("pro-death NON-MITO -- non-circular",
                                            ens_of(prodeath_nomito)),
               "lfc_int", excess_label = "the global attenuation"),
  # the conditional version: is PRO special AMONG MITOCARTA GENES?
  matched_null(ens_of(gate2_genes$gene_symbol),
               "Apoptosis-PRO vs MITOCARTA background", "lfc_int",
               background = mitocarta_bg, excess_label = "the MitoCarta background"))

# the correlation correction Gate 2 never applied
gate2_effective_n <- ct |>
  dplyr::filter(set_name %in% c("MITOCARTA_APOPTOSIS_PRO", "MITOCARTA_APOPTOSIS_ANTI")) |>
  dplyr::transmute(
    set_name, n_genes = camera_ngenes, measured_inter_gene_cor = camera_cor,
    effective_n = camera_ngenes / (1 + (camera_ngenes - 1) * pmax(camera_cor, 0)),
    camera_p, camera_fdr, roast_p, roast_fdr,
    reads_as = paste("Gate 2 used df = n-1 = 22. With the MEASURED inter-gene correlation,",
                     "the effective n is what this column says. CAMERA is the",
                     "correlation-aware version of Gate 2's own question."))

gate2_verdict <- {
  g <- gate2_competitive |> dplyr::filter(grepl("PRO", set))
  e <- gate2_effective_n |> dplyr::filter(set_name == "MITOCARTA_APOPTOSIS_PRO")
  sprintf(paste0(
    "GATE 2's p=0.023 IS THE WRONG TEST. Reproduced: PRO mean interaction LFC = %+.4f, ",
    "t-vs-mu=0 p=%.3f, with only %d of %d genes individually padj<0.05. Two failures. (i) mu=0 ",
    "is not the null: under the GLOBAL attenuation every Myc-induced module has a negative mean ",
    "interaction LFC. Against an expression x dispersion-MATCHED null: observed %+.4f vs null ",
    "%+.4f, z=%.2f, p_emp=%.3f => %s. (ii) It treats %d genes as independent. The MEASURED ",
    "inter-gene correlation is %.2f, giving an effective n of %.1f -- not %d. Gate 2's own ",
    "question, asked correlation-aware, is CAMERA: p=%.3f, FDR=%.3f. Script 17 already ran it ",
    "and found nothing at FDR across the whole apoptosis category. The one positive in the ",
    "death arm does not survive being asked properly."),
    gate2_selfcontained$mean_int_lfc[1], gate2_selfcontained$t_p_mu0[1],
    gate2_selfcontained$n_genes_padj05[1], gate2_selfcontained$n[1],
    g$observed[1], g$null_mean[1], g$z[1], g$p_emp_two_sided[1], g$verdict[1],
    e$n_genes[1], e$measured_inter_gene_cor[1], e$effective_n[1], e$n_genes[1],
    e$camera_p[1], e$camera_fdr[1])
}

# =============================================================================
# PART F: THE NON-CIRCULAR REBUILD -- break the mito <-> priming circularity
# =============================================================================
# pro_comp IS a mitochondrial gene set. "The mito state couples to death priming"
# correlates mito genes with mito genes. Rebuild priming from NON-mito death genes
# and redo everything. If Myc's priming effect survives, it is not an artifact of
# the mito definition; if the COUPLING survives, the bridge deserves another look.
comp_expr <- function(ens) {                     # script 23's comp_expr idiom
  e <- intersect(ens, rownames(vst_mat))
  if (length(e) < 5) return(rep(NA_real_, length(samples)))
  colMeans(t(scale(t(vst_mat[e, , drop = FALSE]))))
}
pro_nomito  <- comp_expr(ens_of(prodeath_nomito))
pro_mito41  <- comp_expr(ens_of(prodeath_mito))
surv_nomito <- comp_expr(ens_of(cdg$mouse_symbol[cdg$effect == "pro-survival" &
                                                 !cdg$is_mitochondrial]))
priming_nomito <- pro_nomito - surv_nomito

noncircular_axes <- list(
  `pro_comp (MITO, 25 genes -- the published axis)` = pro_comp,
  `pro-death NON-MITO z-composite`                 = pro_nomito,
  `pro-death MITO-ONLY z-composite`                = pro_mito41,
  `priming NON-MITO (pro - surv, non-mito)`        = priming_nomito,
  `APOP_BH3_REACTOME (GSVA)`                       = if ("APOP_BH3_REACTOME" %in% rownames(scores))
    scores["APOP_BH3_REACTOME", ] else NULL,
  `CDC_PRODEATH_APOPTOSIS (GSVA)`                  = if ("CDC_PRODEATH_APOPTOSIS" %in% rownames(scores))
    scores["CDC_PRODEATH_APOPTOSIS", ] else NULL)
noncircular_axes <- noncircular_axes[!vapply(noncircular_axes, is.null, logical(1))]

# F1. Does Myc raise the axis, and does it collapse? (genotype + interaction)
fit_axis <- function(y, label) {
  if (all(is.na(y))) return(NULL)
  f <- summary(stats::lm(y ~ ds$timepoint * ds$myc_status))$coefficients
  gr <- tapply(y, ds$group, mean)[group_levels]
  tibble::tibble(
    axis = label,
    geno_beta = f["ds$myc_statuspos", "Estimate"], geno_p = f["ds$myc_statuspos", "Pr(>|t|)"],
    int_beta  = f["ds$timepoint12W:ds$myc_statuspos", "Estimate"],
    int_p     = f["ds$timepoint12W:ds$myc_statuspos", "Pr(>|t|)"],
    gap_6W = gr[["6W_pos"]] - gr[["6W_neg"]], gap_12W = gr[["12W_pos"]] - gr[["12W_neg"]])
}
noncircular_stats <- dplyr::bind_rows(
  lapply(names(noncircular_axes), function(n) fit_axis(noncircular_axes[[n]], n))) |>
  dplyr::mutate(geno_p_bh = stats::p.adjust(geno_p, method = "BH"),
                int_p_bh  = stats::p.adjust(int_p,  method = "BH"),
                geno_sig = geno_p_bh < 0.05, int_sig = int_p_bh < 0.05,
                geno_direction = dplyr::if_else(geno_beta > 0, "Myc RAISES", "Myc LOWERS"))

# F2. Does the IMBALANCE couple to a NON-mito death axis -- and does that survive
# the delta-rho null? This is the bridge's last chance.
noncircular_coupling <- purrr::map_dfr(names(noncircular_axes), function(n) {
  y <- noncircular_axes[[n]]
  if (all(is.na(y))) return(NULL)
  r6  <- suppressWarnings(stats::cor(imbalance[i6],  y[i6],  method = "pearson"))
  r12 <- suppressWarnings(stats::cor(imbalance[i12], y[i12], method = "pearson"))
  fz  <- fisher_z(r6, r12)
  tibble::tibble(death_axis = n, rho_6W = r6, rho_12W = r12, delta = r6 - r12,
                 fisher_z = fz[["z"]], fisher_p = fz[["p"]])
})

# and the delta-rho null for the non-circular axis specifically
delta_rho_null_y <- function(x, y, label, method = "pearson") {
  r6  <- apply(scores[, i6,  drop = FALSE], 1,
               function(s) suppressWarnings(stats::cor(s, x[i6],  method = method)))
  r12 <- apply(scores[, i12, drop = FALSE], 1,
               function(s) suppressWarnings(stats::cor(s, x[i12], method = method)))
  d_s <- (r6 - r12)[is.finite(r6 - r12)]
  d_obs <- suppressWarnings(stats::cor(y[i6], x[i6], method = method)) -
           suppressWarnings(stats::cor(y[i12], x[i12], method = method))
  pp <- mean(abs(d_s) >= abs(d_obs))
  tibble::tibble(death_axis = label, delta_obs = d_obs,
                 null_median_delta = stats::median(d_s),
                 perm_p = pp,
                 # same 3-way ladder as delta_rho_null(); a single >0.10 cut would
                 # print "beats the null" for a p of 0.07, which it does not.
                 verdict = dplyr::case_when(
                   pp > 0.10 ~ "NOT SPECIFIC (any gene set collapses this much)",
                   pp > 0.05 ~ "marginal",
                   TRUE      ~ "beats the null"))
}
noncircular_collapse_null <- dplyr::bind_rows(lapply(
  names(noncircular_axes), function(n) {
    y <- noncircular_axes[[n]]
    if (all(is.na(y))) return(NULL)
    delta_rho_null_y(imbalance, y, n)
  }))

# F3. THE TEST THE NON-CIRCULAR RESULT FORCES. If Myc's "priming" effect is
# carried by MITO membership rather than DEATH membership, then it is script 32's
# content effect wearing a death label: Myc raises nuclear MitoCarta ~19%
# wholesale, so ANY MitoCarta subset rises -- including one labelled APOPTOSIS_PRO.
# This is the `mitocarta-sets-are-membership-loose` trap that already produced the
# retracted "commissioned but unbuilt" claim (Mrpl12 carried MITOCARTA_TRANSCRIPTION).
# The question is CONDITIONAL: is APOPTOSIS_PRO's Myc effect bigger than that of
# expression-matched OTHER MitoCarta genes?
priming_vs_mito_background <- dplyr::bind_rows(
  matched_null(ens_of(pro_mito_syms), "APOPTOSIS_PRO vs ALL genes (unconditional)",
               "lfc6", excess_label = "the transcriptome"),
  matched_null(ens_of(pro_mito_syms), "APOPTOSIS_PRO vs MITOCARTA background",
               "lfc6", background = mitocarta_bg,
               excess_label = "the mito-content effect (script 32: Myc raises MitoCarta ~19%)"),
  matched_null(ens_of(gate2_anti$gene_symbol), "APOPTOSIS_ANTI vs MITOCARTA background",
               "lfc6", background = mitocarta_bg,
               excess_label = "the mito-content effect"),
  matched_null(ens_of(prodeath_nomito), lab("pro-death NON-MITO vs ALL genes",
                                            ens_of(prodeath_nomito)),
               "lfc6", excess_label = "the transcriptome"))

mito_background_verdict <- {
  u <- priming_vs_mito_background |> dplyr::filter(grepl("vs ALL genes \\(unconditional\\)", set))
  c2 <- priming_vs_mito_background |> dplyr::filter(grepl("PRO vs MITOCARTA", set))
  sprintf(paste0(
    "IS 'MYC RAISES APOPTOTIC PRIMING' A DEATH EFFECT OR THE MITO-CONTENT EFFECT? IT IS THE ",
    "CONTENT EFFECT -- AND LESS OF IT THAN AVERAGE. pro_comp is a MitoCarta set, and script 32 ",
    "established that Myc raises nuclear MitoCarta ~19%% WHOLESALE, so ANY MitoCarta subset ",
    "rises whatever its label. UNCONDITIONAL (vs all genes): APOPTOSIS_PRO Myc LFC6 %+.3f vs ",
    "null %+.3f, z=%.2f, p=%.3f => %s. CONDITIONAL (vs expression-matched OTHER MitoCarta ",
    "genes -- the test that separates a death effect from the content effect): observed %+.3f ",
    "vs a MitoCarta background of %+.3f, z=%.2f, p=%.3f => %s. Read the SIGN: the pro-apoptotic ",
    "mito genes rise by HALF what a typical mito gene rises by. So 'Myc raises apoptotic ",
    "priming' (geno p=0.037/0.038) is not a death-specific induction at all -- it is ",
    "BELOW-AVERAGE participation in the wholesale mito-content rise. The composite goes up ",
    "because its genes are MITOCHONDRIAL, not because they are PRO-APOPTOTIC. IMPORTANT -- ",
    "THIS IS NOT NEW, IT IS A CONVERGENCE: script 23 already found exactly this by mitoPPS ",
    "('Myc raises PRO absolutely (+0.20) but RELATIVELY de-prioritises Apoptosis-PRO within the ",
    "mito compartment, mitoPPS diff -0.21, padj 0.022 -- biogenesis crowds it out', ",
    "walkthrough:762-764). Two independent methods agree. The walkthrough reported it as a ",
    "'two-lens nuance' alongside the priming claim; it is better read as the EXPLANATION of the ",
    "priming claim. Same membership-looseness that produced the retracted 'commissioned but ",
    "unbuilt' claim -- see [[mitocarta-sets-are-membership-loose]]. STATISTICAL CEILING ON ",
    "THIS ONE, STATED BECAUSE IT IS THE TRAP THIS SCRIPT IS ABOUT: the matched null samples ",
    "genes INDEPENDENTLY, so p=0.012 is ANTI-CONSERVATIVE (script 21's own AP6.2 caveat -- ",
    "'raw z inflated by inter-gene correlation'). DO NOT QUOTE THE p. What carries this claim ",
    "is (a) the EFFECT SIZE -- half the background rise -- and (b) the INDEPENDENT mitoPPS ",
    "convergence above. Non-significant matched-null results (PARTs D/E) are safe; significant ",
    "ones need corroboration. This one has it."),
    u$observed[1], u$null_mean[1], u$z[1], u$p_emp_two_sided[1], u$verdict[1],
    c2$observed[1], c2$null_mean[1], c2$z[1], c2$p_emp_two_sided[1], c2$verdict[1])
}

noncircular_verdict <- {
  a <- noncircular_stats |> dplyr::filter(grepl("^pro_comp", axis))
  b <- noncircular_stats |> dplyr::filter(grepl("NON-MITO", axis), grepl("z-composite", axis))
  nc <- noncircular_collapse_null
  # every axis that is NOT the mito-defined published one
  nnm <- noncircular_stats |>
    dplyr::filter(!grepl("^pro_comp|MITO-ONLY", axis)) |>
    dplyr::summarise(n_axes = dplyr::n(), n_neg = sum(geno_beta < 0),
                     n_sig = sum(geno_beta < 0 & geno_p_bh < 0.05))
  sprintf(paste0(
    "BREAKING THE CIRCULARITY -- AND THIS IS WHERE THE DEATH ARM ACTUALLY DIES. pro_comp IS a ",
    "mitochondrial gene set (MITOCARTA_APOPTOSIS_PRO, 25 genes), and so is the imbalance -- so ",
    "'the mitochondrial state couples to death priming' has been correlating mito genes with ",
    "mito genes. Rebuilt on the %d NON-mitochondrial pro-death genes of ",
    "cell_death_genes_consolidated, the Myc genotype effect is %+.3f (p=%.3f) -- against the ",
    "MITO axis's %+.3f (p=%.3f). The effect does NOT merely weaken; on non-mito death genes it ",
    "does not reproduce. Interaction: p=%.3f (non-mito) vs %.3f (mito). And the marginal ",
    "collapse of PART B evaporates the moment the shared mito membership is removed: the ",
    "circular mito axis gives perm p=%.3f, the non-mito pro-death axis p=%.3f, the non-mito ",
    "priming axis p=%.3f, CDC_PRODEATH_APOPTOSIS p=%.3f. The ONLY death axis whose coupling to ",
    "the imbalance collapses more than an arbitrary gene set's is the one that SHARES GENES ",
    "with the imbalance. That is what a circularity looks like when you test it. AND THE ",
    "DIRECTION IS THE OPPOSITE OF THE NARRATIVE: ALL %d non-mito/general death axes tested ",
    "have a NEGATIVE Myc genotype effect (%d of them significant after BH) -- APOP_BH3_REACTOME ",
    "%+.3f (p=%.3f) and CDC_PRODEATH_APOPTOSIS %+.3f (p=%.3f). On every death axis that is not ",
    "mito-defined, Myc LOWERS the death programme at 6W rather than raising it. READ THAT WITH ",
    "THE SURVIVOR BIAS, NOT AROUND IT: we sequence the cells that did NOT die, so the most ",
    "death-primed cells are missing by construction, and a suppressed-looking survivor pool is ",
    "exactly what oncogene-induced apoptosis would leave behind. The walkthrough calls survivor ",
    "bias 'a conservative floor' on Myc's death engagement; these numbers suggest it can INVERT ",
    "the sign, not merely shrink it. Either way the transcriptome cannot be cited as showing ",
    "Myc primes for death. See PART F3 (`mito_background_verdict`)."),
    length(ens_of(prodeath_nomito)),
    b$geno_beta[1], b$geno_p[1], a$geno_beta[1], a$geno_p[1], b$int_p[1], a$int_p[1],
    nc$perm_p[grepl("^pro_comp", nc$death_axis)][1],
    nc$perm_p[grepl("pro-death NON-MITO", nc$death_axis)][1],
    nc$perm_p[grepl("priming NON-MITO", nc$death_axis)][1],
    nc$perm_p[grepl("CDC_PRODEATH_APOPTOSIS", nc$death_axis)][1],
    nnm$n_axes[1], nnm$n_sig[1],
    noncircular_stats$geno_beta[grepl("BH3", noncircular_stats$axis)][1],
    noncircular_stats$geno_p[grepl("BH3", noncircular_stats$axis)][1],
    noncircular_stats$geno_beta[grepl("CDC_PRODEATH", noncircular_stats$axis)][1],
    noncircular_stats$geno_p[grepl("CDC_PRODEATH", noncircular_stats$axis)][1])
}

# =============================================================================
# PART G: A COMPETITIVE NULL FOR THE BRANCH-1 BINOMIAL
# =============================================================================
# 16:133 runs binom.test(supporting, total, p = 0.5) on the sign of
# dLFC = myc_6W_log2FC - myc_12W_log2FC. Under the GLOBAL attenuation (Issue #4),
# EVERY Myc-induced gene has dLFC > 0, so 50% is NOT the null. The right null is
# the genome-wide supporting fraction among genes matched on Myc direction and
# |LFC| at 6W. Note this can FLIP the read: 52.7% may sit BELOW the baseline.
supporting_frac <- function(ens, effect) {
  g <- gene_tbl |> dplyr::filter(ensembl %in% ens, abs(lfc6 - lfc12) > 0.2)
  if (!nrow(g)) return(NULL)
  d <- g$lfc6 - g$lfc12
  sup <- if (effect == "pro-death") d > 0 else d < 0
  tibble::tibble(effect = effect, n = length(sup), pct_supporting = 100 * mean(sup),
                 binom_p_vs_50 = stats::binom.test(sum(sup), length(sup), p = 0.5)$p.value)
}
genomewide_baseline <- {
  g <- gene_tbl |> dplyr::filter(abs(lfc6 - lfc12) > 0.2)
  d <- g$lfc6 - g$lfc12
  tibble::tibble(
    universe = "all expressed genes (|dLFC|>0.2)", n = length(d),
    pct_dLFC_positive = 100 * mean(d > 0),
    pct_dLFC_positive_among_myc_induced = 100 * mean((d > 0)[g$lfc6 > 0]),
    pct_dLFC_positive_among_myc_repressed = 100 * mean((d > 0)[g$lfc6 < 0]),
    reads_as = paste("THIS is the null the binomial should have used, not 50%.",
                     "Under a global attenuation a Myc-INDUCED gene has dLFC>0 by",
                     "construction, so 50% tests nothing about death."))
}
binomial_recheck <- dplyr::bind_rows(
  supporting_frac(ens_of(prodeath_all), "pro-death"),
  supporting_frac(ens_of(prodeath_nomito), "pro-death"),
  supporting_frac(ens_of(prosurv_all), "pro-survival")) |>
  dplyr::mutate(set = c(lab("pro-death ALL", ens_of(prodeath_all)),
                        lab("pro-death NON-MITO", ens_of(prodeath_nomito)),
                        lab("pro-survival ALL", ens_of(prosurv_all))), .before = 1) |>
  dplyr::mutate(
    genomewide_pct = genomewide_baseline$pct_dLFC_positive,
    excess_over_genomewide = pct_supporting - genomewide_baseline$pct_dLFC_positive,
    verdict = dplyr::if_else(excess_over_genomewide > 0,
                             "above the genome-wide baseline",
                             "AT OR BELOW the genome-wide baseline"))

binomial_verdict <- {
  pd <- binomial_recheck |> dplyr::filter(grepl("pro-death ALL", set))
  sprintf(paste0(
    "THE BRANCH-1 BINOMIAL'S NULL WAS NEVER 50%%, AND CORRECTING IT FLIPS THE READ. `16:133` ",
    "runs binom.test(supporting, total, p = 0.5) on the sign of dLFC = myc_6W - myc_12W. But ",
    "the Myc programme attenuates GLOBALLY, so dLFC > 0 is what a Myc-INDUCED gene does BY ",
    "CONSTRUCTION. Measured genome-wide: %.1f%% of expressed genes have dLFC>0 -- and among ",
    "Myc-INDUCED genes, %.1f%%. THAT is the null, not 50%%. Consequence: the %d pro-death genes ",
    "are %.1f%% 'supporting', which vs 50%% looks significant (binom p=%.3f) but sits %+.1f ",
    "points against the genome-wide baseline of %.1f%% => %s. So the pro-death set does not ",
    "front-load at 6W any more than an average gene does; it front-loads slightly LESS. Script ",
    "16's published 52.7%%/p=0.27 'clean negative' reached the right conclusion through the ",
    "wrong null -- and had the set scored a few points higher it would have reported a ",
    "significant death result that was nothing but the global attenuation."),
    genomewide_baseline$pct_dLFC_positive,
    genomewide_baseline$pct_dLFC_positive_among_myc_induced,
    pd$n[1], pd$pct_supporting[1], pd$binom_p_vs_50[1],
    pd$excess_over_genomewide[1], genomewide_baseline$pct_dLFC_positive, pd$verdict[1])
}

# =============================================================================
# PART H: VERDICTS
# =============================================================================
message("\n", paste(strwrap(ambient_verdict,    width = 88), collapse = "\n"))
message("\n", paste(strwrap(collapse_verdict,   width = 88), collapse = "\n"))
message("\n", paste(strwrap(priming_verdict,    width = 88), collapse = "\n"))
message("\n", paste(strwrap(atten_verdict,      width = 88), collapse = "\n"))
message("\n", paste(strwrap(tang_verdict,       width = 88), collapse = "\n"))
message("\n", paste(strwrap(gate2_verdict,      width = 88), collapse = "\n"))
message("\n", paste(strwrap(noncircular_verdict, width = 88), collapse = "\n"))
message("\n", paste(strwrap(mito_background_verdict, width = 88), collapse = "\n"))
message("\n", paste(strwrap(binomial_verdict, width = 88), collapse = "\n"), "\n")

# =============================================================================
# PART I: FIGURES
# =============================================================================

# A: the ambient diagnostic -- does EVERY coupling collapse?
amb_df <- ambient_by_timepoint |>
  dplyr::filter(method == "pearson") |>
  dplyr::mutate(window = factor(window, levels = c("6W", "12W")))
p_a <- ggplot2::ggplot(amb_df, ggplot2::aes(x = window, y = null_median_abs)) +
  ggplot2::geom_col(fill = "grey70", width = 0.6) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", null_median_abs)),
                     vjust = -0.4, size = 3.2) +
  ggplot2::labs(
    title = "Ambient coupling level: does EVERY gene set decouple from 6W to 12W?",
    subtitle = paste("Median |r| of the mitonuclear imbalance against ALL 884 library sets.",
                     "\nIf 12W is much lower, the published r 0.75 -> 0.07 collapse is ambient",
                     "structure, not a death-specific decoupling."),
    x = NULL, y = "median |r| vs all 884 sets") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "A_ambient_by_timepoint.pdf"), p_a, width = 5.5, height = 4)

# B: THE DELTA-RHO NULL -- the figure that carries the argument
obs_b <- collapse_null |> dplyr::filter(method == "pearson")
p_b <- ggplot2::ggplot(collapse_distributions, ggplot2::aes(x = delta)) +
  ggplot2::geom_histogram(bins = 50, fill = "grey75", colour = NA) +
  ggplot2::geom_vline(data = obs_b, ggplot2::aes(xintercept = delta_obs),
                      colour = "#D73027", linewidth = 0.9) +
  ggplot2::geom_vline(data = obs_b, ggplot2::aes(xintercept = null_median_delta),
                      colour = "#4575B4", linewidth = 0.6, linetype = 2) +
  ggplot2::facet_wrap(~ axis, ncol = 1, scales = "free_y") +
  ggplot2::labs(
    title = "Does the death coupling COLLAPSE more than an arbitrary gene set's does?",
    subtitle = paste("Grey = that axis's delta (r_6W - r_12W) against ALL 884 library sets.",
                     "RED = its published collapse vs the pro-death composite.",
                     "\nBLUE dashed = null median. Red inside grey => the collapse is what any",
                     "gene set does => it is not evidence of a death-specific decoupling."),
    x = "delta r  (r at 6W - r at 12W)", y = "number of gene sets") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "B_delta_rho_null.pdf"), p_b, width = 8, height = 7)

# D: the mechanism discriminator -- excess over matched attenuation
p_d <- atten_excess |>
  dplyr::mutate(set = stats::reorder(set, z),
                sig = ifelse(p_emp_two_sided < 0.05, "in excess (p<0.05)", "not in excess")) |>
  ggplot2::ggplot(ggplot2::aes(x = z, y = set, fill = sig)) +
  ggplot2::geom_col() +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  ggplot2::scale_fill_manual(values = c("in excess (p<0.05)" = "#D73027",
                                        "not in excess" = "grey70")) +
  ggplot2::labs(
    title = "Is the death attenuation MORE than the Myc programme fading?",
    subtitle = paste("z vs an expression x dispersion-matched null on the aligned per-gene",
                     "attenuation sign(LFC6)*(LFC6-LFC12).\nNull genes get the SAME alignment,",
                     "so the regression-to-the-mean bias cancels. z~0 => the death sets",
                     "attenuate exactly as much as\nmatched genes => priming fades because the",
                     "MYC programme fades, with no death-specific gate."),
    x = "z vs matched null", y = NULL, fill = NULL) +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "D_attenuation_excess.pdf"), p_d, width = 9, height = 4.5)

message("Figures written to ", out_dir)

# =============================================================================
# PART J: SAVE
# =============================================================================
death_out <- list(
  pc1_by_timepoint          = pc1_by_timepoint,
  ambient_by_timepoint      = ambient_by_timepoint,
  ambient_verdict           = ambient_verdict,
  collapse_null             = collapse_null,
  collapse_distributions    = collapse_distributions,
  fisher_z_collapse         = fisher_z_collapse,
  collapse_verdict          = collapse_verdict,
  priming_own_stats         = priming_own_stats,
  priming_group_means       = priming_group_means,
  interaction_percentile    = interaction_percentile,
  death_set_percentile      = death_set_percentile,
  method_concordance        = method_concordance,
  cross_sectional_limit     = cross_sectional_limit,
  priming_verdict           = priming_verdict,
  gene_tbl_params           = list(n_genes = nrow(gene_tbl), B = B,
                                   bins = "10 baseMean x 10 dispersion"),
  atten_excess              = atten_excess,
  death_decomp              = death_decomp,
  tang_shared_vs_specific   = tang_shared_vs_specific,
  tang_verdict              = tang_verdict,
  atten_verdict             = atten_verdict,
  gate2_selfcontained       = gate2_selfcontained,
  gate2_competitive         = gate2_competitive,
  gate2_effective_n         = gate2_effective_n,
  gate2_verdict             = gate2_verdict,
  noncircular_stats         = noncircular_stats,
  noncircular_coupling      = noncircular_coupling,
  noncircular_collapse_null = noncircular_collapse_null,
  noncircular_verdict       = noncircular_verdict,
  priming_vs_mito_background = priming_vs_mito_background,
  mito_background_verdict   = mito_background_verdict,
  genomewide_baseline       = genomewide_baseline,
  binomial_recheck          = binomial_recheck,
  binomial_verdict          = binomial_verdict,
  notes = paste(
    "Block B, author challenge 2026-07-17: 'the main question is whether an increased priming",
    "at 6W is indeed collapsing at 12W, and whether it correlates with mitochondrial",
    "reprioritisation'. Script 33's null was computed WITHIN 6W only, so it tests a",
    "CROSS-SECTIONAL claim and cannot speak to either. This script builds the longitudinal",
    "nulls. PART A = the ambient diagnostic (PC1 within 12W was never computed; if the ambient",
    "coupling level collapses, so does every set's). PART B = the DELTA-RHO null, the null the",
    "collapse claim never had. PART C surfaces the contrast that was FITTED, SAVED and NEVER",
    "READ (23$h1$state_tests: PRO int_p=0.177, priming int_p=0.263) while the narrative quoted",
    "group means ('~3x more at 6W'). PART D is the MECHANISM DISCRIMINATOR: is the death",
    "attenuation in EXCESS of the matched global attenuation? If not, priming fades because the",
    "MYC programme fades (canonical Evan/Lowe BH3-only biology) and no mitochondrial-state gate",
    "is needed. PART E replaces Gate 2's t-vs-mu=0 (its sole positive, p=0.023) with a matched",
    "null + the correlation-corrected effective n. PART F breaks the mito<->priming CIRCULARITY",
    "(pro_comp IS MitoCarta). PART G replaces the binomial's p=0.5 with the genome-wide",
    "supporting fraction. WHAT IS EXPECTED TO SURVIVE: the GENOTYPE main effect (Myc raises",
    "priming, p=0.037/0.038) -- genotype-based, so the genotype-independent IEG axis cannot",
    "bias it -- and the external death PHENOTYPE. CEILINGS: n=6/group, so a null interaction is",
    "UNSUPPORTED not REFUTED; the design is CROSS-SECTIONAL so no within-animal delta exists;",
    "the 884 library sets are mutually correlated so empirical ps are indicative, not exact.",
    "See docs/2026-07-17_evidence_audit_and_narrative.md."))
saveRDS(death_out, here::here("results", "death_priming_reassessment.rds"))
message("Saved results/death_priming_reassessment.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  dp <- readRDS(here::here("results", "death_priming_reassessment.rds"))

  # --- The seven headlines ---
  cat(strwrap(dp$ambient_verdict, 88), sep = "\n")
  cat(strwrap(dp$collapse_verdict, 88), sep = "\n")
  cat(strwrap(dp$priming_verdict, 88), sep = "\n")
  cat(strwrap(dp$atten_verdict, 88), sep = "\n")
  cat(strwrap(dp$tang_verdict, 88), sep = "\n")
  cat(strwrap(dp$gate2_verdict, 88), sep = "\n")
  cat(strwrap(dp$noncircular_verdict, 88), sep = "\n")

  # --- PART A: is the collapse just the ambient correlation dropping? ---
  dp$pc1_by_timepoint |> print()        # PC1 within 6W vs within 12W (43% vs ?)
  dp$ambient_by_timepoint |> print()    # median |r| to all 884 sets, per timepoint

  # --- PART B: the delta-rho null -- does the collapse beat an arbitrary set? ---
  dp$collapse_null |> as.data.frame() |> print()
  dp$fisher_z_collapse |> print()       # and is r6 vs r12 a difference at all?

  # --- PART C: the contrast that was on disk all along ---
  dp$priming_own_stats |> print()       # int_p 0.177 / 0.263 -- never read
  dp$priming_group_means |> print()     # what the narrative quoted instead
  dp$death_set_percentile |> as.data.frame() |> print()
  dp$interaction_percentile |> print()
  dp$method_concordance |> as.data.frame() |> print()
  cat(strwrap(dp$cross_sectional_limit$why, 88), sep = "\n")

  # --- PART D: THE MECHANISM DISCRIMINATOR ---
  dp$atten_excess |> as.data.frame() |> print()   # z ~ 0 => no death-specific attenuation
  dp$death_decomp |> as.data.frame() |> print()   # shared-developmental vs Myc-specific
  dp$tang_shared_vs_specific |> print()           # temporal_pos vs temporal_neg

  # --- PART E: Gate 2 asked properly ---
  dp$gate2_selfcontained |> print()     # the original p=0.023, reproduced
  dp$gate2_competitive |> as.data.frame() |> print()   # vs a matched null
  dp$gate2_effective_n |> as.data.frame() |> print()   # 23 genes -> effective n?

  # --- PART F: break the circularity (pro_comp IS MitoCarta) ---
  dp$noncircular_stats |> as.data.frame() |> print()
  dp$noncircular_coupling |> print()
  dp$noncircular_collapse_null |> print()   # only the CIRCULAR axis collapses
  # F3: is "Myc raises priming" a death effect, or script 32's mito-content effect?
  dp$priming_vs_mito_background |> as.data.frame() |> print()
  cat(strwrap(dp$mito_background_verdict, 88), sep = "\n")

  # --- PART G: the binomial's null should never have been 50% ---
  dp$genomewide_baseline |> as.data.frame() |> print()   # 59.6%, not 50% -- and 94% if Myc-induced
  dp$binomial_recheck |> as.data.frame() |> print()
  cat(strwrap(dp$binomial_verdict, 88), sep = "\n")

  list.files(here::here("outputs", "death_priming_reassessment"), pattern = "\\.pdf$")
}
