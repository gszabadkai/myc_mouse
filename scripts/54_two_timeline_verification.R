# =============================================================================
# 54_two_timeline_verification.R
# -----------------------------------------------------------------------------
# THE VERIFICATION PASS BEHIND FIGURE 1's CLOSING PARAGRAPH, AND THE HOME FOR THE
# NUMBERS IT QUOTES.
#
# Almost nothing here is a new analysis. PARTS A-C re-read fits that scripts 03,
# 40, 42 and 43 already made, put each against a rule fixed BEFORE the numbers
# were retrieved, and add the two matched nulls the two-timeline plane needs and
# nobody had computed. PART D is the one piece of new fitting, and it is a
# re-scaling: the coupling slopes on ONE scale, so that an adjusted and an
# unadjusted fit can be compared at all.
#
# WHY IT EXISTS. The closing paragraph rules out three alternative explanations
# for MYC no longer killing the adult gland, and leaves one hypothesis standing.
# An audit on 2026-09-21 found that four of its numbers had never been checked
# against a rule, and that two of its sentences rested on a choice nobody had
# made deliberately:
#
#   * THE COUPLING PANEL MIXED TWO FITS. Fig. 2I drew UNADJUSTED per-genotype
#     lines beside the ADJUSTED model's interaction p. The two fits disagree
#     about the wild-type slope -- and the -4.80 it drew and the -0.01 the
#     handoff quotes are not the same predictor (mitoPPS against `ox_rel`), not
#     the same response scale (z-scored against raw log2), and not the same
#     model. PART D puts every version on one scale so the choice is visible.
#   * THE RATIO PANEL IMPORTED ITS LINE from another estimator: script 44's
#     gene-level 0.487, against script 40's pathway-level 0.552, against the 0.55
#     constant hard-coded in scripts 42 and 44. PART B replaces it with a line
#     fitted INSIDE the ratio set -- the eight ratios that are not the one the
#     sentence is about.
#
# THE RULES, fixed 2026-09-21 before any number was retrieved and NOT adjusted
# afterwards. They are saved into the object as `rules`, so the panels and the
# dated note quote them from one place rather than re-typing them:
#
#   ON THE DIAGONAL   |interaction log2FC| < 0.20 AND raw interaction p > 0.05.
#                     RAW p throughout, not adjusted.
#   MYC-SPECIFIC      interaction < 0 and |interaction| > 0.20, whatever the sign
#                     of either arm. (Ruling 1: the first version keyed on the
#                     sign of y, which is wrong on a plane where the readable
#                     quantity is the distance BELOW the diagonal.)
#   CHECK 1           Bax and Bcl2l1 both on the diagonal -> "only the PUMA to
#                     BCL-XL balance behaved differently" stands. Either off it
#                     with Bbc3's sign -> that word is FALSE, flagged plainly.
#                     Off it with the opposite sign -> report the coordinates.
#   CHECK 2           Foxo3 MYC-specific -> "rose in the normal gland and did not
#                     rise under MYC".
#   CHECK 3           IQR of the six-week MYC effect across ALL NINE ratios
#                     >= 0.30 log2 -> the panel stays a scatter; below it, the
#                     panel becomes a retention plot.
#   ARMS              the same 0.20 magnitude, applied as declared. No
#                     arm-specific threshold was invented after the values were
#                     seen, and both rulers are reported (ruling 5).
#
# THE LICENCE, GENE BY GENE (ruling 2), because it is not the same for the four:
#   Bbc3    PRE-SPECIFIED. Named in advance from the PGC1a westerns; script 44's
#           `defs$pre_specified_genes`.
#   Bcl2l1  the fixed DENOMINATOR of the pre-specified PUMA:BCL-XL pair (42).
#   Bax     EXPLORATORY. Position is descriptive; no p-value for it in the text.
#   Foxo3   EXPLORATORY. Not pre-specified -- it was found in script 44's scan.
#           Position is descriptive; no p-value for it in the text.
#
# WHAT THIS CANNOT DO, said before the results. BATCH = TIMEPOINT (CLAUDE.md):
# both axes of the plane are the batch-confounded contrast and are DESCRIBED, not
# claimed; the distance from the diagonal is the batch-clean quantity, because
# genotype is balanced within each extraction batch. n = 6 per cell. PART D is a
# correlation among 24 animals and is RANKING, not confirmatory inference.
#
# ONE INPUT IS FLAGGED AND THE FLAG IS DELIBERATE. PART C resolves set membership
# through functions/reconcile_gene_symbols.R, which reads results/ortholog_table.rds
# -- written 2026-02-07, before 01_load_data.R's commits of 2026-02-14 and
# 2026-02-16, so the mtime rule flags it. Neither commit changed the biomaRt query
# (9f63f29 moved it behind a cache; 41af03c added Hallmark loading), it is the
# pinned mapping every downstream object on disk was built with, and re-querying
# biomaRt today would return a different Ensembl release and silently move set
# membership across scripts 15-53. PART C therefore proves its arms against script
# 43's saved values to 1e-9 rather than regenerating anything.
#
# Reads (every other input is FRESH by the mtime-vs-last-commit rule, 2026-09-21):
#   results/interaction_results.rds            (03) five contrasts, raw MLE, IHW padj
#   results/priming_arm_teb.rds                (42) $priming, $machinery,
#                                                   $exclusions$puma_inputs, $axis_scores
#   results/substrate_specificity_tradeoff.rds (43) $comparator, $wt_null, $defs,
#                                                   $tradeoff, $tradeoff_perm
#   results/collapse_module_ownership.rds      (44) $collapse_genes, $defs
#   results/background_vs_myc.rds              (40) $ruler, both rulers x 144 pathways
#   results/gate_model_verification.rds        (48) $scores, $simple_slopes
#   results/gsva_scores.rds                    (15) $expr_mat, to reconcile Fig. 2I's
#                                                   drawn numbers with this scale
#   data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt
# Writes: results/two_timeline_verification.rds
#
# RUNTIME: about two minutes. PART C draws 2000 matched sets for each of ten arms
# on three contrasts; PART D runs 5000 within-timepoint permutations twice.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

set.seed(54)
NSET    <- 2000L   # matched-random-set draws, script 43's count
NBIN    <- 20L     # baseMean bins, script 43's
NPERM   <- 5000L   # within-timepoint permutations, scripts 42/43's count
DIAG    <- 0.20    # the declared |interaction| threshold
ALPHA   <- 0.05    # the declared raw-p threshold
IQR_MIN <- 0.30    # check 3's threshold, in log2 units
GENES   <- c("Bbc3", "Foxo3", "Bax", "Bcl2l1")
TARGET  <- "Bbc3:Bcl2l1"          # the ratio the sentence is about
PCT_TOL <- 5                      # Monte Carlo tolerance, percentile points

rules <- c(
  on_diagonal  = "|interaction log2FC| < 0.20 AND raw interaction p > 0.05",
  myc_specific = "interaction < 0 and |interaction| > 0.20, whatever the sign of either arm",
  check1       = "Bax and Bcl2l1 both on the diagonal -> 'only' stands; either off it with Bbc3's sign -> 'only' is FALSE",
  check2       = "Foxo3 MYC-specific -> 'rose in the normal gland and did not rise under MYC'",
  check3       = "IQR of the six-week MYC effect across all nine ratios >= 0.30 log2 -> keep the scatter",
  arms         = "the same 0.20 magnitude; both rulers reported, mitoPPS drawn",
  p_kind       = "RAW p throughout. The adjusted p on disk is IHW (weighted Benjamini-Hochberg), not plain BH")

licence <- c(
  Bbc3   = "pre-specified (PGC1a westerns; script 44 defs$pre_specified_genes)",
  Bcl2l1 = "the fixed denominator of the pre-specified PUMA:BCL-XL pair (script 42)",
  Bax    = "exploratory -- position descriptive, no p-value quoted in the manuscript",
  Foxo3  = "exploratory -- not pre-specified, found in script 44's scan; position descriptive, no p-value quoted")

# =============================================================================
# PART 0: LOAD, ALIGN, AND REPRODUCE THE RECORD
# -----------------------------------------------------------------------------
# Every part below is checked against a number some other script already wrote.
# If a control fails, the load is wrong and nothing underneath it is readable.
# =============================================================================
message("54 PART 0: load")

ir  <- readRDS(here::here("results", "interaction_results.rds"))
pa  <- readRDS(here::here("results", "priming_arm_teb.rds"))
ss  <- readRDS(here::here("results", "substrate_specificity_tradeoff.rds"))
cmo <- readRDS(here::here("results", "collapse_module_ownership.rds"))
bv  <- readRDS(here::here("results", "background_vs_myc.rds"))
gm  <- readRDS(here::here("results", "gate_model_verification.rds"))
gs  <- readRDS(here::here("results", "gsva_scores.rds"))
gmt <- fgsea::gmtPathways(here::here("data", "genesets_from_library",
                                     "mammary_mito_myc_metab_v1_mouse.gmt"))

CONTRASTS <- c("myc_6W_raw", "myc_12W_raw", "timepoint_neg_raw", "timepoint_pos_raw",
               "interaction_raw")
D <- lapply(ir[CONTRASTS], as.data.frame)          # DESeqResults are S4; coerce once
universe_all <- rownames(D$interaction_raw)
stopifnot(all(vapply(D, function(x) identical(rownames(x), universe_all), logical(1))))

V  <- function(k, col = "log2FoldChange") stats::setNames(D[[k]][[col]], universe_all)
m6 <- V("myc_6W_raw");        m12 <- V("myc_12W_raw")
tn <- V("timepoint_neg_raw"); tp  <- V("timepoint_pos_raw")
it <- V("interaction_raw")
bm <- V("myc_6W_raw", "baseMean")

# --- CONTROL 1: the plane's geometry is exact, gene by gene -------------------
# The vertical distance from the diagonal IS the interaction. That is an identity
# of the saturated fit, not an approximation, and every panel on this plane
# depends on it, so it is measured over the whole transcriptome rather than
# asserted from the design.
identity_tbl <- tibble::tibble(
  identity = c("6>12W_myc = 6>12W_wt + interaction",
               "myc_12W = myc_6W + interaction"),
  n_genes  = c(sum(is.finite(tp - (tn + it))), sum(is.finite(m12 - (m6 + it)))),
  max_abs_dev = c(max(abs(tp - (tn + it)), na.rm = TRUE),
                  max(abs(m12 - (m6 + it)), na.rm = TRUE)))
stopifnot(all(identity_tbl$max_abs_dev < 1e-10),
          all(identity_tbl$n_genes == length(universe_all)))
message(sprintf("54 PART 0: the plane's identity holds to %.1e over %d genes",
                max(identity_tbl$max_abs_dev), length(universe_all)))

# --- the four genes, mapped twice --------------------------------------------
# Route 1 is the reconciler (CLAUDE.md's mandatory route for set membership);
# route 2 is script 44's saved symbol column, which came from the annotation
# table at ITS runtime. Two independent routes must agree or the mapping is not
# safe to quote.
ens_of_one <- function(g) {
  e <- recon_to_ensembl(g, universe_all)
  stopifnot(length(e) == 1L)
  e
}
ens <- vapply(GENES, ens_of_one, character(1))
cg  <- as.data.frame(cmo$collapse_genes)
stopifnot(identical(unname(ens), cg$ens[match(GENES, cg$gene)]))

# --- CONTROL 2: script 42's tables carry the same numbers ---------------------
# Script 42 pulled these four through a third symbol route. Same object, same
# contrast, so any disagreement is a mapping error rather than a rounding one.
s42 <- dplyr::bind_rows(
  as.data.frame(pa$machinery)[, c("gene", "lfc_wt_time", "lfc_myc_time", "lfc_interaction")],
  as.data.frame(pa$exclusions$puma_inputs)[, c("gene", "lfc_wt_time", "lfc_myc_time", "lfc_interaction")])
s42 <- s42[s42$gene %in% GENES, ]
stopifnot(setequal(s42$gene, GENES),
          max(abs(s42$lfc_wt_time     - tn[ens[s42$gene]])) < 1e-12,
          max(abs(s42$lfc_myc_time    - tp[ens[s42$gene]])) < 1e-12,
          max(abs(s42$lfc_interaction - it[ens[s42$gene]])) < 1e-12)

# =============================================================================
# PART A: THE FOUR GENES ON THE PLANE -- checks 1 and 2
# -----------------------------------------------------------------------------
# x is what the wild-type gland does across the window, y what the Myc+ gland
# does, and the interaction is the vertical distance from the diagonal. Raw
# (unshrunken) MLE log2 fold changes with their standard errors and RAW Wald p;
# `padj` is carried too and is IHW -- "Weighted BH adjusted p-values" is the
# object's own description of the column -- so it is named IHW wherever it is
# printed. Four panels and two documents had called it Benjamini-Hochberg.
# =============================================================================
message("54 PART A: the four genes")

pull_gene <- function(g) {
  e <- ens[[g]]
  a <- D$timepoint_neg_raw[e, ]; b <- D$timepoint_pos_raw[e, ]; i <- D$interaction_raw[e, ]
  tibble::tibble(
    gene = g, ensembl = e, baseMean = i$baseMean,
    wt_lfc  = a$log2FoldChange, wt_se  = a$lfcSE, wt_p  = a$pvalue, wt_padj_ihw  = a$padj,
    myc_lfc = b$log2FoldChange, myc_se = b$lfcSE, myc_p = b$pvalue, myc_padj_ihw = b$padj,
    int_lfc = i$log2FoldChange, int_se = i$lfcSE, int_p = i$pvalue, int_padj_ihw = i$padj,
    identity_dev = b$log2FoldChange - (a$log2FoldChange + i$log2FoldChange),
    licence = unname(licence[g]))
}
four_genes <- dplyr::bind_rows(lapply(GENES, pull_gene)) |>
  dplyr::mutate(
    on_diagonal  = abs(int_lfc) < DIAG & int_p > ALPHA,
    myc_specific = int_lfc < 0 & abs(int_lfc) > DIAG,
    verdict = dplyr::case_when(
      on_diagonal  ~ "on the diagonal",
      myc_specific ~ "below the diagonal, MYC-specific",
      int_lfc < 0  ~ "below the diagonal but inside the 0.20 magnitude; off the rule on p",
      TRUE         ~ "above the diagonal"))
stopifnot(max(abs(four_genes$identity_dev)) < 1e-12)

# --- CHECK 1, read against its rule and not against the draft -----------------
bbc3_sign <- sign(four_genes$int_lfc[four_genes$gene == "Bbc3"])
c1 <- four_genes[four_genes$gene %in% c("Bax", "Bcl2l1"), ]
c1_off_same <- c1$gene[!c1$on_diagonal & sign(c1$int_lfc) == bbc3_sign]
c1_off_opp  <- c1$gene[!c1$on_diagonal & sign(c1$int_lfc) != bbc3_sign]
check1 <- list(
  rule    = unname(rules[["check1"]]),
  on_diagonal = c1$gene[c1$on_diagonal],
  off_same_sign_as_Bbc3 = c1_off_same,
  off_opposite_sign     = c1_off_opp,
  verdict = if (length(c1_off_same))
    sprintf(paste("FAILS: %s sits off the diagonal with Bbc3's sign, so \"only the PUMA to",
                  "BCL-XL balance behaved differently\" is FALSE as written."),
            paste(c1_off_same, collapse = " and "))
  else if (length(c1_off_opp))
    sprintf("off the diagonal with the OPPOSITE sign (%s) -- coordinates reported, author decides",
            paste(c1_off_opp, collapse = " and "))
  else "PASSES: both on the diagonal, panel and sentence unchanged")

# --- CHECK 2 ------------------------------------------------------------------
# The pair of contrasts that both miss 0.05 is the GENOTYPE pair, not the two arms
# of this plane. Both pairs are carried, so the note cannot confuse them again.
f3   <- four_genes[four_genes$gene == "Foxo3", ]
g6f  <- D$myc_6W_raw[ens[["Foxo3"]], ]
g12f <- D$myc_12W_raw[ens[["Foxo3"]], ]
check2 <- list(
  rule    = unname(rules[["check2"]]),
  wt_lfc  = f3$wt_lfc, wt_p = f3$wt_p, myc_lfc = f3$myc_lfc, myc_p = f3$myc_p,
  int_lfc = f3$int_lfc, int_p = f3$int_p,
  geno_6W_lfc = g6f$log2FoldChange, geno_6W_p = g6f$pvalue,
  geno_12W_lfc = g12f$log2FoldChange, geno_12W_p = g12f$pvalue,
  myc_specific = f3$myc_specific,
  verdict = if (f3$myc_specific && f3$wt_lfc > 0)
    "PASSES on the amended rule: interaction below -0.20, wild-type arm positive -- 'rose in the normal gland and did not rise under MYC'"
  else if (f3$myc_specific)
    "MYC-specific by the interaction, but the wild-type arm is not positive -- the sentence needs rewording"
  else "NOT MYC-specific by the declared rule",
  note = sprintf(paste("The Myc+ arm is flat (%+.3f, p %.2f), not falling. The pair that both miss 0.05",
                       "is the GENOTYPE pair (6W %+.3f p %.3f, 12W %+.3f p %.3f): the significance",
                       "prohibition binds those, not the temporal arms (ruling 1)."),
                 f3$myc_lfc, f3$myc_p, g6f$log2FoldChange, g6f$pvalue,
                 g12f$log2FoldChange, g12f$pvalue))

# =============================================================================
# PART B: THE NINE RATIOS -- check 3, and the line fitted inside the set
# -----------------------------------------------------------------------------
# Script 42 fits each pro:anti pair per animal as log2(pro) - log2(anti) and takes
# the genotype coefficient within each age, so `d6` and `d12` are the MYC effect
# on that ratio at six and twelve weeks. Those are the saved values; nothing is
# refitted here.
#
# WHY THE LINE IS FITTED RATHER THAN IMPORTED (ruling 4). Three "programme-wide
# retention rates" exist -- 0.487 (script 44, gene-level DESeq2 fold changes),
# 0.552 (script 40, the mitochondrial compartment on the content ruler) and the
# 0.55 CONSTANT hard-coded in scripts 42 and 44 -- and the ratio retentions come
# from none of them: they are per-animal OLS coefficients. Comparing a ratio's
# retention with any of the three compares two estimators. A line fitted to the
# OTHER EIGHT ratios is internal to this estimator, and it makes the comparator
# the ratio set itself.
#
# NO CRITERION FOR "OFF THE LINE" WAS DECLARED IN ADVANCE, so this part reports
# the residual, the residual in units of the fit's own scatter, the rank of that
# residual among the nine, and a 95% prediction interval. It does not rule.
# =============================================================================
message("54 PART B: the nine ratios")

pr <- as.data.frame(pa$priming)
stopifnot(nrow(pr) == 9L, !anyDuplicated(pr$pair), TARGET %in% pr$pair)

# Bax's 0.55 is COMPUTED, not compared into existence: script 42 stores
# retention = d12 / d6 for every pair. The 0.55 constant enters that script only
# as a label threshold (`vs_global`, +/- 0.15) and as the target of the
# matched-pair null control -- never as the value of a retention.
stopifnot(max(abs(pr$retention - pr$d12 / pr$d6)) < 1e-12)
retention_provenance <- tibble::tibble(
  pair = pr$pair, d6 = pr$d6, d12 = pr$d12,
  retention_saved = pr$retention, retention_recomputed = pr$d12 / pr$d6,
  computed_not_compared = TRUE)

ratio_range <- tibble::tibble(
  n_ratios  = nrow(pr),
  median_d6 = stats::median(pr$d6),
  q1_d6     = unname(stats::quantile(pr$d6, 0.25)),
  q3_d6     = unname(stats::quantile(pr$d6, 0.75)),
  iqr_d6    = stats::IQR(pr$d6),
  threshold = IQR_MIN,
  verdict   = ifelse(stats::IQR(pr$d6) >= IQR_MIN,
                     "keep the scatter", "switch to a retention plot"),
  d6_Bax_BclxL  = pr$d6[pr$pair == "Bax:Bcl2l1"],
  d6_PUMA_BclxL = pr$d6[pr$pair == TARGET],
  # the two the old panel did not draw, and why
  not_drawn_before = paste(pr$pair[pr$anti == "Mcl1"], collapse = ", "),
  reason_not_drawn = "a second denominator needs a second encoding (fig2_priming_ratios.R:75-77)")

fit_origin <- function(d) stats::lm(d12 ~ 0 + d6, data = d)
others <- pr[pr$pair != TARGET, ]
tgt    <- pr[pr$pair == TARGET, ]
m8     <- fit_origin(others); s8 <- summary(m8)
pi8    <- stats::predict(m8, newdata = tgt, interval = "prediction", level = 0.95)

# every ratio's distance from that line, the target's out of sample
ratio_resid <- pr |>
  dplyr::mutate(fitted   = unname(stats::coef(m8)[1]) * d6,
                residual = d12 - fitted,
                in_fit   = pair != TARGET) |>
  dplyr::mutate(abs_rank = rank(-abs(residual)))

# the line with its 95% prediction band on a grid, for the panel to draw. The band
# is what the target is read against, so it is computed here, where the fit is,
# rather than re-derived in the figure layer.
xg <- seq(min(0, min(pr$d6)), max(pr$d6) * 1.04, length.out = 101)
ratio_band <- tibble::as_tibble(as.data.frame(
  stats::predict(m8, newdata = data.frame(d6 = xg), interval = "prediction", level = 0.95))) |>
  dplyr::mutate(d6 = xg, .before = 1)

# the comparator set contains ONE other PUMA ratio (Bbc3:Mcl1). Reported as a
# sensitivity, NOT drawn: the ruling names the other eight, and choosing a
# smaller comparator after seeing the answer is the move the pre-declaration
# exists to prevent.
o7  <- others[others$pro != "Bbc3", ]
m7  <- fit_origin(o7)
pi7 <- stats::predict(m7, newdata = tgt, interval = "prediction", level = 0.95)

ratio_line <- tibble::tibble(
  comparator   = c("the other eight ratios", "the seven that are not a PUMA ratio"),
  n_ratios     = c(nrow(others), nrow(o7)),
  slope        = c(unname(stats::coef(m8)[1]), unname(stats::coef(m7)[1])),
  slope_se     = c(s8$coefficients[1, 2], summary(m7)$coefficients[1, 2]),
  sigma        = c(s8$sigma, summary(m7)$sigma),
  df_resid     = c(m8$df.residual, m7$df.residual),
  # R2 as the SQUARED PEARSON CORRELATION, the definition Fig. 1G fixed for this
  # project; a through-origin lm reports an uncentred R2 that is not comparable
  # with anything else, and it is carried alongside rather than quoted.
  r2_pearson   = c(stats::cor(others$d6, others$d12)^2, stats::cor(o7$d6, o7$d12)^2),
  r2_uncentred = c(s8$r.squared, summary(m7)$r.squared),
  target_d6    = tgt$d6, target_d12 = tgt$d12,
  target_pred  = c(pi8[1, "fit"], pi7[1, "fit"]),
  target_resid = c(tgt$d12 - pi8[1, "fit"], tgt$d12 - pi7[1, "fit"]),
  target_resid_in_sigma = c((tgt$d12 - pi8[1, "fit"]) / s8$sigma,
                            (tgt$d12 - pi7[1, "fit"]) / summary(m7)$sigma),
  pi_lwr = c(pi8[1, "lwr"], pi7[1, "lwr"]), pi_upr = c(pi8[1, "upr"], pi7[1, "upr"]),
  target_inside_pi = c(tgt$d12 > pi8[1, "lwr"] & tgt$d12 < pi8[1, "upr"],
                       tgt$d12 > pi7[1, "lwr"] & tgt$d12 < pi7[1, "upr"]))

# =============================================================================
# PART C: THE ARMS, ON BOTH RULERS, WITH MATCHED NULLS ON EVERY CONTRAST THE
#         PLANE DRAWS
# -----------------------------------------------------------------------------
# The arms and their membership are script 43's, READ from its saved `defs` and
# resolved through the same reconciler, then proved against its saved values. The
# new part is the null: script 43 nulled the WILD-TYPE contrast only, so the
# plane's vertical axis and its diagonal distance had none. The three nulls are
# drawn from the SAME 2000 matched sets, so the x, y and interaction percentiles
# for an arm describe one set of draws rather than three.
#
# BOTH RULERS (ruling 5). The content ruler is a set-average log2 fold change, the
# same quantity the genes carry, so the arms and the four genes can share axes.
# mitoPPS is a within-compartment pairwise ratio and exists only for MitoCarta
# pathways: no proliferative arm, no TEB arm, and no expression-matched null.
# =============================================================================
message("54 PART C: the arms")

arms_def    <- as.data.frame(ss$defs$arms)         # arm, set, ruler_pathway
prolif_sets <- ss$defs$prolif_sets
stopifnot(all(arms_def$set %in% names(gmt)), all(prolif_sets %in% names(gmt)),
          length(prolif_sets) >= 5)

ens_set  <- function(syms) { e <- recon_to_ensembl(syms, universe_all); e[!is.na(e)] }
set_mean <- function(v, e) if (length(e)) mean(v[e], na.rm = TRUE) else NA_real_

arm_ens <- c(stats::setNames(lapply(arms_def$set, function(s) ens_set(gmt[[s]])), arms_def$arm),
             list("PROLIF_* pooled" = ens_set(unique(unlist(gmt[prolif_sets])))))

arms_content <- dplyr::bind_rows(lapply(names(arm_ens), function(a) {
  e <- arm_ens[[a]]
  tibble::tibble(arm = a, n_genes = length(e),
                 c_wt = set_mean(tn, e), c_myc = set_mean(tp, e), c_int = set_mean(it, e),
                 c_myc_6W = set_mean(m6, e), c_myc_12W = set_mean(m12, e))
}))

# --- CONTROL 3: script 43's comparator, to the digit --------------------------
cmp <- as.data.frame(ss$comparator)
k43 <- match(arms_content$arm, cmp$arm)
stopifnot(!anyNA(k43),
          all(arms_content$n_genes == cmp$n_genes[k43]),
          max(abs(arms_content$c_wt  - cmp$c_wt_time[k43]))  < 1e-9,
          max(abs(arms_content$c_myc - cmp$c_myc_time[k43])) < 1e-9,
          # the arm's distance from the diagonal is its mean interaction
          max(abs(arms_content$c_int - (arms_content$c_myc - arms_content$c_wt))) < 1e-9)

# --- the matched null, script 43's idiom, on three contrasts at once ----------
expressed <- names(bm)[is.finite(bm) & bm > 0 & is.finite(tn[names(bm)])]
bin_of <- cut(rank(bm[expressed], ties.method = "first"), breaks = NBIN, labels = FALSE)
by_bin <- split(expressed, bin_of)
draw_matched <- function(e) {
  b <- bin_of[match(e, expressed)]
  b <- b[!is.na(b)]
  unlist(lapply(split(b, b), function(k)
    sample(by_bin[[as.character(k[1])]], length(k), replace = TRUE)), use.names = FALSE)
}

arm_null <- dplyr::bind_rows(lapply(names(arm_ens), function(a) {
  e   <- arm_ens[[a]][arm_ens[[a]] %in% expressed]
  obs <- c(set_mean(tn, e), set_mean(tp, e), set_mean(it, e))
  nul <- vapply(seq_len(NSET), function(i) {
    r <- draw_matched(e)
    c(set_mean(tn, r), set_mean(tp, r), set_mean(it, r))
  }, numeric(3))
  tibble::tibble(
    arm = a, n_matched = length(e),
    obs_wt  = obs[1], null_med_wt  = stats::median(nul[1, ]), pct_wt  = 100 * mean(nul[1, ] < obs[1]),
    obs_myc = obs[2], null_med_myc = stats::median(nul[2, ]), pct_myc = 100 * mean(nul[2, ] < obs[2]),
    obs_int = obs[3], null_med_int = stats::median(nul[3, ]), pct_int = 100 * mean(nul[3, ] < obs[3]),
    null_int_lo = unname(stats::quantile(nul[3, ], 0.025)),
    null_int_hi = unname(stats::quantile(nul[3, ], 0.975)))
}))

# --- CONTROL 4: the wild-type half must reproduce script 43 -------------------
# Different draws from the same null, so the observed values must be identical
# and the percentiles must agree within Monte Carlo error.
w43 <- as.data.frame(ss$wt_null)
kw  <- match(arm_null$arm, w43$arm)
stopifnot(!anyNA(kw),
          max(abs(arm_null$obs_wt - w43$observed_wt[kw])) < 1e-9,
          max(abs(arm_null$pct_wt - w43$percentile[kw])) <= PCT_TOL)
message(sprintf("54 PART C: wild-type percentiles reproduce script 43 to %.2f points",
                max(abs(arm_null$pct_wt - w43$percentile[kw]))))

# --- the mitoPPS ruler, where the arm is a MitoPathway ------------------------
rl <- as.data.frame(bv$ruler)
mito_arms <- arms_def[!is.na(arms_def$ruler_pathway), ]
arms_mitopps <- dplyr::bind_rows(lapply(seq_len(nrow(mito_arms)), function(i) {
  r <- rl[rl$pathway == mito_arms$ruler_pathway[i], ]
  stopifnot(nrow(r) == 1L)
  tibble::tibble(arm = mito_arms$arm[i], pathway = r$pathway, n_genes = r$n_genes,
                 p_wt = r$p_tn, p_myc = r$p_tp, p_int = r$p_int,
                 c_wt_ruler = r$c_tn, c_myc_ruler = r$c_tp, c_int_ruler = r$c_int)
}))
# the same identity on the mitoPPS ruler; script 43's own mitoPPS read of these
# arms (CONTROL 4b); and script 40's content values must agree with the
# membership used here on the arm the sentence is about
cp43 <- as.data.frame(ss$comparator_priority)
kp   <- match(arms_mitopps$arm, cp43$arm)
stopifnot(max(abs(arms_mitopps$p_myc - (arms_mitopps$p_wt + arms_mitopps$p_int))) < 1e-9,
          max(abs(arms_mitopps$c_myc_ruler - (arms_mitopps$c_wt_ruler + arms_mitopps$c_int_ruler))) < 1e-9,
          !anyNA(kp),
          max(abs(arms_mitopps$p_wt  - cp43$prio_wt_time[kp]))  < 1e-12,
          max(abs(arms_mitopps$p_myc - cp43$prio_myc_time[kp])) < 1e-12,
          abs(arms_mitopps$c_wt_ruler[arms_mitopps$arm == "OXPHOS subunits"] -
              arms_content$c_wt[arms_content$arm == "OXPHOS subunits"]) < 1e-6)

# --- the declared 0.20, applied to the arms on both rulers --------------------
# "OXPHOS (all)" is the union of the two OXPHOS arms already on the plane, so it
# is carried in the table and NOT drawn -- the same genes would be plotted twice.
arm_diagonal <- dplyr::full_join(
  arms_content |> dplyr::select(arm, content_int = c_int),
  arms_mitopps |> dplyr::select(arm, mitopps_int = p_int), by = "arm") |>
  dplyr::left_join(arm_null |> dplyr::select(arm, content_int_pct = pct_int), by = "arm") |>
  dplyr::mutate(
    content_on_diagonal = abs(content_int) < DIAG,
    mitopps_on_diagonal = ifelse(is.na(mitopps_int), NA, abs(mitopps_int) < DIAG),
    content_margin = DIAG - abs(content_int),
    mitopps_margin = DIAG - abs(mitopps_int),
    draw = arm != "OXPHOS (all)")

# the shared axis limits for the two panels that must be on identical axes: the
# arms on the content ruler (with their null medians, which that panel draws) and
# the four genes. Computed once, here, so neither panel can pick its own.
drawn_arm <- function(a) a != "OXPHOS (all)"
plane_values  <- c(arms_content$c_wt[drawn_arm(arms_content$arm)],
                   arms_content$c_myc[drawn_arm(arms_content$arm)],
                   arm_null$null_med_wt[drawn_arm(arm_null$arm)],
                   arm_null$null_med_myc[drawn_arm(arm_null$arm)],
                   four_genes$wt_lfc, four_genes$myc_lfc)
plane_limits  <- c(-1, 1) * max(abs(plane_values)) * 1.08
mitopps_values <- c(arms_mitopps$p_wt[drawn_arm(arms_mitopps$arm)],
                    arms_mitopps$p_myc[drawn_arm(arms_mitopps$arm)])
plane_limits_mitopps <- c(-1, 1) * max(abs(mitopps_values)) * 1.08

# =============================================================================
# PART D: THE COUPLING, UNADJUSTED AND ADJUSTED, ON ONE SCALE
# -----------------------------------------------------------------------------
# THE PROBLEM THIS PART SOLVES. Three numbers were on record for the same
# relationship and no two of them were comparable:
#   -4.80   Fig. 2I's drawn wild-type line: SD of the VST ratio per UNIT of
#           mitoPPS, no covariates, 12 animals.
#   -2.19   the same fit on the unscaled log2 ratio (script 48 `simple_slopes`).
#   -0.01   a DIFFERENT predictor (`ox_rel`) on the unscaled log2 ratio, adjusted
#           for epithelial and immune composition, 12 animals.
# Here everything is on ONE scale -- SD of the ratio per SD of the predictor,
# both taken over all 24 animals -- and three fits are reported for each axis:
#   U  unadjusted, within genotype (identical to the pooled interaction model
#      without covariates)
#   W  epi + imm within genotype, each genotype with its own covariate
#      coefficients: script 48's form, 8 residual df
#   P  pooled y ~ myc * a + epi + imm, covariate coefficients SHARED: the
#      pre-specified model, 18 residual df, and the one whose interaction the
#      text quotes
# The panel draws U and P, because the difference of P's two slopes IS the
# interaction the sentence is built on (ruling 7).
# =============================================================================
message("54 PART D: the coupling fits")

S   <- gm$scores
ax  <- as.data.frame(pa$axis_scores)
pur <- as.data.frame(pa$purity)
# VALUES are compared with attributes ignored, and ORDER is checked separately.
# Script 48 saved its composites as vectors NAMED by sample; script 42 saved them
# as tibble columns, and once tibble is loaded as.data.frame() drops those names --
# so a plain all.equal() fails on "names for target but not for current" while the
# numbers are identical. (It passed in a session without tibble attached, which is
# how this slipped through; caught on the author's first run, 2026-09-21.)
same_values <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b),
                                               check.attributes = FALSE))
stopifnot(identical(as.character(ax$sample), colnames(gs$expr_mat)),
          identical(as.character(pur$sample), as.character(ax$sample)),
          identical(names(S$epi), as.character(ax$sample)),      # the order, explicitly
          same_values(S$X$ox_ppd, ax$oxphos_ppd),
          same_values(S$epi, pur$epithelial),
          same_values(S$imm, pur$immune),
          identical(levels(S$myc), c("neg", "pos")))

Z     <- function(x) as.numeric(scale(x))
y_log <- S$Y$PUMA                     # log2(normalised count + 1) ratio, raw
y_z   <- Z(y_log)
AXES  <- list(ox_ppd = S$X$ox_ppd, ox_rel = S$X$ox_rel, redox_ppd = ax$redox_ppd)
n_geno <- c(neg = sum(S$myc == "neg"), pos = sum(S$myc == "pos"))
stopifnot(all(n_geno == 12L))

# --- CONTROL 5: script 43's interaction, refitted from these vectors ----------
tr43 <- as.data.frame(ss$tradeoff)
rec  <- tr43[grepl("^Bbc3:Bcl2l1", tr43$outcome) & tr43$axis == "oxphos_ppd", ]
stopifnot(nrow(rec) == 1L)
d43 <- data.frame(y = y_z, m = S$myc, a = S$X$ox_ppd, epi = S$epi, imm = S$imm)
c43 <- summary(stats::lm(y ~ m * a + epi + imm, d43))$coefficients
stopifnot(abs(c43["mpos:a", 1] - rec$myc_x_axis) < 1e-9,
          abs(c43["mpos:a", 4] - rec$p) < 1e-9)
# carried into the object so the panel quotes script 43's timepoint-term p from
# here rather than typing it
tradeoff_43 <- tibble::as_tibble(rec)
message(sprintf("54 PART D: script 43's interaction reproduces at %+.4f (p %.5f)",
                c43["mpos:a", 1], c43["mpos:a", 4]))

# --- CONTROL 6: script 48's simple slopes, in their own units -----------------
simple48 <- as.data.frame(gm$simple_slopes)
rebuilt48 <- dplyr::bind_rows(lapply(seq_len(nrow(simple48)), function(i) {
  an <- simple48$axis[i]; cv <- simple48$covariates[i]
  d  <- data.frame(y = y_log, a = AXES[[an]], epi = S$epi, imm = S$imm, m = S$myc)
  f  <- if (cv == "none") "y ~ a" else "y ~ a + epi + imm"
  w  <- summary(stats::lm(stats::as.formula(f), d[d$m == "neg", ]))$coefficients
  p  <- summary(stats::lm(stats::as.formula(f), d[d$m == "pos", ]))$coefficients
  tibble::tibble(axis = an, covariates = cv, slope_wt = w["a", 1], slope_myc = p["a", 1])
}))
stopifnot(max(abs(rebuilt48$slope_wt  - simple48$slope_wt))  < 1e-9,
          max(abs(rebuilt48$slope_myc - simple48$slope_myc)) < 1e-9)

# --- CONTROL 7: the numbers Fig. 2I actually drew -----------------------------
# Rebuilt from the VST matrix, which is what that panel used, so the old record
# and this table can be read against each other.
E     <- gs$expr_mat
y_vst <- Z(as.numeric(E["Bbc3", ] - E["Bcl2l1", ]))
old_fig2I <- dplyr::bind_rows(lapply(c("oxphos_ppd", "redox_ppd"), function(an) {
  a <- if (an == "oxphos_ppd") S$X$ox_ppd else ax$redox_ppd
  dplyr::bind_rows(lapply(c("neg", "pos"), function(g) {
    k <- S$myc == g
    s <- summary(stats::lm(y_vst[k] ~ a[k]))
    tibble::tibble(axis = an, genotype = g, slope = s$coefficients[2, 1],
                   se = s$coefficients[2, 2], p = s$coefficients[2, 4],
                   r2 = s$r.squared, n = sum(k))
  }))
}))
OLD_DRAWN <- c(neg = -4.80, pos = 4.35)     # as PANELS.md and the Fig. 2I legend record them
stopifnot(abs(old_fig2I$slope[old_fig2I$axis == "oxphos_ppd" & old_fig2I$genotype == "neg"] -
              OLD_DRAWN[["neg"]]) < 0.01,
          abs(old_fig2I$slope[old_fig2I$axis == "oxphos_ppd" & old_fig2I$genotype == "pos"] -
              OLD_DRAWN[["pos"]]) < 0.01)

# --- the fit table, one scale -------------------------------------------------
lincom <- function(mdl, which_terms) {
  cf <- stats::coef(mdl); V <- stats::vcov(mdl)
  L  <- as.numeric(names(cf) %in% which_terms)
  est <- sum(L * cf); se <- sqrt(drop(t(L) %*% V %*% L))
  c(est = est, se = se, p = 2 * stats::pt(-abs(est / se), df = mdl$df.residual))
}

coupling_fits <- dplyr::bind_rows(lapply(names(AXES), function(an) {
  az <- Z(AXES[[an]])
  within <- dplyr::bind_rows(lapply(c("neg", "pos"), function(g) {
    k <- S$myc == g
    u <- summary(stats::lm(y_z[k] ~ az[k]))
    w <- summary(stats::lm(y_z[k] ~ az[k] + S$epi[k] + S$imm[k]))
    dplyr::bind_rows(
      tibble::tibble(axis = an, fit = "U unadjusted, within genotype", genotype = g,
                     slope = u$coefficients[2, 1], se = u$coefficients[2, 2],
                     p = u$coefficients[2, 4], r2 = u$r.squared,
                     n = sum(k), df_resid = u$df[2]),
      tibble::tibble(axis = an, fit = "W epi + imm, within genotype", genotype = g,
                     slope = w$coefficients[2, 1], se = w$coefficients[2, 2],
                     p = w$coefficients[2, 4], r2 = w$r.squared,
                     n = sum(k), df_resid = w$df[2]))
  }))
  d <- data.frame(y = y_z, m = S$myc, a = az, epi = S$epi, imm = S$imm)
  pooled <- dplyr::bind_rows(lapply(c("none", "epi + imm"), function(cv) {
    mdl <- if (cv == "none") stats::lm(y ~ m * a, d) else stats::lm(y ~ m * a + epi + imm, d)
    lw <- lincom(mdl, "a"); lp <- lincom(mdl, c("a", "mpos:a"))
    tibble::tibble(axis = an,
                   fit = paste0(ifelse(cv == "none", "U pooled", "P pooled"),
                                ", shared covariates: ", cv),
                   genotype = c("neg", "pos"),
                   slope = c(lw[["est"]], lp[["est"]]), se = c(lw[["se"]], lp[["se"]]),
                   p = c(lw[["p"]], lp[["p"]]), r2 = summary(mdl)$r.squared,
                   n = 12L, df_resid = mdl$df.residual)
  }))
  dplyr::bind_rows(within, pooled)
}))

# the pooled model with no covariates must reproduce the within-genotype slopes
u_within <- coupling_fits[coupling_fits$fit == "U unadjusted, within genotype", ]
u_pooled <- coupling_fits[grepl("^U pooled", coupling_fits$fit), ]
stopifnot(max(abs(u_within$slope[order(u_within$axis, u_within$genotype)] -
                  u_pooled$slope[order(u_pooled$axis, u_pooled$genotype)])) < 1e-10)

coupling_interaction <- dplyr::bind_rows(lapply(names(AXES), function(an) {
  az <- Z(AXES[[an]])
  d  <- data.frame(y = y_z, m = S$myc, a = az, epi = S$epi, imm = S$imm)
  dplyr::bind_rows(lapply(c("none", "epi + imm"), function(cv) {
    mdl <- if (cv == "none") stats::lm(y ~ m * a, d) else stats::lm(y ~ m * a + epi + imm, d)
    s <- summary(mdl)$coefficients
    tibble::tibble(axis = an, covariates = cv, interaction = s["mpos:a", 1],
                   se = s["mpos:a", 2], p = s["mpos:a", 4], df_resid = mdl$df.residual)
  }))
}))

# --- WHY the wild-type slope moves: the omitted-variable identity -------------
# b_unadjusted - b_adjusted = sum_k (k's coefficient in the adjusted fit) x
# (k's own slope on the axis). It is an identity, so it does not explain the
# move so much as say which covariate carries it.
coupling_decomp <- dplyr::bind_rows(lapply(c("ox_ppd", "ox_rel", "redox_ppd"), function(an) {
  az <- Z(AXES[[an]])
  dplyr::bind_rows(lapply(c("neg", "pos"), function(g) {
    k <- S$myc == g
    bu <- unname(stats::coef(stats::lm(y_z[k] ~ az[k]))[2])
    fa <- stats::coef(stats::lm(y_z[k] ~ az[k] + S$epi[k] + S$imm[k]))
    de <- unname(stats::coef(stats::lm(S$epi[k] ~ az[k]))[2])
    di <- unname(stats::coef(stats::lm(S$imm[k] ~ az[k]))[2])
    tibble::tibble(axis = an, genotype = g, slope_unadjusted = bu,
                   slope_adjusted = unname(fa[2]), difference = bu - unname(fa[2]),
                   via_epi = unname(fa[3]) * de, via_imm = unname(fa[4]) * di)
  }))
})) |>
  dplyr::mutate(check = difference - (via_epi + via_imm))
stopifnot(max(abs(coupling_decomp$check)) < 1e-10)

coupling_cor <- dplyr::bind_rows(lapply(c("neg", "pos"), function(g) {
  k <- S$myc == g
  M <- cbind(ratio = y_log[k], ox_ppd = S$X$ox_ppd[k], ox_rel = S$X$ox_rel[k],
             redox_ppd = ax$redox_ppd[k], epi = S$epi[k], imm = S$imm[k],
             tp12 = as.numeric(S$tp[k] == "12W"))
  cc <- stats::cor(M)
  tibble::tibble(genotype = g, pair = paste(rownames(cc)[row(cc)[upper.tri(cc)]],
                                            colnames(cc)[col(cc)[upper.tri(cc)]], sep = " ~ "),
                 r = cc[upper.tri(cc)])
}))

# --- the within-timepoint permutation null, both specifications ---------------
# Shuffling the axis INSIDE each timepoint keeps the design, the genotype split
# and the cohort contrast, and breaks only the animal-to-animal pairing -- which
# is what the coupling asserts. Its median is not zero. (Fig. 2I's legend said
# script 43 permutes the GENOTYPE LABELS; it does not, and that is corrected
# wherever it is repeated.)
perm_int <- function(az, covars) {
  d <- data.frame(y = y_z, m = S$myc, a = az, epi = S$epi, imm = S$imm)
  f <- if (covars) y ~ m * a + epi + imm else y ~ m * a
  obs <- summary(stats::lm(f, d))$coefficients["mpos:a", 1]
  nul <- vapply(seq_len(NPERM), function(i) {
    dd <- d
    for (tv in levels(S$tp)) { k <- which(S$tp == tv); dd$a[k] <- sample(dd$a[k]) }
    summary(stats::lm(f, dd))$coefficients["mpos:a", 1]
  }, numeric(1))
  c(observed = obs, null_median = stats::median(nul),
    percentile = 100 * mean(nul < obs), p_emp = mean(nul >= obs))
}
coupling_perm <- dplyr::bind_rows(lapply(c("none", "epi + imm"), function(cv) {
  v <- perm_int(Z(AXES$ox_ppd), cv != "none")
  tibble::tibble(axis = "ox_ppd", covariates = cv, n_perm = NPERM,
                 observed = v[["observed"]], null_median = v[["null_median"]],
                 percentile = v[["percentile"]], p_emp = v[["p_emp"]])
}))
# the adjusted one is script 43's own null, recomputed: Monte Carlo agreement only
p43 <- as.data.frame(ss$tradeoff_perm)
p43 <- p43[grepl("^Bbc3:Bcl2l1", p43$outcome) & p43$axis == "oxphos_ppd", ]
stopifnot(abs(coupling_perm$percentile[coupling_perm$covariates == "epi + imm"] -
              p43$percentile) <= PCT_TOL)

# --- what the panel draws -----------------------------------------------------
# The adjusted sub-panel draws the POOLED model, because the difference between
# its two slopes IS the interaction the sentence is built on. Its points are
# partial residuals: the response with the shared covariate terms removed, the
# covariates centred so the panel keeps the response's own units. Within each
# genotype the ordinary least-squares line through those points IS the model's
# slope for that genotype -- an orthogonality identity, asserted below, and the
# reason this is a drawing of the fit rather than a second fit.
az_ox <- Z(AXES$ox_ppd)
d_ox  <- data.frame(y = y_z, m = S$myc, a = az_ox, epi = S$epi, imm = S$imm)
mP    <- stats::lm(y ~ m * a + epi + imm, d_ox)
bP    <- stats::coef(mP)
y_adj <- y_z - bP[["epi"]] * (S$epi - mean(S$epi)) - bP[["imm"]] * (S$imm - mean(S$imm))

coupling_panel <- tibble::tibble(
  sample = as.character(ax$sample), group = as.character(S$group),
  genotype = as.character(S$myc), timepoint = as.character(S$tp),
  x = az_ox, y_unadjusted = y_z, y_adjusted = y_adj)

coupling_lines <- dplyr::bind_rows(lapply(c("unadjusted", "adjusted"), function(fit) {
  dplyr::bind_rows(lapply(c("neg", "pos"), function(g) {
    k <- S$myc == g
    yy <- if (fit == "unadjusted") y_z[k] else y_adj[k]
    cf <- stats::coef(stats::lm(yy ~ az_ox[k]))
    tibble::tibble(fit = fit, genotype = g, intercept = unname(cf[1]), slope = unname(cf[2]))
  }))
}))
# the drawn adjusted lines ARE the pooled model's slopes
pooled_slopes <- c(neg = unname(bP[["a"]]), pos = unname(bP[["a"]] + bP[["mpos:a"]]))
drawn_adj <- coupling_lines$slope[coupling_lines$fit == "adjusted"]
names(drawn_adj) <- coupling_lines$genotype[coupling_lines$fit == "adjusted"]
stopifnot(max(abs(drawn_adj[c("neg", "pos")] - pooled_slopes[c("neg", "pos")])) < 1e-10)

# =============================================================================
# PART E: THE OPEN ITEM, RECORDED AND NOT PURSUED
# -----------------------------------------------------------------------------
# The unadjusted Myc+ slope is largely a difference between the two ages rather
# than animal-to-animal variation inside them: adjusted for timepoint alone it
# shrinks by more than half and loses significance. That bears on whether this
# coupling is a within-animal relationship or a group difference, and the
# author's instruction (2026-09-21) is to record it here and stop. The two
# composites are NOT standing in for timepoint: the table says how little of
# either one timepoint explains within a genotype.
# =============================================================================
message("54 PART E: the timepoint reading, recorded")

coupling_timepoint <- dplyr::bind_rows(lapply(c("ox_ppd", "ox_rel"), function(an) {
  az <- Z(AXES[[an]])
  dplyr::bind_rows(lapply(c("neg", "pos"), function(g) {
    k <- S$myc == g; tpk <- droplevels(S$tp[k])
    t1 <- summary(stats::lm(y_z[k] ~ az[k] + tpk))$coefficients
    t2 <- summary(stats::lm(y_z[k] ~ az[k] + tpk + S$epi[k] + S$imm[k]))$coefficients
    tibble::tibble(axis = an, genotype = g,
                   slope_with_tp = t1[2, 1], p_with_tp = t1[2, 4],
                   slope_with_tp_epi_imm = t2[2, 1], p_with_tp_epi_imm = t2[2, 4],
                   r2_epi_on_tp = summary(stats::lm(S$epi[k] ~ tpk))$r.squared,
                   r2_imm_on_tp = summary(stats::lm(S$imm[k] ~ tpk))$r.squared)
  }))
}))

# =============================================================================
# ASSERTS -- the claims this object is allowed to carry
# =============================================================================
stopifnot(
  nrow(four_genes) == 4L, !anyNA(four_genes$int_p),
  nrow(pr) == 9L, nrow(ratio_line) == 2L,
  nrow(arms_content) == length(arm_ens), nrow(arm_null) == length(arm_ens),
  nrow(coupling_panel) == 24L, nrow(coupling_lines) == 4L,
  # the three fits must be present for the axis the panel draws
  sum(coupling_fits$axis == "ox_ppd") == 8L,
  # and the control axis must be there, or the panel has no comparator
  "redox_ppd" %in% coupling_fits$axis)

# Every number in the notes is formatted from the object, never typed: a number
# copied into prose is never checked, and this object exists to stop that.
fx  <- function(axis, fit, g) coupling_fits$slope[coupling_fits$axis == axis &
                                                  coupling_fits$fit == fit &
                                                  coupling_fits$genotype == g]
FU  <- "U unadjusted, within genotype"
FP  <- "P pooled, shared covariates: epi + imm"
ixn <- function(axis, cv) coupling_interaction[coupling_interaction$axis == axis &
                                                 coupling_interaction$covariates == cv, ]
dwt <- coupling_decomp[coupling_decomp$axis == "ox_ppd" & coupling_decomp$genotype == "neg", ]
ccr <- function(g, pr_) coupling_cor$r[coupling_cor$genotype == g & coupling_cor$pair == pr_]
bx  <- four_genes[four_genes$gene == "Bax", ]
l8  <- ratio_line[1, ]

notes <- c(
  "WHAT CHANGED ON 2026-09-21, and none of it is a new analysis:",
  sprintf(" 1. CHECK 1, read against its rule and not against the draft: %s", check1$verdict),
  sprintf("    Bax: interaction %+.3f (SE %.3f, raw p %.3f) against the declared 0.20. Bax is",
          bx$int_lfc, bx$int_se, bx$int_p),
  "    EXPLORATORY -- no p-value for it in the text -- so what the panel shows is a POSITION.",
  sprintf(" 2. CHECK 2: %s", check2$verdict),
  sprintf("    %s", check2$note),
  sprintf(" 3. CHECK 3: IQR of the six-week effect over all %d ratios is %.3f against %.2f -> %s.",
          ratio_range$n_ratios, ratio_range$iqr_d6, IQR_MIN, ratio_range$verdict),
  sprintf(paste(" 4. THE RATIO LINE IS INTERNAL: slope %.3f through the origin over the other %d",
                "ratios. PUMA:BCL-XL sits %+.3f from it (%.2f of the line's residual SD, rank %d",
                "of 9 by size) and %s its 95%% prediction interval [%+.3f, %+.3f]. No criterion",
                "for 'off the line' was declared, and none is invented here."),
          l8$slope, l8$n_ratios, l8$target_resid, l8$target_resid_in_sigma,
          as.integer(ratio_resid$abs_rank[ratio_resid$pair == TARGET]),
          ifelse(l8$target_inside_pi, "INSIDE", "OUTSIDE"), l8$pi_lwr, l8$pi_upr),
  sprintf(paste(" 5. THE COUPLING IS BUILT ON THE INTERACTION (ruling 7): %+.3f unadjusted (p %.4f),",
                "%+.3f adjusted (p %.4f), SD of the ratio per SD of the mitoPPS axis. The per-genotype",
                "slopes are a drawing choice, never specified either way: wild type %+.3f -> %+.3f,",
                "Myc+ %+.3f -> %+.3f. The adjustment moves BOTH."),
          ixn("ox_ppd", "none")$interaction, ixn("ox_ppd", "none")$p,
          ixn("ox_ppd", "epi + imm")$interaction, ixn("ox_ppd", "epi + imm")$p,
          fx("ox_ppd", FU, "neg"), fx("ox_ppd", FP, "neg"),
          fx("ox_ppd", FU, "pos"), fx("ox_ppd", FP, "pos")),
  "",
  "THE THREE NUMBERS THAT WERE NEVER COMPARABLE (see `old_fig2I` and control 6):",
  sprintf("  %+.2f was SD of the VST ratio per UNIT of mitoPPS, unadjusted (Fig. 2I, as drawn);",
          old_fig2I$slope[old_fig2I$axis == "oxphos_ppd" & old_fig2I$genotype == "neg"]),
  sprintf("  %+.2f is that fit on the unscaled log2 ratio; %+.3f is a DIFFERENT predictor",
          simple48$slope_wt[simple48$axis == "ox_ppd" & simple48$covariates == "none"],
          simple48$slope_wt[simple48$axis == "ox_rel" & simple48$covariates == "epi + imm"]),
  "  (ox_rel), unscaled response, adjusted. The factor of several hundred between the first",
  "  and the last was units and predictor, not adjustment.",
  "",
  sprintf(paste("WHAT CARRIES THE WILD-TYPE MOVE: of the %+.3f change on mitoPPS, %+.3f runs through",
                "`imm` and %+.3f through `epi` (`coupling_decomp`). Within the wild-type animals the",
                "ratio correlates %+.2f with `imm` and the axis %+.2f; `epi` and `imm` correlate",
                "%+.2f (wild type) and %+.2f (Myc+). Each within-genotype adjusted fit has 8 residual",
                "df, the covariates are RNA surrogates from the same count matrix as the response,",
                "and NO POSITIVE CONTROL was carried through the adjustment -- the rule that now",
                "requires one (handoff section 5 item 7) postdates script 48 by five days."),
          dwt$difference, dwt$via_imm, dwt$via_epi,
          ccr("neg", "ratio ~ imm"), ccr("neg", "ox_ppd ~ imm"),
          ccr("neg", "epi ~ imm"), ccr("pos", "epi ~ imm")),
  "",
  "THE ADJUSTED p ON DISK IS IHW, not plain Benjamini-Hochberg: 'Weighted BH adjusted",
  "p-values' is the column's own description. Four panel legends and two documents had",
  "called it BH. The rules here use RAW p throughout, so no verdict depends on it.",
  "",
  "SCOPE. batch = timepoint, so both axes of the plane are DESCRIBED and the distance from",
  "the diagonal is the batch-clean quantity. n = 6 per cell. PART D is 24 correlated mice",
  "and is ranking, not confirmatory inference. The arms' nulls answer 'more than",
  "comparably expressed genes', never 'more than chance'.")

out <- list(
  rules = rules, licence = licence,
  identity = identity_tbl,
  four_genes = four_genes, check1 = check1, check2 = check2,
  ratios = pr, ratio_range = ratio_range, ratio_line = ratio_line,
  ratio_resid = ratio_resid, ratio_band = ratio_band,
  retention_provenance = retention_provenance,
  arms_content = arms_content, arm_null = arm_null, arms_mitopps = arms_mitopps,
  arm_diagonal = arm_diagonal,
  plane_limits = plane_limits, plane_limits_mitopps = plane_limits_mitopps,
  coupling_fits = coupling_fits, coupling_interaction = coupling_interaction,
  coupling_decomp = coupling_decomp, coupling_cor = coupling_cor,
  coupling_perm = coupling_perm, coupling_panel = coupling_panel,
  coupling_lines = coupling_lines, coupling_timepoint = coupling_timepoint,
  old_fig2I = old_fig2I, tradeoff_43 = tradeoff_43,
  params = list(seed = 54, NSET = NSET, NBIN = NBIN, NPERM = NPERM,
                diagonal_threshold = DIAG, alpha = ALPHA, iqr_min = IQR_MIN,
                pct_tolerance = PCT_TOL, genes = GENES, target_ratio = TARGET),
  analysis_date = Sys.Date(), notes = notes)

saveRDS(out, here::here("results", "two_timeline_verification.rds"))
message("54: wrote results/two_timeline_verification.rds")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  res <- readRDS(here::here("results", "two_timeline_verification.rds"))

  ## --- the two checks that decide sentences --------------------------------
  res$four_genes |>
    dplyr::select(gene, wt_lfc, wt_p, myc_lfc, myc_p, int_lfc, int_se, int_p,
                  int_padj_ihw, verdict, licence) |>
    print(width = Inf)
  res$check1$verdict
  res$check2$verdict

  ## --- the ratio set: the dynamic range, then the internal line ------------
  res$ratio_range |> print(width = Inf)
  res$ratio_line  |> print(width = Inf)
  res$ratio_resid |> dplyr::arrange(residual) |> print(n = 9)

  ## --- the arms, both rulers, and the declared 0.20 ------------------------
  res$arm_diagonal |> print(n = 20)
  ## where each arm sits against comparably expressed genes, on all three
  ## contrasts. pct_wt reproduces script 43; pct_myc and pct_int are new.
  res$arm_null |>
    dplyr::select(arm, n_matched, obs_wt, pct_wt, obs_myc, pct_myc, obs_int, pct_int) |>
    print(n = 20)

  ## --- the coupling, one scale ---------------------------------------------
  res$coupling_fits |> dplyr::filter(axis == "ox_ppd") |> print(n = 20)
  res$coupling_interaction |> print(n = 10)
  ## the control axis: the two redox slopes do not differ, on either fit
  res$coupling_fits |> dplyr::filter(axis == "redox_ppd") |> print(n = 20)
  ## which covariate carries the wild-type move
  res$coupling_decomp |> print(n = 12)
  ## the open item, recorded and not pursued
  res$coupling_timepoint |> print(n = 8)

  ## --- the old numbers, for reading the record against this table ----------
  res$old_fig2I |> print(n = 8)
}
