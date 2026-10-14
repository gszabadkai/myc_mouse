# =============================================================================
# 56_window_and_chain_statistics.R
# -----------------------------------------------------------------------------
# TWO STATISTICS THE 2026-10-14 PROPOSAL PANELS DRAW THAT NO SAVED OBJECT HELD.
# The Nat Metab figure layer computes no test (figures/natmetab/_style.R); the
# second-round proposals computed these two inside their panel scripts, marked
# provisional. This script makes them a record the panels read and assert.
#
# PART A -- Fig. 1f with ED Fig. 2f folded in (author, 2026-10-14): THE
#   EXPRESSION-MATCHED NULL FOR THE TWO ROWS SCRIPT 54 HAS NO ARM FOR. Fig. 1f
#   draws the carnitine shuttle and fatty acid oxidation for their priority gain.
#   The merged panel marks every expression point against 2,000 random gene sets
#   of the same size from the same expression ventiles. Script 54 drew those sets
#   for script 43's arms only, and these two MitoCarta pathways are not among them.
#   THE RECIPE is script 54 PART C, line for line: 20 baseMean bins over the
#   expressed genes, 2,000 draws with replacement within each member's bin, and
#   the WT window, the MYC+ window and their difference read from the same draws.
#   It is proved before it is used:
#     CONTROL A1  three of script 54's own arms, re-derived here, land within
#                 script 54's Monte Carlo tolerance (5 percentile points), with
#                 the same membership and the same observed values.
#     CONTROL A2  the two new sets ARE the pathways Fig. 1f draws: their set
#                 averages equal script 40's ruler.
#
# PART B -- ED Fig. 3a, the per-animal form of Fig. 2H+I (alt) (author,
#   2026-10-14): THE TESTS FOR THE TWO COMPOSITES. ED 3a draws, animal by animal,
#   OXPHOS priority (script 42's axis score, which is script 08's OXPHOS-subunit
#   mitoPPS), the PUMA:BCL-XL ratio (script 42's log2(PUMA) - log2(BCL-XL)), Bax
#   and Bmf. Each facet brackets 12W vs 6W in each genotype and the interaction.
#   The two transcripts have DESeq2's record. The two composites have no saved
#   within-genotype contrast, and mitoPPS has no saved interaction. They get
#   SCRIPT 42's MODEL, lm(value ~ age * genotype) over the 24 animals, and the
#   within-genotype contrasts are linear combinations of its coefficients, so the
#   three numbers in a facet come from one fit.
#     CONTROL B1  for the ratio the fit IS script 42's: its saved 6W effect, 12W
#                 effect and interaction, to 1e-8.
#     CONTROL B2  the OXPHOS score IS script 08's: its window differences equal
#                 script 08's, to 1e-6. Script 08's Welch P is stored beside the
#                 model's because Fig. 1f's legend quotes it. The two differ
#                 (pooled against separate variances), so each legend names the
#                 test it draws.
#   The DESeq2 record for Bax and Bmf is copied into the same table (raw P), so the
#   panel reads one object; the copy is asserted against interaction_results.rds.
#
# BATCH = TIMEPOINT: every within-genotype contrast here coincides with the
# extraction batch. The interaction does not, and the matched sets share the
# batch with the arm they are compared with. The legends say so.
#
# SCOPE: n = 6 per group; RAW P throughout (script 54's rule). The matched-set
# percentiles are Monte Carlo (2,000 draws), so another seed moves them by about a
# point. Nothing here is a new claim: each part puts a number on a comparison the
# panels already draw.
#
# Reads:  results/interaction_results.rds          (script 03: raw LFCs, baseMean, the record)
#         results/dds_int_run.rds                  (script 03: normalised counts)
#         results/two_timeline_verification.rds    (script 54: $arm_null, $params)
#         results/substrate_specificity_tradeoff.rds (script 43: $defs$arms)
#         results/background_vs_myc.rds            (script 40: $ruler)
#         results/priming_arm_teb.rds              (script 42: $axis_scores, $priming)
#         results/mitopps_scores.rds               (script 08: $mitopps_pairwise)
#         data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt
# Writes: results/window_and_chain_statistics.rds
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

set.seed(56)

# the two Fig. 1f rows (script 40's pathway names) and their library sets
NEW_ARMS     <- c("Carnitine shuttle"    = "MITOCARTA_CARNITINE_SHUTTLE",
                  "Fatty acid oxidation" = "MITOCARTA_FATTY_ACID_OXIDATION")
# script 54 arms re-derived as the control: mid-range percentiles, where a recipe
# error would show (an arm at the 0th percentile cannot fail the comparison)
CONTROL_ARMS <- c("OXPHOS assembly", "TCA cycle", "nucleotide metabolism")
GENES        <- c("Bbc3", "Bcl2l1", "Bax", "Bmf")

# =============================================================================
# PART 0: LOAD
# =============================================================================
message("56 PART 0: load")
ir  <- readRDS(here::here("results", "interaction_results.rds"))
dds <- readRDS(here::here("results", "dds_int_run.rds"))
tv  <- readRDS(here::here("results", "two_timeline_verification.rds"))
ss  <- readRDS(here::here("results", "substrate_specificity_tradeoff.rds"))
bv  <- readRDS(here::here("results", "background_vs_myc.rds"))
pa  <- readRDS(here::here("results", "priming_arm_teb.rds"))
mp  <- readRDS(here::here("results", "mitopps_scores.rds"))
gmt <- fgsea::gmtPathways(here::here("data", "genesets_from_library",
                                     "mammary_mito_myc_metab_v1_mouse.gmt"))

NSET    <- tv$params$NSET              # script 54's (and script 43's) draw count
NBIN    <- tv$params$NBIN              # script 54's baseMean bins
PCT_TOL <- tv$params$pct_tolerance     # script 54's Monte Carlo tolerance
stopifnot(NSET == 2000L, NBIN == 20L, PCT_TOL == 5,
          all(NEW_ARMS %in% names(gmt)))

# =============================================================================
# PART A: THE MATCHED NULL FOR FIG. 1F'S TWO EXTRA ROWS
# =============================================================================
message("56 PART A: the expression-matched null, script 54's recipe")

CONTRASTS <- c("myc_6W_raw", "myc_12W_raw", "timepoint_neg_raw", "timepoint_pos_raw",
               "interaction_raw")
D <- lapply(ir[CONTRASTS], as.data.frame)
universe_all <- rownames(D$interaction_raw)
stopifnot(all(vapply(D, function(x) identical(rownames(x), universe_all), logical(1))))
V  <- function(k, col = "log2FoldChange") stats::setNames(D[[k]][[col]], universe_all)
tn <- V("timepoint_neg_raw"); tp <- V("timepoint_pos_raw"); it <- V("interaction_raw")
bm <- V("myc_6W_raw", "baseMean")

ens_set  <- function(syms) { e <- recon_to_ensembl(syms, universe_all); e[!is.na(e)] }
set_mean <- function(v, e) if (length(e)) mean(v[e], na.rm = TRUE) else NA_real_

# --- script 54 PART C, verbatim -------------------------------------------------
expressed <- names(bm)[is.finite(bm) & bm > 0 & is.finite(tn[names(bm)])]
bin_of <- cut(rank(bm[expressed], ties.method = "first"), breaks = NBIN, labels = FALSE)
by_bin <- split(expressed, bin_of)
draw_matched <- function(e) {
  b <- bin_of[match(e, expressed)]
  b <- b[!is.na(b)]
  unlist(lapply(split(b, b), function(k)
    sample(by_bin[[as.character(k[1])]], length(k), replace = TRUE)), use.names = FALSE)
}
null_row <- function(a, e) {
  e   <- e[e %in% expressed]
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
}

# --- CONTROL A1: the recipe reproduces script 54 ----------------------------------
arms_def <- as.data.frame(ss$defs$arms)
an54     <- as.data.frame(tv$arm_null)
stopifnot(all(CONTROL_ARMS %in% arms_def$arm), all(CONTROL_ARMS %in% an54$arm))
ctrl <- dplyr::bind_rows(lapply(CONTROL_ARMS, function(a)
  null_row(a, ens_set(gmt[[arms_def$set[arms_def$arm == a]]]))))
k54  <- match(ctrl$arm, an54$arm)
window_null_control <- tibble::tibble(
  arm = ctrl$arm, n_matched = ctrl$n_matched, n_matched_54 = an54$n_matched[k54],
  obs_wt = ctrl$obs_wt, obs_wt_54 = an54$obs_wt[k54],
  pct_wt = ctrl$pct_wt, pct_wt_54 = an54$pct_wt[k54],
  pct_myc = ctrl$pct_myc, pct_myc_54 = an54$pct_myc[k54],
  pct_int = ctrl$pct_int, pct_int_54 = an54$pct_int[k54])
stopifnot(all(ctrl$n_matched == an54$n_matched[k54]),
          max(abs(ctrl$obs_wt  - an54$obs_wt[k54]))  < 1e-12,
          max(abs(ctrl$obs_myc - an54$obs_myc[k54])) < 1e-12,
          max(abs(ctrl$pct_wt  - an54$pct_wt[k54]))  <= PCT_TOL,
          max(abs(ctrl$pct_myc - an54$pct_myc[k54])) <= PCT_TOL,
          max(abs(ctrl$pct_int - an54$pct_int[k54])) <= PCT_TOL)
message(sprintf("56 PART A: script 54's arms reproduce to %.2f percentile points",
                max(abs(c(ctrl$pct_wt - an54$pct_wt[k54], ctrl$pct_myc - an54$pct_myc[k54],
                          ctrl$pct_int - an54$pct_int[k54])))))

# --- the two new arms -----------------------------------------------------------------
window_null <- dplyr::bind_rows(lapply(names(NEW_ARMS), function(a)
  dplyr::mutate(null_row(a, ens_set(gmt[[NEW_ARMS[[a]]]])), set = NEW_ARMS[[a]], .after = arm)))
window_null <- dplyr::mutate(window_null,
  beyond_wt  = pct_wt  <= 2.5 | pct_wt  >= 97.5,     # outside the central 95%, two-sided
  beyond_myc = pct_myc <= 2.5 | pct_myc >= 97.5)

# --- CONTROL A2: they are the pathways Fig. 1f draws -----------------------------------
rl <- as.data.frame(bv$ruler)
kr <- match(window_null$arm, rl$pathway)
stopifnot(!anyNA(kr),
          max(abs(window_null$obs_wt  - rl$c_tn[kr])) < 1e-8,
          max(abs(window_null$obs_myc - rl$c_tp[kr])) < 1e-8)

# =============================================================================
# PART B: THE PER-ANIMAL TESTS FOR ED FIG. 3A
# =============================================================================
message("56 PART B: the per-animal tests")

cd  <- as.data.frame(SummarizedExperiment::colData(dds))
smp <- colnames(dds)
age <- stats::relevel(factor(cd$timepoint), "6W")
gen <- stats::relevel(factor(cd$myc_status), "neg")
stopifnot(all(table(age, gen) == 6L))

L   <- log2(DESeq2::counts(dds, normalized = TRUE) + 1)          # script 42's scale
ens <- vapply(GENES, function(g) recon_to_ensembl(g, rownames(dds)), character(1))
stopifnot(!anyNA(ens), all(ens %in% universe_all))
ax  <- as.data.frame(pa$axis_scores)
stopifnot(setequal(ax$sample, smp))

values <- tibble::tibble(
  sample = smp, timepoint = as.character(age), myc_status = as.character(gen),
  group = paste(as.character(age), as.character(gen), sep = "_"),
  oxphos_mitopps = ax$oxphos_ppd[match(smp, ax$sample)],
  puma_bclxl     = as.numeric(L[ens[["Bbc3"]], ] - L[ens[["Bcl2l1"]], ]),
  Bax            = as.numeric(L[ens[["Bax"]], ]),
  Bmf            = as.numeric(L[ens[["Bmf"]], ]))
stopifnot(!anyNA(values$oxphos_mitopps))

# the five contrasts, as weights on lm(value ~ age * genotype)'s coefficients
W <- rbind(wt_12v6       = c(0, 1, 0, 0),
           myc_12v6      = c(0, 1, 0, 1),
           interaction   = c(0, 0, 0, 1),
           myc_vs_wt_6W  = c(0, 0, 1, 0),
           myc_vs_wt_12W = c(0, 0, 1, 1))
lm_tests <- function(v, quantity, instrument) {
  f  <- stats::lm(v ~ age * gen)
  cf <- stats::coef(f)
  stopifnot(identical(names(cf), c("(Intercept)", "age12W", "genpos", "age12W:genpos")))
  est <- drop(W %*% cf); se <- sqrt(diag(W %*% stats::vcov(f) %*% t(W)))
  tibble::tibble(quantity = quantity, instrument = instrument, contrast = rownames(W),
                 estimate = unname(est), se = unname(se), df = f$df.residual,
                 p = unname(2 * stats::pt(-abs(est / se), f$df.residual)))
}
deseq_tests <- function(g) {
  slot <- c(wt_12v6 = "timepoint_neg_raw", myc_12v6 = "timepoint_pos_raw",
            interaction = "interaction_raw", myc_vs_wt_6W = "myc_6W_raw", myc_vs_wt_12W = "myc_12W_raw")
  r <- lapply(slot, function(s) D[[s]][ens[[g]], c("log2FoldChange", "lfcSE", "pvalue")])
  tibble::tibble(quantity = g, instrument = "DESeq2 Wald (script 03), raw P", contrast = names(slot),
                 estimate = vapply(r, `[[`, numeric(1), "log2FoldChange"),
                 se = vapply(r, `[[`, numeric(1), "lfcSE"), df = NA_real_,
                 p = vapply(r, `[[`, numeric(1), "pvalue"))
}
tests <- dplyr::bind_rows(
  lm_tests(values$oxphos_mitopps, "oxphos_mitopps", "lm(value ~ age * genotype), 24 animals"),
  lm_tests(values$puma_bclxl,     "puma_bclxl",     "lm(value ~ age * genotype), 24 animals (script 42's model)"),
  deseq_tests("Bax"), deseq_tests("Bmf"))
tget <- function(q, ct, col = "p") tests[[col]][tests$quantity == q & tests$contrast == ct]

# --- CONTROL B1: the ratio's fit is script 42's ------------------------------------------
pr42 <- as.data.frame(pa$priming)
pr42 <- pr42[pr42$pair == "Bbc3:Bcl2l1", ]
stopifnot(nrow(pr42) == 1L,
          abs(tget("puma_bclxl", "myc_vs_wt_6W", "estimate")  - pr42$d6)    < 1e-8,
          abs(tget("puma_bclxl", "myc_vs_wt_6W")              - pr42$p6)    < 1e-8,
          abs(tget("puma_bclxl", "myc_vs_wt_12W", "estimate") - pr42$d12)   < 1e-8,
          abs(tget("puma_bclxl", "interaction", "estimate")   - pr42$int)   < 1e-8,
          abs(tget("puma_bclxl", "interaction") - pr42$int_p) < 1e-8)

# --- CONTROL B2: the OXPHOS score is script 08's, and its Welch test kept beside ---------
pw   <- as.data.frame(mp$mitopps_pairwise)
gm   <- function(v, g) mean(v[values$group == g])
w08  <- function(ct) pw[pw$contrast == ct & pw$pathway == "OXPHOS subunits", ]
stopifnot(nrow(w08("Temporal_Myc-")) == 1L, nrow(w08("Temporal_Myc+")) == 1L,
          abs((gm(values$oxphos_mitopps, "12W_neg") - gm(values$oxphos_mitopps, "6W_neg")) -
                w08("Temporal_Myc-")$diff) < 1e-6,
          abs((gm(values$oxphos_mitopps, "12W_pos") - gm(values$oxphos_mitopps, "6W_pos")) -
                w08("Temporal_Myc+")$diff) < 1e-6,
          abs(tget("oxphos_mitopps", "wt_12v6", "estimate")  - w08("Temporal_Myc-")$diff) < 1e-6,
          abs(tget("oxphos_mitopps", "myc_12v6", "estimate") - w08("Temporal_Myc+")$diff) < 1e-6)
welch_08 <- tibble::tibble(
  contrast = c("wt_12v6", "myc_12v6"), script08_contrast = c("Temporal_Myc-", "Temporal_Myc+"),
  diff = c(w08("Temporal_Myc-")$diff, w08("Temporal_Myc+")$diff),
  p_welch = c(w08("Temporal_Myc-")$p_value, w08("Temporal_Myc+")$p_value),
  padj_welch_bh = c(w08("Temporal_Myc-")$padj, w08("Temporal_Myc+")$padj),
  p_lm = c(tget("oxphos_mitopps", "wt_12v6"), tget("oxphos_mitopps", "myc_12v6")))

# --- the copy of the DESeq2 record is the record -----------------------------------------
for (g in c("Bax", "Bmf")) {
  stopifnot(abs(tget(g, "interaction") - as.data.frame(ir$interaction_raw)[ens[[g]], "pvalue"]) < 1e-12,
            abs(tget(g, "myc_12v6", "estimate") -
                  as.data.frame(ir$timepoint_pos_raw)[ens[[g]], "log2FoldChange"]) < 1e-12)
}

# --- the reading, asserted (what ED 3a's legend and the text may say) ---------------------
stopifnot(
  # PUMA:BCL-XL falls across the window only under MYC, and the two windows differ
  tget("puma_bclxl", "myc_12v6", "estimate") < 0, tget("puma_bclxl", "myc_12v6") < 0.05,
  tget("puma_bclxl", "wt_12v6") > 0.05, tget("puma_bclxl", "interaction") < 0.05,
  # Bmf rises in both genotypes and MYC does not change that: the developmental control
  tget("Bmf", "wt_12v6", "estimate") > 0, tget("Bmf", "wt_12v6") < 0.05,
  tget("Bmf", "myc_12v6", "estimate") > 0, tget("Bmf", "myc_12v6") < 0.05,
  tget("Bmf", "interaction") > 0.05,
  # Bax falls under MYC, and the windows do not differ by genotype: it fades
  tget("Bax", "myc_12v6", "estimate") < 0, tget("Bax", "myc_12v6") < 0.05,
  tget("Bax", "interaction") > 0.05,
  # OXPHOS priority falls in the MYC+ window; in WT it falls by as much but does not
  # clear 0.05; the windows do not differ. "Falls on both timelines" is NOT licensed.
  tget("oxphos_mitopps", "myc_12v6", "estimate") < 0, tget("oxphos_mitopps", "myc_12v6") < 0.05,
  tget("oxphos_mitopps", "wt_12v6", "estimate") < 0, tget("oxphos_mitopps", "wt_12v6") > 0.05,
  tget("oxphos_mitopps", "interaction") > 0.05)

# =============================================================================
# PART C: NOTES, AND SAVE
# =============================================================================
fmt <- function(q, ct) sprintf("%+.3f (P %.2g)", tget(q, ct, "estimate"), tget(q, ct))
notes <- c(
  sprintf("PART A: script 54's recipe reproduces its arms to %.2f percentile points (tolerance %g).",
          max(abs(c(window_null_control$pct_wt - window_null_control$pct_wt_54,
                    window_null_control$pct_myc - window_null_control$pct_myc_54,
                    window_null_control$pct_int - window_null_control$pct_int_54))), PCT_TOL),
  sprintf("  %-22s %3d genes: WT %+.3f (pct %.1f), MYC+ %+.3f (pct %.1f)", window_null$arm,
          window_null$n_matched, window_null$obs_wt, window_null$pct_wt,
          window_null$obs_myc, window_null$pct_myc),
  "PART B: 12W vs 6W in WT | in MYC+ | interaction (raw P)",
  sprintf("  %-15s %s | %s | %s", c("OXPHOS mitoPPS", "PUMA:BCL-XL", "Bax", "Bmf"),
          vapply(c("oxphos_mitopps", "puma_bclxl", "Bax", "Bmf"), fmt, character(1), ct = "wt_12v6"),
          vapply(c("oxphos_mitopps", "puma_bclxl", "Bax", "Bmf"), fmt, character(1), ct = "myc_12v6"),
          vapply(c("oxphos_mitopps", "puma_bclxl", "Bax", "Bmf"), fmt, character(1), ct = "interaction")),
  sprintf("  OXPHOS mitoPPS by script 08's Welch test (Fig. 1f's legend): WT P %.2g, MYC+ P %.2g (BH %.2g, %.2g)",
          welch_08$p_welch[1], welch_08$p_welch[2], welch_08$padj_welch_bh[1], welch_08$padj_welch_bh[2]),
  "  The composites' estimates are in their own units (mitoPPS; log2 ratio); the genes' in log2 fold change.",
  "SCOPE: n = 6 per group, raw P. Within-genotype contrasts coincide with the extraction batch;",
  "  the interaction does not. OXPHOS priority does not clear 0.05 in the WT window.")
cat(paste(notes, collapse = "\n"), "\n")

out <- list(
  window_null = window_null, window_null_control = window_null_control,
  values = values, tests = tests, welch_08 = welch_08,
  defs = list(new_arms = NEW_ARMS, control_arms = CONTROL_ARMS, genes = GENES,
              NSET = NSET, NBIN = NBIN, pct_tolerance = PCT_TOL, seed = 56L,
              scale = "log2(DESeq2 normalised counts + 1), script 42's",
              contrasts = W),
  analysis_date = Sys.Date(), notes = notes)

saveRDS(out, here::here("results", "window_and_chain_statistics.rds"))
message("56: wrote results/window_and_chain_statistics.rds")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  res <- readRDS(here::here("results", "window_and_chain_statistics.rds"))

  ## PART A: the two new rows, and the control against script 54
  res$window_null |> dplyr::select(arm, n_matched, obs_wt, pct_wt, obs_myc, pct_myc) |> print()
  res$window_null_control |> print()

  ## PART B: every test ED 3a draws, and the Welch comparison for mitoPPS
  res$tests |> dplyr::filter(contrast %in% c("wt_12v6", "myc_12v6", "interaction")) |>
    dplyr::select(quantity, contrast, estimate, p) |> print(n = 12)
  res$welch_08 |> print()

  ## the per-animal values ED 3a draws (it z-scores each column over the 24)
  res$values |> dplyr::group_by(group) |>
    dplyr::summarise(dplyr::across(c(oxphos_mitopps, puma_bclxl, Bax, Bmf), mean)) |> print()
}
