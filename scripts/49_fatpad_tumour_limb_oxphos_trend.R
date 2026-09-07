# =============================================================================
# 49 -- THE 12W-TO-TUMOUR LIMB: DOES RELATIVE OXPHOS SHARE RISE AGAIN?
# -----------------------------------------------------------------------------
# Cohort: the fat-pad progression timeline (`data/fatpad_timeline/`), NOT the
# purified-MEC dataset that scripts 00-48 analyse. Read that folder's README before
# this file: whole mammary fat pad, adipose-majority throughout, DESeq-normalised
# upstream, no 12W negative control.
#
# THE QUESTION. The MEC cohort establishes that between 6W and 12W the gland
# withdraws nuclear-encoded OXPHOS subunits RELATIVELY, while the compartment gains
# content and the mitoribosome rises -- a reprioritisation WITHIN the compartment,
# not a change in its size. What happens AFTER 12W is unknown, and this is the only
# dataset that can address it. The human arm says MYC-high tumours are
# respiration-high, so the trend must reverse somewhere. The claim tested here is
# that relative OXPHOS share falls to 12W and RISES again as tumours establish.
#
# THIS IS A BETWEEN-GROUP ORDERED TREND, and it is a different estimand from the
# pre-check's (`docs/2026-09-07_fatpad_timeline_oxphos_precheck.md`). That check
# asked whether `ox_rel` could support a WITHIN-TUMOUR correlation at n = 15 and
# correctly returned FAIL. It is not reopened here, not rescued, and no ruler is
# swapped for one that happened to pass. A within-tumour correlation in this cohort
# stays closed.
#
# WHY THE TEST IS ADMISSIBLE DESPITE THE CONFOUND -- fixed in advance, before any
# number below was seen. The pre-check measured adipose loading at +0.538 on the
# nuclear OXPHOS subunits and +0.326 on the rest of MitoCarta, so adipose pushes
# `ox_rel` UP; and adipose FALLS across this series (Adipoq 59,190 -> 31,595 ->
# 26,415 on group medians). COMPOSITION ALONE THEREFORE PREDICTS `ox_rel` FALLING
# from 6W to large tumour. The prediction is that it RISES after 12W. A positive
# result is credible because it runs against the confound; a NULL IS UNINTERPRETABLE
# and must not be reported as evidence of no change. PART D turns this from an
# argument into a measurement (T3) and extends it to the primary endpoint, whose
# confound direction the pre-check never established.
#
# WHY THE LIMB IS RESTRICTED TO THREE GROUPS -- also fixed in advance. Composition
# moves very unevenly: Adipoq falls 1.87x over 6W_POS -> 12W but only 1.20x over
# 12W -> LARGE. The developmental limb is confound-aligned and uninterpretable here,
# and is already clean in the MEC cohort, so nothing is lost by dropping it. The
# 12W-to-tumour limb is where composition is most stable and is the part available
# nowhere else. The 6W groups are scored and reported for completeness and CARRY NO
# CLAIM.
#
# DECISION RULES, FIXED BEFORE THE FIT:
#   T1 PRIMARY   `ox_nuc_mtrib` increases monotonically across the three groups,
#                Jonckheere-Terpstra one-sided p < 0.05.  PASS / FAIL.
#   T2 CONCORD   `ox_rel` agrees in DIRECTION with T1. Concordance, not a test.
#   T3 HARD REQ  the adipose markers must move in the direction that would OPPOSE
#                the observed trend. If they move WITH it, T1 is VOID whatever its
#                p-value.
#   T4 REPORTED  `mito_total` and `ox_lvl`: does the compartment as a whole move, or
#                only its internal proportions? A share change on a flat total is
#                the reprioritisation claim; both moving together is a content claim
#                and a weaker, different finding.
#   T5 OBSERVED  `Bcl2l1` trend direction against the `ox_nuc_mtrib` trend. No
#                threshold: n = 20, and this is an observation.
#
# ADDENDUM (PART G) -- THE ADJUSTED TREND, WHICH CLOSES THIS DATASET. T3 voided T1
# with a SINGLE-marker calculation: Adipoq's predicted push was -0.249 x 0.576 =
# -0.144 against an observed -0.147, so the residual is about zero. But the residual's
# SIGN depends on WHICH adipose marker carries the adjustment -- Cidec loads +0.776 on
# the same endpoint against Adipoq's +0.576, and on a similar trend its push would flip
# the residual positive. Four single-marker adjustments give four incompatible answers.
# PART G replaces them with one number and an interval.
#   NO PASS/FAIL. It returns an interval, not a verdict. The limb was declared
#   uninterpretable and stays so unless the joint interval excludes zero, which is NOT
#   expected: four collinear markers explain a large share of the endpoint at n = 20.
#   A WIDE INTERVAL SPANNING ZERO IS UNINFORMATIVE and must not be presented otherwise.
#   If the joint interval excludes zero POSITIVELY that is a suggestion requiring the
#   orthotopic series, never support for the reversal. If the four single-marker
#   adjustments DISAGREE IN SIGN, that disagreement IS the primary reported result.
#   Mki67 is the NEGATIVE CONTROL and is read FIRST: proliferation genuinely rises on
#   this limb, so if adjustment kills that too, the adjustment is removing signal rather
#   than confound and nothing else in PART G may be read.
#
# WHAT THIS CANNOT ESTABLISH: nothing about 6W->12W (that belongs to the MEC
# cohort); no within-tumour correlation (closed by the pre-check); no causation
# (nothing is manipulated here); and NO CROSS-COHORT COMPARISON -- the branch's
# standing rule holds, these scores are not comparable with the MEC cohort's or the
# orthotopic series'. Directions and orderings travel; values never do.
#
# Reads : data/fatpad_timeline/Normalised_DESseq_counts.csv (DESeq-normalised, LINEAR)
#         data/Mouse.MitoCarta3.0.xls  Sheet 4  (script 08's authoritative splits)
#         data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt (ox_rel only)
#         functions/reconcile_gene_symbols.R (MANDATORY -- vintage-aware membership)
# Writes: results/fatpad_tumour_limb_trend.rds
#         outputs/fatpad_limb/49_per_sample_trends.pdf
#
# RUNTIME: about one minute. The permutation nulls are 20,000 label shuffles of 20
# samples for each endpoint.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

set.seed(1)
NPERM  <- 20000L    # label shuffles for the J-T null
NBOOT  <- 5000L     # bootstrap resamples for the tau interval
MINBM  <- 10        # smallest mean normalised count worth carrying
TREND  <- c("12WK_POS", "SMALL_TUMOUR", "LARGE_TUMOUR")   # the limb, in order
ALLGRP <- c("6WK_NEG", "6WK_POS", TREND)

# =============================================================================
# PART 0: LOAD
# =============================================================================
message("49 PART 0: load")

raw <- utils::read.csv(
  here::here("data", "fatpad_timeline", "Normalised_DESseq_counts.csv"),
  row.names = 1, check.names = FALSE)
NC  <- as.matrix(raw)
grp <- factor(sub("_R[0-9]+$", "", colnames(NC)), levels = ALLGRP)
stopifnot(ncol(NC) == 30L, !anyNA(grp), all(grepl("^ENSMUSG", rownames(NC))))
stopifnot(identical(as.integer(table(grp)), c(5L, 5L, 5L, 6L, 9L)))

keep <- rowMeans(NC) >= MINBM & matrixStats::rowSds(NC) > 0
L    <- log2(NC[keep, , drop = FALSE] + 1)      # z-score constructions live here
message(sprintf("49: %d of %d genes at mean normalised count >= %d",
                nrow(L), nrow(NC), MINBM))

universe_all <- rownames(L)
zrow    <- function(m) t(scale(t(m)))
comp_e  <- function(e) { e <- e[!is.na(e) & e %in% rownames(L)]
                         stopifnot(length(e) >= 3); colMeans(zrow(L[e, , drop = FALSE])) }
ens_set <- function(syms) { e <- recon_to_ensembl(syms, universe_all); e[!is.na(e)] }
one_ens <- function(sym) { e <- ens_set(sym); stopifnot(length(e) >= 1); e[1] }

# =============================================================================
# PART A: THE SETS -- SCRIPT 08's AUTHORITATIVE SPLITS, REBUILT FROM SHEET 4
# -----------------------------------------------------------------------------
# Not the library GMT for the primary endpoint: the session specifies script 08's
# splits, and 08 builds them from MitoCarta Sheet 4 with `cSplit`, then STRIPS every
# mt-* gene out of every pathway into a synthetic mtDNA pathway. Rebuilt here the
# same way rather than read from `mitopps_scores.rds`, whose gene lists were already
# filtered to the MEC matrix -- a cohort-dependent filter that must not travel.
# `ox_rel` keeps script 48's GMT recipe exactly, so T2 is a true concordance check
# between two independently-constructed rulers.
# =============================================================================
message("49 PART A: sets")

s4  <- readxl::read_xls(here::here("data", "Mouse.MitoCarta3.0.xls"), sheet = 4)
s4  <- na.omit(dplyr::select(s4, MitoPathway, Genes))
# cSplit emits one "'as.is' should be specified by the caller" per column from its
# internal type.convert -- 464 of them here, which buries any warning that matters.
# Suppressed at this call only; nothing numeric changes.
g2p <- suppressWarnings(splitstackshape::cSplit(s4, "Genes", ","))
g2p <- tibble::column_to_rownames(g2p, "MitoPathway")
g2p <- as.data.frame(t(g2p))
g2p <- tidyr::pivot_longer(g2p, cols = colnames(g2p),
                           names_to = "Pathway", values_to = "Gene")
g2p <- dplyr::mutate(na.omit(g2p), Gene = as.character(Gene))

is_mt   <- grepl("^[Mm][Tt]-", g2p$Gene)          # 08's mtDNA rule, verbatim
mt_syms <- unique(g2p$Gene[is_mt])
g2p_nuc <- g2p[!is_mt, ]                          # every pathway, mtDNA stripped
path_of <- function(p) { s <- unique(g2p_nuc$Gene[g2p_nuc$Pathway == p])
                         stopifnot(length(s) > 0); s }

ox_sub  <- ens_set(path_of("OXPHOS subunits"))          # nuclear by construction
mtrib   <- ens_set(path_of("Mitochondrial ribosome"))   # nuclear by construction
mito_08 <- ens_set(unique(g2p_nuc$Gene))                # whole nuclear compartment
mt_ens  <- ens_set(mt_syms)

# the separation the primary endpoint depends on, asserted rather than assumed
stopifnot(length(intersect(ox_sub, mt_ens)) == 0L,
          length(intersect(mtrib,  mt_ens)) == 0L,
          length(intersect(ox_sub, mtrib))  == 0L,
          length(mt_ens) == 13L)

# ox_rel keeps script 48's GMT construction, mtDNA genes in the DENOMINATOR
gmt       <- fgsea::gmtPathways(here::here("data", "genesets_from_library",
                                           "mammary_mito_myc_metab_v1_mouse.gmt"))
set_e     <- function(s) { stopifnot(!is.null(gmt[[s]])); e <- ens_set(gmt[[s]]); e[e %in% rownames(L)] }
mito_all  <- unique(unlist(lapply(grep("^MITOCARTA_", names(gmt), value = TRUE), set_e)))
ox_gmt    <- set_e("MITOCARTA_OXPHOS_SUBUNITS")
rest_mito <- setdiff(mito_all, ox_gmt)
stopifnot(all(mt_ens %in% rest_mito), length(intersect(ox_gmt, mt_ens)) == 0L)

set_sizes <- tibble::tibble(
  source = c("sheet4_08", "sheet4_08", "sheet4_08", "sheet4_08", "gmt_48", "gmt_48"),
  set    = c("OXPHOS subunits", "Mitochondrial ribosome", "nuclear MitoCarta",
             "mtDNA-encoded", "MITOCARTA_OXPHOS_SUBUNITS", "rest of MitoCarta"),
  n_mapped = c(length(ox_sub), length(mtrib), length(mito_08), length(mt_ens),
               length(ox_gmt), length(rest_mito)))

# =============================================================================
# PART B: SCORES
# -----------------------------------------------------------------------------
# All z-score constructions on log2(counts + 1). mitoPPS is NOT run: it takes the
# LINEAR counts and never shares an input object with these, and no decision rule
# depends on it.
# =============================================================================
message("49 PART B: scores")

S <- list(
  # --- primary ---
  ox_nuc_mtrib = comp_e(ox_sub) - comp_e(mtrib),      # the reprioritisation claim
  ox_rel       = comp_e(ox_gmt) - comp_e(rest_mito),  # script 48's recipe
  # --- secondary: level vs share, and the compartment as a whole ---
  ox_lvl       = comp_e(ox_sub),
  mtrib_lvl    = comp_e(mtrib),
  mito_total   = comp_e(mito_08),
  ox_mt        = comp_e(mt_ens))

MARK <- c(Adipoq = "Adipoq", Fabp4 = "Fabp4", Plin1 = "Plin1", Cidec = "Cidec",
          Krt8 = "Krt8", Epcam = "Epcam", Ptprc = "Ptprc", Mki67 = "Mki67")
Mk <- lapply(MARK, function(s) as.numeric(L[one_ens(s), ]))

APO <- c(Bcl2l1 = "Bcl2l1", Bbc3 = "Bbc3", Mcl1 = "Mcl1", Bcl2l11 = "Bcl2l11", Myc = "Myc")
Ap  <- lapply(APO, function(s) as.numeric(L[one_ens(s), ]))
Ap$puma_bclxl <- Ap$Bbc3 - Ap$Bcl2l1      # log2 ratio, both already log2
Ap$buffer     <- Ap$Mcl1 - Ap$Bcl2l1

ALL <- c(S, Mk, Ap)

# =============================================================================
# PART C: THE ORDERED-TREND TEST
# -----------------------------------------------------------------------------
# Jonckheere-Terpstra, implemented here rather than taken from a package so the
# script carries no dependency 00_setup_packages.R does not already load. J is the
# sum over ordered group pairs of Mann-Whitney counts, ties at 0.5. The null shuffles
# GROUP LABELS among the 20 samples, preserving group sizes -- the exact null for
# "no association with order". The normal approximation is reported beside it as a
# cross-check, never as the p-value.
# Effect size is Kendall's tau-b between the group index and the score, with a
# percentile bootstrap interval, because at n = 20 an estimate with an interval is
# the honest object and a p-value alone is not.
# =============================================================================
message("49 PART C: ordered trend on the 12W-to-tumour limb")

jt_stat <- function(y, g) {
  lv <- levels(g); J <- 0
  for (i in seq_len(length(lv) - 1L)) for (j in (i + 1L):length(lv)) {
    a <- y[g == lv[i]]; b <- y[g == lv[j]]
    J <- J + sum(outer(a, b, function(u, v) (v > u) + 0.5 * (v == u)))
  }
  J
}
jt_norm <- function(y, g) {           # standardised J under the null, ties ignored
  n  <- length(y); ni <- as.numeric(table(g))
  mu <- (n^2 - sum(ni^2)) / 4
  sg <- sqrt((n^2 * (2 * n + 3) - sum(ni^2 * (2 * ni + 3))) / 72)
  (jt_stat(y, g) - mu) / sg
}
jt_test <- function(y, g, direction = c("increasing", "decreasing"), nperm = NPERM) {
  direction <- match.arg(direction)
  obs  <- jt_stat(y, g)
  lab  <- as.character(g)
  null <- vapply(seq_len(nperm), function(i)
    jt_stat(y, factor(sample(lab), levels = levels(g))), 0)
  p <- if (direction == "increasing") (1 + sum(null >= obs)) / (nperm + 1) else
                                      (1 + sum(null <= obs)) / (nperm + 1)
  list(J = obs, p_one_sided = p, z_normal = jt_norm(y, g),
       null_median = stats::median(null), direction = direction)
}
tau_ci <- function(y, gi, nboot = NBOOT) {
  t0 <- stats::cor(y, gi, method = "kendall")
  bs <- vapply(seq_len(nboot), function(i) {
    k <- sample.int(length(y), replace = TRUE)
    if (length(unique(gi[k])) < 2L) return(NA_real_)
    stats::cor(y[k], gi[k], method = "kendall") }, 0)
  c(tau = t0, lo = unname(stats::quantile(bs, 0.025, na.rm = TRUE)),
              hi = unname(stats::quantile(bs, 0.975, na.rm = TRUE)))
}

lk  <- grp %in% TREND                      # the limb, n = 20
gl  <- factor(as.character(grp[lk]), levels = TREND)
gi  <- as.integer(gl)
stopifnot(sum(lk) == 20L)

trend_of <- function(v, dir = "increasing") {
  y  <- v[lk]
  tt <- jt_test(y, gl, dir)
  tc <- tau_ci(y, gi)
  tibble::tibble(J = tt$J, p_one_sided = tt$p_one_sided, z_normal = tt$z_normal,
                 kruskal_p = stats::kruskal.test(y ~ gl)$p.value,
                 tau = tc[["tau"]], tau_lo = tc[["lo"]], tau_hi = tc[["hi"]],
                 direction_tested = dir)
}

# every endpoint is tested for an INCREASING trend, so the sign of tau -- not the
# choice of tail -- carries the direction, and one table reads consistently
trend <- dplyr::bind_rows(lapply(names(ALL), function(n)
  dplyr::bind_cols(tibble::tibble(endpoint = n), trend_of(ALL[[n]], "increasing"))))

# group summaries, BOTH statistics: the pre-check found them to diverge materially
# at these group sizes, so a claim must name which one it rests on
summ <- dplyr::bind_rows(lapply(names(ALL), function(n) {
  v <- ALL[[n]]
  dplyr::bind_rows(
    tibble::tibble(endpoint = n, stat = "median",
                   !!!stats::setNames(as.list(round(tapply(v, grp, stats::median), 4)), ALLGRP)),
    tibble::tibble(endpoint = n, stat = "mean",
                   !!!stats::setNames(as.list(round(tapply(v, grp, mean), 4)), ALLGRP)))
}))

# =============================================================================
# PART D: T3 -- DOES THE CONFOUND OPPOSE THE TREND, OR RUN WITH IT?
# -----------------------------------------------------------------------------
# The admissibility argument was written for `ox_rel`, whose adipose loading the
# pre-check measured. The PRIMARY endpoint is `ox_nuc_mtrib`, and ITS confound
# direction was never established: if adipose loads more on the mitoribosome than on
# the OXPHOS subunits, then FALLING adipose pushes `ox_nuc_mtrib` UP and the confound
# runs WITH the prediction, which voids T1. That is measured here, not assumed.
# =============================================================================
message("49 PART D: confound direction (T3)")

sp <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman"))
load_dir <- dplyr::bind_rows(lapply(names(Mk), function(mn)
  tibble::tibble(marker = mn,
                 rho_ox_sub  = sp(comp_e(ox_sub)[lk], Mk[[mn]][lk]),
                 rho_mtrib   = sp(comp_e(mtrib)[lk],  Mk[[mn]][lk]),
                 rho_primary = sp(S$ox_nuc_mtrib[lk], Mk[[mn]][lk]),
                 rho_ox_rel  = sp(S$ox_rel[lk],       Mk[[mn]][lk]))))

adipose  <- c("Adipoq", "Fabp4", "Plin1", "Cidec")
adi_tau  <- trend$tau[match(adipose, trend$endpoint)]
prim_tau <- trend$tau[trend$endpoint == "ox_nuc_mtrib"]
adi_load <- load_dir$rho_primary[match(adipose, load_dir$marker)]

# the confound OPPOSES the trend when the adipose marker's own trend, multiplied by
# its loading on the endpoint, points AGAINST the endpoint's observed trend
confound_push <- adi_tau * adi_load                 # predicted push on the endpoint
t3_opposes    <- all(sign(confound_push) != sign(prim_tau) | confound_push == 0)

# -----------------------------------------------------------------------------
# LEAVE-ONE-OUT. At n = 20 a single animal can carry a trend or a loading, and this
# limb contains one that invites the worry: LARGE_TUMOUR_R8 sits at Adipoq 10.7
# where every other sample is 13.9-15.6. Both the primary trend and the T3 loading
# are recomputed dropping each sample in turn; the RANGE is the honest statement.
# -----------------------------------------------------------------------------
loo <- dplyr::bind_rows(lapply(seq_len(sum(lk)), function(i) {
  y <- S$ox_nuc_mtrib[lk][-i]; a <- Mk$Adipoq[lk][-i]; gg <- gi[-i]
  tibble::tibble(dropped = colnames(L)[lk][i],
                 tau_primary = stats::cor(y, gg, method = "kendall"),
                 rho_primary_adipoq = sp(y, a))
}))
loo_range <- tibble::tibble(
  quantity = c("tau primary", "rho(primary, Adipoq)"),
  full     = c(gv_tau <- stats::cor(S$ox_nuc_mtrib[lk], gi, method = "kendall"),
               sp(S$ox_nuc_mtrib[lk], Mk$Adipoq[lk])),
  loo_min  = c(min(loo$tau_primary), min(loo$rho_primary_adipoq)),
  loo_max  = c(max(loo$tau_primary), max(loo$rho_primary_adipoq)))

# =============================================================================
# PART E: T4/T5 AND THE PER-SAMPLE RECORD
# =============================================================================
message("49 PART E: per-sample record and figure")

per_sample <- tibble::tibble(sample = colnames(L), group = as.character(grp),
                             in_limb = lk)
for (n in names(ALL)) per_sample[[n]] <- round(ALL[[n]], 4)

dir.create(here::here("outputs", "fatpad_limb"), showWarnings = FALSE, recursive = TRUE)
plot_endpoints <- c("ox_nuc_mtrib", "ox_rel", "ox_lvl", "mtrib_lvl",
                    "mito_total", "Adipoq", "Bcl2l1", "puma_bclxl")
pd <- do.call(rbind, lapply(plot_endpoints, function(n)
  data.frame(endpoint = n, group = grp, value = ALL[[n]],
             limb = ifelse(lk, "12W -> tumour limb", "6W (no claim)"))))
pd$endpoint <- factor(pd$endpoint, levels = plot_endpoints)

p <- ggplot2::ggplot(pd, ggplot2::aes(group, value)) +
  ggplot2::stat_summary(fun = stats::median, geom = "crossbar",
                        width = 0.55, linewidth = 0.3, colour = "grey35") +
  ggplot2::geom_point(ggplot2::aes(shape = limb), size = 1.7,
                      position = ggplot2::position_jitter(width = 0.12, height = 0),
                      alpha = 0.85) +
  ggplot2::scale_shape_manual(values = c("12W -> tumour limb" = 16, "6W (no claim)" = 1)) +
  ggplot2::facet_wrap(~ endpoint, scales = "free_y", ncol = 4) +
  ggplot2::labs(x = NULL, y = "score (z-composite) or log2 counts", shape = NULL) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 legend.position = "top", panel.grid.minor = ggplot2::element_blank())
ggplot2::ggsave(here::here("outputs", "fatpad_limb", "49_per_sample_trends.pdf"),
                p, width = 9, height = 5)

# -----------------------------------------------------------------------------
# 6W COMPOSITION QC -- not a claim of this script, which excludes the 6W groups.
# The per-sample record surfaced it and it bears on a claim made ELSEWHERE (the
# pre-check note's 6W genotype panel), so it is recorded here rather than left in a
# session transcript. Two of five 6WK_POS samples carry almost no epithelial signal;
# `Adipoq` is matched between the 6W groups but `Krt8`/`Epcam` are not, and adipose
# matching is not composition matching.
# -----------------------------------------------------------------------------
krt8 <- Mk$Krt8; epc <- Mk$Epcam
six_week_qc <- list(
  krt8_spread = dplyr::bind_rows(lapply(ALLGRP, function(g) {
    v <- krt8[grp == g]
    tibble::tibble(group = g, n = length(v), krt8_min = min(v), krt8_max = max(v),
                   krt8_range = diff(range(v)), epcam_min = min(epc[grp == g])) })),
  depleted = colnames(L)[grp == "6WK_POS" & krt8 < 11],
  panel = {
    k6 <- grp %in% c("6WK_NEG", "6WK_POS")
    d  <- data.frame(s = colnames(L)[k6], g = droplevels(grp[k6]),
                     r = Ap$puma_bclxl[k6], k = krt8[k6])
    f  <- function(dd, lab) tibble::tibble(
      subset = lab, n_neg = sum(dd$g == "6WK_NEG"), n_pos = sum(dd$g == "6WK_POS"),
      delta_median = stats::median(dd$r[dd$g == "6WK_POS"]) - stats::median(dd$r[dd$g == "6WK_NEG"]),
      delta_mean   = mean(dd$r[dd$g == "6WK_POS"]) - mean(dd$r[dd$g == "6WK_NEG"]),
      wilcox_p = suppressWarnings(stats::wilcox.test(dd$r[dd$g == "6WK_POS"],
                                                     dd$r[dd$g == "6WK_NEG"])$p.value))
    dplyr::bind_rows(f(d, "all 5 vs 5"), f(d[d$k >= 11, ], "Krt8-low dropped")) })

# =============================================================================
# PART G: THE ADJUSTED TREND -- ONE NUMBER AND AN INTERVAL, NOT FOUR ANSWERS
# -----------------------------------------------------------------------------
# Residualise the endpoint on adipose markers within the limb, then re-test the
# ordered trend on the residual. The bootstrap REFITS the adjustment inside each
# resample, so the interval carries the adjustment's own uncertainty rather than
# treating the residuals as fixed data. Read in the order D -> A -> B -> C.
# =============================================================================
message("49 PART G: adjusted trend")

ADJ <- c("Adipoq", "Cidec", "Fabp4", "Plin1")
AF  <- as.data.frame(lapply(Mk[ADJ], function(v) v[lk]))   # limb only, log2 scale

adj_fit <- function(y, cols, rows = seq_along(y)) {
  d <- cbind(data.frame(y = y[rows]), AF[rows, cols, drop = FALSE])
  stats::lm(y ~ ., d)
}
adj_trend <- function(v, cols, label, endpoint) {
  y  <- v[lk]
  f0 <- adj_fit(y, cols)
  r0 <- stats::resid(f0)
  tt <- jt_test(r0, gl, "increasing")
  bs <- vapply(seq_len(NBOOT), function(i) {
    k <- sample.int(length(y), replace = TRUE)
    if (length(unique(gi[k])) < 2L) return(NA_real_)
    fit <- try(adj_fit(y, cols, k), silent = TRUE)
    if (inherits(fit, "try-error") || anyNA(stats::coef(fit))) return(NA_real_)
    stats::cor(stats::resid(fit), gi[k], method = "kendall")
  }, 0)
  tibble::tibble(endpoint = endpoint, adjustment = label,
                 tau_adj = stats::cor(r0, gi, method = "kendall"),
                 tau_lo = unname(stats::quantile(bs, 0.025, na.rm = TRUE)),
                 tau_hi = unname(stats::quantile(bs, 0.975, na.rm = TRUE)),
                 p_one_sided = tt$p_one_sided,
                 adj_r2 = summary(f0)$adj.r.squared,
                 n_boot_ok = sum(!is.na(bs)))
}
raw_tau <- function(v) stats::cor(v[lk], gi, method = "kendall")

# --- D FIRST: the negative control. If proliferation's rise does not survive, the
# adjustment is over-fitted at n = 20 and nothing below may be read. ------------
adj_control <- dplyr::bind_rows(
  tibble::tibble(endpoint = "Mki67", adjustment = "none", tau_adj = raw_tau(Mk$Mki67),
                 tau_lo = NA_real_, tau_hi = NA_real_,
                 p_one_sided = trend$p_one_sided[trend$endpoint == "Mki67"],
                 adj_r2 = NA_real_, n_boot_ok = NA_integer_),
  adj_trend(Mk$Mki67, ADJ, "joint (4 markers)", "Mki67"))

# --- A: four single-marker adjustments, to show the spread --------------------
adj_single <- dplyr::bind_rows(lapply(ADJ, function(m) {
  out <- adj_trend(S$ox_nuc_mtrib, m, paste0("single: ", m), "ox_nuc_mtrib")
  out$marker_tau     <- trend$tau[trend$endpoint == m]
  out$marker_loading <- load_dir$rho_primary[load_dir$marker == m]
  out$predicted_push <- out$marker_tau * out$marker_loading
  out
}))
signs_disagree <- length(unique(sign(adj_single$tau_adj))) > 1L

# --- B and C: the joint adjustment, both rulers -------------------------------
adj_joint <- dplyr::bind_rows(
  adj_trend(S$ox_nuc_mtrib, ADJ, "joint (4 markers)", "ox_nuc_mtrib"),
  adj_trend(S$ox_rel,       ADJ, "joint (4 markers)", "ox_rel"))

# collinearity: the interval is wide partly because of this, which is worth stating
adj_vif <- tibble::tibble(marker = ADJ, vif = vapply(ADJ, function(m) {
  o <- setdiff(ADJ, m)
  1 / (1 - summary(stats::lm(stats::as.formula(paste(m, "~", paste(o, collapse = " + "))),
                             AF))$r.squared) }, 0))

# leave-one-out on the joint estimate, for completeness (R8 stays in)
adj_loo <- vapply(seq_len(sum(lk)), function(i) {
  y <- S$ox_nuc_mtrib[lk]
  stats::cor(stats::resid(adj_fit(y, ADJ, setdiff(seq_along(y), i))), gi[-i],
             method = "kendall") }, 0)
adj_loo_range <- c(full = adj_joint$tau_adj[1], lo = min(adj_loo), hi = max(adj_loo))

adjusted <- dplyr::bind_rows(adj_control, adj_single, adj_joint)

# =============================================================================
# PART F: THE VERDICT, ON THE RULES FIXED IN THE HEADER
# =============================================================================
message("49 PART F: verdict")

gv  <- function(n, col) trend[[col]][trend$endpoint == n]
t1_pass <- gv("ox_nuc_mtrib", "p_one_sided") < 0.05 && gv("ox_nuc_mtrib", "tau") > 0
t2_agree <- sign(gv("ox_rel", "tau")) == sign(gv("ox_nuc_mtrib", "tau"))

verdict <- tibble::tibble(
  test = c("T1 ox_nuc_mtrib rises (PRIMARY)", "T2 ox_rel agrees in direction",
           "T3 confound opposes the trend (HARD)", "T4 total vs share",
           "T5 Bcl2l1 direction vs primary"),
  result = c(
    if (t1_pass) "PASS" else "FAIL",
    if (t2_agree) "concordant" else "discordant",
    if (t3_opposes) "opposes -- T1 admissible" else "RUNS WITH THE TREND -- T1 VOID",
    sprintf("mito_total tau %+.3f (p %.3f) | ox_lvl tau %+.3f (p %.3f)",
            gv("mito_total", "tau"), gv("mito_total", "p_one_sided"),
            gv("ox_lvl", "tau"), gv("ox_lvl", "p_one_sided")),
    sprintf("Bcl2l1 tau %+.3f (p %.3f) vs primary tau %+.3f",
            gv("Bcl2l1", "tau"), gv("Bcl2l1", "p_one_sided"), prim_tau)),
  detail = c(
    sprintf("J %.0f, one-sided p %.4f, tau %+.3f [%.3f, %.3f], KW p %.3f",
            gv("ox_nuc_mtrib", "J"), gv("ox_nuc_mtrib", "p_one_sided"),
            gv("ox_nuc_mtrib", "tau"), gv("ox_nuc_mtrib", "tau_lo"),
            gv("ox_nuc_mtrib", "tau_hi"), gv("ox_nuc_mtrib", "kruskal_p")),
    sprintf("tau %+.3f [%.3f, %.3f], one-sided p %.4f",
            gv("ox_rel", "tau"), gv("ox_rel", "tau_lo"), gv("ox_rel", "tau_hi"),
            gv("ox_rel", "p_one_sided")),
    paste(sprintf("%s: tau %+.2f x loading %+.2f = push %+.2f",
                  adipose, adi_tau, adi_load, confound_push), collapse = " | "),
    "a share change on a flat total is the reprioritisation claim; both moving is a content claim",
    "no threshold: n = 20, observation not test"))

NOTES <- c(
  "ESTIMAND: a between-group ORDERED TREND across 12WK_POS < SMALL < LARGE (n = 20).",
  "  NOT the within-tumour correlation, which the pre-check closed at n = 15 and which",
  "  this script does not reopen. The 6W groups are scored for completeness only.",
  "ADMISSIBILITY (fixed in advance): adipose loads positively on the respiratory arm",
  "  and FALLS across the series, so composition alone predicts the score FALLING.",
  "  A rise runs against the confound and is credible; a NULL IS UNINTERPRETABLE.",
  "  PART D measures this for the PRIMARY endpoint, which the pre-check never did.",
  "SETS: script 08's splits rebuilt from MitoCarta Sheet 4 with every mt-* gene",
  "  stripped from every pathway -- not read from mitopps_scores.rds, whose lists",
  "  were already filtered to the MEC matrix.",
  "SCALE: z-composites on log2(counts + 1). mitoPPS is not run; it needs linear",
  "  counts and its own input object.",
  "PART G returns an INTERVAL, not a verdict. Read Mki67 (the negative control)",
  "  first: if adjustment kills a trend that is genuinely there, the adjustment is",
  "  over-fitted at n = 20 and nothing else in PART G may be read. A wide interval",
  "  spanning zero is UNINFORMATIVE, not negative.",
  "STANDING RULE: these scores are cohort-relative. They are never pooled, plotted",
  "  or differenced against the MEC cohort's or the orthotopic series'. Directions",
  "  and orderings travel; values never do.")

res <- list(set_sizes = set_sizes, trend = trend, group_summary = summ,
            load_direction = load_dir, loo = loo, loo_range = loo_range,
            six_week_qc = six_week_qc,
            adjusted = adjusted, adj_single = adj_single, adj_joint = adj_joint,
            adj_control = adj_control, adj_vif = adj_vif,
            adj_loo_range = adj_loo_range, signs_disagree = signs_disagree,
            per_sample = per_sample,
            verdict = verdict,
            t1_pass = t1_pass, t2_agree = t2_agree, t3_opposes = t3_opposes,
            confound_push = stats::setNames(confound_push, adipose),
            groups_in_limb = TREND, n_limb = sum(lk),
            params = list(NPERM = NPERM, NBOOT = NBOOT, MINBM = MINBM, seed = 1),
            analysis_date = Sys.Date(), notes = NOTES)

stopifnot(nrow(trend) == length(ALL), sum(lk) == 20L,
          all(c("ox_nuc_mtrib", "ox_rel") %in% trend$endpoint))

saveRDS(res, here::here("results", "fatpad_tumour_limb_trend.rds"))
message("49: wrote results/fatpad_tumour_limb_trend.rds")

cat("\n=============== VERDICT ===============\n")
print(as.data.frame(verdict[, c("test", "result")]))
cat("\n"); print(as.data.frame(trend[, c("endpoint", "tau", "tau_lo", "tau_hi",
                                         "p_one_sided", "kruskal_p")]), digits = 3)

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  res <- readRDS(here::here("results", "fatpad_tumour_limb_trend.rds"))

  ## READ IN THIS ORDER. T3 FIRST: if the confound runs with the trend, T1 is void
  ## and nothing else on the page means anything.
  res$load_direction |> print()
  res$confound_push  |> print()
  res$verdict |> as.data.frame() |> print()

  ## T1 and T2. tau with its bootstrap interval is the estimate; the p-value is the
  ## smaller half of the story at n = 20.
  res$trend |> print(n = 30)

  ## Is any of it one animal? LARGE_TUMOUR_R8 is the adipose-depleted sample.
  ## If rho(primary, Adipoq) stays high across every leave-one-out, T3 is structural.
  res$loo_range |> print()

  ## T4 -- the distinction the whole claim turns on. If mito_total is flat while
  ## ox_nuc_mtrib rises, the compartment is being REPRIORITISED. If both rise, this
  ## is a content claim and a weaker, different finding.
  subset(res$trend, endpoint %in% c("ox_nuc_mtrib", "ox_lvl", "mtrib_lvl", "mito_total"))

  ## Both statistics, because they diverge at these group sizes.
  res$group_summary |> print(n = 40)

  ## The per-sample values. At n = 20 look at the points before the p-value:
  ## outputs/fatpad_limb/49_per_sample_trends.pdf
  res$per_sample |> print(n = 30)

  ## NOT this script's claim, but recorded because the per-sample table surfaced it:
  ## two of five 6WK_POS samples have almost no epithelium (Krt8 102x and 23x below
  ## that group's own median), and they carry its two HIGHEST PUMA:BCL-XL values.
  ## The pre-check's 6W panel needs the caveat -- direction survives, p does not.
  res$six_week_qc |> print()

  ## T5 -- the guardian against the respiratory trend. If the share rises INTO
  ## tumours and Bcl2l1 tracks it, that is the human configuration appearing as a
  ## mouse trajectory. Observation, not test.
  subset(res$trend, endpoint %in% c("Bcl2l1", "Bbc3", "puma_bclxl", "buffer", "Myc"))

  ## --- PART G, THE ADDENDUM THAT CLOSES THIS DATASET ----------------------
  ## READ THE NEGATIVE CONTROL FIRST. Mki67 genuinely rises on this limb. If the
  ## joint adjustment removes that too, the adjustment is eating signal at n = 20
  ## and nothing else below is readable.
  res$adj_control |> print()

  ## A -- four single-marker adjustments. The SPREAD is the point: if they disagree
  ## in sign, no single-marker adjustment is authoritative and that disagreement is
  ## the result. `predicted_push` is marker_tau x marker_loading, i.e. T3's arithmetic.
  res$adj_single |> print()
  res$signs_disagree |> print()

  ## B and C -- one number and an interval, both rulers. adj_r2 says how much of the
  ## endpoint the four markers absorb; the VIFs say why the interval is wide.
  res$adj_joint |> print()
  res$adj_vif |> print()
  res$adj_loo_range |> print()

  cat(res$notes, sep = "\n")
}
