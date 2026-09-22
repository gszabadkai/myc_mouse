# =============================================================================
# 50 -- THE ORTHOTOPIC VECTOR SERIES: LICENSING, CEILING, AND WHICH SIDE OF THE
#       REVERSAL MOUSE TUMOURS FALL ON
# -----------------------------------------------------------------------------
# Cohort: orthotopic implants of MMTV-Myc tumour-derived (MYAZ) cells,
# `data/orthotopic_series/`. Read that folder's README before this file. NOT the
# 6W/12W purified-MEC dataset, NOT the fat-pad timeline.
#
# WHAT THIS DATASET CAN DO, and one thing it cannot.
#
# IT CANNOT ESTIMATE THE MYC x OXPHOS INTERACTION, and that is said first. Every
# arm is MYAZ-derived and MYC-high; there is no MYC-low arm, so the interaction is
# not identifiable here at any n. It lives in the iMMEC rtTA-MYC +/-dox x +/-PGC1a
# design, on the death readout. This script does not fit it, and no group term is
# allowed to stand in for it.
#
# C1 -- THE LICENSING RELATIONSHIP, BY INTERVENTION. PGC1a overexpression alone,
#   with no Bcl-xL construct present, raises `Bcl2l1`. That is a positive
#   association between the guardian and respiration PRODUCED by manipulating
#   respiration, not observed by correlating it -- the argument that the human
#   BCL2L1 coefficient reflects a causal link rather than a cross-sectional one.
#   Directional, one-sided, n = 7 v 7.
#
# C2 -- THE MYC CEILING. A claim about the UPPER TAIL, not the mean: no PGC1a
#   tumour exceeds a `Myc` level that EV and BclxL tumours reach. A shifted median
#   and a ceiling are different statements and the medians do not distinguish them.
#   Both the group MAX and the 80th-percentile contrast are reported, side by side,
#   under one exact permutation null, so it is visible whether the ceiling rests on
#   a single animal. AT n = 7 THIS IS A FIGURE AND AN OBSERVATION, NOT A TEST, and
#   it is written that way.
#
# C3 -- WHICH SIDE OF THE REVERSAL. `Mcl1:Bcl2l1` and both members separately
#   against `ox_rel` across the vector series. In the gland the same endpoint on
#   the same ruler is +0.351; in human tumours it is -0.31 to -0.46. These are
#   mouse TUMOURS. NEGATIVE means they behave like human tumours and the abstract's
#   "mouse tumours resemble human tumours rather than the gland they arose from"
#   survives; POSITIVE means they behave like the gland and THAT SENTENCE MUST BE
#   WITHDRAWN. Both readings are pre-specified as informative. The reference values
#   are quoted for direction only and are NEVER pooled with these scores.
#
# THE SCORING SET IS FIXED AT n = 21 (EV + BclxL + Pgc1a), decided before any
# number was computed. `ox_rel` is a relative score and every value depends on who
# is in the cohort. The p21 series is out of scope. MYAZ / FaMY / BoMY are a
# separate series on the file evidence in the README and are excluded; the wider
# set appears only as a named sensitivity in PART B.
#
# GATE, READ BEFORE ANY CLAIM (PART B). The fat-pad failure's transferable rule:
# before testing anything, check whether the composition axis and the biological
# axis are separable AT ALL, computed on the exact subset the tests run on. The
# POSITIVE CONTROL is `Ppargc1a` itself -- a transgene that is definitionally
# present in the Pgc1a arm. If the composition adjustment attenuates THAT, the
# adjustment is removing signal and PART F cannot be read.
#
# N3. These are transcript associations. The word "primed" is not written of a
# transcript anywhere in this script, its comments, its labels or its outputs.
#
# LIMITATION, STATED NOT GLOSSED. No transcript-level quantification exists on
# disk -- no quant.sf, no salmon directories, no tx2gene. `tximport` is
# unavailable, so counts go through DESeqDataSetFromMatrix and the transcript
# length offset is lost. It applies to every number this script produces.
#
# Reads : data/orthotopic_series/salmon.merged.gene_counts.tsv (salmon merged,
#           gene level, NON-INTEGER estimated counts)
#         data/Mouse.MitoCarta3.0.xls  Sheet 4  (script 08's authoritative splits)
#         functions/reconcile_gene_symbols.R (MANDATORY -- vintage-aware membership)
# Writes: results/orthotopic_vector_series.rds
#
# RUNTIME: about a minute. The C2 null is EXACT -- all 3432 splits of 14 samples.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

set.seed(1)
NBOOT   <- 5000L   # bootstrap resamples for every interval
MINBM   <- 10      # smallest mean normalised count carried into a composite
CUT_MAT <- 0.5     # |rho| at which a composition axis counts as MATERIAL (PART B)
VEC     <- c("EV", "BclxL", "Pgc1a")            # the scoring set, fixed
ALLGRP  <- c(VEC, "NTEV", "KOEV", "NTPgc1a", "KOPgc1a", "MYAZ", "FaMY", "BoMY")

# --- the load control, from the analysis plan's section 2. If these do not come
# back the load or the group parsing has drifted and nothing below is readable.
REF <- data.frame(
  group    = ALLGRP,
  n        = c(7L, 7L, 7L, 7L, 7L, 7L, 7L, 6L, 6L, 6L),
  Bcl2l1   = c(71, 326, 193, 82, 62, 79, 79, 59, 59, 51),
  Ppargc1a = c(0.0, 0.3, 302, 0.2, 0.1, 129, 0.2, 0.0, 0.0, 0.0),
  Myc      = c(1915, 2033, 1440, 2208, 2368, 2286, 2176, 1370, 1280, 1293),
  stringsAsFactors = FALSE)
CTRL_TOL_REL <- 0.02      # 2 pct, the rounding in the reference table

# =============================================================================
# PART 0: LOAD, PARSE, AND PROVE THE LOAD REPRODUCES THE REFERENCE TABLE
# =============================================================================
message("50 PART 0: load")

raw <- data.table::fread(
  here::here("data", "orthotopic_series", "salmon.merged.gene_counts.tsv"),
  data.table = FALSE, check.names = FALSE)
stopifnot(identical(colnames(raw)[1:2], c("gene_id", "gene_name")))
CNT <- as.matrix(raw[, -(1:2)])
rownames(CNT) <- raw$gene_id
sym_of_id <- stats::setNames(raw$gene_name, raw$gene_id)
stopifnot(ncol(CNT) == 67L, all(grepl("^ENSMUSG", rownames(CNT))))

# GROUP PARSING. The tokens are not written the same way twice -- `EV_BRNS` and
# `EVBRNS`, `Pgc1a` and `Pgc`, `BoMY1` and `BoMY`. Match LONGEST-FIRST: `EV` and
# `Pgc` must be matched AFTER `NTEV`/`KOEV`/`NTPgc`/`KOPgc` or four groups
# silently collapse into two. Splitting on `_` or `-` does not work.
PAT <- c("NTEV", "KOEV", "NTPgc1a", "NTPgc", "KOPgc1a", "KOPgc", "BclxL",
         "BoMY1", "BoMY", "FaMY1", "FaMY", "MYAZ", "Pgc1a", "Pgc", "EV")
CANON <- c(NTEV = "NTEV", KOEV = "KOEV", NTPgc1a = "NTPgc1a", NTPgc = "NTPgc1a",
           KOPgc1a = "KOPgc1a", KOPgc = "KOPgc1a", BclxL = "BclxL",
           BoMY1 = "BoMY", BoMY = "BoMY", FaMY1 = "FaMY", FaMY = "FaMY",
           MYAZ = "MYAZ", Pgc1a = "Pgc1a", Pgc = "Pgc1a", EV = "EV")
core <- sub("^[A-H][0-9]-", "", colnames(CNT))
tok  <- vapply(core, function(s) {
  h <- PAT[vapply(PAT, function(p) grepl(paste0("^", p), s), TRUE)]
  if (length(h)) h[1] else NA_character_ }, "")
grp  <- factor(unname(CANON[tok]), levels = ALLGRP)
stopifnot(!anyNA(grp), identical(as.integer(table(grp)), REF$n))

CPM <- t(t(CNT) / colSums(CNT)) * 1e6
cpm_of <- function(sym) {
  r <- which(raw$gene_name == sym)
  stopifnot(length(r) >= 1L)
  if (length(r) == 1L) CPM[r, ] else colSums(CPM[r, , drop = FALSE])
}
load_control <- data.frame(
  group = ALLGRP, n = as.integer(table(grp)),
  Bcl2l1   = round(as.numeric(tapply(cpm_of("Bcl2l1"),   grp, stats::median)), 1),
  Ppargc1a = round(as.numeric(tapply(cpm_of("Ppargc1a"), grp, stats::median)), 1),
  Myc      = round(as.numeric(tapply(cpm_of("Myc"),      grp, stats::median)), 1),
  stringsAsFactors = FALSE)
dev <- max(abs(unlist(load_control[, 3:5]) - unlist(REF[, 3:5])) /
           pmax(unlist(REF[, 3:5]), 1))
message(sprintf("50 PART 0: load control, max relative deviation %.3f%% (tol %.1f%%)",
                100 * dev, 100 * CTRL_TOL_REL))
stopifnot(dev < CTRL_TOL_REL)

# =============================================================================
# PART A: THE COHORT, THE TWO INPUT OBJECTS, AND THE RULERS
# -----------------------------------------------------------------------------
# The scoring set is fixed at the 21 vector-series samples. Counts are ROUNDED for
# DESeq2 because salmon merged counts are non-integer; the transcript-length
# offset is already lost (see the header), and rounding adds nothing to that.
#
# TWO OBJECTS THAT NEVER MIX: `NC` is LINEAR normalised counts (what mitoPPS would
# take); `L` is log2(NC + 1) and feeds every z-composite. mitoPPS and GSVA are NOT
# run -- no pre-specified claim needs either, and computing a score nothing rests
# on is how a cohort-relative number ends up quoted by accident.
# =============================================================================
message("50 PART A: cohort and rulers")

vk  <- grp %in% VEC
sm  <- data.frame(sample = colnames(CNT)[vk],
                  group  = factor(as.character(grp[vk]), levels = VEC))
stopifnot(nrow(sm) == 21L, identical(as.integer(table(sm$group)), c(7L, 7L, 7L)))

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(CNT[, vk, drop = FALSE]), colData = sm, design = ~ group)
dds <- DESeq2::estimateSizeFactors(dds)
NCv <- DESeq2::counts(dds, normalized = TRUE)          # LINEAR -- mitoPPS only
keep <- rowMeans(NCv) >= MINBM & matrixStats::rowSds(NCv) > 0
L    <- log2(NCv[keep, , drop = FALSE] + 1)            # LOG -- every composite
message(sprintf("50 PART A: %d of %d genes at mean normalised count >= %d",
                nrow(L), nrow(NCv), MINBM))

universe_all <- rownames(L)
zrow    <- function(m) t(scale(t(m)))
comp_e  <- function(e) { e <- e[!is.na(e) & e %in% rownames(L)]
                         stopifnot(length(e) >= 3); colMeans(zrow(L[e, , drop = FALSE])) }
ens_set <- function(syms) { e <- recon_to_ensembl(syms, universe_all); e[!is.na(e)] }
one_ens <- function(sym) { e <- ens_set(sym)
                           if (!length(e)) NA_character_ else e[1] }
lg <- function(sym) { e <- one_ens(sym); stopifnot(!is.na(e)); as.numeric(L[e, ]) }

# --- MitoCarta, script 08's authoritative split, rebuilt from Sheet 4 ---------
# Exact filename, never a glob (CLAUDE.md's data trap). Every mt-* gene is
# stripped from every pathway as 08 does, so the OXPHOS numerator is nuclear by
# construction and the 13 mtDNA genes land in the ox_rel DENOMINATOR.
s4  <- readxl::read_xls(here::here("data", "Mouse.MitoCarta3.0.xls"), sheet = 4)
s4  <- na.omit(dplyr::select(s4, MitoPathway, Genes))
g2p <- suppressWarnings(splitstackshape::cSplit(s4, "Genes", ","))
g2p <- as.data.frame(t(tibble::column_to_rownames(g2p, "MitoPathway")))
g2p <- tidyr::pivot_longer(g2p, cols = colnames(g2p),
                           names_to = "Pathway", values_to = "Gene")
g2p <- dplyr::mutate(na.omit(g2p), Gene = as.character(Gene))

is_mt   <- grepl("^[Mm][Tt]-", g2p$Gene)
mt_syms <- unique(g2p$Gene[is_mt])
g2p_nuc <- g2p[!is_mt, ]
ox_sub    <- ens_set(unique(g2p_nuc$Gene[g2p_nuc$Pathway == "OXPHOS subunits"]))
mito_nuc  <- ens_set(unique(g2p_nuc$Gene))
rest_mito <- setdiff(mito_nuc, ox_sub)
mt_ens    <- ens_set(mt_syms)

stopifnot(length(intersect(ox_sub, mt_ens)) == 0L,   # numerator nuclear-only
          length(mt_ens) >= 10L)
# the plan says 88 nuclear subunits; report what ACTUALLY maps here rather than
# forcing the number, and flag a shortfall instead of hiding it
set_sizes <- tibble::tibble(
  set = c("OXPHOS subunits (nuclear)", "rest of nuclear MitoCarta",
          "mtDNA-encoded", "genes carried"),
  n   = c(length(ox_sub), length(rest_mito), length(mt_ens), nrow(L)),
  note = c(sprintf("plan states 88; %d map in this cohort", length(ox_sub)),
           "the ox_rel denominator", "in the denominator, never the numerator", ""))
message(sprintf("50 PART A: ox_sub %d (plan says 88) | rest %d | mtDNA %d",
                length(ox_sub), length(rest_mito), length(mt_ens)))

X <- list(ox_rel = comp_e(ox_sub) - comp_e(rest_mito),
          ox_lvl = comp_e(ox_sub),
          ox_mt  = comp_e(mt_ens))

# --- the endpoints -----------------------------------------------------------
G <- list(Bcl2l1 = lg("Bcl2l1"), Mcl1 = lg("Mcl1"), Bbc3 = lg("Bbc3"),
          Bcl2l11 = lg("Bcl2l11"), Myc = lg("Myc"), Ppargc1a = lg("Ppargc1a"))
G$guardian_ratio <- G$Mcl1 - G$Bcl2l1          # log2(Mcl1/Bcl2l1)
cpm_v <- lapply(c(Bcl2l1 = "Bcl2l1", Ppargc1a = "Ppargc1a", Myc = "Myc"),
                function(s) cpm_of(s)[vk])     # CPM, for continuity with the plan

# =============================================================================
# PART B: THE COMPOSITION GATE -- READ THIS BEFORE ANY CLAIM
# -----------------------------------------------------------------------------
# Six compartment composites plus their single markers, against all three rulers,
# INSIDE the n = 21 scoring set. The fat-pad lesson is that a loading measured on
# a different subset does not transfer, so nothing here is inherited.
# A compartment counts as MATERIAL if |rho| with ox_rel >= CUT_MAT; those and only
# those enter the PART F adjustment. The rule is fixed here, before the values.
# =============================================================================
message("50 PART B: composition gate")

MARK <- list(
  epithelial    = c("Krt8", "Krt18", "Epcam", "Cdh1", "Krt5", "Krt14"),
  stromal       = c("Col1a1", "Pdgfrb", "Acta2", "Thy1"),
  adipose       = c("Adipoq", "Plin1", "Fabp4", "Cidec"),
  immune        = c("Ptprc", "Cd52", "Lyz2", "Cd74"),
  endothelial   = c("Pecam1", "Cdh5", "Cldn5", "Emcn", "Tek"),
  proliferation = c("Mki67", "Top2a", "Ccnb1"))
MARK_e <- lapply(MARK, function(v) { e <- ens_set(v); e[e %in% rownames(L)] })
comp_present <- vapply(MARK_e, length, 0L)
# a "composite" of one or two genes is a gene, so the >= 3 floor in comp_e stands
# and any compartment that cannot reach it is DROPPED and named, not quietly
# rescued by lowering the bar
dropped_comp <- names(comp_present)[comp_present < 3L]
if (length(dropped_comp))
  message(sprintf("50 PART B: dropped (fewer than 3 markers detected): %s",
                  paste(dropped_comp, collapse = ", ")))
K <- lapply(MARK_e[comp_present >= 3L], comp_e)   # compartment composites

sp <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman"))
loading <- dplyr::bind_rows(lapply(names(K), function(cn)
  dplyr::bind_rows(lapply(names(X), function(xn) tibble::tibble(
    compartment = cn, n_genes = comp_present[[cn]], ruler = xn,
    rho_all21 = sp(X[[xn]], K[[cn]]),
    rho_EV    = sp(X[[xn]][sm$group == "EV"],    K[[cn]][sm$group == "EV"]),
    rho_BclxL = sp(X[[xn]][sm$group == "BclxL"], K[[cn]][sm$group == "BclxL"]),
    rho_Pgc1a = sp(X[[xn]][sm$group == "Pgc1a"], K[[cn]][sm$group == "Pgc1a"]))))))

material <- loading$compartment[loading$ruler == "ox_rel" &
                                abs(loading$rho_all21) >= CUT_MAT]
message(sprintf("50 PART B: material composition axes (|rho| >= %.1f): %s",
                CUT_MAT, if (length(material)) paste(material, collapse = ", ") else "NONE"))

# --- the SENSITIVITY: the same table over all 67, to show what the cohort
# choice changes. Scored in its own run -- these values are NOT comparable with
# the ones above and are never mixed with them.
NCa  <- DESeq2::counts(DESeq2::estimateSizeFactors(DESeq2::DESeqDataSetFromMatrix(
          round(CNT), data.frame(group = grp), ~ 1)), normalized = TRUE)
ka   <- rowMeans(NCa) >= MINBM & matrixStats::rowSds(NCa) > 0
La   <- log2(NCa[ka, , drop = FALSE] + 1)
comp_a <- function(e) { e <- e[e %in% rownames(La)]; stopifnot(length(e) >= 3)
                        colMeans(zrow(La[e, , drop = FALSE])) }
oxrel_a <- comp_a(ox_sub) - comp_a(rest_mito)
loading_all67 <- dplyr::bind_rows(lapply(names(MARK_e), function(cn) {
  e <- MARK_e[[cn]][MARK_e[[cn]] %in% rownames(La)]
  if (length(e) < 3L) return(NULL)
  tibble::tibble(compartment = cn, rho_ox_rel_all67 = sp(oxrel_a, comp_a(e))) }))

# =============================================================================
# PART C: C1 -- THE LICENSING RELATIONSHIP
# -----------------------------------------------------------------------------
# One-sided by pre-specification: Bcl2l1 HIGHER in Pgc1a than in EV. BclxL is
# reported as the third arm and never pooled with either.
# =============================================================================
message("50 PART C: C1")

boot_diff <- function(a, b, f = stats::median) {
  d <- vapply(seq_len(NBOOT), function(i)
    f(sample(a, replace = TRUE)) - f(sample(b, replace = TRUE)), 0)
  stats::quantile(d, c(0.025, 0.975), na.rm = TRUE) }
two_group <- function(v, hi, lo, one_sided = TRUE) {
  a <- v[sm$group == hi]; b <- v[sm$group == lo]
  ci <- boot_diff(a, b)
  tibble::tibble(hi = hi, lo = lo, n_hi = length(a), n_lo = length(b),
                 median_hi = stats::median(a), median_lo = stats::median(b),
                 diff_median = stats::median(a) - stats::median(b),
                 ci_lo = ci[[1]], ci_hi = ci[[2]],
                 wilcox_p = suppressWarnings(stats::wilcox.test(
                   a, b, alternative = if (one_sided) "greater" else "two.sided")$p.value)) }

c1 <- dplyr::bind_rows(
  dplyr::bind_cols(tibble::tibble(measure = "Bcl2l1 (CPM)"),
                   two_group(cpm_v$Bcl2l1, "Pgc1a", "EV")),
  dplyr::bind_cols(tibble::tibble(measure = "Ppargc1a (CPM, manipulation check)"),
                   two_group(cpm_v$Ppargc1a, "Pgc1a", "EV")),
  dplyr::bind_cols(tibble::tibble(measure = "ox_lvl (OXPHOS level)"),
                   two_group(X$ox_lvl, "Pgc1a", "EV")),
  dplyr::bind_cols(tibble::tibble(measure = "ox_rel (respiratory share)"),
                   two_group(X$ox_rel, "Pgc1a", "EV")),
  # the third arm, reported not pooled
  dplyr::bind_cols(tibble::tibble(measure = "Bcl2l1 (CPM) -- BclxL arm, construct present"),
                   two_group(cpm_v$Bcl2l1, "BclxL", "EV")))

# =============================================================================
# PART D: C2 -- THE MYC CEILING, AS AN UPPER-TAIL STATEMENT
# -----------------------------------------------------------------------------
# TWO statistics side by side because they fail differently: the MAX is the
# literal ceiling and is hostage to one animal; the 80th percentile is not, and
# at n = 7 sits between the two highest samples. The null is EXACT -- every one of
# the choose(14,7) = 3432 relabellings of the two groups being contrasted.
# A CEILING AT n = 7 IS A FIGURE AND AN OBSERVATION, NOT A TEST.
# =============================================================================
message("50 PART D: C2, the ceiling")

exact_null <- function(v, hi, lo, stat) {
  a <- v[sm$group == hi]; b <- v[sm$group == lo]
  pool <- c(a, b); n <- length(a)
  idx  <- utils::combn(length(pool), n)
  null <- apply(idx, 2, function(k) stat(pool[k]) - stat(pool[-k]))
  obs  <- stat(a) - stat(b)
  list(obs = obs, p_lower = mean(null <= obs), n_perm = ncol(idx),
       null_median = stats::median(null)) }
p80 <- function(z) unname(stats::quantile(z, 0.8, type = 7))

c2 <- dplyr::bind_rows(lapply(list(c("Pgc1a", "EV"), c("Pgc1a", "BclxL")), function(pr)
  dplyr::bind_rows(lapply(list(max = max, pct80 = p80), function(st) {
    r <- exact_null(cpm_v$Myc, pr[1], pr[2], st)
    tibble::tibble(hi = pr[1], lo = pr[2], stat_obs = r$obs,
                   p_one_sided_lower = r$p_lower, n_perm = r$n_perm,
                   null_median = r$null_median) }), .id = "statistic")))
# the companion ceiling/floor reading: how tightly each arm converges
c2_spread <- dplyr::bind_rows(lapply(VEC, function(g) {
  b <- cpm_v$Bcl2l1[sm$group == g]; m <- cpm_v$Myc[sm$group == g]
  tibble::tibble(group = g, Bcl2l1_min = min(b), Bcl2l1_max = max(b),
                 Bcl2l1_fold_spread = max(b) / min(b),
                 Myc_max = max(m), Myc_n_over_2000 = sum(m > 2000)) }))

# =============================================================================
# PART E: C3 -- WHICH SIDE OF THE REVERSAL
# -----------------------------------------------------------------------------
# BOTH readings are pre-specified as informative. Reported two ways, because they
# answer different questions and the note must say which it is quoting:
#   pooled       -- across all 21, which INCLUDES the between-group differences
#                   the constructs created (Pgc1a raises ox_lvl by design)
#   group-adjusted -- ox_rel and the endpoint centred WITHIN group, so it is the
#                   within-arm association and carries no construct effect
# Reference directions, never pooled with these values: gland +0.351, human
# tumours -0.31 to -0.46.
# =============================================================================
message("50 PART E: C3, the reversal")

ctr_g <- function(v) { o <- v; for (g in VEC) { k <- sm$group == g
                       o[k] <- v[k] - mean(v[k]) }; o }
assoc <- function(y, x, lab, endpoint) {
  r  <- sp(y, x)
  bs <- vapply(seq_len(NBOOT), function(i) {
    k <- sample.int(length(y), replace = TRUE)
    if (length(unique(x[k])) < 3L) return(NA_real_)
    sp(y[k], x[k]) }, 0)
  lo <- vapply(seq_along(y), function(i) sp(y[-i], x[-i]), 0)
  tibble::tibble(endpoint = endpoint, reading = lab, rho = r,
                 ci_lo = unname(stats::quantile(bs, 0.025, na.rm = TRUE)),
                 ci_hi = unname(stats::quantile(bs, 0.975, na.rm = TRUE)),
                 loo_min = min(lo), loo_max = max(lo),
                 same_sign_loo = all(sign(lo) == sign(r))) }

EP <- list(guardian_ratio = G$guardian_ratio, Mcl1 = G$Mcl1, Bcl2l1 = G$Bcl2l1)
c3 <- dplyr::bind_rows(lapply(names(EP), function(en) dplyr::bind_rows(
  assoc(EP[[en]],        X$ox_rel,        "pooled (n = 21)",   en),
  assoc(ctr_g(EP[[en]]), ctr_g(X$ox_rel), "group-adjusted",    en))))
c3_ref <- c(gland_within_timepoint = 0.351, human_tumours_lo = -0.46, human_tumours_hi = -0.31)

# =============================================================================
# PART F: THE FAILURE READING -- DOES ANY OF IT SURVIVE COMPOSITION?
# -----------------------------------------------------------------------------
# C1 and C3 repeated with the MATERIAL composition axes from PART B as
# covariates, AND `Ppargc1a` carried through the identical adjustment as the
# POSITIVE CONTROL. Three outcomes, distinguished before looking:
#   control holds, claim holds     -> the claim is not a cellularity artefact
#   control holds, claim reverses  -> the group contrast IS confounded with
#                                     cellularity; the dataset contributes an
#                                     EXCLUSION, not a result
#   control collapses              -> the adjustment is removing signal (the
#                                     fat-pad Mki67 lesson) and NOTHING in this
#                                     part may be read, in either direction
# =============================================================================
message("50 PART F: composition adjustment and its positive control")

adj_cols <- if (length(material)) unique(material) else character(0)
AD <- if (length(adj_cols)) as.data.frame(K[adj_cols]) else NULL

adj_group <- function(v, lab) {
  d <- data.frame(y = v, g = sm$group)
  if (!is.null(AD)) d <- cbind(d, AD)
  f <- stats::lm(y ~ ., d)
  s <- summary(f)$coefficients
  rn <- grep("^gPgc1a$", rownames(s), value = TRUE)
  tibble::tibble(measure = lab, adjusted_for = if (length(adj_cols)) paste(adj_cols, collapse = "+") else "nothing (no material axis)",
                 beta_Pgc1a_vs_EV = if (length(rn)) s[rn, 1] else NA_real_,
                 p = if (length(rn)) s[rn, 4] else NA_real_) }
adj_assoc <- function(y, lab) {
  d <- data.frame(y = y, x = X$ox_rel)
  if (!is.null(AD)) d <- cbind(d, AD)
  s <- summary(stats::lm(y ~ ., d))$coefficients
  tibble::tibble(endpoint = lab, beta_ox_rel = s["x", 1], p = s["x", 4]) }

adjusted <- list(
  c1 = dplyr::bind_rows(
    adj_group(cpm_v$Bcl2l1,   "Bcl2l1 (CPM)"),
    adj_group(cpm_v$Ppargc1a, "Ppargc1a (CPM) -- POSITIVE CONTROL"),
    adj_group(X$ox_lvl,       "ox_lvl")),
  c3 = dplyr::bind_rows(lapply(names(EP), function(en) adj_assoc(EP[[en]], en))))

# =============================================================================
# PART G: VERDICT
# =============================================================================
message("50 PART G: verdict")

gate_ok    <- length(material) == 0L
ctrl_row   <- adjusted$c1[grepl("POSITIVE CONTROL", adjusted$c1$measure), ]
ctrl_holds <- is.na(ctrl_row$p[1]) || (ctrl_row$beta_Pgc1a_vs_EV[1] > 0 & ctrl_row$p[1] < 0.05)
c1_row     <- c1[c1$measure == "Bcl2l1 (CPM)", ]
c3_row     <- c3[c3$endpoint == "guardian_ratio" & c3$reading == "pooled (n = 21)", ]
# THE READING RULE, three-way and fixed here rather than inferred from the sign.
# An interval spanning zero means the dataset does not place these tumours on
# either side, and the abstract sentence is then NEITHER supported NOR withdrawn.
c3_excludes_zero <- c3_row$ci_lo > 0 | c3_row$ci_hi < 0
c3_signs <- c(pooled = sign(c3_row$rho),
              group_adj = sign(c3$rho[c3$endpoint == "guardian_ratio" &
                                      c3$reading == "group-adjusted"]),
              comp_adj = sign(adjusted$c3$beta_ox_rel[adjusted$c3$endpoint == "guardian_ratio"]))
c3_sign_stable <- length(unique(c3_signs)) == 1L
c3_call <- if (!c3_excludes_zero) {
  "UNINFORMATIVE -- interval spans zero; the sentence is neither supported nor withdrawn"
} else if (c3_row$rho < 0) {
  "HUMAN-LIKE -- abstract sentence survives"
} else {
  "GLAND-LIKE -- abstract sentence must be withdrawn"
}

verdict <- tibble::tibble(
  item = c("GATE  composition separable?",
           "C1    Bcl2l1 higher in Pgc1a than EV",
           "C2    MYC ceiling (observation, not a test)",
           "C3    which side of the reversal",
           "F     positive control survives adjustment"),
  result = c(
    if (gate_ok) sprintf("clean -- no compartment reaches |rho| %.1f with ox_rel", CUT_MAT)
    else sprintf("MATERIAL: %s -- PART F decides what survives", paste(material, collapse = ", ")),
    sprintf("%+.1f CPM [%.1f, %.1f], one-sided p %.4f",
            c1_row$diff_median, c1_row$ci_lo, c1_row$ci_hi, c1_row$wilcox_p),
    sprintf("max %+.0f (p %.3f) | 80th pct %+.0f (p %.3f), Pgc1a vs EV",
            c2$stat_obs[c2$statistic == "max" & c2$lo == "EV"],
            c2$p_one_sided_lower[c2$statistic == "max" & c2$lo == "EV"],
            c2$stat_obs[c2$statistic == "pct80" & c2$lo == "EV"],
            c2$p_one_sided_lower[c2$statistic == "pct80" & c2$lo == "EV"]),
    sprintf("rho %+.3f [%.3f, %.3f]; sign across the three readings %s -- %s",
            c3_row$rho, c3_row$ci_lo, c3_row$ci_hi,
            if (c3_sign_stable) "stable" else "FLIPS", c3_call),
    if (!length(adj_cols)) "not applicable -- no material axis to adjust for"
    else if (ctrl_holds) "control holds; the claims above are readable after adjustment"
    else "CONTROL COLLAPSES -- the adjustment removes signal; PART F unreadable"))

NOTES <- c(
  "COHORT: the vector series alone (EV/BclxL/Pgc1a, n = 21), FIXED before any number",
  "  was computed. ox_rel is relative and GSVA is cohort-relative, so this choice sets",
  "  every value. p21 series out of scope; MYAZ/FaMY/BoMY excluded as a separate series",
  "  on the file evidence in data/orthotopic_series/README.md.",
  "NOT ESTIMABLE HERE: the MYC x OXPHOS interaction. Every arm is MYAZ-derived and",
  "  MYC-high; there is no MYC-low arm. It lives in the iMMEC rtTA-MYC +/-dox x",
  "  +/-PGC1a design on the death readout. No group term stands in for it.",
  "LIMITATION: no transcript-level quantification exists on disk, so tximport is",
  "  unavailable and the transcript-length offset is lost. Applies to every number here.",
  "C2 IS AN OBSERVATION, NOT A TEST. n = 7 per arm. Both the max and the 80th",
  "  percentile are reported so it is visible whether the ceiling rests on one animal.",
  "C3 READS BOTH WAYS BY PRE-SPECIFICATION. Negative = mouse tumours resemble human",
  "  tumours and the abstract sentence survives; positive = they resemble the gland and",
  "  it must be withdrawn. Gland +0.351 and human -0.31/-0.46 are DIRECTIONS only and",
  "  are never pooled with these cohort-relative scores.",
  "N3: these are transcript associations. The word 'primed' is not used of a transcript.",
  "SPECIES = COHORT: never pool these scores with the 6W/12W timeline or with human.")

stopifnot(nrow(sm) == 21L, nrow(c1) >= 4L, nrow(c3) == 6L,
          length(c3_signs) == 3L,
          all(c2$n_perm == choose(14, 7)))

res <- list(load_control = load_control, reference_table = REF,
            set_sizes = set_sizes, scoring_set = VEC, sample_table = sm,
            loading = loading, loading_all67 = loading_all67,
            marker_panels = tibble::tibble(compartment = names(comp_present),
                                           n_detected = as.integer(comp_present),
                                           used = comp_present >= 3L),
            material_axes = material, cut_material = CUT_MAT,
            c1 = c1, c2 = c2, c2_spread = c2_spread,
            c3 = c3, c3_reference_directions = c3_ref,
            c3_call = c3_call, c3_sign_stable = c3_sign_stable,
            c3_signs = c3_signs, c3_excludes_zero = c3_excludes_zero,
            adjusted = adjusted, verdict = verdict,
            per_sample = dplyr::bind_cols(
              sm, tibble::as_tibble(X), tibble::as_tibble(G),
              tibble::as_tibble(stats::setNames(cpm_v, paste0(names(cpm_v), "_CPM")))),
            params = list(NBOOT = NBOOT, MINBM = MINBM, CUT_MAT = CUT_MAT, seed = 1),
            analysis_date = Sys.Date(), notes = NOTES)

saveRDS(res, here::here("results", "orthotopic_vector_series.rds"))
message("50: wrote results/orthotopic_vector_series.rds")

cat("\n================ VERDICT ================\n")
print(as.data.frame(verdict), right = FALSE)

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  res <- readRDS(here::here("results", "orthotopic_vector_series.rds"))

  ## 1. THE LOAD CONTROL. If this is not the analysis plan's section 2 table, the
  ## load or the group parsing has drifted and nothing below is readable.
  res$load_control |> print()

  ## 2. THE GATE, BEFORE ANY CLAIM. Is the composition axis separable from the
  ## respiratory axis at all, inside the 21 samples the tests actually run on?
  res$loading |> print(n = 30)
  res$material_axes |> print()
  ## what the cohort choice changes -- NOT comparable with the table above
  res$loading_all67 |> print()

  ## 3. C1 -- the licensing relationship. Read Ppargc1a and ox_lvl alongside: they
  ## show the manipulation worked at the level it claims to. BclxL is the third
  ## arm, reported and never pooled.
  res$c1 |> print()

  ## 4. C2 -- the ceiling. An OBSERVATION. Compare max against the 80th percentile:
  ## if only the max moves, the ceiling is one animal.
  res$c2 |> print()
  res$c2_spread |> print()

  ## 5. C3 -- which side of the reversal, and it can break a written sentence.
  ## READ `c3_call` AND `c3_sign_stable` FIRST. An interval spanning zero places
  ## these tumours on NEITHER side, and a sign that flips between the pooled,
  ## group-adjusted and composition-adjusted readings is the finding, not noise
  ## to be resolved by picking one.
  res$c3_call |> print(); res$c3_signs |> print(); res$c3_sign_stable |> print()
  ## `pooled` includes the construct effect; `group-adjusted` is the within-arm
  ## association. Say which one a sentence is quoting.
  res$c3 |> print()
  res$c3_reference_directions |> print()

  ## 6. THE FAILURE READING. Read the POSITIVE CONTROL row first: if Ppargc1a --
  ## a transgene that is definitionally present -- does not survive the
  ## adjustment, the adjustment is eating signal and neither result below counts.
  res$adjusted$c1 |> print()
  res$adjusted$c3 |> print()

  res$verdict |> as.data.frame() |> print(right = FALSE)
  cat(res$notes, sep = "\n")
}
