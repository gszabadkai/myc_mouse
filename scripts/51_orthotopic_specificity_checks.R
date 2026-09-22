# =============================================================================
# 51 -- ORTHOTOPIC VECTOR SERIES: SPECIFICITY, NORMALISATION, AND ONE
#       CORRECTION TO HOW C3 WAS CONSTRUCTED
# -----------------------------------------------------------------------------
# Cohort: the same 21 vector-series samples as script 50 -- EV / BclxL / Pgc1a,
# 7 each, `data/orthotopic_series/`. Read that README and
# `docs/2026-09-08_orthotopic_vector_series.md` before this file.
#
# THIS IS NOT A FOURTH CLAIM, AND THAT IS SAID FIRST. Script 50 settled three:
# C1 HOLDS, C2 HOLDS AS AN OBSERVATION, C3 IS UNINFORMATIVE. Nothing here
# reopens any of them. This script adds specificity checks on results already
# obtained, one normalisation robustness check, and one correction to how C3 was
# constructed. It qualifies C1's WORDING and corrects C3's CONSTRUCTION. No
# result here becomes a manuscript claim.
#
# C3 STAYS UNINFORMATIVE. PART D.3 and PART F are note-quality work. If a sign
# stabilises, the interval still spans zero and the verdict is unchanged. The
# abstract sentence stays withdrawn.
#
# IT CANNOT ESTIMATE THE MYC x OXPHOS INTERACTION. Every arm is MYAZ-derived and
# MYC-high; there is no MYC-low arm, so it is not identifiable here at any n. It
# lives in the iMMEC rtTA-MYC +/-dox x +/-PGC1a design on the death readout. No
# group term is allowed to stand in for it.
#
# N3. These are transcript associations. The word "primed" is not written of a
# transcript anywhere in this script, its comments, its labels or its outputs.
#
# =============================================================================
# GATE 0, ANSWERED BEFORE THIS FILE WAS WRITTEN
# -----------------------------------------------------------------------------
# `/Users/gs/G/data/MK_myc_2022` was searched exhaustively, no depth limit, for
# *length_scaled*, *transcript_counts*, *gene_lengths*, quant.sf and tx2gene*.
# NONE EXIST. The only salmon artefact anywhere is salmon.merged.gene_counts.tsv
# in three byte-identical copies, MD5 fba3f8fb41de4d4a2c2b76091f933c69:
#   MK_paper_2026/OneDrive_1_15-03-2026/Bulkseq/  (the delivery)
#   orth_tumour_data/                             (working copy)
#   myc_mouse/data/orthotopic_series/             (this repo's snapshot)
#
# STRONGER THAN "NOT FOUND": the delivery folder is a RESULTS folder, not a
# pipeline tree -- the gene matrix plus thirteen downstream MSigDB CSVs and two
# empty directories. No star_salmon/, no per-sample directories, no .rds. There
# is no truncated pipeline output to return to; THESE FILES WERE NEVER
# DELIVERED. Re-obtaining them means going back to whoever ran the pipeline.
#
# TWO CONSEQUENCES, BOTH REAL:
#
#  (1) THE LENGTH OFFSET IS LOST. tximport is unavailable, counts go through
#      DESeqDataSetFromMatrix, and gene-level counts without the offset misstate
#      genes whose isoform usage differs between arms. Applies to every number
#      this script produces, as it did to 50.
#
#  (2) Bcl-xL AND Bcl-xS CANNOT BE SEPARATED, and for this dataset that is the
#      SHARPER limitation. They arise from alternative 5' splice site usage at
#      the same locus and are functionally OPPOSITE -- long anti-apoptotic,
#      short pro-apoptotic. Gene-level `Bcl2l1` pools them.
#      C1, the one clean causal result, rests entirely on `Bcl2l1`. So
#      "PGC1a raises the guardian" CANNOT BE SEPARATED FROM "PGC1a shifts
#      splicing toward the short pro-apoptotic isoform" in these data, and no
#      analysis choice here can separate them. It also qualifies the BclxL arm:
#      the construct is Bcl-xL cDNA, but the readout cannot distinguish
#      construct from endogenous Bcl-xS.
#      THIS IS A LIMITATION ON C1's INTERPRETATION, NOT ON ITS MEASUREMENT, and
#      it is the one a referee finds first. It goes in the note beside the
#      length offset.
#
# =============================================================================
# READING RULES -- FIXED HERE, BEFORE ANY NUMBER IN THIS SCRIPT IS SEEN
# -----------------------------------------------------------------------------
# R1  `Mcl1` in C1. Flat, or moving far less than `Bcl2l1` -> C1 is a SELECTIVE
#     event and the guardian-dissociation sentence is available. Moving with
#     `Bcl2l1` at comparable magnitude -> C1 is a general anti-apoptotic
#     response, the specificity sentence comes out, and C1 is weaker.
#
# R2  `Bbc3` in C1. Flat or falling while `Bcl2l1` rises -> the
#     guardian-sensitiser gap widens BY INTERVENTION and C1 strengthens. Rising
#     proportionally -> the gap is unchanged and C1 is a general BCL2-family
#     response.
#
# R3  The asymmetry, BclxL vs EV. Both rulers flat -> the asymmetry holds and
#     the directed-edge argument is available. `ox_lvl` up with `ox_rel` flat ->
#     ruler-specific; wording restricted to respiratory SHARE, not respiration.
#     `ox_lvl` up -> no asymmetry and the argument comes out.
#     IN ALL THREE CASES this must be reconciled in the note against the
#     standing experimental finding that Bcl-xL overexpressing tumours have
#     higher OXPHOS -- and the note must FLAG WHETHER THAT FINDING WAS
#     TRANSCRIPTIONAL OR FUNCTIONAL, because a functional-versus-transcript
#     discordance is a different statement from an asymmetry.
#
# R4  C2 under median-of-ratios. Survives -> the ceiling is biology and C2
#     stands as an observation. Attenuates materially -> CPM deflation was
#     carrying it and C2 weakens or comes out. STATE WHICH BEFORE QUOTING THE
#     NUMBER ANYWHERE.
#
# R5  PART E. Wide `Ppargc1a` spread with flat `Bcl2l1` -> threshold reading,
#     consistent with selection. `Bcl2l1` tracking `Ppargc1a` -> dose-responsive
#     regulatory reading. NARROW `Ppargc1a` SPREAD -> UNINFORMATIVE, and say so
#     rather than reading the flatness.
#
# =============================================================================
# WHY THE COHORT IS REBUILT FROM SOURCE RATHER THAN READ FROM 50'S RDS
# -----------------------------------------------------------------------------
# PART 0 below is script 50's PART 0 and PART A COPIED VERBATIM -- the load, the
# longest-first group parsing, the two input objects, the MitoCarta split and
# the rulers. It is copied rather than refactored into a shared file because
# script 50 is run and signed off and must not be touched. The copy is then
# PROVEN against 50's published numbers before anything new is computed. If the
# copy has drifted, PART 0 stops the script and nothing below is readable.
#
# Species = cohort. No pooling with the 6W/12W timeline or with human. The gland
# +0.351 and human -0.31/-0.46 are DIRECTIONS ONLY and are never mixed with
# these cohort-relative scores.
#
# Reads : data/orthotopic_series/salmon.merged.gene_counts.tsv
#         data/Mouse.MitoCarta3.0.xls  Sheet 4  (exact filename, never a glob)
#         functions/reconcile_gene_symbols.R (MANDATORY)
# Writes: results/orthotopic_specificity_checks.rds
#         outputs/orthotopic/51_arm_geometry.pdf
#
# RUNTIME: a couple of minutes. The C2 nulls are EXACT -- all 3432 splits.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

set.seed(1)
NBOOT   <- 5000L   # bootstrap resamples for every interval
MINBM   <- 10      # smallest mean normalised count carried into a composite
CUT_MAT <- 0.5     # |rho| at which a composition axis counts as MATERIAL
VEC     <- c("EV", "BclxL", "Pgc1a")            # the scoring set, fixed
ALLGRP  <- c(VEC, "NTEV", "KOEV", "NTPgc1a", "KOPgc1a", "MYAZ", "FaMY", "BoMY")

# --- the load control, from the analysis plan's section 2 --------------------
REF <- data.frame(
  group    = ALLGRP,
  n        = c(7L, 7L, 7L, 7L, 7L, 7L, 7L, 6L, 6L, 6L),
  Bcl2l1   = c(71, 326, 193, 82, 62, 79, 79, 59, 59, 51),
  Ppargc1a = c(0.0, 0.3, 302, 0.2, 0.1, 129, 0.2, 0.0, 0.0, 0.0),
  Myc      = c(1915, 2033, 1440, 2208, 2368, 2286, 2176, 1370, 1280, 1293),
  stringsAsFactors = FALSE)
CTRL_TOL_REL <- 0.02      # 2 pct, the rounding in the reference table

# --- script 50's C1 table, transcribed from docs/2026-09-08_..._series.md -----
# POINT ESTIMATES ONLY. The bootstrap intervals are NOT asserted: they depend on
# the RNG position in the stream, and 51 does not consume randomness in the same
# order as 50. The point estimates and the Wilcoxon p are deterministic, so they
# are what the control tests. The intervals are recomputed and REPORTED beside
# the published ones for eyeballing, never asserted.
REF_C1 <- data.frame(
  measure   = c("Bcl2l1 (CPM)", "Ppargc1a (CPM, manipulation check)",
                "ox_lvl (OXPHOS level)", "ox_rel (respiratory share)",
                "Bcl2l1 (CPM) -- BclxL arm, construct present"),
  hi        = c("Pgc1a", "Pgc1a", "Pgc1a", "Pgc1a", "BclxL"),
  ref_diff  = c(121.8, 301.6, 1.925, 0.915, 255.4),
  ref_ci_lo = c(98.1, 268.6, 1.663, 0.761, 181.9),
  ref_ci_hi = c(151.2, 356.4, 2.296, 1.007, 292.7),
  stringsAsFactors = FALSE)
C1_TOL_REL <- 0.01        # 1 pct on the point estimate
P_FLOOR    <- 1 / 3432    # the one-sided floor for 7 v 7; 50 reports 0.00029

NOTES <- character(0)
note <- function(...) NOTES <<- c(NOTES, paste0(...))

# =============================================================================
# PART 0: LOAD, REBUILD THE COHORT, AND PROVE IT REPRODUCES SCRIPT 50
# -----------------------------------------------------------------------------
# Copied verbatim from 50 PART 0 + PART A, then proven. Nothing new is computed
# until every assertion in this part has passed.
# =============================================================================
message("51 PART 0: load and reproduce script 50")

raw <- data.table::fread(
  here::here("data", "orthotopic_series", "salmon.merged.gene_counts.tsv"),
  data.table = FALSE, check.names = FALSE)
stopifnot(identical(colnames(raw)[1:2], c("gene_id", "gene_name")))
CNT <- as.matrix(raw[, -(1:2)])
rownames(CNT) <- raw$gene_id
stopifnot(ncol(CNT) == 67L, all(grepl("^ENSMUSG", rownames(CNT))))

# GROUP PARSING, longest-first. `EV` and `Pgc` must be matched AFTER
# `NTEV`/`KOEV`/`NTPgc`/`KOPgc` or four groups silently collapse into two.
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
  if (length(r) == 1L) CPM[r, ] else colSums(CPM[r, , drop = FALSE]) }
load_control <- data.frame(
  group = ALLGRP, n = as.integer(table(grp)),
  Bcl2l1   = round(as.numeric(tapply(cpm_of("Bcl2l1"),   grp, stats::median)), 1),
  Ppargc1a = round(as.numeric(tapply(cpm_of("Ppargc1a"), grp, stats::median)), 1),
  Myc      = round(as.numeric(tapply(cpm_of("Myc"),      grp, stats::median)), 1),
  stringsAsFactors = FALSE)
dev <- max(abs(unlist(load_control[, 3:5]) - unlist(REF[, 3:5])) /
           pmax(unlist(REF[, 3:5]), 1))
message(sprintf("51 PART 0: median-CPM control, max relative deviation %.3f%% (tol %.1f%%)",
                100 * dev, 100 * CTRL_TOL_REL))
stopifnot(dev < CTRL_TOL_REL)

# --- the cohort, the two input objects, and the rulers (50 PART A) -----------
# TWO OBJECTS THAT NEVER MIX: `NCv` is LINEAR normalised counts; `L` is
# log2(NCv + 1) and feeds every z-composite. mitoPPS and GSVA are NOT run.
vk  <- grp %in% VEC
sm  <- data.frame(sample = colnames(CNT)[vk],
                  group  = factor(as.character(grp[vk]), levels = VEC))
stopifnot(nrow(sm) == 21L, identical(as.integer(table(sm$group)), c(7L, 7L, 7L)))

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(CNT[, vk, drop = FALSE]), colData = sm, design = ~ group)
dds <- DESeq2::estimateSizeFactors(dds)
NCv <- DESeq2::counts(dds, normalized = TRUE)          # LINEAR
keep <- rowMeans(NCv) >= MINBM & matrixStats::rowSds(NCv) > 0
L    <- log2(NCv[keep, , drop = FALSE] + 1)            # LOG
message(sprintf("51 PART 0: %d of %d genes at mean normalised count >= %d",
                nrow(L), nrow(NCv), MINBM))

universe_all <- rownames(L)
zrow    <- function(m) t(scale(t(m)))
comp_e  <- function(e) { e <- e[!is.na(e) & e %in% rownames(L)]
                         stopifnot(length(e) >= 3); colMeans(zrow(L[e, , drop = FALSE])) }
ens_set <- function(syms) { e <- recon_to_ensembl(syms, universe_all); e[!is.na(e)] }
one_ens <- function(sym) { e <- ens_set(sym)
                           if (!length(e)) NA_character_ else e[1] }
lg <- function(sym) { e <- one_ens(sym); stopifnot(!is.na(e)); as.numeric(L[e, ]) }

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
stopifnot(length(intersect(ox_sub, mt_ens)) == 0L, length(mt_ens) >= 10L)

set_sizes <- tibble::tibble(
  set = c("OXPHOS subunits (nuclear)", "rest of nuclear MitoCarta",
          "mtDNA-encoded", "genes carried"),
  n   = c(length(ox_sub), length(rest_mito), length(mt_ens), nrow(L)))
message(sprintf("51 PART 0: ox_sub %d | rest %d | mtDNA %d",
                length(ox_sub), length(rest_mito), length(mt_ens)))
# 50 reported 87 / 880 / 13 with 15721 genes carried. Asserted, because a
# different ruler would silently change every score below.
stopifnot(length(ox_sub) == 87L, length(rest_mito) == 880L,
          length(mt_ens) == 13L, nrow(L) == 15721L)

X <- list(ox_rel = comp_e(ox_sub) - comp_e(rest_mito),
          ox_lvl = comp_e(ox_sub),
          ox_mt  = comp_e(mt_ens))
# the two halves of ox_rel kept separately -- PART C decomposes the move
X$ox_den <- comp_e(rest_mito)

G <- list(Bcl2l1 = lg("Bcl2l1"), Mcl1 = lg("Mcl1"), Bbc3 = lg("Bbc3"),
          Bcl2l11 = lg("Bcl2l11"), Myc = lg("Myc"), Ppargc1a = lg("Ppargc1a"))
G$guardian_ratio <- G$Mcl1 - G$Bcl2l1          # log2(Mcl1/Bcl2l1)
cpm_v <- lapply(c(Bcl2l1 = "Bcl2l1", Ppargc1a = "Ppargc1a", Myc = "Myc"),
                function(s) cpm_of(s)[vk])

# --- the reproduction control against script 50's published C1 ---------------
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

c1_repro <- dplyr::bind_rows(
  dplyr::bind_cols(tibble::tibble(measure = REF_C1$measure[1]), two_group(cpm_v$Bcl2l1,   "Pgc1a", "EV")),
  dplyr::bind_cols(tibble::tibble(measure = REF_C1$measure[2]), two_group(cpm_v$Ppargc1a, "Pgc1a", "EV")),
  dplyr::bind_cols(tibble::tibble(measure = REF_C1$measure[3]), two_group(X$ox_lvl,       "Pgc1a", "EV")),
  dplyr::bind_cols(tibble::tibble(measure = REF_C1$measure[4]), two_group(X$ox_rel,       "Pgc1a", "EV")),
  dplyr::bind_cols(tibble::tibble(measure = REF_C1$measure[5]), two_group(cpm_v$Bcl2l1,   "BclxL", "EV")))
c1_repro <- dplyr::left_join(c1_repro, REF_C1[, c("measure", "ref_diff", "ref_ci_lo", "ref_ci_hi")],
                             by = "measure")
c1_repro$rel_dev <- abs(c1_repro$diff_median - c1_repro$ref_diff) / abs(c1_repro$ref_diff)
message(sprintf("51 PART 0: C1 reproduction, max relative deviation %.4f%% (tol %.1f%%)",
                100 * max(c1_repro$rel_dev), 100 * C1_TOL_REL))
stopifnot(max(c1_repro$rel_dev) < C1_TOL_REL,
          all(c1_repro$wilcox_p <= P_FLOOR * 1.001))
note("PART 0: the rebuilt cohort reproduces script 50's five C1 point estimates ",
     sprintf("to %.4f%% and all five sit at the 7v7 one-sided floor p = %.5f. ",
             100 * max(c1_repro$rel_dev), P_FLOOR),
     "Bootstrap intervals are reported beside the published ones but NOT ",
     "asserted: they depend on RNG position and 51 does not consume randomness ",
     "in 50's order.")

load_control_out <- list(median_cpm = load_control, reference = REF,
                         max_rel_dev = dev, c1 = c1_repro, set_sizes = set_sizes)

# =============================================================================
# PART A: ARM GEOMETRY -- IS THE NEAR-ZERO POOLED C3 rho A CANCELLING TRIANGLE?
# -----------------------------------------------------------------------------
# ONE specific explanation, confirmed or refuted: the three arm-means form a
# NON-MONOTONE triangle in the (ox_rel, Bcl2l1) plane -- EV low-low, Pgc1a
# high-mid, BclxL low-high -- so the EV->Pgc1a limb and the BclxL arm cancel.
#
# The pooled rho is decomposed into a BETWEEN-arm part (the three arm-means) and
# a WITHIN-arm part (the mean of the three within-arm rhos), which turns "the
# arms cancel" from a picture into a number. The between-arm rho is on n = 3 and
# is DESCRIPTIVE ONLY -- it has no interval and none is invented.
# =============================================================================
message("51 PART A: arm geometry")

sp <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman"))

geom_vars <- list(ox_rel = X$ox_rel, ox_lvl = X$ox_lvl,
                  Bcl2l1 = G$Bcl2l1, Mcl1 = G$Mcl1,
                  guardian_ratio = G$guardian_ratio)
arm_means <- dplyr::bind_rows(lapply(VEC, function(g) {
  k <- sm$group == g
  tibble::tibble(group = g, n = sum(k),
                 !!!lapply(geom_vars, function(v) mean(v[k]))) }))
arm_sds <- dplyr::bind_rows(lapply(VEC, function(g) {
  k <- sm$group == g
  tibble::tibble(group = g,
                 !!!lapply(geom_vars, function(v) stats::sd(v[k]))) }))

# is the (ox_rel, endpoint) arm-mean sequence monotone once ordered by ox_rel?
mono_of <- function(ep) {
  o <- order(arm_means$ox_rel)
  v <- arm_means[[ep]][o]
  list(order_by_ox_rel = paste(arm_means$group[o], collapse = " < "),
       values = round(v, 3),
       monotone = all(diff(v) > 0) || all(diff(v) < 0)) }

decomp <- dplyr::bind_rows(lapply(c("Bcl2l1", "Mcl1", "guardian_ratio"), function(ep) {
  within <- vapply(VEC, function(g) { k <- sm$group == g
                                      sp(geom_vars[[ep]][k], X$ox_rel[k]) }, 0)
  m <- mono_of(ep)
  tibble::tibble(endpoint = ep,
                 rho_pooled_21 = sp(geom_vars[[ep]], X$ox_rel),
                 rho_between_arms_n3 = sp(arm_means$ox_rel, arm_means[[ep]]),
                 rho_within_EV = within[["EV"]], rho_within_BclxL = within[["BclxL"]],
                 rho_within_Pgc1a = within[["Pgc1a"]],
                 rho_within_mean = mean(within),
                 arm_order_by_ox_rel = m$order_by_ox_rel,
                 arm_means_monotone = m$monotone) }))

arm_geometry <- list(arm_means = arm_means, arm_sds = arm_sds,
                     decomposition = decomp,
                     per_sample = dplyr::bind_cols(
                       sm, tibble::as_tibble(geom_vars)))
tri <- decomp$arm_means_monotone[decomp$endpoint == "Bcl2l1"]
note("PART A: the (ox_rel, Bcl2l1) arm-means are ",
     if (tri) { "MONOTONE -- the cancelling-triangle explanation is REFUTED."
     } else { "NON-MONOTONE -- the cancelling-triangle explanation is CONFIRMED." },
     " Between-arm and within-arm rhos are reported separately in $decomposition.")

# --- the figure --------------------------------------------------------------
dir.create(here::here("outputs", "orthotopic"), showWarnings = FALSE, recursive = TRUE)
pg <- arm_geometry$per_sample
pl <- function(ep, ttl) {
  ggplot2::ggplot(pg, ggplot2::aes(x = ox_rel, y = .data[[ep]], colour = group)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::geom_point(data = arm_means, size = 5, shape = 17) +
    ggplot2::geom_path(data = arm_means[order(arm_means$ox_rel), ],
                       ggplot2::aes(group = 1), colour = "grey40",
                       linetype = "dashed", inherit.aes = TRUE) +
    ggplot2::labs(title = ttl, x = "ox_rel (respiratory share)", y = ep) +
    ggplot2::theme_bw() }
p_all <- patchwork::wrap_plots(
  pl("Bcl2l1", "Bcl2l1 vs ox_rel: triangles are arm means"),
  pl("guardian_ratio", "log2(Mcl1/Bcl2l1) vs ox_rel"), ncol = 2) +
  patchwork::plot_annotation(
    title = "51 PART A: arm geometry in the vector series (n = 21)",
    subtitle = "transcript associations; arm means joined in ox_rel order")
ggplot2::ggsave(here::here("outputs", "orthotopic", "51_arm_geometry.pdf"),
                p_all, width = 11, height = 5)
message("51 PART A: wrote outputs/orthotopic/51_arm_geometry.pdf")

# =============================================================================
# PART B: THE C1 PANEL EXTENDED -- SPECIFICITY (R1, R2)
# -----------------------------------------------------------------------------
# Script 50's machinery, Pgc1a vs EV, 7 v 7, one-sided by pre-specification,
# applied to Mcl1, Bbc3 and the Mcl1:Bcl2l1 group contrast that 50 never ran,
# with Bcl2l1 repeated as the ANCHOR. BclxL is reported as the third arm and
# NEVER pooled, exactly as 50 did.
#
# Bcl2a1 is a NON-ONE-TO-ONE MOUSE PARALOGUE. Bcl2a1b and Bcl2a1d are reported
# SEPARATELY and are never summed into a single "Bcl2a1".
# NO PER-GENE FDR across the panel. Estimates with intervals.
# =============================================================================
message("51 PART B: C1 extended")

FAMILY <- c("Bcl2l1", "Mcl1", "Bbc3", "Bad", "Bcl2", "Bcl2a1b", "Bcl2a1d",
            "Bcl2l11", "Bcl2l2", "Bid", "Bik", "Bmf", "Pmaip1")
fam_e  <- vapply(FAMILY, one_ens, "")
fam_ok <- FAMILY[!is.na(fam_e) & fam_e != ""]
fam_missing <- setdiff(FAMILY, fam_ok)
if (length(fam_missing))
  message(sprintf("51 PART B: not detected, dropped and named: %s",
                  paste(fam_missing, collapse = ", ")))

# the three the reading rules turn on, on the log2 scale the ratio needs
c1_panel <- dplyr::bind_rows(
  dplyr::bind_cols(tibble::tibble(measure = "Bcl2l1 (log2 norm) -- ANCHOR", role = "anchor"),
                   two_group(G$Bcl2l1, "Pgc1a", "EV")),
  dplyr::bind_cols(tibble::tibble(measure = "Mcl1 (log2 norm) -- R1", role = "rule"),
                   two_group(G$Mcl1, "Pgc1a", "EV")),
  dplyr::bind_cols(tibble::tibble(measure = "Bbc3 (log2 norm) -- R2", role = "rule"),
                   two_group(G$Bbc3, "Pgc1a", "EV")),
  # the group contrast on the ratio, which 50 never ran. TWO-SIDED: no direction
  # was pre-specified for it, and inventing one now would be the grid trap.
  dplyr::bind_cols(tibble::tibble(measure = "guardian_ratio log2(Mcl1/Bcl2l1) -- two-sided", role = "ratio"),
                   two_group(G$guardian_ratio, "Pgc1a", "EV", one_sided = FALSE)),
  # the third arm, reported and never pooled
  dplyr::bind_cols(tibble::tibble(measure = "Bcl2l1 (log2 norm) -- BclxL arm", role = "third arm"),
                   two_group(G$Bcl2l1, "BclxL", "EV")),
  dplyr::bind_cols(tibble::tibble(measure = "Mcl1 (log2 norm) -- BclxL arm", role = "third arm"),
                   two_group(G$Mcl1, "BclxL", "EV")))

# the whole roster alongside, as a panel and not as a set of tests
fam_panel <- dplyr::bind_rows(lapply(fam_ok, function(s) {
  v <- as.numeric(L[fam_e[[s]], ])
  dplyr::bind_rows(
    dplyr::bind_cols(tibble::tibble(gene = s, contrast = "Pgc1a vs EV"),
                     two_group(v, "Pgc1a", "EV", one_sided = FALSE)),
    dplyr::bind_cols(tibble::tibble(gene = s, contrast = "BclxL vs EV"),
                     two_group(v, "BclxL", "EV", one_sided = FALSE))) }))

# R1 and R2 applied mechanically, on the log2 scale so magnitudes are comparable
d_bcl <- c1_panel$diff_median[c1_panel$measure == "Bcl2l1 (log2 norm) -- ANCHOR"]
d_mcl <- c1_panel$diff_median[c1_panel$measure == "Mcl1 (log2 norm) -- R1"]
d_bbc <- c1_panel$diff_median[c1_panel$measure == "Bbc3 (log2 norm) -- R2"]
r1 <- {
  if (abs(d_mcl) < 0.5 * abs(d_bcl)) {
    "SELECTIVE -- Mcl1 moves less than half as far as Bcl2l1; the guardian-dissociation sentence is available"
  } else {
    "GENERAL -- Mcl1 moves comparably to Bcl2l1; the specificity sentence comes out and C1 is weaker"
  } }
r2 <- {
  if (d_bbc <= 0) {
    "GAP WIDENS -- Bbc3 flat or falling while Bcl2l1 rises; C1 strengthens"
  } else if (d_bbc < 0.5 * d_bcl) {
    "GAP WIDENS (partially) -- Bbc3 rises but far less than Bcl2l1"
  } else {
    "GENERAL -- Bbc3 rises proportionally; the gap is unchanged"
  } }
note("R1: ", r1, ". R2: ", r2, ".")

# =============================================================================
# PART C: THE ASYMMETRY CHECK (R3), AND WHAT THE PGC1a MOVE IS MADE OF
# -----------------------------------------------------------------------------
# BclxL vs EV on both rulers. Plus the PGC1a ox_lvl move DECOMPOSED into its
# numerator (nuclear OXPHOS subunits) and its denominator (the rest of nuclear
# MitoCarta), so how OXPHOS-selective the intervention was is VISIBLE rather
# than inferred from the gap between +1.925 and +0.915.
# =============================================================================
message("51 PART C: asymmetry and decomposition")

asym <- dplyr::bind_rows(
  dplyr::bind_cols(tibble::tibble(measure = "ox_lvl", arm = "BclxL vs EV"),
                   two_group(X$ox_lvl, "BclxL", "EV", one_sided = FALSE)),
  dplyr::bind_cols(tibble::tibble(measure = "ox_rel", arm = "BclxL vs EV"),
                   two_group(X$ox_rel, "BclxL", "EV", one_sided = FALSE)),
  dplyr::bind_cols(tibble::tibble(measure = "ox_mt (mtDNA-encoded)", arm = "BclxL vs EV"),
                   two_group(X$ox_mt, "BclxL", "EV", one_sided = FALSE)))

# the numerator/denominator decomposition, both arms so they can be compared
decomp_ruler <- dplyr::bind_rows(lapply(c("Pgc1a", "BclxL"), function(g)
  dplyr::bind_rows(
    dplyr::bind_cols(tibble::tibble(component = "numerator: OXPHOS subunits", arm = g),
                     two_group(X$ox_lvl, g, "EV", one_sided = FALSE)),
    dplyr::bind_cols(tibble::tibble(component = "denominator: rest of nuclear MitoCarta", arm = g),
                     two_group(X$ox_den, g, "EV", one_sided = FALSE)),
    dplyr::bind_cols(tibble::tibble(component = "ox_rel = numerator - denominator", arm = g),
                     two_group(X$ox_rel, g, "EV", one_sided = FALSE)))))

lvl_b <- asym$diff_median[asym$measure == "ox_lvl"]
rel_b <- asym$diff_median[asym$measure == "ox_rel"]
lo_b  <- asym$ci_lo[asym$measure == "ox_lvl"]; hi_b <- asym$ci_hi[asym$measure == "ox_lvl"]
lo_r  <- asym$ci_lo[asym$measure == "ox_rel"]; hi_r <- asym$ci_hi[asym$measure == "ox_rel"]
r3 <- {
  if (lo_b <= 0 && hi_b >= 0 && lo_r <= 0 && hi_r >= 0) {
    "ASYMMETRY HOLDS -- both rulers flat in BclxL vs EV; the directed-edge argument is available"
  } else if ((lo_b > 0 || hi_b < 0) && lo_r <= 0 && hi_r >= 0) {
    "RULER-SPECIFIC -- ox_lvl moves, ox_rel flat; wording restricted to respiratory SHARE"
  } else {
    "NO ASYMMETRY -- ox_lvl moves in BclxL vs EV; the directed-edge argument comes out"
  } }
note("R3: ", r3, ". THE NOTE MUST RECONCILE THIS against the standing ",
     "experimental finding that Bcl-xL tumours have higher OXPHOS, and must ",
     "flag whether that finding was TRANSCRIPTIONAL or FUNCTIONAL -- a ",
     "functional-versus-transcript discordance is a different statement from ",
     "an asymmetry.")

# =============================================================================
# PART D: NORMALISATION
# -----------------------------------------------------------------------------
# D.1 diagnostic -- size factors against library size, and the MitoCarta share
#     of each library, both summarised BY ARM. This MEASURES the CPM deflation
#     instead of arguing about it.
# D.2 C2 recomputed on median-of-ratios normalised counts instead of CPM, same
#     EXACT permutation null over all 3432 relabellings. This is R4.
# D.3 the composition-adjusted C3 fit under both normalisations.
#
#     A CORRECTION TO THE BRIEF, STATED RATHER THAN SILENTLY APPLIED. The brief
#     asked for 50's PART F C3 fit "recomputed on median-of-ratios rather than
#     CPM". Script 50's PART F ALREADY RUNS ON MEDIAN-OF-RATIOS: its endpoints
#     are lg() = log2(counts(dds, normalized = TRUE) + 1). The direction in the
#     brief is inverted. The informative check is therefore the MIRROR -- the
#     same fit on CPM -- and BOTH are computed here so the comparison exists
#     whichever way round the question was meant.
# =============================================================================
message("51 PART D: normalisation")

sf  <- DESeq2::sizeFactors(dds)
libM <- colSums(CNT[, vk, drop = FALSE]) / 1e6
# MitoCarta share of each library, on RAW counts and over the UNFILTERED matrix,
# so the share is of everything sequenced and not of what survived MINBM
uni_cnt   <- rownames(CNT)
ens_cnt   <- function(syms) { e <- recon_to_ensembl(syms, uni_cnt); e[!is.na(e)] }
mito_all_cnt <- unique(c(ens_cnt(unique(g2p_nuc$Gene)), ens_cnt(mt_syms)))
mito_share   <- colSums(CNT[intersect(mito_all_cnt, uni_cnt), vk, drop = FALSE]) /
                colSums(CNT[, vk, drop = FALSE])

norm_diagnostic <- list(
  per_sample = tibble::tibble(sample = sm$sample, group = sm$group,
                              size_factor = as.numeric(sf),
                              lib_millions = as.numeric(libM),
                              sf_over_libM = as.numeric(sf) / as.numeric(libM),
                              mitocarta_share = as.numeric(mito_share)),
  by_arm = NULL)
norm_diagnostic$by_arm <- norm_diagnostic$per_sample |>
  dplyr::group_by(group) |>
  dplyr::summarise(n = dplyr::n(),
                   sf_median = stats::median(size_factor),
                   libM_median = stats::median(lib_millions),
                   sf_over_libM_median = stats::median(sf_over_libM),
                   mitocarta_share_median = stats::median(mitocarta_share),
                   .groups = "drop")

# --- D.2: C2 under median-of-ratios ------------------------------------------
nc_of <- function(sym) {
  r <- which(raw$gene_name == sym)
  stopifnot(length(r) >= 1L)
  m <- NCv[rownames(CNT)[r], , drop = FALSE]
  if (nrow(m) == 1L) as.numeric(m) else colSums(m) }
myc_mor <- nc_of("Myc")
stopifnot(length(myc_mor) == 21L)

exact_null <- function(v, hi, lo, stat) {
  a <- v[sm$group == hi]; b <- v[sm$group == lo]
  pool <- c(a, b); n <- length(a)
  idx  <- utils::combn(length(pool), n)
  null <- apply(idx, 2, function(k) stat(pool[k]) - stat(pool[-k]))
  obs  <- stat(a) - stat(b)
  list(obs = obs, p_lower = mean(null <= obs), n_perm = ncol(idx),
       null_median = stats::median(null)) }
p80 <- function(z) unname(stats::quantile(z, 0.8, type = 7))

c2_one <- function(v, scale_lab) dplyr::bind_rows(lapply(
  list(c("Pgc1a", "EV"), c("Pgc1a", "BclxL")), function(pr)
    dplyr::bind_rows(lapply(list(max = max, pct80 = p80), function(st) {
      r <- exact_null(v, pr[1], pr[2], st)
      tibble::tibble(scale = scale_lab, hi = pr[1], lo = pr[2], stat_obs = r$obs,
                     p_one_sided_lower = r$p_lower, n_perm = r$n_perm,
                     null_median = r$null_median) }), .id = "statistic")))
c2_mor <- dplyr::bind_rows(c2_one(cpm_v$Myc, "CPM (script 50)"),
                           c2_one(myc_mor,   "median-of-ratios"))

cmp <- c2_mor |>
  dplyr::select(statistic, hi, lo, scale, p_one_sided_lower) |>
  tidyr::pivot_wider(names_from = scale, values_from = p_one_sided_lower)
key <- cmp$`median-of-ratios`[cmp$statistic == "pct80" & cmp$lo == "EV"]
r4 <- {
  if (key <= 0.01) {
    "SURVIVES -- the 80th-percentile contrast against EV keeps p <= 0.01 under median-of-ratios; the ceiling is not CPM deflation and C2 stands as an observation"
  } else if (key <= 0.05) {
    "WEAKENS -- p rises above 0.01 but stays under 0.05; C2 stands as an observation with the normalisation dependence stated"
  } else {
    "ATTENUATES MATERIALLY -- CPM deflation was carrying C2; it weakens or comes out"
  } }
note("R4: ", r4, sprintf(" (80th-pct vs EV, median-of-ratios p = %.4f).", key))

# --- D.3: the composition gate and the adjusted C3 under both normalisations --
MARK <- list(
  epithelial    = c("Krt8", "Krt18", "Epcam", "Cdh1", "Krt5", "Krt14"),
  stromal       = c("Col1a1", "Pdgfrb", "Acta2", "Thy1"),
  adipose       = c("Adipoq", "Plin1", "Fabp4", "Cidec"),
  immune        = c("Ptprc", "Cd52", "Lyz2", "Cd74"),
  endothelial   = c("Pecam1", "Cdh5", "Cldn5", "Emcn", "Tek"),
  proliferation = c("Mki67", "Top2a", "Ccnb1"))
MARK_e <- lapply(MARK, function(v) { e <- ens_set(v); e[e %in% rownames(L)] })
comp_present <- vapply(MARK_e, length, 0L)
K <- lapply(MARK_e[comp_present >= 3L], comp_e)
material <- names(K)[vapply(K, function(k) abs(sp(X$ox_rel, k)) >= CUT_MAT, TRUE)]
message(sprintf("51 PART D: material composition axes: %s",
                if (length(material)) paste(material, collapse = ", ") else "NONE"))

AD <- if (length(material)) as.data.frame(K[material]) else NULL
adj_assoc <- function(y, x, lab, scale_lab) {
  d <- data.frame(y = y, x = x)
  if (!is.null(AD)) d <- cbind(d, AD)
  s <- summary(stats::lm(y ~ ., d))$coefficients
  tibble::tibble(endpoint = lab, scale = scale_lab,
                 beta_ox_rel = s["x", 1], p = s["x", 4]) }

# CPM counterparts of the three endpoints, for the mirror fit
cpm_l2 <- function(sym) log2(cpm_of(sym)[vk] + 1)
EP_MOR <- list(guardian_ratio = G$guardian_ratio, Mcl1 = G$Mcl1, Bcl2l1 = G$Bcl2l1)
EP_CPM <- list(guardian_ratio = cpm_l2("Mcl1") - cpm_l2("Bcl2l1"),
               Mcl1 = cpm_l2("Mcl1"), Bcl2l1 = cpm_l2("Bcl2l1"))
c3_partf_mor <- dplyr::bind_rows(
  dplyr::bind_rows(lapply(names(EP_MOR), function(en)
    adj_assoc(EP_MOR[[en]], X$ox_rel, en, "median-of-ratios (script 50)"))),
  dplyr::bind_rows(lapply(names(EP_CPM), function(en)
    adj_assoc(EP_CPM[[en]], X$ox_rel, en, "CPM (the mirror)"))))
note("PART D.3: script 50's PART F ALREADY ran on median-of-ratios -- its ",
     "endpoints are log2(counts(dds, normalized = TRUE) + 1). The brief's ",
     "direction was inverted, so BOTH normalisations are fitted here and the ",
     "comparison exists whichever way the question was meant.")

# =============================================================================
# PART E: THRESHOLD OR DOSE-RESPONSE (R5)
# -----------------------------------------------------------------------------
# WITHIN the Pgc1a arm only, n = 7. Does Bcl2l1 track Ppargc1a? At n = 7 the
# interval is very wide and the SPREAD of Ppargc1a decides whether the question
# is answerable at all -- a narrow spread makes flatness uninformative rather
# than evidence of a threshold.
# =============================================================================
message("51 PART E: dose-response within the Pgc1a arm")

pk  <- sm$group == "Pgc1a"
pg_ppar <- cpm_v$Ppargc1a[pk]; pg_bcl <- cpm_v$Bcl2l1[pk]
rho_dose <- sp(pg_ppar, pg_bcl)
bs_dose <- vapply(seq_len(NBOOT), function(i) {
  k <- sample.int(length(pg_ppar), replace = TRUE)
  if (length(unique(pg_ppar[k])) < 3L) return(NA_real_)
  sp(pg_ppar[k], pg_bcl[k]) }, 0)
dose <- list(
  per_sample = tibble::tibble(sample = sm$sample[pk], Ppargc1a_CPM = pg_ppar,
                              Bcl2l1_CPM = pg_bcl),
  spread = tibble::tibble(
    measure = c("Ppargc1a", "Bcl2l1"),
    min = c(min(pg_ppar), min(pg_bcl)), max = c(max(pg_ppar), max(pg_bcl)),
    fold_spread = c(max(pg_ppar) / min(pg_ppar), max(pg_bcl) / min(pg_bcl)),
    cv = c(stats::sd(pg_ppar) / mean(pg_ppar), stats::sd(pg_bcl) / mean(pg_bcl))),
  rho = tibble::tibble(rho = rho_dose,
                       ci_lo = unname(stats::quantile(bs_dose, 0.025, na.rm = TRUE)),
                       ci_hi = unname(stats::quantile(bs_dose, 0.975, na.rm = TRUE)),
                       n = sum(pk)))
fold_ppar <- dose$spread$fold_spread[1]
r5 <- {
  if (fold_ppar < 2) {
    "UNINFORMATIVE -- Ppargc1a spans under 2-fold within the arm, so flatness of Bcl2l1 cannot discriminate a threshold from a dose-response"
  } else if (dose$rho$ci_lo > 0) {
    "DOSE-RESPONSIVE -- Bcl2l1 tracks Ppargc1a within the arm"
  } else {
    "THRESHOLD -- Ppargc1a varies materially while Bcl2l1 does not track it, consistent with selection"
  } }
note("R5: ", r5, sprintf(" (Ppargc1a fold spread %.2f, rho = %+.3f [%+.3f, %+.3f], n = 7).",
                         fold_ppar, dose$rho$rho, dose$rho$ci_lo, dose$rho$ci_hi))

# --- the external dose point, HEDGED HARD ------------------------------------
# RAW CPM ONLY. Never scored, never tested, never in a claim. NTPgc1a and NTEV
# are a DIFFERENT GENETIC BACKGROUND and this does not touch the p21 series'
# out-of-scope status. It is in the deposited matrix and a reader will find it,
# which is the only reason it is recorded.
ext_groups <- c("NTEV", "NTPgc1a")
dose$external_raw_cpm <- tibble::tibble(
  group = ext_groups,
  n = as.integer(table(grp)[ext_groups]),
  Ppargc1a_CPM_median = as.numeric(tapply(cpm_of("Ppargc1a"), grp, stats::median)[ext_groups]),
  Bcl2l1_CPM_median   = as.numeric(tapply(cpm_of("Bcl2l1"),   grp, stats::median)[ext_groups]),
  status = "RAW CPM OBSERVATION ONLY -- not scored, not tested, different genetic background")

# =============================================================================
# PART F: THE C3 CONSTRUCTION CORRECTION
# -----------------------------------------------------------------------------
# DESCRIPTIVE. NO NEW VERDICT. C3 stays UNINFORMATIVE.
#
# C1's spec said BclxL is reported and NEVER POOLED. C3's spec POOLED it -- and
# seven of the 21 Bcl2l1 values are therefore set by an EXPRESSION CONSTRUCT
# rather than by the tumour. Same endpoint, opposite handling, caught only
# because C1 and C3 were pre-specified as separate claims.
#
# The pooled C3 readings are recomputed on EV + Pgc1a alone to show they RESTATE
# C1 rather than adding to it. THE SCORES ARE NOT REBUILT: ox_rel is
# cohort-relative and was fixed on the 21, so the 14 samples are SUBSET from the
# 21-cohort scores and never re-scored.
# =============================================================================
message("51 PART F: the C3 construction correction")

assoc <- function(y, x, lab, endpoint, idx = seq_along(y)) {
  y <- y[idx]; x <- x[idx]
  r  <- sp(y, x)
  bs <- vapply(seq_len(NBOOT), function(i) {
    k <- sample.int(length(y), replace = TRUE)
    if (length(unique(x[k])) < 3L) return(NA_real_)
    sp(y[k], x[k]) }, 0)
  lo <- vapply(seq_along(y), function(i) sp(y[-i], x[-i]), 0)
  tibble::tibble(endpoint = endpoint, reading = lab, n = length(y), rho = r,
                 ci_lo = unname(stats::quantile(bs, 0.025, na.rm = TRUE)),
                 ci_hi = unname(stats::quantile(bs, 0.975, na.rm = TRUE)),
                 loo_min = min(lo), loo_max = max(lo),
                 same_sign_loo = all(sign(lo) == sign(r)),
                 excludes_zero = (unname(stats::quantile(bs, 0.025, na.rm = TRUE)) > 0) ||
                                 (unname(stats::quantile(bs, 0.975, na.rm = TRUE)) < 0)) }

k14 <- which(sm$group %in% c("EV", "Pgc1a"))
c3_ev_pgc1a <- dplyr::bind_rows(lapply(names(EP_MOR), function(en) dplyr::bind_rows(
  assoc(EP_MOR[[en]], X$ox_rel, "pooled, all 21 (script 50)",  en),
  assoc(EP_MOR[[en]], X$ox_rel, "EV + Pgc1a only (n = 14)",    en, idx = k14))))
c3_construction <- tibble::tibble(
  point = c("C1 spec", "C3 spec", "consequence"),
  statement = c(
    "BclxL is reported as the third arm and NEVER pooled",
    "BclxL WAS pooled into the n = 21 correlation",
    paste0("7 of 21 Bcl2l1 values are set by an expression construct rather ",
           "than by the tumour; on EV + Pgc1a alone the reading restates C1")))
note("PART F: descriptive only. C3 REMAINS UNINFORMATIVE whatever the sign does ",
     "on 14 samples -- the interval spans zero and the abstract sentence stays ",
     "withdrawn.")

# =============================================================================
# PART G: VERDICT
# =============================================================================
message("51 PART G: verdict")

verdict <- tibble::tibble(
  rule = c("R1 Mcl1 specificity", "R2 Bbc3 gap", "R3 asymmetry",
           "R4 C2 under median-of-ratios", "R5 threshold or dose",
           "PART A cancelling triangle", "PART F C3 construction"),
  reading = c(r1, r2, r3, r4, r5,
              if (tri) { "REFUTED -- arm means are monotone in (ox_rel, Bcl2l1)"
              } else { "CONFIRMED -- arm means are non-monotone; the arms cancel" },
              "CORRECTION RECORDED -- C3 pooled the arm C1 excluded; C3 stays UNINFORMATIVE"))

res <- list(
  load_control = load_control_out,
  arm_geometry = arm_geometry,
  c1_panel = list(headline = c1_panel, family_panel = fam_panel,
                  family_missing = fam_missing),
  asymmetry = list(bclxl_vs_ev = asym, ruler_decomposition = decomp_ruler),
  norm_diagnostic = norm_diagnostic,
  c2_mor = c2_mor,
  c3_partf_mor = c3_partf_mor,
  dose = dose,
  c3_ev_pgc1a = list(readings = c3_ev_pgc1a, construction = c3_construction),
  verdict = verdict,
  material_axes = material,
  scoring_set = VEC,
  sample_table = sm,
  params = list(NBOOT = NBOOT, MINBM = MINBM, CUT_MAT = CUT_MAT,
                C1_TOL_REL = C1_TOL_REL, CTRL_TOL_REL = CTRL_TOL_REL, seed = 1),
  gate0 = list(
    transcript_quant_available = FALSE,
    searched = "/Users/gs/G/data/MK_myc_2022, no depth limit",
    targets = c("*length_scaled*", "*transcript_counts*", "*gene_lengths*",
                "quant.sf", "tx2gene*"),
    found = "salmon.merged.gene_counts.tsv only, 3 byte-identical copies",
    md5 = "fba3f8fb41de4d4a2c2b76091f933c69",
    consequence_1 = "transcript-length offset lost; applies to every number here",
    consequence_2 = paste0("Bcl-xL and Bcl-xS cannot be separated. They arise ",
                           "from alternative 5' splice site usage at the same ",
                           "locus and are functionally opposite; gene-level ",
                           "Bcl2l1 pools them. C1 rests on Bcl2l1, so 'PGC1a ",
                           "raises the guardian' cannot be separated from ",
                           "'PGC1a shifts splicing toward the short isoform'. ",
                           "A limitation on C1's INTERPRETATION, not its ",
                           "measurement.")),
  analysis_date = Sys.Date(),
  notes = NOTES)

saveRDS(res, here::here("results", "orthotopic_specificity_checks.rds"))
message("51: wrote results/orthotopic_specificity_checks.rds")

cat("\n================ VERDICT ================\n")
print(as.data.frame(verdict), right = FALSE)
cat("\n---------------- NOTES ------------------\n")
cat(paste0("- ", NOTES, collapse = "\n"), "\n")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  res <- readRDS(here::here("results", "orthotopic_specificity_checks.rds"))

  ## 1. THE GATE. If the C1 reproduction did not pass, this script did not run
  ## at all -- PART 0 stops it. Read the deviation and the set sizes anyway.
  res$load_control$c1 |> print()
  res$load_control$max_rel_dev
  res$load_control$set_sizes |> print()

  ## 1b. GATE 0, and the limitation that matters more than the length offset.
  res$gate0$consequence_2

  ## 2. PART A -- is the near-zero pooled C3 rho a cancelling triangle? Read the
  ## decomposition, not the picture: between-arm on n = 3 is DESCRIPTIVE ONLY.
  res$arm_geometry$arm_means |> print()
  res$arm_geometry$decomposition |> print()
  ## outputs/orthotopic/51_arm_geometry.pdf

  ## 3. PART B -- R1 and R2. Mcl1 and Bbc3 against the Bcl2l1 anchor, all on the
  ## log2 scale so the magnitudes are comparable. BclxL is the third arm.
  res$c1_panel$headline |> print()
  ## the roster alongside. NO FDR -- this is a panel, not a set of tests, and
  ## Bcl2a1b / Bcl2a1d are separate paralogues that are never summed.
  res$c1_panel$family_panel |> print(n = 30)
  res$c1_panel$family_missing

  ## 4. PART C -- R3. Read ox_lvl and ox_rel together; then the decomposition,
  ## which shows whether PGC1a moved the numerator or shrank the denominator.
  res$asymmetry$bclxl_vs_ev |> print()
  res$asymmetry$ruler_decomposition |> print(n = 20)

  ## 5. PART D -- normalisation. The by-arm diagnostic first: if sf_over_libM
  ## differs by arm, CPM and median-of-ratios disagree BY ARM and that is the
  ## mechanism, not an artefact to argue about.
  res$norm_diagnostic$by_arm |> print()
  ## R4: compare the two scales side by side on the same exact null.
  res$c2_mor |> print()
  ## D.3: script 50's PART F was ALREADY median-of-ratios; the CPM row is the
  ## mirror that the brief actually needed.
  res$c3_partf_mor |> print()

  ## 6. PART E -- R5. Read the SPREAD first. A narrow Ppargc1a spread makes the
  ## flatness of Bcl2l1 uninformative rather than evidence of a threshold.
  res$dose$spread |> print()
  res$dose$rho |> print()
  res$dose$per_sample |> print()
  ## the external point: RAW CPM OBSERVATION ONLY, different background.
  res$dose$external_raw_cpm |> print()

  ## 7. PART F -- the construction correction. DESCRIPTIVE. C3 stays
  ## UNINFORMATIVE and the abstract sentence stays withdrawn.
  res$c3_ev_pgc1a$construction |> print()
  res$c3_ev_pgc1a$readings |> print()

  ## 8. Everything, on the rules fixed in the header before any number.
  res$verdict |> print()
  cat(paste0("- ", res$notes, collapse = "\n"), "\n")

}
