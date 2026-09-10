# =============================================================================
# 52 -- THE ORTHOTOPIC COHORT AS A DOSE-OF-ESCAPE SERIES
# -----------------------------------------------------------------------------
# Cohort: all 49 vector + CRISPR samples of `data/orthotopic_series/` -- EV,
# BclxL, Pgc1a_BclxL, NTEV, NTPgc1a, KOEV, KOPgc1a, 7 each. Read that README and
# `docs/2026-09-09_orthotopic_identity_correction.md` BEFORE this file.
#
# THIS IS A SELECTION EXPERIMENT, NOT A MANIPULATION EXPERIMENT, and that governs
# every sentence this script emits. PGC1a kills MYAZ cells, so a PGC1a-expressing
# line is selected before it is ever injected. Matched EV handling removes the
# protocol asymmetry but not this one. Nothing here may be worded as "PGC1a
# regulates X"; the correct form is "tumours that survived PGC1a expression show
# X". The four PGC1a arms differ in how much of the death pathway was broken and
# therefore in how much PGC1a survived:
#   transgene silenced (MYAZ + PGC1a, no RNA-seq) -> transgene LOST (KOPgc1a,
#   Ppargc1a 0.2 CPM) -> partial dose (NTPgc1a, 128.6) -> full dose under a
#   downstream buffer (Pgc1a_BclxL, 301.6).
#
# WHAT CHANGED, AND WHY 50 AND 51 ARE NOT REPRODUCED HERE. The arm those scripts
# call `Pgc1a` is PGC1a + Bcl-xL. Their C1 read `Pgc1a` vs `EV` as a single-factor
# PGC1a contrast; it is two-factor, so C1 is VOID. Their verdicts are deliberately
# NOT re-derived in PART 0 -- re-deriving them would re-enter the wrong premise.
# What PART 0 does reproduce is the raw median-CPM table, which is a property of
# the matrix and not of the premise.
#
# SCORING SET. All 49 in ONE run. `ox_rel`, GSVA and mitoPPS are cohort-relative,
# so a value from this run is NEVER quotable beside a value from the 21-sample run
# behind scripts 50 and 51, which stays untouched as the record. The prohibition
# is repeated in the saved object as `$scoring_set`.
#
# CONTRAST RULE. Clean contrasts are WITHIN SERIES -- `Pgc1a_BclxL` vs `BclxL`,
# `NTPgc1a` vs `NTEV`, `KOPgc1a` vs `KOEV`. NT and KO are CRISPR clones and the
# vector arms a polyclonal pool, so cross-series contrasts carry a derivation
# difference and are DESCRIPTIVE ONLY, labelled as such in every table. Two
# backgrounds and one manipulation is close to a replication, and is why the
# CRISPR series was returned to scope.
#
# =============================================================================
# READING RULES -- FIXED HERE, BEFORE ANY NUMBER IS COMPUTED
# -----------------------------------------------------------------------------
# R1 -- NTPgc1a's RESPIRATORY RESPONSE. The dose-response reading in the
#   correction document section 3 says NTPgc1a shows no MYC cap because
#   respiration barely rose. SMALL `ox_lvl` in NTPgc1a vs NTEV -> the reading
#   HOLDS. SUBSTANTIAL `ox_lvl` with no `Myc` effect -> the dose-response reading
#   FAILS and section 3's table must be rewritten. State which before quoting
#   anything from section 3.
#
# R2 -- `Bbc3` IN THE CONSTRUCT-FREE ARM. Flat or RISING in NTPgc1a vs NTEV ->
#   the fall in Pgc1a_BclxL requires sustained high dose under a buffer and is an
#   ESCAPE-ROUTE SIGNATURE. FALLING here too -> it tracks PGC1a at any dose and
#   the p53 explanation in section 2 is INSUFFICIENT.
#
# R3 -- THE p53 PANEL. `Cdkn2a`, `Cdkn1a`, `Mdm2`, `Ccng1`, `Zmat3` moving down
#   TOGETHER with `Bbc3` in Pgc1a_BclxL vs BclxL -> the escape-route explanation
#   is SUPPORTED. `Bbc3` moving ALONE -> it is not, and the fall needs a
#   different account.
#
# R4 -- ENDOGENOUS `Bcl2l1`. NTPgc1a vs NTEV is the construct-free test of what
#   C1 claimed, and the only place in the cohort where `Bcl2l1` means endogenous
#   BCL-XL. The medians are 79.3 vs 81.6, so a null is EXPECTED; the fitted
#   contrast and its interval are reported either way. A null here CLOSES C1
#   rather than leaving it merely void.
#
# =============================================================================
# WHAT THIS SCRIPT WILL NOT DO
# -----------------------------------------------------------------------------
# * `Bcl2l1` is UNINTERPRETABLE in the 14 construct-carrying samples (BclxL and
#   Pgc1a_BclxL). Gene-level counts cannot separate construct from endogenous
#   transcript, and BCL-XL protein is equal between those two arms by western, so
#   the transcript difference is not differential construct expression either. It
#   is reported for those arms ONLY behind an explicit uninterpretable flag.
# * Bcl-xL vs Bcl-xS is UNRESOLVABLE anywhere in this cohort -- no transcript-level
#   quantification exists on disk, so `tximport` is unavailable and the
#   transcript-length offset is lost for every number here.
# * C3, THE GUARDIAN RATIO AND THE GUARDIAN-SENSITISER GAP ARE CLOSED. They are
#   not computed, in any arm, in any form.
# * NO `MYC x OXPHOS` INTERACTION. Every arm is MYAZ-derived and MYC-high; there
#   is no MYC-low arm, so it is not identifiable at any n. It lives in the iMMEC
#   rtTA-MYC +/-dox x +/-PGC1a design on the death readout. No group term is
#   allowed to stand in for it.
# * SPECIES = COHORT. No pooling with the 6W/12W timeline or with human.
# * N3. These are transcript associations. The word "primed" is not written of a
#   transcript anywhere in this script, its comments, its labels or its outputs.
#
# Reads : data/orthotopic_series/salmon.merged.gene_counts.tsv (salmon merged,
#           gene level, NON-INTEGER estimated counts)
#         data/Mouse.MitoCarta3.0.xls  Sheet 4  (script 08's authoritative splits)
#         functions/reconcile_gene_symbols.R (MANDATORY -- vintage-aware membership)
# Writes: results/orthotopic_escape_series.rds
#
# RUNTIME: about three minutes. PART F runs an EXACT null -- all choose(14,7) =
# 3432 relabellings -- for two statistics x two scales x three contrasts.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

set.seed(1)
NBOOT   <- 5000L   # bootstrap resamples for every interval
MINBM   <- 10      # smallest mean normalised count carried into a composite
CUT_MAT <- 0.5     # |rho| at which a composition axis counts as MATERIAL
LOWCPM  <- 1       # below this median CPM a gene is reported, never fitted

VEC    <- c("EV", "BclxL", "Pgc1a_BclxL")
CRISPR <- c("NTEV", "NTPgc1a", "KOEV", "KOPgc1a")
SET    <- c(VEC, CRISPR)                       # the 49
EXCL   <- c("MYAZ", "FaMY", "BoMY")            # 18, a separate series
ALLGRP <- c("EV", "BclxL", "Pgc1a_BclxL", "NTEV", "KOEV", "NTPgc1a", "KOPgc1a", EXCL)

NOTES <- character(0)
note  <- function(...) NOTES <<- c(NOTES, paste0(...))

# --- the load control. A property of the MATRIX, not of the premise, so it is
# still the right gate even though 50's and 51's verdicts are void.
REF <- data.frame(
  group    = ALLGRP,
  n        = c(7L, 7L, 7L, 7L, 7L, 7L, 7L, 6L, 6L, 6L),
  Bcl2l1   = c(71, 326, 193, 82, 62, 79, 79, 59, 59, 51),
  Ppargc1a = c(0.0, 0.3, 302, 0.2, 0.1, 129, 0.2, 0.0, 0.0, 0.0),
  Myc      = c(1915, 2033, 1440, 2208, 2368, 2286, 2176, 1370, 1280, 1293),
  stringsAsFactors = FALSE)
CTRL_TOL_REL <- 0.02

# =============================================================================
# PART 0: LOAD, RELABEL, AND ASSERT THE COUNTS
# =============================================================================
message("52 PART 0: load")

raw <- data.table::fread(
  here::here("data", "orthotopic_series", "salmon.merged.gene_counts.tsv"),
  data.table = FALSE, check.names = FALSE)
stopifnot(identical(colnames(raw)[1:2], c("gene_id", "gene_name")))
CNT <- as.matrix(raw[, -(1:2)]); rownames(CNT) <- raw$gene_id
stopifnot(ncol(CNT) == 67L, all(grepl("^ENSMUSG", rownames(CNT))))

# GROUP PARSING, longest-token-first. `EV` and `Pgc` must be matched AFTER
# `NTEV`/`KOEV`/`NTPgc`/`KOPgc` or four groups silently collapse into two.
# THE RELABEL: the matrix says `Pgc1a`; the arm is PGC1a + Bcl-xL, so it is
# canonicalised to `Pgc1a_BclxL` here and never called `Pgc1a` again.
PAT <- c("NTEV", "KOEV", "NTPgc1a", "NTPgc", "KOPgc1a", "KOPgc", "BclxL",
         "BoMY1", "BoMY", "FaMY1", "FaMY", "MYAZ", "Pgc1a", "Pgc", "EV")
CANON <- c(NTEV = "NTEV", KOEV = "KOEV", NTPgc1a = "NTPgc1a", NTPgc = "NTPgc1a",
           KOPgc1a = "KOPgc1a", KOPgc = "KOPgc1a", BclxL = "BclxL",
           BoMY1 = "BoMY", BoMY = "BoMY", FaMY1 = "FaMY", FaMY = "FaMY",
           MYAZ = "MYAZ", Pgc1a = "Pgc1a_BclxL", Pgc = "Pgc1a_BclxL", EV = "EV")
core <- sub("^[A-H][0-9]-", "", colnames(CNT))
tok  <- vapply(core, function(s) {
  h <- PAT[vapply(PAT, function(p) grepl(paste0("^", p), s), TRUE)]
  if (length(h)) h[1] else NA_character_ }, "")
grp  <- factor(unname(CANON[tok]), levels = ALLGRP)

# THE COUNTS, ASSERTED RATHER THAN DESCRIBED. The spec carried "35 vector +
# CRISPR" through five documents before anyone multiplied 7 by 7; a number copied
# between prose files is never checked, and one in a stopifnot is checked every
# run. 7 groups x 7 = 49, and 49 + 18 = 67.
stopifnot(!anyNA(grp),
          identical(as.integer(table(grp)), REF$n),
          sum(grp %in% SET)  == 49L,
          sum(grp %in% EXCL) == 18L,
          sum(grp %in% SET) + sum(grp %in% EXCL) == ncol(CNT),
          length(SET) == 7L, all(table(grp)[SET] == 7L))
message(sprintf("52 PART 0: %d in the scoring set (7 x 7), %d excluded, %d total",
                sum(grp %in% SET), sum(grp %in% EXCL), ncol(CNT)))

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
message(sprintf("52 PART 0: median-CPM control, max relative deviation %.3f%% (tol %.1f%%)",
                100 * dev, 100 * CTRL_TOL_REL))
stopifnot(dev < CTRL_TOL_REL)
note("PART 0: the median-CPM table reproduces to ", sprintf("%.3f%%", 100 * dev),
     ". 50's and 51's VERDICTS are deliberately not re-derived -- they are void, ",
     "and re-deriving them would re-enter the wrong premise.")

# =============================================================================
# PART 0b: THE COHORT AND ITS TWO INPUT OBJECTS
# -----------------------------------------------------------------------------
# TWO OBJECTS THAT NEVER MIX. `NCv` is LINEAR normalised counts and is what
# mitoPPS takes; `L` is log2(NCv + 1) and feeds every z-composite. Counts are
# rounded for DESeq2 because salmon merged counts are non-integer.
# =============================================================================
vk  <- grp %in% SET
sm  <- data.frame(sample = colnames(CNT)[vk],
                  group  = factor(as.character(grp[vk]), levels = SET),
                  series = factor(ifelse(as.character(grp[vk]) %in% VEC, "vector",
                                  ifelse(as.character(grp[vk]) %in% c("NTEV", "NTPgc1a"),
                                         "CRISPR_NT", "CRISPR_KO")),
                                  levels = c("vector", "CRISPR_NT", "CRISPR_KO")))
stopifnot(nrow(sm) == 49L, all(table(sm$group) == 7L))

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(CNT[, vk, drop = FALSE]), colData = sm, design = ~ group)
dds  <- DESeq2::estimateSizeFactors(dds)
NCv  <- DESeq2::counts(dds, normalized = TRUE)          # LINEAR -- mitoPPS only
keep <- rowMeans(NCv) >= MINBM & matrixStats::rowSds(NCv) > 0
L    <- log2(NCv[keep, , drop = FALSE] + 1)             # LOG -- every composite
message(sprintf("52 PART 0b: %d of %d genes at mean normalised count >= %d",
                nrow(L), nrow(NCv), MINBM))

universe_all <- rownames(L)
zrow    <- function(m) t(scale(t(m)))
comp_e  <- function(e) { e <- e[!is.na(e) & e %in% rownames(L)]
                         stopifnot(length(e) >= 3); colMeans(zrow(L[e, , drop = FALSE])) }
ens_set <- function(syms) { e <- recon_to_ensembl(syms, universe_all); e[!is.na(e)] }
one_ens <- function(sym) { e <- ens_set(sym); if (!length(e)) NA_character_ else e[1] }
lg      <- function(sym) { e <- one_ens(sym); stopifnot(!is.na(e)); as.numeric(L[e, ]) }
nc_of   <- function(sym) { r <- which(raw$gene_name == sym); stopifnot(length(r) >= 1L)
                           m <- NCv[rownames(CNT)[r], , drop = FALSE]
                           if (nrow(m) == 1L) as.numeric(m) else colSums(m) }

# --- MitoCarta, script 08's split, exact filename, mt-* stripped -------------
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
stopifnot(length(intersect(ox_sub, mt_ens)) == 0L)

# THE RULER IS COHORT-DEPENDENT AND THAT IS ASSERTED, NOT ASSUMED. The expression
# filter runs over the 49, so a different gene set passes than over the 21: this
# cohort carries 86 nuclear OXPHOS subunits, not 50's and 51's 87. The one that
# drops is `Cox6b2`, at mean normalised count 10.16 over the 21 and 8.17 over the
# 49 against a floor of 10 -- a tissue-restricted paralog that CLAUDE.md already
# flags among the five lowest expressers, so nothing respiratory rests on it. The
# numbers are asserted because a silently different ruler changes every score.
set_sizes <- tibble::tibble(
  set = c("OXPHOS subunits (nuclear)", "rest of nuclear MitoCarta",
          "mtDNA-encoded", "genes carried"),
  n   = c(length(ox_sub), length(rest_mito), length(mt_ens), nrow(L)),
  note = c("86 here, 87 in the 21-sample run; Cox6b2 falls below the floor",
           "the ox_rel denominator", "in the denominator, never the numerator",
           "15788 here, 15721 in the 21-sample run"))
message(sprintf("52 PART 0b: ox_sub %d | rest %d | mtDNA %d | genes %d",
                length(ox_sub), length(rest_mito), length(mt_ens), nrow(L)))
stopifnot(length(ox_sub) == 86L, length(rest_mito) == 880L,
          length(mt_ens) == 13L, nrow(L) == 15788L)
note("PART 0b: the ruler is COHORT-DEPENDENT. 86 nuclear OXPHOS subunits here ",
     "against 87 over the 21 samples -- Cox6b2 sits on the expression floor ",
     "(10.16 -> 8.17). It is a tissue-restricted paralog and nothing rests on it, ",
     "but it is one more reason a score from this run is not comparable with one ",
     "from the 21-sample run.")

X <- list(ox_rel = comp_e(ox_sub) - comp_e(rest_mito),
          ox_lvl = comp_e(ox_sub),
          ox_mt  = comp_e(mt_ens),
          ox_den = comp_e(rest_mito))   # kept separately: PART B decomposes

# =============================================================================
# PART A: THE COMPOSITION GATE -- READ BEFORE ANY CLAIM
# -----------------------------------------------------------------------------
# The fat-pad rule: a loading belongs to the subset the test runs on, so it is
# computed inside the 49. The 09-08 note's lesson applies in BOTH directions -- a
# pooled |rho| can look material while the per-group loadings disagree, and can
# sit below threshold while every arm agrees strongly -- so the per-group
# breakdown is reported regardless of the pooled verdict.
#
# THE POSITIVE CONTROL IS `Ppargc1a`, read FIRST. It is a transgene, definitionally
# present in the arms that carry it. If the composition adjustment attenuates
# THAT, the adjustment is removing signal and nothing adjusted may be read.
# =============================================================================
message("52 PART A: composition gate")

MARK <- list(
  epithelial    = c("Krt8", "Krt18", "Epcam", "Cdh1", "Krt5", "Krt14"),
  stromal       = c("Col1a1", "Pdgfrb", "Acta2", "Thy1"),
  adipose       = c("Adipoq", "Plin1", "Fabp4", "Cidec"),
  immune        = c("Ptprc", "Cd52", "Lyz2", "Cd74"),
  endothelial   = c("Pecam1", "Cdh5", "Cldn5", "Emcn", "Tek"),
  proliferation = c("Mki67", "Top2a", "Ccnb1"))
MARK_e <- lapply(MARK, function(v) { e <- ens_set(v); e[e %in% rownames(L)] })
comp_present <- vapply(MARK_e, length, 0L)
dropped_comp <- names(comp_present)[comp_present < 3L]
if (length(dropped_comp))
  message(sprintf("52 PART A: dropped (fewer than 3 markers): %s",
                  paste(dropped_comp, collapse = ", ")))
K <- lapply(MARK_e[comp_present >= 3L], comp_e)

sp <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman"))
gate <- dplyr::bind_rows(lapply(names(K), function(cn)
  dplyr::bind_rows(lapply(c("ox_rel", "ox_lvl", "ox_mt"), function(xn) {
    per <- vapply(SET, function(g) sp(X[[xn]][sm$group == g], K[[cn]][sm$group == g]), 0)
    tibble::tibble(compartment = cn, n_genes = comp_present[[cn]], ruler = xn,
                   rho_pooled = sp(X[[xn]], K[[cn]]),
                   !!!stats::setNames(as.list(per), paste0("rho_", SET)),
                   n_arms_same_sign = max(sum(per > 0), sum(per < 0))) }))))

material <- unique(gate$compartment[gate$ruler == "ox_rel" &
                                    abs(gate$rho_pooled) >= CUT_MAT])
message(sprintf("52 PART A: material composition axes (|rho| >= %.1f): %s",
                CUT_MAT, if (length(material)) paste(material, collapse = ", ") else "NONE"))
AD <- if (length(material)) as.data.frame(K[material]) else NULL

# --- the positive control, read FIRST ----------------------------------------
adj_group <- function(v, hi, lo, lab) {
  k <- sm$group %in% c(hi, lo)
  d <- data.frame(y = v[k], g = factor(as.character(sm$group[k]), levels = c(lo, hi)))
  if (!is.null(AD)) d <- cbind(d, AD[k, , drop = FALSE])
  s <- summary(stats::lm(y ~ ., d))$coefficients
  rn <- grep(paste0("^g", hi, "$"), rownames(s), value = TRUE)
  tibble::tibble(measure = lab, hi = hi, lo = lo,
                 adjusted_for = if (length(material)) paste(material, collapse = "+") else "nothing",
                 beta = if (length(rn)) s[rn, 1] else NA_real_,
                 p    = if (length(rn)) s[rn, 4] else NA_real_) }

ppargc1a_cpm <- cpm_of("Ppargc1a")[vk]
gate_control <- dplyr::bind_rows(
  adj_group(ppargc1a_cpm, "Pgc1a_BclxL", "BclxL", "Ppargc1a CPM -- POSITIVE CONTROL (vector)"),
  adj_group(ppargc1a_cpm, "NTPgc1a",     "NTEV",  "Ppargc1a CPM -- POSITIVE CONTROL (NT)"))
gate_ok <- all(gate_control$beta > 0 & gate_control$p < 0.05, na.rm = TRUE)
note("PART A: positive control (Ppargc1a, a transgene) ",
     if (gate_ok) "SURVIVES the adjustment -- adjusted results below are readable."
     else "DOES NOT survive -- the adjustment is removing signal and NOTHING adjusted may be read.")

# =============================================================================
# PART B: THE RESPIRATORY RULERS ACROSS ALL SEVEN ARMS
# -----------------------------------------------------------------------------
# `ox_rel` is decomposed into its numerator and denominator per arm, so a rise by
# MASS EXPANSION (both halves up) is separable from a rise by TILT (numerator up,
# denominator flat). mitoPPS is computed here on its OWN LINEAR object.
# `NTPgc1a` vs `NTEV` is the contrast that decides R1.
# =============================================================================
message("52 PART B: rulers")

boot_diff <- function(a, b, f = stats::median) {
  d <- vapply(seq_len(NBOOT), function(i)
    f(sample(a, replace = TRUE)) - f(sample(b, replace = TRUE)), 0)
  stats::quantile(d, c(0.025, 0.975), na.rm = TRUE) }
two_group <- function(v, hi, lo, one_sided = FALSE) {
  a <- v[sm$group == hi]; b <- v[sm$group == lo]
  ci <- boot_diff(a, b)
  tibble::tibble(hi = hi, lo = lo, n_hi = length(a), n_lo = length(b),
                 median_hi = stats::median(a), median_lo = stats::median(b),
                 diff_median = stats::median(a) - stats::median(b),
                 ci_lo = ci[[1]], ci_hi = ci[[2]],
                 wilcox_p = suppressWarnings(stats::wilcox.test(
                   a, b, alternative = if (one_sided) "greater" else "two.sided")$p.value)) }

# the contrast table, with its cleanliness stated in the object rather than a comment
CONTRASTS <- tibble::tribble(
  ~hi,            ~lo,      ~series,       ~status,
  "Pgc1a_BclxL",  "BclxL",  "vector",      "CLEAN -- Bcl-xL held constant, same split pool",
  "NTPgc1a",      "NTEV",   "CRISPR_NT",   "CLEAN -- construct-free PGC1a contrast",
  "KOPgc1a",      "KOEV",   "CRISPR_KO",   "CLEAN -- but the transgene is lost in KOPgc1a",
  "Pgc1a_BclxL",  "EV",     "vector",      "TWO-FACTOR -- not a PGC1a contrast; this is the void C1",
  "BclxL",        "EV",     "vector",      "Bcl-xL alone",
  "NTEV",         "EV",     "cross",       "DESCRIPTIVE ONLY -- clone against polyclonal pool")
WITHIN <- CONTRASTS[startsWith(CONTRASTS$status, "CLEAN"), ]

RULERS <- list(ox_lvl = X$ox_lvl, ox_rel = X$ox_rel, ox_mt = X$ox_mt,
               ox_den = X$ox_den)

# --- mitoPPS, on the LINEAR object, its own input, never shared ---------------
# Monzel's algorithm as script 08 implements it: per-sample pairwise ratios
# between pathway scores, each ratio divided by its cross-sample average, then
# averaged per pathway. Written as matrix algebra here; the loop-based
# equivalence is asserted below on a subset rather than assumed.
pw_syms  <- split(g2p_nuc$Gene, g2p_nuc$Pathway)
pw_syms[["mtDNA-encoded OXPHOS subunits"]] <- mt_syms
sym_row  <- function(syms) { e <- ens_set(syms); e <- e[e %in% rownames(NCv)]; e }
pw_ens   <- lapply(pw_syms, sym_row)
pw_ens   <- pw_ens[vapply(pw_ens, length, 0L) >= 3L]
S <- t(vapply(pw_ens, function(e) colMeans(NCv[e, , drop = FALSE]), numeric(ncol(NCv))))
S <- S[rowSums(!is.finite(S)) == 0 & matrixStats::rowMins(S) > 0, , drop = FALSE]

mitopps_of <- function(S) {
  P <- nrow(S)
  R <- lapply(seq_len(ncol(S)), function(s) outer(S[, s], 1 / S[, s]))
  A <- Reduce(`+`, R) / length(R)
  out <- vapply(seq_along(R), function(s) {
    C <- R[[s]] / A; diag(C) <- NA_real_; rowMeans(C, na.rm = TRUE) }, numeric(P))
  dimnames(out) <- list(rownames(S), colnames(S)); out }
MP <- mitopps_of(S)
# equivalence check against the definition, on one sample and three pathways
{ s <- 1L; ii <- head(seq_len(nrow(S)), 3L)
  A_chk <- sapply(seq_len(nrow(S)), function(j)
            sapply(seq_len(nrow(S)), function(i) mean(S[i, ] / S[j, ])))
  brute <- vapply(ii, function(i) mean(vapply(setdiff(seq_len(nrow(S)), i),
              function(j) (S[i, s] / S[j, s]) / A_chk[i, j], 0)), 0)
  stopifnot(max(abs(brute - MP[ii, s])) < 1e-9) }
stopifnot("OXPHOS subunits" %in% rownames(MP))
RULERS$mitopps_oxphos <- as.numeric(MP["OXPHOS subunits", ])
note("PART B: mitoPPS computed on the LINEAR normalised counts in its own object, ",
     sprintf("%d pathways after the >=3-gene filter. ", nrow(S)),
     "The vectorised implementation is asserted equal to the definition, not assumed.")

rulers <- dplyr::bind_rows(lapply(names(RULERS), function(rn)
  dplyr::bind_rows(lapply(seq_len(nrow(CONTRASTS)), function(i)
    dplyr::bind_cols(tibble::tibble(ruler = rn, series = CONTRASTS$series[i],
                                    status = CONTRASTS$status[i]),
                     two_group(RULERS[[rn]], CONTRASTS$hi[i], CONTRASTS$lo[i]))))))

# R1 is decided here and nowhere else
r1_row <- rulers[rulers$ruler == "ox_lvl" & rulers$hi == "NTPgc1a" & rulers$lo == "NTEV", ]
r1_small <- abs(r1_row$diff_median[1]) < 0.5 | (r1_row$ci_lo[1] < 0 & r1_row$ci_hi[1] > 0)

# =============================================================================
# PART C: THE ESCAPE TABLE -- the display item
# =============================================================================
message("52 PART C: escape table")

p80 <- function(z) unname(stats::quantile(z, 0.8, type = 7))
myc_cpm <- cpm_of("Myc")[vk]
escape_table <- dplyr::bind_rows(lapply(SET, function(g) {
  k <- sm$group == g
  tibble::tibble(
    arm = g,
    series = as.character(sm$series[k][1]),
    n = sum(k),
    Ppargc1a_CPM = stats::median(ppargc1a_cpm[k]),
    ox_lvl = stats::median(X$ox_lvl[k]),
    ox_rel = stats::median(X$ox_rel[k]),
    ox_mt  = stats::median(X$ox_mt[k]),
    mitopps_oxphos = stats::median(RULERS$mitopps_oxphos[k]),
    Myc_median = stats::median(myc_cpm[k]),
    Myc_p80    = p80(myc_cpm[k]),
    Myc_max    = max(myc_cpm[k])) }))
escape_table$arm <- factor(escape_table$arm, levels = SET)
escape_table <- escape_table[order(escape_table$arm), ]

# =============================================================================
# PART D: THE p53 TARGET PANEL -- a PRE-DECLARED test of the correction's §2
# -----------------------------------------------------------------------------
# Section 2 explains `Bbc3`'s fall in Pgc1a_BclxL as a readout of WHICH ESCAPE
# ROUTE that arm took: MYAZ is p53 proficient, the PGC1a survivors have p19^ARF
# and p21 down by western, and Bbc3 is a canonical p53 target coming down with the
# rest of the programme. R3 tests that: the panel should move down TOGETHER with
# Bbc3, not Bbc3 alone.
#
# MOUSE `Cdkn2a` POOLS p16^INK4a AND p19^ARF from alternative first exons, and
# gene-level counts cannot separate them. The western separates what these counts
# cannot; a flat `Cdkn2a` here does not exclude p19^ARF loss.
# =============================================================================
message("52 PART D: p53 panel")

P53 <- c("Cdkn2a", "Cdkn1a", "Mdm2", "Ccng1", "Zmat3", "Bax", "Trp53", "Bbc3")
p53_panel <- dplyr::bind_rows(lapply(P53, function(s) {
  v <- lg(s)
  dplyr::bind_rows(lapply(seq_len(nrow(WITHIN)), function(i)
    dplyr::bind_cols(tibble::tibble(gene = s, series = WITHIN$series[i]),
                     two_group(v, WITHIN$hi[i], WITHIN$lo[i])))) }))
p53_panel$direction <- ifelse(p53_panel$ci_hi < 0, "down",
                       ifelse(p53_panel$ci_lo > 0, "up", "spans zero"))
note("PART D: mouse Cdkn2a pools p16INK4a and p19ARF from alternative first ",
     "exons and gene-level counts cannot separate them -- the western separates ",
     "what these counts cannot, and a flat Cdkn2a here does not exclude ARF loss.")

# R3 is decided here
vec_p53 <- p53_panel[p53_panel$series == "vector" &
                     p53_panel$gene %in% c("Cdkn2a", "Cdkn1a", "Mdm2", "Ccng1", "Zmat3"), ]
bbc3_vec <- p53_panel[p53_panel$series == "vector" & p53_panel$gene == "Bbc3", ]
n_down   <- sum(vec_p53$direction == "down")

# =============================================================================
# PART E: THE CONSTRUCT-FREE TWELVE
# -----------------------------------------------------------------------------
# `NTPgc1a` vs `NTEV` is the ONLY place in this cohort where `Bcl2l1` means
# endogenous BCL-XL. Reported for the construct-carrying arms only behind an
# explicit uninterpretable flag. No per-gene FDR: estimates with intervals.
# `Bik` and `Hrk` sit below 1 CPM and are reported, never fitted -- a bootstrap
# interval on 0.03 CPM is noise wearing an interval.
# =============================================================================
message("52 PART E: the construct-free twelve")

TWELVE <- c("Bcl2l1", "Mcl1", "Bcl2", "Bcl2l2", "Bcl2a1a", "Bcl2l11", "Bbc3",
            "Pmaip1", "Bid", "Bad", "Bik", "Bmf", "Hrk", "Bax", "Bak1", "Foxo3")
med_cpm <- vapply(TWELVE, function(s) stats::median(cpm_of(s)[vk]), 0)
fitable <- TWELVE[med_cpm >= LOWCPM]
lowflag <- TWELVE[med_cpm <  LOWCPM]
if (length(lowflag))
  message(sprintf("52 PART E: below %g CPM, reported not fitted: %s",
                  LOWCPM, paste(lowflag, collapse = ", ")))

twelve_constructfree <- dplyr::bind_rows(lapply(fitable, function(s)
  dplyr::bind_cols(tibble::tibble(gene = s, median_CPM = med_cpm[[s]]),
                   two_group(lg(s), "NTPgc1a", "NTEV"))))
twelve_constructfree$direction <- ifelse(twelve_constructfree$ci_hi < 0, "down",
                                  ifelse(twelve_constructfree$ci_lo > 0, "up", "spans zero"))
twelve_lowexpr <- tibble::tibble(gene = lowflag, median_CPM = med_cpm[lowflag],
                                 note = "below the 1 CPM floor -- reported, not fitted")

# `Bcl2l1` in the construct-carrying arms, behind its flag
bcl2l1_constructarms <- dplyr::bind_cols(
  tibble::tibble(flag = "UNINTERPRETABLE -- construct and endogenous transcript are not separable at gene level; BCL-XL protein is equal between these arms by western"),
  two_group(lg("Bcl2l1"), "Pgc1a_BclxL", "BclxL"))

# R2 and R4 are decided here
r2_row <- twelve_constructfree[twelve_constructfree$gene == "Bbc3", ]
r4_row <- twelve_constructfree[twelve_constructfree$gene == "Bcl2l1", ]

# =============================================================================
# PART F: C2 ON THE ESCAPE SERIES -- AN OBSERVATION, NOT A TEST
# -----------------------------------------------------------------------------
# `Myc` upper tail on every WITHIN-SERIES contrast, both statistics, both
# normalisations, exact null over all choose(14,7) = 3432 relabellings. The by-arm
# size factor and MitoCarta share sit beside them because script 51 showed CPM
# deflation moves this quantity: an arm that nearly doubles its mitochondrial
# share deflates every other gene's CPM without any gene changing.
# =============================================================================
message("52 PART F: C2 on the escape series")

exact_null <- function(v, hi, lo, stat) {
  a <- v[sm$group == hi]; b <- v[sm$group == lo]
  pool <- c(a, b); n <- length(a)
  idx  <- utils::combn(length(pool), n)
  null <- apply(idx, 2, function(k) stat(pool[k]) - stat(pool[-k]))
  obs  <- stat(a) - stat(b)
  list(obs = obs, p_lower = mean(null <= obs), n_perm = ncol(idx),
       null_median = stats::median(null)) }

myc_mor <- nc_of("Myc")
c2 <- dplyr::bind_rows(lapply(seq_len(nrow(WITHIN)), function(i)
  dplyr::bind_rows(lapply(list(`CPM` = myc_cpm, `median-of-ratios` = myc_mor), function(v)
    dplyr::bind_rows(lapply(list(max = max, pct80 = p80), function(st) {
      r <- exact_null(v, WITHIN$hi[i], WITHIN$lo[i], st)
      tibble::tibble(series = WITHIN$series[i], hi = WITHIN$hi[i], lo = WITHIN$lo[i],
                     stat_obs = r$obs, p_one_sided_lower = r$p_lower,
                     n_perm = r$n_perm, null_median = r$null_median) }),
      .id = "statistic")), .id = "scale")))
stopifnot(all(c2$n_perm == choose(14, 7)))

mito_all_ens <- unique(c(ox_sub, rest_mito, mt_ens))
c2_diag <- dplyr::bind_rows(lapply(SET, function(g) {
  k <- sm$group == g
  tibble::tibble(arm = g,
                 size_factor_median = stats::median(DESeq2::sizeFactors(dds)[k]),
                 mitocarta_share = stats::median(
                   colSums(NCv[rownames(NCv) %in% mito_all_ens, k, drop = FALSE]) /
                   colSums(NCv[, k, drop = FALSE]))) }))

# =============================================================================
# PART G: VERDICT ON THE HEADER'S READING RULES
# =============================================================================
message("52 PART G: verdict")

verdict <- tibble::tibble(
  rule = c("GATE  composition separable / control holds",
           "R1    NTPgc1a respiratory response",
           "R2    Bbc3 in the construct-free arm",
           "R3    the p53 panel moves with Bbc3",
           "R4    endogenous Bcl2l1"),
  result = c(
    sprintf("%s | material: %s",
            if (gate_ok) "control HOLDS -- adjusted results readable"
            else "CONTROL FAILS -- nothing adjusted may be read",
            if (length(material)) paste(material, collapse = ", ") else "NONE"),
    sprintf("ox_lvl %+.3f [%.3f, %.3f] -- %s", r1_row$diff_median[1], r1_row$ci_lo[1],
            r1_row$ci_hi[1],
            if (r1_small) "SMALL: section 3's dose-response reading HOLDS"
            else "SUBSTANTIAL: check Myc; if flat, section 3's table must be REWRITTEN"),
    sprintf("%+.3f [%.3f, %.3f] (%s) -- %s", r2_row$diff_median[1], r2_row$ci_lo[1],
            r2_row$ci_hi[1], r2_row$direction[1],
            if (identical(r2_row$direction[1], "down"))
              "FALLS here too: the section 2 p53 account is INSUFFICIENT"
            else "flat or rising: the Pgc1a_BclxL fall is an ESCAPE-ROUTE SIGNATURE"),
    sprintf("%d of 5 panel genes down with Bbc3 %s -- %s", n_down, bbc3_vec$direction[1],
            if (n_down >= 3 && identical(bbc3_vec$direction[1], "down"))
              "moving TOGETHER: escape-route explanation SUPPORTED"
            else "Bbc3 is not moving with the panel: the fall needs a different account"),
    sprintf("%+.3f [%.3f, %.3f] (%s) -- %s", r4_row$diff_median[1], r4_row$ci_lo[1],
            r4_row$ci_hi[1], r4_row$direction[1],
            if (identical(r4_row$direction[1], "spans zero"))
              "NULL: this CLOSES C1 rather than leaving it merely void"
            else "not null: report it, and say which direction")))

note("SELECTION, NOT MANIPULATION: every arm is a survivor population. Nothing ",
     "here may be worded as 'PGC1a regulates X'; the form is 'tumours that ",
     "survived PGC1a expression show X'.")
note("SCORING SET: all 49 in one run. NEVER quote a value from this run beside a ",
     "value from the 21-sample run behind scripts 50 and 51.")
note("C3, the guardian ratio and the guardian-sensitiser gap are CLOSED and are ",
     "not computed here in any form. No MYC x OXPHOS interaction is estimable.")

stopifnot(nrow(sm) == 49L, nrow(escape_table) == 7L, nrow(WITHIN) == 3L,
          all(c("ox_lvl", "ox_rel", "ox_mt", "mitopps_oxphos") %in% rulers$ruler))

res <- list(
  load_control = load_control, reference_table = REF, set_sizes = set_sizes,
  scoring_set = list(
    groups = SET, n = 49L, excluded = EXCL, n_excluded = 18L,
    rule = paste("All 49 vector + CRISPR samples scored in ONE run. ox_rel, GSVA",
                 "and mitoPPS are cohort-relative: a value from this run is NEVER",
                 "quotable beside a value from the 21-sample run behind scripts 50",
                 "and 51, which stays untouched as the record.")),
  sample_table = sm, contrasts = CONTRASTS, within_series = WITHIN,
  gate = gate, gate_control = gate_control, material_axes = material,
  gate_ok = gate_ok,
  rulers = rulers, mitopps_matrix = MP,
  escape_table = escape_table,
  p53_panel = p53_panel,
  twelve_constructfree = twelve_constructfree, twelve_lowexpr = twelve_lowexpr,
  bcl2l1_constructarms = bcl2l1_constructarms,
  c2 = c2, c2_diagnostic = c2_diag,
  verdict = verdict,
  params = list(NBOOT = NBOOT, MINBM = MINBM, CUT_MAT = CUT_MAT, LOWCPM = LOWCPM, seed = 1),
  analysis_date = Sys.Date(), notes = NOTES)

saveRDS(res, here::here("results", "orthotopic_escape_series.rds"))
message("52: wrote results/orthotopic_escape_series.rds")

cat("\n================ VERDICT ================\n")
print(as.data.frame(verdict), right = FALSE)

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  res <- readRDS(here::here("results", "orthotopic_escape_series.rds"))

  ## 1. THE GATE, BEFORE ANY CLAIM. Read the POSITIVE CONTROL row first: if
  ## Ppargc1a -- a transgene, definitionally present -- does not survive the
  ## adjustment, the adjustment is eating signal and nothing adjusted counts.
  res$gate_control |> print()
  res$material_axes |> print()
  ## per-group as well as pooled, because the pooled rho misleads BOTH ways
  subset(res$gate, ruler == "ox_rel") |> print(n = 30)

  ## 2. R1 -- does NTPgc1a's respiration move? This decides whether the
  ## correction document's section 3 dose-response table stands as written.
  subset(res$rulers, hi == "NTPgc1a" & lo == "NTEV") |> print()

  ## 3. THE ESCAPE TABLE -- the display item. Read DOWN the Ppargc1a column and
  ## across to Myc: the claim is that the dose a tumour kept is set by how much of
  ## the death pathway it broke.
  res$escape_table |> print(n = 10)

  ## 4. R2 and R4 -- the construct-free arm, the only place Bcl2l1 means
  ## endogenous BCL-XL. Bik and Hrk are below the floor and are not fitted.
  res$twelve_constructfree |> print(n = 20)
  res$twelve_lowexpr |> print()
  ## and the construct-carrying arms, behind their flag
  res$bcl2l1_constructarms |> print()

  ## 5. R3 -- does the p53 panel move WITH Bbc3, or is Bbc3 alone?
  subset(res$p53_panel, series == "vector") |> print(n = 10)
  res$p53_panel |> print(n = 30)

  ## 6. C2 as an OBSERVATION. Compare max against the 80th percentile, and CPM
  ## against median-of-ratios; the diagnostic says whether an arm's mitochondrial
  ## share is deflating everything else's CPM.
  res$c2 |> print(n = 20)
  res$c2_diagnostic |> print()

  res$verdict |> as.data.frame() |> print(right = FALSE)
  res$scoring_set$rule |> print()
  cat(res$notes, sep = "\n")
}
