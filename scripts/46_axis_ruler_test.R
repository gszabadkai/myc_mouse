# =============================================================================
# 46_axis_ruler_test.R -- does Fig. 2I's interaction survive a change of RULER?
# =============================================================================
#
# THE QUESTION. Fig. 2I draws the animal's OXPHOS-subunit **mitoPPS** score against
# its PUMA:Bcl-xL ratio, and the statistic is the difference between the two
# genotypes' slopes (script 43 PART B: +6.09, p = 0.0052). mitoPPS is a PAIRWISE
# RATIO inside the mitochondrial compartment, so a high value means the compartment
# spends more of its budget on the respiratory chain -- **not that it respires
# more**. The panel's own legend already has to say so.
#
# The closing hypothesis the panel supports is that a HIGH RESPIRATORY STATE is
# required for MYC-PUMA mediated death. That is a statement about how much
# respiratory chain there is, and the ruler that speaks to it is the ABSOLUTE
# LEVEL -- summed DESeq2-normalised counts over the same 87 nuclear OXPHOS
# subunits, which script 45 PART G already computes per animal.
#
# So: refit the identical model with the identical outcome and covariates, change
# ONLY the axis, and report both.
#
# THIS IS A TEST, NOT A FIX, AND IT CAN FAIL. A ratio cancels the global
# common-mode axis; an absolute level does not, and script 36 measured
# cor(mito_oxphos, global mean) = 0.968 in score space. The level ruler's axis may
# therefore BE the global expression axis rather than the respiratory chain
# specifically. PART E measures that directly so a null result can be read
# correctly instead of being over-interpreted in either direction.
#
# THE DECISION RULE IS FIXED BEFORE THE RESULT IS SEEN (see VERDICT_RULE below,
# and the notes block). Whichever way it comes out, BOTH rulers are reported:
# CLAUDE.md's standing rule is that a sentence quoting an OXPHOS abundance number
# must name its ruler, and a ruler-dependence is a finding, not a failure.
#
# WHAT THIS SCRIPT IS NOT. It does not refit DESeq2, does not rebuild any gene set,
# and does not touch the contrasts. It is script 43 PART B's model with a second
# axis pair, plus the diagnostics needed to read the comparison.
#
# Reads:
#   results/dds_int_run.rds                    (script 03) -- script 43's L, sm, covariates
#   results/dds_group_run.rds                  (script 03) -- size factors for the levels
#   results/count_matrix.rds                   (script 01) -- UNfiltered counts
#   results/combined_df_annotated.rds          (script 03) -- symbol -> Ensembl
#   results/mitopps_scores.rds                 (script 08) -- mitoPPS axes + MitoCarta partition
#   results/substrate_specificity_tradeoff.rds (script 43) -- the numbers to reproduce
#   results/state_readings.rds                 (script 45) -- the saved OXPHOS level, asserted
#   data/genesets_from_library/...gmt          -- PROLIF_E2F_HALLMARK only
# Writes: results/axis_ruler_test.rds
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

set.seed(1)
NPERM <- 5000L   # within-timepoint axis shuffles, script 43 PART B's scheme

# THE DECISION RULE, written down before the result is looked at.
VERDICT_RULE <- c(
  reproduces = "same sign as the mitoPPS ruler AND nominal p < 0.05 on model m1",
  partial    = "same sign AND nominal p < 0.20",
  fails      = "opposite sign OR p >= 0.20")
P_REPRO <- 0.05
P_PART  <- 0.20
RESID_FRAC_GATE <- 0.10   # script 36's readability gate: below this the axis IS
                          # the global factor and is not respiratory-chain-specific

# =============================================================================
# PART 0: LOAD, AND ALIGN EVERY OBJECT TO ONE SAMPLE ORDER
# =============================================================================
# Three different DESeq2 objects and a raw count matrix are involved and their
# column orders are not guaranteed to agree. Everything is aligned to colnames(L),
# which is script 43's order, and the alignment is asserted rather than assumed.
dds <- readRDS(here::here("results", "dds_int_run.rds"))
dg  <- readRDS(here::here("results", "dds_group_run.rds"))
cts <- readRDS(here::here("results", "count_matrix.rds"))
cdf <- readRDS(here::here("results", "combined_df_annotated.rds"))
mp  <- readRDS(here::here("results", "mitopps_scores.rds"))
ss  <- readRDS(here::here("results", "substrate_specificity_tradeoff.rds"))
sr  <- readRDS(here::here("results", "state_readings.rds"))
gmt <- fgsea::gmtPathways(here::here("data", "genesets_from_library",
                                     "mammary_mito_myc_metab_v1_mouse.gmt"))

nc <- DESeq2::counts(dds, normalized = TRUE)
L  <- log2(nc + 1)
sm <- as.data.frame(SummarizedExperiment::colData(dds))
sm$tp  <- stats::relevel(as.factor(sm$timepoint),  "6W")
sm$myc <- stats::relevel(as.factor(sm$myc_status), "neg")
samples <- colnames(L)
stopifnot(length(samples) == 24L, all(table(sm$tp, sm$myc) == 6L))

universe_all <- rownames(L)
ens_set <- function(syms) { e <- recon_to_ensembl(syms, universe_all); e[!is.na(e)] }
sym2ens <- cdf[!is.na(cdf$mgi_symbol) & !duplicated(cdf$mgi_symbol), c("mgi_symbol", "gene")]
ens_of  <- function(s) sym2ens$gene[match(s, sym2ens$mgi_symbol)]
zrow    <- function(m) t(scale(t(m)))
comp_of <- function(syms) {
  e <- ens_of(syms); e <- e[!is.na(e) & e %in% rownames(L)]
  stopifnot(length(e) >= 3)
  colMeans(zrow(L[e, , drop = FALSE]))
}
set_score <- function(set, min_n = 5L) {
  if (is.null(gmt[[set]])) return(NULL)
  e <- ens_set(gmt[[set]]); e <- e[e %in% rownames(L)]
  if (length(e) < min_n) return(NULL)
  colMeans(zrow(L[e, , drop = FALSE]))
}

# script 43:100-101 verbatim -- the composition covariates the model adjusts for
epi_comp <- comp_of(c("Epcam", "Krt8", "Krt18", "Krt5", "Krt14", "Cdh1", "Krt17"))
imm_comp <- comp_of(c("Ptprc", "Cd52", "Cd3e", "Lyz2", "Cd74", "H2-Aa", "Itgam", "Ms4a1"))

# --- the mitoPPS side, aligned ------------------------------------------------
mps <- mp$mitopps_scores
mps <- mps[match(samples, mps$sample), , drop = FALSE]
stopifnot(identical(as.character(mps$sample), samples))
OX_P <- "OXPHOS subunits"
RX_P <- "ROS and glutathione metabolism"
stopifnot(all(c(OX_P, RX_P, "Apoptosis-PRO", "Apoptosis-ANTI") %in% names(mps)))

message(sprintf("PART 0: %d samples aligned across dds_int_run, dds_group_run, counts and mitoPPS",
                length(samples)))

# =============================================================================
# PART A: POSITIVE CONTROL -- reproduce script 43 PART B before changing anything
# =============================================================================
# If the rebuilt model does not return script 43's recorded statistics, the
# comparison in PART C is between two things that differ in more than the ruler,
# and nothing downstream is readable. Same discipline as script 45 PART C/D.
ratio_of <- function(a, b) {
  ea <- ens_of(a); eb <- ens_of(b)
  stopifnot(!is.na(ea), !is.na(eb), ea %in% rownames(L), eb %in% rownames(L))
  as.numeric(L[ea, ] - L[eb, ])
}
pr_lvl <- comp_of(c("Mki67", "Top2a", "Ccnb1", "Ccna2", "Bub1", "Plk1", "Aurka",
                    "Cdk1", "Pcna", "Mcm2", "Rrm2", "Tk1"))

outcomes <- list()
outcomes[["priming_ppd (PRO-ANTI)"]] <-
  as.numeric(mps[["Apoptosis-PRO"]]) - as.numeric(mps[["Apoptosis-ANTI"]])
outcomes[["Bbc3:Bcl2l1 (PUMA priming)"]] <- ratio_of("Bbc3", "Bcl2l1")
outcomes[["Bax:Bcl2l1"]]                 <- ratio_of("Bax",  "Bcl2l1")
outcomes[["proliferation (markers)"]]    <- pr_lvl
pp <- set_score("PROLIF_E2F_HALLMARK")
if (!is.null(pp)) outcomes[["proliferation (E2F hallmark)"]] <- pp

HEADLINE <- "Bbc3:Bcl2l1 (PUMA priming)"
stopifnot(HEADLINE %in% names(outcomes))

Z <- function(x) as.numeric(scale(x))
resid_design <- function(y)
  stats::residuals(stats::lm(y ~ tp * myc + epi + imm,
                             data.frame(y = y, tp = sm$tp, myc = sm$myc,
                                        epi = epi_comp, imm = imm_comp)))

# script 43:322-337's fit, extracted so both rulers go through one code path
fit_one <- function(y_raw, a, ry = NULL) {
  y <- Z(y_raw)
  d  <- data.frame(y = y, tp = sm$tp, myc = sm$myc, epi = epi_comp,
                   imm = imm_comp, a = a)
  m1 <- summary(stats::lm(y ~ myc * a + epi + imm, d))$coefficients
  m2 <- summary(stats::lm(y ~ tp * myc + myc:a + a + epi + imm, d))$coefficients
  if (is.null(ry)) ry <- resid_design(y_raw)
  tibble::tibble(
    rho_raw            = stats::cor(a, y_raw, method = "spearman"),
    rho_adj            = stats::cor(resid_design(a), ry, method = "spearman"),
    myc_x_axis         = m1["mycpos:a", 1], se = m1["mycpos:a", 2],
    p                  = m1["mycpos:a", 4],
    myc_x_axis_with_tp = m2["mycpos:a", 1], p_with_tp = m2["mycpos:a", 4])
}

axes_ppd <- list(oxphos_ppd = as.numeric(mps[[OX_P]]),
                 redox_ppd  = as.numeric(mps[[RX_P]]))

# dplyr::bind_cols, NOT cbind: cbind() on two tibbles dispatches to
# cbind.data.frame and silently returns a data.frame, and a data.frame reaching
# print(n = ) reads `n` as `na.print` (R_CODING_INSTRUCTIONS rule 1 -- the trap
# that killed script 45's first run). Everything saved here stays a tibble.
control <- dplyr::bind_rows(lapply(names(outcomes), function(on) {
  ry <- resid_design(outcomes[[on]])
  dplyr::bind_rows(lapply(names(axes_ppd), function(an) {
    dplyr::bind_cols(tibble::tibble(outcome = on, axis = an),
                     fit_one(outcomes[[on]], axes_ppd[[an]], ry))
  }))
}))
stopifnot(tibble::is_tibble(control))

rec <- as.data.frame(ss$tradeoff)
key <- paste(control$outcome, control$axis)
mm  <- match(key, paste(rec$outcome, rec$axis))
stopifnot(!anyNA(mm))
dev_beta <- max(abs(control$myc_x_axis - rec$myc_x_axis[mm]))
dev_p    <- max(abs(control$p          - rec$p[mm]))
message(sprintf("PART A: reproduces script 43 -- max |dbeta| = %.2e, max |dp| = %.2e over %d fits",
                dev_beta, dev_p, nrow(control)))
stopifnot(dev_beta < 1e-8, dev_p < 1e-8)

# =============================================================================
# PART B: THE SECOND RULER -- absolute levels, script 45 PART G's construction
# =============================================================================
# Size factors from the fitted GROUP model applied to the FULL count matrix, NOT
# to counts(dg): the fitted object carries script 03's expression filter, which
# would silently under-count the compartment (45:405-411).
sf <- DESeq2::sizeFactors(dg)
stopifnot(!is.null(sf), all(samples %in% names(sf)), all(samples %in% colnames(cts)))
nrm <- sweep(cts[, samples, drop = FALSE], 2, sf[samples], "/")

gp   <- mp$gene_to_pathway
pcol <- names(gp)[1]; gcol <- names(gp)[2]
paths_all <- setdiff(unique(gp[[pcol]]), mp$mtdna_pathway_name)

path_ens_lvl <- lapply(stats::setNames(paths_all, paths_all), function(p) {
  e <- ens_set(unique(gp[[gcol]][gp[[pcol]] == p]))
  e[e %in% rownames(nrm)]
})
path_ens_lvl <- path_ens_lvl[vapply(path_ens_lvl, length, integer(1)) >= 3L]
# A pathway summing to zero in any animal would give -Inf and poison the ambient
# and the global factor. None is expected; drop and report rather than propagate.
pos_sum <- vapply(path_ens_lvl, function(e)
  min(colSums(nrm[e, , drop = FALSE])) > 0, logical(1))
if (any(!pos_sum))
  message(sprintf("PART B: dropped %d pathway(s) with a zero sum in some animal", sum(!pos_sum)))
path_ens_lvl <- path_ens_lvl[pos_sum]
lvl_of <- function(p) log2(colSums(nrm[path_ens_lvl[[p]], , drop = FALSE]))

stopifnot(all(c(OX_P, RX_P) %in% names(path_ens_lvl)))
axes_lvl <- list(oxphos_lvl = unname(lvl_of(OX_P)),
                 redox_lvl  = unname(lvl_of(RX_P)))

# Cross-check against script 45's saved levels: same genes, same construction, so
# the OXPHOS row must agree exactly. A mismatch means the membership resolution
# has drifted and the two scripts are no longer talking about the same 87 genes.
lv45 <- as.data.frame(sr$levels)
o45  <- lv45[lv45$group_set == OX_P, ]
o45  <- o45[match(samples, as.character(o45$sample)), ]
stopifnot(nrow(o45) == 24L, !anyNA(o45$norm_sum))
dev_lvl <- max(abs(axes_lvl$oxphos_lvl - log2(o45$norm_sum)))
n_ox <- length(path_ens_lvl[[OX_P]]); n_rx <- length(path_ens_lvl[[RX_P]])
message(sprintf("PART B: %s = %d genes, %s = %d genes; max |dlevel| vs script 45 = %.2e",
                OX_P, n_ox, RX_P, n_rx, dev_lvl))
stopifnot(dev_lvl < 1e-9, n_ox == unique(o45$n_genes))

# =============================================================================
# PART C: THE HEAD-TO-HEAD -- one model, four axes, five outcomes
# =============================================================================
# The axis is z-scored so the interaction COEFFICIENTS are comparable between two
# rulers on completely different scales. A linear rescale of `a` leaves t and p
# exactly unchanged, so the p-values compare without any standardisation at all --
# asserted below, because the whole test rests on it.
axes <- c(axes_ppd, axes_lvl)
RULER <- c(oxphos_ppd = "mitoPPS", redox_ppd = "mitoPPS",
           oxphos_lvl = "levels",  redox_lvl = "levels")
ARM   <- c(oxphos_ppd = "OXPHOS subunits", redox_ppd = "redox",
           oxphos_lvl = "OXPHOS subunits", redox_lvl = "redox")

headhead <- dplyr::bind_rows(lapply(names(outcomes), function(on) {
  ry <- resid_design(outcomes[[on]])
  dplyr::bind_rows(lapply(names(axes), function(an) {
    f <- fit_one(outcomes[[on]], Z(axes[[an]]), ry)
    dplyr::bind_cols(tibble::tibble(outcome = on, axis = an,
                                    ruler = unname(RULER[an]),
                                    arm = unname(ARM[an])), f)
  }))
}))
headhead$padj_grid <- stats::p.adjust(headhead$p, method = "BH")
stopifnot(tibble::is_tibble(headhead))

# scale-invariance of the p-value, on the headline cell
chk <- fit_one(outcomes[[HEADLINE]], axes$oxphos_ppd)
chz <- headhead[headhead$outcome == HEADLINE & headhead$axis == "oxphos_ppd", ]
stopifnot(abs(chk$p - chz$p) < 1e-12)
message(sprintf("PART C: %d fits (%d outcomes x %d axes); p is scale-invariant (checked)",
                nrow(headhead), length(outcomes), length(axes)))

# =============================================================================
# PART D: THE PERMUTATION NULL -- script 43 PART B's scheme, all four axes
# =============================================================================
# Shuffle the axis WITHIN timepoint, refit, keep the interaction estimate. The
# model matrix is built by hand so 100,000 fits run in seconds; it is asserted
# against lm() on the observed data before the loop, because a hand-built design
# matrix with the columns in the wrong order would fail silently.
myc_d <- as.numeric(sm$myc == "pos")
Xof   <- function(a) cbind(1, myc_d, a, epi_comp, imm_comp, myc_d * a)
int_fast <- function(y, a) stats::.lm.fit(Xof(a), y)$coefficients[6L]

y_chk <- Z(outcomes[[HEADLINE]])
stopifnot(abs(int_fast(y_chk, axes$oxphos_ppd) -
              summary(stats::lm(y_chk ~ myc * a + epi + imm,
                                data.frame(y = y_chk, myc = sm$myc, epi = epi_comp,
                                           imm = imm_comp,
                                           a = axes$oxphos_ppd)))$coefficients["mycpos:a", 1]) < 1e-10)

tp_idx <- lapply(levels(sm$tp), function(tv) which(sm$tp == tv))
shuffle_within_tp <- function(a) { for (k in tp_idx) a[k] <- sample(a[k]); a }

set.seed(1)
perm <- dplyr::bind_rows(lapply(names(outcomes), function(on) {
  y <- Z(outcomes[[on]])
  dplyr::bind_rows(lapply(names(axes), function(an) {
    a   <- Z(axes[[an]])
    obs <- int_fast(y, a)
    nul <- vapply(seq_len(NPERM), function(i) int_fast(y, shuffle_within_tp(a)), numeric(1))
    tibble::tibble(outcome = on, axis = an, ruler = unname(RULER[an]),
                   observed = obs, null_median = stats::median(nul),
                   null_sd = stats::sd(nul),
                   percentile = 100 * mean(nul < obs), p_emp = mean(nul >= obs))
  }))
}))
message(sprintf("PART D: %d x %d permutations done", nrow(perm), NPERM))

# The mitoPPS arm cannot reproduce script 43's null EXACTLY -- the RNG stream
# differs -- but it must agree within Monte Carlo error, or the scheme has drifted.
rp  <- as.data.frame(ss$tradeoff_perm)
mm2 <- match(paste(perm$outcome, perm$axis), paste(rp$outcome, rp$axis))
ok  <- !is.na(mm2)
mc_dev <- max(abs(perm$p_emp[ok] - rp$p_emp[mm2[ok]]))
message(sprintf("PART D: mitoPPS arm agrees with script 43's null to max |dp_emp| = %.3f (MC error, not identity)",
                mc_dev))
stopifnot(mc_dev < 0.05)

# =============================================================================
# PART E: HOW TO READ THE ANSWER -- what each ruler's axis actually is
# =============================================================================
# Three diagnostics, because a null on the level ruler has two very different
# explanations and they must be told apart.
#
#   1. Are the two rulers even measuring the same thing? A share and a magnitude
#      need not agree, and if they do not then BOTH results can be true.
#   2. The AMBIENT: what does this axis correlate with by default on its own
#      ruler? script 43 measured the mitoPPS ambient; the level ruler has never
#      been measured this way.
#   3. The GLOBAL COMMON-MODE FACTOR (script 36). If the axis IS the global
#      factor, then it is not respiratory-chain-specific, whatever it returns.
zc <- function(x) as.numeric(scale(x))

M_ppd <- as.matrix(mps[, vapply(mps, is.numeric, logical(1)), drop = FALSE])
M_lvl <- vapply(names(path_ens_lvl), function(p) unname(lvl_of(p)), numeric(24L))
stopifnot(nrow(M_lvl) == 24L)

gm_ppd <- rowMeans(scale(M_ppd))
gm_lvl <- rowMeans(scale(M_lvl))

ambient_of <- function(a, M, drop_cols) {
  keep <- setdiff(colnames(M), drop_cols)
  rr <- abs(suppressWarnings(stats::cor(a, M[, keep, drop = FALSE], method = "spearman")))
  c(n = length(keep), med = stats::median(rr, na.rm = TRUE),
    q90 = unname(stats::quantile(rr, 0.90, na.rm = TRUE)))
}

diagnostics <- dplyr::bind_rows(lapply(names(axes), function(an) {
  is_ppd <- RULER[an] == "mitoPPS"
  M  <- if (is_ppd) M_ppd else M_lvl
  gm <- if (is_ppd) gm_ppd else gm_lvl
  amb <- ambient_of(axes[[an]], M, c(OX_P, RX_P))
  r_gm <- stats::cor(axes[[an]], gm)
  tibble::tibble(axis = an, ruler = unname(RULER[an]), arm = unname(ARM[an]),
                 n_pathways = unname(amb["n"]),
                 ambient_median_abs_rho = unname(amb["med"]),
                 ambient_q90 = unname(amb["q90"]),
                 r_global_factor = r_gm,
                 resid_frac = 1 - r_gm^2,
                 is_the_global_factor = (1 - r_gm^2) < RESID_FRAC_GATE)
}))

ruler_agreement <- tibble::tibble(
  arm = c("OXPHOS subunits", "redox"),
  r_pearson  = c(stats::cor(axes$oxphos_ppd, axes$oxphos_lvl),
                 stats::cor(axes$redox_ppd,  axes$redox_lvl)),
  rho_spearman = c(stats::cor(axes$oxphos_ppd, axes$oxphos_lvl, method = "spearman"),
                   stats::cor(axes$redox_ppd,  axes$redox_lvl,  method = "spearman")),
  rho_design_adj = c(stats::cor(resid_design(axes$oxphos_ppd),
                                resid_design(axes$oxphos_lvl), method = "spearman"),
                     stats::cor(resid_design(axes$redox_ppd),
                                resid_design(axes$redox_lvl),  method = "spearman")))

message(sprintf("PART E: the two OXPHOS rulers correlate r = %+.3f (Spearman %+.3f)",
                ruler_agreement$r_pearson[1], ruler_agreement$rho_spearman[1]))
for (i in seq_len(nrow(diagnostics)))
  message(sprintf("PART E: %-11s ambient %.2f | r(global) %+.3f | resid_frac %.3f%s",
                  diagnostics$axis[i], diagnostics$ambient_median_abs_rho[i],
                  diagnostics$r_global_factor[i], diagnostics$resid_frac[i],
                  ifelse(diagnostics$is_the_global_factor[i], "  <- IS the global factor", "")))

# =============================================================================
# PART F: THE VERDICT, against the rule fixed at the top of this script
# =============================================================================
hl <- function(ax) headhead[headhead$outcome == HEADLINE & headhead$axis == ax, ]
pm <- function(ax) perm[perm$outcome == HEADLINE & perm$axis == ax, ]

ref_sign <- sign(hl("oxphos_ppd")$myc_x_axis)
verdict_of <- function(ax) {
  r <- hl(ax)
  if (sign(r$myc_x_axis) != ref_sign) return("fails (opposite sign)")
  if (r$p < P_REPRO) return("reproduces")
  if (r$p < P_PART)  return("partial")
  "fails (p >= 0.20)"
}

verdict <- dplyr::bind_rows(lapply(names(axes), function(ax) {
  r <- hl(ax); q <- pm(ax); d <- diagnostics[diagnostics$axis == ax, ]
  tibble::tibble(
    axis = ax, ruler = unname(RULER[ax]), arm = unname(ARM[ax]),
    myc_x_axis = r$myc_x_axis, p = r$p, padj_grid = r$padj_grid,
    p_with_tp = r$p_with_tp, p_emp = q$p_emp, percentile = q$percentile,
    resid_frac = d$resid_frac,
    verdict = if (ax == "oxphos_ppd") "reference (script 43)" else verdict_of(ax),
    specificity = ifelse(d$is_the_global_factor,
                         "NOT arm-specific: the axis is the global expression factor",
                         "separable from the global factor"))
}))

message("\n", strrep("=", 78))
message("PART F: THE VERDICT on the headline outcome -- ", HEADLINE)
message(strrep("=", 78))
for (i in seq_len(nrow(verdict)))
  message(sprintf("  %-11s (%-7s) beta %+7.3f  p %.4f  p_emp %.3f  resid_frac %.3f  -> %s",
                  verdict$axis[i], verdict$ruler[i], verdict$myc_x_axis[i],
                  verdict$p[i], verdict$p_emp[i], verdict$resid_frac[i],
                  verdict$verdict[i]))

lvl_v   <- verdict$verdict[verdict$axis == "oxphos_lvl"]
lvl_sep <- verdict$resid_frac[verdict$axis == "oxphos_lvl"] >= RESID_FRAC_GATE
ppd_sep <- verdict$resid_frac[verdict$axis == "oxphos_ppd"] >= RESID_FRAC_GATE
agrees  <- grepl("^reproduces|^partial", lvl_v)

message("\nRECOMMENDATION FOR Fig. 2I's x axis:")
if (grepl("^reproduces", lvl_v) && lvl_sep) {
  message("  DRAW ABSOLUTE LEVELS. The interaction holds on a second quantifier whose")
  message("  units match the words in the hypothesis, and that axis has variance")
  message("  separable from the global expression factor. Clearer AND corroborated.")
} else if (agrees && !lvl_sep) {
  message("  KEEP mitoPPS -- and the reason is now MEASURED, not inherited.")
  message(sprintf("  The level ruler agrees in DIRECTION (beta %+.3f against %+.3f) but only",
                  verdict$myc_x_axis[verdict$axis == "oxphos_lvl"],
                  verdict$myc_x_axis[verdict$axis == "oxphos_ppd"]))
  message(sprintf("  %.1f%% of its variance is separable from the global expression factor",
                  100 * verdict$resid_frac[verdict$axis == "oxphos_lvl"]))
  message(sprintf("  (r = %+.3f), against %.1f%% for mitoPPS. On the absolute ruler the axis IS",
                  diagnostics$r_global_factor[diagnostics$axis == "oxphos_lvl"],
                  100 * verdict$resid_frac[verdict$axis == "oxphos_ppd"]))
  message("  'overall expression', not 'respiratory chain', so it cannot carry an")
  message("  ARM-SPECIFIC claim however it comes out. mitoPPS is the arm-specific ruler.")
  message("  CONSEQUENCE FOR THE TEXT: the transcriptomic claim is about respiratory")
  message("  PRIORITY (allocation), not about respiratory CAPACITY. A hypothesis")
  message("  sentence saying 'high respiratory state' is about a rate that no ruler")
  message("  here measures -- that is the perturbation experiment's job.")
} else if (agrees && lvl_sep) {
  message("  REPORT BOTH. The level ruler agrees in direction and is separable from the")
  message("  global factor, but is weaker. Either axis is defensible; name the ruler.")
} else {
  message("  KEEP mitoPPS. The interaction is RULER-DEPENDENT in SIGN, which belongs in")
  message("  the bounds -- the claim is about resource allocation, not about how much")
  message("  respiratory chain there is. The hypothesis sentence must be reworded.")
}
if (!ppd_sep)
  message("  WARNING: the mitoPPS axis is ALSO below the separability gate -- neither",
          "\n  ruler supports an arm-specific reading.")
message(strrep("=", 78), "\n")

# =============================================================================
# NOTES
# =============================================================================
notes <- c(
  "THE QUESTION. Fig. 2I's x axis is the OXPHOS-subunit mitoPPS score, a budget",
  "SHARE inside the mitochondrial compartment. The closing hypothesis is about a",
  "high respiratory STATE, which is a magnitude. This script refits the identical",
  "model with the identical outcome and covariates and changes ONLY the axis.",
  "",
  "THE RULE WAS FIXED BEFORE THE RESULT WAS SEEN:",
  paste0("  reproduces -- ", VERDICT_RULE[["reproduces"]]),
  paste0("  partial    -- ", VERDICT_RULE[["partial"]]),
  paste0("  fails      -- ", VERDICT_RULE[["fails"]]),
  "and BOTH rulers are reported whichever way it comes out. A ruler-dependence is",
  "a finding: CLAUDE.md already requires that a sentence quoting an OXPHOS",
  "abundance number names its ruler (the same 87 subunits give +0.061 unweighted,",
  "+0.201 expression-weighted and +0.226 as a summed-count ratio on the diagonal).",
  "",
  "WHY A NULL ON THE LEVEL RULER WOULD NOT BE A REFUTATION. mitoPPS is a pairwise",
  "ratio, so the global common-mode axis that defeats every score-space coupling",
  "largely cancels; an absolute level does not cancel it. PART E measures",
  "resid_frac = 1 - r(axis, global factor)^2 on each ruler. Below",
  paste0("  ", RESID_FRAC_GATE, " the axis IS the global expression factor (script 36's"),
  "readability gate) and whatever it returns is not a statement about the",
  "respiratory chain specifically. Read the verdict and the specificity column",
  "together, never the verdict alone.",
  "",
  "WHAT IS HELD FIXED. Outcome, covariates (epithelial and immune composition),",
  "model form, sample set, and the within-timepoint permutation scheme are script",
  "43 PART B's, reproduced to 1e-8 in PART A before anything is changed. The five",
  "outcomes are script 43's five, so the multiplicity context is preserved: the",
  "bounds on Fig. 2I say the two axes were tested against five outcomes and only",
  "one reached p < 0.05.",
  "",
  "ONE ASYMMETRY THAT CANNOT BE REMOVED. The `priming_ppd (PRO-ANTI)` outcome is",
  "itself built from mitoPPS scores, so on the mitoPPS ruler that pairing shares a",
  "normalisation the level ruler does not. The HEADLINE outcome does not: the",
  "PUMA:Bcl-xL ratio is two log counts and is ruler-independent, which is why the",
  "verdict is read off it and not off the grid.",
  "",
  "INFERENCE WEIGHT. Ranking, not confirmation. n = 6 per cell, and an interaction",
  "is a difference of two 6-versus-6 contrasts -- the least powered quantity in the",
  "design. The reference result itself (p = 0.0052) sits at the 91.7th percentile",
  "of its own permutation null. Nothing here converts a lead into a result; it",
  "establishes whether the lead is a property of the compartment's ALLOCATION or of",
  "its SIZE, which is what the figure has to decide.",
  "",
  "SCOPE. No DESeq2 refit, no gene set rebuilt, no contrast touched.")

out <- list(
  control         = control,          # PART A, the positive control
  headhead        = headhead,         # PART C, 5 outcomes x 4 axes
  perm            = perm,             # PART D, the permutation nulls
  diagnostics     = diagnostics,      # PART E, ambient + global factor per axis
  ruler_agreement = ruler_agreement,  # PART E, do the two rulers agree at all
  verdict         = verdict,          # PART F
  axes            = as.data.frame(axes),
  defs = list(n_perm = NPERM, headline = HEADLINE, verdict_rule = VERDICT_RULE,
              p_reproduces = P_REPRO, p_partial = P_PART,
              resid_frac_gate = RESID_FRAC_GATE,
              n_genes = c(oxphos = n_ox, redox = n_rx),
              reproduction = c(max_abs_dbeta = dev_beta, max_abs_dp = dev_p,
                               max_abs_dlevel = dev_lvl, perm_mc_dev = mc_dev)),
  analysis_date = Sys.Date(),
  notes = notes)

saveRDS(out, here::here("results", "axis_ruler_test.rds"))
message("46: wrote results/axis_ruler_test.rds")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  res <- readRDS(here::here("results", "axis_ruler_test.rds"))

  ## THE ANSWER, in one table
  res$verdict %>% print(width = Inf)

  ## the headline outcome on both rulers, side by side
  res$headhead[res$headhead$outcome == res$defs$headline, ] %>% print(width = Inf)

  ## are the two rulers even the same variable?
  res$ruler_agreement %>% print()

  ## what each axis IS -- read this before believing any verdict
  res$diagnostics %>% print(width = Inf)

  ## the full grid, ranked. Nothing here is expected to survive BH; the ordering is
  ## the readable quantity, as it is everywhere else in this corpus.
  res$headhead[order(res$headhead$p), c("outcome", "axis", "ruler",
                                        "myc_x_axis", "p", "padj_grid")] %>%
    print(n = 20)

  ## the two OXPHOS axes drawn against each other -- a share against a magnitude.
  ## If this is a tight line the rulers are interchangeable; if it is a cloud they
  ## are measuring different things and both results can be true.
  with(res$axes, plot(oxphos_ppd, oxphos_lvl,
                      xlab = "OXPHOS subunits, mitoPPS (budget share)",
                      ylab = "OXPHOS subunits, log2 summed normalised counts",
                      pch = 19))

  ## the picture Fig. 2I would become, on the level ruler
  dds <- readRDS(here::here("results", "dds_int_run.rds"))
  smx <- as.data.frame(SummarizedExperiment::colData(dds))
  L   <- log2(DESeq2::counts(dds, normalized = TRUE) + 1)
  cdf <- readRDS(here::here("results", "combined_df_annotated.rds"))
  eo  <- function(s) cdf$gene[match(s, cdf$mgi_symbol)]
  ry  <- as.numeric(L[eo("Bbc3"), ] - L[eo("Bcl2l1"), ])
  cl  <- ifelse(smx$myc_status == "pos", "#D55E00", "#0072B2")
  plot(res$axes$oxphos_lvl, ry, col = cl, pch = 19,
       xlab = "OXPHOS subunits, log2 summed normalised counts",
       ylab = "Bbc3 - Bcl2l1 (log2)")
  for (g in c("neg", "pos")) {
    k <- smx$myc_status == g
    abline(stats::lm(ry[k] ~ res$axes$oxphos_lvl[k]),
           col = ifelse(g == "pos", "#D55E00", "#0072B2"))
  }

  ## the permutation null for the headline cell, both rulers
  res$perm[res$perm$outcome == res$defs$headline, ] %>% print(width = Inf)

  ## and the reproduction evidence, if PART A's assertions are ever questioned
  res$defs$reproduction %>% print()
}
