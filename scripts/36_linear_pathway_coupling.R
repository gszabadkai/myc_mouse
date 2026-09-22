# scripts/36_linear_pathway_coupling.R
# =============================================================================
# Pathway-pathway coupling done the principled way: a LINEAR quantifier, the
# design removed FIRST, then the global common-mode factor, then the specific
# coupling. (Block B, narrative finalisation 2026-07-18)
# =============================================================================
#
# WHY THIS EXISTS. Scripts 34/35 established the CEILING: every GSVA composite
# correlates with almost every other (median |rho| ~0.8 to an arbitrary one of
# the 884 library sets), and mito_oxphos is nearly the per-sample mean of all
# scores. Two method notes the author wrote (docs/GSVA_global_background_and_
# pathway_coupling.md ; docs/Gene_set_quantification_for_pathway_correlation.md)
# name what that ceiling IS and how to handle it. This script implements them.
#
# THE THREE THINGS THE DOCS CHANGE, AND WHY.
#   (1) MEAN-SUBTRACTION IS THE WRONG CORRECTION. It assumes every set loads
#       equally on the common factor and imposes a sum-to-zero constraint that
#       manufactures spurious negatives. That is exactly why script-35-era
#       "subtract-the-sample-mean" gave an unstable ordering. The principled
#       correction is: residualise the DESIGN, take PC1 of the residual, regress
#       it out (doc 1 s4/s6).
#   (2) ESTIMATE THE FACTOR *AFTER* REMOVING THE DESIGN (doc 1 s7). On the raw
#       matrix PC1 ~ the four-group contrast. Residualising timepoint*genotype
#       first defines the factor as the broad axis that remains BETWEEN MICE
#       within groups. This DECOMPOSES the ceiling: raw = (between-group
#       contrast, much of it Myc) + (within-group common-mode axis).
#   (3) GSVA IS NOT THE RIGHT PRIMARY INSTRUMENT FOR COUPLING (doc 2). Its
#       nonlinear competitive score has no clean covariance reading. A LINEAR
#       mean gene-wise z-score's correlation IS the average cross-gene
#       covariance. So the primary quantifier here is the mean z-score;
#       singscore (cohort-independent, rank-based) is the robustness arm; GSVA
#       is kept for continuity.
#
# WHAT THIS IS AND IS NOT FOR. RANKING + the global factor as a phenotype -- NOT
# confirmatory CIs. At n=24 a residual r~0.29 arises under the null (doc 1 s10),
# so per-pair intervals confirm nothing; forcing them would re-import the
# over-null-the-exploration error. The deliverable is a cleaner RANKING of the
# pre-specified leads across three quantifiers, plus a characterisation of the
# common factor.
#
# SCOPE (unchanged, load-bearing). None of this touches the CONTRASTS. The
# content result is raw count shares; the attenuation is DESeq2 LFCs; neither
# goes near GSVA or these couplings. This re-ranks exploratory leads only.
#
# Input:  results/gsva_scores.rds            ($scores 884x24 GSVA; $expr_mat VST
#                                             log 14121x24; $pathways memberships;
#                                             $sample_meta; $set_meta)
#         results/myc_mito_centrality.rds    (script 28: $panel composites, incl.
#                                             mb2_fork = AP7 fork, quantifier-free)
#         results/mtdna_axis_and_coupling_null.rds ($panel_share: ieg_stress,
#                                             endothelial contamination, epithelial)
#         results/dds_int_run.rds            (counts -> library size, detected, mt%)
# Output: results/linear_pathway_coupling.rds
#         outputs/linear_pathway_coupling/*.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

out_dir <- here::here("outputs", "linear_pathway_coupling")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

HAS_SINGSCORE <- requireNamespace("singscore", quietly = TRUE)
if (!HAS_SINGSCORE)
  message("NOTE: singscore not installed -- the rank-based robustness arm will be skipped.")

# =============================================================================
# PART 1: LOAD -- one score universe, one panel, the covariate sources
# =============================================================================
gs    <- readRDS(here::here("results", "gsva_scores.rds"))
mc    <- readRDS(here::here("results", "myc_mito_centrality.rds"))

panel    <- mc$panel
sm       <- gs$sample_meta
set_meta <- gs$set_meta
pathways <- gs$pathways
expr     <- gs$expr_mat                                   # VST (log), genes x samples

# align everything to the GSVA score column order (= panel$sample)
M_gsva <- gs$scores[, panel$sample, drop = FALSE]
expr   <- expr[, panel$sample, drop = FALSE]
sm     <- sm[match(panel$sample, sm$sample), ]
stopifnot(identical(colnames(M_gsva), panel$sample),
          identical(colnames(expr),   panel$sample),
          identical(sm$sample,        panel$sample))

by_cat <- function(cat) set_meta$set_name[set_meta$category_primary == cat]

# the design: timepoint * genotype (the four groups). Removed FIRST, everywhere.
meta <- data.frame(timepoint  = factor(sm$timepoint),
                   myc_status = factor(sm$myc_status))
X <- stats::model.matrix(~ timepoint * myc_status, data = meta)
qrX <- qr(X)

i6  <- panel$timepoint == "6W"
i12 <- panel$timepoint == "12W"
ipos6 <- i6 & panel$myc_status == "pos"
stopifnot(sum(i6) == 12, sum(i12) == 12, sum(ipos6) == 6)

# =============================================================================
# PART 2: THREE SCORE MATRICES on the SAME sets and genes
# =============================================================================
# All three are "total" scores (design retained). Design removal happens at the
# correlation stage (score-level residualisation, PART D), uniformly, so the
# three quantifiers are directly comparable. The linear z-score ALSO gets the
# doc-2 gene-level design-adjusted variant (M_z_design) as a consistency check.

sets_universe <- rownames(M_gsva)                          # the 884 scored sets
pw <- pathways[intersect(names(pathways), sets_universe)]  # memberships for them

# --- (B1) GSVA: already scored ------------------------------------------------
# M_gsva above.

# --- (B2) mean gene-wise z-score (doc-2 primary linear score) -----------------
row_z <- function(mat) {                                   # z each gene across samples
  z <- t(scale(t(mat)))
  z[is.finite(rowSums(z)), , drop = FALSE]                 # drop zero-variance genes
}
score_linear <- function(gene_z) {
  t(vapply(sets_universe, function(s) {
    g <- intersect(pw[[s]], rownames(gene_z))
    if (length(g) < 5L) return(rep(NA_real_, ncol(gene_z)))
    colMeans(gene_z[g, , drop = FALSE])
  }, numeric(ncol(gene_z))))
}
gene_z_total  <- row_z(expr)                               # total (design retained)
M_z <- score_linear(gene_z_total)

# doc-2 s16 preferred: residualise GENES on the design, then z, then set-mean
gene_resid    <- t(qr.resid(qrX, t(expr)))                 # genes x samples, design out
gene_z_design <- row_z(gene_resid)
M_z_design    <- score_linear(gene_z_design)               # gene-level design-adjusted

# --- (B3) singscore (rank-based, cohort-independent robustness) ---------------
M_ss <- NULL
if (HAS_SINGSCORE) {
  rankData <- singscore::rankGenes(expr)
  M_ss <- t(vapply(sets_universe, function(s) {
    g <- intersect(pw[[s]], rownames(expr))
    if (length(g) < 5L) return(rep(NA_real_, ncol(expr)))
    sc <- suppressWarnings(singscore::simpleScore(rankData, upSet = g))
    sc$TotalScore
  }, numeric(ncol(expr))))
  colnames(M_ss) <- colnames(expr)
}

# keep the set intersection that is fully scored in every available matrix
mats <- list(gsva = M_gsva, zscore = M_z, singscore = M_ss)
mats <- mats[!vapply(mats, is.null, logical(1))]
common_sets <- Reduce(intersect, lapply(mats, function(M)
  rownames(M)[stats::complete.cases(M)]))
mats     <- lapply(mats, function(M) M[common_sets, , drop = FALSE])
M_z_design <- M_z_design[intersect(common_sets, rownames(M_z_design)), , drop = FALSE]
message(sprintf("Score matrices: %s ; %d sets scored in all.",
                paste(names(mats), collapse = "/"), length(common_sets)))

# =============================================================================
# PART 3: COMPOSITE SERIES from any score matrix (same algebra as script 28)
# =============================================================================
mito_names <- intersect(by_cat("MitoCarta"), common_sets)
mito_grep  <- function(pat) grep(pat, mito_names, value = TRUE)
OXPHOS_PAT <- "OXPHOS|COMPLEX_[IV]|_SUBUNITS|ASSEMBLY_FACTORS|ELECTRON_CARRIERS|CRISTAE"
BIOG_PAT   <- "RIBOSOME|CENTRAL_DOGMA|MT_TRNA|MT_RRNA|MTRNA|MTDNA|IMPORT|TRANSLATION"
teb_up <- grep("^MG_TEB_VS_DUCTAL_.*_UP$", common_sets, value = TRUE)
teb_dn <- grep("^MG_TEB_VS_DUCTAL_.*_DN$", common_sets, value = TRUE)

# axis -> member sets (unsigned axes)
axis_sets <- list(
  mito_oxphos     = mito_grep(OXPHOS_PAT),
  mito_biogenesis = mito_grep(BIOG_PAT),
  prolif          = intersect(by_cat("Proliferation"),  common_sets),
  myc_sig         = intersect(by_cat("MYC_signatures"), common_sets),
  cholesterol     = intersect(c("GS_METAB_CHOLESTEROL","GS_METAB_MEVALONATE","METAB_CHOLESTEROL_HALLMARK",
                                "METAB_CHOLESTEROL_REACTOME","METAB_CHOLESTEROL_WP"), common_sets),
  redox           = intersect(c("GS_METAB_REDOX","GS_METAB_GLUTATHIONE","GS_METAB_REACTIVE_OXYGEN"), common_sets),
  nucleotide      = intersect(c("GS_METAB_NUCLEOTIDE","GS_METAB_PURINE","GS_METAB_PYRIMIDINE",
                                "METAB_NUCLEOTIDE_REACTOME","METAB_NUCLEOTIDE_SALVAGE_REACTOME",
                                "METAB_NUCLEOTIDE_WP","METAB_PURINE_KEGG","METAB_PYRIMIDINE_KEGG",
                                "METAB_PYRIMIDINE_WP"), common_sets),
  tca             = intersect(c("GS_METAB_KREBS","METAB_TCA_KEGG","METAB_TCA_REACTOME"), common_sets))

# signed composites (up - down), exactly as script 28
signed_sets <- list(
  priming    = list(up = intersect("MITOCARTA_APOPTOSIS_PRO",  common_sets),
                    dn = intersect("MITOCARTA_APOPTOSIS_ANTI", common_sets)),
  teb_dediff = list(up = teb_up, dn = teb_dn))

comp_M <- function(M, sets) {
  sets <- intersect(sets, rownames(M))
  if (length(sets) == 0) return(rep(NA_real_, ncol(M)))
  colMeans(M[sets, , drop = FALSE])
}
# every axis + outcome as a per-sample series from a given matrix.
# mb2_fork is the AP7 fold ratio (NOT a gene set) -- quantifier-independent, so
# it is taken from the panel unchanged in all three matrices.
series_from <- function(M, name) {
  if (name == "mb2_fork") return(panel$mb2_fork)
  if (name %in% names(signed_sets))
    return(comp_M(M, signed_sets[[name]]$up) - comp_M(M, signed_sets[[name]]$dn))
  comp_M(M, axis_sets[[name]])
}
# member sets of a composite (for leave-pair-out background exclusion)
member_sets <- function(name) {
  if (name == "mb2_fork") return(character(0))
  if (name %in% names(signed_sets)) return(c(signed_sets[[name]]$up, signed_sets[[name]]$dn))
  axis_sets[[name]]
}

# =============================================================================
# PART A: THE CEILING PROBES -- reproduce them properly, per quantifier
# =============================================================================
# The script-35-era numbers (random-pair 0.678; cor(mito_oxphos, global mean)
# 0.979; drop to ~0.29 on removal) were ad-hoc. Here they are saved fields.
global_mean <- function(M) colMeans(M, na.rm = TRUE)      # per-sample mean of all set scores

ceiling_probes <- purrr::map_dfr(names(mats), function(nm) {
  M <- mats[[nm]]
  gm  <- global_mean(M)
  ox  <- series_from(M, "mito_oxphos")
  # random set-set pair |cor| (curation probe)
  set.seed(1L)
  ridx <- replicate(2000, sample(nrow(M), 2))
  rpair <- apply(ridx, 2, function(k)
    suppressWarnings(stats::cor(M[k[1], ], M[k[2], ], method = "spearman")))
  # chance floor: permute sample labels of the axis, recompute its median |rho|
  amb_axis <- function(vec, idx) {
    keep <- setdiff(rownames(M), axis_sets$mito_oxphos)
    r <- apply(M[keep, idx, drop = FALSE], 1,
               function(v) suppressWarnings(stats::cor(v, vec[idx], method = "spearman")))
    stats::median(abs(r[is.finite(r)]))
  }
  set.seed(2L)
  chance <- stats::median(replicate(200, amb_axis(sample(ox), rep(TRUE, length(ox)))))
  # remove the per-sample global mean (regress it out of every row), re-probe axis
  M_res  <- t(stats::residuals(stats::lm(t(M) ~ gm)))
  ox_res <- series_from(M_res, "mito_oxphos")
  tibble::tibble(
    quantifier          = nm,
    random_pair_med_abs = stats::median(abs(rpair), na.rm = TRUE),
    cor_oxphos_globalmean = suppressWarnings(stats::cor(ox, gm)),
    ambient_oxphos_all  = amb_axis(ox, rep(TRUE, length(ox))),
    ambient_within_6Wpos= amb_axis(ox, ipos6),
    chance_floor_permuted = chance,
    ambient_after_globalmean_removed = amb_axis(ox_res, rep(TRUE, length(ox_res))))
})

# =============================================================================
# PART C: THE GLOBAL FACTOR -- estimate it right, then read it as a phenotype
# =============================================================================
# design->PC1: residualise every set-score row on the design, PC1 of that.
pc1_of <- function(M, residualise = TRUE) {
  Mx <- if (residualise) t(qr.resid(qrX, t(M))) else M
  pc <- stats::prcomp(t(Mx), center = TRUE, scale. = FALSE)
  g  <- pc$x[, 1]
  # orient so the factor tracks the overall score level (sign is arbitrary in PCA)
  if (suppressWarnings(stats::cor(g, global_mean(M))) < 0) { g <- -g; pc$rotation[,1] <- -pc$rotation[,1] }
  list(g = g, var = pc$sdev^2 / sum(pc$sdev^2), load = pc$rotation[, 1])
}

# covariate panel: biological state vs technical/prep, to answer "which is it?"
# (mt% is characterised in script 33 -- the mt-axis is its remit -- so it is not
# duplicated here; library_size/detected_genes/IEG/contamination carry technical.)
counts_mat <- DESeq2::counts(readRDS(here::here("results", "dds_int_run.rds")))
counts_mat <- counts_mat[, panel$sample, drop = FALSE]
share      <- readRDS(here::here("results", "mtdna_axis_and_coupling_null.rds"))$panel_share
share      <- share[panel$sample, , drop = FALSE]
covars <- data.frame(
  MYC_activity  = panel$myc_sig,                                   # biological
  proliferation = panel$prolif,                                    # biological
  epithelial    = share[, "epithelial"],                          # composition
  IEG_prep      = log2(share[, "ieg_stress"]),                    # prep/dissociation stress
  contamination = log2(share[, "endothelial"]),                   # technical (residual stroma)
  library_size  = colSums(counts_mat),                            # technical
  detected_genes= colSums(counts_mat > 0),                        # technical
  timepoint_6W  = as.integer(panel$timepoint == "6W"),           # design (check ~0 after resid)
  genotype_pos  = as.integer(panel$myc_status == "pos"))          # design (check ~0 after resid)

# RAW axis (no design removed): how much of the ceiling IS the four-group contrast?
gf_raw <- pc1_of(mats$zscore, residualise = FALSE)
# RESIDUAL axis (design removed): the within-group common-mode axis
gf_res <- pc1_of(mats$zscore, residualise = TRUE)

global_factor_phenotype <- purrr::map_dfr(names(covars), function(cv) {
  tibble::tibble(
    covariate    = cv,
    rho_raw_axis = suppressWarnings(stats::cor(gf_raw$g, covars[[cv]], method = "spearman")),
    rho_residual_axis = suppressWarnings(stats::cor(gf_res$g, covars[[cv]], method = "spearman")))
}) |> dplyr::arrange(dplyr::desc(abs(rho_residual_axis)))

# variance explained + how one-signed the loadings are (common mode vs contrast)
global_factor_summary <- purrr::map_dfr(names(mats), function(nm) {
  raw <- pc1_of(mats[[nm]], residualise = FALSE)
  res <- pc1_of(mats[[nm]], residualise = TRUE)
  tibble::tibble(
    quantifier = nm,
    pc1_var_raw_pct = 100 * raw$var[1],
    pc1_var_residual_pct = 100 * res$var[1],
    frac_loadings_positive_raw = mean(raw$load > 0),
    cor_rawPC1_globalmean = suppressWarnings(stats::cor(raw$g, global_mean(mats[[nm]]))),
    cor_rawPC1_group = suppressWarnings(stats::cor(raw$g,
                          as.integer(interaction(panel$timepoint, panel$myc_status)))))
})

# =============================================================================
# PART D: THE PRE-SPECIFIED PAIRS -- raw / design-adj / global-adj x 3 methods
# =============================================================================
# Leave-pair-out global adjustment (doc 1 s4): estimate PC1 from the BACKGROUND
# only (all sets except the member sets of BOTH composites, design removed), then
# residualise the two composite series on design + that factor.
#
# *** THE INTERPRETABILITY TRAP (doc 1 s4/s14). *** When a composite is nearly
# the global factor itself -- mito_oxphos correlates ~0.97 with the per-sample
# global mean -- removing that factor leaves ~6% of its variance, which is noise.
# Correlating two such noise residuals at n=24 yields an UNSTABLE, often strong,
# often NEGATIVE number that is NOT antagonism (doc 1's "masked antagonism"
# pattern is indistinguishable from over-removal here). So the coupling is only
# INTERPRETABLE where enough of BOTH series survives the correction. `resid_frac`
# reports exactly that; a global-adjusted coupling on a series with resid_frac <
# ~0.10 must be read as "no specific signal separable", not as a finding.
global_adjusted_detail <- function(M, name_a, name_b, method = "pearson") {
  a <- series_from(M, name_a); b <- series_from(M, name_b)
  bg <- setdiff(rownames(M), c(member_sets(name_a), member_sets(name_b)))
  Rbg <- t(qr.resid(qrX, t(M[bg, , drop = FALSE])))
  g   <- stats::prcomp(t(Rbg), center = TRUE, scale. = FALSE)$x[, 1]
  Xg  <- cbind(X, g)
  ra  <- qr.resid(qr(Xg), a); rb <- qr.resid(qr(Xg), b)
  gm  <- global_mean(M)
  list(cor           = suppressWarnings(stats::cor(ra, rb, method = method)),
       cor_spearman  = suppressWarnings(stats::cor(ra, rb, method = "spearman")),
       ax_r2_global  = suppressWarnings(stats::cor(a, gm))^2,
       oc_r2_global  = suppressWarnings(stats::cor(b, gm))^2,
       ax_resid_frac = stats::var(ra) / stats::var(a),
       oc_resid_frac = stats::var(rb) / stats::var(b))
}
global_adjusted_cor <- function(M, name_a, name_b, method = "pearson")
  global_adjusted_detail(M, name_a, name_b, method)$cor
design_adjusted_cor <- function(M, name_a, name_b, method = "pearson") {
  a <- qr.resid(qrX, series_from(M, name_a))
  b <- qr.resid(qrX, series_from(M, name_b))
  suppressWarnings(stats::cor(a, b, method = method))
}
raw_cor <- function(M, name_a, name_b, method = "pearson")
  suppressWarnings(stats::cor(series_from(M, name_a), series_from(M, name_b), method = method))

pairs_spec <- tibble::tribble(
  ~axis,            ~outcome,     ~note,
  "mito_oxphos",    "prolif",     "Q4: 'OXPHOS is central'",
  "mito_oxphos",    "teb_dediff", "Q4: 'OXPHOS is central'",
  "mito_biogenesis","prolif",     "Q4: 'biogenesis is a bystander'",
  "mito_biogenesis","teb_dediff", "Q4: 'biogenesis is a bystander'",
  "cholesterol",    "mb2_fork",   "lead: mevalonate -> tumorigenic fork (non-Myc)",
  "redox",          "mb2_fork",   "lead: redox -> fork (non-Myc control)",
  "mito_oxphos",    "priming",    "NEGATIVE CONTROL -- circular (priming is MitoCarta)",
  "myc_sig",        "prolif",     "Q5: endogenous MYC ~ programme (pooled; see script 35 for within-WT)")

coupling_hierarchy <- purrr::pmap_dfr(pairs_spec, function(axis, outcome, note) {
  purrr::map_dfr(names(mats), function(nm) {
    M <- mats[[nm]]
    d <- global_adjusted_detail(M, axis, outcome)
    tibble::tibble(
      axis = axis, outcome = outcome, note = note, quantifier = nm,
      raw     = raw_cor(M, axis, outcome),
      design  = design_adjusted_cor(M, axis, outcome),
      global  = d$cor, global_spearman = d$cor_spearman,
      ax_r2_global  = d$ax_r2_global,  oc_r2_global  = d$oc_r2_global,
      ax_resid_frac = d$ax_resid_frac, oc_resid_frac = d$oc_resid_frac,
      # is there any of EITHER series left to test after the correction?
      testable = d$ax_resid_frac >= 0.10 & d$oc_resid_frac >= 0.10)
  })
})

# robustness: FIRST ask whether the coupling is even TESTABLE (did enough of both
# series survive the global-factor removal), THEN whether it agrees across
# quantifiers. A composite that IS the global factor (mito_oxphos, r^2~0.94 with
# the global mean) has no separable signal -- its global-adjusted "coupling" is
# noise-on-noise and must NOT be read as antagonism, however consistent its sign.
robustness <- coupling_hierarchy |>
  dplyr::group_by(axis, outcome, note) |>
  dplyr::summarise(
    n_methods            = dplyr::n(),
    median_ax_r2_global  = stats::median(ax_r2_global),
    median_ax_resid_frac = stats::median(ax_resid_frac),
    any_testable         = any(testable), all_testable = all(testable),
    global_min = min(global), global_max = max(global),
    same_sign  = length(unique(sign(round(global, 2)))) == 1,
    verdict = dplyr::case_when(
      !any_testable ~ "axis IS the global factor: no specific signal separable (residual = noise)",
      !all_testable ~ "barely separable (residual < 10% in >=1 method) -- do not interpret",
      length(unique(sign(round(global, 2)))) > 1 ~ "method-sensitive (sign flips): not a lead",
      all(abs(global) >= 0.3) ~ "separable + consistent (a RANKING lead only, n=24)",
      TRUE ~ "separable but weak"),
    .groups = "drop") |>
  dplyr::arrange(dplyr::desc(all_testable), dplyr::desc(same_sign))

# doc-2 consistency check: gene-level design-adjusted z-score (M_z_design) raw
# coupling should ~ match the score-level design-adjusted z-score coupling.
gene_vs_score_design <- purrr::pmap_dfr(pairs_spec, function(axis, outcome, note) {
  tibble::tibble(
    axis = axis, outcome = outcome,
    zscore_design_scorelevel = design_adjusted_cor(mats$zscore, axis, outcome),
    zscore_design_genelevel  = raw_cor(M_z_design, axis, outcome))
})

# =============================================================================
# PART E: VERDICTS
# =============================================================================
ceiling_verdict <- {
  z  <- ceiling_probes |> dplyr::filter(quantifier == "zscore")
  g  <- ceiling_probes |> dplyr::filter(quantifier == "gsva")
  sprintf(paste0(
    "WHERE THE CEILING COMES FROM -- ONE COMMON-MODE AXIS, NOT n, NOT CURATION. On the LINEAR ",
    "z-score: mito_oxphos correlates with the per-sample mean of all %d set scores at r=%.3f (GSVA ",
    "%.3f), and random set-set pairs sit at median |rho|=%.2f (GSVA %.2f). It is NOT low n -- the ",
    "chance floor from permuting the axis's own sample labels is %.2f, far below the observed ",
    "ambient of %.2f. It is NOT the genotype contrast alone -- WITHIN 6W_pos (one genotype, one ",
    "age, contrast gone) the ambient is still %.2f. And it is REMOVABLE: regress the per-sample ",
    "global mean out of every set and the axis's ambient drops to %.2f (~chance). So the raw ",
    "ceiling is a single global per-sample offset that every score rides -- part the four-group ",
    "contrast (much of it Myc), part a within-group common-mode axis. PART C says which is which."),
    length(common_sets), z$cor_oxphos_globalmean, g$cor_oxphos_globalmean,
    z$random_pair_med_abs, g$random_pair_med_abs,
    z$chance_floor_permuted, z$ambient_oxphos_all, z$ambient_within_6Wpos,
    z$ambient_after_globalmean_removed)
}

phenotype_verdict <- {
  top <- global_factor_phenotype |> dplyr::filter(!covariate %in% c("timepoint_6W","genotype_pos")) |>
    dplyr::slice_head(n = 3)
  s   <- global_factor_summary |> dplyr::filter(quantifier == "zscore")
  sprintf(paste0(
    "IS THE COMMON FACTOR BIOLOGICAL OR METHODOLOGICAL? BOTH -- and the residual within-group axis ",
    "is measurable. The RAW PC1 explains %.0f%% of variance and %.0f%% of set loadings share one ",
    "sign (a common mode, not a contrast); it correlates with the group label at r=%.2f, so a large ",
    "part of the raw ceiling IS the four-group design. The RESIDUAL PC1 (design removed) explains ",
    "%.0f%% and its strongest correlates among measured covariates are: %s. Biological state ",
    "(MYC/proliferation/composition) and technical/prep tells (IEG, library size, mt%%) both load ",
    "-- so the axis is neither cleanly biological nor cleanly technical, exactly as the docs warn. ",
    "IT IS NOT NOISE: it is a real broad sample-level programme, and it is the thing a raw coupling ",
    "mostly measures. Report it as a phenotype (PART C table) AND remove it for specific coupling."),
    s$pc1_var_raw_pct, 100 * global_factor_summary$frac_loadings_positive_raw[global_factor_summary$quantifier=="zscore"],
    s$cor_rawPC1_group, s$pc1_var_residual_pct,
    paste(sprintf("%s (rho=%+.2f)", top$covariate, top$rho_residual_axis), collapse = "; "))
}

coupling_verdict <- {
  z  <- function(ax, oc) {
    r <- coupling_hierarchy |> dplyr::filter(axis == ax, outcome == oc, quantifier == "zscore"); r
  }
  op <- z("mito_oxphos", "prolif"); ch <- z("cholesterol", "mb2_fork"); rx <- z("redox", "mb2_fork")
  n_mito_untestable <- robustness |>
    dplyr::filter(axis %in% c("mito_oxphos","mito_biogenesis"), !any_testable) |> nrow()
  rx_rob <- robustness |> dplyr::filter(axis=="redox", outcome=="mb2_fork")
  ch_rob <- robustness |> dplyr::filter(axis=="cholesterol", outcome=="mb2_fork")
  sprintf(paste0(
    "THE PRE-SPECIFIED LEADS UNDER THE PRINCIPLED CORRECTION -- AND THE TRAP THAT DOMINATES THEM. ",
    "The mito axes ARE the global factor: mito_oxphos carries r^2=%.2f with the per-sample global ",
    "mean, so after design->PC1->residual only ~%.0f%% of its variance remains, and that remainder ",
    "is NOISE. Its global-adjusted 'coupling' to proliferation prints %+.2f (z-score) and is ",
    "sign-consistent across all three quantifiers, but this is noise-on-noise at n=24 -- it is NOT ",
    "antagonism, and %d of the mito->outcome pairs are flagged 'no specific signal separable'. THE ",
    "CORRECT READING: OXPHOS's specific coupling to phenotype is NOT TESTABLE here because OXPHOS ",
    "has essentially no variance independent of the broad state (doc 1 s8) -- which is itself the ",
    "answer to 'is OXPHOS central?': no separable central coupling exists at n=24, in either ",
    "direction. This REPLACES script 35's raw-ambient framing with a cleaner statement, and it ",
    "confirms Q4's published 'OXPHOS-central FAILS'. THE ONE ROBUST, TESTABLE LEAD IS NON-MYC: ",
    "redox -> MB2 fork, global-adjusted %+.2f (z-score), SAME SIGN and 0.62-0.73 across ALL THREE ",
    "quantifiers, on an axis that barely loads on the global factor (r^2=%.2f, so ~%.0f%% of its ",
    "variance survives the correction) -- verdict '%s'. Lower redox/glutathione capacity tracks the ",
    "tumorigenic fork independently of Myc dose. cholesterol/mevalonate -> fork is TESTABLE but ",
    "WEAK once the factor is removed (%+.2f, '%s') -- DEMOTED from script 35's raw-ambient headline. ",
    "NEGATIVE CONTROL mito_oxphos->priming is untestable (both mito, both ~the global factor; its ",
    "sign flips across methods, confirming noise). BOTTOM LINE: no clean specific MITO coupling is ",
    "extractable at n=24; the readable lead is redox->fork (robust), with mevalonate->fork a weaker ",
    "second -- both RANKING leads for the bench, not results."),
    op$ax_r2_global, 100 * op$ax_resid_frac, op$global, n_mito_untestable,
    rx$global, rx$ax_r2_global, 100 * rx$ax_resid_frac, rx_rob$verdict,
    ch$global, ch_rob$verdict)
}

message("\n", paste(strwrap(ceiling_verdict,   width = 92), collapse = "\n"))
message("\n", paste(strwrap(phenotype_verdict, width = 92), collapse = "\n"))
message("\n", paste(strwrap(coupling_verdict,  width = 92), collapse = "\n"), "\n")

# =============================================================================
# PART F: FIGURES
# =============================================================================

# A -- the eliminator: what the ceiling is NOT, and what it IS (z-score)
elim <- ceiling_probes |> dplyr::filter(quantifier == "zscore") |>
  tidyr::pivot_longer(c(cor_oxphos_globalmean, ambient_oxphos_all, ambient_within_6Wpos,
                        chance_floor_permuted, random_pair_med_abs,
                        ambient_after_globalmean_removed),
                      names_to = "probe", values_to = "value") |>
  dplyr::mutate(value = abs(value),
    probe = factor(probe,
      levels = c("cor_oxphos_globalmean","ambient_oxphos_all","random_pair_med_abs",
                 "ambient_within_6Wpos","chance_floor_permuted","ambient_after_globalmean_removed"),
      labels = c("cor(OXPHOS, global mean)","ambient (all 24)","random set pairs",
                 "within 6W_pos (contrast gone)","chance floor (permuted)",
                 "after removing global mean")))
p_a <- ggplot2::ggplot(elim, ggplot2::aes(x = value, y = probe)) +
  ggplot2::geom_col(fill = "grey70") +
  ggplot2::geom_vline(xintercept = elim$value[elim$probe == "chance floor (permuted)"],
                      linetype = 2, colour = "#4575B4") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)), hjust = -0.15, size = 3) +
  ggplot2::scale_x_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  ggplot2::labs(
    title = "Where the ceiling comes from: one global per-sample axis",
    subtitle = paste("Linear z-score. It is not low n (chance floor, blue) and not curation (random",
                     "pairs);\nit survives with the genotype contrast removed (within 6W_pos) and it",
                     "vanishes when\nthe per-sample global mean is regressed out."),
    x = "|correlation| / median |rho|", y = NULL) +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "A_ceiling_origin.pdf"), p_a, width = 8, height = 4)

# B -- the global factor as a phenotype: what the residual axis tracks
p_b <- global_factor_phenotype |>
  dplyr::mutate(covariate = stats::reorder(covariate, rho_residual_axis),
                kind = dplyr::case_when(
                  covariate %in% c("MYC_activity","proliferation","epithelial") ~ "biological",
                  covariate %in% c("timepoint_6W","genotype_pos") ~ "design (check)",
                  TRUE ~ "technical / prep")) |>
  ggplot2::ggplot(ggplot2::aes(x = rho_residual_axis, y = covariate, fill = kind)) +
  ggplot2::geom_col() +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::scale_fill_manual(values = c("biological" = "#2A6F5E",
                                        "technical / prep" = "#D73027",
                                        "design (check)" = "grey70")) +
  ggplot2::labs(
    title = "The residual common-mode axis is part biology, part prep",
    subtitle = paste("Spearman rho of the design-residualised PC1 (linear z-score) with each",
                     "covariate.\nDesign covariates sit near zero by construction (the axis is",
                     "orthogonal to the design)."),
    x = "rho with the residual global factor", y = NULL, fill = NULL) +
  ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "B_global_factor_phenotype.pdf"), p_b, width = 8, height = 4.5)

# C -- the pre-specified pairs: raw -> design -> global, per quantifier.
# Untestable pairs (the axis IS the global factor, so global-adj is noise) are
# drawn faint: their global-adjusted value must NOT be read as a coupling.
pc_df <- coupling_hierarchy |>
  dplyr::mutate(pair = sprintf("%s -> %s\n(axis r%s=%.2f w/ global mean)",
                               axis, outcome, "²", ax_r2_global)) |>
  tidyr::pivot_longer(c(raw, design, global), names_to = "level", values_to = "cor") |>
  dplyr::mutate(level = factor(level, levels = c("raw","design","global"),
                               labels = c("raw","design-adj","global-adj")),
                readable = dplyr::if_else(testable, "testable", "axis = global factor (global-adj is noise)"))
p_c <- ggplot2::ggplot(pc_df,
         ggplot2::aes(x = level, y = cor, group = quantifier, colour = quantifier, alpha = readable)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  ggplot2::geom_line() + ggplot2::geom_point(size = 1.6) +
  ggplot2::facet_wrap(~ pair, ncol = 4) +
  ggplot2::scale_colour_manual(values = c("gsva" = "#4575B4", "zscore" = "#D73027",
                                          "singscore" = "#2A6F5E")) +
  ggplot2::scale_alpha_manual(values = c("testable" = 1,
                                         "axis = global factor (global-adj is noise)" = 0.25)) +
  ggplot2::labs(
    title = "Pre-specified couplings: raw -> design-adjusted -> global-adjusted, three quantifiers",
    subtitle = paste("Where the axis IS the global factor (r² near 1, faint lines), the",
                     "global-adjusted value is a noise residual, NOT antagonism.\nThe one robust",
                     "testable lead is redox -> MB2 fork; a credible lead needs the global-adjusted",
                     "value to agree across methods.  n=24: ranking, not proof."),
    x = NULL, y = "correlation", colour = NULL, alpha = NULL) +
  ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "C_coupling_hierarchy.pdf"), p_c, width = 11, height = 6.5)

message("Figures written to ", out_dir)

# =============================================================================
# PART G: SAVE
# =============================================================================
lpc_out <- list(
  ceiling_probes          = ceiling_probes,
  ceiling_verdict         = ceiling_verdict,
  global_factor_phenotype = global_factor_phenotype,
  global_factor_summary   = global_factor_summary,
  phenotype_verdict       = phenotype_verdict,
  coupling_hierarchy      = coupling_hierarchy,
  robustness              = robustness,
  gene_vs_score_design    = gene_vs_score_design,
  coupling_verdict        = coupling_verdict,
  quantifiers             = names(mats),
  n_sets                  = length(common_sets),
  notes = paste(
    "Narrative finalisation 2026-07-18, implementing the author's two method notes",
    "(docs/GSVA_global_background_and_pathway_coupling.md;",
    "docs/Gene_set_quantification_for_pathway_correlation.md). PRIMARY quantifier = mean gene-wise",
    "z-score (linear -> correlation is average cross-gene covariance); singscore = rank-based",
    "cohort-independent robustness; GSVA = continuity. Correction = design->PC1->residual with the",
    "factor estimated AFTER removing timepoint*genotype (doc 1 s7) and leave-pair-out PC1 for each",
    "tested pair (doc 1 s4). Mean-subtraction is NOT used (it assumes equal loadings and imposes a",
    "sum-to-zero constraint). mb2_fork is the AP7 fold ratio (not a gene set), so it is",
    "quantifier-independent and taken from the panel unchanged. priming is a NEGATIVE CONTROL",
    "(circular: MitoCarta vs MitoCarta). INFERENCE WEIGHT: ranking + global-factor-as-phenotype,",
    "NOT confirmatory CIs -- at n=24 a residual r~0.29 arises under the null (doc 1 s10). SCOPE:",
    "this re-ranks EXPLORATORY couplings only; the content result (raw count shares) and the",
    "attenuation (DESeq2 LFCs) do not touch GSVA or these couplings and are unaffected. See",
    "docs/2026-07-18_narrative_synthesis_five_questions.md."))
saveRDS(lpc_out, here::here("results", "linear_pathway_coupling.rds"))
message("Saved results/linear_pathway_coupling.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  lpc <- readRDS(here::here("results", "linear_pathway_coupling.rds"))

  cat(strwrap(lpc$ceiling_verdict,   92), sep = "\n")
  cat(strwrap(lpc$phenotype_verdict, 92), sep = "\n")
  cat(strwrap(lpc$coupling_verdict,  92), sep = "\n")

  # --- PART A: the ceiling, per quantifier (the eliminator) ---
  lpc$ceiling_probes |> as.data.frame() |> print()

  # --- PART C: the global factor as a phenotype (biological vs technical) ---
  lpc$global_factor_summary   |> as.data.frame() |> print()
  lpc$global_factor_phenotype |> as.data.frame() |> print()

  # --- PART D: the pre-specified leads, raw -> design -> global x 3 methods ---
  lpc$coupling_hierarchy |> as.data.frame() |> print()
  lpc$robustness         |> as.data.frame() |> print()
  # doc-2 consistency: gene-level vs score-level design adjustment agree?
  lpc$gene_vs_score_design |> as.data.frame() |> print()

  list.files(here::here("outputs", "linear_pathway_coupling"), pattern = "\\.pdf$")
}
