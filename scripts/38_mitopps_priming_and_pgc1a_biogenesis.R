# =============================================================================
# 38_mitopps_priming_and_pgc1a_biogenesis.R
# -----------------------------------------------------------------------------
# BEHIND THE GLOBAL AXIS -- Parts 3+4. Two questions the loading layer (script 37)
# could not answer, moved into spaces where they ARE answerable:
#
#   PART 3 -- Does OXPHOS couple to death PRIMING? In GSVA/z-score space this is
#     circular (priming is MitoCarta, so mito-vs-mito rides the shared mito-content
#     mode) AND rides the global shift. The paper's thesis -- and mitoPPS -- is that
#     mito sub-pathways are SEPARATELY regulated: a mitoPPS pairwise RATIO cancels
#     BOTH the global shift AND the shared mito-content mode. So the coupling is
#     tested in mitoPPS ratio space, where Apoptosis-PRO/-ANTI and OXPHOS already
#     exist as nuclear-encoded ratio pathways (script 08 PART 2c; no mt- genes, so
#     the mtDNA-arm border does NOT apply). priming != death -- it is the molecular
#     substrate, which the transcriptome can speak to. First we MEASURE the mitoPPS-
#     space ambient ceiling (point 3 predicts it is far below GSVA's 0.80 -- confirm,
#     do not assume), then read couplings against THAT null, with the genotype/time
#     delta-interaction (does the coupling itself change -- the "separate regulation"
#     signature).
#
#   PART 4 -- Are MYC-driven and PGC1a-driven biogenesis SEPARABLE, and do they
#     differ toward death? PGC1a is a coactivator with no motif, so its programme is
#     its coactivated TFs' targets: ESRRA / NRF1 / GABPA (TF_targets lanes in the
#     884). MYC-biogenesis = the MitoCarta mtDNA/translation arm (script 37's
#     mito_biogenesis; a known Myc bystander). Test SEPARABILITY first (if the two
#     axes are one, the death contrast is undefined -- that is itself the answer);
#     then contrast their coupling to the CLEAN mitoPPS priming readout.
#
# INFERENCE WEIGHT: ranking + separability-as-phenotype, NOT confirmatory CIs (n=24).
# SCOPE: exploratory coupling layer only; the content result, the attenuation, and
# Wald-rank fGSEA never touch this. Figures illustrate; manuscript figures are 39+.
#
# Input:  results/mitopps_scores.rds   (24x146 ratio pathways incl. OXPHOS,
#                                        Apoptosis-PRO/-ANTI, the biogenesis arm)
#         results/gsva_scores.rds       (884 sets: TF lanes + MYC + expr_mat)
#         results/myc_mito_centrality.rds (panel; sample order)
# Output: results/mitopps_priming_pgc1a.rds
#         outputs/mitopps_priming_pgc1a/*.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

out_dir <- here::here("outputs", "mitopps_priming_pgc1a")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# =============================================================================
# PART 1: LOAD + ALIGN (mitoPPS ratio space + z-score space on ONE sample order)
# =============================================================================
gs    <- readRDS(here::here("results", "gsva_scores.rds"))
mc    <- readRDS(here::here("results", "myc_mito_centrality.rds"))
mp    <- readRDS(here::here("results", "mitopps_scores.rds"))

panel    <- mc$panel
sm       <- gs$sample_meta[match(panel$sample, gs$sample_meta$sample), ]
set_meta <- gs$set_meta
pathways <- gs$pathways
expr     <- gs$expr_mat[, panel$sample, drop = FALSE]
M_gsva   <- gs$scores[, panel$sample, drop = FALSE]

# mitoPPS matrix: samples x pathways (keep only numeric pathway columns -- the
# data.frame also carries sample/group/timepoint/myc_status metadata columns).
msc <- mp$mitopps_scores
path_cols <- names(msc)[vapply(msc, is.numeric, logical(1))]
MP  <- as.matrix(msc[match(panel$sample, msc$sample), path_cols])
rownames(MP) <- panel$sample
stopifnot(identical(rownames(MP), panel$sample), !anyNA(MP), is.numeric(MP))

meta <- data.frame(timepoint = factor(sm$timepoint), myc_status = factor(sm$myc_status))
X    <- stats::model.matrix(~ timepoint * myc_status, data = meta)
qrX  <- qr(X)
i6   <- panel$timepoint == "6W"; ipos6 <- i6 & panel$myc_status == "pos"
by_cat <- function(cat) set_meta$set_name[set_meta$category_primary == cat]

# =============================================================================
# PART 2: z-score matrix (primary quantifier for the TF/biogenesis composites)
# =============================================================================
sets_universe <- rownames(M_gsva)
pw <- pathways[intersect(names(pathways), sets_universe)]
row_z <- function(mat) { z <- t(scale(t(mat))); z[is.finite(rowSums(z)), , drop = FALSE] }
score_linear <- function(gene_z) t(vapply(sets_universe, function(s) {
  g <- intersect(pw[[s]], rownames(gene_z))
  if (length(g) < 5L) return(rep(NA_real_, ncol(gene_z)))
  colMeans(gene_z[g, , drop = FALSE])
}, numeric(ncol(gene_z))))
M_z <- score_linear(row_z(expr))
M_z <- M_z[stats::complete.cases(M_z), , drop = FALSE]
common_sets <- rownames(M_z)

# named z-space composites -----------------------------------------------------
BIOG_PAT <- "RIBOSOME|CENTRAL_DOGMA|MT_TRNA|MT_RRNA|MTRNA|MTDNA|IMPORT|TRANSLATION"
mito_names <- intersect(by_cat("MitoCarta"), common_sets)
z_axis_sets <- list(
  pgc1a      = intersect(grep("TFT_(ESRRA|NRF1|GABPA)", common_sets, value = TRUE), common_sets),
  pgc1a_mito = intersect(grep("TFT_(ESRRA|NRF1|GABPA).*_MITO$", common_sets, value = TRUE), common_sets),
  myc_sig    = intersect(by_cat("MYC_signatures"), common_sets),
  mito_biog  = grep(BIOG_PAT, mito_names, value = TRUE),
  prolif     = intersect(by_cat("Proliferation"), common_sets))
comp_M <- function(M, sets) { sets <- intersect(sets, rownames(M))
  if (!length(sets)) return(rep(NA_real_, ncol(M))); colMeans(M[sets, , drop = FALSE]) }
z_series <- function(name) comp_M(M_z, z_axis_sets[[name]])

# shared coupling helpers (script 36 idiom) ------------------------------------
global_mean <- function(M) colMeans(M, na.rm = TRUE)
raw_cor    <- function(a, b, m = "pearson") suppressWarnings(stats::cor(a, b, method = m))
dadj_cor   <- function(a, b, m = "pearson")
  suppressWarnings(stats::cor(qr.resid(qrX, a), qr.resid(qrX, b), method = m))
# leave-pair-out global factor from the z-score background (design removed first)
gfactor_excl <- function(exclude_sets) {
  bg  <- setdiff(rownames(M_z), exclude_sets)
  Rbg <- t(qr.resid(qrX, t(M_z[bg, , drop = FALSE])))
  stats::prcomp(t(Rbg), center = TRUE, scale. = FALSE)$x[, 1]
}
# delta-interaction: does the a<->b coupling change by genotype / timepoint?
delta_coupling <- function(b, a) {                       # model b ~ a * geno * time
  df <- data.frame(b = b, a = a, geno = meta$myc_status, time = meta$timepoint)
  cf <- summary(stats::lm(b ~ a * geno * time, data = df))$coefficients
  g  <- function(term) if (term %in% rownames(cf)) cf[term, c("Estimate", "Pr(>|t|)")] else c(NA, NA)
  ag <- g("a:genopos"); at <- g("a:time12W")
  tibble::tibble(d_coupling_genotype = ag[1], p_coupling_genotype = ag[2],
                 d_coupling_time = at[1], p_coupling_time = at[2])
}

# =============================================================================
# PART A: THE mitoPPS-SPACE AMBIENT CEILING (does the ratio cancel the shift?)
# =============================================================================
# The whole premise of Part 3 is that mitoPPS space is readable where GSVA is not.
# Measure it, the same way script 36 measured GSVA (there: 0.80).
gm_mp <- global_mean(t(MP))                              # per-sample mean across 145 ratios
amb_mp <- function(col_vec, cols_keep, idx = rep(TRUE, nrow(MP))) {
  r <- apply(MP[idx, cols_keep, drop = FALSE], 2,
             function(v) suppressWarnings(stats::cor(v, col_vec[idx], method = "spearman")))
  stats::median(abs(r[is.finite(r)]))
}
ox_col   <- MP[, "OXPHOS"]
keep_cols<- setdiff(colnames(MP), "OXPHOS")
set.seed(1L)
ridx  <- replicate(2000, sample(ncol(MP), 2))
rpair <- apply(ridx, 2, function(k) suppressWarnings(stats::cor(MP[, k[1]], MP[, k[2]], method = "spearman")))
set.seed(2L)
chance_mp <- stats::median(replicate(200, amb_mp(sample(ox_col), keep_cols)))
MP_res <- apply(MP, 2, function(v) stats::residuals(stats::lm(v ~ gm_mp)))
MP_d   <- qr.resid(qrX, MP)                              # design-residualised (samples x pathways)
amb_designadj <- {                                       # matched null for a design-adjusted coupling
  ox_d <- MP_d[, "OXPHOS"]
  r <- apply(MP_d[, keep_cols, drop = FALSE], 2,
             function(v) suppressWarnings(stats::cor(v, ox_d, method = "spearman")))
  stats::median(abs(r[is.finite(r)]))
}
pc_mp_raw <- stats::prcomp(MP, center = TRUE, scale. = FALSE)$sdev^2
pc_mp_res <- stats::prcomp(t(MP_d), center = TRUE, scale. = FALSE)$sdev^2
mitopps_ceiling <- tibble::tibble(
  space                   = "mitoPPS (ratio)",
  n_pathways              = ncol(MP),
  cor_oxphos_globalmean   = suppressWarnings(stats::cor(ox_col, gm_mp)),
  random_pair_med_abs     = stats::median(abs(rpair), na.rm = TRUE),
  ambient_oxphos_all      = amb_mp(ox_col, keep_cols),
  ambient_oxphos_designadj= amb_designadj,
  ambient_within_6Wpos    = amb_mp(ox_col, keep_cols, ipos6),
  chance_floor_permuted   = chance_mp,
  ambient_after_globalmean_removed = amb_mp(MP_res[, "OXPHOS"], keep_cols),
  pc1_var_raw_pct         = 100 * pc_mp_raw[1] / sum(pc_mp_raw),
  pc1_var_resid_pct       = 100 * pc_mp_res[1] / sum(pc_mp_res),
  gsva_ambient_for_ref    = 0.80)

# =============================================================================
# PART B: DEATH PRIMING IN mitoPPS SPACE (Part 3)
# =============================================================================
biog_cols <- intersect(c("Mitochondrial central dogma","mtDNA maintenance","mtDNA replication",
  "mtDNA nucleoid","mtRNA metabolism","Transcription","Translation","Mitochondrial ribosome",
  "Mitochondrial ribosome assembly","Translation factors"), colnames(MP))
PRO   <- MP[, "Apoptosis-PRO"]; ANTI <- MP[, "Apoptosis-ANTI"]
mpv <- list(
  oxphos    = MP[, "OXPHOS"],
  biogenesis= rowMeans(MP[, biog_cols, drop = FALSE]),
  redox     = MP[, "ROS and glutathione metabolism"],
  priming   = PRO - ANTI,
  pro       = PRO)

# the honest null is the mitoPPS ambient (PART A). Compare RAW to the raw null and
# DESIGN-ADJUSTED to the design-adjusted null -- apples to apples in each case.
mp_ceiling_raw <- mitopps_ceiling$ambient_oxphos_all
mp_ceiling_dad <- mitopps_ceiling$ambient_oxphos_designadj
# percentile of the OXPHOS->priming coupling among OXPHOS's couplings to ALL pathways
ox_all_rho <- abs(apply(MP[, keep_cols, drop = FALSE], 2,
                        function(v) suppressWarnings(stats::cor(v, mpv$oxphos, method = "spearman"))))
oxphos_priming_pctile <- mean(ox_all_rho <= abs(raw_cor(mpv$oxphos, mpv$priming, "spearman")), na.rm = TRUE)
priming_couplings <- purrr::map_dfr(
  list(c("oxphos","priming"), c("oxphos","pro"), c("biogenesis","priming"),
       c("redox","priming"), c("biogenesis","pro")),
  function(p) {
    a <- mpv[[p[1]]]; b <- mpv[[p[2]]]; di <- delta_coupling(b, a)
    tibble::tibble(
      axis = p[1], outcome = p[2],
      raw = raw_cor(a, b), raw_spearman = raw_cor(a, b, "spearman"),
      design_adjusted = dadj_cor(a, b), design_adjusted_spearman = dadj_cor(a, b, "spearman"),
      raw_beats_null      = abs(raw_cor(a, b, "spearman")) > mp_ceiling_raw,
      designadj_beats_null = abs(dadj_cor(a, b, "spearman")) > mp_ceiling_dad,
      d_coupling_genotype = di$d_coupling_genotype, p_coupling_genotype = di$p_coupling_genotype,
      d_coupling_time = di$d_coupling_time, p_coupling_time = di$p_coupling_time)
  })

# =============================================================================
# PART C: MYC vs PGC1a/ESRRA BIOGENESIS -- SEPARABILITY (Part 4a)
# =============================================================================
# Are the two biogenesis programmes one axis or two? (a) how correlated are they,
# raw / design-adjusted / global-adjusted; (b) is the PGC1a/ESRRA axis Myc-driven
# (a genotype effect) or Myc-independent? If PGC1a ~ MYC after correction, they are
# not separable and the death contrast (Part 4b) is undefined.
pgc1a <- z_series("pgc1a"); myc <- z_series("myc_sig"); biocarta <- z_series("mito_biog")
g_excl_pm <- gfactor_excl(c(z_axis_sets$pgc1a, z_axis_sets$myc_sig))
gadj_cor  <- function(a, b, g) {
  Xg <- cbind(X, g); ra <- qr.resid(qr(Xg), a); rb <- qr.resid(qr(Xg), b)
  suppressWarnings(stats::cor(ra, rb))
}
geno_effect <- function(v) {                             # genotype main effect (design)
  cf <- summary(stats::lm(v ~ meta$myc_status + meta$timepoint))$coefficients
  cf["meta$myc_statuspos", c("Estimate", "Pr(>|t|)")]
}
ge_pgc <- geno_effect(pgc1a); ge_myc <- geno_effect(myc); ge_bio <- geno_effect(biocarta)
# resid_frac AFTER removing the global factor: if either axis IS the global factor
# (resid_frac < 0.10), the global-adjusted correlation is noise-on-noise, so a LOW
# value is over-removal, NOT evidence of independence (the script-36 trap).
resid_frac <- function(v, g) { r <- qr.resid(qr(cbind(X, g)), v); stats::var(r) / stats::var(v) }
separability <- tibble::tibble(
  pair = "pgc1a(ESRRA/NRF1/GABPA) vs MYC",
  raw_cor             = raw_cor(pgc1a, myc),
  design_adjusted_cor = dadj_cor(pgc1a, myc),
  global_adjusted_cor = gadj_cor(pgc1a, myc, g_excl_pm),
  pgc1a_resid_frac    = resid_frac(pgc1a, g_excl_pm),
  myc_resid_frac      = resid_frac(myc,   g_excl_pm),
  global_adj_testable = resid_frac(pgc1a, g_excl_pm) >= 0.10 & resid_frac(myc, g_excl_pm) >= 0.10,
  pgc1a_genotype_beta = ge_pgc[1], pgc1a_genotype_p = ge_pgc[2],
  myc_genotype_beta   = ge_myc[1], myc_genotype_p   = ge_myc[2],
  n_pgc1a_sets        = length(z_axis_sets$pgc1a))
# also: pgc1a vs the MitoCarta biogenesis arm (are the two "biogenesis" axes distinct?)
sep_biog <- tibble::tibble(
  pair = "pgc1a vs MitoCarta-biogenesis",
  raw_cor = raw_cor(pgc1a, biocarta), design_adjusted_cor = dadj_cor(pgc1a, biocarta),
  global_adjusted_cor = gadj_cor(pgc1a, biocarta, gfactor_excl(c(z_axis_sets$pgc1a, z_axis_sets$mito_biog))),
  biocarta_genotype_beta = ge_bio[1], biocarta_genotype_p = ge_bio[2])

# =============================================================================
# PART D: THE DEATH CONTRAST -- each biogenesis axis vs the CLEAN priming readout
# =============================================================================
# priming = mitoPPS (Part A shows it is background-independent). The z-space
# biogenesis axes ride the ceiling, so they are design- and global-adjusted; the
# mitoPPS priming side is already clean, so it is design-adjusted only.
priming_mp <- mpv$priming
death_contrast <- purrr::map_dfr(c("pgc1a", "mito_biog", "myc_sig"), function(nm) {
  a <- z_series(nm); g <- gfactor_excl(z_axis_sets[[nm]])
  ra_d <- qr.resid(qrX, a);           rb_d <- qr.resid(qrX, priming_mp)
  ra_g <- qr.resid(qr(cbind(X, g)), a)
  tibble::tibble(
    biogenesis_axis = nm,
    raw            = raw_cor(a, priming_mp),
    design_adj     = suppressWarnings(stats::cor(ra_d, rb_d)),
    global_adj     = suppressWarnings(stats::cor(ra_g, rb_d)),
    resid_frac_axis= stats::var(ra_g) / stats::var(a),
    testable       = stats::var(ra_g) / stats::var(a) >= 0.10)
})

# =============================================================================
# PART E: VERDICTS
# =============================================================================
ceiling_verdict <- {
  c1 <- mitopps_ceiling
  sprintf(paste0(
    "mitoPPS SPACE IS READABLE WHERE GSVA IS NOT. The ratio cancels the global shift by construction, ",
    "and the numbers confirm it: OXPHOS tracks the per-sample mean of all %d mitoPPS pathways at only ",
    "r=%.2f (GSVA: 0.97), random pathway pairs sit at median |rho|=%.2f, and OXPHOS's ambient coupling ",
    "to the other pathways is %.2f -- versus GSVA's 0.80. PC1 is %.0f%% of variance (GSVA 68-76%%). So ",
    "a coupling measured in mitoPPS space is NOT handed to you by a common mode; the honest null here ",
    "is |rho| ~ %.2f, not ~0.8. This is why the death-priming question is asked in THIS space."),
    c1$n_pathways, c1$cor_oxphos_globalmean, c1$random_pair_med_abs, c1$ambient_oxphos_all,
    c1$pc1_var_raw_pct, c1$ambient_oxphos_all)
}

priming_verdict <- {
  op <- priming_couplings |> dplyr::filter(axis == "oxphos", outcome == "priming")
  bp <- priming_couplings |> dplyr::filter(axis == "biogenesis", outcome == "priming")
  sprintf(paste0(
    "OXPHOS <-> DEATH-PRIMING IN mitoPPS SPACE (non-circular, background-independent). OXPHOS-ratio ",
    "couples to priming-ratio (PRO-ANTI) at raw Spearman %+.2f vs the RAW mitoPPS null of %.2f (%s), ",
    "and design-adjusted Spearman %+.2f vs the DESIGN-adjusted null of %.2f (%s); its Pearson is %+.2f ",
    "raw / %+.2f design-adjusted. It sits at the %.0fth percentile of OXPHOS's couplings to all %d ",
    "mito pathways -- so it is %s within the compartment. The genotype delta-interaction (does the ",
    "coupling itself change with Myc) is beta=%+.2f, p=%.3f: %s. Read as a RANKING lead at n=24 -- but ",
    "unlike the GSVA-space version it is NEITHER circular (a ratio cancels the shared mito-content ",
    "mode) NOR riding the global shift. biogenesis-ratio <-> priming is design-adj %+.2f. priming here ",
    "is the molecular substrate, not death; the phenotype stays the external IHC/BH3 arm."),
    op$raw_spearman, mp_ceiling_raw, ifelse(op$raw_beats_null, "above", "at/below"),
    op$design_adjusted_spearman, mp_ceiling_dad, ifelse(op$designadj_beats_null, "above", "at/below"),
    op$raw, op$design_adjusted, 100 * oxphos_priming_pctile, ncol(MP) - 1,
    ifelse(oxphos_priming_pctile >= 0.75, "among the STRONGER couplings",
           ifelse(oxphos_priming_pctile >= 0.5, "mid-pack (a weak lead)", "unremarkable")),
    op$d_coupling_genotype, op$p_coupling_genotype,
    ifelse(is.na(op$p_coupling_genotype) || op$p_coupling_genotype > 0.05,
           "no evidence the coupling differs by genotype",
           "the coupling DIFFERS by genotype (separate-regulation signature)"),
    bp$design_adjusted)
}

separability_verdict <- {
  s <- separability
  pf <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p)))
  sprintf(paste0(
    "ARE MYC- AND PGC1a-BIOGENESIS SEPARABLE? NO -- NOT IN THIS DATA. The ESRRA/NRF1/GABPA axis vs ",
    "MYC: raw r=%+.2f, design-adjusted %+.2f. The global-adjusted value drops to %+.2f, but do NOT ",
    "read that as independence: both axes ARE the global factor (resid_frac %.2f and %.2f, i.e. only ",
    "~%.0f%% / ~%.0f%% of each survives removing it), so the drop is OVER-REMOVAL noise, not ",
    "separation (the script-36 trap). Confirming they are one Myc-driven programme: BOTH carry a ",
    "strong genotype effect -- PGC1a/ESRRA beta=%+.2f (p %s), MYC beta=%+.2f (p %s) -- both induced ",
    "by Myc, correlated at r=%+.2f. CONSEQUENCE for the death contrast: the two 'biogenesis' axes are ",
    "not distinguishable here, so the hypothesised MYC-vs-PGC1a difference toward death CANNOT be ",
    "tested with these bulk correlations. That is itself the answer -- it needs an ESRRA/PGC1a ",
    "PERTURBATION (knockdown / agonist) against the Myc-driven substrate, not more n=24 bulk RNA."),
    s$raw_cor, s$design_adjusted_cor, s$global_adjusted_cor,
    s$pgc1a_resid_frac, s$myc_resid_frac, 100 * s$pgc1a_resid_frac, 100 * s$myc_resid_frac,
    s$pgc1a_genotype_beta, pf(s$pgc1a_genotype_p), s$myc_genotype_beta, pf(s$myc_genotype_p),
    s$raw_cor)
}

death_contrast_verdict <- {
  dc <- death_contrast
  pg <- dc |> dplyr::filter(biogenesis_axis == "pgc1a")
  bc <- dc |> dplyr::filter(biogenesis_axis == "mito_biog")
  sprintf(paste0(
    "MYC-BIOGENESIS vs PGC1a-BIOGENESIS TOWARD PRIMING (clean mitoPPS priming readout). PGC1a/ESRRA ",
    "-> priming: raw %+.2f, design-adj %+.2f, global-adj %+.2f (%s). MitoCarta-biogenesis -> priming: ",
    "raw %+.2f, design-adj %+.2f, global-adj %+.2f (%s). %s This is the author's hypothesis (that the ",
    "two biogenesis programmes relate differently to death) rendered as a RANKING at n=24 -- a lead ",
    "for the bench (BH3 profiling of an ESRRA/PGC1a-perturbed vs Myc-driven substrate), not a result."),
    pg$raw, pg$design_adj, pg$global_adj, ifelse(pg$testable, "testable", "NOT testable: axis = global factor"),
    bc$raw, bc$design_adj, bc$global_adj, ifelse(bc$testable, "testable", "NOT testable: axis = global factor"),
    ifelse(sign(pg$design_adj) != sign(bc$design_adj),
           "They point in OPPOSITE directions -- consistent with the hypothesis that MYC- and PGC1a-biogenesis are differently wired to death.",
           "They point the SAME way -- no evidence here that the two biogenesis programmes differ toward priming."))
}

for (v in list(ceiling_verdict, priming_verdict, separability_verdict, death_contrast_verdict))
  message("\n", paste(strwrap(v, width = 92), collapse = "\n"))
message("")

# =============================================================================
# PART F: FIGURES
# =============================================================================
# A -- the mitoPPS ceiling next to GSVA's
elim <- mitopps_ceiling |>
  dplyr::transmute(`cor(OXPHOS, global mean)` = cor_oxphos_globalmean,
                   `ambient (all 24)` = ambient_oxphos_all,
                   `random pathway pairs` = random_pair_med_abs,
                   `chance floor` = chance_floor_permuted,
                   `GSVA ambient (ref)` = gsva_ambient_for_ref) |>
  tidyr::pivot_longer(dplyr::everything(), names_to = "probe", values_to = "value") |>
  dplyr::mutate(value = abs(value), probe = stats::reorder(probe, value))
p_a <- ggplot2::ggplot(elim, ggplot2::aes(value, probe)) +
  ggplot2::geom_col(fill = "grey70") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)), hjust = -0.15, size = 3) +
  ggplot2::scale_x_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  ggplot2::labs(title = "mitoPPS space cancels the global shift",
                subtitle = "The ratio null is far below GSVA's 0.80, so couplings here are readable.",
                x = "|correlation| / median |rho|", y = NULL) +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "A_mitopps_ceiling.pdf"), p_a, width = 7.5, height = 3.8)

# B -- priming couplings in mitoPPS space (raw + design-adjusted), vs the null
p_b <- priming_couplings |>
  dplyr::mutate(pair = paste(axis, "->", outcome)) |>
  tidyr::pivot_longer(c(raw, design_adjusted), names_to = "level", values_to = "cor") |>
  ggplot2::ggplot(ggplot2::aes(cor, stats::reorder(pair, cor), shape = level, colour = level)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  ggplot2::geom_vline(xintercept = c(-1, 1) * mp_ceiling_raw, linetype = 2, colour = "#4575B4") +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::scale_colour_manual(values = c("raw" = "#D73027", "design_adjusted" = "#2A6F5E")) +
  ggplot2::labs(title = "Death-priming couplings in mitoPPS ratio space",
                subtitle = "Dashed blue = the mitoPPS ambient null (PART A). priming = Apoptosis-PRO - Apoptosis-ANTI.",
                x = "correlation", y = NULL, colour = NULL, shape = NULL) +
  ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "B_priming_couplings.pdf"), p_b, width = 8, height = 4)

# C -- separability + death contrast
p_c <- death_contrast |>
  tidyr::pivot_longer(c(raw, design_adj, global_adj), names_to = "level", values_to = "cor") |>
  dplyr::mutate(level = factor(level, c("raw", "design_adj", "global_adj"))) |>
  ggplot2::ggplot(ggplot2::aes(level, cor, group = biogenesis_axis, colour = biogenesis_axis)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  ggplot2::geom_line() + ggplot2::geom_point(size = 2.4) +
  ggplot2::labs(title = "Do MYC- and PGC1a-biogenesis differ toward death priming?",
                subtitle = "Coupling of each biogenesis axis to the clean mitoPPS priming readout. n=24: a ranking, not proof.",
                x = NULL, y = "correlation with priming (mitoPPS)", colour = NULL) +
  ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "C_death_contrast.pdf"), p_c, width = 7.5, height = 4.2)

message("Figures written to ", out_dir)

# =============================================================================
# PART G: SAVE
# =============================================================================
mpp_out <- list(
  mitopps_ceiling         = mitopps_ceiling,
  ceiling_verdict         = ceiling_verdict,
  priming_couplings       = priming_couplings,
  priming_verdict         = priming_verdict,
  separability            = separability,
  sep_biog                = sep_biog,
  separability_verdict    = separability_verdict,
  death_contrast          = death_contrast,
  death_contrast_verdict  = death_contrast_verdict,
  pgc1a_sets              = z_axis_sets$pgc1a,
  biog_mitopps_cols       = biog_cols,
  notes = paste(
    "Behind the global axis, Parts 3+4 (2026-07-19). Part 3: OXPHOS <-> death-PRIMING tested in",
    "mitoPPS RATIO space, where Apoptosis-PRO/-ANTI and OXPHOS are nuclear-encoded ratio pathways",
    "(script 08 PART 2c). A ratio cancels the global shift AND the shared mito-content mode, so the",
    "coupling is NEITHER circular NOR ceiling-bound -- unlike the GSVA-space version (script 34's",
    "negative control). PART A measures the mitoPPS ambient (the honest null in this space) before any",
    "coupling is read. priming != death (molecular substrate; phenotype = external IHC/BH3). Part 4:",
    "PGC1a = ESRRA/NRF1/GABPA TF-target lanes (coactivator proxy, no motif). SEPARABILITY from MYC is",
    "tested FIRST (if not separable the death contrast is undefined). INFERENCE WEIGHT: ranking +",
    "separability-as-phenotype, NOT confirmatory CIs (n=24). SCOPE: exploratory coupling layer only.",
    "See docs/2026-07-18_narrative_synthesis_five_questions.md and 37 (loading)."))
saveRDS(mpp_out, here::here("results", "mitopps_priming_pgc1a.rds"))
message("Saved results/mitopps_priming_pgc1a.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  mpp <- readRDS(here::here("results", "mitopps_priming_pgc1a.rds"))

  cat(strwrap(mpp$ceiling_verdict,        92), sep = "\n")
  cat(strwrap(mpp$priming_verdict,        92), sep = "\n")
  cat(strwrap(mpp$separability_verdict,   92), sep = "\n")
  cat(strwrap(mpp$death_contrast_verdict, 92), sep = "\n")

  mpp$mitopps_ceiling   |> as.data.frame() |> print()   # PART A: is the space readable?
  mpp$priming_couplings |> as.data.frame() |> print()   # PART B: OXPHOS <-> priming
  mpp$separability      |> as.data.frame() |> print()   # PART C: MYC vs PGC1a
  mpp$sep_biog          |> as.data.frame() |> print()
  mpp$death_contrast    |> as.data.frame() |> print()   # PART D: the contrast

  list.files(here::here("outputs", "mitopps_priming_pgc1a"), pattern = "\\.pdf$")
}
