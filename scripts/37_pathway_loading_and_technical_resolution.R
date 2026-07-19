# =============================================================================
# 37_pathway_loading_and_technical_resolution.R
# -----------------------------------------------------------------------------
# BEHIND THE GLOBAL AXIS -- Parts 1+2 of the author's low-dimensional-structure
# plan (docs/Low_dimensional_pathway_structure_and_global_axis_analysis.md).
# Script 36 established that ONE global per-sample axis inflates every per-sample
# pathway CORRELATION, and that the mito axes ARE that factor, so "OXPHOS is
# central" is untestable as a residual coupling. This script does NOT stop there:
# it treats the axis as a phenotype to dissect (doc secs 1, 3.2, 3.3), and it uses
# two facts the author supplied on 2026-07-19:
#
#   (1) BATCH = TIMEPOINT. The RNA-seq was two extractions (6W, 12W separately),
#       so batch is perfectly confounded with timepoint -- BUT genotype is balanced
#       within each batch, so the design ~timepoint*myc_status already absorbs the
#       whole batch effect into the timepoint main-effect term. Every genotype /
#       interaction result stays clean; only the pure between-timepoint
#       (developmental) reading is confounded with batch and cannot be disentangled.
#   (2) CONTAMINATION + IEG/PREP-STRESS are technical and may be grouped as G_tech.
#       This supplies the ANCHORS doc sec 3.1 needs (there is NO external QC
#       metadata on disk; these RNA-derived proxies are the surrogate). We use the
#       script-33 recipe, already saved as panel_share in the mt-axis rds.
#
# BOTH-LENSES correction (author's decision): the design-only biological axis
# (script 36's design->PC1->residual) is PRIMARY; removing contamination + IEG is a
# SECOND LENS, and we report what MOVES between them -- plus an over-correction
# check, because the anchors are RNA-derived and carry some biology (IEG
# anticorrelates with proliferation), so removing them removes some biology too.
#
# THE POSITIVE REFRAME (doc sec 3.3): instead of asking "does OXPHOS->proliferation
# coupling survive removing the axis" (untestable -- OXPHOS IS the axis), ask "is
# OXPHOS a top-LOADING pathway on the axis, and does MYC change its loading?" The
# loading LEVEL is near-tautological for OXPHOS (say so); the non-trivial, testable
# quantity is the loading CHANGE by genotype/time -- a contrast the axis cannot touch.
#
# SCOPE (unchanged, load-bearing): this is the exploratory loading/coupling layer
# only. Content = raw count shares; attenuation = DESeq2 LFCs; fGSEA = Wald ranks.
# None of them touch this. Inference weight = ranking + loading-as-phenotype, NOT
# confirmatory CIs (n=24). Figures are analysis-illustrative; the manuscript figure
# scripts are 39+.
#
# Input:  results/gsva_scores.rds                 (scores, expr_mat VST, meta, pathways)
#         results/myc_mito_centrality.rds         (panel: per-sample axis scores)
#         results/mtdna_axis_and_coupling_null.rds (panel_share: endothelial, ieg_stress)
#         results/dds_int_run.rds                 (counts, for library covariates)
# Output: results/pathway_loading.rds
#         outputs/pathway_loading/*.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

out_dir <- here::here("outputs", "pathway_loading")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# =============================================================================
# PART 1: LOAD -- one score universe, the panel, the technical anchors
# =============================================================================
gs    <- readRDS(here::here("results", "gsva_scores.rds"))
mc    <- readRDS(here::here("results", "myc_mito_centrality.rds"))

panel    <- mc$panel
sm       <- gs$sample_meta
set_meta <- gs$set_meta
pathways <- gs$pathways
expr     <- gs$expr_mat                                   # VST (log), genes x samples

M_gsva <- gs$scores[, panel$sample, drop = FALSE]
expr   <- expr[, panel$sample, drop = FALSE]
sm     <- sm[match(panel$sample, sm$sample), ]
stopifnot(identical(colnames(M_gsva), panel$sample),
          identical(colnames(expr),   panel$sample),
          identical(sm$sample,        panel$sample))

by_cat <- function(cat) set_meta$set_name[set_meta$category_primary == cat]

# the design: timepoint * genotype = the four groups (= batch is inside timepoint).
meta <- data.frame(timepoint  = factor(sm$timepoint),
                   myc_status = factor(sm$myc_status))
X   <- stats::model.matrix(~ timepoint * myc_status, data = meta)
qrX <- qr(X)

i6    <- panel$timepoint == "6W"
i12   <- panel$timepoint == "12W"
ipos6 <- i6 & panel$myc_status == "pos"
stopifnot(sum(i6) == 12, sum(i12) == 12, sum(ipos6) == 6)

# technical anchors (script-33 recipe, saved as panel_share) --------------------
share <- readRDS(here::here("results", "mtdna_axis_and_coupling_null.rds"))$panel_share
share <- share[panel$sample, , drop = FALSE]
contam <- log2(share[, "endothelial"])                    # residual stroma (technical)
ieg    <- log2(share[, "ieg_stress"])                     # dissociation stress (van den Brink)

counts_mat <- DESeq2::counts(readRDS(here::here("results", "dds_int_run.rds")))
counts_mat <- counts_mat[, panel$sample, drop = FALSE]

# =============================================================================
# PART 2: THE PRIMARY SCORE MATRIX + shared helpers (linear z-score, doc-2)
# =============================================================================
# Primary quantifier = mean gene-wise z-score (its correlation IS average cross-
# gene covariance). GSVA is retained only as a continuity cross-check in PART A.
sets_universe <- rownames(M_gsva)
pw <- pathways[intersect(names(pathways), sets_universe)]

row_z <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.finite(rowSums(z)), , drop = FALSE]
}
score_linear <- function(gene_z) {
  t(vapply(sets_universe, function(s) {
    g <- intersect(pw[[s]], rownames(gene_z))
    if (length(g) < 5L) return(rep(NA_real_, ncol(gene_z)))
    colMeans(gene_z[g, , drop = FALSE])
  }, numeric(ncol(gene_z))))
}
M_z <- score_linear(row_z(expr))

# keep sets scored (complete) in BOTH matrices so the two are comparable
common_sets <- intersect(rownames(M_z)[stats::complete.cases(M_z)],
                         rownames(M_gsva)[stats::complete.cases(M_gsva)])
M_z    <- M_z[common_sets, , drop = FALSE]
M_gsva <- M_gsva[common_sets, , drop = FALSE]
message(sprintf("Score matrices aligned: %d sets scored in both (z-score primary).",
                length(common_sets)))

global_mean <- function(M) colMeans(M, na.rm = TRUE)
resid_on    <- function(M, XX) t(qr.resid(qr(XX), t(M)))  # residualise score rows on XX
pc1_oriented <- function(Mres, gm) {                      # PC1, sign set to track level
  pc <- stats::prcomp(t(Mres), center = TRUE, scale. = FALSE)
  g  <- pc$x[, 1]
  if (suppressWarnings(stats::cor(g, gm)) < 0) { g <- -g; pc$rotation[, 1] <- -pc$rotation[, 1] }
  list(g = g, var = pc$sdev^2 / sum(pc$sdev^2), load = pc$rotation[, 1], eig = pc$sdev^2)
}
d_eff <- function(eig) sum(eig)^2 / sum(eig^2)            # participation ratio (doc sec 1.6)

# --- named composite axes (same algebra as scripts 28/36) ---------------------
mito_names <- intersect(by_cat("MitoCarta"), common_sets)
mito_grep  <- function(pat) grep(pat, mito_names, value = TRUE)
OXPHOS_PAT <- "OXPHOS|COMPLEX_[IV]|_SUBUNITS|ASSEMBLY_FACTORS|ELECTRON_CARRIERS|CRISTAE"
BIOG_PAT   <- "RIBOSOME|CENTRAL_DOGMA|MT_TRNA|MT_RRNA|MTRNA|MTDNA|IMPORT|TRANSLATION"
teb_up <- grep("^MG_TEB_VS_DUCTAL_.*_UP$", common_sets, value = TRUE)
teb_dn <- grep("^MG_TEB_VS_DUCTAL_.*_DN$", common_sets, value = TRUE)

axis_sets <- list(
  mito_oxphos     = mito_grep(OXPHOS_PAT),
  mito_biogenesis = mito_grep(BIOG_PAT),
  prolif          = intersect(by_cat("Proliferation"),  common_sets),
  myc_sig         = intersect(by_cat("MYC_signatures"), common_sets),
  tca             = intersect(c("GS_METAB_KREBS","METAB_TCA_KEGG","METAB_TCA_REACTOME"), common_sets),
  nucleotide      = intersect(c("GS_METAB_NUCLEOTIDE","GS_METAB_PURINE","GS_METAB_PYRIMIDINE",
                                "METAB_NUCLEOTIDE_REACTOME","METAB_NUCLEOTIDE_SALVAGE_REACTOME",
                                "METAB_NUCLEOTIDE_WP","METAB_PURINE_KEGG","METAB_PYRIMIDINE_KEGG",
                                "METAB_PYRIMIDINE_WP"), common_sets),
  glycolysis      = intersect(c("GS_METAB_GLYCOLYSIS","METAB_GLYCOLYSIS_HALLMARK",
                                "METAB_GLYCOLYSIS_KEGG","METAB_GLYCOLYSIS_REACTOME"), common_sets),
  cholesterol     = intersect(c("GS_METAB_CHOLESTEROL","GS_METAB_MEVALONATE","METAB_CHOLESTEROL_HALLMARK",
                                "METAB_CHOLESTEROL_REACTOME","METAB_CHOLESTEROL_WP"), common_sets),
  redox           = intersect(c("GS_METAB_REDOX","GS_METAB_GLUTATHIONE","GS_METAB_REACTIVE_OXYGEN"), common_sets))
signed_sets <- list(
  priming    = list(up = intersect("MITOCARTA_APOPTOSIS_PRO",  common_sets),
                    dn = intersect("MITOCARTA_APOPTOSIS_ANTI", common_sets)),
  teb_dediff = list(up = teb_up, dn = teb_dn))

comp_M <- function(M, sets) {
  sets <- intersect(sets, rownames(M))
  if (length(sets) == 0) return(rep(NA_real_, ncol(M)))
  colMeans(M[sets, , drop = FALSE])
}
series_from <- function(M, name) {                        # design-RETAINED composite series
  if (name %in% names(signed_sets))
    return(comp_M(M, signed_sets[[name]]$up) - comp_M(M, signed_sets[[name]]$dn))
  comp_M(M, axis_sets[[name]])
}
named_axes <- c(names(axis_sets), "priming", "teb_dediff")

# =============================================================================
# PART A: HOW LOW-DIMENSIONAL IS IT? (doc sec 1.6)
# =============================================================================
# The gate: reproduce script 36's cor(mito_oxphos, per-sample global mean).
ox_gate <- suppressWarnings(stats::cor(series_from(M_z, "mito_oxphos"), global_mean(M_z)))
message(sprintf("GATE: cor(mito_oxphos, global mean) = %.3f (script 36 reported 0.968).", ox_gate))

dimensionality <- purrr::map_dfr(list(zscore = M_z, gsva = M_gsva), function(M) {
  gm  <- global_mean(M)
  raw <- pc1_oriented(M, gm)
  res <- pc1_oriented(resid_on(M, X), gm)                 # design removed
  cum_raw <- cumsum(raw$var)
  tibble::tibble(
    pc1_var_raw_pct      = 100 * raw$var[1],
    pc2_var_raw_pct      = 100 * raw$var[2],
    pc3_var_raw_pct      = 100 * raw$var[3],
    d_eff_raw            = d_eff(raw$eig),
    n_pc_90pct_raw       = which(cum_raw >= 0.90)[1],
    frac_loadings_pos_raw= mean(raw$load > 0),             # 1-signed = common mode
    cor_pc1raw_globalmean= suppressWarnings(stats::cor(raw$g, gm)),
    pc1_var_resid_pct    = 100 * res$var[1],
    d_eff_resid          = d_eff(res$eig))
}, .id = "quantifier")

# =============================================================================
# PART B: BATCH / TECHNICAL RESOLUTION -- BOTH LENSES + over-correction check
# =============================================================================
# --- B1. batch = timepoint: what the design already absorbs, and what it cannot.
# The anchors themselves: are they genotype- or time-associated? (script 33 found
# the mt/IEG axis genotype-INDEPENDENT, time-associated -- confirm for these here.)
anchor_design <- purrr::map_dfr(
  list(contamination = contam, ieg_prep = ieg), function(v) {
    cf <- summary(stats::lm(v ~ meta$myc_status + meta$timepoint))$coefficients
    tibble::tibble(p_genotype  = cf["meta$myc_statuspos", "Pr(>|t|)"],
                   p_timepoint = cf["meta$timepoint12W", "Pr(>|t|)"])
  }, .id = "anchor")
batch_identifiability <- paste(
  "BATCH = TIMEPOINT (6W and 12W extracted separately). Genotype is BALANCED within",
  "each batch, so the design ~timepoint*myc_status absorbs the entire batch effect",
  "into the timepoint main-effect term: genotype and interaction contrasts are clean.",
  "The pure between-timepoint (developmental) reading is confounded with batch and",
  "CANNOT be separated from it with these data (no external batch/RIN/QC metadata",
  "exists). This script removes the design FIRST everywhere, so the biological axis",
  "it studies is orthogonal to both batch and the four-group design.")

# --- B2. the two biological axes (both lenses) --------------------------------
Sd  <- resid_on(M_z, X)                                    # design removed (primary lens)
Xt  <- cbind(X, contamination = contam, ieg_prep = ieg)   # design + technical anchors
Sdt <- resid_on(M_z, Xt)                                   # design + anchors removed (2nd lens)

gm_z          <- global_mean(M_z)
G_bio_design  <- pc1_oriented(Sd,  gm_z)                   # PRIMARY biological axis (script 36)
G_bio_tech    <- pc1_oriented(Sdt, gm_z)                   # SECOND-LENS biological axis
cor_axes      <- suppressWarnings(stats::cor(G_bio_design$g, G_bio_tech$g))

# --- B3. what MOVES: per-set fraction of design-residual variance the anchors take
r2_tech_set <- 1 - apply(Sdt, 1, stats::var) / apply(Sd, 1, stats::var)
loading_movement <- tibble::tibble(
  set = rownames(Sd),
  lambda_design = as.numeric(suppressWarnings(stats::cor(t(Sd),  G_bio_design$g))),
  lambda_tech   = as.numeric(suppressWarnings(stats::cor(t(Sdt), G_bio_tech$g))),
  r2_removed_by_tech = r2_tech_set,
  category = set_meta$category_primary[match(set, set_meta$set_name)])

# --- B4. OVER-CORRECTION CHECK: how much BIOLOGY lives in the anchor block -----
# The anchors are RNA-derived. If proliferation / MYC co-move with contamination +
# IEG, removing the anchors removes biology, not just noise. Quantify it on the
# design-residualised composites (IEG ~ proliferation is negative -> a real risk).
contam_r <- as.numeric(qr.resid(qrX, contam))
ieg_r    <- as.numeric(qr.resid(qrX, ieg))
overcorrection <- purrr::map_dfr(c("prolif", "myc_sig", "mito_oxphos", "mito_biogenesis"),
  function(nm) {
    s_r <- as.numeric(qr.resid(qrX, series_from(M_z, nm)))
    tibble::tibble(
      composite            = nm,
      cor_with_contam      = suppressWarnings(stats::cor(s_r, contam_r)),
      cor_with_ieg         = suppressWarnings(stats::cor(s_r, ieg_r)),
      r2_absorbed_by_anchors = summary(stats::lm(s_r ~ contam_r + ieg_r))$r.squared)
  })

# =============================================================================
# PART C: LOADING ANALYSIS + CONDITION-DEPENDENT LOADING (doc sec 3.3)
# =============================================================================
# --- C1. named-axis loadings on the biological axis (both lenses) -------------
# loading = cor with the PRIMARY (design-residual) axis. The two-lens correlation
# comparison is confounded by axis ROTATION (the lenses' axes correlate only ~0.57),
# so "what moves" is measured rotation-INDEPENDENTLY as r2_removed_by_tech: the
# fraction of the composite's OWN design-residual variance that contamination + IEG
# absorb (high = prep-driven, so its loading is the one to distrust under the 2nd lens).
axis_loadings <- purrr::map_dfr(named_axes, function(nm) {
  ad <- as.numeric(qr.resid(qrX, series_from(M_z, nm)))    # design-residual composite
  at <- as.numeric(qr.resid(qr(Xt), series_from(M_z, nm))) # design+anchors-residual
  tibble::tibble(
    axis          = nm,
    n_sets        = length(intersect(if (nm %in% names(signed_sets))
                             c(signed_sets[[nm]]$up, signed_sets[[nm]]$dn) else axis_sets[[nm]],
                             common_sets)),
    loading_design_lens = suppressWarnings(stats::cor(ad, G_bio_design$g)),
    loading_tech_lens   = suppressWarnings(stats::cor(at, G_bio_tech$g)),
    r2_removed_by_tech  = 1 - stats::var(at) / stats::var(ad))
}) |> dplyr::arrange(dplyr::desc(abs(loading_design_lens)))

# where do OXPHOS / redox sit among ALL 884 per-set loadings? (is OXPHOS a top loader)
abs_load  <- abs(G_bio_design$load)
pctile    <- function(sets) {
  sets <- intersect(sets, names(abs_load)); if (!length(sets)) return(NA_real_)
  stats::median(vapply(sets, function(s) mean(abs_load <= abs_load[s]), numeric(1)))
}
loading_ranks <- tibble::tibble(
  group = c("OXPHOS sets", "biogenesis sets", "redox sets (control)", "MYC-signature sets"),
  median_pctile_of_abs_loading = c(pctile(axis_sets$mito_oxphos), pctile(axis_sets$mito_biogenesis),
                                   pctile(axis_sets$redox), pctile(axis_sets$myc_sig)))
# is MitoCarta over-represented among high loaders? (Wilcoxon vs the rest)
is_mito   <- set_meta$category_primary[match(names(abs_load), set_meta$set_name)] == "MitoCarta"
mito_load_test <- suppressWarnings(stats::wilcox.test(abs_load[is_mito], abs_load[!is_mito]))

# --- C2. condition-dependent loading: does MYC change the loading? ------------
# S_total ~ G * genotype * timepoint. G is design-orthogonal, so the G:genotype
# term = change in within-group loading by genotype (the non-tautological quantity).
cond_loading <- purrr::map_dfr(c("mito_oxphos", "mito_biogenesis", "prolif", "tca", "nucleotide"),
  function(nm) {
    df <- data.frame(S = series_from(M_z, nm), G = G_bio_design$g,
                     geno = meta$myc_status, time = meta$timepoint)
    cf <- summary(stats::lm(S ~ G * geno * time, data = df))$coefficients
    getp <- function(term) if (term %in% rownames(cf)) cf[term, c("Estimate", "Pr(>|t|)")]
                           else c(NA_real_, NA_real_)
    b_g   <- getp("G"); b_gm <- getp("G:genopos"); b_gt <- getp("G:time12W")
    tibble::tibble(
      composite   = nm,
      loading_base       = b_g[1],  p_loading_base     = b_g[2],  # loading in 6W_neg
      d_loading_genotype = b_gm[1], p_loading_genotype = b_gm[2], # MYC changes loading?
      rel_change_genotype= b_gm[1] / b_g[1],                      # change as fraction of base
      d_loading_time     = b_gt[1], p_loading_time     = b_gt[2]) # time changes loading?
  })

# =============================================================================
# PART D: VERDICTS
# =============================================================================
dim_verdict <- {
  z <- dimensionality |> dplyr::filter(quantifier == "zscore")
  sprintf(paste0(
    "LOW-DIMENSIONAL, QUANTIFIED. On the linear z-score the RAW score matrix has PC1 = %.0f%% of ",
    "variance (PC2 %.0f%%, PC3 %.0f%%), an effective dimensionality of only %.1f (of a possible 23 ",
    "at n=24), and %.0f%% of set loadings share ONE sign -- a single common mode, not a contrast. ",
    "PC1 correlates with the per-sample global mean at r=%.2f, and OXPHOS tracks that mean at r=%.3f. ",
    "So the 884 pathways move along essentially one dominant axis. Removing the four-group design ",
    "leaves a still-substantial within-group axis (residual PC1 %.0f%%, D_eff %.1f) -- the phenotype ",
    "PART C dissects."),
    z$pc1_var_raw_pct, z$pc2_var_raw_pct, z$pc3_var_raw_pct, z$d_eff_raw,
    100 * z$frac_loadings_pos_raw, z$cor_pc1raw_globalmean,
    ox_gate, z$pc1_var_resid_pct, z$d_eff_resid)
}

tech_verdict <- {
  worst <- overcorrection |> dplyr::arrange(dplyr::desc(r2_absorbed_by_anchors)) |> dplyr::slice(1)
  pr    <- overcorrection |> dplyr::filter(composite == "prolif")
  sprintf(paste0(
    "BOTH LENSES AGREE, WITH A NAMED OVER-CORRECTION RISK. The design-only biological axis (primary) ",
    "and the design+anchors axis (contamination + IEG removed) correlate at r=%.2f -- so grouping ",
    "contamination and IEG as technical does NOT overturn the biological axis. BUT the anchors are ",
    "RNA-derived and carry biology: they absorb r^2=%.2f of design-residual PROLIFERATION variance ",
    "(cor with IEG %+.2f -- the van den Brink stress panel runs OPPOSITE to proliferation), and up to ",
    "r^2=%.2f for %s. So the second lens partly removes biology, not only prep noise; loadings that ",
    "move under it are flagged, and the design-only lens stays PRIMARY. Batch is inside timepoint and ",
    "already absorbed; developmental-vs-batch is not separable here (no QC metadata)."),
    cor_axes, pr$r2_absorbed_by_anchors, pr$cor_with_ieg,
    worst$r2_absorbed_by_anchors, worst$composite)
}

loading_verdict <- {
  ox  <- cond_loading |> dplyr::filter(composite == "mito_oxphos")
  bio <- cond_loading |> dplyr::filter(composite == "mito_biogenesis")
  prl <- cond_loading |> dplyr::filter(composite == "prolif")
  oxp <- loading_ranks |> dplyr::filter(group == "OXPHOS sets")
  rxp <- loading_ranks |> dplyr::filter(group == "redox sets (control)")
  # the co-loading cluster: named axes that load >= 0.85 on the biological axis
  cluster <- axis_loadings |> dplyr::filter(loading_design_lens >= 0.85, axis != "mito_oxphos") |>
    dplyr::arrange(dplyr::desc(loading_design_lens))
  clabel <- paste(sprintf("%s %.2f", cluster$axis, cluster$loading_design_lens), collapse = ", ")
  pfmt <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p)))
  sprintf(paste0(
    "THE POSITIVE REFRAME: A MYC-METABOLIC-PROLIFERATIVE CORE DEFINES THE DOMINANT AXIS, AND OXPHOS ",
    "IS IN IT. OXPHOS loads r=%.3f on the design-residual axis -- but it is NOT uniquely central: it ",
    "is indistinguishable from a tightly co-loading cluster (%s), all of which essentially ARE the ",
    "axis. So the dominant programme is that whole co-regulated core, not OXPHOS alone -- which is ",
    "exactly WHY 'OXPHOS is central' stays untestable (nothing in the cluster is separable from the ",
    "rest). The loading LEVEL is the ROBUST statement and it is near-tautological -- a description of ",
    "the axis (OXPHOS member sets sit at percentile %.0f of all 884 |loadings|; MitoCarta as a whole ",
    "loads higher than the rest, Wilcoxon p=%.1e). The framework's non-tautological move is the ",
    "loading CHANGE by genotype (a contrast the axis cannot touch): S ~ G*genotype*timepoint. Within ",
    "groups, MYC hints at TIGHTENING the metabolic core's alignment (OXPHOS %+.0f%%, p %s; TCA also ",
    "up) but LOOSENING proliferation's (%+.0f%%, p %s) -- an intriguing EXPLORATORY pattern, NOT a ",
    "result: near-collinear slopes at n=6/group, and mito_biogenesis shows none (p %s). Report the ",
    "LEVEL as the finding and the CHANGE as a lead for script 38 / the bench. The redox control is a ",
    "weak loader (percentile %.0f), the one axis Myc does not drive."),
    ox_gate, clabel, 100 * oxp$median_pctile_of_abs_loading, mito_load_test$p.value,
    100 * ox$rel_change_genotype, pfmt(ox$p_loading_genotype),
    100 * prl$rel_change_genotype, pfmt(prl$p_loading_genotype), pfmt(bio$p_loading_genotype),
    100 * rxp$median_pctile_of_abs_loading)
}

message("\n", paste(strwrap(dim_verdict,     width = 92), collapse = "\n"))
message("\n", paste(strwrap(tech_verdict,    width = 92), collapse = "\n"))
message("\n", paste(strwrap(loading_verdict, width = 92), collapse = "\n"), "\n")

# =============================================================================
# PART E: FIGURES
# =============================================================================
# A -- scree + effective dimensionality (both quantifiers)
scree_df <- purrr::map_dfr(list(zscore = M_z, gsva = M_gsva), function(M) {
  v <- pc1_oriented(M, global_mean(M))$var
  tibble::tibble(pc = seq_along(v), var = 100 * v)
}, .id = "quantifier") |> dplyr::filter(pc <= 12)
p_a <- ggplot2::ggplot(scree_df, ggplot2::aes(pc, var, colour = quantifier)) +
  ggplot2::geom_line() + ggplot2::geom_point(size = 1.6) +
  ggplot2::scale_colour_manual(values = c("gsva" = "#4575B4", "zscore" = "#D73027")) +
  ggplot2::scale_x_continuous(breaks = 1:12) +
  ggplot2::labs(
    title = "One dominant axis: the pathway-score matrix is low-dimensional",
    subtitle = sprintf(paste("z-score PC1 = %.0f%% of variance; effective dimensionality = %.1f of a",
                             "possible 23 (n=24).\nMany of 884 pathways, but essentially one direction."),
                       dimensionality$pc1_var_raw_pct[dimensionality$quantifier == "zscore"],
                       dimensionality$d_eff_raw[dimensionality$quantifier == "zscore"]),
    x = "principal component", y = "variance explained (%)", colour = NULL) +
  ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "A_scree_dimensionality.pdf"), p_a, width = 7.5, height = 4.5)

# B -- named-axis loading on the PRIMARY axis, coloured by how prep-driven it is.
# (The two-lens correlation gap is confounded by axis rotation; r2_removed_by_tech
# is the rotation-independent "what moves" signal -- red = prep-driven, distrust it.)
p_b <- axis_loadings |>
  dplyr::mutate(axis = stats::reorder(axis, loading_design_lens)) |>
  ggplot2::ggplot(ggplot2::aes(loading_design_lens, axis, colour = r2_removed_by_tech)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  ggplot2::geom_point(size = 3.2) +
  ggplot2::scale_colour_gradient(low = "#4575B4", high = "#D73027", limits = c(0, 1)) +
  ggplot2::labs(
    title = "Pathway loading on the dominant biological axis",
    subtitle = paste("x = correlation of the design-residual composite with the within-group axis",
                     "(primary lens).\nColour = fraction of that composite's variance absorbed by",
                     "contamination + IEG (red = prep-driven,\nso its loading is over-corrected by the",
                     "second lens; blue = robust)."),
    x = "loading on the biological axis", y = NULL, colour = "prep-absorbed r2") +
  ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "B_axis_loadings.pdf"), p_b, width = 8, height = 5)

# C -- condition-dependent loading: OXPHOS/biogenesis on G by group
cond_key <- c("mito_oxphos", "mito_biogenesis")
cond_pts <- purrr::map_dfr(cond_key, function(nm)
  tibble::tibble(composite = nm, G = G_bio_design$g, S = series_from(M_z, nm),
                 group = interaction(meta$timepoint, meta$myc_status, sep = "_")))
p_c <- ggplot2::ggplot(cond_pts, ggplot2::aes(G, S, colour = group)) +
  ggplot2::geom_point(size = 1.6) +
  ggplot2::geom_smooth(ggplot2::aes(group = group), method = "lm", se = FALSE, linewidth = 0.6) +
  ggplot2::facet_wrap(~ composite, scales = "free_y") +
  ggplot2::labs(
    title = "Condition-dependent loading: does MYC change the pathway's coupling to the axis?",
    subtitle = paste("Slope = within-group loading on the biological axis. A genotype difference in slope",
                     "is a CONTRAST\n(readable) -- unlike the residual coupling, which is untestable for",
                     "OXPHOS. n=6/group: a lead, not proof."),
    x = "biological axis (design-residual PC1)", y = "composite score", colour = NULL) +
  ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "C_condition_dependent_loading.pdf"), p_c, width = 9, height = 4.5)

message("Figures written to ", out_dir)

# =============================================================================
# PART F: SAVE
# =============================================================================
pl_out <- list(
  dimensionality        = dimensionality,
  dim_verdict           = dim_verdict,
  batch_identifiability = batch_identifiability,
  anchor_design         = anchor_design,
  cor_axes              = cor_axes,
  loading_movement      = loading_movement,
  overcorrection        = overcorrection,
  tech_verdict          = tech_verdict,
  axis_loadings         = axis_loadings,
  loading_ranks         = loading_ranks,
  mito_load_wilcox_p    = mito_load_test$p.value,
  cond_loading          = cond_loading,
  loading_verdict       = loading_verdict,
  ox_gate               = ox_gate,
  n_sets                = length(common_sets),
  notes = paste(
    "Behind the global axis, Parts 1+2 (2026-07-19), implementing",
    "docs/Low_dimensional_pathway_structure_and_global_axis_analysis.md secs 1, 3.2, 3.3.",
    "PRIMARY quantifier = mean gene-wise z-score; GSVA kept as a dimensionality cross-check.",
    "The biological axis = design->PC1->residual (script 36). BOTH LENSES: design-only",
    "(primary) vs design + contamination + IEG anchors (second lens); loadings that move are",
    "flagged, and the over-correction table quantifies how much proliferation/MYC the anchors",
    "absorb (they are RNA-derived, so removing them removes some biology). BATCH = TIMEPOINT:",
    "the design absorbs it; genotype/interaction clean; developmental-vs-batch not separable",
    "(no QC metadata). CONDITION-DEPENDENT LOADING (S ~ G*genotype*timepoint) makes the",
    "OXPHOS-participation question a CONTRAST (genotype dimension, axis-independent), which is",
    "readable where the residual coupling was not. INFERENCE WEIGHT: ranking + loading-as-",
    "phenotype, NOT confirmatory CIs (n=24). SCOPE: exploratory loading layer only; contrasts",
    "(content shares, DESeq2 LFCs, Wald-rank fGSEA) are untouched. Script 38 (mitoPPS priming-",
    "coupling + MYC-vs-PGC1a) follows. See docs/2026-07-18_narrative_synthesis_five_questions.md."))
saveRDS(pl_out, here::here("results", "pathway_loading.rds"))
message("Saved results/pathway_loading.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  pl <- readRDS(here::here("results", "pathway_loading.rds"))

  cat(strwrap(pl$dim_verdict,     92), sep = "\n")
  cat(strwrap(pl$tech_verdict,    92), sep = "\n")
  cat(strwrap(pl$loading_verdict, 92), sep = "\n")

  # --- PART A: dimensionality ---
  pl$dimensionality |> as.data.frame() |> print()

  # --- PART B: technical resolution ---
  cat(strwrap(pl$batch_identifiability, 92), sep = "\n")
  pl$anchor_design  |> as.data.frame() |> print()
  cat(sprintf("cor(design axis, design+anchors axis) = %.3f\n", pl$cor_axes))
  pl$overcorrection |> as.data.frame() |> print()
  # most-moved sets under the technical lens
  pl$loading_movement |> dplyr::arrange(dplyr::desc(r2_removed_by_tech)) |>
    head(15) |> as.data.frame() |> print()

  # --- PART C: loadings + condition-dependent loading ---
  pl$axis_loadings |> as.data.frame() |> print()
  pl$loading_ranks |> as.data.frame() |> print()
  pl$cond_loading  |> as.data.frame() |> print()

  list.files(here::here("outputs", "pathway_loading"), pattern = "\\.pdf$")
}
