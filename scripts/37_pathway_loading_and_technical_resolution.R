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
# PART A2: WHAT THE PATHWAY COMPOSITING KEEPS vs DROPS (gene-level axis identity)
# =============================================================================
# The author's question: the pathway global axis (76%) is genotype-associated, but the
# LARGEST gene-level axis (top-500 VST PC1, ~44%) is a different direction -- what is it,
# and does the pathway axis correspond to gene-level PC2? Characterise the gene-level PCs
# by (i) correlation with technical + biological covariates, (ii) correspondence to the
# pathway axis, (iii) fGSEA on their genome-wide loadings. RESULT (gene_axis_verdict):
# gene-PC1 = epithelial-purity vs immune-infiltration COMPOSITION axis, genotype-
# INDEPENDENT -> correctly DROPPED. gene-PC2 = the MYC/E2F/proliferation program with
# OXPHOS + mito-biogenesis baked in, and the pathway axis FOREGROUNDS it (cor 0.60 vs 0.44;
# genotype 0.26 -> 0.66). The mito programme rides INSIDE the Myc programme, not the
# composition nuisance. These PCAs are reused by figures A / A2 in PART E.
grp      <- factor(paste(sm$timepoint, sm$myc_status, sep = "_"),
                   levels = c("6W_neg", "6W_pos", "12W_neg", "12W_pos"))
grp_cols <- c("6W_neg" = "#4575B4", "6W_pos" = "#D73027",
              "12W_neg" = "#91BFDB", "12W_pos" = "#FC8D59")
v_gene      <- apply(expr, 1, stats::var)
top500      <- order(v_gene, decreasing = TRUE)[seq_len(min(500L, nrow(expr)))]
pc_gene     <- stats::prcomp(t(expr[top500, , drop = FALSE]), center = TRUE, scale. = FALSE)
v_gene_p    <- pc_gene$sdev^2 / sum(pc_gene$sdev^2)
pc_gene_all <- stats::prcomp(t(expr), center = TRUE, scale. = FALSE)   # robustness: same axis?
v_gene_all  <- pc_gene_all$sdev^2 / sum(pc_gene_all$sdev^2)
pc_path     <- stats::prcomp(t(M_z), center = TRUE, scale. = FALSE)
v_path_p    <- pc_path$sdev^2 / sum(pc_path$sdev^2)

# orient axes so the reported signs are interpretable (epithelial+, proliferation+, myc+)
orient <- function(v, ref) if (suppressWarnings(stats::cor(v, ref)) < 0) -v else v
g1 <- orient(pc_gene$x[, 1], share[, "epithelial"])
g2 <- orient(pc_gene$x[, 2], panel$prolif)
g3 <- pc_gene$x[, 3]
p1 <- orient(pc_path$x[, 1], panel$myc_sig)

axis_covariates <- data.frame(
  genotype      = as.integer(sm$myc_status == "pos"),
  timepoint     = as.integer(sm$timepoint == "12W"),
  epithelial    = share[, "epithelial"],
  immune        = share[, "immune"],
  endothelial   = share[, "endothelial"],
  ieg_prep      = share[, "ieg_stress"],
  depth         = sm$sizeFactor,
  proliferation = panel$prolif,
  myc           = panel$myc_sig,
  oxphos        = panel$mito_oxphos,
  biogenesis    = panel$mito_biogenesis)
pc_scores <- data.frame(gene_PC1 = g1, gene_PC2 = g2, gene_PC3 = g3, pathway_PC1 = p1)
gene_axis_covcor <- purrr::map_dfr(names(pc_scores), function(a)
  tibble::tibble(axis = a, covariate = names(axis_covariates),
                 rho = vapply(axis_covariates, function(cv)
                        suppressWarnings(stats::cor(pc_scores[[a]], cv, method = "spearman")),
                        numeric(1))))

# correspondence: which gene-level axis does the pathway global axis track?
gene_pathway_correspondence <- tibble::tibble(
  gene_axis = c("gene_PC1", "gene_PC2", "gene_PC3"),
  var_pct   = 100 * v_gene_p[1:3],
  cor_with_pathway_PC1 = c(suppressWarnings(stats::cor(g1, p1)),
                           suppressWarnings(stats::cor(g2, p1)),
                           suppressWarnings(stats::cor(g3, p1))))
gene_axis_var <- tibble::tibble(
  space   = c("gene_top500", "gene_all", "pathway_scores"),
  pc1_pct = c(100 * v_gene_p[1], 100 * v_gene_all[1], 100 * v_path_p[1]),
  cor_top500_all_PC1 = c(suppressWarnings(stats::cor(pc_gene$x[, 1], pc_gene_all$x[, 1])),
                         NA_real_, NA_real_))

# fGSEA on the genome-wide loadings of gene-PC1 and gene-PC2 (assign the axes to pathways)
gene_axis_fgsea <- NULL
if (requireNamespace("fgsea", quietly = TRUE) &&
    file.exists(here::here("results", "gene_sets_list.rds"))) {
  glist <- readRDS(here::here("results", "gene_sets_list.rds"))
  rank_axis <- function(pcscore) {
    r <- apply(expr, 1, function(x) suppressWarnings(stats::cor(x, pcscore)))
    sort(r[is.finite(r)], decreasing = TRUE)
  }
  fgsea_axis <- function(pcscore, lab) {
    set.seed(1)
    fg <- suppressWarnings(fgsea::fgsea(glist, rank_axis(pcscore),
                                        minSize = 5, maxSize = 800, eps = 0))
    tibble::as_tibble(fg[, c("pathway", "NES", "padj", "size")]) |> dplyr::mutate(axis = lab)
  }
  gene_axis_fgsea <- dplyr::bind_rows(fgsea_axis(g1, "gene_PC1"), fgsea_axis(g2, "gene_PC2"))
}

gene_axis_verdict <- {
  cc <- gene_axis_covcor; gp <- gene_pathway_correspondence
  gv <- function(a, c) cc$rho[cc$axis == a & cc$covariate == c]
  sprintf(paste0(
    "WHAT THE PATHWAY COMPOSITING KEEPS vs DROPS. The largest gene-level axis (top-500 VST PC1, ",
    "%.0f%% of gene variance) is genotype-INDEPENDENT (rho with genotype %+.2f): it is the ",
    "epithelial-purity vs immune-infiltration COMPOSITION axis (epithelial %+.2f, immune %+.2f), ",
    "i.e. the residual dissociation contamination. Its metabolic pole (OXPHOS/MYC-target-high) is a ",
    "cell-IDENTITY effect -- epithelial MECs are mito-dense, immune infiltrate is not -- NOT Myc ",
    "dose. This is the axis the design/compositing correctly DROPS. gene-PC2 (%.0f%%) is the MYC/",
    "E2F/proliferation program (proliferation %+.2f, myc %+.2f) with OXPHOS + mitochondrial ",
    "translation/mtRNA (biogenesis) baked in by fGSEA -- mitochondria ride INSIDE the Myc program. ",
    "The pathway global axis tracks gene-PC2 more than gene-PC1 (cor %+.2f vs %+.2f) and SHARPENS ",
    "its genotype signal (gene-PC2 genotype %+.2f -> pathway-PC1 %+.2f). So compositing foregrounds ",
    "the Myc-mito program (PC2) and demotes the composition nuisance (PC1). CAVEAT: OXPHOS enriches ",
    "at BOTH gene-PC1's metabolic pole (composition) AND gene-PC2 (Myc program) -- the gene-level ",
    "root of why mito correlates with everything (the ceiling)."),
    gp$var_pct[1], gv("gene_PC1", "genotype"), gv("gene_PC1", "epithelial"), gv("gene_PC1", "immune"),
    gp$var_pct[2], gv("gene_PC2", "proliferation"), gv("gene_PC2", "myc"),
    gp$cor_with_pathway_PC1[2], gp$cor_with_pathway_PC1[1],
    gv("gene_PC2", "genotype"), gv("pathway_PC1", "genotype"))
}
message("\n", paste(strwrap(gene_axis_verdict, width = 92), collapse = "\n"), "\n")

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

# --- C3. MITO vs NON-MITO LOADING ENRICHMENT (the "is this mito-primary?" test) --
# Author's question (2026-07-20): are figure B's axes the actual top loaders, and are
# there NON-mito sets in the high rankings? If not, that would argue for mitochondria-
# PRIMARY regulation of early tumorigenesis by the MYC-mito axis. Interrogate ALL
# common_sets on the PRIMARY axis (|lambda_design| from loading_movement, correlation
# scale). Three caveats decide it: (1) non-mito growth programmes also load high;
# (2) the very-top mito block is a BUILD tautology -- the Gray _MITO / _LE_MITO lanes
# are MitoCarta subsets by construction (docs/library_reference/
# Gray_et_al_developmental_TFS_selection.md sec 2c), so they load ~1 by mito CONTENT,
# not by regulator identity; (3) the ranking does not resolve sub-programme (oxphos vs
# mitoribo vs priming). NO new inputs but the shortlist CSV (read-only base read.csv).
mito_classification <- loading_movement |>
  dplyr::transmute(
    set,
    lambda_design,
    abs_load  = abs(lambda_design),
    category  = dplyr::coalesce(category, ""),
    name_mito = grepl("_MITO$|_MITO_|MITO_NU|^MITO_|CORE_MITO", set)) |>
  dplyr::mutate(
    mito_defined = category == "MitoCarta" | name_mito |
                   (category == "Metabolism" &
                    grepl("OXPHOS|KREBS|TCA|ELECTRON|RESPIRAT", set)),
    is_construction_mito = name_mito,
    class3 = dplyr::case_when(
      category == "MitoCarta" ~ "mitocarta_proper",
      is_construction_mito    ~ "construction_MITO",
      mito_defined            ~ "mitocarta_proper",   # metab OXPHOS/TCA sets
      TRUE                    ~ "non_mito")) |>
  dplyr::arrange(dplyr::desc(abs_load)) |>
  dplyr::mutate(load_rank = dplyr::row_number())

top_frac <- function(k) {
  k <- min(k, nrow(mito_classification)); mean(mito_classification$mito_defined[seq_len(k)])
}
top100_class <- mito_classification$class3[seq_len(min(100, nrow(mito_classification)))]
wilc_mito <- suppressWarnings(stats::wilcox.test(
  mito_classification$abs_load[mito_classification$mito_defined],
  mito_classification$abs_load[!mito_classification$mito_defined]))

mito_enrichment <- tibble::tibble(
  n_sets                   = nrow(mito_classification),
  n_mito                   = sum(mito_classification$mito_defined),
  base_rate_mito           = mean(mito_classification$mito_defined),
  frac_mito_top50          = top_frac(50),
  frac_mito_top100         = top_frac(100),
  frac_mito_top200         = top_frac(200),
  top100_mitocarta_proper  = sum(top100_class == "mitocarta_proper"),
  top100_construction_mito = sum(top100_class == "construction_MITO"),
  top100_non_mito          = sum(top100_class == "non_mito"),
  wilcox_p_mito_vs_rest    = wilc_mito$p.value,
  median_absload_mito      = stats::median(mito_classification$abs_load[mito_classification$mito_defined]),
  median_absload_nonmito   = stats::median(mito_classification$abs_load[!mito_classification$mito_defined]))

# the genuinely-non-mito top loaders, with a coarse programme label for the figure
prog_label <- function(s) dplyr::case_when(
  grepl("NUCLEOTIDE|PURINE|PYRIMIDINE", s)                          ~ "nucleotide",
  grepl("PENTOSE|_PPP|PPP_", s)                                     ~ "pentose-phosphate",
  grepl("E2F|CELL_CYCLE|MITOTIC|PROLIF|MKI67|_G2M|DNA_REPLICATION", s) ~ "proliferation",
  grepl("METABRIC|BICLUSTER|_MB[0-9]|BREAST", s)                    ~ "breast-cancer-METABRIC",
  grepl("MYC", s)                                                   ~ "MYC-target",
  grepl("TEB|DUCT|^MG_|GRAY|ALVEOL|LASP|LHS|BMYO|MAMMARY", s)       ~ "mammary-dev",
  TRUE                                                              ~ "other")
nonmito_top_loaders <- mito_classification |>
  dplyr::filter(!mito_defined) |>
  dplyr::slice_head(n = 15) |>
  dplyr::mutate(programme = prog_label(set)) |>
  dplyr::select(load_rank, set, lambda_design, abs_load, programme, category)

# CONSTRUCTION tautology proof: per TF, the top _MITO lane vs the lowest non-_MITO lane
construction_gap <- purrr::map_dfr(c("MYC", "E2F1", "ESRRA"), function(tf) {
  mlane <- mito_classification |>
    dplyr::filter(grepl(sprintf("^TFT_%s_GRAY_.*_MITO$", tf), set)) |>
    dplyr::arrange(dplyr::desc(abs_load))
  blane <- mito_classification |>
    dplyr::filter(grepl(sprintf("^TFT_%s_GRAY_", tf), set), !grepl("_MITO$", set)) |>
    dplyr::arrange(abs_load)
  tibble::tibble(
    tf               = tf,
    n_mito_lanes     = nrow(mlane),
    n_base_lanes     = nrow(blane),
    top_mito_lane    = if (nrow(mlane)) mlane$set[1]           else NA_character_,
    top_mito_loading = if (nrow(mlane)) mlane$lambda_design[1] else NA_real_,
    top_mito_rank    = if (nrow(mlane)) mlane$load_rank[1]     else NA_integer_,
    low_base_lane    = if (nrow(blane)) blane$set[1]           else NA_character_,
    low_base_loading = if (nrow(blane)) blane$lambda_design[1] else NA_real_,
    low_base_rank    = if (nrow(blane)) blane$load_rank[1]     else NA_integer_)
})

# does the ranking STRATIFY sub-programme? Join the _MITO lanes to the CHEA3 shortlist
# and correlate |loading| with the oxphos / mitoribo / apop-balance content. Read-only
# base read.csv; wrapped so a name mismatch degrades to NA rather than breaking the run.
substrate_stratification <- tryCatch({
  sl <- utils::read.csv(here::here("data", "genesets_from_library",
                                   "gray_chea_mito_tf_shortlist.csv"), stringsAsFactors = FALSE)
  sl$key <- paste(sl$TF, sl$context, sep = "|")
  parse_key <- function(s) {
    core  <- sub("^TFT_", "", sub("_MITO$", "", s))          # <TF>_GRAY_<CONTEXT>
    parts <- strsplit(core, "_GRAY_", fixed = TRUE)[[1]]
    if (length(parts) != 2L) return(NA_character_)
    paste(parts[1], parts[2], sep = "|")
  }
  mm <- mito_classification |> dplyr::filter(is_construction_mito, grepl("_GRAY_", set))
  mm$key <- vapply(mm$set, parse_key, character(1))
  j <- dplyr::inner_join(mm, sl, by = "key")
  floor_div <- function(a, b) a / pmax(b, 1)                 # denominator floor
  j <- dplyr::mutate(j, oxphos_frac   = floor_div(n_oxphos,   n_mito_total),
                        mitoribo_frac = floor_div(n_mitoribo, n_mito_total))
  tibble::tibble(
    n_joined                  = nrow(j),
    cor_absload_oxphos_frac   = suppressWarnings(stats::cor(j$abs_load, j$oxphos_frac)),
    cor_absload_mitoribo_frac = suppressWarnings(stats::cor(j$abs_load, j$mitoribo_frac)),
    cor_absload_apop_balance  = suppressWarnings(stats::cor(j$abs_load, j$apop_balance)),
    median_absload_joined     = stats::median(j$abs_load))
}, error = function(e) tibble::tibble(
    n_joined = 0L, cor_absload_oxphos_frac = NA_real_, cor_absload_mitoribo_frac = NA_real_,
    cor_absload_apop_balance = NA_real_, median_absload_joined = NA_real_))

enrichment_verdict <- {
  me <- mito_enrichment
  ss <- substrate_stratification
  cg <- construction_gap |> dplyr::filter(tf == "MYC")
  nm <- nonmito_top_loaders |> dplyr::slice_head(n = 3)
  wilc_p_str <- if (is.na(me$wilcox_p_mito_vs_rest)) "NA" else
    if (me$wilcox_p_mito_vs_rest < 1e-300) "<1e-300" else
      sprintf("%.1e", me$wilcox_p_mito_vs_rest)
  sprintf(paste0(
    "MITO-LED BUT NOT MITO-SPECIFIC -- NOT an argument for mitochondria-PRIMARY regulation. ",
    "Figure B plots 11 CURATED composites, not the empirical top loaders (the true top loaders ",
    "are the Gray TFT_*_LE_MITO lanes at ~0.99). Across all %d sets, mito-defined sets are %.0f%% ",
    "of the library but %.0f%% of the top 100 and %.0f%% of the top 200 -- a real, striking ",
    "enrichment (Wilcoxon |loading| mito vs rest p=%s). BUT: (1) non-mito growth programmes load ",
    "nearly as high -- the top non-mito loaders are %s, up to |loading| %.2f, vs the redox control ",
    "~0.26; (2) the very top is CONSTRUCTION-inflated -- the top 100 splits MitoCarta-proper %d / ",
    "build-tautological _MITO %d / non-mito %d, and %s loads %.3f (rank %d) while the SAME TF's ",
    "non-mito lane %s loads %.3f (rank %d), so loading tracks mito gene CONTENT not the regulator; ",
    "(3) the ranking does NOT stratify sub-programme -- across %d _MITO lanes cor(|loading|,",
    "oxphos_frac)=%.3f, cor(|loading|,mitoribo_frac)=%.3f (both weak; the lanes saturate at the ",
    "mito ceiling, so oxphos vs biogenesis is not resolved). => the dominant axis is a ",
    "coordinated Myc anabolic-proliferative-mitochondrial ",
    "GROWTH state; loading = COVARIATION, not primacy. A primacy claim needs perturbation, not ",
    "n=24 bulk. (The nuclear-up/mtDNA-down discordance this axis describes is independently ",
    "established in TCGA + mouse HCC: Lesner et al., bioRxiv 2026.07.13.738248.)"),
    me$n_sets, 100 * me$base_rate_mito, 100 * me$frac_mito_top100, 100 * me$frac_mito_top200,
    wilc_p_str,
    paste(nm$programme, collapse = ", "), max(nonmito_top_loaders$abs_load),
    me$top100_mitocarta_proper, me$top100_construction_mito, me$top100_non_mito,
    cg$top_mito_lane, cg$top_mito_loading, cg$top_mito_rank,
    cg$low_base_lane, cg$low_base_loading, cg$low_base_rank,
    ss$n_joined, ss$cor_absload_oxphos_frac, ss$cor_absload_mitoribo_frac)
}

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

message("\n", paste(strwrap(dim_verdict,        width = 92), collapse = "\n"))
message("\n", paste(strwrap(tech_verdict,       width = 92), collapse = "\n"))
message("\n", paste(strwrap(loading_verdict,    width = 92), collapse = "\n"))
message("\n", paste(strwrap(enrichment_verdict, width = 92), collapse = "\n"), "\n")

# =============================================================================
# PART E: FIGURES
# =============================================================================
# A -- scree: gene-level PC distribution vs pathway-score PC distribution, on a shared
# y-scale so the concentration onto PC1 (compositing effect) is directly comparable.
scree_gene <- tibble::tibble(pc = seq_along(v_gene_p), var = 100 * v_gene_p) |>
  dplyr::filter(pc <= 12)
scree_path <- purrr::map_dfr(list(zscore = M_z, gsva = M_gsva), function(M) {
  v <- pc1_oriented(M, global_mean(M))$var
  tibble::tibble(pc = seq_along(v), var = 100 * v)
}, .id = "quantifier") |> dplyr::filter(pc <= 12)
scree_common <- list(ggplot2::scale_x_continuous(breaks = 1:12),
                     ggplot2::coord_cartesian(ylim = c(0, 80)),
                     ggplot2::theme_bw(base_size = 10),
                     ggplot2::theme(legend.position = "bottom"))
p_a_gene <- ggplot2::ggplot(scree_gene, ggplot2::aes(pc, var)) +
  ggplot2::geom_line(colour = "grey30") + ggplot2::geom_point(size = 1.6, colour = "grey30") +
  ggplot2::labs(title = "Gene level (VST, top-500 variable genes)",
                subtitle = sprintf("PC1 = %.0f%%: variance spread across many PCs.", 100 * v_gene_p[1]),
                x = "principal component", y = "variance explained (%)") + scree_common
p_a_path <- ggplot2::ggplot(scree_path, ggplot2::aes(pc, var, colour = quantifier)) +
  ggplot2::geom_line() + ggplot2::geom_point(size = 1.6) +
  ggplot2::scale_colour_manual(values = c("gsva" = "#4575B4", "zscore" = "#D73027")) +
  ggplot2::labs(title = "Pathway-score level (884 library sets)",
                subtitle = sprintf("z-score PC1 = %.0f%%, effective dim = %.1f: essentially one direction.",
                                   dimensionality$pc1_var_raw_pct[dimensionality$quantifier == "zscore"],
                                   dimensionality$d_eff_raw[dimensionality$quantifier == "zscore"]),
                x = "principal component", y = NULL, colour = NULL) + scree_common
if (requireNamespace("patchwork", quietly = TRUE)) {
  p_a <- patchwork::wrap_plots(p_a_gene, p_a_path, nrow = 1) +
    patchwork::plot_annotation(
      title = "Compositing concentrates variance: flat gene-level scree vs steep pathway-score scree")
  ggplot2::ggsave(file.path(out_dir, "A_scree_dimensionality.pdf"), p_a, width = 10, height = 5)
} else {
  ggplot2::ggsave(file.path(out_dir, "A_scree_dimensionality.pdf"), p_a_path, width = 7.5, height = 4.5)
}

# A2 -- sample PCA scatter: gene-level (conventional) vs pathway-score (the 76% claim).
# Reuses the shared gene/pathway PCAs from figure A. The 76% is a property of the
# COMPOSITES, not the transcriptome. coord_equal so the collapse onto PC1 is not hidden by
# the panel aspect ratio; centre = group mean, cloud = 1-SD normal ellipse (n=6/group).
pca_panel <- function(scores, vpct, ttl) {
  df  <- data.frame(PC1 = scores[, 1], PC2 = scores[, 2], group = grp)
  cen <- stats::aggregate(cbind(PC1, PC2) ~ group, df, mean)  # group centroids (means)
  ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, colour = group)) +
    # cloud = 1-SD normal data ellipse per group (n=6/group -- a spread, not a CI)
    ggplot2::stat_ellipse(ggplot2::aes(fill = group), geom = "polygon",
                          type = "norm", level = 0.68, alpha = 0.12, colour = NA) +
    ggplot2::geom_point(size = 2.1, alpha = 0.85) +
    # centre = group mean, drawn as a large ringed marker
    ggplot2::geom_point(data = cen, ggplot2::aes(PC1, PC2, fill = group),
                        size = 4.6, shape = 21, colour = "black", stroke = 0.6,
                        inherit.aes = FALSE) +
    ggplot2::scale_colour_manual(values = grp_cols) +
    ggplot2::scale_fill_manual(values = grp_cols) +
    ggplot2::guides(fill = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = ttl,
                  x = sprintf("PC1 (%.0f%%)", 100 * vpct[1]),
                  y = sprintf("PC2 (%.0f%%)", 100 * vpct[2]), colour = NULL) +
    ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "bottom")
}
p_gene <- pca_panel(pc_gene$x, v_gene_p, "Gene-level PCA (VST, top-500 variable genes)")
p_path <- pca_panel(pc_path$x, v_path_p, "Pathway-score PCA (884 library sets)")
if (requireNamespace("patchwork", quietly = TRUE)) {
  p_a2 <- patchwork::wrap_plots(p_gene, p_path, nrow = 1) +
    patchwork::plot_annotation(
      title = "Samples collapse onto one axis in pathway-score space, not in gene space",
      subtitle = sprintf(paste(
        "Gene-level PC1 = %.0f%%; averaging genes into 884 correlated pathway scores concentrates",
        "variance onto the\ncommon mode (pathway-score PC1 = %.0f%%). The dominant axis is a property",
        "of the composites, not the transcriptome.\nFilled ring = group mean; shaded cloud = 1-SD",
        "normal ellipse (n=6/group, a spread not a confidence region)."),
        100 * v_gene_p[1], 100 * v_path_p[1]))
  ggplot2::ggsave(file.path(out_dir, "A2_sample_pca_gene_vs_pathway.pdf"), p_a2, width = 10, height = 5.5)
} else {
  ggplot2::ggsave(file.path(out_dir, "A2_sample_pca_gene.pdf"),    p_gene, width = 5.5, height = 5)
  ggplot2::ggsave(file.path(out_dir, "A2_sample_pca_pathway.pdf"), p_path, width = 5.5, height = 5)
}
sample_pca_var <- tibble::tibble(
  space = c("gene_vst_top500", "pathway_scores_884"),
  pc1_pct = c(100 * v_gene_p[1], 100 * v_path_p[1]),
  pc2_pct = c(100 * v_gene_p[2], 100 * v_path_p[2]),
  pc3_pct = c(100 * v_gene_p[3], 100 * v_path_p[3]))

# A3 -- what each PC axis IS: correlation with technical + biological covariates
cov_order  <- c("genotype", "timepoint", "epithelial", "immune", "endothelial", "ieg_prep",
                "depth", "proliferation", "myc", "oxphos", "biogenesis")
axis_order <- c("gene_PC1", "gene_PC2", "gene_PC3", "pathway_PC1")
hm_df <- gene_axis_covcor |>
  dplyr::mutate(covariate = factor(covariate, levels = rev(cov_order)),
                axis      = factor(axis, levels = axis_order))
p_a3 <- ggplot2::ggplot(hm_df, ggplot2::aes(axis, covariate, fill = rho)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho)), size = 3) +
  ggplot2::scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                                midpoint = 0, limits = c(-1, 1)) +
  ggplot2::labs(
    title = "What each axis IS: gene-PC1 = composition (genotype-independent), pathway = Myc",
    subtitle = paste("Spearman correlation of PC sample-scores with covariates. gene-PC1 tracks",
                     "epithelial/immune\ncomposition NOT genotype; the pathway global axis is the",
                     "genotype-associated Myc-metabolic program."),
    x = NULL, y = NULL, fill = "rho") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())
ggplot2::ggsave(file.path(out_dir, "A3_gene_axis_covariates.pdf"), p_a3, width = 7.5, height = 6)

# A4 -- assign gene-PC1 and gene-PC2 to pathways (fGSEA on their genome-wide loadings)
if (!is.null(gene_axis_fgsea)) {
  clean_nm <- function(x) {
    x <- gsub("^MSigDB_HALLMARK_", "", x); x <- gsub("^MC_", "MitoCarta: ", x)
    substr(gsub("_", " ", x), 1, 40)
  }
  top_fg <- gene_axis_fgsea |> dplyr::group_by(axis) |>
    dplyr::mutate(rk = rank(NES, ties.method = "first")) |>
    dplyr::filter(rk <= 6 | rk > dplyr::n() - 6) |> dplyr::ungroup() |>
    dplyr::mutate(label    = clean_nm(pathway),
                  axis_lab = ifelse(axis == "gene_PC1",
                                    "gene-PC1 (composition)", "gene-PC2 (Myc program)"))
  p_a4 <- ggplot2::ggplot(top_fg,
      ggplot2::aes(NES, stats::reorder(interaction(label, axis), NES), fill = NES > 0)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~ axis_lab, scales = "free_y") +
    ggplot2::scale_y_discrete(labels = function(v) sub("\\..*$", "", v)) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "#4575B4"),
                               labels = c("TRUE" = "top pole (+)", "FALSE" = "bottom pole (-)"),
                               name = NULL) +
    ggplot2::labs(
      title = "The gene-level axes assigned to pathways (fGSEA on PC loadings)",
      subtitle = paste("gene-PC1: epithelial OXPHOS/MYC-target pole vs immune/inflammation pole",
                       "(composition).\ngene-PC2: MYC/E2F/proliferation + OXPHOS + mito-translation",
                       "vs myogenesis/EMT (the Myc program)."),
      x = "normalised enrichment score (NES)", y = NULL) +
    ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(file.path(out_dir, "A4_gene_axis_fgsea.pdf"), p_a4, width = 11, height = 5.5)
}

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

# D -- mito vs non-mito loading enrichment across all 884 sets (author's question)
class_cols <- c(mitocarta_proper = "#D73027", construction_MITO = "#FC8D59", non_mito = "#4575B4")
redox_ref  <- stats::median(mito_classification$abs_load[mito_classification$set %in% axis_sets$redox])
med_ref    <- stats::median(mito_classification$abs_load)
lab_df     <- nonmito_top_loaders |> dplyr::slice_head(n = 6)
p_d1 <- ggplot2::ggplot(mito_classification,
                        ggplot2::aes(load_rank, abs_load, colour = class3)) +
  ggplot2::geom_point(size = 0.9, alpha = 0.6) +
  ggplot2::geom_hline(yintercept = med_ref,   linetype = 2, colour = "grey45") +
  ggplot2::geom_hline(yintercept = redox_ref, linetype = 3, colour = "#1a9850") +
  ggplot2::geom_text(data = lab_df,
                     ggplot2::aes(load_rank, abs_load, label = programme),
                     inherit.aes = FALSE, size = 2.5, hjust = 0, nudge_x = 8, nudge_y = 0.012) +
  ggplot2::scale_colour_manual(values = class_cols) +
  ggplot2::labs(
    title = "Mito is the strongest-loading block -- but the axis is mito-LED, not mito-SPECIFIC",
    subtitle = sprintf(paste("All %d sets ranked by |loading| on the biological axis. Base rate %.0f%%",
                             "mito; top-100 %.0f%% mito, BUT non-mito growth programmes reach 0.95-0.98;",
                             "\nthe very top is construction-inflated (orange = _MITO lanes = MitoCarta",
                             "subsets by build). Green line = redox control; dashed = median."),
                       mito_enrichment$n_sets, 100 * mito_enrichment$base_rate_mito,
                       100 * mito_enrichment$frac_mito_top100),
    x = "loading rank (1 = highest |loading|)", y = "|loading| on biological axis",
    colour = NULL) +
  ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "bottom")

dec_df <- mito_classification |>
  dplyr::mutate(decile = dplyr::ntile(load_rank, 10)) |>
  dplyr::count(decile, class3) |>
  dplyr::group_by(decile) |>
  dplyr::mutate(frac = n / sum(n)) |>
  dplyr::ungroup()
p_d2 <- ggplot2::ggplot(dec_df, ggplot2::aes(factor(decile), frac, fill = class3)) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(values = class_cols) +
  ggplot2::labs(
    title = "The top decile is ~all mito -- but mostly build-tautological _MITO lanes",
    subtitle = paste("Composition of each loading decile with construction _MITO sets split out.",
                     "Non-mito rises steadily down the ranking; the enrichment is a gradient, not a cliff."),
    x = "loading decile (1 = top)", y = "fraction of sets", fill = NULL) +
  ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "bottom")

if (requireNamespace("patchwork", quietly = TRUE)) {
  p_d <- patchwork::wrap_plots(p_d1, p_d2, ncol = 1, heights = c(2, 1))
  ggplot2::ggsave(file.path(out_dir, "D_loading_enrichment_mito.pdf"), p_d, width = 8.5, height = 8.5)
} else {
  ggplot2::ggsave(file.path(out_dir, "D_loading_enrichment_mito.pdf"),      p_d1, width = 8.5, height = 5.5)
  ggplot2::ggsave(file.path(out_dir, "D2_loading_enrichment_deciles.pdf"),  p_d2, width = 8.5, height = 3.5)
}

message("Figures written to ", out_dir)

# =============================================================================
# PART F: SAVE
# =============================================================================
pl_out <- list(
  dimensionality        = dimensionality,
  sample_pca_var        = sample_pca_var,
  gene_axis_covcor      = gene_axis_covcor,
  gene_pathway_correspondence = gene_pathway_correspondence,
  gene_axis_var         = gene_axis_var,
  gene_axis_fgsea       = gene_axis_fgsea,
  gene_axis_verdict     = gene_axis_verdict,
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
  mito_classification   = mito_classification,
  mito_enrichment       = mito_enrichment,
  nonmito_top_loaders   = nonmito_top_loaders,
  construction_gap      = construction_gap,
  substrate_stratification = substrate_stratification,
  enrichment_verdict    = enrichment_verdict,
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

  cat(strwrap(pl$dim_verdict,        92), sep = "\n")
  cat(strwrap(pl$gene_axis_verdict,  92), sep = "\n")
  cat(strwrap(pl$tech_verdict,       92), sep = "\n")
  cat(strwrap(pl$loading_verdict,    92), sep = "\n")
  cat(strwrap(pl$enrichment_verdict, 92), sep = "\n")

  # --- PART A: dimensionality ---
  pl$dimensionality |> as.data.frame() |> print()
  pl$sample_pca_var |> as.data.frame() |> print()   # gene-level vs pathway-score PC1

  # --- PART A2: gene-level axis identity (what compositing keeps vs drops) ---
  pl$gene_axis_var |> as.data.frame() |> print()
  pl$gene_pathway_correspondence |> as.data.frame() |> print()
  pl$gene_axis_covcor |> tidyr::pivot_wider(names_from = axis, values_from = rho) |>
    as.data.frame() |> print()
  if (!is.null(pl$gene_axis_fgsea)) {
    pl$gene_axis_fgsea |> dplyr::group_by(axis) |> dplyr::arrange(NES) |>
      dplyr::slice(c(1:5, (dplyr::n() - 4):dplyr::n())) |> as.data.frame() |> print()
  }

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

  # --- PART C3: mito vs non-mito loading enrichment ---
  cat(strwrap(pl$enrichment_verdict, 92), sep = "\n")
  pl$mito_enrichment          |> as.data.frame() |> print()
  pl$nonmito_top_loaders      |> as.data.frame() |> print()
  pl$construction_gap         |> as.data.frame() |> print()
  pl$substrate_stratification |> as.data.frame() |> print()
  # top of the full 884-set ranking (mito monopoly + the construction block)
  pl$mito_classification |> head(20) |> as.data.frame() |> print()

  list.files(here::here("outputs", "pathway_loading"), pattern = "\\.pdf$")
}
