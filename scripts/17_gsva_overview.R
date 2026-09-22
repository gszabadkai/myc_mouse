# scripts/17_gsva_overview.R
# =============================================================================
# GSVA overview + trajectory visualisation (Block A, Day 3 pre-step)
# =============================================================================
#
# Motivation:
#   Script 15 produced results/gsva_scores.rds -- 902 GSVA-tagged library sets
#   scored per-sample (VST, cohort-relative). Before the Day-3 action-point
#   scripts consume that matrix, we read it: per gene-set category, WHICH
#   programs move across the four conditions, and -- for the trajectory-tied
#   categories -- in WHICH direction each genotype moves.
#
#   This is the H1-vs-H3 discriminator. Gate 1 found the Myc effect ~halves from
#   6W to 12W (a negative interaction). Two mechanisms give the same attenuation:
#     H1 front-loaded footprint -- Myc+ program high at 6W, falls by 12W; WT flat.
#     H3 WT developmental catch-up -- Myc+ holds; WT (Myc-) rises to meet it.
#   The interaction term equals d12 - d6 (the attenuation) so we RANK on it, but
#   it is mechanism-blind. The two-line trajectory profile shows which genotype
#   line moved, which is what names the mechanism.
#
# Statistical basis -- lm(score ~ timepoint * myc_status), treatment coding,
#   reference 6W / neg. Every coefficient is a difference of group means:
#     (Intercept)              mean(6W_neg)                    WT baseline at 6W
#     timepoint (12W-6W)       mean(12W_neg) - mean(6W_neg)    Myc- (WT) trajectory
#     myc_status (pos-neg)     mean(6W_pos)  - mean(6W_neg)    d6 (Myc gap at 6W)
#     timepoint:myc_status     d12 - d6                        the Gate-1 attenuation
#   Derived: d12 = myc_status + interaction; Myc+ slope = timepoint + interaction.
#   We fit all 902 sets at once by ordinary least squares on the shared design
#   matrix (identical to per-set lm(); no moderation).
#
# Statistical caveat + validation:
#   The per-set lm interaction is underpowered (n=6/group; interactions are the
#   lowest-power contrast) and nothing survives BH within category. The library
#   sets also overlap heavily (set algebra), so per-set "consistency" is inflated
#   by non-independence. PART 6 therefore adds two correlation-aware / independent
#   checks on the SAME data before any of this is trusted for figure text:
#     (a) CAMERA (competitive, inter-gene-correlation aware) + ROAST (self-
#         contained rotation) on the gene-level interaction contrast, per set and
#         per category (union meta-set = one powered test per category); and
#     (b) a cross-check correlating each set's GSVA beta_int (d12 - d6) against the
#         mean DESeq2 interaction Wald stat of its member genes -- convergence of
#         the GSVA trajectory with the independent count-level interaction.
#
# Input:
#   - results/gsva_scores.rds : list(scores [set x 24], set_meta, sample_meta,
#                                     expr_mat [symbol x 24 VST], pathways)
#   - results/interaction_results.rds : $interaction_raw (DESeq2 interaction Wald)
#   - results/ortholog_table.rds : ENSMUSG -> mouse symbol (for the stat mapping)
#
# Output:
#   - results/gsva_overview.rds : list(coef_table [+ camera/roast/mean_int_stat],
#                                     cond_means, cat_camera, crosscheck, ...)
#   - outputs/gsva_overview/<category>_landscape.pdf        (all categories)
#   - outputs/gsva_overview/<category>_profiles.pdf         (trajectory cats)
#   - outputs/gsva_overview/<category>_dumbbell.pdf         (trajectory cats)
#   - outputs/gsva_overview/crosscheck_beta_int_vs_gene_interaction.pdf
#   - outputs/gsva_overview/perm_contrast_biogenesis_vs_mammary.pdf
#   - outputs/gsva_overview/mechanism_biogenesis_vs_development.pdf
#
# PART 6b (between-category permutation) is the powered test: are program CLASSES
#   on DIFFERENT trajectories (biogenesis/Myc decline vs lineage rise)? Gene-level,
#   pre-defined categories, sample-label rotation null (correlation-aware).
# PART 6c (per-sample coupling) asks whether developmental re-differentiation
#   ABSORBS the Myc biogenesis advantage (Myc drive is constant, so this is about
#   the substrate, not inhibition of Myc activity).
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# =============================================================================
# PART 1: LOAD GSVA SCORES + ALIGN SAMPLE METADATA
# =============================================================================

gsva_out    <- readRDS(here::here("results", "gsva_scores.rds"))
scores      <- gsva_out$scores                       # set x sample (24)
set_meta    <- gsva_out$set_meta                     # set_name / category_primary / ...
sample_meta <- as.data.frame(gsva_out$sample_meta)
expr_mat    <- gsva_out$expr_mat                     # symbol x 24 VST (full universe)
pathways    <- gsva_out$pathways                     # mouse-symbol gene-set lists

# Gene-level validation inputs (PART 6). Fail early with a clear message if an
# older gsva_scores.rds (without expr_mat/pathways) is present -- re-run 15.
if (is.null(expr_mat) || is.null(pathways)) {
  stop("gsva_scores.rds lacks expr_mat/pathways; re-run scripts/15_gsva_scoring.R ",
       "(it now stores them for the gene-level tests in script 17 PART 6).")
}
interaction_raw <- readRDS(here::here("results", "interaction_results.rds"))$interaction_raw
ortholog        <- readRDS(here::here("results", "ortholog_table.rds"))

stopifnot(all(c("group", "myc_status", "timepoint") %in% colnames(sample_meta)))

# Align sample_meta rows to the score matrix column order.
sample_meta <- sample_meta[colnames(scores), , drop = FALSE]
stopifnot(identical(rownames(sample_meta), colnames(scores)))

# Fix factor coding explicitly (do not depend on how colData was stored).
sample_meta$timepoint  <- stats::relevel(as.factor(sample_meta$timepoint),  "6W")
sample_meta$myc_status <- stats::relevel(as.factor(sample_meta$myc_status), "neg")
sample_meta$group      <- factor(sample_meta$group,
                                 levels = c("6W_neg", "6W_pos", "12W_neg", "12W_pos"))
stopifnot(!anyNA(sample_meta$group))

message(sprintf("Loaded %d sets x %d samples; groups: %s",
                nrow(scores), ncol(scores),
                paste(levels(sample_meta$group), collapse = ", ")))

# =============================================================================
# PART 2: CONDITION-MEAN MATRIX + TRAJECTORY POINT ESTIMATES
# =============================================================================
# Every point estimate we report is a cell-mean difference, so we compute them
# straight from the group means -- contrast-coding independent, hence robust to
# whatever global options(contrasts=...) the R session happens to have set (a
# non-treatment default would make lm coefficients NOT equal the simple gaps).

group_order <- c("6W_neg", "6W_pos", "12W_neg", "12W_pos")

# Single source of truth for the 4-cell partition: derive it from the model
# factors (timepoint x myc_status) so it cannot drift from the design used for
# the p-value below. Diagnose (do not fail) if a stored 'group' column disagrees.
cell <- factor(paste(sample_meta$timepoint, sample_meta$myc_status, sep = "_"),
               levels = group_order)
stopifnot(!anyNA(cell))
if (!identical(as.character(cell), as.character(sample_meta$group))) {
  message("NOTE: stored 'group' differs from timepoint x myc_status; ",
          "using the factor-derived partition.")
}

cond_means <- vapply(group_order, function(g) {
  cols <- rownames(sample_meta)[cell == g]
  rowMeans(scores[, cols, drop = FALSE])
}, numeric(nrow(scores)))
rownames(cond_means) <- rownames(scores)             # set names
colnames(cond_means) <- group_order

d6        <- cond_means[, "6W_pos"]  - cond_means[, "6W_neg"]   # Myc gap at 6W
d12       <- cond_means[, "12W_pos"] - cond_means[, "12W_neg"]  # Myc gap at 12W
beta_time <- cond_means[, "12W_neg"] - cond_means[, "6W_neg"]   # Myc- (WT) slope
myc_slope <- cond_means[, "12W_pos"] - cond_means[, "6W_pos"]   # Myc+ slope
beta_int  <- d12 - d6                                          # attenuation

# =============================================================================
# PART 3: INTERACTION P-VALUE (vectorised OLS; treatment contrasts forced)
# =============================================================================
# The lm is used ONLY for the interaction p-value -- the significance of the gap
# change (d12 - d6). Force treatment contrasts so the OLS interaction estimate
# equals d12 - d6 exactly (the p-value itself is contrast-invariant for this
# 1-df test). df = n - p = 24 - 4 = 20.

X <- stats::model.matrix(
  ~ timepoint * myc_status, data = sample_meta,
  contrasts.arg = list(timepoint  = "contr.treatment",
                       myc_status = "contr.treatment")
)
cn      <- colnames(X)
int_col <- grep(":", cn)
stopifnot(length(int_col) == 1)

Y        <- t(scores)                                # samples (24) x sets
XtX_inv  <- solve(crossprod(X))
beta     <- XtX_inv %*% crossprod(X, Y)              # coeff x sets
resid    <- Y - X %*% beta
df_resid <- nrow(X) - ncol(X)
sigma2   <- colSums(resid^2) / df_resid
se_int   <- sqrt(XtX_inv[int_col, int_col] * sigma2)
t_int    <- beta[int_col, ] / se_int
p_int    <- 2 * stats::pt(-abs(t_int), df = df_resid)

# Sanity: under treatment contrasts the OLS interaction estimate == d12 - d6.
stopifnot(max(abs(beta[int_col, rownames(scores)] - beta_int)) < 1e-6)
message("Coding check passed: OLS interaction == d12 - d6 (means).")

coef_table <- tibble::tibble(
  set_name  = rownames(scores),
  beta_time = beta_time,                             # WT (Myc-) trajectory
  d6        = d6,                                    # Myc gap at 6W
  d12       = d12,                                   # Myc gap at 12W
  beta_int  = beta_int,                              # d12 - d6 (attenuation)
  myc_slope = myc_slope,                             # Myc+ trajectory
  int_p     = p_int
) |>
  dplyr::left_join(
    set_meta |> dplyr::select(set_name, category_primary),
    by = "set_name"
  ) |>
  # BH-adjust the interaction p WITHIN each category (per-category emphasis)
  dplyr::group_by(category_primary) |>
  dplyr::mutate(int_padj_within_cat = stats::p.adjust(int_p, method = "BH")) |>
  dplyr::ungroup()

# =============================================================================
# PART 4: LANDSCAPE HEATMAP, ONE PER CATEGORY
# =============================================================================
# Rows = the top ~35 sets by interaction p (the movers); row-z-scored across the
# 4 condition means (GSVA scores are small / cohort-relative, so without scaling
# the map shows magnitude not movement). Row annotations carry the WT slope and
# the interaction so the mechanism is legible per row.

out_dir <- here::here("outputs", "gsva_overview")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

top_n_landscape <- 35L
categories <- sort(unique(coef_table$category_primary))

safe_name <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

for (cat in categories) {
  sub <- coef_table |>
    dplyr::filter(category_primary == cat) |>
    dplyr::arrange(int_p) |>
    head(top_n_landscape)
  if (nrow(sub) < 2) next

  cm <- cond_means[sub$set_name, , drop = FALSE]
  # row z-score across the 4 means; drop zero-variance rows defensively
  rsd <- apply(cm, 1, stats::sd)
  keep <- rsd > 0
  cm <- cm[keep, , drop = FALSE]
  sub <- sub[keep, , drop = FALSE]
  if (nrow(cm) < 2) next
  cm_z <- t(scale(t(cm)))

  amax_t <- max(abs(sub$beta_time)); amax_t <- if (amax_t == 0) 1 else amax_t
  amax_i <- max(abs(sub$beta_int));  amax_i <- if (amax_i == 0) 1 else amax_i
  row_ha <- ComplexHeatmap::rowAnnotation(
    `WT slope`    = sub$beta_time,
    `interaction` = sub$beta_int,
    col = list(
      `WT slope`    = circlize::colorRamp2(c(-amax_t, 0, amax_t),
                                           c("#4575B4", "white", "#D73027")),
      `interaction` = circlize::colorRamp2(c(-amax_i, 0, amax_i),
                                           c("#4575B4", "white", "#D73027"))
    ),
    annotation_name_gp = grid::gpar(fontsize = 7)
  )

  n_row <- nrow(cm_z)
  fs_row <- max(5, min(9, 11 - 0.08 * n_row))
  ht <- ComplexHeatmap::Heatmap(
    cm_z,
    name             = "row z\n(cond mean)",
    cluster_rows     = TRUE,
    clustering_distance_rows = "pearson",
    clustering_method_rows   = "ward.D2",
    cluster_columns  = FALSE,
    column_order     = group_order,
    show_row_names   = TRUE,
    row_names_gp     = grid::gpar(fontsize = fs_row),
    column_names_gp  = grid::gpar(fontsize = 8),
    right_annotation = row_ha,
    column_title     = sprintf("%s -- top %d sets by interaction p", cat, n_row),
    column_title_gp  = grid::gpar(fontsize = 9)
  )

  pdf_h <- max(4, min(20, 1.5 + n_row * 0.22))
  grDevices::pdf(file.path(out_dir, paste0(safe_name(cat), "_landscape.pdf")),
                 width = 7.5, height = pdf_h)
  ComplexHeatmap::draw(ht, merge_legend = TRUE)
  grDevices::dev.off()
}
message(sprintf("Wrote %d landscape heatmaps to %s", length(categories), out_dir))

# =============================================================================
# PART 5: TRAJECTORY PANELS (trajectory categories only)
# =============================================================================
# Panel A: 4-point two-line profile -- reads H1 (Myc+ descends) vs H3 (Myc- ascends).
# Panel B: d6 -> d12 dumbbell -- the gap-collapse summary across many sets.

trajectory_categories <- c("Mammary_development", "MYC_signatures",
                           "Biogenesis_discrimination", "Apoptosis",
                           "Biogenesis_apoptosis_intersections")
present <- intersect(trajectory_categories, categories)
missing <- setdiff(trajectory_categories, categories)
if (length(missing) > 0) {
  warning("Trajectory categories not found in set_meta: ",
          paste(missing, collapse = ", "))
}

top_n_profile  <- 15L
top_n_dumbbell <- 30L

# long condition-mean table for profiles/dumbbell construction
cond_long_all <- tibble::as_tibble(cond_means, rownames = "set_name") |>
  tidyr::pivot_longer(dplyr::all_of(group_order),
                      names_to = "group", values_to = "mean_score") |>
  tidyr::separate(group, into = c("timepoint", "myc_status"),
                  sep = "_", remove = FALSE) |>
  dplyr::mutate(
    timepoint  = factor(timepoint, levels = c("6W", "12W")),
    myc_status = factor(myc_status, levels = c("neg", "pos"))
  )

myc_cols <- c(neg = "#4575B4", pos = "#D73027")

for (cat in present) {
  ord <- coef_table |>
    dplyr::filter(category_primary == cat) |>
    dplyr::arrange(int_p)

  # --- Panel A: two-line profiles (top movers) ---
  prof_sets <- head(ord$set_name, top_n_profile)
  prof_df <- cond_long_all |>
    dplyr::filter(set_name %in% prof_sets) |>
    dplyr::mutate(set_name = factor(set_name, levels = prof_sets))

  p_prof <- ggplot2::ggplot(
    prof_df,
    ggplot2::aes(x = timepoint, y = mean_score,
                 colour = myc_status, group = myc_status)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_colour_manual(values = myc_cols, name = "Myc") +
    ggplot2::facet_wrap(~ set_name, scales = "free_y", ncol = 3) +
    ggplot2::labs(
      title = sprintf("%s -- genotype trajectories (top %d by interaction)",
                      cat, length(prof_sets)),
      subtitle = "Myc+ line descending = H1 (front-loaded); Myc- line ascending = H3 (WT catch-up)",
      x = NULL, y = "mean GSVA score") +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 6))

  n_prof_row <- ceiling(length(prof_sets) / 3)
  ggplot2::ggsave(file.path(out_dir, paste0(safe_name(cat), "_profiles.pdf")),
                  p_prof, width = 8, height = max(3, 1.2 + n_prof_row * 1.5))

  # --- Panel B: d6 -> d12 dumbbell (gap collapse) ---
  dumb <- head(ord, top_n_dumbbell) |>
    dplyr::mutate(set_name = factor(set_name, levels = rev(set_name)))

  p_dumb <- ggplot2::ggplot(dumb) +
    ggplot2::geom_segment(
      ggplot2::aes(x = d6, xend = d12, y = set_name, yend = set_name),
      colour = "grey60", linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(x = d6,  y = set_name, colour = "d6 (6W gap)"),  size = 2) +
    ggplot2::geom_point(ggplot2::aes(x = d12, y = set_name, colour = "d12 (12W gap)"), size = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    ggplot2::scale_colour_manual(
      values = c("d6 (6W gap)" = "#D73027", "d12 (12W gap)" = "#4575B4"),
      name = NULL) +
    ggplot2::labs(
      title = sprintf("%s -- Myc gap collapse (top %d by interaction)",
                      cat, nrow(dumb)),
      subtitle = "d6 = 6W_pos - 6W_neg; d12 = 12W_pos - 12W_neg",
      x = "genotype difference (GSVA score)", y = NULL) +
    ggplot2::theme_bw(base_size = 8)

  ggplot2::ggsave(file.path(out_dir, paste0(safe_name(cat), "_dumbbell.pdf")),
                  p_dumb, width = 7, height = max(3, 1 + nrow(dumb) * 0.22))
}
message(sprintf("Wrote trajectory panels for: %s",
                paste(present, collapse = ", ")))

# =============================================================================
# PART 6: GENE-LEVEL VALIDATION (CAMERA / ROAST + cross-check)
# =============================================================================
# The per-set GSVA lm interaction is underpowered and its cross-set "consistency"
# is inflated by set overlap. These gene-level tests are correlation-aware /
# independent, on the SAME VST matrix + design used for the GSVA scores.

# Align the VST matrix columns to the model rows (both are dds order, but assert).
stopifnot(all(rownames(sample_meta) %in% colnames(expr_mat)))
expr_mat <- expr_mat[, rownames(sample_meta), drop = FALSE]

# Restrict the gene-set list to the sets that were actually scored (= coef_table).
pathways_scored <- pathways[names(pathways) %in% coef_table$set_name]
idx <- limma::ids2indices(pathways_scored, rownames(expr_mat))

# --- (a) CAMERA: competitive, inter-gene-correlation aware, interaction contrast
# inter.gene.cor = NA estimates the ACTUAL per-set correlation (vs the fixed 0.01
# preset), so the test genuinely corrects for the set overlap that inflates the
# per-set GSVA "consistency"; the estimate is reported in camera_cor.
cam <- limma::camera(expr_mat, idx, design = X, contrast = int_col,
                     inter.gene.cor = NA)
cam_df <- tibble::as_tibble(cam, rownames = "set_name") |>
  dplyr::select(set_name,
                camera_ngenes = NGenes, camera_cor = Correlation,
                camera_dir = Direction, camera_p = PValue, camera_fdr = FDR)

# --- (a) ROAST: self-contained rotation test (does the set move at all?)
roa <- limma::mroast(expr_mat, idx, design = X, contrast = int_col,
                     nrot = 9999, set.statistic = "mean")
roa_df <- tibble::as_tibble(roa, rownames = "set_name") |>
  dplyr::select(set_name,
                roast_dir = Direction, roast_p = PValue, roast_fdr = FDR)

# --- Category-level CAMERA: union meta-set per category = one powered test each
cats_all  <- sort(unique(coef_table$category_primary))
meta_sets <- lapply(cats_all, function(cc) {
  sets <- coef_table$set_name[coef_table$category_primary == cc]
  unique(unlist(pathways_scored[sets], use.names = FALSE))
})
names(meta_sets) <- cats_all
meta_idx  <- limma::ids2indices(meta_sets, rownames(expr_mat))
cat_camera <- tibble::as_tibble(
  limma::camera(expr_mat, meta_idx, design = X, contrast = int_col,
                inter.gene.cor = NA),
  rownames = "category_primary"
)

# --- (b) Cross-check: per-set mean DESeq2 interaction Wald stat vs GSVA beta_int
int_df <- tibble::as_tibble(as.data.frame(interaction_raw), rownames = "ensembl_gene_id")
ens2sym <- ortholog |>
  dplyr::select(ensembl_gene_id, external_gene_name) |>
  dplyr::filter(!is.na(external_gene_name), external_gene_name != "") |>
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
# Map ENSMUSG -> symbol; collapse duplicate symbols to the highest-baseMean gene
# (matches the most-expressed-transcript rule used to build expr_mat in 15).
gene_stat <- int_df |>
  dplyr::inner_join(ens2sym, by = "ensembl_gene_id") |>
  dplyr::filter(!is.na(stat)) |>
  dplyr::arrange(external_gene_name, dplyr::desc(baseMean)) |>
  dplyr::distinct(external_gene_name, .keep_all = TRUE)
stat_by_symbol <- stats::setNames(gene_stat$stat, gene_stat$external_gene_name)

mean_int_stat <- vapply(pathways_scored, function(g) {
  s <- stat_by_symbol[intersect(g, names(stat_by_symbol))]
  if (length(s) == 0) NA_real_ else mean(s)
}, numeric(1))
set_mean_stat <- tibble::tibble(set_name = names(mean_int_stat),
                                mean_int_stat = unname(mean_int_stat))

# Fold the gene-level results into coef_table; BH camera within category too.
coef_table <- coef_table |>
  dplyr::left_join(cam_df, by = "set_name") |>
  dplyr::left_join(roa_df, by = "set_name") |>
  dplyr::left_join(set_mean_stat, by = "set_name") |>
  dplyr::group_by(category_primary) |>
  dplyr::mutate(camera_padj_within_cat = stats::p.adjust(camera_p, method = "BH")) |>
  dplyr::ungroup()

# Convergence: GSVA beta_int (d12 - d6) should track the count-level interaction.
cc_df <- coef_table |> dplyr::filter(!is.na(mean_int_stat), !is.na(beta_int))
overall_cc <- suppressWarnings(
  stats::cor.test(cc_df$beta_int, cc_df$mean_int_stat, method = "spearman"))
per_cat_cc <- cc_df |>
  dplyr::group_by(category_primary) |>
  dplyr::summarise(
    n   = dplyr::n(),
    rho = suppressWarnings(stats::cor(beta_int, mean_int_stat, method = "spearman")),
    .groups = "drop"
  )

p_cc <- ggplot2::ggplot(cc_df,
    ggplot2::aes(x = mean_int_stat, y = beta_int, colour = category_primary)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  ggplot2::geom_point(alpha = 0.6, size = 1) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "black",
                       linewidth = 0.5, formula = y ~ x) +
  ggplot2::labs(
    title = "GSVA trajectory vs count-level interaction (convergence check)",
    subtitle = sprintf("Spearman rho = %.3f, p = %.2g (all %d sets)",
                       overall_cc$estimate, overall_cc$p.value, nrow(cc_df)),
    x = "mean DESeq2 interaction Wald stat (member genes)",
    y = "GSVA beta_int (d12 - d6)") +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(
  file.path(out_dir, "crosscheck_beta_int_vs_gene_interaction.pdf"),
  p_cc, width = 8, height = 5.5)

message(sprintf("CAMERA/ROAST done; cross-check Spearman rho = %.3f (p = %.2g)",
                overall_cc$estimate, overall_cc$p.value))

# =============================================================================
# PART 6b: BETWEEN-CATEGORY TRAJECTORY CONTRAST (correlation-aware permutation)
# =============================================================================
# The per-set "vs zero" interaction is underpowered; the powered, decision-
# relevant question is whether program CLASSES follow DIFFERENT trajectories
# (biogenesis/Myc DECLINE vs lineage RISE). This is a large effect (mean
# interaction-t ~ -0.8 vs +0.9), tested at the gene level on PRE-DEFINED library
# categories (not data-selected sets -> not circular).
#
# Metric per gene: the VST-based interaction t (same substrate as the GSVA
# scores), by fast OLS with the treatment-contrast design X (int_col). Null:
# permute SAMPLE labels (rotate the design) so gene-gene correlation is preserved
# -- gene-label shuffling would be anticonservative under the 0.25-0.36 inter-
# gene correlation we measured. The statistic is a DIFFERENCE of category means,
# so any main-effect leakage under permutation cancels between the two groups.

perm_trajectory_contrast <- function(genesA, genesB, B = 4999, seed = 1) {
  gA <- setdiff(intersect(genesA, rownames(expr_mat)), genesB)  # exclusive to A
  gB <- setdiff(intersect(genesB, rownames(expr_mat)), genesA)  # exclusive to B
  Ysub  <- t(expr_mat[c(gA, gB), , drop = FALSE])               # samples x genes
  lab_A <- seq_along(gA); lab_B <- length(gA) + seq_along(gB)
  XtXi  <- solve(crossprod(X))          # invariant to row permutation of X
  vdiag <- XtXi[int_col, int_col]
  df_r  <- nrow(X) - ncol(X)
  tvec  <- function(design) {
    b   <- XtXi %*% crossprod(design, Ysub)
    res <- Ysub - design %*% b
    s2  <- colSums(res^2) / df_r
    b[int_col, ] / sqrt(vdiag * s2)
  }
  del   <- function(tv) mean(tv[lab_A]) - mean(tv[lab_B])
  d_obs <- del(tvec(X))
  set.seed(seed)
  d_null <- vapply(seq_len(B),
                   function(i) del(tvec(X[sample.int(nrow(X)), , drop = FALSE])),
                   numeric(1))
  list(nA = length(gA), nB = length(gB), delta_obs = d_obs,
       p_perm = (1 + sum(abs(d_null) >= abs(d_obs))) / (B + 1), d_null = d_null)
}

onco_metab <- unique(c(meta_sets[["Biogenesis_discrimination"]],
                       meta_sets[["MYC_signatures"]],
                       meta_sets[["MitoCarta"]]))
contr_specs <- list(
  "Biogenesis_discrimination vs Mammary_development" =
    list(meta_sets[["Biogenesis_discrimination"]], meta_sets[["Mammary_development"]]),
  "MYC_signatures vs Mammary_development" =
    list(meta_sets[["MYC_signatures"]], meta_sets[["Mammary_development"]]),
  "Oncogenic-metabolic vs Mammary_development" =
    list(onco_metab, meta_sets[["Mammary_development"]])
)
perm_res <- lapply(contr_specs, function(s) perm_trajectory_contrast(s[[1]], s[[2]]))
perm_contrasts <- tibble::tibble(
  contrast  = names(perm_res),
  nA        = vapply(perm_res, `[[`, integer(1), "nA"),
  nB        = vapply(perm_res, `[[`, integer(1), "nB"),
  delta_obs = vapply(perm_res, `[[`, numeric(1), "delta_obs"),
  p_perm    = vapply(perm_res, `[[`, numeric(1), "p_perm")
)

prim <- perm_res[[1]]
p_perm_plot <- ggplot2::ggplot(data.frame(d = prim$d_null), ggplot2::aes(x = d)) +
  ggplot2::geom_histogram(bins = 60, fill = "grey80", colour = "grey60") +
  ggplot2::geom_vline(xintercept = prim$delta_obs, colour = "#D73027", linewidth = 1) +
  ggplot2::labs(
    title = "Between-category trajectory contrast (permutation null)",
    subtitle = sprintf("Biogenesis vs Mammary-development: observed delta = %.3f, p = %.4f (%d vs %d genes)",
                       prim$delta_obs, prim$p_perm, prim$nA, prim$nB),
    x = "mean interaction-t difference (biogenesis - mammary) under label rotation",
    y = "count") +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "perm_contrast_biogenesis_vs_mammary.pdf"),
                p_perm_plot, width = 7.5, height = 4.5)

message("Permutation contrasts:")
print(perm_contrasts)

# =============================================================================
# PART 6c: IS THE DEVELOPMENTAL TRAJECTORY THE 'CULPRIT'? (per-sample coupling)
# =============================================================================
# Constant Myc drive (MYC_signatures stable) means development does not inhibit
# Myc ACTIVITY. Test instead whether developmental re-differentiation ABSORBS the
# Myc biogenesis advantage: (i) across samples, is a biogenesis composite anti-
# correlated with a developmental composite; (ii) does adjusting for the
# developmental score REMOVE the 6W->12W biogenesis decline within Myc+ (i.e.
# development statistically accounts for the window closing)?

cat_of   <- stats::setNames(set_meta$category_primary, set_meta$set_name)
sets_bio <- rownames(scores)[cat_of[rownames(scores)] == "Biogenesis_discrimination"]
sets_dev <- rownames(scores)[cat_of[rownames(scores)] == "Mammary_development"]
bio_score <- colMeans(scores[sets_bio, , drop = FALSE])   # per-sample composite
dev_score <- colMeans(scores[sets_dev, , drop = FALSE])

mech_df <- tibble::tibble(
  sample     = colnames(scores),
  group      = sample_meta$group,
  myc_status = sample_meta$myc_status,
  timepoint  = sample_meta$timepoint,
  bio        = bio_score,
  dev        = dev_score
)

mech_cor <- mech_df |>
  dplyr::group_by(myc_status) |>
  dplyr::summarise(
    n           = dplyr::n(),
    rho_bio_dev = suppressWarnings(stats::cor(bio, dev, method = "spearman")),
    .groups     = "drop"
  )

# Within Myc+: does the developmental composite account for the biogenesis drop?
posd <- mech_df |> dplyr::filter(myc_status == "pos")
tc   <- list(timepoint = "contr.treatment")   # force treatment coding for coef name
m1   <- stats::lm(bio ~ timepoint, data = posd, contrasts = tc)
m2   <- stats::lm(bio ~ timepoint + dev, data = posd, contrasts = tc)
tp_effect <- c(unadjusted   = unname(stats::coef(m1)["timepoint12W"]),
               adj_for_dev  = unname(stats::coef(m2)["timepoint12W"]))

p_mech <- ggplot2::ggplot(mech_df, ggplot2::aes(x = dev, y = bio, colour = group)) +
  ggplot2::geom_point(size = 2) +
  ggplot2::labs(
    title = "Biogenesis vs developmental program (per sample)",
    subtitle = "Absorption predicts anti-correlation; dev covariate should shrink the Myc+ 6W->12W biogenesis drop",
    x = "Mammary_development composite GSVA", y = "Biogenesis_discrimination composite GSVA") +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "mechanism_biogenesis_vs_development.pdf"),
                p_mech, width = 7, height = 5)

message(sprintf("Myc+ biogenesis 6W->12W effect: unadjusted %.3f, adj-for-dev %.3f",
                tp_effect["unadjusted"], tp_effect["adj_for_dev"]))

# =============================================================================
# PART 7: SAVE OVERVIEW OBJECT
# =============================================================================

overview_out <- list(
  coef_table            = coef_table,
  cond_means            = cond_means,
  group_order           = group_order,
  trajectory_categories = present,
  cat_camera            = cat_camera,
  crosscheck            = list(overall = overall_cc, per_category = per_cat_cc),
  perm_contrasts        = perm_contrasts,
  mechanism             = list(cor = mech_cor, tp_effect = tp_effect, data = mech_df),
  n_sets                = nrow(scores),
  notes                 = paste(
    "coef_table: per-set lm(score ~ timepoint*myc_status) OLS (beta_int = d12 - d6,",
    "the Gate-1 attenuation / rank metric; beta_time = WT slope; myc_slope = Myc+",
    "slope; int_padj_within_cat = BH within category). PLUS gene-level validation:",
    "camera_* (competitive, correlation-aware) + roast_* (self-contained) on the",
    "interaction contrast; mean_int_stat = mean DESeq2 interaction Wald of member",
    "genes; cat_camera = category union meta-set CAMERA; crosscheck = Spearman of",
    "beta_int vs mean_int_stat. Gene-level tests are the trustworthy layer; the",
    "per-set GSVA lm is underpowered / overlap-inflated and directional only.")
)
saveRDS(overview_out, here::here("results", "gsva_overview.rds"))
message("Saved results/gsva_overview.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  ov <- readRDS(here::here("results", "gsva_overview.rds"))
  ct <- ov$coef_table

  # Coding sanity: d6 equals the 6W group-mean gap (already stopifnot'd above)
  cm <- ov$cond_means
  d6_check <- cm[ct$set_name, "6W_pos"] - cm[ct$set_name, "6W_neg"]
  summary(unname(d6_check) - ct$d6)          # should be ~0

  # Per trajectory category: top movers and which mechanism dominates.
  # H1 signature: beta_int < 0 with Myc+ slope (myc_slope) < 0 and beta_time ~ 0.
  # H3 signature: beta_int < 0 with beta_time > 0 and myc_slope ~ 0.
  for (cat in ov$trajectory_categories) {
    cat("\n====", cat, "====\n")
    ct |>
      dplyr::filter(category_primary == cat) |>
      dplyr::arrange(int_p) |>
      dplyr::mutate(
        mech = dplyr::case_when(
          beta_int < 0 & myc_slope < 0 & abs(beta_time) < abs(myc_slope) ~ "H1 (Myc+ down)",
          beta_int < 0 & beta_time > 0 & abs(myc_slope) < abs(beta_time) ~ "H3 (WT up)",
          beta_int < 0 ~ "attenuates (mixed)",
          TRUE ~ "grows/other"
        )
      ) |>
      dplyr::select(set_name, d6, d12, beta_time, myc_slope, beta_int,
                    int_p, camera_p, camera_fdr, roast_p, mean_int_stat, mech) |>
      head(12) |>
      print()
  }

  # Landscape row counts after the top-35 cut (readable?)
  ct |> dplyr::count(category_primary) |> print(n = Inf)

  # --- Gene-level validation (the trustworthy layer) ---
  # Category-level CAMERA: one powered, correlation-aware test per category.
  ov$cat_camera |> dplyr::arrange(PValue) |> print(n = Inf)

  # How many individual sets survive CAMERA FDR < 0.05 (vs zero for the GSVA lm)?
  ct |> dplyr::summarise(
    camera_sig   = sum(camera_fdr < 0.05, na.rm = TRUE),
    roast_sig    = sum(roast_fdr  < 0.05, na.rm = TRUE),
    lm_int_sig   = sum(int_padj_within_cat < 0.05, na.rm = TRUE)
  ) |> print()

  # Convergence: does the GSVA trajectory track the count-level interaction?
  ov$crosscheck$overall                      # overall Spearman rho + p
  ov$crosscheck$per_category |> dplyr::arrange(dplyr::desc(abs(rho))) |> print(n = Inf)

  # Do CAMERA and the GSVA lm agree on the movers? (rank concordance)
  suppressWarnings(stats::cor(ct$camera_p, ct$int_p, method = "spearman",
                              use = "complete.obs"))

  # --- The powered claim: do program CLASSES move oppositely? (permutation) ---
  ov$perm_contrasts |> print()             # delta_obs strongly negative + p_perm

  # --- Is development the 'culprit'? (per-sample coupling) ---
  ov$mechanism$cor                         # bio~dev Spearman per genotype (expect neg, esp pos)
  ov$mechanism$tp_effect                   # Myc+ biogenesis 6W->12W: unadjusted vs adj-for-dev
  # If |adj_for_dev| << |unadjusted|, the developmental composite statistically
  # accounts for the Myc+ biogenesis window closing (absorption).

  # Confirm output PDFs exist (incl. cross-check, permutation, mechanism)
  list.files(here::here("outputs", "gsva_overview"), pattern = "\\.pdf$")
}
