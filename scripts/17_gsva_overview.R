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
# Input:
#   - results/gsva_scores.rds : list(scores [set x 24], set_meta, sample_meta, ...)
#
# Output:
#   - results/gsva_overview.rds : list(coef_table, cond_means, trajectory_categories)
#   - outputs/gsva_overview/<category>_landscape.pdf        (all categories)
#   - outputs/gsva_overview/<category>_profiles.pdf         (trajectory cats)
#   - outputs/gsva_overview/<category>_dumbbell.pdf         (trajectory cats)
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# =============================================================================
# PART 1: LOAD GSVA SCORES + ALIGN SAMPLE METADATA
# =============================================================================

gsva_out    <- readRDS(here::here("results", "gsva_scores.rds"))
scores      <- gsva_out$scores                       # set x sample (24)
set_meta    <- gsva_out$set_meta                     # set_name / category_primary / ...
sample_meta <- as.data.frame(gsva_out$sample_meta)

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
# PART 6: SAVE OVERVIEW OBJECT
# =============================================================================

overview_out <- list(
  coef_table            = coef_table,
  cond_means            = cond_means,
  group_order           = group_order,
  trajectory_categories = present,
  n_sets                = nrow(scores),
  notes                 = paste(
    "coef_table: per-set lm(score ~ timepoint*myc_status) OLS.",
    "beta_int = d12 - d6 (Gate-1 attenuation, the rank metric);",
    "beta_time = Myc- (WT) slope; myc_slope = Myc+ slope.",
    "int_padj_within_cat = BH within category_primary.")
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
                    int_p, int_padj_within_cat, mech) |>
      head(12) |>
      print()
  }

  # Landscape row counts after the top-35 cut (readable?)
  ct |> dplyr::count(category_primary) |> print(n = Inf)

  # Confirm output PDFs exist
  list.files(here::here("outputs", "gsva_overview"), pattern = "\\.pdf$")
}
