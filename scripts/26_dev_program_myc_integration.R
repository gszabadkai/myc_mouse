# scripts/26_dev_program_myc_integration.R
# =============================================================================
# How does Myc integrate with / modify the background developmental program?
# (Block A revision, Issue #1 -- docs/2026-07-08_BlockA_revision_plan.md)
# =============================================================================
#
# A REFRAME on already-scored developmental sets (no new GSVA/fGSEA/DESeq run).
# Supersedes the superficial answer in 19_dev_composition.R, which (a) never asked
# the integration question, (b) collapsed ~179 signatures to 5 lineage composites,
# (c) reduced differentiation to a single HR-LP scalar, and (d) averaged the Myc
# effect over time as (d6+d12)/2 -- discarding the timepoint resolution.
#
# Here every MG_* developmental set is resolved INDIVIDUALLY, labelled by source +
# putative state + modality, across the FOUR-CONTRAST trajectory, dual lens:
#   GSVA (per-sample cohort-relative STATE) + fGSEA (importance vs transcriptome).
# Canonical contrasts (both lenses):
#   WT_6W->12W  = beta_time / timepoint_neg   (the background developmental axis)
#   Myc+_6W->12W= myc_slope / timepoint_pos   (where the Myc+ gland travels)
#   Myc@6W      = d6        / myc_6W          (Myc displacement already at 6W)
#   Myc@12W     = d12       / myc_12W         (where it lands)
#   interaction = beta_int  / interaction
#
# The headline reframe question: where is the WT developmental trajectory, how is
# Myc+ ALREADY DIVERGED at 6W, and where does 6W_pos -> 12W_pos land -- does Myc
# ACCELERATE (parallel to the WT axis), DIVERT (orthogonal), or REVERSE (against) it?
#
# Parts:
#   0  set annotation from the curated data/dev_mec_annotation.csv (final state /
#      state_fine subcategory / lineage / modality) -- all 179 sets.
#   A  four-contrast dual-lens matrix + paired nets (UP-DN, OPEN-CLOSED).
#   B  WT-axis + Myc-divergence geometry (main-state convergence; other for observation).
#   C  three-state MEC profiles (BMYO/LASP/LHS) + other-subgroup profiles.
#   D  per-source (main) + per-subgroup (other) trajectories, heatmaps, convergence,
#      pair nets.
#
# State uses the consensus MEC nomenclature (Gray, Kessenbrock, Khaled et al.,
# Dev Cell 2025; docs/MEC_types.pdf): three phenotypically discontinuous resting-
# gland epithelial types -- BMYO (basal-myoepithelial), LASP (luminal adaptive
# secretory precursor), LHS (luminal hormone-sensing). The annotation is the author's
# hand curation (data/dev_mec_annotation.csv): state = final MEC type or 'other';
# state_fine = analysis SUBCATEGORY (source/lineage for the mains; SIG/MATRIX/IMMUNE/
# PAL2017/HENRY/SCHEELE_pub/CHUNG_C6/FETAL and the GRAY directional families for
# 'other'); lineage = the consensus lineage a set represents even when demoted to
# 'other' (used for the paired nets).
#
# Ceiling (notes): per-set GSVA contrasts are the powered layer (per-sample, all 24);
# fGSEA NES adds vs-transcriptome importance; n=6/group -> within-set interaction is
# directional. Association, not causation.
#
# Input:  results/gsva_scores.rds      (per-sample GSVA matrix + sample_meta)
#         results/gsva_overview.rds     (coef_table: beta_time/d6/d12/beta_int/myc_slope)
#         results/fgsea_percategory.rds (per-set NES + padj, 5 contrasts)
#         data/dev_mec_annotation.csv   (author-curated MEC annotation)
# Output: results/dev_program_myc_integration.rds
#         outputs/dev_program_myc_integration/*.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

out_dir <- here::here("outputs", "dev_program_myc_integration")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# =============================================================================
# PART 1: LOAD + contrast crosswalk
# =============================================================================

gsva   <- readRDS(here::here("results", "gsva_scores.rds"))
gov    <- readRDS(here::here("results", "gsva_overview.rds"))
fgp    <- readRDS(here::here("results", "fgsea_percategory.rds"))$fgsea

scores      <- gsva$scores
sample_meta <- as.data.frame(gsva$sample_meta)
sample_meta <- sample_meta[colnames(scores), , drop = FALSE]
stopifnot(identical(rownames(sample_meta), colnames(scores)))
sample_meta$timepoint  <- factor(sample_meta$timepoint,  levels = c("6W", "12W"))
sample_meta$myc_status <- factor(sample_meta$myc_status, levels = c("neg", "pos"))
sample_meta$group      <- factor(sample_meta$group,
                                 levels = c("6W_neg", "6W_pos", "12W_neg", "12W_pos"))

coef_tbl <- gov$coef_table
dev_sets <- grep("^MG_", rownames(scores), value = TRUE)   # all 179 developmental
message(sprintf("Developmental MG_* sets scored (GSVA): %d", length(dev_sets)))

# fGSEA carries the 5 contrasts as 'ranking'; canonicalise to the shared vocabulary.
canon_contrast <- function(x) {
  dplyr::case_when(
    x %in% c("myc_6W")        ~ "Myc@6W",
    x %in% c("myc_12W")       ~ "Myc@12W",
    x %in% c("timepoint_neg") ~ "WT_6W->12W",
    x %in% c("timepoint_pos") ~ "Myc+_6W->12W",
    x %in% c("interaction")   ~ "interaction",
    TRUE ~ NA_character_
  )
}
contrast_levels <- c("WT_6W->12W", "Myc@6W", "Myc@12W", "Myc+_6W->12W", "interaction")

# =============================================================================
# PART 0: SET ANNOTATION -- author-curated MEC consensus states + subcategories
# =============================================================================
# Authoritative annotation curated by the author (data/dev_mec_annotation.csv;
# seeded from an algorithmic first pass, then hand-corrected against the catalog and
# the Gray 2025 consensus). Columns:
#   state_fine : ANALYSIS SUBCATEGORY. Main states -> source/lineage tag; 'other' ->
#     the meaningful subgroup (SIG, MATRIX, IMMUNE, PAL2017, HENRY, SCHEELE_pub,
#     CHUNG_C6, FETAL, and the GRAY directional families X_HEVSLE / X_TEB_VS_DUCTAL).
#   state      : FINAL consensus MEC type (BMYO/LASP/LHS) or 'other' (= new_state
#     where the author demoted/promoted a set, else the collapsed consensus).
#   lineage    : the consensus lineage a set represents EVEN IF demoted to 'other'
#     (so the directional UP-DN and CHUNG OPEN-CLOSED pairs keep BMYO/LASP/LHS).
#   modality   : measurement type (state-marker / trajectory-cluster / directional-DE
#     / ATAC / pubertal / embryonic-fetal / contamination-control).

anno_path    <- here::here("data", "dev_mec_annotation.csv")
stopifnot(file.exists(anno_path))
main_states  <- c("BMYO", "LASP", "LHS")
state_levels <- c(main_states, "other")

annot <- readr::read_csv(anno_path, show_col_types = FALSE) |>
  dplyr::mutate(
    new_state = dplyr::na_if(as.character(new_state), ""),
    lineage   = dplyr::case_when(               # uses the pre-coalesce 'state'
      new_state %in% main_states ~ new_state,
      state     %in% main_states ~ state,
      TRUE ~ "other"),
    state     = dplyr::coalesce(new_state, state)) |>   # FINAL state
  dplyr::select(set, source, state_fine, modality, state, lineage, note)

# Coverage: every scored dev set must be annotated (and vice versa).
miss  <- setdiff(dev_sets, annot$set)
extra <- setdiff(annot$set, dev_sets)
if (length(miss))  stop("Scored but unannotated: ", paste(miss, collapse = ", "))
if (length(extra)) message(sprintf("Annotated but not scored (dropped): %s",
                                    paste(extra, collapse = ", ")))
annot <- annot |> dplyr::filter(set %in% dev_sets)
annot$state   <- factor(annot$state,   levels = state_levels)
annot$lineage <- factor(annot$lineage, levels = state_levels)

readr::write_csv(annot, file.path(out_dir, "annot_full.csv"))
message(sprintf("Annotation: %d sets | %s", nrow(annot),
                paste(sprintf("%s=%d", levels(annot$state), as.integer(table(annot$state))),
                      collapse = " ")))

# =============================================================================
# PART A: FOUR-CONTRAST DUAL-LENS MATRIX  (replaces (d6+d12)/2)
# =============================================================================

# GSVA lens: coef_table columns -> canonical contrasts (long).
gsva_long <- coef_tbl |>
  dplyr::filter(set_name %in% dev_sets) |>
  dplyr::transmute(set = set_name,
                   `WT_6W->12W`   = beta_time,
                   `Myc+_6W->12W` = myc_slope,
                   `Myc@6W`       = d6,
                   `Myc@12W`      = d12,
                   `interaction`  = beta_int) |>
  tidyr::pivot_longer(-set, names_to = "contrast", values_to = "gsva_effect")

# fGSEA lens: per-set NES + padj (167 sets; the ~12 *_GSVA sets have no NES).
fgsea_long <- fgp |>
  dplyr::filter(category == "03_mammary_development", grepl("^MG_", pathway)) |>
  dplyr::transmute(set = pathway, contrast = canon_contrast(ranking),
                   fgsea_NES = NES, fgsea_padj = padj_within_category) |>
  dplyr::filter(!is.na(contrast))

dual_matrix <- gsva_long |>
  dplyr::left_join(fgsea_long, by = c("set", "contrast")) |>
  dplyr::left_join(annot, by = "set") |>
  dplyr::mutate(contrast = factor(contrast, levels = contrast_levels)) |>
  dplyr::arrange(source, state, set, contrast)

# Paired-set nets (author request): sets that come as a +/- pair collapse to ONE net
# per contrast per lens = pos - neg, tagged by lineage (the mains BMYO/LASP/LHS, kept
# even where the pair is demoted to state 'other'). Difference not ratio (NES can be
# negative).
#   directional_pairs = GRAY UP - DN  (toward the _UP pole of the lineage contrast)
#   chung_pairs       = CHUNG ATAC OPEN - CLOSED  (toward the open-chromatin pole)
directional_pairs <- dual_matrix |>
  dplyr::filter(modality == "directional-DE", grepl("_(UP|DN)$", set)) |>
  dplyr::mutate(dir  = ifelse(grepl("_UP$", set), "UP", "DN"),
                base = sub("_(UP|DN)$", "", set)) |>
  dplyr::select(base, source, lineage, contrast, dir, gsva_effect, fgsea_NES) |>
  tidyr::pivot_wider(names_from = dir,
                     values_from = c(gsva_effect, fgsea_NES)) |>
  dplyr::mutate(gsva_net  = gsva_effect_UP - gsva_effect_DN,
                fgsea_net = fgsea_NES_UP  - fgsea_NES_DN) |>
  dplyr::filter(!is.na(gsva_net)) |>
  dplyr::arrange(source, lineage, base, contrast)

chung_pairs <- dual_matrix |>
  dplyr::filter(modality == "ATAC", state != "other", grepl("_(OPEN|CLOSED)_", set)) |>
  dplyr::mutate(dir  = ifelse(grepl("_OPEN_", set), "OPEN", "CLOSED"),
                base = sub("_(OPEN|CLOSED)_", "_", set)) |>
  dplyr::select(base, source, lineage, contrast, dir, gsva_effect, fgsea_NES) |>
  tidyr::pivot_wider(names_from = dir,
                     values_from = c(gsva_effect, fgsea_NES)) |>
  dplyr::mutate(gsva_net  = gsva_effect_OPEN - gsva_effect_CLOSED,
                fgsea_net = fgsea_NES_OPEN  - fgsea_NES_CLOSED) |>
  dplyr::filter(!is.na(gsva_net)) |>
  dplyr::arrange(source, lineage, base, contrast)

# Combined tidy net table (for the pair-net figure)
pair_nets <- dplyr::bind_rows(
  directional_pairs |> dplyr::transmute(base, source, lineage, contrast,
                                        gsva_net, fgsea_net, kind = "directional (UP-DN)"),
  chung_pairs       |> dplyr::transmute(base, source, lineage, contrast,
                                        gsva_net, fgsea_net, kind = "ATAC (OPEN-CLOSED)"))

# =============================================================================
# PART B: WT-AXIS + MYC-DIVERGENCE GEOMETRY
# =============================================================================
# Per-set: WT axis (beta_time) vs Myc displacement (d6 at 6W, d12 at 12W). With a
# single set the only relations are ACCELERATE (Myc same sign as WT axis) or REVERSE
# (opposite); ORTHOGONAL/DIVERT is an ENSEMBLE property (correlation of the Myc-effect
# vector against the WT-axis vector across sets) -- computed below.

eps <- 0.02   # WT-axis magnitude below which the developmental direction is ~flat
geom_perset <- coef_tbl |>
  dplyr::filter(set_name %in% dev_sets) |>
  dplyr::transmute(set = set_name, wt_axis = beta_time,
                   myc_6W = d6, myc_12W = d12,
                   mycpos_traj = myc_slope, int = beta_int) |>
  dplyr::left_join(annot, by = "set") |>
  dplyr::mutate(
    class_6W = dplyr::case_when(
      abs(wt_axis) < eps                    ~ "WT-flat",
      sign(myc_6W) == sign(wt_axis)         ~ "accelerate",
      sign(myc_6W) == -sign(wt_axis)        ~ "reverse",
      TRUE ~ "none"),
    class_12W = dplyr::case_when(
      abs(wt_axis) < eps                    ~ "WT-flat",
      sign(myc_12W) == sign(wt_axis)        ~ "accelerate",
      sign(myc_12W) == -sign(wt_axis)       ~ "reverse",
      TRUE ~ "none"),
    # does Myc's displacement grow or shrink from 6W to 12W (re-converge vs diverge)?
    div_change = abs(myc_12W) - abs(myc_6W))

# per-set class tally (by state and by state_fine subcategory), 6W and 12W
class_tally <- dplyr::bind_rows(
  geom_perset |> dplyr::count(state, class_6W)  |> dplyr::rename(class = class_6W)  |> dplyr::mutate(window = "6W"),
  geom_perset |> dplyr::count(state, class_12W) |> dplyr::rename(class = class_12W) |> dplyr::mutate(window = "12W"))
class_tally_fine <- dplyr::bind_rows(
  geom_perset |> dplyr::count(state, state_fine, class_6W)  |> dplyr::rename(class = class_6W)  |> dplyr::mutate(window = "6W"),
  geom_perset |> dplyr::count(state, state_fine, class_12W) |> dplyr::rename(class = class_12W) |> dplyr::mutate(window = "12W"))

# ENSEMBLE convergence vector (honest AP3): correlate Myc effect vs WT axis across sets.
conv_vec <- function(df) {
  ok <- stats::complete.cases(df$wt_axis, df$myc_6W)
  r6  <- suppressWarnings(stats::cor(df$wt_axis[ok], df$myc_6W[ok],  method = "spearman"))
  ok2 <- stats::complete.cases(df$wt_axis, df$myc_12W)
  r12 <- suppressWarnings(stats::cor(df$wt_axis[ok2], df$myc_12W[ok2], method = "spearman"))
  ok3 <- stats::complete.cases(df$wt_axis, df$mycpos_traj)
  rtraj <- suppressWarnings(stats::cor(df$wt_axis[ok3], df$mycpos_traj[ok3], method = "spearman"))
  tibble::tibble(n = sum(ok), rho_myc6_vs_wt = r6, rho_myc12_vs_wt = r12,
                 rho_mycpostraj_vs_wt = rtraj)
}
# CALCULATION excludes 'other' (only the three main MEC states); 'other' shown in the
# plot for observation only.
geom_main      <- geom_perset |> dplyr::filter(state != "other")
conv_overall   <- conv_vec(geom_main)
conv_by_source <- geom_main |> dplyr::group_by(source) |>
  dplyr::group_modify(~ conv_vec(.x)) |> dplyr::ungroup()
conv_by_state  <- geom_main |> dplyr::group_by(state) |>
  dplyr::group_modify(~ conv_vec(.x)) |> dplyr::ungroup()
# observation: per state_fine subcategory (incl. 'other' subgroups)
conv_by_state_fine <- geom_perset |> dplyr::group_by(state, state_fine) |>
  dplyr::group_modify(~ conv_vec(.x)) |> dplyr::ungroup()

# =============================================================================
# PART C: THREE-STATE MEC PROFILES  (BMYO/LASP/LHS; replaces HR-LP scalar)
# =============================================================================
# Per-sample GSVA of every dev set -> mean by (group, state) = the MEC-composition
# PROFILE of each of the 4 conditions across the three consensus types. Read as a
# profile over discontinuous states (NOT an ordered trajectory, per Gray 2025);
# flag cross-source (dis)agreement per state.

scores_long <- as.data.frame(scores[dev_sets, , drop = FALSE]) |>
  tibble::rownames_to_column("set") |>
  tidyr::pivot_longer(-set, names_to = "sample", values_to = "gsva") |>
  dplyr::left_join(annot, by = "set") |>
  dplyr::mutate(group      = sample_meta[sample, "group"],
                timepoint  = sample_meta[sample, "timepoint"],
                myc_status = sample_meta[sample, "myc_status"])

# state x group profile (mean over sets-in-state and samples-in-group)
state_profile <- scores_long |>
  dplyr::filter(state != "other") |>
  dplyr::group_by(group, timepoint, myc_status, state) |>
  dplyr::summarise(mean_gsva = mean(gsva), .groups = "drop")

# cross-source agreement: per source x state x group mean (to eyeball dispersion)
source_state_profile <- scores_long |>
  dplyr::filter(state != "other") |>
  dplyr::group_by(source, state, group) |>
  dplyr::summarise(mean_gsva = mean(gsva), n_sets = dplyr::n_distinct(set), .groups = "drop")

# 'other' subdivided by state_fine subcategory: profile per condition (observation)
subgroup_profile <- scores_long |>
  dplyr::filter(state == "other") |>
  dplyr::group_by(group, timepoint, myc_status, state_fine) |>
  dplyr::summarise(mean_gsva = mean(gsva), n_sets = dplyr::n_distinct(set), .groups = "drop")

# per-set per-group means + se (for the trajectory small-multiples)
set_group_means <- scores_long |>
  dplyr::group_by(set, source, state, state_fine, lineage, modality,
                  group, timepoint, myc_status) |>
  dplyr::summarise(mean_gsva = mean(gsva),
                   se = stats::sd(gsva) / sqrt(dplyr::n()), .groups = "drop")

# =============================================================================
# PART D: FIGURES
# =============================================================================
geno_cols <- c(neg = "#4575B4", pos = "#D73027")

# --- D1: small-multiple trajectories ------------------------------------------
# Shared plotter: facet by a caption column, one PDF per facet-group value.
traj_pdf <- function(d, facet_by, title, fname) {
  if (nrow(d) == 0) return(invisible(NULL))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = timepoint, y = mean_gsva,
                                       colour = myc_status, group = myc_status)) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_gsva - se, ymax = mean_gsva + se),
                           width = 0.12, linewidth = 0.3) +
    ggplot2::facet_wrap(stats::as.formula(paste0("~ ", facet_by, " + set")),
                        scales = "free_y",
                        labeller = ggplot2::labeller(.multi_line = FALSE)) +
    ggplot2::scale_colour_manual(values = geno_cols) +
    ggplot2::labs(title = title,
      subtitle = "per-set GSVA, 4-group; blue = Myc-, red = Myc+ (mean +- se)",
      x = NULL, y = "GSVA score") +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 5.5),
                   legend.position = "top")
  h <- max(4, ceiling(dplyr::n_distinct(d$set) / 4) * 1.6)
  ggplot2::ggsave(file.path(out_dir, fname), p, width = 10, height = h, limitsize = FALSE)
}

# D1a: per-source trajectories -- MAIN states only (BMYO/LASP/LHS).
invisible(lapply(sort(unique(annot$source)), function(src)
  traj_pdf(set_group_means |> dplyr::filter(source == src, state != "other"),
           "state", sprintf("Developmental trajectories (main states) -- %s", src),
           sprintf("trajectories_main_%s.pdf", src))))

# D1b: per-state_fine trajectories -- the 'other' subgroups, analysed separately.
other_subgroups <- sort(unique(as.character(
  set_group_means$state_fine[set_group_means$state == "other"])))
invisible(lapply(other_subgroups, function(sf) {
  safe <- gsub("[^A-Za-z0-9]+", "_", sf)
  traj_pdf(set_group_means |> dplyr::filter(state == "other", state_fine == sf),
           "lineage", sprintf("'other' subgroup trajectories -- %s", sf),
           sprintf("trajectories_other_%s.pdf", safe))
}))

# --- D2: master dual-lens heatmaps (GSVA effect; fGSEA NES) --------------------
set_order <- dual_matrix |> dplyr::distinct(set, state, state_fine, source) |>
  dplyr::arrange(state, state_fine, source, set) |> dplyr::pull(set)
dm <- dual_matrix |> dplyr::mutate(set = factor(set, levels = rev(set_order)))

p_gsva <- ggplot2::ggplot(dm, ggplot2::aes(x = contrast, y = set, fill = gsva_effect)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                                midpoint = 0, name = "GSVA effect") +
  ggplot2::labs(title = "Dev sets x contrast -- GSVA effect (per-sample state)",
    subtitle = "rows grouped by source, state; all 179 MG_* sets", x = NULL, y = NULL) +
  ggplot2::theme_bw(base_size = 6) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                 axis.text.y = ggplot2::element_text(size = 3.5))
ggplot2::ggsave(file.path(out_dir, "heatmap_gsva_effect.pdf"), p_gsva,
                width = 7, height = 20, limitsize = FALSE)

p_fg <- ggplot2::ggplot(dm |> dplyr::filter(!is.na(fgsea_NES)),
                        ggplot2::aes(x = contrast, y = set, fill = fgsea_NES)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                                midpoint = 0, name = "fGSEA NES") +
  ggplot2::labs(title = "Dev sets x contrast -- fGSEA NES (vs transcriptome)",
    subtitle = "167 sets with NES; rows grouped by source, state", x = NULL, y = NULL) +
  ggplot2::theme_bw(base_size = 6) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                 axis.text.y = ggplot2::element_text(size = 3.5))
ggplot2::ggsave(file.path(out_dir, "heatmap_fgsea_nes.pdf"), p_fg,
                width = 7, height = 20, limitsize = FALSE)

# --- D3: three-state MEC profile per group (BMYO/LASP/LHS) --------------------
p_cont <- ggplot2::ggplot(state_profile,
    ggplot2::aes(x = state, y = mean_gsva, colour = myc_status, group = group)) +
  ggplot2::geom_line(ggplot2::aes(linetype = timepoint), linewidth = 0.6) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_colour_manual(values = geno_cols) +
  ggplot2::labs(title = "Three-state MEC profile by condition (Gray 2025 consensus)",
    subtitle = "mean GSVA over sets per MEC type; replaces the single HR-LP scalar",
    x = "consensus MEC type (BMYO / LASP / LHS)", y = "mean GSVA") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "continuum_profile.pdf"), p_cont,
                width = 7.5, height = 5)

# --- D4: ensemble convergence vector (Myc effect vs WT axis, per set) ----------
# rho computed on the three MAIN states only; 'other' shown for observation (squares).
gp <- geom_perset |>
  dplyr::mutate(shape_grp = ifelse(state == "other", "other (excl. from rho)", "main state"))
p_conv <- ggplot2::ggplot(gp,
    ggplot2::aes(x = wt_axis, y = myc_6W, colour = state, shape = shape_grp)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  ggplot2::geom_point(size = 1.8, alpha = 0.8) +
  ggplot2::scale_shape_manual(values = c(`main state` = 16, `other (excl. from rho)` = 15),
                              name = NULL) +
  ggplot2::labs(title = "Myc 6W displacement vs WT developmental axis (per set)",
    subtitle = sprintf("Spearman rho(Myc@6W, WT axis) = %.2f on the 3 main states; dashed = pure acceleration",
                       conv_overall$rho_myc6_vs_wt),
    x = "WT axis (beta_time = 6W_neg -> 12W_neg)", y = "Myc@6W displacement (d6)") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "convergence_vector_myc6.pdf"), p_conv,
                width = 7.5, height = 6)

# --- D5: paired-set net direction (directional UP-DN + CHUNG OPEN-CLOSED) ------
if (nrow(pair_nets) > 0) {
  p_pair <- ggplot2::ggplot(pair_nets,
      ggplot2::aes(x = contrast, y = gsva_net, colour = base, group = base)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_grid(kind ~ lineage) +
    ggplot2::labs(title = "Paired-set net direction (GSVA): UP-DN and OPEN-CLOSED",
      subtitle = "positive = toward the UP / OPEN pole of that lineage; by consensus lineage",
      x = NULL, y = "net GSVA effect (pos - neg)") +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                   legend.position = "right")
  ggplot2::ggsave(file.path(out_dir, "pair_nets.pdf"), p_pair, width = 9, height = 5.5)
}

# =============================================================================
# PART E: SAVE
# =============================================================================
dev_out <- list(
  annot           = annot,
  dual_matrix     = dual_matrix,
  geom_perset     = geom_perset,
  class_tally     = class_tally,
  class_tally_fine = class_tally_fine,
  convergence     = list(overall = conv_overall, by_source = conv_by_source,
                         by_state = conv_by_state, by_state_fine = conv_by_state_fine),
  state_profile   = state_profile,
  source_state_profile = source_state_profile,
  subgroup_profile = subgroup_profile,
  set_group_means = set_group_means,
  directional_pairs = directional_pairs,
  chung_pairs     = chung_pairs,
  pair_nets       = pair_nets,
  notes = paste(
    "Issue #1 reframe: all 179 MG_* dev sets resolved individually. Annotation is the",
    "author's hand curation (data/dev_mec_annotation.csv). state = final consensus MEC",
    "type (Gray/Kessenbrock/Khaled, Dev Cell 2025; docs/MEC_types.pdf): BMYO / LASP /",
    "LHS, or 'other'. state_fine = analysis SUBCATEGORY (mains: source/lineage; other:",
    "SIG/MATRIX/IMMUNE/PAL2017/HENRY/SCHEELE_pub/CHUNG_C6/FETAL + GRAY directional",
    "families). lineage = consensus lineage kept even when demoted to other (for the",
    "nets). Analysed across the 4-contrast trajectory (WT_6W->12W, Myc@6W, Myc@12W,",
    "Myc+_6W->12W + interaction), dual lens GSVA (per-sample state) + fGSEA (vs",
    "transcriptome). Replaces 19_dev_composition's 5-composite / HR-LP scalar /",
    "(d6+d12)/2 collapse. Geometry: per-set accelerate/reverse vs WT axis + ensemble",
    "Spearman convergence (Myc effect vs WT axis) on the 3 MAIN states only ('other'",
    "plotted for observation). directional_pairs = GRAY UP-DN net; chung_pairs =",
    "CHUNG ATAC OPEN-CLOSED net (BASAL/LP/ML -> BMYO/LASP/LHS); both per contrast per",
    "lens, lineage-tagged. 'other' subdivided by state_fine (subgroup_profile + per-",
    "subgroup trajectory PDFs). Ceiling: GSVA per-set contrasts powered (24 samples);",
    "fGSEA adds importance; n=6/group -> interaction directional; association not causation.")
)
saveRDS(dev_out, here::here("results", "dev_program_myc_integration.rds"))
message("Saved results/dev_program_myc_integration.rds")
message(sprintf("Ensemble convergence: rho(Myc@6W,WT)=%.2f  rho(Myc@12W,WT)=%.2f  rho(Myc+traj,WT)=%.2f",
                conv_overall$rho_myc6_vs_wt, conv_overall$rho_myc12_vs_wt,
                conv_overall$rho_mycpostraj_vs_wt))

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  dp <- readRDS(here::here("results", "dev_program_myc_integration.rds"))

  # 0. Annotation: final state counts (expect BMYO 26 / LASP 29 / LHS 29 / other 95),
  #    and the state_fine subcategories within 'other'.
  dp$annot |> dplyr::count(state) |> print(n = Inf)
  dp$annot |> dplyr::filter(state == "other") |> dplyr::count(state_fine) |> print(n = Inf)
  dp$annot |> dplyr::count(state, state_fine) |> print(n = Inf)

  # A. dual-lens matrix: 179 GSVA sets x 5 contrasts; fGSEA present for 167.
  dp$dual_matrix |> dplyr::summarise(
    n_rows = dplyr::n(), n_sets = dplyr::n_distinct(set),
    n_fgsea_sets = dplyr::n_distinct(set[!is.na(fgsea_NES)])) |> print()
  dp$dual_matrix |> dplyr::filter(contrast == "WT_6W->12W", !is.na(fgsea_NES)) |>
    dplyr::summarise(rho = stats::cor(gsva_effect, fgsea_NES, method = "spearman")) |> print()

  # B. geometry (MAIN states only in the rho): 6W orthogonal, 12W oppositional?
  dp$convergence$overall |> print()
  dp$convergence$by_state |> print(n = Inf)
  dp$convergence$by_state_fine |> print(n = Inf)   # incl. 'other' subgroups (observation)
  dp$class_tally |> tidyr::pivot_wider(names_from = class, values_from = n,
                                       values_fill = 0) |> print(n = Inf)

  # C. three-state MEC profile + the 'other' subgroup profiles
  dp$state_profile |> tidyr::pivot_wider(names_from = state, values_from = mean_gsva) |>
    print(n = Inf)
  dp$subgroup_profile |> tidyr::pivot_wider(names_from = state_fine, values_from = mean_gsva) |>
    print(width = Inf)

  # C2. paired nets by lineage (Myc-effect contrasts): UP-DN and OPEN-CLOSED
  dp$directional_pairs |> dplyr::filter(contrast %in% c("Myc@6W", "Myc@12W")) |>
    dplyr::select(base, lineage, contrast, gsva_net, fgsea_net) |> print(n = Inf)
  dp$chung_pairs |>
    dplyr::select(base, lineage, contrast, gsva_net, fgsea_net) |> print(n = Inf)

  # D. outputs (per-source main trajectories + per-subgroup 'other' + heatmaps + nets)
  list.files(here::here("outputs", "dev_program_myc_integration"), pattern = "\\.pdf$")
}
