# scripts/22_reframe_mtdna_abund_retention.R
# =============================================================================
# Reframe: mtDNA / abundance / retention -> Supp 3 inputs (Block A, Day 4-5)
# =============================================================================
#
# The final piece of the action-point menu. A pure REFRAME: reads existing
# .rds (mitoPPS, per-category fGSEA), computes no new model. Three lenses on the
# mitochondrial signal, each feeding a Supp-3 panel:
#
#   AP-mtDNA  -- mtDNA-encoded vs nuclear OXPHOS split. mitoPPS already isolates
#     the 13 mtDNA-encoded genes into the synthetic pathway
#     "mtDNA-encoded OXPHOS subunits" (script 08), with mt-* stripped from the
#     nuclear complexes. Here: per-group trajectory of the synthetic mtDNA
#     pathway vs the nuclear-OXPHOS composite (corroborates script 19 AP1 -- WT
#     mtDNA-OXPHOS UP, nuclear subunits DOWN across 6W->12W), plus a per-complex
#     mtDNA-vs-nuclear correlation table (per-complex detail is cut-line 5, so
#     table not figure).
#
#   AP-abund  -- fGSEA (rank-relative, ABSOLUTE enrichment) vs mitoPPS
#     (within-compartment REPRIORITISATION) are complementary lenses that can
#     DISAGREE. The answer is TIME-DEPENDENT: at 6W the Myc effect is a selective
#     reprioritisation (mitoPPS diffs widely spread across mito pathways), by 12W
#     it is more uniform biogenesis (diffs narrower / shifted together). We
#     quantify the reprioritisation dispersion per timepoint and annotate with
#     the fGSEA compartment NES. Do NOT collapse to one word.
#
#   AP-retention -- how much of the WT developmental enrichment is "retained"
#     under Myc across time. retention_abs = NES_pos - NES_neg (the meaningful,
#     absolute one) and retention_ratio = NES_pos / NES_neg on the temporal
#     contrasts (timepoint_pos = Myc+ 6W->12W, timepoint_neg = WT 6W->12W), per
#     category. FELSHER GUARDRAIL: a declining Myc-signature NES over time reads
#     as "developmental decline of the substrate + a STABLE Myc offset", NOT loss
#     of Myc activity -- Myc dose is constant (mRNA/blot/literature, CLAUDE.md).
#
# Input:  results/mitopps_scores.rds (raw_pathway_scores, mitopps_scores,
#           mitopps_pairwise, pathway_tier1_map, mtdna_pathway_name),
#         results/fgsea_percategory.rds (script 20 tidy fGSEA),
#         results/fgsea_results.rds + results/fgsea_xs_results.rds (legacy pooled,
#           optional comparator -- loaded defensively).
# Output: results/reframe_supp3.rds; outputs/reframe_supp3/*.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

out_dir <- here::here("outputs", "reframe_supp3")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# =============================================================================
# PART 1: LOAD
# =============================================================================

mps  <- readRDS(here::here("results", "mitopps_scores.rds"))
fgp  <- readRDS(here::here("results", "fgsea_percategory.rds"))$fgsea

# Legacy pooled fGSEA -- comparator only; structure not relied on. Load
# defensively so an unexpected shape never breaks the reframe.
legacy_pooled <- tryCatch(readRDS(here::here("results", "fgsea_results.rds")),
                          error = function(e) NULL)
legacy_xs     <- tryCatch(readRDS(here::here("results", "fgsea_xs_results.rds")),
                          error = function(e) NULL)

group_levels <- c("6W_neg", "6W_pos", "12W_neg", "12W_pos")

# Per-sample score frames (samples x pathways + group/timepoint/myc_status).
# Pathway columns are exactly the numeric columns (metadata cols are character).
raw_scores <- mps$raw_pathway_scores
pps_scores <- mps$mitopps_scores
num_cols   <- function(df) names(df)[vapply(df, is.numeric, logical(1))]
raw_paths  <- num_cols(raw_scores)
pps_paths  <- num_cols(pps_scores)

tier1     <- mps$pathway_tier1_map           # named vector: pathway -> Tier1
mtdna_nm  <- mps$mtdna_pathway_name           # "mtDNA-encoded OXPHOS subunits"
stopifnot(!is.null(mtdna_nm), mtdna_nm %in% pps_paths)

# OXPHOS-tier pathways present as columns; nuclear = OXPHOS-tier minus synthetic mtDNA
oxphos_all     <- intersect(names(tier1)[tier1 == "OXPHOS"], pps_paths)
oxphos_nuclear <- setdiff(oxphos_all, mtdna_nm)
message(sprintf("OXPHOS-tier pathways: %d (nuclear %d + synthetic mtDNA)",
                length(oxphos_all), length(oxphos_nuclear)))
stopifnot(length(oxphos_nuclear) >= 4)

# =============================================================================
# PART 2: AP-mtDNA -- mtDNA-encoded vs nuclear OXPHOS split
# =============================================================================

# Per-group means of the synthetic mtDNA pathway and the nuclear-OXPHOS composite
# (mean across nuclear OXPHOS pathways), on BOTH the raw score and the mitoPPS.
grp_mean <- function(df, cols, comp_name) {
  df |>
    dplyr::mutate(.comp = rowMeans(dplyr::across(dplyr::all_of(cols)), na.rm = TRUE)) |>
    dplyr::group_by(group) |>
    dplyr::summarise(mean_score = mean(.comp), .groups = "drop") |>
    dplyr::mutate(component = comp_name)
}

mtdna_traj <- dplyr::bind_rows(
  grp_mean(pps_scores, mtdna_nm,       "mtDNA-encoded (mitoPPS)"),
  grp_mean(pps_scores, oxphos_nuclear, "nuclear OXPHOS (mitoPPS)"),
  grp_mean(raw_scores, mtdna_nm,       "mtDNA-encoded (raw)"),
  grp_mean(raw_scores, oxphos_nuclear, "nuclear OXPHOS (raw)")
) |>
  tidyr::separate(group, into = c("timepoint", "myc_status"), sep = "_", remove = FALSE) |>
  dplyr::mutate(group      = factor(group, levels = group_levels),
                timepoint  = factor(timepoint, levels = c("6W", "12W")),
                myc_status = factor(myc_status, levels = c("neg", "pos")))

# Per-complex mtDNA-vs-nuclear correlation across the 24 samples (and within WT --
# the substrate the acute pulse hits). Each nuclear OXPHOS pathway vs the mtDNA
# pathway; mitoPPS scale (within-compartment reallocation).
mtdna_vec_all <- pps_scores[[mtdna_nm]]
wt_idx        <- pps_scores$myc_status == "neg"
per_complex_cor <- purrr::map_dfr(oxphos_nuclear, function(pw) {
  v <- pps_scores[[pw]]
  ct_all <- suppressWarnings(stats::cor.test(mtdna_vec_all, v, method = "pearson"))
  ct_wt  <- suppressWarnings(stats::cor.test(mtdna_vec_all[wt_idx], v[wt_idx],
                                             method = "pearson"))
  tibble::tibble(nuclear_pathway = pw,
                 r_all = unname(ct_all$estimate), p_all = ct_all$p.value,
                 r_wt  = unname(ct_wt$estimate),  p_wt  = ct_wt$p.value)
}) |>
  dplyr::arrange(r_all)

# Plot: mtDNA vs nuclear OXPHOS trajectory, mitoPPS lens (the reallocation view)
p_mtdna <- mtdna_traj |>
  dplyr::filter(grepl("mitoPPS", component)) |>
  ggplot2::ggplot(ggplot2::aes(x = timepoint, y = mean_score,
                               colour = component, group = component)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::facet_wrap(~ myc_status, labeller = ggplot2::labeller(
    myc_status = c(neg = "Myc- (WT substrate)", pos = "Myc+"))) +
  ggplot2::scale_colour_manual(values = c(
    "mtDNA-encoded (mitoPPS)"  = "#B2182B",
    "nuclear OXPHOS (mitoPPS)" = "#2166AC")) +
  ggplot2::labs(
    title = "AP-mtDNA: mtDNA-encoded vs nuclear OXPHOS reprioritisation (mitoPPS)",
    subtitle = "WT 6W->12W: mtDNA-OXPHOS up, nuclear subunits down (cf. script 19 AP1)",
    x = NULL, y = "mean mitoPPS", colour = NULL) +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "mtdna_vs_nuclear_trajectory.pdf"),
                p_mtdna, width = 8, height = 4.5)

# =============================================================================
# PART 3: AP-abund -- fGSEA (rank-relative) vs mitoPPS (within-compartment) lens
# =============================================================================
# The two lenses answer different questions; their agreement is TIME-DEPENDENT.
# mitoPPS reprioritisation dispersion (spread of Myc-effect diffs across mito
# pathways) is the operational "selective (6W) vs uniform (12W)" measure; the
# fGSEA compartment NES is the absolute-enrichment annotation.

pw_pairs <- mps$mitopps_pairwise         # pathway, contrast, diff, padj, ...
mito_pps_paths <- intersect(oxphos_all,  pps_paths)   # OXPHOS-tier for the compartment view

abund_disp <- pw_pairs |>
  dplyr::filter(contrast %in% c("Myc_effect_6W", "Myc_effect_12W"),
                pathway %in% mito_pps_paths) |>
  dplyr::mutate(timepoint = ifelse(contrast == "Myc_effect_6W", "6W", "12W")) |>
  dplyr::group_by(timepoint) |>
  dplyr::summarise(
    n_pathways   = dplyr::n(),
    diff_sd      = stats::sd(diff, na.rm = TRUE),        # reprioritisation dispersion
    diff_iqr     = stats::IQR(diff, na.rm = TRUE),
    diff_median  = stats::median(diff, na.rm = TRUE),    # net direction
    frac_sig     = mean(padj < 0.1, na.rm = TRUE),
    .groups = "drop") |>
  dplyr::mutate(timepoint = factor(timepoint, levels = c("6W", "12W")))

# fGSEA compartment NES per timepoint (MitoCarta category, Myc genotype rankings)
fgsea_mito_nes <- fgp |>
  dplyr::filter(category == "01_mitocarta",
                ranking %in% c("myc_6W", "myc_12W")) |>
  dplyr::mutate(timepoint = ifelse(ranking == "myc_6W", "6W", "12W")) |>
  dplyr::group_by(timepoint) |>
  dplyr::summarise(fgsea_median_NES = stats::median(NES, na.rm = TRUE),
                   fgsea_n_sig      = sum(padj_within_category < 0.05, na.rm = TRUE),
                   .groups = "drop")

abund_summary <- dplyr::left_join(abund_disp, fgsea_mito_nes, by = "timepoint")

# Plot: mitoPPS Myc-effect diff distributions per timepoint (the lens contrast)
abund_plot_df <- pw_pairs |>
  dplyr::filter(contrast %in% c("Myc_effect_6W", "Myc_effect_12W"),
                pathway %in% mito_pps_paths) |>
  dplyr::mutate(timepoint = factor(ifelse(contrast == "Myc_effect_6W", "6W", "12W"),
                                   levels = c("6W", "12W")))
p_abund <- ggplot2::ggplot(abund_plot_df,
                           ggplot2::aes(x = timepoint, y = diff, fill = timepoint)) +
  ggplot2::geom_violin(alpha = 0.5, colour = NA) +
  ggplot2::geom_jitter(width = 0.12, size = 0.8, alpha = 0.5) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::scale_fill_manual(values = c("6W" = "#E41A1C", "12W" = "#FF7F00")) +
  ggplot2::labs(
    title = "AP-abund: mitoPPS Myc-effect reprioritisation dispersion, 6W vs 12W",
    subtitle = "Wider spread at 6W (selective) -> narrower/uniform at 12W; fGSEA NES annotates absolute enrichment",
    x = NULL, y = "mitoPPS diff (Myc+ - Myc-), OXPHOS-tier pathways") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "none")
ggplot2::ggsave(file.path(out_dir, "abund_lens_dispersion.pdf"),
                p_abund, width = 7, height = 5)

# =============================================================================
# PART 4: AP-retention -- retention index on the temporal contrasts
# =============================================================================
# retention_abs = NES_pos - NES_neg (absolute, meaningful); retention_ratio =
# NES_pos / NES_neg. Per category, median NES on the two temporal rankings.
# Felsher guardrail: a Myc-signature NES that declines over time reflects the
# developmental substrate decline + a stable Myc offset, NOT loss of Myc activity.

retention <- fgp |>
  dplyr::filter(ranking %in% c("timepoint_pos", "timepoint_neg")) |>
  dplyr::group_by(category, ranking) |>
  dplyr::summarise(median_NES = stats::median(NES, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = ranking, values_from = median_NES) |>
  dplyr::rename(NES_pos = timepoint_pos, NES_neg = timepoint_neg) |>
  dplyr::mutate(
    retention_abs   = NES_pos - NES_neg,
    retention_ratio = NES_pos / NES_neg) |>
  dplyr::arrange(retention_abs)

p_ret <- retention |>
  dplyr::mutate(category = stats::reorder(category, retention_abs)) |>
  ggplot2::ggplot(ggplot2::aes(x = retention_abs, y = category,
                               fill = retention_abs > 0)) +
  ggplot2::geom_col() +
  ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
  ggplot2::scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                             guide = "none") +
  ggplot2::labs(
    title = "AP-retention: temporal NES retained under Myc (NES_pos - NES_neg)",
    subtitle = "Per category, temporal contrasts. >0 = Myc+ amplifies the WT temporal shift; <0 = attenuates",
    x = "retention_abs = median NES(Myc+ 6W->12W) - median NES(WT 6W->12W)", y = NULL) +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "retention_index_by_category.pdf"),
                p_ret, width = 8, height = 5)

# =============================================================================
# PART 5: SAVE
# =============================================================================

reframe_out <- list(
  mtdna = list(
    trajectory      = mtdna_traj,
    per_complex_cor = per_complex_cor,
    mtdna_pathway   = mtdna_nm,
    nuclear_oxphos  = oxphos_nuclear),
  abund = list(
    summary   = abund_summary,        # dispersion + fGSEA NES per timepoint
    diffs     = abund_plot_df),
  retention = retention,
  params = list(
    group_levels   = group_levels,
    mito_compartment = mito_pps_paths,
    legacy_loaded  = c(pooled = !is.null(legacy_pooled), xs = !is.null(legacy_xs))),
  notes = paste(
    "Reframe only (no model fitting). AP-mtDNA: synthetic 'mtDNA-encoded OXPHOS",
    "subunits' vs nuclear OXPHOS-tier composite, mitoPPS + raw, per group;",
    "corroborates script 19 AP1 (WT mtDNA up / nuclear down). AP-abund: fGSEA",
    "(rank-relative, absolute) vs mitoPPS (within-compartment reprioritisation)",
    "are complementary; agreement is TIME-DEPENDENT -- 6W selective (wide diff",
    "spread) -> 12W uniform (narrower). Not one word. AP-retention: NES_pos -",
    "NES_neg per category on temporal contrasts. FELSHER GUARDRAIL: declining",
    "Myc-signature NES = developmental substrate decline + stable Myc offset,",
    "NOT loss of Myc activity (dose constant, CLAUDE.md).")
)
saveRDS(reframe_out, here::here("results", "reframe_supp3.rds"))
message("Saved results/reframe_supp3.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  ref <- readRDS(here::here("results", "reframe_supp3.rds"))

  # AP-mtDNA: does the WT substrate reallocate mtDNA-up / nuclear-down 6W->12W?
  ref$mtdna$trajectory |>
    dplyr::filter(grepl("mitoPPS", component)) |>
    dplyr::arrange(component, group) |> print(n = Inf)
  # Per-complex correlations (mtDNA vs each nuclear complex), all samples + WT
  ref$mtdna$per_complex_cor |> head(15) |> print()

  # AP-abund: reprioritisation dispersion narrows 6W->12W? + fGSEA NES annotation
  ref$abund$summary |> print()

  # AP-retention: which categories does Myc amplify vs attenuate over time?
  ref$retention |> print(n = Inf)
  # Sign convention: MitoCarta/biogenesis expected retention_abs < 0 (attenuation
  # via convergence, cf. script 20); Myc-signatures carry the Felsher guardrail.

  list.files(here::here("outputs", "reframe_supp3"), pattern = "\\.pdf$")
}
