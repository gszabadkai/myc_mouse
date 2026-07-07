# scripts/24_biogenesis_discrimination.R
# =============================================================================
# Biogenesis discrimination: Myc vs ER/PGC1a, the per-pathway detail
# (Block A, Day 5 -- resolving the review discussion, docs/BlockA_review_discussion.md)
# =============================================================================
#
# A REFRAME on already-scored sets (no new GSVA/fGSEA computation). Answers three
# of the author's review points:
#
#   Part A (point 1) -- the WITHIN-MITO imbalance. CORRECTION (author, 2026-07-07):
#     the imbalance is a within-mitochondrial quantity, so it is read on the
#     mitoPPS lens ALONE (relative prioritisation inside the fixed mito budget).
#     fGSEA is a SEPARATE axis -- importance vs the WHOLE transcriptome -- kept as
#     annotation, NOT merged into the imbalance. We:
#       (a) compare mitoPPS across pathways and Level2 pathway-GROUPS to find the
#           largest within-mito priority differences (the imbalance structure);
#       (b) foreground the headline mitonuclear imbalance -- mtDNA-encoded vs
#           nuclear-encoded OXPHOS, which ARE separated in the mitoPPS set (13 mt-*
#           genes in a synthetic pathway; nuclear complexes have mt-* stripped);
#       (c) ANCHOR on the WT 6W->12W temporal contrast -- the developmental
#           substrate the acute Myc-ER pulse hits -- then overlay the Myc effect
#           and the Myc+ window.
#     The bulk-vs-within-mito DISSOCIATION (fGSEA sign disagreeing with mitoPPS
#     sign) is a SEPARATE, secondary observation (the script-22 AP-abund finding),
#     reported as annotation but explicitly NOT called the imbalance. Candidate
#     ISR/UPR^mt frame (proteostasis retention) is tested against the Level2
#     groups, not asserted. Scope note: the mtDNA-vs-nuclear composite trajectory
#     is in script 22 and its developmental decomposition in script 25 Part C --
#     here Part A is the fuller within-mito MAP (all Level2 groups), not a redo.
#
#   Part B (point 2) -- WHICH biogenesis program is active at each stage/condition.
#     Category-7 discrimination sets grouped into lanes (MYC / ER-PGC1a / shared /
#     fork) x GSVA trajectory (coef_table) + fGSEA NES per contrast. Stated as "the
#     data are CONSISTENT WITH" -- in-silico set memberships, not a functional
#     attribution of which biogenesis "kills".
#
#   Part C (point 4d) -- where the external "ER/PGC1a kills" result lands.
#     Category-9 biogenesis x death intersections, grouped by biogenesis lane x
#     death modality: GSVA + fGSEA trajectory, and per-sample coupling of each
#     intersection to the death-priming composite (from script 23). Maps coupling,
#     not causation.
#
#   Part D (point 2) -- the TF drivers. Category-6 TF-target sets: the MYC TF lane
#     vs the ESRRA/NRF1/GABPA/ESR1 (PGC1a-axis) TF lanes across contrasts, annotated
#     with the Gray/CHEA mito-TF shortlist -- which TF program is active when.
#
# Input:  results/fgsea_percategory.rds        (per-pathway NES, 5 contrasts)
#         results/mitopps_fgsea_comparison.rds  (fGSEA x mitoPPS join, 4 contrasts)
#         results/interaction_fgsea_mitopps.rds (interaction-contrast join)
#         results/mitopps_scores.rds            (raw_pairwise = abundance lens;
#                                                pathway_levels; mtdna_pathway_name)
#         results/gsva_scores.rds               (per-sample GSVA matrix)
#         results/gsva_overview.rds             (coef_table trajectory betas)
#         results/death_timing_substrate.rds    (per-sample death priming pro_comp)
#         data/genesets_from_library/gray_chea_mito_tf_shortlist.csv
# Output: results/biogenesis_discrimination.rds; outputs/biogenesis_discrimination/*.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

out_dir <- here::here("outputs", "biogenesis_discrimination")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# =============================================================================
# PART 1: LOAD + contrast crosswalk
# =============================================================================

fgp      <- readRDS(here::here("results", "fgsea_percategory.rds"))$fgsea
cmp      <- readRDS(here::here("results", "mitopps_fgsea_comparison.rds"))
int_cmp  <- readRDS(here::here("results", "interaction_fgsea_mitopps.rds"))
mps      <- readRDS(here::here("results", "mitopps_scores.rds"))
gsva     <- readRDS(here::here("results", "gsva_scores.rds"))
gov      <- readRDS(here::here("results", "gsva_overview.rds"))
death    <- readRDS(here::here("results", "death_timing_substrate.rds"))

coef_tbl <- gov$coef_table
comparison <- tibble::as_tibble(cmp$comparison)   # pathway x contrast: NES + mitopps_diff
raw_pw     <- tibble::as_tibble(mps$raw_pairwise)  # pathway x contrast: abundance diff
lvl        <- tibble::as_tibble(mps$pathway_levels)
mtdna_nm   <- mps$mtdna_pathway_name

# The three objects label the same biology with different contrast strings.
# Canonicalise everything to one vocabulary so the lenses align cleanly.
#   Myc effect (pos vs neg) at each timepoint, and the temporal (6W->12W) contrast
#   within each genotype. fGSEA also carries the interaction ranking.
canon_contrast <- function(x) {
  dplyr::case_when(
    x %in% c("Myc_effect_6W",  "myc_6W")                    ~ "Myc@6W",
    x %in% c("Myc_effect_12W", "myc_12W")                   ~ "Myc@12W",
    x %in% c("Temporal_Myc_neg", "Temporal_Myc-", "timepoint_neg") ~ "WT_6W->12W",
    x %in% c("Temporal_Myc_pos", "Temporal_Myc+", "timepoint_pos") ~ "Myc+_6W->12W",
    x %in% c("interaction")                                 ~ "interaction",
    TRUE ~ NA_character_
  )
}
contrast_levels <- c("Myc@6W", "Myc@12W", "WT_6W->12W", "Myc+_6W->12W", "interaction")

message(sprintf("Loaded: comparison %d rows, raw_pairwise %d rows, coef_table %d sets",
                nrow(comparison), nrow(raw_pw), nrow(coef_tbl)))

# =============================================================================
# PART A: WITHIN-MITO IMBALANCE (review point 1)
# =============================================================================
# Imbalance = mitoPPS structure (within-mito reprioritisation). fGSEA is kept as
# a SEPARATE importance axis (annotation only). Anchor on the WT 6W->12W substrate.

# ---- the axes, kept separate ----
mitopps_diff_tbl <- comparison |>                     # IMBALANCE axis
  dplyr::transmute(pathway, contrast = canon_contrast(contrast),
                   mitopps_diff, padj_mitopps, tier1) |>
  dplyr::filter(!is.na(contrast))
fgsea_imp_tbl <- comparison |>                        # IMPORTANCE axis (annotation)
  dplyr::transmute(pathway, contrast = canon_contrast(contrast),
                   fgsea_NES = NES, padj_fgsea) |>
  dplyr::filter(!is.na(contrast))
abund_tbl <- raw_pw |>                                # within-mito level (support)
  dplyr::transmute(pathway, contrast = canon_contrast(contrast),
                   abund_diff = diff) |>
  dplyr::filter(!is.na(contrast))

# per-group mitoPPS level + the four contrast shifts, with Level2 grouping
group_levels_mp <- c("6W_neg", "6W_pos", "12W_neg", "12W_pos")
lvl2 <- lvl |>
  dplyr::transmute(pathway = Pathway, L1 = Pathway_Level1, L2 = Pathway_Level2,
                   grp = dplyr::coalesce(Pathway_Level2, Pathway_Level1)) |>
  # The synthetic mtDNA-encoded pathway co-files under the Level2 "OXPHOS subunits"
  # node with the NUCLEAR subunits -- averaging them would MASK the mitonuclear
  # imbalance at group level (the very thing we are mapping). Split it out so the
  # nuclear "OXPHOS subunits" group stays purely nuclear and mtDNA is its own group.
  dplyr::mutate(grp = ifelse(pathway == mtdna_nm, "mtDNA-encoded OXPHOS", grp))
mp_wide <- tibble::as_tibble(mps$mitopps_group_means) |>
  dplyr::filter(group %in% group_levels_mp) |>
  tidyr::pivot_wider(id_cols = pathway, names_from = group,
                     values_from = mean_mitopps) |>
  dplyr::left_join(lvl2, by = "pathway") |>
  dplyr::mutate(`WT_6W->12W`   = `12W_neg` - `6W_neg`,   # developmental substrate (PRIMARY)
                `Myc+_6W->12W` = `12W_pos` - `6W_pos`,
                `Myc@6W`       = `6W_pos`  - `6W_neg`,
                `Myc@12W`      = `12W_pos` - `12W_neg`)

# --- A1: imbalance MAP at Level2-group resolution (ranked by the WT substrate) ---
imbalance_group_map <- mp_wide |>
  dplyr::filter(!is.na(grp)) |>
  dplyr::group_by(L1, grp) |>
  dplyr::summarise(n = dplyr::n(),
                   `WT_6W->12W`   = median(`WT_6W->12W`,   na.rm = TRUE),
                   `Myc@6W`       = median(`Myc@6W`,       na.rm = TRUE),
                   `Myc@12W`      = median(`Myc@12W`,      na.rm = TRUE),
                   `Myc+_6W->12W` = median(`Myc+_6W->12W`, na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::arrange(`WT_6W->12W`)

# --- A1b: largest INDIVIDUAL-pathway priority differences in the WT substrate ---
wt_pathway_rank <- mp_wide |>
  dplyr::filter(!is.na(grp)) |>
  dplyr::select(pathway, L1, grp, `6W_neg`, `12W_neg`, `WT_6W->12W`) |>
  dplyr::arrange(`WT_6W->12W`)

# --- A2: headline mitonuclear imbalance -- mtDNA-encoded vs nuclear OXPHOS ---
nuclear_oxphos_subunits <- c("OXPHOS subunits", "CI subunits", "CII subunits",
                             "CIII subunits", "CIV subunits", "CV subunits")
mtnuc_traj <- tibble::as_tibble(mps$mitopps_group_means) |>
  dplyr::filter(pathway %in% c(mtdna_nm, nuclear_oxphos_subunits),
                group %in% group_levels_mp) |>
  tidyr::separate(group, into = c("timepoint", "myc_status"), sep = "_",
                  remove = FALSE) |>
  dplyr::mutate(encoding   = ifelse(pathway == mtdna_nm, "mtDNA-encoded",
                                    "nuclear-encoded"),
                group      = factor(group, levels = group_levels_mp),
                timepoint  = factor(timepoint, levels = c("6W", "12W")),
                myc_status = factor(myc_status, levels = c("neg", "pos")))
# imbalance index per group = mean(nuclear-OXPHOS-subunit mitoPPS) - mtDNA-encoded
# (positive = nuclear-prioritised / mtDNA-lagging = imbalanced)
mtnuc_index <- mtnuc_traj |>
  dplyr::group_by(group, timepoint, myc_status, encoding) |>
  dplyr::summarise(mp = mean(mean_mitopps), .groups = "drop") |>
  tidyr::pivot_wider(names_from = encoding, values_from = mp) |>
  dplyr::mutate(imbalance = `nuclear-encoded` - `mtDNA-encoded`)

# --- A3 (SEPARATE axis): bulk-vs-within-mito DISSOCIATION -- NOT the imbalance ---
# Where transcriptome-relative importance (fGSEA) disagrees with within-mito
# priority (mitoPPS): the script-22 AP-abund observation, reported as a distinct
# annotation. Explicitly not the imbalance.
dissociation <- mitopps_diff_tbl |>
  dplyr::left_join(fgsea_imp_tbl, by = c("pathway", "contrast")) |>
  dplyr::left_join(abund_tbl,     by = c("pathway", "contrast")) |>
  dplyr::mutate(
    dissociation = dplyr::case_when(
      fgsea_NES < 0 & mitopps_diff > 0 ~ "down_vs_transcriptome__up_within_mito",
      fgsea_NES > 0 & mitopps_diff < 0 ~ "up_vs_transcriptome__down_within_mito",
      fgsea_NES > 0 & mitopps_diff > 0 ~ "concordant_up",
      fgsea_NES < 0 & mitopps_diff < 0 ~ "concordant_down",
      TRUE ~ "mixed_or_zero"),
    contrast = factor(contrast, levels = contrast_levels))
dissociation_tally <- dissociation |>
  dplyr::filter(contrast != "interaction") |>
  dplyr::count(contrast, dissociation, name = "n_pathways")

# ---- PLOTS ----
# A1: imbalance map (Level2 group x contrast), ordered by the WT shift
grp_ord <- imbalance_group_map$grp
imb_long <- imbalance_group_map |>
  tidyr::pivot_longer(c(`WT_6W->12W`, `Myc@6W`, `Myc@12W`, `Myc+_6W->12W`),
                      names_to = "contrast", values_to = "mitopps_shift") |>
  dplyr::mutate(grp = factor(grp, levels = grp_ord),
                contrast = factor(contrast, levels = c("WT_6W->12W", "Myc@6W",
                                                       "Myc@12W", "Myc+_6W->12W")))
p_A1 <- ggplot2::ggplot(imb_long,
                        ggplot2::aes(x = contrast, y = grp, fill = mitopps_shift)) +
  ggplot2::geom_tile(colour = "grey85") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.2f", mitopps_shift)), size = 2) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                midpoint = 0) +
  ggplot2::labs(title = "Within-mito imbalance map (mitoPPS shift by pathway group)",
                subtitle = "ordered by the WT 6W->12W substrate shift; + = gains within-mito priority",
                x = NULL, y = NULL, fill = "mitoPPS\nshift") +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
ggplot2::ggsave(file.path(out_dir, "A1_imbalance_group_map.pdf"),
                p_A1, width = 7.5, height = 8)

# A2: headline mitonuclear imbalance trajectory (mtDNA-encoded vs nuclear OXPHOS)
p_A2 <- ggplot2::ggplot(mtnuc_traj,
                        ggplot2::aes(x = timepoint, y = mean_mitopps,
                                     group = pathway, colour = encoding)) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 1.6) +
  ggplot2::facet_wrap(~ myc_status, labeller = ggplot2::label_both) +
  ggplot2::scale_colour_manual(values = c("mtDNA-encoded" = "#B2182B",
                                          "nuclear-encoded" = "#2166AC")) +
  ggplot2::labs(title = "Mitonuclear imbalance (headline): mtDNA-encoded vs nuclear OXPHOS",
                subtitle = "mitoPPS priority 6W->12W; the substrate rebalances as mtDNA catches up",
                x = NULL, y = "mitoPPS (within-mito priority)", colour = "OXPHOS encoding") +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "A2_mitonuclear_imbalance.pdf"),
                p_A2, width = 8, height = 4.5)

# A3: WT substrate -- ranked group shift (what changes 6W->12W in WT)
p_A3 <- imbalance_group_map |>
  dplyr::mutate(grp = factor(grp, levels = grp_ord)) |>
  ggplot2::ggplot(ggplot2::aes(x = `WT_6W->12W`, y = grp, fill = `WT_6W->12W`)) +
  ggplot2::geom_col() +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                midpoint = 0) +
  ggplot2::labs(title = "WT developmental substrate: within-mito priority shift 6W->12W",
                subtitle = "the state the acute Myc-ER pulse acts on",
                x = "mitoPPS shift (12W_neg - 6W_neg)", y = NULL, fill = "shift") +
  ggplot2::theme_bw(base_size = 8)
ggplot2::ggsave(file.path(out_dir, "A3_wt_substrate_shift.pdf"),
                p_A3, width = 7, height = 8)

message("Part A: within-mito imbalance (mitoPPS) map + mitonuclear headline + WT substrate built")

# =============================================================================
# PART B: Myc vs ER/PGC1a biogenesis lanes (review point 2)
# =============================================================================
# WHICH biogenesis program is active when. Category-7 discrimination sets, grouped
# into lanes. GSVA trajectory (coef_table) + fGSEA NES per contrast.

cat7_lane <- function(set_name) {
  dplyr::case_when(
    set_name %in% c("MYC_MITO", "MYC_SPECIFIC_MITO", "MYC_AND_CORE_MITO",
                    "MYC_AND_DEVELOPMENTAL_MITO")                 ~ "MYC",
    set_name %in% c("ESR1_AND_CORE_MITO", "ESR1_NOT_CORE_MITO",
                    "ESRRA_MITO", "NRF1_MITO", "GABPA_MITO")      ~ "ER/PGC1a",
    set_name %in% c("CORE_MITO", "DEVELOPMENTAL_MITO", "E2F1_MITO") ~ "shared",
    grepl("^MYC_TARGETS_AND_METABRIC", set_name)                 ~ "MB_fork",
    TRUE ~ NA_character_
  )
}

# GSVA trajectory betas for Category-7 sets
cat7_gsva <- coef_tbl |>
  dplyr::filter(category_primary == "Biogenesis_discrimination") |>
  dplyr::transmute(set_name, lane = cat7_lane(set_name),
                   d6, d12, beta_int, int_p, beta_time, myc_slope) |>
  dplyr::filter(!is.na(lane))

# fGSEA NES per contrast for Category-7 sets (category 07)
cat7_fgsea <- fgp |>
  dplyr::filter(category == "07_biogenesis_discrimination") |>
  dplyr::transmute(set_name = pathway,
                   contrast = canon_contrast(ranking),
                   fgsea_NES = NES, padj = padj_within_category) |>
  dplyr::filter(!is.na(contrast)) |>
  dplyr::mutate(lane = cat7_lane(set_name)) |>
  dplyr::filter(!is.na(lane))

# Per-lane summary: GSVA Myc effect at each timepoint + interaction
cat7_lane_summary <- cat7_gsva |>
  dplyr::group_by(lane) |>
  dplyr::summarise(n_sets      = dplyr::n(),
                   med_d6      = median(d6,       na.rm = TRUE),
                   med_d12     = median(d12,      na.rm = TRUE),
                   med_beta_int= median(beta_int, na.rm = TRUE),
                   min_int_p   = min(int_p,       na.rm = TRUE),
                   .groups = "drop")

cat7_fgsea_summary <- cat7_fgsea |>
  dplyr::group_by(lane, contrast) |>
  dplyr::summarise(med_NES = median(fgsea_NES, na.rm = TRUE),
                   n_sig   = sum(padj < 0.05, na.rm = TRUE),
                   n_sets  = dplyr::n(), .groups = "drop") |>
  dplyr::mutate(contrast = factor(contrast, levels = contrast_levels))

# --- PLOT B1: Cat-7 lane fGSEA NES across contrasts ---
p_B1 <- ggplot2::ggplot(cat7_fgsea_summary,
                        ggplot2::aes(x = contrast, y = med_NES,
                                     colour = lane, group = lane)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2) +
  ggplot2::labs(title = "Biogenesis lanes: fGSEA NES across contrasts (Cat-7)",
                subtitle = "which biogenesis program is active at each stage/condition (median NES per lane)",
                x = NULL, y = "median fGSEA NES", colour = "lane") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
ggplot2::ggsave(file.path(out_dir, "B1_cat7_lane_fgsea.pdf"),
                p_B1, width = 7, height = 4.5)

# --- PLOT B2: per-set GSVA Myc effect (d6 vs d12) coloured by lane ---
p_B2 <- ggplot2::ggplot(cat7_gsva,
                        ggplot2::aes(x = d6, y = d12, colour = lane, label = set_name)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70") +
  ggplot2::geom_point(size = 2.4) +
  ggrepel::geom_text_repel(size = 2.3, max.overlaps = 20, show.legend = FALSE) +
  ggplot2::labs(title = "Cat-7 sets: Myc GSVA effect at 6W (d6) vs 12W (d12)",
                subtitle = "below the diagonal = Myc effect attenuates with age; by lane",
                x = "d6 (Myc effect @6W)", y = "d12 (Myc effect @12W)", colour = "lane") +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "B2_cat7_gsva_d6_d12.pdf"),
                p_B2, width = 7.5, height = 6)

message("Part B: Cat-7 Myc-vs-ER lanes built")

# =============================================================================
# PART C: biogenesis x death coupling (review point 4d)
# =============================================================================
# Category-9 intersections, grouped by biogenesis lane x death modality. GSVA +
# fGSEA trajectory, and per-sample coupling to the death-priming composite.

cat9_parse <- function(set_name) {
  bio <- dplyr::case_when(
    grepl("^MYC_SPECIFIC_MITO_AND_", set_name) ~ "MYC_SPECIFIC",
    grepl("^MYC_MITO_AND_",          set_name) ~ "MYC",
    grepl("^CORE_MITO_AND_",         set_name) ~ "CORE",
    grepl("^DEVELOPMENTAL_MITO_AND_",set_name) ~ "DEVELOPMENTAL",
    grepl("^MITO_NU_AND_",           set_name) ~ "MITO_NU",
    TRUE ~ "other")
  death <- sub("^.*_AND_", "", set_name)
  tibble::tibble(bio_lane = bio, death_modality = death)
}

cat9_sets <- coef_tbl |>
  dplyr::filter(category_primary == "Biogenesis_apoptosis_intersections") |>
  dplyr::pull(set_name)
cat9_map <- dplyr::bind_cols(tibble::tibble(set_name = cat9_sets),
                             purrr::map_dfr(cat9_sets, cat9_parse))

cat9_gsva <- coef_tbl |>
  dplyr::filter(category_primary == "Biogenesis_apoptosis_intersections") |>
  dplyr::transmute(set_name, d6, d12, beta_int, int_p, beta_time) |>
  dplyr::left_join(cat9_map, by = "set_name")

cat9_fgsea <- fgp |>
  dplyr::filter(category == "09_biogenesis_apoptosis_intersections") |>
  dplyr::transmute(set_name = pathway,
                   contrast = canon_contrast(ranking),
                   fgsea_NES = NES, padj = padj_within_category) |>
  dplyr::filter(!is.na(contrast)) |>
  dplyr::left_join(cat9_map, by = "set_name") |>
  dplyr::mutate(contrast = factor(contrast, levels = contrast_levels))

# Per-sample coupling to death priming. GSVA scores (sets x samples) -> align to
# the per-sample death composite (pro_comp) from script 23; correlate within
# timepoint. Map sample columns via gsva sample_meta.
gsva_scores <- gsva$scores
smeta       <- gsva$sample_meta
death_ps    <- death$h2$per_sample     # sample, group, mitonuclear_imbalance, bio_comp, pro_comp

# Build a per-sample tidy frame: one row per (sample, cat9 set) with GSVA + pro_comp
cat9_present <- intersect(cat9_map$set_name, rownames(gsva_scores))
gsva_c9 <- gsva_scores[cat9_present, , drop = FALSE]
# sample_meta rownames are the sample codes = colnames of gsva_scores
smeta_df <- tibble::as_tibble(smeta) |>
  dplyr::mutate(.col = rownames(smeta))
# death per-sample keyed by `sample`; align to smeta by `sample`
stopifnot(all(c("sample", "pro_comp", "mitonuclear_imbalance") %in% names(death_ps)))
samp_join <- smeta_df |>
  dplyr::select(.col, sample, group, timepoint) |>
  dplyr::left_join(dplyr::select(death_ps, sample, pro_comp, mitonuclear_imbalance),
                   by = "sample")
if (any(is.na(samp_join$pro_comp))) {
  warning("Part C: some samples did not match death per_sample by `sample`; ",
          "check ID alignment (n unmatched = ",
          sum(is.na(samp_join$pro_comp)), ")")
}

coupling_within_tp <- function(set_row_name) {
  v <- gsva_c9[set_row_name, samp_join$.col]
  df <- samp_join |> dplyr::mutate(gsva = as.numeric(v))
  df |>
    dplyr::group_by(timepoint) |>
    dplyr::summarise(
      r_pro       = suppressWarnings(stats::cor(gsva, pro_comp,
                                                use = "complete.obs")),
      r_imbalance = suppressWarnings(stats::cor(gsva, mitonuclear_imbalance,
                                                use = "complete.obs")),
      .groups = "drop") |>
    dplyr::mutate(set_name = set_row_name)
}
cat9_coupling <- purrr::map_dfr(cat9_present, coupling_within_tp) |>
  dplyr::left_join(cat9_map, by = "set_name") |>
  dplyr::mutate(timepoint = factor(timepoint, levels = c("6W", "12W")))

# Per-lane median coupling (does the balanced CORE/ER intersection couple to death
# differently from the unbalanced MYC_SPECIFIC one, and how does it move 6W->12W?)
cat9_coupling_summary <- cat9_coupling |>
  dplyr::group_by(bio_lane, timepoint) |>
  dplyr::summarise(n_sets       = dplyr::n(),
                   med_r_pro    = median(r_pro,       na.rm = TRUE),
                   med_r_imbal  = median(r_imbalance, na.rm = TRUE),
                   .groups = "drop")

# --- PLOT C1: Cat-9 death coupling by biogenesis lane, 6W vs 12W ---
p_C1 <- ggplot2::ggplot(cat9_coupling_summary,
                        ggplot2::aes(x = timepoint, y = med_r_pro,
                                     colour = bio_lane, group = bio_lane)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::labs(title = "Biogenesis x death coupling (Cat-9), by lane",
                subtitle = "per-sample corr of intersection GSVA with death priming (pro_comp), within timepoint",
                x = NULL, y = "median r (GSVA vs pro_comp)", colour = "biogenesis lane") +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "C1_cat9_death_coupling.pdf"),
                p_C1, width = 7, height = 4.5)

# --- PLOT C2: Cat-9 fGSEA NES across contrasts, faceted by biogenesis lane ---
cat9_fgsea_summary <- cat9_fgsea |>
  dplyr::group_by(bio_lane, contrast) |>
  dplyr::summarise(med_NES = median(fgsea_NES, na.rm = TRUE),
                   n_sets  = dplyr::n(), .groups = "drop")
p_C2 <- ggplot2::ggplot(cat9_fgsea_summary,
                        ggplot2::aes(x = contrast, y = med_NES,
                                     colour = bio_lane, group = bio_lane)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::labs(title = "Cat-9 intersection fGSEA NES across contrasts, by lane",
                x = NULL, y = "median fGSEA NES", colour = "biogenesis lane") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
ggplot2::ggsave(file.path(out_dir, "C2_cat9_fgsea_by_lane.pdf"),
                p_C2, width = 7.5, height = 4.5)

message("Part C: Cat-9 biogenesis x death coupling built")

# =============================================================================
# PART D: TF driver lanes (review point 2)
# =============================================================================
# Category-6 TF-target sets: the MYC TF lane vs the ESRRA/NRF1/GABPA/ESR1
# (PGC1a-axis) TF lanes. Annotated with the Gray/CHEA mito-TF shortlist.

tf_lane <- function(set_name) {
  s <- toupper(set_name)
  dplyr::case_when(
    grepl("MYC", s)                                  ~ "MYC",
    grepl("ESRRA|NRF1|GABPA|PPARGC|PGC1", s)          ~ "PGC1a-axis",
    grepl("ESR1|ESR2|\\bER\\b|ESTROGEN", s)           ~ "ER",
    grepl("E2F", s)                                   ~ "E2F",
    TRUE ~ "other")
}

tf_gsva <- coef_tbl |>
  dplyr::filter(category_primary == "TF_targets") |>
  dplyr::transmute(set_name, lane = tf_lane(set_name),
                   d6, d12, beta_int, int_p, beta_time)

tf_fgsea <- fgp |>
  dplyr::filter(category == "06_tf_targets") |>
  dplyr::transmute(set_name = pathway,
                   contrast = canon_contrast(ranking),
                   fgsea_NES = NES, padj = padj_within_category) |>
  dplyr::filter(!is.na(contrast)) |>
  dplyr::mutate(lane = tf_lane(set_name),
                contrast = factor(contrast, levels = contrast_levels))

# Focus lanes of interest for the Myc-vs-PGC1a-axis attribution
tf_focus_lanes <- c("MYC", "PGC1a-axis", "ER", "E2F")
tf_fgsea_summary <- tf_fgsea |>
  dplyr::filter(lane %in% tf_focus_lanes) |>
  dplyr::group_by(lane, contrast) |>
  dplyr::summarise(med_NES = median(fgsea_NES, na.rm = TRUE),
                   n_sig   = sum(padj < 0.05, na.rm = TRUE),
                   n_sets  = dplyr::n(), .groups = "drop")

# Gray/CHEA mito-TF shortlist -- which TFs are the mito-biogenesis enriched drivers
# (context annotation; loaded defensively).
gray_chea <- tryCatch(
  readr::read_csv(here::here("data", "genesets_from_library",
                             "gray_chea_mito_tf_shortlist.csv"),
                  show_col_types = FALSE),
  error = function(e) NULL)

p_D1 <- ggplot2::ggplot(tf_fgsea_summary,
                        ggplot2::aes(x = contrast, y = med_NES,
                                     colour = lane, group = lane)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::labs(title = "TF driver lanes: fGSEA NES across contrasts (Cat-6)",
                subtitle = "which TF program (MYC vs PGC1a-axis vs ER vs E2F) is active when",
                x = NULL, y = "median fGSEA NES", colour = "TF lane") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
ggplot2::ggsave(file.path(out_dir, "D1_tf_lane_fgsea.pdf"),
                p_D1, width = 7, height = 4.5)

message("Part D: TF driver lanes built")

# =============================================================================
# SAVE
# =============================================================================

biog_out <- list(
  # Part A -- within-mito imbalance (mitoPPS); fGSEA kept SEPARATE as annotation
  imbalance_group_map = imbalance_group_map,
  wt_pathway_rank     = wt_pathway_rank,
  mtnuc_traj          = mtnuc_traj,
  mtnuc_index         = mtnuc_index,
  fgsea_importance    = fgsea_imp_tbl,
  dissociation        = dissociation,
  dissociation_tally  = dissociation_tally,
  # Part B
  cat7_gsva         = cat7_gsva,
  cat7_fgsea        = cat7_fgsea,
  cat7_lane_summary = cat7_lane_summary,
  cat7_fgsea_summary= cat7_fgsea_summary,
  # Part C
  cat9_map          = cat9_map,
  cat9_gsva         = cat9_gsva,
  cat9_fgsea        = cat9_fgsea,
  cat9_coupling     = cat9_coupling,
  cat9_coupling_summary = cat9_coupling_summary,
  cat9_fgsea_summary= cat9_fgsea_summary,
  # Part D
  tf_gsva           = tf_gsva,
  tf_fgsea          = tf_fgsea,
  tf_fgsea_summary  = tf_fgsea_summary,
  gray_chea         = gray_chea,
  contrast_levels   = contrast_levels,
  notes = paste(
    "Reframe on already-scored sets (no new GSVA/fGSEA). Part A: WITHIN-MITO",
    "IMBALANCE read on mitoPPS ALONE (fGSEA is a SEPARATE importance axis, kept",
    "as annotation -- author correction 2026-07-07). imbalance_group_map = median",
    "mitoPPS shift per Level2 pathway-group across the four contrasts, ordered by",
    "the WT 6W->12W substrate; mtnuc_traj/mtnuc_index = headline mtDNA-encoded vs",
    "nuclear-OXPHOS imbalance (the largest within-mito shift: WT mtDNA +0.52 vs",
    "nuclear complexes -0.13/-0.18 across 6W->12W -> the substrate rebalances as",
    "mtDNA catches up). dissociation = the bulk-vs-within-mito divergence",
    "(script-22 AP-abund), reported as annotation, explicitly NOT the imbalance.",
    "Scope: mtDNA/nuclear composite is in script 22 and its decomposition in 25",
    "Part C; here Part A is the fuller all-group MAP. Part B: Cat-7 lanes (MYC /",
    "ER-PGC1a / shared / MB_fork) x GSVA",
    "trajectory + fGSEA NES = which biogenesis active when (CONSISTENT WITH, not",
    "causal). Part C: Cat-9 biogenesis x death intersections, per-sample coupling",
    "of intersection GSVA to death priming (pro_comp, script 23) within timepoint",
    "-- maps coupling not causation; places the external ER/PGC1a-kills result.",
    "Part D: Cat-6 TF lanes (MYC vs PGC1a-axis vs ER vs E2F) + Gray/CHEA shortlist.",
    "CEILING: n=6 directional, inter-gene correlation, survivor bias -> substrate",
    "association not causation."),
  analysis_date = Sys.Date()
)
saveRDS(biog_out, here::here("results", "biogenesis_discrimination.rds"))
message("Saved results/biogenesis_discrimination.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  bd <- readRDS(here::here("results", "biogenesis_discrimination.rds"))

  # --- Part A: the WITHIN-MITO imbalance (mitoPPS), anchored on the WT substrate ---
  # A1: imbalance map -- which Level2 groups gain/lose within-mito priority, ranked
  #     by the WT 6W->12W substrate (the state the acute Myc-ER pulse hits).
  bd$imbalance_group_map |> print(n = Inf)

  # A2: the headline mitonuclear imbalance -- nuclear-OXPHOS minus mtDNA-encoded
  #     per group (positive = nuclear-prioritised / mtDNA-lagging = imbalanced).
  #     Expect the imbalance to peak young (6W) and rebalance by 12W as mtDNA rises.
  bd$mtnuc_index |> dplyr::arrange(group) |> print()

  # A1b: largest INDIVIDUAL-pathway priority shifts in the WT substrate
  bd$wt_pathway_rank |> head(12) |> print()      # most DE-prioritised with age
  bd$wt_pathway_rank |> dplyr::arrange(dplyr::desc(`WT_6W->12W`)) |>
    head(12) |> print()                          # most UP-prioritised with age
  # Proteostasis check for the ISR/UPR^mt candidate: are Chaperones / Proteases /
  # Protein homeostasis retained while OXPHOS complexes fall? (read A1 rows)

  # A3 (SEPARATE, secondary): the bulk-vs-within-mito DISSOCIATION -- NOT imbalance
  bd$dissociation_tally |> print(n = Inf)

  # --- Part B: which biogenesis lane is active when? ---
  bd$cat7_lane_summary |> print()
  bd$cat7_fgsea_summary |> dplyr::arrange(lane, contrast) |> print(n = Inf)
  # MYC lane vs ER/PGC1a lane: direction + timing (6W vs 12W, WT vs Myc window)

  # --- Part C: does the balanced (CORE/ER) intersection couple to death
  #     differently from the unbalanced MYC_SPECIFIC one, and how 6W->12W? ---
  bd$cat9_coupling_summary |> dplyr::arrange(bio_lane, timepoint) |> print(n = Inf)
  bd$cat9_fgsea_summary |> dplyr::arrange(bio_lane, contrast) |> print(n = Inf)

  # --- Part D: MYC TF lane vs PGC1a-axis TF lane across contrasts ---
  bd$tf_fgsea_summary |> dplyr::arrange(lane, contrast) |> print(n = Inf)
  if (!is.null(bd$gray_chea)) bd$gray_chea |> head(10) |> print()

  list.files(here::here("outputs", "biogenesis_discrimination"), pattern = "\\.pdf$")
}
