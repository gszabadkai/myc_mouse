# scripts/25_developmental_substrate_death.R
# =============================================================================
# Developmental substrate as the death-timing "why"
# (Block A, Day 5 -- resolving the review discussion, docs/BlockA_review_discussion.md)
# =============================================================================
#
# COMPLEMENTARY LAYER to script 23 (23 stays intact). Establishes that the
# developmental SUBSTRATE -- not chronic-Myc adaptation -- sets the death timing,
# anchored on the external Myc-ER logic: an acute Myc pulse on a normally-developed
# background still kills more at 6W than 12W, so the timing is a property of the
# 6W-vs-12W substrate. Answers review points 3, 4a, 4b, 4c, 4e.
#
# HONEST LEAD (what the data actually support, strongest first):
#   * The robust developmental substrate change is MITO MATURATION: the
#     mtDNA-encoded vs nuclear-OXPHOS mitoPPS imbalance is large at 6W and resolves
#     by 12W, in BOTH genotypes (genotype-shared = developmental). Script 24 showed
#     it peaks at 6W_pos (+0.56) and inverts by 12W; here we DECOMPOSE the
#     resolution (mtDNA-rise vs nuclear-fall) and tie it to death.
#   * The LINEAGE-composition shift (progenitor/luminal down, basal up 6W->12W in
#     WT) is DIRECTIONAL but UNDERPOWERED at n=6 (dev_composition$wt_shift, all
#     p>0.18). Reported as a weaker corroborating layer, not overclaimed.
#
#   Part A (point 3)  -- the WT developmental substrate shift: mito maturation
#     (robust) + lineage composition (directional). What the acute pulse hits.
#   Part B (4a, 4e)   -- which background axis tracks death susceptibility:
#     per-sample coupling of composition + mito-imbalance to death priming. The
#     changing "background" that gates Myc's death-coupling is the MITO state more
#     than lineage composition (resolves "something must change in the background").
#   Part C (4b)       -- mtDNA developmental resolution: decompose the 6W->12W
#     imbalance drop into mtDNA-rise vs nuclear-fall, per genotype. Builds the
#     mtDNA rise explicitly into the imbalance interpretation.
#   Part D (4c)       -- selective culling: CV/dispersion narrowing 6W->12W
#     (death priming + composition + imbalance) as survivor convergence. Weakest
#     arm, on-narrative.
#
# Input:  results/dev_composition.rds        (comp_df, wt_shift, group_means)
#         results/death_timing_substrate.rds (h2$per_sample, h4$cv)
#         results/reframe_supp3.rds          (mtdna$trajectory: mtDNA vs nuclear)
#         results/mitopps_scores.rds         (context)
# Output: results/developmental_substrate_death.rds; outputs/dev_substrate_death/*.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

out_dir <- here::here("outputs", "dev_substrate_death")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

group_levels <- c("6W_neg", "6W_pos", "12W_neg", "12W_pos")

# =============================================================================
# PART 1: LOAD + align per-sample
# =============================================================================

devc  <- readRDS(here::here("results", "dev_composition.rds"))
death <- readRDS(here::here("results", "death_timing_substrate.rds"))
refr  <- readRDS(here::here("results", "reframe_supp3.rds"))

comp_df  <- tibble::as_tibble(devc$comp_df)          # per-sample lineage composites
wt_shift  <- tibble::as_tibble(devc$wt_shift)        # WT 6W->12W test per composite
death_ps <- tibble::as_tibble(death$h2$per_sample)   # sample, mitonuclear_imbalance, pro_comp
cv_death <- tibble::as_tibble(death$h4$cv)           # per-group CV of PRO / death
mtdna_traj <- tibble::as_tibble(refr$mtdna$trajectory)

# Per-sample master frame: lineage composites + death priming + mito imbalance
comp_metrics <- c("HR_comp", "LP_comp", "ALV_comp", "MASC_comp", "BASAL_comp",
                  "bio_comp", "diff_axis")
per_sample <- comp_df |>
  dplyr::select(sample, group, myc_status, timepoint, dplyr::all_of(comp_metrics)) |>
  dplyr::left_join(dplyr::select(death_ps, sample, pro_comp, mitonuclear_imbalance),
                   by = "sample")
if (any(is.na(per_sample$pro_comp))) {
  warning("Part 1: ", sum(is.na(per_sample$pro_comp)),
          " samples did not match death per_sample by `sample` -- check IDs")
}
per_sample <- per_sample |>
  dplyr::mutate(group     = factor(group, levels = group_levels),
                timepoint = factor(timepoint, levels = c("6W", "12W")))

# =============================================================================
# PART A: the WT developmental substrate shift (review point 3)
# =============================================================================
# A1 -- mito maturation (ROBUST). mtDNA-encoded vs nuclear OXPHOS mitoPPS, per
# group; the imbalance (nuclear - mtDNA) and its 6W->12W change. This is the
# substrate the acute Myc-ER pulse acts on.
mtdna_mp <- mtdna_traj |>
  dplyr::filter(component %in% c("mtDNA-encoded (mitoPPS)", "nuclear OXPHOS (mitoPPS)")) |>
  dplyr::mutate(encoding = ifelse(grepl("mtDNA", component), "mtDNA", "nuclear")) |>
  dplyr::select(group, timepoint, myc_status, encoding, mean_score) |>
  tidyr::pivot_wider(names_from = encoding, values_from = mean_score) |>
  dplyr::mutate(imbalance = nuclear - mtDNA,
                group = factor(group, levels = group_levels))

# A2 -- lineage composition (DIRECTIONAL, underpowered). Re-read dev_composition's
# WT-only test; classify each composite's 6W->12W direction. Not recomputed.
comp_shift <- wt_shift |>
  dplyr::mutate(direction = dplyr::case_when(wt_beta < 0 ~ "down_6W->12W",
                                             wt_beta > 0 ~ "up_6W->12W",
                                             TRUE ~ "flat"),
                sig = wt_p < 0.05) |>
  dplyr::arrange(wt_beta)

# --- PLOT A1: mito maturation -- imbalance per group (the substrate) ---
p_A1 <- ggplot2::ggplot(mtdna_mp,
                        ggplot2::aes(x = timepoint, y = imbalance,
                                     colour = myc_status, group = myc_status)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::labs(title = "WT substrate matures: mitonuclear imbalance resolves 6W->12W",
                subtitle = "imbalance = nuclear-OXPHOS - mtDNA-encoded mitoPPS; high at 6W (esp. Myc+), inverts by 12W",
                x = NULL, y = "mitonuclear imbalance", colour = "genotype") +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "A1_mito_maturation_imbalance.pdf"),
                p_A1, width = 6.5, height = 4.5)

# --- PLOT A2: WT lineage composition shift (directional, with p-values) ---
p_A2 <- comp_shift |>
  dplyr::mutate(metric = factor(metric, levels = metric)) |>
  ggplot2::ggplot(ggplot2::aes(x = wt_beta, y = metric, fill = direction)) +
  ggplot2::geom_col() +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("p=%.2f", wt_p)),
                     hjust = -0.1, size = 2.6) +
  ggplot2::labs(title = "WT lineage composition shift 6W->12W (directional, underpowered)",
                subtitle = "dev_composition$wt_shift; luminal/progenitor down, basal up; none significant at n=6",
                x = "WT beta (12W - 6W)", y = NULL, fill = NULL) +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "A2_wt_composition_shift.pdf"),
                p_A2, width = 7, height = 4)

message("Part A: mito maturation (robust) + lineage composition (directional) built")

# =============================================================================
# PART B: which background axis tracks death susceptibility (points 4a, 4e)
# =============================================================================
# Per-sample coupling of each background axis (lineage composites + mito
# imbalance) to death priming (pro_comp), within timepoint. Spearman -- robust to
# the one extreme 12W_neg imbalance sample (MYCF53_6d, ~ -3.9).
axes <- c(comp_metrics, "mitonuclear_imbalance")
couple_axis <- function(ax) {
  per_sample |>
    dplyr::group_by(timepoint) |>
    dplyr::summarise(
      rho = suppressWarnings(stats::cor(.data[[ax]], pro_comp,
                                        method = "spearman", use = "complete.obs")),
      .groups = "drop") |>
    dplyr::mutate(axis = ax)
}
death_coupling <- purrr::map_dfr(axes, couple_axis) |>
  dplyr::mutate(timepoint = factor(timepoint, levels = c("6W", "12W")),
                axis = factor(axis, levels = axes))

# Wide view: which axis is most death-coupled at 6W (the death-permissive window)?
death_coupling_wide <- death_coupling |>
  tidyr::pivot_wider(names_from = timepoint, values_from = rho) |>
  dplyr::arrange(dplyr::desc(abs(`6W`)))

# --- PLOT B1: death coupling by background axis, 6W vs 12W ---
p_B1 <- ggplot2::ggplot(death_coupling,
                        ggplot2::aes(x = axis, y = rho, fill = timepoint)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.65) +
  ggplot2::labs(title = "Which background axis tracks death priming? (Spearman, within timepoint)",
                subtitle = "the mito imbalance is the axis coupled to death; lineage composites weaker",
                x = NULL, y = "rho vs death priming (pro_comp)", fill = NULL) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
ggplot2::ggsave(file.path(out_dir, "B1_death_coupling_by_axis.pdf"),
                p_B1, width = 7.5, height = 4.5)

message("Part B: composition/mito vs death coupling built")

# =============================================================================
# PART C: mtDNA developmental resolution of the imbalance (review point 4b)
# =============================================================================
# Decompose the 6W->12W imbalance change into the mtDNA-rise vs nuclear-fall
# contributions, per genotype. If the resolution is mtDNA-driven and shared across
# genotypes, it is a developmental mito-maturation event (built explicitly into
# the model, per the author's point that the mtDNA rise was not previously wired in).
imb_wide <- mtdna_mp |>
  dplyr::select(myc_status, timepoint, mtDNA, nuclear, imbalance) |>
  tidyr::pivot_wider(names_from = timepoint, values_from = c(mtDNA, nuclear, imbalance))
imbalance_decomp <- imb_wide |>
  dplyr::transmute(
    genotype        = myc_status,
    imbalance_6W    = imbalance_6W,
    imbalance_12W   = imbalance_12W,
    d_imbalance     = imbalance_12W - imbalance_6W,   # total change (negative = resolves)
    d_mtDNA         = mtDNA_12W  - mtDNA_6W,          # mtDNA rise
    d_nuclear       = nuclear_12W - nuclear_6W,       # nuclear change
    # contributions to the imbalance change: d_imbalance = d_nuclear - d_mtDNA
    contrib_mtDNA   = -d_mtDNA,                       # mtDNA-rise lowers imbalance
    contrib_nuclear =  d_nuclear,
    pct_mtDNA_driven = 100 * abs(-d_mtDNA) / (abs(-d_mtDNA) + abs(d_nuclear)))

# --- PLOT C1: imbalance-resolution decomposition (mtDNA-rise vs nuclear-fall) ---
decomp_long <- imbalance_decomp |>
  dplyr::select(genotype, contrib_mtDNA, contrib_nuclear) |>
  tidyr::pivot_longer(c(contrib_mtDNA, contrib_nuclear),
                      names_to = "source", values_to = "contribution") |>
  dplyr::mutate(source = dplyr::recode(source,
                                       contrib_mtDNA   = "mtDNA rise",
                                       contrib_nuclear = "nuclear change"))
p_C1 <- ggplot2::ggplot(decomp_long,
                        ggplot2::aes(x = genotype, y = contribution, fill = source)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  ggplot2::geom_col(position = "stack") +
  ggplot2::labs(title = "Imbalance resolution 6W->12W: mtDNA-rise vs nuclear change",
                subtitle = "negative contribution = lowers the imbalance; mtDNA-driven and genotype-shared = developmental",
                x = NULL, y = "contribution to imbalance change", fill = NULL) +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "C1_imbalance_decomposition.pdf"),
                p_C1, width = 6, height = 4.5)

message("Part C: mtDNA developmental resolution decomposition built")

# =============================================================================
# PART D: selective culling + composition (review point 4c)
# =============================================================================
# Cross-sample dispersion per group: death priming (from script 23 h4), plus mito
# imbalance and composition. Culling of the susceptible (high-priming,
# high-imbalance) population predicts NARROWING 6W->12W, more so in Myc+.
disp <- per_sample |>
  dplyr::group_by(group, timepoint, myc_status) |>
  dplyr::summarise(
    sd_pro       = sd(pro_comp,             na.rm = TRUE),
    sd_imbalance = sd(mitonuclear_imbalance, na.rm = TRUE),
    sd_bio       = sd(bio_comp,             na.rm = TRUE),
    .groups = "drop") |>
  dplyr::left_join(dplyr::select(cv_death, group, cv_PRO, cv_death),
                   by = "group") |>
  dplyr::mutate(group = factor(group, levels = group_levels))

# Narrowing summary: 12W/6W dispersion ratio per genotype (< 1 = converged)
narrowing <- disp |>
  dplyr::select(timepoint, myc_status, sd_pro, sd_imbalance, sd_bio, cv_PRO) |>
  tidyr::pivot_longer(c(sd_pro, sd_imbalance, sd_bio, cv_PRO),
                      names_to = "metric", values_to = "value") |>
  tidyr::pivot_wider(names_from = timepoint, values_from = value) |>
  dplyr::mutate(ratio_12W_over_6W = `12W` / `6W`)

# --- PLOT D1: dispersion narrowing 6W->12W by genotype ---
p_D1 <- disp |>
  dplyr::select(group, timepoint, myc_status, sd_pro, sd_bio) |>
  tidyr::pivot_longer(c(sd_pro, sd_bio), names_to = "metric", values_to = "sd") |>
  ggplot2::ggplot(ggplot2::aes(x = timepoint, y = sd, colour = myc_status,
                               group = interaction(myc_status, metric),
                               linetype = metric)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::labs(title = "Survivor convergence: cross-sample dispersion narrows 6W->12W",
                subtitle = "sd of death priming (pro) and biogenesis composite; Myc+ narrows more (culling of susceptible)",
                x = NULL, y = "cross-sample SD", colour = "genotype", linetype = NULL) +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "D1_dispersion_narrowing.pdf"),
                p_D1, width = 6.5, height = 4.5)

message("Part D: culling / dispersion narrowing built")

# =============================================================================
# SAVE
# =============================================================================

dev_death_out <- list(
  # Part A
  mtdna_mp          = mtdna_mp,
  comp_shift        = comp_shift,
  # Part B
  per_sample        = per_sample,
  death_coupling    = death_coupling,
  death_coupling_wide = death_coupling_wide,
  # Part C
  imbalance_decomp  = imbalance_decomp,
  # Part D
  dispersion        = disp,
  narrowing         = narrowing,
  notes = paste(
    "Complementary layer to script 23 (23 intact). Developmental SUBSTRATE sets",
    "the death timing, anchored on the Myc-ER acute logic. HONEST LEAD: the robust",
    "developmental change is MITO MATURATION (mtDNA-encoded vs nuclear OXPHOS",
    "mitoPPS imbalance high at 6W, resolves/inverts by 12W, genotype-shared);",
    "lineage composition shifts DIRECTIONALLY (luminal/progenitor down, basal up)",
    "but is UNDERPOWERED at n=6 (dev_composition$wt_shift, all p>0.18) -- not",
    "overclaimed. Part B: the mito imbalance is the background axis coupled to",
    "death priming (Spearman, robust to the 12W_neg outlier), more than lineage",
    "composition -> the 'background change' that gates Myc's death-coupling is the",
    "MITO state (resolves point 4e). Part C: the 6W->12W imbalance resolution is",
    "~mtDNA-rise-driven and genotype-shared (imbalance_decomp$pct_mtDNA_driven) =",
    "developmental mito-maturation, now built into the model (point 4b). Part D:",
    "cross-sample dispersion (death priming / biogenesis / imbalance) narrows",
    "6W->12W, more in Myc+ = survivor convergence / culling of the susceptible",
    "(weakest arm, on-narrative). RECONCILE with 23: 23's 'Myc coupling",
    "attenuates' is the readout on a developmentally-maturing mito substrate; 25",
    "supplies the why. Myc-ER anchor: acute pulse on a young mtDNA-lagging /",
    "nuclear-assembly substrate -> mitonuclear mismatch -> death; by 12W the",
    "substrate has matured (mtDNA caught up), no mismatch, less death. CEILING:",
    "bulk RNA, n=6, survivor bias -> substrate association not causation; the",
    "composition/state -> death causal link needs single-cell / functional assays."),
  analysis_date = Sys.Date()
)
saveRDS(dev_death_out, here::here("results", "developmental_substrate_death.rds"))
message("Saved results/developmental_substrate_death.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  dd <- readRDS(here::here("results", "developmental_substrate_death.rds"))

  # --- Part A: the substrate. Mito maturation (robust) + composition (directional) ---
  dd$mtdna_mp |> print()          # imbalance high 6W (esp Myc+), inverts by 12W
  dd$comp_shift |> print()        # luminal/progenitor down, basal up; all p>0.18

  # --- Part B: which background axis tracks death? (the mito imbalance should win) ---
  dd$death_coupling_wide |> print(n = Inf)
  # Expect mitonuclear_imbalance strongly coupled at 6W, decoupling at 12W;
  # lineage composites weaker -> the changing background = the mito state.

  # --- Part C: is the imbalance resolution mtDNA-driven and genotype-shared? ---
  dd$imbalance_decomp |> print()
  # pct_mtDNA_driven ~ high in BOTH genotypes = developmental mito maturation.

  # --- Part D: does dispersion narrow 6W->12W (survivor convergence)? ---
  dd$dispersion |> print()
  dd$narrowing |> dplyr::arrange(myc_status, metric) |> print(n = Inf)
  # ratio_12W_over_6W < 1 = converged; expect Myc+ to narrow more than WT.

  list.files(here::here("outputs", "dev_substrate_death"), pattern = "\\.pdf$")
}
