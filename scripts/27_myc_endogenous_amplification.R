# scripts/27_myc_endogenous_amplification.R
# =============================================================================
# Endogenous Myc & the pubertal TEB/proliferative phenotype -- does the Myc+
# transgene AMPLIFY the endogenous axis? (Block A revision, Issue #2)
# =============================================================================
#
# Reported (literature): endogenous Myc is instrumental in the WT for establishing
# the proliferative phenotype of the pubertal gland (esp. terminal end buds / TEBs).
# This is a quick quantification, a reframe on already-scored GSVA data (no re-run),
# showing (i) endogenous Myc marks that program in the WT and (ii) the Myc+ transgene
# amplifies the SAME axis, most strongly at 6W.
#
# Design: decompose each program's variation into
#   ENDOGENOUS = WT (Myc-) 6W->12W temporal shift (pubertal -> adult), and
#   TRANSGENE  = Myc genotype MAIN effect (Myc+ - WT, pooled over time),
# and show Myc activity / proliferation / TEB form ONE coupled axis (same in both
# genotypes = amplification, not a new axis).
#
# Programs (all GSVA-scored per sample):
#   myc         = MYC_signatures composite (17 sets: Felsher, Hallmark_V2, Dang core,
#                 Coller/Schuhmacher/Yu/... targets)   [Hallmark_V1 is NOT scored]
#   felsher     = MYC_felsher_integrative_signature    (reported individually)
#   hallmark_v2 = MYC_HALLMARK_MYC_TARGETS_V2          (reported individually)
#   prolif      = Proliferation composite (14 sets)     -- the functional phenotype
#   myc_in_teb  = TFT_MYC_GRAY_*_TEB composite          -- MYC targets in the TEB context
#   teb_ductal  = mean(MG_TEB_VS_DUCTAL_*_UP) - mean(*_DN) -- the TEB->ductal phenotype
#
# Ceiling (notes): the TRANSGENE genotype effect is POWERED (main effect, p~1e-5); the
# ENDOGENOUS WT-temporal shift is DIRECTIONAL (ns at n=6/timepoint); composite p's and
# correlations use correlated sets -> indicative. "Endogenous Myc establishes the TEB
# phenotype" is a LITERATURE claim -- our data show WT co-variation consistent with it
# and that the transgene amplifies the same axis; not a causal proof.
#
# Input:  results/gsva_scores.rds (scores + set_meta + sample_meta)
# Output: results/myc_endogenous_amplification.rds
#         outputs/myc_endogenous_amplification/*.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

out_dir <- here::here("outputs", "myc_endogenous_amplification")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# =============================================================================
# PART 1: LOAD
# =============================================================================
gsva     <- readRDS(here::here("results", "gsva_scores.rds"))
scores   <- gsva$scores
set_meta <- gsva$set_meta
sm       <- as.data.frame(gsva$sample_meta)
sm       <- sm[colnames(scores), , drop = FALSE]
stopifnot(identical(rownames(sm), colnames(scores)))
sm$timepoint  <- factor(sm$timepoint,  levels = c("6W", "12W"))
sm$myc_status <- factor(sm$myc_status, levels = c("neg", "pos"))
sm$group      <- factor(sm$group, levels = c("6W_neg", "6W_pos", "12W_neg", "12W_pos"))

# =============================================================================
# PART 0: PROGRAM PANELS (per-sample composites)
# =============================================================================
comp <- function(sets) {
  sets <- intersect(sets, rownames(scores))
  stopifnot(length(sets) > 0)
  colMeans(scores[sets, , drop = FALSE])
}
myc_sets     <- set_meta$set_name[set_meta$category_primary == "MYC_signatures"]
prol_sets    <- set_meta$set_name[set_meta$category_primary == "Proliferation"]
teb_myc_sets <- grep("^TFT_MYC_GRAY_.*TEB", rownames(scores), value = TRUE)
teb_up       <- grep("^MG_TEB_VS_DUCTAL_.*_UP$", rownames(scores), value = TRUE)
teb_dn       <- grep("^MG_TEB_VS_DUCTAL_.*_DN$", rownames(scores), value = TRUE)
stopifnot("MYC_felsher_integrative_signature" %in% rownames(scores),
          "MYC_HALLMARK_MYC_TARGETS_V2"       %in% rownames(scores),
          length(teb_myc_sets) > 0, length(teb_up) > 0, length(teb_dn) > 0)
message(sprintf("Panels: MYC_sig %d | Prolif %d | MYC-in-TEB %d | TEB UP/DN %d/%d",
                length(intersect(myc_sets, rownames(scores))),
                length(intersect(prol_sets, rownames(scores))),
                length(teb_myc_sets), length(teb_up), length(teb_dn)))

prog_wide <- tibble::tibble(
  sample     = colnames(scores),
  group      = sm$group, timepoint = sm$timepoint, myc_status = sm$myc_status,
  myc         = comp(myc_sets),
  felsher     = scores["MYC_felsher_integrative_signature", ],
  hallmark_v2 = scores["MYC_HALLMARK_MYC_TARGETS_V2", ],
  prolif      = comp(prol_sets),
  myc_in_teb  = comp(teb_myc_sets),
  teb_ductal  = comp(teb_up) - comp(teb_dn))

prog_names  <- c("myc", "felsher", "hallmark_v2", "prolif", "myc_in_teb", "teb_ductal")
prog_labels <- c(myc = "MYC signatures (17)", felsher = "Felsher",
                 hallmark_v2 = "Hallmark MYC V2", prolif = "Proliferation (14)",
                 myc_in_teb = "MYC-in-TEB (Gray)", teb_ductal = "TEB - ductal")

# =============================================================================
# PART A: ENDOGENOUS vs TRANSGENE DECOMPOSITION  (prog_stats)
# =============================================================================
# Reuses the state_stats lm pattern from script 26: WT-only temporal (endogenous),
# genotype MAIN effect (transgene, powered), interaction; effect sizes vs within-group SD.
prog_stat_one <- function(nm) {
  d   <- data.frame(y = prog_wide[[nm]], myc = prog_wide$myc_status,
                    tp = prog_wide$timepoint, grp = prog_wide$group)
  gm  <- tapply(d$y, d$grp, mean)
  wsd <- sqrt(mean(tapply(d$y, d$grp, stats::var)))
  wt  <- summary(stats::lm(y ~ tp, dplyr::filter(d, myc == "neg")))$coefficients["tp12W", ]
  mp  <- summary(stats::lm(y ~ tp, dplyr::filter(d, myc == "pos")))$coefficients["tp12W", ]
  ma  <- summary(stats::lm(y ~ myc + tp, d))$coefficients["mycpos", ]
  mi  <- summary(stats::lm(y ~ tp * myc, d))$coefficients["tp12W:mycpos", ]
  tibble::tibble(
    program = nm, within_sd = wsd,
    wt_shift = unname(wt["Estimate"]), wt_d = unname(wt["Estimate"]) / wsd,
    wt_p = unname(wt["Pr(>|t|)"]),
    mycpos_shift = unname(mp["Estimate"]), mycpos_p = unname(mp["Pr(>|t|)"]),
    geno_beta = unname(ma["Estimate"]), geno_d = unname(ma["Estimate"]) / wsd,
    geno_p = unname(ma["Pr(>|t|)"]),
    int_beta = unname(mi["Estimate"]), int_p = unname(mi["Pr(>|t|)"]),
    m6_neg = unname(gm["6W_neg"]), m6_pos = unname(gm["6W_pos"]),
    m12_neg = unname(gm["12W_neg"]), m12_pos = unname(gm["12W_pos"]))
}
prog_stats <- dplyr::bind_rows(lapply(prog_names, prog_stat_one)) |>
  dplyr::mutate(
    transgene_verdict = dplyr::case_when(
      geno_p < 0.001 ~ "powered (strong)",
      geno_p < 0.05  ~ "powered",
      TRUE           ~ "weak / ns"),
    # endogenous pubertal expectation = 6W > 12W in WT (wt_shift < 0)
    endogenous_dir = dplyr::case_when(
      wt_shift < 0 & wt_p < 0.05 ~ "6W>12W (sig)",
      wt_shift < 0               ~ "6W>12W (directional, ns)",
      wt_shift > 0               ~ "12W>6W",
      TRUE                       ~ "flat"),
    label = prog_labels[program])

# =============================================================================
# PART B: THE COUPLED AXIS  (Myc -> proliferation -> TEB)
# =============================================================================
prog_mat <- as.matrix(prog_wide[, prog_names])
cor_overall <- stats::cor(prog_mat, method = "spearman")

# key pairs, overall + within genotype + within timepoint
cor_group <- function(rows, tag) {
  m <- stats::cor(prog_mat[rows, , drop = FALSE], method = "spearman")
  tibble::tibble(
    subset = tag,
    myc_prolif      = m["myc", "prolif"],
    myc_myc_in_teb  = m["myc", "myc_in_teb"],
    myc_teb_ductal  = m["myc", "teb_ductal"],
    prolif_teb      = m["prolif", "teb_ductal"])
}
cor_key <- dplyr::bind_rows(
  cor_group(rep(TRUE, nrow(prog_wide)), "all"),
  cor_group(prog_wide$myc_status == "neg", "WT (Myc-)"),
  cor_group(prog_wide$myc_status == "pos", "Myc+"),
  cor_group(prog_wide$timepoint == "6W", "6W"),
  cor_group(prog_wide$timepoint == "12W", "12W"))

# same-axis check: slope of the phenotype on Myc activity within each genotype
axis_slopes <- dplyr::bind_rows(lapply(c("neg", "pos"), function(gt) {
  d <- dplyr::filter(prog_wide, myc_status == gt)
  tibble::tibble(genotype = gt,
    slope_teb_on_myc    = stats::coef(stats::lm(teb_ductal ~ myc, d))["myc"],
    slope_prolif_on_myc = stats::coef(stats::lm(prolif ~ myc, d))["myc"],
    mean_myc = mean(d$myc))
}))

# =============================================================================
# PART C: FIGURES
# =============================================================================
geno_cols <- c(neg = "#4575B4", pos = "#D73027")
prog_long <- prog_wide |>
  tidyr::pivot_longer(dplyr::all_of(prog_names), names_to = "program", values_to = "score") |>
  dplyr::mutate(program = factor(prog_labels[program], levels = unname(prog_labels)))

# C1: 6 programs x 4 conditions (amplification + WT 6W>12W)
gm_long <- prog_long |>
  dplyr::group_by(program, group, timepoint, myc_status) |>
  dplyr::summarise(mean = mean(score), se = stats::sd(score) / sqrt(dplyr::n()), .groups = "drop")
p_panel <- ggplot2::ggplot(gm_long, ggplot2::aes(x = timepoint, y = mean,
    colour = myc_status, group = myc_status)) +
  ggplot2::geom_line(linewidth = 0.6) + ggplot2::geom_point(size = 2) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), width = 0.12) +
  ggplot2::facet_wrap(~ program, scales = "free_y") +
  ggplot2::scale_colour_manual(values = geno_cols) +
  ggplot2::labs(title = "Endogenous Myc & the pubertal proliferative/TEB program",
    subtitle = "blue = WT (endogenous), red = Myc+ (transgene amplifies); mean +- se",
    x = NULL, y = "GSVA score") +
  ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "top")
ggplot2::ggsave(file.path(out_dir, "program_panels.pdf"), p_panel, width = 8.5, height = 5.5)

# C2: endogenous (WT) vs transgene (genotype) shift per program, with p's + d
eff_long <- prog_stats |>
  dplyr::transmute(label, `endogenous (WT 6->12)` = wt_shift, `transgene (Myc+ - WT)` = geno_beta) |>
  tidyr::pivot_longer(-label, names_to = "effect", values_to = "shift")
eff_lab <- prog_stats |>
  dplyr::transmute(label,
    lab = sprintf("transgene p=%.1e (d=%.2f)\nWT p=%.2g (%s)", geno_p, geno_d, wt_p, endogenous_dir))
p_eff <- ggplot2::ggplot(eff_long, ggplot2::aes(x = effect, y = shift, fill = effect)) +
  ggplot2::geom_col() + ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
  ggplot2::facet_wrap(~ label, scales = "free_y") +
  ggplot2::geom_text(data = eff_lab, ggplot2::aes(x = 1.5, y = Inf, label = lab),
                     inherit.aes = FALSE, vjust = 1.2, size = 2.4, lineheight = 0.9) +
  ggplot2::scale_fill_manual(values = c(`endogenous (WT 6->12)` = "#4575B4",
                                        `transgene (Myc+ - WT)` = "#D73027")) +
  ggplot2::labs(title = "Endogenous (WT temporal) vs transgene (genotype) effect",
    subtitle = "transgene = powered main effect; endogenous = directional (n=6/tp)",
    x = NULL, y = "GSVA shift") +
  ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
ggplot2::ggsave(file.path(out_dir, "endogenous_vs_transgene.pdf"), p_eff, width = 8.5, height = 5.5)

# C3: the axis -- Myc activity vs TEB and vs proliferation, coloured by group
axis_long <- prog_wide |>
  dplyr::select(myc, group, teb_ductal, prolif) |>
  tidyr::pivot_longer(c(teb_ductal, prolif), names_to = "phenotype", values_to = "y") |>
  dplyr::mutate(phenotype = dplyr::recode(phenotype, teb_ductal = "TEB - ductal",
                                          prolif = "Proliferation"))
p_axis <- ggplot2::ggplot(axis_long, ggplot2::aes(x = myc, y = y, colour = group)) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "grey50", linewidth = 0.5,
                       formula = y ~ x) +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::facet_wrap(~ phenotype, scales = "free_y") +
  ggplot2::labs(title = "One coupled axis: Myc activity -> proliferation / TEB",
    subtitle = "Myc+ groups (esp. 6W_pos) shifted further along the same axis",
    x = "MYC signatures (GSVA)", y = NULL) +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "myc_axis_scatter.pdf"), p_axis, width = 8.5, height = 4.5)

# C4: correlation heatmap of the 6 programs
cor_df <- as.data.frame(as.table(cor_overall))
names(cor_df) <- c("p1", "p2", "rho")
cor_df$p1 <- factor(prog_labels[as.character(cor_df$p1)], levels = unname(prog_labels))
cor_df$p2 <- factor(prog_labels[as.character(cor_df$p2)], levels = unname(prog_labels))
p_cor <- ggplot2::ggplot(cor_df, ggplot2::aes(x = p1, y = p2, fill = rho)) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho)), size = 3) +
  ggplot2::scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                                midpoint = 0, limits = c(-1, 1)) +
  ggplot2::labs(title = "Program coupling (per-sample Spearman)", x = NULL, y = NULL) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
ggplot2::ggsave(file.path(out_dir, "program_correlation.pdf"), p_cor, width = 6.5, height = 5.5)

# =============================================================================
# PART D: SAVE
# =============================================================================
group_means <- prog_long |>
  dplyr::group_by(program, group) |>
  dplyr::summarise(mean = mean(score), .groups = "drop") |>
  tidyr::pivot_wider(names_from = group, values_from = mean)

myc_out <- list(
  prog_wide    = prog_wide,
  prog_stats   = prog_stats,
  cor_overall  = cor_overall,
  cor_key      = cor_key,
  axis_slopes  = axis_slopes,
  group_means  = group_means,
  panels = list(myc_sets = myc_sets, prol_sets = prol_sets,
                teb_myc_sets = teb_myc_sets, teb_up = teb_up, teb_dn = teb_dn),
  notes = paste(
    "Issue #2: endogenous Myc & the pubertal TEB/proliferative phenotype. Reframe on",
    "already-scored GSVA (no re-run). Programs: myc (MYC_signatures 17-set composite),",
    "felsher, hallmark_v2 (individual), prolif (Proliferation 14), myc_in_teb",
    "(TFT_MYC_GRAY_*_TEB = MYC targets in the Gray TEB context), teb_ductal",
    "(MG_TEB_VS_DUCTAL UP-DN). Decomposition: ENDOGENOUS = WT 6W->12W temporal;",
    "TRANSGENE = Myc genotype MAIN effect. Result: transgene amplification POWERED",
    "(MYC_sig geno p~1e-5, ~+0.5 GSVA, Myc+>>WT both ages); endogenous WT is 6W>12W",
    "for Myc activity/proliferation/TEB (DIRECTIONAL, ns at n=6/tp); the three form",
    "ONE coupled axis (per-sample rho 0.58-0.86) and Myc+ is shifted to high Myc +",
    "high TEB/proliferation (6W_pos most TEB/proliferative). NB the within-Myc+",
    "TEB/prolif-on-Myc slope is SHALLOWER than WT (1.65->0.83; likely saturation at",
    "high Myc), so 'amplifies the same axis' holds BETWEEN groups, not as an identical",
    "within-group slope. Interaction ns -> the transgene boost is ~additive over time.",
    "Ceiling: transgene powered (main effect); endogenous directional; composite p's/",
    "correlations use correlated sets -> indicative; endogenous-Myc-establishes-TEB is",
    "a LITERATURE claim, our data show consistent WT co-variation + amplification, not",
    "causation. Hallmark_V1 not scored (V2 + Felsher + 15 targets used).")
)
saveRDS(myc_out, here::here("results", "myc_endogenous_amplification.rds"))
message("Saved results/myc_endogenous_amplification.rds")
message(sprintf("Transgene (MYC_sig genotype): beta=%+.3f p=%.2g | WT endogenous 12W-6W=%+.3f p=%.2g",
                prog_stats$geno_beta[prog_stats$program == "myc"],
                prog_stats$geno_p[prog_stats$program == "myc"],
                prog_stats$wt_shift[prog_stats$program == "myc"],
                prog_stats$wt_p[prog_stats$program == "myc"]))

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  mp <- readRDS(here::here("results", "myc_endogenous_amplification.rds"))

  # A. endogenous vs transgene, per program (the amplification quantification)
  mp$prog_stats |> dplyr::select(program, m6_neg, m6_pos, m12_neg, m12_pos,
    wt_shift, wt_p, endogenous_dir, geno_beta, geno_d, geno_p, transgene_verdict) |>
    print(width = Inf)

  # B. the coupled axis: correlations overall + within genotype/timepoint
  mp$cor_key |> print()
  round(mp$cor_overall, 2)
  mp$axis_slopes |> print()   # same slope in WT vs Myc+ = amplification, not new axis

  # group means (6 programs x 4 conditions)
  mp$group_means |> print(width = Inf)

  # C. outputs
  list.files(here::here("outputs", "myc_endogenous_amplification"), pattern = "\\.pdf$")
}
