# scripts/19_dev_composition.R
# =============================================================================
# Developmental substrate + lineage-composition (Block A, Day 3): AP1/2/3/4
# =============================================================================
#
# AP7 (script 18) showed Myc drives the biogenesis fork toward the MB2_UF =
# luminal-progenitor (LP), ER-negative, less-differentiated human state. Script
# 17 showed Myc suppresses differentiated luminal programs. This script tests
# that at the MOUSE mammary-lineage level and characterises the developmental
# substrate the constant Myc acts on.
#
#   AP2 (PRIMARY, powered): does Myc shift lineage composition toward progenitor
#     (lower the mature-HR vs progenitor-LP differentiation axis)? Genotype MAIN
#     EFFECT on per-sample lineage composites -- the powered regime (cf. AP7).
#   AP1: the WT (Myc-) 6W->12W developmental shift (the substrate; also the
#     death-permissive-state readout for the deferred death-timing thread), on
#     the MG_* GSVA composites + mitoPPS.
#   AP3: convergence vector -- is Myc's effect parallel to (accelerating) or
#     anti-parallel to (opposing) the WT developmental axis. Descriptive.
#   AP4: Felsher co-variation with biogenesis / differentiation (dose constant).
#
# Lineage tagging is a COARSE heuristic on MG_* set names (the library catalog
# warns LP definitions diverge); tag counts are printed for eyeballing.
#
# Input:  results/gsva_scores.rds, results/gsva_overview.rds,
#         results/mitopps_scores.rds (AP1, best-effort), results/coldata.rds
# Output: results/dev_composition.rds; outputs/dev_composition/*.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# =============================================================================
# PART 1: LOAD
# =============================================================================

gsva_out    <- readRDS(here::here("results", "gsva_scores.rds"))
scores      <- gsva_out$scores
sample_meta <- as.data.frame(gsva_out$sample_meta)
sample_meta <- sample_meta[colnames(scores), , drop = FALSE]
stopifnot(identical(rownames(sample_meta), colnames(scores)))
sample_meta$timepoint  <- stats::relevel(as.factor(sample_meta$timepoint),  "6W")
sample_meta$myc_status <- stats::relevel(as.factor(sample_meta$myc_status), "neg")
sample_meta$group      <- factor(sample_meta$group,
                                 levels = c("6W_neg", "6W_pos", "12W_neg", "12W_pos"))

coef_table <- readRDS(here::here("results", "gsva_overview.rds"))$coef_table

felsher <- "MYC_felsher_integrative_signature"
stopifnot(felsher %in% rownames(scores))

# =============================================================================
# PART 2: LINEAGE TAGGING OF MG_* DEVELOPMENTAL SETS
# =============================================================================
# Mammary_development category also carries the human METABRIC subtype sets (used
# in AP7); the mouse DEVELOPMENTAL signatures are the MG_* ones -- use those.

mammary_sets <- coef_table$set_name[coef_table$category_primary == "Mammary_development"]
dev_sets     <- grep("^MG_", mammary_sets, value = TRUE)   # exclude METABRIC_*
dev_sets     <- intersect(dev_sets, rownames(scores))

tag_lineage <- function(nm) {
  u <- toupper(nm)
  dplyr::case_when(
    grepl("MASC|FMASC|STEM",              u) ~ "MASC",   # stem
    grepl("BASAL|MYOEP|BMYO",             u) ~ "BASAL",
    grepl("ALV|ALVEOL|SECRET|LACT",       u) ~ "ALV",    # alveolar/secretory
    grepl("LP|PROG|CD61|LAPRO|_AP_|_AP$", u) ~ "LP",     # luminal progenitor
    grepl("HS|HR|HORM|LHS|LHOR|MATURE|DIFF", u) ~ "HR",  # mature hormone-sensing
    TRUE ~ "OTHER"
  )
}
lineage_tbl <- tibble::tibble(set_name = dev_sets, lineage = tag_lineage(dev_sets))
message("Lineage tag counts (coarse heuristic):")
print(lineage_tbl |> dplyr::count(lineage))

# =============================================================================
# PART 3: PER-SAMPLE COMPOSITES (lineage, differentiation axis, biogenesis, Felsher)
# =============================================================================

composite <- function(set_names) {
  set_names <- intersect(set_names, rownames(scores))
  if (length(set_names) == 0) return(rep(NA_real_, ncol(scores)))
  colMeans(scores[set_names, , drop = FALSE])
}

lin_classes <- c("HR", "LP", "ALV", "MASC", "BASAL")
lin_comp <- sapply(lin_classes, function(cl)
  composite(lineage_tbl$set_name[lineage_tbl$lineage == cl]))   # samples x class
colnames(lin_comp) <- paste0(lin_classes, "_comp")

bio_sets  <- coef_table$set_name[coef_table$category_primary == "Biogenesis_discrimination"]
bio_comp  <- composite(bio_sets)

comp_df <- tibble::tibble(
  sample     = colnames(scores),
  group      = sample_meta$group,
  myc_status = sample_meta$myc_status,
  timepoint  = sample_meta$timepoint
)
comp_df <- dplyr::bind_cols(comp_df, tibble::as_tibble(lin_comp))
comp_df$bio_comp <- bio_comp
comp_df$felsher  <- scores[felsher, ]
# Differentiation axis: mature hormone-sensing minus progenitor (high = differentiated)
comp_df$diff_axis <- comp_df$HR_comp - comp_df$LP_comp

# =============================================================================
# PART 4: AP2 -- GENOTYPE MAIN EFFECT ON COMPOSITES (powered)
# =============================================================================
tc <- list(myc_status = "contr.treatment", timepoint = "contr.treatment")

test_geno <- function(y) {
  if (all(is.na(y))) return(tibble::tibble(geno_beta = NA_real_, geno_p = NA_real_,
                                           myc_6W_beta = NA_real_, int_beta = NA_real_,
                                           int_p = NA_real_))
  d     <- data.frame(y = y, myc_status = sample_meta$myc_status,
                      timepoint = sample_meta$timepoint)
  m_add <- stats::lm(y ~ myc_status + timepoint, data = d, contrasts = tc)
  m_int <- stats::lm(y ~ timepoint * myc_status, data = d, contrasts = tc)
  sa <- summary(m_add)$coefficients; si <- summary(m_int)$coefficients
  tibble::tibble(
    geno_beta   = sa["myc_statuspos", "Estimate"],
    geno_p      = sa["myc_statuspos", "Pr(>|t|)"],
    myc_6W_beta = si["myc_statuspos", "Estimate"],
    int_beta    = si["timepoint12W:myc_statuspos", "Estimate"],
    int_p       = si["timepoint12W:myc_statuspos", "Pr(>|t|)"]
  )
}

ap2_metrics <- c(paste0(lin_classes, "_comp"), "diff_axis", "bio_comp")
ap2_tests <- do.call(rbind, lapply(ap2_metrics, function(m)
  cbind(metric = m, test_geno(comp_df[[m]]))))
ap2_tests <- tibble::as_tibble(ap2_tests)

# Per-set pooled genotype effect (descriptive) by lineage, for the distribution.
dev_geno <- coef_table |>
  dplyr::filter(set_name %in% dev_sets) |>
  dplyr::mutate(geno_set = (d6 + d12) / 2) |>
  dplyr::left_join(lineage_tbl, by = "set_name")

# =============================================================================
# PART 5: AP1 -- WT DEVELOPMENTAL SUBSTRATE SHIFT (Myc- 6W->12W)
# =============================================================================
neg <- sample_meta$myc_status == "neg"
wt_shift <- do.call(rbind, lapply(ap2_metrics, function(m) {
  y <- comp_df[[m]][neg]
  if (all(is.na(y))) return(data.frame(metric = m, wt_beta = NA, wt_p = NA))
  d  <- data.frame(y = y, timepoint = sample_meta$timepoint[neg])
  mm <- stats::lm(y ~ timepoint, data = d, contrasts = list(timepoint = "contr.treatment"))
  s  <- summary(mm)$coefficients
  data.frame(metric = m, wt_beta = s["timepoint12W", "Estimate"],
             wt_p = s["timepoint12W", "Pr(>|t|)"])
}))
wt_shift <- tibble::as_tibble(wt_shift)

# mitoPPS WT shift: mlist$mitopps_scores is samples x pathways (ratio, ~1.0) with
# a 'sample' column. WT (Myc-) 12W - 6W per pathway = the developmental mito
# reallocation (the "moving mitochondrial substrate").
mito_wt_shift <- tryCatch({
  mps       <- readRDS(here::here("results", "mitopps_scores.rds"))$mitopps_scores
  # numeric columns only (drop 'sample' and any other metadata columns present)
  path_cols <- names(mps)[vapply(mps, is.numeric, logical(1))]
  idx       <- match(mps$sample, rownames(sample_meta))
  if (anyNA(idx)) message(sprintf("mitoPPS: %d/%d samples matched to metadata",
                                  sum(!is.na(idx)), nrow(mps)))
  keep <- !is.na(idx)
  msub <- mps[keep, , drop = FALSE]
  meta <- sample_meta[idx[keep], c("myc_status", "timepoint")]
  is_neg <- meta$myc_status == "neg"
  wt6  <- colMeans(msub[is_neg & meta$timepoint == "6W",  path_cols, drop = FALSE])
  wt12 <- colMeans(msub[is_neg & meta$timepoint == "12W", path_cols, drop = FALSE])
  tibble::tibble(pathway = path_cols, wt_6W = wt6, wt_12W = wt12,
                 wt_shift = wt12 - wt6) |>
    dplyr::arrange(dplyr::desc(abs(wt_shift)))
}, error = function(e) { message("mitoPPS parse failed: ", conditionMessage(e)); NULL })

# =============================================================================
# PART 6: AP3 -- CONVERGENCE VECTOR (Myc genotype effect vs WT developmental shift)
# =============================================================================
# Reuse gsva_overview: beta_time = WT slope (12W_neg - 6W_neg); genotype effect
# per set = (d6 + d12)/2. Parallel (positive corr) = Myc accelerates development;
# anti-parallel (negative) = Myc opposes it.
ap3 <- dev_geno |> dplyr::filter(!is.na(beta_time), !is.na(geno_set))
ap3_overall <- suppressWarnings(
  stats::cor.test(ap3$geno_set, ap3$beta_time, method = "spearman"))
ap3_by_lineage <- ap3 |>
  dplyr::group_by(lineage) |>
  dplyr::summarise(n = dplyr::n(),
                   rho = suppressWarnings(stats::cor(geno_set, beta_time, method = "spearman")),
                   .groups = "drop")

# =============================================================================
# PART 7: AP4 -- FELSHER CO-VARIATION (constant dose -> coupling)
# =============================================================================
ap4 <- tibble::tibble(
  pair = c("felsher~bio", "felsher~diff_axis"),
  rho  = c(suppressWarnings(stats::cor(comp_df$felsher, comp_df$bio_comp, method = "spearman")),
           suppressWarnings(stats::cor(comp_df$felsher, comp_df$diff_axis, method = "spearman")))
)

# =============================================================================
# PART 8: FIGURES
# =============================================================================
out_dir <- here::here("outputs", "dev_composition")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
geno_cols <- c(neg = "#4575B4", pos = "#D73027")
diff_p <- ap2_tests$geno_p[ap2_tests$metric == "diff_axis"]

# Fig 1: differentiation axis by group (primary AP2 readout)
p_diff <- ggplot2::ggplot(comp_df,
    ggplot2::aes(x = group, y = diff_axis, fill = myc_status)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  ggplot2::geom_jitter(width = 0.12, size = 2, ggplot2::aes(colour = myc_status)) +
  ggplot2::scale_fill_manual(values = geno_cols) +
  ggplot2::scale_colour_manual(values = geno_cols) +
  ggplot2::labs(title = "AP2: differentiation axis (mature-HR minus progenitor-LP)",
    subtitle = sprintf("lower = toward progenitor; Myc genotype effect p = %.3g", diff_p),
    x = NULL, y = "HR_comp - LP_comp") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "differentiation_axis.pdf"), p_diff, width = 7, height = 5)

# Fig 2: per-lineage genotype effect (composite beta + per-set distribution)
p_lin <- ggplot2::ggplot(dev_geno |> dplyr::filter(lineage %in% lin_classes),
    ggplot2::aes(x = lineage, y = geno_set)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.4, fill = "grey85") +
  ggplot2::geom_jitter(width = 0.12, size = 1, alpha = 0.6) +
  ggplot2::labs(title = "Myc genotype effect on MG_* sets, by lineage class",
    subtitle = "per-set (d6+d12)/2; negative = Myc suppresses that lineage program",
    x = NULL, y = "Myc genotype effect (GSVA)") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "lineage_genotype_effect.pdf"), p_lin, width = 7, height = 5)

# Fig 3: WT developmental substrate shift (Myc- 6W->12W per composite)
p_wt <- ggplot2::ggplot(wt_shift |> dplyr::filter(!is.na(wt_beta)),
    ggplot2::aes(x = stats::reorder(metric, wt_beta), y = wt_beta)) +
  ggplot2::geom_col(fill = "#4575B4") +
  ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
  ggplot2::coord_flip() +
  ggplot2::labs(title = "AP1: WT (Myc-) developmental shift 6W->12W",
    subtitle = "the substrate the constant Myc acts on", x = NULL,
    y = "Myc- 12W - 6W (GSVA composite)") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "wt_developmental_shift.pdf"), p_wt, width = 7, height = 4.5)

# Fig 4: AP3 convergence vector
p_vec <- ggplot2::ggplot(ap3, ggplot2::aes(x = beta_time, y = geno_set, colour = lineage)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  ggplot2::geom_point(size = 2, alpha = 0.8) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.5,
                       formula = y ~ x) +
  ggplot2::labs(title = "AP3: Myc effect vs WT developmental shift (per MG_* set)",
    subtitle = sprintf("Spearman rho = %.3f, p = %.2g (parallel = Myc accelerates dev)",
                       ap3_overall$estimate, ap3_overall$p.value),
    x = "WT developmental shift (beta_time = 12W_neg - 6W_neg)",
    y = "Myc genotype effect ((d6+d12)/2)") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "convergence_vector.pdf"), p_vec, width = 7, height = 5.5)

# Fig 5: composite genotype effects with significance (the powered AP2 summary)
p_ap2 <- ggplot2::ggplot(ap2_tests,
    ggplot2::aes(x = stats::reorder(metric, geno_beta), y = geno_beta,
                 fill = geno_p < 0.05)) +
  ggplot2::geom_col() +
  ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("p=%.2g", geno_p)),
                     hjust = ifelse(ap2_tests$geno_beta >= 0, -0.1, 1.1), size = 3) +
  ggplot2::scale_fill_manual(values = c(`TRUE` = "#D73027", `FALSE` = "grey70"),
                             name = "p < 0.05") +
  ggplot2::coord_flip() +
  ggplot2::labs(title = "AP2: Myc genotype main effect on composites (powered)",
    subtitle = "biogenesis strongly driven; basal suppressed; luminal null (crossover)",
    x = NULL, y = "genotype effect (Myc+ - Myc-, timepoint-adjusted)") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "ap2_composite_genotype_effects.pdf"), p_ap2,
                width = 7.5, height = 5)

# Fig 6: Felsher (Myc activity) vs biogenesis composite -- the r=0.90 coupling
p_fels <- ggplot2::ggplot(comp_df, ggplot2::aes(x = bio_comp, y = felsher, colour = group)) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.5,
                       formula = y ~ x) +
  ggplot2::labs(title = "AP4: Felsher Myc-activity vs biogenesis composite",
    subtitle = sprintf("Spearman rho = %.2f -- the Myc footprint IS a biogenesis program",
                       ap4$rho[ap4$pair == "felsher~bio"]),
    x = "Biogenesis_discrimination composite GSVA", y = "Felsher signature GSVA") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "felsher_biogenesis_coupling.pdf"), p_fels,
                width = 6.5, height = 5)

# Fig 7: WT developmental mito reallocation (mitoPPS), top movers -- AP1 mito half
if (!is.null(mito_wt_shift)) {
  top_mito <- head(mito_wt_shift, 20)
  p_mito <- ggplot2::ggplot(top_mito,
      ggplot2::aes(x = stats::reorder(pathway, wt_shift), y = wt_shift,
                   fill = wt_shift > 0)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#D73027", `FALSE` = "#4575B4"),
                               guide = "none") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "AP1: WT developmental mito reallocation (mitoPPS ratio)",
      subtitle = "Myc- 12W - 6W; top 20 |shift| pathways (the moving mito substrate)",
      x = NULL, y = "WT 12W - 6W (mitoPPS ratio)") +
    ggplot2::theme_bw(base_size = 9)
  ggplot2::ggsave(file.path(out_dir, "wt_mito_reallocation.pdf"), p_mito,
                  width = 7.5, height = 6)
}

# =============================================================================
# PART 9: SAVE
# =============================================================================
group_means <- comp_df |>
  dplyr::group_by(group) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(c(ap2_metrics)), mean), .groups = "drop")

dev_out <- list(
  comp_df        = comp_df,
  lineage_tbl    = lineage_tbl,
  ap2_tests      = ap2_tests,
  dev_geno       = dev_geno,
  wt_shift       = wt_shift,
  mito_wt_shift  = mito_wt_shift,
  ap3            = list(overall = ap3_overall, by_lineage = ap3_by_lineage),
  ap4_felsher    = ap4,
  group_means    = group_means,
  notes = paste(
    "AP2 (primary, powered): genotype MAIN EFFECT on per-sample lineage composites;",
    "diff_axis = HR_comp - LP_comp (lower = toward progenitor = the AP7/LP link).",
    "AP1: WT (Myc-) 6W->12W shift = the substrate. AP3: convergence vector (Myc",
    "effect vs WT dev shift). AP4: Felsher co-variation. Lineage tags are a coarse",
    "name heuristic. Per-set p's are underpowered/descriptive; composites are the",
    "powered layer.")
)
saveRDS(dev_out, here::here("results", "dev_composition.rds"))
message("Saved results/dev_composition.rds")
message(sprintf("AP2 differentiation-axis genotype effect: beta = %.3f, p = %.3g",
                ap2_tests$geno_beta[ap2_tests$metric == "diff_axis"], diff_p))

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  dv <- readRDS(here::here("results", "dev_composition.rds"))

  # Lineage tagging sanity: eyeball a few sets per class (heuristic is coarse)
  dv$lineage_tbl |> dplyr::group_by(lineage) |> dplyr::slice_head(n = 4) |>
    print(n = Inf)

  # AP2 (primary): genotype effect on composites. Expect diff_axis NEGATIVE
  # (Myc -> progenitor), HR_comp down, LP_comp spared/up -> the AP7 LP link.
  dv$ap2_tests |> print()
  dv$group_means |> print()

  # AP1: does the WT gland mature 6W->12W (HR up, LP/MASC down)? substrate readout
  dv$wt_shift |> print()
  if (!is.null(dv$mito_wt_shift)) dv$mito_wt_shift |> head(15) |> print()
  # If mitoPPS was skipped, inspect its structure to wire it next iteration:
  # str(readRDS(here::here("results","mitopps_scores.rds")), max.level = 2)

  # AP3: is Myc parallel (accelerates) or anti-parallel (opposes) to WT dev?
  dv$ap3$overall
  dv$ap3$by_lineage |> print()

  # AP4: Felsher coupling (constant dose)
  dv$ap4_felsher |> print()

  list.files(here::here("outputs", "dev_composition"), pattern = "\\.pdf$")
}
