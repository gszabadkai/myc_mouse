# scripts/28_myc_mito_centrality.R
# =============================================================================
# Block A revision -- Issue #3: is the MYC effect REALLY "primarily mitochondrial"?
# =============================================================================
#
# The synthesis leads with "Myc preferentially amplifies mitochondrial biogenesis"
# (theme B / model lead A), resting on AP6.1 (script 20, MitoCarta NES 2.18) and
# AP6.2 (script 21, MitoCarta the top permutation-null z). The critique: mito NES is
# high, but the MYC signature is ALSO high and spans many arms (ribosome biogenesis,
# proliferation, metabolism) -- does mito have independent substance, or is
# "MYC effect = mitochondrial" just restating "MYC effect = MYC targets"?
#
# The one word "primarily" does two jobs; this script separates them.
#
#   Q1  PREFERENTIAL?  Is mito the single most-altered compartment? -- the overclaim.
#       AP6.2 ranks compartments by a permutation-null z, but z is SIZE-CONFOUNDED
#       (a larger set -> tighter matched null -> bigger z: MitoCarta n=1041 z=19.3 vs
#       OXPHOS n=185 z=10.9 at equal effect). By the size-fair metrics -- effect
#       MAGNITUDE (mean|LFC|) and housekeeping-corrected RELATIVE EXCESS over the
#       matched null ((obs-null_mean)/null_mean) -- the OXPHOS/biogenesis mito CORE
#       is top-tier (co-leading with the MYC-target core), the mito compartment
#       WHOLESALE is mid-pack, and MYC-targets OVERTAKE at 12W. Part A reframes the
#       existing null (no recompute) + adds a powered per-sample Cohen's-d ranking
#       across all arms incl. central metabolism. Verdict: soften "the single most
#       preferential" to "the OXPHOS/biogenesis core is top-tier, co-equal with the
#       MYC-target core".
#
#   Q2  CENTRAL?  Is the mito change a CORE, phenotype-coupled arm, or a bystander
#       that merely rides MYC dose? We cannot prove causation transcriptomically (no
#       clean mito-only KO of MYC's effect exists), but we can test whether mito
#       BEHAVES like a core arm. Part B couples each mito / metabolic axis to the
#       phenotypic OUTCOMES we already measured -- the MB2 tumorigenic fork (AP7,
#       script 18), mitochondrial death priming (script 23), TEB/dedifferentiation
#       (Issue #1/#2), proliferation -- and, critically, tests whether that coupling
#       SURVIVES removing the generic MYC/proliferation axis (PARTIAL correlation vs
#       the Felsher core, which is the 67-gene MYC phenotype core with the mito
#       signature stripped: 61 genes, 8% mito, 0% OXPHOS). An axis that tracks the
#       phenotype only through MYC dose is a bystander; one that tracks it BEYOND MYC
#       dose is central-by-coupling.
#
# Reframe on already-scored data (GSVA per-sample) + saved composites; no re-score.
#
# Input:  results/gsva_scores.rds          (scores + set_meta + sample_meta)
#         results/ap6_permutation_null.rds  (matched-null table, for the Q1 reframe)
#         results/ap7_mb_fork.rds           (MB2_over_MB1_UF per sample)
#         results/fgsea_percategory.rds     (per-category NES, context)
# Output: results/myc_mito_centrality.rds
#         outputs/myc_mito_centrality/*.pdf
#
# Ceiling: Q1 = powered (main effects / 24 samples). Q2 coupling/partial correlation
# are per-sample associations at n=24 (n=6/group) -> indicative, correlated sets, and
# co-regulation is NOT causation. The definitive centrality test is genetic.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# =============================================================================
# PART 1: LOAD + SAMPLE METADATA
# =============================================================================

gsva_out    <- readRDS(here::here("results", "gsva_scores.rds"))
scores      <- gsva_out$scores
set_meta    <- gsva_out$set_meta
sm          <- as.data.frame(gsva_out$sample_meta)
sm          <- sm[colnames(scores), , drop = FALSE]
stopifnot(identical(rownames(sm), colnames(scores)))
sm$timepoint  <- stats::relevel(as.factor(sm$timepoint),  "6W")
sm$myc_status <- stats::relevel(as.factor(sm$myc_status), "neg")
sm$group      <- factor(sm$group, levels = c("6W_neg", "6W_pos", "12W_neg", "12W_pos"))

ap6 <- readRDS(here::here("results", "ap6_permutation_null.rds"))
ap7 <- readRDS(here::here("results", "ap7_mb_fork.rds"))$fork_df
stopifnot(all(ap7$sample %in% colnames(scores)))
ap7 <- ap7[match(colnames(scores), ap7$sample), ]

# =============================================================================
# PART 2: COMPARTMENT + METABOLIC-AXIS PER-SAMPLE COMPOSITES
# =============================================================================
comp <- function(sets) {
  sets <- intersect(sets, rownames(scores))
  stopifnot(length(sets) > 0)
  colMeans(scores[sets, , drop = FALSE])
}
by_cat  <- function(cat) set_meta$set_name[set_meta$category_primary == cat]
mito_names <- by_cat("MitoCarta")
mito_grep  <- function(pat) grep(pat, mito_names, value = TRUE)

# --- MYC identity arms ---
myc_sig     <- comp(by_cat("MYC_signatures"))                 # 17-set MYC composite
felsher     <- scores["MYC_felsher_integrative_signature", ]  # mito-stripped MYC core
hallmark_v2 <- scores["MYC_HALLMARK_MYC_TARGETS_V2", ]

# --- Mitochondrial arms (whole vs the OXPHOS/biogenesis core) ---
mito_all         <- if ("MITOCARTA_ALL" %in% rownames(scores)) scores["MITOCARTA_ALL", ] else comp(mito_names)
mito_oxphos      <- comp(mito_grep("OXPHOS|COMPLEX_[IV]|_SUBUNITS|ASSEMBLY_FACTORS|ELECTRON_CARRIERS|CRISTAE"))
mito_biogenesis  <- comp(mito_grep("RIBOSOME|CENTRAL_DOGMA|MT_TRNA|MT_RRNA|MTRNA|MTDNA|IMPORT|TRANSLATION"))

# --- Proliferation ---
prolif <- comp(by_cat("Proliferation"))                       # 14-set

# --- Central metabolic axes (cancer-metabolism remit) ---
metab_axes <- list(
  glycolysis  = c("GS_METAB_GLYCOLYSIS", "METAB_GLYCOLYSIS_HALLMARK", "METAB_GLYCOLYSIS_KEGG", "METAB_GLYCOLYSIS_WP"),
  tca         = c("GS_METAB_KREBS", "METAB_TCA_KEGG", "METAB_TCA_REACTOME"),
  oxphos_met  = c("METAB_OXPHOS_HALLMARK", "METAB_OXPHOS_KEGG", "METAB_OXPHOS_REACTOME",
                  "GS_METAB_COMPLEX_I", "GS_METAB_COMPLEX_III", "GS_METAB_COMPLEX_IV",
                  "GS_METAB_PROTON_TRANSPORT", "GS_METAB_UBIQUINONE"),
  ppp         = c("GS_METAB_PENTOSE_PHOSPHATE", "METAB_PPP_KEGG", "METAB_PPP_REACTOME"),
  nucleotide  = c("GS_METAB_NUCLEOTIDE", "GS_METAB_PURINE", "GS_METAB_PYRIMIDINE",
                  "METAB_NUCLEOTIDE_REACTOME", "METAB_NUCLEOTIDE_SALVAGE_REACTOME",
                  "METAB_NUCLEOTIDE_WP", "METAB_PURINE_KEGG", "METAB_PYRIMIDINE_KEGG", "METAB_PYRIMIDINE_WP"),
  one_carbon  = c("METAB_ONE_CARBON_KEGG", "METAB_SER_GLY_KEGG", "GS_METAB_FOLATE",
                  "METAB_FOLATE_WP", "GS_METAB_METHIONINE", "METAB_MET_CYS_KEGG"),
  glutamine   = c("GS_METAB_GLUTAMATE", "METAB_GLN_GLU_KEGG", "METAB_GLN_GLU_REACTOME"),
  fao         = c("METAB_FA_BETAOX_REACTOME", "METAB_FA_HALLMARK", "METAB_FA_KEGG", "GS_METAB_FATTY_ACID"),
  fa_synth    = c("METAB_FA_SYNTHESIS_REACTOME", "METAB_FA_SYNTHESIS_WP"),
  amino_acid  = c("GS_METAB_AMINO_ACID", "GS_METAB_BCAAS", "METAB_BCAA_KEGG", "METAB_BCAA_REACTOME", "METAB_BCAA_WP"),
  cholesterol = c("GS_METAB_CHOLESTEROL", "GS_METAB_MEVALONATE", "METAB_CHOLESTEROL_HALLMARK",
                  "METAB_CHOLESTEROL_REACTOME", "METAB_CHOLESTEROL_WP"),
  redox       = c("GS_METAB_REDOX", "GS_METAB_GLUTATHIONE", "GS_METAB_REACTIVE_OXYGEN"))
metab_mat <- vapply(metab_axes, comp, numeric(ncol(scores)))  # samples x axes

# --- Phenotypic OUTCOME anchors (Q2 targets) ---
mb2_fork   <- ap7$MB2_over_MB1_UF                              # tumorigenic Myc fork (AP7)
priming    <- scores["MITOCARTA_APOPTOSIS_PRO", ] - scores["MITOCARTA_APOPTOSIS_ANTI", ]  # mito death priming
teb_up     <- grep("^MG_TEB_VS_DUCTAL_.*_UP$", rownames(scores), value = TRUE)
teb_dn     <- grep("^MG_TEB_VS_DUCTAL_.*_DN$", rownames(scores), value = TRUE)
teb_dediff <- comp(teb_up) - comp(teb_dn)                     # TEB/dedifferentiation (Issue #1/#2)

# One tidy panel of every per-sample series
panel <- data.frame(
  sample = colnames(scores), group = sm$group,
  timepoint = sm$timepoint, myc_status = sm$myc_status,
  myc_sig = myc_sig, felsher = felsher, hallmark_v2 = hallmark_v2,
  mito_all = mito_all, mito_oxphos = mito_oxphos, mito_biogenesis = mito_biogenesis,
  prolif = prolif,
  mb2_fork = mb2_fork, priming = priming, teb_dediff = teb_dediff,
  metab_mat, check.names = FALSE)

message(sprintf("Panel built: %d samples x %d series (mito arms 3, MYC 3, prolif 1, metab %d, outcomes 3)",
                nrow(panel), ncol(panel) - 4L, length(metab_axes)))

# =============================================================================
# PART A1: Q1 -- reframe the AP6.2 matched null (size-fair, no recompute)
# =============================================================================
# z is size-confounded; report magnitude + housekeeping-corrected relative excess.
null_reframe <- ap6$null_table |>
  dplyr::mutate(
    rel_excess = (obs - null_mean) / null_mean,   # size-robust % over matched null
    magnitude  = obs) |>
  dplyr::group_by(contrast) |>
  dplyr::mutate(
    rank_by_z         = rank(-z),
    rank_by_excess    = rank(-rel_excess),
    rank_by_magnitude = rank(-magnitude)) |>
  dplyr::ungroup() |>
  dplyr::arrange(contrast, dplyr::desc(rel_excess)) |>
  dplyr::select(compartment, contrast, n, obs, null_mean, z,
                rel_excess, rank_by_z, rank_by_excess, rank_by_magnitude)

# =============================================================================
# PART A2: Q1 -- powered per-sample Cohen's-d ranking (all arms + metabolism)
# =============================================================================
# Genotype MAIN effect (Myc+ - WT, timepoint-adjusted) / within-group SD. Same lm
# pattern as script 27 prog_stat_one. Size-confound-free (per-sample enrichment).
arm_series <- c("myc_sig", "felsher", "hallmark_v2",
                "mito_all", "mito_oxphos", "mito_biogenesis", "prolif",
                names(metab_axes))

geno_effect_one <- function(nm) {
  d   <- data.frame(y = panel[[nm]], myc = panel$myc_status, tp = panel$timepoint, grp = panel$group)
  wsd <- sqrt(mean(tapply(d$y, d$grp, stats::var)))
  ma  <- summary(stats::lm(y ~ myc + tp, d))$coefficients["mycpos", ]
  mi  <- summary(stats::lm(y ~ tp * myc, d))$coefficients["tp12W:mycpos", ]
  gm  <- tapply(d$y, d$grp, mean)
  tibble::tibble(
    arm = nm, within_sd = wsd,
    geno_beta = unname(ma["Estimate"]), geno_d = unname(ma["Estimate"]) / wsd,
    geno_p = unname(ma["Pr(>|t|)"]),
    int_beta = unname(mi["Estimate"]), int_p = unname(mi["Pr(>|t|)"]),
    m6_neg = unname(gm["6W_neg"]), m6_pos = unname(gm["6W_pos"]),
    m12_neg = unname(gm["12W_neg"]), m12_pos = unname(gm["12W_pos"]))
}
compartment_d <- dplyr::bind_rows(lapply(arm_series, geno_effect_one)) |>
  dplyr::mutate(
    arm_type = dplyr::case_when(
      arm %in% c("myc_sig", "felsher", "hallmark_v2")                 ~ "MYC identity",
      arm %in% c("mito_all", "mito_oxphos", "mito_biogenesis")        ~ "Mito arm",
      arm == "prolif"                                                 ~ "Proliferation",
      TRUE                                                            ~ "Metabolic axis"),
    verdict = dplyr::case_when(geno_p < 0.001 ~ "powered (strong)",
                               geno_p < 0.05  ~ "powered",
                               TRUE           ~ "weak / ns")) |>
  dplyr::arrange(dplyr::desc(geno_d))

# =============================================================================
# PART B: Q2 -- centrality by coupling to phenotypic outcomes (+ partial)
# =============================================================================
outcomes <- c(mb2_fork = "MB2 tumorigenic fork", priming = "mito death priming",
              teb_dediff = "TEB / dedifferentiation", prolif = "proliferation")
axes     <- c("mito_all", "mito_oxphos", "mito_biogenesis",
              names(metab_axes), "myc_sig", "hallmark_v2")

scor <- function(a, b, idx = TRUE) suppressWarnings(
  stats::cor(panel[[a]][idx], panel[[b]][idx], method = "spearman"))

# B1 -- raw coupling: axis x outcome, overall + within genotype + within timepoint
subsets <- list(
  all    = rep(TRUE, nrow(panel)),
  wt     = panel$myc_status == "neg",
  mycpos = panel$myc_status == "pos",
  at6    = panel$timepoint == "6W",
  at12   = panel$timepoint == "12W")
coupling <- purrr::map_dfr(names(outcomes), function(oc)
  purrr::map_dfr(axes, function(ax)
    purrr::map_dfr(names(subsets), function(sub)
      tibble::tibble(outcome = oc, axis = ax, subset = sub,
                     rho = scor(ax, oc, subsets[[sub]])))))

# B2 -- PARTIAL correlation: does axis track outcome BEYOND the MYC axis? partial
# rho(axis, outcome | z) = spearman of the rank-residuals after regressing each on z.
# Deconfounders, weakest -> strongest MYC proxy:
#   rho_pcore = | Felsher 67-gene core (mito-stripped MYC phenotype signature)
#   rho_pmyc  = | full 17-set MYC composite (the strongest MYC-dose proxy) -- the
#               robustness that attacks "it was just MYC dose the Felsher proxy missed"
#   rho_pprol = | proliferation (is the coupling just generic proliferation?)
# CAUTION: surviving = the phenotype signal is ORTHOGONAL to the MYC-signature axis,
# NOT that it is MYC-INDEPENDENT or causal (proxy under-capture / non-linear MYC route
# / common cause all survive too). Direction is unidentified. See notes + revision log.
partial_rho <- function(a, b, z) {
  if (identical(a, z) || identical(b, z)) return(NA_real_)
  rx <- stats::residuals(stats::lm(rank(panel[[a]]) ~ rank(panel[[z]])))
  ry <- stats::residuals(stats::lm(rank(panel[[b]]) ~ rank(panel[[z]])))
  suppressWarnings(stats::cor(rx, ry, method = "pearson"))  # pearson on ranks = spearman-partial
}
partial <- purrr::map_dfr(names(outcomes), function(oc)
  purrr::map_dfr(setdiff(axes, oc), function(ax)
    tibble::tibble(
      outcome    = oc, axis = ax,
      rho_raw    = scor(ax, oc),
      rho_pcore  = partial_rho(ax, oc, "felsher"),   # beyond MYC phenotype core
      rho_pmyc   = partial_rho(ax, oc, "myc_sig"),    # beyond full MYC composite (robust)
      rho_pprol  = if (oc == "prolif") NA_real_ else partial_rho(ax, oc, "prolif"),
      # within-WT raw coupling (no transgene varying) = the cleanest "beyond transgene"
      rho_wt     = scor(ax, oc, panel$myc_status == "neg")))) |>
  dplyr::mutate(
    survives_core = abs(rho_pcore) >= 0.30 & sign(rho_pcore) == sign(rho_raw),
    survives_myc  = !is.na(rho_pmyc) & abs(rho_pmyc) >= 0.30 & sign(rho_pmyc) == sign(rho_raw),
    centrality = dplyr::case_when(
      abs(rho_raw) < 0.30           ~ "not coupled",
      survives_core & survives_myc  ~ "central (robust: beyond full MYC)",
      survives_core                 ~ "central (beyond MYC core)",
      TRUE                          ~ "rides MYC core (bystander)"))

# B3 -- metabolic-axis centrality summary: powered effect (Part A2 d) x max
# outcome-coupling-beyond-core -> which metabolic axes are BOTH moved and central.
metab_summary <- compartment_d |>
  dplyr::filter(arm_type == "Metabolic axis") |>
  dplyr::select(axis = arm, geno_d, geno_p) |>
  dplyr::left_join(
    partial |> dplyr::group_by(axis) |>
      dplyr::summarise(max_rho_raw = rho_raw[which.max(abs(rho_raw))],
                       max_rho_core = rho_pcore[which.max(abs(rho_pcore))],
                       max_rho_myc  = rho_pmyc[which.max(abs(rho_pcore))],
                       top_outcome = outcome[which.max(abs(rho_pcore))], .groups = "drop"),
    by = "axis") |>
  dplyr::arrange(dplyr::desc(geno_d))

# =============================================================================
# PART C: FIGURES
# =============================================================================
out_dir <- here::here("outputs", "myc_mito_centrality")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
mito_pal <- c("MYC identity" = "#4575B4", "Mito arm" = "#D73027",
              "Proliferation" = "#7B3294", "Metabolic axis" = "#1A9850")

# C1 -- Q1 reframe: matched-null relative excess (size-fair), both contrasts
p_c1 <- null_reframe |>
  dplyr::mutate(compartment = stats::reorder(compartment, rel_excess)) |>
  ggplot2::ggplot(ggplot2::aes(x = rel_excess, y = compartment, colour = contrast)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
  ggplot2::scale_colour_manual(values = c(myc_6W = "#D73027", myc_12W = "#4575B4")) +
  ggplot2::labs(
    title = "Q1: size-fair preferentiality (housekeeping-corrected excess over matched null)",
    subtitle = "OXPHOS/MYC-target core lead; MitoCarta-wholesale mid-pack; MYC-targets overtake at 12W",
    x = "relative excess (obs - null_mean) / null_mean", y = NULL, colour = "contrast") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "q1_preferentiality_reframed.pdf"), p_c1, width = 8, height = 4.5)

# C2 -- Q1 powered: genotype Cohen's d across all arms
p_c2 <- compartment_d |>
  dplyr::mutate(arm = stats::reorder(arm, geno_d)) |>
  ggplot2::ggplot(ggplot2::aes(x = geno_d, y = arm, fill = arm_type)) +
  ggplot2::geom_col() +
  ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
  ggplot2::scale_fill_manual(values = mito_pal) +
  ggplot2::labs(
    title = "Q1 powered: Myc genotype main effect (Cohen's d) across arms",
    subtitle = "per-sample, size-confound-free; is mito the biggest arm, and which metabolic axes move?",
    x = "genotype-effect Cohen's d (Myc+ - WT, timepoint-adjusted)", y = NULL, fill = NULL) +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "q1_genotype_effect_ranking.pdf"), p_c2, width = 8, height = 7)

# C3 -- Q2 coupling heatmap: axis x outcome (overall rho)
p_c3 <- coupling |>
  dplyr::filter(subset == "all") |>
  dplyr::mutate(outcome = factor(outcome, levels = names(outcomes), labels = outcomes)) |>
  ggplot2::ggplot(ggplot2::aes(x = outcome, y = axis, fill = rho)) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho)), size = 2.6) +
  ggplot2::scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                                midpoint = 0, limits = c(-1, 1)) +
  ggplot2::labs(title = "Q2: axis x phenotypic-outcome coupling (Spearman, all samples)",
                x = NULL, y = NULL) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
ggplot2::ggsave(file.path(out_dir, "q2_coupling_heatmap.pdf"), p_c3, width = 7.5, height = 7)

# C4 -- Q2 partial: raw vs beyond-MYC-core coupling (does it survive?). Colour = the
# FULL 4-tier verdict (robust tier must be in the scale, else it silently greys out).
cent_levels <- c("central (robust: beyond full MYC)", "central (beyond MYC core)",
                 "rides MYC core (bystander)", "not coupled")
cent_cols <- c("central (robust: beyond full MYC)" = "#1A9850",   # dark green
               "central (beyond MYC core)"         = "#A6D96A",   # light green
               "rides MYC core (bystander)"        = "#D73027",   # red
               "not coupled"                       = "grey60")
p_c4 <- partial |>
  dplyr::mutate(outcome = factor(outcome, levels = names(outcomes), labels = outcomes),
                centrality = factor(centrality, levels = cent_levels)) |>
  ggplot2::ggplot(ggplot2::aes(x = rho_raw, y = rho_pcore, colour = centrality)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  ggplot2::geom_hline(yintercept = c(-0.3, 0.3), linetype = "dotted", colour = "grey70") +
  ggplot2::geom_point(size = 2.4) +
  ggrepel::geom_text_repel(ggplot2::aes(label = axis), size = 2.3, max.overlaps = 20) +
  ggplot2::facet_wrap(~ outcome) +
  ggplot2::scale_colour_manual(values = cent_cols, drop = FALSE) +
  ggplot2::labs(
    title = "Q2 partial correlation: coupling raw vs beyond the Felsher (MYC) core",
    subtitle = paste("on-diagonal = coupling survives (central); collapse toward 0 = MYC-dose bystander.",
                     "Dark green also survives the FULL MYC composite. Orthogonal to MYC-signature != causal."),
    x = "raw Spearman rho", y = "partial rho | Felsher core", colour = NULL) +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "q2_partial_correlation.pdf"), p_c4, width = 10, height = 7.5)

# C5 -- metabolic centrality: effect size x coupling. Sign-aware: a strong NEGATIVE
# partial is coupling too (inverse tracking), so plot |partial| and mark direction +
# robustness (survives the full MYC composite), not just a signed upper-right quadrant.
p_c5 <- metab_summary |>
  dplyr::mutate(
    robust    = !is.na(max_rho_myc) & abs(max_rho_myc) >= 0.30,
    direction = ifelse(max_rho_core >= 0, "positive", "inverse"),
    label     = sprintf("%s (%s%s)", axis, top_outcome, ifelse(direction == "inverse", ", -", ""))) |>
  ggplot2::ggplot(ggplot2::aes(x = geno_d, y = abs(max_rho_core),
                               colour = robust, shape = direction)) +
  ggplot2::geom_hline(yintercept = 0.3, linetype = "dotted", colour = "grey70") +
  ggplot2::geom_point(size = 3) +
  ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 2.4, max.overlaps = 20) +
  ggplot2::scale_colour_manual(values = c(`TRUE` = "#1A9850", `FALSE` = "grey55"),
                               labels = c(`TRUE` = "survives full MYC", `FALSE` = "Felsher-only / weak")) +
  ggplot2::scale_shape_manual(values = c(positive = 16, inverse = 17)) +
  ggplot2::labs(
    title = "Central metabolic axes: Myc effect size x phenotype coupling (magnitude)",
    subtitle = "upper band = coupled beyond MYC dose (either sign); green = also survives the full MYC composite",
    x = "genotype-effect Cohen's d", y = "|max partial rho to a phenotype| | Felsher core",
    colour = NULL, shape = "coupling sign") +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "metabolic_centrality.pdf"), p_c5, width = 8, height = 6)

# =============================================================================
# PART D: SAVE
# =============================================================================
mito_centrality <- list(
  panel         = tibble::as_tibble(panel),
  null_reframe  = null_reframe,
  compartment_d = compartment_d,
  coupling      = coupling,
  partial       = partial,
  metab_summary = metab_summary,
  defs = list(
    mito_oxphos_sets     = mito_grep("OXPHOS|COMPLEX_[IV]|_SUBUNITS|ASSEMBLY_FACTORS|ELECTRON_CARRIERS|CRISTAE"),
    mito_biogenesis_sets = mito_grep("RIBOSOME|CENTRAL_DOGMA|MT_TRNA|MT_RRNA|MTRNA|MTDNA|IMPORT|TRANSLATION"),
    metab_axes           = metab_axes,
    outcomes             = outcomes),
  notes = paste(
    "Issue #3. Q1 (PREFERENTIAL): AP6.2's mito #1 rank is a SIZE artifact of the",
    "permutation-null z; by size-fair relative-excess/magnitude the OXPHOS/biogenesis",
    "mito CORE is top-tier (co-leading the MYC-target core), MitoCarta-wholesale is",
    "mid-pack, MYC-targets overtake at 12W -> soften 'the single most preferential'.",
    "compartment_d = powered per-sample genotype Cohen's d across arms + metabolism.",
    "Q2 (CENTRAL): coupling = axis x outcome Spearman (MB2 fork / mito priming / TEB",
    "dedifferentiation / proliferation), overall + by genotype + by timepoint; partial",
    "deconfounds against Felsher core (rho_pcore), the full 17-set MYC composite",
    "(rho_pmyc, strongest MYC proxy), and proliferation (rho_pprol); rho_wt = within-WT",
    "raw coupling (no transgene varying). Tiered verdict: 'central (robust)' survives",
    "BOTH Felsher and full-MYC; 'central (beyond MYC core)' survives Felsher only.",
    "INTERPRETATION CEILING (critical): surviving the partial means the phenotype signal",
    "is ORTHOGONAL to the MYC-signature axis -- NOT that it is MYC-INDEPENDENT or causal.",
    "Proxy under-capture (MYC dose the signature misses), a non-linear MYC route, and a",
    "common cause all survive too; direction is unidentified. rho_pmyc + within-WT",
    "coupling tighten the 'beyond MYC dose' reading but cannot establish independence or",
    "mediation. Co-variation says NOTHING about necessity (a MYC-readout arm may still be",
    "required). Per-sample n=6/group, correlated sets, priming mito-defined (mito<->priming",
    "partly circular; nucleotide<->priming is the non-circular corroborator). Definitive",
    "test is genetic. Association layer -> Block B. See docs/2026-07-08_BlockA_revision_plan.md."))
saveRDS(mito_centrality, here::here("results", "myc_mito_centrality.rds"))
message("Saved results/myc_mito_centrality.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  mc <- readRDS(here::here("results", "myc_mito_centrality.rds"))

  # --- Q1a: the size confound made explicit -- z rank vs size-fair excess rank ---
  mc$null_reframe |> dplyr::filter(contrast == "myc_6W") |>
    dplyr::select(compartment, n, obs, rel_excess, rank_by_z, rank_by_excess) |> print()
  mc$null_reframe |> dplyr::filter(contrast == "myc_12W") |>
    dplyr::arrange(dplyr::desc(rel_excess)) |>
    dplyr::select(compartment, obs, rel_excess, rank_by_excess) |> print()

  # --- Q1b: powered arm ranking -- is mito the biggest? where do MYC-targets sit? ---
  mc$compartment_d |>
    dplyr::select(arm, arm_type, geno_d, geno_p, verdict) |> print(n = Inf)

  # --- Q2: which axes are CENTRAL, and do they survive the FULL MYC composite? ---
  # The robust tier (survives both Felsher AND the 17-set MYC composite) is the one
  # that resists "it was just MYC dose the Felsher proxy missed".
  mc$partial |> dplyr::filter(grepl("^central", centrality)) |>
    dplyr::arrange(dplyr::desc(abs(rho_pmyc))) |>
    dplyr::select(outcome, axis, rho_raw, rho_pcore, rho_pmyc, rho_wt, centrality) |>
    print(n = Inf)
  # mito arms specifically, per outcome (raw -> |Felsher -> |fullMYC -> within-WT):
  mc$partial |> dplyr::filter(axis %in% c("mito_all", "mito_oxphos", "mito_biogenesis")) |>
    dplyr::select(outcome, axis, rho_raw, rho_pcore, rho_pmyc, rho_wt, centrality) |>
    print(n = Inf)
  # How many 'central' rows survive the full-MYC deconfound (robust) vs Felsher-only?
  mc$partial |> dplyr::filter(grepl("^central", centrality)) |>
    dplyr::count(centrality) |> print()

  # --- Coupling tightest at 6W_pos (the phenotype-strong condition)? ---
  mc$coupling |> dplyr::filter(axis == "mito_oxphos") |>
    tidyr::pivot_wider(names_from = subset, values_from = rho) |> print()

  # --- Metabolic centrality: which metabolic axes are moved AND central? ---
  mc$metab_summary |> print(n = Inf)

  # Sanity: mito_oxphos vs mito_biogenesis vs mito_all correlation (arms coherent?)
  cat(sprintf("mito_all~oxphos %.2f  all~biog %.2f  oxphos~biog %.2f\n",
              stats::cor(mc$panel$mito_all, mc$panel$mito_oxphos),
              stats::cor(mc$panel$mito_all, mc$panel$mito_biogenesis),
              stats::cor(mc$panel$mito_oxphos, mc$panel$mito_biogenesis)))

  list.files(here::here("outputs", "myc_mito_centrality"), pattern = "\\.pdf$")
}
