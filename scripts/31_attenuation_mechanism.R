# scripts/31_attenuation_mechanism.R
# =============================================================================
# Block A revision -- Issue #6: a CITABLE MECHANISM for the OXPHOS attenuation
# =============================================================================
#
# Issues #4/#5 fixed WHAT the attenuation is (magnitude compression, not rank; Myc raises
# absolute nuclear OXPHOS) and that no developmental/TF axis ABSORBS the genotype x time
# interaction. What was still missing is a CITABLE mechanism -- "compression" only renames the
# attenuation. Issue #6 delivers one, using an exact identity in the interaction model.
#
# THE IDENTITY (verified, cor = 1.0000): the per-gene attenuation is
#   myc_12W - myc_6W  ==  timepoint_pos - timepoint_neg
# i.e. the change in the genotype gap 6W->12W equals (Myc+ temporal change) - (WT temporal
# change). Aligning each gene to Myc's induction direction d = sign(myc_6W), the gap-shrink
# splits cleanly into TWO mechanisms:
#   WT-convergence  =  d * timepoint_neg   (wild-type tissue matures TOWARD Myc's state --
#                                            the "moving background")
#   Myc-fade        =  d * timepoint_pos   (the oncogenic program retreats on the aging
#                                            substrate; negative = retreat)
#   attenuation     =  aligned_gap6 - aligned_gap12  =  WT-convergence - Myc-fade
#
# PART A finding (preview, all 2777 6W-divergent genes): the attenuation is ~66% Myc-fade,
# ~34% WT-convergence, and it is PROGRAM-SPECIFIC -- WT-convergence is real for the
# biosynthetic arm (amino-acid 37%, lipid 41%, nucleotide 31%) but NEGATIVE for OXPHOS (WT
# moves AWAY; the gap closes purely by fade) and for the MYC-target core (which is ALSO the
# least attenuated -- protected). This is the POWERED, gene-level resolution of the old
# H1 (front-loaded fade) vs H3 (WT catch-up) question that Gate 1 / script 17 left directional.
#
# The dominant Myc-fade term is the ambiguous one: per-cell weakening OR compositional dilution
# (the Myc-responsive proliferative/TEB compartment shrinking with maturation). PART B (a
# broadened-TF absorption panel, extending Issue #5) both completes the "normalise for other TF
# activity" check AND diagnoses the compositional hypothesis: if the E2F/proliferation TF axis
# co-declines with the OXPHOS fade, that SUPPORTS (does not prove) dilution. PART C (reference
# deconvolution) is the definitive test and is DEFERRED (no sc reference / deconv pkg on disk).
#
# Reframe on already-fitted data (DESeq2 contrasts, GSVA, VST); no DESeq/GSVA re-run.
#
# Input:  results/interaction_results.rds        (raw LFC contrasts -- the identity)
#         results/combined_df_annotated.rds       (mgi_symbol <-> ensembl)
#         data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt  (986 pathways)
#         results/gsva_scores.rds                 (TF composites, Part B)
#         results/dds_int_run.rds                 (VST for the oxphos_abs outcome)
#         results/attenuation_decomposition.rds   ($defs: roster, nuclear_oxphos_ens, luminal_sets,
#                                                   mito_ox_gsva_sets)
#         data/genesets_from_library/gray_chea_mito_tf_shortlist.csv  (data-driven mito-TF panel)
# Output: results/attenuation_mechanism.rds
#         outputs/attenuation_mechanism/*.pdf
#
# CEILING: Part A is POWERED (contrast-level, all genes; the identity is exact) -> the
# convergence/fade split + H1/H3 resolution are solid but DESCRIPTIVE of the transcriptional
# change (composition still confounds the fade term; not a per-cell causal claim). Part B is
# BOUNDED exactly as Issue #5: endogenous TFs (Myc drives them), n=6/group, directional
# baseline, collinear axes -> absorption bounds; a TF co-declining with the fade SUPPORTS but
# does not prove dilution. Definitive = deconvolution/single-cell (Part C, deferred).
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# =============================================================================
# PART 1: LOAD + IDENTIFIER MAP + HELPERS
# =============================================================================
ir   <- readRDS(here::here("results", "interaction_results.rds"))
cdf  <- readRDS(here::here("results", "combined_df_annotated.rds"))
gmt  <- fgsea::gmtPathways(
  here::here("data", "genesets_from_library", "mammary_mito_myc_metab_v1_mouse.gmt"))
ad   <- readRDS(here::here("results", "attenuation_decomposition.rds"))
defs <- ad$defs

sym2ens <- cdf[!is.na(cdf$mgi_symbol) & !duplicated(cdf$mgi_symbol), c("mgi_symbol", "gene")]
universe_all <- rownames(ir$myc_6W_raw)
# Vintage-aware symbol -> Ensembl for gene-set membership (recovers renamed genes
# such as ATP synthase; a plain current-symbol match drops ~12% of nuclear OXPHOS).
source(here::here("functions", "reconcile_gene_symbols.R"))
ens_of <- function(syms) recon_to_ensembl(syms, universe_all)
pool_ens <- function(paths) ens_of(unique(unlist(gmt[intersect(paths, names(gmt))])))

# raw (unshrunken) LFC vectors indexed by Ensembl (the interaction identity uses raw LFCs)
V <- function(k) stats::setNames(ir[[k]]$log2FoldChange, rownames(ir[[k]]))
m6 <- V("myc_6W_raw"); m12 <- V("myc_12W_raw")
tn <- V("timepoint_neg_raw"); tp <- V("timepoint_pos_raw")
p6 <- stats::setNames(ir$myc_6W_raw$padj, rownames(ir$myc_6W_raw))

# =============================================================================
# PART A: CONVERGENCE / FADE DECOMPOSITION (POWERED, contrast-level)
# =============================================================================
# identity check: gap-change == Myc+ temporal - WT temporal
g_all <- Reduce(intersect, list(names(m6), names(m12), names(tn), names(tp)))
ident <- stats::cor((m12[g_all] - m6[g_all]), (tp[g_all] - tn[g_all]), use = "complete.obs")
message(sprintf("Identity check cor((m12-m6),(tp-tn)) = %.4f", ident))
stopifnot(ident > 0.999)

# per-gene aligned decomposition on a chosen gene set within a chosen universe
decomp_one <- function(genes, label, arm, universe_genes) {
  e <- intersect(genes, universe_genes)
  e <- e[!is.na(m6[e]) & !is.na(m12[e]) & !is.na(tn[e]) & !is.na(tp[e]) & m6[e] != 0]
  if (length(e) < 5) return(NULL)
  d  <- sign(m6[e])                                   # Myc induction direction
  gap6  <- d * m6[e];  gap12 <- d * m12[e]            # aligned gap (positive = induction)
  wc <- d * tn[e];     mf <- d * tp[e]                # aligned WT temporal / Myc+ temporal
  atten <- mean(gap6) - mean(gap12)                   # = mean(wc) - mean(mf) (the identity)
  tibble::tibble(
    program = label, arm = arm, n = length(e),
    gap6 = mean(gap6), gap12 = mean(gap12), atten = atten,
    wt_conv = mean(wc), myc_fade = mean(mf),
    conv_pct = 100 * mean(wc) / atten, fade_pct = 100 * (-mean(mf)) / atten,
    # correlated genes -> report robust spread + sign fraction, NOT an independence-p
    wt_conv_median = stats::median(wc), wt_conv_iqr = stats::IQR(wc),
    frac_wt_toward = mean(wc > 0), frac_myc_retreat = mean(mf < 0))
}

# program roster: the script 29 mito/MYC arms + two non-mito representatives (what ELSE
# converges): pooled Proliferation + pooled mammary-luminal (the WT-declining dev axis).
roster <- tibble::tribble(
  ~program,                            ~arm,                                          ~paths,
  "MITOCARTA_OXPHOS_SUBUNITS",         "OXPHOS core",                                 "MITOCARTA_OXPHOS_SUBUNITS",
  "MITOCARTA_OXPHOS",                  "OXPHOS core",                                 "MITOCARTA_OXPHOS",
  "MITOCARTA_TCA_CYCLE",               "TCA",                                         "MITOCARTA_TCA_CYCLE",
  "MITOCARTA_NUCLEOTIDE_METABOLISM",   "nucleotide",                                  "MITOCARTA_NUCLEOTIDE_METABOLISM",
  "MITOCARTA_AMINO_ACID_METABOLISM",   "biosynthetic",                                "MITOCARTA_AMINO_ACID_METABOLISM",
  "MITOCARTA_LIPID_METABOLISM",        "biosynthetic",                                "MITOCARTA_LIPID_METABOLISM",
  "MITOCARTA_MITOCHONDRIAL_RIBOSOME",  "biogenesis/translation",                      "MITOCARTA_MITOCHONDRIAL_RIBOSOME",
  "MYC_HALLMARK_MYC_TARGETS_V2",       "MYC-target core",                             "MYC_HALLMARK_MYC_TARGETS_V2",
  "MYC_felsher_integrative_signature", "MYC-target core",                             "MYC_felsher_integrative_signature")
prog_ens <- lapply(stats::setNames(roster$program, roster$program), function(p) ens_of(gmt[[p]]))
# two pooled non-mito representatives
prog_ens[["PROLIFERATION_pooled"]]   <- pool_ens(grep("^PROLIF_", names(gmt), value = TRUE))
prog_ens[["MAMMARY_LUMINAL_pooled"]] <- pool_ens(defs$luminal_sets)
roster <- dplyr::bind_rows(roster, tibble::tibble(
  program = c("PROLIFERATION_pooled", "MAMMARY_LUMINAL_pooled"),
  arm = c("proliferation", "mammary-dev"), paths = NA_character_))

# two universes: significance-selected (divergent) and effect-based (selection-independent)
universes <- list(
  divergent = names(m6)[!is.na(p6) & p6 < 0.1],       # padj<0.1 @6W (matches Issue #4)
  effect    = names(m6)[!is.na(m6) & abs(m6) > 0.5])  # |LFC6|>0.5, significance-independent
message(sprintf("Universes: divergent n=%d ; effect(|LFC6|>0.5) n=%d",
                length(universes$divergent), length(universes$effect)))

decomp_conv <- dplyr::bind_rows(lapply(names(universes), function(u) {
  rows <- list(decomp_one(universes[[u]], "ALL", "global", universes[[u]]))   # global
  rows <- c(rows, lapply(seq_len(nrow(roster)), function(i)
    decomp_one(prog_ens[[roster$program[i]]], roster$program[i], roster$arm[i], universes[[u]])))
  dplyr::bind_rows(rows) |> dplyr::mutate(universe = u, .before = 1)
}))

# per-gene aligned frame (divergent set) for the 2D map, tagged by first matching roster arm
div <- universes$divergent
div <- div[!is.na(m6[div]) & m6[div] != 0 & !is.na(tn[div]) & !is.na(tp[div])]
tag_of <- rep("other", length(div)); names(tag_of) <- div
for (i in rev(seq_len(nrow(roster)))) {         # rev so earlier (mito) arms win ties
  hit <- intersect(div, prog_ens[[roster$program[i]]]); tag_of[hit] <- roster$arm[i]
}
d_div <- sign(m6[div])
conv_fade_genes <- tibble::tibble(
  gene = div, arm = unname(tag_of[div]),
  wt_conv  = d_div * tn[div],                   # aligned WT temporal
  myc_fade = d_div * tp[div],                   # aligned Myc+ temporal
  aligned_gap6 = d_div * m6[div])

# =============================================================================
# PART B: BROADENED-TF ABSORPTION (extends Issue #5; BOUNDED, per-sample GSVA)
# =============================================================================
dds <- readRDS(here::here("results", "dds_int_run.rds"))
sm  <- as.data.frame(SummarizedExperiment::colData(dds))
sm$timepoint  <- stats::relevel(as.factor(sm$timepoint),  "6W")
sm$myc_status <- stats::relevel(as.factor(sm$myc_status), "neg")
sm$group      <- factor(sm$group, levels = c("6W_neg", "6W_pos", "12W_neg", "12W_pos"))
samples  <- colnames(dds)
gsva_out <- readRDS(here::here("results", "gsva_scores.rds"))
scores   <- gsva_out$scores[, samples, drop = FALSE]
vst_mat  <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))[, samples, drop = FALSE]

zrow      <- function(m) t(scale(t(m)))
composite <- function(m, genes) colMeans(zrow(m[intersect(genes, rownames(m)), , drop = FALSE]))
comp_gsva <- function(sets) {
  s <- intersect(sets, rownames(scores)); stopifnot(length(s) > 0)
  colMeans(scores[s, , drop = FALSE])
}
grep_sets <- function(pattern, extra = character(0))
  unique(c(grep(pattern, rownames(scores), value = TRUE), intersect(extra, rownames(scores))))

# data-driven mito-TF panel: shortlist LE-context TFs (strong mito-target enrichment), EXCLUDING
# MYC (self-defeating for a MYC-independent test) and the E2F/biogenesis/ESR1/lipogenic TFs
# (already their own axes), mapped to available TFT_*_GRAY_*_MITO GSVA sets.
sl <- utils::read.csv(here::here("data", "genesets_from_library", "gray_chea_mito_tf_shortlist.csv"),
                      stringsAsFactors = FALSE)
excl_tf <- c("MYC", "E2F1", "E2F4", "E2F6", "E2F7", "E2F8",
             "ESRRA", "GABPA", "NRF1", "ESR1", "SREBF1", "SREBF2", "MLX")
sl_le <- sl[grepl("LE", sl$context) & sl$bh_fdr < 1e-10 & !(sl$TF %in% excl_tf), ]
panel_tf <- utils::head(unique(sl_le$TF[order(sl_le$bh_fdr)]), 20)
mito_panel_sets <- unique(unlist(lapply(panel_tf, function(tf)
  grep(paste0("^TFT_", tf, "_GRAY_.*_MITO$"), rownames(scores), value = TRUE))))

tf_axes <- list(
  tf_biogenesis = grep_sets("^(ESRRA|GABPA|NRF1)_MITO$"),
  tf_e2f        = grep_sets("^TFT_E2F[0-9]_GRAY_.*_MITO$", "PROLIF_E2F_HALLMARK"),
  tf_esr1       = grep_sets("^TFT_ESR1_GRAY_", "ESR1_AND_CORE_MITO"),
  tf_lipogenic  = grep_sets("^TFT_(SREBF1|SREBF2|MLX)_GRAY_.*_MITO$"),
  tf_mito_panel = mito_panel_sets)
stopifnot(all(vapply(tf_axes, length, integer(1)) > 0))
message(sprintf("TF axes (n sets): %s",
                paste(sprintf("%s=%d", names(tf_axes), lengths(tf_axes)), collapse = " ")))

# per-sample moderation frame (outcomes + scaled TF axes)
fr <- data.frame(
  sample = samples, group = sm$group, tp = sm$timepoint, myc = sm$myc_status,
  oxphos_abs  = composite(vst_mat, defs$nuclear_oxphos_ens),
  oxphos_gsva = comp_gsva(defs$mito_ox_gsva_sets),
  stringsAsFactors = FALSE)
for (ax in names(tf_axes)) fr[[ax]] <- as.numeric(scale(comp_gsva(tf_axes[[ax]])))
stopifnot(!anyNA(fr[, c("oxphos_abs", "oxphos_gsva", names(tf_axes))]))

outcomes <- c("oxphos_abs", "oxphos_gsva")
INT <- "tp12W:mycpos"
adj_formula <- function(ax) stats::as.formula(paste("y ~ tp*myc +", ax, "+ tp:", ax))
b_int_of <- function(fit) { co <- stats::coef(fit); stopifnot(INT %in% names(co)); unname(co[INT]) }

absorb_one <- function(oc, ax) {
  d    <- data.frame(y = fr[[oc]], tp = fr$tp, myc = fr$myc, ax = fr[[ax]])
  names(d)[4] <- ax
  base <- stats::lm(y ~ tp * myc, d); adj <- stats::lm(adj_formula(ax), d)
  b0 <- b_int_of(base); b1 <- b_int_of(adj); av <- stats::anova(base, adj)
  tibble::tibble(
    outcome = oc, axis = ax, b_int_base = b0, b_int_adj = b1,
    delta_b_int = b1 - b0, absorption_frac = 1 - b1 / b0,
    int_survives_p = summary(adj)$coefficients[INT, "Pr(>|t|)"],
    delta_r2 = summary(adj)$r.squared - summary(base)$r.squared,
    lrt_F = av$F[2], lrt_p = av$`Pr(>F)`[2])
}
tf_absorption <- dplyr::bind_rows(
  lapply(outcomes, function(oc) dplyr::bind_rows(lapply(names(tf_axes), function(ax) absorb_one(oc, ax)))))

# stratified case bootstrap CI on delta_b_int + absorption_frac (reuse Issue #5 approach)
grid <- expand.grid(outcome = outcomes, axis = names(tf_axes), stringsAsFactors = FALSE)
metric_vec <- function(d) unlist(lapply(seq_len(nrow(grid)), function(i) {
  oc <- grid$outcome[i]; ax <- grid$axis[i]
  dd <- data.frame(y = d[[oc]], tp = d$tp, myc = d$myc, ax = d[[ax]]); names(dd)[4] <- ax
  b0 <- b_int_of(stats::lm(y ~ tp * myc, dd)); b1 <- b_int_of(stats::lm(adj_formula(ax), dd))
  stats::setNames(c(b1 - b0, 1 - b1 / b0), paste(oc, ax, c("delta", "frac"), sep = "__"))
}))
set.seed(1)
boot_out <- boot::boot(fr, function(data, idx) metric_vec(data[idx, ]), R = 2000, strata = fr$group)
tf_boot_ci <- tibble::tibble(
  key = names(boot_out$t0), t0 = unname(boot_out$t0),
  lo = vapply(seq_along(boot_out$t0), function(j) stats::quantile(boot_out$t[, j], 0.025, na.rm = TRUE, names = FALSE), numeric(1)),
  hi = vapply(seq_along(boot_out$t0), function(j) stats::quantile(boot_out$t[, j], 0.975, na.rm = TRUE, names = FALSE), numeric(1))) |>
  tidyr::separate(key, into = c("outcome", "axis", "metric"), sep = "__")

# bridge to deferred Part C: does each TF axis co-decline with the OXPHOS fade? Report the
# within-Myc+ and within-WT 6W->12W slope of every axis + the OXPHOS outcomes as reference.
temporal_one <- function(v, label) {
  d  <- data.frame(y = v, tp = sm$timepoint, myc = sm$myc_status)
  mp <- summary(stats::lm(y ~ tp, subset(d, myc == "pos")))$coefficients["tp12W", ]
  wt <- summary(stats::lm(y ~ tp, subset(d, myc == "neg")))$coefficients["tp12W", ]
  tibble::tibble(series = label,
    mycpos_temporal = unname(mp["Estimate"]), mycpos_p = unname(mp["Pr(>|t|)"]),
    wt_temporal = unname(wt["Estimate"]), wt_p = unname(wt["Pr(>|t|)"]))
}
tf_temporal <- dplyr::bind_rows(
  temporal_one(scale(fr$oxphos_abs)[, 1],  "oxphos_abs (fade reference)"),
  temporal_one(scale(fr$oxphos_gsva)[, 1], "oxphos_gsva (fade reference)"),
  dplyr::bind_rows(lapply(names(tf_axes), function(ax) temporal_one(fr[[ax]], ax))))

# =============================================================================
# PART C: DEFERRED -- documented only
# =============================================================================
settle_it <- tibble::tribble(
  ~step,                              ~what_it_would_show,
  "Reference deconvolution (bulk)",   "estimate proliferative/TEB/progenitor FRACTION per sample from a mouse-mammary sc atlas (e.g. Bach 2017 GSE106273) via MuSiC/Bisque; test whether the Myc-fade tracks the compartment shrinking. NOT on disk (no sc ref / deconv pkg) -> next script.",
  "Single-cell / snRNA-seq (4 grps)", "separate a composition shift from a per-cell program change directly -- the definitive test bulk cannot substitute for.",
  "Perturb E2F / candidate TFs",      "genetic necessity of the proliferation/TF axis for the OXPHOS level -- turns the co-decline into causation.")

# =============================================================================
# PART D: FIGURES
# =============================================================================
out_dir <- here::here("outputs", "attenuation_mechanism")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
arm_pal <- c("global" = "grey40", "OXPHOS core" = "#D73027", "TCA" = "#FDAE61",
             "nucleotide" = "#F46D43", "biosynthetic" = "#1A9850",
             "biogenesis/translation" = "#7B3294", "MYC-target core" = "#4575B4",
             "proliferation" = "#E7298A", "mammary-dev" = "#00A0B0", "other" = "grey75")

# A1 -- per-program attenuation split into WT-conv + Myc-fade (divergent set)
a1_df <- decomp_conv |> dplyr::filter(universe == "divergent") |>
  dplyr::select(program, arm, wt_conv, myc_fade) |>
  tidyr::pivot_longer(c(wt_conv, myc_fade), names_to = "mechanism", values_to = "value") |>
  dplyr::mutate(mechanism = factor(ifelse(mechanism == "wt_conv", "WT-convergence (background)",
                                          "Myc-fade (retreat)"),
                                   levels = c("WT-convergence (background)", "Myc-fade (retreat)")))
ord <- decomp_conv |> dplyr::filter(universe == "divergent") |> dplyr::arrange(conv_pct) |> dplyr::pull(program)
p_a1 <- a1_df |> dplyr::mutate(program = factor(program, levels = ord)) |>
  ggplot2::ggplot(ggplot2::aes(x = value, y = program, fill = mechanism)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
  ggplot2::scale_fill_manual(values = c("WT-convergence (background)" = "#00A0B0",
                                        "Myc-fade (retreat)" = "#D73027")) +
  ggplot2::labs(title = "Part A: the attenuation decomposed -- WT-convergence vs Myc-fade (aligned to Myc direction)",
                subtitle = "6W-divergent genes. WT-conv>0 = wild-type matures toward Myc (moving background); Myc-fade<0 = oncogenic program retreats",
                x = "aligned mean LFC 6W->12W", y = NULL, fill = NULL) +
  ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "A1_convergence_fade_by_program.pdf"), p_a1, width = 9, height = 5)

# A2 -- 2D gene map: aligned WT temporal vs aligned Myc+ temporal, coloured by arm
p_a2 <- conv_fade_genes |>
  dplyr::mutate(arm = factor(arm, levels = names(arm_pal))) |>
  ggplot2::ggplot(ggplot2::aes(x = wt_conv, y = myc_fade)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
  ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
  ggplot2::geom_point(data = ~dplyr::filter(.x, arm == "other"),
                      colour = "grey80", size = 0.5, alpha = 0.3) +
  ggplot2::geom_point(data = ~dplyr::filter(.x, arm != "other"),
                      ggplot2::aes(colour = arm), size = 1.1, alpha = 0.8) +
  ggplot2::scale_colour_manual(values = arm_pal, drop = TRUE, name = NULL) +
  ggplot2::labs(title = "Part A: convergence-fade gene map (aligned to Myc induction direction)",
                subtitle = "x>0 = WT moves toward Myc (convergence); y<0 = Myc+ retreats (fade). Biosynthetic converges; OXPHOS/MYC-core fade",
                x = "WT-convergence  (aligned WT 6W->12W)", y = "Myc-fade  (aligned Myc+ 6W->12W)") +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "A2_conv_fade_gene_map.pdf"), p_a2, width = 8, height = 6.5)

# A3 -- convergence vs fade SHARE by program, both universes (robustness)
p_a3 <- decomp_conv |> dplyr::filter(program != "ALL") |>
  ggplot2::ggplot(ggplot2::aes(x = conv_pct, y = stats::reorder(program, conv_pct),
                               colour = universe, shape = universe)) +
  ggplot2::geom_vline(xintercept = c(0, 50, 100), linetype = "dotted", colour = "grey70") +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::facet_grid(arm ~ ., scales = "free_y", space = "free_y") +
  ggplot2::labs(title = "Part A: WT-convergence share of the attenuation, by program (both universes)",
                subtitle = "conv% = 100 x WT-conv / attenuation. >50 = convergence-led; <0 = WT diverges (fade-only). Robust across selection",
                x = "convergence share (%)", y = NULL) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0), legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "A3_conv_fade_shares.pdf"), p_a3, width = 8.5, height = 6)

# B1 -- TF absorption: b_int base->adj per axis (does any absorb?)
p_b1 <- tf_absorption |>
  dplyr::mutate(axis = factor(axis, levels = names(tf_axes))) |>
  ggplot2::ggplot() +
  ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
  ggplot2::geom_segment(ggplot2::aes(x = b_int_base, xend = b_int_adj, y = axis, yend = axis),
                        colour = "grey70", arrow = grid::arrow(length = grid::unit(0.10, "cm"))) +
  ggplot2::geom_point(ggplot2::aes(x = b_int_base, y = axis, shape = "baseline"), size = 2.6, colour = "#4575B4") +
  ggplot2::geom_point(ggplot2::aes(x = b_int_adj, y = axis, shape = "adjusted"), size = 2.6, colour = "#D73027") +
  ggplot2::facet_wrap(~ outcome, scales = "free_x") +
  ggplot2::scale_shape_manual(values = c(baseline = 16, adjusted = 17)) +
  ggplot2::labs(title = "Part B: does any broadened TF axis absorb the attenuation?",
                subtitle = "b_int (genotype x time) before -> after adjusting for TF activity. Toward 0 = absorbs (expect: none, Issue #5 ceiling)",
                x = "interaction coefficient b_int", y = NULL, shape = NULL) +
  ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "B1_tf_bint_before_after.pdf"), p_b1, width = 10, height = 4.5)

# B2 -- delta_b_int with bootstrap CI per TF axis (stable metric)
p_b2 <- tf_boot_ci |> dplyr::filter(metric == "delta") |>
  dplyr::mutate(axis = factor(axis, levels = names(tf_axes))) |>
  ggplot2::ggplot(ggplot2::aes(x = t0, y = axis)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi), orientation = "y", width = 0.25, colour = "grey55") +
  ggplot2::geom_point(size = 2.6, colour = "#1A9850") +
  ggplot2::facet_wrap(~ outcome, scales = "free_x") +
  ggplot2::labs(title = "Part B: change in the attenuation, Delta b_int, per TF axis (bootstrap 95% CI)",
                subtitle = "stable metric (no division). CI crossing 0 = the TF axis does not reliably absorb the attenuation",
                x = "Delta b_int (adjusted - baseline)", y = NULL) +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "B2_tf_delta_bint_bootCI.pdf"), p_b2, width = 10, height = 4.5)

# B3 -- compositional diagnostic: each TF axis's Myc+ 6W->12W slope vs the OXPHOS fade
ox_fade_ref <- tf_temporal$mycpos_temporal[tf_temporal$series == "oxphos_abs (fade reference)"]
p_b3 <- tf_temporal |> dplyr::filter(grepl("^tf_", series)) |>
  ggplot2::ggplot(ggplot2::aes(x = mycpos_temporal, y = stats::reorder(series, mycpos_temporal))) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
  ggplot2::geom_vline(xintercept = ox_fade_ref, linetype = "dashed", colour = "#D73027") +
  ggplot2::geom_col(fill = "#7B3294", width = 0.6) +
  ggplot2::annotate("text", x = ox_fade_ref, y = 0.6, label = "OXPHOS fade", colour = "#D73027",
                    size = 3, hjust = -0.05) +
  ggplot2::labs(title = "Part B bridge (compositional diagnostic): TF axis Myc+ 6W->12W slope vs the OXPHOS fade",
                subtitle = "a TF axis that co-declines with OXPHOS (near the dashed line) SUPPORTS -- does not prove -- proliferative-compartment dilution",
                x = "Myc+ temporal slope (6W->12W, scaled)", y = NULL) +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "B3_tf_temporal_vs_oxphos_fade.pdf"), p_b3, width = 8.5, height = 4)

# =============================================================================
# PART E: SAVE
# =============================================================================
mechanism <- list(
  decomp_conv     = decomp_conv,
  conv_fade_genes = conv_fade_genes,
  tf_absorption   = tf_absorption,
  tf_boot_ci      = tf_boot_ci,
  tf_temporal     = tf_temporal,
  settle_it       = settle_it,
  defs = list(roster = roster, tf_axes = tf_axes, mito_panel_tf = panel_tf,
              universes = lapply(universes, length), identity_cor = ident),
  notes = paste(
    "Issue #6: a citable MECHANISM for the attenuation via the exact DESeq2 identity",
    "myc_12W-myc_6W == timepoint_pos-timepoint_neg (cor 1.0). Aligned to d=sign(myc_6W), the",
    "attenuation = WT-convergence (d*timepoint_neg, WT matures toward Myc) - Myc-fade",
    "(d*timepoint_pos, oncogenic retreat). PART A (POWERED, decomp_conv): global ~66% fade /",
    "34% convergence; PROGRAM-SPECIFIC -- biosynthetic arm converges (WT rises toward Myc),",
    "OXPHOS + MYC-target core do NOT (WT diverges; gap closes by fade), MYC-core also LEAST",
    "attenuated (protected). Powered gene-level resolution of H1(fade) vs H3(WT catch-up),",
    "robust across the divergent and effect(|LFC6|>0.5) universes. PART B (tf_absorption/",
    "tf_boot_ci, BOUNDED as Issue #5): broadened TF panel (biogenesis/E2F/ESR1/lipogenic/",
    "data-driven mito) -- delta_b_int primary, absorption_frac caveated; expect none absorb.",
    "tf_temporal = compositional DIAGNOSTIC: each TF axis Myc+ 6W->12W slope vs the OXPHOS",
    "fade -- an axis co-declining with OXPHOS SUPPORTS (not proves) proliferative-compartment",
    "dilution. PART C (deconvolution/single-cell, settle_it) is the definitive test, DEFERRED",
    "(no sc reference / deconv pkg on disk). CEILING: Part A descriptive of transcriptional",
    "change (composition still confounds the fade); Part B endogenous TFs + n=6/group ->",
    "bounds. See docs/2026-07-08_BlockA_revision_plan.md."))
saveRDS(mechanism, here::here("results", "attenuation_mechanism.rds"))
message("Saved results/attenuation_mechanism.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  am <- readRDS(here::here("results", "attenuation_mechanism.rds"))

  # --- Part A: the mechanism split per program (divergent set). conv%>50 = convergence-led;
  #     conv%<0 = WT diverges (fade-only). OXPHOS/MYC-core should be fade; biosynthetic convergence.
  am$decomp_conv |> dplyr::filter(universe == "divergent") |>
    dplyr::select(program, arm, n, atten, wt_conv, myc_fade, conv_pct, fade_pct) |> print(n = Inf)

  # --- robustness: same split on the effect-based (selection-independent) universe ---
  am$decomp_conv |> dplyr::filter(universe == "effect") |>
    dplyr::select(program, arm, n, wt_conv, myc_fade, conv_pct) |> print(n = Inf)

  # --- Part B: does any broadened TF axis absorb? delta_b_int primary; lrt_p = adds variance? ---
  am$tf_absorption |>
    dplyr::select(outcome, axis, b_int_base, b_int_adj, delta_b_int, absorption_frac, lrt_p) |> print(n = Inf)
  am$tf_boot_ci |> dplyr::filter(metric == "delta") |> print(n = Inf)

  # --- Part B bridge: which TF axis co-declines with the OXPHOS fade (compositional hint)? ---
  am$tf_temporal |> print(n = Inf)

  am$settle_it |> print()
  list.files(here::here("outputs", "attenuation_mechanism"), pattern = "\\.pdf$")
}
