# =============================================================================
# 39_lesner_mechanism_test.R
# -----------------------------------------------------------------------------
# TESTING LESNER et al. (bioRxiv 2026.07.13.738248) IN THE MMTV-Myc MAMMARY MODEL.
# Lesner report, in MYC-driven HCC, a hypoxia-linked TRANSCRIPTIONAL DOWNREGULATION
# of MYC activity plus MYC-driven mitochondrial turnover (mitophagy). Two proposed
# hypoxia mechanisms:
#   (M1) HIF-1a -> MXI1; MXI1 binds MAX -> inhibits MYC transcriptional activity.
#   (M2) HIF-2 (EPAS1) -> inhibits JUN -> reduces induction of MYC and TFAM.
# Plus mitophagy/autophagy activation (DRP1/Dnm1l, FUNDC1, ULK1, fission/fusion).
# See docs/2026-07-20_lesner_mtDNA_external_validation.md and
# docs/2026-07-19_oxphos_axis_biology_and_mtdna_priming.md.
#
# INSTRUMENT (load-bearing): the CLEAN, powered readout here is the DESeq2 CONTRAST
# (Wald stat / raw LFC), NOT per-sample GSVA couplings (the correlation ceiling makes
# n=24 per-sample correlations unreliable -- scripts 35/36). BATCH = TIMEPOINT, so the
# pure between-timepoint contrasts (timepoint_pos / timepoint_neg) are batch-confounded;
# the INTERACTION term (Myc-specific time change) is clean. Claims lean on the
# interaction and on the (clean) genotype contrasts where possible.
#
# Contrasts (results/interaction_results.rds):
#   myc_6W / myc_12W        = Myc-vs-WT genotype gap at each timepoint (CLEAN)
#   timepoint_pos           = 12W vs 6W within Myc+  (time change; batch-confounded)
#   timepoint_neg           = 12W vs 6W within WT    (WT timecourse; batch-confounded)
#   interaction             = Myc-specific time change (CLEAN, batch-orthogonal)
#
# Input:  results/interaction_results.rds        (per-contrast DESeq2 results)
#         results/combined_df_annotated_raw.rds  (gene<->symbol map, baseMean, raw LFCs)
#         results/gene_sets_list.rds             (Hallmark + MitoCarta MC_ sets)
#         results/gsva_scores.rds                (MITOCARTA_* + TFT_JUN sets, expr_mat)
# Output: results/lesner_mechanism_test.rds
#         outputs/lesner/*.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

out_dir <- here::here("outputs", "lesner")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# =============================================================================
# PART 1: LOAD + symbol map + contrast pullers
# =============================================================================
ir  <- readRDS(here::here("results", "interaction_results.rds"))
cd  <- readRDS(here::here("results", "combined_df_annotated_raw.rds"))
gl  <- readRDS(here::here("results", "gene_sets_list.rds"))
gsv <- readRDS(here::here("results", "gsva_scores.rds"))

map <- cd |> dplyr::select(gene, mgi_symbol, baseMean) |> dplyr::distinct()
cons <- c(myc_6W = "myc_6W_raw", myc_12W = "myc_12W_raw",
          timepoint_pos = "timepoint_pos_raw", timepoint_neg = "timepoint_neg_raw",
          interaction = "interaction_raw")

# per-contrast tidy table (symbol, lfc, stat, padj)
pull_contrast <- function(nm) {
  r <- as.data.frame(ir[[nm]]); r$gene <- rownames(r)
  dplyr::left_join(r, map, by = "gene") |>
    dplyr::transmute(mgi_symbol, lfc = log2FoldChange, stat, padj)
}
con_tabs <- lapply(cons, pull_contrast)

# =============================================================================
# PART 2: GENE-LEVEL TEST of the Lesner mechanism genes (LFC per contrast)
# =============================================================================
gene_groups <- list(
  hypoxia_HIF    = c("Hif1a", "Epas1", "Arnt", "Vhl"),
  myc_repressor  = c("Mxi1", "Max", "Mnt", "Mga"),
  AP1_JUN        = c("Jun", "Junb", "Jund", "Fos"),
  biogenesis_TF  = c("Tfam"),
  fission        = c("Dnm1l", "Fis1", "Mff"),
  fusion         = c("Mfn1", "Mfn2", "Opa1"),
  mitophagy      = c("Fundc1", "Bnip3", "Bnip3l", "Pink1", "Prkn", "Ulk1"),
  autophagy      = c("Sqstm1", "Gabarapl1", "Becn1", "Atg7"),
  hypoxia_target = c("Vegfa", "Slc2a1", "Pgk1", "Ldha"))
mech_map <- tibble::tibble(category = rep(names(gene_groups), lengths(gene_groups)),
                           gene     = unlist(gene_groups, use.names = FALSE))

gene_contrast_tbl <- purrr::imap_dfr(con_tabs, function(tb, cn) {
  tb |> dplyr::filter(mgi_symbol %in% mech_map$gene) |>
    dplyr::transmute(gene = mgi_symbol, contrast = cn, lfc, stat, padj)
}) |>
  dplyr::left_join(mech_map, by = "gene") |>
  dplyr::mutate(contrast = factor(contrast, levels = names(cons)),
                category = factor(category, levels = names(gene_groups)))

# =============================================================================
# PART 3: PATHWAY fGSEA across the five contrasts (the powered readout)
# =============================================================================
pick <- function(nm, src) if (nm %in% names(src)) src[[nm]] else NULL
focus_sets <- list(
  HYPOXIA         = pick("MSigDB_HALLMARK_HYPOXIA", gl),
  GLYCOLYSIS      = pick("MSigDB_HALLMARK_GLYCOLYSIS", gl),
  MYC_TARGETS     = pick("MSigDB_HALLMARK_MYC_TARGETS_V1", gl),
  OXPHOS_hallmark = pick("MSigDB_HALLMARK_OXIDATIVE_PHOSPHORYLATION", gl),
  OXPHOS_subunits = pick("MC_OXPHOS > OXPHOS subunits", gl),
  MITOPHAGY       = pick("MITOCARTA_MITOPHAGY", gsv$pathways),
  AUTOPHAGY       = pick("MITOCARTA_AUTOPHAGY", gsv$pathways),
  FISSION         = pick("MITOCARTA_FISSION", gsv$pathways),
  FUSION          = pick("MITOCARTA_FUSION", gsv$pathways),
  DYNAMICS        = pick("MITOCARTA_MITOCHONDRIAL_DYNAMICS_AND_SURVEILLANCE", gsv$pathways),
  JUN_targets     = unique(unlist(gsv$pathways[grep("^TFT_JUN_GRAY", names(gsv$pathways))])))
focus_sets <- focus_sets[!vapply(focus_sets, is.null, logical(1))]

rank_contrast <- function(nm) {                       # symbol-collapsed Wald-stat ranks
  r <- as.data.frame(ir[[nm]]); r$gene <- rownames(r)
  r <- dplyr::left_join(r, map, by = "gene") |>
    dplyr::filter(!is.na(mgi_symbol), is.finite(stat)) |>
    dplyr::group_by(mgi_symbol) |>
    dplyr::summarise(stat = stat[which.max(abs(stat))], .groups = "drop")
  sort(stats::setNames(r$stat, r$mgi_symbol), decreasing = TRUE)
}
fgsea_contrasts <- purrr::imap_dfr(cons, function(raw_nm, cn) {
  set.seed(1)
  fg <- suppressWarnings(fgsea::fgsea(focus_sets, rank_contrast(raw_nm),
                                      minSize = 3, maxSize = 800, eps = 0))
  tibble::as_tibble(fg[, c("pathway", "NES", "padj", "size")]) |>
    dplyr::mutate(contrast = cn)
})

# =============================================================================
# PART 4: IS THE MITOPHAGY/DYNAMICS ENRICHMENT MYC-SPECIFIC, OR JUST MITO CONTENT?
# =============================================================================
# Myc raises MitoCarta ~20% wholesale (script 32), so any MitoCarta subset rises with
# it (the death-priming trap, script 34). Test each mito-dynamics set's genotype (myc_6W)
# effect against an EXPRESSION-MATCHED background drawn from the MitoCarta universe.
# NOTE (script 34 rail): independent-draw null is ANTI-CONSERVATIVE for co-regulated
# sets, so a NON-significant result is safe; a significant one is not safe on p alone.
mito_universe <- unique(unlist(gsv$pathways[grep("^MITOCARTA_", names(gsv$pathways))]))
# one value per symbol (max|LFC|, matching the fGSEA collapse) + one baseMean per symbol
myc6_sym <- con_tabs$myc_6W |> dplyr::filter(!is.na(mgi_symbol), is.finite(lfc)) |>
  dplyr::group_by(mgi_symbol) |>
  dplyr::summarise(lfc = lfc[which.max(abs(lfc))], .groups = "drop")
base_sym <- map |> dplyr::filter(!is.na(mgi_symbol), is.finite(baseMean)) |>
  dplyr::group_by(mgi_symbol) |> dplyr::summarise(baseMean = max(baseMean), .groups = "drop")
uni <- dplyr::inner_join(myc6_sym, base_sym, by = "mgi_symbol") |>
  dplyr::filter(mgi_symbol %in% mito_universe)
uni$bin <- dplyr::ntile(uni$baseMean, 10)

matched_null <- function(target_syms, nperm = 5000) {
  tgt <- uni |> dplyr::filter(mgi_symbol %in% target_syms)
  if (nrow(tgt) < 3) return(tibble::tibble(n = nrow(tgt), obs = NA, null_mean = NA, z = NA, p = NA))
  obs  <- mean(tgt$lfc)
  null <- replicate(nperm, mean(vapply(tgt$bin, function(b) {
    pool <- uni$lfc[uni$bin == b]; pool[sample.int(length(pool), 1L)]
  }, numeric(1))))
  tibble::tibble(n = nrow(tgt), obs = obs, null_mean = mean(null),
                 z = (obs - mean(null)) / stats::sd(null), p = mean(null >= obs))
}
set.seed(1)
mito_specificity <- purrr::map_dfr(
  c("MITOPHAGY", "FISSION", "FUSION", "DYNAMICS", "AUTOPHAGY"),
  function(s) matched_null(focus_sets[[s]]) |> dplyr::mutate(set = s, .before = 1))

# =============================================================================
# PART 5: VERDICTS
# =============================================================================
nes <- fgsea_contrasts |> dplyr::select(pathway, contrast, NES) |>
  tidyr::pivot_wider(names_from = contrast, values_from = NES)
pad <- fgsea_contrasts |> dplyr::select(pathway, contrast, padj) |>
  tidyr::pivot_wider(names_from = contrast, values_from = padj)
gv <- function(p, cc) { v <- nes[[cc]][nes$pathway == p]; if (length(v)) v else NA_real_ }
gp <- function(p, cc) { v <- pad[[cc]][pad$pathway == p]; if (length(v)) v else NA_real_ }
gene_lfc <- function(g, cc) {
  v <- gene_contrast_tbl$lfc[gene_contrast_tbl$gene == g & gene_contrast_tbl$contrast == cc]
  if (length(v)) v[1] else NA_real_
}

hypoxia_verdict <- sprintf(paste0(
  "M1/M2 (hypoxia -> MYC-off) do NOT reproduce in the mammary model. HYPOXIA is LOWER in ",
  "Myc+ (myc_6W NES %+.2f, padj %.1g) and does NOT rise over time in Myc+ (timepoint_pos ",
  "%+.2f, padj %.1g); the positive interaction (%+.2f) is because hypoxia falls LESS in Myc+ ",
  "than WT (timepoint_neg %+.2f), not because it rises. Gene level: Hif1a myc_6W %+.2f (low in ",
  "Myc+); Mxi1 myc_6W %+.2f but flat over time (timepoint_pos %+.2f); Tfam ~flat (myc_6W %+.2f). ",
  "The one thread is Epas1/HIF-2 rising over time in Myc+ (interaction %+.2f) -- but its predicted ",
  "downstream FAILS: Jun rises over time in Myc+ (timepoint_pos %+.2f), not down, and JUN targets ",
  "are not suppressed (interaction NES %+.2f). => the HIF/MXI1/JUN axis is liver-specific, not ",
  "transferable here."),
  gv("HYPOXIA","myc_6W"), gp("HYPOXIA","myc_6W"), gv("HYPOXIA","timepoint_pos"),
  gp("HYPOXIA","timepoint_pos"), gv("HYPOXIA","interaction"), gv("HYPOXIA","timepoint_neg"),
  gene_lfc("Hif1a","myc_6W"), gene_lfc("Mxi1","myc_6W"), gene_lfc("Mxi1","timepoint_pos"),
  gene_lfc("Tfam","myc_6W"), gene_lfc("Epas1","interaction"), gene_lfc("Jun","timepoint_pos"),
  gv("JUN_targets","interaction"))

oxphos_verdict <- sprintf(paste0(
  "The nuclear-OXPHOS time-decline in Myc+ is REAL and MYC-SPECIFIC, but NOT hypoxia-linked. ",
  "Hallmark OXPHOS: Myc raises it (myc_6W %+.2f, myc_12W %+.2f) then it DECLINES over time in ",
  "Myc+ (timepoint_pos %+.2f, padj %.1g) with the WT timecourse NOT significant (timepoint_neg ",
  "%+.2f, padj %.1g) and a clean Myc-specific interaction (%+.2f, padj %.1g). It falls in lockstep ",
  "with MYC targets (interaction %+.2f) while HYPOXIA does not rise -- so it is the general Myc-",
  "PROGRAMME fade (canonical attenuation), not a hypoxia switch. (Nuance: nuclear OXPHOS subunits ",
  "decline in WT too, timepoint_neg %+.2f -- the Issue #6 convergence component; the Hallmark set ",
  "gives the cleaner Myc-specific signal.)"),
  gv("OXPHOS_hallmark","myc_6W"), gv("OXPHOS_hallmark","myc_12W"),
  gv("OXPHOS_hallmark","timepoint_pos"), gp("OXPHOS_hallmark","timepoint_pos"),
  gv("OXPHOS_hallmark","timepoint_neg"), gp("OXPHOS_hallmark","timepoint_neg"),
  gv("OXPHOS_hallmark","interaction"), gp("OXPHOS_hallmark","interaction"),
  gv("MYC_TARGETS","interaction"), gv("OXPHOS_subunits","timepoint_neg"))

mitophagy_verdict <- {
  mp <- mito_specificity |> dplyr::filter(set == "MITOPHAGY")
  dy <- mito_specificity |> dplyr::filter(set == "DYNAMICS")
  sprintf(paste0(
    "Lesner's DOWNSTREAM observation transfers: Myc ACTIVATES the mito dynamics/mitophagy ",
    "machinery at the genotype level (DYNAMICS myc_6W NES %+.2f padj %.1g; FUSION %+.2f; ",
    "MITOPHAGY %+.2f; FISSION %+.2f; AUTOPHAGY %+.2f), and it attenuates over time like the rest ",
    "of the programme (interaction DYNAMICS %+.2f). TWO CAVEATS: (a) at the gene level DRP1/Dnm1l ",
    "(myc_6W %+.2f) and Fundc1 (%+.2f) are only weakly Myc-induced, and the big mitophagy movers ",
    "Bnip3/Bnip3l/Prkn rise EQUALLY in WT (timepoint_pos ~ timepoint_neg -- shared/batch, not ",
    "Myc-specific); (b) vs an EXPRESSION-MATCHED MitoCarta background, the genotype rise is NOT ",
    "mitophagy-specific -- MITOPHAGY obs %+.2f vs mito background %+.2f (z %+.2f, p %.2f), DYNAMICS ",
    "obs %+.2f vs %+.2f (z %+.2f, p %.2f). So Myc raises mito-turnover genes largely BECAUSE it ",
    "raises mito content ~20%% wholesale (script 32/34), not a mitophagy-specific programme. The ",
    "matched null is anti-conservative, so this NON-specific verdict is safe."),
    gv("DYNAMICS","myc_6W"), gp("DYNAMICS","myc_6W"), gv("FUSION","myc_6W"),
    gv("MITOPHAGY","myc_6W"), gv("FISSION","myc_6W"), gv("AUTOPHAGY","myc_6W"),
    gv("DYNAMICS","interaction"), gene_lfc("Dnm1l","myc_6W"), gene_lfc("Fundc1","myc_6W"),
    mp$obs, mp$null_mean, mp$z, mp$p, dy$obs, dy$null_mean, dy$z, dy$p)
}

for (v in list(hypoxia_verdict, oxphos_verdict, mitophagy_verdict))
  message("\n", paste(strwrap(v, width = 92), collapse = "\n"))
message("")

# =============================================================================
# PART 6: FIGURES
# =============================================================================
# A -- gene-level LFC heatmap for the Lesner mechanism genes, grouped by category
gene_ord <- mech_map$gene[order(mech_map$category)]
p_a <- gene_contrast_tbl |>
  dplyr::mutate(gene = factor(gene, levels = rev(gene_ord))) |>
  ggplot2::ggplot(ggplot2::aes(contrast, gene, fill = lfc)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", lfc)), size = 2.6) +
  ggplot2::facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  ggplot2::scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                                midpoint = 0, limits = c(-1, 1), oob = scales::squish) +
  ggplot2::labs(
    title = "Lesner mechanism genes across the DESeq2 contrasts",
    subtitle = paste("myc_6W/12W = genotype gap; timepoint_pos/neg = time within genotype;",
                     "interaction = Myc-specific time change.\nHypoxia/HIF/MXI1 not elevated in Myc+;",
                     "mito-turnover genes weakly Myc-induced; Bnip3/Prkn rise shared with WT."),
    x = NULL, y = NULL, fill = "raw LFC") +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                 strip.text.y.left = ggplot2::element_text(angle = 0))
ggplot2::ggsave(file.path(out_dir, "A_mechanism_genes.pdf"), p_a, width = 7.5, height = 9)

# B -- fGSEA NES heatmap across contrasts (the powered readout)
set_ord <- c("HYPOXIA", "JUN_targets", "GLYCOLYSIS", "MYC_TARGETS", "OXPHOS_hallmark",
             "OXPHOS_subunits", "DYNAMICS", "FISSION", "FUSION", "MITOPHAGY", "AUTOPHAGY")
p_b <- fgsea_contrasts |>
  dplyr::mutate(pathway  = factor(pathway, levels = rev(set_ord)),
                contrast = factor(contrast, levels = names(cons)),
                sig      = ifelse(padj < 0.05, "*", "")) |>
  ggplot2::ggplot(ggplot2::aes(contrast, pathway, fill = NES)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%s", NES, sig)), size = 2.9) +
  ggplot2::scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                                midpoint = 0, limits = c(-3.5, 3.5), oob = scales::squish) +
  ggplot2::labs(
    title = "Pathway enrichment across contrasts (fGSEA, Wald-stat ranks; * padj < 0.05)",
    subtitle = paste("OXPHOS/MYC up in Myc+ then decline Myc-specifically (interaction);",
                     "HYPOXIA down in Myc+, no rise over time;\nMyc raises mito-dynamics/mitophagy",
                     "(genotype) but see the matched-null panel for specificity."),
    x = NULL, y = NULL, fill = "NES") +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
ggplot2::ggsave(file.path(out_dir, "B_fgsea_contrasts.pdf"), p_b, width = 7.5, height = 6)

# C -- is the mito-turnover genotype rise mitophagy-specific? (expression-matched null)
p_c <- mito_specificity |> dplyr::filter(!is.na(z)) |>
  dplyr::mutate(set = factor(set, levels = c("DYNAMICS","FISSION","FUSION","MITOPHAGY","AUTOPHAGY"))) |>
  ggplot2::ggplot(ggplot2::aes(z, set)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = z, yend = set), linewidth = 0.5, colour = "grey50") +
  ggplot2::geom_point(ggplot2::aes(colour = p < 0.05), size = 3.5) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("p=%.2f", p)), nudge_y = 0.28, size = 3) +
  ggplot2::scale_colour_manual(values = c("TRUE" = "#D73027", "FALSE" = "grey40"),
                               labels = c("TRUE" = "p<0.05", "FALSE" = "n.s."), name = NULL) +
  ggplot2::labs(
    title = "Is Myc's mito-turnover rise SPECIFIC, or the general mito-content effect?",
    subtitle = paste("z of the set's genotype (myc_6W) LFC vs an expression-matched MitoCarta",
                     "background.\nz<=0 => the set rises no more than an average mito gene",
                     "(content, not mitophagy-specific)."),
    x = "z vs expression-matched MitoCarta background", y = NULL) +
  ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "C_mitophagy_specificity_null.pdf"), p_c, width = 7.5, height = 4.5)

message("Figures written to ", out_dir)

# =============================================================================
# PART 7: SAVE + SANDBOX
# =============================================================================
lesner_out <- list(
  gene_contrast_tbl = gene_contrast_tbl,
  fgsea_contrasts   = fgsea_contrasts,
  nes_matrix        = nes,
  padj_matrix       = pad,
  mito_specificity  = mito_specificity,
  hypoxia_verdict   = hypoxia_verdict,
  oxphos_verdict    = oxphos_verdict,
  mitophagy_verdict = mitophagy_verdict,
  notes = paste(
    "Test of Lesner et al. (bioRxiv 2026.07.13.738248) in the MMTV-Myc mammary model",
    "(2026-07-20). INSTRUMENT = DESeq2 contrasts (clean genotype gaps + Myc-specific",
    "interaction), NOT per-sample GSVA couplings (correlation ceiling). BATCH=TIMEPOINT:",
    "timepoint_pos/neg are batch-confounded; the interaction is clean. FINDINGS: (1) the",
    "hypoxia->MXI1/JUN Myc-suppression mechanisms do NOT reproduce (hypoxia LOWER in Myc+,",
    "no rise over time; HIF1A low; MXI1 flat; TFAM flat; JUN not suppressed) -- liver-specific.",
    "(2) nuclear OXPHOS declines over time Myc-SPECIFICALLY (interaction), but NOT hypoxia-",
    "linked -- it is the general Myc-programme fade. (3) Myc raises the mito dynamics/mitophagy",
    "machinery (genotype) but NOT mitophagy-specifically -- vs an expression-matched MitoCarta",
    "background the rise is ~average for a mito gene (the script-32/34 content effect); the big",
    "mitophagy movers (Bnip3/Prkn) rise shared with WT. SCOPE: exploratory contrast layer;",
    "single-gene LFCs are illustration, set-level fGSEA + the matched null are the claims.",
    "The matched null is anti-conservative (script 34 rail) -> the NON-specific verdict is safe."))
saveRDS(lesner_out, here::here("results", "lesner_mechanism_test.rds"))
message("Saved results/lesner_mechanism_test.rds")

if (FALSE) {

  le <- readRDS(here::here("results", "lesner_mechanism_test.rds"))
  cat(strwrap(le$hypoxia_verdict,   92), sep = "\n")
  cat(strwrap(le$oxphos_verdict,    92), sep = "\n")
  cat(strwrap(le$mitophagy_verdict, 92), sep = "\n")

  le$nes_matrix        |> as.data.frame() |> print()
  le$padj_matrix       |> as.data.frame() |> print()
  le$mito_specificity  |> as.data.frame() |> print()
  le$gene_contrast_tbl |> dplyr::filter(contrast %in% c("myc_6W", "timepoint_pos", "interaction")) |>
    tidyr::pivot_wider(id_cols = c(category, gene), names_from = contrast, values_from = lfc) |>
    as.data.frame() |> print()

  list.files(here::here("outputs", "lesner"), pattern = "\\.pdf$")
}
