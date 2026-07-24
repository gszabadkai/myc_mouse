# scripts/29_attenuation_decomposition.R
# =============================================================================
# Block A revision -- Issue #4: what IS the 6W->12W attenuation? (the "NES paradox")
# =============================================================================
#
# Gate 1 (script 13) showed the genotype (Myc+/-) effect SHRINKS ~2-fold 6W->12W and
# that it is biological (effect-size, not power). "Attenuation" stayed fuzzy because two
# different rulers get read as one, and the author (Issue #4) asked to resolve it, verify
# it in ABSOLUTE mRNA, and connect it to Issue #3.
#
# THE "NES PARADOX" (resolved): the fGSEA NES of the Myc-vs-WT contrast is FLAT across age
# (MitoCarta stays top-enriched at both timepoints) yet the divergence shrinks. Not a
# paradox -- two quantities:
#   - NES is RANK-BASED / magnitude-blind: it asks "are these genes at the TOP of this
#     contrast?" -> yes at both ages. It measures identity/priority, NOT size. Using NES
#     flatness to argue "no attenuation" is a category error.
#   - What SHRINKS is MAGNITUDE: DE count 2777->239, |LFC| slope ~0.46 (Gate 1).
#   So SHRINKS = the absolute WT<->Myc+ gap; DOES NOT shrink = (i) transgene level, (ii)
#   program rank (NES). Myc points at the same programs at both ages; the distance shrinks.
#
# CORRECTION (2026-07-11, author caught it) -- the analysis is PATHWAY-RESOLVED, never an
# all-MitoCarta aggregate. The tempting aggregate "WT mito UP with age +0.82" is a MEDIAN
# over 63 heterogeneous MitoCarta sets, carried by the BIOSYNTHETIC arm (amino-acid metab
# +1.86, lipid +1.51). The OXPHOS pathways FALL in WT with age (MITOCARTA_OXPHOS
# timepoint_neg -1.94, OXPHOS_SUBUNITS -2.47) and fall FURTHER in Myc+ (-2.61). There is
# NO WT developmental OXPHOS rise; OXPHOS declines in BOTH genotypes, Myc+ faster. Any
# gap-shrink at OXPHOS is Myc+ falling toward WT, not WT rising. So the Issue #3 non-MYC
# OXPHOS axis stays GENUINELY OPEN (Part D probes it; it is not pre-answered here).
#
# mtDNA is contained: the fGSEA MITOCARTA_OXPHOS set is already nuclear (zero mt-* genes);
# mtDNA bias lives only in the GSVA MITOCARTA_ALL composite and mitoPPS (which isolates the
# 13 mt genes). Part B excludes mt-* from nuclear OXPHOS and reports mtDNA separately.
#
# METHOD VALIDITY (author's Q): NES for a magnitude claim = INVALID (rank-based) and NES on
# an all-MitoCarta aggregate hides opposite-direction pathways. mitoPPS "OXPHOS down" = a
# within-compartment REPRIORITISATION ratio -> can read "down" as reallocation even if
# absolute mRNA is flat/up. The missing check = ABSOLUTE normalized OXPHOS levels,
# pathway-resolved. That gap is the core of Issue #4 (Part B).
#
# Reframe on already-fitted data (DESeq2 contrasts, GSVA scores, mitoPPS, fGSEA); no
# DESeq/GSVA re-run. VST/log-CPM are deterministic transforms of the fitted dds.
#
# Input:  results/interaction_results.rds  (raw DESeqResults per contrast; use _raw)
#         results/dds_int_run.rds          (DESeqDataSet; normalized counts + VST source)
#         results/fgsea_percategory.rds     ($fgsea: NES per ranking x pathway)
#         results/mitopps_scores.rds        (reprioritisation OXPHOS, nuclear vs mtDNA)
#         results/gsva_scores.rds           (developmental + TF composites, Part D)
#         results/combined_df_annotated.rds (ENSMUSG <-> mgi_symbol map)
#         data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt  (per-pathway genes)
#         data/dev_mec_annotation.csv       (LASP/LHS luminal sets, Part D dev axis)
# Output: results/attenuation_decomposition.rds
#         outputs/attenuation_decomposition/*.pdf
#
# TPM SKIPPED -- for a STATISTICAL reason, not a data gap. Lengths are obtainable (biomaRt;
# counts are RSEM/tximport-style). But Part B compares each gene ACROSS the 4 groups
# (per-gene z-scores); gene length is a per-gene constant that cancels under z-scoring, so
# TPM and log-CPM give an IDENTICAL between-group picture. TPM only matters for cross-gene
# absolute-abundance comparisons, which Part B does not make.
#
# CEILING: attenuation magnitude = powered (Gate 1, 24 samples). Part D within-WT
# regressions = n=12, correlated composites -> HYPOTHESIS-GENERATING, association not
# causation; "the developmental program regulates OXPHOS" is a candidate, not proof (TF
# causality needs ATAC/ChIP). Normalization comparison is a robustness check, not new
# biology. mitoPPS = relative; absolute = levels -- report both, do not collapse.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# =============================================================================
# PART 1: LOAD + SAMPLE METADATA + IDENTIFIER MAP
# =============================================================================
dds <- readRDS(here::here("results", "dds_int_run.rds"))
sm  <- as.data.frame(SummarizedExperiment::colData(dds))
sm$timepoint  <- stats::relevel(as.factor(sm$timepoint),  "6W")
sm$myc_status <- stats::relevel(as.factor(sm$myc_status), "neg")
sm$group      <- factor(sm$group, levels = c("6W_neg", "6W_pos", "12W_neg", "12W_pos"))
samples       <- colnames(dds)                                   # master order

ir       <- readRDS(here::here("results", "interaction_results.rds"))
fg       <- readRDS(here::here("results", "fgsea_percategory.rds"))$fgsea
mp       <- readRDS(here::here("results", "mitopps_scores.rds"))
gsva_out <- readRDS(here::here("results", "gsva_scores.rds"))
scores   <- gsva_out$scores[, samples, drop = FALSE]             # align to master
set_meta <- gsva_out$set_meta
cdf      <- readRDS(here::here("results", "combined_df_annotated.rds"))

# ENSMUSG <-> mgi_symbol (sym2ens retained for ensembl -> symbol LABELLING only).
sym2ens <- cdf[!is.na(cdf$mgi_symbol) & !duplicated(cdf$mgi_symbol),
               c("mgi_symbol", "gene")]
# Vintage-aware symbol -> Ensembl for gene-set membership: a plain current-symbol
# match drops renamed genes (ATP synthase / ~12% of nuclear OXPHOS). See the helper.
source(here::here("functions", "reconcile_gene_symbols.R"))
ens_of  <- function(syms, universe = rownames(dds)) recon_to_ensembl(syms, universe)

gmt <- fgsea::gmtPathways(
  here::here("data", "genesets_from_library", "mammary_mito_myc_metab_v1_mouse.gmt"))

# per-contrast raw (unshrunken) LFC + padj vectors, indexed by Ensembl
lfc  <- lapply(c(myc_6W = "myc_6W_raw", myc_12W = "myc_12W_raw",
                 timepoint_neg = "timepoint_neg_raw", timepoint_pos = "timepoint_pos_raw"),
               function(k) stats::setNames(ir[[k]]$log2FoldChange, rownames(ir[[k]])))
padj <- lapply(c(myc_6W = "myc_6W_raw", myc_12W = "myc_12W_raw"),
               function(k) stats::setNames(ir[[k]]$padj, rownames(ir[[k]])))

# =============================================================================
# PART 2: PATHWAY ROSTER -- resolved individually, tagged by arm (never aggregated)
# =============================================================================
roster <- tibble::tribble(
  ~pathway,                            ~arm,
  "MITOCARTA_OXPHOS",                  "OXPHOS core (Issue#3 detached)",
  "MITOCARTA_OXPHOS_SUBUNITS",         "OXPHOS core (Issue#3 detached)",
  "MITOCARTA_COMPLEX_I",               "OXPHOS core (Issue#3 detached)",
  "MITOCARTA_NUCLEOTIDE_METABOLISM",   "nucleotide (Issue#3 detached)",
  "MITOCARTA_TCA_CYCLE",               "TCA (Issue#3 detached)",
  "MITOCARTA_AMINO_ACID_METABOLISM",   "biosynthetic (aggregate-inflating)",
  "MITOCARTA_LIPID_METABOLISM",        "biosynthetic (aggregate-inflating)",
  "MITOCARTA_MITOCHONDRIAL_RIBOSOME",  "biogenesis/translation (MYC-dose bystander)",
  "MYC_HALLMARK_MYC_TARGETS_V2",       "MYC-target core",
  "MYC_felsher_integrative_signature", "MYC-target core")
stopifnot(all(roster$pathway %in% names(gmt)))
path_ens <- lapply(stats::setNames(roster$pathway, roster$pathway),
                   function(p) ens_of(gmt[[p]]))

# =============================================================================
# PART A: ATTENUATION DECOMPOSITION -- pathway-resolved, two rulers side by side
# =============================================================================

# --- A1: RANK layer (fGSEA NES) -- genotype (flat/high at both ages) + temporal
nes_wide <- fg |>
  dplyr::filter(pathway %in% roster$pathway,
                ranking %in% c("myc_6W", "myc_12W", "timepoint_neg", "timepoint_pos")) |>
  dplyr::select(pathway, ranking, NES) |>
  tidyr::pivot_wider(names_from = ranking, values_from = NES) |>
  dplyr::rename(nes_geno_6W = myc_6W, nes_geno_12W = myc_12W,
                nes_wt_temporal = timepoint_neg, nes_myc_temporal = timepoint_pos)

# --- A2: MAGNITUDE layer (raw |LFC|) -- the ruler that actually attenuates
mag_one <- function(p) {
  e   <- path_ens[[p]]
  a6  <- abs(lfc$myc_6W[e]);  a12 <- abs(lfc$myc_12W[e])
  # Gate-1 style fixed set: genes 6W-divergent (padj<0.1 at 6W), scored at both ages
  div <- e[!is.na(padj$myc_6W[e]) & padj$myc_6W[e] < 0.1]
  f6  <- abs(lfc$myc_6W[div]); f12 <- abs(lfc$myc_12W[div])
  tibble::tibble(
    pathway = p, n_genes = length(e), n_div6 = length(div),
    mabs_geno_6W  = mean(a6,  na.rm = TRUE), mabs_geno_12W = mean(a12, na.rm = TRUE),
    attn_ratio_all   = mean(a12, na.rm = TRUE) / mean(a6, na.rm = TRUE),
    mabs_fixed_6W = mean(f6, na.rm = TRUE), mabs_fixed_12W = mean(f12, na.rm = TRUE),
    attn_ratio_fixed = mean(f12, na.rm = TRUE) / mean(f6, na.rm = TRUE),
    de6 = sum(padj$myc_6W[e]  < 0.1, na.rm = TRUE),
    de12 = sum(padj$myc_12W[e] < 0.1, na.rm = TRUE),
    # SIGNED temporal LFC: does the pathway rise or fall with age, per genotype?
    lfc_wt_temporal  = mean(lfc$timepoint_neg[e], na.rm = TRUE),
    lfc_myc_temporal = mean(lfc$timepoint_pos[e], na.rm = TRUE))
}
magnitude <- dplyr::bind_rows(lapply(roster$pathway, mag_one))

decomp <- roster |>
  dplyr::left_join(nes_wide, by = "pathway") |>
  dplyr::left_join(magnitude, by = "pathway") |>
  dplyr::arrange(arm, pathway)

# =============================================================================
# PART B: ABSOLUTE nuclear-OXPHOS levels, normalization-robust (the un-done check)
# =============================================================================
ox_ens_all <- ens_of(gmt[["MITOCARTA_OXPHOS_SUBUNITS"]])
mt_syms    <- mp$mtdna_genes_separated                            # 13 mtDNA-encoded
mt_ens     <- ens_of(mt_syms)
ox_ens_nuc <- setdiff(ox_ens_all, mt_ens)                         # nuclear OXPHOS subunits
message(sprintf("Nuclear OXPHOS subunits: %d genes; mtDNA-encoded: %d genes",
                length(ox_ens_nuc), length(mt_ens)))

# THREE normalizations of the same fitted dds
raw_counts <- DESeq2::counts(dds, normalized = FALSE)
norm_mats <- list(
  VST      = SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE)),
  normcnt  = log2(DESeq2::counts(dds, normalized = TRUE) + 1),
  logCPM   = log2(sweep(raw_counts, 2, colSums(raw_counts), "/") * 1e6 + 1))
norm_mats <- lapply(norm_mats, function(m) m[, samples, drop = FALSE])

# per-sample composite (mean of per-gene z-scores) for a gene set, per normalization
zrow      <- function(m) t(scale(t(m)))                           # z-score each gene (row)
composite <- function(m, genes) colMeans(zrow(m[intersect(genes, rownames(m)), , drop = FALSE]))

# lm stats reused from script 27/28 prog_stat_one, extended with per-genotype temporal
level_stat_one <- function(y, label) {
  d   <- data.frame(y = y, myc = sm$myc_status, tp = sm$timepoint, grp = sm$group)
  wsd <- sqrt(mean(tapply(d$y, d$grp, stats::var)))
  ma  <- summary(stats::lm(y ~ myc + tp, d))$coefficients["mycpos", ]
  mi  <- summary(stats::lm(y ~ tp * myc, d))$coefficients["tp12W:mycpos", ]
  wt  <- summary(stats::lm(y ~ tp, subset(d, myc == "neg")))$coefficients["tp12W", ]
  mc  <- summary(stats::lm(y ~ tp, subset(d, myc == "pos")))$coefficients["tp12W", ]
  gm  <- tapply(d$y, d$grp, mean)
  tibble::tibble(
    metric = label, within_sd = wsd,
    geno_beta = unname(ma["Estimate"]), geno_d = unname(ma["Estimate"]) / wsd,
    geno_p = unname(ma["Pr(>|t|)"]),
    int_beta = unname(mi["Estimate"]), int_p = unname(mi["Pr(>|t|)"]),
    wt_temporal_beta = unname(wt["Estimate"]), wt_temporal_p = unname(wt["Pr(>|t|)"]),
    myc_temporal_beta = unname(mc["Estimate"]), myc_temporal_p = unname(mc["Pr(>|t|)"]),
    m6_neg = unname(gm["6W_neg"]), m6_pos = unname(gm["6W_pos"]),
    m12_neg = unname(gm["12W_neg"]), m12_pos = unname(gm["12W_pos"]))
}

# composites + stats per normalization, for nuclear OXPHOS and mtDNA separately
absolute_stats <- dplyr::bind_rows(lapply(names(norm_mats), function(nm) {
  dplyr::bind_rows(
    level_stat_one(composite(norm_mats[[nm]], ox_ens_nuc), "nuclear_OXPHOS") |>
      dplyr::mutate(norm = nm, .before = 1),
    level_stat_one(composite(norm_mats[[nm]], mt_ens),     "mtDNA_OXPHOS") |>
      dplyr::mutate(norm = nm, .before = 1))
}))

# per-group per-gene means (VST) for the heatmap + a normalization-agreement check
group_gene_means <- function(m, genes) {
  g   <- intersect(genes, rownames(m))
  sub <- m[g, , drop = FALSE]
  sapply(levels(sm$group), function(gr) rowMeans(sub[, sm$group == gr, drop = FALSE]))
}
gm_by_norm <- lapply(norm_mats, group_gene_means, genes = ox_ens_nuc)
norm_agreement <- utils::combn(names(gm_by_norm), 2, simplify = FALSE) |>
  lapply(function(pr) tibble::tibble(
    norm_a = pr[1], norm_b = pr[2],
    spearman_group_means = suppressWarnings(stats::cor(
      as.vector(gm_by_norm[[pr[1]]]), as.vector(gm_by_norm[[pr[2]]]), method = "spearman")))) |>
  dplyr::bind_rows()

# =============================================================================
# PART C: mitoPPS (reprioritisation) vs ABSOLUTE (level) reconciliation
# =============================================================================
grp_of <- sm$group[match(mp$mitopps_scores$sample, samples)]     # align mitoPPS -> group
mitopps_traj <- tibble::tibble(
  pathway = rep(c("OXPHOS subunits (nuclear)", "mtDNA-encoded OXPHOS subunits"), each = 4),
  group   = factor(rep(levels(sm$group), 2), levels = levels(sm$group)),
  lens    = "mitoPPS (reprioritisation ratio)",
  value   = c(tapply(mp$mitopps_scores$`OXPHOS subunits`, grp_of, mean)[levels(sm$group)],
              tapply(mp$mitopps_scores$`mtDNA-encoded OXPHOS subunits`, grp_of, mean)[levels(sm$group)]))
absolute_traj <- tibble::tibble(
  pathway = rep(c("OXPHOS subunits (nuclear)", "mtDNA-encoded OXPHOS subunits"), each = 4),
  group   = factor(rep(levels(sm$group), 2), levels = levels(sm$group)),
  lens    = "absolute mRNA (VST composite)",
  value   = c(tapply(composite(norm_mats$VST, ox_ens_nuc), sm$group, mean)[levels(sm$group)],
              tapply(composite(norm_mats$VST, mt_ens),     sm$group, mean)[levels(sm$group)]))
# z-score each (lens x pathway) trajectory across the 4 groups so shapes overlay
recon <- dplyr::bind_rows(mitopps_traj, absolute_traj) |>
  dplyr::group_by(pathway, lens) |>
  dplyr::mutate(z = as.numeric(scale(value))) |>
  dplyr::ungroup()

# =============================================================================
# PART D: FIRST-PASS regulator regression (within WT; the OPEN end)
# =============================================================================
# Candidate MYC-INDEPENDENT drivers of the WT OXPHOS trajectory (a DECLINE, per the
# correction -- not a rise). Within WT the transgene is fixed, so this asks what tracks
# OXPHOS once Myc dose is held out. HYPOTHESIS-GENERATING (n=12, correlated composites).
comp_gsva <- function(sets) {
  s <- intersect(sets, rownames(scores)); stopifnot(length(s) > 0)
  colMeans(scores[s, , drop = FALSE])
}
# developmental maturation axis = luminal (LASP/LHS) state composite (Issue #1: declines
# in WT with age). Rebuilt from the author-curated annotation (final state coalesce).
anno <- utils::read.csv(here::here("data", "dev_mec_annotation.csv"),
                        stringsAsFactors = FALSE)
anno$fstate <- ifelse(!is.na(anno$new_state) & anno$new_state != "",
                      anno$new_state, anno$state)
luminal_sets <- intersect(anno$set[anno$fstate %in% c("LASP", "LHS")], rownames(scores))
mito_ox_sets <- grep("OXPHOS|COMPLEX_[IV]|_SUBUNITS|ASSEMBLY_FACTORS|ELECTRON_CARRIERS|CRISTAE",
                     set_meta$set_name[set_meta$category_primary == "MitoCarta"], value = TRUE)

drivers <- data.frame(
  sample       = samples,
  myc_status   = sm$myc_status, timepoint = sm$timepoint,
  dev_luminal  = comp_gsva(luminal_sets),                          # dev maturation axis
  tf_biogenesis = comp_gsva(c("ESRRA_MITO", "GABPA_MITO", "NRF1_MITO")),  # ER/PGC1a TFs
  mtdna_reprior = mp$mitopps_scores$`mtDNA-encoded OXPHOS subunits`[
                    match(samples, mp$mitopps_scores$sample)],     # mtDNA shift
  oxphos_abs   = composite(norm_mats$VST, ox_ens_nuc),             # absolute level outcome
  oxphos_gsva  = comp_gsva(mito_ox_sets))                          # GSVA enrichment outcome

wt <- drivers[drivers$myc_status == "neg", ]
fit_within_wt <- function(outcome) {
  d <- data.frame(y = scale(wt[[outcome]]),
                  dev = scale(wt$dev_luminal), tf = scale(wt$tf_biogenesis),
                  mtdna = scale(wt$mtdna_reprior))
  full <- stats::lm(y ~ dev + tf + mtdna, d)
  co   <- summary(full)$coefficients
  # zero-order + partial (within WT) for each driver
  pcor <- sapply(c("dev", "tf", "mtdna"), function(v) {
    others <- setdiff(c("dev", "tf", "mtdna"), v)
    rx <- stats::residuals(stats::lm(stats::reformulate(others, v), d))
    ry <- stats::residuals(stats::lm(stats::reformulate(others, "y"), d))
    stats::cor(rx, ry)
  })
  tibble::tibble(
    outcome = outcome, driver = c("dev_luminal", "tf_biogenesis", "mtdna_reprior"),
    zero_order_r = c(stats::cor(d$y, d$dev), stats::cor(d$y, d$tf), stats::cor(d$y, d$mtdna)),
    partial_r = unname(pcor),
    beta = co[c("dev", "tf", "mtdna"), "Estimate"],
    p = co[c("dev", "tf", "mtdna"), "Pr(>|t|)"],
    model_r2 = summary(full)$r.squared)
}
regulator_wt <- dplyr::bind_rows(fit_within_wt("oxphos_abs"), fit_within_wt("oxphos_gsva"))

# =============================================================================
# PART E: FIGURES
# =============================================================================
out_dir <- here::here("outputs", "attenuation_decomposition")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
arm_pal <- c("OXPHOS core (Issue#3 detached)"              = "#D73027",
             "nucleotide (Issue#3 detached)"               = "#F46D43",
             "TCA (Issue#3 detached)"                      = "#FDAE61",
             "biosynthetic (aggregate-inflating)"          = "#1A9850",
             "biogenesis/translation (MYC-dose bystander)" = "#7B3294",
             "MYC-target core"                             = "#4575B4")

# E1 -- the NES paradox in one figure: rank (flat/high) next to magnitude (shrinks)
e1a <- decomp |>
  ggplot2::ggplot(ggplot2::aes(x = nes_geno_6W, y = nes_geno_12W, colour = arm)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  ggplot2::geom_point(size = 3) +
  ggrepel::geom_text_repel(ggplot2::aes(label = pathway), size = 2.2, max.overlaps = 20) +
  ggplot2::scale_colour_manual(values = arm_pal, guide = "none") +
  ggplot2::labs(title = "RANK ruler (fGSEA NES): flat/high at both ages",
                subtitle = "on the diagonal = same rank enrichment 6W vs 12W (magnitude-blind)",
                x = "genotype NES @ 6W", y = "genotype NES @ 12W") +
  ggplot2::theme_bw(base_size = 9)
e1b <- decomp |>
  dplyr::mutate(pathway = stats::reorder(pathway, attn_ratio_fixed)) |>
  ggplot2::ggplot(ggplot2::aes(x = attn_ratio_fixed, y = pathway, fill = arm)) +
  ggplot2::geom_col() +
  ggplot2::geom_vline(xintercept = 1, colour = "grey40") +
  ggplot2::scale_fill_manual(values = arm_pal) +
  ggplot2::labs(title = "MAGNITUDE ruler (raw |LFC|): what actually shrinks",
                subtitle = "12W/6W |LFC| ratio on the 6W-divergent set; <1 = attenuates",
                x = "attenuation ratio (mean|LFC| 12W / 6W)", y = NULL, fill = NULL) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(legend.position = "bottom")
p_e1 <- gridExtra::arrangeGrob(e1a, e1b, ncol = 2, widths = c(1, 1.15))
ggplot2::ggsave(file.path(out_dir, "A_nes_paradox_two_rulers.pdf"), p_e1, width = 12, height = 6)

# E2 -- convergence, pathway-resolved: signed temporal LFC (WT vs Myc+). Shows OXPHOS
# falls in BOTH (no WT rise); the biosynthetic arm rises in WT (the aggregate artifact).
p_e2 <- decomp |>
  dplyr::mutate(pathway = stats::reorder(pathway, lfc_wt_temporal)) |>
  ggplot2::ggplot() +
  ggplot2::geom_segment(ggplot2::aes(x = lfc_wt_temporal, xend = lfc_myc_temporal,
                                     y = pathway, yend = pathway), colour = "grey70") +
  ggplot2::geom_point(ggplot2::aes(x = lfc_wt_temporal, y = pathway, shape = "WT 6->12"),
                      size = 3, colour = "#4575B4") +
  ggplot2::geom_point(ggplot2::aes(x = lfc_myc_temporal, y = pathway, shape = "Myc+ 6->12"),
                      size = 3, colour = "#D73027") +
  ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
  ggplot2::facet_grid(arm ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_shape_manual(values = c("WT 6->12" = 16, "Myc+ 6->12" = 17)) +
  ggplot2::labs(title = "Temporal direction, pathway-resolved (signed mean raw LFC 6W->12W)",
                subtitle = "OXPHOS/nucleotide/TCA fall in BOTH genotypes (no WT rise); biosynthetic rises in WT = the aggregate artifact",
                x = "mean raw LFC 6W->12W", y = NULL, shape = NULL) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0), legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "A_convergence_temporal_by_pathway.pdf"), p_e2, width = 9, height = 7)

# E3 -- ABSOLUTE nuclear-OXPHOS heatmap (VST z-scores x 4 group means, row-clustered)
zmat <- t(scale(t(gm_by_norm$VST)))
zmat <- zmat[stats::order.dendrogram(stats::as.dendrogram(
  stats::hclust(stats::dist(zmat)))), , drop = FALSE]
sym_lab <- sym2ens$mgi_symbol[match(rownames(zmat), sym2ens$gene)]
heat_df <- tibble::tibble(gene = factor(rownames(zmat), levels = rownames(zmat)),
                          symbol = sym_lab)[rep(1:nrow(zmat), 4), ]
heat_df$group <- factor(rep(colnames(zmat), each = nrow(zmat)), levels = levels(sm$group))
heat_df$z     <- as.vector(zmat)
p_e3 <- ggplot2::ggplot(heat_df, ggplot2::aes(x = group, y = symbol, fill = z)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
  ggplot2::labs(title = "Absolute nuclear-OXPHOS subunit mRNA (VST z-score, group means)",
                subtitle = sprintf("%d nuclear subunits; is the level reduced in absolute mRNA, or only in mitoPPS?",
                                   length(ox_ens_nuc)),
                x = NULL, y = NULL, fill = "z") +
  ggplot2::theme_bw(base_size = 6) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 4))
ggplot2::ggsave(file.path(out_dir, "B_nuclear_oxphos_heatmap_vst.pdf"), p_e3, width = 5.5, height = 9)

# E4 -- absolute composite trajectory across 4 groups, all 3 normalizations agree
comp_traj <- dplyr::bind_rows(lapply(names(norm_mats), function(nm) tibble::tibble(
  norm = nm, group = levels(sm$group),
  nuclear_OXPHOS = tapply(composite(norm_mats[[nm]], ox_ens_nuc), sm$group, mean)[levels(sm$group)],
  mtDNA_OXPHOS   = tapply(composite(norm_mats[[nm]], mt_ens),     sm$group, mean)[levels(sm$group)]))) |>
  tidyr::pivot_longer(c(nuclear_OXPHOS, mtDNA_OXPHOS), names_to = "compartment", values_to = "composite") |>
  dplyr::mutate(group = factor(group, levels = levels(sm$group)))
p_e4 <- ggplot2::ggplot(comp_traj, ggplot2::aes(x = group, y = composite, colour = norm, group = norm)) +
  ggplot2::geom_line() + ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~ compartment, scales = "free_y") +
  ggplot2::labs(title = "Absolute OXPHOS composite across groups -- 3 normalizations",
                subtitle = "normalization-robustness: shapes should agree (see norm_agreement Spearman)",
                x = NULL, y = "mean per-gene z composite", colour = "normalization") +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "B_absolute_composite_trajectory.pdf"), p_e4, width = 9, height = 4.5)

# E5 -- mitoPPS vs ABSOLUTE reconciliation (z-scored trajectories overlaid)
p_e5 <- ggplot2::ggplot(recon, ggplot2::aes(x = group, y = z, colour = lens, group = lens)) +
  ggplot2::geom_line() + ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey70", linetype = "dotted") +
  ggplot2::facet_wrap(~ pathway) +
  ggplot2::scale_colour_manual(values = c("mitoPPS (reprioritisation ratio)" = "#D73027",
                                          "absolute mRNA (VST composite)"     = "#4575B4")) +
  ggplot2::labs(title = "Part C: mitoPPS reprioritisation vs absolute mRNA (z across groups)",
                subtitle = "agree = real level change; diverge = pure reprioritisation. Is mitoPPS 'OXPHOS down' in absolute mRNA?",
                x = NULL, y = "z-scored trajectory", colour = NULL) +
  ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "C_mitopps_vs_absolute.pdf"), p_e5, width = 9, height = 5)

# E6 -- Part D within-WT regulators: partial r of each candidate to the OXPHOS trajectory
p_e6 <- regulator_wt |>
  dplyr::mutate(driver = stats::reorder(driver, partial_r)) |>
  ggplot2::ggplot(ggplot2::aes(x = partial_r, y = driver, fill = outcome)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
  ggplot2::labs(title = "Part D (within WT, n=12): candidate MYC-independent drivers of OXPHOS",
                subtitle = "partial r (holding the other two) to the WT OXPHOS trajectory. HYPOTHESIS-GENERATING, not causal.",
                x = "partial correlation to OXPHOS | other drivers", y = NULL, fill = "OXPHOS outcome") +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "D_regulator_partials_within_wt.pdf"), p_e6, width = 8, height = 4)

# =============================================================================
# PART F: SAVE
# =============================================================================
attenuation <- list(
  decomp          = decomp,
  absolute_stats  = absolute_stats,
  norm_agreement  = norm_agreement,
  recon           = recon,
  regulator_wt    = regulator_wt,
  defs = list(
    roster            = roster,
    nuclear_oxphos_ens = ox_ens_nuc,
    mtdna_oxphos_ens  = mt_ens,
    luminal_sets      = luminal_sets,
    tf_sets           = c("ESRRA_MITO", "GABPA_MITO", "NRF1_MITO"),
    mito_ox_gsva_sets = mito_ox_sets,
    normalizations    = names(norm_mats)),
  notes = paste(
    "Issue #4. The 'NES paradox' = two rulers: fGSEA NES is RANK-based (flat/high at both",
    "ages, decomp$nes_geno_*) while MAGNITUDE (raw |LFC|, decomp$attn_ratio_*) is what",
    "attenuates. PATHWAY-RESOLVED throughout (author correction 2026-07-11): the",
    "all-MitoCarta aggregate 'WT mito +0.82' is a MEDIAN artifact of averaging the",
    "biosynthetic arm (rises in WT) with OXPHOS (FALLS in WT); see decomp$lfc_wt_temporal",
    "vs lfc_myc_temporal -- OXPHOS/nucleotide/TCA fall in BOTH genotypes (Myc+ faster), so",
    "there is NO WT developmental OXPHOS rise and any OXPHOS gap-shrink is Myc+ falling",
    "toward WT. Part B: ABSOLUTE nuclear-OXPHOS subunit mRNA across 4 groups on 3",
    "normalizations (VST / normalized-counts / log-CPM; TPM skipped as redundant for a",
    "within-gene between-group z-score design -- length cancels). absolute_stats gives",
    "genotype d + WT/Myc+ temporal betas; norm_agreement = Spearman of per-gene group",
    "means across normalizations (robustness). mtDNA (13 mt-*) reported SEPARATELY. Part C:",
    "mitoPPS reprioritisation ratio vs absolute mRNA per group (recon) -- agree=real level,",
    "diverge=pure reprioritisation; resolves whether mitoPPS 'OXPHOS down' is in absolute",
    "mRNA. Part D (regulator_wt): FIRST-PASS within-WT regression of the OXPHOS trajectory",
    "(a DECLINE) on candidate MYC-independent drivers -- luminal/developmental axis,",
    "ER/PGC1a TF activity (ESRRA/NRF1/GABPA), mtDNA reprioritisation. CEILING: n=12,",
    "correlated composites -> HYPOTHESIS-GENERATING, association not causation; the Issue #3",
    "non-MYC OXPHOS axis stays OPEN. See docs/2026-07-08_BlockA_revision_plan.md."))
saveRDS(attenuation, here::here("results", "attenuation_decomposition.rds"))
message("Saved results/attenuation_decomposition.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  ad <- readRDS(here::here("results", "attenuation_decomposition.rds"))

  # --- The two rulers side by side: NES flat/high vs magnitude attenuates ---
  ad$decomp |>
    dplyr::select(arm, pathway, nes_geno_6W, nes_geno_12W,
                  mabs_geno_6W, mabs_geno_12W, attn_ratio_fixed, de6, de12) |>
    print(n = Inf)

  # --- CORRECTION made concrete: temporal direction per pathway (WT vs Myc+) ---
  # OXPHOS/nucleotide/TCA should be NEGATIVE in BOTH columns (fall in both genotypes);
  # biosynthetic (amino-acid/lipid) POSITIVE in WT = the median artifact.
  ad$decomp |>
    dplyr::select(arm, pathway, lfc_wt_temporal, lfc_myc_temporal) |>
    dplyr::arrange(lfc_wt_temporal) |> print(n = Inf)

  # --- Part B: is nuclear OXPHOS reduced in ABSOLUTE mRNA? (per normalization) ---
  ad$absolute_stats |>
    dplyr::filter(metric == "nuclear_OXPHOS") |>
    dplyr::select(norm, geno_d, geno_p, wt_temporal_beta, wt_temporal_p,
                  myc_temporal_beta, myc_temporal_p, int_p) |> print()
  ad$absolute_stats |> dplyr::filter(metric == "mtDNA_OXPHOS") |>
    dplyr::select(norm, wt_temporal_beta, myc_temporal_beta, geno_d) |> print()
  ad$norm_agreement |> print()   # should be ~1 (robust across normalizations)

  # --- Part C: does mitoPPS 'OXPHOS down' show up in absolute mRNA? ---
  ad$recon |> tidyr::pivot_wider(names_from = lens, values_from = c(value, z)) |> print(n = Inf)

  # --- Part D: what tracks the WT OXPHOS trajectory beyond the other drivers? ---
  ad$regulator_wt |> print(n = Inf)

  list.files(here::here("outputs", "attenuation_decomposition"), pattern = "\\.pdf$")
}
