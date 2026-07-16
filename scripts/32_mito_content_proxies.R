# scripts/32_mito_content_proxies.R
# =============================================================================
# Does MITOCHONDRIAL CONTENT PER CELL change? -- transcriptomic content proxies
# (Block B, author question 2026-07-16)
# =============================================================================
#
# THE QUESTION (a reviewer will ask it). The manuscript says Myc drives a "robust
# increase in mitochondrial biogenesis". Is that mitochondrial CONTENT/mass per cell,
# or only a transcriptional PROGRAM? Nothing in Block A answers this, and three
# objects that sound like they do, do not:
#
#   - reframe_supp3$abund (script 22 "AP-abund") is NOT an abundance measure. It is
#     the SD/IQR of mitoPPS Myc-effect diffs across 19 OXPHOS pathways (0.140 -> 0.100)
#     = dispersion of REALLOCATION on ratio-normalised scores. The finalisation plan
#     glosses it as "uniform biogenesis (an abundance-LIKE change)" -- "-like" is
#     load-bearing.
#   - mitoPPS is abundance-BLIND by construction (08:404 "cancels out both total
#     mitochondrial content..."). So mtnuc_index (script 24), the mitonuclear-imbalance
#     headline, is a PRIORITY RATIO: a mitochondrion with 2x of everything scores
#     identically.
#   - total_mito_score (08:1153, subtitle "reflects mito content") is plotted to a PDF
#     and never saved or tested. It is also rowSums over 142 OVERLAPPING MitoCarta
#     pathways (double-counts genes) and ~68% mtDNA-driven. PART 5 rebuilds it.
#
#   And bio_comp (19:91) -- the composite behind "the Myc footprint IS a biogenesis
#   program, rho=0.90" -- is a mean of GSVA scores: rank-based, measures PROGRAM.
#
# The one real abundance result is attenuation_decomposition$absolute_stats (29 PART B):
# Myc raises nuclear-OXPHOS transcript (d=1.27, p=0.005), does nothing to mtDNA-encoded
# (p=0.94). Never generalised beyond OXPHOS subunits. This script generalises it.
#
# THE RESULT THAT CAME OUT OF GENERALISING IT ("commissioned but unbuilt"). Myc raises
# every NUCLEAR mito arm (mass markers, mitoribosome, import, OXPHOS_NU) AND raises the
# mtDNA-handling machinery hardest of all -- nucleoid d=2.32 p=9e-6, mt-transcription
# (Tfam/Polrmt/Tfb2m/Tefm) +43% p=8e-5 -- while mtDNA-encoded OUTPUT stays flat (p=0.83).
# So the mitonuclear imbalance is NOT Myc neglecting the mtDNA arm: Myc commissions the
# machinery to replicate and transcribe mtDNA, and the mtDNA-encoded output does not
# follow. Bulk cannot say why (copy number? mt-transcription rate? turnover?) -- mtDNA
# qPCR is the one-experiment test.
#
# WHAT THIS SCRIPT DOES. Pure reframe on fitted/saved data -- no DESeq re-fit, no
# GSVA/fGSEA re-run, scripts 26-31 untouched. Computes per-sample COMPARTMENT SHARES
# (share of transcriptome held by a mito gene panel) on raw counts, fits the genotype x
# time model, gates every claim through a library-depth/complexity QC block, and
# reconciles against the abundance-blind lenses already on disk.
#
# WHAT IT CANNOT DO (three ceilings, restated in PART 6 and in the saved notes):
#   (i)   Bulk polyA with no spike-ins and no cell counts CANNOT measure per-cell content
#         in absolute terms. DESeq2 median-of-ratios normalisation removes exactly that
#         scale. Every number here is a SHARE OF TRANSCRIPTOME.
#   (ii)  Myc globally amplifies total RNA per cell, so a CONSTANT share already implies
#         more mitochondria per cell => the share effect is a LOWER BOUND on the
#         Myc-driven content increase, and bulk cannot recover the true one.
#   (iii) Transcript share is not protein and not organelle volume. Blot (TOMM20 / VDAC /
#         CS / HSP60), mtDNA qPCR, or EM SETTLES it. n = 6/group.
#
# THE CONFOUND THAT GATES THE TIME AXIS (PART 4). The mt-* share ranges 3.4%-40% across
# samples and tracks library depth (rho ~ -0.47), which is PERFECTLY CONFOUNDED with
# timepoint: 6W = 18.8-29.0M reads from animals MYCF62-65/MYBS10x; 12W = 9.9-16.1M from
# MYCF52-56, a different cohort. Within timepoint, depth is genotype-BALANCED and both
# genotypes occur inside single litters (MYCF62, MYCF52, MYCF56). So the GENOTYPE axis is
# clean and every TEMPORAL claim is cohort/depth-exposed -- including script 29's
# "mtDNA absolute z -1.14 -> +1.23", which inherits the same exposure. There is no RIN or
# batch column on disk. This script FLAGS rather than corrects: the confound is a design
# fact, not a bug to regress away.
#
# Input:  results/dds_int_run.rds             (fitted DESeq2, sample metadata)
#         results/count_matrix.rds            (raw counts, 54838 x 24, Ensembl)
#         results/combined_df_annotated.rds   (mgi_symbol <-> Ensembl map)
#         results/mitopps_scores.rds          (mtdna_genes_separated; raw_pathway_scores)
#         results/attenuation_decomposition.rds (absolute_stats -- reconciliation)
#         results/biogenesis_discrimination.rds (mtnuc_index -- reconciliation)
#         results/reframe_supp3.rds           (mtdna$trajectory -- unread raw rows)
#         data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt
# Output: results/mito_content_proxies.rds
#         outputs/mito_content_proxies/*.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

out_dir <- here::here("outputs", "mito_content_proxies")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

group_levels <- c("6W_neg", "6W_pos", "12W_neg", "12W_pos")
geno_cols    <- c(neg = "#4575B4", pos = "#D73027")

# =============================================================================
# PART 1: LOAD + SAMPLE METADATA + IDENTIFIER MAP
# =============================================================================
# Identifier idiom mirrors script 29 (29:70-96) so the two scripts resolve the same
# gene sets to the same Ensembl IDs and their numbers are directly comparable.
dds <- readRDS(here::here("results", "dds_int_run.rds"))
sm  <- as.data.frame(SummarizedExperiment::colData(dds))
sm$timepoint  <- stats::relevel(as.factor(sm$timepoint),  "6W")
sm$myc_status <- stats::relevel(as.factor(sm$myc_status), "neg")
sm$group      <- factor(sm$group, levels = group_levels)
samples       <- colnames(dds)                                   # master order

cts <- readRDS(here::here("results", "count_matrix.rds"))
cts <- cts[, samples, drop = FALSE]                              # align to master
stopifnot(identical(colnames(cts), samples))

cdf     <- readRDS(here::here("results", "combined_df_annotated.rds"))
sym2ens <- cdf[!is.na(cdf$mgi_symbol) & !duplicated(cdf$mgi_symbol),
               c("mgi_symbol", "gene")]
ens_of  <- function(syms, universe = rownames(cts)) {
  e <- sym2ens$gene[match(intersect(syms, sym2ens$mgi_symbol), sym2ens$mgi_symbol)]
  intersect(e, universe)
}

gmt <- fgsea::gmtPathways(
  here::here("data", "genesets_from_library", "mammary_mito_myc_metab_v1_mouse.gmt"))

mp  <- readRDS(here::here("results", "mitopps_scores.rds"))
ad  <- readRDS(here::here("results", "attenuation_decomposition.rds"))
bd  <- readRDS(here::here("results", "biogenesis_discrimination.rds"))
ref <- readRDS(here::here("results", "reframe_supp3.rds"))

# The 13 mtDNA-encoded genes -- taken from the mitoPPS object (the project's single
# source of truth for the split), NOT re-derived by a regex.
mt_syms <- mp$mtdna_genes_separated
mt_ens  <- ens_of(mt_syms)
message(sprintf("mtDNA-encoded genes resolved: %d of %d", length(mt_ens), length(mt_syms)))
stopifnot(length(mt_ens) >= 10)

# =============================================================================
# PART 2: PANEL ROSTER + COMPARTMENT SHARES
# =============================================================================

# --- 2a. The MASS_MARKERS panel -----------------------------------------------
# DELIBERATE, FLAGGED EXCEPTION to CLAUDE.md's "do not rebuild gene sets": this is not
# a library set being reconstructed -- the library has NO mass-marker equivalent. It is
# the in-silico stand-in for the TOMM20/VDAC/CS/HSP60 blot we do not have.
#
# Split into three arms so no single route carries the claim. Hspa9 (Grp75/mortalin) and
# Hspd1 (Hsp60) are standard mito-mass markers AND direct MYC transactivation targets --
# so the chaperone arm has a route to elevation that the structural arm does not. The
# claim-bearing composite therefore EXCLUDES the chaperones (mass_nochap); the full panel
# is corroboration. Tspo is reported individually: it is the lone non-mover.
mass_roster <- tibble::tribble(
  ~gene,      ~arm,
  "Vdac1",    "OMM_structural",
  "Vdac2",    "OMM_structural",
  "Vdac3",    "OMM_structural",
  "Tomm20",   "OMM_structural",
  "Tomm22",   "OMM_structural",
  "Tomm40",   "OMM_structural",
  "Tspo",     "OMM_structural",
  "Cs",       "matrix_IMM",
  "Immt",     "matrix_IMM",
  "Slc25a3",  "matrix_IMM",
  "Timm23",   "matrix_IMM",
  "Timm44",   "matrix_IMM",
  "Hspa9",    "chaperone",
  "Hspd1",    "chaperone")

mass_all    <- mass_roster$gene
mass_nochap <- mass_roster$gene[mass_roster$arm != "chaperone"]

# --- 2b. The panel roster ------------------------------------------------------
panel_roster <- tibble::tribble(
  ~panel,                     ~source,   ~tag,
  "MITOCARTA_ALL",            "gmt",     "compartment",
  "MITOCARTA_NUCLEAR_ENCODED","gmt",     "compartment",
  "MITOCARTA_MTDNA_ENCODED",  "gmt",     "compartment",
  "MASS_MARKERS",             "roster",  "mass proxy",
  "MASS_MARKERS_NOCHAP",      "roster",  "mass proxy (claim-bearing)",
  "MASS_OMM_structural",      "roster",  "mass proxy (arm)",
  "MASS_matrix_IMM",          "roster",  "mass proxy (arm)",
  "MASS_chaperone",           "roster",  "mass proxy (arm)",
  "MITOCARTA_MITOCHONDRIAL_RIBOSOME",   "gmt", "biogenesis machinery",
  "MITOCARTA_PROTEIN_IMPORT_AND_SORTING","gmt", "biogenesis machinery",
  "MITOCARTA_MTDNA_REPLICATION","gmt",   "mtDNA machinery",
  "MITOCARTA_MTDNA_NUCLEOID", "gmt",     "mtDNA machinery",
  "MITOCARTA_TRANSCRIPTION",  "gmt",     "mtDNA machinery",
  "MITOCARTA_FISSION",        "gmt",     "dynamics",
  "MITOCARTA_FUSION",         "gmt",     "dynamics",
  "MITOCARTA_OXPHOS_NU",      "gmt",     "OXPHOS (reconciles script 29)",
  "MITOCARTA_OXPHOS_MT",      "gmt",     "OXPHOS (reconciles script 29)")

gmt_panels <- panel_roster$panel[panel_roster$source == "gmt"]
stopifnot(all(gmt_panels %in% names(gmt)))

panel_syms <- c(
  stats::setNames(lapply(gmt_panels, function(p) gmt[[p]]), gmt_panels),
  list(MASS_MARKERS        = mass_all,
       MASS_MARKERS_NOCHAP = mass_nochap,
       MASS_OMM_structural = mass_roster$gene[mass_roster$arm == "OMM_structural"],
       MASS_matrix_IMM     = mass_roster$gene[mass_roster$arm == "matrix_IMM"],
       MASS_chaperone      = mass_roster$gene[mass_roster$arm == "chaperone"]))
panel_syms <- panel_syms[panel_roster$panel]                     # roster order
panel_ens  <- lapply(panel_syms, ens_of)

panel_roster$n_genes <- vapply(panel_ens, length, integer(1))
message("Panels resolved:")
panel_roster |> print(n = nrow(panel_roster))
stopifnot(all(panel_roster$n_genes > 0))

# --- 2c. Shares ---------------------------------------------------------------
# Computed on RAW counts. Gene length is constant across samples, so share-of-counts is
# valid BETWEEN samples for a fixed gene set (it is not valid BETWEEN sets -- never
# compare two panels' absolute share to each other, only each panel across groups).
#
# TWO denominators:
#   share_all  = share of the whole transcriptome. Interpretable, but the 13 mt-* genes
#                are 7-16% of the library and swing wildly, so they move the denominator.
#   share_nomt = denominator EXCLUDES the 13 mt-* genes. Immune to mt inflation; this is
#                the denominator every NUCLEAR panel is read on.
den_all  <- colSums(cts)
den_nomt <- colSums(cts[setdiff(rownames(cts), mt_ens), , drop = FALSE])

share_of <- function(ens, den) 100 * colSums(cts[ens, , drop = FALSE]) / den

shares <- purrr::map_dfr(names(panel_ens), function(p) {
  tibble::tibble(
    panel      = p,
    sample     = samples,
    group      = sm$group,
    timepoint  = sm$timepoint,
    myc_status = sm$myc_status,
    share_all  = share_of(panel_ens[[p]], den_all),
    share_nomt = share_of(panel_ens[[p]], den_nomt))
})

share_group_means <- shares |>
  tidyr::pivot_longer(c(share_all, share_nomt),
                      names_to = "denominator", values_to = "share") |>
  dplyr::group_by(panel, denominator, group) |>
  dplyr::summarise(mean_share = mean(share), sd_share = stats::sd(share),
                   .groups = "drop") |>
  tidyr::pivot_wider(names_from = group, values_from = c(mean_share, sd_share))

# --- 2d. Per-gene mass panel (never hide the members inside a composite) --------
mass_per_gene <- purrr::map_dfr(seq_len(nrow(mass_roster)), function(i) {
  g <- ens_of(mass_roster$gene[i])
  if (length(g) == 0) return(NULL)
  f <- share_of(g, den_nomt)
  d <- data.frame(y = log2(f), myc = sm$myc_status, tp = sm$timepoint)
  co <- summary(stats::lm(y ~ myc + tp, d))$coefficients["mycpos", ]
  gm <- tapply(f, sm$group, mean)
  tibble::tibble(gene = mass_roster$gene[i], arm = mass_roster$arm[i],
                 myc_l2fc = unname(co["Estimate"]), myc_p = unname(co["Pr(>|t|)"]),
                 m6_neg = unname(gm["6W_neg"]),  m6_pos  = unname(gm["6W_pos"]),
                 m12_neg = unname(gm["12W_neg"]), m12_pos = unname(gm["12W_pos"]))
}) |>
  dplyr::arrange(dplyr::desc(myc_l2fc))

mass_coherence <- list(
  n_genes    = nrow(mass_per_gene),
  n_myc_up   = sum(mass_per_gene$myc_l2fc > 0),
  n_myc_sig  = sum(mass_per_gene$myc_p < 0.05),
  non_movers = mass_per_gene$gene[mass_per_gene$myc_p >= 0.05])

# =============================================================================
# PART 3: GENOTYPE x TIME MODEL ON THE SHARES
# =============================================================================
# share_stat_one mirrors script 29's level_stat_one (29:186-203) column-for-column, so
# the two tables can be read side by side. Fitted on log2(share) => betas are log2 fold
# changes, geno_d is a within-group-SD-standardised effect size on that scale.
share_stat_one <- function(y, label, denominator, meta = sm) {
  d   <- data.frame(y = y, myc = meta$myc_status, tp = meta$timepoint, grp = meta$group)
  wsd <- sqrt(mean(tapply(d$y, d$grp, stats::var)))
  ma  <- summary(stats::lm(y ~ myc + tp, d))$coefficients["mycpos", ]
  mi  <- summary(stats::lm(y ~ tp * myc, d))$coefficients["tp12W:mycpos", ]
  wt  <- summary(stats::lm(y ~ tp, subset(d, myc == "neg")))$coefficients["tp12W", ]
  mc  <- summary(stats::lm(y ~ tp, subset(d, myc == "pos")))$coefficients["tp12W", ]
  gm  <- tapply(d$y, d$grp, mean)
  tibble::tibble(
    panel = label, denominator = denominator, within_sd = wsd,
    geno_beta = unname(ma["Estimate"]), geno_d = unname(ma["Estimate"]) / wsd,
    geno_p = unname(ma["Pr(>|t|)"]),
    int_beta = unname(mi["Estimate"]), int_p = unname(mi["Pr(>|t|)"]),
    wt_temporal_beta = unname(wt["Estimate"]), wt_temporal_p = unname(wt["Pr(>|t|)"]),
    myc_temporal_beta = unname(mc["Estimate"]), myc_temporal_p = unname(mc["Pr(>|t|)"]),
    m6_neg = unname(gm["6W_neg"]),   m6_pos  = unname(gm["6W_pos"]),
    m12_neg = unname(gm["12W_neg"]), m12_pos = unname(gm["12W_pos"]))
}

share_stats <- purrr::map_dfr(names(panel_ens), function(p) {
  s <- shares[shares$panel == p, ]
  dplyr::bind_rows(
    share_stat_one(log2(s$share_all),  p, "share_all"),
    share_stat_one(log2(s$share_nomt), p, "share_nomt"))
}) |>
  dplyr::left_join(panel_roster[, c("panel", "tag", "n_genes")], by = "panel")

# =============================================================================
# PART 4: THE QC GATE -- which claims survive the cohort/depth confound
# =============================================================================
# This block does NOT correct the confound. Timepoint is perfectly confounded with
# sequencing cohort, so there is no within-design contrast that separates them: any
# "correction" would be regressing the time effect on itself. It DIAGNOSES and FLAGS.

# --- 4a. Per-sample diagnostics -----------------------------------------------
nomt_rows <- setdiff(rownames(cts), mt_ens)
qc_sample <- tibble::tibble(
  sample     = samples,
  group      = sm$group,
  timepoint  = sm$timepoint,
  myc_status = sm$myc_status,
  depth_M    = den_all / 1e6,
  n_detected = colSums(cts > 0),
  # library complexity: share of NON-mt counts held by that sample's top 100 genes.
  # High = a few genes dominate = the classic low-complexity / degradation signature.
  top100_pct = apply(cts[nomt_rows, , drop = FALSE], 2,
                     function(x) 100 * sum(sort(x, decreasing = TRUE)[1:100]) / sum(x)),
  mt_pct     = share_of(mt_ens, den_all))

# --- 4b. Depth balance: clean on genotype, confounded on time ------------------
depth_balance <- qc_sample |>
  dplyr::group_by(group) |>
  dplyr::summarise(n = dplyr::n(), mean_depth_M = mean(depth_M),
                   min_depth_M = min(depth_M), max_depth_M = max(depth_M),
                   mean_mt_pct = mean(mt_pct), sd_mt_pct = stats::sd(mt_pct),
                   .groups = "drop")

depth_geno_p <- summary(stats::lm(depth_M ~ myc_status + timepoint,
                                  qc_sample))$coefficients["myc_statuspos", "Pr(>|t|)"]
depth_time_p <- summary(stats::lm(depth_M ~ myc_status + timepoint,
                                  qc_sample))$coefficients["timepoint12W", "Pr(>|t|)"]

# --- 4c. Does each panel's share track depth/complexity? ------------------------
cor_or_na <- function(a, b) if (length(a) < 3) NA_real_ else
  suppressWarnings(stats::cor(a, b, method = "spearman"))

qc_coupling <- purrr::map_dfr(names(panel_ens), function(p) {
  s <- shares[shares$panel == p, ]
  s <- s[match(qc_sample$sample, s$sample), ]
  i6 <- qc_sample$timepoint == "6W"; i12 <- qc_sample$timepoint == "12W"
  tibble::tibble(
    panel        = p,
    rho_depth    = cor_or_na(s$share_nomt, qc_sample$depth_M),
    rho_depth_6W = cor_or_na(s$share_nomt[i6],  qc_sample$depth_M[i6]),
    rho_depth_12W= cor_or_na(s$share_nomt[i12], qc_sample$depth_M[i12]),
    rho_top100   = cor_or_na(s$share_nomt, qc_sample$top100_pct))
})

# --- 4d. Leave-one-out sensitivity on the most extreme mt sample ---------------
# The confound's worst case: one sample with an mt share far outside the rest. If a
# panel's genotype effect depends on it, the panel is not reportable.
loo_sample <- qc_sample$sample[which.max(qc_sample$mt_pct)]
keep       <- samples != loo_sample
message(sprintf("Leave-one-out sample (max mt%%): %s (mt %.1f%%, depth %.1fM, %d genes)",
                loo_sample, max(qc_sample$mt_pct),
                qc_sample$depth_M[qc_sample$sample == loo_sample],
                qc_sample$n_detected[qc_sample$sample == loo_sample]))

loo_sensitivity <- purrr::map_dfr(names(panel_ens), function(p) {
  s    <- shares[shares$panel == p, ]
  full <- share_stat_one(log2(s$share_nomt), p, "share_nomt")
  drop <- share_stat_one(log2(s$share_nomt[keep]), p, "share_nomt", meta = sm[keep, ])
  tibble::tibble(
    panel = p,
    geno_beta_full = full$geno_beta, geno_p_full = full$geno_p,
    geno_beta_loo  = drop$geno_beta, geno_p_loo  = drop$geno_p,
    sign_stable    = sign(full$geno_beta) == sign(drop$geno_beta),
    sig_stable     = (full$geno_p < 0.05) == (drop$geno_p < 0.05))
})

# --- 4e. The flag every downstream statement inherits ---------------------------
# GENOTYPE contrasts: depth-balanced within timepoint, litter-controlled -> CLEAN,
#   conditional on surviving leave-one-out.
# TEMPORAL contrasts: timepoint == sequencing cohort here -> TIME-EXPOSED, always.
qc_flags <- loo_sensitivity |>
  dplyr::left_join(qc_coupling, by = "panel") |>
  dplyr::mutate(
    geno_flag = dplyr::case_when(
      !sign_stable | !sig_stable ~ "GENOTYPE: LOO-FRAGILE",
      abs(rho_depth_6W) > 0.6 | abs(rho_depth_12W) > 0.6 ~ "GENOTYPE: depth-coupled within tp",
      TRUE ~ "GENOTYPE: CLEAN"),
    time_flag = "TIME-EXPOSED (timepoint confounded with sequencing cohort)") |>
  dplyr::select(panel, geno_flag, time_flag, sign_stable, sig_stable,
                rho_depth, rho_depth_6W, rho_depth_12W)

# =============================================================================
# PART 5: RECONCILIATION + REBUILT total_mito_score
# =============================================================================

# --- 5a. Rebuild script 08's total_mito_score correctly -------------------------
# The original (08:1153-1171) is rowSums over 142 OVERLAPPING MitoCarta pathways, so a
# gene in k pathways is counted k times, and the sum is ~68% the mt-* pathway. Rebuilt:
# MITOCARTA_ALL summed ONCE, reported with and without mt-*.
total_mito_rebuilt <- dplyr::bind_rows(
  share_stat_one(log2(shares$share_all[shares$panel == "MITOCARTA_ALL"]),
                 "MITOCARTA_ALL (rebuilt total)", "share_all"),
  share_stat_one(log2(shares$share_nomt[shares$panel == "MITOCARTA_NUCLEAR_ENCODED"]),
                 "MITOCARTA_NUCLEAR_ENCODED (rebuilt total, mt-free)", "share_nomt"))

mt_dominance <- {
  a <- shares$share_all[shares$panel == "MITOCARTA_MTDNA_ENCODED"]
  b <- shares$share_all[shares$panel == "MITOCARTA_ALL"]
  tibble::tibble(mean_mtdna_share = mean(a), mean_mitocarta_share = mean(b),
                 mtdna_frac_of_mitocarta = mean(a / b))
}

# --- 5b. Against script 29's absolute_stats (VST z-composite, OXPHOS only) ------
recon_29 <- ad$absolute_stats |>
  dplyr::filter(norm == "VST") |>
  dplyr::transmute(lens = "29 absolute (VST z-composite)",
                   metric = metric, geno_d, geno_p,
                   wt_temporal_beta, myc_temporal_beta) |>
  dplyr::bind_rows(
    share_stats |>
      dplyr::filter(panel %in% c("MITOCARTA_OXPHOS_NU", "MITOCARTA_OXPHOS_MT"),
                    denominator == "share_nomt") |>
      dplyr::transmute(lens = "32 share (log2 % of non-mt transcriptome)",
                       metric = ifelse(panel == "MITOCARTA_OXPHOS_NU",
                                       "nuclear_OXPHOS", "mtDNA_OXPHOS"),
                       geno_d, geno_p, wt_temporal_beta, myc_temporal_beta))

# --- 5c. Against script 24's mtnuc_index (mitoPPS ratio) -----------------------
# Do the abundance-BLIND and abundance-CAPABLE lenses agree in SHAPE across groups?
# Agreement is the argument; divergence gets named, not smoothed.
z <- function(x) as.numeric(scale(x))
share_grp <- function(p) {
  s <- shares[shares$panel == p, ]
  vapply(group_levels, function(g) mean(s$share_nomt[s$group == g]), numeric(1))
}
mtnuc <- as.data.frame(bd$mtnuc_index)
mtnuc <- mtnuc[match(group_levels, mtnuc$group), ]
recon_24 <- tibble::tibble(
  group          = factor(group_levels, levels = group_levels),
  mitopps_mtdna  = z(mtnuc[["mtDNA-encoded"]]),
  share_mtdna    = z(share_grp("MITOCARTA_MTDNA_ENCODED")),
  mitopps_nuc    = z(mtnuc[["nuclear-encoded"]]),
  share_nuc      = z(share_grp("MITOCARTA_OXPHOS_NU")))
recon_24_agreement <- tibble::tibble(
  component = c("mtDNA-encoded", "nuclear-encoded"),
  spearman_mitopps_vs_share = c(
    cor_or_na(recon_24$mitopps_mtdna, recon_24$share_mtdna),
    cor_or_na(recon_24$mitopps_nuc,   recon_24$share_nuc)))

# --- 5d. The CHAPERONE reconciliation ------------------------------------------
# Script 24 reported "Chaperones fall -> UPR^mt REFUTED" -- a mitoPPS result (within-
# budget PRIORITISATION, over time). In absolute share the chaperones are Myc-ELEVATED
# and age-FLAT. NOT a contradiction: they deprioritise inside a reallocating compartment
# while their absolute share holds. A second worked example of the lens split, and it
# protects the script-24 result from looking overturned.
recon_chaperone <- share_stats |>
  dplyr::filter(panel %in% c("MASS_chaperone", "MASS_OMM_structural", "MASS_matrix_IMM"),
                denominator == "share_nomt") |>
  dplyr::select(panel, geno_beta, geno_p, wt_temporal_beta, wt_temporal_p,
                myc_temporal_beta, myc_temporal_p, m6_neg, m6_pos, m12_neg, m12_pos)

# --- 5e. The unread raw rows already sitting in reframe_supp3 -------------------
# Script 22 computed raw-abundance trajectories (22:96-99) and then plotted only the
# mitoPPS ones -- so these rows have never been read. Surfaced here for the record.
recon_22_raw <- ref$mtdna$trajectory |>
  dplyr::filter(grepl("raw", component)) |>
  dplyr::select(group, component, mean_score)

# =============================================================================
# PART 6: VERDICT
# =============================================================================
stat_of <- function(p, col) {
  v <- share_stats[[col]][share_stats$panel == p & share_stats$denominator == "share_nomt"]
  if (length(v) == 0) NA_real_ else v[1]
}
pct_of <- function(p) 100 * (2^stat_of(p, "geno_beta") - 1)   # log2 beta -> % change

content_verdict <- sprintf(paste0(
  "MITOCHONDRIAL CONTENT: Myc RAISES the mitochondrial share of the transcriptome -- ",
  "mass markers (chaperone-free) %+.0f%% (p=%.2g), nuclear MitoCarta %+.0f%% (p=%.2g), ",
  "mitoribosome %+.0f%% (p=%.2g); panel coherent (%d/%d member genes up with Myc, %d at ",
  "p<0.05). But Myc does NOT scale mtDNA-encoded output with it (%+.0f%%, p=%.2g) -- the ",
  "mitonuclear imbalance of script 24, now in ABSOLUTE SHARE rather than a mitoPPS ratio. ",
  "COMMISSIONED-BUT-UNBUILT: Myc DOES raise the mtDNA machinery -- nucleoid %+.0f%% ",
  "(p=%.2g), mt-transcription %+.0f%% (p=%.2g), mtDNA-replication %+.0f%% (p=%.2g) -- so ",
  "the imbalance is not neglect of the mtDNA arm; the machinery to make and read mtDNA is ",
  "built and the mtDNA-encoded OUTPUT does not follow. The genotype axis is %s. TIME is ",
  "NOT interpretable here: timepoint is perfectly confounded with sequencing cohort ",
  "(6W %.0f-%.0fM reads, 12W %.0f-%.0fM; depth ~ timepoint p=%.2g vs ~ genotype p=%.2g) ",
  "-- every temporal number, including script 29's mtDNA absolute z rise, inherits this. ",
  "CEILING: shares, not per-cell content; Myc's global RNA amplification makes the share ",
  "a LOWER BOUND; blot/qPCR/EM settles it."),
  pct_of("MASS_MARKERS_NOCHAP"),               stat_of("MASS_MARKERS_NOCHAP", "geno_p"),
  pct_of("MITOCARTA_NUCLEAR_ENCODED"),         stat_of("MITOCARTA_NUCLEAR_ENCODED", "geno_p"),
  pct_of("MITOCARTA_MITOCHONDRIAL_RIBOSOME"),  stat_of("MITOCARTA_MITOCHONDRIAL_RIBOSOME", "geno_p"),
  mass_coherence$n_myc_up, mass_coherence$n_genes, mass_coherence$n_myc_sig,
  pct_of("MITOCARTA_MTDNA_ENCODED"),           stat_of("MITOCARTA_MTDNA_ENCODED", "geno_p"),
  pct_of("MITOCARTA_MTDNA_NUCLEOID"),          stat_of("MITOCARTA_MTDNA_NUCLEOID", "geno_p"),
  pct_of("MITOCARTA_TRANSCRIPTION"),           stat_of("MITOCARTA_TRANSCRIPTION", "geno_p"),
  pct_of("MITOCARTA_MTDNA_REPLICATION"),       stat_of("MITOCARTA_MTDNA_REPLICATION", "geno_p"),
  qc_flags$geno_flag[qc_flags$panel == "MASS_MARKERS_NOCHAP"],
  min(qc_sample$depth_M[qc_sample$timepoint == "6W"]),
  max(qc_sample$depth_M[qc_sample$timepoint == "6W"]),
  min(qc_sample$depth_M[qc_sample$timepoint == "12W"]),
  max(qc_sample$depth_M[qc_sample$timepoint == "12W"]),
  depth_time_p, depth_geno_p)

message("\n", strwrap(content_verdict, width = 88) |> paste(collapse = "\n"), "\n")

# =============================================================================
# PART 7: FIGURES
# =============================================================================

# --- A: compartment shares across the 4 groups, faceted by lens tag -------------
plot_df_a <- shares |>
  dplyr::left_join(panel_roster[, c("panel", "tag")], by = "panel") |>
  dplyr::filter(!panel %in% c("MASS_MARKERS", "MITOCARTA_ALL")) |>
  dplyr::mutate(group = factor(group, levels = group_levels))

p_a <- ggplot2::ggplot(plot_df_a,
                       ggplot2::aes(x = group, y = share_nomt, colour = myc_status)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  ggplot2::geom_jitter(width = 0.15, size = 1.2, alpha = 0.8) +
  ggplot2::facet_wrap(~ panel, scales = "free_y", ncol = 4) +
  ggplot2::scale_colour_manual(values = geno_cols) +
  ggplot2::labs(
    title = "Mitochondrial compartment shares across the four groups",
    subtitle = paste("% of the non-mtDNA transcriptome (raw counts).",
                     "Genotype axis is depth-balanced; the TIME axis is cohort-confounded (see C)."),
    x = NULL, y = "share of non-mt transcriptome (%)", colour = "Myc") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggplot2::ggsave(file.path(out_dir, "A_compartment_shares.pdf"), p_a,
                width = 12, height = 8)

# --- B: the headline -- Myc raises the nuclear arms, not mtDNA ------------------
plot_df_b <- share_stats |>
  dplyr::filter(denominator == "share_nomt",
                !panel %in% c("MASS_MARKERS", "MITOCARTA_ALL")) |>
  dplyr::mutate(panel = stats::reorder(panel, geno_beta),
                sig = ifelse(geno_p < 0.05, "p < 0.05", "ns"))

p_b <- ggplot2::ggplot(plot_df_b,
                       ggplot2::aes(x = geno_beta, y = panel, fill = sig)) +
  ggplot2::geom_col() +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  ggplot2::facet_grid(tag ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_fill_manual(values = c("p < 0.05" = "#D73027", "ns" = "grey70")) +
  ggplot2::labs(
    title = "Myc raises the nuclear mitochondrial arms but not mtDNA-encoded output",
    subtitle = paste("Genotype main effect on log2(share of non-mt transcriptome).",
                     "The mitonuclear imbalance in ABSOLUTE share, not a mitoPPS ratio."),
    x = "genotype beta (log2 fold-change in share)", y = NULL, fill = NULL) +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "B_myc_vs_age_axes.pdf"), p_b,
                width = 9, height = 7)

# --- C: the QC gate, made visible rather than buried ----------------------------
p_c <- ggplot2::ggplot(qc_sample,
                       ggplot2::aes(x = depth_M, y = mt_pct, colour = myc_status,
                                    shape = timepoint)) +
  ggplot2::geom_point(size = 3, alpha = 0.85) +
  ggplot2::geom_text(ggplot2::aes(label = ifelse(sample == loo_sample, sample, "")),
                     hjust = -0.15, size = 2.5, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = geno_cols) +
  ggplot2::labs(
    title = "QC gate: the mtDNA share tracks library depth, and depth tracks TIMEPOINT",
    subtitle = sprintf(paste("6W and 12W are different sequencing cohorts (depth ~ timepoint",
                             "p=%.2g; ~ genotype p=%.2g).\nSo genotype contrasts are clean and",
                             "every temporal claim is exposed. Labelled = leave-one-out sample."),
                       depth_time_p, depth_geno_p),
    x = "library depth (M reads)", y = "mtDNA-encoded share of library (%)",
    colour = "Myc", shape = NULL) +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "C_qc_gate.pdf"), p_c, width = 8, height = 5.5)

message("Figures written to ", out_dir)

# =============================================================================
# PART 8: SAVE
# =============================================================================
content_out <- list(
  panels            = panel_roster,
  mass_roster       = mass_roster,
  shares            = shares,
  share_group_means = share_group_means,
  share_stats       = share_stats,
  mass_per_gene     = mass_per_gene,
  mass_coherence    = mass_coherence,
  qc = list(
    per_sample    = qc_sample,
    depth_balance = depth_balance,
    depth_geno_p  = depth_geno_p,
    depth_time_p  = depth_time_p,
    coupling      = qc_coupling,
    flags         = qc_flags,
    loo_sample    = loo_sample),
  loo_sensitivity   = loo_sensitivity,
  total_mito_rebuilt = total_mito_rebuilt,
  mt_dominance      = mt_dominance,
  reconciliation = list(
    vs_29_absolute   = recon_29,
    vs_24_mitopps    = recon_24,
    vs_24_agreement  = recon_24_agreement,
    chaperone        = recon_chaperone,
    raw_22_unread    = recon_22_raw),
  verdict = content_verdict,
  notes = paste(
    "Block B, author question 2026-07-16: is the 'robust increase in mitochondrial",
    "biogenesis' a CONTENT change or only a PROGRAM? Nothing in Block A answered it:",
    "reframe_supp3$abund ('AP-abund') is mitoPPS diff DISPERSION, not abundance; mitoPPS",
    "is abundance-blind by construction (08:404) so mtnuc_index is a priority ratio;",
    "total_mito_score (08:1153, captioned 'reflects mito content') was never saved or",
    "tested AND double-counts genes across 142 overlapping pathways (~68% mt-driven) --",
    "rebuilt here as total_mito_rebuilt. bio_comp (19:91) is a GSVA composite = program,",
    "not organelle. RESULT: Myc raises the nuclear mito share (mass markers, nuclear",
    "MitoCarta, mitoribosome) and does NOT scale mtDNA-encoded output with it = the",
    "script-24 mitonuclear imbalance restated in ABSOLUTE SHARE. COMMISSIONED-BUT-UNBUILT:",
    "Myc raises the mtDNA MACHINERY hardest of all (nucleoid d=2.32; mt-transcription",
    "Tfam/Polrmt/Tfb2m +43%) while mtDNA-encoded OUTPUT is flat -- so the imbalance is not",
    "neglect of the mtDNA arm, it is machinery built and output not following. Bulk cannot",
    "say why (copy number vs transcription rate vs turnover); mtDNA qPCR is the one",
    "experiment that separates them. MASS_MARKERS is a",
    "flagged exception to the 'do not rebuild gene sets' rule -- the library has no",
    "mass-marker set; it is the in-silico stand-in for the TOMM20/VDAC/CS/HSP60 blot we",
    "do not have. Hspa9/Hspd1 are standard mass markers AND direct MYC targets, so the",
    "CLAIM-BEARING composite is MASS_MARKERS_NOCHAP and the full panel is corroboration.",
    "THREE CEILINGS: (i) bulk polyA, no spike-ins, no cell counts -> per-cell content is",
    "NOT measurable; median-of-ratios normalisation removes that scale; every number is a",
    "SHARE OF TRANSCRIPTOME. (ii) Myc globally amplifies total RNA per cell, so a constant",
    "share already means more mitochondria per cell -> the share effect is a LOWER BOUND.",
    "(iii) transcript share is not protein and not organelle volume -- blot / mtDNA qPCR /",
    "EM SETTLES it; n=6/group. QC GATE: timepoint is PERFECTLY confounded with sequencing",
    "cohort (6W 18.8-29.0M reads from MYCF62-65/MYBS10x; 12W 9.9-16.1M from MYCF52-56) and",
    "the mt share tracks depth (rho ~ -0.47), so ALL temporal claims are TIME-EXPOSED --",
    "including script 29's mtDNA absolute z -1.14 -> +1.23. Genotype is depth-balanced",
    "within timepoint and litter-controlled (both genotypes inside MYCF62/MYCF52/MYCF56)",
    "-> CLEAN. The gate FLAGS rather than corrects: with timepoint == cohort there is no",
    "contrast that separates them. See docs/2026-07-13_BlockA_revision_walkthrough_and_",
    "intro_alignment.md."))
saveRDS(content_out, here::here("results", "mito_content_proxies.rds"))
message("Saved results/mito_content_proxies.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  mc <- readRDS(here::here("results", "mito_content_proxies.rds"))

  # --- The headline: does Myc raise the mitochondrial share? ---
  cat(strwrap(mc$verdict, width = 88), sep = "\n")

  mc$share_stats |>
    dplyr::filter(denominator == "share_nomt") |>
    dplyr::mutate(pct_change = 100 * (2^geno_beta - 1)) |>
    dplyr::select(panel, tag, n_genes, pct_change, geno_d, geno_p, int_p) |>
    dplyr::arrange(dplyr::desc(pct_change)) |>
    print(n = 20)

  # --- The mass panel gene by gene: is it coherent, or carried by the chaperones? ---
  mc$mass_per_gene |> print(n = 20)
  mc$mass_coherence

  # --- PART 2: group means, both denominators ---
  mc$share_group_means |>
    dplyr::filter(denominator == "share_nomt") |> print(n = 20)

  # --- PART 4: the QC gate. Which panels survive, and what is time-exposed? ---
  mc$qc$per_sample |> print(n = 24)
  mc$qc$depth_balance |> print()
  c(depth_vs_genotype_p = mc$qc$depth_geno_p, depth_vs_timepoint_p = mc$qc$depth_time_p)
  mc$qc$flags |> print(n = 20)
  mc$loo_sensitivity |> print(n = 20)

  # --- The "commissioned but unbuilt" dissociation: machinery up, output flat ---
  mc$share_stats |>
    dplyr::filter(denominator == "share_nomt",
                  panel %in% c("MITOCARTA_MTDNA_NUCLEOID", "MITOCARTA_TRANSCRIPTION",
                               "MITOCARTA_MTDNA_REPLICATION", "MITOCARTA_MTDNA_ENCODED")) |>
    dplyr::mutate(pct_change = 100 * (2^geno_beta - 1)) |>
    dplyr::select(panel, n_genes, pct_change, geno_d, geno_p) |> print()

  # --- PART 5a: why script 08's total_mito_score was null ---
  mc$mt_dominance                 # what fraction of MitoCarta counts are just the 13 mt-*
  mc$total_mito_rebuilt |> print()

  # --- PART 5b/c: do the abundance-blind and abundance-capable lenses agree? ---
  mc$reconciliation$vs_29_absolute |> print()
  mc$reconciliation$vs_24_mitopps |> print()
  mc$reconciliation$vs_24_agreement |> print()

  # --- PART 5d: the chaperone reconciliation (script 24's UPR^mt result is safe) ---
  mc$reconciliation$chaperone |> print()

  # --- PART 5e: the raw-abundance rows script 22 computed and never read ---
  mc$reconciliation$raw_22_unread |> print(n = 8)

  list.files(here::here("outputs", "mito_content_proxies"), pattern = "\\.pdf$")
}
