# scripts/33_mtdna_axis_and_coupling_null.R
# =============================================================================
# What IS the mt-transcript axis, and do the per-sample death couplings beat a
# null? (Block B, author challenge 2026-07-16)
# =============================================================================
#
# WHY THIS EXISTS. The author asked two questions that this corpus never tested:
#   (1) "These are purified MECs (enzymatic dissociation) -- any other cell type is
#       contamination. But why
#       would only mtDNA-coded genes contaminate? Shouldn't the nuclear be there too?"
#   (2) Is the mt-transcript signal a stress artifact -- and if so, does the
#       mitonuclear imbalance survive?
# Both questions are sharper than the analyses they interrogate. Answering (1)
# refutes contamination and exposes what the mt axis actually is; answering (2)
# turns out to require a NULL that the death spine never had.
#
# THE THREE FINDINGS THIS SCRIPT ESTABLISHES (read-only reframe; nothing re-fitted):
#
#   A. mt% IS NOT CONTAMINATION, and it tracks an IEG-marked dominant axis. The
#      contaminating fraction is only ~3-12% of cells, so to lift mt% from 3.4% to
#      40.1% by MIXING the contaminant would need an mt-share of 126-737% --
#      impossible. The author's asymmetry is the tell: nuclear MEC genes track
#      contamination NEGATIVELY (dilution, as mixing predicts) while mt% tracks it
#      POSITIVELY. Mixing cannot produce opposite signs. What mt% DOES track is a
#      dominant sample-level axis marked by the canonical dissociation IEG signature
#      (rho +0.87 within 6W); adding it lifts the mt% R2 from 0.26 to 0.72.
#
#      *** WITHDRAWN (2026-07-17, author's catch) -- the "direction test". ***
#      An earlier version argued: "if mt% measured mitochondrial content, nuclear
#      mito genes would RISE with it; they FALL (rho -0.65); therefore mt% is not a
#      mito readout." THAT IS CIRCULAR. The mitonuclear imbalance IS the claim that
#      mt-encoded and nuclear-encoded dissociate -- so an anticorrelation is what the
#      imbalance PREDICTS, not evidence against a mitochondrial reading. The argument
#      assumed the phenomenon's absence to prove the measurement broken. So THIS
#      SCRIPT DOES NOT ESTABLISH THAT mt% IS NON-MITOCHONDRIAL, and PART A does NOT
#      refute the imbalance. mt% may be the mtDNA arm of a real imbalance whose
#      between-sample variation happens to track the IEG axis; bulk cannot separate
#      those. The load-bearing result is PART C, which is INDEPENDENT of this
#      question by construction -- that is now doing more work than intended.
#
#   B. THE MITONUCLEAR IMBALANCE IS MOSTLY THAT AXIS. imbalance ~ mt%: r = -0.85,
#      R^2 = 0.73. Script 24's mtnuc_index is ~three-quarters the mt% metric.
#
#   C. THE DEATH COUPLING IS NOT SPECIFIC -- THE CORE RESULT. The published
#      "imbalance couples to pro-apoptotic priming at 6W (r=0.75, p=0.005)" sits at
#      the ~61st percentile of the imbalance's couplings to ALL 884 library gene
#      sets (median +0.62); permutation-style p ~ 0.49. ANY gene set gives r~0.6.
#      Without a null, r=0.75 / p=0.005 looks compelling and means nothing. This is
#      exactly the AP6 logic script 21 applied to fGSEA -- never applied to the
#      per-sample couplings. PART C is that null.
#
# WHAT THIS DOES *NOT* CLAIM. The IEG/mt axis is genotype-INDEPENDENT (p=0.49), so
# it CANNOT bias any Myc contrast -- every genotype result in the corpus, including
# script 32's content claim, is untouched. It IS time-associated (p=0.0003), so it
# entangles the time axis. Whether the axis is TECHNICAL or BIOLOGICAL is NOT
# resolved: the samples are ENZYMATICALLY DISSOCIATED (author 2026-07-17 -- NOT
# FACS-sorted, an earlier version of this header said FACS; warm collagenase digest
# is precisely what van den Brink 2017 characterised, so the IEG reading is if
# anything strengthened, while "purified" is weaker -- consistent with the 3-12%
# residual stroma). No dissociation-batch, viability or RIN metadata is on disk, so
# the IEG signature is an INFERRED covariate, not a measured one. PART D states both
# readings. The specificity failure in PART C holds under EITHER -- it is a property
# of the coupling, not of the axis's origin. That is why PART C, not PART A, is the
# load-bearing part.
#
#   D2 (NEW 2026-07-17, from the author's question "why is mtDNA so closely
#   correlated with all the others, when those signals look so different?"). Because
#   THEY ARE NOT DIFFERENT SIGNALS. Within 6W, PC1 of the non-mt transcriptome
#   explains ~43% of the variance among 12 samples and EVERYTHING loads on it: IEG
#   -0.71, proliferation +0.62, contamination -0.62, MYC activity +0.59, mt% -0.51.
#   That is the STRUCTURAL reason the PART C null is flat -- with one axis at 43% and
#   n=12, pairwise couplings between composites carry almost no information about
#   specific mechanisms. It is not a quirk of the death sets. Note the axis is MIXED:
#   contamination (which CANNOT be biological -- it is residual stroma in the prep)
#   loads on it alongside proliferation (which can). So it is neither cleanly
#   technical nor cleanly biological, and this script does not pretend otherwise.
#
#   D3 (NEW). Are IEGs simply part of the pubertal TEB/proliferative programme (a
#   good alternative the author raised)? NO -- refuted on direction AND timing. A
#   growth-factor/TEB reading predicts IEGs POSITIVELY coupled to proliferation and
#   HIGHER at 6W (the TEB-rich timepoint). Observed: IEG ~ proliferation rho -0.62
#   (all) / -0.76 (6W); IEG ~ TEB sets -0.57; IEG ~ MYC activity -0.41; and IEG is
#   LOWER at 6W (0.72/0.60) than 12W (1.03/1.02). Both predictions fail. IEGs here
#   anticorrelate with the proliferative biology, as an overlay would.
#
# Input:  results/mito_content_proxies.rds        (script 32 shares + QC)
#         results/mitopps_scores.rds              (mitoPPS; mtdna_genes_separated)
#         results/developmental_substrate_death.rds (script 25 per-sample composites)
#         results/gsva_scores.rds                 (884 scored sets = the null universe)
#         results/count_matrix.rds, results/dds_int_run.rds
#         results/combined_df_annotated.rds       (symbol <-> Ensembl)
# Output: results/mtdna_axis_and_coupling_null.rds
#         outputs/mtdna_axis_null/*.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

out_dir <- here::here("outputs", "mtdna_axis_null")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

group_levels <- c("6W_neg", "6W_pos", "12W_neg", "12W_pos")
geno_cols    <- c(neg = "#4575B4", pos = "#D73027")

# =============================================================================
# PART 1: LOAD + IDENTIFIER MAP
# =============================================================================
mc  <- readRDS(here::here("results", "mito_content_proxies.rds"))
mp  <- readRDS(here::here("results", "mitopps_scores.rds"))
ds  <- readRDS(here::here("results", "developmental_substrate_death.rds"))$per_sample
gs  <- readRDS(here::here("results", "gsva_scores.rds"))
cts <- readRDS(here::here("results", "count_matrix.rds"))
cdf <- readRDS(here::here("results", "combined_df_annotated.rds"))

qc      <- mc$qc$per_sample
samples <- qc$sample                                        # master order
cts     <- cts[, samples, drop = FALSE]
scores  <- gs$scores[, samples, drop = FALSE]
ds      <- ds[match(samples, ds$sample), ]
stopifnot(identical(ds$sample, samples), identical(colnames(scores), samples))

sym2ens <- cdf[!is.na(cdf$mgi_symbol) & !duplicated(cdf$mgi_symbol),
               c("mgi_symbol", "gene")]
ens_of  <- function(syms, universe = rownames(cts)) {
  e <- sym2ens$gene[match(intersect(syms, sym2ens$mgi_symbol), sym2ens$mgi_symbol)]
  intersect(e, universe)
}

mt_ens   <- ens_of(mp$mtdna_genes_separated)
den_nomt <- colSums(cts[setdiff(rownames(cts), mt_ens), , drop = FALSE])
shr      <- function(syms) 100 * colSums(cts[ens_of(syms), , drop = FALSE]) / den_nomt

i6  <- qc$timepoint == "6W"                                 # the death spine's window
i12 <- qc$timepoint == "12W"

# =============================================================================
# PART 2: MARKER PANELS -- the axes the author's question forces us to separate
# =============================================================================
# CONTAMINATION: in purified MECs these are ~absent, so any signal is foreign cells.
# IEG/STRESS: the canonical dissociation-artifact panel (van den Brink 2017) -- the
#   transcriptional response to warm enzymatic digest -- exactly this protocol. NOT death.
panels <- list(
  endothelial = c("Pecam1", "Cdh5", "Ptprb", "Adgrl4", "Mmrn2", "Sox17", "Tie1",
                  "Kdr", "Egfl7", "Cldn5", "Emcn"),
  adipocyte   = c("Adipoq", "Lep", "Plin1", "Cidec", "Retn", "Lipe", "Pparg"),
  fibroblast  = c("Col1a1", "Col1a2", "Pdgfra", "Dcn", "Lum", "Fbn1"),
  immune      = c("Ptprc", "Cd68", "Adgre1", "Cd3e", "Cd19"),
  epithelial  = c("Krt8", "Krt18", "Krt5", "Krt14", "Epcam", "Cdh1"),
  ieg_stress  = c("Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Ier2", "Ier3",
                  "Atf3", "Hspa1a", "Hspa1b", "Socs3", "Zfp36"))
panel_share <- vapply(panels, shr, numeric(length(samples)))
contam      <- log2(panel_share[, "endothelial"])           # contamination proxy
stress      <- log2(panel_share[, "ieg_stress"])            # inferred prep-quality axis
mt_pct      <- qc$mt_pct

cor_p <- function(a, b, m = "spearman") {
  ct <- suppressWarnings(stats::cor.test(a, b, method = m))
  c(rho = unname(ct$estimate), p = ct$p.value)
}

# =============================================================================
# PART A: WHAT IS THE mt AXIS? (three independent refutations of contamination)
# =============================================================================

# --- A1. The author's asymmetry: mixing cannot give opposite signs -------------
# Real contamination DILUTES MEC transcripts -> nuclear MEC genes fall. If mt% rose
# from the same mixing it would have to fall too. It rises. Different mechanism.
asymmetry <- purrr::map_dfr(
  c("MITOCARTA_NUCLEAR_ENCODED", "MITOCARTA_OXPHOS_NU", "MASS_MARKERS_NOCHAP",
    "MITOCARTA_MITOCHONDRIAL_RIBOSOME"),
  function(p) {
    s <- mc$shares |> dplyr::filter(panel == p)
    s <- s[match(samples, s$sample), ]
    a <- cor_p(s$share_nomt, contam); b <- cor_p(s$share_nomt[i6], mt_pct[i6])
    tibble::tibble(measure = p,
                   rho_vs_contamination = a[["rho"]], p_vs_contamination = a[["p"]],
                   rho_vs_mtpct_6W = b[["rho"]],      p_vs_mtpct_6W = b[["p"]])
  }) |>
  dplyr::bind_rows(tibble::tibble(
    measure = "mt_pct (for contrast)",
    rho_vs_contamination = cor_p(mt_pct, contam)[["rho"]],
    p_vs_contamination   = cor_p(mt_pct, contam)[["p"]],
    rho_vs_mtpct_6W = NA_real_, p_vs_mtpct_6W = NA_real_))

# --- A2. The arithmetic: contamination is far too small to do this -------------
# In PURE endothelium this marker panel runs ~1-3% of the transcriptome, so
# implied contaminant fraction f ~ observed_share / 2%. For bulk mt% to go from
# min to max by MIXING, the contaminant needs mt-share m_c solving
#   mt_max = (1-f)*mt_min + f*m_c  =>  m_c = (mt_max - (1-f)*mt_min)/f.
# If m_c > 100% the mixing explanation is arithmetically impossible.
contam_frac <- panel_share[, "endothelial"] / 2            # crude but bounded
arithmetic_refutation <- purrr::map_dfr(c(0.05, 0.10, 0.20, 0.30), function(f) {
  m_c <- (max(mt_pct) / 100 - (1 - f) * min(mt_pct) / 100) / f
  tibble::tibble(assumed_contaminant_fraction = f,
                 required_contaminant_mt_share_pct = 100 * m_c,
                 possible = 100 * m_c <= 100)
})

# --- A3. The nuclear-vs-mt anticorrelation -- REPORTED, NOT INTERPRETED ---------
# WITHDRAWN AS EVIDENCE (author's catch, 2026-07-17). This was published here as a
# "direction test": nuclear mito genes FALL with mt%, so mt% is not a mito readout.
# That is CIRCULAR -- the mitonuclear imbalance IS the claim that these two arms
# dissociate, so the anticorrelation is equally what a REAL imbalance predicts. The
# number is kept because it is a fact about the data; the inference is not.
mt_vs_nuclear_anticorrelation <- asymmetry |>
  dplyr::filter(!is.na(rho_vs_mtpct_6W)) |>
  dplyr::summarise(n_measures = dplyr::n(),
                   median_rho = stats::median(rho_vs_mtpct_6W),
                   reads_as = paste(
                     "nuclear mito arms anticorrelate with mt% (median rho",
                     sprintf("%+.2f)", stats::median(rho_vs_mtpct_6W)),
                     "-- CONSISTENT WITH EITHER a broken mt metric OR a genuine",
                     "mitonuclear imbalance. NOT evidence for either. Do not cite",
                     "as a direction test (withdrawn 2026-07-17)."))

# --- A4. What the mt axis DOES track ------------------------------------------
axis_tracking <- purrr::map_dfr(names(panels), function(n) {
  a <- cor_p(panel_share[, n], mt_pct); b <- cor_p(panel_share[i6, n], mt_pct[i6])
  tibble::tibble(panel = n, rho_all = a[["rho"]], p_all = a[["p"]],
                 rho_6W = b[["rho"]], p_6W = b[["p"]])
}) |>
  dplyr::arrange(dplyr::desc(rho_6W))

mt_variance_model <- tibble::tibble(
  model = c("group only", "group + contamination + IEG"),
  r2 = c(summary(stats::lm(log2(mt_pct) ~ qc$group))$r.squared,
         summary(stats::lm(log2(mt_pct) ~ qc$group + contam + stress))$r.squared))

# --- A5. Is the axis genotype- or time-associated? (what it can and cannot bias)
axis_design <- purrr::map_dfr(
  list(contamination = contam, ieg_stress = stress, mt_pct = log2(mt_pct)),
  function(v) {
    m <- summary(stats::lm(v ~ qc$myc_status + qc$timepoint))$coefficients
    tibble::tibble(p_genotype = m["qc$myc_statuspos", "Pr(>|t|)"],
                   p_timepoint = m["qc$timepoint12W", "Pr(>|t|)"])
  }, .id = "axis") |>
  dplyr::mutate(
    can_bias_genotype = p_genotype < 0.05,
    can_bias_time     = p_timepoint < 0.05)

# =============================================================================
# PART B: THE MITONUCLEAR IMBALANCE IS MOSTLY THE mt AXIS
# =============================================================================
# Rebuild script 24's mtnuc_index PER SAMPLE (24 saves only the 4 group means), so
# it can be regressed. Group means are checked against 24 in the sandbox.
pps    <- mp$mitopps_scores
rownames(pps) <- pps$sample
pps    <- pps[samples, ]
num    <- names(pps)[vapply(pps, is.numeric, logical(1))]
mtnm   <- mp$mtdna_pathway_name
ox_nuc <- setdiff(intersect(names(mp$pathway_tier1_map)[mp$pathway_tier1_map == "OXPHOS"],
                            num), mtnm)
imbalance <- rowMeans(pps[, ox_nuc, drop = FALSE]) - pps[[mtnm]]

imbalance_decomposition <- tibble::tibble(
  term = c("imbalance ~ mt_pct", "imbalance ~ IEG/stress", "imbalance ~ contamination"),
  rho  = c(cor_p(imbalance, log2(mt_pct))[["rho"]], cor_p(imbalance, stress)[["rho"]],
           cor_p(imbalance, contam)[["rho"]]),
  p    = c(cor_p(imbalance, log2(mt_pct))[["p"]], cor_p(imbalance, stress)[["p"]],
           cor_p(imbalance, contam)[["p"]]),
  r2_linear = c(summary(stats::lm(imbalance ~ log2(mt_pct)))$r.squared,
                summary(stats::lm(imbalance ~ stress))$r.squared,
                summary(stats::lm(imbalance ~ contam))$r.squared))

imbalance_group_means <- tapply(imbalance, qc$group, mean)[group_levels]

# =============================================================================
# PART C: THE COUPLING NULL -- the load-bearing part
# =============================================================================
# THE QUESTION script 25 never asked. It reports per-sample couplings to the
# pro-death composite at 6W (imbalance r=0.75 p=0.005; bio_comp 0.83; MASC 0.78) as
# evidence for a death-permissive substrate. But with n=12 and one dominant
# sample-level axis, ANY two composites correlate. The test is not "is r big" but
# "is r BIGGER THAN what an arbitrary gene set gives" -- the AP6 logic of script 21,
# applied to sample-level couplings instead of fGSEA.
#
# Null universe = all 884 GSVA-scored library sets. For a candidate axis X, compute
# rho(X, every set) and locate rho(X, pro_comp) in that distribution.
null_for_axis <- function(x, label, idx = i6) {
  r_sets <- apply(scores[, idx, drop = FALSE], 1,
                  function(s) suppressWarnings(stats::cor(s, x[idx], method = "spearman")))
  r_sets <- r_sets[is.finite(r_sets)]
  r_obs  <- suppressWarnings(stats::cor(ds$pro_comp[idx], x[idx], method = "spearman"))
  tibble::tibble(
    axis = label, rho_vs_pro_comp = r_obs,
    null_median = stats::median(r_sets), null_q90 = stats::quantile(r_sets, 0.90),
    n_sets = length(r_sets),
    frac_sets_above_0.5 = mean(abs(r_sets) > 0.5),
    empirical_percentile = 100 * mean(r_sets < r_obs),
    perm_p = mean(abs(r_sets) >= abs(r_obs)),      # two-sided empirical p
    verdict = dplyr::case_when(
      mean(abs(r_sets) >= abs(r_obs)) > 0.10 ~ "NOT SPECIFIC (any gene set does this)",
      mean(abs(r_sets) >= abs(r_obs)) > 0.05 ~ "marginal",
      TRUE ~ "beats the null"))
}

coupling_null <- dplyr::bind_rows(
  null_for_axis(imbalance,                     "mitonuclear_imbalance (script 24/25)"),
  null_for_axis(ds$bio_comp,                   "bio_comp (script 25 Part B)"),
  null_for_axis(ds$MASC_comp,                  "MASC_comp (script 25 Part B)"),
  null_for_axis(log2(mt_pct),                  "mt_pct (the raw axis)"),
  null_for_axis(stress,                        "IEG/stress (the inferred prep axis)"))

# The full null distributions, kept for the figure
null_distributions <- purrr::map_dfr(
  list(mitonuclear_imbalance = imbalance, bio_comp = ds$bio_comp,
       MASC_comp = ds$MASC_comp),
  function(x) {
    r <- apply(scores[, i6, drop = FALSE], 1,
               function(s) suppressWarnings(stats::cor(s, x[i6], method = "spearman")))
    tibble::tibble(rho = r[is.finite(r)])
  }, .id = "axis")

# =============================================================================
# PART D: WHAT SURVIVES -- and the technical-vs-biological adjudication
# =============================================================================

# --- D1. Every genotype claim is untouched (the axis is genotype-independent) ---
genotype_untouched <- purrr::map_dfr(
  c("MASS_MARKERS_NOCHAP", "MITOCARTA_NUCLEAR_ENCODED", "MITOCARTA_OXPHOS_NU",
    "MITOCARTA_MITOCHONDRIAL_RIBOSOME", "MITOCARTA_MTDNA_ENCODED"),
  function(p) {
    s <- mc$shares |> dplyr::filter(panel == p); s <- s[match(samples, s$sample), ]
    y <- log2(s$share_nomt)
    r <- summary(stats::lm(y ~ qc$myc_status + qc$timepoint))$coefficients["qc$myc_statuspos", ]
    a <- summary(stats::lm(y ~ qc$myc_status + qc$timepoint + stress + contam))$coefficients["qc$myc_statuspos", ]
    tibble::tibble(panel = p,
                   raw_pct = 100 * (2^r[["Estimate"]] - 1), raw_p = r[["Pr(>|t|)"]],
                   adj_pct = 100 * (2^a[["Estimate"]] - 1), adj_p = a[["Pr(>|t|)"]],
                   pct_retained = 100 * a[["Estimate"]] / r[["Estimate"]])
  })

# --- D2. Over-adjustment check: adjusting for a covariate ON the causal path
# would destroy real biology. It cannot here -- the axis is genotype-independent,
# so it carries no Myc signal to remove. Quantified rather than asserted.
overadjust_guard <- tibble::tibble(
  check = c("IEG axis ~ genotype", "contamination ~ genotype"),
  p = c(axis_design$p_genotype[axis_design$axis == "ieg_stress"],
        axis_design$p_genotype[axis_design$axis == "contamination"]),
  reads_as = "genotype-independent => adjusting cannot remove Myc biology => genotype claims safe")

# --- D2. WHY does everything correlate with everything? The structural answer ----
# Author's question (2026-07-17): "I'm surprised mtDNA expression is so closely
# correlated with all the others, since the other signals seem so different." They
# are not different signals. One axis dominates and every measure is a projection of
# it -- which is WHY the PART C null is flat. Not a quirk of the death sets.
vst_mat  <- SummarizedExperiment::assay(
  DESeq2::vst(readRDS(here::here("results", "dds_int_run.rds")), blind = TRUE))
vst_mat  <- vst_mat[setdiff(rownames(vst_mat), mt_ens), samples, drop = FALSE]
v6       <- vst_mat[order(-matrixStats::rowVars(vst_mat[, i6, drop = FALSE]))[1:2000],
                    i6, drop = FALSE]
pca6     <- stats::prcomp(t(v6), scale. = TRUE)
pc_var   <- 100 * summary(pca6)$importance[2, 1:3]
prolif_c <- colMeans(scores[grep("^PROLIF_", rownames(scores)), , drop = FALSE])
myc_c    <- colMeans(scores[grep("^MYC_",    rownames(scores)), , drop = FALSE])
teb_c    <- colMeans(scores[grep("TEB",      rownames(scores)), , drop = FALSE])

dominant_axis <- tibble::tibble(
  loads_on_PC1 = c("IEG/stress", "proliferation", "contamination", "MYC activity", "mt_pct"),
  rho = c(cor_p(pca6$x[, 1], stress[i6])[["rho"]],
          cor_p(pca6$x[, 1], prolif_c[i6])[["rho"]],
          cor_p(pca6$x[, 1], contam[i6])[["rho"]],
          cor_p(pca6$x[, 1], myc_c[i6])[["rho"]],
          cor_p(pca6$x[, 1], log2(mt_pct[i6]))[["rho"]]),
  p   = c(cor_p(pca6$x[, 1], stress[i6])[["p"]],
          cor_p(pca6$x[, 1], prolif_c[i6])[["p"]],
          cor_p(pca6$x[, 1], contam[i6])[["p"]],
          cor_p(pca6$x[, 1], myc_c[i6])[["p"]],
          cor_p(pca6$x[, 1], log2(mt_pct[i6]))[["p"]])) |>
  dplyr::mutate(pc1_variance_pct = unname(pc_var[1]),
                pc2_variance_pct = unname(pc_var[2]),
                reads_as = paste0(
                  "PC1 = ", sprintf("%.0f%%", pc_var[1]), " of within-6W variance and EVERY measure ",
                  "loads on it => these are not independent signals => pairwise composite ",
                  "couplings are structurally uninformative at n=12. NOTE the axis is MIXED: ",
                  "contamination (cannot be biological) loads alongside proliferation (can) -- ",
                  "so it is neither cleanly technical nor cleanly biological."))

# --- D3. Are IEGs just the pubertal TEB / proliferative programme? --------------
# The author's alternative (2026-07-17). A growth-factor/TEB reading predicts IEGs
# POSITIVELY coupled to proliferation and HIGHER at 6W. Both predictions fail.
ieg_teb_test <- tibble::tibble(
  test = c("IEG ~ proliferation (all)", "IEG ~ proliferation (6W)",
           "IEG ~ TEB sets (all)", "IEG ~ MYC activity (all)"),
  rho = c(cor_p(stress, prolif_c)[["rho"]], cor_p(stress[i6], prolif_c[i6])[["rho"]],
          cor_p(stress, teb_c)[["rho"]], cor_p(stress, myc_c)[["rho"]]),
  teb_story_predicts = "POSITIVE") |>
  dplyr::mutate(consistent_with_teb = rho > 0)

ieg_group_means <- tapply(2^stress, qc$group, mean)[group_levels]
ieg_teb_verdict <- sprintf(paste0(
  "IEGs are NOT the TEB/proliferative programme: they anticorrelate with ",
  "proliferation (rho %+.2f all / %+.2f within 6W) and with TEB sets (%+.2f), and ",
  "they are LOWER at 6W (%.2f/%.2f) than 12W (%.2f/%.2f) -- the TEB story predicts ",
  "the opposite on BOTH direction and timing. IEGs here anticorrelate with the ",
  "proliferative biology, as an overlay would."),
  ieg_teb_test$rho[1], ieg_teb_test$rho[2], ieg_teb_test$rho[3],
  ieg_group_means[["6W_neg"]], ieg_group_means[["6W_pos"]],
  ieg_group_means[["12W_neg"]], ieg_group_means[["12W_pos"]])

# --- D4. The adjudication we CANNOT make, stated as such -----------------------
technical_vs_biological <- tibble::tribble(
  ~evidence,                                            ~favours,      ~note,
  "IEG panel is the canonical dissociation signature",  "technical",   "van den Brink 2017; warm enzymatic digest induces it in LIVE cells -- exactly this protocol",
  "mt% high = the scRNA-seq stressed/permeabilised tell","technical",   "cytoplasmic mRNA lost, mt transcripts retained in organelle",
  "axis is genotype-independent (p~0.49)",              "either",      "rules out a Myc-driven biological axis; consistent with prep",
  "axis is time-associated (p~0.0003)",                 "either",      "6W/12W are different animal cohorts = different prep sessions",
  "no dissociation-batch / viability / RIN metadata",    "unresolvable","enzymatic dissociation, NOT FACS (author 2026-07-17); IEG panel is INFERRED, not measured",
  "PART C specificity failure holds under both",        "n/a",         "the coupling is unremarkable whatever the axis IS -- that is the point")

# =============================================================================
# PART E: VERDICT
# =============================================================================
cn <- function(a, col) coupling_null[[col]][coupling_null$axis == a]
imb_lab <- "mitonuclear_imbalance (script 24/25)"

axis_verdict <- sprintf(paste0(
  "THE mt AXIS IS NOT CONTAMINATION -- but whether it is MITOCHONDRIAL is NOT settled ",
  "here. (a) Arithmetic: contamination is ~3-12%% of cells, so mixing would need a ",
  "contaminant mt-share of %.0f-%.0f%% to span the observed 3.4-40.1%% -- IMPOSSIBLE. ",
  "(b) Asymmetry (the author's point): nuclear MEC genes track contamination ",
  "NEGATIVELY (dilution, as mixing predicts) while mt%% tracks it POSITIVELY -- mixing ",
  "cannot give opposite signs. What mt%% DOES track is a dominant sample-level axis ",
  "marked by the dissociation IEG signature (rho %+.2f within 6W); adding it lifts the ",
  "mt%% R2 from %.2f to %.2f, and that axis is NOT the TEB/proliferative programme ",
  "(IEG ~ proliferation rho %+.2f within 6W; IEGs LOWER at 6W -- both TEB predictions ",
  "fail). *** WITHDRAWN: the 'direction test'. *** An earlier verdict argued nuclear ",
  "mito genes FALL with mt%% (median rho %+.2f) therefore mt%% is not mitochondrial. ",
  "That is CIRCULAR -- the mitonuclear imbalance IS the claim that these arms ",
  "dissociate, so the anticorrelation is equally what a REAL imbalance predicts. PART A ",
  "therefore does NOT refute the imbalance; mt%% may be its mtDNA arm with ",
  "between-sample variation that happens to track the IEG axis. Bulk cannot separate ",
  "them. STRUCTURAL CONTEXT: within 6W, PC1 = %.0f%% of variance and EVERY measure ",
  "loads on it (IEG %+.2f, proliferation %+.2f, contamination %+.2f, MYC %+.2f) -- ",
  "these are not independent signals, which is WHY the PART C null is flat. The axis is ",
  "MIXED (contamination cannot be biological; proliferation can), so it is neither ",
  "cleanly technical nor cleanly biological."),
  min(arithmetic_refutation$required_contaminant_mt_share_pct),
  max(arithmetic_refutation$required_contaminant_mt_share_pct),
  axis_tracking$rho_6W[axis_tracking$panel == "ieg_stress"],
  mt_variance_model$r2[1], mt_variance_model$r2[2],
  ieg_teb_test$rho[2], mt_vs_nuclear_anticorrelation$median_rho,
  dominant_axis$pc1_variance_pct[1],
  dominant_axis$rho[dominant_axis$loads_on_PC1 == "IEG/stress"],
  dominant_axis$rho[dominant_axis$loads_on_PC1 == "proliferation"],
  dominant_axis$rho[dominant_axis$loads_on_PC1 == "contamination"],
  dominant_axis$rho[dominant_axis$loads_on_PC1 == "MYC activity"])

coupling_verdict <- sprintf(paste0(
  "THE DEATH COUPLING IS NOT SPECIFIC. The mitonuclear imbalance is %.0f%% the mt%% ",
  "metric (r=%+.2f, R2=%.2f). Its published coupling to pro-death priming at 6W ",
  "(rho=%+.2f) sits at the %.0fth PERCENTILE of its couplings to all %d library gene ",
  "sets (median %+.2f; %.0f%% of sets exceed |rho|=0.5) -- permutation p=%.2f => %s. ",
  "Any gene set gives ~%.2f. Without a null, r=0.75/p=0.005 looks compelling and ",
  "means nothing; this is the AP6 logic of script 21, never applied to the ",
  "per-sample couplings. Comparators: bio_comp p=%.2f (%s); MASC_comp p=%.2f (%s). ",
  "CRITICALLY: this holds whether the axis is technical or biological -- specificity ",
  "is a property of the coupling, not of the axis's origin. GENOTYPE CLAIMS ARE ",
  "UNTOUCHED: the axis is genotype-independent (p=%.2f), so it cannot bias any Myc ",
  "contrast; script 32's content claim retains %.0f%% of its effect after adjustment."),
  100 * imbalance_decomposition$r2_linear[1],
  imbalance_decomposition$rho[1], imbalance_decomposition$r2_linear[1],
  cn(imb_lab, "rho_vs_pro_comp"), cn(imb_lab, "empirical_percentile"),
  cn(imb_lab, "n_sets"), cn(imb_lab, "null_median"),
  100 * cn(imb_lab, "frac_sets_above_0.5"), cn(imb_lab, "perm_p"),
  cn(imb_lab, "verdict"), cn(imb_lab, "null_median"),
  cn("bio_comp (script 25 Part B)", "perm_p"), cn("bio_comp (script 25 Part B)", "verdict"),
  cn("MASC_comp (script 25 Part B)", "perm_p"), cn("MASC_comp (script 25 Part B)", "verdict"),
  overadjust_guard$p[1],
  genotype_untouched$pct_retained[genotype_untouched$panel == "MASS_MARKERS_NOCHAP"])

message("\n", paste(strwrap(axis_verdict, width = 88), collapse = "\n"))
message("\n", paste(strwrap(coupling_verdict, width = 88), collapse = "\n"), "\n")

# =============================================================================
# PART F: FIGURES
# =============================================================================

# --- A: the author's asymmetry, in one panel ----------------------------------
asym_df <- dplyr::bind_rows(
  tibble::tibble(sample = samples, axis = "nuclear MEC mito genes\n(share of non-mt)",
                 x = log2(panel_share[, "endothelial"]),
                 y = as.numeric(scale(log2((mc$shares |> dplyr::filter(panel == "MITOCARTA_NUCLEAR_ENCODED"))$share_nomt[match(samples, (mc$shares |> dplyr::filter(panel == "MITOCARTA_NUCLEAR_ENCODED"))$sample)])))),
  tibble::tibble(sample = samples, axis = "mtDNA-encoded share\n(mt%)",
                 x = log2(panel_share[, "endothelial"]),
                 y = as.numeric(scale(log2(mt_pct)))))
asym_df$group <- qc$group[match(asym_df$sample, samples)]

p_a <- ggplot2::ggplot(asym_df, ggplot2::aes(x = x, y = y, colour = group)) +
  ggplot2::geom_point(size = 2.5, alpha = 0.85) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "grey40",
                       linewidth = 0.6, formula = y ~ x) +
  ggplot2::facet_wrap(~ axis) +
  ggplot2::labs(
    title = "Contamination cannot explain the mt signal: nuclear and mt move OPPOSITE ways",
    subtitle = paste("Purified MECs. Real contamination DILUTES MEC transcripts (left, negative).",
                     "\nmt% RISES with it (right) -- mixing cannot produce opposite signs."),
    x = "log2 endothelial marker share (contamination proxy)",
    y = "z-scored measure", colour = NULL) +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "A_contamination_asymmetry.pdf"), p_a,
                width = 9, height = 4.5)

# --- B: THE NULL -- the figure that carries the argument -----------------------
obs_pts <- coupling_null |>
  dplyr::filter(axis %in% c("mitonuclear_imbalance (script 24/25)",
                            "bio_comp (script 25 Part B)", "MASC_comp (script 25 Part B)")) |>
  dplyr::mutate(axis = sub(" \\(.*", "", axis))

p_b <- ggplot2::ggplot(null_distributions, ggplot2::aes(x = rho)) +
  ggplot2::geom_histogram(bins = 50, fill = "grey75", colour = NA) +
  ggplot2::geom_vline(data = obs_pts, ggplot2::aes(xintercept = rho_vs_pro_comp),
                      colour = "#D73027", linewidth = 0.9) +
  ggplot2::geom_vline(data = obs_pts, ggplot2::aes(xintercept = null_median),
                      colour = "#4575B4", linewidth = 0.6, linetype = 2) +
  ggplot2::facet_wrap(~ axis, ncol = 1, scales = "free_y") +
  ggplot2::labs(
    title = "Do the per-sample death couplings beat a null? (the test script 25 never ran)",
    subtitle = paste("Grey = that axis's coupling to ALL 884 library gene sets at 6W.",
                     "RED = its published coupling to the pro-death composite.",
                     "\nBLUE dashed = null median. Red sitting inside grey => ANY gene set gives",
                     "the same r => the coupling is not evidence of a death-specific mechanism."),
    x = "Spearman rho vs the axis (within 6W, n=12)", y = "number of gene sets") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "B_coupling_null.pdf"), p_b, width = 8, height = 7)

# --- C: what the mt axis tracks -----------------------------------------------
p_c <- axis_tracking |>
  dplyr::mutate(panel = stats::reorder(panel, rho_6W),
                sig = ifelse(p_6W < 0.05, "p < 0.05", "ns")) |>
  ggplot2::ggplot(ggplot2::aes(x = rho_6W, y = panel, fill = sig)) +
  ggplot2::geom_col() +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  ggplot2::scale_fill_manual(values = c("p < 0.05" = "#D73027", "ns" = "grey70")) +
  ggplot2::labs(
    title = "What the mt-transcript share actually tracks (within 6W)",
    subtitle = paste("The IEG/dissociation panel leads. Epithelial markers fall.",
                     "In enzymatically dissociated MECs this axis is prep-linked.\nWhether it is",
                     "technical or biological is NOT settled -- and the nuclear/mt anticorrelation",
                     "\nis NOT evidence either way (a real imbalance predicts it too)."),
    x = "Spearman rho vs mt% (within 6W)", y = NULL, fill = NULL) +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "C_what_mt_tracks.pdf"), p_c, width = 7.5, height = 4)

message("Figures written to ", out_dir)

# =============================================================================
# PART G: SAVE
# =============================================================================
axis_out <- list(
  panels                  = panels,
  panel_share             = panel_share,
  asymmetry               = asymmetry,
  arithmetic_refutation   = arithmetic_refutation,
  mt_vs_nuclear_anticorrelation = mt_vs_nuclear_anticorrelation,  # NOT a direction test (withdrawn)
  axis_tracking           = axis_tracking,
  mt_variance_model       = mt_variance_model,
  axis_design             = axis_design,
  imbalance               = tibble::tibble(sample = samples, group = qc$group,
                                           imbalance = imbalance),
  imbalance_decomposition = imbalance_decomposition,
  imbalance_group_means   = imbalance_group_means,
  coupling_null           = coupling_null,
  null_distributions      = null_distributions,
  genotype_untouched      = genotype_untouched,
  overadjust_guard        = overadjust_guard,
  dominant_axis           = dominant_axis,
  ieg_teb_test            = ieg_teb_test,
  ieg_group_means         = ieg_group_means,
  ieg_teb_verdict         = ieg_teb_verdict,
  technical_vs_biological = technical_vs_biological,
  axis_verdict            = axis_verdict,
  coupling_verdict        = coupling_verdict,
  notes = paste(
    "Block B, author challenge 2026-07-16, two questions: (1) these are PURIFIED",
    "MECs so other cell types are CONTAMINATION -- but why would only mtDNA-coded",
    "genes contaminate, shouldn't nuclear be there too? (2) is the mt signal a",
    "stress artifact, and does the mitonuclear imbalance survive? PART A answers",
    "(1) THREE ways: direction (nuclear mito genes FALL with mt%, so mt% is not a",
    "content readout), arithmetic (contamination is ~3-12% of cells; mixing would",
    "need a contaminant mt-share of 187-737% -- impossible), and the author's",
    "ASYMMETRY (nuclear MEC genes track contamination NEGATIVELY = dilution, mt%",
    "POSITIVELY -- mixing cannot give opposite signs). mt% instead tracks a dominant",
    "sample-level axis marked by the canonical dissociation IEG panel (rho +0.87 at",
    "6W; R2 0.23 -> 0.72). PART B: script 24's mtnuc_index is ~73% that axis",
    "(r=-0.85). PART C is THE LOAD-BEARING PART and the test script 25 never ran:",
    "the published 'imbalance couples to pro-death priming at 6W (r=0.75 p=0.005)'",
    "sits at the ~61st percentile of the imbalance's couplings to ALL 884 library",
    "sets (median ~+0.62), permutation p~0.49 -- ANY gene set gives r~0.6. This is",
    "the AP6 logic of script 21 applied to per-sample couplings. CRITICAL: the",
    "specificity failure holds whether the axis is TECHNICAL or BIOLOGICAL -- it is",
    "a property of the coupling, not of the axis's origin, so it does not depend on",
    "the unresolved adjudication in technical_vs_biological. WHAT IS UNTOUCHED: the",
    "axis is genotype-INDEPENDENT (p~0.49) so it cannot bias ANY Myc contrast --",
    "script 32's content claim retains ~95% of its effect after adjusting for",
    "stress+contamination, and Issues #1/#2/#3/#5/#6 are genotype/interaction-based.",
    "WHAT IS EXPOSED: the axis IS time-associated (p~0.0003), so script 24's",
    "mitonuclear narrative, script 22/29's mtDNA time story, and script 25's",
    "per-sample death couplings all need re-examination. CEILINGS: samples are",
    "ENZYMATICALLY DISSOCIATED (not FACS; author 2026-07-17) with NO batch/viability/RIN",
    "metadata on disk, so the IEG panel",
    "is an INFERRED covariate, not a measured one -- technical vs biological is NOT",
    "resolved. n=12 within 6W. The null universe is the 884 GSVA-scored library sets,",
    "which are correlated with each other, so the empirical p is indicative, not",
    "exact -- but the observed coupling sitting at the MEDIAN needs no precision to",
    "read. See docs/2026-07-13_BlockA_revision_walkthrough_and_intro_alignment.md."))
saveRDS(axis_out, here::here("results", "mtdna_axis_and_coupling_null.rds"))
message("Saved results/mtdna_axis_and_coupling_null.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  ax <- readRDS(here::here("results", "mtdna_axis_and_coupling_null.rds"))

  # --- The two headlines ---
  cat(strwrap(ax$axis_verdict, width = 88), sep = "\n")
  cat(strwrap(ax$coupling_verdict, width = 88), sep = "\n")

  # --- PART A: the author's question, answered three ways ---
  ax$asymmetry |> print()              # opposite signs => not mixing
  ax$arithmetic_refutation |> print()  # required contaminant mt-share > 100% => impossible
  ax$mt_vs_nuclear_anticorrelation |> print()   # a FACT, not a direction test (withdrawn)
  ax$dominant_axis |> print()          # WHY everything correlates: PC1 ~43%, all measures load
  ax$ieg_teb_test |> print(); cat(strwrap(ax$ieg_teb_verdict, 88), sep="\n")  # IEGs are not TEB
  ax$axis_tracking |> print()          # what mt% DOES track: the IEG panel leads
  ax$mt_variance_model |> print()      # R2 0.23 -> 0.72 once the axis is added

  # --- PART B: how much of the mitonuclear imbalance is just mt%? ---
  ax$imbalance_decomposition |> print()
  round(ax$imbalance_group_means, 3)   # should reproduce script 24's mtnuc_index

  # --- PART C: THE NULL. Does any published coupling beat an arbitrary gene set? ---
  ax$coupling_null |> as.data.frame() |> print()
  # read `perm_p` and `empirical_percentile`: ~0.5 / ~50th pct = not specific

  # --- PART D: what survives, and the adjudication we cannot make ---
  ax$genotype_untouched |> print()     # Myc content claim intact after adjustment
  ax$axis_design |> print()            # genotype-independent, time-associated
  ax$overadjust_guard |> print()
  ax$technical_vs_biological |> print()

  list.files(here::here("outputs", "mtdna_axis_null"), pattern = "\\.pdf$")
}
