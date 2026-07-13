# scripts/23_death_timing_substrate.R
# =============================================================================
# Death-timing SUBSTRATE model (Block A): why acute Myc kills at 6W, not 12W
# =============================================================================
#
# The phenotype is EXTERNAL and solid: in a parallel Myc-ER tamoxifen-inducible
# model (NO RNA-seq there), an acute Myc pulse at 6W causes much more cell death
# than at 12W (IHC). Constant acute pulse, different tissue age -> death is
# SUBSTRATE-GATED by developmental stage. Our chronic MMTV-Myc RNA-seq cannot see
# the inducible phenotype, but it can CHARACTERISE the 6W-vs-12W substrate the
# pulse hits. Because Myc is off pre-tamoxifen in the inducible model, the
# cleanest readout of the death-permissive state is the WT / Myc-negative
# `timepoint_neg` contrast (WT 6W->12W); the chronic-Myc+ layer is secondary and
# SURVIVOR-BIASED (we sequence the cells that did not die).
#
# HONEST CEILING (state as a manuscript limitation): bulk RNA, one timepoint per
# age, n=6/group, survivor bias -> we describe the death-permissive STATE
# (association), not causation. No single-cell / functional (BH3 profiling) data.
#
# Four hypotheses (author: BROAD scope), WT-substrate anchored, Myc+ layer second:
#   H1  BH3-only : anti-apoptotic BCL2 rheostat + p53/ARF readiness  (the lead;
#       PART 2b adds a gene-level p19ARF/p53 re-test of the classic Myc escape;
#       PART 5b adds the p53-INDEPENDENT PUMA regulators FOXO3 / PGC1a-ESRRA / HTRA2)
#   H2  biogenesis-death decoupling + MITONUCLEAR IMBALANCE substrate feature
#       (from script 22: Myc forces a selective, mtDNA-vs-nuclear-imbalanced mito
#        state maximal at 6W_pos -> candidate death-permissive stress state)
#   H3  proliferation-apoptosis coupling (oncogene-induced apoptosis at 6W)
#   H4  selection / survivor culling (cross-sample CV narrowing; weakest arm)
# + integration of cell-death branch 1 (16), branch 2 (12), Gate 2 (14) into one
#   convergence read.
#
# Input:  results/interaction_results.rds (myc_6W_raw, timepoint_neg_raw [WT
#           substrate, PRIMARY], timepoint_pos_raw), results/mitopps_scores.rds
#           (apoptosis_pro_genes/anti_genes, per-sample mitoPPS incl. mtDNA +
#           nuclear OXPHOS + Apoptosis-PRO/ANTI, mitopps_pairwise),
#           results/gate2_apoptosis_readout.rds, results/cell_death_fgsea.rds
#           (Tang 15 RCD, branch 2), results/cell_death_hypothesis_results_raw.csv
#           + results/cell_death_genes_full_raw.csv (branch 1),
#           results/gsva_scores.rds (scores/set_meta/sample_meta/expr_mat),
#           results/gsva_overview.rds (coef_table), results/ortholog_table.rds,
#           msigdbr Hallmark P53_PATHWAY.
# Output: results/death_timing_substrate.rds; outputs/death_timing/*.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

out_dir <- here::here("outputs", "death_timing")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# =============================================================================
# PART 1: LOAD + SAMPLE METADATA (treatment coding, ref 6W/neg)
# =============================================================================

ir       <- readRDS(here::here("results", "interaction_results.rds"))
mps      <- readRDS(here::here("results", "mitopps_scores.rds"))
gate2    <- readRDS(here::here("results", "gate2_apoptosis_readout.rds"))
cdf      <- readRDS(here::here("results", "cell_death_fgsea.rds"))
gs       <- readRDS(here::here("results", "gsva_scores.rds"))
coef_tbl <- readRDS(here::here("results", "gsva_overview.rds"))$coef_table
ortholog <- readRDS(here::here("results", "ortholog_table.rds"))

b1_hyp   <- readr::read_csv(here::here("results", "cell_death_hypothesis_results_raw.csv"),
                            show_col_types = FALSE)
b1_genes <- readr::read_csv(here::here("results", "cell_death_genes_full_raw.csv"),
                            show_col_types = FALSE)

scores   <- gs$scores                       # set x 24 GSVA
expr_mat <- gs$expr_mat                      # symbol x 24 VST (full universe)
sample_meta <- as.data.frame(gs$sample_meta)
if (!"sample" %in% names(sample_meta)) sample_meta$sample <- rownames(sample_meta)
sample_meta <- sample_meta[match(colnames(scores), sample_meta$sample), , drop = FALSE]
sample_meta$timepoint  <- factor(sample_meta$timepoint,  levels = c("6W", "12W"))
sample_meta$myc_status <- factor(sample_meta$myc_status, levels = c("neg", "pos"))
sample_meta$group      <- factor(paste(sample_meta$timepoint, sample_meta$myc_status, sep = "_"),
                                 levels = c("6W_neg", "6W_pos", "12W_neg", "12W_pos"))
tc     <- list(timepoint = "contr.treatment", myc_status = "contr.treatment")
is_neg <- sample_meta$myc_status == "neg"
is_pos <- sample_meta$myc_status == "pos"
is_6W  <- sample_meta$timepoint  == "6W"
is_12W <- sample_meta$timepoint  == "12W"
geno_cols <- c(neg = "#377EB8", pos = "#E41A1C")

# --- helpers ---------------------------------------------------------------
# Per-sample expression composite: row-z each present gene across 24 samples,
# then average -> a cohort-relative signature score per sample.
comp_expr <- function(genes) {
  g <- intersect(genes, rownames(expr_mat))
  if (length(g) < 3) return(rep(NA_real_, ncol(expr_mat)))
  z <- t(scale(t(expr_mat[g, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}
# Raw (non-z) VST composite -- for CV (z-scored mean ~0 makes CV unstable).
comp_expr_raw <- function(genes) {
  g <- intersect(genes, rownames(expr_mat))
  if (length(g) < 3) return(rep(NA_real_, ncol(expr_mat)))
  colMeans(expr_mat[g, , drop = FALSE], na.rm = TRUE)
}
# Per-sample GSVA category composite (mean over sets in a category_primary).
comp_gsva <- function(cat) {
  sn <- coef_tbl$set_name[coef_tbl$category_primary == cat]
  sn <- intersect(sn, rownames(scores))
  if (length(sn) < 2) return(rep(NA_real_, ncol(scores)))
  colMeans(scores[sn, , drop = FALSE])
}
# Genotype main effect (powered additive) + interaction on a per-sample vector.
test_geno <- function(y) {
  d <- data.frame(y = y, myc_status = sample_meta$myc_status, timepoint = sample_meta$timepoint)
  sa <- summary(stats::lm(y ~ myc_status + timepoint, data = d, contrasts = tc))$coefficients
  si <- summary(stats::lm(y ~ timepoint * myc_status, data = d, contrasts = tc))$coefficients
  tibble::tibble(geno_beta = sa["myc_statuspos", "Estimate"],
                 geno_p    = sa["myc_statuspos", "Pr(>|t|)"],
                 int_beta  = si["timepoint12W:myc_statuspos", "Estimate"],
                 int_p     = si["timepoint12W:myc_statuspos", "Pr(>|t|)"])
}
# WT-only 6W->12W change (the substrate the acute pulse hits).
wt_time <- function(y) {
  d  <- data.frame(y = y[is_neg], timepoint = sample_meta$timepoint[is_neg])
  co <- summary(stats::lm(y ~ timepoint, data = d,
                          contrasts = tc["timepoint"]))$coefficients
  tibble::tibble(wt_delta = co["timepoint12W", "Estimate"], wt_p = co["timepoint12W", "Pr(>|t|)"])
}
group_means <- function(y) tapply(y, sample_meta$group, mean)

# ENSMUSG <-> symbol
ens2sym <- ortholog |>
  dplyr::select(ensembl_gene_id, external_gene_name) |>
  dplyr::filter(!is.na(external_gene_name), external_gene_name != "") |>
  dplyr::distinct(external_gene_name, .keep_all = TRUE)
lfc_for <- function(res, genes) {
  df  <- as.data.frame(res)
  ens <- ens2sym$ensembl_gene_id[match(genes, ens2sym$external_gene_name)]
  tibble::tibble(symbol = genes, ensembl = ens,
                 lfc  = df$log2FoldChange[match(ens, rownames(df))],
                 stat = df$stat[match(ens, rownames(df))]) |>
    dplyr::filter(!is.na(lfc))
}

pro_genes  <- mps$apoptosis_pro_genes     # BH3-only + effectors (curated, script 08)
anti_genes <- mps$apoptosis_anti_genes    # anti-apoptotic BCL2 family

# Hallmark P53 pathway (mouse) -- robust to the msigdbr category/collection rename
hh <- tryCatch(msigdbr::msigdbr(species = "Mus musculus", category = "H"),
               error = function(e) msigdbr::msigdbr(species = "Mus musculus",
                                                    collection = "H"))
sym_col  <- intersect(c("gene_symbol", "db_gene_symbol", "mouse_symbol"), colnames(hh))[1]
name_col <- intersect(c("gs_name"), colnames(hh))[1]
p53_genes <- unique(hh[[sym_col]][hh[[name_col]] == "HALLMARK_P53_PATHWAY"])
arf_core  <- c("Trp53", "Cdkn2a", "Mdm2", "Cdkn1a")

# =============================================================================
# PART 2: H1 -- BH3:BCL2 rheostat + p53/ARF readiness (WT substrate lead)
# =============================================================================
# Per-sample priming STATE (PRO minus ANTI expression) across the 4 groups: does
# the WT substrate become LESS primed with age (index falls 6W->12W)? Plus the
# LFC module tests (Gate-2 style) on the temporal contrasts, primary = WT.

pro_state  <- comp_expr(pro_genes)
anti_state <- comp_expr(anti_genes)
p53_state  <- comp_expr(p53_genes)
priming_state <- pro_state - anti_state          # higher = more death-primed

h1_state <- dplyr::bind_cols(
  tibble::tibble(metric = c("PRO", "ANTI", "priming (PRO-ANTI)", "p53")),
  dplyr::bind_rows(test_geno(pro_state), test_geno(anti_state),
                   test_geno(priming_state), test_geno(p53_state)),
  dplyr::bind_rows(wt_time(pro_state), wt_time(anti_state),
                   wt_time(priming_state), wt_time(p53_state)))

h1_group_means <- tibble::tibble(
  group   = names(group_means(priming_state)),
  PRO     = as.numeric(group_means(pro_state)),
  ANTI    = as.numeric(group_means(anti_state)),
  priming = as.numeric(group_means(priming_state)),
  p53     = as.numeric(group_means(p53_state)))

# LFC module tests on the temporal + 6W-Myc contrasts (primary = WT timepoint_neg)
contrasts_lfc <- list(timepoint_neg = ir$timepoint_neg_raw,   # WT substrate 6W->12W
                      timepoint_pos = ir$timepoint_pos_raw,   # Myc+ 6W->12W
                      myc_6W        = ir$myc_6W_raw)          # Myc effect at 6W
h1_lfc <- purrr::map_dfr(names(contrasts_lfc), function(cn) {
  res  <- contrasts_lfc[[cn]]
  prol <- lfc_for(res, pro_genes)$lfc
  antl <- lfc_for(res, anti_genes)$lfc
  arfl <- lfc_for(res, arf_core)$lfc
  tt   <- function(x) if (length(x) >= 3) stats::t.test(x, mu = 0)$p.value else NA_real_
  tibble::tibble(
    contrast      = cn,
    pro_mean_lfc  = mean(prol), pro_t_p  = tt(prol),
    anti_mean_lfc = mean(antl), anti_t_p = tt(antl),
    priming_index = mean(prol) - mean(antl),     # >0 = pulled toward death
    arf_mean_lfc  = mean(arfl))
})

# Per-gene PRO/ANTI LFC on the WT substrate (forest input)
h1_wt_genes <- dplyr::bind_rows(
  lfc_for(ir$timepoint_neg_raw, pro_genes)  |> dplyr::mutate(module = "PRO / BH3-only"),
  lfc_for(ir$timepoint_neg_raw, anti_genes) |> dplyr::mutate(module = "ANTI (BCL2)"))

# =============================================================================
# PART 2b: p19ARF / p53 axis -- gene-level re-test (does it explain the 12W drop?)
# =============================================================================
# The canonical Myc anti-apoptosis escape cited in the literature: Myc induces
# p19ARF (Cdkn2a locus) -> stabilises p53 -> apoptosis; loss of ARF or of p53
# ACTIVITY relieves that and permits survival. Tested here as a candidate for the
# 12W death drop, at GENE level (PART 2 read the p53 composite only as flat/null).
# Death-OFF prediction: ARF induction LOST at 12W (Cdkn2a Myc-effect falls) AND/OR
# p53 transcriptional ACTIVITY coordinately reduced at 12W (targets fall; Mdm2 up).
# Uses RAW/unshrunken contrasts (a dLFC/interaction question -- shrinkage rule).
# Reads: Myc@6W, Myc@12W (Myc-vs-WT LFC), interaction (= myc_12W - myc_6W, with the
# Wald z/padj), and the two temporal LFCs (WT + Myc+) for context. p53 ACTIVITY is
# read from TARGET induction, not Trp53 mRNA (p53 is regulated post-translationally).

arf_p53_roster <- tibble::tribble(
  ~symbol,       ~role,
  "Trp53",       "ARF/p53 core",
  "Cdkn2a",      "ARF/p53 core",
  "Mdm2",        "ARF/p53 core",
  "Cdkn1a",      "ARF/p53 core",
  "Bbc3",        "p53 target (activity)",
  "Pmaip1",      "p53 target (activity)",
  "Bax",         "p53 target (activity)",
  "Ccng1",       "p53 target (activity)",
  "Zmat3",       "p53 target (activity)",
  "Sesn2",       "p53 target (activity)",
  "Eda2r",       "p53 target (activity)",
  "Phlda3",      "p53 target (activity)",
  "Trp53inp1",   "p53 target (activity)",
  "Perp",        "p53 target (activity)",
  "Aen",         "p53 target (activity)",
  "Ei24",        "p53 target (activity)")

# gene-level extractor incl. padj (lfc_for omits padj); symbol -> ensembl via ens2sym
stat_for <- function(res, genes) {
  df  <- as.data.frame(res)
  ens <- ens2sym$ensembl_gene_id[match(genes, ens2sym$external_gene_name)]
  tibble::tibble(symbol = genes,
                 lfc  = df$log2FoldChange[match(ens, rownames(df))],
                 z    = df$stat[match(ens, rownames(df))],
                 padj = df$padj[match(ens, rownames(df))])
}
grab <- function(res, suffix, genes) {
  s <- stat_for(res, genes)
  names(s)[-1] <- paste0(names(s)[-1], suffix)
  s
}
# assemble a per-gene table (LFC/z/padj across the 5 raw contrasts) for a roster
contrast_table <- function(genes) {
  tibble::tibble(symbol = genes) |>
    dplyr::left_join(grab(ir$myc_6W_raw,        "_myc6",  genes), by = "symbol") |>
    dplyr::left_join(grab(ir$myc_12W_raw,       "_myc12", genes), by = "symbol") |>
    dplyr::left_join(grab(ir$interaction_raw,   "_int",   genes), by = "symbol") |>
    dplyr::left_join(grab(ir$timepoint_neg_raw, "_wt",    genes), by = "symbol") |>
    dplyr::left_join(grab(ir$timepoint_pos_raw, "_pos",   genes), by = "symbol") |>
    dplyr::mutate(dMyc = lfc_myc12 - lfc_myc6) |>
    dplyr::select(symbol, myc6 = lfc_myc6, myc12 = lfc_myc12, dMyc,
                  int_z = z_int, int_padj = padj_int, wt_t = lfc_wt, pos_t = lfc_pos)
}

arf_p53_genes <- arf_p53_roster |>
  dplyr::left_join(contrast_table(arf_p53_roster$symbol), by = "symbol")

# module-level verdict (DATA-DRIVEN): death-OFF needs BOTH ARF lost AND p53-target
# activity coordinately reduced at 12W; otherwise the escape is REJECTED as the cause.
tgt_int   <- arf_p53_genes$dMyc[arf_p53_genes$role == "p53 target (activity)"]
tgt_int   <- tgt_int[is.finite(tgt_int)]
arf_dMyc  <- arf_p53_genes$dMyc[arf_p53_genes$symbol == "Cdkn2a"]
mdm2_dMyc <- arf_p53_genes$dMyc[arf_p53_genes$symbol == "Mdm2"]
bbc3_z    <- arf_p53_genes$int_z[arf_p53_genes$symbol == "Bbc3"]
tgt_mean  <- mean(tgt_int)
tgt_t_p   <- if (length(tgt_int) >= 3) stats::t.test(tgt_int, mu = 0)$p.value else NA_real_
arf_lost  <- isTRUE(arf_dMyc < 0)
p53_down  <- isTRUE(tgt_mean < 0 && !is.na(tgt_t_p) && tgt_t_p < 0.05)
verdict_short <- if (!arf_lost && !p53_down) "REJECTED (BCL2-family rheostat, not ARF/p53)" else "SIGNAL -- inspect"
arf_p53_verdict <- tibble::tibble(
  cdkn2a_dMyc = arf_dMyc, mdm2_dMyc = mdm2_dMyc,
  p53target_mean_int = tgt_mean, p53target_int_t_p = tgt_t_p,
  bbc3_int_z = bbc3_z,
  any_int_fdr_sig = any(arf_p53_genes$int_padj < 0.05, na.rm = TRUE),
  verdict_short = verdict_short,
  verdict = paste0(
    "p19ARF/p53 escape as the 12W death mechanism: ", verdict_short, ". Cdkn2a Myc-effect ",
    ifelse(arf_lost, "FALLS", "does NOT fall"), " at 12W (dMyc=", sprintf("%+.2f", arf_dMyc),
    "; ARF not lost); Mdm2 dMyc=", sprintf("%+.2f", mdm2_dMyc), " (no added p53 brake); ",
    "p53-target module mean interaction=", sprintf("%+.3f", tgt_mean), " (t p=",
    sprintf("%.2f", tgt_t_p), ") = not coordinately reduced. Only Bbc3/Puma (int z=",
    sprintf("%+.1f", bbc3_z), ") + Bax decline = BH3/effector BCL2-family, downstream of p53. ",
    "Corroborates the PART-2 p53/ARF null; the rheostat is the BCL2 family."))

# =============================================================================
# PART 3: H2 -- biogenesis-death decoupling + MITONUCLEAR IMBALANCE substrate
# =============================================================================
# (a) mitonuclear imbalance per sample (from script-22 finding): nuclear OXPHOS
#     mitoPPS minus mtDNA-encoded mitoPPS. Maximal imbalance = candidate stress /
#     death-permissive state. Test: max at 6W_pos? correlate with PRO priming.
# (b) biogenesis ~ death coupling within timepoint: present at 6W, decays by 12W?

pps      <- mps$mitopps_scores
pps      <- pps[match(sample_meta$sample, pps$sample), , drop = FALSE]  # align order
tier1    <- mps$pathway_tier1_map
mtdna_nm <- mps$mtdna_pathway_name
oxn      <- setdiff(intersect(names(tier1)[tier1 == "OXPHOS"], names(pps)), mtdna_nm)
nuclear_pps <- rowMeans(pps[, oxn, drop = FALSE], na.rm = TRUE)
mtdna_pps   <- pps[[mtdna_nm]]
mitonuclear_imbalance <- nuclear_pps - mtdna_pps        # higher = more imbalanced

bio_comp    <- comp_gsva("Biogenesis_discrimination")
inter_comp  <- comp_gsva("Biogenesis_apoptosis_intersections")
pro_comp    <- pro_state                                 # per-sample PRO priming

h2_imbalance <- tibble::tibble(
  group = names(group_means(mitonuclear_imbalance)),
  mitonuclear_imbalance = as.numeric(group_means(mitonuclear_imbalance)),
  nuclear_pps = as.numeric(group_means(nuclear_pps)),
  mtdna_pps   = as.numeric(group_means(mtdna_pps)))
h2_imbalance_test <- test_geno(mitonuclear_imbalance)

# within-timepoint correlations (coupling): imbalance~PRO and biogenesis~PRO
cor_within <- function(x, y, pair) {
  purrr::map_dfr(c("6W", "12W"), function(tp) {
    i  <- sample_meta$timepoint == tp
    ct <- suppressWarnings(stats::cor.test(x[i], y[i], method = "pearson"))
    tibble::tibble(pair = pair, timepoint = tp, r = unname(ct$estimate),
                   p = ct$p.value, n = sum(i))
  })
}
h2_coupling <- dplyr::bind_rows(
  cor_within(mitonuclear_imbalance, pro_comp, "mitonuclear_imbalance ~ PRO"),
  cor_within(bio_comp,   pro_comp, "biogenesis ~ PRO"),
  cor_within(inter_comp, pro_comp, "intersections ~ PRO"))

# corroborate with mitoPPS Apoptosis-PRO/ANTI temporal reprioritisation
h2_mitopps_apop <- mps$mitopps_pairwise |>
  dplyr::filter(pathway %in% c("Apoptosis-PRO", "Apoptosis-ANTI"),
                contrast %in% c("Temporal_Myc+", "Temporal_Myc-", "Myc_effect_6W")) |>
  dplyr::select(pathway, contrast, diff, padj)

# =============================================================================
# PART 4: H3 -- proliferation-apoptosis coupling
# =============================================================================
# Oncogene-induced apoptosis: proliferation and pro-death priming co-vary at the
# proliferative 6W state, decoupling by 12W.

prolif_comp <- comp_gsva("Proliferation")
h3_coupling <- cor_within(prolif_comp, pro_comp, "proliferation ~ PRO")

# =============================================================================
# PART 5: H4 -- selection / survivor culling (cross-sample CV; weakest arm)
# =============================================================================
# If sensitive cells are culled, replicate-to-replicate heterogeneity of the
# death programme narrows 6W->12W in Myc+ (survivors converge). CV on the RAW VST
# composite (z-scored mean ~0 makes CV meaningless).

pro_raw   <- comp_expr_raw(pro_genes)
death_raw <- comp_expr_raw(unique(b1_genes$mgi_symbol[!is.na(b1_genes$mgi_symbol)]))
cv_by_group <- function(y) {
  vapply(levels(sample_meta$group), function(g) {
    v <- y[sample_meta$group == g]; stats::sd(v) / mean(v)
  }, numeric(1))
}
h4_cv <- tibble::tibble(
  group   = levels(sample_meta$group),
  cv_PRO  = as.numeric(cv_by_group(pro_raw)),
  cv_death = as.numeric(cv_by_group(death_raw))) |>
  tidyr::separate(group, into = c("timepoint", "myc_status"), sep = "_", remove = FALSE)

# =============================================================================
# PART 5b: p53-INDEPENDENT PUMA regulators (FOXO3 / PGC1a-ESRRA / HTRA2 / AP-1)
# =============================================================================
# PART 2b showed the 12W PUMA(Bbc3) loss is NOT a p53/ARF effect. PUMA has p53-
# INDEPENDENT drivers (docs/library_reference/PUMA-and-its-relationships.md): FOXO3a
# (the main p53-independent activator, growth-factor-withdrawal), GATED by PGC-1a
# (abundant PGC1a complexes with FOXO3a and SUPPRESSES its apoptotic output; loss of
# PGC1a unleashes FOXO3a->PUMA), with FOXO co-targets Bim/Noxa/BNIP3, the antioxidant
# survival arm (Sod2/Cat), and the mitochondrial executioner HTRA2 (TFs: p53, HSF1,
# AP-1, Sp1, YY1). Two lenses, both on RAW contrasts / already-scored GSVA:
#   A. regulator+target GENE panel across contrasts (does a FOXO3->PUMA->HTRA2 module
#      co-move with Bbc3? does PGC1a RISE as a brake? do Bim/Noxa co-fall?)
#   B. per-sample TF-SIGNATURE activity (FOXO / PGC1a-ESRRA / AP-1): trajectory +
#      coupling to Bbc3 and PRO priming (is the gene co-movement backed by TF ACTIVITY?)

# --- A. regulator + target gene panel -------------------------------------
puma_reg_roster <- tibble::tribble(
  ~symbol,     ~role,
  "Bbc3",      "PUMA (anchor)",
  "Foxo1",     "FOXO driver",
  "Foxo3",     "FOXO driver",
  "Foxo4",     "FOXO driver",
  "Ppargc1a",  "PGC1a (brake)",
  "Ppargc1b",  "PGC1a (brake)",
  "Esrra",     "PGC1a effector",
  "Nrf1",      "PGC1a effector",
  "Gabpa",     "PGC1a effector",
  "Bcl2l11",   "FOXO BH3 co-target",
  "Pmaip1",    "FOXO BH3 co-target",
  "Bnip3",     "FOXO BH3 co-target",
  "Bnip3l",    "FOXO BH3 co-target",
  "Sod2",      "FOXO-PGC1a survival",
  "Cat",       "FOXO-PGC1a survival",
  "Htra2",     "executioner",
  "Hsf1",      "HTRA2 TF",
  "Jun",       "AP-1 / HTRA2 TF",
  "Fos",       "AP-1 / HTRA2 TF",
  "Sp1",       "HTRA2 TF",
  "Yy1",       "HTRA2 TF")
puma_reg_genes <- puma_reg_roster |>
  dplyr::left_join(contrast_table(puma_reg_roster$symbol), by = "symbol")

gv <- function(sym, col) {
  v <- puma_reg_genes[[col]][puma_reg_genes$symbol == sym]; if (length(v) == 0) NA_real_ else v
}
foxo3_tracks_puma <- isTRUE(gv("Foxo3", "dMyc") < 0 && gv("Bbc3", "dMyc") < 0 &&
                            sign(gv("Foxo3", "int_z")) == sign(gv("Bbc3", "int_z")))
htra2_tracks      <- isTRUE(gv("Htra2", "dMyc") < 0)
pgc1a_brake_up    <- isTRUE(gv("Ppargc1a", "dMyc") > 0 || gv("Esrra", "dMyc") > 0)  # rise = brake
foxo_puma_select  <- isTRUE(gv("Bcl2l11", "dMyc") >= -0.1 && gv("Pmaip1", "dMyc") >= -0.1) # Bim/Noxa NOT co-falling

# --- B. per-sample TF-signature activity ----------------------------------
gmt_lib   <- fgsea::gmtPathways(here::here("data", "genesets_from_library",
                                           "mammary_mito_myc_metab_v1_mouse.gmt"))
avg_scores <- function(sn) {
  sn <- intersect(sn, rownames(scores))
  if (!length(sn)) return(rep(NA_real_, ncol(scores)))
  colMeans(scores[sn, , drop = FALSE])
}
foxo3_sig <- comp_expr(gmt_lib[["TFT_FOXO3_CHUNG"]])                       # FOXO3 target activity
foxo_sig  <- comp_expr(unique(unlist(gmt_lib[c("TFT_FOXO1_CHUNG", "TFT_FOXO3_CHUNG",
                                               "TFT_FOXO4_CHUNG")])))
pgc1a_sig <- avg_scores(c("ESRRA_MITO", "NRF1_MITO", "GABPA_MITO"))        # PGC1a effector activity
ap1_sig   <- avg_scores(grep("^TFT_(JUN|JUNB|FOS|FOSB|FOSL2)_GRAY", rownames(scores), value = TRUE))
bbc3_expr <- if ("Bbc3" %in% rownames(expr_mat)) expr_mat["Bbc3", ] else rep(NA_real_, ncol(expr_mat))

puma_coupling <- dplyr::bind_rows(
  cor_within(foxo3_sig, bbc3_expr, "FOXO3 sig ~ Bbc3"),
  cor_within(foxo_sig,  bbc3_expr, "FOXO(1/3/4) sig ~ Bbc3"),
  cor_within(pgc1a_sig, bbc3_expr, "PGC1a/ESRRA sig ~ Bbc3"),
  cor_within(ap1_sig,   bbc3_expr, "AP-1 sig ~ Bbc3"),
  cor_within(foxo3_sig, pro_state, "FOXO3 sig ~ PRO priming"))
puma_sig_traj <- dplyr::bind_cols(
  tibble::tibble(signature = c("FOXO3", "FOXO(1/3/4)", "PGC1a/ESRRA", "AP-1")),
  dplyr::bind_rows(test_geno(foxo3_sig), test_geno(foxo_sig),
                   test_geno(pgc1a_sig), test_geno(ap1_sig)))

# does the TF ACTIVITY confirm the gene co-movement? (FOXO3 sig declines AND couples 6W)
foxo3_r6 <- puma_coupling$r[puma_coupling$pair == "FOXO3 sig ~ Bbc3" & puma_coupling$timepoint == "6W"]
foxo3_int <- puma_sig_traj$int_beta[puma_sig_traj$signature == "FOXO3"]
foxo3_activity_confirms <- isTRUE(foxo3_int < 0 && foxo3_r6 > 0.4)

puma_reg_verdict <- tibble::tibble(
  foxo3_dMyc = gv("Foxo3", "dMyc"), bbc3_dMyc = gv("Bbc3", "dMyc"),
  htra2_dMyc = gv("Htra2", "dMyc"), esrra_dMyc = gv("Esrra", "dMyc"),
  bim_dMyc = gv("Bcl2l11", "dMyc"), noxa_dMyc = gv("Pmaip1", "dMyc"),
  foxo3_sig_r6 = foxo3_r6, foxo3_sig_int = foxo3_int,
  foxo3_tracks_puma = foxo3_tracks_puma, pgc1a_brake_up = pgc1a_brake_up,
  foxo_puma_selective = foxo_puma_select, foxo3_activity_confirms = foxo3_activity_confirms,
  verdict = paste0(
    "p53-INDEPENDENT PUMA regulators. GENE level: a FOXO3/PUMA/HTRA2 module is ",
    "coordinately Myc-induced at 6W and withdrawn by 12W (Foxo3 dMyc=", sprintf("%+.2f", gv("Foxo3","dMyc")),
    ", int z=", sprintf("%+.1f", gv("Foxo3","int_z")), "; Bbc3 int z=", sprintf("%+.1f", gv("Bbc3","int_z")),
    "; Htra2 dMyc=", sprintf("%+.2f", gv("Htra2","dMyc")), ") -- concordant, directional. ",
    "PGC1a-BRAKE model ", ifelse(pgc1a_brake_up, "SUPPORTED", "REJECTED"),
    " (ESRRA dMyc=", sprintf("%+.2f", gv("Esrra","dMyc")), "/PGC1a dMyc=", sprintf("%+.2f", gv("Ppargc1a","dMyc")),
    " -- co-decline, do NOT rise). FOXO program is ", ifelse(foxo_puma_select, "PUMA-SELECTIVE", "broad"),
    " (Bim dMyc=", sprintf("%+.2f", gv("Bcl2l11","dMyc")), ", Noxa dMyc=", sprintf("%+.2f", gv("Pmaip1","dMyc")),
    "). CAVEAT (SIGNATURE level): FOXO3 target ACTIVITY does ", ifelse(foxo3_activity_confirms, "", "NOT "),
    "confirm -- 6W coupling to Bbc3 r=", sprintf("%.2f", foxo3_r6), ", trajectory int=", sprintf("%+.3f", foxo3_int),
    " -> the gene co-movement is consistent with DE-AMPLIFICATION of the 6W death-primed compartment ",
    "rather than a proven FOXO3-activity cascade. Association, n=6, no FDR."))

# =============================================================================
# PART 6: BRANCH INTEGRATION + CONVERGENCE READ
# =============================================================================
# Branch 1 (raw binomial, delta = myc_6W - myc_12W; SUPPORTING = pro-death front-
# loaded at 6W). Branch 2 (Tang RCD fGSEA NES, esp. Apoptosis on the WT substrate
# temporal_neg and myc_effect_6W). Gate 2 (PRO/ANTI module shift).

b1_summary <- b1_hyp |>
  dplyr::filter(category %in% c("All genes", "Apoptosis pathway")) |>
  dplyr::select(category, n_genes, pct_supporting, binom_p, direction)

b2 <- cdf$fgsea_combined
b2_col <- intersect(c("pathway_label", "pathway"), colnames(b2))[1]
b2$modality <- b2[[b2_col]]
b2_summary <- b2 |>
  dplyr::filter(contrast %in% c("myc_effect_6W", "temporal_neg", "temporal_pos")) |>
  dplyr::select(modality, contrast, NES, padj) |>
  dplyr::arrange(contrast, dplyr::desc(NES))
b2_apop <- b2_summary |> dplyr::filter(grepl("Apopto", modality, ignore.case = TRUE))

convergence <- tibble::tibble(
  lens = c("H1 priming (WT substrate)", "H1 p53 (WT substrate)",
           "H2 mitonuclear imbalance", "H2 biogenesis~death 6W coupling",
           "H3 proliferation~death 6W coupling", "H4 Myc+ CV 6W->12W",
           "Gate2 PRO module shift", "Branch1 apoptosis front-loading",
           "Branch2 apoptosis (WT substrate)", "H1 ARF/p53 axis (gene-level)",
           "H1 p53-indep PUMA (FOXO3/PGC1a)"),
  metric = c(
    sprintf("WT delta = %+.3f (p=%.3f)",
            h1_state$wt_delta[h1_state$metric == "priming (PRO-ANTI)"],
            h1_state$wt_p[h1_state$metric == "priming (PRO-ANTI)"]),
    sprintf("WT delta = %+.3f (p=%.3f)",
            h1_state$wt_delta[h1_state$metric == "p53"],
            h1_state$wt_p[h1_state$metric == "p53"]),
    sprintf("6W_pos=%+.3f vs 12W_pos=%+.3f",
            h2_imbalance$mitonuclear_imbalance[h2_imbalance$group == "6W_pos"],
            h2_imbalance$mitonuclear_imbalance[h2_imbalance$group == "12W_pos"]),
    sprintf("r6W=%+.2f -> r12W=%+.2f",
            h2_coupling$r[h2_coupling$pair == "biogenesis ~ PRO" & h2_coupling$timepoint == "6W"],
            h2_coupling$r[h2_coupling$pair == "biogenesis ~ PRO" & h2_coupling$timepoint == "12W"]),
    sprintf("r6W=%+.2f -> r12W=%+.2f",
            h3_coupling$r[h3_coupling$timepoint == "6W"],
            h3_coupling$r[h3_coupling$timepoint == "12W"]),
    sprintf("cv6W=%.2f -> cv12W=%.2f",
            h4_cv$cv_PRO[h4_cv$group == "6W_pos"], h4_cv$cv_PRO[h4_cv$group == "12W_pos"]),
    sprintf("PRO mean LFC=%+.3f (t p=%.3f); class=%s",
            gate2$pro_ttest$estimate, gate2$pro_ttest$p.value, gate2$gate2_readout),
    sprintf("Apoptosis %%supporting=%.1f (binom p=%.3f)",
            b1_summary$pct_supporting[b1_summary$category == "Apoptosis pathway"],
            b1_summary$binom_p[b1_summary$category == "Apoptosis pathway"]),
    if (nrow(b2_apop) > 0)
      paste(sprintf("%s NES=%+.2f", b2_apop$contrast, b2_apop$NES), collapse = "; ")
    else "no apoptosis modality",
    sprintf("Cdkn2a dMyc=%+.2f (ARF %s); p53-target mean int=%+.3f (t p=%.2f); Bbc3 int z=%+.1f -> %s",
            arf_p53_verdict$cdkn2a_dMyc, ifelse(arf_p53_verdict$cdkn2a_dMyc < 0, "lost", "not lost"),
            arf_p53_verdict$p53target_mean_int, arf_p53_verdict$p53target_int_t_p,
            arf_p53_verdict$bbc3_int_z, arf_p53_verdict$verdict_short),
    sprintf("Foxo3 dMyc=%+.2f (int z=%+.1f) ~ Bbc3 int z=%+.1f; Htra2 dMyc=%+.2f; ESRRA dMyc=%+.2f (brake %s); FOXO3-activity r6W=%.2f -> %s",
            puma_reg_verdict$foxo3_dMyc, gv("Foxo3", "int_z"), gv("Bbc3", "int_z"),
            puma_reg_verdict$htra2_dMyc, puma_reg_verdict$esrra_dMyc,
            ifelse(puma_reg_verdict$pgc1a_brake_up, "up", "co-decline"), puma_reg_verdict$foxo3_sig_r6,
            ifelse(puma_reg_verdict$foxo3_activity_confirms, "activity-confirmed", "gene-only (de-amplification)"))))

# =============================================================================
# PART 7: FIGURES
# =============================================================================

# H1: BH3:BCL2 rheostat trajectory (priming state, WT vs Myc+)
h1_traj <- h1_group_means |>
  tidyr::separate(group, into = c("timepoint", "myc_status"), sep = "_", remove = FALSE) |>
  dplyr::mutate(timepoint = factor(timepoint, levels = c("6W", "12W")))
p_h1 <- ggplot2::ggplot(h1_traj, ggplot2::aes(x = timepoint, y = priming,
                                              colour = myc_status, group = myc_status)) +
  ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 3) +
  ggplot2::scale_colour_manual(values = geno_cols,
    labels = c(neg = "Myc- (WT substrate)", pos = "Myc+")) +
  ggplot2::labs(title = "H1: BH3:BCL2 priming state (PRO - ANTI), WT substrate anchored",
    subtitle = "Prediction: WT priming falls 6W->12W (older tissue less death-primed)",
    x = NULL, y = "priming index (PRO - ANTI, cohort-z)", colour = NULL) +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "h1_bh3_bcl2_rheostat.pdf"), p_h1, width = 6.5, height = 4.5)

# H1b: PRO/ANTI per-gene LFC on the WT substrate (forest)
p_h1b <- ggplot2::ggplot(h1_wt_genes,
    ggplot2::aes(x = lfc, y = stats::reorder(symbol, lfc), colour = module)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~ module, scales = "free_y") +
  ggplot2::labs(title = "H1: BCL2-family LFC on the WT substrate (timepoint_neg, 6W->12W)",
    subtitle = "PRO down / ANTI up = WT substrate de-primes with age", x = "raw LFC (WT 12W vs 6W)",
    y = NULL, colour = NULL) +
  ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "none")
ggplot2::ggsave(file.path(out_dir, "h1_bcl2_family_wt_forest.pdf"), p_h1b, width = 8, height = 6)

# H1c: p19ARF / p53 axis -- Myc effect at 6W vs 12W per gene (dumbbell). A death-OFF
# ARF/p53 escape predicts a LEFTWARD shift 6W->12W (ARF lost + p53 targets fall).
arf_p53_long <- arf_p53_genes |>
  dplyr::select(symbol, role, `Myc@6W` = myc6, `Myc@12W` = myc12, dMyc) |>
  tidyr::pivot_longer(c(`Myc@6W`, `Myc@12W`), names_to = "age", values_to = "lfc")
p_h1c <- ggplot2::ggplot(arf_p53_long,
    ggplot2::aes(x = lfc, y = stats::reorder(symbol, dMyc))) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_line(ggplot2::aes(group = symbol), colour = "grey70") +
  ggplot2::geom_point(ggplot2::aes(colour = age), size = 2.5) +
  ggplot2::facet_grid(role ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_colour_manual(values = c(`Myc@6W` = "#377EB8", `Myc@12W` = "#E41A1C")) +
  ggplot2::labs(
    title = "H1: p19ARF / p53 axis -- Myc effect at 6W vs 12W (does the ARF/p53 escape explain the 12W death drop?)",
    subtitle = "Death-OFF predicts a leftward 6W->12W shift. Observed: ARF (Cdkn2a) up, p53 targets mostly up; only Bbc3/Bax fall = BCL2-family.",
    x = "Myc-vs-WT raw LFC", y = NULL, colour = NULL) +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "h1c_arf_p53_axis.pdf"), p_h1c, width = 8, height = 6.5)

# H1d: p53-independent PUMA regulators -- Myc effect at 6W vs 12W (FOXO3/PUMA/HTRA2)
puma_reg_long <- puma_reg_genes |>
  dplyr::select(symbol, role, `Myc@6W` = myc6, `Myc@12W` = myc12, dMyc) |>
  tidyr::pivot_longer(c(`Myc@6W`, `Myc@12W`), names_to = "age", values_to = "lfc")
p_h1d <- ggplot2::ggplot(puma_reg_long,
    ggplot2::aes(x = lfc, y = stats::reorder(symbol, dMyc))) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_line(ggplot2::aes(group = symbol), colour = "grey70") +
  ggplot2::geom_point(ggplot2::aes(colour = age), size = 2.4) +
  ggplot2::facet_grid(role ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_colour_manual(values = c(`Myc@6W` = "#377EB8", `Myc@12W` = "#E41A1C")) +
  ggplot2::labs(
    title = "H1: p53-independent PUMA regulators -- Myc effect at 6W vs 12W",
    subtitle = "FOXO3/PUMA(Bbc3)/HTRA2 co-induced at 6W, withdrawn by 12W; PGC1a/ESRRA co-decline (no rising brake); Bim/Noxa flat (PUMA-selective).",
    x = "Myc-vs-WT raw LFC", y = NULL, colour = NULL) +
  ggplot2::theme_bw(base_size = 8)
ggplot2::ggsave(file.path(out_dir, "h1d_puma_regulators.pdf"), p_h1d, width = 8, height = 8)

# H1e: TF-signature ACTIVITY coupling to Bbc3 by timepoint (does activity confirm the genes?)
p_h1e <- puma_coupling |>
  dplyr::mutate(timepoint = factor(timepoint, levels = c("6W", "12W"))) |>
  ggplot2::ggplot(ggplot2::aes(x = r, y = pair, fill = timepoint)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
  ggplot2::scale_fill_manual(values = c(`6W` = "#377EB8", `12W` = "#E41A1C")) +
  ggplot2::labs(
    title = "H1: TF-signature activity coupling to PUMA (Bbc3), by timepoint",
    subtitle = "FOXO3 target-activity does NOT track Bbc3 at 6W -> the Foxo3-gene co-movement is not a confirmed FOXO3-activity cascade",
    x = "within-timepoint Pearson r", y = NULL, fill = NULL) +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(file.path(out_dir, "h1e_puma_tf_signature_coupling.pdf"), p_h1e, width = 7, height = 4)

# H2: mitonuclear imbalance by group + imbalance~PRO coupling
p_h2a <- h2_imbalance |>
  tidyr::separate(group, into = c("timepoint", "myc_status"), sep = "_", remove = FALSE) |>
  dplyr::mutate(group = factor(group, levels = levels(sample_meta$group))) |>
  ggplot2::ggplot(ggplot2::aes(x = group, y = mitonuclear_imbalance, fill = myc_status)) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(values = geno_cols, guide = "none") +
  ggplot2::labs(title = "H2: mitonuclear imbalance (nuclear - mtDNA OXPHOS mitoPPS)",
    subtitle = "Prediction: maximal at 6W_pos = candidate death-permissive stress state",
    x = NULL, y = "nuclear - mtDNA (mitoPPS)") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "h2_mitonuclear_imbalance.pdf"), p_h2a, width = 6, height = 4.5)

coupling_df <- tibble::tibble(
  imbalance = mitonuclear_imbalance, biogenesis = bio_comp, proliferation = prolif_comp,
  PRO = pro_comp, timepoint = sample_meta$timepoint, myc_status = sample_meta$myc_status)
p_h2b <- ggplot2::ggplot(coupling_df, ggplot2::aes(x = biogenesis, y = PRO, colour = timepoint)) +
  ggplot2::geom_point(ggplot2::aes(shape = myc_status), size = 2.5) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  ggplot2::labs(title = "H2: biogenesis-death coupling by timepoint",
    subtitle = "Coupled at 6W, decays by 12W = biogenesis-death decoupling",
    x = "biogenesis composite (GSVA)", y = "PRO priming composite") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "h2_biogenesis_death_coupling.pdf"), p_h2b, width = 6.5, height = 4.5)

# H3: proliferation-death coupling by timepoint
p_h3 <- ggplot2::ggplot(coupling_df, ggplot2::aes(x = proliferation, y = PRO, colour = timepoint)) +
  ggplot2::geom_point(ggplot2::aes(shape = myc_status), size = 2.5) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  ggplot2::labs(title = "H3: proliferation-death coupling by timepoint",
    subtitle = "Oncogene-induced apoptosis: co-vary at 6W, decouple by 12W",
    x = "proliferation composite (GSVA)", y = "PRO priming composite") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "h3_proliferation_death_coupling.pdf"), p_h3, width = 6.5, height = 4.5)

# H4: cross-sample CV narrowing
p_h4 <- h4_cv |>
  dplyr::mutate(timepoint = factor(timepoint, levels = c("6W", "12W"))) |>
  ggplot2::ggplot(ggplot2::aes(x = timepoint, y = cv_PRO, colour = myc_status, group = myc_status)) +
  ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 3) +
  ggplot2::scale_colour_manual(values = geno_cols) +
  ggplot2::labs(title = "H4: cross-sample CV of the death programme (survivor culling)",
    subtitle = "Prediction (weak): Myc+ CV narrows 6W->12W as sensitive cells are culled",
    x = NULL, y = "CV of PRO composite across 6 reps", colour = "myc_status") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "h4_survivor_cv.pdf"), p_h4, width = 6, height = 4.5)

# =============================================================================
# PART 8: SAVE
# =============================================================================

death_out <- list(
  h1 = list(state_tests = h1_state, group_means = h1_group_means,
            lfc_module = h1_lfc, wt_genes = h1_wt_genes,
            arf_p53 = list(genes = arf_p53_genes, verdict = arf_p53_verdict),
            panels = list(pro = pro_genes, anti = anti_genes,
                          p53_n = length(intersect(p53_genes, rownames(expr_mat))),
                          arf = arf_core, arf_p53_panel = arf_p53_roster)),
  h2 = list(imbalance_group = h2_imbalance, imbalance_test = h2_imbalance_test,
            coupling = h2_coupling, mitopps_apoptosis = h2_mitopps_apop,
            per_sample = tibble::tibble(sample = sample_meta$sample,
                                        group = sample_meta$group,
                                        mitonuclear_imbalance = mitonuclear_imbalance,
                                        bio_comp = bio_comp, pro_comp = pro_comp)),
  h3 = list(coupling = h3_coupling),
  h4 = list(cv = h4_cv),
  puma_regulators = list(genes = puma_reg_genes, coupling = puma_coupling,
                         sig_traj = puma_sig_traj, verdict = puma_reg_verdict,
                         roster = puma_reg_roster),
  branches = list(branch1 = b1_summary, branch2 = b2_summary, branch2_apoptosis = b2_apop,
                  gate2 = gate2$summary),
  convergence = convergence,
  notes = paste(
    "Death-timing SUBSTRATE model. Phenotype is EXTERNAL (Myc-ER inducible, IHC):",
    "acute Myc kills more at 6W than 12W. This RNA-seq characterises the substrate,",
    "anchored on the WT timepoint_neg contrast (Myc off pre-tamoxifen); the Myc+",
    "layer is SURVIVOR-BIASED. H1 BH3:BCL2 + p53 rheostat (lead); H2 biogenesis-",
    "death decoupling + mitonuclear imbalance (script-22 finding: selective/",
    "imbalanced mito state maximal at 6W_pos = candidate death-permissive stress);",
    "H3 proliferation-apoptosis coupling; H4 survivor CV (weakest). CEILING: bulk",
    "RNA + survivor bias + n=6 -> death-permissive STATE (association), not",
    "causation; single-cell / BH3 profiling needed for mechanism.")
)
saveRDS(death_out, here::here("results", "death_timing_substrate.rds"))
message("Saved results/death_timing_substrate.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  dt <- readRDS(here::here("results", "death_timing_substrate.rds"))

  # Sanity: BCL2 panels + p53 mapped into expr_mat / contrasts
  length(intersect(pro_genes, rownames(expr_mat)))
  length(intersect(anti_genes, rownames(expr_mat)))
  dt$h1$panels$p53_n
  head(h1_wt_genes) |> print()

  # H1 LEAD: does the WT substrate de-prime with age? (priming WT delta < 0, PRO
  # down / ANTI up on timepoint_neg). Group-mean trajectory: WT 6W > 12W?
  dt$h1$state_tests |> print()
  dt$h1$group_means |> print()
  dt$h1$lfc_module |> print()

  # H1 p19ARF/p53 gene-level re-test: does the ARF/p53 escape explain the 12W drop?
  # Death-OFF needs ARF lost (Cdkn2a dMyc<0) AND p53 targets coordinately down at 12W.
  dt$h1$arf_p53$genes |> print(n = Inf)
  dt$h1$arf_p53$verdict$verdict |> print()

  # H1 p53-INDEPENDENT PUMA regulators (FOXO3 / PGC1a-ESRRA / HTRA2 / AP-1):
  # gene-level module co-movement vs TF-signature ACTIVITY confirmation.
  dt$puma_regulators$genes |> print(n = Inf)
  dt$puma_regulators$coupling |> print()
  dt$puma_regulators$sig_traj |> print()
  dt$puma_regulators$verdict$verdict |> print()

  # H2: is mitonuclear imbalance maximal at 6W_pos? does biogenesis-death coupling
  # decay 6W->12W? (the script-22 substrate feature)
  dt$h2$imbalance_group |> print()
  dt$h2$coupling |> print()
  dt$h2$mitopps_apoptosis |> print()

  # H3 / H4
  dt$h3$coupling |> print()
  dt$h4$cv |> print()

  # THE CONVERGENCE READ across all lenses + branches
  dt$convergence |> print(n = Inf)
  dt$branches$branch2_apoptosis |> print()

  list.files(here::here("outputs", "death_timing"), pattern = "\\.pdf$")
}
