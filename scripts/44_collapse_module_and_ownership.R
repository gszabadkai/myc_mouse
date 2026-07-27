# =============================================================================
# 44_collapse_module_and_ownership.R
# -----------------------------------------------------------------------------
# CAN WE DEFINE AN OXPHOS/APOPTOSIS MODULE? The proposed title claims "a
# PGC-1a/Nrf1 regulated common OXPHOS/apoptosis module". This script asks whether
# such a module exists, in three places, and answers a fourth question the first
# three raise.
#
# PART A -- OWNERSHIP. Which mitochondrial regulon actually carries the
#   pro-apoptotic genes? Hypergeometric overlap of MITOCARTA_APOPTOSIS_PRO/_ANTI
#   with the PGC1a axis (CORE_MITO), MYC_MITO, MYC_SPECIFIC_MITO and
#   DEVELOPMENTAL_MITO. Plus the base-rate control the Gray/ChEA shortlist needs:
#   with 25 PRO genes in ~1140 MitoCarta genes, ANY regulon picks up 1-2 by
#   chance, so "the biogenesis TFs also touch pro-apoptotic effectors" has to be
#   read against that rate before it can mean anything.
#
# PART B -- THE WILD-TYPE EXPRESSION TEST, on BOTH rulers. Does the apoptotic arm
#   de-prioritise together with OXPHOS across 6W->12W in the wild-type gland? The
#   content ruler (set-mean raw LFC, matched-random-set null) and the priority
#   ruler (mitoPPS, from script 40's table) can disagree, and here they must both
#   be shown. Plus the CORE_MITO decomposition: the regulon contains 66% of the
#   OXPHOS subunits, so if the PGC1a axis were withdrawing as a unit the whole
#   300-gene regulon would move, not only its respiratory part.
#
# PART B2 -- THE LINEAGE-SUPPRESSION CARRIER. Gray 2023 (p.7) defines LE as
#   "LOW-EXPRESSING" -- a lineage-program-suppressed state -- NOT low oestrogen.
#   (Our provenance_table.csv glosses MG_HEVSLE_* as "HE vs LE (AP) | adult" and
#   the catalog expands it to "High- vs low-estrogen"; that gloss is WRONG and is
#   recorded here so it cannot propagate again.) ESRRA's and GABPA's mitochondrial
#   programmes are detected in AP_LE -- the lineage-suppressed alveolar-progenitor
#   state the pubescent gland is rich in. So the live alternative to "PGC1a
#   activity declines" is "the cell state that CARRIES the PGC1a-driven OXPHOS
#   programme becomes less abundant". Directional prediction, stated before the
#   test: LE-marker sets (_DN) fall, HE-marker sets (_UP) rise. AP is the
#   prediction; HS and BA are the specificity controls.
#
# PART C -- THE COLLAPSE SCAN. The centrepiece. Since the wild-type gland does not
#   move its apoptotic transcripts at all, a module cannot be defined by
#   co-expression -- so define it from the PHENOMENON instead: which genes lose
#   their Myc-inducibility between 6W and 12W the way PUMA does, and are they a
#   set or is PUMA a solo? Under the global x0.55 attenuation EVERY Myc-induced
#   gene has a negative interaction (script 34 PART G's trap), so the raw
#   interaction is the wrong statistic; rank on the RESIDUAL FROM THE GLOBAL RATE.
#
# PART D -- MECHANISM CANDIDATES. For the routes by which PUMA priming could
#   depend on the respiratory state: does each candidate group collapse like PUMA?
#   Ranking only.
#
# PRE-SPECIFICATION, and why it is not circularity. `Bbc3` (PUMA) and `Bcl2l11`
# (BIM) are named in advance because the PGC1a perturbation induces both at
# protein level (author's westerns). Those experiments are INDEPENDENT of this
# RNA-seq, so using them to fix which genes and sets are tested is
# pre-registration, not circularity -- and at n=24 it is the strongest available
# position. `Bcl2l11` is read as a pre-specified gene DESPITE its non-significant
# 6W effect (p 0.46), with that weakness stated: script 42's correction #1 still
# applies, you cannot lose an effect you never had, so its retention is reported
# but never headlined.
#
# WHAT THIS CANNOT DO, stated before the results. `batch = timepoint` (CLAUDE.md):
# every wild-type temporal statement is DESCRIBED, not claimed. PART C is a
# genome-wide RANKING at n=24 -- a collapse module here is a hypothesis about
# co-dependence, not a demonstration of co-regulation, and it needs the PGC1a
# perturbation transcriptome to become one. PART A is set membership in a curated
# resource: it constrains what a TITLE can claim, it does not establish regulation
# in these mice.
#
# Reads : results/interaction_results.rds, results/dds_int_run.rds,
#         results/combined_df_annotated.rds, results/gene_sets_list.rds,
#         results/background_vs_myc.rds        (priority ruler + global rate),
#         results/attenuation_decomposition.rds (positive control),
#         data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt,
#         data/genesets_from_library/gray_chea_mito_tf_shortlist.csv
#           (the TRACKED copy -- script 42 reads the untracked docs/ duplicate;
#            both are md5-identical, this one is the library snapshot),
#         functions/reconcile_gene_symbols.R (MANDATORY -- vintage-aware membership)
# Writes: results/collapse_module_ownership.rds
#
# RUNTIME NOTE: PART C reconciles ~900 pathways to Ensembl one set at a time so
# that renamed symbols are not silently dropped (the 2026-07-24 fix). That costs a
# couple of minutes and prints progress.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

set.seed(1)
NSET        <- 2000L   # matched-random-set draws
NBIN        <- 20L     # baseMean bins for the matched null
GLOBAL_RATE <- 0.55    # script 40's rescaling slope; PART C re-fits it independently
PADJ6       <- 0.1     # "a real 6W Myc effect"
LFC6_FLOOR  <- 0.2     # script 42's floor: a retention is a ratio of noisy quantities
BM_FLOOR    <- 20      # script 42's expression floor
OX_WT_REF   <- -0.2548 # Issue #4's OXPHOS-subunit wild-type value (positive control)

# =============================================================================
# PART 0: LOAD
# =============================================================================
ir   <- readRDS(here::here("results", "interaction_results.rds"))
dds  <- readRDS(here::here("results", "dds_int_run.rds"))
cdf  <- readRDS(here::here("results", "combined_df_annotated.rds"))
ap6  <- readRDS(here::here("results", "gene_sets_list.rds"))
bgm  <- readRDS(here::here("results", "background_vs_myc.rds"))
gmt  <- fgsea::gmtPathways(here::here("data", "genesets_from_library",
                                      "mammary_mito_myc_metab_v1_mouse.gmt"))
shortlist <- readr::read_csv(
  here::here("data", "genesets_from_library", "gray_chea_mito_tf_shortlist.csv"),
  show_col_types = FALSE)

CONTRASTS <- c("myc_6W_raw", "myc_12W_raw", "timepoint_neg_raw", "timepoint_pos_raw",
               "interaction_raw")
D  <- lapply(ir[CONTRASTS], as.data.frame)      # DESeqResults are S4; coerce once
nc <- DESeq2::counts(dds, normalized = TRUE)
L  <- log2(nc + 1)
sm <- as.data.frame(SummarizedExperiment::colData(dds))
sm$tp  <- stats::relevel(as.factor(sm$timepoint),  "6W")
sm$myc <- stats::relevel(as.factor(sm$myc_status), "neg")

universe_all <- rownames(D$myc_6W_raw)
V   <- function(k) stats::setNames(D[[k]]$log2FoldChange, rownames(D[[k]]))
m6  <- V("myc_6W_raw");        m12 <- V("myc_12W_raw")
tn  <- V("timepoint_neg_raw"); tpz <- V("timepoint_pos_raw")
bm  <- stats::setNames(D$myc_6W_raw$baseMean, rownames(D$myc_6W_raw))

# set membership ALWAYS through the reconciler (CLAUDE.md / the 2026-07-24 fix)
ens_set  <- function(syms) recon_to_ensembl(syms, universe_all)
set_mean <- function(v, e) if (length(e)) mean(v[e], na.rm = TRUE) else NA_real_

sym2ens <- cdf[!is.na(cdf$mgi_symbol) & !duplicated(cdf$mgi_symbol), c("mgi_symbol", "gene")]
ens_of  <- function(s) sym2ens$gene[match(s, sym2ens$mgi_symbol)]
zrow    <- function(m) t(scale(t(m)))
comp_of <- function(syms) {
  e <- ens_of(syms); e <- e[!is.na(e) & e %in% rownames(L)]
  stopifnot(length(e) >= 3)
  colMeans(zrow(L[e, , drop = FALSE]))
}
set_score <- function(set, min_n = 5L) {
  if (is.null(gmt[[set]])) return(NULL)
  e <- ens_set(gmt[[set]]); e <- e[e %in% rownames(L)]
  if (length(e) < min_n) return(NULL)
  colMeans(zrow(L[e, , drop = FALSE]))
}

# GUARD: every named set must exist before anything uses it. Script 42's failure
# mode was a wrong name giving a silent all-NaN composite.
NEEDED <- c("MITOCARTA_OXPHOS_SUBUNITS", "MITOCARTA_OXPHOS",
            "MITOCARTA_OXPHOS_ASSEMBLY_FACTORS", "MITOCARTA_APOPTOSIS_PRO",
            "MITOCARTA_APOPTOSIS_ANTI", "CORE_MITO", "MYC_MITO",
            "MYC_SPECIFIC_MITO", "DEVELOPMENTAL_MITO",
            "MG_AP_GRAY", "MG_HS_GRAY", "MG_BASAL_GRAY",
            "MG_HEVSLE_AP_GRAY_UP", "MG_HEVSLE_AP_GRAY_DN",
            "MG_HEVSLE_HS_GRAY_UP", "MG_HEVSLE_HS_GRAY_DN",
            "MG_HEVSLE_BA_GRAY_UP", "MG_HEVSLE_BA_GRAY_DN",
            "MG_TEB_VS_DUCTAL_HS_GRAY_UP", "MG_TEB_VS_DUCTAL_AP_GRAY_UP",
            "TFT_ESRRA_GRAY_AP_LE", "TFT_GABPA_GRAY_AP_LE", "TFT_NRF1_GRAY_AP_TEB")
missing_sets <- NEEDED[!NEEDED %in% names(gmt)]
if (length(missing_sets))
  stop("44: set names absent from the GMT -> ", paste(missing_sets, collapse = ", "))

# the matched-random-set null, script 43's idiom -------------------------------
expressed <- names(bm)[is.finite(bm) & bm > 0 & is.finite(tn[names(bm)])]
bin_of    <- cut(rank(bm[expressed], ties.method = "first"), breaks = NBIN, labels = FALSE)
names(bin_of) <- expressed
by_bin    <- split(expressed, bin_of)
draw_matched <- function(e) {
  b <- bin_of[e]; b <- b[!is.na(b)]
  unlist(lapply(split(b, b), function(k)
    sample(by_bin[[as.character(k[1])]], length(k), replace = TRUE)), use.names = FALSE)
}
null_pct <- function(v, e, n = NSET) {
  e <- e[e %in% expressed]
  if (length(e) < 3) return(c(observed = NA_real_, null_median = NA_real_,
                              percentile = NA_real_))
  obs <- set_mean(v, e)
  nul <- vapply(seq_len(n), function(i) set_mean(v, draw_matched(e)), numeric(1))
  c(observed = obs, null_median = stats::median(nul), percentile = 100 * mean(nul < obs))
}

# =============================================================================
# PART A: OWNERSHIP -- which regulon carries the pro-apoptotic genes?
# -----------------------------------------------------------------------------
# CORE_MITO = intersect(union(ESRRA, NRF1, GABPA)_targets, MitoCarta), i.e. the
# PGC1a axis at the DNA level (PGC1a is a coactivator with no motif of its own).
# MYC_MITO = intersect(MYC_targets, MitoCarta). Both come from EXTERNAL regulons
# (DoRothEA / ChIP-Atlas), not from this experiment, so neither is circular with
# our differential expression.
# =============================================================================
mito_universe <- unique(unlist(gmt[grep("^MITOCARTA_", names(gmt))], use.names = FALSE))
apop_sets     <- c("MITOCARTA_APOPTOSIS_PRO", "MITOCARTA_APOPTOSIS_ANTI")
prog_sets     <- c("CORE_MITO", "MYC_MITO", "MYC_SPECIFIC_MITO", "DEVELOPMENTAL_MITO")

ownership <- purrr::map_dfr(prog_sets, function(pg) {
  s <- intersect(gmt[[pg]], mito_universe)
  purrr::map_dfr(apop_sets, function(ap) {
    a <- intersect(gmt[[ap]], mito_universe)
    k <- length(intersect(s, a))
    e <- length(s) * length(a) / length(mito_universe)
    tibble::tibble(
      programme = pg, apoptosis_set = ap,
      n_programme = length(s), n_apoptosis = length(a), overlap = k,
      expected = e, fold = k / e,
      p_enrich  = stats::phyper(k - 1, length(a), length(mito_universe) - length(a),
                                length(s), lower.tail = FALSE),
      p_deplete = stats::phyper(k, length(a), length(mito_universe) - length(a),
                                length(s), lower.tail = TRUE),
      genes = paste(sort(intersect(s, a)), collapse = ", "),
      has_PUMA = "Bbc3" %in% s, has_BIM = "Bcl2l11" %in% s)
  })
})

# --- the base-rate control the shortlist claim needs --------------------------
# With 25 PRO genes among ~1140 MitoCarta genes, the expected PRO fraction of any
# regulon's mito overlap is ~0.022. If the observed rate across 741 TFs matches
# that, "the biogenesis TFs also touch pro-apoptotic effectors" is arithmetic.
n_pro_cat  <- length(intersect(gmt[["MITOCARTA_APOPTOSIS_PRO"]],  mito_universe))
n_anti_cat <- length(intersect(gmt[["MITOCARTA_APOPTOSIS_ANTI"]], mito_universe))
shortlist_baserate <- tibble::tibble(
  n_rows              = nrow(shortlist),
  n_TFs               = dplyr::n_distinct(shortlist$TF),
  pro_per_mito_gene   = sum(shortlist$n_apoptosis_pro, na.rm = TRUE) /
                        sum(shortlist$n_mito_total,    na.rm = TRUE),
  catalog_expectation = n_pro_cat / length(mito_universe),
  cor_size_vs_pro     = stats::cor(shortlist$k, shortlist$n_apoptosis_pro,
                                   use = "pairwise.complete.obs"),
  rows_balance_zero   = sum(shortlist$apop_balance == 0, na.rm = TRUE),
  rows_balance_plus   = sum(shortlist$apop_balance >  0, na.rm = TRUE),
  mean_apop_balance   = mean(shortlist$apop_balance, na.rm = TRUE),
  n_pro_catalog       = n_pro_cat, n_anti_catalog = n_anti_cat,
  pro_anti_prior      = n_pro_cat / n_anti_cat)

# The roster rows, located against the resource's OWN promotion gate.
GATE_FDR <- 0.01; GATE_MITO <- 15L
roster <- c("PPARGC1A", "ESRRA", "NRF1", "GABPA", "MYC", "E2F1", "ESR1")
shortlist_roster <- shortlist |>
  dplyr::filter(TF %in% roster) |>
  dplyr::transmute(TF, context, n_mito_total, n_oxphos,
                   n_apoptosis_pro, n_apoptosis_anti, apop_balance,
                   apoptosis_pro_genes, k, hyperg_p, bh_fdr,
                   passes_promotion_gate = bh_fdr < GATE_FDR & n_mito_total >= GATE_MITO) |>
  dplyr::arrange(TF, dplyr::desc(n_mito_total))

# =============================================================================
# PART B: THE WILD-TYPE EXPRESSION TEST, BOTH RULERS
# -----------------------------------------------------------------------------
# The content ruler is a set-mean raw LFC; the priority ruler is mitoPPS, where a
# uniform dose scaling cancels in the pairwise ratio. They fail in different ways,
# so the apoptosis-versus-OXPHOS question is asked on both.
# =============================================================================
wt_arms <- tibble::tribble(
  ~arm,                    ~set,                                ~ruler_pathway,
  "OXPHOS subunits",       "MITOCARTA_OXPHOS_SUBUNITS",         "OXPHOS subunits",
  "OXPHOS (all)",          "MITOCARTA_OXPHOS",                  "OXPHOS",
  "OXPHOS assembly",       "MITOCARTA_OXPHOS_ASSEMBLY_FACTORS", "OXPHOS assembly factors",
  "Apoptosis-PRO",         "MITOCARTA_APOPTOSIS_PRO",           "Apoptosis-PRO",
  "Apoptosis-ANTI",        "MITOCARTA_APOPTOSIS_ANTI",          "Apoptosis-ANTI",
  "PGC1a axis (CORE_MITO)","CORE_MITO",                         NA_character_,
  "MYC_MITO",              "MYC_MITO",                          NA_character_,
  "MYC_SPECIFIC_MITO",     "MYC_SPECIFIC_MITO",                 NA_character_,
  "DEVELOPMENTAL_MITO",    "DEVELOPMENTAL_MITO",                NA_character_)

ruler_tab <- bgm$ruler
wt_content <- purrr::pmap_dfr(wt_arms, function(arm, set, ruler_pathway) {
  e <- ens_set(gmt[[set]])
  n <- null_pct(tn, e)
  r <- if (!is.na(ruler_pathway) && ruler_pathway %in% ruler_tab$pathway)
         ruler_tab[match(ruler_pathway, ruler_tab$pathway), ] else NULL
  tibble::tibble(
    arm = arm, n_genes = length(e),
    c_wt_time   = set_mean(tn,  e),  c_mycpos_time = set_mean(tpz, e),
    c_myc_6W    = set_mean(m6,  e),  c_myc_12W     = set_mean(m12, e),
    wt_null_median = n[["null_median"]], wt_null_pct = n[["percentile"]],
    p_wt_time   = if (is.null(r)) NA_real_ else r$p_tn,
    p_myc_6W    = if (is.null(r)) NA_real_ else r$p_m6,
    p_myc_6W_padj = if (is.null(r)) NA_real_ else r$p_m6_padj)
})

# POSITIVE CONTROL: OXPHOS subunits on the wild-type axis must reproduce Issue #4.
ox_wt <- wt_content$c_wt_time[wt_content$arm == "OXPHOS subunits"]
if (!is.finite(ox_wt) || abs(ox_wt - OX_WT_REF) > 0.02)
  warning(sprintf("44 PART B: OXPHOS-subunit WT %.4f departs from Issue #4's %.4f",
                  ox_wt, OX_WT_REF))

# --- gene level, so a set-level flat is not hiding a member that moves ---------
wt_genes <- purrr::map_dfr(apop_sets, function(s) {
  e <- ens_set(gmt[[s]])
  i <- match(e, rownames(D$timepoint_neg_raw))
  tibble::tibble(
    gene = cdf$mgi_symbol[match(e, cdf$gene)], arm = sub("MITOCARTA_APOPTOSIS_", "", s),
    baseMean = as.numeric(bm[e]), wt_time = as.numeric(tn[e]),
    padj_wt = D$timepoint_neg_raw$padj[i],
    myc_6W = as.numeric(m6[e]), myc_12W = as.numeric(m12[e]))
}) |> dplyr::arrange(wt_time)

# --- CORE_MITO decomposed: is the fall regulon-wide, or only its OXPHOS part? --
core_e <- ens_set(gmt[["CORE_MITO"]])
oxs_e  <- ens_set(gmt[["MITOCARTA_OXPHOS_SUBUNITS"]])
core_parts <- list(`CORE_MITO n OXPHOS subunits` = intersect(core_e, oxs_e),
                   `CORE_MITO, rest`             = setdiff(core_e, oxs_e),
                   `CORE_MITO, all`              = core_e)
core_decomp <- purrr::imap_dfr(core_parts, function(e, nm) {
  n <- null_pct(tn, e)
  tibble::tibble(part = nm, n_genes = length(e),
                 c_wt_time = n[["observed"]], wt_null_pct = n[["percentile"]],
                 c_myc_6W = set_mean(m6, e))
})

# =============================================================================
# PART B2: THE LINEAGE-SUPPRESSION CARRIER
# -----------------------------------------------------------------------------
# Gray 2023 p.7: "the AP4, HS4, and BA2 clusters constituted a 'low-expressing'
# (LE) state within their respective lineages compared to their high-expressing
# (HE) counterparts". LE = LINEAGE-PROGRAM SUPPRESSION. So MG_HEVSLE_*_UP are
# HE-markers and MG_HEVSLE_*_DN are LE-markers.
#
# PREDICTION, stated before the test: if the maturing gland loses LE abundance,
# _DN falls and _UP rises on the wild-type axis. AP is the prediction (ESRRA and
# GABPA's mito programmes are detected in AP_LE); HS and BA are the SPECIFICITY
# CONTROLS -- a shift in all three lineages is global lineage suppression, not an
# AP-carrier effect.
#
# NUANCE from the paper that bounds the claim: AP4/HS4/BA2 are "broadly
# distributed among samples of all ages", so what a shift here would show is a
# change in LE ABUNDANCE, not a state appearing or vanishing. The pubertally
# restricted LE clusters are AP5 and HS5.
# =============================================================================
le_sets <- tibble::tribble(
  ~lineage, ~role,        ~set,
  "AP",     "HE marker",  "MG_HEVSLE_AP_GRAY_UP",
  "AP",     "LE marker",  "MG_HEVSLE_AP_GRAY_DN",
  "HS",     "HE marker",  "MG_HEVSLE_HS_GRAY_UP",
  "HS",     "LE marker",  "MG_HEVSLE_HS_GRAY_DN",
  "BA",     "HE marker",  "MG_HEVSLE_BA_GRAY_UP",
  "BA",     "LE marker",  "MG_HEVSLE_BA_GRAY_DN",
  "AP",     "lineage",    "MG_AP_GRAY",
  "HS",     "lineage",    "MG_HS_GRAY",
  "BA",     "lineage",    "MG_BASAL_GRAY",
  "HS",     "TEB vs duct","MG_TEB_VS_DUCTAL_HS_GRAY_UP",
  "AP",     "TEB vs duct","MG_TEB_VS_DUCTAL_AP_GRAY_UP")

# CIRCULARITY GUARD, and it is not hypothetical: 25 of the 89 OXPHOS subunits ARE
# AP LE-marker genes (28% of the set). So "OXPHOS falls" and "AP LE markers fall"
# are partly the same measurement, and a lineage claim built on the raw sets would
# be the OXPHOS claim restated. Every row therefore carries a MITO-REMOVED value;
# the lineage reading must survive on `c_wt_time_nonmito`, not on `c_wt_time`.
mito_ens <- ens_set(mito_universe)
le_content <- purrr::pmap_dfr(le_sets, function(lineage, role, set) {
  e  <- ens_set(gmt[[set]])
  nm <- setdiff(e, mito_ens)
  n  <- null_pct(tn, e); n2 <- null_pct(tn, nm)
  tibble::tibble(lineage = lineage, role = role, set = set, n_genes = length(e),
                 c_wt_time = n[["observed"]], null_median = n[["null_median"]],
                 wt_null_pct = n[["percentile"]],
                 n_mito_in_set = length(e) - length(nm),
                 n_genes_nonmito = length(nm),
                 c_wt_time_nonmito = n2[["observed"]],
                 wt_null_pct_nonmito = n2[["percentile"]],
                 c_mycpos_time = set_mean(tpz, e), c_myc_6W = set_mean(m6, e))
})

# Is the PGC1a-axis programme actually carried by the LE / AP state? Overlap, so
# "lineage-restricted" is a number and not a reading of the label.
carrier_overlap <- purrr::map_dfr(
  c("CORE_MITO", "TFT_ESRRA_GRAY_AP_LE", "TFT_GABPA_GRAY_AP_LE",
    "TFT_NRF1_GRAY_AP_TEB", "MITOCARTA_OXPHOS_SUBUNITS"),
  function(pg) purrr::map_dfr(le_sets$set, function(ls) {
    a <- gmt[[pg]]; b <- gmt[[ls]]
    tibble::tibble(programme = pg, state_set = ls, n_programme = length(a),
                   n_state = length(b), overlap = length(intersect(a, b)),
                   frac_of_programme = length(intersect(a, b)) / length(a))
  }))

# Descriptive only (n = 12 wild-types, batch = timepoint): does the OXPHOS score
# track the LE score across the wild-type mice? This CANNOT separate a lineage
# shift from developmental regulation -- deconvolution is the design that could.
# The adjustment MUST use a DISJOINT score, or it is self-adjustment: the raw AP-LE
# score shares 25 genes with the OXPHOS subunits, so regressing one on the other
# would remove the outcome from itself and the beta would collapse for arithmetic
# reasons. `ap_le_nm` is the AP LE-marker set with every MitoCarta gene removed.
wt_i <- which(sm$myc_status == "neg")
score_from <- function(ens, min_n = 5L) {
  e <- ens[ens %in% rownames(L)]
  if (length(e) < min_n) return(NULL)
  colMeans(zrow(L[e, , drop = FALSE]))
}
le_within <- local({
  ox_e    <- ens_set(gmt[["MITOCARTA_OXPHOS_SUBUNITS"]])
  ap_le_e <- ens_set(gmt[["MG_HEVSLE_AP_GRAY_DN"]])
  ox      <- score_from(ox_e)
  ap_le   <- score_from(ap_le_e)
  ap_le_nm <- score_from(setdiff(ap_le_e, mito_ens))
  ap_he_nm <- score_from(setdiff(ens_set(gmt[["MG_HEVSLE_AP_GRAY_UP"]]), mito_ens))
  if (is.null(ox) || is.null(ap_le_nm)) return(NULL)
  b <- function(f) unname(stats::coef(f)[2])
  tibble::tibble(
    n_wt                  = length(wt_i),
    n_shared_ox_apLE      = length(intersect(ox_e, ap_le_e)),
    cor_ox_apLE_raw       = stats::cor(ox[wt_i], ap_le[wt_i]),
    cor_ox_apLE_nonmito   = stats::cor(ox[wt_i], ap_le_nm[wt_i]),
    cor_ox_apHE_nonmito   = stats::cor(ox[wt_i], ap_he_nm[wt_i]),
    ox_beta_time_raw      = b(stats::lm(ox[wt_i] ~ sm$tp[wt_i])),
    ox_beta_time_adj_nm   = b(stats::lm(ox[wt_i] ~ sm$tp[wt_i] + ap_le_nm[wt_i])),
    ox_beta_time_adj_self = b(stats::lm(ox[wt_i] ~ sm$tp[wt_i] + ap_le[wt_i])))
})

# =============================================================================
# PART C: THE COLLAPSE SCAN
# -----------------------------------------------------------------------------
# Statistic. The interaction is LFC12 - LFC6, and under a global rescaling every
# Myc-responsive gene has one, so "has a negative interaction" is not a finding.
# What is a finding is departing from the GLOBAL RATE:
#
#     resid_aligned = sign(LFC6) * (LFC12 - RATE * LFC6)
#
# Sign-aligning matters: for a Myc-INDUCED gene a collapse makes the residual
# negative, for a REPRESSED gene it makes it positive, so without the alignment
# the two halves of the transcriptome cancel. Negative = collapsed more than the
# global rate; positive = retained more.
#
# LFC6 and LFC12 are estimated from DISJOINT sample sets (6W and 12W), so their
# sampling errors are independent and the residual's standard error is
# sqrt(SE12^2 + RATE^2 * SE6^2). Ranking on the standardised residual is therefore
# the statistically right choice; the unstandardised one and the retention ratio
# are carried alongside as robustness.
#
# The scan is restricted to genes with a REAL 6W effect (script 42's correction
# #1: you cannot lose an effect you never had), which also changes the question
# fgsea answers -- it becomes "among Myc-responsive genes, which programmes
# collapse", which is the question we want.
# =============================================================================
se6  <- stats::setNames(D$myc_6W_raw$lfcSE,  rownames(D$myc_6W_raw))
se12 <- stats::setNames(D$myc_12W_raw$lfcSE, rownames(D$myc_12W_raw))
padj6 <- stats::setNames(D$myc_6W_raw$padj,  rownames(D$myc_6W_raw))

# TWO gene sets, deliberately. The RANKING set carries the padj filter, so the
# fgsea question stays "among genes with a demonstrable 6W Myc effect, which
# programmes collapse". The REPORTING set drops it, because the pre-specified
# genes fail it -- Bbc3 padj 0.22 and Bcl2l11 padj 0.46 at 6W -- and excluding
# them would remove exactly what the scan was built to locate. PUMA's finding was
# never in its level (that is script 42's PUMA:Bcl-xL ratio); reporting it here
# without letting it into the ranking is the honest way to have both.
keep_rep <- universe_all[abs(m6) >= LFC6_FLOOR & bm >= BM_FLOOR & is.finite(m12)]
keep6    <- keep_rep[!is.na(padj6[keep_rep]) & padj6[keep_rep] < PADJ6]

# ASSERTION: the interaction really is LFC12 - LFC6 with the sign we assume.
int_raw   <- V("interaction_raw")
int_check <- stats::cor(int_raw[keep6], (m12 - m6)[keep6])
if (!is.finite(int_check) || int_check < 0.99)
  stop(sprintf("44 PART C: interaction_raw is not (myc_12W - myc_6W); cor = %.3f", int_check))

# Re-fit the global rate on this gene set rather than assuming script 40's 0.55.
rate_fit  <- stats::lm(m12[keep6] ~ 0 + m6[keep6])
RATE_HAT  <- unname(stats::coef(rate_fit)[1])

stat_for <- function(g) {
  ra <- sign(m6[g]) * (m12[g] - RATE_HAT * m6[g])
  sr <- sqrt(se12[g]^2 + RATE_HAT^2 * se6[g]^2)
  list(resid_aligned = as.numeric(ra), se_resid = as.numeric(sr),
       z_resid = as.numeric(ra / sr), retention = as.numeric(m12[g] / m6[g]))
}
S <- stat_for(keep_rep)

collapse_genes <- tibble::tibble(
  ens = keep_rep, gene = cdf$mgi_symbol[match(keep_rep, cdf$gene)],
  baseMean = as.numeric(bm[keep_rep]),
  lfc_6W = as.numeric(m6[keep_rep]), lfc_12W = as.numeric(m12[keep_rep]),
  padj_6W = as.numeric(padj6[keep_rep]),
  resid_aligned = S$resid_aligned, se_resid = S$se_resid, z_resid = S$z_resid,
  retention = S$retention, ret_minus_rate = S$retention - RATE_HAT,
  myc_induced = as.numeric(m6[keep_rep]) > 0,
  in_ranking = keep_rep %in% keep6) |>
  dplyr::mutate(pct_z = 100 * rank(z_resid) / dplyr::n()) |>
  dplyr::arrange(z_resid)

# ROBUSTNESS GATE: the two statistics must agree, or neither is usable.
rank_agreement <- stats::cor(collapse_genes$z_resid, collapse_genes$ret_minus_rate,
                             method = "spearman", use = "pairwise.complete.obs")

# --- the pre-specified genes --------------------------------------------------
# `in_ranking = FALSE` marks a gene reported but kept out of the fgsea ranking.
# `pct_z` still locates it in the full reporting distribution, which is the number
# to read: it says where this gene sits among all Myc-responsive genes.
PRESPEC <- c("Bbc3", "Bcl2l11")
BH3     <- c("Bbc3", "Bcl2l11", "Bid", "Bik", "Bad", "Bmf", "Pmaip1", "Bnip3",
             "Bnip3l", "Hrk", "Bok", "Bax", "Bak1")
prespec_position <- collapse_genes |>
  dplyr::filter(gene %in% union(PRESPEC, BH3)) |>
  dplyr::mutate(pre_specified = gene %in% PRESPEC) |>
  dplyr::arrange(z_resid)

# --- fgsea on the residual ranking -------------------------------------------
# Pathways are reconciled to Ensembl one set at a time so that renamed symbols are
# not silently dropped; the ranking is over Ensembl IDs for the same reason.
all_paths <- c(gmt, ap6)
all_paths <- all_paths[!duplicated(names(all_paths))]
message("44 PART C: reconciling ", length(all_paths), " pathways to Ensembl ...")
paths_ens <- vector("list", length(all_paths)); names(paths_ens) <- names(all_paths)
for (i in seq_along(all_paths)) {
  paths_ens[[i]] <- ens_set(all_paths[[i]])
  if (i %% 200 == 0) message("   ... ", i, " / ", length(all_paths))
}
paths_ens <- paths_ens[lengths(paths_ens) >= 8]
# The library has no ATF4/ISR set, and grepping for one matched only Hallmark UPR,
# which is a different question. Add the named target list explicitly so question
# (d) is actually asked rather than approximated.
paths_ens[["ISR_ATF4_TARGETS_ADHOC"]] <- ens_set(c(
  "Atf4", "Ddit3", "Trib3", "Chac1", "Sesn2", "Asns", "Atf3", "Nupr1", "Eif2ak3",
  "Ppp1r15a", "Cebpb", "Slc7a11", "Aldh18a1", "Mthfd2", "Psat1", "Phgdh", "Shmt2",
  "Cth", "Wars1", "Gars1", "Sars1", "Cebpg", "Xbp1", "Herpud1"))

rank_in <- collapse_genes |> dplyr::filter(in_ranking)
rank_z  <- sort(stats::setNames(rank_in$z_resid, rank_in$ens), decreasing = TRUE)
set.seed(1)
collapse_fgsea <- suppressWarnings(fgsea::fgsea(
  pathways = paths_ens, stats = rank_z, minSize = 8, maxSize = 800, eps = 0)) |>
  tibble::as_tibble() |>
  dplyr::select(pathway, NES, pval, padj, size) |>
  dplyr::arrange(NES)

# The pre-registered questions, reported whether they pass or fail. Nothing here
# is a post-hoc pick: the list is fixed before the ranking is looked at.
PREREG <- list(
  `a. OXPHOS / MitoCarta`     = c("MITOCARTA_OXPHOS_SUBUNITS", "MITOCARTA_OXPHOS",
                                  "MITOCARTA_OXPHOS_ASSEMBLY_FACTORS"),
  `b. apoptosis`              = c("MITOCARTA_APOPTOSIS_PRO", "MITOCARTA_APOPTOSIS_ANTI",
                                  "TANG_APOPTOSIS", "HALLMARK_APOPTOSIS",
                                  "CDC_PRODEATH_APOPTOSIS"),
  `c. AP_TEB accessible`      = grep("_AP_TEB$", names(paths_ens), value = TRUE),
  `d. ISR / ATF4 targets`     = c("ISR_ATF4_TARGETS_ADHOC",
                                  grep("UNFOLDED_PROTEIN", names(paths_ens),
                                       value = TRUE, ignore.case = TRUE)),
  `e. TEB vs ductal`          = grep("^MG_TEB_VS_DUCTAL", names(paths_ens), value = TRUE),
  `f. NRF1 / PGC1a axis`      = c(grep("^TFT_NRF1", names(paths_ens), value = TRUE),
                                  "CORE_MITO", "TFT_ESRRA_GRAY_AP_LE",
                                  "TFT_GABPA_GRAY_AP_LE"))
prereg_results <- purrr::imap_dfr(PREREG, function(sets, q)
  collapse_fgsea |>
    dplyr::filter(pathway %in% sets) |>
    dplyr::mutate(question = q, .before = 1))

# =============================================================================
# PART D: MECHANISM CANDIDATES
# -----------------------------------------------------------------------------
# Routes by which PUMA priming could depend on the respiratory state. Each group
# is scored on the same residual and tested against a DOUBLY matched null (script
# 42's idiom: matched on baseMean AND |LFC6| decile, which is what makes it a test
# of EXCESS collapse rather than of effect size). Ranking only; nothing here is
# offered as the mechanism.
#
# The ISR row is the test script 42 did NOT do: it asked whether ATF4 and friends
# TRACK PUMA per sample, not whether they COLLAPSE like PUMA -- and ATF4 is
# translationally regulated, so its mRNA is a poor readout of the pathway either
# way.
# =============================================================================
MECH <- list(
  `ISR / ATF4 targets`   = c("Atf4", "Ddit3", "Trib3", "Chac1", "Sesn2", "Asns",
                             "Atf3", "Nupr1", "Eif2ak3"),
  `cardiolipin`          = c("Crls1", "Taz", "Ptpmt1", "Pla2g6", "Pgs1"),
  `cristae / OPA1`       = c("Opa1", "Immt", "Chchd3", "Chchd6", "Apoo", "Apool"),
  `IMS payload`          = c("Cycs", "Aifm1", "Aifm2", "Endog", "Diablo", "Htra2"),
  `BH3-only, prespecified` = PRESPEC,
  `BH3-only, the rest`   = setdiff(BH3, PRESPEC),
  `FOXO`                 = c("Foxo1", "Foxo3", "Foxo4"))

# Deciles, not ventiles, and over the REPORTING set: the double stratification on
# the ranking set alone left 5-13 genes per cell, and a percentile out of seven
# genes is not a measurement. Deciles over ~keep_rep give order-100 per cell.
q10 <- function(x) cut(x, unique(stats::quantile(x, 0:10 / 10)), include.lowest = TRUE)
db  <- q10(bm[keep_rep]); dl <- q10(abs(m6[keep_rep]))
names(db) <- keep_rep; names(dl) <- keep_rep
res_named <- stats::setNames(collapse_genes$resid_aligned, collapse_genes$ens)

mech_genes <- purrr::imap_dfr(MECH, function(syms, grp) {
  e <- ens_of(syms); ok <- !is.na(e) & e %in% keep_rep
  if (!any(ok)) return(NULL)
  purrr::map_dfr(which(ok), function(j) {
    g <- e[j]
    sel <- keep_rep[db[keep_rep] == db[g] & dl[keep_rep] == dl[g]]
    tibble::tibble(group = grp, gene = syms[j],
                   lfc_6W = as.numeric(m6[g]), lfc_12W = as.numeric(m12[g]),
                   resid_aligned = as.numeric(res_named[g]),
                   retention = as.numeric(m12[g] / m6[g]),
                   matched_n = length(sel),
                   matched_median = stats::median(res_named[sel], na.rm = TRUE),
                   percentile = 100 * mean(res_named[sel] < res_named[g], na.rm = TRUE))
  })
})

# min/max percentile alongside the median, because a small group whose members sit
# at OPPOSITE ends has a median that describes neither of them -- which is exactly
# what the pre-specified pair does (PUMA collapses, BIM does not).
mech_groups <- mech_genes |>
  dplyr::group_by(group) |>
  dplyr::summarise(n_genes = dplyr::n(),
                   median_resid = stats::median(resid_aligned, na.rm = TRUE),
                   median_retention = stats::median(retention, na.rm = TRUE),
                   median_percentile = stats::median(percentile, na.rm = TRUE),
                   min_percentile = min(percentile, na.rm = TRUE),
                   max_percentile = max(percentile, na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::arrange(median_percentile)

# =============================================================================
# NOTES -- written before the numbers are read, so they bound the reading
# =============================================================================
pf <- function(p) if (is.na(p)) "NA" else if (p < 1e-4) "<1e-4" else sprintf("%.3f", p)
own_core <- ownership |> dplyr::filter(programme == "CORE_MITO",
                                       apoptosis_set == "MITOCARTA_APOPTOSIS_PRO")
own_myc  <- ownership |> dplyr::filter(programme == "MYC_MITO",
                                       apoptosis_set == "MITOCARTA_APOPTOSIS_PRO")

notes <- c(
  "44 -- CAN WE DEFINE AN OXPHOS/APOPTOSIS MODULE?",
  "",
  sprintf(paste("PART A OWNERSHIP. The PGC1a axis (CORE_MITO, n=%d) overlaps the 25",
                "pro-apoptotic MitoCarta genes in %d (expected %.1f, fold %.2f, p_enrich %s):",
                "%s. MYC_MITO (n=%d) overlaps in %d (expected %.1f, fold %.2f, p_enrich %s).",
                "Neither contains Bbc3. This is SET MEMBERSHIP in a curated resource -- it",
                "constrains what a TITLE may claim, it does not establish regulation here."),
          own_core$n_programme, own_core$overlap, own_core$expected, own_core$fold,
          pf(own_core$p_enrich), own_core$genes,
          own_myc$n_programme, own_myc$overlap, own_myc$expected, own_myc$fold,
          pf(own_myc$p_enrich)),
  "",
  sprintf(paste("  Base rate: across %d shortlist rows / %d TFs the pro-apoptotic capture is",
                "%.4f per mito gene against a catalog expectation of %.4f, and",
                "cor(regulon size, n_pro) = %.2f. So a biogenesis TF 'touching' a",
                "pro-apoptotic effector is arithmetic, not signal, and the SHORTLIST cannot",
                "be cited for the NRF1->Bbc3 axis. The axis itself rests on the author's",
                "PGC1a westerns (PUMA and BIM), which are independent of this dataset."),
          shortlist_baserate$n_rows, shortlist_baserate$n_TFs,
          shortlist_baserate$pro_per_mito_gene, shortlist_baserate$catalog_expectation,
          shortlist_baserate$cor_size_vs_pro),
  "",
  paste("PART B. Read c_wt_time and wt_null_pct together with p_wt_time. If OXPHOS",
        "subunits sit at the bottom on BOTH rulers while the apoptosis arms sit mid-pack,",
        "there is no common OXPHOS/apoptosis de-prioritisation in the wild-type gland, and",
        "a module cannot be defined by co-expression. core_decomp then says whether the",
        "PGC1a axis withdraws as a unit or only through its respiratory part."),
  "",
  sprintf(paste("PART B2. PREDICTION STATED IN ADVANCE: LE-marker sets (_DN) fall and",
                "HE-marker sets (_UP) rise if the gland loses LE abundance. AP is the",
                "prediction; HS and BA are specificity controls -- movement in all three is",
                "global lineage suppression, not an AP carrier. CIRCULARITY GUARD: %d of the",
                "%d OXPHOS subunits are AP LE-marker genes, so read c_wt_time_nonmito, not",
                "c_wt_time -- the raw column is partly the OXPHOS result restated. Likewise",
                "le_within reports ox_beta_time_adj_nm (disjoint) next to",
                "ox_beta_time_adj_self (shared genes), and only the first is a test.",
                "le_within is DESCRIPTIVE ONLY (n=12, batch = timepoint) and cannot separate",
                "a lineage shift from developmental regulation."),
          length(intersect(ens_set(gmt[["MITOCARTA_OXPHOS_SUBUNITS"]]),
                           ens_set(gmt[["MG_HEVSLE_AP_GRAY_DN"]]))),
          length(ens_set(gmt[["MITOCARTA_OXPHOS_SUBUNITS"]]))),
  "",
  sprintf(paste("PART C. Two gene sets: %d genes are REPORTED (baseMean >= %d, |LFC6| >=",
                "%.1f) and %d of those enter the fgsea RANKING (6W padj < %.2f). The",
                "pre-specified genes fail the padj filter -- Bbc3 and Bcl2l11 -- so they are",
                "reported with in_ranking = FALSE rather than dropped; PUMA's finding was",
                "never in its level. Global rate re-fitted on the ranking set: %.3f (script",
                "40's value 0.55). Rank agreement between the standardised residual and the",
                "retention ratio: Spearman %.3f -- BELOW ~0.7 THIS IS A STOP, not a footnote,",
                "because a ratio of noisy quantities cannot stand alone."),
          length(keep_rep), BM_FLOOR, LFC6_FLOOR, length(keep6), PADJ6,
          RATE_HAT, rank_agreement),
  "",
  paste("  A collapse module here is a HYPOTHESIS ABOUT CO-DEPENDENCE, not a demonstration",
        "of co-regulation: n=24, batch = timepoint, and the ranking is over Myc-responsive",
        "genes only. It needs the PGC1a perturbation transcriptome to become one. Bbc3 and",
        "Bcl2l11 are PRE-SPECIFIED from the westerns, not scan hits; Bcl2l11's 6W effect is",
        "non-significant (p 0.46), so its retention is reported and never headlined."),
  "",
  local({
    p <- prespec_position
    g <- function(x) p[match(x, p$gene), ]
    sprintf(paste("  Where the two PRE-SPECIFIED genes landed, both named in advance from the",
                  "PGC1a westerns: Bbc3 at the %.2f percentile (retention %+.2f) and Bcl2l11",
                  "at the %.1f percentile (retention %+.2f). They are at OPPOSITE ends. The",
                  "in-vivo collapse is PUMA-specific; BIM does not collapse, so the western's",
                  "PUMA+BIM pair does NOT transfer as a pair. Report both -- half a",
                  "pre-specification failing is the part that makes the other half worth",
                  "anything."),
            g("Bbc3")$pct_z, g("Bbc3")$retention,
            g("Bcl2l11")$pct_z, g("Bcl2l11")$retention)
  }),
  "",
  paste("PART D. Ranking only. The ISR row is the test script 42 did not do -- it asked",
        "whether ATF4 and friends TRACK PUMA per sample, not whether they COLLAPSE like it.",
        "Read min_percentile and max_percentile, not the median: these groups are small and",
        "a split group's median describes neither member."),
  "",
  paste("LIBRARY METADATA CORRECTION recorded here: Gray 2023 p.7 defines LE as",
        "'LOW-EXPRESSING' (lineage-program suppression), NOT low oestrogen. Our",
        "provenance_table.csv glosses MG_HEVSLE_* as 'HE vs LE (AP) | adult' and the",
        "catalog expands it to 'High- vs low-estrogen'. That gloss is wrong. AP =",
        "alveolar progenitor = LASP in the consensus nomenclature."))

out <- list(
  ownership          = ownership,
  shortlist_baserate = shortlist_baserate,
  shortlist_roster   = shortlist_roster,
  wt_content         = wt_content,
  wt_genes           = wt_genes,
  core_decomp        = core_decomp,
  le_content         = le_content,
  carrier_overlap    = carrier_overlap,
  le_within          = le_within,
  collapse_genes     = collapse_genes,
  prespec_position   = prespec_position,
  collapse_fgsea     = collapse_fgsea,
  prereg_results     = prereg_results,
  mech_genes         = mech_genes,
  mech_groups        = mech_groups,
  defs = list(global_rate_assumed = GLOBAL_RATE, global_rate_fitted = RATE_HAT,
              rank_agreement_spearman = rank_agreement,
              n_reported_genes = length(keep_rep), n_ranking_genes = length(keep6),
              padj6 = PADJ6,
              lfc6_floor = LFC6_FLOOR, basemean_floor = BM_FLOOR,
              n_set_draws = NSET, n_bins = NBIN,
              oxphos_wt_reference = OX_WT_REF, oxphos_wt_observed = ox_wt,
              interaction_sign_check = int_check,
              pre_specified_genes = PRESPEC,
              promotion_gate = c(bh_fdr = GATE_FDR, n_mito = GATE_MITO)),
  analysis_date = Sys.Date(),
  notes = notes)

saveRDS(out, here::here("results", "collapse_module_ownership.rds"))
message("44: wrote results/collapse_module_ownership.rds")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  res <- readRDS(here::here("results", "collapse_module_ownership.rds"))

  ## PART A -- the title question. Read `fold` and `genes`. If CORE_MITO is at or
  ## below 1.0 while MYC_MITO is above it, the pro-apoptotic genes belong to Myc's
  ## regulon and not to the PGC1a axis, and "a PGC-1a/Nrf1 regulated OXPHOS/
  ## apoptosis module" cannot be the title.
  res$ownership |> dplyr::filter(apoptosis_set == "MITOCARTA_APOPTOSIS_PRO") |>
    dplyr::select(programme, n_programme, overlap, expected, fold, p_enrich, genes) |>
    print()

  ## The base rate that disqualifies the shortlist as a citation. pro_per_mito_gene
  ## should sit on top of catalog_expectation.
  res$shortlist_baserate |> print()

  ## The roster rows against the resource's OWN gate. NRF1/AP_TEB should fail it.
  res$shortlist_roster |> print(n = 30)

  ## PART B -- both rulers side by side. c_wt_time + wt_null_pct is the content
  ## ruler; p_wt_time is mitoPPS priority. The claim under test is whether
  ## apoptosis moves WITH OXPHOS. Opposite signs on the priority ruler would mean
  ## it does not.
  res$wt_content |> print(n = 12)

  ## Nothing hiding at gene level? Only significant padj_wt entries matter.
  res$wt_genes |> dplyr::filter(!is.na(padj_wt) & padj_wt < 0.05) |> print(n = 20)

  ## Does the PGC1a axis withdraw as a unit, or only through its OXPHOS part?
  res$core_decomp |> print()

  ## PART B2 -- the prediction. READ THE NONMITO COLUMNS: n_mito_in_set says how
  ## much of each set is shared with the mitochondrial result, and c_wt_time is
  ## partly that result restated. For AP the "LE marker" row should be NEGATIVE and
  ## the "HE marker" row POSITIVE. Then check HS and BA: if they move the same way,
  ## it is global lineage suppression rather than an AP-specific carrier.
  res$le_content |> dplyr::arrange(lineage, role) |>
    dplyr::select(lineage, role, n_genes, c_wt_time, wt_null_pct,
                  n_mito_in_set, c_wt_time_nonmito, wt_null_pct_nonmito) |> print(n = 12)

  ## Is the programme actually carried by that state? frac_of_programme is the
  ## number that replaces the label reading.
  res$carrier_overlap |> dplyr::filter(programme == "CORE_MITO") |> print(n = 12)

  ## Descriptive only, and compare the two adjustments: ox_beta_time_adj_self shares
  ## n_shared_ox_apLE genes with the outcome and will collapse for arithmetic
  ## reasons. Only ox_beta_time_adj_nm is a test -- and at n=12 with batch =
  ## timepoint it ranks the idea, it does not settle it.
  res$le_within |> print()

  ## PART C -- THE GATE FIRST. Below ~0.7 stop and do not read anything else here.
  res$defs$rank_agreement_spearman
  res$defs$global_rate_fitted        # next to 0.55

  ## Where do the pre-specified genes sit? pct_z is the percentile in the REPORTED
  ## distribution: low = collapsed far more than the global rate. in_ranking = FALSE
  ## on Bbc3/Bcl2l11 is expected -- they fail the 6W padj filter and are reported,
  ## not ranked. Read them against Bax, which retains at the global rate.
  res$prespec_position |> print(n = 15)

  ## The collapse tail itself. Is it a set, or is PUMA a solo?
  res$collapse_genes |> head(30) |> print()

  ## The mirror -- what KEEPS its Myc-inducibility. Without this the tail is not
  ## interpretable.
  res$collapse_genes |> dplyr::arrange(dplyr::desc(z_resid)) |> head(30) |> print()

  ## THE MODULE QUESTION. Negative NES = collapses faster than the global rate.
  ## Anything with padj < 0.05 at the negative end is a candidate module.
  res$collapse_fgsea |> head(25) |> print()

  ## The pre-registered questions, pass or fail. Report these whatever they say --
  ## they were fixed before the ranking was looked at.
  res$prereg_results |> print(n = 40)

  ## PART D -- mechanism candidates, ranking only. median_percentile below ~25 means
  ## the group collapses more than expression-matched genes.
  res$mech_groups |> print()
  res$mech_genes |> dplyr::arrange(percentile) |> print(n = 40)

  ## Sanity: the positive control, within 0.02 of -0.2548.
  res$defs$oxphos_wt_observed

  cat(res$notes, sep = "\n")
}
