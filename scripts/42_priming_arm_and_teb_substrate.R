# scripts/42_priming_arm_and_teb_substrate.R
# =============================================================================
# Block B -- DOSE vs COMPETENCE: which apoptotic arm does Myc lose by 12W, and
#            what does the gland lose underneath it?
# =============================================================================
#
# WHY THIS SCRIPT EXISTS. Two blots reframed the whole attenuation question.
#
#   (1) MYC PROTEIN falls ~50% 6W->12W at CONSTANT transcript (figS8). That single
#       fact explains the entire transcriptome attenuation: script 40 measured a
#       uniform x0.55 rescaling with no reshaping, script 34 found no programme in
#       excess, and the 2026-07-25 scans found no repressor and no moderating axis.
#       Output per unit MYC is unchanged. THE MMTV-Myc TIMELINE IS A DOSE EXPERIMENT.
#
#   (2) The MYC-ER mouse is NOT. Acute induction at 12W gives ~80% less death and
#       only ~30% less proliferation at EQUAL or HIGHER MYC protein. Those animals
#       had never seen MYC before induction, so the change cannot be MYC-driven
#       selection -- it is a property of the maturing gland. THAT IS A COMPETENCE
#       EXPERIMENT, and it is the reason a substrate story is needed at all.
#
# This script asks what our transcriptome can say about the competence half.
#
# THE ANSWER, in one line: Myc builds a death-ready mitochondrion at BOTH ages --
# Htra2, Bax up, Bcl-xL down, biogenesis up -- and what it loses by 12W is the
# TRIGGER. Bax:Bcl-xL priming fades at exactly the global rate; PUMA:Bcl-xL priming
# collapses to zero. PUMA is precisely the BH3-only that the PGC1a cell experiments
# single out, so the tested pair was pre-specified by external perturbation data.
#
# THREE NEGATIVES TRAVEL WITH IT AND ARE REPORTED AS PROMINENTLY (PART C, PART D):
#   - the cross-sectional OXPHOS<->PUMA coupling FAILS (at/below the ambient), so the
#     "OXPHOS-PUMA bundle" is carried by the cell perturbations, NEVER by in-vivo
#     correlation;
#   - Bbc3 is NOT in the TEB-vs-ductal signature -- the TEB claim is about TF-binding
#     CONTEXT (Gray/ChEA AP_TEB regulons), not about PUMA being a TEB transcript;
#   - the MEC state composites (LASP/LHS/BMYO/LP) are not significant; only markers move.
#
# CEILING. n=6/group, exploratory. The 6W-vs-12W axis is BATCH-CONFOUNDED (batch =
# timepoint), so every WT temporal statement in PART D is DESCRIBED, not claimed --
# with the one mitigation that the protein blots replicate its direction off the RNA
# batch entirely. The batch-clean quantities are the genotype effects at each age and
# their interaction. Script 34's circularity caveat still binds COMPOSITE death
# scores (MITOCARTA_APOPTOSIS_PRO contains Bbc3, Bax, Bcl2l11 by construction); this
# script therefore makes NO composite-correlation claim -- it tests named genes with
# genome-wide FDR and per-mouse ratios against expression-matched nulls.
#
# Input:  results/interaction_results.rds   (raw/unshrunken LFC + p + padj)
#         results/dds_int_run.rds           (normalised counts, colData)
#         results/combined_df_annotated.rds (mgi_symbol <-> ensembl)
#         data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt
#         docs/library_reference/gray_chea_mito_tf_shortlist.csv  (read, not re-derived)
# Output: results/priming_arm_teb.rds
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

set.seed(1)
NPAIR       <- 4000L    # random pro:anti pairs for the matched-pair null
NAMB        <- 3000L    # arbitrary genes for the ambient-coupling null
NPERM       <- 5000L    # within-timepoint shuffles for PART H's axis null
GLOBAL_RATE <- 0.55     # script 40's rescaling slope; PART B reproduces it independently

# =============================================================================
# PART 0: LOAD
# =============================================================================
ir  <- readRDS(here::here("results", "interaction_results.rds"))
dds <- readRDS(here::here("results", "dds_int_run.rds"))
cdf <- readRDS(here::here("results", "combined_df_annotated.rds"))
gmt <- fgsea::gmtPathways(here::here("data", "genesets_from_library",
                                     "mammary_mito_myc_metab_v1_mouse.gmt"))

# DESeqResults are S4; coerce once so rownames() and [ , ] behave as data.frames.
CONTRASTS <- c("myc_6W_raw", "myc_12W_raw", "timepoint_neg_raw", "timepoint_pos_raw",
               "interaction_raw")
D  <- lapply(ir[CONTRASTS], as.data.frame)
nc <- DESeq2::counts(dds, normalized = TRUE)
L  <- log2(nc + 1)
sm <- as.data.frame(SummarizedExperiment::colData(dds))
sm$group <- factor(sm$group, levels = c("6W_neg", "6W_pos", "12W_neg", "12W_pos"))
sm$tp    <- stats::relevel(as.factor(sm$timepoint),  "6W")
sm$myc   <- stats::relevel(as.factor(sm$myc_status), "neg")

sym2ens <- cdf[!is.na(cdf$mgi_symbol) & !duplicated(cdf$mgi_symbol), c("mgi_symbol", "gene")]
ens_of  <- function(s) sym2ens$gene[match(s, sym2ens$mgi_symbol)]
zrow    <- function(m) t(scale(t(m)))
comp_of <- function(syms) {                       # per-sample z-mean over a symbol list
  e <- ens_of(syms); e <- e[!is.na(e) & e %in% rownames(L)]
  stopifnot(length(e) >= 3)
  colMeans(zrow(L[e, , drop = FALSE]))
}
set_score <- function(set, min_n = 5L) {          # per-sample z-mean over a GMT set
  if (is.null(gmt[[set]])) return(NULL)
  e <- recon_to_ensembl(gmt[[set]], rownames(L)); e <- e[!is.na(e)]
  if (length(e) < min_n) return(NULL)
  colMeans(zrow(L[e, , drop = FALSE]))
}

# The two purity covariates. There is no external QC metadata (CLAUDE.md), so these
# are RNA surrogates for composition -- adequate for asking "is the effect merely
# contamination?", not for measuring purity.
epi_comp <- comp_of(c("Epcam", "Krt8", "Krt18", "Krt5", "Krt14", "Cdh1", "Krt17"))
imm_comp <- comp_of(c("Ptprc", "Cd52", "Cd3e", "Lyz2", "Cd74", "H2-Aa", "Itgam", "Ms4a1"))

# =============================================================================
# PART A: THE DEATH-READY MITOCHONDRION
# -----------------------------------------------------------------------------
# What Myc does to the execution machinery and to biogenesis at EACH age, on the
# clean genotype axis. `retention` = LFC12 / LFC6, to be read against GLOBAL_RATE:
# at 0.55 a gene is simply following the dose drop; far below it is losing something
# extra. Retention is unstable when LFC6 is near zero, so it is suppressed there.
# =============================================================================
roster <- tibble::tribble(
  ~gene,      ~arm,
  "Htra2",    "execution (IMS protease)",
  "Diablo",   "execution (IAP antagonist)",
  "Cycs",     "execution (apoptosome)",
  "Apaf1",    "execution (apoptosome)",
  "Casp9",    "execution (apoptosome)",
  "Casp3",    "execution (effector caspase)",
  "Bax",      "effector (pro)",
  "Bak1",     "effector (pro)",
  "Bbc3",     "BH3-only trigger (PUMA)",
  "Bcl2l11",  "BH3-only trigger (BIM)",
  "Pmaip1",   "BH3-only trigger (NOXA)",
  "Bid",      "BH3-only trigger",
  "Bmf",      "BH3-only trigger",
  "Bcl2l1",   "brake (Bcl-xL)",
  "Mcl1",     "brake",
  "Bcl2",     "brake",
  "Xiap",     "brake (IAP)",
  "Birc5",    "brake (IAP, survivin)",
  "Tomm20",   "biogenesis / import",
  "Tomm22",   "biogenesis / import",
  "Tomm40",   "biogenesis / import",
  "Timm23",   "biogenesis / import",
  "Timm44",   "biogenesis / import",
  "Hspd1",    "biogenesis / import",
  "Immt",     "biogenesis / import",
  "Opa1",     "biogenesis / import",
  "Atp5f1a",  "OXPHOS subunit (blot: Atp5a)",
  "Uqcrc2",   "OXPHOS subunit (blot)",
  "Ndufb8",   "OXPHOS subunit (blot)",
  "Sdhb",     "OXPHOS subunit (blot)",
  "Cox4i1",   "OXPHOS subunit",
  "mt-Co1",   "OXPHOS subunit (mtDNA, blot: Mt-Co1)",
  "Vdac1",    "outer membrane",
  "Esrra",    "biogenesis TF (PGC1a axis)",
  "Nrf1",     "biogenesis TF (PGC1a axis)",
  "Gabpa",    "biogenesis TF (PGC1a axis)",
  "Ppargc1a", "biogenesis coactivator (PGC1a -- see note, baseMean too low)",
  "Ppargc1b", "biogenesis coactivator",
  "Tfam",     "mtDNA machinery")

pull <- function(genes) {
  e  <- ens_of(genes)
  ok <- !is.na(e) & e %in% rownames(D$myc_6W_raw)
  g  <- function(contrast, col) D[[contrast]][match(e[ok], rownames(D[[contrast]])), col]
  tibble::tibble(
    gene            = genes[ok],
    baseMean        = round(g("myc_6W_raw", "baseMean")),
    lfc_myc_6W      = g("myc_6W_raw",        "log2FoldChange"),
    padj_myc_6W     = g("myc_6W_raw",        "padj"),
    lfc_myc_12W     = g("myc_12W_raw",       "log2FoldChange"),
    padj_myc_12W    = g("myc_12W_raw",       "padj"),
    lfc_wt_time     = g("timepoint_neg_raw", "log2FoldChange"),
    padj_wt_time    = g("timepoint_neg_raw", "padj"),
    lfc_myc_time    = g("timepoint_pos_raw", "log2FoldChange"),
    padj_myc_time   = g("timepoint_pos_raw", "padj"),
    lfc_interaction = g("interaction_raw",   "log2FoldChange"),
    padj_interaction= g("interaction_raw",   "padj"))
}

machinery <- pull(roster$gene) |>
  dplyr::left_join(roster, by = "gene") |>
  dplyr::mutate(
    # only interpretable where Myc actually does something at 6W
    retention = ifelse(abs(lfc_myc_6W) >= 0.2, round(lfc_myc_12W / lfc_myc_6W, 2), NA_real_),
    vs_global = ifelse(is.na(retention), NA_character_,
                       ifelse(retention < GLOBAL_RATE - 0.15, "below global",
                       ifelse(retention > GLOBAL_RATE + 0.15, "above global", "at global"))))
unresolved_roster <- setdiff(roster$gene, machinery$gene)
if (length(unresolved_roster))
  message("42 PART A: unresolved -> ", paste(unresolved_roster, collapse = ", "))

# =============================================================================
# PART B: ARM SPECIFICITY -- the headline
# -----------------------------------------------------------------------------
# Per-mouse log2 priming ratio r = log2(pro) - log2(anti). This is the quantity a
# BH3 profile approximates and it is the quantity the Bcl-xL rescue experiment
# manipulates directly, so it is mechanistically meaningful, not just convenient.
#
#   lm(r ~ tp * myc)            -> the genotype effect at 6W and its change (interaction)
#   lm(r ~ tp * myc + epi + imm) -> the same, asking whether composition explains it
#
# The primary read is the CONTRAST between arms, not any single p-value: the panel
# shares a denominator, so if BAX priming retains the global rate while PUMA priming
# collapses, "it is just the global attenuation" is answered from inside the panel.
# =============================================================================
ratio_of <- function(pro, anti) {
  ep <- ens_of(pro); ea <- ens_of(anti)
  stopifnot(!is.na(ep), !is.na(ea), ep %in% rownames(L), ea %in% rownames(L))
  as.numeric(L[ep, ]) - as.numeric(L[ea, ])
}
fit_ratio <- function(r, adjust = FALSE) {
  d  <- data.frame(r = r, tp = sm$tp, myc = sm$myc, epi = epi_comp, imm = imm_comp)
  fo <- if (adjust) r ~ tp * myc + epi + imm else r ~ tp * myc
  m  <- summary(stats::lm(fo, d))$coefficients
  d6 <- m["mycpos", "Estimate"]; i <- m["tp12W:mycpos", "Estimate"]
  c(d6 = d6, p6 = m["mycpos", "Pr(>|t|)"],
    d12 = d6 + i, int = i, int_p = m["tp12W:mycpos", "Pr(>|t|)"])
}

pair_panel <- tibble::tribble(
  ~pro,      ~anti,
  "Bax",     "Bcl2l1",
  "Bbc3",    "Bcl2l1",
  "Bcl2l11", "Bcl2l1",
  "Bak1",    "Bcl2l1",
  "Bid",     "Bcl2l1",
  "Bmf",     "Bcl2l1",
  "Pmaip1",  "Bcl2l1",
  "Bax",     "Mcl1",
  "Bbc3",    "Mcl1")

priming <- dplyr::bind_rows(lapply(seq_len(nrow(pair_panel)), function(i) {
  p <- pair_panel$pro[i]; a <- pair_panel$anti[i]
  r <- ratio_of(p, a); f0 <- fit_ratio(r, FALSE); f1 <- fit_ratio(r, TRUE)
  tibble::tibble(
    pair = paste0(p, ":", a), pro = p, anti = a,
    d6 = f0[["d6"]], p6 = f0[["p6"]], d12 = f0[["d12"]],
    retention = f0[["d12"]] / f0[["d6"]], int = f0[["int"]], int_p = f0[["int_p"]],
    d6_adj = f1[["d6"]], p6_adj = f1[["p6"]], int_adj = f1[["int"]], int_p_adj = f1[["int_p"]])
})) |>
  dplyr::mutate(int_p_bh = stats::p.adjust(int_p, "BH"))

# --- the matched-pair null ----------------------------------------------------
# Random (pro-like, anti-like) gene pairs matched on baseMean within 2-fold of the
# real members. Two readings:
#   unconditional -- how extreme is this interaction among arbitrary pairs;
#   CONDITIONAL on a 6W effect at least as large -- among pairs that start as high,
#   how many retain more at 12W. The conditional one is the real test: it asks about
#   SELECTIVE LOSS rather than about effect size.
bm_all <- rowMeans(nc)
near_pool <- function(g, fold = 2) {
  b <- bm_all[ens_of(g)]
  names(bm_all)[bm_all > b / fold & bm_all < b * fold & bm_all >= 20]
}
null_for <- function(pro, anti, n = NPAIR) {
  pa <- near_pool(pro); pb <- near_pool(anti)
  stopifnot(length(pa) >= 50, length(pb) >= 50)
  draws <- t(vapply(seq_len(n), function(i)
    fit_ratio(as.numeric(L[sample(pa, 1), ]) - as.numeric(L[sample(pb, 1), ])), numeric(5)))
  list(draws = draws, n_pro_pool = length(pa), n_anti_pool = length(pb))
}
pair_null <- dplyr::bind_rows(lapply(which(pair_panel$anti == "Bcl2l1"), function(i) {
  p <- pair_panel$pro[i]; a <- pair_panel$anti[i]
  obs <- fit_ratio(ratio_of(p, a)); nl <- null_for(p, a)
  dr  <- nl$draws
  ret_obs  <- obs[["d12"]] / obs[["d6"]]
  cond     <- dr[dr[, "d6"] >= obs[["d6"]], , drop = FALSE]
  ret_cond <- if (nrow(cond)) cond[, "d12"] / cond[, "d6"] else NA_real_
  tibble::tibble(
    pair = paste0(p, ":", a),
    retention = ret_obs,
    null_pairs = nrow(dr),
    pct_int_uncond = 100 * mean(abs(dr[, "int"]) < abs(obs[["int"]])),
    n_conditional = nrow(cond),
    null_retention_median = stats::median(ret_cond),
    pct_retention_cond = 100 * mean(ret_cond > ret_obs),
    p_emp_cond = mean(ret_cond <= ret_obs))
}))

# POSITIVE CONTROL, asserted below: the conditional null's own median retention must
# land near GLOBAL_RATE. It is an independent reconstruction of script 40's x0.55
# from random gene pairs -- if it does not, the null is not measuring what it should.
ctrl_median <- stats::median(pair_null$null_retention_median, na.rm = TRUE)

# --- Bbc3 as a single gene, against a doubly-matched null ---------------------
# script 34's aligned attenuation: sign(lfc6) * (lfc6 - lfc12). Matching on baseMean
# AND |lfc6| decile is what makes it a test of EXCESS attenuation.
keep    <- rownames(D$myc_6W_raw)[!is.na(D$myc_6W_raw$padj) & D$myc_6W_raw$baseMean >= 20]
a6      <- D$myc_6W_raw[keep, "log2FoldChange"]
a12     <- D$myc_12W_raw[match(keep, rownames(D$myc_12W_raw)), "log2FoldChange"]
att_all <- sign(a6) * (a6 - a12)
bm_keep <- D$myc_6W_raw[keep, "baseMean"]
# unique() on the breaks: tied quantiles would otherwise make cut() error out
q20     <- function(x) cut(x, unique(stats::quantile(x, 0:20 / 20)), include.lowest = TRUE)
db      <- q20(bm_keep); dl <- q20(abs(a6))
single_gene_null <- dplyr::bind_rows(lapply(c("Bbc3", "Bax", "Bcl2l1", "Htra2"), function(g) {
  i <- match(ens_of(g), keep)
  if (is.na(i)) return(NULL)
  sel <- which(db == db[i] & dl == dl[i])
  tibble::tibble(gene = g, lfc6 = a6[i], lfc12 = a12[i], aligned_attenuation = att_all[i],
                 matched_n = length(sel), matched_median = stats::median(att_all[sel]),
                 percentile = 100 * mean(att_all[sel] < att_all[i]),
                 p_emp = mean(att_all[sel] >= att_all[i]))
}))

# =============================================================================
# PART C: DOES THE PGC1a / BIOGENESIS AXIS REACH PUMA?
# -----------------------------------------------------------------------------
# Three lines, of decreasing strength, reported with their weight attached.
#   C1  the Gray/ChEA shortlist -- WHICH TF x context regulons contain Bbc3, and which
#       biogenesis TFs carry pro-apoptotic effectors. READ from the library CSV, not
#       re-derived. This is the in-silico support for the co-regulation hypothesis.
#   C2  the ESRRA/NRF1/GABPA GSVA lanes -- is the biogenesis-TF programme Myc-driven,
#       and does it attenuate faster than the global rate?
#   C3  the per-sample OXPHOS <-> PUMA coupling with an ambient null. THIS IS A
#       NEGATIVE. It is computed and kept so the paper cannot quietly claim it.
# Ppargc1a itself has baseMean ~30 in this tissue and is NOT readable; ESRRA/NRF1
# lanes are the activity proxy (as in script 38 PART C).
# =============================================================================
# the library snapshot is the canonical copy (CLAUDE.md); the docs/library_reference
# copy is byte-identical but untracked, so read the tracked one
shortlist_path <- here::here("data", "genesets_from_library", "gray_chea_mito_tf_shortlist.csv")
tf_puma <- NULL; tf_biogenesis <- NULL
if (file.exists(shortlist_path)) {
  sl <- utils::read.csv(shortlist_path, stringsAsFactors = FALSE)
  tf_puma <- sl[grepl("Bbc3", sl$apoptosis_pro_genes, fixed = TRUE),
                c("TF", "context", "n_mito_total", "n_oxphos", "n_apoptosis_pro",
                  "apoptosis_pro_genes", "bh_fdr")]
  # READ THIS AS A CONTEXT STATEMENT, NOT A TF STATEMENT. Bbc3 turns up in ~72 TF
  # regulons and every one of them is AP_TEB -- it appears in NO other context. That
  # is what a broadly-bound promoter in one cell context looks like, not evidence that
  # any particular TF targets PUMA. The specific claim that does survive is in
  # tf_biogenesis: MYC/E2F1 reach Bbc3 ONLY in AP_TEB, and reach Bid/Bok/Bak1/Pmaip1
  # in every other context -- so which BH3-only MYC can touch is context-dependent.
  tf_puma_contexts <- sort(unique(tf_puma$context))
  tf_biogenesis <- sl[sl$TF %in% c("ESRRA", "ESRRB", "ESRRG", "NRF1", "GABPA", "MYC", "E2F1") &
                      sl$n_apoptosis_pro > 0,
                      c("TF", "context", "n_mito_total", "n_oxphos", "n_apoptosis_pro",
                        "apoptosis_pro_genes", "bh_fdr")]
} else {
  tf_puma_contexts <- character(0)
  message("42 PART C: shortlist CSV not found in data/genesets_from_library -- C1 skipped")
}

myc_effect_by_age <- function(s) {          # genotype effect within each timepoint
  d <- data.frame(s = s, myc = sm$myc, epi = epi_comp, imm = imm_comp)
  f <- function(k) summary(stats::lm(s ~ myc + epi + imm, d[k, ]))$coefficients["mycpos", ]
  m6 <- f(sm$tp == "6W"); m12 <- f(sm$tp == "12W")
  c(myc6 = m6[["Estimate"]], p6 = m6[["Pr(>|t|)"]],
    myc12 = m12[["Estimate"]], p12 = m12[["Pr(>|t|)"]])
}
bio_lanes <- grep("^TFT_(ESRRA|NRF1|GABPA|MYC|E2F1)_GRAY", names(gmt), value = TRUE)
tf_lane_effects <- dplyr::bind_rows(lapply(bio_lanes, function(s) {
  v <- set_score(s); if (is.null(v)) return(NULL)
  f <- myc_effect_by_age(v)
  tibble::tibble(lane = s, n_genes = length(gmt[[s]]),
                 myc_6W = f[["myc6"]], p_6W = f[["p6"]],
                 myc_12W = f[["myc12"]], p_12W = f[["p12"]],
                 retention = ifelse(abs(f[["myc6"]]) >= 0.2, f[["myc12"]] / f[["myc6"]], NA_real_))
})) |> dplyr::arrange(dplyr::desc(myc_6W))

# C3 -- the negative. OXPHOS composite vs PUMA / BAX, per timepoint, against the
# ambient distribution of couplings to arbitrary expressed genes (script 35 logic).
ox_ens  <- recon_to_ensembl(gmt[["MITOCARTA_OXPHOS_NU"]], rownames(L))
ox_ens  <- ox_ens[!is.na(ox_ens)]
ox_comp <- colMeans(zrow(L[ox_ens, , drop = FALSE]))
amb_pool <- rownames(L)[rowMeans(nc) >= 20]
coupling_null <- dplyr::bind_rows(lapply(c("6W", "12W"), function(tpv) {
  k   <- sm$tp == tpv
  amb <- vapply(sample(amb_pool, NAMB), function(g)
    stats::cor(ox_comp[k], as.numeric(L[g, k]), method = "spearman"), numeric(1))
  dplyr::bind_rows(lapply(c("Bbc3", "Bax", "Bcl2l1"), function(g) {
    r <- stats::cor(ox_comp[k], as.numeric(L[ens_of(g), k]), method = "spearman")
    tibble::tibble(timepoint = tpv, gene = g, rho = r,
                   ambient_median_abs = stats::median(abs(amb)),
                   percentile = 100 * mean(amb < r))
  }))
}))

# =============================================================================
# PART D: THE DEVELOPMENTAL SUBSTRATE -- does the TEB regress?
# -----------------------------------------------------------------------------
# The MYC-ER experiment demands a substrate change that is DEVELOPMENTAL (those mice
# never saw Myc). The library documents that the HS_TEB context is apoptosis-loaded
# ("TEB lumen clearance by apoptosis during pubertal elongation",
# docs/library_reference/Gray_et_al_developmental_TFS_selection.md), so the compartment
# to look at is the terminal end bud. Purity-adjusted throughout.
#
# BATCH CAVEAT, restated because this PART is the one it bites: WT 6->12W is the
# batch-confounded axis. These are DESCRIBED observations. Their one mitigation is
# that the mito and PUMA protein blots move the same way, off the RNA batch entirely.
# =============================================================================
score_and_fit <- function(sets) {
  dplyr::bind_rows(lapply(sets, function(s) {
    v <- set_score(s); if (is.null(v)) return(NULL)
    d <- data.frame(s = v, tp = sm$tp, myc = sm$myc, epi = epi_comp, imm = imm_comp)
    m <- summary(stats::lm(s ~ tp * myc + epi + imm, d))$coefficients
    tibble::tibble(set = s, n_genes = length(gmt[[s]]),
                   wt_time = m["tp12W", "Estimate"],   wt_p = m["tp12W", "Pr(>|t|)"],
                   myc_6W  = m["mycpos", "Estimate"],  myc6_p = m["mycpos", "Pr(>|t|)"],
                   interaction = m["tp12W:mycpos", "Estimate"])
  }))
}
teb_sets   <- grep("TEB", names(gmt), value = TRUE)
state_sets <- c("MG_LASP_CONSENSUS", "MG_LHS_CONSENSUS", "MG_BMYO_CONSENSUS",
                "MG_LP_PAL2017_CIII", "MG_LUMINAL_PROGENITOR_GARCIASOLA", "POMMIER_LP_UP",
                "MG_LUMINAL_ALVPROG_GARCIASOLA", "MG_BASAL_GRAY", "MG_LUMINAL_PAL2017",
                "MG_BASAL_MYOEPITHELIAL_GARCIASOLA", "MG_LUMINAL_HORMSENSDIF_GARCIASOLA")
teb_signatures <- score_and_fit(teb_sets)   |> dplyr::arrange(wt_time)
mec_states     <- score_and_fit(state_sets) |> dplyr::arrange(wt_time)

# Is PUMA a TEB transcript? (Expected: no. The TEB claim is about TF-binding context.)
puma_in_teb <- tibble::tibble(
  set = grep("TEB_VS_DUCTAL", names(gmt), value = TRUE)) |>
  dplyr::mutate(bbc3_member = vapply(set, function(s) "Bbc3" %in% gmt[[s]], logical(1)))

substrate_markers <- pull(c("Aldh1a3", "Nrg1", "Areg", "Esr1", "Sox9", "Cited1", "Wnt4",
                            "Oxtr", "Myh11", "Trp63", "Acta2", "Krt5", "Krt14", "Krt8",
                            "Krt18", "Elf5", "Prlr", "Epcam"))

# =============================================================================
# PART E: THE EXCLUSIONS -- what is NOT the mechanism
# -----------------------------------------------------------------------------
# Each of these narrows the claim. The buffer does not rise, p53 does not move, none
# of PUMA's known transcriptional inputs tracks it, and bulk RNA cannot see the death
# at all -- which also means the 6W pro-apoptotic induction is UNDER-reported (the
# dying cells are not in the library), so the PUMA collapse is a LOWER BOUND. That is
# the in-vivo analogue of why Bcl-xL was needed in vitro to reveal the full induction.
# =============================================================================
exclusions <- list(
  buffer          = pull(c("Bcl2l1", "Mcl1", "Bcl2", "Bcl2l2", "Bcl2a1b", "Xiap",
                           "Birc2", "Birc3", "Birc5")),
  p53_axis        = pull(c("Trp53", "Mdm2", "Cdkn1a", "Cdkn2a", "Trp53inp1", "Zmat3",
                           "Eda2r", "Phlda3", "Ccng1")),
  puma_inputs     = pull(c("E2f1", "Trp73", "Atf4", "Ddit3", "Trib3", "Chac1",
                           "Foxo1", "Foxo3", "Foxo4", "Sesn2", "Eif2ak3", "Nupr1")),
  replication     = pull(c("Atr", "Chek1", "Chek2", "Rad51", "Fancd2", "Brca1", "Brca2",
                           "Rpa1", "Cdc25a", "Wee1", "Clspn", "Exo1", "Blm", "Mre11a",
                           "Atm", "Parp1")),
  efferocytosis   = pull(c("Mertk", "Gas6", "Axl", "Trem2", "C1qa", "C1qb", "C1qc",
                           "Cd68", "Adgre1", "Itgav", "Timd4", "Stab1", "Msr1", "Csf1r")),
  proliferation   = pull(c("Mki67", "Top2a", "Ccnb1", "Ccna2", "Bub1", "Plk1", "Aurka",
                           "Cdk1", "Pcna", "Mcm2", "Rrm2", "Tk1")))
exclusions <- lapply(exclusions, function(x)
  dplyr::mutate(x, retention = ifelse(abs(lfc_myc_6W) >= 0.2,
                                      round(lfc_myc_12W / lfc_myc_6W, 2), NA_real_)))

# =============================================================================
# PART G: THE TWO-RULER OXPHOS DECLINE -- evidence strand 1
# -----------------------------------------------------------------------------
# figS6's reading, tabulated. On the PRIORITY ruler (mitoPPS) OXPHOS is the most
# de-prioritised Level-1 tier on BOTH temporal axes, and within it the STRUCTURAL
# SUBUNITS carry it while the assembly factors do not move. The subunits are the
# PGC1a / NRF1 / ERRalpha output, so what loses priority is exactly the arm PGC1a
# builds.
#
# WHY THIS STRAND IS WORTH HAVING, structurally: mitoPPS is a PAIRWISE RATIO within
# the mitochondrial compartment, so a uniform x0.55 scaling of the whole programme
# CANCELS in it. The content ruler is dose-dominated; the priority ruler is not. An
# OXPHOS decline that survives on the priority ruler is therefore NOT the dose effect.
#
# The apoptosis pathways are put on the same table deliberately: their composite
# balance is FLAT over the timeline (PRO and ANTI move together), which is why
# script 34's composite analysis found nothing and why the PUMA result had to be
# gene-level. That is also what the cell work predicts -- PGC1a induces PUMA "but not
# other BH3-only proteins", so a 25-gene PRO composite is the wrong instrument.
# =============================================================================
bvm_path  <- here::here("results", "background_vs_myc.rds")
two_ruler_tier <- NULL; two_ruler_pathway <- NULL
if (file.exists(bvm_path)) {
  ruler <- tibble::as_tibble(readRDS(bvm_path)$ruler)
  two_ruler_tier <- ruler |>
    dplyr::filter(!is_mtdna) |>
    dplyr::group_by(tier) |>
    dplyr::summarise(
      n_pathways      = dplyr::n(),
      content_myc_6W  = stats::median(c_m6),  content_wt_time  = stats::median(c_tn),
      content_myc_time= stats::median(c_tp),
      prio_myc_6W     = stats::median(p_m6),  prio_wt_time     = stats::median(p_tn),
      prio_myc_time   = stats::median(p_tp),  .groups = "drop") |>
    dplyr::arrange(prio_myc_time)
  focus <- c("OXPHOS", "OXPHOS subunits", "OXPHOS assembly factors",
             "mtDNA-encoded OXPHOS subunits", "Apoptosis", "Apoptosis-PRO",
             "Apoptosis-ANTI")
  two_ruler_pathway <- ruler |>
    dplyr::filter(pathway %in% focus) |>
    dplyr::select(pathway, tier, n_genes,
                  c_m6, c_m12, c_tn, c_tp, p_m6, p_m12, p_tn, p_tp)
} else {
  message("42 PART G: results/background_vs_myc.rds absent -- run script 40 first")
}

# =============================================================================
# PART H: THE COINCIDENCE MODEL, IN THE SPACE WHERE IT IS READABLE
# -----------------------------------------------------------------------------
# The model (the author's, 2026-07-25): MYC-driven death needs TWO inputs -- MYC, and
# a competent mitochondrial substrate. The substrate input falls developmentally, so
# the SAME MYC no longer clears the threshold. This is not a dose model: the MYC-ER
# experiment holds MYC fixed and varies the substrate, so it TESTS this model rather
# than excluding it.
#
# It also earns its place on three counts. It explains the death >> proliferation
# asymmetry with no extra assumption (proliferation needs one input, death needs two,
# so death falls superlinearly); a THRESHOLD on the product explains why PUMA priming
# collapses THROUGH ZERO rather than merely attenuating; and it is a literal statement
# of the manuscript's title thesis -- mitochondria INTEGRATING an oncogenic and a
# metabolic input.
#
# Four candidate second inputs, plus a control:
#   oxphos_ppd  mitoPPS OXPHOS-subunit ratio  -- content-blind and dose-cancelling
#   oxphos_lvl  log-expression OXPHOS composite -- the LEVEL, not the priority
#   teb         MG_TEB_VS_DUCTAL_HS_GRAY_UP    -- the compositional/context reading
#   redox_ppd   mitoPPS ROS-and-glutathione    -- NEGATIVE CONTROL (fails in script 38)
#
# WHAT THIS CANNOT DO, stated before the results: it tests the biogenesis LEVEL and
# PRIORITY, not PGC1a ACTIVITY. Ppargc1a is baseMean ~30 here and unreadable, and a
# coactivator's activity is post-translational. A null here bounds the level version
# and leaves the activity version untouched -- for which the cell perturbations are
# the evidence, and they are perturbations, which outrank any correlation below.
# EVERYTHING IN THIS PART IS HYPOTHESIS-GENERATING.
# =============================================================================
mp_path <- here::here("results", "mitopps_scores.rds")
axes <- list(oxphos_lvl = ox_comp, teb = set_score("MG_TEB_VS_DUCTAL_HS_GRAY_UP"))
if (file.exists(mp_path)) {
  mp  <- readRDS(mp_path)
  mps <- mp$mitopps_scores
  mps <- mps[match(colnames(L), mps$sample), , drop = FALSE]
  stopifnot(identical(as.character(mps$sample), colnames(L)))
  # exact column names (script 08 PART 2c); never grep -- "OXPHOS" also prefixes others
  need <- c("OXPHOS subunits", "ROS and glutathione metabolism",
            "Apoptosis-PRO", "Apoptosis-ANTI")
  missing_mps <- setdiff(need, names(mps))
  if (length(missing_mps)) {
    message("42 PART H: mitoPPS columns absent -> ", paste(missing_mps, collapse = ", "))
  } else {
    axes$oxphos_ppd <- as.numeric(mps[["OXPHOS subunits"]])
    axes$redox_ppd  <- as.numeric(mps[["ROS and glutathione metabolism"]])
    # the mitoPPS-space priming ratio, for continuity with script 38 PART B
    axes_priming_ppd <- as.numeric(mps[["Apoptosis-PRO"]]) - as.numeric(mps[["Apoptosis-ANTI"]])
  }
} else {
  message("42 PART H: results/mitopps_scores.rds absent -- run script 08 first")
}
axes <- axes[!vapply(axes, is.null, logical(1))]

ratios <- list(`Bbc3:Bcl2l1` = ratio_of("Bbc3", "Bcl2l1"),
               `Bax:Bcl2l1`  = ratio_of("Bax",  "Bcl2l1"))

# H1 -- absorption, in script 30's idiom (30_attenuation_moderation.R:174 absorb_one).
# CAVEAT CARRIED FROM SCRIPT 30 VERBATIM: the mediator here is ENDOGENOUS and itself
# Myc-driven, so this is biased. Absorption BOUNDS a claim; it does not prove mediation.
base_int <- function(r) {
  d <- data.frame(r = r, tp = sm$tp, myc = sm$myc, epi = epi_comp, imm = imm_comp)
  summary(stats::lm(r ~ tp * myc + epi + imm, d))$coefficients["tp12W:mycpos", c(1, 4)]
}
coincidence_absorption <- dplyr::bind_rows(lapply(names(ratios), function(rn) {
  r <- ratios[[rn]]; b <- base_int(r)
  dplyr::bind_rows(lapply(names(axes), function(an) {
    d <- data.frame(r = r, tp = sm$tp, myc = sm$myc, epi = epi_comp, imm = imm_comp,
                    a = axes[[an]])
    m <- summary(stats::lm(r ~ tp * myc + epi + imm + a + tp:a, d))$coefficients
    tibble::tibble(ratio = rn, axis = an,
                   int_base = b[[1]], p_base = b[[2]],
                   int_adj = m["tp12W:mycpos", 1], p_adj = m["tp12W:mycpos", 4],
                   absorbed_frac = 1 - m["tp12W:mycpos", 1] / b[[1]])
  }))
}))

# H2 -- the coincidence test proper: does Myc raise priming ONLY where the substrate
# is competent? H3 -- and does that survive the timepoint it is confounded with?
coincidence_fit <- dplyr::bind_rows(lapply(names(ratios), function(rn) {
  r <- ratios[[rn]]
  dplyr::bind_rows(lapply(names(axes), function(an) {
    d  <- data.frame(r = r, tp = sm$tp, myc = sm$myc, epi = epi_comp, imm = imm_comp,
                     a = axes[[an]])
    m1 <- summary(stats::lm(r ~ myc * a + epi + imm, d))$coefficients
    m2 <- summary(stats::lm(r ~ tp * myc + myc:a + a + epi + imm, d))$coefficients
    tibble::tibble(ratio = rn, axis = an,
                   myc_x_axis = m1["mycpos:a", 1], p = m1["mycpos:a", 4],
                   myc_x_axis_with_tp = m2["mycpos:a", 1], p_with_tp = m2["mycpos:a", 4],
                   myc_x_tp_with_axis = m2["tp12W:mycpos", 1])
  }))
}))

# H3b -- within-timepoint permutation of the axis score. Shuffling INSIDE each
# timepoint preserves the timepoint structure and breaks only the mouse-to-mouse link,
# which is the thing being claimed.
perm_axis <- dplyr::bind_rows(lapply(names(ratios), function(rn) {
  r <- ratios[[rn]]
  dplyr::bind_rows(lapply(names(axes), function(an) {
    d   <- data.frame(r = r, myc = sm$myc, epi = epi_comp, imm = imm_comp, a = axes[[an]])
    obs <- summary(stats::lm(r ~ myc * a + epi + imm, d))$coefficients["mycpos:a", 1]
    nul <- vapply(seq_len(NPERM), function(i) {
      dd <- d
      for (tv in levels(sm$tp)) { k <- which(sm$tp == tv); dd$a[k] <- sample(dd$a[k]) }
      summary(stats::lm(r ~ myc * a + epi + imm, dd))$coefficients["mycpos:a", 1]
    }, numeric(1))
    tibble::tibble(ratio = rn, axis = an, observed = obs,
                   null_median = stats::median(nul),
                   percentile = 100 * mean(nul < obs), p_emp = mean(nul >= obs))
  }))
}))

# H4 -- substrate magnitude. Is the developmental fall big enough to matter next to
# the Myc effect? On RNA it is not obviously so; this table is where that shows.
substrate_magnitude <- dplyr::bind_rows(lapply(names(axes), function(an) {
  v <- axes[[an]]; g <- function(k) mean(v[sm$group == k])
  tibble::tibble(axis = an,
                 m_6W_neg = g("6W_neg"), m_6W_pos = g("6W_pos"),
                 m_12W_neg = g("12W_neg"), m_12W_pos = g("12W_pos"),
                 wt_developmental_fall = g("12W_neg") - g("6W_neg"),
                 myc_effect_6W = g("6W_pos") - g("6W_neg"),
                 fall_vs_myc = (g("12W_neg") - g("6W_neg")) / (g("6W_pos") - g("6W_neg")))
}))

# =============================================================================
# ASSERTS
# =============================================================================
stopifnot(
  nrow(priming) == nrow(pair_panel),                       # panel complete
  all(c("Bbc3", "Bax", "Bcl2l1", "Htra2") %in% machinery$gene),
  "Bbc3" %in% single_gene_null$gene,
  # the two arms the headline compares must both have a conditional null
  all(!is.na(pair_null$null_retention_median[pair_null$pair %in%
        c("Bax:Bcl2l1", "Bbc3:Bcl2l1")])))
# a pair with a very large 6W effect can leave the conditional pool empty; that is a
# missing comparison, not an error -- name it rather than dropping it silently
if (any(is.na(pair_null$null_retention_median)))
  message("42 PART B: no conditional null (6W effect too extreme) for -> ",
          paste(pair_null$pair[is.na(pair_null$null_retention_median)], collapse = ", "))
stopifnot(
  nrow(coincidence_fit) == length(ratios) * length(axes),
  all(is.finite(perm_axis$p_emp)))
# the control axis must be present, or PART H has no negative control and its
# positives cannot be read
if (!"redox_ppd" %in% names(axes))
  warning("42 PART H: redox_ppd control axis missing -- read the myc:axis terms with care")
if (is.null(two_ruler_tier))
  warning("42 PART G: no two-ruler table (script 40 output absent)")
# the built-in positive control: random matched pairs must reproduce script 40's rate
if (!is.finite(ctrl_median) || abs(ctrl_median - GLOBAL_RATE) > 0.20)
  warning(sprintf(paste("42 PART B: matched-pair null median retention %.2f is far from the",
                        "global rate %.2f -- the null may not be measuring the attenuation."),
                  ctrl_median, GLOBAL_RATE))

# =============================================================================
# SAVE
# =============================================================================
notes <- c(
  "DOSE vs COMPETENCE. The MMTV-Myc 6W->12W timeline is a DOSE experiment: MYC protein",
  "falls ~50% at constant transcript, and the transcriptome rescales x0.55 with no",
  "reshaping (script 40), no programme in excess (script 34) and no repressor (the",
  "2026-07-25 scans). Output per unit MYC is unchanged. The MYC-ER mouse is a COMPETENCE",
  "experiment: equal or higher MYC protein, ~80% less death, ~30% less proliferation --",
  "and since those animals never saw MYC before induction, the substrate change cannot be",
  "MYC-driven selection. It is developmental. The two are not rival explanations; they",
  "answer different questions.",
  "",
  "PART B is the headline and it is a CONTRAST, not a p-value: within one panel sharing",
  "the Bcl-xL denominator, Bax:Bcl-xL priming retains the global rate while PUMA:Bcl-xL",
  "priming collapses. PUMA is the BH3-only the PGC1a cell experiments single out, so the",
  "pair was pre-specified externally -- but the 9-pair panel is exploratory context and",
  "int_p_bh is reported for it.",
  "",
  "THREE NEGATIVES, to be quoted as prominently as the positives:",
  " 1. PART C3 -- the per-sample OXPHOS<->PUMA coupling sits at or below the ambient. The",
  "    OXPHOS-PUMA bundle is carried by the CELL PERTURBATIONS, never by in-vivo",
  "    correlation. Do not claim it from these data.",
  " 2. PART D -- Bbc3 is NOT in the TEB-vs-ductal signature. The TEB claim is about",
  "    TF-binding CONTEXT (Gray/ChEA AP_TEB regulons contain Bbc3), not about PUMA being",
  "    a TEB transcript.",
  " 3. PART D -- the MEC state composites are not significant. Only markers move, so the",
  "    substrate change is documented at marker level, not as a compartment shift.",
  "",
  "CEILING. n=6/group, exploratory, ranking not confirmatory. WT temporal statements are",
  "on the BATCH-CONFOUNDED axis (batch = timepoint) and are DESCRIBED, not claimed; their",
  "mitigation is that the protein blots move the same way off the RNA batch entirely.",
  "Script 34's circularity caveat still binds COMPOSITE death scores (Bbc3, Bax and",
  "Bcl2l11 are all MitoCarta APOPTOSIS members) -- this script makes no composite",
  "correlation claim. Ppargc1a is baseMean ~30 here and is NOT readable; ESRRA/NRF1/GABPA",
  "lanes are the activity proxy, as in script 38 PART C.",
  "",
  "PART E bounds the whole enterprise: there is no efferocytosis signature at 6W, so bulk",
  "RNA cannot see the death. It also UNDER-reports the 6W pro-apoptotic induction, because",
  "the cells that died are not in the library -- making the PUMA collapse a lower bound.",
  "",
  "THE EPISTEMIC CONTRACT for PARTS G and H. The in-vivo transcriptome GENERATES the",
  "hypothesis; the cell perturbations (PGC1a gain sensitises; SS and passaged cells lose",
  "biogenesis and resist) PROVE it. No correlation here is asked to carry a causal claim.",
  "",
  "PART G -- why the priority ruler earns its place. mitoPPS is a PAIRWISE RATIO within the",
  "mitochondrial compartment, so a uniform x0.55 dose scaling CANCELS in it. The content",
  "ruler is dose-dominated; the priority ruler is not. OXPHOS is nevertheless the most",
  "de-prioritised Level-1 tier on BOTH temporal axes, carried by the structural SUBUNITS",
  "while the assembly factors do not move -- and the subunits are the PGC1a/NRF1/ERRalpha",
  "output. So the arm PGC1a builds loses priority, and that is NOT the dose effect.",
  "Note in the same table that the apoptosis COMPOSITE balance is flat over the timeline",
  "(PRO and ANTI move together): the death finding is gene-specific to PUMA, which is why",
  "script 34's composite analysis found nothing and is exactly what the cell work predicts",
  "(PGC1a induces PUMA 'but not other BH3-only proteins').",
  "",
  "PART H -- the coincidence model. Death needs MYC AND a competent substrate; the substrate",
  "input falls developmentally. This is NOT a dose model, and MYC-ER (fixed MYC, varied",
  "substrate) TESTS it rather than excluding it. Three things it buys: the death >>",
  "proliferation asymmetry follows with no extra assumption (one input vs two); a THRESHOLD",
  "on the product explains why PUMA priming collapses THROUGH ZERO rather than attenuating;",
  "and it is the manuscript's title thesis stated literally.",
  "WHAT PART H CANNOT DO: it tests the biogenesis LEVEL and PRIORITY, never PGC1a ACTIVITY.",
  "Ppargc1a is baseMean ~30 here; a coactivator's activity is post-translational. A null",
  "bounds the level version and leaves the activity version untouched.",
  "Two findings to carry forward with their weight attached: (a) the biogenesis LEVEL does",
  "NOT mediate the PUMA collapse, and adjusting for it can make the interaction stronger --",
  "whereas Bax priming IS largely absorbed by it, so in vivo BAX tracks mitochondrial mass",
  "and PUMA does not, which is the opposite of the cell result and must be said; (b) nothing",
  "in PART H separates from timepoint at n=24 -- the axis x Myc terms shrink once tp*myc is",
  "in the model. Read PART H as RANKING candidate second inputs, nothing more.",
  "",
  "HOW THE FIRST RUN (2026-07-25) READ, so the object is not over-interpreted later:",
  " - PART B control passed: matched pairs retain 0.55, script 40's global rate exactly.",
  " - PART B specificity is GRADED, not binary. Bmf and Bcl2l11 have MORE extreme",
  "   retentions than PUMA but NON-SIGNIFICANT 6W effects (p6 0.41 and 0.46), and a",
  "   retention is a ratio of two noisy quantities -- you cannot lose an effect you never",
  "   had. Among the pairs with a real 6W effect (Bax, Bbc3, Bak1, Bid, Bax:Mcl1), PUMA is",
  "   the ONLY one that reverses sign, and the only conditional p_emp below 0.1.",
  " - PART C1 is a CONTEXT statement: Bbc3 sits in ~72 AP_TEB regulons and in NO other",
  "   context. That is promoter accessibility in one cell context, not TF specificity.",
  " - PART H ABSORPTION IS UNINFORMATIVE -- THE CONTROL FAILED. redox_ppd 'absorbs' 41% of",
  "   the PUMA interaction, more than any real axis. Adding any axis with the right noise",
  "   structure moves the interaction, so absorption here measures nothing. Do not cite it.",
  " - PART H COINCIDENCE IS THE POSITIVE, and its control WORKED: myc x oxphos_ppd on PUMA",
  "   priming is +2.79 (p 0.005), the largest of any axis; it partially survives tp*myc",
  "   (+2.06, p 0.088) where the TEB axis does not (1.15 -> 0.53); the within-timepoint",
  "   permutation puts it at the 92nd percentile (p_emp 0.081); and redox_ppd is null",
  "   throughout (p 0.72, permutation 51st percentile). The mitoPPS OXPHOS-priority axis is",
  "   the best-supported candidate second input.",
  " - PART H4 answers the size objection: on the LEVEL ruler the WT developmental fall is",
  "   26% of the Myc effect, but on the PRIORITY ruler it is 97% of it. In the space where",
  "   dose cancels, the developmental decline is as large as the oncogene's own effect.",
  " - Htra2 is a useful internal control in single_gene_null: the strongest Myc-induced",
  "   apoptotic gene attenuates LESS than matched genes (37th pct), so the PUMA result is",
  "   not a property of Myc-induced apoptotic genes in general.")

out <- list(
  machinery         = machinery,
  priming           = priming,
  pair_null         = pair_null,
  single_gene_null  = single_gene_null,
  tf_puma           = tf_puma,
  tf_puma_contexts  = tf_puma_contexts,
  tf_biogenesis     = tf_biogenesis,
  tf_lane_effects   = tf_lane_effects,
  coupling_null     = coupling_null,
  teb_signatures    = teb_signatures,
  mec_states        = mec_states,
  puma_in_teb       = puma_in_teb,
  substrate_markers = substrate_markers,
  exclusions        = exclusions,
  two_ruler_tier    = two_ruler_tier,
  two_ruler_pathway = two_ruler_pathway,
  coincidence_absorption = coincidence_absorption,
  coincidence_fit   = coincidence_fit,
  perm_axis         = perm_axis,
  substrate_magnitude = substrate_magnitude,
  axis_scores       = tibble::as_tibble(c(list(sample = colnames(L),
                                               group = as.character(sm$group)), axes)),
  purity            = tibble::tibble(sample = colnames(L), group = sm$group,
                                     epithelial = epi_comp, immune = imm_comp),
  params            = list(NPAIR = NPAIR, NAMB = NAMB, NPERM = NPERM,
                           GLOBAL_RATE = GLOBAL_RATE,
                           null_control_median = ctrl_median,
                           axes_used = names(axes),
                           unresolved = unresolved_roster),
  analysis_date     = Sys.Date(),
  notes             = notes)

saveRDS(out, here::here("results", "priming_arm_teb.rds"))
message("Saved results/priming_arm_teb.rds")
message(sprintf("42: matched-pair null control median retention = %.2f (target ~%.2f)",
                ctrl_median, GLOBAL_RATE))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  ## PART A -- what Myc builds, and how much of it survives to 12W
  machinery |>
    dplyr::select(gene, arm, baseMean, lfc_myc_6W, padj_myc_6W, lfc_myc_12W,
                  retention, vs_global) |>
    dplyr::arrange(arm, dplyr::desc(lfc_myc_6W)) |>
    print(n = 40)

  ## PART B -- THE HEADLINE. Read the two Bcl-xL rows against each other.
  priming |>
    dplyr::select(pair, d6, p6, d12, retention, int_p, int_p_bh, d6_adj, p6_adj, int_p_adj) |>
    print(n = 12)
  pair_null |> print(n = 12)
  single_gene_null |> print()
  cat(sprintf("\nnull control: matched pairs retain %.2f (script 40 global rate %.2f)\n",
              ctrl_median, GLOBAL_RATE))

  ## PART C -- does the biogenesis axis reach PUMA?
  tf_puma |> print()            # expected: MYC / MYCN / E2F1 / KDM5B / MITF, all AP_TEB
  tf_biogenesis |> print()      # expected: NRF1 -> Bbc3, GABPA -> Bak1, ESRRA -> Aifm2
  tf_lane_effects |> print(n = 30)
  coupling_null |> print(n = 12)   # THE NEGATIVE -- at or below ambient

  ## PART D -- the developmental substrate
  teb_signatures |> dplyr::filter(wt_p < 0.1) |> print(n = 30)
  mec_states |> print(n = 12)      # expected: nothing significant
  puma_in_teb |> print()           # expected: all FALSE
  substrate_markers |>
    dplyr::select(gene, baseMean, lfc_wt_time, padj_wt_time, lfc_myc_time, lfc_myc_6W) |>
    dplyr::arrange(lfc_wt_time) |> print(n = 20)

  ## PART E -- the exclusions
  lapply(exclusions, function(x)
    dplyr::select(x, gene, lfc_myc_6W, padj_myc_6W, lfc_myc_12W, retention,
                  lfc_wt_time, padj_wt_time))

  ## PART G -- the two-ruler decline. Read prio_wt_time / prio_myc_time: OXPHOS should
  ## be the most negative tier on both, and a ratio ruler cannot show the dose effect.
  two_ruler_tier |> print(n = 10)
  two_ruler_pathway |> print(n = 10)   # subunits move, assembly does not; PRO ~ ANTI

  ## PART H -- the coincidence model. Read the redox_ppd rows as the control: if they
  ## behave like the OXPHOS rows, the OXPHOS rows say nothing.
  coincidence_absorption |> print(n = 20)
  coincidence_fit        |> print(n = 20)
  perm_axis              |> print(n = 20)
  substrate_magnitude    |> print(n = 10)

  cat(notes, sep = "\n")
}
