# =============================================================================
# 43_substrate_specificity_and_tradeoff.R
# -----------------------------------------------------------------------------
# THE TRADE-OFF MODEL. Mitochondrial respiratory capacity is both an ASSET and a
# LIABILITY to a MYC-driven cell: it supports growth, and it primes death. The
# cell arm shows the sign reversal directly -- PGC1a KILLS MYAZ cells, while
# PGC1a + Bcl-xL grows BETTER than Bcl-xL alone. Which limb is expressed depends
# entirely on whether death can be executed. This script asks the two in-vivo
# questions that the model needs, and it is scoped to exactly those two.
#
# PART A -- SUBSTRATE SPECIFICITY. The adult gland withdraws from the respiratory
#   arm. Is that withdrawal separable from a PROLIFERATIVE withdrawal, or is it
#   just the gland leaving the pubertal/TEB proliferative state? This matters for
#   the two-input model: if respiration falls only because proliferation falls,
#   the model collapses into "less proliferation, less death" and the 80%-vs-30%
#   asymmetry has no explanation. Read-only reconciliation (2026-07-26) says the
#   arms dissociate; PART A puts a NULL under that so it is a result and not a
#   description.
#
# PART B -- THE TRADE-OFF ASYMMETRY. In mitoPPS ratio space -- the one space this
#   project has found readable (ambient ~0.42, not the ~0.8 of GSVA space; script
#   38 PART A) -- does respiratory PRIORITY track death competence more tightly
#   than it tracks proliferation? That asymmetry is what makes respiration a
#   trade-off rather than a general growth input.
#
# WHAT THIS CANNOT DO, stated before the results. `batch = timepoint` (CLAUDE.md):
# every wild-type temporal statement here is DESCRIBED, not claimed, and its only
# mitigation is that the protein blots move the same way off the RNA batch. And
# nothing at n=24 separates a candidate axis from the timepoint it is confounded
# with -- script 42 PART H already showed the best term dropping p 0.005 -> 0.088
# once `tp*myc` is in the model. Read both parts as RANKING, under the standing
# epistemic contract: the in-vivo transcriptome GENERATES the hypothesis, the cell
# perturbations PROVE it.
#
# Reads : results/interaction_results.rds, results/dds_int_run.rds,
#         results/combined_df_annotated.rds, results/mitopps_scores.rds,
#         results/attenuation_decomposition.rds (positive control),
#         results/background_vs_myc.rds (priority ruler),
#         data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt,
#         functions/reconcile_gene_symbols.R (MANDATORY -- vintage-aware membership)
# Writes: results/substrate_specificity_tradeoff.rds
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

set.seed(1)
NSET  <- 2000L   # matched-random-set draws for PART A's null
NPERM <- 5000L   # within-timepoint shuffles for PART B's axis null
NBIN  <- 20L     # baseMean bins for the matched null

# =============================================================================
# PART 0: LOAD
# =============================================================================
ir  <- readRDS(here::here("results", "interaction_results.rds"))
dds <- readRDS(here::here("results", "dds_int_run.rds"))
cdf <- readRDS(here::here("results", "combined_df_annotated.rds"))
gmt <- fgsea::gmtPathways(here::here("data", "genesets_from_library",
                                     "mammary_mito_myc_metab_v1_mouse.gmt"))

CONTRASTS <- c("myc_6W_raw", "myc_12W_raw", "timepoint_neg_raw", "timepoint_pos_raw",
               "interaction_raw")
D  <- lapply(ir[CONTRASTS], as.data.frame)          # DESeqResults are S4; coerce once
nc <- DESeq2::counts(dds, normalized = TRUE)
L  <- log2(nc + 1)
sm <- as.data.frame(SummarizedExperiment::colData(dds))
sm$group <- factor(sm$group, levels = c("6W_neg", "6W_pos", "12W_neg", "12W_pos"))
sm$tp    <- stats::relevel(as.factor(sm$timepoint),  "6W")
sm$myc   <- stats::relevel(as.factor(sm$myc_status), "neg")

universe_all <- rownames(D$myc_6W_raw)
V   <- function(k) stats::setNames(D[[k]]$log2FoldChange, rownames(D[[k]]))
m6  <- V("myc_6W_raw");        m12 <- V("myc_12W_raw")
tn  <- V("timepoint_neg_raw"); tp  <- V("timepoint_pos_raw")
bm  <- stats::setNames(D$myc_6W_raw$baseMean, rownames(D$myc_6W_raw))

# set membership ALWAYS through the reconciler (CLAUDE.md / the 2026-07-24 fix)
ens_set <- function(syms) {
  e <- recon_to_ensembl(syms, universe_all)
  e[!is.na(e)]
}
set_mean <- function(v, e) if (length(e)) mean(v[e], na.rm = TRUE) else NA_real_

# single-gene lookup keeps script 42's symbol table (cheap, and it is only for genes)
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

epi_comp <- comp_of(c("Epcam", "Krt8", "Krt18", "Krt5", "Krt14", "Cdh1", "Krt17"))
imm_comp <- comp_of(c("Ptprc", "Cd52", "Cd3e", "Lyz2", "Cd74", "H2-Aa", "Itgam", "Ms4a1"))

# =============================================================================
# PART A: SUBSTRATE SPECIFICITY -- does the gland de-respire without de-proliferating?
# -----------------------------------------------------------------------------
# The comparator arms are chosen so that the alternative hypothesis has its best
# possible representatives: NUCLEOTIDE metabolism and the MITORIBOSOME are the two
# mitochondrial arms most tightly tied to growth, and PROLIF_* pooled is the
# proliferative programme itself. If the wild-type respiratory withdrawal were a
# growth withdrawal, those three should move with OXPHOS. The TEB signature is
# included because the competing (and compatible) reading is a change of cell
# IDENTITY rather than of proliferative rate.
# =============================================================================
arms <- tibble::tribble(
  ~arm,                    ~set,                                 ~ruler_pathway,
  "OXPHOS subunits",       "MITOCARTA_OXPHOS_SUBUNITS",          "OXPHOS subunits",
  "OXPHOS (all)",          "MITOCARTA_OXPHOS",                   "OXPHOS",
  "OXPHOS assembly",       "MITOCARTA_OXPHOS_ASSEMBLY_FACTORS",  "OXPHOS assembly factors",
  "nucleotide metabolism", "MITOCARTA_NUCLEOTIDE_METABOLISM",    "Nucleotide metabolism",
  "mitoribosome",          "MITOCARTA_MITOCHONDRIAL_RIBOSOME",   "Mitochondrial ribosome",
  "TCA cycle",             "MITOCARTA_TCA_CYCLE",                "TCA cycle",
  "amino-acid metabolism", "MITOCARTA_AMINO_ACID_METABOLISM",    "Amino acid metabolism",
  "lipid metabolism",      "MITOCARTA_LIPID_METABOLISM",         "Lipid metabolism",
  "TEB vs ductal (HS)",    "MG_TEB_VS_DUCTAL_HS_GRAY_UP",        NA_character_)

# THE SCRIPT-42 FAILURE MODE: a wrong set name gives a silent all-NaN composite.
missing_sets <- arms$set[!arms$set %in% names(gmt)]
if (length(missing_sets))
  stop("43 PART A: set names absent from the GMT -> ", paste(missing_sets, collapse = ", "))

# PROLIF_* is pooled: no single set is "proliferation", and the union is the fair
# representative of the alternative hypothesis.
prolif_sets <- grep("^PROLIF_", names(gmt), value = TRUE)
stopifnot(length(prolif_sets) >= 5)
arm_ens <- c(
  stats::setNames(lapply(arms$set, function(s) ens_set(gmt[[s]])), arms$arm),
  list("PROLIF_* pooled" = ens_set(unique(unlist(gmt[prolif_sets])))))

comparator <- dplyr::bind_rows(lapply(names(arm_ens), function(a) {
  e <- arm_ens[[a]]
  tibble::tibble(arm = a, n_genes = length(e),
                 c_myc_6W = set_mean(m6, e), c_myc_12W = set_mean(m12, e),
                 c_wt_time = set_mean(tn, e), c_myc_time = set_mean(tp, e))
})) |>
  dplyr::mutate(wt_pct_of_myc = round(100 * c_wt_time / c_myc_time)) |>
  dplyr::arrange(c_wt_time)

# --- the priority ruler alongside, where the arm is a MitoPathway ---------------
bvm_path <- here::here("results", "background_vs_myc.rds")
comparator_priority <- NULL
if (file.exists(bvm_path)) {
  rl <- tibble::as_tibble(readRDS(bvm_path)$ruler)
  comparator_priority <- arms |>
    dplyr::filter(!is.na(ruler_pathway)) |>
    dplyr::left_join(dplyr::select(rl, ruler_pathway = pathway,
                                   p_m6, p_m12, p_tn, p_tp), by = "ruler_pathway") |>
    dplyr::select(arm, ruler_pathway, prio_myc_6W = p_m6, prio_myc_12W = p_m12,
                  prio_wt_time = p_tn, prio_myc_time = p_tp)
  unmatched <- comparator_priority$arm[is.na(comparator_priority$prio_wt_time)]
  if (length(unmatched))
    message("43 PART A: no priority row for -> ", paste(unmatched, collapse = ", "))
} else {
  message("43 PART A: results/background_vs_myc.rds absent -- run script 40 first")
}

# --- matched-random-set null on the WILD-TYPE temporal contrast -----------------
# The claim is "the normal gland withdraws from the respiratory arm more than from
# an arbitrary set of comparably expressed genes". Matching is on baseMean, which
# is the variable that drives set-mean LFC precision.
expressed <- names(bm)[is.finite(bm) & bm > 0 & is.finite(tn[names(bm)])]
bin_of <- cut(rank(bm[expressed], ties.method = "first"),
              breaks = NBIN, labels = FALSE)
by_bin <- split(expressed, bin_of)
draw_matched <- function(e) {
  b <- bin_of[match(e, expressed)]
  b <- b[!is.na(b)]
  unlist(lapply(split(b, b), function(k)
    sample(by_bin[[as.character(k[1])]], length(k), replace = TRUE)), use.names = FALSE)
}
null_wt <- function(e, n = NSET)
  vapply(seq_len(n), function(i) set_mean(tn, draw_matched(e)), numeric(1))

wt_null <- dplyr::bind_rows(lapply(names(arm_ens), function(a) {
  e   <- arm_ens[[a]][arm_ens[[a]] %in% expressed]
  obs <- set_mean(tn, e)
  nul <- null_wt(e)
  tibble::tibble(arm = a, n_matched = length(e), observed_wt = obs,
                 null_median = stats::median(nul),
                 percentile  = 100 * mean(nul < obs),
                 p_emp_lower = mean(nul <= obs))
})) |> dplyr::arrange(percentile)

# --- the PAIRED null: OXPHOS-versus-proliferation is the actual claim -----------
# A one-set null answers "is OXPHOS extreme?". The claim is comparative, so the
# null has to be on the DIFFERENCE, with both sets redrawn together.
pair_contrasts <- list(
  c("OXPHOS subunits", "PROLIF_* pooled"),
  c("OXPHOS subunits", "nucleotide metabolism"),
  c("OXPHOS subunits", "mitoribosome"))
paired_null <- dplyr::bind_rows(lapply(pair_contrasts, function(p) {
  ea <- arm_ens[[p[1]]][arm_ens[[p[1]]] %in% expressed]
  eb <- arm_ens[[p[2]]][arm_ens[[p[2]]] %in% expressed]
  obs <- set_mean(tn, ea) - set_mean(tn, eb)
  nul <- vapply(seq_len(NSET), function(i)
    set_mean(tn, draw_matched(ea)) - set_mean(tn, draw_matched(eb)), numeric(1))
  tibble::tibble(arm_a = p[1], arm_b = p[2], observed_diff = obs,
                 null_median = stats::median(nul),
                 percentile  = 100 * mean(nul < obs),
                 p_emp_lower = mean(nul <= obs))
}))

# --- per-sample, within the 12 wild-type mice ----------------------------------
# DESCRIPTIVE ONLY: `batch = timepoint`, so a surviving `tp` term is not evidence
# of a developmental effect. What it can show is whether the respiratory fall is
# carried by the proliferation score -- i.e. collinearity, not causation.
wt_i <- which(sm$myc == "neg")
wt_within <- NULL
ox_lvl <- set_score("MITOCARTA_OXPHOS_SUBUNITS")
pr_lvl <- comp_of(c("Mki67", "Top2a", "Ccnb1", "Ccna2", "Bub1", "Plk1", "Aurka",
                    "Cdk1", "Pcna", "Mcm2", "Rrm2", "Tk1"))
if (!is.null(ox_lvl)) {
  dw <- data.frame(ox = ox_lvl[wt_i], pr = pr_lvl[wt_i], tp = droplevels(sm$tp[wt_i]),
                   epi = epi_comp[wt_i], imm = imm_comp[wt_i])
  f1 <- summary(stats::lm(ox ~ tp + epi + imm, dw))$coefficients
  f2 <- summary(stats::lm(ox ~ tp + pr + epi + imm, dw))$coefficients
  wt_within <- tibble::tibble(
    model     = c("ox ~ tp", "ox ~ tp + proliferation"),
    tp_beta   = c(f1["tp12W", 1], f2["tp12W", 1]),
    tp_p      = c(f1["tp12W", 4], f2["tp12W", 4]),
    prolif_beta = c(NA_real_, f2["pr", 1]),
    prolif_p    = c(NA_real_, f2["pr", 4]),
    cor_ox_prolif = stats::cor(dw$ox, dw$pr))
}

# --- the anti-apoptotic buffer over the timeline -------------------------------
# The gland's route through the trade-off is DE-PRIORITISATION, not buffering. If
# the buffer rose developmentally that would be a rival mechanism, so it is tested
# rather than assumed (death-narrative section 6 asserts it is flat).
pull <- function(genes) {
  e  <- ens_of(genes)
  ok <- !is.na(e) & e %in% rownames(D$myc_6W_raw)
  g  <- function(k, col) D[[k]][match(e[ok], rownames(D[[k]])), col]
  tibble::tibble(gene = genes[ok], baseMean = round(g("myc_6W_raw", "baseMean")),
                 lfc_myc_6W = g("myc_6W_raw", "log2FoldChange"),
                 padj_myc_6W = g("myc_6W_raw", "padj"),
                 lfc_wt_time = g("timepoint_neg_raw", "log2FoldChange"),
                 padj_wt_time = g("timepoint_neg_raw", "padj"),
                 lfc_myc_time = g("timepoint_pos_raw", "log2FoldChange"),
                 padj_myc_time = g("timepoint_pos_raw", "padj"))
}
buffer <- pull(c("Bcl2", "Bcl2l1", "Mcl1", "Bcl2l2", "Bcl2a1b", "Xiap", "Birc2",
                 "Birc3", "Birc5"))

# =============================================================================
# PART B: THE TRADE-OFF ASYMMETRY, IN mitoPPS RATIO SPACE
# -----------------------------------------------------------------------------
# mitoPPS is a PAIRWISE RATIO within the mitochondrial compartment, so the global
# common-mode axis that defeats every GSVA-space coupling (scripts 35/36) largely
# cancels: script 38 PART A measured the mitoPPS ambient at ~0.42 against GSVA's
# ~0.80. That ambient is recomputed here so the number is self-contained.
#
# ONE METHODOLOGICAL POINT THAT MATTERS. The outcomes are on different scales -- a
# log2 priming ratio and a z-mean proliferation composite are not comparable -- so
# every outcome is STANDARDISED before fitting. The `myc x axis` terms below are
# therefore comparable ACROSS outcomes, which is the whole point of the part, but
# they are NOT numerically comparable with script 42 PART H's +2.79, which was
# fitted on the unstandardised ratio. PART H's priming fit is repeated here on the
# standardised scale so the head-to-head is internal.
# =============================================================================
mp_path <- here::here("results", "mitopps_scores.rds")
axes <- list(); outcomes <- list(); ambient <- NULL
if (file.exists(mp_path)) {
  mp  <- readRDS(mp_path)
  mps <- mp$mitopps_scores
  mps <- mps[match(colnames(L), mps$sample), , drop = FALSE]
  stopifnot(identical(as.character(mps$sample), colnames(L)))
  need <- c("OXPHOS subunits", "ROS and glutathione metabolism",
            "Apoptosis-PRO", "Apoptosis-ANTI")
  missing_mps <- setdiff(need, names(mps))
  if (length(missing_mps)) {
    message("43 PART B: mitoPPS columns absent -> ", paste(missing_mps, collapse = ", "))
  } else {
    axes$oxphos_ppd <- as.numeric(mps[["OXPHOS subunits"]])
    axes$redox_ppd  <- as.numeric(mps[["ROS and glutathione metabolism"]])
    outcomes[["priming_ppd (PRO-ANTI)"]] <-
      as.numeric(mps[["Apoptosis-PRO"]]) - as.numeric(mps[["Apoptosis-ANTI"]])
    # the ambient: |rho| of each axis to every OTHER mitoPPS pathway score
    num <- vapply(mps, is.numeric, logical(1))
    M   <- as.matrix(mps[, num, drop = FALSE])
    ambient <- dplyr::bind_rows(lapply(names(axes), function(an) {
      keep <- setdiff(colnames(M), c("OXPHOS subunits", "ROS and glutathione metabolism"))
      rr <- abs(suppressWarnings(
        stats::cor(axes[[an]], M[, keep, drop = FALSE], method = "spearman")))
      tibble::tibble(axis = an, n_pathways = length(keep),
                     ambient_median_abs_rho = stats::median(rr, na.rm = TRUE),
                     ambient_q90 = stats::quantile(rr, 0.90, na.rm = TRUE, names = FALSE))
    }))
  }
} else {
  message("43 PART B: results/mitopps_scores.rds absent -- run script 08 first")
}

ratio_of <- function(a, b) {
  ea <- ens_of(a); eb <- ens_of(b)
  stopifnot(!is.na(ea), !is.na(eb), ea %in% rownames(L), eb %in% rownames(L))
  as.numeric(L[ea, ] - L[eb, ])
}
outcomes[["Bbc3:Bcl2l1 (PUMA priming)"]] <- ratio_of("Bbc3", "Bcl2l1")
outcomes[["Bax:Bcl2l1"]]                 <- ratio_of("Bax",  "Bcl2l1")
outcomes[["proliferation (markers)"]]    <- pr_lvl
pp <- set_score("PROLIF_E2F_HALLMARK")
if (!is.null(pp)) outcomes[["proliferation (E2F hallmark)"]] <- pp

Z <- function(x) as.numeric(scale(x))
resid_design <- function(y)
  stats::residuals(stats::lm(y ~ tp * myc + epi + imm,
                             data.frame(y = y, tp = sm$tp, myc = sm$myc,
                                        epi = epi_comp, imm = imm_comp)))

tradeoff <- NULL; tradeoff_perm <- NULL
if (length(axes)) {
  tradeoff <- dplyr::bind_rows(lapply(names(outcomes), function(on) {
    y <- Z(outcomes[[on]]); ry <- resid_design(outcomes[[on]])
    dplyr::bind_rows(lapply(names(axes), function(an) {
      a  <- axes[[an]]
      d  <- data.frame(y = y, tp = sm$tp, myc = sm$myc, epi = epi_comp,
                       imm = imm_comp, a = a)
      m1 <- summary(stats::lm(y ~ myc * a + epi + imm, d))$coefficients
      m2 <- summary(stats::lm(y ~ tp * myc + myc:a + a + epi + imm, d))$coefficients
      tibble::tibble(
        outcome = on, axis = an,
        rho_raw   = stats::cor(a, outcomes[[on]], method = "spearman"),
        rho_adj   = stats::cor(resid_design(a), ry, method = "spearman"),
        myc_x_axis = m1["mycpos:a", 1], p = m1["mycpos:a", 4],
        myc_x_axis_with_tp = m2["mycpos:a", 1], p_with_tp = m2["mycpos:a", 4])
    }))
  }))

  tradeoff_perm <- dplyr::bind_rows(lapply(names(outcomes), function(on) {
    y <- Z(outcomes[[on]])
    dplyr::bind_rows(lapply(names(axes), function(an) {
      d   <- data.frame(y = y, myc = sm$myc, epi = epi_comp, imm = imm_comp,
                        a = axes[[an]])
      obs <- summary(stats::lm(y ~ myc * a + epi + imm, d))$coefficients["mycpos:a", 1]
      nul <- vapply(seq_len(NPERM), function(i) {
        dd <- d
        for (tv in levels(sm$tp)) { k <- which(sm$tp == tv); dd$a[k] <- sample(dd$a[k]) }
        summary(stats::lm(y ~ myc * a + epi + imm, dd))$coefficients["mycpos:a", 1]
      }, numeric(1))
      tibble::tibble(outcome = on, axis = an, observed = obs,
                     null_median = stats::median(nul),
                     percentile = 100 * mean(nul < obs), p_emp = mean(nul >= obs))
    }))
  }))
}

# =============================================================================
# ASSERTS
# =============================================================================
# PART A's built-in positive control: the comparator must reproduce Issue #4's
# wild-type OXPHOS-subunit value. If it does not, the membership resolution has
# drifted and nothing downstream is comparable with the committed numbers.
OX_WT_REF <- -0.2548
ox_wt <- comparator$c_wt_time[comparator$arm == "OXPHOS subunits"]
stopifnot(length(ox_wt) == 1L)
if (!is.finite(ox_wt) || abs(ox_wt - OX_WT_REF) > 0.02)
  warning(sprintf(paste("43 PART A: OXPHOS-subunit wild-type temporal %.4f does not match",
                        "attenuation_decomposition.rds (%.4f) -- membership has drifted."),
                  ox_wt, OX_WT_REF))
stopifnot(
  nrow(comparator) == length(arm_ens),
  all(comparator$n_genes >= 15),                       # no silently empty composite
  all(is.finite(wt_null$p_emp_lower)),
  nrow(paired_null) == length(pair_contrasts))
if (!is.null(tradeoff)) {
  stopifnot(nrow(tradeoff) == length(outcomes) * length(axes),
            all(is.finite(tradeoff_perm$p_emp)))
  if (!"redox_ppd" %in% names(axes))
    warning("43 PART B: redox_ppd control axis missing -- the positives cannot be read")
} else {
  warning("43 PART B: not run (mitoPPS scores absent)")
}

# =============================================================================
# SAVE
# =============================================================================
notes <- c(
  "THE TRADE-OFF MODEL. Mitochondrial respiratory capacity is simultaneously the",
  "strongest pro-growth and the strongest pro-death input to a MYC-driven cell, and",
  "which limb is expressed depends on whether death can be executed. The cell arm",
  "shows the sign reversal outright: PGC1a KILLS MYAZ cells; PGC1a + Bcl-xL grows",
  "BETTER than Bcl-xL alone. Three solutions to one trade-off appear in the paper --",
  "development gives the escape away free (the window), selection finds it expensively",
  "(SS and passaged cells lose the asset with the liability), and anti-apoptotic",
  "buffering solves it properly (respiration becomes pure asset). The closing claim:",
  "what changes across progression is not what the tumour NEEDS from its mitochondria",
  "but what it can AFFORD.",
  "",
  "PART A -- the alternative hypothesis is that the adult gland's respiratory",
  "withdrawal is just its exit from the pubertal/TEB proliferative state. The",
  "comparator is built to give that hypothesis its best representatives: NUCLEOTIDE",
  "metabolism and the MITORIBOSOME are the two mitochondrial arms most tied to growth,",
  "and PROLIF_* pooled is the proliferative programme itself. The matched null makes",
  "'OXPHOS withdraws more than expected' a percentile; the PAIRED null is the one that",
  "matches the claim, because the claim is comparative.",
  "",
  "PART A IS THE LOAD-BEARING ONE FOR THE MODEL. If respiration fell only because",
  "proliferation fell, the two-input model would collapse into 'less proliferation,",
  "less death' and the MYC-ER asymmetry (~80% less death vs ~30% less proliferation)",
  "would have no explanation. A dissociation is what the model needs.",
  "",
  "PART A's per-sample wild-type fit is DESCRIPTIVE ONLY. batch = timepoint, so a",
  "surviving tp term is not evidence of a developmental effect; what it can show is",
  "whether the respiratory fall is collinear with the proliferation score.",
  "",
  "PART B -- the outcomes are STANDARDISED before fitting, so the myc x axis terms",
  "are comparable ACROSS outcomes. They are NOT comparable with script 42 PART H's",
  "+2.79, which was fitted on the unstandardised ratio; the priming fit is repeated",
  "here on the standardised scale so the head-to-head is internal. redox_ppd is the",
  "NEGATIVE CONTROL and must be null for the positives to be readable (it was null",
  "throughout in script 42 PART H: p 0.72, permutation 51st percentile).",
  "",
  "SCOPE. Nothing here separates from timepoint at n=24 -- script 42 PART H already",
  "showed the best term dropping p 0.005 -> 0.088 once tp*myc is in the model. Both",
  "parts are RANKING. The epistemic contract stands: the in-vivo transcriptome",
  "GENERATES the hypothesis, the cell perturbations PROVE it. And the death phenotype",
  "is external (IHC, MYC-ER counts) -- priming is not death.")

out <- list(
  comparator          = comparator,
  comparator_priority = comparator_priority,
  wt_null             = wt_null,
  paired_null         = paired_null,
  wt_within           = wt_within,
  buffer              = buffer,
  ambient             = ambient,
  tradeoff            = tradeoff,
  tradeoff_perm       = tradeoff_perm,
  defs = list(arms = arms, prolif_sets = prolif_sets,
              n_set_draws = NSET, n_perm = NPERM, n_bins = NBIN,
              oxphos_wt_reference = OX_WT_REF),
  analysis_date = Sys.Date(),
  notes = notes)

saveRDS(out, here::here("results", "substrate_specificity_tradeoff.rds"))
message("43: wrote results/substrate_specificity_tradeoff.rds")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  res <- readRDS(here::here("results", "substrate_specificity_tradeoff.rds"))

  ## PART A -- the headline. Read c_wt_time down the column: OXPHOS subunits should
  ## sit at the top (most negative) and the growth-coupled arms near zero. If
  ## nucleotide metabolism or the mitoribosome tracked OXPHOS, the proliferation
  ## explanation would be alive and the two-input model would be in trouble.
  res$comparator |> print(n = 12)

  ## The same arms on the dose-cancelling priority ruler.
  res$comparator_priority |> print(n = 12)

  ## Is the wild-type respiratory withdrawal extreme against matched random genes?
  ## Low percentile = extreme in the negative direction.
  res$wt_null |> print(n = 12)

  ## THE ACTUAL CLAIM, which is comparative: OXPHOS versus each growth-coupled arm,
  ## with both sets redrawn together.
  res$paired_null |> print()

  ## Descriptive only (batch = timepoint). Does the proliferation score absorb the
  ## respiratory fall within the 12 wild-type mice? cor_ox_prolif says how collinear
  ## the question even is.
  res$wt_within |> print()

  ## The rival mechanism: does the gland buffer instead of de-prioritising? Every
  ## padj_wt_time should be non-significant (death narrative section 6).
  res$buffer |> print(n = 10)

  ## PART B -- the ambient first. mitoPPS should be ~0.42, not GSVA's ~0.80; the
  ## couplings below are only readable against this.
  res$ambient |> print()

  ## The asymmetry. Compare myc_x_axis for the priming outcomes against the
  ## proliferation outcomes ON THE SAME AXIS (oxphos_ppd). The trade-off model
  ## predicts respiratory priority tracks death competence more tightly than
  ## proliferation. Read every redox_ppd row as the control: if redox is not null,
  ## nothing in the oxphos_ppd rows can be read.
  res$tradeoff |> dplyr::arrange(axis, dplyr::desc(abs(myc_x_axis))) |> print(n = 20)

  ## And the within-timepoint permutation, which is the honest null for the above.
  res$tradeoff_perm |> dplyr::arrange(axis, dplyr::desc(percentile)) |> print(n = 20)

  ## Sanity: the positive control. Should be within 0.02 of -0.2548.
  res$comparator$c_wt_time[res$comparator$arm == "OXPHOS subunits"]

  cat(res$notes, sep = "\n")
}
