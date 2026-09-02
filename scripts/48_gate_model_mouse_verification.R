# =============================================================================
# 48_gate_model_mouse_verification.R
# -----------------------------------------------------------------------------
# THE GATE MODEL: a quantitative statement of MYC x OXPHOS -> apoptotic priming
# that can be carried to human tumour cohorts, and its verification in the mouse.
#
# WHY THIS SCRIPT EXISTS. The human arm (docs/2026-08-27_human_validation_plan.md)
# proposes to test `PRIME ~ MYC * OXPHOS` in TCGA / METABRIC / SCAN-B. That is the
# right shape, but three of its measurement choices are assumptions, not results,
# and all three are decidable HERE, in the mouse, before a single human sample is
# downloaded:
#   (1) the human M is a CONTINUOUS MYC activity score; the mouse M that produced
#       the published interaction is a BINARY genotype. Does the interaction
#       survive the substitution?
#   (2) the human plan makes the OXPHOS LEVEL primary and mitoPPS secondary. The
#       published mouse interaction is on mitoPPS. Which ruler transfers?
#   (3) the endpoint is `log2(BBC3) - log2(BCL2L1)`. Which HALF of that ratio
#       carries the interaction -- the sensor or the guardian?
#
# THE MODEL, stated before it is fitted. Write M for MYC dose, X for the
# respiratory state of the mitochondrial compartment, T for the apoptotic trigger
# read as a BH3-sensor : guardian log-ratio. The claim is not that M raises T, and
# not that X raises T. It is that the two are INTEGRATED -- T responds to X only
# in proportion to M:
#
#       T  =  b0  +  bM * M  +  bX * X  +  bMX * (M x X)  +  g' C          (1)
#
#   H-gate    bMX > 0                    the AND-gate; the whole claim
#   H-null    bX  = 0 at M = 0           respiration alone does not prime
#   H-spec    bMX = 0 when X is any      the gate is respiratory, not
#             other mitochondrial arm    mitochondrial-in-general
#   H-endp    bMX = 0 for sensors that   the gate is not a general apoptotic
#             the model does not name    shift
#
# If H-null holds, (1) collapses to a ONE-PARAMETER model with no independent
# OXPHOS term -- T = b0 + bM*M + bMX*(M x X) -- which is the manuscript's title
# thesis written as an equation, and a sharper thing to test than (1).
#
# The gate has a scale-free summary that transfers across species and platforms:
#
#       M*  =  -bX / bMX                                                   (2)
#
# the MYC level, in SD units, at which dT/dX changes sign -- below it respiration
# is neutral-to-protective, above it respiration primes. M* is the "asset becomes
# liability" point, and it is a number, not a metaphor.
#
# WHAT THE MOUSE CAN AND CANNOT VERIFY. It can verify (1) and (2): both inputs
# vary, the endpoint is measured, and the design is not the tumour. It CANNOT
# verify the two layers that only exist in tumours -- that the lethal corner is
# under-occupied because its occupants died (survivor bias), and that the survivors
# bought their way in by buffering. Those are human-only by construction, which is
# the point of the human arm.
#
# PART A -- THE MEASUREMENT BRIDGE. Three OXPHOS rulers (absolute level;
#   compartment-relative level; mitoPPS priority), four MYC estimators (genotype;
#   MYC transcript; MSigDB target signature; DoRothEA regulon -- the last two
#   MitoCarta-stripped, as the human plan requires), and the separability table
#   that says how far each ruler is from the other input. An interaction test needs
#   its two inputs to be distinguishable; this part measures whether they are.
#
# PART B -- THE GATE, GENOTYPE FORM. Equation (1) for every axis, with the
#   within-timepoint permutation null of scripts 42/43 (shuffling X inside each
#   timepoint preserves the design and breaks only the mouse-to-mouse link, which
#   is the thing being claimed).
#
# PART C -- THE GATE, CONTINUOUS-MYC FORM. The transfer test. Each estimator's
#   coefficient is rescaled by its own genotype gap so the four are comparable:
#   an estimator that recovers half the genotype effect is a half-power instrument
#   for the human study, and the human study must be powered on ITS number, not on
#   the genotype number.
#
# PART D -- MODEL FORM. Additive vs interaction vs gate-only (drop the X main
#   effect); the crossover M*; the simple slopes dT/dX within each genotype.
#
# PART E -- WHICH HALF OF THE RATIO. Numerator and denominator fitted separately,
#   each against an expression-matched null and against all expressed genes. A
#   ratio endpoint whose signal sits in its denominator is not a statement about
#   the numerator, and the human pre-specification depends on knowing which it is.
#
# PART F -- THE BUFFERING ARM. MCL1 vs BCL2L1 under the same gate. If the mouse
#   already switches guardian dependence as M x X rises, the human buffering
#   hypothesis (H1) has a mouse-side prior and a drug-shaped prediction.
#
# PART G -- WITHIN-TIMEPOINT ONLY. Every variable centred inside its timepoint, so
#   the between-cohort contrast -- which is also the batch contrast -- is gone
#   entirely. This is the only regime a cross-sectional human cohort has, so it is
#   the honest transfer estimate even though it is the weakest.
#
# PART H -- POWER TRANSFER. What n a human cohort needs to see the mouse effect,
#   and a third of it, on each estimator.
#
# PART I -- VERDICT, on a rule fixed before the tests (script 46/47 pattern).
#
# SCOPE, said once. BATCH = TIMEPOINT (CLAUDE.md), so PART B's full-cohort fits
# carry the cohort contrast; PART G removes it and is the number to quote when the
# claim is cross-sectional. n = 24: this is RANKING, a set of negatives, and a
# measurement-design decision for the human arm -- not confirmatory inference. And
# every model here is a CORRELATION among 24 mice; the causal claims in this paper
# belong to the perturbations, not to this script.
#
# Reads : results/dds_int_run.rds, results/combined_df_annotated.rds,
#         results/mitopps_scores.rds            (script 08),
#         data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt,
#         functions/reconcile_gene_symbols.R (MANDATORY -- vintage-aware membership)
# Writes: results/gate_model_verification.rds
#
# SELF-CONTAINED: those three .rds are the only state it needs, so a cold R session
# is the right way to run it. The positive control is REFITTED here from dds +
# mitopps_scores rather than read back from script 42's or 43's saved object -- that
# is the point of it, since reading the answer would not test the load block.
#
# RUNTIME: about two minutes. PART B runs 5000 within-timepoint permutations for
# each of 8 axes; PART E fits ~17k gene-wise models by a single lm.fit.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

set.seed(1)
NPERM <- 5000L    # within-timepoint permutations (scripts 42/43 idiom)
NBIN  <- 20L      # baseMean bins for the expression-matched null
MINBM <- 10       # smallest mean normalised count worth a gene-wise fit

# POSITIVE CONTROL -- script 43 PART B, row `Bbc3:Bcl2l1 x oxphos_ppd`, refitted
# here from raw objects. Same numbers or the load block is wrong.
GATE_PPD_REF <- 6.089      # myc x oxphos_ppd on the standardised PUMA ratio
GATE_PPD_P   <- 0.00523
CTRL_TOL     <- 5e-3

# VERDICT RULE, fixed before the tests.
VERDICT_RULE <- c(
  transfers      = "same sign, p < 0.05 in PART B, and >= 50 pct of the genotype effect recovered in PART C",
  transfers_weak = "same sign and p < 0.05 in PART B, but < 50 pct recovered in PART C",
  fails          = "sign reverses, or PART B p >= 0.05")

# =============================================================================
# PART 0: LOAD, ALIGN, AND PROVE THE MACHINERY REPRODUCES SCRIPT 43
# =============================================================================
message("48 PART 0: load")

dds <- readRDS(here::here("results", "dds_int_run.rds"))
cdf <- readRDS(here::here("results", "combined_df_annotated.rds"))
gmt <- fgsea::gmtPathways(here::here("data", "genesets_from_library",
                                     "mammary_mito_myc_metab_v1_mouse.gmt"))

NC <- DESeq2::counts(dds, normalized = TRUE)
L  <- log2(NC + 1)
sm <- as.data.frame(SummarizedExperiment::colData(dds))
sm$group <- factor(sm$group, levels = c("6W_neg", "6W_pos", "12W_neg", "12W_pos"))
sm$tp    <- stats::relevel(as.factor(sm$timepoint),  "6W")
sm$myc   <- stats::relevel(as.factor(sm$myc_status), "neg")
universe_all <- rownames(L)

sym2ens <- cdf[!is.na(cdf$mgi_symbol) & !duplicated(cdf$mgi_symbol),
               c("mgi_symbol", "gene")]
ens_of  <- function(s) sym2ens$gene[match(s, sym2ens$mgi_symbol)]
zrow    <- function(m) t(scale(t(m)))
comp_e  <- function(e) {
  e <- e[!is.na(e) & e %in% rownames(L)]
  stopifnot(length(e) >= 3)
  colMeans(zrow(L[e, , drop = FALSE]))
}
comp_of <- function(syms) comp_e(ens_of(syms))
# set membership ALWAYS through the reconciler (CLAUDE.md / the 2026-07-24 fix)
ens_set <- function(syms) { e <- recon_to_ensembl(syms, universe_all); e[!is.na(e)] }
set_e   <- function(s) {
  stopifnot(!is.null(gmt[[s]]))
  e <- ens_set(gmt[[s]]); e[e %in% rownames(L)]
}
Z     <- function(x) as.numeric(scale(x))
ratio <- function(a, b) as.numeric(L[ens_of(a), ] - L[ens_of(b), ])
ixof  <- function(s) grep(":a$", rownames(s), value = TRUE)[1]   # the M x X row
# centre inside each timepoint: removes the cohort/batch contrast entirely
ctr <- function(v) {
  o <- v
  for (tv in levels(sm$tp)) { k <- sm$tp == tv; o[k] <- v[k] - mean(v[k]) }
  o
}

# covariates, script 42/43 idiom verbatim
epi <- comp_of(c("Epcam", "Krt8", "Krt18", "Krt5", "Krt14", "Cdh1", "Krt17"))
imm <- comp_of(c("Ptprc", "Cd52", "Cd3e", "Lyz2", "Cd74", "H2-Aa", "Itgam", "Ms4a1"))
pro <- comp_of(c("Mki67", "Top2a", "Ccnb1", "Plk1", "Aurka", "Bub1", "Ccna2", "Rrm2"))

mp <- readRDS(here::here("results", "mitopps_scores.rds"))$mitopps_scores
mp <- mp[match(colnames(L), mp$sample), , drop = FALSE]
stopifnot(identical(as.character(mp$sample), colnames(L)))

# --- POSITIVE CONTROL ---------------------------------------------------------
d0 <- data.frame(y = Z(ratio("Bbc3", "Bcl2l1")), m = sm$myc,
                 a = as.numeric(mp[["OXPHOS subunits"]]), epi = epi, imm = imm)
s0 <- summary(stats::lm(y ~ m * a + epi + imm, d0))$coefficients
ctrl <- c(observed = s0[ixof(s0), 1], reference = GATE_PPD_REF,
          p_observed = s0[ixof(s0), 4], p_reference = GATE_PPD_P)
message(sprintf("48 PART 0: positive control %+.3f (ref %+.3f), p %.5f (ref %.5f)",
                ctrl[1], ctrl[2], ctrl[3], ctrl[4]))
stopifnot(abs(ctrl[["observed"]] - GATE_PPD_REF) < CTRL_TOL,
          abs(ctrl[["p_observed"]] - GATE_PPD_P) < CTRL_TOL)

# =============================================================================
# PART A: THE MEASUREMENT BRIDGE
# -----------------------------------------------------------------------------
# Three rulers for X. `ox_lvl` is the human plan's primary: the absolute level of
# the nuclear OXPHOS subunits. `ox_ppd` is mitoPPS, the ruler the published mouse
# interaction was fitted on -- pairwise-ratio, so the compartment-wide common mode
# cancels, but composition-dependent and NEVER comparable across cohorts. `ox_rel`
# is the bridge this script proposes: the respiratory arm's level MINUS the level
# of the rest of the mitochondrial compartment. It is a within-compartment share,
# like mitoPPS, but it is a difference of two z-composites -- so it needs nothing
# but MitoCarta membership and an expression matrix, and it is computed the same
# way in any cohort of any species.
#
# Four estimators for M, ordered from the experimental fact to the most derived.
# The two signature scores are MitoCarta-stripped, as the human plan requires:
# an unstripped MYC signature shares genes with X and the interaction would be
# partly a set overlap.
# =============================================================================
message("48 PART A: rulers and estimators")

mito_sets <- grep("^MITOCARTA_", names(gmt), value = TRUE)
mito_all  <- unique(unlist(lapply(mito_sets, set_e)))
ox_sub    <- set_e("MITOCARTA_OXPHOS_SUBUNITS")       # nuclear only by construction
rest_mito <- setdiff(mito_all, ox_sub)

relify <- function(s) { e <- set_e(s); comp_e(e) - comp_e(setdiff(mito_all, e)) }

X <- list(
  ox_rel      = comp_e(ox_sub) - comp_e(rest_mito),   # proposed bridge ruler
  ox_ppd      = as.numeric(mp[["OXPHOS subunits"]]),  # published mouse ruler
  ox_lvl      = comp_e(ox_sub),                       # human plan's primary
  oxasm_rel   = relify("MITOCARTA_OXPHOS_ASSEMBLY_FACTORS"),
  mitorib_rel = relify("MITOCARTA_MITOCHONDRIAL_RIBOSOME"),
  tca_rel     = relify("MITOCARTA_TCA_CYCLE"),
  fao_rel     = relify("MITOCARTA_FATTY_ACID_OXIDATION"),
  redox_rel   = relify("MITOCARTA_ROS_AND_GLUTATHIONE_METABOLISM"))

strip <- function(e) setdiff(e, mito_all)
M <- list(
  geno = sm$myc,
  mRNA = Z(as.numeric(L[ens_of("Myc"), ])),
  msig = Z(comp_e(strip(set_e("TFT_MYC_MSIGDB")))),
  dor  = Z(comp_e(strip(set_e("TFT_MYC_DOROTHEA_ABC")))))

Y <- list(PUMA  = ratio("Bbc3", "Bcl2l1"),    BIM   = ratio("Bcl2l11", "Bcl2l1"),
          BAX   = ratio("Bax",  "Bcl2l1"),    BID   = ratio("Bid",  "Bcl2l1"),
          BAK1  = ratio("Bak1", "Bcl2l1"),    NOXA  = ratio("Pmaip1", "Bcl2l1"),
          BMF   = ratio("Bmf",  "Bcl2l1"),    BUFFER = ratio("Mcl1", "Bcl2l1"),
          PUMA_over_MCL1 = ratio("Bbc3", "Mcl1"))

# an interaction test needs its two inputs separable; this is that measurement
separability <- dplyr::bind_rows(lapply(names(X), function(an) tibble::tibble(
  axis = an,
  r_myc_mRNA = stats::cor(X[[an]], M$mRNA),
  r_myc_sig  = stats::cor(X[[an]], M$msig),
  r_myc_dor  = stats::cor(X[[an]], M$dor),
  r_prolif   = stats::cor(X[[an]], pro),
  r_epi      = stats::cor(X[[an]], epi))))

# the arms are not independent of each other either; the head-to-head in PART B
# only means something next to how correlated the two arms it separates are
axis_cor <- stats::cor(do.call(cbind, X))

estimator_gap <- dplyr::bind_rows(lapply(c("mRNA", "msig", "dor"), function(n)
  tibble::tibble(estimator = n,
                 genotype_gap_sd = mean(M[[n]][sm$myc == "pos"]) -
                                   mean(M[[n]][sm$myc == "neg"]),
                 r_to_myc_mRNA = stats::cor(M[[n]], M$mRNA))))

set_sizes <- tibble::tibble(
  n_mitocarta = length(mito_all), n_oxphos_subunits = length(ox_sub),
  n_rest_compartment = length(rest_mito),
  n_myc_msigdb_stripped = length(strip(set_e("TFT_MYC_MSIGDB"))),
  n_myc_dorothea_stripped = length(strip(set_e("TFT_MYC_DOROTHEA_ABC"))))

# =============================================================================
# PART B: THE GATE, GENOTYPE FORM, WITH A WITHIN-TIMEPOINT PERMUTATION NULL
# -----------------------------------------------------------------------------
# The null shuffles X inside each timepoint. It therefore keeps the design, the
# genotype split and the cohort contrast, and destroys only the animal-to-animal
# pairing of X with the endpoint -- which is exactly what the gate asserts. Its
# median is NOT zero (scripts 33/35: at this n everything correlates with
# everything), and the median is the honest reference, not zero.
# =============================================================================
message("48 PART B: the gate, genotype form")

gate_fit <- function(y, m, a, extra = NULL) {
  d <- data.frame(y = Z(y), m = m, a = Z(a), epi = epi, imm = imm,
                  tp = sm$tp, pro = pro)
  f <- "y ~ m * a + epi + imm"
  if (!is.null(extra)) f <- paste(f, "+", extra)
  s <- summary(stats::lm(stats::as.formula(f), d))$coefficients
  list(bM = s[grep("^m", rownames(s))[1], 1], bX = s["a", 1],
       bMX = s[ixof(s), 1], p = s[ixof(s), 4], se = s[ixof(s), 2])
}

perm_gate <- function(y, a) {
  d <- data.frame(y = Z(y), m = sm$myc, a = Z(a), epi = epi, imm = imm)
  o <- summary(stats::lm(y ~ m * a + epi + imm, d))$coefficients
  obs <- o[ixof(o), 1]
  nul <- vapply(seq_len(NPERM), function(i) {
    dd <- d
    for (tv in levels(sm$tp)) { k <- which(sm$tp == tv); dd$a[k] <- sample(dd$a[k]) }
    ss <- summary(stats::lm(y ~ m * a + epi + imm, dd))$coefficients
    ss[ixof(ss), 1]
  }, numeric(1))
  c(observed = obs, null_median = stats::median(nul),
    percentile = 100 * mean(nul < obs), p_emp = mean(nul >= obs))
}

gate_axes <- dplyr::bind_rows(lapply(names(X), function(an) {
  f0 <- gate_fit(Y$PUMA, sm$myc, X[[an]])
  f1 <- gate_fit(Y$PUMA, sm$myc, X[[an]], "tp")
  f2 <- gate_fit(Y$PUMA, sm$myc, X[[an]], "tp + pro")
  pp <- perm_gate(Y$PUMA, X[[an]])
  tibble::tibble(axis = an, bX = f0$bX, bMX = f0$bMX, p = f0$p,
                 bMX_tp = f1$bMX, p_tp = f1$p,
                 bMX_tp_prolif = f2$bMX, p_tp_prolif = f2$p,
                 perm_null_median = pp[["null_median"]],
                 perm_pct = pp[["percentile"]], perm_p = pp[["p_emp"]])
}))

# the mitoribosome is the one control arm that can rival OXPHOS here; it is also
# 0.75-correlated with it, so the head-to-head is the test that separates them
d_hh <- data.frame(y = Z(Y$PUMA), m = sm$myc, ox = Z(X$ox_rel),
                   mr = Z(X$mitorib_rel), epi = epi, imm = imm)
head_to_head <- as.data.frame(summary(
  stats::lm(y ~ m * ox + m * mr + epi + imm, d_hh))$coefficients)
head_to_head$term <- rownames(head_to_head)

# =============================================================================
# PART C: THE GATE, CONTINUOUS-MYC FORM -- THE TRANSFER TEST
# -----------------------------------------------------------------------------
# The human cohorts have no genotype. Every estimator's coefficient is per SD of
# that estimator, so the four are not comparable as printed; multiplying by the
# estimator's own genotype gap puts them all on the genotype scale and answers the
# only question that matters -- how much of the experimentally-created effect does
# this instrument see?
# =============================================================================
message("48 PART C: the transfer test")

geno_ref <- gate_fit(Y$PUMA, sm$myc, X$ox_rel)$bMX

gate_estimators <- dplyr::bind_rows(lapply(c("ox_rel", "ox_ppd", "ox_lvl"), function(an)
  dplyr::bind_rows(lapply(c("mRNA", "msig", "dor"), function(mn) {
    f   <- gate_fit(Y$PUMA, M[[mn]], X[[an]])
    gap <- mean(M[[mn]][sm$myc == "pos"]) - mean(M[[mn]][sm$myc == "neg"])
    tibble::tibble(axis = an, estimator = mn,
                   bMX_per_sd = f$bMX, p = f$p, genotype_gap_sd = gap,
                   bMX_in_genotype_units = f$bMX * gap,
                   frac_of_genotype = f$bMX * gap /
                     gate_fit(Y$PUMA, sm$myc, X[[an]])$bMX)
  }))))

# =============================================================================
# PART D: MODEL FORM, AND THE CROSSOVER
# -----------------------------------------------------------------------------
# H-null says X has no effect of its own. With a binary M that is testable by
# dropping the X main effect: if the gate-only model fits as well, equation (1)
# has one fewer parameter and reads as a pure product.
# =============================================================================
message("48 PART D: model form")

# BOTH co-primary endpoints: the collapse to (1') is a property of the ENDPOINT,
# not of the model. PUMA drops its X main effect for free; BUFFER does not, and a
# single row would have hidden that.
model_form <- dplyr::bind_rows(lapply(c("PUMA", "BUFFER"), function(ep)
                dplyr::bind_rows(lapply(c("geno", "mRNA", "msig"), function(mn) {
  mv <- if (mn == "geno") ifelse(sm$myc == "pos", 1, 0) else M[[mn]]
  d  <- data.frame(y = Z(Y[[ep]]), m = mv, a = Z(X$ox_rel), epi = epi, imm = imm)
  m_add  <- stats::lm(y ~ m + a + epi + imm, d)
  m_full <- stats::lm(y ~ m * a + epi + imm, d)
  m_gate <- stats::lm(y ~ m + I(m * a) + epi + imm, d)
  s <- summary(m_full)$coefficients
  tibble::tibble(endpoint = ep, estimator = mn,
                 bX = s["a", 1], p_bX = s["a", 4],
                 r2_additive = summary(m_add)$r.squared,
                 r2_interaction = summary(m_full)$r.squared,
                 r2_gate_only = summary(m_gate)$r.squared,
                 aic_additive = stats::AIC(m_add), aic_interaction = stats::AIC(m_full),
                 p_interaction_vs_additive = stats::anova(m_add, m_full)$`Pr(>F)`[2],
                 p_drop_X_main = stats::anova(m_gate, m_full)$`Pr(>F)`[2],
                 crossover_Mstar = -s["a", 1] / s[ixof(s), 1])
}))))

simple_slopes <- dplyr::bind_rows(lapply(c("ox_rel", "ox_ppd"), function(an)
  dplyr::bind_rows(lapply(c("none", "epi + imm"), function(cv) {
    d <- data.frame(y = Y$PUMA, a = X[[an]], epi = epi, imm = imm, m = sm$myc)
    f <- if (cv == "none") "y ~ a" else "y ~ a + epi + imm"
    w <- summary(stats::lm(stats::as.formula(f), d[d$m == "neg", ]))$coefficients
    p <- summary(stats::lm(stats::as.formula(f), d[d$m == "pos", ]))$coefficients
    tibble::tibble(axis = an, covariates = cv,
                   slope_wt = w["a", 1], p_wt = w["a", 4],
                   r2_wt = summary(stats::lm(stats::as.formula(f),
                                             d[d$m == "neg", ]))$r.squared,
                   slope_myc = p["a", 1], p_myc = p["a", 4],
                   r2_myc = summary(stats::lm(stats::as.formula(f),
                                              d[d$m == "pos", ]))$r.squared)
  }))))

# =============================================================================
# PART E: WHICH HALF OF THE RATIO CARRIES THE GATE
# -----------------------------------------------------------------------------
# Every endpoint in the menu shares the BCL2L1 denominator. If BCL2L1 itself moves
# under the gate, then several "specific" endpoints are one moving denominator
# seen through different numerators, and the pre-specification has to say so.
# The genome-wide fit gives the only null that means anything at this n: where
# does a named gene sit among all expressed genes, and among genes of its own
# expression level?
# =============================================================================
message("48 PART E: numerator, denominator, and the genome-wide calibration")

single_gene <- function(sym) {
  d <- data.frame(y = Z(as.numeric(L[ens_of(sym), ])), m = sm$myc,
                  a = Z(X$ox_rel), epi = epi, imm = imm)
  s <- summary(stats::lm(y ~ m * a + epi + imm, d))$coefficients
  dw <- data.frame(y = ctr(as.numeric(L[ens_of(sym), ])), m = sm$myc,
                   a = Z(ctr(X$ox_rel)), epi = ctr(epi), imm = ctr(imm))
  sw <- summary(stats::lm(y ~ m * a + epi + imm, dw))$coefficients
  tibble::tibble(gene = sym, bX_wt = s["a", 1], bX_myc = s["a", 1] + s[ixof(s), 1],
                 bMX = s[ixof(s), 1], p = s[ixof(s), 4],
                 bMX_within_tp = sw[ixof(sw), 1], p_within_tp = sw[ixof(sw), 4])
}
ROSTER <- c("Bbc3", "Bcl2l1", "Mcl1", "Bcl2l11", "Bax", "Bak1", "Bid", "Bmf",
            "Pmaip1", "Bcl2", "Bcl2l2", "Xiap", "Birc5", "Myc")
gene_gate <- dplyr::bind_rows(lapply(ROSTER, single_gene))

keep <- rowMeans(NC) >= MINBM
set_sizes$n_genes_tested <- sum(keep)   # the null PART E's percentiles are against
Lk   <- L[keep, , drop = FALSE]
bmk  <- rowMeans(NC)[keep]
mm   <- stats::model.matrix(~ m * a + epi + imm,
                            data.frame(m = sm$myc, a = Z(X$ox_rel), epi = epi, imm = imm))
ixc  <- grep(":a$", colnames(mm), value = TRUE)[1]
bvec <- stats::lm.fit(mm, t(zrow(Lk)))$coefficients[ixc, ]
names(bvec) <- rownames(Lk)
bins <- cut(log10(bmk + 1),
            breaks = stats::quantile(log10(bmk + 1), seq(0, 1, length.out = NBIN + 1)),
            include.lowest = TRUE)
calib <- dplyr::bind_rows(lapply(ROSTER, function(g) {
  e <- ens_of(g)
  if (is.na(e) || !(e %in% names(bvec)))
    return(tibble::tibble(gene = g, beta = NA_real_, pct_all_genes = NA_real_,
                          pct_expression_matched = NA_real_, n_matched = NA_integer_))
  b    <- bins[match(e, names(bvec))]
  pool <- setdiff(names(bvec)[bins == b], e)
  tibble::tibble(gene = g, beta = bvec[[e]],
                 pct_all_genes = 100 * mean(bvec < bvec[[e]]),
                 pct_expression_matched = 100 * mean(bvec[pool] < bvec[[e]]),
                 n_matched = length(pool))
}))

endpoint_menu <- dplyr::bind_rows(lapply(names(Y), function(n) {
  d  <- data.frame(y = Z(Y[[n]]), m = sm$myc, a = Z(X$ox_rel), epi = epi, imm = imm)
  s  <- summary(stats::lm(y ~ m * a + epi + imm, d))$coefficients
  dw <- data.frame(y = ctr(Y[[n]]), m = sm$myc, a = Z(ctr(X$ox_rel)),
                   epi = ctr(epi), imm = ctr(imm))
  sw <- summary(stats::lm(y ~ m * a + epi + imm, dw))$coefficients
  wt <- stats::coef(stats::lm(y ~ a + epi + imm, d[d$m == "neg", ]))[["a"]]
  pz <- stats::coef(stats::lm(y ~ a + epi + imm, d[d$m == "pos", ]))[["a"]]
  tibble::tibble(endpoint = n, bMX = s[ixof(s), 1], p = s[ixof(s), 4],
                 bMX_within_tp = sw[ixof(sw), 1], p_within_tp = sw[ixof(sw), 4],
                 slope_wt = wt, slope_myc = pz)
}))

# =============================================================================
# PART F: THE BUFFERING ARM
# -----------------------------------------------------------------------------
# The human model's escape clause is that tumours occupying the lethal corner
# bought their way in by buffering. That cannot be tested here -- there are no
# tumours. What CAN be tested is whether the guardian pair already moves under the
# gate in normal epithelium, which is what decides whether the human buffering
# test starts from a mouse-side prior or from nothing.
# =============================================================================
message("48 PART F: the buffering arm")

d_buf <- data.frame(y = Z(Y$BUFFER), m = sm$myc, a = Z(X$ox_rel), epi = epi, imm = imm)
buffer_fit <- as.data.frame(summary(stats::lm(y ~ m * a + epi + imm, d_buf))$coefficients)
buffer_fit$term <- rownames(buffer_fit)
buffer_perm <- perm_gate(Y$BUFFER, X$ox_rel)

# =============================================================================
# PART G: WITHIN-TIMEPOINT ONLY -- THE CROSS-SECTIONAL ANALOGUE
# =============================================================================
message("48 PART G: within-timepoint only")

within_tp <- dplyr::bind_rows(lapply(c("ox_rel", "ox_ppd", "ox_lvl"), function(an) {
  d <- data.frame(y = ctr(Y$PUMA), m = sm$myc, a = Z(ctr(X[[an]])),
                  epi = ctr(epi), imm = ctr(imm))
  s <- summary(stats::lm(y ~ m * a + epi + imm, d))$coefficients
  db <- data.frame(y = ctr(Y$BUFFER), m = sm$myc, a = Z(ctr(X[[an]])),
                   epi = ctr(epi), imm = ctr(imm))
  sb <- summary(stats::lm(y ~ m * a + epi + imm, db))$coefficients
  tibble::tibble(axis = an,
                 bMX_puma = s[ixof(s), 1], p_puma = s[ixof(s), 4],
                 bMX_buffer = sb[ixof(sb), 1], p_buffer = sb[ixof(sb), 4])
}))

# =============================================================================
# PART H: POWER TRANSFER
# -----------------------------------------------------------------------------
# SE of an OLS coefficient scales as 1/sqrt(n) at fixed design correlation, so the
# mouse SE at n = 24 fixes the whole curve. 2.8 SE is the usual 80 pct / two-sided
# 5 pct rule of thumb. The human effect will be ATTENUATED by survivor bias -- the
# tumours that expressed the gate most strongly are the ones that are not in the
# cohort -- so the third-of-the-effect row is the one to budget on.
# =============================================================================
message("48 PART H: power transfer")

n_for <- function(se, n0, eff) ceiling((2.8 * se * sqrt(n0) / eff)^2)
power_transfer <- dplyr::bind_rows(lapply(c("mRNA", "msig", "dor"), function(mn) {
  f <- gate_fit(Y$PUMA, M[[mn]], X$ox_rel)
  tibble::tibble(estimator = mn, bMX_per_sd = f$bMX, se = f$se, n_mouse = nrow(sm),
                 n_for_full_effect  = n_for(f$se, nrow(sm), f$bMX),
                 n_for_half_effect  = n_for(f$se, nrow(sm), f$bMX / 2),
                 n_for_third_effect = n_for(f$se, nrow(sm), f$bMX / 3))
}))

# =============================================================================
# PART I: VERDICT
# =============================================================================
message("48 PART I: verdict")

classify <- function(axis) {
  b <- gate_axes[gate_axes$axis == axis, ]
  e <- gate_estimators[gate_estimators$axis == axis &
                       gate_estimators$estimator == "msig", ]
  if (nrow(e) == 0 || b$bMX <= 0 || b$p >= 0.05) return("fails")
  if (e$frac_of_genotype >= 0.5) "transfers" else "transfers_weak"
}
verdict <- dplyr::bind_rows(lapply(c("ox_rel", "ox_ppd", "ox_lvl"), function(an)
  tibble::tibble(axis = an, verdict = classify(an),
                 bMX = gate_axes$bMX[gate_axes$axis == an],
                 p = gate_axes$p[gate_axes$axis == an],
                 perm_p = gate_axes$perm_p[gate_axes$axis == an],
                 frac_recovered_by_signature =
                   gate_estimators$frac_of_genotype[
                     gate_estimators$axis == an & gate_estimators$estimator == "msig"])))

# =============================================================================
# ASSERTS
# =============================================================================
stopifnot(
  nrow(gate_axes) == length(X),
  all(c("ox_rel", "ox_ppd", "ox_lvl") %in% gate_axes$axis),
  all(ROSTER %in% gene_gate$gene),
  # the two halves of the pre-specified ratio must both have been fitted alone,
  # or PART E cannot answer the question it exists for
  all(c("Bbc3", "Bcl2l1") %in% calib$gene[!is.na(calib$beta)]),
  nrow(power_transfer) == 3L,
  # both co-primaries must be present, or the gate-collapse claim is unscoped
  all(c("PUMA", "BUFFER") %in% model_form$endpoint))

NOTES <- c(
  "MODEL: T = b0 + bM*M + bX*X + bMX*(M x X) + covariates, T a BH3-sensor:guardian",
  "  log-ratio, M MYC dose, X the respiratory state of the mitochondrial compartment.",
  "  H-gate bMX > 0; H-null bX = 0 at M = 0; crossover M* = -bX/bMX in SD of M.",
  "SCOPE: n = 24, batch = timepoint. PART B carries the cohort contrast; PART G",
  "  removes it and is the number to quote for a cross-sectional claim. Ranking and",
  "  measurement design, not confirmatory inference. Causal claims belong to the",
  "  perturbations, not to 24 correlated mice.",
  "TRANSFER: the three questions this script exists to answer are (1) PART C, does a",
  "  continuous MYC estimator recover the genotype effect; (2) PART B, which OXPHOS",
  "  ruler; (3) PART E, which half of the ratio. Read those three before anything else.",
  "PART F is the only mouse-side prior the human buffering hypothesis has. It is a",
  "  correlation in normal epithelium, not a statement about tumours.")

res <- list(
  control = ctrl, separability = separability, estimator_gap = estimator_gap,
  set_sizes = set_sizes, axis_cor = axis_cor,
  gate_axes = gate_axes, head_to_head = head_to_head,
  gate_estimators = gate_estimators, model_form = model_form,
  simple_slopes = simple_slopes, gene_gate = gene_gate, calibration = calib,
  endpoint_menu = endpoint_menu, buffer_fit = buffer_fit, buffer_perm = buffer_perm,
  within_tp = within_tp, power_transfer = power_transfer, verdict = verdict,
  verdict_rule = VERDICT_RULE, geno_ref = geno_ref,
  scores = list(X = X, M = M, Y = Y, epi = epi, imm = imm, pro = pro,
                group = sm$group, tp = sm$tp, myc = sm$myc),
  params = list(NPERM = NPERM, NBIN = NBIN, MINBM = MINBM, seed = 1),
  analysis_date = Sys.Date(), notes = NOTES)

saveRDS(res, here::here("results", "gate_model_verification.rds"))
message("48: wrote results/gate_model_verification.rds")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  res <- readRDS(here::here("results", "gate_model_verification.rds"))

  ## The positive control first. If this is not script 43's +6.089 / 0.00523 the
  ## load block has drifted and nothing below is readable.
  res$control |> print()

  ## --- THE THREE ANSWERS THE HUMAN ARM NEEDS -------------------------------

  ## (1) WHICH RULER. Read bMX, p and perm_p together: the permutation median is
  ## not zero, so the excess over it is the claim, not the coefficient itself.
  ## Expect ox_rel >= ox_ppd > ox_lvl. If ox_lvl is the weakest, the human plan's
  ## primary axis measure is the wrong one and section 7.2 needs rewriting.
  res$gate_axes |> print(n = 20)

  ## The one control arm that can rival OXPHOS is the mitoribosome, and it is
  ## ~0.75 correlated with it. Only the head-to-head separates them: read the
  ## `mpos:ox` and `mpos:mr` rows.
  res$head_to_head |> print()

  ## (2) DOES A CONTINUOUS MYC ESTIMATOR TRANSFER. frac_of_genotype is the
  ## instrument's efficiency. Anything well under 1 means the human study must be
  ## powered on the estimator's own effect size, not on the mouse genotype effect.
  res$gate_estimators |> print(n = 20)
  res$estimator_gap |> print()

  ## Why: an interaction needs separable inputs. This is the diagnostic.
  res$separability |> print()

  ## (3) WHICH HALF OF THE RATIO. Fit the numerator and the denominator alone.
  ## If BCL2L1 carries as much as BBC3, the endpoint is a BALANCE and must be
  ## described as one -- and the several "specific" endpoints that share the
  ## denominator are not independent.
  res$gene_gate |> print(n = 20)
  res$calibration |> print(n = 20)
  res$endpoint_menu |> print(n = 20)

  ## --- MODEL FORM ----------------------------------------------------------
  ## p_drop_X_main is the H-null test. A large p means the OXPHOS main effect can
  ## be dropped: the model is a pure gate and has one parameter fewer.
  res$model_form |> print()

  ## The crossover is only meaningful if it falls INSIDE the observed MYC range.
  ## Check it against the simple slopes: an unadjusted negative wild-type slope
  ## that goes to zero on covariates is a composition effect, not a crossover.
  res$simple_slopes |> print()

  ## --- THE BUFFERING ARM ---------------------------------------------------
  res$buffer_fit |> print()
  res$buffer_perm |> print()

  ## --- THE HONEST TRANSFER ESTIMATE ----------------------------------------
  ## Cross-sectional human data has only within-cohort variation. This is that
  ## regime, and it is the weakest. Quote it when the claim is cross-sectional.
  res$within_tp |> print()

  ## --- WHAT n THE HUMAN COHORTS NEED ---------------------------------------
  res$power_transfer |> print()

  ## --- VERDICT -------------------------------------------------------------
  res$verdict |> print()
  res$verdict_rule |> print()
  cat(res$notes, sep = "\n")
}
