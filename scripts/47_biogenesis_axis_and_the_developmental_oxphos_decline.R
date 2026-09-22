# =============================================================================
# 47_biogenesis_axis_and_the_developmental_oxphos_decline.R
# -----------------------------------------------------------------------------
# IS THE WILD-TYPE 6W->12W OXPHOS DECLINE A PGC1a-AXIS CHANGE? And if it is not,
# what is the TF-activity change, and what licenses the PGC1a/NRF1 rescue in the
# experimental part?
#
# The question was raised repeatedly across Block A and never closed. The reason
# was not power. It was that BOTH obvious instruments are broken for this contrast
# and nobody had measured how badly:
#   (i)  the `_MITO` TF lanes are BUILT from MitoCarta genes, so a lane whose
#        overlap is OXPHOS-heavy falls because OXPHOS falls -- by construction;
#   (ii) the `_GRAY_<context>` lanes report a CELL STATE. Within one context every
#        TF moves together, so the layer ranks contexts, not factors.
# PART C measures both, and fixes the reading rule that follows.
#
# PART A -- THE ARM MAP. Every MitoCarta arm on the content ruler (set-mean raw
#   LFC, matched-random-set null), wild-type timeline, labelled chain_structural /
#   chain_assembly / organelle_build / metabolic / other. The claim under test is
#   that the decline is arm-selective. A biogenesis-axis withdrawal predicts the
#   structural subunits AND the assembly factors AND the mitoribosome AND the
#   import channels moving together; arm-selectivity is the alternative.
#
# PART B -- THE WITHIN-REGULON SPLIT, the decisive test. Script 44 did this for
#   CORE_MITO alone (`core_decomp`) and found the OXPHOS part at the 0.05th
#   percentile against the rest at the 99.55th. Here it is generalised to seven
#   regulons. A change in a factor's ACTIVITY acts on its REGULON; it cannot act on
#   a functional subset of the regulon. If every regulon splits the same way, the
#   split is a property of the OXPHOS genes and not of any factor.
#
# PART C -- WHAT THE TF LAYER CAN AND CANNOT SAY (see above). Ambient over all 421
#   TF lanes; adjusted-R2 of context vs TF identity over the 352 Gray lanes;
#   within-context spread against between-context spread; and the demonstration
#   rows -- a general transcription factor and a MICOS structural protein ranking
#   BELOW ESRRA in the same context. Reading rule fixed BEFORE the numbers are
#   read: a lane may be read as TF activity only if it beats BOTH the layer ambient
#   AND its own context median.
#
# PART D -- THE FACTORS THEMSELVES. Roster of biogenesis TFs, coactivators and
#   corepressors: baseMean, all four contrasts, interaction padj, and an
#   expression-matched-null percentile for the wild-type LFC. Carries the fact that
#   decides how the rescue must be described: Ppargc1a is at the floor of
#   expression in MEC, and the only PGC-1 family member expressed at a real level
#   is Pprc1 (PRC), which is also the only one that moves.
#
# PART E -- WITHIN-TIMEPOINT PER-SAMPLE COUPLINGS. Timepoint means removed, so the
#   between-cohort contrast -- which is also the batch contrast -- is gone
#   ENTIRELY. What does the OXPHOS composite track among animals of the same age?
#   Percentile against all expressed genes, which is the only null that means
#   anything at this n (scripts 33/35: everything correlates with everything).
#
# PART F -- THE CARRIER, AND THE LIMIT THAT STOPS IT BEING A RESULT. Reproduces
#   script 44's `le_within`. The adjusted betas are reported as a BOUND, not as a
#   mediation estimate: carrier and cargo correlate at ~0.93, and regressing one on
#   the other at that collinearity is not identified. Naming the design that would
#   settle it is the deliverable here, not a p-value.
#
# PART G -- VERDICT, on a rule fixed before the tests (script 46's pattern).
#
# PART H -- WHAT LICENSES THE PGC1a/NRF1 RESCUE. If the axis did not cause the
#   decline, the experimental part needs its justification stated in numbers, not
#   in prose: the axis's REACH over the genes that actually fell; the CONVERGENCE
#   of the two interventions on the respiratory arm; the fact that the death
#   effectors are NOT in the PGC1a regulon (so a restored death phenotype cannot be
#   direct transcription of the machinery); the CALIBRATION target the mouse sets
#   for the rescue; and the OVERSHOOT that must be pre-specified because it will
#   happen.
#
# SCOPE, stated once and not repeated. BATCH = TIMEPOINT (CLAUDE.md), so every
# between-age number here is confounded with cohort. What survives that is (a) the
# ARM-SELECTIVITY -- a shared batch effect cannot move one arm of one regulon to
# the 0.05th percentile while moving the rest of the same regulon to the 99.55th --
# and (b) PART E, which removes the timepoint means entirely. What does NOT survive
# is the MAGNITUDE of the developmental decline. n = 24 (n = 12 wild-type): this is
# RANKING plus a set of negatives, not confirmatory inference.
#
# Reads : results/interaction_results.rds, results/dds_int_run.rds,
#         results/combined_df_annotated.rds,
#         results/fgsea_percategory.rds          (script 20, the NES layer),
#         results/collapse_module_ownership.rds  (script 44, POSITIVE CONTROL),
#         results/biogenesis_discrimination.rds  (script 24, TF lane summary),
#         data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt,
#         functions/reconcile_gene_symbols.R (MANDATORY -- vintage-aware membership)
# Writes: results/biogenesis_axis_developmental.rds
#
# RUNTIME: a few minutes. PART A draws a matched null for ~65 arms and PART B for
# ~21 regulon parts, 2000 draws each; PART E correlates ~17k genes against a
# composite twice. Progress is printed.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

set.seed(1)
NSET   <- 2000L    # matched-random-set draws (script 43/44 idiom)
NBIN   <- 20L      # baseMean bins for the matched null
MIN_N  <- 8L       # smallest arm worth a null

# POSITIVE CONTROLS -- fixed constants, checked before anything new is computed.
OX_WT_REF        <- -0.2548   # script 44 wt_content, OXPHOS subunits, WT timeline
CORE_OX_REF      <- -0.2159   # script 44 core_decomp, CORE_MITO n OXPHOS subunits
CORE_REST_REF    <-  0.0689   # script 44 core_decomp, CORE_MITO rest
NES_OX_REF       <- -2.670    # script 20 fgsea_percategory, timepoint_neg
CTRL_TOL_CONTENT <- 1e-6      # same code path, same object -> must be exact
CTRL_TOL_NES     <- 5e-3      # read back from a stored table

# PART C reading rule, fixed before the numbers are looked at. TWO-SIDED on
# purpose: this contrast is negative-going as a whole (ambient NES is about -1.5),
# so a one-sided rule would score every falling lane as a signal and miss every
# rising one -- and the rising lanes are where the only readable TF statement is.
LANE_RULE <- paste(
  "`_MITO` lanes are NEVER read as TF activity: they are built from MitoCarta",
  "genes, so an OXPHOS-heavy overlap falls because OXPHOS falls. A Gray context",
  "lane is readable only if it departs from its OWN context median by more than",
  "that context's IQR, in the SAME direction as its departure from the layer",
  "ambient. A non-Gray lane (ChIP-Atlas / DoRothEA / MSigDB / Chung) carries no",
  "context confound and is readable if it departs from the layer ambient by more",
  "than the layer IQR. Everything else reports cell state or gene content.")

# PART G verdict rule, fixed before the tests.
VERDICT_RULE <- c(
  axis_change   = "regulon moves as a unit: whole-regulon null percentile outside [5, 95]",
  arm_selective = "OXPHOS part below the 5th pct AND the rest above the 50th, same regulon",
  inconclusive  = "neither")

# =============================================================================
# PART 0: LOAD, ALIGN, AND PROVE THE MACHINERY REPRODUCES SCRIPT 44
# =============================================================================
message("47 PART 0: load")

ir  <- readRDS(here::here("results", "interaction_results.rds"))
dds <- readRDS(here::here("results", "dds_int_run.rds"))
cdf <- readRDS(here::here("results", "combined_df_annotated.rds"))
fgc <- readRDS(here::here("results", "fgsea_percategory.rds"))
c44 <- readRDS(here::here("results", "collapse_module_ownership.rds"))
b24 <- readRDS(here::here("results", "biogenesis_discrimination.rds"))
gmt <- fgsea::gmtPathways(here::here("data", "genesets_from_library",
                                     "mammary_mito_myc_metab_v1_mouse.gmt"))

CONTRASTS <- c("myc_6W_raw", "myc_12W_raw", "timepoint_neg_raw",
               "timepoint_pos_raw", "interaction_raw")
D  <- lapply(ir[CONTRASTS], as.data.frame)
nc <- DESeq2::counts(dds, normalized = TRUE)
L  <- log2(nc + 1)
sm <- as.data.frame(SummarizedExperiment::colData(dds))
sm$tp  <- stats::relevel(as.factor(sm$timepoint),  "6W")
sm$myc <- stats::relevel(as.factor(sm$myc_status), "neg")
stopifnot(identical(rownames(sm), colnames(L)))

universe_all <- rownames(D$myc_6W_raw)
V   <- function(k) stats::setNames(D[[k]]$log2FoldChange, rownames(D[[k]]))
m6  <- V("myc_6W_raw");        m12 <- V("myc_12W_raw")
tn  <- V("timepoint_neg_raw"); tpz <- V("timepoint_pos_raw")
bm  <- stats::setNames(D$myc_6W_raw$baseMean, rownames(D$myc_6W_raw))

# Membership ALWAYS through the reconciler (CLAUDE.md; the 2026-07-24 fix). A naive
# symbol match loses renamed genes silently -- for MITOCARTA_OXPHOS_SUBUNITS it
# recovers 69 of 89 where the reconciler recovers 87.
ens_set  <- function(syms) recon_to_ensembl(syms, universe_all)
set_mean <- function(v, e) if (length(e)) mean(v[e], na.rm = TRUE) else NA_real_
zrow     <- function(m) t(scale(t(m)))

sym2ens <- cdf[!is.na(cdf$mgi_symbol) & !duplicated(cdf$mgi_symbol),
               c("mgi_symbol", "gene")]
ens_of  <- function(s) sym2ens$gene[match(s, sym2ens$mgi_symbol)]
sym_of  <- function(e) cdf$mgi_symbol[match(e, cdf$gene)]

# the matched-random-set null, script 43/44 idiom -----------------------------
expressed <- names(bm)[is.finite(bm) & bm > 0 & is.finite(tn[names(bm)])]
bin_of    <- cut(rank(bm[expressed], ties.method = "first"),
                 breaks = NBIN, labels = FALSE)
names(bin_of) <- expressed
by_bin    <- split(expressed, bin_of)
draw_matched <- function(e) {
  b <- bin_of[e]; b <- b[!is.na(b)]
  unlist(lapply(split(b, b), function(k)
    sample(by_bin[[as.character(k[1])]], length(k), replace = TRUE)),
    use.names = FALSE)
}
null_pct <- function(v, e, n = NSET) {
  e <- e[e %in% expressed]
  if (length(e) < 3) return(c(observed = NA_real_, null_median = NA_real_,
                              percentile = NA_real_))
  obs <- set_mean(v, e)
  nul <- vapply(seq_len(n), function(i) set_mean(v, draw_matched(e)), numeric(1))
  c(observed = obs, null_median = stats::median(nul), percentile = 100 * mean(nul < obs))
}

mito_universe <- unique(unlist(gmt[grep("^MITOCARTA_", names(gmt))], use.names = FALSE))
mito_ens      <- ens_set(mito_universe)

# GUARD: every named set must exist before anything uses it (script 42's failure
# mode was a wrong name giving a silent all-NaN composite).
NEEDED <- c("MITOCARTA_OXPHOS_SUBUNITS", "MITOCARTA_OXPHOS",
            "MITOCARTA_OXPHOS_ASSEMBLY_FACTORS", "MITOCARTA_MITOCHONDRIAL_RIBOSOME",
            "MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA", "MITOCARTA_MTDNA_MAINTENANCE",
            "MITOCARTA_APOPTOSIS_PRO", "MITOCARTA_APOPTOSIS_ANTI",
            "CORE_MITO", "ESRRA_MITO", "NRF1_MITO", "GABPA_MITO", "E2F1_MITO",
            "MYC_MITO", "MYC_SPECIFIC_MITO", "DEVELOPMENTAL_MITO",
            "MG_HEVSLE_AP_GRAY_DN", "MG_HEVSLE_AP_GRAY_UP",
            "MG_TEB_VS_DUCTAL_HS_GRAY_UP", "PROLIF_CELL_CYCLE_REACTOME")
missing_sets <- NEEDED[!NEEDED %in% names(gmt)]
if (length(missing_sets))
  stop("47: set names absent from the GMT -> ", paste(missing_sets, collapse = ", "))

# --- POSITIVE CONTROLS -------------------------------------------------------
# Nothing new is computed until the machinery reproduces what is already on disk.
ox_e   <- ens_set(gmt[["MITOCARTA_OXPHOS_SUBUNITS"]])
core_e <- ens_set(gmt[["CORE_MITO"]])
ctrl <- c(
  ox_wt        = set_mean(tn, ox_e),
  core_ox      = set_mean(tn, intersect(core_e, ox_e)),
  core_rest    = set_mean(tn, setdiff(core_e, ox_e)))
ctrl_ref <- c(ox_wt = OX_WT_REF, core_ox = CORE_OX_REF, core_rest = CORE_REST_REF)
ctrl_dev <- abs(ctrl - ctrl_ref)

nes_ox <- fgc$fgsea$NES[fgc$fgsea$ranking == "timepoint_neg" &
                        fgc$fgsea$pathway == "MITOCARTA_OXPHOS_SUBUNITS"]
stopifnot(length(nes_ox) == 1L)

# script 44's own stored values, so the check is against the object and not only
# against a transcribed constant
c44_ox   <- c44$wt_content$c_wt_time[c44$wt_content$arm == "OXPHOS subunits"]
c44_cox  <- c44$core_decomp$c_wt_time[c44$core_decomp$part == "CORE_MITO n OXPHOS subunits"]
c44_crst <- c44$core_decomp$c_wt_time[c44$core_decomp$part == "CORE_MITO, rest"]

message(sprintf("47 PART 0 controls: OXPHOS %.4f (44: %.4f, ref %.4f) | CORE n OX %.4f | CORE rest %+.4f",
                ctrl[["ox_wt"]], c44_ox, OX_WT_REF, ctrl[["core_ox"]], ctrl[["core_rest"]]))
message(sprintf("47 PART 0 controls: OXPHOS NES %.3f (ref %.3f)", nes_ox, NES_OX_REF))

if (max(abs(c(ctrl[["ox_wt"]] - c44_ox, ctrl[["core_ox"]] - c44_cox,
              ctrl[["core_rest"]] - c44_crst))) > CTRL_TOL_CONTENT)
  stop("47 PART 0: the content ruler does not reproduce script 44 -- stop and fix ",
       "before reading anything below.")
if (max(ctrl_dev) > 0.02)
  stop("47 PART 0: content ruler departs from the recorded references.")
if (abs(nes_ox - NES_OX_REF) > CTRL_TOL_NES)
  stop("47 PART 0: fGSEA NES does not reproduce script 20.")
message("47 PART 0: positive controls PASS")

# =============================================================================
# PART A: THE ARM MAP -- is the decline arm-selective?
# -----------------------------------------------------------------------------
# One ruler, every arm. The labels are the corpus's own "run vs build" split
# (docs/2026-07-19_oxphos_axis_biology_and_mtdna_priming.md section 3), made
# explicit here and stored in defs so the assignment is auditable rather than
# asserted. chain_assembly is separated from organelle_build on purpose: the
# assembly factors are the one arm that builds the chain itself, so they are the
# sharpest control on "the gland stopped building respiratory complexes".
# =============================================================================
message("47 PART A: arm map (matched nulls, this takes a minute)")

arm_class <- function(nm) {
  s <- sub("^MITOCARTA_", "", nm)
  if (grepl("ASSEMBLY_FACTORS", s))                                  return("chain_assembly")
  # the 13 mtDNA-encoded subunits get their OWN class so the confounded arm can
  # never blend into the clean one when a class median is taken (CLAUDE.md mt-axis)
  if (s %in% c("MTDNA_ENCODED", "OXPHOS_MT"))                        return("chain_mtDNA")
  if (grepl("^(OXPHOS|OXPHOS_NU|COMPLEX_|C(I|II|III|IV|V)_SUBUNITS)", s) ||
      s %in% c("OXPHOS_SUBUNITS", "ELECTRON_CARRIERS"))              return("chain_structural")
  if (grepl("^(MITOCHONDRIAL_RIBOSOME|TRANSLATION|MT_TRNA|MT_RRNA|MTRNA|MITOCHONDRIAL_CENTRAL_DOGMA|MTDNA_MAINTENANCE|PROTEIN_IMPORT|CHAPERONES|PROTEIN_HOMEOSTASIS|CARDIOLIPIN)", s))
    return("organelle_build")
  if (grepl("METABOLISM|^TCA|FATTY_ACID|NUCLEOTIDE|FOLATE|HEME|VITAMIN|PYRUVATE|LYSINE|GLYCINE|BRANCHED_CHAIN|CARBOHYDRATE|AMINO_ACID|LIPID|XENOBIOTIC", s))
    return("metabolic")
  "other"
}

mito_arms <- grep("^MITOCARTA_", names(gmt), value = TRUE)
nes_wt <- fgc$fgsea |>
  dplyr::filter(ranking == "timepoint_neg") |>
  dplyr::select(pathway, NES, pval, padj_within_category, size)

armmap <- purrr::map_dfr(mito_arms, function(a) {
  e <- ens_set(gmt[[a]])
  if (length(e) < MIN_N) return(NULL)
  n <- null_pct(tn, e)
  tibble::tibble(
    arm = a, class = arm_class(a),
    n_symbols = length(gmt[[a]]), n_genes = length(e),
    c_wt_time = n[["observed"]], null_median = n[["null_median"]],
    wt_null_pct = n[["percentile"]],
    c_mycpos_time = set_mean(tpz, e),
    c_myc_6W = set_mean(m6, e), c_myc_12W = set_mean(m12, e))
}) |>
  dplyr::left_join(nes_wt, by = c("arm" = "pathway")) |>
  dplyr::arrange(c_wt_time)

arm_class_summary <- armmap |>
  dplyr::group_by(class) |>
  dplyr::summarise(n_arms = dplyr::n(),
                   median_content = stats::median(c_wt_time, na.rm = TRUE),
                   median_pct     = stats::median(wt_null_pct, na.rm = TRUE),
                   median_NES     = stats::median(NES, na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::arrange(median_content)

# The mtDNA-encoded arm is reported, and nothing rests on it: CLAUDE.md flags the
# mt axis as time-associated (p = 0.0003) and three-way confounded (content x
# proliferation denominator x prep leak). Recorded here so the confound travels
# with the number instead of being rediscovered.
mt_genes <- grep("^mt-", cdf$mgi_symbol, value = TRUE)
mt_arm <- {
  e <- ens_of(mt_genes); e <- e[!is.na(e) & e %in% expressed]
  tibble::tibble(arm = "mtDNA-encoded (13)", n_genes = length(e),
                 c_wt_time = set_mean(tn, e), c_mycpos_time = set_mean(tpz, e),
                 note = "CONFOUNDED -- see CLAUDE.md mt-axis note; no argument rests on this")
}

# =============================================================================
# PART B: THE WITHIN-REGULON SPLIT -- the decisive test
# -----------------------------------------------------------------------------
# Generalises script 44's core_decomp from one regulon to seven. The logic is the
# whole argument: a change in a factor's ACTIVITY acts on its REGULON. If the
# OXPHOS members of a regulon fall while the other members do not, no activity
# change of that factor can be the explanation -- whatever moved the OXPHOS genes
# moved them as OXPHOS genes, not as targets.
# =============================================================================
message("47 PART B: within-regulon split")

REGULONS <- c("CORE_MITO", "ESRRA_MITO", "NRF1_MITO", "GABPA_MITO",
              "E2F1_MITO", "MYC_MITO", "DEVELOPMENTAL_MITO")

split_one <- function(rg) {
  re <- ens_set(gmt[[rg]])
  parts <- list(`n OXPHOS subunits` = intersect(re, ox_e),
                `rest`              = setdiff(re, ox_e),
                `all`               = re)
  purrr::imap_dfr(parts, function(e, nm) {
    if (length(e) < 3) return(NULL)
    n <- null_pct(tn, e)
    tibble::tibble(regulon = rg, part = nm, n_genes = length(e),
                   c_wt_time = n[["observed"]], null_median = n[["null_median"]],
                   wt_null_pct = n[["percentile"]],
                   c_mycpos_time = set_mean(tpz, e), c_myc_6W = set_mean(m6, e))
  })
}
split <- purrr::map_dfr(REGULONS, split_one)

verdict_of <- function(d) {
  p_all  <- d$wt_null_pct[d$part == "all"]
  p_ox   <- d$wt_null_pct[d$part == "n OXPHOS subunits"]
  p_rest <- d$wt_null_pct[d$part == "rest"]
  if (length(p_ox) && length(p_rest) && !is.na(p_ox) && !is.na(p_rest) &&
      p_ox < 5 && p_rest > 50) return("arm_selective")
  if (length(p_all) && !is.na(p_all) && (p_all < 5 || p_all > 95)) return("axis_change")
  "inconclusive"
}
split_verdict <- purrr::map_dfr(REGULONS, function(rg) {
  d <- split[split$regulon == rg, ]
  tibble::tibble(regulon = rg,
                 n_ox   = d$n_genes[d$part == "n OXPHOS subunits"][1],
                 n_rest = d$n_genes[d$part == "rest"][1],
                 pct_ox   = d$wt_null_pct[d$part == "n OXPHOS subunits"][1],
                 pct_rest = d$wt_null_pct[d$part == "rest"][1],
                 pct_all  = d$wt_null_pct[d$part == "all"][1],
                 verdict  = verdict_of(d))
})

# =============================================================================
# PART C: WHAT THE TF LAYER CAN AND CANNOT SAY
# -----------------------------------------------------------------------------
# Two failure modes, measured rather than assumed. The output is a reading rule,
# not a ranking of factors -- the point of this part is that the layer CANNOT rank
# factors in this contrast, which is why the question stayed open.
# =============================================================================
message("47 PART C: TF layer")

tf_wt <- fgc$fgsea |> dplyr::filter(ranking == "timepoint_neg",
                                    category == "06_tf_targets")
tf_ambient <- stats::median(tf_wt$NES, na.rm = TRUE)

gray <- tf_wt |> dplyr::filter(grepl("_GRAY_", pathway))
gb   <- sub("^TFT_", "", sub("_MITO$", "", gray$pathway))
gsp  <- strsplit(gb, "_GRAY_", fixed = TRUE)
gray$tf   <- vapply(gsp, `[`, character(1), 1L)
gray$ctx  <- vapply(gsp, `[`, character(1), 2L)
gray$mito <- grepl("_MITO$", gray$pathway)

adjr2 <- function(f) summary(stats::lm(f, data = gray))$adj.r.squared
tf_layer <- tibble::tibble(
  n_lanes_total   = nrow(tf_wt),
  n_lanes_gray    = nrow(gray),
  n_tf            = dplyr::n_distinct(gray$tf),
  n_context       = dplyr::n_distinct(gray$ctx),
  ambient_NES     = tf_ambient,
  ambient_iqr     = stats::IQR(tf_wt$NES, na.rm = TRUE),
  ambient_mito    = stats::median(tf_wt$NES[grepl("_MITO$", tf_wt$pathway)], na.rm = TRUE),
  ambient_nonmito = stats::median(tf_wt$NES[!grepl("_MITO$", tf_wt$pathway)], na.rm = TRUE),
  adj_r2_context  = adjr2(NES ~ ctx),
  adj_r2_tf       = adjr2(NES ~ tf),
  adj_r2_both     = adjr2(NES ~ ctx + tf),
  ctx_median_range = diff(range(tapply(gray$NES, gray$ctx, stats::median))))

tf_context <- gray |>
  dplyr::group_by(ctx, mito) |>
  dplyr::summarise(n = dplyr::n(), median_NES = stats::median(NES),
                   iqr_NES = stats::IQR(NES), range_NES = diff(range(NES)),
                   .groups = "drop") |>
  dplyr::arrange(median_NES)

# The demonstration: inside ONE context, rank the lanes. If a general transcription
# factor and a structural mitochondrial protein outrank the biogenesis factors, the
# lane is reporting content and state, not activity.
tf_demo <- gray |>
  dplyr::filter(mito, ctx == "AP_LE") |>
  dplyr::arrange(NES) |>
  dplyr::mutate(rank_in_context = dplyr::row_number()) |>
  dplyr::select(rank_in_context, tf, ctx, NES, size)

tf_layer_iqr <- stats::IQR(tf_wt$NES, na.rm = TRUE)

# apply LANE_RULE to EVERY lane, not only the roster, so the rule's own hit rate
# is visible: if it passes half the layer it is not a rule.
lane_rule_apply <- function(d) {
  d |>
    dplyr::mutate(
      is_mito_lane = grepl("_MITO$", pathway),
      is_gray      = grepl("_GRAY_", pathway),
      ctx = ifelse(is_gray,
                   sub("^TFT_[^_]+_GRAY_", "", sub("_MITO$", "", pathway)),
                   NA_character_)) |>
    dplyr::left_join(tf_context |> dplyr::filter(!mito) |>
                       dplyr::select(ctx, ctx_median = median_NES, ctx_iqr = iqr_NES),
                     by = "ctx") |>
    dplyr::mutate(
      dev_ambient = NES - tf_ambient,
      dev_context = NES - ctx_median,
      readable_as_activity = dplyr::case_when(
        is_mito_lane ~ FALSE,
        is_gray      ~ !is.na(dev_context) &
                       sign(dev_ambient) == sign(dev_context) &
                       abs(dev_context) > ctx_iqr,
        TRUE         ~ abs(dev_ambient) > tf_layer_iqr))
}

tf_all_lanes <- lane_rule_apply(tf_wt)
lane_rule_hits <- tf_all_lanes |>
  dplyr::group_by(is_mito_lane, is_gray) |>
  dplyr::summarise(n = dplyr::n(), n_readable = sum(readable_as_activity),
                   frac = mean(readable_as_activity), .groups = "drop")

# The context-free sub-layer is the ONLY place a TF statement can be made here:
# these lanes are generic regulons, so the Gray context artifact does not apply.
tf_contextfree <- tf_all_lanes |>
  dplyr::filter(!is_gray, !is_mito_lane) |>
  dplyr::arrange(dplyr::desc(NES)) |>
  dplyr::select(pathway, NES, pval, padj_within_category, size,
                dev_ambient, readable_as_activity)

ROSTER_TF <- c("ESRRA", "NRF1", "GABPA", "YY1", "MYC", "E2F1", "ESR1", "FOXO3",
               "PPARG", "SREBF1", "SREBF2", "CREB1", "CHCHD3", "MTERF3", "GTF3A")
tf_roster_lanes <- tf_all_lanes |>
  dplyr::mutate(tf = sub("^TFT_", "", sub("_(GRAY|CHIPATLAS|DOROTHEA|MSIGDB|CHUNG).*$", "",
                                          pathway))) |>
  dplyr::filter(tf %in% ROSTER_TF) |>
  dplyr::arrange(NES) |>
  dplyr::select(pathway, tf, ctx, NES, pval, padj_within_category, size,
                is_mito_lane, is_gray, ctx_median, ctx_iqr,
                dev_ambient, dev_context, readable_as_activity)

# script 24's lane-level summary, joined so the PGC1a axis is read against the
# layer ambient rather than against zero
tf_lane_summary_24 <- b24$tf_fgsea_summary |>
  dplyr::filter(contrast == "WT_6W->12W") |>
  dplyr::mutate(ambient_NES = tf_ambient,
                above_ambient = med_NES > tf_ambient)

# =============================================================================
# PART D: THE FACTORS THEMSELVES
# -----------------------------------------------------------------------------
# An expression-matched percentile per gene, because at this n a lone LFC says
# nothing: the question is whether the factor moves MORE than genes of the same
# abundance move. Ppargc1a's baseMean is the number that decides how the rescue
# has to be described in the paper.
# =============================================================================
message("47 PART D: factor roster")

ROSTER <- c(
  # PGC-1 family coactivators
  "Ppargc1a", "Ppargc1b", "Pprc1",
  # the axis TFs
  "Esrra", "Esrrb", "Esrrg", "Nrf1", "Gabpa", "Gabpb1", "Gabpb2", "Yy1",
  # mtDNA transcription/replication machinery
  "Tfam", "Tfb1m", "Tfb2m", "Polrmt", "Mterf3", "Mterf4",
  # coregulators
  "Nrip1", "Ncor1", "Ncor2", "Ncoa1", "Ncoa2", "Ncoa3", "Med1",
  # nutrient/quiescence arm
  "Sirt1", "Sirt3", "Prkaa1", "Prkaa2", "Mtor", "Foxo1", "Foxo3", "Foxo4",
  "Bnip3", "Bnip3l",
  # growth arm and controls
  "Myc", "Mycn", "Max", "Mxi1", "Mybl2", "Foxm1", "E2f1", "Tfdp1",
  # other candidates raised during Block A
  "Nfe2l1", "Nfe2l2", "Ppara", "Ppard", "Pparg", "Srebf1", "Srebf2",
  "Hif1a", "Epas1", "Tfeb", "Tfe3", "Esr1", "Pgr", "Gata3", "Elf5")

gene_pct <- function(e, v) {
  # percentile of this gene's LFC among expressed genes in its own baseMean bin
  if (is.na(e) || !e %in% expressed) return(NA_real_)
  b    <- bin_of[[e]]
  peer <- by_bin[[as.character(b)]]
  100 * mean(v[peer] < v[[e]], na.rm = TRUE)
}
padj_of <- function(k, e) if (is.na(e)) NA_real_ else D[[k]]$padj[match(e, rownames(D[[k]]))]

roster <- purrr::map_dfr(ROSTER, function(s) {
  e <- ens_of(s)
  tibble::tibble(
    sym = s, gene = e,
    baseMean   = if (is.na(e)) NA_real_ else as.numeric(bm[e]),
    wt_lfc     = if (is.na(e)) NA_real_ else as.numeric(tn[e]),
    wt_padj    = padj_of("timepoint_neg_raw", e),
    wt_pct     = gene_pct(e, tn),
    mycpos_lfc = if (is.na(e)) NA_real_ else as.numeric(tpz[e]),
    mycpos_padj = padj_of("timepoint_pos_raw", e),
    myc6_lfc   = if (is.na(e)) NA_real_ else as.numeric(m6[e]),
    myc6_padj  = padj_of("myc_6W_raw", e),
    int_padj   = padj_of("interaction_raw", e))
}) |> dplyr::arrange(wt_lfc)

pgc1_family <- roster |>
  dplyr::filter(sym %in% c("Ppargc1a", "Ppargc1b", "Pprc1")) |>
  dplyr::select(sym, baseMean, wt_lfc, wt_padj, wt_pct, mycpos_lfc, mycpos_padj, int_padj)

# =============================================================================
# PART E: WITHIN-TIMEPOINT PER-SAMPLE COUPLINGS
# -----------------------------------------------------------------------------
# Timepoint means removed within the genotype, so the between-cohort (= batch)
# contrast is gone. The percentile against every expressed gene is the null that
# matters: at n = 12 a raw rho of 0.8 is what the window hands you (scripts 33/35).
# =============================================================================
message("47 PART E: within-timepoint couplings")

couple_in <- function(which_myc) {
  idx <- which(sm$myc_status == which_myc)
  tpv <- droplevels(sm$tp[idx])
  M   <- L[, idx, drop = FALSE]
  M   <- M[rowMeans(nc[, idx, drop = FALSE]) >= 10, , drop = FALSE]
  # residualise on timepoint (removes the batch contrast entirely)
  R <- M
  for (lv in levels(tpv)) {
    k <- which(tpv == lv)
    R[, k] <- M[, k, drop = FALSE] - rowMeans(M[, k, drop = FALSE])
  }
  Z <- zrow(R)
  Z <- Z[is.finite(rowSums(Z)), , drop = FALSE]
  comp <- function(e) {
    e <- intersect(e, rownames(Z))
    if (length(e) < 5) return(NULL)
    colMeans(Z[e, , drop = FALSE])
  }
  oxc <- comp(ox_e)
  stopifnot(!is.null(oxc))
  allr <- as.numeric(stats::cor(t(Z), oxc))
  names(allr) <- rownames(Z)
  pct <- function(r) if (is.na(r)) NA_real_ else 100 * mean(allr < r, na.rm = TRUE)

  set_rows <- purrr::map_dfr(
    c("PROLIF_CELL_CYCLE_REACTOME", "MITOCARTA_OXPHOS_ASSEMBLY_FACTORS",
      "MITOCARTA_MITOCHONDRIAL_RIBOSOME", "CORE_MITO",
      "MG_HEVSLE_AP_GRAY_DN", "MG_TEB_VS_DUCTAL_HS_GRAY_UP"),
    function(s) {
      v <- comp(ens_set(gmt[[s]]))
      if (is.null(v)) return(NULL)
      tibble::tibble(genotype = which_myc, kind = "set", label = s,
                     r = stats::cor(v, oxc), pct = NA_real_)
    })
  gene_rows <- purrr::map_dfr(ROSTER, function(s) {
    e <- ens_of(s)
    if (is.na(e) || !e %in% rownames(Z)) return(NULL)
    tibble::tibble(genotype = which_myc, kind = "gene", label = s,
                   r = allr[[e]], pct = pct(allr[[e]]))
  })
  list(rows = dplyr::bind_rows(set_rows, gene_rows),
       meta = tibble::tibble(genotype = which_myc, n_samples = length(idx),
                             n_genes = nrow(Z), n_oxphos = length(intersect(ox_e, rownames(Z)))))
}
cpl_neg <- couple_in("neg"); cpl_pos <- couple_in("pos")
couplings <- dplyr::bind_rows(cpl_neg$rows, cpl_pos$rows) |>
  dplyr::arrange(genotype, dplyr::desc(r))
couplings_meta <- dplyr::bind_rows(cpl_neg$meta, cpl_pos$meta)

# =============================================================================
# PART F: THE CARRIER, AND WHY IT CANNOT BE SETTLED HERE
# -----------------------------------------------------------------------------
# Reproduces script 44's le_within and then states the limit that stops it being a
# mediation result. At cor ~ 0.93 between carrier and cargo the adjusted beta is
# not identified: a covariate that shares 87% of its variance with the outcome
# removes the outcome from itself. The adjusted numbers are a BOUND on how much of
# the decline COULD be carried, not an estimate of how much IS.
# =============================================================================
message("47 PART F: carrier")

wt_i <- which(sm$myc_status == "neg")
score_from <- function(ens, min_n = 5L) {
  e <- ens[ens %in% rownames(L)]
  if (length(e) < min_n) return(NULL)
  colMeans(zrow(L[e, , drop = FALSE]))
}
carrier <- local({
  ap_le_e  <- ens_set(gmt[["MG_HEVSLE_AP_GRAY_DN"]])
  teb_e    <- ens_set(gmt[["MG_TEB_VS_DUCTAL_HS_GRAY_UP"]])
  ox       <- score_from(ox_e)
  ap_le    <- score_from(ap_le_e)
  ap_le_nm <- score_from(setdiff(ap_le_e, mito_ens))
  ap_he_nm <- score_from(setdiff(ens_set(gmt[["MG_HEVSLE_AP_GRAY_UP"]]), mito_ens))
  teb_nm   <- score_from(setdiff(teb_e, mito_ens))
  b   <- function(f) unname(stats::coef(f)[2])
  vif <- function(x, y) 1 / (1 - stats::cor(x, y)^2)
  tibble::tibble(
    n_wt                  = length(wt_i),
    n_shared_ox_apLE      = length(intersect(ox_e, ap_le_e)),
    cor_ox_apLE_raw       = stats::cor(ox[wt_i], ap_le[wt_i]),
    cor_ox_apLE_nonmito   = stats::cor(ox[wt_i], ap_le_nm[wt_i]),
    cor_ox_apHE_nonmito   = stats::cor(ox[wt_i], ap_he_nm[wt_i]),
    cor_ox_teb_nonmito    = stats::cor(ox[wt_i], teb_nm[wt_i]),
    vif_apLE_nonmito      = vif(ox[wt_i], ap_le_nm[wt_i]),
    ox_beta_time_raw      = b(stats::lm(ox[wt_i] ~ sm$tp[wt_i])),
    ox_beta_time_adj_nm   = b(stats::lm(ox[wt_i] ~ sm$tp[wt_i] + ap_le_nm[wt_i])),
    ox_beta_time_adj_teb  = b(stats::lm(ox[wt_i] ~ sm$tp[wt_i] + teb_nm[wt_i])),
    ox_beta_time_adj_self = b(stats::lm(ox[wt_i] ~ sm$tp[wt_i] + ap_le[wt_i])),
    identified            = FALSE,
    limit_note = paste("cor(carrier, cargo) ~0.93 -> the adjusted betas are a BOUND,",
                       "not a mediation estimate. Separating them needs cell-resolved",
                       "measurement (sorted AP/HS/BA or single-cell), not more bulk."))
})

# =============================================================================
# PART G: VERDICT
# =============================================================================
verdict <- tibble::tibble(
  question = c("Does the PGC1a axis withdraw as a unit in the wild-type gland?",
               "Is the decline arm-selective within the regulon?",
               "Do the axis factors themselves move?",
               "Can the TF layer rank factors in this contrast?",
               "What does OXPHOS track among animals of the same age?"),
  answer = c(
    paste0("NO -- CORE_MITO whole-regulon percentile ",
           sprintf("%.1f", split_verdict$pct_all[split_verdict$regulon == "CORE_MITO"])),
    paste0("YES -- ", sum(split_verdict$verdict == "arm_selective"), " of ",
           nrow(split_verdict), " regulons split arm-selectively"),
    paste0("NO -- see roster; Ppargc1a baseMean ",
           sprintf("%.0f", roster$baseMean[roster$sym == "Ppargc1a"]),
           ", Pprc1 the only PGC-1 member that moves"),
    paste0("NO -- adj-R2 context ", sprintf("%.3f", tf_layer$adj_r2_context),
           " vs TF ", sprintf("%.3f", tf_layer$adj_r2_tf)),
    "see couplings: proliferation composite, not the biogenesis factors"),
  rule = c(VERDICT_RULE[["axis_change"]], VERDICT_RULE[["arm_selective"]],
           "expression-matched percentile outside [5, 95]", LANE_RULE,
           "percentile against all expressed genes"))

# =============================================================================
# PART H: WHAT LICENSES THE PGC1a/NRF1 RESCUE
# -----------------------------------------------------------------------------
# The experimental part reverts the OXPHOS decline with PGC1a and NRF1. If the axis
# did not CAUSE the decline, that choice needs a justification, and the justification
# should be numbers rather than prose. Five of them:
#
#   H1 REACH        -- the axis did not lower these genes but it reaches them.
#                      Reverting a deficit needs reach over the affected genes, not
#                      authorship of the deficit.
#   H2 CONVERGENCE  -- ESRRA and NRF1 regulons overlap little, and their INTERSECTION
#                      is where the respiratory subunits concentrate. Two interventions
#                      with different reach converging on one phenotype localise the
#                      effect to their shared part. This is why two factors is a
#                      control and not a redundancy.
#   H3 OWNERSHIP    -- the death effectors are NOT in the PGC1a regulon (they are in
#                      MYC's). So a restored death phenotype cannot be direct
#                      transcription of the machinery; it has to run through the
#                      respiratory state, which is the claim under test. The single
#                      exception, Cycs, is itself a respiratory carrier.
#   H4 CALIBRATION  -- the mouse sets the target effect size for the rescue. Since
#                      Ppargc1a is at the floor of expression in MEC, the transgene
#                      level is not a calibration; the OUTPUT is.
#   H5 OVERSHOOT    -- the developmental change was arm-selective and the rescue will
#                      not be. Pre-specify it: the arms that never fell are still
#                      inside the regulon, so they will rise. Naming this before the
#                      experiment is what keeps the claim to "raising respiratory
#                      capacity restores death competence" rather than "reverting the
#                      developmental change restores it".
# =============================================================================
message("47 PART H: rescue pre-specification")

ox_syms <- gmt[["MITOCARTA_OXPHOS_SUBUNITS"]]

reach <- purrr::map_dfr(REGULONS, function(rg) {
  s <- gmt[[rg]]; k <- intersect(ox_syms, s)
  tibble::tibble(regulon = rg, n_regulon = length(s),
                 covers_ox = length(k), n_ox_total = length(ox_syms),
                 frac_of_declining = length(k) / length(ox_syms),
                 ox_frac_of_regulon = length(k) / length(s))
}) |> dplyr::arrange(dplyr::desc(frac_of_declining))

convergence <- local({
  e <- gmt[["ESRRA_MITO"]]; n <- gmt[["NRF1_MITO"]]
  parts <- list(`ESRRA only`  = setdiff(e, n),
                `NRF1 only`   = setdiff(n, e),
                `shared`      = intersect(e, n),
                `union`       = union(e, n))
  purrr::imap_dfr(parts, function(s, nm) {
    k <- length(intersect(ox_syms, s))
    tibble::tibble(part = nm, n_genes = length(s), n_oxphos = k,
                   ox_fraction = k / length(s))
  })
})
# is the shared part enriched for the respiratory subunits relative to the union?
convergence_test <- local({
  u <- union(gmt[["ESRRA_MITO"]], gmt[["NRF1_MITO"]])
  s <- intersect(gmt[["ESRRA_MITO"]], gmt[["NRF1_MITO"]])
  a <- intersect(ox_syms, u); k <- length(intersect(ox_syms, s))
  tibble::tibble(
    n_union = length(u), n_shared = length(s), n_ox_in_union = length(a),
    n_ox_in_shared = k, expected = length(s) * length(a) / length(u),
    fold = k / (length(s) * length(a) / length(u)),
    p_enrich = stats::phyper(k - 1, length(a), length(u) - length(a), length(s),
                             lower.tail = FALSE))
})

ownership <- purrr::map_dfr(c("CORE_MITO", "ESRRA_MITO", "NRF1_MITO", "MYC_MITO"),
  function(pg) {
    s <- intersect(gmt[[pg]], mito_universe)
    purrr::map_dfr(c("MITOCARTA_APOPTOSIS_PRO", "MITOCARTA_APOPTOSIS_ANTI"), function(ap) {
      a <- intersect(gmt[[ap]], mito_universe)
      k <- length(intersect(s, a))
      e <- length(s) * length(a) / length(mito_universe)
      tibble::tibble(
        programme = pg, apoptosis_set = ap, n_programme = length(s),
        n_apoptosis = length(a), overlap = k, expected = e, fold = k / e,
        p_enrich = stats::phyper(k - 1, length(a), length(mito_universe) - length(a),
                                 length(s), lower.tail = FALSE),
        genes = paste(sort(intersect(s, a)), collapse = ", "),
        has_PUMA = "Bbc3" %in% intersect(s, a),
        has_BAX  = "Bax"  %in% intersect(s, a))
    })
  })

calibration <- local({
  d_wt  <- armmap$c_wt_time[armmap$arm == "MITOCARTA_OXPHOS_SUBUNITS"]
  d_myc <- armmap$c_mycpos_time[armmap$arm == "MITOCARTA_OXPHOS_SUBUNITS"]
  tibble::tibble(
    context = c("wild-type 6W->12W", "Myc+ 6W->12W"),
    deficit_log2 = c(d_wt, d_myc),
    remaining_fraction = 2 ^ c(d_wt, d_myc),
    percent_lost = 100 * (1 - 2 ^ c(d_wt, d_myc)),
    restoration_fold_needed = 1 / 2 ^ c(d_wt, d_myc),
    note = "target for the rescue, measured at OXPHOS-subunit level, not transgene level")
})

WATCH <- c("MITOCARTA_OXPHOS_ASSEMBLY_FACTORS", "MITOCARTA_MITOCHONDRIAL_RIBOSOME",
           "MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA", "MITOCARTA_MTDNA_MAINTENANCE")
overshoot <- purrr::map_dfr(WATCH, function(a) {
  s <- gmt[[a]]
  row <- armmap[armmap$arm == a, ]
  tibble::tibble(
    arm = a, n_genes = length(s),
    wt_content = if (nrow(row)) row$c_wt_time else NA_real_,
    wt_null_pct = if (nrow(row)) row$wt_null_pct else NA_real_,
    in_CORE_MITO  = length(intersect(s, gmt[["CORE_MITO"]])) / length(s),
    in_ESRRA_MITO = length(intersect(s, gmt[["ESRRA_MITO"]])) / length(s),
    in_NRF1_MITO  = length(intersect(s, gmt[["NRF1_MITO"]])) / length(s),
    prediction = "did NOT fall developmentally but IS inside the regulon -> expect it to rise")
})

rescue <- list(reach = reach, convergence = convergence,
               convergence_test = convergence_test, ownership = ownership,
               calibration = calibration, overshoot = overshoot)

# =============================================================================
# PART I: THE FOXO3 ARM -- the one readable TF signal, tested rather than told
# -----------------------------------------------------------------------------
# PART C leaves exactly one sub-layer in which a TF statement is possible, and
# FOXO3 is at the top of it. That makes it the strongest READABLE signal in this
# contrast, which is not the same as the strongest hypothesis. Four tests decide
# which it is, and they are the same tests PART H applies to PGC1a:
#
#   I1 CONTENT   -- is the lane's rise an OXPHOS artifact? (it is a rise, so the
#                   artifact would have to be a rising mito arm)
#   I2 REACH     -- does FOXO3 reach the genes that fell? A factor with no reach
#                   over the affected genes cannot be the proximal cause, however
#                   clean its own signal is. This is the test that decides
#                   "driver" against "marker", and it is the same one that
#                   disqualified the PGC1a axis as the cause.
#   I3 TARGETS   -- does the programme move gene by gene, and WHICH arm of it?
#                   FOXO3 has an arrest arm, an atrophy/turnover arm and an
#                   apoptotic arm, and they do not have to move together.
#   I4 SEPARABLE -- is the FOXO3 axis distinguishable from the carrier state and
#                   from OXPHOS itself, or is it another reading of the same
#                   variable? Percentile against all expressed genes, since at
#                   n = 12 a raw rho is the ceiling (scripts 33/35).
#
# The apoptotic arm matters beyond this script: FOXO3 is the PUMA regulator, so
# "the adult gland raises FOXO3 and arms PUMA" is the story this section is most
# likely to be read as supporting. I3 tests it explicitly rather than leaving it
# to be assumed.
# =============================================================================
message("47 PART I: the FOXO3 arm")

FOXO_SETS <- c("TFT_FOXO3_CHUNG", "TFT_FOXO1_CHUNG", "TFT_FOXO4_CHUNG")
stopifnot(all(FOXO_SETS %in% names(gmt)))

# --- I1 + I2: content and reach, on the same footing as PART H's reach table ---
foxo_reach <- purrr::map_dfr(FOXO_SETS, function(s) {
  m <- gmt[[s]]
  tibble::tibble(
    set = s, n_genes = length(m),
    n_mitocarta = length(intersect(m, mito_universe)),
    mito_members = paste(sort(intersect(m, mito_universe)), collapse = ", "),
    covers_ox = length(intersect(m, ox_syms)),
    n_ox_total = length(ox_syms),
    frac_of_declining = length(intersect(m, ox_syms)) / length(ox_syms))
})

# --- I3: the lane across all five rankings, and gene by gene ------------------
foxo_lanes <- fgc$fgsea |>
  dplyr::filter(pathway %in% FOXO_SETS) |>
  dplyr::select(pathway, ranking, NES, pval, padj_within_category, size) |>
  dplyr::arrange(pathway, ranking)

# The three arms of the FOXO3 programme, named BEFORE the values are read, so
# "which arm moved" is a result and not a description of whatever moved.
FOXO_ARMS <- tibble::tribble(
  ~arm,                  ~gene,
  "arrest",              "Cdkn1a",   "arrest",   "Cdkn1b",  "arrest", "Gadd45a",
  "arrest",              "Ccng2",    "arrest",   "Rbl2",
  "atrophy_turnover",    "Fbxo32",   "atrophy_turnover", "Bnip3",
  "atrophy_turnover",    "Bnip3l",   "atrophy_turnover", "Pink1",
  "atrophy_turnover",    "Gabarapl1","atrophy_turnover", "Sirt1",
  "antioxidant",         "Sod2",     "antioxidant", "Cat",  "antioxidant", "Prdx3",
  "antioxidant",         "Txnip",    "antioxidant", "Sesn1",
  "apoptotic",           "Bcl2l11",  "apoptotic", "Bbc3",   "apoptotic", "Fasl",
  "apoptotic",           "Tnfsf10")

foxo_targets <- purrr::pmap_dfr(FOXO_ARMS, function(arm, gene) {
  e <- ens_of(gene)
  tibble::tibble(
    arm = arm, sym = gene,
    in_chung_set = gene %in% gmt[["TFT_FOXO3_CHUNG"]],
    baseMean   = if (is.na(e)) NA_real_ else as.numeric(bm[e]),
    wt_lfc     = if (is.na(e)) NA_real_ else as.numeric(tn[e]),
    wt_padj    = padj_of("timepoint_neg_raw", e),
    wt_pct     = gene_pct(e, tn),
    mycpos_lfc = if (is.na(e)) NA_real_ else as.numeric(tpz[e]),
    myc6_lfc   = if (is.na(e)) NA_real_ else as.numeric(m6[e]),
    int_padj   = padj_of("interaction_raw", e))
}) |> dplyr::arrange(arm, dplyr::desc(wt_lfc))

foxo_arm_summary <- foxo_targets |>
  dplyr::group_by(arm) |>
  dplyr::summarise(n = dplyr::n(),
                   n_padj05 = sum(wt_padj < 0.05, na.rm = TRUE),
                   median_wt_lfc = stats::median(wt_lfc, na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::arrange(dplyr::desc(median_wt_lfc))

# --- I4: separability, on the non-mito part of the programme ------------------
# The 6 MitoCarta members are stripped first, or the coupling to a mitochondrial
# composite would be partly self-correlation.
foxo_separability <- purrr::map_dfr(c("neg", "pos"), function(gt) {
  idx <- which(sm$myc_status == gt)
  tpv <- droplevels(sm$tp[idx])
  M <- L[, idx, drop = FALSE]
  M <- M[rowMeans(nc[, idx, drop = FALSE]) >= 10, , drop = FALSE]
  R <- M
  for (lv in levels(tpv)) {
    k <- which(tpv == lv); R[, k] <- M[, k, drop = FALSE] - rowMeans(M[, k, drop = FALSE])
  }
  Z <- zrow(R); Z <- Z[is.finite(rowSums(Z)), , drop = FALSE]
  cmp <- function(syms) {
    e <- intersect(ens_set(syms), rownames(Z))
    if (length(e) < 5) return(NULL)
    colMeans(Z[e, , drop = FALSE])
  }
  f3 <- cmp(setdiff(gmt[["TFT_FOXO3_CHUNG"]], mito_universe))
  ox <- cmp(gmt[["MITOCARTA_OXPHOS_SUBUNITS"]])
  ap <- cmp(setdiff(gmt[["MG_HEVSLE_AP_GRAY_DN"]], mito_universe))
  allr <- as.numeric(stats::cor(t(Z), ox))
  r_f3 <- stats::cor(f3, ox)
  tibble::tibble(genotype = gt,
                 n_foxo3_nonmito = length(intersect(
                   ens_set(setdiff(gmt[["TFT_FOXO3_CHUNG"]], mito_universe)), rownames(Z))),
                 r_foxo3_oxphos = r_f3,
                 pct_vs_all_genes = 100 * mean(allr < r_f3, na.rm = TRUE),
                 r_foxo3_carrier = stats::cor(f3, ap))
})

foxo_verdict <- tibble::tibble(
  claim = c(
    "FOXO3 is the strongest READABLE TF signal in the wild-type timeline",
    "The rise is an OXPHOS-content artifact",
    "FOXO3 is the DRIVER of the OXPHOS decline",
    "The rise is the arrest arm",
    "The adult wild-type gland arms PUMA through FOXO3",
    "Myc BLOCKS the developmental FOXO3 rise",
    "Myc suppresses the FOXO3 programme"),
  verdict = c(
    "YES -- top of the 69 context-free lanes, passes the lane rule",
    "NO -- zero OXPHOS subunits in the set",
    "NO -- reach is 0 of 89 declining subunits (CORE_MITO reaches 59)",
    "NO -- see foxo_arm_summary; the atrophy/turnover arm carries it",
    "NO -- Bbc3 is flat in the wild-type timeline",
    "NO -- the programme rises in BOTH genotypes; interaction is ns",
    "YES -- significant at both ages, with FOXO1 null as the control"),
  read_from = c("tf_contextfree", "foxo_reach", "foxo_reach",
                "foxo_arm_summary", "foxo_targets", "foxo_lanes", "foxo_lanes"))

# =============================================================================
# SAVE
# =============================================================================
notes <- c(
  "47: is the wild-type 6W->12W OXPHOS decline a PGC1a-axis change, and what",
  "licenses the PGC1a/NRF1 rescue.",
  "",
  "SCOPE. BATCH = TIMEPOINT: every between-age number is confounded with cohort.",
  "What survives is the ARM-SELECTIVITY (PART B -- a shared batch effect cannot",
  "move one arm of a regulon to the 0.05th percentile and the rest to the 99.55th)",
  "and PART E (timepoint means removed entirely). What does NOT survive is the",
  "MAGNITUDE of the decline. Ranking plus negatives, not confirmatory inference.",
  "",
  "PART A. Arm-selectivity on one ruler. chain_structural vs chain_assembly is the",
  "sharp control: the assembly factors build the very complexes whose subunits fell.",
  "",
  "PART B. The decisive test. Read pct_ox against pct_rest WITHIN a regulon. An",
  "activity change acts on the regulon; it cannot act on a functional subset.",
  "",
  "PART C. The layer cannot rank factors here. Read LANE_RULE before quoting any",
  "TF lane; `_MITO` lanes are never TF activity.",
  "",
  "PART D. wt_pct is the expression-matched percentile -- a lone LFC says nothing",
  "at this n. Ppargc1a's baseMean is the number that decides how the rescue is",
  "described in the paper.",
  "",
  "PART E. Percentile against all expressed genes is the only null that means",
  "anything here (scripts 33/35: at n = 12 everything correlates at 0.6-0.8).",
  "",
  "PART F. NOT identified. The adjusted betas are a bound. Cell-resolved",
  "measurement is the design that settles it.",
  "",
  "PART H. The rescue's justification in numbers: reach, convergence, ownership,",
  "calibration, overshoot. H3 is the one a reviewer will test -- if the death",
  "effectors were in the PGC1a regulon the experiment would be circular. They are",
  "not; they are in MYC's.")

out <- list(
  armmap            = armmap,
  arm_class_summary = arm_class_summary,
  mt_arm            = mt_arm,
  split             = split,
  split_verdict     = split_verdict,
  tf_layer          = tf_layer,
  tf_context        = tf_context,
  tf_demo           = tf_demo,
  tf_all_lanes      = tf_all_lanes,
  lane_rule_hits    = lane_rule_hits,
  tf_contextfree    = tf_contextfree,
  tf_roster_lanes   = tf_roster_lanes,
  tf_lane_summary_24 = tf_lane_summary_24,
  roster            = roster,
  pgc1_family       = pgc1_family,
  couplings         = couplings,
  couplings_meta    = couplings_meta,
  carrier           = carrier,
  verdict           = verdict,
  foxo_reach        = foxo_reach,
  foxo_lanes        = foxo_lanes,
  foxo_targets      = foxo_targets,
  foxo_arm_summary  = foxo_arm_summary,
  foxo_separability = foxo_separability,
  foxo_verdict      = foxo_verdict,
  rescue            = rescue,
  defs = list(
    n_set_draws = NSET, n_bins = NBIN, min_arm_n = MIN_N,
    controls_observed = ctrl, controls_reference = ctrl_ref,
    control_nes_observed = nes_ox, control_nes_reference = NES_OX_REF,
    arm_class_rule = "see arm_class(); stored so the assignment is auditable",
    lane_rule = LANE_RULE, lane_layer_iqr = tf_layer_iqr,
    verdict_rule = VERDICT_RULE,
    regulons = REGULONS, roster_genes = ROSTER, roster_tf_lanes = ROSTER_TF,
    watch_arms = WATCH),
  analysis_date = Sys.Date(),
  notes = notes)

saveRDS(out, here::here("results", "biogenesis_axis_developmental.rds"))
message("47: wrote results/biogenesis_axis_developmental.rds")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  res <- readRDS(here::here("results", "biogenesis_axis_developmental.rds"))

  ## PART 0 -- the controls. Must reproduce script 44 exactly, or nothing below is
  ## readable. controls_observed vs controls_reference.
  res$defs$controls_observed
  res$defs$controls_reference

  ## PART A -- IS THE DECLINE ARM-SELECTIVE? Read class first, then the arms. The
  ## claim holds if chain_structural is at the bottom while chain_assembly and
  ## organelle_build sit near the middle of their nulls (pct ~50).
  res$arm_class_summary |> print()

  ## The arms themselves, most negative first. wt_null_pct is the number to read --
  ## a content value without its null is not interpretable.
  res$armmap |>
    dplyr::select(arm, class, n_genes, c_wt_time, wt_null_pct, NES) |>
    print(n = 30)

  ## The assembly factors are the sharp control: same complexes, opposite result.
  res$armmap |> dplyr::filter(class %in% c("chain_structural", "chain_assembly")) |>
    dplyr::select(arm, class, n_genes, c_wt_time, wt_null_pct) |> print(n = 25)

  ## The mtDNA arm, reported and NOT used. Read the note with it.
  res$mt_arm |> print()

  ## PART B -- THE DECISIVE TEST. pct_ox should be at the floor and pct_rest above
  ## the middle, in EVERY regulon. If it is, the split belongs to the OXPHOS genes
  ## and not to any factor, and no axis-activity change can explain it.
  res$split_verdict |> print()

  ## The same thing gene-count by gene-count, if the summary looks too tidy.
  res$split |> print(n = 25)

  ## PART C -- CAN THE LAYER RANK FACTORS? adj_r2_context vs adj_r2_tf. If context
  ## wins, the lanes report a cell state and the question cannot be answered here.
  res$tf_layer |> print()

  ## Between-context spread against within-context spread. ctx_median_range is the
  ## between; iqr_NES is the within.
  res$tf_context |> print(n = 20)

  ## THE DEMONSTRATION. Inside AP_LE, where do a general transcription factor
  ## (GTF3A) and a MICOS structural protein (CHCHD3) rank relative to ESRRA? If
  ## they outrank it, the lane is measuring gene content, not activity.
  res$tf_demo |> head(12) |> print()
  res$tf_demo |> dplyr::filter(tf %in% c("ESRRA", "GABPA", "NRF1", "GTF3A", "CHCHD3")) |>
    print()

  ## Does the reading rule actually discriminate? If it passed half the layer it
  ## would not be a rule. Read frac by lane kind.
  res$lane_rule_hits |> print()

  ## The roster lanes against the reading rule. readable_as_activity should be
  ## FALSE almost everywhere -- that is the result, not a failure.
  res$tf_roster_lanes |> dplyr::filter(!is_mito_lane) |> print(n = 40)

  ## THE ONLY SUB-LAYER WHERE A TF STATEMENT IS POSSIBLE: the generic regulons,
  ## which carry no Gray context. Top and bottom. This is where FOXO3 sits, and
  ## where the biogenesis factors visibly do NOT fall.
  res$tf_contextfree |> head(10) |> print()
  res$tf_contextfree |> dplyr::arrange(NES) |> head(10) |> print()

  ## Caution to carry with the bottom of that list: the strongest context-free
  ## fallers are the NF-kB / AP-1 immediate-early lanes, which are the
  ## dissociation-stress axis (CLAUDE.md, van den Brink) -- prep signal, not
  ## necessarily developmental biology.
  res$tf_contextfree |> dplyr::filter(grepl("RELA|NFKB1|REL_|FOS|EGR1|JUN", pathway)) |>
    print()

  ## script 24's lane medians read against the LAYER AMBIENT rather than zero. The
  ## PGC1a axis being ABOVE ambient is the point.
  res$tf_lane_summary_24 |> print(n = 20)

  ## PART D -- the factors. Nothing should clear the expression-matched percentile
  ## except Pprc1. wt_pct outside [5, 95] is the bar.
  res$roster |>
    dplyr::select(sym, baseMean, wt_lfc, wt_padj, wt_pct, mycpos_lfc, int_padj) |>
    print(n = 60)

  ## THE NUMBER FOR THE PAPER: Ppargc1a's baseMean against Pprc1's.
  res$pgc1_family |> print()

  ## PART E -- what does OXPHOS track among animals of the same age? Sets first.
  res$couplings |> dplyr::filter(kind == "set") |> print(n = 15)

  ## Then genes, by percentile. The top and the bottom are both informative: the
  ## growth machinery at the top, the quiescence arm at the bottom, and the whole
  ## biogenesis axis in the middle.
  res$couplings |> dplyr::filter(kind == "gene", genotype == "neg") |>
    dplyr::arrange(dplyr::desc(r)) |> print(n = 20)
  res$couplings |> dplyr::filter(kind == "gene", genotype == "neg") |>
    dplyr::arrange(r) |> print(n = 12)

  ## PART F -- the carrier. Read cor_ox_apLE_nonmito and vif_apLE_nonmito BEFORE
  ## the betas: at that collinearity the adjusted beta is a bound, not a result.
  res$carrier |> print()
  cat(res$carrier$limit_note, "\n")

  ## PART G -- the verdict table, on a rule fixed before the tests.
  res$verdict |> print()

  ## PART H -- what licenses the rescue.
  ## H1 REACH: does the axis reach the genes that fell?
  res$rescue$reach |> print()

  ## H2 CONVERGENCE: ESRRA and NRF1 share little, and the shared part is where the
  ## respiratory subunits concentrate. That is why using both is a control.
  res$rescue$convergence |> print()
  res$rescue$convergence_test |> print()

  ## H3 OWNERSHIP: the death effectors must NOT be in the PGC1a regulon, or the
  ## rescue is circular. CORE_MITO fold should sit at or below 1 while MYC_MITO is
  ## above it. Read the `genes` column -- Cycs is the one respiratory/death gene.
  res$rescue$ownership |>
    dplyr::filter(apoptosis_set == "MITOCARTA_APOPTOSIS_PRO") |>
    dplyr::select(programme, overlap, expected, fold, p_enrich, genes) |> print()

  ## H4 CALIBRATION: the target effect size the mouse sets for the rescue.
  res$rescue$calibration |> print()

  ## H5 OVERSHOOT: the arms that did not fall but are inside the regulon. Pre-specify
  ## that these will rise; the claim then stays "raising respiratory capacity",
  ## not "reverting the developmental change".
  res$rescue$overshoot |> print()

  ## PART I -- THE FOXO3 ARM. Read the verdict table first; every row names the
  ## object it is read from, so nothing here has to be taken on trust.
  res$foxo_verdict |> print()

  ## I2 is the one that decides driver against marker. covers_ox = 0 means FOXO3
  ## does not touch a single gene that fell, so it cannot be the proximal cause --
  ## the same test that disqualified the PGC1a axis, applied consistently.
  res$foxo_reach |> dplyr::select(set, n_genes, n_mitocarta, covers_ox,
                                  frac_of_declining) |> print()
  cat(res$foxo_reach$mito_members[res$foxo_reach$set == "TFT_FOXO3_CHUNG"], "\n")

  ## I3 -- the lane across all five rankings. The developmental rise is in BOTH
  ## genotypes and the interaction is ns; the Myc SUPPRESSION is significant at
  ## both ages. FOXO1 should be null throughout (the specificity control).
  res$foxo_lanes |> print(n = 15)

  ## WARNING before quoting the two genotype rows as "it does not attenuate":
  ## NES is scale-free and CANNOT see an amplitude fade (the four-list result,
  ## PANELS.md under Fig. 1H). Equal NES at both ages is not equal effect size.

  ## Which arm of the programme moved? Named before the values were read.
  res$foxo_arm_summary |> print()
  res$foxo_targets |> dplyr::select(arm, sym, wt_lfc, wt_padj, wt_pct, mycpos_lfc) |>
    print(n = 25)

  ## The story this section is most likely to be misread as supporting, tested:
  ## Bbc3 in the wild-type timeline. Flat means the adult gland does NOT arm PUMA.
  res$foxo_targets |> dplyr::filter(sym == "Bbc3") |> print()

  ## I4 -- separability. pct_vs_all_genes is the number: mid-range means the
  ## FOXO3 axis is not a standout inverse of OXPHOS, and r_foxo3_carrier says it
  ## is not separable from the cell state either.
  res$foxo_separability |> print()

  cat(res$notes, sep = "\n")
}
