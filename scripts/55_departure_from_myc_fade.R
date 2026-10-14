# =============================================================================
# 55_departure_from_myc_fade.R
# -----------------------------------------------------------------------------
# HOW FAR EACH GENE'S 12-WEEK MYC EFFECT DEPARTS FROM MYC'S OVERALL FADE, FOR
# EVERY EXPRESSED GENE -- the record behind the Fig. 1 apoptosis panels:
# ED Fig. 3b (the BCL-2 family and downstream genes, ranked) and ED Fig. 3c (the
# PUMA regulators), and the numbers their legends and the text quote.
#
# WHY IT EXISTS (author, 2026-10-14). Script 44 PART C computed the departure
# only for its REPORTING SET -- the genes MYC moved at 6W (|LFC6| >= 0.2,
# baseMean >= 20) -- because its question was "which MYC programmes collapse",
# and you cannot lose an effect you never had. The apoptosis panels need the same
# number for the genes MYC did NOT move at 6W as well: for them the fade predicts
# about zero at 12W, so their departure is simply a 12-week effect, and it still
# moves the balance (Mcl1, Pmaip1, most of the p53 arm). The 2026-10-13 proposal
# panels computed it inside the panel script; this script makes it a record the
# panels read and assert against.
#
# THE STATISTIC is script 44's, with script 44's fitted rate, left UNSIGNED:
#   departure   d  = LFC12 - RATE x LFC6      RATE = 0.4872, fitted by script 44
#                                             through the origin over the genes MYC
#                                             moves (the rate of record for genes)
#   its SE         = sqrt(SE12^2 + RATE^2 x SE6^2)   the 6W and 12W contrasts come
#                                                     from disjoint animals
#   z              = d / SE
#   p_norm         = 2 x pnorm(-|z|): nominal, two-sided, a Wald-type test of d = 0
# Negative: the 12W MYC effect is LOWER than the fade predicts. Script 44 signs it
# by the direction of the 6W effect instead (z_resid = sign(LFC6) x z); for its
# reporting set the two are asserted identical up to that sign.
#
# WHAT IT IS NOT: THE INTERACTION. The interaction (LFC12 - LFC6) is the whole
# change in MYC's effect between the ages. It mixes the common fade,
# -(1 - RATE) x LFC6, with the departure; both shares are stored, and they sum to
# the interaction (asserted). Htra2 and Bbc3 have the same interaction (-0.55,
# -0.54) for opposite reasons: Htra2's is almost all fade, Bbc3's is mostly
# departure (2026-10-14).
#
# BATCH = TIMEPOINT does not bite: d is built from two genotype contrasts, each
# within one extraction batch.
#
# THE REFERENCE BANDS, what an ordinary gene does, kept apart:
#   responsive   |z| over script 44's reporting set (8,774 genes MYC moved at 6W;
#                the set is rebuilt here and asserted identical)
#   unmoved      |z| over the expressed genes MYC did not move (|LFC6| < 0.2)
# Each gene gets the percentile of its |z| within its own reference, and is
# "outside" when |z| exceeds that reference's 95th percentile.
#
# THE ROSTERS, declared here before any number is read:
#   BCL-2 family    the anti-apoptotic members, the effectors and the BH3-only
#                   proteins, with a ROLE for the Fig. 1g heatmap. Activators are
#                   BID, BIM and PUMA (the common modern grouping; PUMA's place is
#                   debated -- the original sensitiser/activator split put it with
#                   the sensitisers). Below the expression floor and dropped: Hrk,
#                   Bcl2a1d, Bcl2l10 (none is read in these libraries).
#   downstream      script 42's execution and IAP genes, its own arm labels
#   PUMA regulators script 42's $exclusions$p53_axis (the p53 arm, 9) and
#                   $exclusions$puma_inputs (12: the ISR, E2f1, Trp73 and the FOXO
#                   family), plus Bbc3 -- the old Fig. S2D (alt) roster. Trp73 and
#                   Cdkn2a fall below the floor (and Cdkn2a cannot separate
#                   p19Arf from p16Ink4a anyway; CLAUDE.md).
# FUNCTION SIGN for the BCL-2 family and downstream: + = toward survival, a
# pro-apoptotic transcript lower, or an anti-apoptotic one higher, than the fade
# predicts. The PUMA regulators get none (Mdm2 restrains p53; the others induce
# PUMA), so their panel draws the unsigned departure.
#
# PART 0 (TEXT) declares the numbers the panels and the draft text quote, and each
# is asserted where it is computed. A stop there means the text and the analysis
# have drifted apart.
#
# Reads:  results/interaction_results.rds      (script 03)
#         results/collapse_module_ownership.rds (script 44: $defs, $collapse_genes)
#         results/priming_arm_teb.rds          (script 42: $machinery, $exclusions)
#         results/combined_df_annotated.rds    (symbols, as script 44 reads them)
# Writes: results/departure_from_myc_fade.rds
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

BM_FLOOR   <- 20     # scripts 42 and 44's expression floor
LFC6_FLOOR <- 0.2    # script 44: "MYC moved it at 6W"
PADJ6      <- 0.1    # script 44's ranking set, used here only to re-derive the rate
ALPHA      <- 0.05   # nominal threshold for the panel labels (raw p)
BAND_Q     <- 0.95   # a gene is outside when |z| exceeds this quantile of its reference

# =============================================================================
# PART 0 (TEXT): the numbers quoted, declared before anything is computed
# =============================================================================
TEXT <- list(
  rate3              = 0.487,                    # the fade rate, to three decimals
  band_responsive    = 1.73,                     # 95th percentile of |z|, two decimals
  bcl2_outside       = "Bbc3",                   # the only BCL-2 family member outside
  downstream_outside = "Apaf1",                  # the only downstream gene outside
  regulators_p05     = c("Bbc3", "Eda2r", "Foxo3", "Zmat3"),  # raw interaction p < 0.05
  regulators_bh_none = TRUE)                     # none survives BH within the roster

# =============================================================================
# PART A: LOAD, AND REPRODUCE SCRIPT 44
# =============================================================================
message("55 PART A: load")
ir  <- readRDS(here::here("results", "interaction_results.rds"))
cmo <- readRDS(here::here("results", "collapse_module_ownership.rds"))
pa  <- readRDS(here::here("results", "priming_arm_teb.rds"))
cdf <- readRDS(here::here("results", "combined_df_annotated.rds"))

D   <- lapply(ir[c("myc_6W_raw", "myc_12W_raw", "interaction_raw")], as.data.frame)
U   <- rownames(D$myc_6W_raw)
stopifnot(identical(rownames(D$myc_12W_raw), U), identical(rownames(D$interaction_raw), U))
col <- function(k, v) stats::setNames(D[[k]][[v]], U)
m6  <- col("myc_6W_raw", "log2FoldChange");  se6  <- col("myc_6W_raw", "lfcSE")
m12 <- col("myc_12W_raw", "log2FoldChange"); se12 <- col("myc_12W_raw", "lfcSE")
p6a <- col("myc_6W_raw", "padj");            bm   <- col("myc_6W_raw", "baseMean")
mi  <- col("interaction_raw", "log2FoldChange")
pi_ <- col("interaction_raw", "pvalue");     pia  <- col("interaction_raw", "padj")

RATE <- cmo$defs$global_rate_fitted
cg   <- as.data.frame(cmo$collapse_genes)

# the rate, re-derived in closed form from script 44's own ranking set
keep_rep <- U[is.finite(m6) & is.finite(m12) & abs(m6) >= LFC6_FLOOR & bm >= BM_FLOOR]
keep6    <- keep_rep[!is.na(p6a[keep_rep]) & p6a[keep_rep] < PADJ6]
stopifnot(abs(RATE - sum(m6[keep6] * m12[keep6]) / sum(m6[keep6]^2)) < 1e-12,
          round(RATE, 3) == TEXT$rate3,
          # the reporting set is script 44's, gene for gene
          setequal(keep_rep, cg$ens), nrow(cg) == 8774L)

# the interaction IS LFC12 - LFC6 (the saturated model's identity)
fin <- is.finite(m6) & is.finite(m12) & is.finite(mi)
stopifnot(max(abs(mi[fin] - (m12[fin] - m6[fin]))) < 1e-6)

# =============================================================================
# PART B: THE DEPARTURE, FOR EVERY EXPRESSED GENE
# =============================================================================
message("55 PART B: the departure for every expressed gene")
expressed <- U[is.finite(m6) & is.finite(m12) & is.finite(se6) & is.finite(se12) & bm >= BM_FLOOR]
sym <- cdf$mgi_symbol[match(expressed, cdf$gene)]

genes <- tibble::tibble(
  ens = expressed, gene = sym, baseMean = unname(bm[expressed]),
  lfc6 = unname(m6[expressed]), se6 = unname(se6[expressed]), padj6 = unname(p6a[expressed]),
  lfc12 = unname(m12[expressed]), se12 = unname(se12[expressed]),
  interaction = unname(mi[expressed]), int_p = unname(pi_[expressed]),
  int_padj_ihw = unname(pia[expressed])) |>
  dplyr::mutate(
    expected12  = RATE * lfc6,
    fade_share  = -(1 - RATE) * lfc6,            # what the common fade alone does
    departure   = lfc12 - RATE * lfc6,           # what changes beyond it
    se_dep      = sqrt(se12^2 + RATE^2 * se6^2),
    z           = departure / se_dep,
    p_norm      = 2 * stats::pnorm(-abs(z)),
    reference   = ifelse(abs(lfc6) >= LFC6_FLOOR, "responsive", "unmoved"))

# the decomposition is exact: fade share + departure = interaction
stopifnot(max(abs(genes$fade_share + genes$departure - genes$interaction), na.rm = TRUE) < 1e-6)

# script 44's statistic, reproduced: z_resid = sign(LFC6) x z on its reporting set
k44 <- match(cg$ens, genes$ens)
stopifnot(!anyNA(k44),
          max(abs(sign(genes$lfc6[k44]) * genes$z[k44] - cg$z_resid)) < 1e-9,
          all(genes$reference[k44] == "responsive"),
          sum(genes$reference == "responsive") == nrow(cg))

# =============================================================================
# PART C: THE REFERENCE BANDS
# =============================================================================
message("55 PART C: the reference bands")
band <- genes |>
  dplyr::group_by(reference) |>
  dplyr::summarise(n = dplyr::n(), q_abs_z = unname(stats::quantile(abs(z), BAND_Q)),
                   q025 = unname(stats::quantile(z, 0.025)), q975 = unname(stats::quantile(z, 0.975)),
                   .groups = "drop")
BAND <- stats::setNames(band$q_abs_z, band$reference)
stopifnot(round(BAND[["responsive"]], 2) == TEXT$band_responsive)

ref_abs <- split(abs(genes$z), genes$reference)
genes <- genes |>
  dplyr::mutate(
    pct_abs_z = vapply(seq_len(dplyr::n()), function(i)
      100 * mean(ref_abs[[reference[i]]] < abs(z[i])), numeric(1)),
    outside_band = abs(z) > BAND[reference])

# =============================================================================
# PART D: THE ROSTERS
# =============================================================================
message("55 PART D: the rosters")
ROSTER_BCL2 <- tibble::tribble(
  ~gene,      ~role,
  "Bax",      "effector",      "Bak1",    "effector",      "Bok",     "effector",
  "Bid",      "activator",     "Bcl2l11", "activator",     "Bbc3",    "activator",
  "Bad",      "sensitiser",    "Bik",     "sensitiser",    "Bmf",     "sensitiser",
  "Pmaip1",   "sensitiser",    "Hrk",     "sensitiser",
  "Bcl2",     "anti-apoptotic", "Bcl2l1", "anti-apoptotic", "Bcl2l2", "anti-apoptotic",
  "Mcl1",     "anti-apoptotic", "Bcl2a1a", "anti-apoptotic", "Bcl2a1b", "anti-apoptotic",
  "Bcl2a1d",  "anti-apoptotic", "Bcl2l10", "anti-apoptotic")

m42  <- as.data.frame(pa$machinery)
down <- m42[grepl("execution|apoptosome|IAP", m42$arm), c("gene", "arm")]
stopifnot(setequal(down$gene, c("Htra2", "Diablo", "Cycs", "Apaf1", "Casp9", "Casp3",
                                "Xiap", "Birc5")))
# Anti-apoptotic downstream = script 42's "brake (IAP...)" arms only. NOT a match on
# "IAP": Diablo's arm is "execution (IAP antagonist)" -- SMAC, pro-apoptotic. The
# first version matched "IAP" and signed Diablo as anti-apoptotic (caught in the
# author's first run, 2026-10-14).
ROSTER_DOWN <- tibble::tibble(gene = down$gene, role = ifelse(grepl("^brake", down$arm),
                                                              "anti-apoptotic", down$arm))
stopifnot(setequal(ROSTER_DOWN$gene[ROSTER_DOWN$role == "anti-apoptotic"], c("Xiap", "Birc5")),
          ROSTER_DOWN$role[ROSTER_DOWN$gene == "Diablo"] == "execution (IAP antagonist)")

p53 <- as.data.frame(pa$exclusions$p53_axis)$gene
pin <- as.data.frame(pa$exclusions$puma_inputs)$gene
stopifnot(length(p53) == 9L, length(pin) == 12L, "Foxo3" %in% pin, "Foxo1" %in% pin)
ROSTER_REG <- tibble::tibble(gene = c("Bbc3", p53, pin),
                             role = c("PUMA", rep("p53 arm", length(p53)),
                                      rep("PUMA inducer", length(pin))))

ens_of <- function(g) {
  e <- tryCatch(recon_to_ensembl(g, U), error = function(e) NA_character_)
  if (length(e) != 1L) NA_character_ else unname(e)
}
build_roster <- function(r, name, signed) {
  r$roster <- name
  r$ens <- vapply(r$gene, ens_of, character(1))
  r$dropped <- ifelse(is.na(r$ens), "not in the count matrix",
                      ifelse(!r$ens %in% genes$ens, "below the expression floor", NA_character_))
  r <- dplyr::left_join(r, dplyr::select(genes, -gene), by = "ens")
  r$fn_sign <- if (signed) ifelse(r$role == "anti-apoptotic", 1, -1) else NA_real_
  r$z_function <- r$fn_sign * r$z                 # + = toward survival
  r
}
rosters <- dplyr::bind_rows(
  build_roster(ROSTER_BCL2, "BCL-2 family", TRUE),
  build_roster(ROSTER_DOWN, "downstream", TRUE),
  build_roster(ROSTER_REG,  "PUMA regulators", FALSE))

# within-roster multiplicity, on the genes that are read
rosters <- rosters |>
  dplyr::group_by(roster) |>
  dplyr::mutate(int_p_bh  = ifelse(is.na(dropped), stats::p.adjust(ifelse(is.na(dropped), int_p, NA), "BH"), NA),
                p_norm_bh = ifelse(is.na(dropped), stats::p.adjust(ifelse(is.na(dropped), p_norm, NA), "BH"), NA)) |>
  dplyr::ungroup()

read <- function(name) rosters[rosters$roster == name & is.na(rosters$dropped), ]
fam <- read("BCL-2 family"); dn <- read("downstream"); rg <- read("PUMA regulators")

# --- the claims, asserted where computed ------------------------------------------
stopifnot(
  identical(sort(rosters$gene[rosters$roster == "BCL-2 family" & !is.na(rosters$dropped)]),
            sort(c("Hrk", "Bcl2a1d", "Bcl2l10"))),
  identical(fam$gene[fam$outside_band], TEXT$bcl2_outside),
  fam$gene[which.max(fam$z_function)] == "Bbc3",
  identical(dn$gene[dn$outside_band], TEXT$downstream_outside),
  identical(sort(rg$gene[rg$int_p < ALPHA]), sort(TEXT$regulators_p05)),
  identical(sort(rg$gene[rg$p_norm < ALPHA]), sort(TEXT$regulators_p05)),
  isTRUE(all(rg$int_p_bh >= ALPHA)) == TEXT$regulators_bh_none,
  identical(sort(rosters$gene[rosters$roster == "PUMA regulators" & !is.na(rosters$dropped)]),
            sort(c("Trp73", "Cdkn2a"))))

# =============================================================================
# PART E: NOTES, AND SAVE
# =============================================================================
g1 <- function(x) genes[genes$ens == ens_of(x), ]
dec <- function(x) { r <- g1(x); sprintf("%-7s interaction %+.3f = fade share %+.3f + departure %+.3f (z %+.2f)",
                                          x, r$interaction, r$fade_share, r$departure, r$z) }
notes <- c(
  sprintf("THE FADE RATE: %.4f (script 44, fitted over the genes MYC moves; re-derived here).", RATE),
  sprintf("THE BANDS: 95th percentile of |z| = %.2f over %d responsive genes, %.2f over %d unmoved ones.",
          BAND[["responsive"]], band$n[band$reference == "responsive"],
          BAND[["unmoved"]], band$n[band$reference == "unmoved"]),
  "THE DECOMPOSITION, for the genes argued over on 2026-10-13/14:",
  paste0("  ", vapply(c("Bbc3", "Htra2", "Bax", "Bcl2l1", "Bmf", "Apaf1"), dec, character(1))),
  sprintf("BCL-2 FAMILY: %s alone outside its band (z %+.2f toward survival).",
          TEXT$bcl2_outside, fam$z_function[fam$gene == "Bbc3"]),
  sprintf("DOWNSTREAM: %s alone outside (z %+.2f toward survival): its 6W repression deepened at 12W instead of fading.",
          TEXT$downstream_outside, dn$z_function[dn$gene == "Apaf1"]),
  sprintf("PUMA REGULATORS (%d read): raw interaction p < %.2f for %s; none survives BH within the roster.",
          nrow(rg), ALPHA, paste(sort(TEXT$regulators_p05), collapse = ", ")),
  "  Bbc3 and Foxo3 depart BELOW the fade; the two p53 targets Zmat3 and Eda2r depart ABOVE it,",
  "  through a fall in the WT gland with the MYC+ gland flat. Label all four (author, 2026-10-14).",
  "SCOPE: n = 6 per group; p_norm is nominal; the bands describe what an ordinary gene does,",
  "they are not a test. Bbc3's licence is pre-specification; Foxo3, Zmat3 and Eda2r are exploratory.")
cat(paste(notes, collapse = "\n"), "\n")

out <- list(
  genes = genes, band = band, rosters = rosters,
  defs = list(rate = RATE, bm_floor = BM_FLOOR, lfc6_floor = LFC6_FLOOR, alpha = ALPHA,
              band_q = BAND_Q, roster_bcl2 = ROSTER_BCL2, roster_down = ROSTER_DOWN,
              roster_regulators = ROSTER_REG, text = TEXT),
  analysis_date = Sys.Date(), notes = notes)

saveRDS(out, here::here("results", "departure_from_myc_fade.rds"))
message("55: wrote results/departure_from_myc_fade.rds")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  res <- readRDS(here::here("results", "departure_from_myc_fade.rds"))

  ## the bands, by reference
  res$band

  ## the three rosters, ranked (+ = toward survival for the first two)
  res$rosters |>
    dplyr::filter(is.na(dropped), roster != "PUMA regulators") |>
    dplyr::arrange(roster, dplyr::desc(z_function)) |>
    dplyr::select(roster, gene, role, lfc6, lfc12, z_function, pct_abs_z, outside_band) |>
    print(n = 30)
  res$rosters |>
    dplyr::filter(is.na(dropped), roster == "PUMA regulators") |>
    dplyr::arrange(z) |>
    dplyr::select(gene, role, lfc6, lfc12, z, p_norm, int_p, int_p_bh) |>
    print(n = 25)

  ## the decomposition for any gene: interaction = fade share + departure
  res$genes |> dplyr::filter(gene %in% c("Bbc3", "Htra2", "Bax", "Bcl2l1")) |>
    dplyr::select(gene, lfc6, lfc12, interaction, fade_share, departure, z) |>
    print()
}
