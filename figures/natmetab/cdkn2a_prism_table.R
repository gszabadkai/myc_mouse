# =============================================================================
# cdkn2a_prism_table.R -- Cdkn2a per animal, four arms, for a Prism column plot
# -----------------------------------------------------------------------------
# Requested 2026-10-29: Cdkn2a CPM per tumour, broken down to individual animals,
# over NTEV, NTPgc1a, BclxL and Pgc1a_BclxL. It draws nothing.
#
# WHY THIS READS THE MATRIX AND NOT THE SAVED OBJECT. Script 52 carries Cdkn2a in
#   `p53_panel`, but that is a CONTRAST -- a log2 difference with a bootstrap
#   interval, on the three within-series pairs. No per-animal CPM is saved for it
#   anywhere. So PART 0 of script 52 is reproduced here exactly: the same file,
#   the same longest-token-first group parse, the same
#   `t(t(CNT)/colSums(CNT))*1e6`, and the same assertions, so a drift in any of
#   them stops the script rather than changing a number quietly. Nothing is
#   re-fitted and no saved interval is recomputed; script 52's own control is
#   re-run below and must still reproduce to within 2%.
#
# Cdkn2a IS NOT p19ARF, AND A COLUMN OF IT MAY NOT BE LABELLED p19ARF.
#   Mouse Cdkn2a encodes p16INK4a and p19ARF from alternative first exons sharing
#   exons 2 and 3, and a gene-level count pools them (script 52 PART D). The
#   western and the exon-1b qPCR separate what these counts cannot. A flat Cdkn2a
#   here does not exclude p19ARF loss, and a fall here is not a fall in p19ARF.
#   Axis title: "Cdkn2a (CPM)", italic gene name, nothing about ARF.
#
# THE FOUR ARMS SPAN TWO SERIES, AND THE COLUMN ORDER KEEPS THE PAIRING.
#   NTEV/NTPgc1a are CRISPR non-targeting clones; BclxL/Pgc1a_BclxL are the
#   polyclonal vector pool. Script 52's contrast rule: clean comparisons are
#   WITHIN series (NTPgc1a vs NTEV, Pgc1a_BclxL vs BclxL); across series they
#   carry a derivation difference and are DESCRIPTIVE ONLY. Four columns side by
#   side invite the across-series read, so the order is control-then-PGC1a twice,
#   and any bracket Prism draws should stay inside a series.
#   The p19ARF arms (KOEV, KOPgc1a) are deliberately absent: the locus is
#   CRISPR-targeted there, so a count is not a read-out of p53 output. That is
#   the same exclusion the specificity strip makes (edfig5a_bbc3_routes.R).
#
# CPM HAS A DENOMINATOR, AND IN THIS COHORT IT MOVES. Script 51 showed that an
#   arm which nearly doubles its mitochondrial share deflates every other gene's
#   CPM without any gene changing, and Pgc1a_BclxL is that arm: MitoCarta is
#   14.7% of its library against 7.8-8.6% elsewhere. A ~7-point shift in the
#   denominator takes ~7% off every non-mito gene's CPM in that arm alone.
#   So the median-of-ratios normalised count is written BESIDE the CPM, from
#   script 52's own size factors. Read the two together before reading a
#   difference between the vector columns; they are printed side by side below.
#
# PER-ANIMAL VALUES, SO PRISM'S TEST IS PRISM'S. These are raw values; anything
#   Prism computes from them is Prism's own statistic and the legend must cite
#   that test, not script 52's interval. Script 52's within-series contrast is
#   printed below for comparison only -- it is a bootstrap of a difference of
#   medians on log2 normalised counts and Prism cannot reproduce it. The split is
#   the one `fig3i_prism_tables.R` sets out; read that header.
#
# THE FIGURE NUMBER IS PROVISIONAL, as for EDFig5a/b (spec section 4, item 8).
#
# Reads : data/orthotopic_series/salmon.merged.gene_counts.tsv
#         results/orthotopic_escape_series.rds   (script 52, for its saved contrast)
# Writes: outputs/natmetab/EDFig5/EDFig5_prism_cdkn2a_per_animal.csv           (CPM)
#         outputs/natmetab/EDFig5/EDFig5_prism_cdkn2a_normcount_per_animal.csv (control)
#         outputs/natmetab/EDFig5/EDFig5_cdkn2a_per_animal_long.csv            (both, tidy)
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

GENE  <- "Cdkn2a"
ARMS  <- c("NTEV", "NTPgc1a", "BclxL", "Pgc1a_BclxL")       # control, PGC1a, twice
LABEL <- c(NTEV        = "NT, EV",
           NTPgc1a     = "NT, PGC-1a",
           BclxL       = "Bcl-xL",
           Pgc1a_BclxL = "PGC-1a + Bcl-xL")
SERIES <- c(NTEV = "CRISPR_NT", NTPgc1a = "CRISPR_NT",
            BclxL = "vector",   Pgc1a_BclxL = "vector")

# =============================================================================
# PART 0 -- script 52's load, parse and CPM, reproduced and asserted
# =============================================================================
# ALLGRP, SET and REF are copied VERBATIM from script 52 lines 103-121. The order
# of ALLGRP is load-bearing: REF is positional, so a reordered level set silently
# compares the wrong groups.
ALLGRP <- c("EV", "BclxL", "Pgc1a_BclxL", "NTEV", "KOEV", "NTPgc1a", "KOPgc1a",
            "MYAZ", "FaMY", "BoMY")
SET    <- c("EV", "BclxL", "Pgc1a_BclxL", "NTEV", "NTPgc1a", "KOEV", "KOPgc1a")
REF <- data.frame(
  group    = ALLGRP,
  n        = c(7L, 7L, 7L, 7L, 7L, 7L, 7L, 6L, 6L, 6L),
  Bcl2l1   = c(71, 326, 193, 82, 62, 79, 79, 59, 59, 51),
  Ppargc1a = c(0.0, 0.3, 302, 0.2, 0.1, 129, 0.2, 0.0, 0.0, 0.0),
  Myc      = c(1915, 2033, 1440, 2208, 2368, 2286, 2176, 1370, 1280, 1293),
  stringsAsFactors = FALSE)

raw <- data.table::fread(
  here::here("data", "orthotopic_series", "salmon.merged.gene_counts.tsv"),
  data.table = FALSE, check.names = FALSE)
stopifnot(identical(colnames(raw)[1:2], c("gene_id", "gene_name")))
CNT <- as.matrix(raw[, -(1:2)]); rownames(CNT) <- raw$gene_id
stopifnot(ncol(CNT) == 67L, all(grepl("^ENSMUSG", rownames(CNT))))

PAT <- c("NTEV", "KOEV", "NTPgc1a", "NTPgc", "KOPgc1a", "KOPgc", "BclxL",
         "BoMY1", "BoMY", "FaMY1", "FaMY", "MYAZ", "Pgc1a", "Pgc", "EV")
CANON <- c(NTEV = "NTEV", KOEV = "KOEV", NTPgc1a = "NTPgc1a", NTPgc = "NTPgc1a",
           KOPgc1a = "KOPgc1a", KOPgc = "KOPgc1a", BclxL = "BclxL",
           BoMY1 = "BoMY", BoMY = "BoMY", FaMY1 = "FaMY", FaMY = "FaMY",
           MYAZ = "MYAZ", Pgc1a = "Pgc1a_BclxL", Pgc = "Pgc1a_BclxL", EV = "EV")
core <- sub("^[A-H][0-9]-", "", colnames(CNT))
tok  <- vapply(core, function(s) {
  h <- PAT[vapply(PAT, function(p) grepl(paste0("^", p), s), TRUE)]
  if (length(h)) h[1] else NA_character_ }, "")
grp  <- factor(unname(CANON[tok]), levels = ALLGRP)
stopifnot(!anyNA(grp), sum(grp %in% SET) == 49L, all(table(grp)[SET] == 7L))

CPM <- t(t(CNT) / colSums(CNT)) * 1e6
cpm_of <- function(sym) {
  r <- which(raw$gene_name == sym)
  stopifnot(length(r) >= 1L)
  if (length(r) == 1L) CPM[r, ] else colSums(CPM[r, , drop = FALSE]) }

# script 52's own median-CPM control, so a changed file cannot pass silently
ctrl <- data.frame(
  group = ALLGRP, n = as.integer(table(grp)),
  Bcl2l1   = round(as.numeric(tapply(cpm_of("Bcl2l1"),   grp, stats::median)), 1),
  Ppargc1a = round(as.numeric(tapply(cpm_of("Ppargc1a"), grp, stats::median)), 1),
  Myc      = round(as.numeric(tapply(cpm_of("Myc"),      grp, stats::median)), 1))
stopifnot(identical(ctrl$n, REF$n))
dev <- max(abs(unlist(ctrl[, 3:5]) - unlist(REF[, 3:5])) / pmax(unlist(REF[, 3:5]), 1))
message(sprintf("PART 0: median-CPM control, max relative deviation %.3f%% (tol 2.0%%)",
                100 * dev))
stopifnot(dev < 0.02)

# =============================================================================
# PART 1 -- the second ruler: median-of-ratios, from script 52's size factors
# -----------------------------------------------------------------------------
# Size factors only. No model is fitted and no contrast is estimated here.
# =============================================================================
vk <- grp %in% SET
sm <- data.frame(sample = colnames(CNT)[vk],
                 group  = factor(as.character(grp[vk]), levels = SET))
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(CNT[, vk, drop = FALSE]), colData = sm, design = ~ group)
dds <- DESeq2::estimateSizeFactors(dds)
NCv <- DESeq2::counts(dds, normalized = TRUE)
nc_of <- function(sym) { r <- which(raw$gene_name == sym); stopifnot(length(r) >= 1L)
                         m <- NCv[rownames(CNT)[r], , drop = FALSE]
                         if (nrow(m) == 1L) as.numeric(m) else colSums(m) }

g_cpm <- cpm_of(GENE)[vk]
g_mor <- nc_of(GENE)
arm   <- as.character(sm$group)
stopifnot(length(g_cpm) == 49L, length(g_mor) == 49L, !anyNA(g_cpm), !anyNA(g_mor))

# =============================================================================
# PART 2 -- the Prism block: one column per arm, seven rows
# =============================================================================
wide <- function(v, digits) {
  m <- vapply(ARMS, function(a) round(v[arm == a], digits), numeric(7))
  colnames(m) <- sprintf("%s (%s)", LABEL[ARMS], ARMS)
  as.data.frame(m) }

cpm_block <- wide(g_cpm, 2)
mor_block <- wide(g_mor, 1)

# the long form, so the series stays attached to every value
long <- do.call(rbind, lapply(ARMS, function(a) data.frame(
  arm = a, label = unname(LABEL[a]), series = unname(SERIES[a]),
  animal = sm$sample[arm == a],
  cdkn2a_cpm = round(g_cpm[arm == a], 2),
  cdkn2a_normcount = round(g_mor[arm == a], 1),
  stringsAsFactors = FALSE)))

d5 <- here::here("outputs", "natmetab", "EDFig5")
dir.create(d5, recursive = TRUE, showWarnings = FALSE)
f <- file.path(d5, "EDFig5_prism_cdkn2a_per_animal.csv")
utils::write.csv(cpm_block, f, row.names = FALSE)
message("wrote ", f)
f2 <- file.path(d5, "EDFig5_prism_cdkn2a_normcount_per_animal.csv")
utils::write.csv(mor_block, f2, row.names = FALSE)
message("wrote ", f2)
f3 <- file.path(d5, "EDFig5_cdkn2a_per_animal_long.csv")
utils::write.csv(long, f3, row.names = FALSE)
message("wrote ", f3)

# =============================================================================
# PRINTED, so the table can be read without opening a file
# =============================================================================
op <- options(width = 130); on.exit(options(op))

cat("\n", strrep("=", 78), "\n",
    "1. Cdkn2a, CPM, one column per arm, seven animals -- PASTE THIS INTO PRISM\n",
    "   control then PGC-1a, twice; the two pairs are different series\n", sep = "")
print(cpm_block, row.names = FALSE)

cat("\n", strrep("=", 78), "\n",
    "2. the same animals as median-of-ratios normalised counts (the CPM control)\n",
    "   Pgc1a_BclxL carries a 14.7% MitoCarta share against 7.8-8.6% elsewhere,\n",
    "   which deflates every other gene's CPM in that arm by roughly 7%\n", sep = "")
print(mor_block, row.names = FALSE)

cat("\n", strrep("=", 78), "\n", "3. medians, both rulers\n", sep = "")
summ <- data.frame(
  arm = ARMS, label = unname(LABEL[ARMS]), series = unname(SERIES[ARMS]),
  median_cpm = round(vapply(ARMS, function(a) stats::median(g_cpm[arm == a]), 0), 1),
  min_cpm    = round(vapply(ARMS, function(a) min(g_cpm[arm == a]), 0), 1),
  max_cpm    = round(vapply(ARMS, function(a) max(g_cpm[arm == a]), 0), 1),
  median_normcount = round(vapply(ARMS, function(a) stats::median(g_mor[arm == a]), 0), 0),
  stringsAsFactors = FALSE)
print(summ, row.names = FALSE)

cat("\n", strrep("=", 78), "\n",
    "4. script 52's SAVED within-series contrast for Cdkn2a, for comparison only\n",
    "   log2 difference of normalised counts with a bootstrap interval; Prism\n",
    "   cannot reproduce it, so a Prism panel must cite Prism's own test\n", sep = "")
p52 <- readRDS(here::here("results", "orthotopic_escape_series.rds"))$p53_panel
print(as.data.frame(p52[p52$gene == GENE, ]), row.names = FALSE)

cat("\n   REMINDER: mouse Cdkn2a pools p16INK4a and p19ARF. This axis says\n",
    "   Cdkn2a (CPM); it does not say p19ARF. The p19ARF arms are not in the\n",
    "   table because the locus is CRISPR-targeted there.\n", sep = "")

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  ## does the CPM denominator explain any vector-arm difference? compare the
  ## two rulers' arm ratios against their own controls
  r_cpm <- c(NT = stats::median(g_cpm[arm == "NTPgc1a"]) / stats::median(g_cpm[arm == "NTEV"]),
             vec = stats::median(g_cpm[arm == "Pgc1a_BclxL"]) / stats::median(g_cpm[arm == "BclxL"]))
  r_mor <- c(NT = stats::median(g_mor[arm == "NTPgc1a"]) / stats::median(g_mor[arm == "NTEV"]),
             vec = stats::median(g_mor[arm == "Pgc1a_BclxL"]) / stats::median(g_mor[arm == "BclxL"]))
  rbind(CPM = r_cpm, normcount = r_mor)
  ## the whole p53 panel, saved, for context
  as.data.frame(p52)
}
