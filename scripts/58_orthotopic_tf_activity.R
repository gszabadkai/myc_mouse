# =============================================================================
# 58 -- ORTHOTOPIC TF ACTIVITY (ED Fig. 5d) AND THE PRE-SPECIFIED INVASION GENES
# -----------------------------------------------------------------------------
# Cohort: the 21 vector-series samples of data/orthotopic_series/, 7 per arm --
# `EV` (EV + EV), `BclxL` (Bcl-xL + EV), `Pgc1a_BclxL` (Bcl-xL + PGC-1a) -- taken
# from script 52's saved sample table. Read data/orthotopic_series/README.md and
# docs/2026-11-07_edfig5d_tf_activity.md before this file.
#
# Contrasts: Pgc1a_BclxL vs BclxL (the PGC-1a contrast on a shared Bcl-xL pool) and
# Pgc1a_BclxL vs EV (two manipulations apart), as on ED 5d.
#
# WHAT ED 5d NOW SHOWS (author, 2026-11-07). The TF charts are DROPPED: TF activity
# does not test the hypothesis -- a PGC-1a effect on invasion separate from its
# mitochondrial role -- the mRNA of the pre-specified genes does. ED 5d is the
# heatmap of PART 5's 35 genes (figures/natmetab/edfig5d_invasion_genes.R), drawn
# from PART 6b: the arm means of z, a star where Bcl-xL + PGC-1a differs from at
# least one control, gene names in bold where the change is in the predicted
# direction (rules in PART 6b). PARTs 2-4 (the network and TF activity) stay as the
# record of why the original TF charts were withdrawn; nothing of them is drawn.
#
# WHAT THIS REPLACES. ED 5d's two bar charts came from scripts/external/Pgc bclxl
# analysis.R. The author chose the corrected analysis (option (b), 2026-11-07) and
# ruled that the original need not be reproduced. Three corrections:
#   1. THE CONSTRUCT TRANSCRIPTS ARE LEFT OUT of the gene statistics. `Ppargc1a` is
#      a CollecTRI target of 50 TFs and `Bcl2l1` of 60, and neither is a regulated
#      target in these arms: one is the PGC-1a transgene (+10 to +12 log2), the
#      other the Bcl-xL construct. Left in, the transgene alone lowers Twist1's
#      score by about 0.6-0.7, because CollecTRI lists Ppargc1a as a target Twist1
#      represses (note section 4.4; PART 4 reports it again on this pipeline).
#      Both stay in the DESeq2 fit, where they are the construct check.
#   2. THE REGULONS ARE THE PINNED SNAPSHOT in data/collectri_mouse/ (OmniPath's
#      mouse CollecTRI, 2026-11-07), not a live get_collectri() call, which fails
#      here (OmnipathR 4.0.0). The 1,077 nodes OmniPath names only by UniProt
#      accession are resolved to gene symbols by the README's four rules; unresolved,
#      5,870 edges could never match a gene.
#   3. ONLY THE VECTOR SERIES IS FITTED. MYAZ, FaMY and BoMY are unidentified and
#      excluded repo-wide, and both contrasts sit inside the vector series.
#
# THE COHORT'S FIRST DESeq2 WALD FIT. Scripts 50-53 and 57 deliberately fitted none
# (57's header). TF activity needs a gene-level statistic, so this fit exists for
# the two contrasts above and nothing else. Adjusted P values are IHW, as in
# script 03 and FIGURE_RULES.
#
# THIS IS A SELECTION EXPERIMENT (script 52's header). PGC-1a kills MYAZ cells, so
# the arm is a survivor population. Every sentence takes the form "tumours that
# survived PGC-1a expression show X", never "PGC-1a regulates X".
#
# -----------------------------------------------------------------------------
# PART 5's GENES, FIXED HERE BEFORE ANY OF THEIR NUMBERS IS READ
# -----------------------------------------------------------------------------
# The author (2026-11-07): the mRNA needs a hypothesis -- genes linked to PGC-1a and
# to migration or invasion, and independent of PGC-1a-driven mitochondrial
# biogenesis. The literature (Consensus and PubMed, 2026-11-07) gives one such
# route. PGC-1a raises ID proteins, which sequester the bHLH E-protein TCF4 (E2-2)
# and so disable TCF4 and its cooperation with TWIST1:
#   - Luo et al. 2016, Nature 537:422 (doi:10.1038/nature19347): in melanoma PGC-1a
#     binds the ID2 promoter; ID2 inactivates TCF4; TCF4's targets, integrins among
#     them, fall. ID2 and TCF4 leave metabolism unchanged, so the axis is separable
#     from PGC-1a's bioenergetic role. Its Methods list a 23-gene "PGC1a-regulated
#     metastasis/invasion signature" -- the genes_of_interest list of the external
#     script, which therefore already followed this paper.
#   - Oh et al. 2021, Cancers 13:1772 (doi:10.3390/cancers13081772) and Fang et al.
#     2022, Int J Mol Sci 23:8247 (doi:10.3390/ijms23158247): in lung cancer PGC-1a
#     raises ID1, which keeps TCF4 and TWIST1 off the CDH1 and CDH2 promoters: CDH1
#     up; CDH2, VIM and PTHLH down. Pgc1a loss promotes bone metastasis.
#   - Nam et al. 2020, Matrix Biol: PGC-1a lowers SNAIL1 and SNAIL2 in renal
#     carcinoma, through miR-29a and collagen-DDR1 signalling.
#
#  module                           genes                                  predicted, Bcl-xL + PGC-1a
#  A  ID proteins (PGC-1a targets)   Id1 Id2 Id3                            up
#  B  TCF4-TWIST1 programme          Twist1 Cdh2 Vim Pthlh                  down
#                                    Cdh1                                   up
#                                    Tcf4                                   none (ID acts on the protein)
#  C  Luo 2016 signature, 22 of 23   Itga3 Itga4 Itga10 Itgb3 Itgb5 Cav1    down, and one set test
#     (IL8 has no mouse orthologue;  Cav2 Actn2 Lama4 Col4a1 Inhba Tgfbi
#     shown as the "Adhesion,        Tgfbr3 Tgfbr2 Smad3 Smad7 Il11 Lef1
#     TGF-beta and Wnt signature")   Tcf7l2 Dkk3 Ppp3ca Sfrp1
#
# MODULE C IS A SIGNATURE AND IS JUDGED ONLY AS A SET (author, 2026-11-08): its
# prediction is for the set, so a single gene moving the predicted way carries no
# weight. Its genes get no arrow and no bold on ED 5d; the set's arrow follows its
# name, and its test of record is the fgsea set test below. Modules A, B, D and E
# are judged gene by gene.
#  D  the text's effectors (ED 5e)   Mmp2 Mmp9                              down
#  E  SNAIL route                    Snai1 Snai2                            down
#
# ALREADY SEEN: Twist1, Mmp2, Mmp9 and Cdh1 were read in the scratch run of
# 2026-11-07 (note sections 4.4-4.5); none of the others was. Every gene is reported
# whatever it shows, and none is added or dropped after the read-out. The set test
# is fgsea on the Wald statistic (script 04's convention), one-sided in reading:
# the prediction is NES < 0.
#
# -----------------------------------------------------------------------------
# WHAT THIS SCRIPT WILL NOT DO
# -----------------------------------------------------------------------------
# * Reproduce the external analysis (author, 2026-11-07).
# * Read `Bcl2l1` in either arm that carries the construct (README): it is reported
#   as the construct check only.
# * Call any result "causal": the westerns (ED 5e) and the mRNA here are
#   associations in survivor tumours.
# * Adjust any gene for composition. These are bulk tumours, and a gene the host's
#   stroma expresses moves with the stromal fraction too. PART 5b reads script 52's
#   compartment markers in the same two contrasts, as context; whether a gene is
#   also to be read adjusted is the author's call (open, 2026-11-07).
#
# DRY RUN. Claude dry-runs this script with the saves sent to a scratch folder:
#   options(myc.s58_rds = "<dir>/orthotopic_tf_activity.rds", myc.s58_outdir = "<dir>")
# The author's run, with neither option set, is the record.
#
# Reads : data/orthotopic_series/salmon.merged.gene_counts.tsv
#         results/orthotopic_escape_series.rds        (script 52: the sample table)
#         data/collectri_mouse/collectri_mouse_omnipath.tsv.gz
#         data/collectri_mouse/uniprot_accessions_2026-11-07.tsv
# Writes: results/orthotopic_tf_activity.rds
#         outputs/orthotopic/58_tf_activity_panel.csv     (the 25 per contrast, drawn on ED 5d)
#         outputs/orthotopic/58_tf_activity_all.csv       (every TF scored)
#         outputs/orthotopic/58_invasion_genes.csv        (PART 5, both contrasts)
#         outputs/orthotopic/58_invasion_genes_cpm_per_animal.csv (for Prism)
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

set.seed(42)

VEC        <- c("EV", "BclxL", "Pgc1a_BclxL")
CONTRASTS  <- list(vs_bclxl = c("Pgc1a_BclxL", "BclxL"),
                   vs_ev    = c("Pgc1a_BclxL", "EV"))
CONTRAST_LABELS <- c(vs_bclxl = "Bcl-xL + PGC-1a vs Bcl-xL + EV",
                     vs_ev    = "Bcl-xL + PGC-1a vs EV + EV")
CONSTRUCTS <- c("Ppargc1a", "Bcl2l1")    # out of the gene statistics (header, 1)
MIN_COUNT  <- 10L                        # the original analysis's prefilter:
MIN_SAMPLES <- 2L                        #   >= 10 counts in >= 2 samples
MINSIZE    <- 5L                         # run_ulm(): targets present, the original's value
N_PANEL    <- 25L                        # TFs drawn per contrast: the 25 largest |score|
FDR        <- 0.05                       # BH across the TFs scored, per contrast
NET_MD5    <- "d138f1f2689d5dee6416ba302a2f847f"   # data/collectri_mouse/README.md
ACC_MD5    <- "54e5e309d44ec6ddca41ca35b95210d2"
ACC_RX     <- "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"

OUT_RDS <- getOption("myc.s58_rds", here::here("results", "orthotopic_tf_activity.rds"))
OUT_DIR <- getOption("myc.s58_outdir", here::here("outputs", "orthotopic"))
dir.create(dirname(OUT_RDS), recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

NOTES <- character(0)
note  <- function(...) NOTES <<- c(NOTES, paste0(...))

# --- the pre-specified genes (header) ------------------------------------------
GENES <- data.frame(
  module = c(rep("A", 3), rep("B", 6), rep("C", 22), rep("D", 2), rep("E", 2)),
  gene = c("Id1", "Id2", "Id3",
           "Twist1", "Cdh2", "Vim", "Pthlh", "Cdh1", "Tcf4",
           "Itga3", "Itga4", "Itga10", "Itgb3", "Itgb5", "Cav1", "Cav2", "Actn2", "Lama4",
           "Col4a1", "Inhba", "Tgfbi", "Tgfbr3", "Tgfbr2", "Smad3", "Smad7", "Il11", "Lef1",
           "Tcf7l2", "Dkk3", "Ppp3ca", "Sfrp1",
           "Mmp2", "Mmp9",
           "Snai1", "Snai2"),
  predicted = c(rep("up", 3),
                "down", "down", "down", "down", "up", "none",
                rep("down", 22),
                "down", "down",
                "down", "down"),
  stringsAsFactors = FALSE)
GENES$seen_before <- GENES$gene %in% c("Twist1", "Mmp2", "Mmp9", "Cdh1")
MODULES <- c(A = "Id proteins (PGC-1a targets)", B = "Tcf4-Twist1 programme",
             C = "Adhesion, TGF-beta and Wnt signature (Luo 2016)",
             D = "MMPs (the effectors named in the text)", E = "Snail route")
MODULE_ORDER <- c("A", "B", "D", "C", "E")   # the panel's order (author, 2026-11-07: MMPs above the signature)
MODULE_SOURCES <- c(A = "Luo 2016 (ID2, ID3); Oh 2021 (ID1)",
                    B = "Oh 2021; Fang 2022",
                    C = "Luo 2016, Methods (IL8 has no mouse orthologue)",
                    D = "ED 5e westerns; Mmp2 is a Twist1 target in CollecTRI",
                    E = "Nam 2020")
stopifnot(nrow(GENES) == 35L, !anyDuplicated(GENES$gene), sum(GENES$module == "C") == 22L)

# =============================================================================
# PART 0: THE 21 SAMPLES, FROM THE RECORD AND FROM AN INDEPENDENT PARSE
# -----------------------------------------------------------------------------
# Script 52's saved sample table is the record of which sample is in which arm.
# The column headers are parsed independently (the token rule of the external
# script, Pgc -> Pgc1a_BclxL after the 2026-09-09 identity correction) and must give
# the same 21 samples in the same arms.
# =============================================================================
message("58 PART 0: the 21 vector-series samples")

res52 <- readRDS(here::here("results", "orthotopic_escape_series.rds"))
sm <- as.data.frame(res52$sample_table)
sm$sample <- as.character(sm$sample)
sm <- sm[as.character(sm$group) %in% VEC, ]
sm$arm <- factor(as.character(sm$group), levels = VEC)
sm <- sm[order(sm$arm, sm$sample), c("sample", "arm", "series")]
rownames(sm) <- NULL
stopifnot(nrow(sm) == 21L, all(table(sm$arm) == 7L), all(as.character(sm$series) == "vector"))

raw <- data.table::fread(
  here::here("data", "orthotopic_series", "salmon.merged.gene_counts.tsv"),
  data.table = FALSE, check.names = FALSE)
stopifnot(identical(colnames(raw)[1:2], c("gene_id", "gene_name")))
CNT <- as.matrix(raw[, -(1:2)]); rownames(CNT) <- raw$gene_id
stopifnot(ncol(CNT) == 67L, !anyDuplicated(rownames(CNT)))

tok <- sub("^[^-]+-", "", colnames(CNT))
tok <- sub("BRNS.*$", "", tok)
tok <- sub("_.*$", "", tok)
tok <- sub("([A-Za-z]+).*", "\\1", tok)
parsed <- c(EV = "EV", BclxL = "BclxL", Pgc = "Pgc1a_BclxL")[tok]
names(parsed) <- colnames(CNT)
parsed <- parsed[!is.na(parsed)]
stopifnot(setequal(names(parsed), sm$sample),
          identical(unname(parsed[sm$sample]), as.character(sm$arm)))

# =============================================================================
# PART 1: DESeq2, AND THE GENE STATISTIC
# -----------------------------------------------------------------------------
# ~ arm over the 21; salmon's merged counts are non-integer, so they are rounded,
# as in script 52. The Wald statistic of the unshrunken fit is the gene statistic
# (CLAUDE.md); adjusted P values are IHW (filterFun = ihw, as script 03).
# Symbols come from org.Mm.eg.db (Ensembl -> symbol, first match), the annotation
# the network's accession names are resolved against in PART 2.
# =============================================================================
message("58 PART 1: DESeq2 on the vector series")

cd  <- data.frame(row.names = sm$sample, arm = sm$arm)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(CNT[, sm$sample]),
                                      colData = cd, design = ~ arm)
keep <- rowSums(DESeq2::counts(dds) >= MIN_COUNT) >= MIN_SAMPLES
dds  <- dds[keep, ]
dds  <- DESeq2::DESeq(dds, quiet = TRUE)
n_fit <- nrow(dds)

sym <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, keys = rownames(dds),
                             column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

de <- lapply(CONTRASTS, function(cn) {
  r <- DESeq2::results(dds, contrast = c("arm", cn), filterFun = IHW::ihw)
  d <- as.data.frame(r)
  d$ensembl <- rownames(d)
  d$symbol  <- unname(sym[d$ensembl])
  d$lo95 <- d$log2FoldChange - 1.96 * d$lfcSE
  d$hi95 <- d$log2FoldChange + 1.96 * d$lfcSE
  rownames(d) <- NULL
  d[, c("ensembl", "symbol", "baseMean", "log2FoldChange", "lfcSE", "lo95", "hi95",
        "stat", "pvalue", "padj")]
})

# one row per symbol: the gene statistic with the largest |stat|, the original
# analysis's rule for duplicated symbols; Cook's outliers (pvalue NA) are out
stat_vec <- function(d, drop = CONSTRUCTS) {
  d <- d[is.finite(d$stat) & !is.na(d$pvalue) & !is.na(d$symbol) & !(d$symbol %in% drop), ]
  d <- d[order(-abs(d$stat)), ]
  d <- d[!duplicated(d$symbol), ]
  stats::setNames(d$stat, d$symbol)
}
sv <- lapply(de, stat_vec)

# the construct check: the transgene must be the largest change in both contrasts
lookup <- function(d, g) {
  r <- d[!is.na(d$symbol) & d$symbol == g, ]
  if (!nrow(r)) return(NULL)
  r[which.max(r$baseMean), ]
}
cons <- do.call(rbind, lapply(names(de), function(k) do.call(rbind, lapply(CONSTRUCTS, function(g) {
  r <- lookup(de[[k]], g)
  data.frame(contrast = k, gene = g, log2FC = r$log2FoldChange, stat = r$stat, padj = r$padj)
}))))
stopifnot(all(cons$log2FC[cons$gene == "Ppargc1a"] > 8),
          all(!(CONSTRUCTS %in% unlist(lapply(sv, names)))))

# =============================================================================
# PART 2: THE NETWORK -- THE PINNED SNAPSHOT, PROCESSED AS get_collectri() DOES,
#         AND ITS ACCESSION NAMES RESOLVED
# -----------------------------------------------------------------------------
# decoupleR::get_collectri(organism = "mouse", split_complexes = FALSE), decoupleR
# 2.17.0, line for line: True/False to 1/0; COMPLEX sources to AP1/NFKB (the mouse
# table has none); one edge per source-target pair; mor +1 for stimulation, -1
# otherwise. Then the README's four rules for accession names. A pair that the
# resolution makes twice with OPPOSITE signs is dropped, being undecidable; a pair
# made twice with the same sign is kept once.
# =============================================================================
message("58 PART 2: the pinned CollecTRI mouse network")

net_gz  <- here::here("data", "collectri_mouse", "collectri_mouse_omnipath.tsv.gz")
net_raw <- memDecompress(readBin(net_gz, "raw", file.info(net_gz)$size), type = "gzip")
tmp <- tempfile(fileext = ".tsv"); writeBin(net_raw, tmp)
stopifnot(unname(tools::md5sum(tmp)) == NET_MD5)
ws <- utils::read.delim(tmp, check.names = FALSE, colClasses = "character", quote = "",
                        na.strings = character(0))
unlink(tmp)
stopifnot(nrow(ws) == 39849L)

tf01 <- function(v) as.integer(toupper(v) %in% c("TRUE", "1"))
ws$is_stimulation <- tf01(ws$is_stimulation)
ws$is_inhibition  <- tf01(ws$is_inhibition)
cols <- c("source_genesymbol", "target_genesymbol", "is_stimulation", "is_inhibition")
ci <- ws[!grepl("COMPLEX", ws$source), cols]
cc <- ws[grepl("COMPLEX", ws$source), cols]
cc$source_genesymbol <- dplyr::case_when(
  grepl("JUN", cc$source_genesymbol) | grepl("FOS", cc$source_genesymbol) ~ "AP1",
  grepl("REL", cc$source_genesymbol) | grepl("NFKB", cc$source_genesymbol) ~ "NFKB")
net0 <- rbind(ci, cc)
net0 <- net0[!duplicated(net0[, c("source_genesymbol", "target_genesymbol")]), ]
net0 <- data.frame(source = net0$source_genesymbol, target = net0$target_genesymbol,
                   mor = ifelse(net0$is_stimulation == 1, 1, -1), stringsAsFactors = FALSE)
net0 <- net0[!is.na(net0$source), ]
stopifnot(nrow(cc) == 0L, nrow(net0) == 39521L)

acc_path <- here::here("data", "collectri_mouse", "uniprot_accessions_2026-11-07.tsv")
stopifnot(unname(tools::md5sum(acc_path)) == ACC_MD5)
acc_tab <- utils::read.delim(acc_path, colClasses = "character", na.strings = character(0))
accs <- sort(unique(c(net0$source, net0$target)))
accs <- accs[grepl(ACC_RX, accs)]
stopifnot(length(accs) == 1077L, setequal(accs, acc_tab$accession))

om <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = accs,
                                             keytype = "UNIPROT", columns = "SYMBOL"))
om <- om[!is.na(om$SYMBOL), ]
org_syms <- split(om$SYMBOL, om$UNIPROT)
valid    <- AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db, keytype = "SYMBOL")
resolve_one <- function(a) {
  o <- unique(org_syms[[a]])
  u <- acc_tab$uniprot_gene_primary[acc_tab$accession == a]
  if (length(o) == 1L) return(c(o, "1 org.Mm.eg.db, one symbol"))
  if (length(o) > 1L && nzchar(u) && u %in% o) return(c(u, "2 UniProt, among org.Mm.eg.db's"))
  if (length(o) == 0L && nzchar(u) && u %in% valid) return(c(u, "3 UniProt, a current symbol"))
  c(a, "4 unresolved")
}
rr <- t(vapply(accs, resolve_one, character(2)))
acc_map <- data.frame(accession = accs, symbol = rr[, 1], rule = rr[, 2],
                      stringsAsFactors = FALSE, row.names = NULL)

relabel <- function(x) ifelse(x %in% accs, acc_map$symbol[match(x, accs)], x)
net1 <- data.frame(source = relabel(net0$source), target = relabel(net0$target), mor = net0$mor,
                   stringsAsFactors = FALSE)
key <- paste(net1$source, net1$target, sep = "\r")
sign_n <- tapply(net1$mor, key, function(m) length(unique(m)))
conflict_keys <- names(sign_n)[sign_n > 1]
net <- net1[!(key %in% conflict_keys), ]
net <- net[!duplicated(paste(net$source, net$target, sep = "\r")), ]
rownames(net) <- NULL

old_sources <- unique(net0$source[!grepl(ACC_RX, net0$source)])
renamed_src <- unique(acc_map$symbol[acc_map$accession %in% net0$source & acc_map$rule != "4 unresolved"])
net_summary <- list(
  edges_snapshot = nrow(net0), sources_snapshot = length(unique(net0$source)),
  accession_nodes = length(accs),
  rules = as.data.frame(table(rule = acc_map$rule), stringsAsFactors = FALSE),
  edges_accession_target_before = sum(grepl(ACC_RX, net0$target)),
  edges_accession_target_after  = sum(grepl(ACC_RX, net$target)),
  edges_accession_source_before = sum(grepl(ACC_RX, net0$source)),
  edges_accession_source_after  = sum(grepl(ACC_RX, net$source)),
  pairs_merged_same_sign = nrow(net1) - length(unique(key)),
  pairs_dropped_conflicting_sign = length(conflict_keys),
  resolved_sources_colliding = sum(renamed_src %in% old_sources),
  edges_final = nrow(net), sources_final = length(unique(net$source)))
if (net_summary$resolved_sources_colliding > 0)
  note(net_summary$resolved_sources_colliding, " resolved TF name(s) already existed as a ",
       "symbol-named TF; their regulons were merged (same-sign duplicates kept once).")

# =============================================================================
# PART 3: TF ACTIVITY (decoupleR ULM) AND THE PANEL TABLE
# -----------------------------------------------------------------------------
# Per TF, run_ulm() regresses the gene statistic on the TF's signed targets (mor,
# 0 for every other gene) and reports the t-value of the slope as the score, with
# its two-sided P. The P treats genes as independent, so it is anti-conservative;
# BH across the TFs scored is the panel's test of record, and the panel's 25 are
# chosen by |score| alone, one rule a legend can state.
# =============================================================================
message("58 PART 3: TF activity")

ulm_one <- function(v) {
  m <- matrix(v, ncol = 1, dimnames = list(names(v), "stat"))
  a <- as.data.frame(decoupleR::run_ulm(mat = m, network = net, .source = "source",
                                        .target = "target", .mor = "mor", minsize = MINSIZE))
  a$padj_bh <- stats::p.adjust(a$p_value, method = "BH")
  present <- net[net$target %in% names(v), ]
  nt <- table(present$source)
  a$n_targets <- as.integer(nt[a$source])
  a$rank_abs  <- rank(-abs(a$score), ties.method = "first")
  a <- a[order(-a$score), c("source", "score", "p_value", "padj_bh", "n_targets", "rank_abs")]
  rownames(a) <- NULL
  a
}
ulm <- lapply(sv, ulm_one)
stopifnot(all(vapply(ulm, nrow, integer(1)) > 500L))

panel <- do.call(rbind, lapply(names(ulm), function(k) {
  a <- ulm[[k]]
  p <- a[a$rank_abs <= N_PANEL, ]
  p$contrast <- k
  p$label    <- CONTRAST_LABELS[[k]]
  p$bh_sig   <- p$padj_bh < FDR
  p[order(-p$score), c("contrast", "label", "source", "score", "p_value", "padj_bh",
                       "bh_sig", "n_targets", "rank_abs")]
}))
rownames(panel) <- NULL
stopifnot(all(table(panel$contrast) == N_PANEL))
left_acc <- unique(panel$source[grepl(ACC_RX, panel$source)])
if (length(left_acc))
  note("A TF still named by accession is among the panel's 25: ", paste(left_acc, collapse = ", "))

# =============================================================================
# PART 4: TWIST1 -- ITS REGULON, ITS DRIVERS, AND WHAT THE CONSTRUCTS DID
# -----------------------------------------------------------------------------
# The text names Twist1 and its targets Mmp2 and Mmp9. Reported: whether each is a
# Twist1 target here; the targets that carry Twist1's score (stat x mor); and the
# score again with each construct transcript put back, the sensitivity behind
# correction 1. The record is the constructs-out score.
# =============================================================================
message("58 PART 4: Twist1")

tw_reg <- net[net$source == "Twist1", ]
mmp_membership <- data.frame(target = c("Mmp2", "Mmp9"),
                             twist1_target = c("Mmp2", "Mmp9") %in% tw_reg$target,
                             mor = tw_reg$mor[match(c("Mmp2", "Mmp9"), tw_reg$target)])
tw_drivers <- lapply(names(sv), function(k) {
  v <- sv[[k]]
  d <- tw_reg[tw_reg$target %in% names(v), ]
  d$stat <- unname(v[d$target])
  d$contribution <- d$stat * d$mor
  d <- d[order(d$contribution), c("target", "mor", "stat", "contribution")]
  rownames(d) <- NULL
  d
}); names(tw_drivers) <- names(sv)

tw_score <- function(a) a[a$source == "Twist1", c("score", "p_value", "padj_bh", "rank_abs")]
tw_sens <- do.call(rbind, lapply(names(de), function(k) {
  variants <- list(`constructs out (record)` = CONSTRUCTS,
                   `Ppargc1a back in`        = "Bcl2l1",
                   `Bcl2l1 back in`          = "Ppargc1a",
                   `both back in`            = character(0))
  do.call(rbind, lapply(names(variants), function(vn) {
    a <- if (vn == "constructs out (record)") ulm[[k]] else ulm_one(stat_vec(de[[k]], drop = variants[[vn]]))
    cbind(data.frame(contrast = k, variant = vn), tw_score(a))
  }))
}))
rownames(tw_sens) <- NULL
stopifnot(nrow(tw_reg) >= MINSIZE, all(table(tw_sens$contrast) == 4L))

# =============================================================================
# PART 5: THE PRE-SPECIFIED GENES (header), BOTH CONTRASTS
# -----------------------------------------------------------------------------
# Per gene: the DESeq2 log2 fold change with its 95% Wald interval, Wald P and IHW-
# adjusted P (genome-wide), and whether the sign agrees with the prediction. A gene
# below the expression filter is reported as such. Where a symbol has several
# Ensembl rows, the one with the largest baseMean is read. Module C also gets its
# one set test: fgsea on the Wald statistic, constructs out, as script 04 runs it.
# =============================================================================
message("58 PART 5: the pre-specified invasion genes")

gene_tab <- do.call(rbind, lapply(names(de), function(k) {
  do.call(rbind, lapply(seq_len(nrow(GENES)), function(i) {
    g <- GENES$gene[i]
    r <- lookup(de[[k]], g)
    base <- data.frame(contrast = k, module = GENES$module[i], gene = g,
                       predicted = GENES$predicted[i], seen_before = GENES$seen_before[i])
    if (is.null(r)) return(cbind(base, ensembl = NA_character_, baseMean = NA_real_,
                                 log2FC = NA_real_, lo95 = NA_real_, hi95 = NA_real_,
                                 pvalue = NA_real_, padj = NA_real_,
                                 status = "below the expression filter"))
    cbind(base, ensembl = r$ensembl, baseMean = r$baseMean, log2FC = r$log2FoldChange,
          lo95 = r$lo95, hi95 = r$hi95, pvalue = r$pvalue, padj = r$padj, status = "fitted")
  }))
}))
rownames(gene_tab) <- NULL
gene_tab$agrees <- ifelse(gene_tab$predicted == "none" | is.na(gene_tab$log2FC), NA,
                          (gene_tab$predicted == "up") == (gene_tab$log2FC > 0))
stopifnot(nrow(gene_tab) == 2L * nrow(GENES))

module_summary <- do.call(rbind, lapply(names(de), function(k) {
  do.call(rbind, lapply(names(MODULES), function(m) {
    t <- gene_tab[gene_tab$contrast == k & gene_tab$module == m, ]
    tp <- t[t$predicted != "none" & t$status == "fitted", ]
    data.frame(contrast = k, module = m, name = MODULES[[m]], genes = nrow(t),
               fitted = sum(t$status == "fitted"), with_prediction = nrow(tp),
               agree_in_sign = sum(tp$agrees, na.rm = TRUE),
               agree_and_padj_05 = sum(tp$agrees & tp$padj < 0.05, na.rm = TRUE),
               against_and_padj_05 = sum(!tp$agrees & tp$padj < 0.05, na.rm = TRUE))
  }))
}))
rownames(module_summary) <- NULL

set_test <- do.call(rbind, lapply(names(sv), function(k) {
  luo <- intersect(GENES$gene[GENES$module == "C"], names(sv[[k]]))
  set.seed(42)
  f <- fgsea::fgseaMultilevel(pathways = list(luo2016_signature = luo), stats = sv[[k]],
                              minSize = 10, maxSize = 500, eps = 0)
  data.frame(contrast = k, set = "Luo 2016 signature (22 mouse genes)", genes_present = length(luo),
             size = f$size, NES = f$NES, pval = f$pval, padj = f$padj,
             leading_edge = paste(unlist(f$leadingEdge), collapse = " "))
}))
rownames(set_test) <- NULL
stopifnot(nrow(set_test) == 2L)

# --- PART 5b: composition, read as context (nothing is adjusted) ---------------
# Script 52's compartment markers, verbatim (52, lines 285-290), in the same two
# contrasts. Fibroblasts express Mmp2 and Twist1, so a lower stromal fraction alone
# would lower both; the read-out says how large that difference is.
COMPARTMENTS <- list(
  epithelial    = c("Krt8", "Krt18", "Epcam", "Cdh1", "Krt5", "Krt14"),
  stromal       = c("Col1a1", "Pdgfrb", "Acta2", "Thy1"),
  adipose       = c("Adipoq", "Plin1", "Fabp4", "Cidec"),
  immune        = c("Ptprc", "Cd52", "Lyz2", "Cd74"),
  endothelial   = c("Pecam1", "Cdh5", "Cldn5", "Emcn", "Tek"),
  proliferation = c("Mki67", "Top2a", "Ccnb1"))
composition <- do.call(rbind, lapply(names(de), function(k) {
  do.call(rbind, lapply(names(COMPARTMENTS), function(cmp) {
    rows <- lapply(COMPARTMENTS[[cmp]], function(g) lookup(de[[k]], g))
    names(rows) <- COMPARTMENTS[[cmp]]
    rows <- rows[!vapply(rows, is.null, logical(1))]
    fc <- vapply(rows, function(r) r$log2FoldChange, numeric(1))
    data.frame(contrast = k, compartment = cmp, markers_found = length(fc),
               median_log2FC = stats::median(fc),
               markers = paste(sprintf("%s %+.2f", names(fc), fc), collapse = ", "))
  }))
}))
rownames(composition) <- NULL
stopifnot(all(composition$markers_found >= 3L))

# =============================================================================
# PART 6: PER-ANIMAL VALUES FOR PRISM (CPM, the orthotopic reporting unit, script 57)
# =============================================================================
message("58 PART 6: per-animal CPM for Prism")

CPM <- t(t(CNT[, sm$sample]) / colSums(CNT[, sm$sample])) * 1e6
cpm_rows <- unique(gene_tab[!is.na(gene_tab$ensembl), c("gene", "ensembl")])
cons_rows <- do.call(rbind, lapply(CONSTRUCTS, function(g) {
  r <- lookup(de[[1]], g); data.frame(gene = g, ensembl = r$ensembl) }))
cpm_rows <- rbind(cons_rows, cpm_rows[!(cpm_rows$gene %in% CONSTRUCTS), ])
cpm_long <- do.call(rbind, lapply(seq_len(nrow(cpm_rows)), function(i) {
  data.frame(gene = cpm_rows$gene[i], ensembl = cpm_rows$ensembl[i],
             arm = as.character(sm$arm), sample = sm$sample,
             cpm = as.numeric(CPM[cpm_rows$ensembl[i], sm$sample]))
}))
stopifnot(!anyNA(cpm_long$cpm))

# =============================================================================
# PART 6b: THE HEATMAP'S VALUES, IN FIG. 1g's FORM
# -----------------------------------------------------------------------------
# VST of this fit (blind = FALSE), z-scored per gene over the 21 tumours, the mean
# per arm: for the 35 genes and for script 52's compartment markers. Descriptive
# values for a panel; the tests are PART 5's Wald/IHW P values, carried alongside
# as one flag per gene and contrast. Set rows (the Luo signature, each compartment)
# are the mean of their genes' arm means.
# =============================================================================
message("58 PART 6b: heatmap values")

vsd <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))
hm_rows <- rbind(
  data.frame(block = "genes", set = GENES$module, gene = GENES$gene),
  do.call(rbind, lapply(names(COMPARTMENTS), function(cmp)
    data.frame(block = "composition", set = cmp, gene = COMPARTMENTS[[cmp]]))))
hm <- do.call(rbind, lapply(seq_len(nrow(hm_rows)), function(i) {
  r <- lookup(de[[1]], hm_rows$gene[i])
  if (is.null(r)) return(NULL)
  z <- as.numeric(scale(vsd[r$ensembl, sm$sample]))
  m <- tapply(z, sm$arm, mean)
  data.frame(hm_rows[i, ], ensembl = r$ensembl, arm = names(m), mean_z = as.numeric(m),
             row.names = NULL)
}))
hm_sets <- stats::aggregate(mean_z ~ block + set + arm, data = hm, FUN = mean)
hm_sig <- data.frame(
  gene = GENES$gene,
  padj_vs_bclxl = gene_tab$padj[gene_tab$contrast == "vs_bclxl"][match(GENES$gene, gene_tab$gene[gene_tab$contrast == "vs_bclxl"])],
  padj_vs_ev    = gene_tab$padj[gene_tab$contrast == "vs_ev"][match(GENES$gene, gene_tab$gene[gene_tab$contrast == "vs_ev"])])
hm_sig$sig_vs_bclxl <- !is.na(hm_sig$padj_vs_bclxl) & hm_sig$padj_vs_bclxl < 0.05
hm_sig$sig_vs_ev    <- !is.na(hm_sig$padj_vs_ev) & hm_sig$padj_vs_ev < 0.05
# THE PANEL'S TWO MARKS (author, 2026-11-07), decided here so the panel only draws:
#   star: Bcl-xL + PGC-1a differs from AT LEAST ONE control (IHW-adjusted P < 0.05);
#   bold: the change is in the PREDICTED direction -- significant (IHW < 0.05) against
#         at least one control in that direction, and not significant against the
#         prediction in the other. Significance is required so that a change of
#         about zero (Tgfbr2, -0.00) cannot read as agreeing. Tcf4 has no prediction,
#         and the signature's genes are not judged one by one (module C, header).
lfc_of <- function(k) gene_tab$log2FC[gene_tab$contrast == k][match(GENES$gene, gene_tab$gene[gene_tab$contrast == k])]
hm_sig$lfc_vs_bclxl <- lfc_of("vs_bclxl")
hm_sig$lfc_vs_ev    <- lfc_of("vs_ev")
hm_sig$predicted    <- GENES$predicted
hm_sig$star <- hm_sig$sig_vs_bclxl | hm_sig$sig_vs_ev
right <- function(lfc) ifelse(hm_sig$predicted == "up", lfc > 0, lfc < 0)
hm_sig$bold <- hm_sig$predicted != "none" & GENES$module != "C" &   # the signature: set test only
  ((hm_sig$sig_vs_bclxl & right(hm_sig$lfc_vs_bclxl)) | (hm_sig$sig_vs_ev & right(hm_sig$lfc_vs_ev))) &
  !(hm_sig$sig_vs_bclxl & !right(hm_sig$lfc_vs_bclxl)) & !(hm_sig$sig_vs_ev & !right(hm_sig$lfc_vs_ev))
stopifnot(all(GENES$gene %in% hm$gene), all(abs(hm$mean_z) < 3), !anyNA(hm_sig$bold),
          all(hm_sig$star[hm_sig$bold]))

# =============================================================================
# PART 7: SAVE, AND THE NUMBERS FOR THE LEGEND
# =============================================================================
message("58 PART 7: save")

out <- list(
  params = list(samples = 21L, arms = VEC, contrasts = CONTRASTS, contrast_labels = CONTRAST_LABELS,
                constructs_left_out = CONSTRUCTS, prefilter = c(min_count = MIN_COUNT,
                min_samples = MIN_SAMPLES), ulm_minsize = MINSIZE, n_panel = N_PANEL, fdr = FDR,
                network = "data/collectri_mouse/collectri_mouse_omnipath.tsv.gz",
                network_md5 = NET_MD5, accession_table_md5 = ACC_MD5),
  sample_table = sm,
  genes_fitted = n_fit,
  genes_into_ulm = vapply(sv, length, integer(1)),
  de = de,
  constructs = cons,
  network = list(summary = net_summary, accession_map = acc_map, edges = net),
  ulm = ulm,
  panel = panel,
  twist1 = list(regulon = tw_reg, mmp_membership = mmp_membership, drivers = tw_drivers,
                sensitivity = tw_sens),
  genes = list(spec = GENES, modules = MODULES, module_order = MODULE_ORDER, module_sources = MODULE_SOURCES,
               readout = gene_tab, module_summary = module_summary, set_test = set_test,
               cpm_per_animal = cpm_long),
  composition = composition,
  heatmap = list(genes = hm, sets = hm_sets, sig = hm_sig,
                 scale = "VST (blind = FALSE), z per gene over the 21 tumours, mean per arm"),
  versions = c(R = R.version.string,
               DESeq2 = as.character(utils::packageVersion("DESeq2")),
               IHW = as.character(utils::packageVersion("IHW")),
               decoupleR = as.character(utils::packageVersion("decoupleR")),
               fgsea = as.character(utils::packageVersion("fgsea")),
               org.Mm.eg.db = as.character(utils::packageVersion("org.Mm.eg.db"))),
  notes = NOTES,
  analysis_date = Sys.Date())
saveRDS(out, OUT_RDS)

utils::write.csv(panel, file.path(OUT_DIR, "58_tf_activity_panel.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, lapply(names(ulm), function(k) cbind(contrast = k, ulm[[k]]))),
                 file.path(OUT_DIR, "58_tf_activity_all.csv"), row.names = FALSE)
utils::write.csv(gene_tab, file.path(OUT_DIR, "58_invasion_genes.csv"), row.names = FALSE)
utils::write.csv(cpm_long, file.path(OUT_DIR, "58_invasion_genes_cpm_per_animal.csv"), row.names = FALSE)
message("wrote ", OUT_RDS, " and four CSVs in ", OUT_DIR)

# --- the read-out ----------------------------------------------------------------
cat("\n58 -- read-out (", format(Sys.Date()), ")\n", sep = "")
cat(sprintf("  samples 21 (7/7/7); genes fitted %d; into ULM: %s\n", n_fit,
            paste(sprintf("%s %d", names(sv), vapply(sv, length, integer(1))), collapse = ", ")))
cat("  construct check (log2FC):\n")
for (i in seq_len(nrow(cons)))
  cat(sprintf("    %-9s %-8s %+.2f (padj %.2g)\n", cons$contrast[i], cons$gene[i], cons$log2FC[i], cons$padj[i]))
s <- net_summary
cat(sprintf("  network: %d edges, %d TFs in the snapshot -> %d edges, %d TFs after resolution\n",
            s$edges_snapshot, s$sources_snapshot, s$edges_final, s$sources_final))
cat(sprintf("    accession-named targets %d -> %d edges; sources %d -> %d edges; same-sign pairs merged %d; conflicting pairs dropped %d\n",
            s$edges_accession_target_before, s$edges_accession_target_after,
            s$edges_accession_source_before, s$edges_accession_source_after,
            s$pairs_merged_same_sign, s$pairs_dropped_conflicting_sign))
for (i in seq_len(nrow(s$rules))) cat(sprintf("    rule %-36s %4d\n", s$rules$rule[i], s$rules$Freq[i]))
for (k in names(ulm)) {
  a <- ulm[[k]]; p <- panel[panel$contrast == k, ]
  cat(sprintf("\n  %s: %d TFs scored, %d at BH < %.2f (%s); the 25th |score| %.2f\n",
              CONTRAST_LABELS[[k]], nrow(a), sum(a$padj_bh < FDR), FDR,
              paste(a$source[a$padj_bh < FDR], collapse = ", "),
              min(abs(p$score))))
  cat("    up:  ", paste(sprintf("%s %+.2f", p$source[p$score > 0], p$score[p$score > 0]), collapse = ", "), "\n")
  cat("    down:", paste(sprintf("%s %+.2f", p$source[p$score < 0], p$score[p$score < 0]), collapse = ", "), "\n")
}
cat(sprintf("\n  Twist1 (score, P, BH across the TFs scored, rank by |score| of %s):\n",
            paste(vapply(ulm, nrow, integer(1)), collapse = " / ")))
for (i in seq_len(nrow(tw_sens)))
  cat(sprintf("    %-9s %-24s %+.2f  P %.3g  BH %.2g  rank %d\n", tw_sens$contrast[i], tw_sens$variant[i],
              tw_sens$score[i], tw_sens$p_value[i], tw_sens$padj_bh[i], tw_sens$rank_abs[i]))
cat(sprintf("    Twist1 targets: %d; Mmp2 %s; Mmp9 %s\n", nrow(tw_reg),
            if (mmp_membership$twist1_target[1]) "is one" else "is not",
            if (mmp_membership$twist1_target[2]) "is one" else "is not"))
cat("\n  the pre-specified genes (log2FC [95% CI], IHW padj; * = sign as predicted):\n")
for (k in names(de)) {
  cat("   ", CONTRAST_LABELS[[k]], "\n")
  t <- gene_tab[gene_tab$contrast == k, ]
  for (i in seq_len(nrow(t))) {
    if (t$status[i] != "fitted") { cat(sprintf("      %s %-7s below the expression filter\n", t$module[i], t$gene[i])); next }
    cat(sprintf("      %s %-7s %+.2f [%+.2f, %+.2f]  padj %-8.2g pred %-4s %s\n", t$module[i], t$gene[i],
                t$log2FC[i], t$lo95[i], t$hi95[i], t$padj[i], t$predicted[i],
                ifelse(is.na(t$agrees[i]), "", ifelse(t$agrees[i], "*", ""))))
  }
}
cat("\n  module summary:\n"); print(module_summary)
cat("\n  set test (fgsea, Wald statistic):\n"); print(set_test[, c("contrast", "genes_present", "NES", "pval")])
cat("\n  composition (script 52's markers; median log2FC, context only):\n")
for (i in seq_len(nrow(composition)))
  cat(sprintf("    %-9s %-13s %+.2f  (%s)\n", composition$contrast[i], composition$compartment[i],
              composition$median_log2FC[i], composition$markers[i]))
cat(sprintf("\n  ED 5d marks: star (differs from >= 1 control) %d of %d genes; bold (predicted direction) %d of %d judged gene by gene\n",
            sum(hm_sig$star), nrow(hm_sig), sum(hm_sig$bold), sum(hm_sig$predicted != "none" & GENES$module != "C")))
cat("    bold:", paste(hm_sig$gene[hm_sig$bold], collapse = ", "), "\n")
cat("    star, not bold:", paste(hm_sig$gene[hm_sig$star & !hm_sig$bold], collapse = ", "), "\n")
if (length(NOTES)) { cat("\n  notes:\n"); cat(paste0("   - ", NOTES), sep = "\n") }

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  o <- readRDS(here::here("results", "orthotopic_tf_activity.rds"))
  o$network$summary
  head(o$ulm$vs_bclxl, 30)
  o$panel[o$panel$contrast == "vs_ev", ]
  o$twist1$sensitivity
  o$twist1$drivers$vs_bclxl
  o$genes$module_summary
  o$genes$set_test
  ## one gene per animal, as Prism would draw it
  subset(o$genes$cpm_per_animal, gene == "Twist1")
  ## any TF still named by accession on the panel
  o$panel$source[grepl("^[A-Z][0-9]", o$panel$source)]
}
