# =============================================================================
# 57 -- SECTION 3 PANEL READ-OUT: PER-ANIMAL EXPORTS, THE ROUTE COMPARISON,
#       AND THE NUMBERS THE LEGENDS AND METHODS QUOTE
# -----------------------------------------------------------------------------
# THIS IS A READ-OUT JOB, NOT AN ANALYSIS. Every contrast is re-read from
# `results/orthotopic_escape_series.rds` (script 52, committed at faf1e8a).
# Nothing script 52 fitted is re-fitted here. What this script DOES compute is
# the per-sample layer script 52 never saved -- it saved arm medians -- plus the
# three read-outs the panel spec asks for that no object holds (the within-arm
# Bbc3/Ppargc1a correlation, the Myc median contrast, the construct check).
# Where a per-sample score is rebuilt, the rebuild is ASSERTED to reproduce
# script 52's saved arm medians before it is used. If an assertion fails, the
# rebuild has drifted and nothing below it may be read.
#
# Spec: docs/2026-10-27_section3_final_narrative_and_panel_spec.md (v2), with the
# author's build list of 2026-10-28, which supersedes it where they differ.
# Limitations: docs/2026-09-23_section3_narrative_and_figure3.md section 2.5.
#
# -----------------------------------------------------------------------------
# THE UNIT, FIXED HERE AND USED THROUGHOUT (author, 2026-10-28)
# -----------------------------------------------------------------------------
# EXPRESSION LEVELS ARE CPM. Script 52, its saved object, the escape-series note
# and the identity-correction note all quote CPM -- 301.6 / 128.6 / 0.17 are CPM.
# DESeq2-normalised counts differ from CPM here by a near-constant 26.4-28.3x
# (size factors run 0.958-1.104), so the ordering and the story are identical and
# the only thing a switch would buy is renumbering four documents. The existing
# Prism axis saying "DESeq2 normalised" is relabelled to CPM instead.
#   ONE EXCEPTION, NAMED RATHER THAN HIDDEN. Script 52's gene CONTRASTS are
#   differences of medians of log2(DESeq2-normalised + 1). A fold change read
#   from this script's CPM columns will therefore not equal the saved contrast
#   exactly. Both are exported: `*_cpm` is the reporting unit, `*_log2norm` is
#   the internal scale the saved contrasts live on, so the two are reconcilable
#   and nobody has to guess which produced a number.
#
# -----------------------------------------------------------------------------
# THE P VALUE FOR FIG 3g (author, 2026-10-28)
# -----------------------------------------------------------------------------
# NO DESeq2 WALD TEST IS FITTED, HERE OR ANYWHERE. Script 52 never ran `DESeq()`
# -- it called `estimateSizeFactors` and took normalised counts -- so there is no
# Wald P to read out, and fitting one now would add a model the record does not
# have. The panel and its legend both report MANN-WHITNEY ON THE EXPORTED
# PER-ANIMAL VALUES, which is what Prism computes from the pasted replicates.
# This script prints that P so Prism can be checked against it.
#   At 7 vs 7 the two-sided floor is 2 / choose(14,7) = 5.83e-04, so "P < 0.0001"
#   is not reachable at this n and the placeholder must go. It is not doing any
#   work: the panel shows a roughly 1,100-fold transgene difference.
#
# -----------------------------------------------------------------------------
# WHAT THIS SCRIPT WILL NOT DO
# -----------------------------------------------------------------------------
# * NO composition-adjusted re-fits. PART 7 REPORTS the screen as script 52 ran
#   it: no axis reached |rho| >= 0.5, so no adjustment was applied to anything.
#   "Adjusted vs unadjusted" has no adjusted side in this cohort, and inventing
#   one now would be a new analysis wearing a read-out's clothes.
# * NO growth panel and NO blot panel (ED b, ED c). Tumour volumes and blot
#   quantification live outside the repo; both are built in Prism.
# * NO `Cdkn2a` in the KO arm. The locus is CRISPR-targeted there, so its +0.587
#   is not a read-out of anything (author, 2026-10-28). Vector and NT only.
# * NO cross-run quotes. Every number here is the 49-sample run. The spec's
#   adipose loadings (+0.621 pooled, +0.036 / -0.393 / +0.714) are the 21-SAMPLE
#   vector run from `docs/2026-09-08_orthotopic_vector_series.md:94`, where the
#   three columns are EV / BclxL / Pgc1a_BclxL -- not the three PGC-1a arms. The
#   author is correcting the spec; PART 7 prints the 49-sample values that
#   replace them.
# * NO causal verb. Every arm is a survivor population, so these are associations
#   across escape routes: "tumours that retained more PGC-1a show less Bbc3",
#   never "PGC-1a represses Bbc3".
# * N3. The word "primed" is not written of a transcript anywhere.
#
# -----------------------------------------------------------------------------
# THE CRISPR ARMS ARE POOLED, NOT CLONES (author, 2026-10-28)
# -----------------------------------------------------------------------------
# The p19ARF knockout was POOLED CRISPR with no single-cell selection, so the arm
# is p19ARF-DISRUPTED, not null. Editing was validated independently: the edited
# fraction was checked and the p19ARF western QC passed. This matters because the
# KO negative carries the closing sentence of the section.
#   Three places on disk still say "clone" and are wrong. `data/orthotopic_series/
#   README.md` is corrected alongside this script. Script 52's header comment and
#   script 53's figure caption are NOT edited: they are the committed analysis
#   behind a saved object and a published PDF, and changing their text would make
#   the files differ from what produced those outputs. The correction is recorded
#   here and in the note this script writes.
#
# Reads : results/orthotopic_escape_series.rds        (script 52 -- the record)
#         data/orthotopic_series/salmon.merged.gene_counts.tsv
#         data/Mouse.MitoCarta3.0.xls  Sheet 4        (script 08's splits)
#         functions/reconcile_gene_symbols.R          (MANDATORY)
#         docs/<facility manifest>                    (OPTIONAL -- see MANIFEST)
# Writes: results/orthotopic_section3_panels.rds
#         outputs/orthotopic/57_per_animal_long.csv
#         outputs/orthotopic/57_fig3g_ppargc1a_prism.csv
#         outputs/orthotopic/57_fig3i_route_comparison.csv
#         outputs/orthotopic/57_fig3i_route_comparison.pdf
#         outputs/orthotopic/57_eda_bbc3_dose.csv
#         outputs/orthotopic/57_eda_bbc3_dose.pdf
#         docs/2026-10-28_section3_panel_readout.md
#
# RUNTIME: under a minute. The only resampling is 5,000-draw bootstraps.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

set.seed(1)
NBOOT  <- 5000L    # script 52's resample count, so intervals are built alike
MINBM  <- 10       # script 52's expression floor -- must match or the ruler moves
CPM_ON <- 10       # a sample counts as transgene-positive above this CPM (PART 3)

VEC    <- c("EV", "BclxL", "Pgc1a_BclxL")
CRISPR <- c("NTEV", "NTPgc1a", "KOEV", "KOPgc1a")
SET    <- c(VEC, CRISPR)

# The facility manifest had not landed when this script was written. If it is
# placed at one of these paths the per-sample verification in PART 3 runs; if it
# is not, PART 3 says so rather than silently skipping. Point MANIFEST at the
# file if it lands under another name.
MANIFEST <- c(here::here("docs", "orthotopic_facility_manifest.csv"),
              here::here("docs", "orthotopic_sample_manifest.csv"),
              here::here("data", "orthotopic_series", "facility_manifest.csv"))

NOTES <- character(0)
note  <- function(...) NOTES <<- c(NOTES, paste0(...))

res_path <- here::here("results", "orthotopic_escape_series.rds")
if (!file.exists(res_path)) stop("run scripts/52_orthotopic_escape_series.R first")
res <- readRDS(res_path)
stopifnot(identical(as.character(res$scoring_set$groups), SET), res$scoring_set$n == 49L)

# =============================================================================
# PART 0: REBUILD SCRIPT 52's SCORING SET, AND PROVE IT IS THE SAME ONE
# -----------------------------------------------------------------------------
# `ox_rel` and the MitoCarta share are cohort-relative composites that script 52
# saved only as arm medians, so the per-animal layer has to be rebuilt. It is
# rebuilt with script 52's recipe line for line -- same floor, same filter, same
# z-construction, same gene sets -- and then asserted against the saved medians.
# The assertion is the point: a composite rebuilt with a different ruler would
# look perfectly reasonable and be silently incomparable with the record.
# =============================================================================
message("57 PART 0: rebuild and verify the 49-sample scoring set")

raw <- data.table::fread(
  here::here("data", "orthotopic_series", "salmon.merged.gene_counts.tsv"),
  data.table = FALSE, check.names = FALSE)
stopifnot(identical(colnames(raw)[1:2], c("gene_id", "gene_name")))
CNT <- as.matrix(raw[, -(1:2)]); rownames(CNT) <- raw$gene_id
stopifnot(ncol(CNT) == 67L)

# the 49 come from the SAVED sample table, not from a re-parse: script 52 is the
# record of which sample is in which arm, and PART 3 checks the labels separately
sm <- as.data.frame(res$sample_table)
sm$sample <- as.character(sm$sample)
sm$group  <- factor(as.character(sm$group), levels = SET)
vk <- match(sm$sample, colnames(CNT))
stopifnot(nrow(sm) == 49L, !anyNA(vk), all(table(sm$group) == 7L))

CPM    <- t(t(CNT) / colSums(CNT)) * 1e6
cpm_of <- function(sym) {
  r <- which(raw$gene_name == sym); stopifnot(length(r) >= 1L)
  if (length(r) == 1L) CPM[r, vk] else colSums(CPM[r, vk, drop = FALSE]) }

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(CNT[, vk, drop = FALSE]), colData = sm, design = ~ group)
dds <- DESeq2::estimateSizeFactors(dds)
NCv  <- DESeq2::counts(dds, normalized = TRUE)          # LINEAR -- mitoPPS scale
keep <- rowMeans(NCv) >= MINBM & matrixStats::rowSds(NCv) > 0
L    <- log2(NCv[keep, , drop = FALSE] + 1)             # LOG -- every composite

universe_all <- rownames(L)
zrow    <- function(m) t(scale(t(m)))
comp_e  <- function(e) { e <- e[!is.na(e) & e %in% rownames(L)]
                         stopifnot(length(e) >= 3); colMeans(zrow(L[e, , drop = FALSE])) }
ens_set <- function(syms) { e <- recon_to_ensembl(syms, universe_all); e[!is.na(e)] }
one_ens <- function(sym) { e <- ens_set(sym); if (!length(e)) NA_character_ else e[1] }
lg      <- function(sym) { e <- one_ens(sym); stopifnot(!is.na(e)); as.numeric(L[e, ]) }

# --- MitoCarta, script 08's split, exact filename, mt-* stripped --------------
s4  <- readxl::read_xls(here::here("data", "Mouse.MitoCarta3.0.xls"), sheet = 4)
s4  <- na.omit(dplyr::select(s4, MitoPathway, Genes))
g2p <- suppressWarnings(splitstackshape::cSplit(s4, "Genes", ","))
g2p <- as.data.frame(t(tibble::column_to_rownames(g2p, "MitoPathway")))
g2p <- tidyr::pivot_longer(g2p, cols = colnames(g2p),
                           names_to = "Pathway", values_to = "Gene")
g2p <- dplyr::mutate(na.omit(g2p), Gene = as.character(Gene))
is_mt   <- grepl("^[Mm][Tt]-", g2p$Gene)
mt_syms <- unique(g2p$Gene[is_mt])
g2p_nuc <- g2p[!is_mt, ]

ox_sub    <- ens_set(unique(g2p_nuc$Gene[g2p_nuc$Pathway == "OXPHOS subunits"]))
mito_nuc  <- ens_set(unique(g2p_nuc$Gene))
rest_mito <- setdiff(mito_nuc, ox_sub)
mt_ens    <- ens_set(mt_syms)

# the ruler is cohort-dependent; these four numbers ARE script 52's ruler
ss <- as.data.frame(res$set_sizes)
stopifnot(length(ox_sub) == ss$n[ss$set == "OXPHOS subunits (nuclear)"],
          length(rest_mito) == ss$n[ss$set == "rest of nuclear MitoCarta"],
          length(mt_ens) == ss$n[ss$set == "mtDNA-encoded"],
          nrow(L) == ss$n[ss$set == "genes carried"])

ox_rel <- comp_e(ox_sub) - comp_e(rest_mito)

# the MitoCarta share of the library: script 52 PART F's denominator, per sample
mito_all_ens <- unique(c(ox_sub, rest_mito, mt_ens))
mito_share   <- colSums(NCv[rownames(NCv) %in% mito_all_ens, , drop = FALSE]) /
                colSums(NCv)

# mitoPPS is NOT recomputed -- script 52 saved the whole matrix, per sample
MP <- res$mitopps_matrix
stopifnot("OXPHOS subunits" %in% rownames(MP), identical(colnames(MP), sm$sample))
mitopps_ox <- as.numeric(MP["OXPHOS subunits", ])

# --- THE ASSERTIONS. Arm medians must reproduce script 52's saved table -------
ET <- as.data.frame(res$escape_table); ET$arm <- as.character(ET$arm)
med_by_arm <- function(v) as.numeric(tapply(v, sm$group, stats::median)[ET$arm])
chk <- data.frame(
  quantity = c("Ppargc1a CPM", "ox_rel", "mitoPPS OXPHOS", "MitoCarta share"),
  max_abs_diff = c(
    max(abs(med_by_arm(cpm_of("Ppargc1a")) - ET$Ppargc1a_CPM)),
    max(abs(med_by_arm(ox_rel)             - ET$ox_rel)),
    max(abs(med_by_arm(mitopps_ox)         - ET$mitopps_oxphos)),
    max(abs(med_by_arm(mito_share) -
            as.data.frame(res$c2_diagnostic)$mitocarta_share[
              match(ET$arm, as.data.frame(res$c2_diagnostic)$arm)]))),
  stringsAsFactors = FALSE)
print(chk, row.names = FALSE, digits = 3)
stopifnot(all(chk$max_abs_diff < 1e-8))
note("PART 0: the per-sample rebuild reproduces script 52's saved arm medians to ",
     sprintf("%.1e", max(chk$max_abs_diff)), " on all four quantities, so the ",
     "per-animal layer below sits on the same ruler as the record.")
message("57 PART 0: rebuild verified against the saved object")

# =============================================================================
# PART 1: THE PER-ANIMAL TABLE AND THE TWO PRISM EXPORTS  (spec Step 2, Fig 3g)
# -----------------------------------------------------------------------------
# Prism computes its own statistics from pasted replicates, so it needs VALUES,
# not summaries. Two files: a long table carrying every measure for all 49, and
# the grouped two-column layout for the one panel that gets pasted.
# =============================================================================
message("57 PART 1: per-animal export")

# animal id = the sample name with its plate prefix and arm token removed, e.g.
# `A1-EV_BRNS336-5c` -> `BRNS336-5c`. The group token is not written the same way
# twice in this cohort (see the series README), so it is stripped by longest match.
TOKENS <- c("NTEV", "KOEV", "NTPgc1a", "NTPgc", "KOPgc1a", "KOPgc", "BclxL",
            "Pgc1a", "Pgc", "EV")
core   <- sub("^[A-H][0-9]-", "", sm$sample)
tok    <- vapply(core, function(s) {
  h <- TOKENS[vapply(TOKENS, function(p) grepl(paste0("^", p), s), TRUE)]
  if (length(h)) h[1] else NA_character_ }, "")
stopifnot(!anyNA(tok))
animal <- unname(mapply(function(c_, t_) sub(paste0("^", t_, "_?"), "", c_), core, tok))
stopifnot(!anyDuplicated(animal))

per_animal <- tibble::tibble(
  animal_id          = animal,
  sample             = sm$sample,
  plate              = sub("^([A-H][0-9])-.*$", "\\1", sm$sample),
  arm                = as.character(sm$group),
  series             = as.character(sm$series),
  Ppargc1a_cpm       = as.numeric(cpm_of("Ppargc1a")),
  ox_rel             = as.numeric(ox_rel),
  mitopps_oxphos     = mitopps_ox,
  mito_share         = as.numeric(mito_share),
  Bbc3_cpm           = as.numeric(cpm_of("Bbc3")),
  Myc_cpm            = as.numeric(cpm_of("Myc")),
  Bcl2l1_cpm         = as.numeric(cpm_of("Bcl2l1")),
  # the internal scale script 52's contrasts live on -- see the header
  Ppargc1a_log2norm  = lg("Ppargc1a"),
  Bbc3_log2norm      = lg("Bbc3"),
  Myc_log2norm       = lg("Myc"))
per_animal <- per_animal[order(factor(per_animal$arm, levels = SET),
                               per_animal$animal_id), ]
stopifnot(nrow(per_animal) == 49L)

dir.create(here::here("outputs", "orthotopic"), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(per_animal,
                 here::here("outputs", "orthotopic", "57_per_animal_long.csv"),
                 row.names = FALSE)

# --- the grouped layout for Fig 3g: NTEV against NTPgc1a, Ppargc1a alone ------
g3g  <- per_animal[per_animal$arm %in% c("NTEV", "NTPgc1a"), ]
wide <- data.frame(NTEV    = g3g$Ppargc1a_cpm[g3g$arm == "NTEV"],
                   NTPgc1a = g3g$Ppargc1a_cpm[g3g$arm == "NTPgc1a"])
stopifnot(nrow(wide) == 7L)
utils::write.csv(wide,
                 here::here("outputs", "orthotopic", "57_fig3g_ppargc1a_prism.csv"),
                 row.names = FALSE)

# the test the panel and its legend both report -- what Prism should reproduce
mw <- stats::wilcox.test(wide$NTPgc1a, wide$NTEV, alternative = "two.sided", exact = TRUE)
mw_floor <- 2 / choose(14, 7)
fig3g <- list(
  test = "Mann-Whitney (Wilcoxon rank sum), two-sided, exact, on the exported CPM",
  p = unname(mw$p.value), floor_at_7v7 = mw_floor,
  at_floor = isTRUE(all.equal(unname(mw$p.value), mw_floor)),
  median_NTPgc1a = stats::median(wide$NTPgc1a), median_NTEV = stats::median(wide$NTEV),
  fold = stats::median(wide$NTPgc1a) / stats::median(wide$NTEV))
note("PART 1 (Fig 3g): Mann-Whitney on the exported per-animal CPM gives P = ",
     signif(fig3g$p, 3), ". The two-sided floor at 7 vs 7 is ", signif(mw_floor, 3),
     ", so 'P < 0.0001' is NOT reachable at this n and comes off the panel. The ",
     "panel shows a ", round(fig3g$fold), "-fold median transgene difference; the ",
     "P value is not carrying the result.")

# =============================================================================
# PART 2: THE SUMMARY NUMBERS TABLE  (spec Step 3)
# -----------------------------------------------------------------------------
# Every row is READ from script 52's object. The last column flags any value that
# differs from what the escape-series note records, so a transcription slip in
# either document shows up as a number rather than as nobody noticing.
# `mito_share` is the DENOMINATOR: it is here so the reader can see where the
# compartment is flat and where it expands.
# =============================================================================
message("57 PART 2: summary numbers")

RU <- as.data.frame(res$rulers)
P5 <- as.data.frame(res$p53_panel)
DG <- as.data.frame(res$c2_diagnostic)
ARMS <- tibble::tribble(
  ~arm,           ~control, ~route,
  "BclxL",        "EV",     "no transgene (zero arm)",
  "Pgc1a_BclxL",  "BclxL",  "effector blocked",
  "NTPgc1a",      "NTEV",   "pathway intact",
  "KOPgc1a",      "KOEV",   "checkpoint removed")

pick <- function(tab, ruler_or_gene, hi, lo, col) {
  k <- if ("ruler" %in% names(tab)) tab$ruler == ruler_or_gene else tab$gene == ruler_or_gene
  r <- tab[k & tab$hi == hi & tab$lo == lo, ]
  if (!nrow(r)) return(NA_real_)
  as.numeric(r[[col]][1]) }

summary_numbers <- dplyr::bind_rows(lapply(seq_len(nrow(ARMS)), function(i) {
  a <- ARMS$arm[i]; c_ <- ARMS$control[i]
  tibble::tibble(
    route = ARMS$route[i], arm = a, control = c_,
    Ppargc1a_cpm_arm   = ET$Ppargc1a_CPM[ET$arm == a],
    Ppargc1a_cpm_ctrl  = ET$Ppargc1a_CPM[ET$arm == c_],
    ox_rel_diff        = pick(RU, "ox_rel", a, c_, "diff_median"),
    ox_rel_lo          = pick(RU, "ox_rel", a, c_, "ci_lo"),
    ox_rel_hi          = pick(RU, "ox_rel", a, c_, "ci_hi"),
    mitopps_diff       = pick(RU, "mitopps_oxphos", a, c_, "diff_median"),
    mitopps_lo         = pick(RU, "mitopps_oxphos", a, c_, "ci_lo"),
    mitopps_hi         = pick(RU, "mitopps_oxphos", a, c_, "ci_hi"),
    mito_share_arm     = DG$mitocarta_share[DG$arm == a],
    mito_share_ctrl    = DG$mitocarta_share[DG$arm == c_],
    Bbc3_diff          = pick(P5, "Bbc3", a, c_, "diff_median"),
    Bbc3_lo            = pick(P5, "Bbc3", a, c_, "ci_lo"),
    Bbc3_hi            = pick(P5, "Bbc3", a, c_, "ci_hi")) }))

# what the note records, transcribed ONCE, here, so the comparison is mechanical.
# Source: docs/2026-09-10_orthotopic_escape_series.md sections 2 and 3.
NOTE_RECORDS <- tibble::tribble(
  ~arm,          ~quantity,  ~value,  ~lo,     ~hi,     ~source,
  "NTPgc1a",     "ox_rel",   0.176,  0.039,   0.348,   "09-10 note section 2",
  "NTPgc1a",     "mitopps",  0.053,  0.030,   0.116,   "09-10 note section 2",
  "Pgc1a_BclxL", "Bbc3",    -0.731, -1.120,  -0.450,   "09-10 note section 3",
  "NTPgc1a",     "Bbc3",    -0.345, -0.890,  -0.010,   "09-10 note section 3",
  "KOPgc1a",     "Bbc3",     0.044, -0.231,   1.180,   "09-10 note section 3")
NOTE_RECORDS$object <- mapply(function(a, q) {
  r <- summary_numbers[summary_numbers$arm == a, ]
  switch(q, ox_rel = r$ox_rel_diff, mitopps = r$mitopps_diff, Bbc3 = r$Bbc3_diff) },
  NOTE_RECORDS$arm, NOTE_RECORDS$quantity)
# the note rounds to three decimals, so anything under 5e-4 is rounding, not drift
NOTE_RECORDS$agrees <- abs(NOTE_RECORDS$object - NOTE_RECORDS$value) < 5e-4
if (all(NOTE_RECORDS$agrees)) {
  note("PART 2: all five values the escape-series note records reproduce from the ",
       "saved object to rounding. No flag.")
} else {
  note("PART 2: *** ", sum(!NOTE_RECORDS$agrees), " of ", nrow(NOTE_RECORDS),
       " recorded values DISAGREE with the object. The object is right and the ",
       "note needs a correction: ",
       paste(sprintf("%s %s (note %.3f, object %.3f)",
                     NOTE_RECORDS$arm[!NOTE_RECORDS$agrees],
                     NOTE_RECORDS$quantity[!NOTE_RECORDS$agrees],
                     NOTE_RECORDS$value[!NOTE_RECORDS$agrees],
                     NOTE_RECORDS$object[!NOTE_RECORDS$agrees]), collapse = "; "))
}

# the sentence the author flagged: a transcript-level tilt beside a flat western
nt <- summary_numbers[summary_numbers$arm == "NTPgc1a", ]
note("PART 2 (the NTPgc1a ox_rel the manuscript sentence turns on): ",
     sprintf("%+.3f [%+.3f, %+.3f]", nt$ox_rel_diff, nt$ox_rel_lo, nt$ox_rel_hi),
     ", excluding zero, at TRANSCRIPT level, in the arm whose western shows ETC ",
     "subunits unchanged. The compartment does not expand with it: MitoCarta ",
     "share ", sprintf("%.1f%% vs %.1f%%", 100 * nt$mito_share_arm,
                       100 * nt$mito_share_ctrl),
     ". That is a TILT, not mass expansion, and it is why the prose says 'OXPHOS ",
     "subunit protein' in the 3g sentence and makes 3i a dose-scaled relation.")

# =============================================================================
# PART 3: SAMPLE IDENTITY -- MANIFEST VERIFICATION AND THE CONSTRUCT CHECK
# -----------------------------------------------------------------------------
# Arm identity in this cohort has no on-disk metadata: the column headers are the
# only record, and the `Pgc1a` -> `Pgc1a_BclxL` correction came from the
# manuscript's Figure 5 legend. The author confirmed the assignment with a
# colleague and a facility manifest is being added. Two passes:
#   (1) if the manifest is present, every sample's arm is checked against it;
#   (2) regardless, a construct self-consistency check on expression.
#
# WHAT THE CONSTRUCT CHECK CAN AND CANNOT DO, stated before it runs.
#   `Ppargc1a` IS a per-sample transgene assay: the carrying arms sit three
#   orders of magnitude above the rest, so a single sample in the wrong arm is
#   visible. It CANNOT distinguish "transgene lost" from "mislabelled" for a
#   WHOLE arm -- which is exactly `KOPgc1a`, a recorded biological result, not a
#   labelling fault. The check therefore reports designed-vs-observed per arm and
#   flags INDIVIDUAL samples against their own arm, never a whole arm.
#   `Bcl2l1` is NOT a per-sample construct assay. Gene-level counts pool the
#   construct with the endogenous transcript, and BCL-XL protein is equal between
#   BclxL and Pgc1a_BclxL by western, so a level difference is not differential
#   construct expression. It is reported as a distribution, not used as a gate.
# =============================================================================
message("57 PART 3: identity and construct checks")

man_path <- MANIFEST[file.exists(MANIFEST)][1]
if (!is.na(man_path)) {
  man <- utils::read.csv(man_path, stringsAsFactors = FALSE)
  cols <- tolower(names(man))
  s_col <- names(man)[cols %in% c("sample", "sample_id", "sample_name", "library")][1]
  a_col <- names(man)[cols %in% c("arm", "group", "cohort", "construct")][1]
  if (is.na(s_col) || is.na(a_col)) {
    manifest_check <- list(status = "PRESENT BUT UNREADABLE", path = man_path,
                           detail = paste("no sample/arm columns found; saw:",
                                          paste(names(man), collapse = ", ")))
  } else {
    m <- data.frame(sample = as.character(man[[s_col]]),
                    manifest_arm = as.character(man[[a_col]]), stringsAsFactors = FALSE)
    cmp <- merge(data.frame(sample = sm$sample, parsed_arm = as.character(sm$group),
                            stringsAsFactors = FALSE), m, by = "sample", all.x = TRUE)
    cmp$agrees <- !is.na(cmp$manifest_arm) & cmp$manifest_arm == cmp$parsed_arm
    manifest_check <- list(status = if (all(cmp$agrees)) "ALL 49 AGREE" else "MISMATCH",
                           path = man_path, n_missing = sum(is.na(cmp$manifest_arm)),
                           mismatches = cmp[!cmp$agrees, ])
  }
} else {
  manifest_check <- list(
    status = "NOT PRESENT",
    detail = paste("The facility manifest had not landed when 57 was run. The arm",
                   "assignment below is script 52's, which comes from the column",
                   "headers. Place the manifest at one of:",
                   paste(basename(MANIFEST), collapse = ", "),
                   "and re-run to verify it per sample."))
}
note("PART 3 (manifest): ", manifest_check$status,
     if (!is.null(manifest_check$detail)) paste0(" -- ", manifest_check$detail) else "")

# --- designed vs observed, per arm -------------------------------------------
CARRIES_PGC   <- c("Pgc1a_BclxL", "NTPgc1a", "KOPgc1a")   # by design
CARRIES_BCLXL <- c("BclxL", "Pgc1a_BclxL")                # by design
construct_by_arm <- dplyr::bind_rows(lapply(SET, function(g) {
  k <- per_animal$arm == g
  tibble::tibble(
    arm = g, n = sum(k),
    designed_Ppargc1a = g %in% CARRIES_PGC,
    Ppargc1a_median_cpm = stats::median(per_animal$Ppargc1a_cpm[k]),
    Ppargc1a_min_cpm = min(per_animal$Ppargc1a_cpm[k]),
    Ppargc1a_max_cpm = max(per_animal$Ppargc1a_cpm[k]),
    n_positive = sum(per_animal$Ppargc1a_cpm[k] >= CPM_ON),
    designed_Bclxl = g %in% CARRIES_BCLXL,
    Bcl2l1_median_cpm = stats::median(per_animal$Bcl2l1_cpm[k]),
    Bcl2l1_min_cpm = min(per_animal$Bcl2l1_cpm[k]),
    Bcl2l1_max_cpm = max(per_animal$Bcl2l1_cpm[k])) }))
construct_by_arm$Ppargc1a_matches_design <-
  construct_by_arm$designed_Ppargc1a == (construct_by_arm$n_positive > 0)

# --- individual samples against their OWN arm --------------------------------
# a robust z on the log scale: median and MAD within the arm, so one odd sample
# does not set the standard it is judged against
rob_z <- function(v) { m <- stats::median(v); s <- stats::mad(v)
                       if (s == 0) rep(0, length(v)) else (v - m) / s }
per_animal$Ppargc1a_z_in_arm <- unsplit(
  lapply(split(log2(per_animal$Ppargc1a_cpm + 1), per_animal$arm), rob_z),
  per_animal$arm)
per_animal$Bcl2l1_z_in_arm <- unsplit(
  lapply(split(log2(per_animal$Bcl2l1_cpm + 1), per_animal$arm), rob_z),
  per_animal$arm)
construct_outliers <- per_animal[
  abs(per_animal$Ppargc1a_z_in_arm) > 5 | abs(per_animal$Bcl2l1_z_in_arm) > 5,
  c("animal_id", "sample", "arm", "Ppargc1a_cpm", "Ppargc1a_z_in_arm",
    "Bcl2l1_cpm", "Bcl2l1_z_in_arm")]

lost <- construct_by_arm[construct_by_arm$designed_Ppargc1a &
                         construct_by_arm$n_positive == 0, ]
note("PART 3 (construct check): Ppargc1a is a per-sample transgene assay and ",
     "Bcl2l1 is not -- gene-level counts pool construct with endogenous ",
     "transcript, so it is reported as a distribution, never used as a gate. ",
     if (nrow(construct_outliers) == 0)
       "No individual sample sits more than 5 robust MADs from its own arm on either gene. "
     else paste0(nrow(construct_outliers), " sample(s) sit beyond 5 robust MADs within their arm: ",
                 paste(construct_outliers$sample, collapse = ", "), ". "),
     if (nrow(lost))
       paste0("Arm(s) designed to carry Ppargc1a with no positive sample: ",
              paste(lost$arm, collapse = ", "),
              " -- this is the recorded transgene loss, not a labelling fault, and ",
              "a whole-arm loss is exactly what this check cannot tell apart from ",
              "a whole-arm mislabel.")
     else "Every arm designed to carry Ppargc1a has at least one positive sample.")

# =============================================================================
# PART 4: FIG 3i -- THE ROUTE COMPARISON  (spec Step 4b)
# -----------------------------------------------------------------------------
# Ordered by ROUTE, not by construct name and not alphabetically: zero arm,
# effector blocked, pathway intact, checkpoint removed. Each arm is labelled by
# WHAT WAS DONE TO THE DEATH PATHWAY, with the construct beneath it -- the spec
# asks for the construct as a subscript, and the CSV carries `route` and
# `construct` in separate columns so Prism and Illustrator can set a true
# subscript. The reference PDF renders it as a second line.
#
# THE DOSE AXIS SAYS "RETAINED". Retention is an outcome of selection, not a
# designed variable: PGC-1a kills MYAZ cells, so each arm is the population that
# survived expressing it.
#
# Two quantities on one figure, so the facets carry free axes -- what the panel
# is for is the two read-outs in the same route order, side by side.
# =============================================================================
message("57 PART 4: Fig 3i route comparison")

CONSTRUCT <- c(BclxL = "Bcl-xL", Pgc1a_BclxL = "PGC-1a + Bcl-xL",
               NTPgc1a = "PGC-1a, NT", KOPgc1a = "PGC-1a, p19ARF disrupted")
fig3i <- summary_numbers
fig3i$construct  <- unname(CONSTRUCT[fig3i$arm])
fig3i$label      <- paste0(fig3i$route, "\n(", fig3i$construct, ")")
fig3i$label      <- factor(fig3i$label, levels = fig3i$label)   # route order, as built
fig3i$ox_excl_0  <- fig3i$ox_rel_hi < 0 | fig3i$ox_rel_lo > 0
utils::write.csv(
  fig3i[, c("route", "arm", "construct", "control", "Ppargc1a_cpm_arm",
            "Ppargc1a_cpm_ctrl", "ox_rel_diff", "ox_rel_lo", "ox_rel_hi",
            "mitopps_diff", "mitopps_lo", "mitopps_hi")],
  here::here("outputs", "orthotopic", "57_fig3i_route_comparison.csv"),
  row.names = FALSE)

p3i_dose <- ggplot2::ggplot(fig3i, ggplot2::aes(label, Ppargc1a_cpm_arm)) +
  ggplot2::geom_col(width = 0.5, fill = "grey35") +
  # %.4g, not %.3g: the record quotes 301.6 and 128.6, and %.3g rounds them to
  # 302 and 129, so a legend copied off this panel would not match the notes
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.4g", Ppargc1a_cpm_arm)),
                     vjust = -0.5, size = 2.5) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
  ggplot2::labs(x = NULL, y = "Ppargc1a retained (CPM, arm median)") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

p3i_ox <- ggplot2::ggplot(fig3i, ggplot2::aes(label, ox_rel_diff)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ox_rel_lo, ymax = ox_rel_hi),
                         width = 0.12, linewidth = 0.4) +
  ggplot2::geom_point(ggplot2::aes(shape = ox_excl_0), size = 2.8, fill = "white") +
  ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
  ggplot2::labs(x = NULL,
                y = "relative OXPHOS share\n(difference vs matched control)") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

CAP_3i <- paste0(
  "Each arm against its OWN matched control, within series: ",
  paste(sprintf("%s vs %s", fig3i$arm, fig3i$control), collapse = "; "),
  ". Arms are ordered by ROUTE, not by construct name. Ppargc1a RETAINED is an ",
  "outcome of selection, not an assigned dose: PGC-1a kills MYAZ cells, so each ",
  "arm is the population that survived expressing it. Filled points exclude ",
  "zero. The vector arms are a polyclonal pool; the CRISPR arms are POOLED ",
  "CRISPR without single-cell selection, so the p19ARF arm is disrupted rather ",
  "than null. n = 7 per arm; intervals are ", format(NBOOT, big.mark = ","),
  "-resample bootstraps of the difference in medians. Transcript associations ",
  "across escape routes, not a regulatory claim.")

fig_3i <- patchwork::wrap_plots(p3i_dose, p3i_ox, ncol = 1, heights = c(1, 1.1)) +
  patchwork::plot_annotation(
    title = "Fig. 3i -- the route comparison: what each escape route retained, and its respiratory share",
    subtitle = "ordered by what was done to the death pathway; the construct is beneath each route",
    caption = paste(strwrap(CAP_3i, width = 140), collapse = "\n"),
    tag_levels = list(c("", ""))) &
  ggplot2::theme(plot.caption = ggplot2::element_text(size = 6.5, hjust = 0,
                                                      colour = "grey25"))
ggplot2::ggsave(here::here("outputs", "orthotopic", "57_fig3i_route_comparison.pdf"),
                fig_3i, width = 8, height = 7.5)

# =============================================================================
# PART 5: DOES Bbc3 TRACK Ppargc1a WITHIN THE Pgc1a_BclxL ARM?  (spec Step 5a)
# -----------------------------------------------------------------------------
# This decides how Extended Data a is DESCRIBED, not whether it is drawn. If Bbc3
# tracks dose across the seven animals of one arm, the panel is a dose-response.
# If it is flat, the panel is a three-point comparison between three backgrounds
# that differ in genotype as well as in dose, and the legend must say that rather
# than implying a within-tumour dose-response.
#
# n = 7. Report it whichever way it comes out, and read the INTERVAL, not the
# point estimate: at this n a correlation CI spans most of the range, so "flat"
# and "strong" are usually both inside it. Both coefficients are given because
# Pearson answers the linear question the panel implies and Spearman is the one
# robust to the single high-dose animal.
# =============================================================================
message("57 PART 5: within-arm Bbc3 vs Ppargc1a")

ka <- per_animal$arm == "Pgc1a_BclxL"
xa <- per_animal$Ppargc1a_log2norm[ka]      # log scale: the panel's linear question
ya <- per_animal$Bbc3_log2norm[ka]
n_a <- sum(ka); stopifnot(n_a == 7L)

pear <- stats::cor.test(xa, ya, method = "pearson")
spr  <- suppressWarnings(stats::cor.test(xa, ya, method = "spearman"))
bs   <- vapply(seq_len(NBOOT), function(i) {
  j <- sample.int(n_a, replace = TRUE)
  if (length(unique(xa[j])) < 3 || length(unique(ya[j])) < 3) return(NA_real_)
  suppressWarnings(stats::cor(xa[j], ya[j], method = "spearman")) }, 0)
spr_ci <- stats::quantile(bs, c(0.025, 0.975), na.rm = TRUE)

within_arm <- list(
  arm = "Pgc1a_BclxL", n = n_a,
  pearson_r = unname(pear$estimate),
  pearson_ci = unname(pear$conf.int), pearson_p = unname(pear$p.value),
  spearman_rho = unname(spr$estimate),
  spearman_boot_ci = unname(spr_ci), spearman_p = unname(spr$p.value),
  scale = "log2(DESeq2-normalised + 1) on both axes")
spans <- within_arm$pearson_ci[1] < 0 && within_arm$pearson_ci[2] > 0
note("PART 5 (spec Step 5a): within Pgc1a_BclxL, n = 7, Bbc3 vs Ppargc1a -- ",
     sprintf("Pearson r %+.3f [%+.3f, %+.3f], p %.3f; Spearman rho %+.3f [%+.3f, %+.3f].",
             within_arm$pearson_r, within_arm$pearson_ci[1], within_arm$pearson_ci[2],
             within_arm$pearson_p, within_arm$spearman_rho,
             within_arm$spearman_boot_ci[1], within_arm$spearman_boot_ci[2]),
     if (spans)
       paste0(" THE INTERVAL SPANS ZERO. Extended Data a is therefore a ",
              "THREE-POINT COMPARISON BETWEEN THREE BACKGROUNDS that differ in ",
              "genotype as well as in retained dose, NOT a within-tumour ",
              "dose-response. The legend must say so. At n = 7 this is weak ",
              "evidence of absence, so the right wording is that the data do not ",
              "resolve a within-arm relation, not that there is none.")
     else
       paste0(" The interval excludes zero, so Bbc3 does track retained dose ",
              "WITHIN one background and Extended Data a may be described as a ",
              "dose-response -- still as selection, not regulation, and still at n = 7."))

# =============================================================================
# PART 6: EXTENDED DATA a -- Bbc3 AGAINST RETAINED DOSE, WITH THE SPECIFICITY
#         STRIP  (spec Step 4c)
# -----------------------------------------------------------------------------
# Bbc3 across the three backgrounds INCLUDING the arm where the transgene was
# lost, ordered by retained dose. The x axis is an ordered factor carrying the
# CPM, as script 53 draws it: the three doses are 0.17, 128.6 and 301.6, so a
# linear axis would stack two points on top of each other and a log axis has no
# zero for an arm that lost the transgene.
#
# THE STRIP IS A SPECIFICITY CONTROL, NOT A RESULT. All four canonical p53
# targets are labelled; Cdkn1a is NOT dropped, although it is also in 3f falling
# in the in vitro surviving fraction. Down in culture, flat in tumours, is what
# the narrative predicts, because those arms differ in buffer status.
#
# Cdkn2a IS SHOWN FOR THE VECTOR AND NT ARMS ONLY (author, 2026-10-28). The locus
# is CRISPR-targeted in the KO arm, so its value there is not a read-out of p53
# output; mouse Cdkn2a also pools p16INK4a and p19ARF from alternative first
# exons, which gene-level counts cannot separate. The legend says so.
# =============================================================================
message("57 PART 6: Extended Data a")

DOSE_ARMS <- c("Pgc1a_BclxL", "NTPgc1a", "KOPgc1a")
dose_of   <- stats::setNames(ET$Ppargc1a_CPM, ET$arm)
ord       <- DOSE_ARMS[order(dose_of[DOSE_ARMS])]
dose_lab  <- stats::setNames(
  sprintf("%s\n%.4g CPM retained", ord, dose_of[ord]), ord)

grab <- function(gene, arms) {
  d <- P5[P5$gene == gene & P5$hi %in% arms, ]
  d <- d[match(arms[arms %in% d$hi], d$hi), ]
  tibble::tibble(gene = gene, arm = d$hi, control = d$lo,
                 diff = d$diff_median, lo = d$ci_lo, hi = d$ci_hi,
                 wilcox_p = d$wilcox_p, direction = d$direction) }

eda_bbc3 <- grab("Bbc3", ord)
eda_bbc3$label <- factor(dose_lab[eda_bbc3$arm], levels = unname(dose_lab))
eda_bbc3$excl  <- eda_bbc3$hi < 0 | eda_bbc3$lo > 0

P53_FOUR <- c("Cdkn1a", "Mdm2", "Ccng1", "Zmat3")
eda_strip <- dplyr::bind_rows(
  dplyr::bind_rows(lapply(P53_FOUR, grab, arms = ord)),
  grab("Cdkn2a", ord[ord != "KOPgc1a"]))          # KO arm omitted: locus targeted
eda_strip$label <- factor(dose_lab[eda_strip$arm], levels = unname(dose_lab))
eda_strip$gene  <- factor(eda_strip$gene, levels = c("Cdkn2a", P53_FOUR))
eda_strip$excl  <- eda_strip$hi < 0 | eda_strip$lo > 0

COLS <- c("gene", "arm", "control", "diff", "lo", "hi", "wilcox_p", "direction")
utils::write.csv(
  dplyr::bind_rows(cbind(panel = "a (Bbc3)",          eda_bbc3[, COLS]),
                   cbind(panel = "b (specificity)",   eda_strip[, COLS])),
  here::here("outputs", "orthotopic", "57_eda_bbc3_dose.csv"), row.names = FALSE)

pE1 <- ggplot2::ggplot(eda_bbc3, ggplot2::aes(label, diff)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.4) +
  ggplot2::geom_point(ggplot2::aes(shape = excl), size = 2.8, fill = "white") +
  ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
  ggplot2::labs(x = "Ppargc1a retained", y = "Bbc3, log2 difference\nvs matched control") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

# the strip repeats the same three arms five times, so it takes the DOSE ALONE as
# its tick label -- the full two-line labels collide across five facets and the
# axis title carries the unit. Panel a keeps the long form.
short_lab <- stats::setNames(sprintf("%.4g", dose_of[ord]), unname(dose_lab))
eda_strip$short <- factor(short_lab[as.character(eda_strip$label)],
                          levels = unname(short_lab))

# the KO arm is OMITTED from Cdkn2a, not measured-and-flat: say so on the panel,
# or an empty slot reads as a null result
omitted <- tibble::tibble(
  gene  = factor("Cdkn2a", levels = levels(eda_strip$gene)),
  short = factor(short_lab[[unname(dose_lab[["KOPgc1a"]])]], levels = levels(eda_strip$short)),
  diff  = 0)

pE2 <- ggplot2::ggplot(eda_strip, ggplot2::aes(short, diff)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  ggplot2::geom_text(data = omitted, ggplot2::aes(short, diff, label = "locus\ntargeted,\nnot shown"),
                     size = 1.9, colour = "grey45", lineheight = 0.95) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.35) +
  ggplot2::geom_point(ggplot2::aes(shape = excl), size = 2, fill = "white") +
  ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
  ggplot2::facet_wrap(~ gene, nrow = 1) +
  ggplot2::labs(x = "Ppargc1a retained (CPM, arm median)",
                y = "log2 difference\nvs matched control") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(size = 6))

CAP_ED <- paste0(
  "Bbc3 against the PGC-1a dose each arm RETAINED, three backgrounds, each ",
  "against its own matched control (", paste(sprintf("%s vs %s", eda_bbc3$arm,
  eda_bbc3$control), collapse = "; "), "). Retention is an outcome of selection, ",
  "not an assigned dose. The strip is a SPECIFICITY control: four canonical p53 ",
  "targets, all labelled. Cdkn2a is shown for the vector and NT arms only -- the ",
  "locus is CRISPR-targeted in the p19ARF arm, so its value there is not a ",
  "read-out of p53 output, and mouse Cdkn2a pools p16INK4a and p19ARF from ",
  "alternative first exons, which gene-level counts cannot separate. Filled ",
  "points exclude zero. n = 7 per arm; intervals are ",
  format(NBOOT, big.mark = ","), "-resample bootstraps. Script 52's R3 found one ",
  "of five panel genes moving with Bbc3, so the fall is unlikely to reflect a ",
  "change in p53 transcriptional output. Association across escape routes. ",
  if (spans)
    paste0("THREE BACKGROUNDS COMPARED AT THREE DOSES, NOT A DOSE-RESPONSE: within ",
           "Pgc1a_BclxL, Bbc3 does not track Ppargc1a across its seven animals ",
           sprintf("(Pearson r %+.3f [%+.3f, %+.3f]), ", within_arm$pearson_r,
                   within_arm$pearson_ci[1], within_arm$pearson_ci[2]),
           "so the three points differ in genotype as well as in retained dose. At ",
           "n = 7 that is weak evidence of absence, not a demonstration of none.")
  else
    paste0("Bbc3 also tracks retained dose WITHIN one background ",
           sprintf("(Pearson r %+.3f [%+.3f, %+.3f] in Pgc1a_BclxL), ",
                   within_arm$pearson_r, within_arm$pearson_ci[1], within_arm$pearson_ci[2]),
           "so the across-arm pattern is not carried by background alone."))

fig_ed <- patchwork::wrap_plots(pE1, pE2, ncol = 1, heights = c(1, 0.95)) +
  patchwork::plot_annotation(
    title = "Extended Data -- Bbc3 is lower in the arms that retained more PGC-1a, and the p53 programme does not follow",
    caption = paste(strwrap(CAP_ED, width = 140), collapse = "\n"),
    tag_levels = list(c("a", "b"))) &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 10, face = "bold"),
                 plot.caption = ggplot2::element_text(size = 6.5, hjust = 0,
                                                      colour = "grey25"))
ggplot2::ggsave(here::here("outputs", "orthotopic", "57_eda_bbc3_dose.pdf"),
                fig_ed, width = 9, height = 7.5)

# =============================================================================
# PART 7: THE COMPOSITION SCREEN, REPORTED AS IT STANDS  (spec Step 5b)
# -----------------------------------------------------------------------------
# NO FITS. Script 52 screened six composition axes against three respiratory
# rulers and adjusted for those reaching |rho| >= 0.5 pooled. NOTHING REACHED IT,
# so `material_axes` is empty, `AD` was NULL, and every estimate in the record --
# including the Ppargc1a positive control, whose `adjusted_for` column reads
# "nothing" -- is UNADJUSTED. There is no adjusted side to compare against.
#
# WHERE THE THRESHOLD WAS PRE-SPECIFIED. `CUT_MAT <- 0.5` sits in script 52's
# parameter block (line 101), above every computation, and is saved as
# `$params$CUT_MAT`. It carried over from the 21-sample vector-series run, where
# `docs/2026-09-08_orthotopic_vector_series.md:100` calls it "the pre-registered
# |rho| >= 0.5". It is NOT in the 2026-09-07 analysis plan.
#   NOT TO BE CONFUSED with the other 0.5 in script 52, at line 415, which turns
#   R1's ox_lvl difference into the word "SMALL". The escape-series note flags
#   THAT one as chosen in code and not from the spec. Different threshold,
#   different purpose.
# =============================================================================
message("57 PART 7: composition screen, report only")

GT <- as.data.frame(res$gate)
gate_ox_rel <- GT[GT$ruler == "ox_rel", ]
strongest <- gate_ox_rel[which.max(abs(gate_ox_rel$rho_pooled)), ]
arm_cols  <- paste0("rho_", SET)

composition_screen <- list(
  threshold = res$params$CUT_MAT,
  threshold_source = paste(
    "scripts/52_orthotopic_escape_series.R line 101 (parameter block, above every",
    "computation); saved as $params$CUT_MAT; called pre-registered at",
    "docs/2026-09-08_orthotopic_vector_series.md:100. Not in the 2026-09-07 plan."),
  material_axes = res$material_axes,
  n_material = length(res$material_axes),
  adjustment_applied = length(res$material_axes) > 0,
  gate_control = as.data.frame(res$gate_control),
  strongest_axis = strongest$compartment,
  strongest_rho_pooled = strongest$rho_pooled,
  strongest_by_arm = stats::setNames(as.numeric(strongest[arm_cols]), SET),
  full_ox_rel_table = gate_ox_rel)

note("PART 7 (spec Step 5b, report only): no composition axis reached ",
     "|rho| >= ", res$params$CUT_MAT, " pooled on ox_rel, so NO ADJUSTMENT WAS ",
     "APPLIED to anything in the record and there is no adjusted estimate for ",
     "ox_rel or Bbc3 to compare against. The strongest axis is ",
     strongest$compartment, ", pooled ", sprintf("%+.3f", strongest$rho_pooled),
     ", by arm ",
     paste(sprintf("%s %+.3f", SET, as.numeric(strongest[arm_cols])), collapse = ", "),
     ". ALL 49-SAMPLE VALUES. The spec's adipose figures (+0.621 pooled, +0.036 / ",
     "-0.393 / +0.714) are the 21-SAMPLE vector run and their three columns are ",
     "EV / BclxL / Pgc1a_BclxL, not the three PGC-1a arms; 49-sample adipose on ",
     "ox_rel is ",
     sprintf("%+.3f pooled, Pgc1a_BclxL %+.3f, NTPgc1a %+.3f, KOPgc1a %+.3f",
             gate_ox_rel$rho_pooled[gate_ox_rel$compartment == "adipose"],
             gate_ox_rel$rho_Pgc1a_BclxL[gate_ox_rel$compartment == "adipose"],
             gate_ox_rel$rho_NTPgc1a[gate_ox_rel$compartment == "adipose"],
             gate_ox_rel$rho_KOPgc1a[gate_ox_rel$compartment == "adipose"]),
     ". Note the trap: NTPgc1a's +0.036 happens to equal EV's value in the ",
     "21-sample table.")

# =============================================================================
# PART 8: Myc IN KOPgc1a AGAINST KOEV, AND Myc ACROSS ALL SEVEN ARMS (Step 5c)
# -----------------------------------------------------------------------------
# NEITHER the 80th-percentile NOR the maximum contrast. Script 52 computed both
# under an exact null, and KOPgc1a carries a Myc outlier that makes its MAXIMUM
# contrast positive and uninterpretable. The author asked for a centre with an
# interval instead, so this is a median difference with a bootstrap interval,
# built with script 52's recipe.
#
# Script 52 never fitted a median Myc contrast, so this is one of the three
# genuinely new numbers in this script -- not a re-fit of something on the record.
# The mean is given beside it because a mean is what gets written by accident,
# and the outlier-dropped mean shows exactly how much that choice costs.
# NO PANEL. The author decides once the number is visible.
# =============================================================================
message("57 PART 8: Myc")

myc_cpm <- as.numeric(cpm_of("Myc"))
boot_diff <- function(a, b, f = stats::median) {
  d <- vapply(seq_len(NBOOT), function(i)
    f(sample(a, replace = TRUE)) - f(sample(b, replace = TRUE)), 0)
  stats::quantile(d, c(0.025, 0.975), na.rm = TRUE) }

a_ko <- myc_cpm[sm$group == "KOPgc1a"]; b_ko <- myc_cpm[sm$group == "KOEV"]
ci_med  <- boot_diff(a_ko, b_ko, stats::median)
ci_mean <- boot_diff(a_ko, b_ko, mean)
out_i   <- which.max(a_ko)
a_no    <- a_ko[-out_i]
ci_mean_no <- boot_diff(a_no, b_ko, mean)

myc_ko <- list(
  contrast = "KOPgc1a vs KOEV, Myc CPM", n_hi = length(a_ko), n_lo = length(b_ko),
  median_diff = stats::median(a_ko) - stats::median(b_ko),
  median_ci = unname(ci_med),
  mean_diff = mean(a_ko) - mean(b_ko), mean_ci = unname(ci_mean),
  outlier_cpm = max(a_ko),
  mean_diff_outlier_dropped = mean(a_no) - mean(b_ko),
  mean_ci_outlier_dropped = unname(ci_mean_no),
  wilcox_p = suppressWarnings(stats::wilcox.test(a_ko, b_ko)$p.value))

myc_by_arm <- dplyr::bind_rows(lapply(SET, function(g) {
  v <- myc_cpm[sm$group == g]
  tibble::tibble(arm = g, n = length(v), Myc_median_cpm = stats::median(v),
                 Myc_mean_cpm = mean(v), Myc_min_cpm = min(v), Myc_max_cpm = max(v),
                 Myc_IQR_cpm = stats::IQR(v)) }))

med_spans <- myc_ko$median_ci[1] < 0 && myc_ko$median_ci[2] > 0
note("PART 8 (spec Step 5c): Myc in KOPgc1a vs KOEV, CPM -- median difference ",
     sprintf("%+.0f [%+.0f, %+.0f]", myc_ko$median_diff, myc_ko$median_ci[1],
             myc_ko$median_ci[2]),
     sprintf("; mean %+.0f [%+.0f, %+.0f]", myc_ko$mean_diff, myc_ko$mean_ci[1],
             myc_ko$mean_ci[2]),
     sprintf("; mean with the %.0f CPM outlier dropped %+.0f [%+.0f, %+.0f].",
             myc_ko$outlier_cpm, myc_ko$mean_diff_outlier_dropped,
             myc_ko$mean_ci_outlier_dropped[1], myc_ko$mean_ci_outlier_dropped[2]),
     if (med_spans)
       paste0(" THE MEDIAN INTERVAL SPANS ZERO, so by the spec's own rule (section 7) ",
              "the Myc ceiling stays what the September note called it -- the weakest ",
              "of the three sequenced results -- and ED d is NOT built.")
     else
       paste0(" The median interval excludes zero, so Myc does come down under ",
              "PGC-1a pressure in the p19ARF-disrupted background as well. By the ",
              "spec's section 7 that makes ED d buildable, at the cost of saying ",
              "CEILING rather than abolition, precisely, because 3c has just ",
              "established Myc unchanged in the in vitro survivors."))

# =============================================================================
# PART 9: SAVE, AND WRITE THE NOTE
# -----------------------------------------------------------------------------
# The note is generated from the objects above rather than transcribed, so the
# numbers in it cannot drift from the numbers in the panels.
# =============================================================================
message("57 PART 9: save and write the note")

out <- list(
  unit = list(levels = "CPM",
              internal_contrast_scale = "log2(DESeq2-normalised + 1)",
              reason = paste("script 52, its object and both orthotopic notes quote CPM;",
                             "DESeq2-normalised differs by a near-constant 26.4-28.3x")),
  rebuild_check = chk,
  per_animal = per_animal,
  fig3g = fig3g, fig3g_prism_wide = wide,
  summary_numbers = summary_numbers, note_records = NOTE_RECORDS,
  manifest_check = manifest_check,
  construct_by_arm = construct_by_arm, construct_outliers = construct_outliers,
  fig3i = fig3i, eda_bbc3 = eda_bbc3, eda_strip = eda_strip,
  within_arm_bbc3_ppargc1a = within_arm,
  composition_screen = composition_screen,
  myc_ko = myc_ko, myc_by_arm = myc_by_arm,
  crispr_derivation = paste(
    "POOLED CRISPR, no single-cell selection, so the p19ARF arm is DISRUPTED,",
    "not null. Edited fraction checked and p19ARF western QC passed (author,",
    "2026-10-28). data/orthotopic_series/README.md corrected; script 52's header",
    "and script 53's caption still say 'clone' and are deliberately not edited --",
    "they are the committed analysis behind a saved object and a published PDF."),
  limitations = c(
    "Every arm is a survivor population, so a PGC-1a effect measured here is an escape-route signature, not a response to a manipulation.",
    "No transcript-level quantification exists, so gene-level Bcl2l1 pools Bcl-xL with Bcl-xS, which are functionally opposite products of one locus; the transcript-length offset is lost for every number.",
    "Bcl2l1 is construct-contaminated in every Bcl-xL-carrying arm, so NTPgc1a vs NTEV is the only readable guardian contrast.",
    "The retained dose is an outcome of selection: label it retained, never assigned.",
    "The intervention moves the composition markers, so pooled loadings mislead; read the per-arm breakdown."),
  source_spec = "docs/2026-10-27_section3_final_narrative_and_panel_spec.md (v2) + author build list 2026-10-28",
  params = list(NBOOT = NBOOT, MINBM = MINBM, CPM_ON = CPM_ON, seed = 1),
  analysis_date = Sys.Date(), notes = NOTES)

saveRDS(out, here::here("results", "orthotopic_section3_panels.rds"))

fmt_ci <- function(v, lo, hi, d = 3) sprintf("%+.*f [%+.*f, %+.*f]", d, v, d, lo, d, hi)
md <- c(
  "---",
  "date: 2026-10-28",
  "tags: [project/myc_mouse, nat-metab, section-3, figure-3, orthotopic, panel-build]",
  "status: generated by scripts/57_section3_panel_readout.R -- do not hand-edit the numbers",
  "relates-to:",
  "  - scripts/57_section3_panel_readout.R (writes this file)",
  "  - results/orthotopic_section3_panels.rds",
  "  - docs/2026-10-27_section3_final_narrative_and_panel_spec.md (v2, plus the author's 2026-10-28 build list)",
  "  - docs/2026-09-23_section3_narrative_and_figure3.md section 2.5 (limitations)",
  "  - scripts/52_orthotopic_escape_series.R at faf1e8a; results/orthotopic_escape_series.rds",
  "---",
  "",
  "# Section 3 panel read-out",
  "",
  sprintf("Generated %s. Every contrast is re-read from script 52's saved object; the", Sys.Date()),
  "per-animal layer is rebuilt with script 52's recipe and asserted against its saved arm",
  sprintf("medians (largest disagreement %.1e).", max(chk$max_abs_diff)),
  "",
  "## 1. The unit",
  "",
  "**CPM for every expression level.** Script 52, its object and both orthotopic notes quote",
  "CPM; DESeq2-normalised counts differ by a near-constant 26.4-28.3x here, so the ordering is",
  "identical and a switch would only renumber four documents. The Prism axis reading \"DESeq2",
  "normalised\" is relabelled to CPM.",
  "",
  "Script 52's gene contrasts are differences of medians of log2(DESeq2-normalised + 1), so a",
  "fold change read off the CPM columns will not equal a saved contrast exactly. Both scales are",
  "exported (`*_cpm`, `*_log2norm`).",
  "",
  "## 2. Arm identity",
  "",
  sprintf("Seven arms, n = 7 each, 49 sequenced. Manifest verification: **%s**.", manifest_check$status),
  if (!is.null(manifest_check$detail)) paste0("", manifest_check$detail) else "",
  "",
  "The CRISPR arms are **pooled CRISPR with no single-cell selection**, so the p19ARF arm is",
  "**disrupted, not null**. Editing was validated independently: the edited fraction was checked",
  "and the p19ARF western QC passed. `data/orthotopic_series/README.md` is corrected; script 52's",
  "header and script 53's caption still say \"clone\" and are deliberately left alone, because they",
  "are the committed analysis behind a saved object and a published PDF.",
  "",
  "## 3. Summary numbers, each arm against its own control",
  "",
  "| route | arm | control | Ppargc1a CPM (arm) | ox_rel vs control | mitoPPS vs control | MitoCarta share, arm vs control | Bbc3 log2 vs control |",
  "|---|---|---|---|---|---|---|---|",
  apply(summary_numbers, 1, function(r) sprintf(
    "| %s | `%s` | `%s` | %.4g | %s | %s | %.1f%% vs %.1f%% | %s |",
    r[["route"]], r[["arm"]], r[["control"]], as.numeric(r[["Ppargc1a_cpm_arm"]]),
    fmt_ci(as.numeric(r[["ox_rel_diff"]]), as.numeric(r[["ox_rel_lo"]]), as.numeric(r[["ox_rel_hi"]])),
    fmt_ci(as.numeric(r[["mitopps_diff"]]), as.numeric(r[["mitopps_lo"]]), as.numeric(r[["mitopps_hi"]])),
    100 * as.numeric(r[["mito_share_arm"]]), 100 * as.numeric(r[["mito_share_ctrl"]]),
    fmt_ci(as.numeric(r[["Bbc3_diff"]]), as.numeric(r[["Bbc3_lo"]]), as.numeric(r[["Bbc3_hi"]])))),
  "",
  sprintf("Intervals are %s-resample bootstraps of the difference in medians, n = 7 per arm.",
          format(NBOOT, big.mark = ",")),
  sprintf("Agreement with the escape-series note: %s.",
          if (all(NOTE_RECORDS$agrees)) "all five recorded values reproduce to rounding"
          else sprintf("**%d of %d DISAGREE -- the object is right**", sum(!NOTE_RECORDS$agrees), nrow(NOTE_RECORDS))),
  "",
  "## 4. Fig. 3g, the P value",
  "",
  sprintf("Mann-Whitney on the exported per-animal CPM: **P = %s**, against a two-sided floor of",
          signif(fig3g$p, 3)),
  sprintf("%s at 7 vs 7. The placeholder \"P < 0.0001\" is not reachable at this n and comes off the",
          signif(mw_floor, 3)),
  sprintf("panel; it is not carrying the result, which is a %d-fold median transgene difference.",
          round(fig3g$fold)),
  "No DESeq2 Wald test is fitted anywhere: script 52 never ran `DESeq()`, so there is none to cite.",
  "",
  "## 5. The three read-outs the spec asked for",
  "",
  "**Within `Pgc1a_BclxL`, does Bbc3 track Ppargc1a? (Step 5a)**",
  sprintf("n = 7. Pearson r %s, P = %.3f. Spearman rho %s.",
          fmt_ci(within_arm$pearson_r, within_arm$pearson_ci[1], within_arm$pearson_ci[2]),
          within_arm$pearson_p,
          fmt_ci(within_arm$spearman_rho, within_arm$spearman_boot_ci[1],
                 within_arm$spearman_boot_ci[2])),
  if (spans)
    paste("The interval spans zero, so Extended Data a is a **three-point comparison between three",
          "backgrounds** that differ in genotype as well as in retained dose, not a within-tumour",
          "dose-response, and the legend must say so. At n = 7 this is weak evidence of absence:",
          "the data do not resolve a within-arm relation, which is not the same as there being none.")
  else
    paste("The interval excludes zero, so Bbc3 tracks retained dose within one background and",
          "Extended Data a may be described as a dose-response -- as selection, not regulation."),
  "",
  "**Composition adjustment (Step 5b, report only)**",
  sprintf("No axis reached |rho| >= %.1f pooled on `ox_rel`, so **no adjustment was applied to anything**",
          composition_screen$threshold),
  "in the record, and there is no adjusted estimate for `ox_rel` or `Bbc3` to compare against.",
  sprintf("Strongest axis: **%s**, pooled %+.3f; by arm %s.",
          composition_screen$strongest_axis, composition_screen$strongest_rho_pooled,
          paste(sprintf("`%s` %+.3f", SET, composition_screen$strongest_by_arm), collapse = ", ")),
  sprintf("Threshold provenance: %s", composition_screen$threshold_source),
  "",
  "All 49-sample values. The spec's adipose figures are the 21-sample vector run and its three",
  "columns are EV / BclxL / Pgc1a_BclxL, not the three PGC-1a arms; 49-sample adipose on `ox_rel`",
  sprintf("is %+.3f pooled, `Pgc1a_BclxL` %+.3f, `NTPgc1a` %+.3f, `KOPgc1a` %+.3f.",
          gate_ox_rel$rho_pooled[gate_ox_rel$compartment == "adipose"],
          gate_ox_rel$rho_Pgc1a_BclxL[gate_ox_rel$compartment == "adipose"],
          gate_ox_rel$rho_NTPgc1a[gate_ox_rel$compartment == "adipose"],
          gate_ox_rel$rho_KOPgc1a[gate_ox_rel$compartment == "adipose"]),
  "",
  "**Myc in `KOPgc1a` vs `KOEV` (Step 5c)**",
  sprintf("Median difference %s CPM; mean %s; mean with the %.0f CPM outlier dropped %s.",
          fmt_ci(myc_ko$median_diff, myc_ko$median_ci[1], myc_ko$median_ci[2], 0),
          fmt_ci(myc_ko$mean_diff, myc_ko$mean_ci[1], myc_ko$mean_ci[2], 0),
          myc_ko$outlier_cpm,
          fmt_ci(myc_ko$mean_diff_outlier_dropped, myc_ko$mean_ci_outlier_dropped[1],
                 myc_ko$mean_ci_outlier_dropped[2], 0)),
  if (med_spans)
    "The median interval spans zero, so by the spec's section 7 the Myc ceiling stays the weakest of the three sequenced results and **ED d is not built**."
  else
    "The median interval excludes zero, so by the spec's section 7 ED d is buildable, at the cost of writing ceiling rather than abolition.",
  "",
  "## 6. For Methods",
  "",
  paste("-", out$crispr_derivation),
  paste("-", out$limitations),
  paste("- Composition markers were screened against three respiratory rulers;",
        sprintf("no axis reached |rho| >= %.1f pooled, so no adjustment was applied.",
                composition_screen$threshold)),
  "",
  "## 7. What this script does not do",
  "",
  "- No DESeq2 Wald fit, here or anywhere in the orthotopic cohort.",
  "- No composition-adjusted re-fits: there is no adjusted side to report.",
  "- No growth panel and no blot panel; both are built in Prism from data outside the repo.",
  "- No `Cdkn2a` in the KO arm: the locus is CRISPR-targeted there.",
  "- No value is quoted across the 21-sample and 49-sample scoring runs.",
  "",
  "## 8. Run notes",
  "",
  paste("-", NOTES))
writeLines(md, here::here("docs", "2026-10-28_section3_panel_readout.md"))

message("57: wrote results/orthotopic_section3_panels.rds")
message("57: wrote docs/2026-10-28_section3_panel_readout.md")
message("57: wrote 4 CSVs and 2 PDFs under outputs/orthotopic/")

cat("\n================ SUMMARY NUMBERS (spec Step 3) ================\n")
print(as.data.frame(summary_numbers), row.names = FALSE, digits = 3)
cat("\n================ FIG 3g ================\n")
cat(sprintf("Mann-Whitney P = %s (two-sided floor at 7v7 = %s); medians %.3g vs %.3g CPM, %d-fold\n",
            signif(fig3g$p, 3), signif(fig3g$floor_at_7v7, 3),
            fig3g$median_NTPgc1a, fig3g$median_NTEV, round(fig3g$fold)))
cat("\n================ CONSTRUCT CHECK ================\n")
print(as.data.frame(construct_by_arm), row.names = FALSE, digits = 4)
cat("\n================ Myc ACROSS THE SEVEN ARMS ================\n")
print(as.data.frame(myc_by_arm), row.names = FALSE, digits = 5)
cat("\n================ NOTES ================\n")
cat(paste(strwrap(NOTES, width = 100, exdent = 2), collapse = "\n"), sep = "\n")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  out <- readRDS(here::here("results", "orthotopic_section3_panels.rds"))

  ## 1. THE REBUILD GATE, before anything else. If these are not ~0 the
  ## per-animal layer is on a different ruler from script 52 and nothing below
  ## it may be compared with the record.
  out$rebuild_check |> print()

  ## 2. THE PRISM EXPORT. The wide one is what gets pasted into Fig 3g; the P
  ## value is what Prism should reproduce on those same numbers.
  out$fig3g_prism_wide |> print()
  out$fig3g |> print()

  ## 3. THE SUMMARY TABLE, and whether it agrees with the escape-series note.
  out$summary_numbers |> as.data.frame() |> print(digits = 3)
  out$note_records |> as.data.frame() |> print(digits = 3)

  ## 4. IDENTITY. The manifest check first -- if it says NOT PRESENT, the arm
  ## assignment still rests on the column headers alone.
  out$manifest_check |> print()
  out$construct_by_arm |> as.data.frame() |> print(digits = 4)
  out$construct_outliers |> as.data.frame() |> print()

  ## 5. THE DECISION THAT DESCRIBES EXTENDED DATA a: dose-response within one
  ## arm, or three backgrounds compared at three doses?
  out$within_arm_bbc3_ppargc1a |> print()

  ## 6. THE COMPOSITION SCREEN. material_axes is expected to be EMPTY -- that is
  ## the result, not a failure, and it is why there is no adjusted estimate.
  out$composition_screen$material_axes |> print()
  out$composition_screen$full_ox_rel_table |> as.data.frame() |> print(digits = 3)

  ## 7. Myc. Read the MEDIAN with its interval; the max and 80th-percentile
  ## contrasts are uninterpretable in KOPgc1a because of the outlier.
  out$myc_ko |> print()
  out$myc_by_arm |> as.data.frame() |> print(digits = 5)

  ## 8. per-animal values, for any panel not built here
  out$per_animal |> print(n = 49)

  cat(out$notes, sep = "\n\n")
}
