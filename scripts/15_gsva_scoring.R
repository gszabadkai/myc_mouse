# scripts/15_gsva_scoring.R
# =============================================================================
# GSVA scoring engine: per-sample gene-set scores for the mammary geneset
# library (Block A, Day 2)
# =============================================================================
#
# Motivation:
#   fGSEA (scripts 04/10/19) reads a ranked contrast and asks "is this set
#   enriched at the top". GSVA asks a different question: for EACH sample,
#   how active is each program? That per-sample matrix is what AP1/AP2/AP3
#   (dev composition, script 18), AP4 (Felsher co-variation), and AP7 (the
#   MB-fork projection, script 17) all build on. This script is the shared
#   engine; downstream scripts subset the matrix by category.
#
# SCALE TRAP (read before touching this file):
#   GSVA takes LOG-SCALE input (here VST) with kcdf = "Gaussian". This is the
#   OPPOSITE of mitoPPS (script 08), which is a ratio method on LINEAR-scale
#   DESeq2 normalised counts. Do not feed linear counts here, and do not feed
#   VST into mitoPPS. See CLAUDE.md "Key methods and analytical conventions".
#
# Scope (decided 2026-07-06): score ALL GSVA-tagged sets (method 'gsva' or
#   'both' in provenance_table.csv == 902 sets) in ONE cohort-relative run.
#   fGSEA (script 19) covers every category as ranked enrichment; this GSVA
#   matrix is the per-sample lens. Scoring the full 902 makes it a reusable
#   superset so later opportunistic uses (08_apoptosis -> H4/AP8,
#   06_tf_targets -> AP3, 04/05 -> Fig 2) need no re-run.
#
# Sets are already MOUSE symbols (by_category/*_mouse.gmt); the VST matrix is
#   keyed on ENSMUSG ensembl IDs. We map ensembl -> mouse symbol via the cached
#   ortholog table's LOCAL symbol columns (ensembl_gene_id, external_gene_name),
#   NOT the human-ortholog columns (no human mapping is needed here).
#
# Input:
#   - results/dds_int_run.rds                          (DESeqDataSet, 24 samples)
#   - results/ortholog_table.rds                       (ensembl -> mouse symbol)
#   - data/genesets_from_library/by_category/*_mouse.gmt   (9 category GMTs)
#   - data/genesets_from_library/provenance_table.csv  (method tag per set)
#
# Output:
#   - results/gsva_scores.rds : list(scores, set_meta, sample_meta,
#                                     n_sets_scored, n_symbols_collapsed)
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# =============================================================================
# PART 1: LOAD INPUTS
# =============================================================================

dds_int   <- readRDS(here::here("results", "dds_int_run.rds"))
ortholog  <- readRDS(here::here("results", "ortholog_table.rds"))
provenance <- readr::read_csv(
  here::here("data", "genesets_from_library", "provenance_table.csv"),
  show_col_types = FALSE
)

gmt_dir  <- here::here("data", "genesets_from_library", "by_category")
gmt_files <- list.files(gmt_dir, pattern = "_mouse\\.gmt$", full.names = TRUE)
stopifnot(length(gmt_files) == 9L)

message(sprintf("Loaded dds (%d genes x %d samples), %d ortholog rows, %d provenance rows",
                nrow(dds_int), ncol(dds_int), nrow(ortholog), nrow(provenance)))

# =============================================================================
# PART 2: VST EXPRESSION MATRIX (log scale) -> mouse-symbol rows
# =============================================================================
# VST is the GSVA-appropriate log-scale transform. blind = FALSE keeps the
# design-aware dispersion trend (we are not using VST for unsupervised QC here).

vsd      <- DESeq2::vst(dds_int, blind = FALSE)
vst_mat  <- SummarizedExperiment::assay(vsd)          # ENSMUSG x 24 (log scale)

# ensembl -> mouse symbol map (local columns; dedup, drop blanks)
ens2sym <- ortholog |>
  dplyr::select(ensembl_gene_id, external_gene_name) |>
  dplyr::filter(!is.na(external_gene_name), external_gene_name != "") |>
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

vst_df <- tibble::as_tibble(vst_mat, rownames = "ensembl_gene_id") |>
  dplyr::inner_join(ens2sym, by = "ensembl_gene_id")

n_before_symbol <- nrow(vst_df)

# Collapse duplicate symbols (several ENSMUSG -> one mouse symbol): keep the
# highest-mean-VST row per symbol (most-expressed transcript wins).
sample_cols <- colnames(vst_mat)
vst_collapsed <- vst_df |>
  dplyr::mutate(row_mean = rowMeans(dplyr::across(dplyr::all_of(sample_cols)))) |>
  dplyr::arrange(external_gene_name, dplyr::desc(row_mean)) |>
  dplyr::distinct(external_gene_name, .keep_all = TRUE)

n_symbols_collapsed <- n_before_symbol - nrow(vst_collapsed)

expr_mat <- vst_collapsed |>
  dplyr::select(dplyr::all_of(sample_cols)) |>
  as.matrix()
rownames(expr_mat) <- vst_collapsed$external_gene_name

message(sprintf("VST matrix: %d ensembl -> %d unique mouse symbols (collapsed %d)",
                n_before_symbol, nrow(expr_mat), n_symbols_collapsed))

# =============================================================================
# PART 3: GENE SETS (GSVA-tagged only) FROM THE LIBRARY GMTs
# =============================================================================

# Union of all 9 category GMTs (mouse symbols). Named list: set_name -> genes.
all_pathways <- unlist(
  lapply(gmt_files, fgsea::gmtPathways),
  recursive = FALSE
)
# Guard against any accidental name collision across category files.
stopifnot(!any(duplicated(names(all_pathways))))

gsva_tags   <- provenance$set_name[provenance$method %in% c("gsva", "both")]
pathways    <- all_pathways[names(all_pathways) %in% gsva_tags]

message(sprintf("Pathways: %d total in GMTs, %d GSVA-tagged (method gsva/both)",
                length(all_pathways), length(pathways)))
stopifnot(length(pathways) > 0)

# =============================================================================
# PART 4: GSVA (one cohort-relative run, all 24 samples)
# =============================================================================
# Size filter: min.sz = 5, NO max cap -- the library's large-set size-filter
# exemption for GSVA applies (see CLAUDE.md / library decisions). GSVA drops
# sets with < min.sz genes present in expr_mat automatically.
#
# GSVA API split: >= 1.50 uses gsvaParam() + gsva(); older takes the legacy
# positional gsva(expr, gset, ...). kcdf = "Gaussian" for log-scale (VST) input
# either way. Detect at run time so this is robust across installs.

min_sz <- 5L

if ("gsvaParam" %in% getNamespaceExports("GSVA")) {
  gp <- GSVA::gsvaParam(
    exprData = expr_mat,
    geneSets = pathways,
    kcdf     = "Gaussian",
    minSize  = min_sz,
    maxSize  = Inf
  )
  gsva_scores <- GSVA::gsva(gp, verbose = TRUE)
} else {
  gsva_scores <- GSVA::gsva(
    expr    = expr_mat,
    gset.idx.list = pathways,
    kcdf    = "Gaussian",
    min.sz  = min_sz,
    max.sz  = Inf,
    verbose = TRUE
  )
}

# Preserve dds sample order for downstream joins.
gsva_scores <- gsva_scores[, colnames(expr_mat), drop = FALSE]

message(sprintf("GSVA scored %d sets x %d samples", nrow(gsva_scores), ncol(gsva_scores)))
stopifnot(!anyNA(gsva_scores))

# =============================================================================
# PART 5: METADATA + SAVE
# =============================================================================

set_meta <- provenance |>
  dplyr::filter(set_name %in% rownames(gsva_scores)) |>
  dplyr::select(set_name, category_primary, category, method, size_mouse) |>
  dplyr::arrange(match(set_name, rownames(gsva_scores)))

sample_meta <- as.data.frame(SummarizedExperiment::colData(dds_int))

gsva_out <- list(
  scores              = gsva_scores,
  set_meta            = set_meta,
  sample_meta         = sample_meta,
  n_sets_scored       = nrow(gsva_scores),
  n_symbols_collapsed = n_symbols_collapsed,
  input_scale         = "VST (log); kcdf=Gaussian; opposite of mitoPPS (linear)",
  # Also stored so downstream gene-level tests (script 17 ROAST/CAMERA +
  # cross-check) use EXACTLY the matrix + set lists that produced the scores.
  expr_mat            = expr_mat,   # symbol x 24 VST (full universe)
  pathways            = pathways    # GSVA-tagged mouse-symbol gene-set lists
)

saveRDS(gsva_out, here::here("results", "gsva_scores.rds"))
message("Saved results/gsva_scores.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  gsva_out <- readRDS(here::here("results", "gsva_scores.rds"))
  sc <- gsva_out$scores

  # Shape + sample order matches dds (24 cols, colData order)
  dim(sc)
  identical(colnames(sc), rownames(gsva_out$sample_meta))

  # No all-NA / degenerate rows; scores are the expected GSVA range (~[-1, 1])
  stopifnot(!anyNA(sc))
  range(sc)
  summary(as.vector(sc))
  # Any zero-variance sets (constant across samples) would be suspect
  sum(apply(sc, 1, function(x) stats::sd(x) == 0))

  # How many sets per category made it through the size filter
  gsva_out$set_meta |> dplyr::count(category_primary) |> print(n = Inf)
  cat("collapsed symbols:", gsva_out$n_symbols_collapsed, "\n")

  # Spot-check a known developmental set (MG_*) and a MYC set: score spread and
  # the expected Myc+ vs Myc- / timepoint separation.
  smeta <- gsva_out$sample_meta
  dev_set <- grep("^MG_", rownames(sc), value = TRUE)[1]
  myc_set <- grep("^MYC_", rownames(sc), value = TRUE)[1]
  data.frame(group = smeta$group,
             dev = sc[dev_set, ],
             myc = sc[myc_set, ]) |>
    dplyr::group_by(group) |>
    dplyr::summarise(dev_mean = mean(dev), myc_mean = mean(myc), .groups = "drop") |>
    print()

  # Felsher integrative signature present (AP4 consumer)?
  "MYC_felsher_integrative_signature" %in% rownames(sc)
  # Biogenesis-discrimination sets present (Fig 2 / AP7 consumers)?
  grep("MITO$|BICLUSTER", rownames(sc), value = TRUE) |> head(20)
}
