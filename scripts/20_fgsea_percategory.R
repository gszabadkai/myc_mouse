# scripts/20_fgsea_percategory.R
# =============================================================================
# Per-category directional fGSEA (Block A, Day 3): AP6.1 + retention/abund inputs
# =============================================================================
#
# The decomposed, directional enrichment layer. Unlike the pooled run in script
# 10 (kept as exploratory-omnibus only), this runs ONE fGSEA per (ranking x
# category), so BH is WITHIN each category -- the scope publication claims cite,
# and what makes the "preferential, not just absolute" mito-focus argument (AP6)
# legible.
#
# fGSEA reads ABSOLUTE transcriptional enrichment per contrast and ranks on the
# unshrunken Wald `stat` (shrinkage-independent) -- a different lens from the
# underpowered interaction tests. The Myc genotype rankings (myc_6W, myc_12W) are
# the informative, well-powered lens (cf. Gate 1); the interaction ranking is
# expected weak (script 10 found no interaction pathways) and is reported as such.
#
# Rankings: 5 unshrunken contrasts from interaction_results.rds, stat-ranked;
#   ENSMUSG -> mouse symbol via inner_join on the cached ortholog table (per
#   CLAUDE.md), duplicate symbols collapsed by max |stat| (matching script 10).
# Pathways: 9 library category GMTs filtered to fgsea_eligible sets (drops the
#   gsva-only large sets), plus fresh MSigDB Hallmark (mouse) as a comparator.
#
# Input:  results/interaction_results.rds, results/ortholog_table.rds,
#         data/genesets_from_library/by_category/*_mouse.gmt,
#         data/genesets_from_library/provenance_table.csv, msigdbr (Hallmark)
# Output: results/fgsea_percategory.rds (tidy: ranking, category, pathway, NES,
#         pval, padj_within_category, size, leadingEdge) + a summary + AP6 preview
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# =============================================================================
# PART 1: LOAD
# =============================================================================

ir       <- readRDS(here::here("results", "interaction_results.rds"))
ortholog <- readRDS(here::here("results", "ortholog_table.rds"))
prov     <- readr::read_csv(
  here::here("data", "genesets_from_library", "provenance_table.csv"),
  show_col_types = FALSE)
fgsea_ok <- prov$set_name[prov$fgsea_eligible]

gmt_dir   <- here::here("data", "genesets_from_library", "by_category")
gmt_files <- list.files(gmt_dir, pattern = "_mouse\\.gmt$", full.names = TRUE)
stopifnot(length(gmt_files) == 9L)

# =============================================================================
# PART 2: RANKINGS (5 unshrunken contrasts, Wald stat, ENSMUSG -> symbol)
# =============================================================================

ens2sym <- ortholog |>
  dplyr::select(ensembl_gene_id, external_gene_name) |>
  dplyr::filter(!is.na(external_gene_name), external_gene_name != "") |>
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

make_ranks <- function(res) {
  df <- tibble::as_tibble(as.data.frame(res), rownames = "ensembl_gene_id") |>
    dplyr::filter(!is.na(stat)) |>
    dplyr::inner_join(ens2sym, by = "ensembl_gene_id") |>
    dplyr::group_by(external_gene_name) |>
    dplyr::slice_max(abs(stat), n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(stat))
  stats::setNames(df$stat, df$external_gene_name)
}

rankings <- list(
  myc_6W        = ir$myc_6W_raw,
  myc_12W       = ir$myc_12W_raw,
  timepoint_neg = ir$timepoint_neg_raw,
  timepoint_pos = ir$timepoint_pos_raw,
  interaction   = ir$interaction_raw
)
ranks_list <- lapply(rankings, make_ranks)
message(sprintf("Rankings built: %s",
  paste(sprintf("%s(%d)", names(ranks_list), vapply(ranks_list, length, integer(1))),
        collapse = ", ")))

# =============================================================================
# PART 3: PATHWAYS (fgsea-eligible library categories + fresh Hallmark)
# =============================================================================

cat_pathways <- lapply(gmt_files, function(f) {
  p <- fgsea::gmtPathways(f)
  p[names(p) %in% fgsea_ok]                         # fgsea-eligible only
})
names(cat_pathways) <- sub("_mouse\\.gmt$", "", basename(gmt_files))

# Fresh MSigDB Hallmark (mouse) -- robust to the msigdbr category/collection rename
hh <- tryCatch(msigdbr::msigdbr(species = "Mus musculus", category = "H"),
               error = function(e) msigdbr::msigdbr(species = "Mus musculus",
                                                    collection = "H"))
sym_col  <- intersect(c("gene_symbol", "db_gene_symbol", "mouse_symbol"), colnames(hh))[1]
name_col <- intersect(c("gs_name"), colnames(hh))[1]
stopifnot(!is.na(sym_col), !is.na(name_col))
cat_pathways[["hallmark_msigdb"]] <- split(hh[[sym_col]], hh[[name_col]])

message(sprintf("Categories: %s",
  paste(sprintf("%s(%d)", names(cat_pathways),
                vapply(cat_pathways, length, integer(1))), collapse = ", ")))

# =============================================================================
# PART 4: fGSEA PER (ranking x category), BH WITHIN EACH CALL
# =============================================================================

min_sz <- 10L; max_sz <- 500L

run_one <- function(ranks, pathways, ranking_name, category_name) {
  if (length(pathways) == 0) return(NULL)
  res <- tryCatch(
    suppressWarnings(fgsea::fgsea(pathways = pathways, stats = ranks,
                                  minSize = min_sz, maxSize = max_sz, eps = 0)),
    error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) return(NULL)
  tibble::as_tibble(res) |>
    dplyr::transmute(ranking = ranking_name, category = category_name,
                     pathway, NES, pval, padj_within_category = padj,
                     size, leadingEdge)
}

fgsea_all <- purrr::map_dfr(names(ranks_list), function(rn)
  purrr::map_dfr(names(cat_pathways), function(cn)
    run_one(ranks_list[[rn]], cat_pathways[[cn]], rn, cn)))

message(sprintf("fGSEA rows: %d across %d ranking x category runs",
                nrow(fgsea_all),
                length(ranks_list) * length(cat_pathways)))

# =============================================================================
# PART 5: SUMMARY + AP6 PREVIEW
# =============================================================================

fgsea_summary <- fgsea_all |>
  dplyr::group_by(ranking, category) |>
  dplyr::summarise(n_sets     = dplyr::n(),
                   n_sig      = sum(padj_within_category < 0.05, na.rm = TRUE),
                   median_NES = stats::median(NES, na.rm = TRUE),
                   .groups    = "drop")

# AP6 preview: on the myc_6W ranking, is MitoCarta preferentially enriched vs the
# metabolism / proliferation / Hallmark comparators? (descriptive; the hard test
# is the expression x dispersion permutation null in script 21).
out_dir <- here::here("outputs", "fgsea_percategory")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

p_ap6 <- fgsea_all |>
  dplyr::filter(ranking == "myc_6W") |>
  ggplot2::ggplot(ggplot2::aes(x = stats::reorder(category, NES, FUN = stats::median),
                               y = NES)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
  ggplot2::geom_boxplot(outlier.size = 0.6, fill = "grey85") +
  ggplot2::coord_flip() +
  ggplot2::labs(title = "AP6 preview: per-category NES on the Myc-at-6W contrast",
    subtitle = "preferential mito focus = MitoCarta shifted vs metabolism/proliferation/Hallmark",
    x = NULL, y = "fGSEA NES (myc_6W)") +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "ap6_nes_by_category_myc6W.pdf"), p_ap6,
                width = 7.5, height = 5)

# =============================================================================
# PART 6: SAVE
# =============================================================================

fgsea_out <- list(
  fgsea      = fgsea_all,
  summary    = fgsea_summary,
  params     = list(minSize = min_sz, maxSize = max_sz, eps = 0,
                    rankings = names(ranks_list),
                    n_genes_ranked = vapply(ranks_list, length, integer(1))),
  notes = paste(
    "One fGSEA per (ranking x category); padj_within_category = BH within that",
    "run (the publication scope, vs script 10's pooled). Ranks = unshrunken Wald",
    "stat (shrinkage-independent). fgsea-eligible library sets + fresh Hallmark.",
    "Myc genotype rankings (myc_6W/12W) are the powered lens; interaction ranking",
    "expected weak. AP6 'preferential' hard test = permutation null in script 21.")
)
saveRDS(fgsea_out, here::here("results", "fgsea_percategory.rds"))
message("Saved results/fgsea_percategory.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  fo <- readRDS(here::here("results", "fgsea_percategory.rds"))

  # Rankings sane? stat ranges; interaction should be ~N(0,1) / weak (script 10).
  lapply(ranks_list, function(r) round(range(r), 2))
  stats::ks.test(ranks_list$interaction, "pnorm", 0, 1)   # expect close to normal
  stats::ks.test(ranks_list$myc_6W, "pnorm", 0, 1)        # expect strong departure

  # Per-category coverage + how many sig per run
  fo$summary |> dplyr::arrange(ranking, dplyr::desc(n_sig)) |> print(n = Inf)

  # AP6: MitoCarta vs comparators on myc_6W (median NES + n_sig)
  fo$summary |> dplyr::filter(ranking == "myc_6W",
      category %in% c("01_mitocarta", "04_metabolism", "05_proliferation",
                      "hallmark_msigdb")) |> print()

  # Top MitoCarta pathways enriched by Myc at 6W
  fo$fgsea |> dplyr::filter(ranking == "myc_6W", category == "01_mitocarta") |>
    dplyr::arrange(padj_within_category) |>
    dplyr::select(pathway, NES, padj_within_category, size) |> head(15) |> print()

  # Interaction ranking: confirm it is weak (few/no sig), unlike the genotype ones
  fo$fgsea |> dplyr::filter(ranking == "interaction",
      padj_within_category < 0.05) |> nrow()

  list.files(here::here("outputs", "fgsea_percategory"), pattern = "\\.pdf$")
}
