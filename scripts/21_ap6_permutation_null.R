# scripts/21_ap6_permutation_null.R
# =============================================================================
# AP6.2 -- expression x dispersion-matched permutation null (Block A, Day 4)
# =============================================================================
#
# fGSEA (AP6.1, script 20) showed MitoCarta strongly enriched by Myc. But mito
# genes are highly expressed and low-dispersion, and such genes carry larger,
# more-reliable effect statistics regardless of biology -- so "enriched" could be
# a housekeeping/expression artifact. This script tests PREFERENTIALITY against a
# null that matches each compartment's baseMean x dispersion profile: if the
# mitochondrial |LFC| exceeds what expression-matched random genes give, the
# alteration is preferential, not just a consequence of high expression.
#
# Metric: mean |raw LFC| (raw per the shrinkage rule -- a magnitude test),
#   primary on the Myc-at-6W contrast, myc_12W as robustness.
# Null: bin all genes by baseMean-decile x dispersion-decile (100 bins); for each
#   member gene draw a random gene from its OWN bin; B matched sets -> null of
#   mean |LFC|; locate the observed value (z, empirical p).
#
# Honest scope: this defeats the EXPRESSION/housekeeping confound (AP6's concern).
#   It does not fully escape inter-gene correlation (genes sampled independently),
#   but for a MAGNITUDE statistic that inflation is milder, and a large observed-
#   vs-matched gap is the point.
#
# Input:  results/interaction_results.rds (myc_6W_raw/myc_12W_raw: baseMean, raw LFC),
#         results/dds_int_run.rds (dispersions), results/ortholog_table.rds,
#         by_category/{01_mitocarta,04_metabolism,05_proliferation,03_mammary_development}
#         *_mouse.gmt, msigdbr Hallmark (MYC_TARGETS_V1, OXPHOS, E2F_TARGETS)
# Output: results/ap6_permutation_null.rds; outputs/ap6_null/preferential_alteration.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# =============================================================================
# PART 1: LOAD
# =============================================================================

ir       <- readRDS(here::here("results", "interaction_results.rds"))
dds_int  <- readRDS(here::here("results", "dds_int_run.rds"))
ortholog <- readRDS(here::here("results", "ortholog_table.rds"))

gmt <- function(f) fgsea::gmtPathways(
  here::here("data", "genesets_from_library", "by_category", f))

# =============================================================================
# PART 2: GENE TABLE (baseMean, dispersion, |LFC|, symbol) + binning
# =============================================================================

r6  <- ir$myc_6W_raw
r12 <- ir$myc_12W_raw
disp <- DESeq2::dispersions(dds_int)
names(disp) <- rownames(dds_int)

ens2sym <- ortholog |>
  dplyr::select(ensembl_gene_id, external_gene_name) |>
  dplyr::filter(!is.na(external_gene_name), external_gene_name != "") |>
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

gene_tbl <- tibble::tibble(
  ensembl    = rownames(r6),
  baseMean   = r6$baseMean,
  lfc        = r6$log2FoldChange,
  lfc_12W    = r12$log2FoldChange[match(rownames(r6), rownames(r12))],
  dispersion = disp[rownames(r6)]
) |>
  dplyr::filter(baseMean > 0, !is.na(dispersion), dispersion > 0,
                !is.na(lfc), !is.na(lfc_12W)) |>
  dplyr::mutate(absLFC = abs(lfc), absLFC_12W = abs(lfc_12W)) |>
  dplyr::left_join(ens2sym, by = c("ensembl" = "ensembl_gene_id")) |>
  dplyr::rename(symbol = external_gene_name)

# 10 x 10 expression x dispersion bins
gene_tbl <- gene_tbl |>
  dplyr::mutate(
    bm_bin   = dplyr::ntile(log10(baseMean), 10),
    disp_bin = dplyr::ntile(log10(dispersion), 10),
    bin      = (bm_bin - 1L) * 10L + disp_bin
  )
message(sprintf("Gene table: %d genes, %d with symbol, %d bins occupied",
                nrow(gene_tbl), sum(!is.na(gene_tbl$symbol)),
                dplyr::n_distinct(gene_tbl$bin)))

# =============================================================================
# PART 3: COMPARTMENTS (query = MitoCarta; comparators)
# =============================================================================

mg_union <- function(f) {
  p <- gmt(f)
  unique(unlist(p, use.names = FALSE))
}
mammary_syms <- {
  p <- gmt("03_mammary_development_mouse.gmt")
  unique(unlist(p[grepl("^MG_", names(p))], use.names = FALSE))  # dev sets, not METABRIC
}

hh <- tryCatch(msigdbr::msigdbr(species = "Mus musculus", category = "H"),
               error = function(e) msigdbr::msigdbr(species = "Mus musculus",
                                                    collection = "H"))
sym_col  <- intersect(c("gene_symbol", "db_gene_symbol", "mouse_symbol"), colnames(hh))[1]
name_col <- intersect(c("gs_name"), colnames(hh))[1]
hh_set <- function(nm) hh[[sym_col]][hh[[name_col]] == nm]

compartments <- list(
  MitoCarta               = mg_union("01_mitocarta_mouse.gmt"),
  Metabolism              = mg_union("04_metabolism_mouse.gmt"),
  Proliferation           = mg_union("05_proliferation_mouse.gmt"),
  Mammary_development     = mammary_syms,
  HALLMARK_MYC_TARGETS_V1 = hh_set("HALLMARK_MYC_TARGETS_V1"),
  HALLMARK_OXPHOS         = hh_set("HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
  HALLMARK_E2F_TARGETS    = hh_set("HALLMARK_E2F_TARGETS")
)
stopifnot(all(vapply(compartments, length, integer(1)) > 0))

# =============================================================================
# PART 4: MATCHED PERMUTATION NULL
# =============================================================================

B <- 5000L
run_null <- function(member_sym, metric_vec, bin_vals) {
  qi <- which(gene_tbl$symbol %in% member_sym & !is.na(gene_tbl$symbol))
  if (length(qi) < 5) return(NULL)
  obs    <- mean(metric_vec[qi])
  bins_q <- gene_tbl$bin[qi]
  samp   <- vapply(bins_q,
                   function(b) sample(bin_vals[[as.character(b)]], B, replace = TRUE),
                   numeric(B))                      # B x length(qi)
  null_means <- rowMeans(samp)
  nm <- mean(null_means); ns <- stats::sd(null_means)
  list(n = length(qi), obs = obs, null_mean = nm, null_sd = ns,
       z = (obs - nm) / ns,
       p_emp = (1 + sum(null_means >= obs)) / (B + 1),
       null = null_means)
}

set.seed(1)
bin_vals_6W  <- split(gene_tbl$absLFC,     gene_tbl$bin)
bin_vals_12W <- split(gene_tbl$absLFC_12W, gene_tbl$bin)

res_6W  <- lapply(compartments, run_null, metric_vec = gene_tbl$absLFC,     bin_vals = bin_vals_6W)
res_12W <- lapply(compartments, run_null, metric_vec = gene_tbl$absLFC_12W, bin_vals = bin_vals_12W)

to_row <- function(res_list, contrast) {
  do.call(rbind, lapply(names(res_list), function(nm) {
    r <- res_list[[nm]]
    if (is.null(r)) return(NULL)
    data.frame(compartment = nm, contrast = contrast, n = r$n, obs = r$obs,
               null_mean = r$null_mean, null_sd = r$null_sd, z = r$z, p_emp = r$p_emp)
  }))
}
null_table <- tibble::as_tibble(rbind(to_row(res_6W, "myc_6W"),
                                      to_row(res_12W, "myc_12W")))
message("AP6.2 matched-null (myc_6W):")
print(null_table |> dplyr::filter(contrast == "myc_6W") |> dplyr::arrange(dplyr::desc(z)))

# =============================================================================
# PART 5: FIGURE (observed vs matched null, myc_6W)
# =============================================================================

out_dir <- here::here("outputs", "ap6_null")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

plt <- null_table |> dplyr::filter(contrast == "myc_6W") |>
  dplyr::mutate(compartment = stats::reorder(compartment, z),
                preferential = p_emp < 0.05)
p_ap6 <- ggplot2::ggplot(plt, ggplot2::aes(y = compartment)) +
  ggplot2::geom_segment(ggplot2::aes(x = null_mean - 2 * null_sd,
                                     xend = null_mean + 2 * null_sd,
                                     yend = compartment), colour = "grey70", linewidth = 3) +
  ggplot2::geom_point(ggplot2::aes(x = null_mean), colour = "grey40", size = 2) +
  ggplot2::geom_point(ggplot2::aes(x = obs, colour = preferential), size = 3.5) +
  ggplot2::scale_colour_manual(values = c(`TRUE` = "#D73027", `FALSE` = "grey50"),
                               name = "p_emp < 0.05") +
  ggplot2::labs(
    title = "AP6.2: preferential alteration vs expression x dispersion-matched null",
    subtitle = "grey band = matched null mean +/- 2SD; red point = observed mean |LFC| (myc_6W)",
    x = "mean |raw LFC|", y = NULL) +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(out_dir, "preferential_alteration.pdf"), p_ap6,
                width = 8, height = 5)

# =============================================================================
# PART 6: SAVE (+ ROAST/CAMERA cut-line stub)
# =============================================================================
# ROAST/CAMERA competitive/self-contained tests are a plan cut-line (section 9);
# the correlation-aware CAMERA read already exists (script 17). Guarded stub only:
if (FALSE) {
  # e.g. limma::camera(vst_expr, idx_mito, design, contrast = myc_6W_col) -- not run.
}

ap6_out <- list(
  null_table = null_table,
  results_6W = res_6W,
  results_12W = res_12W,
  params = list(B = B, n_genes = nrow(gene_tbl), bins = "10 baseMean x 10 dispersion",
                metric = "mean |raw LFC|", primary = "myc_6W"),
  notes = paste(
    "AP6.2 expression x dispersion-matched permutation null. z = (obs - null_mean)/",
    "null_sd; p_emp one-sided (obs >= null). Preferential = obs beyond the matched",
    "null. Defeats the expression/housekeeping confound (AP6 concern); shares the",
    "inter-gene-correlation caveat but milder for a magnitude statistic. Compare",
    "MitoCarta z to generic comparators (metabolism/Hallmark); proliferation may be",
    "co-equal (Myc drives it too).")
)
saveRDS(ap6_out, here::here("results", "ap6_permutation_null.rds"))
message("Saved results/ap6_permutation_null.rds")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  ap6 <- readRDS(here::here("results", "ap6_permutation_null.rds"))

  # First-line: dispersions retrievable + gene-table coverage
  length(DESeq2::dispersions(dds_int)); summary(gene_tbl$dispersion)
  # Bin occupancy (want no near-empty bins)
  gene_tbl |> dplyr::count(bin) |> dplyr::arrange(n) |> head(10) |> print()

  # THE RESULT: is MitoCarta preferential, and beyond the comparators?
  ap6$null_table |> dplyr::filter(contrast == "myc_6W") |>
    dplyr::arrange(dplyr::desc(z)) |> print()
  # Robustness at 12W
  ap6$null_table |> dplyr::filter(contrast == "myc_12W") |>
    dplyr::arrange(dplyr::desc(z)) |> print()

  # Where does MitoCarta's observed sit in its own null?
  m <- ap6$results_6W$MitoCarta
  cat(sprintf("MitoCarta obs %.3f vs null %.3f +/- %.3f (z=%.1f, p=%.4g)\n",
              m$obs, m$null_mean, m$null_sd, m$z, m$p_emp))

  list.files(here::here("outputs", "ap6_null"), pattern = "\\.pdf$")
}
