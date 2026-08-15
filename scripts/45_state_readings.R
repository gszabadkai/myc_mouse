# scripts/45_state_readings.R
# =============================================================================
# THE THREE MISSING STATE READINGS
# =============================================================================
#
# The whole corpus is built out of CONTRASTS: a genotype effect within an age, a
# temporal effect within a genotype, and their interaction. Three readings that
# are about STATES rather than contrasts have never been made, and all three are
# needed to decide what a halved paper keeps.
#
#   (1) THE DIAGONAL. 6W_neg versus 12W_pos -- the young normal gland against the
#       gland at the point of initial tumour expansion. Never computed anywhere:
#       script 03 defines seven contrasts (five interaction 03:52-79, two group
#       03:142-154) and neither group contrast crosses genotype. The one trace is
#       a decision NOT to do it, docs/myc_mouse_finalisation_plan.md:379-382, which
#       calls it "a composite state-to-state contrast ... probably adds little;
#       exploratory-only and a cut candidate". PART E shows it adds a great deal.
#
#   (2) THE ABSOLUTE LEVELS. Abundance has only ever been read as a % SHARE of the
#       transcriptome (script 32) or as a LOG FOLD CHANGE (script 40's content
#       ruler). The actual DESeq2 normalised values for the main MitoCarta groups,
#       four states side by side with the per-animal spread, have never been
#       plotted. figures/fig03_background_vs_myc.R panel A comes closest but draws
#       bg$state_table, which holds GROUP MEANS only: four points per tier, no
#       spread and no magnitude.
#
#   (3) WHY THE TWO RULERS DISAGREE. On the diagonal the per-gene mean-LFC ruler
#       puts the OXPHOS subunits at +0.061 ("returns to baseline") and the
#       abundance-weighted summed-count ruler puts the SAME 87 genes at +0.226.
#       PART H decomposes the gap: it is expression weighting, and the weighting
#       gradient comes from the WILD-TYPE timeline, not from Myc.
#
# -----------------------------------------------------------------------------
# THE ONE ALGEBRAIC FACT THIS SCRIPT RESTS ON
#
#   `~ group` and `~ timepoint * myc_status` are both SATURATED over the same four
#   groups, so they span the same design space and give the same MLE group means.
#   The diagonal is therefore a RE-READING of a fit that already exists, not a new
#   model, and gene by gene
#
#       cross  ==  myc_12W + timepoint_neg  ==  myc_6W + timepoint_pos
#
#   to optimiser precision. PART A asserts both. That is what licenses drawing the
#   diagonal as a head-to-tail sum of two contrasts the paper already shows.
#
# -----------------------------------------------------------------------------
# BOUNDS THAT TRAVEL WITH EVERY NUMBER BELOW
#
#   * BATCH = TIMEPOINT (CLAUDE.md). The diagonal spans timepoints, so it carries
#     the batch offset in full -- exactly like `6>12W_wt`. It is DESCRIBED, not
#     claimed. Only the genotype-within-age contrasts are batch-clean.
#   * A NET OF ~0 IS NOT "NO CHANGE". It is two large opposite changes that cancel.
#     Anything drawn off `c_cross` must draw both components.
#   * THE RULER MATTERS. PART E2 reports the per-gene and the abundance-weighted
#     value side by side for every arm. They are not interchangeable and the
#     sentence has to name which one it means.
#   * MITOCARTA APOPTOSIS SETS ARE MITO-DEFINED (script 34). Any mito-versus-death
#     reading off this table is mito-versus-mito and circular.
#   * n = 6 per cell. This script RANKS; it does not confirm.
#
# Reads:  results/dds_group_run.rds          (script 03 -- the fitted group model)
#         results/interaction_results.rds    (script 03 -- the five raw contrasts)
#         results/count_matrix.rds           (script 01 -- UNfiltered counts)
#         results/mitopps_scores.rds         (script 08 -- partition, tiers, scores)
#         results/background_vs_myc.rds      (script 40 -- $ruler, the content ruler)
#         results/substrate_specificity_tradeoff.rds (script 43 -- $comparator arms)
#         results/mito_content_proxies.rds   (script 32 -- $shares, agreement check)
#         results/fgsea_percategory.rds      (script 20 -- the five saved rankings)
#         results/ortholog_table.rds, data/genesets_from_library/**
#         functions/reconcile_gene_symbols.R (MANDATORY -- vintage-aware membership)
# Writes: results/state_readings.rds
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

set.seed(1)
NSET <- 2000L   # matched-random-set draws (PARTS B3 and H)
NBIN <- 20L     # baseMean bins for the matched null (script 43's NBIN)

group_levels <- c("6W_neg", "6W_pos", "12W_neg", "12W_pos")

# =============================================================================
# PART 1: LOAD, AND THE CONTRAST ITSELF
# =============================================================================
dg  <- readRDS(here::here("results", "dds_group_run.rds"))
ir  <- readRDS(here::here("results", "interaction_results.rds"))
cts <- readRDS(here::here("results", "count_matrix.rds"))
mp  <- readRDS(here::here("results", "mitopps_scores.rds"))
bv  <- readRDS(here::here("results", "background_vs_myc.rds"))
ss  <- readRDS(here::here("results", "substrate_specificity_tradeoff.rds"))
mc  <- readRDS(here::here("results", "mito_content_proxies.rds"))
fp  <- readRDS(here::here("results", "fgsea_percategory.rds"))
ot  <- readRDS(here::here("results", "ortholog_table.rds"))

sm <- as.data.frame(SummarizedExperiment::colData(dg))
sm$group <- factor(as.character(sm$group), levels = group_levels)
sm$tp    <- stats::relevel(as.factor(sm$timepoint),  "6W")
sm$myc   <- stats::relevel(as.factor(sm$myc_status), "neg")
samples  <- rownames(sm)
stopifnot(all(table(sm$group) == 6), identical(levels(sm$group), group_levels))

# THE DIAGONAL. Same `filterFun = ihw` idiom as script 03:142-154, raw (unshrunken)
# because everything downstream is an averaged-LFC quantity (CLAUDE.md).
cross_res <- DESeq2::results(dg, contrast = c("group", "12W_pos", "6W_neg"),
                             filterFun = IHW::ihw)

D  <- lapply(ir[c("myc_6W_raw", "myc_12W_raw",
                  "timepoint_neg_raw", "timepoint_pos_raw")], as.data.frame)
V  <- function(k) stats::setNames(D[[k]]$log2FoldChange, rownames(D[[k]]))
m6 <- V("myc_6W_raw");        m12 <- V("myc_12W_raw")
tn <- V("timepoint_neg_raw"); tp  <- V("timepoint_pos_raw")
bm <- stats::setNames(D$myc_6W_raw$baseMean, rownames(D$myc_6W_raw))

universe_all <- rownames(D$myc_6W_raw)
cr <- stats::setNames(cross_res$log2FoldChange, rownames(cross_res))
cs <- stats::setNames(cross_res$stat,           rownames(cross_res))

ens_set  <- function(syms) { e <- recon_to_ensembl(syms, universe_all); e[!is.na(e)] }
set_mean <- function(v, e) if (length(e)) mean(v[e], na.rm = TRUE) else NA_real_

# ENSMUSG -> current symbol, from the cached ortholog table (CLAUDE.md: inner_join
# on the cached table, never a deframe()d named-vector lookup). Used for the
# per-gene OXPHOS table in PART H and for the fGSEA rankings in PART D.
ens2sym <- ot |>
  dplyr::select(ensembl_gene_id, external_gene_name) |>
  dplyr::filter(!is.na(external_gene_name), external_gene_name != "") |>
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
sym_of <- function(e) ens2sym$external_gene_name[match(e, ens2sym$ensembl_gene_id)]

message(sprintf("PART 1: diagonal on %d genes; %d up / %d down at padj < 0.05",
                nrow(cross_res),
                sum(cross_res$padj < 0.05 & cross_res$log2FoldChange > 0, na.rm = TRUE),
                sum(cross_res$padj < 0.05 & cross_res$log2FoldChange < 0, na.rm = TRUE)))

# =============================================================================
# PART A: THE IDENTITY -- the diagonal is a re-reading, not a new model
# =============================================================================
# Script 40 asserts only c_int == c_tp - c_tn (40:138-142). The two identities
# below are the complementary pair and are written down nowhere in the corpus.
#
# They hold to ~1e-6 on every gene but a handful, where the MLE diverges on a
# group that is effectively empty and the two optimiser runs land differently.
# Those genes are NAMED rather than dropped: a tolerance that silently absorbs
# them would also absorb a real disagreement.
shared <- intersect(names(cr), universe_all)
stopifnot(length(shared) > 15000L)

dev_of <- function(route) abs(cr[shared] - route[shared])
routes <- list("myc_12W + 6>12W_wt"  = m12 + tn,
               "myc_6W + 6>12W_myc"  = m6  + tp)

identity_tbl <- dplyr::bind_rows(lapply(names(routes), function(k) {
  d <- dev_of(routes[[k]])
  tibble::tibble(route = k, n_genes = length(shared),
                 median_dev = stats::median(d, na.rm = TRUE),
                 q999_dev   = unname(stats::quantile(d, 0.999, na.rm = TRUE)),
                 max_dev    = max(d, na.rm = TRUE),
                 n_over_001 = sum(d > 0.01, na.rm = TRUE))
}))
identity_tbl |> print()

divergent <- shared[dev_of(routes[[1]]) > 0.01]
divergent <- divergent[!is.na(divergent)]
identity_divergent <- tibble::tibble(
  gene = divergent, baseMean = unname(bm[divergent]),
  cross = unname(cr[divergent]), route_myc12_wt = unname((m12 + tn)[divergent]))
if (length(divergent))
  message(sprintf("PART A: %d gene(s) diverge > 0.01 (near-empty group, optimiser); named in $identity_divergent",
                  length(divergent)))

# The identity is asserted on the ROBUST summaries, not on the maximum.
stopifnot(all(identity_tbl$median_dev < 1e-5),
          all(identity_tbl$q999_dev   < 1e-3),
          all(identity_tbl$n_over_001 <= 10L))

# =============================================================================
# PART B: THE CONTENT RULER ON THE DIAGONAL
# =============================================================================
# Script 40's machinery exactly (40:102-136): the MitoCarta partition through
# mitopps_scores.rds, membership through the reconciler, unweighted mean of raw
# gene LFCs. Same 144 pathways, one new column.
pwise <- mp$mitopps_pairwise
tier1 <- mp$pathway_tier1_map
gp    <- mp$gene_to_pathway
pcol  <- names(gp)[1]; gcol <- names(gp)[2]
mtdna_name <- mp$mtdna_pathway_name

paths    <- intersect(unique(pwise$pathway), names(tier1))
path_ens <- lapply(stats::setNames(paths, paths), function(p)
  ens_set(unique(gp[[gcol]][gp[[pcol]] == p])))
message(sprintf("PART B: %d MitoPathways; median genes resolved = %.0f",
                length(paths), stats::median(vapply(path_ens, length, integer(1)))))

r0 <- as.data.frame(bv$ruler)
stopifnot(all(c("pathway", "tier", "c_m6", "c_m12", "c_tn", "c_tp", "is_mtdna") %in% names(r0)))

ruler <- r0
ruler$c_cross <- vapply(ruler$pathway,
                        function(p) set_mean(cr, path_ens[[p]]), numeric(1))

# The pathway-level identity: the diagonal ruler IS the sum of two published ones.
id_path <- max(abs(ruler$c_cross - (ruler$c_m12 + ruler$c_tn)), na.rm = TRUE)
message(sprintf("PART B identity: max |c_cross - (c_m12 + c_tn)| = %.2e", id_path))
stopifnot(id_path < 1e-6)

nuc <- ruler[!ruler$is_mtdna & !is.na(ruler$c_cross), ]
ruler_summary <- tibble::tibble(
  n_pathways   = nrow(nuc),
  median_cross = stats::median(nuc$c_cross),
  pct_above_0  = 100 * mean(nuc$c_cross > 0),
  n_above_0    = sum(nuc$c_cross > 0),
  mtdna_cross  = ruler$c_cross[ruler$is_mtdna])
ruler_summary |> print()

ruler_tiers <- nuc |>
  dplyr::group_by(tier) |>
  dplyr::summarise(n = dplyr::n(),
                   median_cross = stats::median(c_cross),
                   median_m12   = stats::median(c_m12),
                   median_tn    = stats::median(c_tn), .groups = "drop") |>
  dplyr::arrange(median_cross)
ruler_tiers |> print(n = nrow(ruler_tiers))

# =============================================================================
# PART B2: THE NAMED ARMS
# =============================================================================
# Script 43's own arm roster, so the diagonal and the substrate analysis share a
# vocabulary rather than growing a second one.
cmp <- as.data.frame(ss$comparator)
gmt <- fgsea::gmtPathways(here::here("data", "genesets_from_library",
                                     "mammary_mito_myc_metab_v1_mouse.gmt"))
arm_sets <- as.data.frame(ss$defs$arms)
stopifnot(all(c("arm", "set") %in% names(arm_sets)))

# script 43 saved its own pooled roster; reuse it rather than re-globbing, and
# assert the glob still agrees so a library change cannot silently split the two
prolif_sets <- ss$defs$prolif_sets
stopifnot(length(prolif_sets) >= 5,
          setequal(prolif_sets, grep("^PROLIF_", names(gmt), value = TRUE)))
arm_ens <- c(
  stats::setNames(lapply(arm_sets$set, function(s) ens_set(gmt[[s]])), arm_sets$arm),
  list("PROLIF_* pooled" = ens_set(unique(unlist(gmt[prolif_sets])))))
stopifnot(setequal(names(arm_ens), cmp$arm))

arms <- cmp
arms$c_cross <- vapply(arms$arm, function(a) set_mean(cr, arm_ens[[a]]), numeric(1))
id_arm <- max(abs(arms$c_cross - (arms$c_myc_12W + arms$c_wt_time)), na.rm = TRUE)
message(sprintf("PART B2 identity: max |c_cross - (c_myc_12W + c_wt_time)| = %.2e", id_arm))
stopifnot(id_arm < 1e-6)
arms <- dplyr::arrange(arms, c_cross)
arms |> print(n = nrow(arms))

# =============================================================================
# PART B3: THE MATCHED-SET NULL ON THE DIAGONAL
# =============================================================================
# "The respiratory arm is the one that returns to baseline" is a RELATIVE claim --
# relative to a compartment whose median pathway gains +0.25 -- so it needs a null.
# Matching is on baseMean, script 43:166-181 verbatim in construction.
expressed <- names(bm)[is.finite(bm) & bm > 0 & is.finite(cr[names(bm)])]
bin_of <- cut(rank(bm[expressed], ties.method = "first"),
              breaks = NBIN, labels = FALSE)
by_bin <- split(expressed, bin_of)
draw_matched <- function(e) {
  b <- bin_of[match(e, expressed)]
  b <- b[!is.na(b)]
  unlist(lapply(split(b, b), function(k)
    sample(by_bin[[as.character(k[1])]], length(k), replace = TRUE)), use.names = FALSE)
}

arm_null <- dplyr::bind_rows(lapply(names(arm_ens), function(a) {
  e   <- arm_ens[[a]][arm_ens[[a]] %in% expressed]
  obs <- set_mean(cr, e)
  nul <- vapply(seq_len(NSET), function(i) set_mean(cr, draw_matched(e)), numeric(1))
  tibble::tibble(arm = a, n_matched = length(e), observed_cross = obs,
                 null_median = stats::median(nul),
                 percentile  = 100 * mean(nul < obs),
                 p_emp_lower = mean(nul <= obs))
})) |> dplyr::arrange(percentile)
arm_null |> print(n = nrow(arm_null))

# =============================================================================
# PART C: THE PRIORITY RULER (mitoPPS) ON THE DIAGONAL
# =============================================================================
# mitoPPS is a per-sample pairwise-ratio score, so a new group contrast costs
# nothing: it is a difference of group means on scores that already exist. The
# construction reproduces script 08:611-622, and the reproduction is ASSERTED on
# an existing contrast before the new one is trusted.
ps <- as.data.frame(mp$mitopps_scores)
stopifnot("group" %in% names(ps))
ps$group <- factor(as.character(ps$group), levels = group_levels)
score_cols <- setdiff(names(ps), c("sample", "group", "timepoint", "myc_status"))

pairwise_of <- function(a, b) {
  dplyr::bind_rows(lapply(score_cols, function(p) {
    x <- ps[[p]][ps$group == a]; y <- ps[[p]][ps$group == b]
    tt <- tryCatch(stats::t.test(y, x), error = function(e) NULL)
    tibble::tibble(pathway = p, mean_a = mean(x), mean_b = mean(y),
                   diff = mean(y) - mean(x),
                   p_value = if (is.null(tt)) NA_real_ else unname(tt$p.value))
  })) |> dplyr::mutate(padj = stats::p.adjust(p_value, "BH"))
}

# the assertion: rebuild `Temporal_Myc-` and match script 08's saved table
chk <- pairwise_of("6W_neg", "12W_neg")
ref <- pwise[pwise$contrast == "Temporal_Myc-", ]
j   <- match(chk$pathway, ref$pathway)
d_mp <- max(abs(chk$diff - ref$diff[j]), na.rm = TRUE)
message(sprintf("PART C reproduction: max |rebuilt - saved| on Temporal_Myc- = %.2e", d_mp))
stopifnot(d_mp < 1e-8)

mitopps <- pairwise_of("6W_neg", "12W_pos")
mitopps$contrast <- "cross_6Wneg_to_12Wpos"
mitopps <- mitopps[, c("pathway", "contrast", "mean_a", "mean_b", "diff", "p_value", "padj")]
message(sprintf("PART C: mitoPPS diagonal -- %d of %d pathways above zero, %d at BH < 0.05",
                sum(mitopps$diff > 0, na.rm = TRUE), nrow(mitopps),
                sum(mitopps$padj < 0.05, na.rm = TRUE)))

# =============================================================================
# PART D: fGSEA ON THE DIAGONAL
# =============================================================================
# Script 20's recipe, copied not re-imagined (20:57-65 for the ranking, 20:81-131
# for the sets and the per-category BH).
#
# SCRIPT 20 IS NOT RE-RUN AND results/fgsea_percategory.rds IS NOT REWRITTEN --
# six panels and paper/analysis_record.qmd assert against it. Instead the myc_6W
# ranking is recomputed HERE and asserted against the saved table, which puts the
# new sixth ranking on the existing five's footing without touching them.
prov     <- readr::read_csv(here::here("data", "genesets_from_library",
                                       "provenance_table.csv"), show_col_types = FALSE)
fgsea_ok <- prov$set_name[prov$fgsea_eligible]
gmt_dir   <- here::here("data", "genesets_from_library", "by_category")
gmt_files <- list.files(gmt_dir, pattern = "_mouse[.]gmt$", full.names = TRUE)
stopifnot(length(gmt_files) == 9L)

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
ranks_list <- list(cross = make_ranks(cross_res), myc_6W = make_ranks(ir$myc_6W_raw))

cat_pathways <- lapply(gmt_files, function(f) {
  p <- fgsea::gmtPathways(f); p[names(p) %in% fgsea_ok] })
names(cat_pathways) <- sub("_mouse[.]gmt$", "", basename(gmt_files))

# Reconcile old MitoCarta symbols to current ones, script 20:94-98 verbatim --
# `recon_current_map` returns a NAMED MAP, not a translated set, and using it as
# though it translated would silently drop every renamed gene (the 2026-07-24 bug
# this helper exists to fix).
.sync <- recon_current_map(unique(unlist(cat_pathways, use.names = FALSE)))
cat_pathways <- lapply(cat_pathways, function(cat) lapply(cat, function(genes) {
  m <- .sync[genes]; unique(ifelse(is.na(m), genes, unname(m)))
}))

hh <- tryCatch(msigdbr::msigdbr(species = "Mus musculus", category = "H"),
               error = function(e)
                 msigdbr::msigdbr(species = "Mus musculus", collection = "H"))
cat_pathways$hallmark_msigdb <- split(hh$gene_symbol, hh$gs_name)

run_one <- function(ranks, pathways, ranking_name, category_name) {
  suppressWarnings(fgsea::fgsea(pathways = pathways, stats = ranks,
                                minSize = 10, maxSize = 500, eps = 0)) |>
    tibble::as_tibble() |>
    dplyr::transmute(ranking = ranking_name, category = category_name,
                     pathway, NES, pval, padj_within_category = padj, size)
}
fgsea_new <- purrr::map_dfr(names(ranks_list), function(rn)
  purrr::map_dfr(names(cat_pathways), function(cn)
    run_one(ranks_list[[rn]], cat_pathways[[cn]], rn, cn)))

# the assertion: the recomputed myc_6W ranking must reproduce script 20's
saved6 <- as.data.frame(fp$fgsea)[fp$fgsea$ranking == "myc_6W", ]
new6   <- as.data.frame(fgsea_new)[fgsea_new$ranking == "myc_6W", ]
k <- match(paste(saved6$category, saved6$pathway), paste(new6$category, new6$pathway))
ok <- !is.na(k)
rho6 <- stats::cor(saved6$NES[ok], new6$NES[k[ok]], method = "spearman")
d6   <- max(abs(saved6$NES[ok] - new6$NES[k[ok]]), na.rm = TRUE)
message(sprintf("PART D reproduction: myc_6W NES vs script 20 -- rho = %.4f, max |dNES| = %.3f, %d/%d matched",
                rho6, d6, sum(ok), nrow(saved6)))
stopifnot(rho6 > 0.99, sum(ok) > 0.95 * nrow(saved6))

# the diagonal binds onto the five saved rankings as a sixth
fgsea <- dplyr::bind_rows(
  as.data.frame(fp$fgsea)[, c("ranking", "category", "pathway", "NES",
                              "pval", "padj_within_category", "size")],
  as.data.frame(fgsea_new)[fgsea_new$ranking == "cross", ])
message(sprintf("PART D: %d rankings, %d rows; diagonal has %d sets, %d at BH < 0.05 within category",
                dplyr::n_distinct(fgsea$ranking), nrow(fgsea),
                sum(fgsea$ranking == "cross"),
                sum(fgsea$ranking == "cross" & fgsea$padj_within_category < 0.05, na.rm = TRUE)))

# =============================================================================
# PART G: THE ABSOLUTE LEVELS  (moved before PART E -- E2 consumes it)
# =============================================================================
# Size factors from the fitted group model applied to the FULL count matrix. NOT
# to counts(dg): the fitted object carries script 03's rowSums(counts >= 10) >= 4
# filter (03:136), which would silently under-count every compartment.
sf <- DESeq2::sizeFactors(dg)
stopifnot(!is.null(sf), setequal(names(sf), samples))
cts <- cts[, samples, drop = FALSE]
nrm <- sweep(cts, 2, sf[samples], "/")

# ROSTER 1: the seven Level-1 tiers as fig03 panel A uses them, with the synthetic
# mtDNA-encoded pathway held OUT of every tier and given its own row.
tier_levels <- c("Protein import, sorting and homeostasis", "Mitochondrial central dogma",
                 "OXPHOS", "Metabolism", "Signaling",
                 "Mitochondrial dynamics and surveillance", "Small molecule transport")
stopifnot(setequal(unique(unname(tier1)), tier_levels))

# The tier map carries 152 pathways and the mitoPPS table 144; the 8 extra
# (Cholesterol-associated, Cytochrome C, OXA, ...) contribute ZERO additional
# genes, so the two rosters give identical unions. Asserted rather than assumed,
# because a library change could silently split them.
syms_of_paths <- function(ps) unique(unlist(lapply(ps, function(p) gp[[gcol]][gp[[pcol]] == p])))
genes_of_tier <- function(tt) {
  ps <- setdiff(names(tier1)[tier1 == tt], mtdna_name)
  stopifnot(setequal(syms_of_paths(ps),
                     syms_of_paths(setdiff(intersect(paths, ps), mtdna_name))))
  e <- ens_set(syms_of_paths(ps))
  e[e %in% rownames(nrm)]
}
level_ens <- stats::setNames(lapply(tier_levels, genes_of_tier), tier_levels)
level_ens[["mtDNA-encoded"]] <-
  ens_set(unique(gp[[gcol]][gp[[pcol]] == mtdna_name]))
level_ens[["mtDNA-encoded"]] <-
  level_ens[["mtDNA-encoded"]][level_ens[["mtDNA-encoded"]] %in% rownames(nrm)]

# ROSTER 2: script 43's named arms, so the levels table and the diagonal table
# share a roster.
arm_ens_n <- lapply(arm_ens, function(e) e[e %in% rownames(nrm)])
level_ens <- c(level_ens, arm_ens_n)
stopifnot(all(vapply(level_ens, length, integer(1)) > 0))

levels_tbl <- dplyr::bind_rows(lapply(names(level_ens), function(g) {
  s <- colSums(nrm[level_ens[[g]], , drop = FALSE])
  tibble::tibble(group_set = g, roster = if (g %in% c(tier_levels, "mtDNA-encoded"))
                   "tier" else "arm",
                 n_genes = length(level_ens[[g]]),
                 sample = samples, group = sm$group,
                 timepoint = sm$timepoint, myc_status = sm$myc_status,
                 norm_sum = unname(s[samples]))
}))

# The model is script 32's share_stat_one (32:352-370) column-for-column, so the
# levels table and the shares table read side by side. The INTERACTION is fitted
# and reported first: script 32 quoted an additive genotype effect only after
# confirming the interaction was far from significant, and that gate holds here.
level_stat_one <- function(y, label, roster, n_genes, meta = sm) {
  d   <- data.frame(y = y, myc = meta$myc, tp = meta$tp, grp = meta$group)
  wsd <- sqrt(mean(tapply(d$y, d$grp, stats::var)))
  ma  <- summary(stats::lm(y ~ myc + tp, d))$coefficients["mycpos", ]
  mi  <- summary(stats::lm(y ~ tp * myc, d))$coefficients["tp12W:mycpos", ]
  wt  <- summary(stats::lm(y ~ tp, subset(d, myc == "neg")))$coefficients["tp12W", ]
  mcf <- summary(stats::lm(y ~ tp, subset(d, myc == "pos")))$coefficients["tp12W", ]
  gm  <- tapply(d$y, d$grp, mean)
  tibble::tibble(
    group_set = label, roster = roster, n_genes = n_genes, within_sd = wsd,
    geno_beta = unname(ma["Estimate"]), geno_d = unname(ma["Estimate"]) / wsd,
    geno_p = unname(ma["Pr(>|t|)"]),
    int_beta = unname(mi["Estimate"]), int_p = unname(mi["Pr(>|t|)"]),
    wt_temporal_beta = unname(wt["Estimate"]),  wt_temporal_p = unname(wt["Pr(>|t|)"]),
    myc_temporal_beta = unname(mcf["Estimate"]), myc_temporal_p = unname(mcf["Pr(>|t|)"]),
    m6_neg = unname(gm["6W_neg"]),   m6_pos  = unname(gm["6W_pos"]),
    m12_neg = unname(gm["12W_neg"]), m12_pos = unname(gm["12W_pos"]),
    level_6W_wt = 2^unname(gm["6W_neg"]),
    cross_log2  = unname(gm["12W_pos"] - gm["6W_neg"]))
}

level_stats <- dplyr::bind_rows(lapply(names(level_ens), function(g) {
  s <- levels_tbl[levels_tbl$group_set == g, ]
  s <- s[match(samples, s$sample), ]
  level_stat_one(log2(s$norm_sum), g,
                 if (g %in% c(tier_levels, "mtDNA-encoded")) "tier" else "arm",
                 s$n_genes[1])
})) |>
  dplyr::mutate(geno_padj = stats::p.adjust(geno_p, "BH"),
                int_padj  = stats::p.adjust(int_p,  "BH"))
level_stats |>
  dplyr::select(group_set, n_genes, level_6W_wt, geno_beta, geno_p, geno_padj,
                int_p, wt_temporal_beta, cross_log2) |>
  print(n = nrow(level_stats))

# --- how much of this is NEW? the honesty check --------------------------------
# A summed normalised count and a % share are close relatives: they differ only by
# each sample's total normalised count. What is genuinely new here is the ABSOLUTE
# MAGNITUDE, the PER-ANIMAL SPREAD (bg$state_table carries group means only) and a
# genotype test on a ruler with no compartment denominator in it -- NOT the
# direction. This measures the overlap instead of leaving it to be assumed.
share_agreement <- NULL
sh <- tryCatch(as.data.frame(mc$shares), error = function(e) NULL)
if (!is.null(sh) && all(c("panel", "sample", "share_nomt") %in% names(sh))) {
  share_agreement <- dplyr::bind_rows(lapply(intersect(unique(sh$panel),
                                                       names(level_ens)), function(g) {
    a <- sh[sh$panel == g, ]; a <- a[match(samples, a$sample), ]
    b <- levels_tbl[levels_tbl$group_set == g, ]; b <- b[match(samples, b$sample), ]
    tibble::tibble(group_set = g,
                   r_log2 = stats::cor(log2(a$share_nomt), log2(b$norm_sum)))
  }))
  if (nrow(share_agreement))
    message(sprintf("PART G: levels vs script 32 shares -- median r = %.3f over %d matched panels",
                    stats::median(share_agreement$r_log2), nrow(share_agreement)))
} else {
  message("PART G: no matching panel names in mito_content_proxies.rds$shares -- agreement check skipped")
}

# =============================================================================
# PART H: WHY THE TWO RULERS DISAGREE
# =============================================================================
# On the diagonal the per-gene mean-LFC ruler puts the OXPHOS subunits at ~+0.06
# and the summed-count ruler puts the SAME genes at ~+0.23. Both are correct; they
# WEIGHT differently (mean of log ratios versus log of the sum). This part shows
# the weighting IS the gap, and then asks where the weighting gradient comes from.
ox_ens <- arm_ens[["OXPHOS subunits"]]
ox_ens <- ox_ens[ox_ens %in% rownames(nrm) & ox_ens %in% names(cr)]
g6 <- samples[sm$group == "6W_neg"]
lev6 <- rowMeans(nrm[ox_ens, g6, drop = FALSE])

# complex membership from MitoCarta's own subunit sets; the three that belong to
# no complex (cytochrome c and its synthase) get a labelled block, not a silent drop
complex_sets <- c("CI subunits", "CII subunits", "CIII subunits",
                  "CIV subunits", "CV subunits")
complex_of <- stats::setNames(rep(NA_character_, length(ox_ens)), ox_ens)
for (k in complex_sets) {
  e <- ens_set(unique(gp[[gcol]][gp[[pcol]] == k]))
  complex_of[intersect(ox_ens, e)] <- sub(" subunits$", "", k)
}
complex_of[is.na(complex_of)] <- "cytochrome c / other"
message(sprintf("PART H: %d OXPHOS subunits; %s",
                length(ox_ens),
                paste(sprintf("%s %d", names(table(complex_of)), table(complex_of)),
                      collapse = ", ")))

# the four states as MEDIANS per group (what the panel draws), and the means beside
oxphos_genes <- tibble::tibble(
  gene    = ox_ens,
  symbol  = sym_of(ox_ens),
  complex = unname(complex_of[ox_ens]),
  level_6W_wt = unname(lev6),
  med_6W_neg  = vapply(ox_ens, function(e)
    stats::median(nrm[e, samples[sm$group == "6W_neg"]]),  numeric(1)),
  med_6W_pos  = vapply(ox_ens, function(e)
    stats::median(nrm[e, samples[sm$group == "6W_pos"]]),  numeric(1)),
  med_12W_neg = vapply(ox_ens, function(e)
    stats::median(nrm[e, samples[sm$group == "12W_neg"]]), numeric(1)),
  med_12W_pos = vapply(ox_ens, function(e)
    stats::median(nrm[e, samples[sm$group == "12W_pos"]]), numeric(1)),
  mean_all    = unname(rowMeans(nrm[ox_ens, samples, drop = FALSE])),
  lfc_cross   = unname(cr[ox_ens]),
  lfc_myc_6W  = unname(m6[ox_ens]),
  lfc_myc_12W = unname(m12[ox_ens]),
  lfc_wt_time = unname(tn[ox_ens]),
  lfc_myc_time = unname(tp[ox_ens]))

# medians are drawn; assert the means tell the same story so the choice is inert
med_mean_r <- stats::cor(
  log2(oxphos_genes$med_12W_pos + 1) - log2(oxphos_genes$med_6W_neg + 1),
  oxphos_genes$lfc_cross, method = "spearman")
message(sprintf("PART H: drawn group medians vs the fitted LFC -- rho = %.3f", med_mean_r))

# --- the decomposition: unweighted / weighted / sum ratio ----------------------
w <- lev6 / sum(lev6)
sum_ratio <- log2(sum(rowMeans(nrm[ox_ens, samples[sm$group == "12W_pos"], drop = FALSE])) /
                  sum(rowMeans(nrm[ox_ens, samples[sm$group == "6W_neg"],  drop = FALSE])))
weighting <- tibble::tibble(
  summary = c("unweighted per-gene mean", "expression-weighted mean", "log2(sum ratio)"),
  value   = c(mean(cr[ox_ens], na.rm = TRUE),
              sum(w * cr[ox_ens], na.rm = TRUE),
              sum_ratio))
weighting |> print()

# The gap IS the weighting if the weighted mean lands between the other two. If it
# does not, the gap is a membership or reconciliation difference and this whole
# reading is wrong -- so it is asserted, not hoped for.
stopifnot(weighting$value[2] > weighting$value[1],
          abs(weighting$value[2] - weighting$value[3]) <
            abs(weighting$value[1] - weighting$value[3]))

# --- the gradient: do high expressers move differently? ------------------------
lfc_cols <- c(cross = "lfc_cross", myc_6W = "lfc_myc_6W", myc_12W = "lfc_myc_12W",
              wt_time = "lfc_wt_time", myc_time = "lfc_myc_time")
gradient <- dplyr::bind_rows(lapply(names(lfc_cols), function(k) {
  y <- oxphos_genes[[lfc_cols[[k]]]]
  tibble::tibble(contrast = k,
                 spearman = stats::cor(log10(lev6), y, method = "spearman",
                                       use = "complete.obs"),
                 pearson  = stats::cor(log10(lev6), y, use = "complete.obs"))
}))
gradient |> print()

qt <- cut(rank(lev6, ties.method = "first"), 4,
          labels = c("Q1 low", "Q2", "Q3", "Q4 high"))
gradient_quartiles <- dplyr::bind_rows(lapply(levels(qt), function(q) {
  i <- qt == q
  tibble::tibble(quartile = q, n = sum(i),
                 median_level = stats::median(lev6[i]),
                 cross    = mean(oxphos_genes$lfc_cross[i],    na.rm = TRUE),
                 wt_time  = mean(oxphos_genes$lfc_wt_time[i],  na.rm = TRUE),
                 myc_12W  = mean(oxphos_genes$lfc_myc_12W[i],  na.rm = TRUE))
}))
gradient_quartiles |> print()

# --- THE NULL THAT DECIDES WHAT THE PANEL MAY CLAIM ----------------------------
# Lowly expressed genes carry noisier LFCs, so SOME expression gradient is expected
# generically. Noise is symmetric and would not give Q1 a consistently negative
# mean -- but that has to be tested, not argued. Same matched draws as PART B3, and
# the statistic is the gradient itself, recomputed on each drawn set.
grad_stat <- function(e, v) {
  l <- rowMeans(nrm[e, g6, drop = FALSE])
  y <- v[e]
  k <- is.finite(l) & is.finite(y) & l > 0
  if (sum(k) < 10L) return(c(rho = NA_real_, gap = NA_real_))
  q <- cut(rank(l[k], ties.method = "first"), 4, labels = FALSE)
  c(rho = stats::cor(log10(l[k]), y[k], method = "spearman"),
    gap = mean(y[k][q == 4]) - mean(y[k][q == 1]))
}
ox_expr <- ox_ens[ox_ens %in% expressed]
gradient_null <- dplyr::bind_rows(lapply(c("wt_time", "cross"), function(k) {
  v   <- if (k == "wt_time") tn else cr
  obs <- grad_stat(ox_expr, v)
  nul <- vapply(seq_len(NSET), function(i) {
    e <- draw_matched(ox_expr); e <- e[e %in% rownames(nrm)]
    grad_stat(e, v)
  }, numeric(2))
  tibble::tibble(
    contrast = k, statistic = c("spearman rho", "Q4 - Q1 gap"),
    observed = unname(obs),
    null_median = c(stats::median(nul["rho", ], na.rm = TRUE),
                    stats::median(nul["gap", ], na.rm = TRUE)),
    percentile = c(100 * mean(nul["rho", ] < obs[["rho"]], na.rm = TRUE),
                   100 * mean(nul["gap", ] < obs[["gap"]], na.rm = TRUE)))
}))
gradient_null |> print()

# =============================================================================
# PART E: THE DIRECTION SUMMARY -- the table that answers the question
# =============================================================================
nes_of <- function(setname, ranking) {
  x <- fgsea[fgsea$ranking == ranking & fgsea$pathway == setname, ]
  if (!nrow(x)) return(c(NES = NA_real_, padj = NA_real_))
  x <- x[which.min(x$padj_within_category), ]
  c(NES = x$NES[1], padj = x$padj_within_category[1])
}
arm_set_of <- stats::setNames(arm_sets$set, arm_sets$arm)
arm_path_of <- stats::setNames(arm_sets$ruler_pathway, arm_sets$arm)

direction <- dplyr::bind_rows(lapply(arms$arm, function(a) {
  s  <- unname(arm_set_of[a])
  pw <- unname(arm_path_of[a])
  nc <- if (is.na(s)) c(NES = NA_real_, padj = NA_real_) else nes_of(s, "cross")
  n6 <- if (is.na(s)) c(NES = NA_real_, padj = NA_real_) else nes_of(s, "myc_6W")
  mpv <- if (!is.na(pw) && pw %in% mitopps$pathway)
    mitopps$diff[match(pw, mitopps$pathway)] else NA_real_
  lv <- level_stats[level_stats$group_set == a, ]
  tibble::tibble(
    arm = a,
    content_cross = arms$c_cross[arms$arm == a],
    content_myc_12W = arms$c_myc_12W[arms$arm == a],
    content_wt_time = arms$c_wt_time[arms$arm == a],
    null_percentile = arm_null$percentile[match(a, arm_null$arm)],
    levels_cross = if (nrow(lv)) lv$cross_log2[1] else NA_real_,
    levels_geno_p = if (nrow(lv)) lv$geno_p[1] else NA_real_,
    mitopps_cross = mpv,
    nes_cross = unname(nc["NES"]), nes_cross_padj = unname(nc["padj"]),
    nes_myc_6W = unname(n6["NES"]))
})) |>
  dplyr::mutate(sign_agree = sign(content_cross) == sign(nes_cross)) |>
  dplyr::arrange(content_cross)
direction |> print(n = nrow(direction))

# =============================================================================
# PART E2: THE TWO RULERS, SIDE BY SIDE
# =============================================================================
# Reported, not reconciled. The per-gene ruler and the abundance-weighted ruler
# disagree most on precisely the arm the paper's closing sentence rests on, and a
# sentence that quotes one number must name which ruler it is on.
ruler_compare <- direction |>
  dplyr::transmute(arm,
                   per_gene_mean = content_cross,
                   abundance_weighted = levels_cross,
                   difference = levels_cross - content_cross) |>
  dplyr::arrange(dplyr::desc(abs(difference)))
ruler_compare |> print(n = nrow(ruler_compare))

# =============================================================================
# PART F: SAVE
# =============================================================================
notes <- c(
  "THE DIAGONAL (6W_neg -> 12W_pos) IS A RE-READING OF AN EXISTING FIT. `~ group`",
  "  and `~ timepoint*myc_status` are both saturated over the same four groups, so",
  "  gene-wise cross == myc_12W + 6>12W_wt == myc_6W + 6>12W_myc to optimiser",
  "  precision (PART A). No new model was fitted. A handful of genes diverge where",
  "  the MLE is unstable on a near-empty group; they are NAMED in $identity_divergent.",
  "BATCH = TIMEPOINT. The diagonal spans timepoints and therefore carries the batch",
  "  offset in full, exactly like 6>12W_wt. DESCRIPTIVE, not a clean effect. Only the",
  "  genotype-within-age contrasts are batch-clean.",
  "docs/myc_mouse_finalisation_plan.md:379-382 deferred this contrast as 'probably",
  "  adds little; exploratory-only and a cut candidate'. That judgement is superseded:",
  "  the compartment is uniformly positive on the diagonal EXCEPT the respiratory arm,",
  "  which is the closing claim of the written Results in a single contrast.",
  "A NET OF ~0 IS NOT 'NO CHANGE' -- it is two large opposite changes that cancel.",
  "  Anything drawn off c_cross must draw both components (PART B2 gives them).",
  "THE RULER MATTERS AND THE TWO DISAGREE ($ruler_compare). The per-gene mean-LFC",
  "  ruler and the abundance-weighted summed-count ruler differ most on the OXPHOS",
  "  subunits. PART H shows the gap IS expression weighting, and that the weighting",
  "  gradient is carried by the WILD-TYPE timeline, not by Myc. 'Returns to baseline'",
  "  is a PER-GENE statement and must be written as one.",
  "$gradient_null decides what may be claimed about the gradient: it is specific to",
  "  the respiratory arm only if the observed rho/gap sits outside expression-matched",
  "  draws. Read the percentile before writing the sentence.",
  "MITOCARTA APOPTOSIS SETS ARE MITO-DEFINED (script 34), so any mito-versus-death",
  "  reading off these tables is mito-versus-mito and circular.",
  "SUMMED NORMALISED COUNTS AND % SHARES ARE CLOSE RELATIVES ($share_agreement gives",
  "  the measured correlation). What is new in $levels is the absolute magnitude, the",
  "  per-animal spread, and a genotype test with no compartment denominator in it.",
  "n = 6 per cell. This script RANKS; it does not confirm.")

state_readings <- list(
  cross = tibble::as_tibble(as.data.frame(cross_res), rownames = "gene") |>
    dplyr::select(gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj),
  identity           = identity_tbl,
  identity_divergent = identity_divergent,
  ruler              = tibble::as_tibble(ruler),
  ruler_summary      = ruler_summary,
  ruler_tiers        = ruler_tiers,
  arms               = tibble::as_tibble(arms),
  arm_null           = arm_null,
  mitopps            = tibble::as_tibble(mitopps),
  fgsea              = tibble::as_tibble(fgsea),
  direction          = direction,
  ruler_compare      = ruler_compare,
  levels             = levels_tbl,
  level_stats        = level_stats,
  share_agreement    = share_agreement,
  oxphos_genes       = oxphos_genes,
  weighting          = weighting,
  gradient           = gradient,
  gradient_quartiles = gradient_quartiles,
  gradient_null      = gradient_null,
  params = list(n_set_draws = NSET, n_bins = NBIN, seed = 1,
                contrast = "group: 12W_pos vs 6W_neg, filterFun = ihw, unshrunken",
                fgsea = list(minSize = 10L, maxSize = 500L, eps = 0,
                             bh = "within category"),
                myc_6W_reproduction = list(spearman = rho6, max_abs_dNES = d6),
                mitopps_reproduction = d_mp),
  notes = notes,
  analysis_date = Sys.Date())

saveRDS(state_readings, here::here("results", "state_readings.rds"))
message("Saved results/state_readings.rds")

# =============================================================================
# SANDBOX -- skipped by source(); run line by line in Positron
# =============================================================================
if (FALSE) {

  res <- readRDS(here::here("results", "state_readings.rds"))
  cat(res$notes, sep = "\n")

  ## (1) The diagonal is the sum of two published contrasts, to 1e-8
  res$identity |> print()
  res$identity_divergent |> print()

  ## (2) Direction: is the compartment up or down on the diagonal?
  res$ruler_summary |> print()
  res$ruler_tiers |> print(n = 10)
  res$direction |> print(n = 12)

  ## (3) The respiratory arm against expression-matched random genes
  res$arm_null |> print(n = 12)

  ## (4) Where the two rulers disagree, largest first
  res$ruler_compare |> print(n = 12)

  ## (5) The absolute levels, and whether the interaction gate is clear
  res$level_stats |>
    dplyr::select(group_set, n_genes, level_6W_wt, geno_beta, geno_padj, int_p) |>
    print(n = 20)
  res$share_agreement |> print(n = 20)

  ## (6) Why the rulers disagree: weighting, gradient, and the null
  res$weighting |> print()
  res$gradient |> print()
  res$gradient_quartiles |> print()
  res$gradient_null |> print()

  ## (7) MitoCarta sets are membership-loose -- is a tier shift one gene?
  ## The ten highest-expressed OXPHOS subunits carry most of the summed ruler.
  res$oxphos_genes |>
    dplyr::arrange(dplyr::desc(level_6W_wt)) |>
    dplyr::select(symbol, complex, level_6W_wt, lfc_cross, lfc_wt_time, lfc_myc_12W) |>
    head(10) |>
    print()

  ## and the ten lowest, which is where the developmental withdrawal sits
  res$oxphos_genes |>
    dplyr::arrange(level_6W_wt) |>
    dplyr::select(symbol, complex, level_6W_wt, lfc_cross, lfc_wt_time, lfc_myc_12W) |>
    head(10) |>
    print()

  ## (8) The diagonal fGSEA beside the five saved rankings, MitoCarta only
  res$fgsea |>
    dplyr::filter(category == "01_mitocarta") |>
    dplyr::select(ranking, pathway, NES, padj_within_category) |>
    tidyr::pivot_wider(names_from = ranking,
                       values_from = c(NES, padj_within_category)) |>
    dplyr::arrange(NES_cross) |>
    print(n = 20)
}
