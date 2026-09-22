# =============================================================================
# figS1b_mitocarta_survey_lfc.R -- LFC survey across ALL main MitoCarta groups
# -----------------------------------------------------------------------------
# Supplementary companion to figS1 (compartment-share survey). Answers the
# author's question: in the WT 6->12W timeline the TOTAL nuclear MitoCarta and
# the biogenesis arm sit at ~0 LFC while OXPHOS is mostly DOWN -- so which groups
# RISE to hold the total flat? This panel reads that off the "Time WT" column's
# direction share, group by group.
#
# The 16 survey groups = the 7 canonical MitoCarta MitoPathway top-level
# categories, with METABOLISM split into its 9 depth-2 children (hierarchy from
# Sheet 4), plus mtDNA-encoded as the flat reference. Every set is already mt-*
# free (library legacy convention), so MITOCARTA_OXPHOS_NU is the nuclear OXPHOS.
#
# Same FOUR comparisons as fig01/fig01b, in per-gene raw-LFC space:
#   Myc@6W  = myc_6W_raw   (genotype at 6W)   Myc@12W = myc_12W_raw  (genotype at 12W)
#   Time WT = timepoint_neg_raw (6W->12W, WT)  Time Myc+ = timepoint_pos_raw (6W->12W, Myc+)
# Genotype columns are the CLEAN axis; the two Time columns are batch-confounded
# (batch = timepoint) -- read as described, not claimed.
#
# Reads (read-only, already on disk -- NO analysis re-run needed):
#   results/interaction_results.rds        -- the four RAW DESeqResults (CLAUDE.md:
#       averaged-LFC visuals use raw, not shrunken).
#   results/combined_df_annotated_raw.rds  -- symbol <-> ensembl dictionary.
#   data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt -- group memberships.
#
# Also emits the "most divergent genes" supplementary table: per group, the N genes
# whose LFC varies most across the four contrasts (spread = max - min). Written to
# outputs/figures/figS1b_diverse_genes.csv under the same save guard.
# =============================================================================

source(here::here("figures", "theme_myc.R"))
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("figS1b needs DESeq2 to coerce the DESeqResults in interaction_results.rds")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("figS1b needs patchwork to compose the direction track over the boxplots")
}

out_dir    <- here::here("outputs", "figures")
N_DIVERSE  <- 15   # genes per group in the divergence table (author: 10-20)

ir       <- readRDS(here::here("results", "interaction_results.rds"))
combined <- readRDS(here::here("results", "combined_df_annotated_raw.rds"))

read_gmt <- function(path) {
  ln <- strsplit(readLines(path), "\t")
  stats::setNames(lapply(ln, function(x) x[-(1:2)]),
                  vapply(ln, `[`, character(1), 1))
}
gmt <- read_gmt(here::here("data", "genesets_from_library",
                           "mammary_mito_myc_metab_v1_mouse.gmt"))

# --- the 16 survey groups: 6 top-level (non-metabolism) + 9 metabolism L2 + ref -
arms <- c(
  # top-level MitoPathways -- OXPHOS split into its nuclear (155) and mtDNA (13,
  # = MITOCARTA_MTDNA_ENCODED) halves; all other sets are already mt-free.
  "MITOCARTA_OXPHOS_NU",
  "MITOCARTA_OXPHOS_MT",
  "MITOCARTA_ELECTRON_CARRIERS",
  "MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA",
  "MITOCARTA_PROTEIN_IMPORT_SORTING_AND_HOMEOSTASIS",
  "MITOCARTA_MITOCHONDRIAL_DYNAMICS_AND_SURVEILLANCE",
  "MITOCARTA_SIGNALING",
  "MITOCARTA_SMALL_MOLECULE_TRANSPORT",
  # Metabolism, split into its depth-2 children (electron carriers grouped with
  # the respiratory chain above, next to OXPHOS)
  "MITOCARTA_AMINO_ACID_METABOLISM",
  "MITOCARTA_CARBOHYDRATE_METABOLISM",
  "MITOCARTA_LIPID_METABOLISM",
  "MITOCARTA_NUCLEOTIDE_METABOLISM",
  "MITOCARTA_VITAMIN_METABOLISM",
  "MITOCARTA_METALS_AND_COFACTORS",
  "MITOCARTA_DETOXIFICATION",
  "MITOCARTA_SULFUR_METABOLISM")
arm_name <- c(
  MITOCARTA_OXPHOS_NU                              = "Nuclear OXPHOS",
  MITOCARTA_OXPHOS_MT                              = "mtDNA OXPHOS",
  MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA            = "Central dogma",
  MITOCARTA_PROTEIN_IMPORT_SORTING_AND_HOMEOSTASIS = "Import / homeostasis",
  MITOCARTA_MITOCHONDRIAL_DYNAMICS_AND_SURVEILLANCE= "Dynamics & surveillance",
  MITOCARTA_SIGNALING                              = "Signaling",
  MITOCARTA_SMALL_MOLECULE_TRANSPORT               = "Small-molecule transport",
  MITOCARTA_AMINO_ACID_METABOLISM                  = "Amino acid metab.",
  MITOCARTA_CARBOHYDRATE_METABOLISM                = "Carbohydrate metab.",
  MITOCARTA_LIPID_METABOLISM                       = "Lipid metab.",
  MITOCARTA_NUCLEOTIDE_METABOLISM                  = "Nucleotide metab.",
  MITOCARTA_VITAMIN_METABOLISM                     = "Vitamin metab.",
  MITOCARTA_METALS_AND_COFACTORS                   = "Metals & cofactors",
  MITOCARTA_DETOXIFICATION                         = "Detoxification",
  MITOCARTA_SULFUR_METABOLISM                      = "Sulfur metab.",
  MITOCARTA_ELECTRON_CARRIERS                      = "Electron carriers")
stopifnot(all(arms %in% names(gmt)))
arm_syms <- stats::setNames(lapply(arms, function(a) gmt[[a]]), arms)

# group sets are by symbol (any vintage); DESeqResults are ensembl. Reconcile via the
# shared helper so renamed genes (esp. ATP synthase) are not dropped. dict is kept
# only to LABEL genes (ensembl -> current symbol) in the divergence table.
source(here::here("functions", "reconcile_gene_symbols.R"))
dict    <- unique(combined[!is.na(combined$mgi_symbol), c("mgi_symbol", "gene")])
de_univ <- rownames(as.data.frame(ir[["myc_6W_raw"]]))
arm_ens <- lapply(arm_syms, function(s) recon_to_ensembl(s, de_univ))

# --- the four contrasts (raw), tidied to gene x {lfc, padj} -------------------
contr <- c(myc_6W_raw        = "Myc@6W",
           myc_12W_raw       = "Myc@12W",
           timepoint_neg_raw = "Time WT",
           timepoint_pos_raw = "Time Myc+")
contrast_cols <- c("Myc@6W"    = "#F4A582", "Myc@12W"   = "#B2182B",   # genotype = red
                   "Time WT"   = "#BDBDBD", "Time Myc+" = "#7B7B7B")   # time = grey
tcs <- lapply(names(contr), function(nm) {
  r <- as.data.frame(ir[[nm]])
  data.frame(ensembl = rownames(r), lfc = r$log2FoldChange, padj = r$padj)
})
names(tcs) <- names(contr)
lfc_vec <- lapply(tcs, function(t) stats::setNames(t$lfc, t$ensembl))

# --- long per-gene table, group x contrast (lfc + padj) ----------------------
long <- do.call(rbind, lapply(arms, function(a) {
  do.call(rbind, lapply(names(contr), function(cc) {
    tc <- tcs[[cc]][tcs[[cc]]$ensembl %in% arm_ens[[a]], ]
    data.frame(arm = a, contrast = contr[[cc]], lfc = tc$lfc, padj = tc$padj,
               stringsAsFactors = FALSE)
  }))
}))
long <- long[!is.na(long$lfc), ]

# set size (genes present in the DE table) for the boxplot facet strip
Ntot    <- vapply(arms, function(a) sum(tcs[["myc_6W_raw"]]$ensembl %in% arm_ens[[a]]),
                  integer(1))
arm_lab <- stats::setNames(sprintf("%s\n(%d genes)", arm_name[arms], Ntot), arms)
long$arm_f    <- factor(long$arm, levels = arms, labels = arm_lab[arms])   # boxplot strips
long$contrast <- factor(long$contrast, levels = unname(contr))

# --- per group x contrast: median LFC + coordinate direction + sig up/down ---
# (16 groups x 4 contrasts does NOT fit the fig01b bar-atop-box composite -- that
# aligns only for a single facet row. The summary is a heatmap; the boxplot grid
# carries the distributions.)
stat <- long |>
  dplyr::group_by(arm, contrast) |>
  dplyr::summarise(
    med_lfc = stats::median(lfc, na.rm = TRUE),
    pct_up  = mean(lfc > 0, na.rm = TRUE) * 100,
    n_up    = sum(padj < 0.05 & lfc > 0, na.rm = TRUE),
    n_dn    = sum(padj < 0.05 & lfc < 0, na.rm = TRUE),
    .groups = "drop")

# =========================== SUMMARY: median-LFC heatmap =====================
# rows = the 16 groups (OXPHOS on top), cols = the four contrasts; fill = median
# raw LFC (diverging at 0); text = count of individually significant genes
# (up in white / down as -n) so power and direction are both visible.
LIM <- 0.6   # fill clipped so a few strong groups don't wash the rest out
hm <- stat
hm$group <- factor(unname(arm_name[hm$arm]), levels = rev(unname(arm_name[arms])))
hm$fillv <- pmax(pmin(hm$med_lfc, LIM), -LIM)
hm$sig   <- ifelse(hm$n_up + hm$n_dn == 0, "",
                   ifelse(hm$n_up >= hm$n_dn, sprintf("%d", hm$n_up),
                          sprintf("-%d", hm$n_dn)))

p_heat <- ggplot2::ggplot(hm, ggplot2::aes(contrast, group, fill = fillv)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
  ggplot2::geom_text(ggplot2::aes(label = sig), size = 2.2, colour = "grey15") +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-LIM, LIM), name = "median\nlog2FC",
    breaks = c(-LIM, 0, LIM), labels = c(paste0("<=-", LIM), "0", paste0(">=", LIM))) +
  ggplot2::scale_x_discrete(
    position = "top",
    labels = c("Myc@6W" = "Myc @6W", "Myc@12W" = "Myc @12W",
               "Time WT" = "WT\n6->12W", "Time Myc+" = "Myc+\n6->12W")) +
  ggplot2::labs(x = NULL, y = NULL) +
  theme_myc(base_size = 8) +
  ggplot2::theme(
    axis.line       = ggplot2::element_blank(),
    axis.ticks      = ggplot2::element_blank(),
    axis.text.x.top = ggplot2::element_text(face = "bold", size = 7.5, lineheight = 0.9),
    axis.text.y     = ggplot2::element_text(size = 7.5),
    legend.position = "right",
    legend.key.size = ggplot2::unit(3.5, "mm"),
    legend.title    = ggplot2::element_text(size = 7))

# =========================== DETAIL: LFC boxplot grid ========================
p_box <- ggplot2::ggplot(long, ggplot2::aes(contrast, lfc)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = contrast), outlier.shape = NA,
                        width = 0.62, alpha = 0.6, colour = "grey30", linewidth = 0.3) +
  ggplot2::facet_wrap(~ arm_f, ncol = 4) +
  ggplot2::coord_cartesian(ylim = c(-1.6, 1.6)) +
  ggplot2::scale_fill_manual(
    values = contrast_cols, name = NULL,
    labels = c("Myc@6W" = "Myc effect at 6W", "Myc@12W" = "Myc effect at 12W",
               "Time WT" = "WT 6->12W", "Time Myc+" = "Myc+ 6->12W")) +
  ggplot2::labs(x = NULL, y = "per-gene raw log2 fold-change") +
  theme_myc(base_size = 8) +
  ggplot2::theme(
    strip.text      = ggplot2::element_text(size = 6.6, lineheight = 0.9),
    axis.text.x     = ggplot2::element_blank(),
    axis.ticks.x    = ggplot2::element_blank(),
    axis.line.x     = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.text     = ggplot2::element_text(size = 7),
    legend.key.size = ggplot2::unit(3.5, "mm")) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, override.aes = list(alpha = 0.6)))

# =========================== COMPOSE (A heatmap / B boxplots) ================
p <- patchwork::wrap_plots(p_heat, p_box, ncol = 1, heights = c(1, 1.45)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Which mitochondrial groups rise while OXPHOS falls? -- a MitoPathways survey",
    subtitle = "A: median per-gene log2FC (n = sig genes padj<0.05, -n = down). B: LFC distributions.",
    caption = paste(
      "16 MitoCarta MitoPathway groups (Metabolism split into its 9 depth-2 children) x the four fig01 contrasts, in LFC space.",
      "RAW LFCs + per-contrast padj from interaction_results. Genotype columns are clean; Time columns are batch-confounded (batch=timepoint).",
      sep = "\n"),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 8.2, colour = "grey20"),
      plot.caption  = ggplot2::element_text(size = 6.3, hjust = 0, colour = "grey30",
                                            lineheight = 1.1)))

# =========================== DIVERGENCE TABLE ================================
# Per group, the N genes whose LFC varies most ACROSS the four contrasts
# (spread = max - min; sd kept alongside) -- surfaces genes that behave
# differently between genotype and time (the context-dependent movers).
diverse_genes <- do.call(rbind, lapply(arms, function(a) {
  e <- arm_ens[[a]]
  m <- vapply(names(contr), function(cc) unname(lfc_vec[[cc]][e]), numeric(length(e)))
  rownames(m) <- e
  spread <- apply(m, 1, function(x) { x <- x[is.finite(x)]
                                      if (length(x) < 2) NA_real_ else max(x) - min(x) })
  sdv    <- apply(m, 1, function(x) { x <- x[is.finite(x)]
                                      if (length(x) < 2) NA_real_ else stats::sd(x) })
  ord  <- order(spread, decreasing = TRUE)
  keep <- utils::head(ord[is.finite(spread[ord])], N_DIVERSE)
  data.frame(
    group          = unname(arm_name[[a]]),
    gene           = dict$mgi_symbol[match(e[keep], dict$gene)],
    lfc_Myc6W      = m[keep, "myc_6W_raw"],
    lfc_Myc12W     = m[keep, "myc_12W_raw"],
    lfc_TimeWT     = m[keep, "timepoint_neg_raw"],
    lfc_TimeMycPos = m[keep, "timepoint_pos_raw"],
    spread         = spread[keep],
    sd             = sdv[keep],
    row.names      = NULL, stringsAsFactors = FALSE)
}))
diverse_genes$group <- factor(diverse_genes$group, levels = unname(arm_name[arms]))

# Guard: sourced only to obtain `p` / `diverse_genes` (e.g. Quarto) when nosave = TRUE.
if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "figS1b_mitocarta_survey_lfc.pdf"),
             width = fig_w[["double"]], height = 165)
  utils::write.csv(diverse_genes, file.path(out_dir, "figS1b_diverse_genes.csv"),
                   row.names = FALSE)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(rbind(members = vapply(arm_syms, length, integer(1)), in_table = Ntot))
  # which groups rise in the WT timeline while OXPHOS falls?
  stat |>
    dplyr::filter(contrast == "Time WT") |>
    dplyr::arrange(dplyr::desc(pct_up)) |>
    as.data.frame() |> print()
  print(subset(diverse_genes, group == "OXPHOS"))
  print(subset(diverse_genes, group == "Lipid metab."))
  print(p_heat); print(p_box); print(p)
}
