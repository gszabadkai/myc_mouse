# =============================================================================
# figS2b_oxphos_complex_lfc.R -- LFC breakdown of the five OXPHOS complexes
# -----------------------------------------------------------------------------
# The granular companion to figS1/figS1b: having seen that NUCLEAR OXPHOS is the
# arm that falls across the WT timeline (and rises with Myc), which of the five
# respiratory complexes (CI-CV) carries it? Same two views as figS1b, in LFC space.
#
# Groups = MITOCARTA_COMPLEX_I..V. All are NUCLEAR-encoded (mt-* removed by the
# library convention): the 13 mtDNA-encoded subunits (mt-Nd*/mt-Cytb/mt-Co*/mt-Atp*)
# are the mtDNA OXPHOS facet of figS1, not here. Complex II is fully nuclear.
# Complex = subunits + assembly factors (CI 57, CII 8, CIII 15, CIV 47, CV 24).
#
# Same FOUR comparisons as fig01/figS1b, per-gene raw log2FC:
#   Myc@6W / Myc@12W (genotype, clean) ; Time WT / Time Myc+ (6W->12W, batch-confounded).
#
# Reads (read-only; NO analysis re-run needed):
#   results/interaction_results.rds        -- the four RAW DESeqResults.
#   results/combined_df_annotated_raw.rds  -- symbol <-> ensembl dictionary.
#   data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt -- complex memberships.
#
# Also emits the divergence table (per complex, the genes whose LFC varies most
# across the four contrasts) -> outputs/figures/figS2b_diverse_genes.csv.
# =============================================================================

source(here::here("figures", "theme_myc.R"))
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("figS2b needs DESeq2 to coerce the DESeqResults in interaction_results.rds")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("figS2b needs patchwork to compose the heatmap over the boxplots")
}
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("figS2b needs readxl to map MitoCarta symbols to Ensembl (see mapping note)")
}

out_dir   <- here::here("outputs", "figures")
N_DIVERSE <- 12   # genes per complex in the divergence table (small sets cap below this)

ir       <- readRDS(here::here("results", "interaction_results.rds"))
combined <- readRDS(here::here("results", "combined_df_annotated_raw.rds"))

read_gmt <- function(path) {
  ln <- strsplit(readLines(path), "\t")
  stats::setNames(lapply(ln, function(x) x[-(1:2)]),
                  vapply(ln, `[`, character(1), 1))
}
gmt <- read_gmt(here::here("data", "genesets_from_library",
                           "mammary_mito_myc_metab_v1_mouse.gmt"))

arms <- c("MITOCARTA_COMPLEX_I", "MITOCARTA_COMPLEX_II", "MITOCARTA_COMPLEX_III",
          "MITOCARTA_COMPLEX_IV", "MITOCARTA_COMPLEX_V")
arm_name <- c(MITOCARTA_COMPLEX_I   = "Complex I",
              MITOCARTA_COMPLEX_II  = "Complex II",
              MITOCARTA_COMPLEX_III = "Complex III",
              MITOCARTA_COMPLEX_IV  = "Complex IV",
              MITOCARTA_COMPLEX_V   = "Complex V")
stopifnot(all(arms %in% names(gmt)))
arm_syms <- stats::setNames(lapply(arms, function(a) gmt[[a]]), arms)

# --- MAPPING NOTE: MitoCarta symbol -> Ensembl (NOT current-symbol match) ------
# MitoCarta's complex sets carry the OLD ATP-synthase names (Atp5a1/b/c1/d/e/g1-3/
# h/j/k/l/o, Atpif1); the DESeq table uses the CURRENT names (Atp5f1a..e, Atp5mc*,
# Atp5pd/pf/po, Atp5if1). A current-symbol match silently drops 17 of Complex V's
# 24 genes. So membership is resolved by MitoCarta's OWN Ensembl IDs (Sheet 2,
# EnsemblGeneID) intersected with the DE table -- recovering CV 24/24. dict below
# is used ONLY to LABEL genes (ensembl -> current symbol) in the divergence table.
dict   <- unique(combined[!is.na(combined$mgi_symbol), c("mgi_symbol", "gene")])
de_ens <- rownames(as.data.frame(ir[["myc_6W_raw"]]))
mc <- as.data.frame(readxl::read_excel(
  here::here("data", "Mouse.MitoCarta3.0.xls"), sheet = 2))
mc <- mc[!is.na(mc$Symbol) & !is.na(mc$EnsemblGeneID) & mc$EnsemblGeneID != "", ]
mc_map <- do.call(rbind, lapply(seq_len(nrow(mc)), function(i)
  data.frame(sym = mc$Symbol[i],
             ens = trimws(strsplit(mc$EnsemblGeneID[i], "[|]")[[1]]),
             stringsAsFactors = FALSE)))
arm_ens <- lapply(arm_syms, function(s)
  intersect(unique(mc_map$ens[mc_map$sym %in% s]), de_ens))

# --- the four contrasts (raw), tidied to gene x {lfc, padj} -------------------
contr <- c(myc_6W_raw        = "Myc@6W",
           myc_12W_raw       = "Myc@12W",
           timepoint_neg_raw = "Time WT",
           timepoint_pos_raw = "Time Myc+")
contrast_cols <- c("Myc@6W" = "#F4A582", "Myc@12W" = "#B2182B",
                   "Time WT" = "#BDBDBD", "Time Myc+" = "#7B7B7B")
tcs <- lapply(names(contr), function(nm) {
  r <- as.data.frame(ir[[nm]])
  data.frame(ensembl = rownames(r), lfc = r$log2FoldChange, padj = r$padj)
})
names(tcs) <- names(contr)
lfc_vec <- lapply(tcs, function(t) stats::setNames(t$lfc, t$ensembl))

# --- long per-gene table, complex x contrast --------------------------------
long <- do.call(rbind, lapply(arms, function(a) {
  do.call(rbind, lapply(names(contr), function(cc) {
    tc <- tcs[[cc]][tcs[[cc]]$ensembl %in% arm_ens[[a]], ]
    data.frame(arm = a, contrast = contr[[cc]], lfc = tc$lfc, padj = tc$padj,
               stringsAsFactors = FALSE)
  }))
}))
long <- long[!is.na(long$lfc), ]

Ntot    <- vapply(arms, function(a) sum(tcs[["myc_6W_raw"]]$ensembl %in% arm_ens[[a]]),
                  integer(1))
arm_lab <- stats::setNames(sprintf("%s\n(%d genes)", arm_name[arms], Ntot), arms)
long$arm_f    <- factor(long$arm, levels = arms, labels = arm_lab[arms])
long$contrast <- factor(long$contrast, levels = unname(contr))

stat <- long |>
  dplyr::group_by(arm, contrast) |>
  dplyr::summarise(
    med_lfc = stats::median(lfc, na.rm = TRUE),
    pct_up  = mean(lfc > 0, na.rm = TRUE) * 100,
    n_up    = sum(padj < 0.05 & lfc > 0, na.rm = TRUE),
    n_dn    = sum(padj < 0.05 & lfc < 0, na.rm = TRUE),
    .groups = "drop")

# =========================== SUMMARY: median-LFC heatmap =====================
LIM <- 0.6
hm <- stat
hm$group <- factor(unname(arm_name[hm$arm]), levels = rev(unname(arm_name[arms])))
hm$fillv <- pmax(pmin(hm$med_lfc, LIM), -LIM)
hm$sig   <- ifelse(hm$n_up + hm$n_dn == 0, "",
                   ifelse(hm$n_up >= hm$n_dn, sprintf("%d", hm$n_up),
                          sprintf("-%d", hm$n_dn)))

p_heat <- ggplot2::ggplot(hm, ggplot2::aes(contrast, group, fill = fillv)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
  ggplot2::geom_text(ggplot2::aes(label = sig), size = 2.3, colour = "grey15") +
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
    axis.text.y     = ggplot2::element_text(size = 8),
    legend.position = "right",
    legend.key.size = ggplot2::unit(3.5, "mm"),
    legend.title    = ggplot2::element_text(size = 7))

# =========================== DETAIL: LFC boxplots ============================
p_box <- ggplot2::ggplot(long, ggplot2::aes(contrast, lfc)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = contrast), outlier.shape = NA,
                        width = 0.62, alpha = 0.6, colour = "grey30", linewidth = 0.3) +
  ggplot2::facet_wrap(~ arm_f, nrow = 1) +
  ggplot2::coord_cartesian(ylim = c(-1.6, 1.6)) +
  ggplot2::scale_fill_manual(
    values = contrast_cols, name = NULL,
    labels = c("Myc@6W" = "Myc effect at 6W", "Myc@12W" = "Myc effect at 12W",
               "Time WT" = "WT 6->12W", "Time Myc+" = "Myc+ 6->12W")) +
  ggplot2::labs(x = NULL, y = "per-gene raw log2 fold-change") +
  theme_myc(base_size = 8) +
  ggplot2::theme(
    strip.text      = ggplot2::element_text(size = 7, lineheight = 0.9),
    axis.text.x     = ggplot2::element_blank(),
    axis.ticks.x    = ggplot2::element_blank(),
    axis.line.x     = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.text     = ggplot2::element_text(size = 7),
    legend.key.size = ggplot2::unit(3.5, "mm")) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, override.aes = list(alpha = 0.6)))

p <- patchwork::wrap_plots(p_heat, p_box, ncol = 1, heights = c(1, 1.5)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "The five OXPHOS complexes: which carries the nuclear-OXPHOS signal?",
    subtitle = "A: median per-gene log2FC (n = sig genes padj<0.05, -n = down). B: LFC distributions.",
    caption = paste(
      "MITOCARTA_COMPLEX_I..V (nuclear-encoded subunits + assembly factors; the 13 mtDNA subunits are the mtDNA OXPHOS facet of figS1).",
      "RAW LFCs + per-contrast padj from interaction_results. Genotype columns clean; Time columns batch-confounded (batch=timepoint).",
      sep = "\n"),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 8.2, colour = "grey20"),
      plot.caption  = ggplot2::element_text(size = 6.3, hjust = 0, colour = "grey30",
                                            lineheight = 1.1)))

# =========================== DIVERGENCE TABLE ===============================
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
    complex        = unname(arm_name[[a]]),
    gene           = dict$mgi_symbol[match(e[keep], dict$gene)],
    lfc_Myc6W      = m[keep, "myc_6W_raw"],
    lfc_Myc12W     = m[keep, "myc_12W_raw"],
    lfc_TimeWT     = m[keep, "timepoint_neg_raw"],
    lfc_TimeMycPos = m[keep, "timepoint_pos_raw"],
    spread         = spread[keep],
    sd             = sdv[keep],
    row.names      = NULL, stringsAsFactors = FALSE)
}))
diverse_genes$complex <- factor(diverse_genes$complex, levels = unname(arm_name[arms]))

if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "figS2b_oxphos_complex_lfc.pdf"),
             width = fig_w[["double"]], height = 118)
  utils::write.csv(diverse_genes, file.path(out_dir, "figS2b_diverse_genes.csv"),
                   row.names = FALSE)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(rbind(members = vapply(arm_syms, length, integer(1)), in_table = Ntot))
  stat |> dplyr::filter(contrast == "Time WT") |>
    dplyr::arrange(dplyr::desc(pct_up)) |> as.data.frame() |> print()
  print(subset(diverse_genes, complex == "Complex I"))
  print(p_heat); print(p_box); print(p)
}
