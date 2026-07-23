# =============================================================================
# fig01b_mito_content_lfc.R -- LFC companion to fig01 (mito content)
# -----------------------------------------------------------------------------
# Rebuts the compositional reading of the share panel (fig01): the share can
# rise because OTHER genes fall. This panel shows the arm's genes go up in
# ABSOLUTE terms (per-gene raw LFC) and COORDINATELY (the whole distribution
# shifts + most genes individually significant), only in the genotype contrasts.
#
# Same FIVE arms and the SAME four comparisons as fig01's brackets, in per-gene
# log2-fold-change space:
#   Myc@6W  = myc_6W_raw        (genotype at 6W)    -- fig01 bracket C1
#   Myc@12W = myc_12W_raw       (genotype at 12W)   -- fig01 bracket C3
#   Time WT = timepoint_neg_raw (6W->12W in WT)     -- fig01 bracket C2
#   Time Myc+= timepoint_pos_raw(6W->12W in Myc+)   -- fig01 bracket C4
# The two Time columns are kept as a deliberate preview of the attenuation /
# mtDNA-over-time stories (developed in later figures).
#
# Reads (read-only, already on disk):
#   results/interaction_results.rds -- the four contrasts as RAW (unshrunken)
#       DESeqResults (per-gene log2FoldChange + padj). RAW is required for
#       averaged-LFC visuals (CLAUDE.md). LFC here == combined_df_raw (cor=1).
#   results/combined_df_annotated_raw.rds -- symbol<->ensembl dictionary.
#   results/mito_content_proxies.rds      -- mass_per_gene, for the mass roster.
#   data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt -- the four
#       MitoCarta arm memberships (identical to fig01; symbol match reproduces
#       script 32's counts: nuclear 1027, oxphos_nu 134, ribo 83, mtDNA 13).
#
# Significance shown = COUNT of individually significant genes (padj<0.05), split
# up (red) / down (blue), over the set size (n) -- concrete and direction-clear.
# The formal set-level competitive test (fGSEA / NES vs a matched background) is
# deferred (used later for the attenuation contrast).
# =============================================================================

source(here::here("figures", "theme_myc.R"))
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("fig01b needs DESeq2 to coerce the DESeqResults in interaction_results.rds")
}

out_dir <- here::here("outputs", "figures")

ir       <- readRDS(here::here("results", "interaction_results.rds"))
combined <- readRDS(here::here("results", "combined_df_annotated_raw.rds"))
content  <- readRDS(here::here("results", "mito_content_proxies.rds"))

# --- arm memberships: 4 from the GMT + the mass roster (chaperone-free) --------
read_gmt <- function(path) {
  ln <- strsplit(readLines(path), "\t")
  stats::setNames(lapply(ln, function(x) x[-(1:2)]),
                  vapply(ln, `[`, character(1), 1))
}
gmt <- read_gmt(here::here("data", "genesets_from_library",
                           "mammary_mito_myc_metab_v1_mouse.gmt"))
mass_nochap <- content$mass_per_gene$gene[content$mass_per_gene$arm != "chaperone"]

arms <- c("MASS_MARKERS_NOCHAP",
          "MITOCARTA_NUCLEAR_ENCODED",
          "MITOCARTA_OXPHOS_NU",
          "MITOCARTA_MITOCHONDRIAL_RIBOSOME",
          "MITOCARTA_MTDNA_ENCODED")
arm_name <- c(MASS_MARKERS_NOCHAP              = "Mass markers",
              MITOCARTA_NUCLEAR_ENCODED        = "Nuclear MitoCarta",
              MITOCARTA_OXPHOS_NU              = "Nuclear OXPHOS",
              MITOCARTA_MITOCHONDRIAL_RIBOSOME = "Mitoribosome",
              MITOCARTA_MTDNA_ENCODED          = "mtDNA-encoded")
arm_syms <- list(
  MASS_MARKERS_NOCHAP              = mass_nochap,
  MITOCARTA_NUCLEAR_ENCODED        = gmt[["MITOCARTA_NUCLEAR_ENCODED"]],
  MITOCARTA_OXPHOS_NU              = gmt[["MITOCARTA_OXPHOS_NU"]],
  MITOCARTA_MITOCHONDRIAL_RIBOSOME = gmt[["MITOCARTA_MITOCHONDRIAL_RIBOSOME"]],
  MITOCARTA_MTDNA_ENCODED          = gmt[["MITOCARTA_MTDNA_ENCODED"]])

# symbol -> ensembl dictionary (arm sets are by symbol; DESeqResults are ensembl)
dict    <- unique(combined[!is.na(combined$mgi_symbol), c("mgi_symbol", "gene")])
arm_ens <- lapply(arm_syms, function(s) unique(dict$gene[dict$mgi_symbol %in% s]))

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

# --- long per-gene table, arm x contrast (lfc + padj) ------------------------
long <- do.call(rbind, lapply(arms, function(a) {
  do.call(rbind, lapply(names(contr), function(cc) {
    tc <- tcs[[cc]][tcs[[cc]]$ensembl %in% arm_ens[[a]], ]
    data.frame(arm = a, contrast = contr[[cc]], lfc = tc$lfc, padj = tc$padj,
               stringsAsFactors = FALSE)
  }))
}))
long <- long[!is.na(long$lfc), ]

# set size (genes present in the DE table) for the facet strip
Ntot <- vapply(arms, function(a) sum(tcs[["myc_6W_raw"]]$ensembl %in% arm_ens[[a]]),
               integer(1))
arm_lab <- sprintf("%s\n(%d genes)", arm_name[arms], Ntot)
long$arm      <- factor(long$arm, levels = arms, labels = arm_lab)
long$contrast <- factor(long$contrast, levels = unname(contr))

# --- per arm x contrast: coordinate direction (% up) + significant up/down ----
#   pct_up  = fraction of set genes with LFC>0  (direction; power-independent --
#             tracks the box shift, so 12W doesn't look weak under FDR power loss)
#   n_up/dn = genes individually padj<0.05, up / down (stringent per-gene lens)
stat <- long |>
  dplyr::group_by(arm, contrast) |>
  dplyr::summarise(
    pct_up = mean(lfc > 0, na.rm = TRUE) * 100,
    n_up   = sum(padj < 0.05 & lfc > 0, na.rm = TRUE),
    n_dn   = sum(padj < 0.05 & lfc < 0, na.rm = TRUE),
    .groups = "drop")

YLIM <- c(-1.5, 3.05)                          # per-gene LFC outliers clipped
stat$pct_lab <- sprintf("%.0f%%", stat$pct_up)
stat$y_pct   <- YLIM[2] - 0.15                 # % up   (grey)  -- coordinate shift
stat$y_up    <- YLIM[2] - 0.52                 # sig up (red)
stat$y_dn    <- YLIM[2] - 0.87                 # sig dn (blue)

p <- ggplot2::ggplot(long, ggplot2::aes(contrast, lfc)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = contrast), outlier.shape = NA,
                        width = 0.62, alpha = 0.6, colour = "grey30", linewidth = 0.3) +
  # annotation stack (top -> bottom): % up (grey), sig up (red), sig down (blue)
  ggplot2::geom_text(data = stat, ggplot2::aes(contrast, y_pct, label = pct_lab),
                     colour = "grey25", size = 1.8, inherit.aes = FALSE) +
  ggplot2::geom_text(data = stat, ggplot2::aes(contrast, y_up, label = n_up),
                     colour = "#D73027", size = 1.8, inherit.aes = FALSE) +
  ggplot2::geom_text(data = stat, ggplot2::aes(contrast, y_dn, label = n_dn),
                     colour = "#2166AC", size = 1.8, inherit.aes = FALSE) +
  ggplot2::facet_wrap(~ arm, nrow = 1) +
  ggplot2::coord_cartesian(ylim = YLIM) +
  ggplot2::scale_fill_manual(
    values = contrast_cols, name = NULL,
    labels = c("Myc@6W" = "Myc effect at 6W", "Myc@12W" = "Myc effect at 12W",
               "Time WT" = "WT 6->12W", "Time Myc+" = "Myc+ 6->12W")) +
  ggplot2::labs(
    x = NULL, y = "per-gene raw log2 fold-change",
    title = "Myc coordinately induces the nuclear arm genes; mtDNA-encoded does not",
    subtitle = "Per contrast top-to-bottom: % genes LFC>0 (grey); padj<0.05 up (red) / down (blue); n per facet.",
    caption = paste(
      "Per-gene RAW (unshrunken) log2FC; box = median/IQR, dashed line = 0. Same four comparisons as fig01's brackets, in LFC space; per-contrast padj from interaction_results (raw).",
      "The two Time contrasts preview the attenuation / mtDNA-over-time stories. Set-level competitive significance (fGSEA / NES vs a matched background) is deferred. LFC outliers beyond the axis are clipped.",
      sep = "\n")) +
  theme_myc(base_size = 9) +
  ggplot2::theme(
    axis.text.x     = ggplot2::element_blank(),
    axis.ticks.x    = ggplot2::element_blank(),
    axis.line.x     = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.text     = ggplot2::element_text(size = 7),
    legend.key.size = ggplot2::unit(3.5, "mm")) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, override.aes = list(alpha = 0.6)))

# Guard: sourced only to obtain `p` (e.g. Quarto) when myc.fig.nosave = TRUE.
if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "fig01b_mito_content_lfc.pdf"),
             width = fig_w[["double"]], height = 92)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(as.data.frame(stat))
  print(rbind(members = vapply(arm_syms, length, integer(1)),
              in_table = Ntot))
  print(p)
}
