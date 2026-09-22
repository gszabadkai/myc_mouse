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
#       script 32's counts: nuclear 1027, oxphos_nu 134, mtDNA 13).
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
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("fig01b needs patchwork to compose the significance-bar track over the boxplots")
}

# NOTE (nesting, kept as-is by author decision 2026-07-23): Nuclear MitoCarta is an
# UMBRELLA -- it fully contains OXPHOS_NU (155/155), the biogenesis arm (314/314) and the
# mass markers (12/12). So those facets overlap the Nuclear MitoCarta facet. mtDNA-encoded
# is disjoint (0/13). A future "Other nuclear MitoCarta" partition is deferred (would also
# need to change fig01 for consistency).
# The two CLAIM-bearing facets are nevertheless effectively disjoint of each other:
# BIOGENESIS_FULL n MITOCARTA_OXPHOS_NU = 7 genes (2% of 314). That is the point of the
# arm -- OXPHOS vs GENERAL biogenesis, which the previous mitoribosome facet (a sub-arm of
# biogenesis translation) could not deliver.

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
          "BIOGENESIS_FULL",
          "MITOCARTA_MTDNA_ENCODED")
arm_name <- c(MASS_MARKERS_NOCHAP       = "Mass markers",
              MITOCARTA_NUCLEAR_ENCODED = "Nuclear MitoCarta",
              MITOCARTA_OXPHOS_NU       = "Nuclear OXPHOS",
              BIOGENESIS_FULL           = "Mito biogenesis",
              MITOCARTA_MTDNA_ENCODED   = "mtDNA-encoded")
# BIOGENESIS_FULL is built by the SAME union expression as script 32's roster, so the
# share panel (fig01) and this LFC panel cannot drift apart.
arm_syms <- list(
  MASS_MARKERS_NOCHAP       = mass_nochap,
  MITOCARTA_NUCLEAR_ENCODED = gmt[["MITOCARTA_NUCLEAR_ENCODED"]],
  MITOCARTA_OXPHOS_NU       = gmt[["MITOCARTA_OXPHOS_NU"]],
  BIOGENESIS_FULL           = union(gmt[["MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA"]],
                                    gmt[["MITOCARTA_PROTEIN_IMPORT_SORTING_AND_HOMEOSTASIS"]]),
  MITOCARTA_MTDNA_ENCODED   = gmt[["MITOCARTA_MTDNA_ENCODED"]])

# arm sets are by symbol (any vintage); DESeqResults are Ensembl. Reconcile via the
# shared helper so renamed genes are not silently dropped -- a plain symbol match
# loses 17 of Complex V's 24 ATP-synthase genes and ~12% of nuclear OXPHOS.
source(here::here("functions", "reconcile_gene_symbols.R"))
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

# =========================== TOP TRACK: significance bars =====================
# Per arm x contrast: a 100% stacked direction bar (up red / down blue by % of
# set genes) with the count of individually significant genes (padj<0.05) inside.
bar <- do.call(rbind, lapply(seq_len(nrow(stat)), function(i) {
  r <- stat[i, ]
  data.frame(arm = r$arm, contrast = r$contrast,
             direction = c("up", "down"),
             frac = c(r$pct_up, 100 - r$pct_up),
             nsig = c(r$n_up, r$n_dn), stringsAsFactors = FALSE)
}))
bar$direction <- factor(bar$direction, levels = c("down", "up"))   # up sits on top

# Label only where it can be read INSIDE its own segment: a 0 carries no information
# and a count printed on a sliver lands outside the fill.
bar$lab <- ifelse(bar$nsig == 0 | bar$frac < 6, "", as.character(bar$nsig))

p_bar <- ggplot2::ggplot(bar, ggplot2::aes(contrast, frac, fill = direction)) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::geom_text(ggplot2::aes(label = lab),
                     position = ggplot2::position_stack(vjust = 0.5),
                     colour = "white", fontface = "bold", size = 2.0) +
  ggplot2::facet_wrap(~ arm, nrow = 1) +
  # muted RdBu mid tones (not the saturated #D73027/#2166AC): softer on the page,
  # still dark enough to carry the bold white counts.
  ggplot2::scale_fill_manual(values = c(up = "#D6604D", down = "#4393C3"), guide = "none") +
  ggplot2::scale_y_continuous(breaks = c(0, 50, 100), expand = ggplot2::expansion(0)) +
  ggplot2::labs(x = NULL, y = "% genes\n(direction)") +
  theme_myc(base_size = 9) +
  ggplot2::theme(
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.line.x  = ggplot2::element_blank())

# =========================== BOTTOM TRACK: LFC boxplots =======================
p_box <- ggplot2::ggplot(long, ggplot2::aes(contrast, lfc)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = contrast), outlier.shape = NA,
                        width = 0.62, alpha = 0.6, colour = "grey30", linewidth = 0.3) +
  ggplot2::facet_wrap(~ arm, nrow = 1) +
  ggplot2::coord_cartesian(ylim = c(-1.4, 1.85)) +
  ggplot2::scale_fill_manual(
    values = contrast_cols, name = NULL,
    labels = c("Myc@6W" = "Myc effect at 6W", "Myc@12W" = "Myc effect at 12W",
               "Time WT" = "WT 6->12W", "Time Myc+" = "Myc+ 6->12W")) +
  ggplot2::labs(x = NULL, y = "per-gene raw log2 fold-change") +
  theme_myc(base_size = 9) +
  ggplot2::theme(
    strip.text      = ggplot2::element_blank(),   # arm labels live on the top track
    axis.text.x     = ggplot2::element_blank(),
    axis.ticks.x    = ggplot2::element_blank(),
    axis.line.x     = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.text     = ggplot2::element_text(size = 7),
    legend.key.size = ggplot2::unit(3.5, "mm")) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, override.aes = list(alpha = 0.6)))

# =========================== COMPOSE (bars atop boxes) ========================
p <- patchwork::wrap_plots(p_bar, p_box, ncol = 1, heights = c(1, 2.9)) +
  patchwork::plot_annotation(
    title = "Myc coordinately induces the nuclear arm genes; mtDNA-encoded does not",
    subtitle = "Top: direction share, up (red) / down (blue); number = genes at padj<0.05. Bottom: per-gene raw log2FC.",
    caption = paste(
      "Same four comparisons as fig01's brackets, in LFC space; RAW LFCs + per-contrast padj from interaction_results.",
      "Time contrasts preview attenuation / mtDNA-over-time. Competitive set-level test (fGSEA/NES vs background) deferred; box outliers clipped.",
      sep = "\n"),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 8.2, colour = "grey20"),
      plot.caption  = ggplot2::element_text(size = 6.3, hjust = 0, colour = "grey30",
                                            lineheight = 1.1)))

# Guard: sourced only to obtain `p` (e.g. Quarto) when myc.fig.nosave = TRUE.
if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "fig01b_mito_content_lfc.pdf"),
             width = fig_w[["double"]], height = 118)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(as.data.frame(stat))
  print(rbind(members = vapply(arm_syms, length, integer(1)),
              in_table = Ntot))
  print(p_bar)   # top track alone
  print(p_box)   # bottom track alone
  print(p)       # composite
}
