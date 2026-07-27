# =============================================================================
# figureS1_compartment_detail.R -- SUPPLEMENTARY FIGURE 1
# "The mitochondrial compartment in detail" -- the controls behind Figure 1.
# -----------------------------------------------------------------------------
# Figure 1 makes three claims about the compartment. This figure is the evidence
# a reader needs to believe them, in the same three orders:
#
#   A  THE SHARE RISES BECAUSE THE GENES RISE. The obvious objection to a
#      compartment SHARE is that a fraction can rise because everything else
#      falls. The same four contrasts in per-gene log2FC space answer it: the
#      genes of every nuclear arm rise ABSOLUTELY and COHERENTLY, and the 13
#      mtDNA-encoded transcripts do not. The share is not arithmetic.
#   B  WHAT RISES WHEN OXPHOS FALLS. The whole nuclear compartment is roughly
#      flat across the wild-type window while OXPHOS falls -- so something must
#      rise to hold the total. The 16-group survey names it: the substrate-
#      metabolism arms. A flat total is REALLOCATION, not stasis.
#   C  THE RESPIRATORY CHAIN MOVES AS ONE BLOCK. If MYC raised one complex and
#      not the others, "MYC raises OXPHOS" would be the wrong sentence. All five
#      move together on the genotype axis, and all five attenuate together.
#      Complex II is the small underpowered exception, and it is the one complex
#      with no mtDNA-encoded subunit.
#
# MEMBERSHIP GOES THROUGH THE RECONCILER, always. The library GMTs carry the
# ORIGINAL MitoCarta symbols (Atp5a1, Atp5b, ...) while the DESeq table carries
# the current ones; a plain symbol match silently drops 17 of Complex V's 24
# genes and ~12% of nuclear OXPHOS. `recon_to_ensembl()` recovers them.
#
# READ THE `% up` LINE, NOT ONLY THE SIGNIFICANT-GENE COUNTS. Counts track set
# size and power; the direction fraction does not. At 12W the per-gene magnitude
# is attenuated and almost nothing clears FDR even where the compartment share is
# firmest -- the counts would otherwise make 12W look artificially weak.
#
# BATCH = TIMEPOINT: the two temporal columns are confounded and are DESCRIBED,
# not claimed. The genotype columns are clean. n = 6 per group.
#
# Reads (read-only; the author runs scripts 01-12 and 32 first):
#   results/interaction_results.rds        -- the four RAW DESeqResults
#   results/combined_df_annotated_raw.rds  -- ensembl <-> current symbol
#   results/mito_content_proxies.rds       -- the mass roster (panel A)
#   data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt
# =============================================================================

source(here::here("figures", "theme_myc.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))
if (!requireNamespace("patchwork", quietly = TRUE)) stop("figureS1 needs patchwork")
if (!requireNamespace("DESeq2", quietly = TRUE))    stop("figureS1 needs DESeq2")
if (!requireNamespace("dplyr", quietly = TRUE))     stop("figureS1 needs dplyr")

out_dir <- here::here("outputs", "figures")

ir       <- readRDS(here::here("results", "interaction_results.rds"))
combined <- readRDS(here::here("results", "combined_df_annotated_raw.rds"))

read_gmt <- function(path) {
  ln <- strsplit(readLines(path), "\t")
  stats::setNames(lapply(ln, function(x) x[-(1:2)]),
                  vapply(ln, `[`, character(1), 1))
}
gmt <- read_gmt(here::here("data", "genesets_from_library",
                           "mammary_mito_myc_metab_v1_mouse.gmt"))

# --- the four raw contrasts, shared by every panel ---------------------------
contr <- c(myc_6W_raw        = "Myc @6W",
           myc_12W_raw       = "Myc @12W",
           timepoint_neg_raw = "WT 6->12W",
           timepoint_pos_raw = "Myc+ 6->12W")
tcs <- lapply(names(contr), function(nm) {
  r <- as.data.frame(ir[[nm]])
  data.frame(ensembl = rownames(r), lfc = r$log2FoldChange, padj = r$padj,
             stringsAsFactors = FALSE)
})
names(tcs) <- names(contr)
de_ens <- tcs[[1]]$ensembl

# --- one tabulator, used by all three panels ---------------------------------
# Returns median LFC, the direction fraction and the FDR counts per set x contrast.
tabulate_arms <- function(sym_list, labels) {
  set_names <- names(sym_list)
  ens <- stats::setNames(lapply(sym_list, function(s)
    recon_to_ensembl(s, de_ens)), set_names)
  long <- do.call(rbind, lapply(set_names, function(a) {
    do.call(rbind, lapply(names(contr), function(cc) {
      t <- tcs[[cc]][tcs[[cc]]$ensembl %in% ens[[a]], ]
      data.frame(arm = a, contrast = contr[[cc]], lfc = t$lfc, padj = t$padj,
                 stringsAsFactors = FALSE)
    }))
  }))
  long <- long[is.finite(long$lfc), ]
  st <- long |>
    dplyr::group_by(arm, contrast) |>
    dplyr::summarise(med_lfc = stats::median(lfc, na.rm = TRUE),
                     pct_up  = 100 * mean(lfc > 0, na.rm = TRUE),
                     n_up    = sum(padj < 0.05 & lfc > 0, na.rm = TRUE),
                     n_dn    = sum(padj < 0.05 & lfc < 0, na.rm = TRUE),
                     n_set   = dplyr::n(), .groups = "drop") |>
    as.data.frame()
  st$label <- sprintf("%s  (%d)", labels[st$arm], st$n_set[match(st$arm, st$arm)])
  st$label <- factor(st$label,
                     levels = rev(unique(st$label[order(match(st$arm, set_names))])))
  st$contrast <- factor(st$contrast, levels = unname(contr))
  st
}

# --- the shared heatmap grammar ----------------------------------------------
LIM <- 0.6
heat_panel <- function(st, title, subtitle, legend = "none") {
  st$fillv <- pmax(pmin(st$med_lfc, LIM), -LIM)
  st$sig   <- ifelse(st$n_up + st$n_dn == 0, "",
                     ifelse(st$n_up >= st$n_dn, sprintf("%d", st$n_up),
                            sprintf("-%d", st$n_dn)))
  ggplot2::ggplot(st, ggplot2::aes(contrast, label, fill = fillv)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.55) +
    ggplot2::geom_text(ggplot2::aes(label = sig), size = 2.1, colour = "grey15") +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-LIM, LIM), name = "median\nlog2FC",
      breaks = c(-LIM, 0, LIM),
      labels = c(paste0("<=-", LIM), "0", paste0(">=", LIM))) +
    ggplot2::scale_x_discrete(position = "top",
                              labels = c("Myc @6W" = "Myc\n6W", "Myc @12W" = "Myc\n12W",
                                         "WT 6->12W" = "WT\ntime",
                                         "Myc+ 6->12W" = "Myc+\ntime")) +
    ggplot2::labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_myc(base_size = 8) +
    ggplot2::theme(
      axis.line       = ggplot2::element_blank(),
      axis.ticks      = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(face = "bold", size = 6.6,
                                              lineheight = 0.9),
      axis.text.y     = ggplot2::element_text(size = 6.6),
      legend.position = legend,
      legend.key.size = ggplot2::unit(3.2, "mm"),
      legend.title    = ggplot2::element_text(size = 6.5),
      legend.text     = ggplot2::element_text(size = 6),
      plot.title      = ggplot2::element_text(face = "bold", size = 8.5),
      plot.subtitle   = ggplot2::element_text(size = 6.2, colour = "grey25",
                                              lineheight = 1.15))
}

# =============================================================================
# PANEL A -- the share rises because the genes rise
# =============================================================================
armsA <- c("MASS_MARKERS_NOCHAP", "MITOCARTA_NUCLEAR_ENCODED", "MITOCARTA_OXPHOS_NU",
           "BIOGENESIS_FULL", "MITOCARTA_MTDNA_ENCODED")
labA <- c(MASS_MARKERS_NOCHAP       = "Mass markers",
          MITOCARTA_NUCLEAR_ENCODED = "Nuclear MitoCarta",
          MITOCARTA_OXPHOS_NU       = "Nuclear OXPHOS",
          BIOGENESIS_FULL           = "Mito biogenesis",
          MITOCARTA_MTDNA_ENCODED   = "mtDNA-encoded")
# Two of these are NOT GMT sets. `MASS_MARKERS_NOCHAP` is script 32's mass roster
# minus the chaperone arm, and `BIOGENESIS_FULL` is the union script 32 builds --
# reproduced by the SAME expression here so the share panel and this one cannot
# drift apart.
content_32  <- readRDS(here::here("results", "mito_content_proxies.rds"))
mass_nochap <- content_32$mass_per_gene$gene[content_32$mass_per_gene$arm != "chaperone"]
from_gmt <- c("MITOCARTA_NUCLEAR_ENCODED", "MITOCARTA_OXPHOS_NU",
              "MITOCARTA_MTDNA_ENCODED", "MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA",
              "MITOCARTA_PROTEIN_IMPORT_SORTING_AND_HOMEOSTASIS")
missA <- from_gmt[!from_gmt %in% names(gmt)]
if (length(missA))
  stop("figureS1 A: set names absent from the GMT -> ", paste(missA, collapse = ", "))
symsA <- list(
  MASS_MARKERS_NOCHAP       = mass_nochap,
  MITOCARTA_NUCLEAR_ENCODED = gmt[["MITOCARTA_NUCLEAR_ENCODED"]],
  MITOCARTA_OXPHOS_NU       = gmt[["MITOCARTA_OXPHOS_NU"]],
  BIOGENESIS_FULL           = union(gmt[["MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA"]],
                                    gmt[["MITOCARTA_PROTEIN_IMPORT_SORTING_AND_HOMEOSTASIS"]]),
  MITOCARTA_MTDNA_ENCODED   = gmt[["MITOCARTA_MTDNA_ENCODED"]])
symsA <- symsA[armsA]
stA <- tabulate_arms(symsA, labA)

# 100% stacked direction bars: the power-INDEPENDENT read. The FDR count rides
# inside the bar, blanked where it is 0 or the segment is too narrow to hold it.
dirA <- rbind(
  data.frame(stA[, c("arm", "contrast", "label")], dir = "up",
             frac = stA$pct_up, n = stA$n_up),
  data.frame(stA[, c("arm", "contrast", "label")], dir = "down",
             frac = 100 - stA$pct_up, n = stA$n_dn))
dirA$dir <- factor(dirA$dir, levels = c("down", "up"))
dirA$lab <- ifelse(dirA$n == 0 | dirA$frac < 12, "", as.character(dirA$n))

pA <- ggplot2::ggplot(dirA, ggplot2::aes(frac, label, fill = dir)) +
  ggplot2::geom_col(width = 0.72, colour = "white", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = lab), position = ggplot2::position_stack(vjust = 0.5),
                     size = 1.95, colour = "grey15") +
  ggplot2::geom_vline(xintercept = 50, colour = "grey35", linewidth = 0.3, linetype = 2) +
  ggplot2::facet_wrap(~ contrast, nrow = 1) +
  ggplot2::scale_fill_manual(values = c(up = "#D6604D", down = "#4393C3"),
                             labels = c(up = "up", down = "down"), name = NULL,
                             breaks = c("up", "down")) +
  ggplot2::scale_x_continuous(breaks = c(0, 100), expand = c(0, 0)) +
  ggplot2::labs(
    x = "% of the set's genes moving in each direction", y = NULL,
    title = "A   The compartment share rises because the genes rise",
    subtitle = "Per-gene raw log2FC, the answer to 'a fraction can rise because everything else falls'. Numbers inside the bars\nare FDR-significant genes -- they track set size and power; the % up line does not. Dashed line = 50%.") +
  theme_myc(base_size = 8) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.key.size = ggplot2::unit(3, "mm"),
    legend.text     = ggplot2::element_text(size = 6.2),
    axis.text.y     = ggplot2::element_text(size = 6.6),
    axis.text.x     = ggplot2::element_text(size = 6),
    panel.spacing.x = ggplot2::unit(3.2, "mm"),
    strip.text      = ggplot2::element_text(face = "bold", size = 7),
    plot.title      = ggplot2::element_text(face = "bold", size = 8.5),
    plot.subtitle   = ggplot2::element_text(size = 6.2, colour = "grey25",
                                            lineheight = 1.15))

# =============================================================================
# PANEL B -- what rises when OXPHOS falls
# =============================================================================
armsB <- c("MITOCARTA_OXPHOS_NU", "MITOCARTA_OXPHOS_MT", "MITOCARTA_ELECTRON_CARRIERS",
           "MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA",
           "MITOCARTA_PROTEIN_IMPORT_SORTING_AND_HOMEOSTASIS",
           "MITOCARTA_MITOCHONDRIAL_DYNAMICS_AND_SURVEILLANCE",
           "MITOCARTA_SIGNALING", "MITOCARTA_SMALL_MOLECULE_TRANSPORT",
           "MITOCARTA_AMINO_ACID_METABOLISM", "MITOCARTA_CARBOHYDRATE_METABOLISM",
           "MITOCARTA_LIPID_METABOLISM", "MITOCARTA_NUCLEOTIDE_METABOLISM",
           "MITOCARTA_VITAMIN_METABOLISM", "MITOCARTA_METALS_AND_COFACTORS",
           "MITOCARTA_DETOXIFICATION", "MITOCARTA_SULFUR_METABOLISM")
labB <- c(MITOCARTA_OXPHOS_NU                              = "Nuclear OXPHOS",
          MITOCARTA_OXPHOS_MT                              = "mtDNA OXPHOS",
          MITOCARTA_ELECTRON_CARRIERS                      = "Electron carriers",
          MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA            = "Central dogma",
          MITOCARTA_PROTEIN_IMPORT_SORTING_AND_HOMEOSTASIS = "Import / homeostasis",
          MITOCARTA_MITOCHONDRIAL_DYNAMICS_AND_SURVEILLANCE = "Dynamics & surveillance",
          MITOCARTA_SIGNALING                              = "Signaling",
          MITOCARTA_SMALL_MOLECULE_TRANSPORT               = "SM transport",
          MITOCARTA_AMINO_ACID_METABOLISM                  = "Amino acid metab.",
          MITOCARTA_CARBOHYDRATE_METABOLISM                = "Carbohydrate metab.",
          MITOCARTA_LIPID_METABOLISM                       = "Lipid metab.",
          MITOCARTA_NUCLEOTIDE_METABOLISM                  = "Nucleotide metab.",
          MITOCARTA_VITAMIN_METABOLISM                     = "Vitamin metab.",
          MITOCARTA_METALS_AND_COFACTORS                   = "Metals & cofactors",
          MITOCARTA_DETOXIFICATION                         = "Detoxification",
          MITOCARTA_SULFUR_METABOLISM                      = "Sulfur metab.")
missB <- armsB[!armsB %in% names(gmt)]
if (length(missB))
  stop("figureS1 B: set names absent from the GMT -> ", paste(missB, collapse = ", "))
stB <- tabulate_arms(stats::setNames(lapply(armsB, function(a) gmt[[a]]), armsB), labB)

ox_wt <- stB$med_lfc[stB$arm == "MITOCARTA_OXPHOS_NU" & stB$contrast == "WT 6->12W"]
n_rise <- sum(stB$med_lfc[stB$contrast == "WT 6->12W"] > 0)

pB <- heat_panel(
  stB,
  "B   What rises when OXPHOS falls",
  sprintf("Read the WT column: nuclear OXPHOS is the arm that falls while %d of %d groups rise\n-- a flat compartment total is REALLOCATION, not stasis. (Figure 2 reads the OXPHOS\nSUBUNITS specifically, which fall further than the arm as a whole.) Numbers = FDR-sig genes.",
          n_rise, length(armsB)),
  legend = "right")

# =============================================================================
# PANEL C -- the respiratory chain moves as one block
# =============================================================================
armsC <- c("MITOCARTA_COMPLEX_I", "MITOCARTA_COMPLEX_II", "MITOCARTA_COMPLEX_III",
           "MITOCARTA_COMPLEX_IV", "MITOCARTA_COMPLEX_V")
labC <- c(MITOCARTA_COMPLEX_I   = "Complex I",
          MITOCARTA_COMPLEX_II  = "Complex II",
          MITOCARTA_COMPLEX_III = "Complex III",
          MITOCARTA_COMPLEX_IV  = "Complex IV",
          MITOCARTA_COMPLEX_V   = "Complex V")
missC <- armsC[!armsC %in% names(gmt)]
if (length(missC))
  stop("figureS1 C: set names absent from the GMT -> ", paste(missC, collapse = ", "))
stC <- tabulate_arms(stats::setNames(lapply(armsC, function(a) gmt[[a]]), armsC), labC)

# GUARD on the reconciler: a plain symbol match recovers only 7 of Complex V's 24
# genes. If this assertion fails the mapping has regressed and every set-level
# number in this figure is under-counted.
cv_n <- unique(stC$n_set[stC$arm == "MITOCARTA_COMPLEX_V"])
if (cv_n < 20)
  stop(sprintf(paste("figureS1 C: Complex V resolved to only %d genes -- the symbol",
                     "reconciler has regressed (expect 24)."), cv_n))

pC <- heat_panel(
  stC,
  "C   OXPHOS moves together",
  "All five complexes rise with MYC and attenuate\ntogether. Complex II is the small underpowered\nexception.")

# =============================================================================
# ASSEMBLY
# =============================================================================
bottom <- patchwork::wrap_plots(pB, pC, nrow = 1, widths = c(1.5, 1))
p <- patchwork::wrap_plots(pA, bottom, ncol = 1, heights = c(1, 1.1)) +
  patchwork::plot_annotation(
    caption = paste(
      "Raw (unshrunken) DESeq2 log2FC throughout; n=6 per group. Set membership resolved with functions/reconcile_gene_symbols.R, which recovers genes renamed since",
      "MitoCarta 3.0 -- a plain symbol match drops 17 of Complex V's 24 subunits and ~12% of nuclear OXPHOS. OXPHOS is shown as its nuclear and mtDNA-encoded halves; set sizes in brackets are the genes resolved in the DE universe.",
      "The two GENOTYPE columns are the batch-clean axis. BATCH = TIMEPOINT, so the two temporal columns are confounded and are DESCRIBED, not claimed.",
      sep = "\n"),
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(size = 5.7, hjust = 0, colour = "grey30",
                                           lineheight = 1.15)))

if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "figureS1_compartment_detail.pdf"),
             width = fig_w[["double"]], height = 165)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  ## A: the compositional answer, as numbers -- % up per arm x contrast
  stA[, c("label", "contrast", "med_lfc", "pct_up", "n_up", "n_dn")] |> print(n = 20)

  ## B: the WT column, sorted -- what actually rises while OXPHOS falls
  b <- stB[stB$contrast == "WT 6->12W", c("label", "med_lfc", "pct_up")]
  b[order(b$med_lfc), ] |> print(n = 16)

  ## C: the reconciler guard -- Complex V must resolve to ~24 genes, not 7
  unique(stC[, c("arm", "n_set")]) |> print()

  print(pA); print(pB); print(pC)
  print(p)
}
