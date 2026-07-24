# =============================================================================
# figS1_mitocarta_survey_share.R -- compartment-share survey, ALL main MitoCarta groups
# -----------------------------------------------------------------------------
# Supplementary companion to figS1b (the LFC survey). The share-of-transcriptome
# view of the same 16 MitoCarta MitoPathway groups across the four groups, so the
# author's question -- "OXPHOS falls in the WT timeline but the total is flat, so
# what rises?" -- can be read in absolute compartment terms as well as in LFC.
#
# 16 groups = 7 top-level MitoPathway categories with METABOLISM split into its 9
# depth-2 children (hierarchy from Sheet 4), plus mtDNA-encoded as the flat
# reference. All sets are already mt-* free, so MITOCARTA_OXPHOS_NU is OXPHOS and
# share_nomt (denominator drops only the 13 mt-* genes) is the consistent scale.
#
# Reads (read-only; NEEDS script 32 re-run so content$shares carries the 12 new
#   survey panels -- see plan / script 32 PART 2b roster additions):
#   results/mito_content_proxies.rds -- $shares (per-sample % of the non-mtDNA
#       transcriptome, n=6/group).
#
# Survey, not a claim panel: no per-contrast brackets. Genotype (Myc+ vs WT within
# a timepoint) is the clean axis; the 6W->12W time axis is cohort/batch-confounded
# (batch = timepoint) and read as context. Per-group genotype/time stats already
# sit in content$share_stats for the Quarto text.
# =============================================================================

source(here::here("figures", "theme_myc.R"))

out_dir <- here::here("outputs", "figures")

content <- readRDS(here::here("results", "mito_content_proxies.rds"))

# --- the 16 survey groups (same identities + labels as figS1b) ----------------
arms <- c(
  "MITOCARTA_OXPHOS_NU",
  "MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA",
  "MITOCARTA_PROTEIN_IMPORT_SORTING_AND_HOMEOSTASIS",
  "MITOCARTA_MITOCHONDRIAL_DYNAMICS_AND_SURVEILLANCE",
  "MITOCARTA_SIGNALING",
  "MITOCARTA_SMALL_MOLECULE_TRANSPORT",
  "MITOCARTA_AMINO_ACID_METABOLISM",
  "MITOCARTA_CARBOHYDRATE_METABOLISM",
  "MITOCARTA_LIPID_METABOLISM",
  "MITOCARTA_NUCLEOTIDE_METABOLISM",
  "MITOCARTA_VITAMIN_METABOLISM",
  "MITOCARTA_METALS_AND_COFACTORS",
  "MITOCARTA_DETOXIFICATION",
  "MITOCARTA_SULFUR_METABOLISM",
  "MITOCARTA_ELECTRON_CARRIERS",
  "MITOCARTA_MTDNA_ENCODED")
arm_name <- c(
  MITOCARTA_OXPHOS_NU                              = "OXPHOS",
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
  MITOCARTA_ELECTRON_CARRIERS                      = "Electron carriers",
  MITOCARTA_MTDNA_ENCODED                          = "mtDNA-encoded")

missing <- setdiff(arms, unique(content$shares$panel))
if (length(missing) > 0) {
  stop("figS1: content$shares is missing ", length(missing), " survey panel(s): ",
       paste(missing, collapse = ", "),
       ".\n  -> re-run scripts/32_mito_content_proxies.R (the roster gained these rows).")
}

df <- content$shares
df <- df[df$panel %in% arms, ]
df$panel <- factor(df$panel, levels = arms, labels = arm_name[arms])
df$group <- factor(df$group, levels = names(group_labels))

n_per <- min(table(df$group[df$panel == arm_name[[arms[1]]]]))

# --- point geom: quasirandom if available, else jitter (n=6 -> show all) -------
pts_layer <- function(dat) {
  if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
    ggbeeswarm::geom_quasirandom(data = dat, width = 0.22, size = 0.9, alpha = 0.9)
  } else {
    ggplot2::geom_jitter(data = dat, width = 0.16, height = 0, size = 0.9, alpha = 0.9)
  }
}

p <- ggplot2::ggplot(df, ggplot2::aes(group, share_nomt, colour = group, fill = group)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.28,
                        colour = "grey35", linewidth = 0.3) +
  pts_layer(df) +
  ggplot2::facet_wrap(~ panel, ncol = 4, scales = "free_y") +
  ggplot2::scale_colour_manual(values = group_cols,
                               breaks = names(group_cols), labels = group_labels) +
  ggplot2::scale_fill_manual(values = group_cols, guide = "none") +
  ggplot2::labs(
    x = NULL, y = "share of the nuclear transcriptome (%)", colour = NULL,
    title = "Mitochondrial compartment shares across all main MitoCarta groups",
    caption = paste(
      sprintf("Points, n=%d/group; box = median/IQR. %% of the non-mtDNA transcriptome (raw counts); free y per group.", n_per),
      "16 MitoCarta MitoPathway groups (Metabolism split into its 9 depth-2 children); mtDNA-encoded = flat reference.",
      "Genotype (Myc+ vs WT within a timepoint) is the clean axis; the 6W-vs-12W time axis is cohort/batch-confounded (batch=timepoint).",
      sep = "\n")) +
  theme_myc(base_size = 8) +
  ggplot2::theme(
    axis.text.x     = ggplot2::element_blank(),
    axis.ticks.x    = ggplot2::element_blank(),
    axis.line.x     = ggplot2::element_blank(),
    strip.text      = ggplot2::element_text(size = 6.8),
    legend.position = "bottom",
    legend.key.size = ggplot2::unit(3.5, "mm")) +
  ggplot2::guides(colour = ggplot2::guide_legend(
    nrow = 1, override.aes = list(size = 2.4, alpha = 1, shape = 16)))

# Guard: sourced only to obtain `p` (e.g. Quarto) when myc.fig.nosave = TRUE.
if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "figS1_mitocarta_survey_share.pdf"),
             width = fig_w[["double"]], height = 150)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(n_per)
  print(table(df$panel))
  print(p)
}
