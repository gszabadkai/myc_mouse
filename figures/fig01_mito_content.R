# =============================================================================
# fig01_mito_content.R -- "Myc raises mitochondrial content"
# -----------------------------------------------------------------------------
# The single strongest new result of Block B, taken to publication quality.
# Genotype (Myc+ vs WT) is the CLEAN axis (depth-balanced, litter-controlled);
# the 6W->12W time axis is cohort/batch-confounded and is shown only as context.
#
# Reads (read-only, already on disk):
#   results/mito_content_proxies.rds   (script 32) -- $shares (per-sample % of the
#       non-mtDNA transcriptome, n=6/group), $share_stats (genotype main effect
#       geno_beta/geno_p per compartment arm).
#   results/mtdna_axis_and_coupling_null.rds (script 33) -- for the adjusted
#       headline (+21%, prep-stress + contamination adjusted) quoted in the caption.
#
# Distils script 32's 14-facet draft (outputs/mito_content_proxies/A) to the four
# claim-bearing arms so "nuclear up, mtDNA flat" reads in one panel.
#
# Design notes (skills: distribution-plots, ggplot2-fundamentals, figure-export):
#   n=6/group -> show EVERY point (quasirandom) over a light box; never bar-of-mean.
#   Compartment arm + its genotype effect (%/p) live in the facet strip.
# =============================================================================

source(here::here("figures", "theme_myc.R"))

out_dir <- here::here("outputs", "figures")

content <- readRDS(here::here("results", "mito_content_proxies.rds"))

# --- the four claim-bearing arms, ordered nuclear-up first, mtDNA (flat) last --
arms <- c("MASS_MARKERS_NOCHAP",
          "MITOCARTA_NUCLEAR_ENCODED",
          "MITOCARTA_MITOCHONDRIAL_RIBOSOME",
          "MITOCARTA_MTDNA_ENCODED")
arm_name <- c(MASS_MARKERS_NOCHAP              = "Mass markers",
              MITOCARTA_NUCLEAR_ENCODED        = "Nuclear MitoCarta",
              MITOCARTA_MITOCHONDRIAL_RIBOSOME = "Mitoribosome",
              MITOCARTA_MTDNA_ENCODED          = "mtDNA-encoded")

# --- genotype effect per arm (from script 32), for the facet-strip annotation --
ss <- content$share_stats[content$share_stats$denominator == "share_nomt" &
                            content$share_stats$panel %in% arms, ]
ss$pct  <- 100 * (2^ss$geno_beta - 1)
strip_lab <- vapply(arms, function(a) {
  r  <- ss[ss$panel == a, ]
  ns <- if (r$geno_p >= 0.05) " (ns)" else ""
  sprintf("%s\nMyc %+.0f%%, p=%s%s",
          arm_name[[a]], r$pct, formatC(r$geno_p, format = "g", digits = 2), ns)
}, character(1))
names(strip_lab) <- arms

# --- per-sample data for the four arms ---------------------------------------
df <- content$shares
df <- df[df$panel %in% arms, ]
df$panel      <- factor(df$panel, levels = arms, labels = strip_lab[arms])
df$group      <- factor(df$group, levels = names(group_labels), labels = group_labels)
df$myc_status <- factor(df$myc_status, levels = c("neg", "pos"))

n_per <- min(table(content$shares$group[content$shares$panel == arms[1]]))

# --- point geom: quasirandom if available, else jitter (n=6 -> show all) -------
pts <- if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
  ggbeeswarm::geom_quasirandom(width = 0.22, size = 1.3, alpha = 0.9)
} else {
  ggplot2::geom_jitter(width = 0.18, height = 0, size = 1.3, alpha = 0.9)
}

p <- ggplot2::ggplot(df, ggplot2::aes(group, share_nomt,
                                      colour = myc_status, fill = myc_status)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.16,
                        linewidth = 0.35) +
  pts +
  ggplot2::facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  ggplot2::scale_colour_manual(values = geno_cols) +
  ggplot2::scale_fill_manual(values = geno_cols) +
  ggplot2::scale_x_discrete(labels = function(x) sub(" ", "\n", x)) +
  ggplot2::labs(
    x = NULL, y = "share of non-mtDNA transcriptome (%)",
    title = "Myc raises the nuclear mitochondrial arms; mtDNA-encoded output stays flat",
    subtitle = sprintf(
      "Genotype is the clean axis (depth-balanced); the 6W->12W time axis is cohort-confounded. n=%d/group.",
      n_per),
    caption = paste(
      "Share of the non-mtDNA transcriptome (raw counts); % and p are the genotype main effect (script 32).",
      "The compartment-wide nuclear increase is +21% (p=0.0006) after adjusting for prep-stress + contamination (script 33).",
      "Transcript share, not per-cell content: with global RNA amplification this is a LOWER BOUND; mtDNA qPCR / blot settles it.",
      sep = "\n")) +
  theme_myc(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7))

# Guard: when this script is sourced only to obtain `p` (e.g. from the Quarto doc),
# set options(myc.fig.nosave = TRUE) to display without re-exporting the canonical PDF.
if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "fig01_mito_content.pdf"),
             width = fig_w[["double"]], height = 90)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(ss[, c("panel", "geno_beta", "pct", "geno_p")])
  print(strip_lab)
  print(n_per)
  print(p)
}
