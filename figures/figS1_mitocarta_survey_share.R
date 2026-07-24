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
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("figS1 needs patchwork to divide the top-level and metabolism blocks")
}

out_dir <- here::here("outputs", "figures")

content <- readRDS(here::here("results", "mito_content_proxies.rds"))

# --- the 16 survey groups (same identities + labels as figS1b) ----------------
# OXPHOS is shown as its two halves -- MITOCARTA_OXPHOS_NU (nuclear, 155, no mt-*
# genes) and MITOCARTA_OXPHOS_MT (the 13 mtDNA-encoded subunits, identical to
# MITOCARTA_MTDNA_ENCODED). So there is no separate "mtDNA reference" facet: the
# mtDNA-encoded compartment IS mtDNA OXPHOS.
arms <- c(
  "MITOCARTA_OXPHOS_NU",
  "MITOCARTA_OXPHOS_MT",
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
  "MITOCARTA_ELECTRON_CARRIERS")
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

# --- p-value brackets, SIGNIFICANT ONLY (author request) ----------------------
# Same four comparisons and the SAME log2-share simple-effect method as fig01
# (script 32's own approach, so the numbers reconcile). Only brackets with
# p<0.05 are drawn -- genotype (clean) in red, 6W-vs-12W time (batch-confounded)
# in grey. Non-significant comparisons are simply omitted to keep 16 facets legible.
comp <- data.frame(
  x1    = c(1, 1, 3, 2),
  x2    = c(2, 3, 4, 4),
  kind  = c("geno", "time", "geno", "time"),
  key   = c("6W", "neg", "12W", "pos"),   # geno -> timepoint subset; time -> myc subset
  clean = c(TRUE, FALSE, TRUE, FALSE),     # x: 6W_neg=1, 6W_pos=2, 12W_neg=3, 12W_pos=4
  stringsAsFactors = FALSE)

pval_for <- function(d, kind, key) {
  d$yv <- log2(d$share_nomt)
  m <- if (kind == "geno") {
    stats::lm(yv ~ myc_status, data = d[d$timepoint == key, ])
  } else {
    stats::lm(yv ~ timepoint, data = d[d$myc_status == key, ])
  }
  summary(m)$coefficients[2, "Pr(>|t|)"]
}
fmt_p <- function(p) if (p < 0.001) "<0.001" else formatC(p, format = "g", digits = 2)

build_brackets <- function(arm_ids) {
  do.call(rbind, lapply(arm_ids, function(a) {
    d <- content$shares[content$shares$panel == a, ]
    b <- comp
    b$p <- vapply(seq_len(nrow(b)), function(i) pval_for(d, b$kind[i], b$key[i]),
                  numeric(1))
    b <- b[b$p < 0.05, , drop = FALSE]              # SIGNIFICANT ONLY
    if (nrow(b) == 0) return(NULL)
    b <- b[order(!b$clean, b$x1), ]                 # genotype first, then by width
    facmax <- max(d$share_nomt)
    data.frame(
      panel = unname(arm_name[[a]]),
      x1 = b$x1, x2 = b$x2, xmid = (b$x1 + b$x2) / 2,
      y  = facmax * (1 + 0.11 * seq_len(nrow(b))),
      tick = facmax * 0.02,
      lab = vapply(b$p, fmt_p, character(1)),
      col = ifelse(b$kind == "geno", "geno_sig", "time_sig"),
      stringsAsFactors = FALSE)
  }))
}

# --- two blocks with a divider between them (author request 2026-07-24) --------
# ncol = 3 fills the 9-child metabolism block exactly. OXPHOS leads the top block
# as its two halves (nuclear / mtDNA); there is no separate mtDNA reference facet.
blk_top <- arms[1:7]    # Nuclear OXPHOS, mtDNA OXPHOS + 5 other top-level categories
blk_met <- arms[8:16]   # 9 Metabolism level-2 children

make_block <- function(arm_ids, subtitle, ylab = NULL) {
  labs <- unname(arm_name[arm_ids])
  d <- df[df$panel %in% labs, ]
  d$panel <- factor(as.character(d$panel), levels = labs)
  brk <- build_brackets(arm_ids)
  # per-facet headroom so the top bracket + label fit (or a touch above the data)
  hr <- do.call(rbind, lapply(arm_ids, function(a) {
    lab   <- unname(arm_name[[a]])
    dmax  <- max(content$shares$share_nomt[content$shares$panel == a])
    ytop  <- if (!is.null(brk) && any(brk$panel == lab)) {
      max(brk$y[brk$panel == lab]) * 1.07
    } else dmax * 1.03
    data.frame(panel = lab, group = names(group_labels)[1], y = ytop,
               stringsAsFactors = FALSE)
  }))
  hr$panel  <- factor(hr$panel, levels = labs)
  hr$group  <- factor(hr$group, levels = names(group_labels))
  brk_layers <- NULL
  if (!is.null(brk)) {
    brk$panel <- factor(brk$panel, levels = labs)
    brk_layers <- list(
      ggplot2::geom_segment(data = brk, inherit.aes = FALSE,
        ggplot2::aes(x = x1, xend = x2, y = y, yend = y, colour = col), linewidth = 0.25),
      ggplot2::geom_segment(data = brk, inherit.aes = FALSE,
        ggplot2::aes(x = x1, xend = x1, y = y, yend = y - tick, colour = col), linewidth = 0.25),
      ggplot2::geom_segment(data = brk, inherit.aes = FALSE,
        ggplot2::aes(x = x2, xend = x2, y = y, yend = y - tick, colour = col), linewidth = 0.25),
      ggplot2::geom_text(data = brk, inherit.aes = FALSE,
        ggplot2::aes(x = xmid, y = y, label = lab, colour = col), vjust = -0.2, size = 1.8))
  }
  ggplot2::ggplot(d, ggplot2::aes(group, share_nomt, colour = group, fill = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.28,
                          colour = "grey35", linewidth = 0.3) +
    pts_layer(d) +
    brk_layers +
    ggplot2::geom_blank(data = hr, ggplot2::aes(x = group, y = y), inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ panel, ncol = 3, scales = "free_y") +
    ggplot2::scale_colour_manual(
      values = c(group_cols, geno_sig = "#E41A1C", time_sig = "grey45"),
      breaks = names(group_cols), labels = group_labels) +
    ggplot2::scale_fill_manual(values = group_cols, guide = "none") +
    ggplot2::labs(x = NULL, y = ylab, colour = NULL, subtitle = subtitle) +
    theme_myc(base_size = 8) +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_blank(),
      axis.ticks.x    = ggplot2::element_blank(),
      axis.line.x     = ggplot2::element_blank(),
      strip.text      = ggplot2::element_text(size = 6.8),
      plot.subtitle   = ggplot2::element_text(face = "bold", size = 8.5, colour = "grey20"),
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(3.5, "mm")) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      nrow = 1, override.aes = list(size = 2.4, alpha = 1, shape = 16)))
}

p_top <- make_block(blk_top, "Top-level MitoPathway categories (OXPHOS split nuclear / mtDNA)",
                    ylab = "share of the nuclear transcriptome (%)")
p_met <- make_block(blk_met, "Metabolism (level-2 children)")

p <- patchwork::wrap_plots(p_top, p_met, ncol = 1, heights = c(3, 3)) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "Mitochondrial compartment shares across all main MitoCarta groups",
    caption = paste(
      sprintf("Points, n=%d/group; box = median/IQR. %% of the non-mtDNA transcriptome (raw counts); free y per group.", n_per),
      "Brackets: genotype (WT vs Myc+, red) and 6W-vs-12W (grey) p-values, drawn ONLY where p<0.05 (log2-share simple-effect models, script 32).",
      "Metabolism divided from the other top-level categories, split into its 9 depth-2 children. Genotype is the clean axis; the time axis is batch-confounded (batch=timepoint).",
      sep = "\n"),
    theme = ggplot2::theme(
      plot.title   = ggplot2::element_text(face = "bold", size = 11),
      plot.caption = ggplot2::element_text(size = 6.3, hjust = 0, colour = "grey30",
                                           lineheight = 1.1))) &
  ggplot2::theme(legend.position = "bottom")

# Guard: sourced only to obtain `p` (e.g. Quarto) when myc.fig.nosave = TRUE.
if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "figS1_mitocarta_survey_share.pdf"),
             width = fig_w[["double"]], height = 170)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(n_per)
  print(table(df$panel))
  print(p)
}
