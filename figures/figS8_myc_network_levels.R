# =============================================================================
# figS8_myc_network_levels.R -- is the Myc DOSE constant? The driver, its proximal
# network, and its output, measured the same way.
# -----------------------------------------------------------------------------
# The project has assumed throughout that Myc dose is stable across the 6W-12W
# window ("constant driver, moving substrate") on the strength of mRNA + blot +
# literature. Everything in sec-attenuation rests on it: if the driver fell, the
# attenuation would need no further explanation. This panel puts the assumption on
# the same footing as the results it supports.
#
#   A  Myc itself, normalised counts, four groups.
#   B  The proximal MYC/MAX/MXD network: MAX (the obligate partner), the competing
#      repressors MXD1/MXD3/MXD4/MXI1/MNT/MGA, and the MLX arm MLX/MLXIP/MLXIPL.
#      (MXI1 is MXD2 -- the family is complete.)
#   C  The three quantities that would register a falling dose -- Myc, Myc:MAX, and
#      Myc relative to the summed repressors -- against Myc's OUTPUT (Hallmark MYC
#      targets V2). The genotype gap in the DRIVER widens 6W->12W while the gap in
#      the OUTPUT narrows: the attenuation is not a dose effect.
#
# Ratios are read only BETWEEN groups. Normalised counts are not length-normalised,
# so the absolute value of Myc/MAX is meaningless, but transcript length is constant
# across samples and cancels in a group-to-group comparison of the same ratio.
#
# CAVEAT: these are reads on the mouse Myc locus. Total Myc mRNA in Myc+ is
# transgene + endogenous and the two CANNOT be separated here (see script 27 for the
# endogenous-vs-transgene decomposition at the programme level). "Dose is constant"
# is therefore a statement about total Myc message, which is what the substrate sees.
#
# Reads (read-only, no re-run):
#   results/dds_int_run.rds            -- normalised counts + colData
#   results/combined_df_annotated.rds  -- mgi_symbol <-> ensembl
#   data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt  (target composite)
#   functions/reconcile_gene_symbols.R
# =============================================================================

source(here::here("figures", "theme_myc.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))
if (!requireNamespace("patchwork", quietly = TRUE)) stop("figS8 needs patchwork")
if (!requireNamespace("DESeq2", quietly = TRUE))    stop("figS8 needs DESeq2 (normalised counts)")

out_dir <- here::here("outputs", "figures")

dds <- readRDS(here::here("results", "dds_int_run.rds"))
cdf <- readRDS(here::here("results", "combined_df_annotated.rds"))
gmt <- fgsea::gmtPathways(
  here::here("data", "genesets_from_library", "mammary_mito_myc_metab_v1_mouse.gmt"))

nc <- DESeq2::counts(dds, normalized = TRUE)
sm <- as.data.frame(SummarizedExperiment::colData(dds))
sm$group     <- factor(sm$group, levels = names(group_labels))
sm$timepoint <- stats::relevel(as.factor(sm$timepoint),  "6W")
sm$myc_status<- stats::relevel(as.factor(sm$myc_status), "neg")

sym2ens <- cdf[!is.na(cdf$mgi_symbol) & !duplicated(cdf$mgi_symbol), c("mgi_symbol", "gene")]
ens_of  <- function(s) sym2ens$gene[match(s, sym2ens$mgi_symbol)]
expr_of <- function(s) { e <- ens_of(s); if (is.na(e) || !e %in% rownames(nc)) NULL else nc[e, ] }

# --- the roster, grouped by role -------------------------------------------------
net <- tibble::tribble(
  ~gene,    ~role,
  "Max",    "obligate partner",
  "Mxd1",   "MXD/MNT repressors", "Mxd3", "MXD/MNT repressors", "Mxd4", "MXD/MNT repressors",
  "Mxi1",   "MXD/MNT repressors", "Mnt",  "MXD/MNT repressors", "Mga",  "MXD/MNT repressors",
  "Mlx",    "MLX arm", "Mlxip", "MLX arm", "Mlxipl", "MLX arm")
missing <- net$gene[vapply(net$gene, function(g) is.null(expr_of(g)), logical(1))]
if (length(missing)) message("figS8: not in the count matrix: ", paste(missing, collapse = ", "))
net <- net[!net$gene %in% missing, ]
stopifnot(!is.null(expr_of("Myc")))

# --- simple-effect p-values on log2 expression (fig01 / script 32 method) --------
pval_for <- function(y, kind, key) {
  d <- data.frame(y = y, timepoint = sm$timepoint, myc_status = sm$myc_status)
  m <- if (kind == "geno") stats::lm(y ~ myc_status, d[d$timepoint == key, ])
       else                stats::lm(y ~ timepoint,  d[d$myc_status == key, ])
  summary(m)$coefficients[2, "Pr(>|t|)"]
}
fmt_p <- function(p) if (p < 0.001) "<0.001" else formatC(p, format = "g", digits = 2)

# comparison geometry, identical to fig01: 6W_neg=1, 6W_pos=2, 12W_neg=3, 12W_pos=4
comp <- data.frame(
  x1 = c(1, 3, 1, 2), x2 = c(2, 4, 3, 4),
  kind = c("geno", "geno", "time", "time"),
  key  = c("6W", "12W", "neg", "pos"),
  level = c(1, 1, 2, 3), stringsAsFactors = FALSE)

# `hr` = the headroom the facet must reserve. A geom_text label has no data extent,
# so free_y scales to the top BRACKET and the top LABEL is then clipped by the strip;
# a geom_blank at `hr` reserves the space explicitly (the fig01 idiom).
brackets_for <- function(y, facet_val, facet_name) {
  b <- comp
  b$p    <- vapply(seq_len(nrow(b)), function(i) pval_for(y, b$kind[i], b$key[i]), numeric(1))
  b$lab  <- vapply(b$p, fmt_p, character(1))
  top    <- max(y); rng <- diff(range(y))
  b$y    <- top + rng * 0.11 * b$level
  b$tick <- rng * 0.02
  b$xmid <- (b$x1 + b$x2) / 2
  b$col  <- ifelse(b$p < 0.05, "sig", "ns")
  b$hr   <- max(b$y) + rng * 0.14          # bracket line + room for its label above it
  b[[facet_name]] <- facet_val
  b
}

# --- the shared panel builder ----------------------------------------------------
group_panel <- function(df, ylab, ttl, sub = NULL, facet = NULL, ncol = 5,
                        brk = NULL, base = 8) {
  p <- ggplot2::ggplot(df, ggplot2::aes(group, y, colour = group, fill = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.28,
                          colour = "grey35", linewidth = 0.3) +
    ggplot2::geom_point(size = 1.1, alpha = 0.95,
                        position = ggplot2::position_jitter(width = 0.16, height = 0, seed = 1))
  if (!is.null(brk)) {
    p <- p +
      ggplot2::geom_segment(data = brk, inherit.aes = FALSE,
                            ggplot2::aes(x = x1, xend = x2, y = y, yend = y, colour = col),
                            linewidth = 0.25) +
      ggplot2::geom_segment(data = brk, inherit.aes = FALSE,
                            ggplot2::aes(x = x1, xend = x1, y = y, yend = y - tick, colour = col),
                            linewidth = 0.25) +
      ggplot2::geom_segment(data = brk, inherit.aes = FALSE,
                            ggplot2::aes(x = x2, xend = x2, y = y, yend = y - tick, colour = col),
                            linewidth = 0.25) +
      ggplot2::geom_text(data = brk, inherit.aes = FALSE,
                         ggplot2::aes(x = xmid, y = y, label = lab, colour = col),
                         vjust = -0.3, size = 1.8) +
      ggplot2::geom_blank(data = unique(brk[, c(setdiff(names(brk), names(comp)), "hr")]),
                          inherit.aes = FALSE,
                          ggplot2::aes(x = 1, y = hr))
  }
  if (!is.null(facet)) p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)),
                                                    ncol = ncol, scales = "free_y")
  p +
    ggplot2::scale_colour_manual(values = c(group_cols, sig = "#E41A1C", ns = "grey45"),
                                 breaks = names(group_cols), labels = group_labels) +
    ggplot2::scale_fill_manual(values = group_cols, guide = "none") +
    ggplot2::labs(x = NULL, y = ylab, colour = NULL, title = ttl, subtitle = sub) +
    theme_myc(base_size = base) +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_blank(),
      axis.ticks.x  = ggplot2::element_blank(),
      axis.line.x   = ggplot2::element_blank(),
      strip.text    = ggplot2::element_text(size = base - 1, face = "bold"),
      plot.title    = ggplot2::element_text(face = "bold", size = base + 1),
      plot.subtitle = ggplot2::element_text(size = base - 1.5, colour = "grey25",
                                            lineheight = 1.15))
}

# =============================================================================
# PANEL A -- Myc itself
# =============================================================================
myc <- expr_of("Myc")
dA  <- data.frame(group = sm$group, y = myc)
bA  <- brackets_for(myc, NA, ".dummy"); bA$.dummy <- NULL
gapA <- c(g6  = mean(log2(myc[sm$group == "6W_pos"]))  - mean(log2(myc[sm$group == "6W_neg"])),
          g12 = mean(log2(myc[sm$group == "12W_pos"])) - mean(log2(myc[sm$group == "12W_neg"])))

pA <- group_panel(
  dA, "Myc, normalised counts", "A   Myc",
  sprintf("genotype gap %+.2f log2 at 6W -> %+.2f at 12W;\nno decline within either genotype", gapA[1], gapA[2]),
  brk = bA, base = 8)

# =============================================================================
# PANEL B -- the proximal network
# =============================================================================
dB <- do.call(rbind, lapply(seq_len(nrow(net)), function(i) {
  data.frame(group = sm$group, y = expr_of(net$gene[i]),
             gene = net$gene[i], role = net$role[i])
}))
dB$gene <- factor(dB$gene, levels = net$gene)
bB <- do.call(rbind, lapply(net$gene, function(g) brackets_for(expr_of(g), g, "gene")))
bB$gene <- factor(bB$gene, levels = net$gene)

pB <- group_panel(
  dB, "normalised counts", "B   The proximal MYC/MAX/MXD network",
  "MAX (partner) | MXD1/3/4, MXI1 (=MXD2), MNT, MGA (competing repressors) | MLX, MLXIP, MLXIPL",
  facet = "gene", ncol = 5, brk = bB, base = 8)

# =============================================================================
# PANEL C -- driver vs output: does anything that would register a falling dose fall?
# =============================================================================
rep_syms <- c("Mxd1", "Mxd3", "Mxd4", "Mxi1", "Mnt", "Mga")
rep_syms <- rep_syms[rep_syms %in% net$gene]
rep_sum  <- Reduce(`+`, lapply(rep_syms, expr_of))
tgt_ens  <- recon_to_ensembl(gmt[["MYC_HALLMARK_MYC_TARGETS_V2"]], rownames(nc))
tgt_ens  <- tgt_ens[!is.na(tgt_ens)]
tgt_z    <- t(scale(t(log2(nc[tgt_ens, , drop = FALSE] + 1))))
tgt_z    <- tgt_z[is.finite(rowSums(tgt_z)), , drop = FALSE]

measures <- list(
  `Myc`                    = log2(myc),
  `Myc : MAX`              = log2(myc / expr_of("Max")),
  `Myc : repressor sum`    = log2(myc / rep_sum),
  `MYC targets (V2)`       = as.numeric(colMeans(tgt_z)))
kind_of <- c(`Myc` = "driver", `Myc : MAX` = "driver", `Myc : repressor sum` = "driver",
             `MYC targets (V2)` = "output")

dC <- do.call(rbind, lapply(names(measures), function(k) {
  y <- measures[[k]]
  data.frame(measure = k, kind = kind_of[[k]],
             timepoint = c("6W", "12W"),
             gap = c(mean(y[sm$group == "6W_pos"])  - mean(y[sm$group == "6W_neg"]),
                     mean(y[sm$group == "12W_pos"]) - mean(y[sm$group == "12W_neg"])))
}))
dC$measure   <- factor(dC$measure, levels = names(measures))
dC$timepoint <- factor(dC$timepoint, levels = c("6W", "12W"))
wideC <- data.frame(
  measure = factor(names(measures), levels = names(measures)),
  kind    = unname(kind_of[names(measures)]),
  g6      = dC$gap[match(names(measures), dC$measure[dC$timepoint == "6W"])],
  g12     = dC$gap[dC$timepoint == "12W"][match(names(measures),
                                                dC$measure[dC$timepoint == "12W"])])
wideC$g6    <- vapply(names(measures), function(k) dC$gap[dC$measure == k & dC$timepoint == "6W"], numeric(1))
wideC$g12   <- vapply(names(measures), function(k) dC$gap[dC$measure == k & dC$timepoint == "12W"], numeric(1))
wideC$delta <- wideC$g12 - wideC$g6

kind_col <- c(driver = "#D73027", output = "#4575B4")
pC <- ggplot2::ggplot(wideC, ggplot2::aes(y = measure)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  ggplot2::geom_segment(ggplot2::aes(x = g6, xend = g12, yend = measure, colour = kind),
                        linewidth = 0.5,
                        arrow = ggplot2::arrow(length = ggplot2::unit(1.6, "mm"), type = "closed")) +
  ggplot2::geom_point(ggplot2::aes(x = g6, colour = kind), shape = 21, fill = "white",
                      size = 2, stroke = 0.5) +
  ggplot2::geom_text(ggplot2::aes(x = pmax(g6, g12), label = sprintf("%+.2f", delta),
                                  colour = kind),
                     hjust = -0.35, size = 2.2, fontface = "bold") +
  ggplot2::scale_colour_manual(values = kind_col, name = NULL,
                               labels = c(driver = "driver", output = "output")) +
  ggplot2::scale_y_discrete(limits = rev(levels(dC$measure))) +
  ggplot2::coord_cartesian(xlim = c(0, max(c(wideC$g6, wideC$g12)) * 1.25)) +
  ggplot2::labs(
    x = "genotype gap (Myc+ - WT), log2", y = NULL,
    title = "C   The driver does not fall while its output does",
    subtitle = "open = 6W, arrow to 12W; label = the change in the gap.\nRatios read between groups only (transcript length cancels).") +
  theme_myc(base_size = 8) +
  ggplot2::theme(legend.position = "bottom",
                 legend.key.size = ggplot2::unit(3, "mm"),
                 legend.text     = ggplot2::element_text(size = 6.5),
                 axis.text.y     = ggplot2::element_text(size = 7),
                 plot.title      = ggplot2::element_text(face = "bold", size = 9),
                 plot.subtitle   = ggplot2::element_text(size = 6.4, colour = "grey25",
                                                         lineheight = 1.15))

# =============================================================================
# ASSEMBLY
# =============================================================================
top <- patchwork::wrap_plots(pA, pC, nrow = 1, widths = c(1, 1.5)) +
  patchwork::plot_layout(guides = "keep")

p <- patchwork::wrap_plots(top, pB, ncol = 1, heights = c(1, 1.25)) +
  patchwork::plot_annotation(
    caption = paste(
      "DESeq2 median-of-ratios normalised counts, n=6/group; points = mice, box = median/IQR. Brackets: genotype (WT vs Myc+) within a timepoint and 6W-vs-12W within a genotype,",
      "simple-effect lm on log2 expression; red = p<0.05. Reads map to the mouse Myc locus: total Myc message in Myc+ is transgene + endogenous and the two are not separable here.",
      "The 6W-vs-12W axis is batch-confounded (batch = timepoint), so a time bracket is described, not claimed -- but the claim being tested is the ABSENCE of a temporal decline.",
      "p-values are UNADJUSTED across 11 genes x 4 comparisons; nothing in panel B clears BH<0.05 over that family, and the genome-wide DESeq2 interaction padj is >0.05 for every network gene.",
      sep = "\n"),
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(size = 5.8, hjust = 0, colour = "grey30",
                                           lineheight = 1.15)))

if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "figS8_myc_network_levels.pdf"),
             width = fig_w[["double"]], height = 150)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  ## group means of every gene in the roster
  round(t(vapply(c("Myc", net$gene), function(g)
    tapply(expr_of(g), sm$group, mean), numeric(4))), 0)

  ## the driver-vs-output table (panel C)
  print(wideC, digits = 3)

  ## every simple-effect p for Myc
  print(bA[, c("kind", "key", "p", "lab")], digits = 3)

  print(pA); print(pB); print(pC); print(p)
}
