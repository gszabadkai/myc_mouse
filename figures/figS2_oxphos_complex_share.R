# =============================================================================
# figS2_oxphos_complex_share.R -- compartment share of the five OXPHOS complexes
# -----------------------------------------------------------------------------
# The granular share companion to figS2b (LFC) and to figS1 (all groups): the
# transcriptome share held by each respiratory complex (CI-CV) across the four
# groups, so the nuclear-OXPHOS content story can be read complex by complex.
#
# Groups = MITOCARTA_COMPLEX_I..V, all NUCLEAR-encoded (mt-* removed): the 13
# mtDNA subunits are the mtDNA OXPHOS facet of figS1, not here. CII is fully nuclear.
#
# SELF-CONTAINED (does NOT route through script 32's roster): the shares are
# computed here, on RAW counts, exactly as script 32 does (share_nomt = 100 *
# panel counts / non-mtDNA counts). Two reasons: (1) it keeps script 32 -- and
# the already-signed-off fig01/figS1 -- untouched, needing no re-run; (2) it lets
# the complex membership use MitoCarta's OWN Ensembl IDs rather than a current-
# symbol match. That MAPPING matters: MitoCarta carries the OLD ATP-synthase names
# (Atp5a1/b/c1/...), the count matrix uses the CURRENT ones (Atp5f1a..), and a
# symbol match silently drops 17 of Complex V's 24 genes. Mapping by Ensembl
# recovers CV 24/24.
#
# Reads (read-only; NO analysis re-run needed):
#   results/count_matrix.rds         -- raw counts (Ensembl rows).
#   results/mito_content_proxies.rds -- sample -> group/timepoint/myc_status metadata
#                                        (script 32's; the buggy complex mapping is NOT used).
#   data/Mouse.MitoCarta3.0.xls      -- Sheet 2, Symbol + EnsemblGeneID (the map).
#   data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt -- complex memberships.
#
# Genotype (Myc+ vs WT within a timepoint) is the clean axis; the 6W->12W time
# axis is cohort/batch-confounded (batch = timepoint). Brackets: the same four
# comparisons and log2-share simple-effect method as figS1/fig01, drawn ONLY
# where p<0.05 (genotype red, 6W-vs-12W grey).
# =============================================================================

source(here::here("figures", "theme_myc.R"))
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("figS2 needs readxl to map MitoCarta symbols to Ensembl (see mapping note)")
}

out_dir <- here::here("outputs", "figures")

cts     <- readRDS(here::here("results", "count_matrix.rds"))
content <- readRDS(here::here("results", "mito_content_proxies.rds"))

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

# --- sample metadata (from script 32's shares; the complex mapping is NOT used) -
meta <- unique(content$shares[, c("sample", "group", "timepoint", "myc_status")])
meta <- meta[match(colnames(cts), meta$sample), ]
stopifnot(!anyNA(meta$group), identical(meta$sample, colnames(cts)))
meta$group     <- factor(meta$group, levels = names(group_labels))
meta$timepoint <- factor(meta$timepoint)
meta$myc_status<- factor(meta$myc_status)

# --- MitoCarta symbol -> Ensembl map (Sheet 2), then complexes + the 13 mt genes
mc <- as.data.frame(readxl::read_excel(
  here::here("data", "Mouse.MitoCarta3.0.xls"), sheet = 2))
mc <- mc[!is.na(mc$Symbol) & !is.na(mc$EnsemblGeneID) & mc$EnsemblGeneID != "", ]
mc_map <- do.call(rbind, lapply(seq_len(nrow(mc)), function(i)
  data.frame(sym = mc$Symbol[i],
             ens = trimws(strsplit(mc$EnsemblGeneID[i], "[|]")[[1]]),
             stringsAsFactors = FALSE)))
ens_of <- function(syms) intersect(unique(mc_map$ens[mc_map$sym %in% syms]), rownames(cts))

mt_ens  <- ens_of(gmt[["MITOCARTA_OXPHOS_MT"]])           # the 13 mtDNA subunits
stopifnot(length(mt_ens) >= 10)
arm_ens <- stats::setNames(lapply(arms, function(a) ens_of(gmt[[a]])), arms)

# --- shares on RAW counts (denominator excludes only the 13 mt-* genes) --------
den_nomt <- colSums(cts[setdiff(rownames(cts), mt_ens), , drop = FALSE])
share_of <- function(ens) 100 * colSums(cts[ens, , drop = FALSE]) / den_nomt

shares <- do.call(rbind, lapply(arms, function(a)
  data.frame(panel = a, sample = colnames(cts), group = meta$group,
             timepoint = meta$timepoint, myc_status = meta$myc_status,
             share_nomt = share_of(arm_ens[[a]]), stringsAsFactors = FALSE)))

df <- shares
df$panel <- factor(df$panel, levels = arms, labels = arm_name[arms])

n_per <- min(table(df$group[df$panel == arm_name[[arms[1]]]]))

# --- point geom: quasirandom if available, else jitter (n=6 -> show all) -------
pts_layer <- function(dat) {
  if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
    ggbeeswarm::geom_quasirandom(data = dat, width = 0.22, size = 1.0, alpha = 0.9)
  } else {
    ggplot2::geom_jitter(data = dat, width = 0.16, height = 0, size = 1.0, alpha = 0.9)
  }
}

# --- p-value brackets, SIGNIFICANT ONLY (same method as figS1/fig01) ----------
comp <- data.frame(
  x1    = c(1, 1, 3, 2),
  x2    = c(2, 3, 4, 4),
  kind  = c("geno", "time", "geno", "time"),
  key   = c("6W", "neg", "12W", "pos"),
  clean = c(TRUE, FALSE, TRUE, FALSE),
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

brk <- do.call(rbind, lapply(arms, function(a) {
  d <- shares[shares$panel == a, ]
  b <- comp
  b$p <- vapply(seq_len(nrow(b)), function(i) pval_for(d, b$kind[i], b$key[i]), numeric(1))
  b <- b[b$p < 0.05, , drop = FALSE]
  if (nrow(b) == 0) return(NULL)
  b <- b[order(!b$clean, b$x1), ]
  facmax <- max(d$share_nomt)
  data.frame(
    panel = factor(arm_name[[a]], levels = arm_name[arms]),
    x1 = b$x1, x2 = b$x2, xmid = (b$x1 + b$x2) / 2,
    y  = facmax * (1 + 0.11 * seq_len(nrow(b))),
    tick = facmax * 0.02,
    lab = vapply(b$p, fmt_p, character(1)),
    col = ifelse(b$kind == "geno", "geno_sig", "time_sig"),
    stringsAsFactors = FALSE)
}))

hr <- do.call(rbind, lapply(arms, function(a) {
  lab  <- unname(arm_name[[a]])
  dmax <- max(shares$share_nomt[shares$panel == a])
  ytop <- if (!is.null(brk) && any(brk$panel == lab)) max(brk$y[brk$panel == lab]) * 1.07
          else dmax * 1.03
  data.frame(panel = factor(lab, levels = arm_name[arms]),
             group = factor("6W_neg", levels = names(group_labels)), y = ytop)
}))

brk_layers <- if (!is.null(brk)) list(
  ggplot2::geom_segment(data = brk, inherit.aes = FALSE,
    ggplot2::aes(x = x1, xend = x2, y = y, yend = y, colour = col), linewidth = 0.25),
  ggplot2::geom_segment(data = brk, inherit.aes = FALSE,
    ggplot2::aes(x = x1, xend = x1, y = y, yend = y - tick, colour = col), linewidth = 0.25),
  ggplot2::geom_segment(data = brk, inherit.aes = FALSE,
    ggplot2::aes(x = x2, xend = x2, y = y, yend = y - tick, colour = col), linewidth = 0.25),
  ggplot2::geom_text(data = brk, inherit.aes = FALSE,
    ggplot2::aes(x = xmid, y = y, label = lab, colour = col), vjust = -0.2, size = 2.0)) else NULL

p <- ggplot2::ggplot(df, ggplot2::aes(group, share_nomt, colour = group, fill = group)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.28,
                        colour = "grey35", linewidth = 0.3) +
  pts_layer(df) +
  brk_layers +
  ggplot2::geom_blank(data = hr, ggplot2::aes(x = group, y = y), inherit.aes = FALSE) +
  ggplot2::facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  ggplot2::scale_colour_manual(
    values = c(group_cols, geno_sig = "#E41A1C", time_sig = "grey45"),
    breaks = names(group_cols), labels = group_labels) +
  ggplot2::scale_fill_manual(values = group_cols, guide = "none") +
  ggplot2::labs(
    x = NULL, y = "share of the nuclear transcriptome (%)", colour = NULL,
    title = "Nuclear-OXPHOS content, complex by complex",
    caption = paste(
      sprintf("Points, n=%d/group; box = median/IQR. %% of the non-mtDNA transcriptome (raw counts); free y per complex.", n_per),
      "Brackets: genotype (WT vs Myc+, red) and 6W-vs-12W (grey), drawn ONLY where p<0.05 (log2-share simple-effect models).",
      "MITOCARTA_COMPLEX_I..V mapped by MitoCarta Ensembl IDs (CI 57 / CII 8 / CIII 15 / CIV 45 / CV 24); nuclear-encoded. Genotype clean; time batch-confounded.",
      sep = "\n")) +
  theme_myc(base_size = 9) +
  ggplot2::theme(
    axis.text.x     = ggplot2::element_blank(),
    axis.ticks.x    = ggplot2::element_blank(),
    axis.line.x     = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.key.size = ggplot2::unit(3.5, "mm")) +
  ggplot2::guides(colour = ggplot2::guide_legend(
    nrow = 1, override.aes = list(size = 2.4, alpha = 1, shape = 16)))

# Guard: sourced only to obtain `p` (e.g. Quarto) when myc.fig.nosave = TRUE.
if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "figS2_oxphos_complex_share.pdf"),
             width = fig_w[["double"]], height = 85)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(vapply(arm_ens, length, integer(1)))   # CI 57 / CII 8 / CIII 15 / CIV 45 / CV 24
  print(n_per)
  print(brk)
  print(p)
}
