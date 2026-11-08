# =============================================================================
# p5d_expression_options.R -- PROPOSAL: how ED Fig. 5d shows the mRNA changes of the
#                             pre-specified invasion genes (script 58, PART 5)
# -----------------------------------------------------------------------------
# The author (2026-11-07): the TF activities do not test the hypothesis -- a
# bioenergetics-independent PGC-1a effect on invasion -- the mRNA changes do. "Would
# a heatmap with the three experimental groups be informative? Labelling the cells
# with asterisk that are significant, using the usual heatmap colour scheme? Or are
# there alternative plans?" Three forms, at printed size, then on the ED 5 page:
#   H1  the author's heatmap, all 35 genes: Fig. 1g's form (VST, z per gene over the
#       21 tumours, the arm mean in each tile; espresso-white-mint), rows grouped by
#       script 58's five modules; a star in the Bcl-xL + Pgc-1a tile where that arm
#       differs from BOTH controls (IHW-adjusted P < 0.05 in each Wald contrast);
#       an arrow after each gene name for the direction the literature predicts.
#   H2  the same, compact: modules A, B, D and E gene by gene, the 22-gene Luo 2016
#       signature as one row (the mean of its genes), and script 52's compartment
#       markers as context rows (no test); tile numbers as Fig. 1g prints them.
#   F   a forest plot instead: each gene's log2 fold change with its 95% Wald
#       interval in both contrasts, filled where IHW-adjusted P < 0.05 (Fig. 1e's
#       rule), the clean contrast (against Bcl-xL + EV) larger (the size rule).
#
# Reads : results/orthotopic_tf_activity.rds (script 58), or the file named by
#         options(natmetab.s58_rds = ...); options(natmetab.ed5_png = ...) for the
#         page mock-ups (the assembled ED 5 rendered to an A4 PNG)
# Writes: outputs/natmetab/proposals/p5d_expr/EDFig5d_<h1|h2|forest>.pdf
#         outputs/natmetab/proposals/p5d_expression_options.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))
suppressPackageStartupMessages(library(grid))

S58  <- getOption("natmetab.s58_rds", here::here("results", "orthotopic_tf_activity.rds"))
ED5  <- getOption("natmetab.ed5_png", NA_character_)
OUTD <- here::here("outputs", "natmetab", "proposals")
o <- readRDS(S58)
stopifnot(!is.null(o$heatmap))

TXT   <- 5
TITLE <- 6
LINE  <- lw_pt(NM_LINE_PRISM[["Fig3"]])
SH    <- "\u00ad"                                     # soft hyphen
ARMS  <- c("EV", "BclxL", "Pgc1a_BclxL")
ARM_LAB <- c(EV = 'atop("EV", "+ EV")',
             BclxL = sprintf('atop("Bcl%sxL", "+ EV")', SH),
             Pgc1a_BclxL = sprintf('atop("Bcl%sxL", "+ Pgc%s1"*alpha)', SH, SH))
MOD_LAB <- c(A = "ID proteins", B = sprintf("TCF4%sTWIST1\nprogramme", SH),
             C = "Luo 2016\nsignature", D = "MMPs", E = "SNAIL",
             comp = "Composition\n(markers)")
COMP_LAB <- c(epithelial = "Epithelium", stromal = "Stroma", endothelial = "Endothelium",
              adipose = "Adipose", immune = "Immune")

G   <- as.data.frame(o$genes$spec)
gt  <- as.data.frame(o$genes$readout)
hm  <- as.data.frame(o$heatmap$genes)
hs  <- as.data.frame(o$heatmap$sets)
sig <- as.data.frame(o$heatmap$sig)
sig$both <- sig$sig_vs_bclxl & sig$sig_vs_ev
fc1 <- gt[gt$contrast == "vs_bclxl", ]

# gene order within a module: the Bcl-xL + Pgc-1a effect against Bcl-xL + EV, largest first
gene_order <- function(mods) {
  g <- G[G$module %in% mods, ]
  g$lfc <- fc1$log2FC[match(g$gene, fc1$gene)]
  g$gene[order(match(g$module, mods), -g$lfc)]
}
LIM <- 1.5
stopifnot(max(abs(hm$mean_z[hm$block == "genes"])) <= LIM + 0.25)

# --- labels: genes italic with the predicted direction as an arrow; set rows plain ---
row_labels <- function(x) {
  lab  <- ifelse(x %in% G$gene, sprintf('italic("%s")', x), sprintf('"%s"', x))
  pred <- G$predicted[match(x, G$gene)]
  arr  <- ifelse(is.na(pred), "",
                 ifelse(pred == "none", '~~phantom("" %up% "")',     # keeps Tcf4 in line
                 ifelse(pred == "up", '~~"" %up% ""', '~~"" %down% ""')))
  parse(text = paste0(lab, arr))
}
SIG_ROW <- "Signature (22 genes)"
MOD_LAB_COMPACT <- replace(MOD_LAB, "C", "Luo 2016")

# --- the heatmap (H1, H2) --------------------------------------------------------
heat <- function(compact, row_mm) {
  if (!compact) {
    rows <- gene_order(c("A", "B", "C", "D", "E"))
    d <- hm[hm$block == "genes" & hm$gene %in% rows, ]
    d$row <- d$gene; d$module <- G$module[match(d$gene, G$gene)]
    mods <- c("A", "B", "C", "D", "E"); labs <- MOD_LAB
  } else {
    rows_g <- gene_order(c("A", "B", "D", "E"))
    d1 <- hm[hm$block == "genes" & hm$gene %in% rows_g, ]
    d1$row <- d1$gene; d1$module <- G$module[match(d1$gene, G$gene)]
    d2 <- hs[hs$block == "genes" & hs$set == "C", ]
    d2$gene <- NA; d2$row <- SIG_ROW; d2$module <- "C"
    d3 <- hs[hs$block == "composition" & hs$set %in% names(COMP_LAB), ]
    d3$gene <- NA; d3$row <- unname(COMP_LAB[d3$set]); d3$module <- "comp"
    keep <- c("row", "gene", "module", "arm", "mean_z")
    d <- rbind(d1[, keep], d2[, keep], d3[, keep])
    mg <- G$module[match(rows_g, G$gene)]
    rows <- c(rows_g[mg %in% c("A", "B")], SIG_ROW, rows_g[mg %in% c("D", "E")], unname(COMP_LAB))
    mods <- c("A", "B", "C", "D", "E", "comp"); labs <- MOD_LAB_COMPACT
  }
  d$module <- factor(d$module, levels = mods, labels = labs[mods])
  d$row <- factor(d$row, levels = rev(rows))
  d$arm <- factor(d$arm, levels = ARMS)
  d$star <- ifelse(d$arm == "Pgc1a_BclxL" & !is.na(d$gene) & d$gene %in% sig$gene[sig$both], "*", "")
  d$lab <- if (compact) sprintf("%+.1f", d$mean_z) else ""
  if (compact) d$lab[d$lab %in% c("+0.0", "-0.0")] <- "0"
  d$lab <- paste0(d$lab, d$star)
  # the column headers, two lines set 2.0 mm apart in the top facet (plotmath's atop()
  # spaces them about twice that, FIGURE_RULES), drawn past the panel's top edge
  n_top <- length(unique(d$row[d$module == levels(d$module)[1]]))
  y2 <- n_top + 0.5 + 0.8 / row_mm
  hdr <- data.frame(module = factor(levels(d$module)[1], levels = levels(d$module)),
                    arm = factor(rep(ARMS, 2), levels = ARMS),
                    y = rep(c(y2 + 2.0 / row_mm, y2), each = 3),
                    lab = c('"EV"', sprintf('"Bcl%sxL"', SH), sprintf('"Bcl%sxL"', SH),
                            '"+ EV"', '"+ EV"', sprintf('"+ Pgc%s1"*alpha', SH)))
  ggplot(d, aes(arm, row)) +
    geom_tile(aes(fill = mean_z), colour = "white", linewidth = 0.3) +
    geom_text(aes(label = lab), size = pt(TXT), family = NM_FONT, colour = "black",
              vjust = if (compact) 0.5 else 0.75) +
    geom_text(data = hdr, aes(x = arm, y = y, label = lab), inherit.aes = FALSE, parse = TRUE,
              vjust = 0, size = pt(TXT), family = NM_FONT, colour = "black") +
    coord_cartesian(clip = "off") +
    facet_grid(module ~ ., scales = "free_y", space = "free_y") +
    scale_fill_div(c(-LIM, LIM), name = sprintf("z%sscore", SH), breaks = c(-LIM, 0, LIM)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(labels = row_labels, expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_colourbar(barwidth = unit(18, "mm"), barheight = unit(1.6, "mm"),
                                  ticks.colour = "grey30", title.position = "left",
                                  title.vjust = 1)) +
    theme_nm(legend = "bottom", txt = TXT, txt_title = TITLE) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          strip.text.y = element_text(size = TITLE, angle = 0, hjust = 0, lineheight = 0.9),
          panel.spacing.y = unit(0.8, "mm"),
          legend.title = element_text(size = TXT),
          plot.margin = margin(5.2, 0.5, 0.5, 0.5, "mm"))
}

# --- the forest (F) --------------------------------------------------------------
forest <- function() {
  rows <- gene_order(c("A", "B", "C", "D", "E"))
  d <- gt[gt$status == "fitted", ]
  d$module <- factor(MOD_LAB[d$module], levels = MOD_LAB[c("A", "B", "C", "D", "E")])
  d$gene <- factor(d$gene, levels = rev(rows))
  d$contrast <- factor(d$contrast, levels = c("vs_bclxl", "vs_ev"))
  d$filled <- !is.na(d$padj) & d$padj < 0.05
  d$col <- ifelse(d$contrast == "vs_bclxl", "#000000", "#8C8C8C")
  d$fill <- ifelse(d$filled, d$col, "white")
  d$size <- ifelse(d$contrast == "vs_bclxl", 1.4, 1.0)
  pd <- position_dodge(width = 0.55)
  lo <- min(d$lo95); hi <- max(d$hi95); pad <- 0.04 * (hi - lo)
  ggplot(d, aes(log2FC, gene, group = contrast)) +
    geom_vline(xintercept = 0, linewidth = LINE, colour = "black") +
    geom_errorbar(aes(xmin = lo95, xmax = hi95, colour = col), width = 0, linewidth = LINE * 0.8,
                  position = pd, orientation = "y") +
    geom_point(aes(fill = fill, colour = col, size = size), shape = 21, stroke = 0.35, position = pd) +
    scale_colour_identity() + scale_fill_identity() + scale_size_identity() +
    facet_grid(module ~ ., scales = "free_y", space = "free_y") +
    scale_x_nm(step = 2, top = hi + pad, bottom = lo - pad, minor = TRUE, labels = lab_signed) +
    scale_y_discrete(labels = row_labels) +
    labs(x = "Expression (log2 fold change)", y = NULL) +
    theme_nm(txt = TXT, txt_title = TITLE, line = LINE) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          strip.text.y = element_text(size = TITLE, angle = 0, hjust = 0, lineheight = 0.9),
          panel.spacing.y = unit(0.8, "mm"),
          plot.margin = margin(0.5, 1.5, 0.5, 0.5, "mm"))
}

# --- save each (every rule checked), capture for the sheet ---------------------------
op <- options(natmetab.capture = TRUE)
P <- list(); S <- list()
add <- function(key, p, w, h) {
  st <- fixed_size(p, w, h)
  sz <- ceiling(natural_size(st) * 2) / 2
  stopifnot(!anyNA(sz))
  save_panel(st, fig = "EDFig5", panel = "d", name = key, width = sz[["width"]],
             height = sz[["height"]], dir = file.path(OUTD, "p5d_expr"))
  P[[key]] <<- get(".nm_captured", envir = globalenv()); S[[key]] <<- sz
}
n_full <- nrow(G); n_compact <- 13 + 1 + length(COMP_LAB)
add("h1_heatmap_all35", heat(FALSE, 2.1), 3 * 7, n_full * 2.1)
add("h2_heatmap_compact", heat(TRUE, 2.6), 3 * 7.5, n_compact * 2.6)
add("forest_all35", forest(), 36, n_full * 2.6)
options(op)

# --- the sheet -----------------------------------------------------------------
mm <- function(v) unit(v, "mm")
hy <- function(x) gsub("(?<=[A-Za-z0-9])-(?=[A-Za-z0-9])", SH, x, perl = TRUE)
PW <- 210; PH <- 297
place <- function(p, x, y, sz) print(p, vp = viewport(x = mm(x), y = mm(PH - y), width = mm(sz[["width"]]),
                                     height = mm(sz[["height"]]), just = c("left", "top")), newpage = FALSE)
txt <- function(label, x, y, size = 7, face = "plain")
  grid.text(hy(label), x = mm(x), y = mm(PH - y), just = c("left", "top"),
            gp = gpar(fontsize = size, fontface = face, fontfamily = NM_FONT))
blurb <- list(
  h1_heatmap_all35 = c("H1  The heatmap as asked: all 35 genes in script 58's five modules.",
                       "Tiles: VST, z per gene over the 21 tumours, the arm mean (Fig. 1g's form), espresso-white-mint +/-1.5.",
                       "* in the Bcl-xL + Pgc-1a tile: differs from both controls (DESeq2 Wald, IHW-adjusted P < 0.05 in each).",
                       "Arrows after the gene names: the direction the literature predicts in Bcl-xL + Pgc-1a tumours (none for Tcf4)."),
  h2_heatmap_compact = c("H2  Compact: modules A, B, D, E gene by gene; the Luo 2016 signature as one row",
                         "(the mean of its 22 genes); script 52's compartment markers as context rows (no test).",
                         "Tile numbers as Fig. 1g prints them; the star and arrows as in H1."),
  forest_all35 = c("F  A forest instead: log2 fold change of Bcl-xL + Pgc-1a against Bcl-xL + EV (black, larger)",
                   "and against EV + EV (grey), with 95% Wald intervals; filled = IHW-adjusted P < 0.05 (Fig. 1e's rule).",
                   "Arrows after the gene names: the predicted direction."))
pdf(file.path(OUTD, "p5d_expression_options.pdf"), width = PW / 25.4, height = PH / 25.4,
    family = NM_FONT, useDingbats = FALSE)
for (key in names(P)) {
  grid.newpage()
  txt("ED Fig. 5d, the invasion genes' mRNA: three forms at printed size", 15, 10, 10, "bold")
  txt(sprintf("Script 58 (%s); 21 tumours, 7 per arm.",
              if (identical(S58, here::here("results", "orthotopic_tf_activity.rds"))) "the record" else "a dry run"),
      15, 16, 6.5)
  for (i in seq_along(blurb[[key]])) txt(blurb[[key]][i], 15, 22 + 3.4 * (i - 1), 7)
  place(P[[key]], 15, 40, S[[key]])
  txt(sprintf("%.1f x %.1f mm", S[[key]][["width"]], S[[key]][["height"]]), 15 + S[[key]][["width"]] + 3, 40, 6.5)
}
if (!is.na(ED5) && file.exists(ED5)) {
  bg <- png::readPNG(ED5, native = TRUE)
  for (key in c("h2_heatmap_compact", "h1_heatmap_all35", "forest_all35")) {
    grid.newpage()
    grid.raster(bg, x = 0, y = 0, width = 1, height = 1, just = c("left", "bottom"))
    grid.rect(x = mm(24.5), y = mm(PH - 75), width = mm(98), height = mm(38.5),
              just = c("left", "top"), gp = gpar(fill = "white", col = NA))
    place(P[[key]], 26, 76, S[[key]])
    txt(sprintf("%s in place of the TF charts on the assembled ED 5 (2026-11-06): %.0f x %.0f mm",
                toupper(sub("_.*", "", key)), S[[key]][["width"]], S[[key]][["height"]]), 15, 8, 8, "bold")
  }
}
invisible(dev.off())
message("wrote ", file.path(OUTD, "p5d_expression_options.pdf"))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(heat(FALSE)); print(heat(TRUE)); print(forest())
  sig[sig$sig_vs_bclxl | sig$sig_vs_ev, ]
  S
}
