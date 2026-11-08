# =============================================================================
# edfig5d_invasion_genes.R -- ED Fig. 5d: the pre-specified invasion genes in the
#                             three orthotopic arms (heatmap)
# -----------------------------------------------------------------------------
# THE AUTHOR'S CHOICE (2026-11-07): form H1 of proposals/p5d_expression_options.R --
# all 35 genes of script 58's five groups, "that spares an additional table and
# anyway makes the point". The TF charts that held this slot are DROPPED: they do
# not test the hypothesis, a PGC-1a effect on invasion separate from its
# mitochondrial role; the mRNA does. The author's changes to H1:
#   - MMPs above the signature;
#   - the signature named for what it holds AND as a signature: "Adhesion, TGF-beta
#     and Wnt signature", after Luo et al.'s own description (genes of focal adhesion
#     and ECM interactions, integrins, and TGF-beta and Wnt signalling); the citation
#     is in the legend. It is judged ONLY AS A SET (author, 2026-11-08): a single gene
#     moving the predicted way carries no weight, so its genes get no arrow and no
#     bold; the set's arrow follows its name, and its fgsea test goes in the legend;
#   - a gene name in BOLD where the change is in the predicted direction;
#   - a STAR where Bcl-xL + Pgc-1a differs from at least one control;
#   - colour limits +/-1.5.
# Both marks are decided in script 58 (PART 6b): star, IHW-adjusted P < 0.05 against
# at least one control (DESeq2 Wald); bold, significant in the predicted direction
# against at least one control and not significant against it in the other -- so a
# change of about zero cannot read as agreeing.
#
# FORM: Fig. 1g's. VST, z-scored per gene over the 21 tumours, the arm mean in each
#   tile; the diverging espresso-white-mint scale, white at 0. No tile numbers.
# ARROWS after each name: the direction the cited studies predict in Bcl-xL +
#   Pgc-1a tumours. Tcf4 has none (ID proteins act on the TCF4 protein), nor do the
#   signature's genes (set only); they carry an invisible arrow so the names stay in
#   line.
# ROWS: groups in the author's order; within a group, the effect against Bcl-xL + EV,
#   largest first.
# HEADERS AND GROUP NAMES are placed by calculation, not by plotmath's atop(), whose
#   two lines sit about twice as far apart as text lines (FIGURE_RULES), and so the
#   group names can carry Greek.
# TYPE: gene names, column headers and key labels 5 pt; group names and the key
#   title 6 pt (as Fig. 1g); stars 7 pt (as ED 2f-g's asterisks).
#
# NO TEST IS COMPUTED HERE: every value and flag is script 58's.
# A DRY RUN'S OBJECT NEVER WRITES INTO THE FIGURE FOLDER: when the object read is not
#   results/orthotopic_tf_activity.rds, the PDF goes to outputs/natmetab/proposals/preview/.
#
# Reads : results/orthotopic_tf_activity.rds (script 58), or options(natmetab.s58_rds = ...)
# Output: outputs/natmetab/EDFig5/EDFig5d_invasion_genes.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

REC <- here::here("results", "orthotopic_tf_activity.rds")
S58 <- getOption("natmetab.s58_rds", REC)
o   <- readRDS(S58)
is_record <- identical(normalizePath(S58), normalizePath(REC, mustWork = FALSE))
OUT <- if (is_record) NULL else here::here("outputs", "natmetab", "proposals", "preview")

TXT    <- 5      # pt; gene names, column headers, key labels
TITLE  <- 6      # pt; group names, key title
STAR   <- 7      # pt
LIM    <- 1.5
TILE_W <- 7      # mm
ROW_H  <- 2.1    # mm
SH     <- "\u00ad"                                  # soft hyphen: base pdf() prints "-" as a minus
ARMS   <- c("EV", "BclxL", "Pgc1a_BclxL")

G   <- as.data.frame(o$genes$spec)
ORD <- o$genes$module_order
sig <- as.data.frame(o$heatmap$sig)
hm  <- as.data.frame(o$heatmap$genes)
hm  <- hm[hm$block == "genes", ]
gt  <- as.data.frame(o$genes$readout)
fc1 <- gt[gt$contrast == "vs_bclxl", ]
stopifnot(identical(ORD, c("A", "B", "D", "C", "E")), nrow(G) == 35L,
          setequal(hm$gene, G$gene), setequal(sig$gene, G$gene),
          all(c("star", "bold") %in% names(sig)))
# the signature is judged only as a set (2026-11-08): an object that bolds any of its
# genes was written by script 58 before that rule, and is stale
if (any(sig$bold[sig$gene %in% G$gene[G$module == "C"]]))
  stop(S58, " was written before the signature rule of 2026-11-08 (it bolds signature genes: ",
       paste(sig$gene[sig$bold & sig$gene %in% G$gene[G$module == "C"]], collapse = ", "),
       "). Re-run scripts/58_orthotopic_tf_activity.R, then this script.", call. = FALSE)

# Mouse proteins in title case (the author's naming, 2026-11-06); MMPs as on the ED 5e
# blot; Wnt as Luo et al. write the pathway.
GROUP_LINES <- list(
  A = '"Id proteins"',
  B = c(sprintf('"Tcf4%sTwist1"', SH), '"programme"'),
  D = '"MMPs"',
  C = c(sprintf('"Adhesion, TGF%s"*beta', SH), '"and Wnt signature"~"" %down% ""'),   # the set's arrow
  E = '"Snail"')

G$lfc <- fc1$log2FC[match(G$gene, fc1$gene)]
G <- G[order(match(G$module, ORD), -G$lfc), ]

d <- hm[, c("gene", "arm", "mean_z")]
d$module <- factor(G$module[match(d$gene, G$gene)], levels = ORD)
d$gene   <- factor(d$gene, levels = rev(G$gene))      # first row on top
d$arm    <- factor(d$arm, levels = ARMS)
d$star   <- ifelse(d$arm == "Pgc1a_BclxL" & as.character(d$gene) %in% sig$gene[sig$star], "*", "")
stopifnot(!anyNA(d$module), !anyNA(d$arm), nrow(d) == 3L * 35L,
          max(abs(d$mean_z)) < LIM)                     # no tile clipped by the colour limits

gene_label <- function(x) {
  bold <- sig$bold[match(x, sig$gene)]
  pred <- G$predicted[match(x, G$gene)]
  pred[G$module[match(x, G$gene)] == "C"] <- "none"    # the signature's genes: no arrow of their own
  nm <- ifelse(bold, sprintf('bolditalic("%s")', x), sprintf('italic("%s")', x))
  ar <- ifelse(pred == "up", '~~"" %up% ""',
        ifelse(pred == "down", '~~"" %down% ""', '~~phantom("" %up% "")'))
  parse(text = paste0(nm, ar))
}

# column headers: two lines 2.0 mm apart, in the top group, just above its tiles
n_in <- table(factor(G$module, levels = ORD))
y2   <- n_in[[ORD[1]]] + 0.5 + 0.8 / ROW_H
hdr  <- data.frame(module = factor(ORD[1], levels = ORD), arm = factor(rep(ARMS, 2), levels = ARMS),
                   y = rep(c(y2 + 2.0 / ROW_H, y2), each = 3),
                   # phantom("g"): every header carries the descender of Pgc's g, so the
                   # bottom-aligned lines share one baseline
                   lab = paste0(c('"EV"', sprintf('"Bcl%sxL"', SH), sprintf('"Bcl%sxL"', SH),
                                  '"+ EV"', '"+ EV"', sprintf('"+ Pgc%s1"*alpha', SH)),
                                '*phantom("g")'))

# group names: right of the tiles, centred on their rows, lines 2.3 mm apart
GAP_MM <- 1.5; LINE_MM <- 2.3
grp <- do.call(rbind, lapply(ORD, function(m) {
  L <- GROUP_LINES[[m]]; k <- length(L)
  off <- (seq_len(k) - (k + 1) / 2) * LINE_MM / ROW_H
  data.frame(module = factor(m, levels = ORD), x = length(ARMS) + 0.5 + GAP_MM / TILE_W,
             y = (n_in[[m]] + 1) / 2 - off, lab = L)
}))
right_mm <- GAP_MM + 3 + max(text_mm(c("Adhesion, TGF-b", "Tcf4-Twist1", "and Wnt signature"), txt = TITLE)[, "w"])

p <- ggplot(d, aes(arm, gene)) +
  geom_tile(aes(fill = mean_z), colour = "white", linewidth = 0.3) +
  geom_text(aes(label = star), size = pt(STAR), family = NM_FONT, colour = "black", vjust = 0.75) +
  geom_text(data = hdr, aes(x = arm, y = y, label = lab), inherit.aes = FALSE, parse = TRUE,
            vjust = 0, size = pt(TXT), family = NM_FONT, colour = "black") +
  geom_text(data = grp, aes(x = x, y = y, label = lab), inherit.aes = FALSE, parse = TRUE,
            hjust = 0, size = pt(TITLE), family = NM_FONT, colour = "black") +
  coord_cartesian(clip = "off") +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_fill_div(c(-LIM, LIM), name = sprintf("z%sscore", SH), breaks = c(-LIM, 0, LIM)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(labels = gene_label, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(barwidth = unit(18, "mm"), barheight = unit(1.6, "mm"),
                                ticks.colour = "grey30", title.position = "left", title.vjust = 1)) +
  theme_nm(legend = "bottom", txt = TXT, txt_title = TITLE) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
        strip.text.y = element_blank(), strip.background = element_blank(),
        panel.spacing.y = unit(0.8, "mm"),
        legend.title = element_text(size = TITLE),
        plot.margin = margin(5.2, right_mm, 0.5, 0.5, "mm"))

panel <- fixed_size(p, length(ARMS) * TILE_W, nrow(G) * ROW_H)
size  <- ceiling(natural_size(panel) * 2) / 2
stopifnot(!anyNA(size))
save_panel(panel, fig = "EDFig5", panel = "d", name = "invasion_genes",
           width = size[["width"]], height = size[["height"]], dir = OUT)

# --- numbers for the legend (printed, never drawn) ----------------------------------
cat(sprintf("\nED 5d -- for the legend (%g x %g mm)%s\n", size[["width"]], size[["height"]],
            if (is_record) "" else "  [FROM A DRY RUN: not the record]"))
cat("  n = 7 tumours per arm (EV + EV, Bcl-xL + EV, Bcl-xL + Pgc-1a); 35 genes in five groups:\n")
for (m in ORD) cat(sprintf("    %-2s %-45s %2d genes\n", m, o$genes$modules[[m]], n_in[[m]]))
cat(sprintf("  tiles: VST, z per gene over 21 tumours, arm mean; colour limits +/-%.1f (largest |z| %.2f)\n",
            LIM, max(abs(d$mean_z))))
cat(sprintf("  star (differs from >= 1 control, IHW-adjusted P < 0.05): %d genes\n", sum(sig$star)))
cat(sprintf("  bold (change in the predicted direction): %d of %d genes judged gene by gene: %s\n",
            sum(sig$bold), sum(G$predicted != "none" & G$module != "C"), paste(sig$gene[sig$bold], collapse = ", ")))
cat("  the signature's set test, the legend's numbers (fgsea on the DESeq2 Wald statistic):\n")
print(o$genes$set_test[, c("contrast", "genes_present", "NES", "pval")])

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p)
  sig[sig$star, c("gene", "padj_vs_bclxl", "padj_vs_ev", "lfc_vs_bclxl", "lfc_vs_ev", "bold")]
  d[d$arm == "Pgc1a_BclxL", ][order(d$mean_z[d$arm == "Pgc1a_BclxL"]), ]
}
