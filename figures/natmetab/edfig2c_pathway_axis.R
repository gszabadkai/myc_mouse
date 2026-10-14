# =============================================================================
# edfig2c_pathway_axis.R -- Extended Data Fig. 2c: the pathway-level PCA, PC1
# broken for the one far animal, beside the ranked loadings
# -----------------------------------------------------------------------------
# PROMOTED from proposals/ on 2026-10-14 (author: "promote them to final panels").
# It replaces the 2026-10-12 narrow version, retired/edfig2c_pathway_axis.R.
#
# Wider, with a broken PC1 axis (author,
# 2026-10-14: "the PC plot is really crammed into a little space with
# unproportionately large fonts ... the x axis can be broken for the outlier point
# and the panel should be wider", as the old Fig. 1C did)
# -----------------------------------------------------------------------------
# The same content and the same asserts as the retired narrow version. Three parts side by
# side: the 24 animals on PC1 and PC2 with the axis broken between the main cloud
# and the one 6W MYC+ animal far out on PC1 (it stays in every number); then the
# ranked absolute PC1 loadings.
#
# UNITS: the two PCA parts share one y scale and their data areas are FIXED in
# millimetres at the same millimetres per unit on PC1 and PC2, so the cloud has
# equal units -- the flatness is the 77% against the 8%. The loadings take the rest
# of the width.
#
# 2026-10-14 (author): lower ("could be done with reduced height"), the classes
# coloured mitochondrial against not (mito_class_cols, purple against grey), and
# more sets named at the top: MYC (Felsher), MYC V1 and proliferation beside E2F.
#
# Reads:  results/gsva_scores.rds     (script 15)
#         results/pathway_loading.rds (script 37)
# Output: outputs/natmetab/EDFig2/EDFig2c_pathway_axis.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

gs <- readRDS(here::here("results", "gsva_scores.rds"))
pl <- readRDS(here::here("results", "pathway_loading.rds"))

# --- the pathway-level PCA (script 37), rebuilt and asserted ---------------------
M_gsva <- gs$scores
expr   <- gs$expr_mat[, colnames(M_gsva), drop = FALSE]          # VST, genes x samples
sets   <- rownames(M_gsva)
pw     <- gs$pathways[intersect(names(gs$pathways), sets)]
stopifnot(identical(as.character(gs$sample_meta$sample), colnames(M_gsva)))
zg <- t(scale(t(expr)))
zg <- zg[is.finite(rowSums(zg)), , drop = FALSE]
M  <- t(vapply(sets, function(s) {
  g <- intersect(pw[[s]], rownames(zg))
  if (length(g) < 5L) rep(NA_real_, ncol(zg)) else colMeans(zg[g, , drop = FALSE])
}, numeric(ncol(zg))))
keep <- intersect(rownames(M)[stats::complete.cases(M)],
                  rownames(M_gsva)[stats::complete.cases(M_gsva)])
M  <- M[keep, , drop = FALSE]
pc <- stats::prcomp(t(M), center = TRUE, scale. = FALSE)
vfr <- pc$sdev^2 / sum(pc$sdev^2)
gm  <- colMeans(M)                                    # per-sample mean of all set scores
sc  <- pc$x
if (stats::cor(sc[, 1], gm) < 0) sc[, 1] <- -sc[, 1]  # PC1 high = high everywhere
lam <- stats::setNames(as.numeric(stats::cor(t(M), sc[, 1])), rownames(M))

ref_v <- pl$sample_pca_var[pl$sample_pca_var$space == "pathway_scores_884", ]
mito  <- intersect(gs$set_meta$set_name[gs$set_meta$category_primary == "MitoCarta"], keep)
ox    <- grep("OXPHOS|COMPLEX_[IV]|_SUBUNITS|ASSEMBLY_FACTORS|ELECTRON_CARRIERS|CRISTAE",
              mito, value = TRUE)
ox_gate <- stats::cor(colMeans(M[ox, , drop = FALSE]), gm)
conc    <- stats::cor(lam[pl$loading_movement$set], pl$loading_movement$lambda_design)
stopifnot(length(keep) == pl$n_sets, nrow(ref_v) == 1L,
          max(abs(100 * vfr[1:3] - c(ref_v$pc1_pct, ref_v$pc2_pct, ref_v$pc3_pct))) < 1e-6,
          abs(ox_gate - pl$ox_gate) < 1e-6,
          conc >= 0.9)                    # script 37 ranks on its design-residual axis

# --- LEFT: the animals, PC1 broken between the cloud and the one far animal ----------
sm  <- gs$sample_meta
pts <- data.frame(PC1 = sc[, 1], PC2 = sc[, 2],
                  group = factor(as.character(sm$group), levels = order_by_age))
stopifnot(nrow(pts) == 24L, !anyNA(pts$group))
cen <- stats::aggregate(cbind(PC1, PC2) ~ group, pts, mean)

far  <- which.max(pts$PC1)
gap  <- sort(pts$PC1, decreasing = TRUE)[1:2]
stopifnot(gap[1] > 35, gap[2] < 25,                          # one animal, alone, far out
          as.character(pts$group[far]) == "6W_pos",
          all(cen$PC1 < gap[2] + 1))                          # every group mean in the cloud
ypad <- 0.06 * diff(range(pts$PC2))
YLO <- min(pts$PC2) - ypad; YHI <- max(pts$PC2) + ypad
MAIN <- c(min(pts$PC1[-far]) - 1.5, gap[2] + 1.5)
OUTR <- c(floor(gap[1]) - 2, floor(gap[1]) + 2)              # a 4-unit window around it

pca_part <- function(d, xlim, main) {
  ggplot(d, aes(PC1, PC2, colour = group)) +
    geom_hline(yintercept = 0, linewidth = NM_LINE, colour = "grey85") +
    (if (main) geom_vline(xintercept = 0, linewidth = NM_LINE, colour = "grey85")) +
    geom_point(size = 1.1, stroke = 0) +
    (if (main) geom_point(data = cen, aes(fill = group), shape = 21, size = 2.4,
                          colour = "black", stroke = 0.4, show.legend = FALSE)) +
    scale_colour_manual(values = group_cols, labels = group_keys, name = NULL, drop = FALSE) +
    scale_fill_manual(values = group_cols, guide = "none") +
    scale_y_nm(step = 10, top = YHI, bottom = YLO, minor = TRUE, labels = lab_signed) +
    coord_cartesian(clip = "off")
}
# the break mark: two short parallel diagonals across the x axis line, drawn in panel
# coordinates (npc) so the axis limits cannot drop them; `at` = 0 (left end) or 1 (right)
brk <- function(at) {
  mm <- function(v) grid::unit(v, "mm")
  g <- grid::segmentsGrob(x0 = grid::unit(at, "npc") + mm(c(-1.1, -0.3)),
                          x1 = grid::unit(at, "npc") + mm(c(-0.3, 0.5)),
                          y0 = grid::unit(0, "npc") - mm(1), y1 = grid::unit(0, "npc") + mm(1),
                          gp = grid::gpar(lwd = NM_LINE * ggplot2::.pt, col = "black"))
  annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
}
p_main <- pca_part(pts[-far, ], MAIN, TRUE) +
  scale_x_nm(step = 10, top = MAIN[2], bottom = MAIN[1], minor = TRUE, labels = lab_signed) +
  labs(x = sprintf("PC1 (%.0f%%)", 100 * vfr[1]), y = sprintf("PC2 (%.0f%%)", 100 * vfr[2])) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.5))) +
  theme_nm(legend = "top") +
  theme(legend.key.spacing.x = unit(1, "mm"), plot.margin = margin(1, 0.5, 1, 1, "mm"))
xr_main <- layer_scales(p_main)$x$get_limits()[2]
p_main <- p_main + brk(1)
p_far <- pca_part(pts[far, , drop = FALSE], OUTR, FALSE) +
  scale_x_continuous(limits = OUTR, breaks = c(OUTR[1], round(gap[1]), OUTR[2]),
                     labels = c("", sprintf("+%d", round(gap[1])), ""), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_nm() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), plot.margin = margin(1, 3, 1, 0.5, "mm"))
# (one break mark, on the main part, is enough)
# Equal units: the data area is PANEL_H_MM tall for the PC2 span, and each PCA part
# is as wide as its PC1 span at the same millimetres per unit.
PANEL_H_MM <- 30
span_y  <- diff(.nm_axis(10, YHI, YLO, TRUE)$limits)
mm_unit <- PANEL_H_MM / span_y
W_MAIN_MM <- diff(layer_scales(p_main)$x$get_limits()) * mm_unit
W_FAR_MM  <- diff(OUTR) * mm_unit
stopifnot(W_MAIN_MM > 40, W_MAIN_MM < 80)

# --- RIGHT: what the axis is made of --------------------------------------------------
cls <- as.data.frame(pl$mito_classification)
ld  <- data.frame(set = names(lam), abs_load = abs(as.numeric(lam)))
ld  <- merge(ld, cls[, c("set", "class3")], by = "set")
ld  <- ld[order(-ld$abs_load), ]
ld$rank   <- seq_len(nrow(ld))
ld$class3 <- factor(ld$class3, levels = names(mito_class_cols))
stopifnot(nrow(ld) == length(keep), !anyNA(ld$class3))

# Named sets (author, 2026-10-14: "label more gene sets, I'd also add MYC (Felsher),
# MYC v1, proliferation (although E2F is already there)"):
#   OXPHOS          MitoCarta's OXPHOS, the mitochondrial programme the text names
#   MYC (Felsher)   the Felsher integrative MYC signature
#   MYC V1          MYC_UP.V1_UP (MSigDB C6). The library has no Hallmark MYC
#                   targets V1, only V2 (rank 334, 0.93; the legend can say so)
#   Cell cycle      Reactome's cell cycle, the proliferation programme beside E2F
#   E2F targets     Hallmark E2F
#   Glutathione     the redox control (the legend says so)
# The non-mitochondrial programme that loads as high as OXPHOS (the human
# breast-cancer MYC bicluster, 0.98) stays in the legend.
# LAYOUT: a staircase under the curve. Each label sits left of its own leader, and
# the further right a set ranks the lower its label, so no label crosses a leader.
# Glutathione's label goes up and right, above the falling end of the curve.
lab <- data.frame(
  set   = c("MITOCARTA_OXPHOS", "MYC_felsher_integrative_signature", "MYC_MYC_UP.V1_UP",
            "PROLIF_CELL_CYCLE_REACTOME", "PROLIF_E2F_HALLMARK", "GS_METAB_GLUTATHIONE"),
  label = c("OXPHOS", "MYC (Felsher)", "MYC V1", "Cell cycle", "E2F targets", "Glutathione"),
  ly    = c(0.78, 0.64, 0.50, 0.36, 0.22, 0.80),
  hj    = c(1, 1, 1, 1, 1, 0),
  stringsAsFactors = FALSE)
lab <- merge(lab, ld, by = "set")
lab <- lab[order(lab$rank), ]
lab$ye <- lab$ly                                  # the leader stops at its label
lab$lx <- ifelse(lab$hj == 1, lab$rank - 8, lab$rank + 8)
stopifnot(nrow(lab) == 6L,
          all(diff(lab$ly[lab$hj == 1]) < 0),     # further right, lower: no crossings
          all(lab$abs_load[lab$hj == 1] > lab$ly[lab$hj == 1] + 0.05))

p_load <- ggplot(ld, aes(rank, abs_load)) +
  geom_point(aes(colour = class3), size = 0.7, stroke = 0) +
  geom_segment(data = lab, aes(x = rank, xend = rank, y = abs_load, yend = ye),
               linewidth = 0.25, colour = "grey45") +
  geom_point(data = lab, size = 1, colour = "black", stroke = 0) +
  geom_text(data = lab, aes(x = lx, y = ly, label = label, hjust = hj), vjust = 0.5,
            lineheight = 0.85, size = pt(NM_TXT), family = NM_FONT) +
  scale_colour_manual(values = mito_class_cols, guide = "none") +
  scale_x_nm(step = 200, top = nrow(ld), bottom = 0, minor = TRUE, labels = lab_num) +   # ends at 900
  scale_y_nm(step = 0.5, top = 1, minor = TRUE, labels = lab_num) +
  labs(x = "Gene sets, ranked", y = "PC1 loading (absolute)") +
  theme_nm()

# THE TEXT'S NUMBERS (suggested 2026-10-12): "Across 885 gene sets, a single axis
# carried 77% of the variance; it separated MYC+ from WT glands at both ages and was
# led by mitochondrial gene sets" -- asserted here, where they are computed.
sep_p <- vapply(c("6W", "12W"), function(a) {
  k <- sm$timepoint == a
  stats::t.test(sc[k & sm$myc_status == "pos", 1], sc[k & sm$myc_status == "neg", 1])$p.value
}, numeric(1))
gm1 <- tapply(sc[, 1], paste(sm$timepoint, sm$myc_status), mean)
stopifnot(nrow(M) == 885L, round(100 * vfr[1]) == 77,
          all(sep_p < 0.05), gm1[["6W pos"]] > gm1[["6W neg"]], gm1[["12W pos"]] > gm1[["12W neg"]],
          all(ld$class3[1:50] != "non_mito"))

save_panel(stack_parts(p_main, p_far, p_load, ncol = 3,
                       widths = grid::unit(c(W_MAIN_MM, W_FAR_MM, 1), c("mm", "mm", "null")),
                       heights = grid::unit(PANEL_H_MM, "mm")),
           fig = "EDFig2", panel = "c", name = "pathway_axis", width = 176, height = 46)

cat("\nED 2c -- named sets (rank of", nrow(ld), "by |PC1 loading|):\n")
for (i in seq_len(nrow(lab))) cat(sprintf("  %-16s %-40s rank %3d  loading %.3f  %s\n", lab$label[i],
                                          lab$set[i], lab$rank[i], lab$abs_load[i], lab$class3[i]))
cat(sprintf("  top 100 by class: %s\n", paste(names(table(ld$class3[1:100])), table(ld$class3[1:100]),
                                                sep = " ", collapse = ", ")))

