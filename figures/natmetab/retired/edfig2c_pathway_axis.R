# RETIRED 2026-10-14 -- replaced by the wide edfig2c_pathway_axis.R. Kept for the record;
# not a panel, and outside the dry-run set.
# =============================================================================
# edfig2c_pathway_axis.R -- Extended Data Fig. 2c: one axis carries most of the
# variance across gene sets; it separates the genotypes and mitochondrial sets
# lead it
# -----------------------------------------------------------------------------
# The text (suggested 2026-10-12): "Across 885 gene sets, a single axis carried
# 77% of the variance; it separated MYC+ from WT glands at both ages and was led
# by mitochondrial gene sets (ED Fig. 2a-c)." From the old Fig. 1C (the PCA) and
# Fig. S1C (its loadings), one PDF, two parts side by side:
#
#   LEFT   the 24 animals on PC1 and PC2 of the pathway-level PCA, equal units on
#          both axes, so the flatness is the 77% against 8%. Ringed: group means.
#   RIGHT  every gene set's absolute loading on PC1, ranked, in the mitochondrial
#          classes of ED 2b (same colours; ED 2b carries the key).
#
# THE PCA is script 37's: the linear mean gene-wise z-score per set (885 sets
# scored in both matrices), PCA across the 24 samples. Rebuilt here and asserted
# against results/pathway_loading.rds (the rebuild of the old layer's
# pathway_axis(), without its freshness check).
#
# WHAT BOUNDS "LED BY MITOCHONDRIAL SETS" (all printed below, none drawn):
#   - loading is co-variation, not cause (script 37's enrichment verdict);
#   - the axis is common-mode: r 0.98 with the mean of all 885 scores, 79% of
#     sets load one way, median |loading| 0.88;
#   - 47 of the top 50 are transcription-factor target sets intersected with
#     MitoCarta, mitochondrial by construction;
#   - non-mitochondrial growth programmes reach 0.98 (the labelled human
#     breast-cancer MYC bicluster), and glutathione, a redox set, loads 0.67.
# PC1 is a genotype axis, not an age axis: batch = age, so this is the reading
# that survives.
#
# Reads:  results/gsva_scores.rds     (script 15: VST matrix, gene sets, sample sheet)
#         results/pathway_loading.rds (script 37: assertions, mitochondrial classes)
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

# --- LEFT: the animals ---------------------------------------------------------------
sm  <- gs$sample_meta
pts <- data.frame(PC1 = sc[, 1], PC2 = sc[, 2],
                  group = factor(as.character(sm$group), levels = order_by_age))
stopifnot(nrow(pts) == 24L, !anyNA(pts$group))
cen <- stats::aggregate(cbind(PC1, PC2) ~ group, pts, mean)

xpad <- 0.03 * diff(range(pts$PC1)); ypad <- 0.06 * diff(range(pts$PC2))
p_pca <- ggplot(pts, aes(PC1, PC2, colour = group)) +
  geom_hline(yintercept = 0, linewidth = NM_LINE, colour = "grey85") +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = "grey85") +
  geom_point(size = 0.9, stroke = 0) +
  geom_point(data = cen, aes(fill = group), shape = 21, size = 2, colour = "black",
             stroke = 0.4, show.legend = FALSE) +
  scale_colour_manual(values = group_cols, labels = group_keys, name = NULL) +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_x_nm(step = 20, top = max(pts$PC1) + xpad, bottom = min(pts$PC1) - xpad,
             minor = TRUE, labels = lab_signed) +
  scale_y_nm(step = 10, top = max(pts$PC2) + ypad, bottom = min(pts$PC2) - ypad,
             minor = TRUE, labels = lab_signed) +
  coord_equal() +
  labs(x = sprintf("PC1 (%.0f%%)", 100 * vfr[1]), y = sprintf("PC2 (%.0f%%)", 100 * vfr[2])) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.5))) +
  theme_nm(legend = "top") +
  theme(legend.key.spacing.x = unit(0.5, "mm"))

# --- RIGHT: what the axis is made of --------------------------------------------------
cls <- as.data.frame(pl$mito_classification)
ld  <- data.frame(set = names(lam), abs_load = abs(as.numeric(lam)))
ld  <- merge(ld, cls[, c("set", "class3")], by = "set")
ld  <- ld[order(-ld$abs_load), ]
ld$rank   <- seq_len(nrow(ld))
ld$class3 <- factor(ld$class3, levels = names(mito_class_cols))
stopifnot(nrow(ld) == length(keep), !anyNA(ld$class3))

# Three named sets: OXPHOS, the E2F targets (a growth programme), and glutathione,
# the redox control (the legend says so). Four did not fit at 30 mm (tried
# 2026-10-12): the non-mitochondrial programme that loads as high as OXPHOS (the
# human breast-cancer MYC bicluster, 0.98) is given in the legend instead.
# Labels go where the curve leaves room: OXPHOS below it, E2F above it, and
# glutathione right-aligned under its own point, so its leader drops onto the end
# of the word and the word stays clear of the OXPHOS leader on the left.
# Positions are data units: lx/ly the label's anchor, hj its justification.
lab <- data.frame(
  set   = c("MITOCARTA_OXPHOS", "PROLIF_E2F_HALLMARK", "GS_METAB_GLUTATHIONE"),
  label = c("OXPHOS", "E2F targets", "Glutathione"),
  ly    = c(0.40, 0.95, 0.54),
  hj    = c(0, 0, 1),
  ye    = c(0.40, 0.95, 0.58),        # where the leader stops: glutathione's on the word
  stringsAsFactors = FALSE)
lab <- merge(lab, ld, by = "set")
lab$lx <- ifelse(lab$hj == 0, lab$rank + 18, lab$rank + 25)
stopifnot(nrow(lab) == 3L)

p_load <- ggplot(ld, aes(rank, abs_load)) +
  geom_point(aes(colour = class3), size = 0.7, stroke = 0) +
  geom_segment(data = lab, aes(x = rank, xend = rank, y = abs_load, yend = ye),
               linewidth = 0.25, colour = "grey45") +
  geom_point(data = lab, size = 1, colour = "black", stroke = 0) +
  geom_text(data = lab, aes(x = lx, y = ly, label = label, hjust = hj), vjust = 0.5,
            lineheight = 0.85, size = pt(NM_TXT), family = NM_FONT) +
  scale_colour_manual(values = mito_class_cols, guide = "none") +
  scale_x_nm(step = 400, top = nrow(ld), bottom = 0, minor = TRUE, labels = lab_num) +
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

save_panel(stack_parts(p_pca, p_load, ncol = 2, widths = c(1.15, 1)),
           fig = "EDFig2", panel = "c", name = "pathway_axis", width = 72, height = 40)

# --- numbers for the legend and the text (printed, never drawn) -------------------
fit <- function(f) stats::coef(summary(stats::lm(f, data = data.frame(
  PC1 = pts$PC1, age = factor(sm$timepoint, levels = c("6W", "12W")),
  genotype = factor(sm$myc_status, levels = c("neg", "pos"))))))
add <- fit(PC1 ~ age + genotype); int <- fit(PC1 ~ age * genotype)
per_age <- vapply(c("6W", "12W"), function(a) {
  k <- sm$timepoint == a
  stats::t.test(pts$PC1[k & sm$myc_status == "pos"], pts$PC1[k & sm$myc_status == "neg"])$p.value
}, numeric(1))
gmean <- tapply(pts$PC1, pts$group, mean)
top50 <- table(ld$class3[1:50])
cat("\nED Fig. 2c -- for the legend\n")
cat(sprintf("  %d gene sets, each the mean gene-wise z-score per sample (script 37); PCA over 24 animals\n", nrow(M)))
cat(sprintf("  PC1 %.1f%%, PC2 %.1f%%, PC3 %.1f%% of the variance\n", 100 * vfr[1], 100 * vfr[2], 100 * vfr[3]))
cat(sprintf("  PC1 ~ age + genotype: genotype P = %.2g, age P = %.2g; interaction P = %.2g\n",
            add["genotypepos", 4], add["age12W", 4], int["age12W:genotypepos", 4]))
cat(sprintf("  MYC+ vs WT within age (Welch): 6W P = %.2g, 12W P = %.2g\n", per_age[["6W"]], per_age[["12W"]]))
cat(sprintf("  group means on PC1: %s\n", paste(sprintf("%s %+.1f", group_keys[names(gmean)], gmean), collapse = ", ")))
cat(sprintf("  the axis is common-mode: r(PC1, mean of all sets) = %.3f; %.0f%% of loadings one-signed; median |loading| %.2f\n",
            stats::cor(pts$PC1, gm), 100 * mean(lam > 0), stats::median(ld$abs_load)))
cat(sprintf("  top 50 by |loading|: %s\n", paste(sprintf("%s %d", names(top50), top50), collapse = ", ")))
cat(sprintf("  median |loading| by class: %s\n",
            paste(sprintf("%s %.2f", levels(ld$class3), tapply(ld$abs_load, ld$class3, stats::median)), collapse = ", ")))
cat("  named: ", paste(sprintf("%s rank %d, %.2f", lab$set, lab$rank, lab$abs_load), collapse = "; "), "\n")
mb2 <- ld[ld$set == "METABRIC_MB2_HI_CV_GROUP1", ]
cat(sprintf("  not named, for the legend: the human breast-cancer MYC bicluster (MB2) loads %.2f\n    (rank %d), as high as OXPHOS\n", mb2$abs_load, mb2$rank))
cat("  loading is co-variation, not cause; batch = age, and PC1 is a genotype axis\n")

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(.nm_drawable(stack_parts(p_pca, p_load, ncol = 2, widths = c(1.15, 1))))
  head(ld, 20)
}
