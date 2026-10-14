# =============================================================================
# fig1g_bcl2_heatmap.R -- Fig. 1g: the BCL-2 family in the four groups (heatmap)
# -----------------------------------------------------------------------------
# PROMOTED from proposals/ on 2026-10-14 (author: "promote them to final panels").
# The author's plan (2026-10-14): "While at 6 weeks MYC shifted the balance in the
# transcriptome toward pro-apoptotic factors, at 12 weeks MYC-driven pro-apoptotic
# pathways were suppressed: expression of both ... declined compared to the
# normal gland, whereas ... remained stable (Fig. 1g)." The form is Fig. 2H+I
# (alt)'s top-right heatmap: each gene standardised over the 24 animals, the group
# mean in each tile.
#
# ROWS: script 55's BCL-2 family roster (16 genes read), grouped by ROLE so that
# the 6-week column reads as it should: the effectors and the activators rise
# under MYC while BCL-XL falls, and the sensitisers (BIK, BMF, NOXA) and BIM fall
# too -- without the grouping a reader sees "pro-apoptotic genes down" and doubts
# the balance. Within a role, rows run by MYC's 6-week effect, largest first.
# Activators are BID, BIM and PUMA (script 55; PUMA's place is debated).
#
# VALUES: VST expression (script 15's matrix), z-scored per gene over all 24
# animals, mean per group. Descriptive only: no test is drawn. The tests are in
# ED Fig. 3a (the animals) and ED Fig. 3b (the departures), and the legend. The two
# ages are two extraction batches, so the WT 6W -> 12W columns include batch.
#
# COLOUR: the diverging espresso-white-mint scale, white at 0, limits +/-2 so that
# the darkest tile (|z| 1.3) stays light enough for black numbers.
#
# Reads: results/gsva_scores.rds (script 15: $expr_mat, VST; $sample_meta)
#        results/departure_from_myc_fade.rds (script 55: the roster and roles)
# Output: outputs/natmetab/Fig1/Fig1g_bcl2_heatmap.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

gs  <- readRDS(here::here("results", "gsva_scores.rds"))
dep <- readRDS(here::here("results", "departure_from_myc_fade.rds"))
sm  <- as.data.frame(gs$sample_meta)
E   <- gs$expr_mat[, as.character(sm$sample), drop = FALSE]
grp <- factor(paste(sm$timepoint, sm$myc_status, sep = "_"), levels = order_by_geno)
stopifnot(!anyNA(grp), all(table(grp) == 6L))

fam <- as.data.frame(dep$rosters)
fam <- fam[fam$roster == "BCL-2 family" & is.na(fam$dropped), ]
stopifnot(nrow(fam) == 16L, all(fam$gene %in% rownames(E)))

ROLES <- c(effector = "Effectors", activator = "Activators", sensitiser = "Sensitisers",
           `anti-apoptotic` = "Anti­apoptotic")
stopifnot(setequal(fam$role, names(ROLES)))

Z  <- t(scale(t(E[fam$gene, , drop = FALSE])))
MZ <- t(apply(Z, 1, function(v) tapply(v, grp, mean)))

d <- do.call(rbind, lapply(seq_len(nrow(fam)), function(i)
  data.frame(gene = fam$gene[i], role = fam$role[i], lfc6 = fam$lfc6[i],
             group = factor(colnames(MZ), levels = order_by_geno), z = MZ[fam$gene[i], ])))
ord <- fam$gene[order(match(fam$role, names(ROLES)), -fam$lfc6)]
d$gene <- factor(d$gene, levels = rev(ord))                    # first row on top
d$role <- factor(ROLES[d$role], levels = ROLES)
d$lab  <- sprintf("%+.1f", d$z)
d$lab[d$lab %in% c("+0.0", "-0.0")] <- "0"

# --- the reading, asserted -------------------------------------------------------------
g <- function(x, grp) d$z[d$gene == x & d$group == grp]
stopifnot(g("Bax", "6W_pos") > g("Bax", "6W_neg"),            # MYC raises the effector
          g("Bcl2l1", "6W_pos") < g("Bcl2l1", "6W_neg"),      # and lowers BCL-XL at 6W
          g("Bbc3", "12W_pos") < g("Bbc3", "12W_neg"),        # PUMA ends below WT at 12W
          g("Bbc3", "12W_pos") < g("Bbc3", "6W_pos"),
          max(abs(d$z)) < 2)                                    # inside the colour limits

p <- ggplot(d, aes(group, gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = lab), size = pt(NM_TXT_MIN), family = NM_FONT, colour = "black") +
  facet_grid(role ~ ., scales = "free_y", space = "free_y") +
  scale_fill_div(c(-2, 2), name = "z score\n(group mean)", breaks = c(-2, -1, 0, 1, 2)) +
  scale_x_discrete(labels = group_labels, position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(barwidth = unit(2, "mm"), barheight = unit(18, "mm"),
                                ticks.colour = "grey30", frame.colour = "grey30")) +
  theme_nm(legend = "right") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(face = "italic"),
        axis.text.x.top = element_text(lineheight = 0.85),
        strip.text.y = element_text(size = NM_TXT, angle = 0, hjust = 0),
        panel.spacing.y = unit(0.8, "mm"),
        legend.title = element_text(size = NM_TXT, lineheight = 0.9))

save_panel(p, fig = "Fig1", panel = "g", name = "bcl2_heatmap", width = 80, height = 56)

cat("\n1g heatmap -- for the legend: VST expression, z-scored per gene over 24 animals; group means\n")
print(round(MZ[ord, ], 2))
