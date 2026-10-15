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
# TYPE (2026-10-15, as Fig. 1d and 1e): gene and group labels, tile numbers and the
# colour bar's labels 5 pt, the role names and the colour bar's title 6 pt.
#
# SIZE (author, 2026-10-15): g, h and i share one row of 180 mm, so g keeps its 56-mm
# height and narrows: the key's title is "z-score" alone, the bar has no frame, and the
# key moves left, to the top of the role names' column beside the two short names
# (Effectors, Activators) rather than past the longest (Anti-apoptotic).
#
# Reads: results/gsva_scores.rds (script 15: $expr_mat, VST; $sample_meta)
#        results/departure_from_myc_fade.rds (script 55: the roster and roles)
# Output: outputs/natmetab/Fig1/Fig1g_bcl2_heatmap.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

gs  <- readRDS(here::here("results", "gsva_scores.rds"))
dep <- readRDS(here::here("results", "departure_from_myc_fade.rds"))

TXT   <- 5      # pt; gene and group labels, tile numbers, the colour bar's labels
TITLE <- 6      # pt; the role names and the colour bar's title
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
  geom_text(aes(label = lab), size = pt(TXT), family = NM_FONT, colour = "black") +
  facet_grid(role ~ ., scales = "free_y", space = "free_y") +
  scale_fill_div(c(-2, 2), name = "z­score", breaks = c(-2, -1, 0, 1, 2)) +   # soft hyphen: pdf() prints "-" as a minus
  scale_x_discrete(labels = group_labels, position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(barwidth = unit(2, "mm"), barheight = unit(18, "mm"),
                                ticks.colour = "grey30")) +                  # no frame (author)
  theme_nm(legend = "right", txt = TXT, txt_title = TITLE) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(face = "italic"),
        axis.text.x.top = element_text(lineheight = 0.85),
        strip.text.y = element_text(size = TITLE, angle = 0, hjust = 0),
        panel.spacing.y = unit(0.8, "mm"),
        legend.title = element_text(size = TITLE, lineheight = 0.9),
        legend.justification.right = "top")

# The key starts 1 mm past the wider of the two names beside it: from the role column's
# right edge (its longest name plus the strip's margin), that is a negative spacing.
# Widths at TITLE pt with the PDF device's Helvetica metrics. The 18-mm bar ends above
# the third name, Sensitisers, which is wider (checked in the PDF, 2026-10-15).
key_spacing <- local({
  grDevices::pdf(NULL, family = NM_FONT)
  on.exit(grDevices::dev.off())
  mm <- function(u) grid::convertWidth(u, "mm", valueOnly = TRUE)
  w  <- vapply(ROLES, function(s) mm(grid::grobWidth(grid::textGrob(s,
          gp = grid::gpar(fontsize = TITLE, fontfamily = NM_FONT)))), numeric(1))
  m  <- ggplot2::calc_element("strip.text.y.right", ggplot2::complete_theme(p$theme))$margin
  max(w[c("effector", "activator")]) + 1 - max(w) - mm(m[2])
})
p <- p + theme(legend.box.spacing = unit(key_spacing, "mm"))

save_panel(p, fig = "Fig1", panel = "g", name = "bcl2_heatmap", width = 66, height = 56)

cat("\n1g heatmap -- for the legend: VST expression, z-scored per gene over 24 animals; group means\n")
print(round(MZ[ord, ], 2))
