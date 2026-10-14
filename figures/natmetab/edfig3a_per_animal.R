# =============================================================================
# edfig3a_per_animal.R -- Extended Data Fig. 3a: OXPHOS priority, the PUMA:BCL-XL
# balance, BAX and BMF in every animal
# -----------------------------------------------------------------------------
# PROMOTED from proposals/ on 2026-10-14 (author: "promote them to final panels").
# The author (2026-10-14), on the five-gene version (p1h_puma_members.R): "This
# kind of just repeats the heatmap data. I would rather use the per animal z values
# and plot them together as in Fig 2H+I (alt): OXPHOS mitoPPS, PUMA:BCL-XL and Bax,
# Bmf values." The form is Fig. 2H+I (alt) part A (figures/panels/
# fig2_puma_chain_alt.R). Each quantity is standardised over the 24 animals, so
# four quantities in three units share one axis. Each group shows a box and every
# animal, genotype-major, so each genotype's 6W -> 12W pair sits side by side.
#
#   OXPHOS (mitoPPS)  script 42's axis score = script 08's OXPHOS-subunit mitoPPS
#   Bbc3:Bcl2l1       log2(PUMA) - log2(BCL-XL), script 42's ratio (log2 of
#                     normalised counts + 1), the ratio of Fig. 1h-i
#   Bax, Bmf          log2 normalised counts + 1, the same scale; Bmf is the
#                     developmental control (it rises in both genotypes)
#
# BRACKETS: lower = 12W vs 6W within each genotype; upper = the difference between
# the two changes (age x genotype interaction), the one contrast the batch design
# leaves clean. RAW P throughout (script 54's rule), all from script 56 PART B:
#   Bax, Bmf          DESeq2 Wald (script 03), the record the text quotes
#   the two composites  lm(value ~ age * genotype), script 42's model; for the
#                     ratio its interaction IS script 42's saved value
# Script 08's Welch P for mitoPPS (Fig. 1f's test) differs a little from the model's
# (pooled against separate variances): the legend names the one drawn, and this
# script prints both.
# NOT LICENSED: "OXPHOS priority falls on both timelines". The WT window falls by as
# much but does not clear 0.05 (script 56 asserts it).
#
# Reads: results/window_and_chain_statistics.rds (script 56: $values, $tests, $welch_08)
# Output: outputs/natmetab/EDFig3/EDFig3a_per_animal.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

s56 <- readRDS(here::here("results", "window_and_chain_statistics.rds"))
val <- as.data.frame(s56$values)
tst <- as.data.frame(s56$tests)
grp <- factor(val$group, levels = order_by_geno)
stopifnot(nrow(val) == 24L, !anyNA(grp), all(table(grp) == 6L))

Q     <- c("oxphos_mitopps", "puma_bclxl", "Bax", "Bmf")
STRIP <- c(oxphos_mitopps = 'OXPHOS~"(mitoPPS)"', puma_bclxl = 'italic("Bbc3:Bcl2l1")',
           Bax = 'italic("Bax")', Bmf = 'italic("Bmf")')
tget  <- function(q, ct, col = "p") tst[[col]][tst$quantity == q & tst$contrast == ct]
stopifnot(all(Q %in% names(val)), all(Q %in% tst$quantity))

# the drawn values move the way their tests say (group means against the estimates)
gmean <- function(v, g) mean(v[grp == g])
for (q in Q) {
  v <- val[[q]]
  stopifnot(sign(gmean(v, "12W_neg") - gmean(v, "6W_neg")) == sign(tget(q, "wt_12v6", "estimate")),
            sign(gmean(v, "12W_pos") - gmean(v, "6W_pos")) == sign(tget(q, "myc_12v6", "estimate")))
}

# --- the drawn values: z over the 24 animals ----------------------------------------------
d <- do.call(rbind, lapply(Q, function(q) data.frame(q = q, group = grp, z = as.numeric(scale(val[[q]])))))
d$q <- factor(d$q, levels = Q, labels = STRIP[Q])

# brackets: the lower pair at one base in every facet; staggered where a P needs an exponent
STEP <- 0.62
BASE <- max(d$z) + 0.3
pos  <- function(g) match(g, order_by_geno)
br <- do.call(rbind, lapply(Q, function(q) {
  p    <- c(tget(q, "wt_12v6"), tget(q, "myc_12v6"), tget(q, "interaction"))
  wide <- any(p[1:2] < 0.001)
  data.frame(q = factor(STRIP[[q]], levels = levels(d$q)),
             group1 = c(pos("6W_neg"), pos("6W_pos"), 1.5), group2 = c(pos("12W_neg"), pos("12W_pos"), 3.5),
             p = p, level = if (wide) c(1, 2, 3) else c(1, 1, 2))
}))
bl  <- brackets(br, base = rep(BASE, nrow(br)), step = STEP)
TOP <- max(attr(bl, "tops"))

p <- ggplot(d, aes(group, z)) +
  geom_hline(yintercept = 0, linewidth = NM_LINE, colour = "grey85") +
  geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.62, colour = box_line,
               linewidth = NM_LINE) +
  ggbeeswarm::geom_quasirandom(aes(colour = group), width = 0.18, size = 0.7, shape = 16) +
  bl +
  facet_wrap(~ q, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = group_cols, labels = group_keys, name = NULL) +
  scale_fill_manual(values = group_fill, guide = "none") +
  scale_y_nm(step = 1, top = TOP, bottom = min(d$z) - 0.1, minor = TRUE, labels = lab_signed) +
  labs(x = NULL, y = "Per animal (z score)") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.5))) +
  coord_cartesian(clip = "off") +
  theme_nm(legend = "bottom") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        # 6 pt and unclipped: "OXPHOS (mitoPPS)" is wider than a quarter of the panel at 7 pt
        strip.text = element_text(size = NM_TXT, margin = margin(0, 0, 0.5, 0, "mm")),
        strip.clip = "off",
        panel.spacing.x = unit(2, "mm"), legend.key.spacing.x = unit(1, "mm"))

save_panel(p, fig = "EDFig3", panel = "a", name = "per_animal", width = 104, height = 60)

cat("\nED 3a -- for the legend: z over the 24 animals; RAW P (script 56); lower brackets 12W vs 6W within a\n")
cat("  genotype, upper bracket the interaction (the difference between the two changes)\n")
for (q in Q) cat(sprintf("  %-15s %-58s WT %+.3f P %.2g | MYC+ %+.3f P %.2g | interaction %+.3f P %.2g\n", q,
  unique(tst$instrument[tst$quantity == q]), tget(q, "wt_12v6", "estimate"), tget(q, "wt_12v6"),
  tget(q, "myc_12v6", "estimate"), tget(q, "myc_12v6"), tget(q, "interaction", "estimate"), tget(q, "interaction")))
w8 <- as.data.frame(s56$welch_08)
cat(sprintf("  OXPHOS mitoPPS by script 08's Welch test (Fig. 1f): WT P %.2g (BH %.2g), MYC+ P %.2g (BH %.2g)\n",
            w8$p_welch[1], w8$padj_welch_bh[1], w8$p_welch[2], w8$padj_welch_bh[2]))
cat("  estimates: the composites in their own units (mitoPPS; log2 ratio), the genes in log2 fold change\n")
