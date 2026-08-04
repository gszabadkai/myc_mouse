# =============================================================================
# figS1_myc_network.R -- the proximal MYC/MAX/MXD network does not move either
# -----------------------------------------------------------------------------
# SLOT: Fig. S1F, the second half of the sentence Fig. S1E opens.
#
#   "This decline occurred despite stable transcript levels; the expression gap
#    remained constant between 6W_myc and 12W_myc for MYC and the proximal
#    MYC/MAX/MXD network (Fig. S1E, F)."
#
# S1E is Myc itself. This is the network around it, and it is a NEGATIVE CONTROL
# on the most obvious alternative to a dose explanation: MYC works as a
# heterodimer with MAX and competes for the same E-boxes with the MXD/MNT family,
# so a rise in the competing repressors, or a fall in MAX, would attenuate Myc's
# output at constant Myc. Neither happens. Nothing in the network moves with
# genotype at either age, and nothing changes between the ages.
#
# THE PANEL DRAWS THE CONTRAST, NOT THE LEVEL, which is the difference from S1E.
# The sentence's claim is about the GAP, so the gap is the quantity on the axis:
# each gene's Myc genotype effect at six weeks and at twelve, joined. A level plot
# of ten genes would need ten facets and 183 mm to say the same thing less
# directly -- that form is figures/figS8_myc_network_levels.R panel B.
#
# WHY THE ERROR BAR IS ON THE SIX-WEEK POINT ONLY. This is a negative, so the
# reader has to be able to tell "measured to be zero" from "unmeasurable", which
# is what the bar is for. Drawing it twice would double the ink for nothing: the
# two ages' standard errors are the same to within 0.024 log2 across all ten
# genes (median ratio 1.01). That equality is worth knowing on its own account --
# twelve weeks is not the noisier cohort, which is what protects every attenuation
# result in the paper from the obvious objection.
#
# NO SIGNIFICANCE ENCODING, because nothing is significant and saying so once in
# the legend is more honest than a key whose filled entry never appears: no gene
# clears padj 0.05 at either age (closest Mxi1 at 0.096 and Mlx at 0.156), and no
# interaction is below padj 0.83.
#
# Reads (read-only, no re-run):
#   results/interaction_results.rds   -- the raw, unshrunken contrasts
#   results/combined_df_annotated.rds -- mgi_symbol <-> ensembl
# Output: outputs/figures/panels/figS1_myc_network.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

# Not gated by require_fresher_than(): these are the Step 1 objects, upstream of
# the 2026-07-24 gene-symbol reconciliation, which was a set-membership fix and
# never touched per-gene DESeq2 results. Same reasoning as Fig. S1E.
ir  <- readRDS(here::here("results", "interaction_results.rds"))
cdf <- readRDS(here::here("results", "combined_df_annotated.rds"))

# --- the roster, by role -----------------------------------------------------
# MXI1 is MXD2, so the MXD family is complete. MLX heterodimerises with MLXIP
# (MondoA) and MLXIPL (ChREBP) rather than with MYC, and is included because it
# competes for MAX-family partners -- the arm script 31's scan left open.
net <- data.frame(
  gene = c("Max",
           "Mxd1", "Mxd3", "Mxd4", "Mxi1", "Mnt", "Mga",
           "Mlx", "Mlxip", "Mlxipl"),
  role = c("obligate partner",
           rep("MXD / MNT repressors", 6),
           rep("MLX arm", 3)),
  stringsAsFactors = FALSE)
net$role <- factor(net$role,
                   levels = c("obligate partner", "MXD / MNT repressors", "MLX arm"))

ens <- cdf$gene[match(net$gene, cdf$mgi_symbol)]
stopifnot(!anyNA(ens), all(ens %in% rownames(as.data.frame(ir[["myc_6W_raw"]]))))

grab <- function(slot) {
  r <- as.data.frame(ir[[slot]])[ens, ]
  data.frame(lfc = r$log2FoldChange, se = r$lfcSE, padj = r$padj,
             baseMean = r$baseMean)
}
g6  <- grab("myc_6W_raw")
g12 <- grab("myc_12W_raw")
gi  <- grab("interaction_raw")
gtp <- grab("timepoint_pos_raw")

d <- cbind(net, m6 = g6$lfc, se6 = g6$se, p6 = g6$padj, baseMean = g6$baseMean,
           m12 = g12$lfc, se12 = g12$se, p12 = g12$padj,
           int = gi$lfc, pint = gi$padj, tp = gtp$lfc, ptp = gtp$padj)

# ASSERTION: the two ages are equally precise, which is what licenses drawing one
# error bar and is a claim in its own right.
SE_TOL <- 0.05
stopifnot(nrow(d) == 10L, max(abs(d$se12 - d$se6)) < SE_TOL)

# rows ordered by the six-week effect within each role block
d <- d[order(d$role, d$m6), ]
d$gene <- factor(d$gene, levels = d$gene)

dl <- rbind(
  data.frame(gene = d$gene, role = d$role, eff = d$m6,  contrast = "myc_6W"),
  data.frame(gene = d$gene, role = d$role, eff = d$m12, contrast = "myc_12W"))
dl$contrast <- factor(dl$contrast, levels = contrast_geno)

p <- ggplot2::ggplot(dl, ggplot2::aes(eff, gene)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  # Connector first, then the error bar OVER it in the six-week colour: the bar
  # belongs to the six-week estimate and has to read that way, or a twelve-week
  # point falling outside it looks like a significant difference between the ages.
  ggplot2::geom_segment(data = d, inherit.aes = FALSE,
                        ggplot2::aes(x = m6, xend = m12, y = gene, yend = gene),
                        colour = "grey45", linewidth = 0.25) +
  ggplot2::geom_segment(data = d, inherit.aes = FALSE,
                        ggplot2::aes(x = m6 - se6, xend = m6 + se6,
                                     y = gene, yend = gene),
                        colour = contrast_cols[["myc_6W"]], alpha = 0.4,
                        linewidth = 0.9, lineend = "round") +
  ggplot2::geom_point(data = dl[dl$contrast == "myc_12W", ],
                      ggplot2::aes(colour = contrast), size = 0.9) +
  ggplot2::geom_point(data = dl[dl$contrast == "myc_6W", ],
                      ggplot2::aes(colour = contrast), size = 1.5) +
  ggplot2::facet_grid(role ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_colour_manual(values = contrast_cols[contrast_geno],
                               breaks = contrast_geno, name = NULL) +
  ggplot2::scale_x_continuous(labels = lab_signed) +
  ggplot2::labs(x = "Myc genotype effect (log2 fold change)", y = NULL) +
  ggplot2::guides(colour = ggplot2::guide_legend(
    override.aes = list(size = 1.5))) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    axis.text.y     = ggplot2::element_text(size = 5.8, face = "italic"),
    axis.line.y     = ggplot2::element_blank(),
    axis.ticks.y    = ggplot2::element_blank(),
    strip.text.y    = ggplot2::element_text(size = 5.2, angle = 0, hjust = 0,
                                            margin = ggplot2::margin(0, 0, 0, 1, "mm")),
    panel.spacing.y = ggplot2::unit(1.2, "mm"),
    # Under the plot, not inside: unlike Figs. 1G and 1H this panel has no empty
    # corner -- the MLX rows run to the right edge and a key there sits on them.
    legend.position = "bottom",
    legend.margin   = ggplot2::margin(-2, 0, 0, 0),
    legend.key.size = ggplot2::unit(2.4, "mm"),
    plot.margin     = ggplot2::margin(1.5, 1.5, 1, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
worst6  <- d[order(d$p6), ][1:2, ]
lowexp  <- d[d$baseMean < 100, ]
movers  <- d[!is.na(d$ptp) & d$ptp < 0.05, ]

LEGEND <- panel_legend(
  slot = "Fig. S1F",
  what = paste0(
    "The proximal MYC/MAX/MXD network. For each gene, the Myc genotype effect at ",
    "six weeks and at twelve, joined; the soft band is one standard error on the ",
    "six-week estimate. The claim the panel supports is a negative: nothing in ",
    "the network that could attenuate Myc's output at constant Myc changes."),
  detail = c(
    "n = 6 animals per group. Values are raw, unshrunken DESeq2 log2 fold changes of the genotype contrast (Myc+ minus wild type) within each age; both contrasts are clean, genotype being balanced within each extraction batch.",
    "MAX is the obligate heterodimerisation partner: MYC cannot bind an E-box without it, so a fall in MAX would attenuate Myc's output at constant Myc. MXD1, MXD3, MXD4, MXI1 (which is MXD2, so the family is complete), MNT and MGA are the competing repressors that dimerise with MAX at the same E-boxes; a rise in them would do the same. MLX, MLXIP (MondoA) and MLXIPL (ChREBP) are the parallel MAX-family arm.",
    sprintf("NOTHING CLEARS SIGNIFICANCE, which is why the panel carries no significance encoding: the two closest at six weeks are %s (%+.2f, padj %.3f) and %s (%+.2f, padj %.2f), and every gene is further from significance at twelve weeks. No interaction is below padj %.2f.",
            worst6$gene[1], worst6$m6[1], worst6$p6[1],
            worst6$gene[2], worst6$m6[2], worst6$p6[2], min(d$pint, na.rm = TRUE)),
    sprintf("The error bar is drawn on the six-week estimate only because the two ages are equally precise: standard errors agree to within %.3f log2 across all ten genes, median ratio %.2f. Twelve weeks is not the noisier cohort, which is the objection every attenuation result in the paper has to survive.",
            max(abs(d$se12 - d$se6)), stats::median(d$se12 / d$se6)),
    sprintf("Across age within the Myc+ gland the only gene that moves at all is %s (%+.2f, padj %.3f), and it is one of the two lowest-expressed on the panel.",
            paste(movers$gene, collapse = ", "), movers$tp[1], movers$ptp[1]),
    "This is the gene-level counterpart of the scan script 31 ran at set level, which found every E-box binder flat and no significant interaction among 111 nuclear receptors, corepressors and DREAM components. The two agree."),
  bounds = c(
    sprintf("TWO GENES ARE TOO LOWLY EXPRESSED TO CARRY A NEGATIVE. %s, with base means of %s, have standard errors of %s -- wide enough that a real effect of half a log2 unit would not be detected. Their flatness is uninformative and the panel shows the bars so that this is visible rather than assumed.",
            paste(lowexp$gene, collapse = " and "),
            paste(round(lowexp$baseMean), collapse = " and "),
            paste(sprintf("%.2f", lowexp$se6), collapse = " and ")),
    "A NEGATIVE AT THE MESSAGE LEVEL IS NOT A NEGATIVE AT THE PROTEIN LEVEL. MAX availability, MXD/MNT competition and MYC itself are all controlled post-translationally; this panel rules out a transcriptional shift in the network and nothing more. The same limit applies to Fig. S1E.",
    "The panel is a control, not evidence for the dose model. What supports that model is the western blot of Fig. S1D; what this panel does is close the most obvious alternative to it.",
    "The two genotype contrasts drawn are clean. The temporal comparison quoted in the detail is exposed to batch = timepoint, as every wild-type or within-genotype temporal statement in this paper is.",
    "MLX, MLXIP and MLXIPL are drawn for completeness of the MAX-family arm; they are not MYC partners, and MLXIP's route to the Myc programme -- competition for MLX release -- was left untestable by script 31."),
  source = c(
    "results/interaction_results.rds (scripts/03) -- the raw, unshrunken contrasts",
    "results/combined_df_annotated.rds -- mgi_symbol to ensembl",
    "The set-level counterpart: scripts/31, the E-box binder and corepressor scan",
    "Earlier double-column form: figures/figS8_myc_network_levels.R panel B, which draws the LEVELS in ten facets"))

save_panel_p(p, "figS1_myc_network", height = 55)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## every gene, every contrast
  d[, c("gene", "role", "baseMean", "m6", "se6", "p6", "m12", "p12",
        "tp", "ptp", "int", "pint")] |> print(row.names = FALSE, digits = 3)

  ## the two ages' standard errors, which is what licenses one error bar
  data.frame(gene = d$gene, se6 = d$se6, se12 = d$se12,
             diff = d$se12 - d$se6, ratio = d$se12 / d$se6) |>
    print(row.names = FALSE, digits = 3)

  ## the LEVELS form, ten facets -- the double-column version this replaces
  ## source(here::here("figures", "figS8_myc_network_levels.R"))
}
