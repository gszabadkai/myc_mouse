# =============================================================================
# figS1_myc_transcript_stable.R -- the protein halves; the message does not
# -----------------------------------------------------------------------------
# SLOT: Fig. S1E.
#
#   "Mechanistically, this attenuation was driven by a ~50% reduction in MYC
#    protein levels at 12W compared to 6W (Fig. S1D). This decline occurred
#    despite stable transcript levels; the expression gap remained constant
#    between 6W_myc and 12W_myc for MYC and the proximal MYC/MAX/MXD network
#    (Fig. S1E, F)."
#
# S1D is the western blot, the author's bench data. This panel is the control the
# blot needs: if the message fell too, the attenuation would need no further
# explanation and the whole dose argument would collapse into "less transgene".
# It does not fall. Two things are drawn and both matter:
#
#   the GAP, Myc+ against wild type, is +1.68 at six weeks and +1.82 at twelve --
#     constant, if anything wider;
#   the LEVEL within the Myc+ gland does not decline across the window
#     (-0.10 log2, padj 0.76).
#
# Measured as the ratio of its twelve-week to its six-week effect, Myc is the
# single most retained gene in the transcriptome -- retention 1.085, the 99.99th
# percentile of 8,774 (script 44). The protein is what falls.
#
# GROUPS ARE DRAWN AGE-MAJOR here, the opposite of Fig. 1E, and the reason is the
# same in both cases: the order follows the TEST. Fig. 1E draws one pooled
# genotype main effect, so both wild-type boxes sit together and one bracket spans
# the contrast. Here the tests are the two WITHIN-AGE gaps, so each age's pair has
# to be adjacent.
#
# THE Y AXIS IS log2 OF NORMALISED COUNTS, labelled in counts. On that scale the
# claim "the gap remained constant" is a constant vertical distance and can be
# read off the panel; on a linear count axis it would not be.
#
# THE BRACKETS CARRY DESeq2, NOT A t-TEST ON THE DRAWN POINTS. The effects and
# adjusted p-values come from the raw (unshrunken) interaction model, which is the
# analysis of record and what the manuscript quotes; the script asserts that the
# empirical log2 gap between the drawn groups reproduces the DESeq2 log2 fold
# change to within 0.1, so the picture and the statistic are the same thing.
#
# Reads (read-only, no re-run):
#   results/dds_int_run.rds           -- normalised counts + colData
#   results/combined_df_annotated.rds -- mgi_symbol <-> ensembl
#   results/interaction_results.rds   -- the four raw contrasts (effects, padj)
#   results/collapse_module_ownership.rds -- $collapse_genes, for the retention
#                                            percentile quoted in the legend
# Output: outputs/figures/panels/figS1_myc_transcript_stable.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))
if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Fig. S1E needs DESeq2 for normalised counts")

# NOT gated by require_fresher_than(). dds_int_run.rds and interaction_results.rds
# are the Step 1 objects: they predate the 2026-07-24 gene-symbol reconciliation
# because the reconciliation is a SET-membership fix and never touched counts or
# per-gene DESeq2 results (docs/2026-07-24_symbol_reconciliation.md). The guard is
# for objects downstream of gsva_scores.rds; these are upstream of everything.
dds <- readRDS(here::here("results", "dds_int_run.rds"))
cdf <- readRDS(here::here("results", "combined_df_annotated.rds"))
ir  <- readRDS(here::here("results", "interaction_results.rds"))

nc <- DESeq2::counts(dds, normalized = TRUE)
sm <- as.data.frame(SummarizedExperiment::colData(dds))

ens <- cdf$gene[match("Myc", cdf$mgi_symbol)]
stopifnot(length(ens) == 1L, !is.na(ens), ens %in% rownames(nc))

# AGE-MAJOR: 6W_wt, 6W_myc | 12W_wt, 12W_myc -- so each age's pair is adjacent.
grp_levels <- c("6W_neg", "6W_pos", "12W_neg", "12W_pos")
dat <- data.frame(group = factor(as.character(sm$group), levels = grp_levels),
                  y = log2(nc[ens, ]))
stopifnot(nrow(dat) == 24L, !anyNA(dat$group), all(table(dat$group) == 6L))

# --- the four contrasts, from DESeq2 -----------------------------------------
slot_of <- c(myc_6W = "myc_6W_raw", myc_12W = "myc_12W_raw",
             `6>12W_wt` = "timepoint_neg_raw", `6>12W_myc` = "timepoint_pos_raw")
res <- do.call(rbind, lapply(names(slot_of), function(k) {
  r <- as.data.frame(ir[[slot_of[[k]]]])[ens, ]
  data.frame(contrast = k, lfc = r$log2FoldChange, se = r$lfcSE, padj = r$padj)
}))
rownames(res) <- res$contrast
int <- as.data.frame(ir[["interaction_raw"]])[ens, ]

# ASSERTION: the drawn groups reproduce the model. Group log2 means are of
# normalised counts and DESeq2 fits a negative-binomial GLM, so they agree to
# about a hundredth of a log2 unit rather than exactly.
gm  <- tapply(dat$y, dat$group, mean)
emp <- c(myc_6W    = unname(gm[["6W_pos"]]  - gm[["6W_neg"]]),
         myc_12W   = unname(gm[["12W_pos"]] - gm[["12W_neg"]]),
         `6>12W_wt`  = unname(gm[["12W_neg"]] - gm[["6W_neg"]]),
         `6>12W_myc` = unname(gm[["12W_pos"]] - gm[["6W_pos"]]))
stopifnot(max(abs(emp[res$contrast] - res$lfc)) < 0.1)

# --- brackets ----------------------------------------------------------------
# x positions follow the age-major order above: 6W_wt 1, 6W_myc 2, 12W_wt 3,
# 12W_myc 4. The two genotype gaps are adjacent pairs; the two trajectories span
# them and stack above.
fmt_p <- function(p) if (is.na(p)) "n.a." else if (p < 0.001) "<0.001" else
  formatC(p, format = "g", digits = 2)

comps <- data.frame(
  contrast = c("myc_6W", "myc_12W", "6>12W_myc", "6>12W_wt"),
  x1 = c(1, 3, 2, 1), x2 = c(2, 4, 4, 3), level = c(1, 1, 2, 3),
  stringsAsFactors = FALSE)
comps$lab <- vapply(res[comps$contrast, "padj"], fmt_p, character(1))
comps$col <- ifelse(!is.na(res[comps$contrast, "padj"]) &
                      res[comps$contrast, "padj"] < 0.05, "sig", "ns")
brk <- bracket_frame(comps, range(dat$y))

pts_layer <- if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
  ggbeeswarm::geom_quasirandom(width = 0.22, size = 1.0, alpha = 0.95)
} else {
  ggplot2::geom_jitter(width = 0.16, height = 0, size = 1.0, alpha = 0.95)
}

# counts on a log2 axis: the breaks are chosen in count units and placed by log2,
# so a constant genotype gap is a constant vertical distance on the page
brk_counts <- c(2000, 3000, 4000, 6000, 8000, 12000)

p <- ggplot2::ggplot(dat, ggplot2::aes(group, y, colour = group, fill = group)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.28,
                        colour = "grey35", linewidth = 0.25) +
  pts_layer +
  bracket_layers(brk) +
  ggplot2::geom_blank(inherit.aes = FALSE,
                      data = data.frame(group = factor("6W_neg", levels = grp_levels),
                                        y = attr(brk, "headroom")),
                      ggplot2::aes(x = group, y = y)) +
  # key order follows the DRAWN order, which is age-major here, so the reader can
  # match key to box left to right without re-sorting
  ggplot2::scale_colour_manual(values = c(group_cols, sig_cols),
                               breaks = grp_levels, labels = group_labels[grp_levels]) +
  ggplot2::scale_fill_manual(values = group_cols, guide = "none") +
  ggplot2::scale_y_continuous(breaks = log2(brk_counts),
                              labels = format(brk_counts, big.mark = ",", trim = TRUE),
                              expand = ggplot2::expansion(mult = c(0.06, 0.02))) +
  ggplot2::labs(x = NULL, y = "Myc, normalised counts (log2 scale)", colour = NULL) +
  ggplot2::guides(colour = ggplot2::guide_legend(
    nrow = 2, override.aes = list(size = 1.5, alpha = 1, shape = 16))) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    axis.text.x     = ggplot2::element_blank(),
    axis.ticks.x    = ggplot2::element_blank(),
    axis.line.x     = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.margin   = ggplot2::margin(-2, 0, 0, 0),
    legend.key.size = ggplot2::unit(2.4, "mm"),
    plot.margin     = ggplot2::margin(1.5, 1.5, 1, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
cg  <- as.data.frame(readRDS(here::here("results",
                                        "collapse_module_ownership.rds"))$collapse_genes)
myc <- cg[which(cg$gene == "Myc"), ]
f   <- function(k, w) res[k, w]

LEGEND <- panel_legend(
  slot = "Fig. S1E",
  what = paste0(
    "Myc transcript, normalised counts, one point per animal. The control the ",
    "western blot of Fig. S1D needs: the protein falls by about half across this ",
    "window and the message does not."),
  detail = c(
    "n = 6 animals per group, n = 24. Boxes are median and quartiles, whiskers 1.5x the interquartile range. Groups are ordered by age so that each age's genotype pair is adjacent, which is the comparison the two lower brackets make. The x axis is blank because the colour key carries the group names.",
    "The y axis is log2 of DESeq2 normalised counts, labelled in counts. On that scale a constant genotype gap is a constant vertical distance, which is the claim.",
    sprintf("Brackets carry the raw (unshrunken) DESeq2 contrasts and their Benjamini-Hochberg adjusted p-values, red where padj < 0.05: the genotype gap is %+.2f log2 at six weeks (padj %.2g) and %+.2f at twelve (padj %.2g); across age the wild-type gland changes by %+.2f (padj %.2f) and the Myc+ gland by %+.2f (padj %.2f). The interaction is %+.2f, padj %.2g.",
            f("myc_6W","lfc"), f("myc_6W","padj"), f("myc_12W","lfc"), f("myc_12W","padj"),
            f("6>12W_wt","lfc"), f("6>12W_wt","padj"),
            f("6>12W_myc","lfc"), f("6>12W_myc","padj"),
            int$log2FoldChange, int$padj),
    sprintf("The gap does not narrow; it widens slightly. Measured as the ratio of its twelve-week effect to its six-week one, Myc is the single most retained gene in the transcriptome -- retention %.3f, the %.3fth percentile of the %s genes script 44 ranks.",
            myc$retention, myc$pct_z, format(nrow(cg), big.mark = ",")),
    sprintf("The empirical log2 difference between the drawn group means reproduces each DESeq2 fold change to within %.3f, so the picture and the statistic are the same quantity.",
            max(abs(emp[res$contrast] - res$lfc)))),
  bounds = c(
    "THESE ARE READS ON THE MOUSE Myc LOCUS. Total Myc message in a Myc+ animal is transgene plus endogenous and the two cannot be separated here; the MMTV-Myc construct carries the mouse coding sequence. \"Dose is constant\" is therefore a statement about total Myc message, which is what the substrate sees, and script 27 has the endogenous-versus-transgene decomposition at the programme level.",
    "The two genotype brackets are clean -- genotype is balanced within each extraction batch. The two 6W-versus-12W brackets are not: batch is aligned with timepoint by design. Both of those are non-significant, and the panel's claim is that nothing falls, so the batch exposure works against the claim rather than for it.",
    "A stable message with a halved protein is consistent with translational or post-translational control and this panel cannot say which. Script 41 scanned the MYC post-translational stability machinery and found it is not transcriptionally down -- what moves in it predicts a MORE stable protein, not a less stable one -- so the mechanism of the protein fall is unresolved and is not claimed.",
    "The blot is the primary evidence for the dose statement. This panel is a negative control on the obvious alternative, not independent support: it rules out \"the transgene was silenced\" and nothing more.",
    sprintf("The manuscript sentence reads \"the expression gap remained constant between 6W_myc and 12W_myc\". A gap is between genotypes, not between ages. Both readings are true here and the panel draws both -- the genotype gap %+.2f to %+.2f, and no decline within the Myc+ gland (%+.2f, padj %.2f) -- but the sentence should say which one it means.",
            f("myc_6W","lfc"), f("myc_12W","lfc"),
            f("6>12W_myc","lfc"), f("6>12W_myc","padj"))),
  source = c(
    "results/dds_int_run.rds (scripts/03) -- DESeq2 normalised counts and colData",
    "results/interaction_results.rds (scripts/03) -- the raw, unshrunken contrasts",
    "results/combined_df_annotated.rds -- mgi_symbol to ensembl",
    "results/collapse_module_ownership.rds (scripts/44) -- $collapse_genes, the retention percentile",
    "Earlier double-column form: figures/figS8_myc_network_levels.R panel A"))

save_panel_p(p, "figS1_myc_transcript_stable", width = 55, height = 55)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the four contrasts and the interaction, as drawn
  rbind(res, data.frame(contrast = "interaction", lfc = int$log2FoldChange,
                        se = int$lfcSE, padj = int$padj)) |>
    print(row.names = FALSE, digits = 3)

  ## group means in counts, and the empirical gaps the assertion checks
  round(2^gm) |> print()
  round(emp, 3) |> print()

  ## Myc against the rest of the transcriptome on the departure-from-dose ranking
  cg[order(-cg$pct_z), c("gene", "lfc_6W", "lfc_12W", "retention", "pct_z")] |>
    head(8) |> print(row.names = FALSE, digits = 3)
}
