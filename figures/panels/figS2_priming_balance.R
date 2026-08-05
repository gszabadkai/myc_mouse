# =============================================================================
# figS2_priming_balance.R -- the death machinery is not what the window changes
# -----------------------------------------------------------------------------
# SLOT: Fig. S2C.
#
#   "The overall apoptotic priming remains stable in the WT timeline."
#
# THE JOB, AND WHY IT IS A NEGATIVE WORTH DRAWING. Section 2's argument is that
# the wild-type gland changes its respiratory STATE between six and twelve weeks
# (Fig. 2F) while the apparatus that executes death does not. If the machinery
# moved too, the death phenotype would need no further explanation and the whole
# competence argument would collapse into "the gland dismantled its apoptosome".
# It does not: of 37 transcripts, one moves, and it moves the wrong way for that
# reading -- Bnip3 RISES.
#
# THE ROSTER is figures/fig04_substrate_specificity.R panel D's: the 25 pro- and
# 7 anti-apoptotic MitoCarta transcripts, plus five more anti-apoptotic genes that
# MitoCarta does not contain and the apoptosis sets therefore miss -- the caspase
# inhibitors XIAP, cIAP1/2 (Birc2/3) and survivin (Birc5), and a Bcl2a1 paralog.
# They get their own row because they answer a different question: "does the gland
# BUFFER against death?" rather than "does it move its apoptotic transcripts?".
# The row is labelled for what they ARE (IAPs & Bcl2a1) rather than for how they
# were excluded, and the script asserts the membership so the label stays true.
#
# "OVERALL PRIMING" IS A BALANCE, so the composite is quoted in the legend on both
# rulers of Fig. 2F: across the window the pro arm moves +0.025 and the anti arm
# +0.045 on content, +0.054 and +0.046 on priority. The two arms move TOGETHER,
# so the difference -- which is what priming means -- is -0.020 and +0.008. For
# comparison the Myc genotype effect at six weeks moves them APART (+0.156 pro,
# -0.024 anti). The balance is something Myc changes and the window does not.
#
# THE POWER CONTROL, because a negative at n = 6 needs one: the same transcripts,
# in the same libraries, at the same n, DO move under the genotype contrast. The
# count is computed below rather than asserted.
#
# ONE TRAP THIS PANEL DOES NOT FALL INTO. `MITOCARTA_APOPTOSIS_PRO`/`_ANTI` are
# MitoCarta sets, so any coupling between "priming" and a mitochondrial axis is
# mito-vs-mito and circular (script 34). THAT WARNING DOES NOT APPLY HERE: this
# is a temporal contrast on the transcripts themselves, not a coupling, and
# nothing on this panel is correlated with a mitochondrial score.
#
# BATCH = TIMEPOINT (CLAUDE.md): the wild-type temporal axis is confounded with
# extraction batch, so this reads as "no detectable movement at n = 6 on a
# confounded axis", not as "no movement".
#
# Reads (read-only, no re-run):
#   results/collapse_module_ownership.rds      (script 44) -- $wt_genes, the 32
#                                                 MitoCarta pro/anti transcripts
#   results/substrate_specificity_tradeoff.rds (script 43) -- $buffer, the brakes
#   results/background_vs_myc.rds              (script 40) -- $ruler, the arm
#                                                 composites on both rulers
#   results/combined_df_annotated_raw.rds      (script 03) -- symbol -> Ensembl
#   results/interaction_results.rds            (script 03) -- padj per contrast
# Output: outputs/figures/panels/figS2_priming_balance.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))
if (!requireNamespace("DESeq2", quietly = TRUE))
  stop("figS2C needs DESeq2 to coerce the DESeqResults in interaction_results.rds")

cmo <- readRDS(here::here("results", "collapse_module_ownership.rds"))
ss  <- readRDS(here::here("results", "substrate_specificity_tradeoff.rds"))
bv  <- readRDS(here::here("results", "background_vs_myc.rds"))
cdf <- readRDS(here::here("results", "combined_df_annotated_raw.rds"))
ir  <- readRDS(here::here("results", "interaction_results.rds"))

# =============================================================================
# the roster -- fig04 panel D's, rebuilt from the same two objects
# =============================================================================
wg <- as.data.frame(cmo$wt_genes)
bf <- as.data.frame(ss$buffer)
stopifnot(all(c("gene", "arm", "wt_time", "padj_wt", "myc_6W") %in% names(wg)))

d <- data.frame(
  gene = wg$gene, baseMean = wg$baseMean,
  class = ifelse(wg$arm == "PRO", "pro-apoptotic", "anti-apoptotic"),
  lfc = wg$wt_time, padj = wg$padj_wt, lfc_myc = wg$myc_6W,
  stringsAsFactors = FALSE)
# The third row is the OTHER anti-apoptotic arm: the caspase inhibitors (XIAP and
# the cIAPs, plus survivin) and a Bcl2a1 paralog. They are not MitoCarta members,
# so they are absent from the apoptosis sets the first two rows come from -- and
# they answer a different question, "does the gland BUFFER?", rather than "does it
# move its apoptotic transcripts?". fig04's rule, verbatim: whatever of $buffer is
# not already present.
#
# THE ROW IS NAMED FOR ITS MEMBERS, so the members are checked. If $buffer ever
# gains a gene that is neither an IAP nor a Bcl2a1 paralog, the label stops being
# true and this stops the panel rather than mislabelling it.
add <- bf[!bf$gene %in% d$gene, ]
stopifnot(setequal(add$gene, c("Xiap", "Birc2", "Birc3", "Birc5", "Bcl2a1b")))
d <- rbind(d, data.frame(
  gene = add$gene, baseMean = add$baseMean, class = "IAPs & Bcl2a1",
  lfc = add$lfc_wt_time, padj = add$padj_wt_time, lfc_myc = add$lfc_myc_6W,
  stringsAsFactors = FALSE))
d <- d[is.finite(d$lfc), ]

n_class <- table(d$class)
stopifnot(nrow(d) == 37L, n_class[["pro-apoptotic"]] == 25L,
          n_class[["anti-apoptotic"]] == 7L,
          n_class[["IAPs & Bcl2a1"]] == 5L)

# Two of the 32 apoptosis-set members have no current MGI symbol in the
# annotation table; they stay in the distribution and simply go unlabelled.
d$gene[is.na(d$gene)] <- ""
d$sig <- !is.na(d$padj) & d$padj < 0.05

# =============================================================================
# the power control -- the same transcripts under the genotype contrast
# =============================================================================
# A negative at n = 6 is only worth drawing if the same measurement can detect
# movement. Adjusted p-values for the genotype contrast are keyed by Ensembl, so
# the symbols are mapped through the annotation table and the mapping is CHECKED
# against baseMean rather than trusted.
ann <- as.data.frame(cdf)
ens <- ann$gene[match(d$gene, ann$mgi_symbol)]
bm_ann <- ann$baseMean[match(d$gene, ann$mgi_symbol)]
mapped <- !is.na(ens) & nzchar(d$gene)
stopifnot(sum(mapped) >= 34L,
          max(abs(bm_ann[mapped] - d$baseMean[mapped])) < 1)   # same gene, same object

m6 <- as.data.frame(ir$myc_6W_raw)
d$padj_myc <- NA_real_
d$padj_myc[mapped] <- m6$padj[match(ens[mapped], rownames(m6))]
n_sig_myc <- sum(!is.na(d$padj_myc) & d$padj_myc < 0.05)
n_sig_wt  <- sum(d$sig)

# =============================================================================
# the composite -- "overall priming" is a balance, and it is on script 40's ruler
# =============================================================================
r <- as.data.frame(bv$ruler)
arm_of <- function(pw, col) {
  v <- r[[col]][r$pathway == pw]
  as.numeric(unname(v))
}
comp <- data.frame(
  arm      = c("Apoptosis-PRO", "Apoptosis-ANTI", "Apoptosis"),
  content  = vapply(c("Apoptosis-PRO", "Apoptosis-ANTI", "Apoptosis"),
                    arm_of, numeric(1), "c_tn"),
  priority = vapply(c("Apoptosis-PRO", "Apoptosis-ANTI", "Apoptosis"),
                    arm_of, numeric(1), "p_tn"),
  myc_content = vapply(c("Apoptosis-PRO", "Apoptosis-ANTI", "Apoptosis"),
                       arm_of, numeric(1), "c_m6"),
  row.names = NULL, stringsAsFactors = FALSE)
prime_wt  <- comp$content[1]     - comp$content[2]
prime_p   <- comp$priority[1]    - comp$priority[2]
prime_myc <- comp$myc_content[1] - comp$myc_content[2]
stopifnot(nrow(comp) == 3L, !anyNA(comp$content))

# =============================================================================
# the panel
# =============================================================================
LEV <- c("IAPs & Bcl2a1", "anti-apoptotic", "pro-apoptotic")
d$row <- factor(sprintf("%s (%d)", d$class, n_class[d$class]),
                levels = sprintf("%s (%d)", LEV, n_class[LEV]))

# DETERMINISTIC VERTICAL OFFSETS, not jitter. Within a row the genes are ordered
# by fold change and the offsets cycle through a fixed ladder, so neighbours in x
# are separated in y by construction -- better than random jitter at avoiding
# overplot, reproducible without a seed, and it lets the labels below be placed at
# coordinates that are KNOWN rather than guessed at from a render.
LADDER <- c(0, 0.20, -0.20, 0.10, -0.10, 0.30, -0.30)
d <- d[order(d$row, d$lfc), ]
d$dy <- unlist(lapply(split(seq_len(nrow(d)), d$row),
                      function(i) LADDER[(seq_along(i) - 1L) %% length(LADDER) + 1L]),
               use.names = FALSE)
d$y <- as.integer(d$row) + d$dy

# The three the later panels turn on: the one mover, PUMA (Fig. 2H) and Bcl-xL,
# the denominator of the ratio Fig. 2G reports.
# LABELS ARE PLACED, and every leader leaves the text on the text's own line --
# the convention Fig. 2F settled. Two go into a lane above the top row, which the
# y limit reserves; Bcl-xL goes out to the left, into the half of the anti row
# that is empty (nothing in it falls below -0.15). The regions are checked below
# rather than eyeballed.
SHOW <- c("Bnip3", "Bbc3", "Bcl2l1")
lab <- d[match(SHOW, d$gene), ]
stopifnot(!anyNA(lab$lfc))
TOP_LANE <- length(levels(d$row)) + 0.55
LEFT_X   <- -0.30
lab$mode  <- c("top", "top", "left")
lab$x_lab <- ifelse(lab$mode == "top", lab$lfc, LEFT_X)
lab$y_lab <- ifelse(lab$mode == "top", TOP_LANE, lab$y)
# where the leader leaves the text: just under it for the top lane, just past its
# right-hand end for the left-hand one
lab$x0 <- ifelse(lab$mode == "top", lab$x_lab, lab$x_lab + 0.012)
lab$y0 <- ifelse(lab$mode == "top", lab$y_lab - 0.05, lab$y_lab)
stopifnot(
  # the two top-lane labels are far enough apart in x not to touch
  abs(diff(lab$lfc[lab$mode == "top"])) > 0.3,
  # and the left-hand lane is empty of data at that height
  !any(d$lfc < LEFT_X + 0.06 & abs(d$y - lab$y[lab$mode == "left"]) < 0.35))

XR <- range(d$lfc) + c(-1, 1) * diff(range(d$lfc)) * 0.10

p <- ggplot2::ggplot(d, ggplot2::aes(lfc, y)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey45") +
  ggplot2::geom_point(ggplot2::aes(shape = sig), size = 1.3, stroke = 0.32,
                      colour = "grey20", fill = "grey20") +
  ggplot2::geom_segment(data = lab, inherit.aes = FALSE,
                        ggplot2::aes(x = x0, xend = lfc, y = y0, yend = y),
                        linewidth = 0.2, colour = "grey55") +
  ggplot2::geom_text(data = lab, inherit.aes = FALSE,
                     ggplot2::aes(x = x_lab, y = y_lab, label = gene,
                                  hjust = ifelse(mode == "top", 0.5, 1),
                                  vjust = ifelse(mode == "top", 0, 0.5)),
                     size = 1.7, colour = "grey15") +
  ggplot2::scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1),
                              breaks = c(TRUE, FALSE),
                              labels = c("padj < 0.05", "n.s."), name = NULL) +
  ggplot2::scale_x_continuous(limits = XR, labels = lab_signed,
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::scale_y_continuous(breaks = seq_along(levels(d$row)),
                              labels = levels(d$row),
                              limits = c(0.4, length(levels(d$row)) + 0.95),
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::labs(x = "wild-type 6>12W  (raw log2FC)", y = NULL) +
  ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 1.4))) +
  theme_panel(base_size = 6) +
  # Key inside, top left: every transcript that moves at all moves right, so the
  # left of the top row is the empty corner.
  ggplot2::theme(
    legend.position        = "inside",
    legend.position.inside = c(0.005, 0.99),
    legend.justification   = c(0, 1),
    legend.background      = ggplot2::element_blank(),
    legend.margin          = ggplot2::margin(0, 0, 0, 0),
    legend.key.size        = ggplot2::unit(2.4, "mm"),
    legend.spacing.y       = ggplot2::unit(0.3, "mm"),
    axis.text.y            = ggplot2::element_text(size = 6),
    axis.ticks.y           = ggplot2::element_blank(),
    plot.margin            = ggplot2::margin(1.5, 2.5, 1, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
one <- d[d$sig, ]
LEGEND <- panel_legend(
  slot = "Fig. S2C",
  what = paste0(
    "Every transcript of the mitochondrial death apparatus across the wild-type ",
    "6 to 12 week window: the 25 pro- and 7 anti-apoptotic MitoCarta genes, plus ",
    "the five non-MitoCarta brakes that are not in those sets. One point per ",
    "gene, positioned by its raw log2 fold change; filled points are significant."),
  detail = c(
    sprintf("n = 6 wild-type animals per timepoint, %d transcripts. DESeq2 raw (unshrunken) log2 fold changes on the `6>12W_wt` contrast, Benjamini-Hochberg adjusted genome-wide.",
            nrow(d)),
    sprintf("%d of %d transcripts reaches padj < 0.05, and it moves the WRONG WAY for a loss of priming: %s %+.2f (padj %.3f), a pro-apoptotic gene going UP. Nothing else clears 0.05 on either arm or among the brakes.",
            n_sig_wt, nrow(d), one$gene[1], one$lfc[1], one$padj[1]),
    sprintf("OVERALL PRIMING IS A BALANCE, and both arms move the same way, so the balance does not move. On the content ruler the pro arm shifts %+.4f and the anti arm %+.4f across the window, a difference of %+.4f; on the content-blind mitoPPS ruler %+.4f and %+.4f, a difference of %+.4f. (Script 40's ruler, the same instrument as Fig. 2F.)",
            comp$content[1], comp$content[2], prime_wt,
            comp$priority[1], comp$priority[2], prime_p),
    sprintf("FOR CONTRAST, MYC MOVES THE TWO ARMS APART: at six weeks the genotype effect is %+.4f on the pro arm and %+.4f on the anti arm, a difference of %+.4f -- an order of magnitude larger than anything the window does, and in the opposite geometry (apart rather than together). The balance is something Myc changes and development does not.",
            comp$myc_content[1], comp$myc_content[2], prime_myc),
    sprintf("THE POWER CONTROL, because a negative at n = 6 needs one: on these SAME %d transcripts, in the same libraries and at the same n, the genotype contrast at six weeks reaches padj < 0.05 for %d of them, against %d across the window. The measurement can see movement in these genes; there is none to see here.",
            sum(!is.na(d$padj_myc)), n_sig_myc, n_sig_wt),
    sprintf("Three genes are labelled because later panels turn on them: %s, the only mover; %s (PUMA) at %+.3f, padj %.2f, which is flat across the window and is the subject of Figs. 2G and 2H; and %s (Bcl-xL) at %+.3f, padj %.2f, the denominator of the ratio Fig. 2G reports.",
            "Bnip3", "Bbc3", d$lfc[d$gene == "Bbc3"], d$padj[d$gene == "Bbc3"],
            "Bcl2l1", d$lfc[d$gene == "Bcl2l1"], d$padj[d$gene == "Bcl2l1"]),
    "The third row is the OTHER anti-apoptotic arm: the caspase inhibitors XIAP, cIAP1 and cIAP2 (Birc2, Birc3) and survivin (Birc5), plus the Bcl2a1b paralog. None is a MitoCarta gene, so none is in the apoptosis sets the first two rows are drawn from, and they answer a different question -- whether the gland BUFFERS against death rather than whether it moves its apoptotic transcripts. Their membership is asserted in the script, so the row label cannot drift away from what it names. Nothing in this row moves either (padj 0.32 to 0.88)."),
  bounds = c(
    "BATCH = TIMEPOINT. The 6W and 12W cohorts were extracted as two separate batches, so this panel reads as \"no detectable movement at n = 6 on a confounded axis\", not as \"no movement\". The power control above is what makes the first reading worth having.",
    "A NEGATIVE ON TRANSCRIPTS IS NOT A NEGATIVE ON PRIMING. Apoptotic priming is a property of the protein complement and of how close the mitochondrion sits to the threshold; transcript levels are a substrate for it, not a measurement of it. The measurement is BH3 profiling, and this panel is the reason to do it rather than a substitute.",
    "THE CIRCULARITY WARNING OF SCRIPT 34 DOES NOT APPLY HERE, and it should not be imported. `MITOCARTA_APOPTOSIS_PRO`/`_ANTI` are MitoCarta sets, so a coupling between a priming composite and a mitochondrial axis is mito-vs-mito and circular. This panel is a temporal contrast on the transcripts themselves; nothing here is correlated with a mitochondrial score.",
    "Two of the 32 apoptosis-set members have no current MGI symbol in the annotation table. They are in the distribution, unlabelled, and are excluded from the genotype-contrast power count, which is why that count is over the mapped transcripts.",
    "Vertical offsets within a row are deterministic (genes ordered by fold change, offsets cycling through a fixed ladder), not random jitter. They carry no information: only the horizontal position is data."),
  source = c(
    "results/collapse_module_ownership.rds (scripts/44_collapse_module_and_ownership.R) -- $wt_genes, the 25 pro- and 7 anti-apoptotic MitoCarta transcripts on the wild-type temporal contrast",
    "results/substrate_specificity_tradeoff.rds (scripts/43_substrate_specificity_and_tradeoff.R) -- $buffer, the five non-MitoCarta brakes",
    "results/background_vs_myc.rds (scripts/40_background_vs_myc_decomposition.R) -- $ruler, the Apoptosis / Apoptosis-PRO / Apoptosis-ANTI composites on both rulers",
    "results/interaction_results.rds and results/combined_df_annotated_raw.rds (scripts/03) -- adjusted p-values for the genotype contrast, and the symbol-to-Ensembl mapping the power control needs",
    "The panel is figures/fig04_substrate_specificity.R panel D ported into this directory's idiom"))

save_panel_p(p, "figS2_priming_balance", height = 38)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## every transcript, both contrasts, ranked by what the window does
  d[order(d$lfc), c("gene", "class", "baseMean", "lfc", "padj", "lfc_myc", "padj_myc")] |>
    print(row.names = FALSE, digits = 3)

  ## the ones that move under the genotype contrast but not across the window --
  ## the power control, gene by gene
  d[!is.na(d$padj_myc) & d$padj_myc < 0.05,
    c("gene", "class", "lfc", "padj", "lfc_myc", "padj_myc")] |>
    print(row.names = FALSE, digits = 3)

  ## the composites the legend quotes, on both rulers
  comp |> print(row.names = FALSE, digits = 4)
  c(priming_window_content = prime_wt, priming_window_priority = prime_p,
    priming_myc_content = prime_myc) |> round(4) |> print()

  ## the whole dynamics-and-surveillance tier for context: apoptosis is one of the
  ## quietest rows in it across the window
  rr <- as.data.frame(bv$ruler)
  rr[rr$tier == "Mitochondrial dynamics and surveillance",
     c("pathway", "n_genes", "c_tn", "p_tn")] |>
    print(row.names = FALSE, digits = 3)
}
