# =============================================================================
# fig2_departure_from_dose.R -- the two genes that leave the dose line together
# -----------------------------------------------------------------------------
# SLOT: Fig. 2H.
#
#   "... closely tracking Bbc3 was its known p53-independent activator, Foxo3."
#
# THE SCAN. Fig. 1G showed that the Myc effect at twelve weeks is the six-week
# effect rescaled by a global factor. Script 44 turns that into a per-gene
# residual: for each of 8,774 genes, how far its twelve-week effect differs from
# what the global Myc scaling predicts, divided by its own standard error, signed so
# that NEGATIVE means "collapsed further than the dose explains".
#
#   Bbc3   z = -2.49, 44 of 8,774 genes below it       (percentile 0.51)
#   Foxo3  z = -2.53, 39 below it                      (percentile 0.46)
#
# and ONLY FOUR GENES LIE BETWEEN THEM. That is what "closely tracking" means
# here, and it is the panel.
#
# THE PANEL ALSO CARRIES THE HONEST PART, which is why the third mark is drawn.
# The scan was pre-registered with TWO genes (script 44's `pre_specified_genes`:
# Bbc3 and Bcl2l11) and they SPLIT: Bbc3 collapses at the 0.51st percentile,
# Bcl2l11 sits at the 91.6th, which is the ordinary middle of the distribution.
# So the pre-registration half-succeeded, and drawing only the half that worked
# would be the wrong panel. WHICH gene was pre-specified is not drawn -- all three
# marks are identical -- because that is a fact about the analysis rather than
# about a gene's position. It is in the legend block.
#
# AND FOXO3 WAS NOT PRE-SPECIFIED. It was found in this scan. What makes it a
# lead rather than one of forty-four names is the prior -- FOXO3 is PUMA's
# canonical p53-independent activator, with direct ChIP evidence -- not its
# position, and the legend block says so in as many words. The forty-three other
# genes below the 0.5th percentile are mostly unrelated (Fut9, Cntnap2, Inhba,
# Wfdc2, Vtcn1), which is exactly why adjacency on its own proves nothing.
#
# THE MECHANISM IS IN THE TWO GENES' OWN NUMBERS, not on the panel:
#   Foxo3  Myc effect +0.243 at 6W, -0.242 at 12W; the WILD-TYPE gland raises it
#          with age (+0.473) and the Myc+ gland does not (-0.013)
#   Bbc3   +0.258 to -0.283; wild-type flat (+0.062), Myc+ falls (-0.479)
# so the same interaction is reached by two different routes, which is a
# distinction the text should keep.
#
# Reads (read-only, no re-run):
#   results/collapse_module_ownership.rds (script 44) -- $collapse_genes, the scan
#                                            itself; $defs for the rate and the
#                                            pre-specified roster
#   results/interaction_results.rds       (script 03) -- the interaction p-values
#   results/combined_df_annotated_raw.rds (script 03) -- symbol -> Ensembl
# Output: outputs/figures/panels/fig2_departure_from_dose.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))
if (!requireNamespace("DESeq2", quietly = TRUE))
  stop("fig2H needs DESeq2 to coerce the DESeqResults in interaction_results.rds")

cm_path <- here::here("results", "collapse_module_ownership.rds")
require_fresher_than(cm_path)
cmo <- readRDS(cm_path)
ir  <- readRDS(here::here("results", "interaction_results.rds"))
ann <- as.data.frame(readRDS(here::here("results", "combined_df_annotated_raw.rds")))

cg <- as.data.frame(cmo$collapse_genes)
stopifnot(nrow(cg) == 8774L,
          all(c("gene", "z_resid", "pct_z", "lfc_6W", "lfc_12W") %in% names(cg)),
          identical(sort(cmo$defs$pre_specified_genes), sort(c("Bbc3", "Bcl2l11"))))

# =============================================================================
# the three marks
# =============================================================================
# Bbc3 and Foxo3 are the sentence; Bcl2l11 is the pre-registration's other half
# and is drawn because leaving it out would show only the gene that worked.
MARK <- data.frame(
  gene = c("Foxo3", "Bbc3", "Bcl2l11"),
  role = c("found in the scan", "pre-specified", "pre-specified"),
  stringsAsFactors = FALSE)
MARK$z   <- cg$z_resid[match(MARK$gene, cg$gene)]
MARK$pct <- cg$pct_z[match(MARK$gene, cg$gene)]
MARK$below <- vapply(MARK$z, function(v) sum(cg$z_resid < v), integer(1))
stopifnot(!anyNA(MARK$z),
          # the claim, asserted: the two are adjacent and both in the far tail
          sum(cg$z_resid > MARK$z[MARK$gene == "Foxo3"] &
              cg$z_resid < MARK$z[MARK$gene == "Bbc3"]) == 4L,
          MARK$pct[MARK$gene == "Bbc3"]  < 1,
          MARK$pct[MARK$gene == "Foxo3"] < 1,
          MARK$pct[MARK$gene == "Bcl2l11"] > 50)

# interaction p-values, keyed by Ensembl; the mapping is checked against the
# scan's own fold changes rather than trusted
int <- as.data.frame(ir$interaction_raw)
MARK$ens <- ann$gene[match(MARK$gene, ann$mgi_symbol)]
MARK$int_p    <- int$pvalue[match(MARK$ens, rownames(int))]
MARK$int_padj <- int$padj[match(MARK$ens, rownames(int))]
stopifnot(!anyNA(MARK$int_p),
          max(abs(ann$myc_6W_log2FC_raw[match(MARK$gene, ann$mgi_symbol)] -
                  cg$lfc_6W[match(MARK$gene, cg$gene)])) < 1e-6)

# =============================================================================
# the panel
# =============================================================================
# Windowed for display: a handful of genes sit past z = +4 and are named in the
# legend -- the right tail is real (these are genes Myc RETAINS better than the
# dose predicts) but it is empty enough to waste a third of the panel.
# Every number is computed on the complete 8,774.
XLO <- min(cg$z_resid) - 0.2; XHI <- 4
out <- sum(cg$z_resid > XHI)
dz  <- cg$z_resid[cg$z_resid <= XHI]

den <- stats::density(cg$z_resid, n = 1024)
dd  <- data.frame(x = den$x, y = den$y)
dd  <- dd[dd$x >= XLO & dd$x <= XHI, ]
HMAX <- max(dd$y)

# The marks stand above the curve so they are visible where it is flat, and the
# two adjacent ones take different heights so their labels can sit on their own
# lines -- they are 0.04 apart in z, which is less than a character.
MARK$h   <- c(0.62, 0.42, 0.42) * HMAX          # Foxo3 higher than Bbc3
MARK$lab <- MARK$gene
tail_rug <- cg[cg$z_resid <= stats::quantile(cg$z_resid, 0.01), ]

p <- ggplot2::ggplot(dd, ggplot2::aes(x, y)) +
  ggplot2::geom_area(fill = "grey88") +
  ggplot2::geom_line(linewidth = 0.3, colour = "grey35") +
  # every gene in the bottom percentile, so the tail's sparseness is visible
  ggplot2::geom_rug(data = tail_rug, inherit.aes = FALSE,
                    ggplot2::aes(x = z_resid), sides = "b",
                    length = ggplot2::unit(1.4, "mm"),
                    linewidth = 0.18, colour = "grey55") +
  ggplot2::geom_segment(data = MARK, inherit.aes = FALSE,
                        ggplot2::aes(x = z, xend = z, y = 0, yend = h),
                        linewidth = 0.35, colour = "grey15") +
  ggplot2::geom_point(data = MARK, inherit.aes = FALSE,
                      ggplot2::aes(x = z, y = h),
                      size = 1.5, stroke = 0.35, colour = "grey15") +
  ggplot2::geom_text(data = MARK, inherit.aes = FALSE,
                     ggplot2::aes(x = z, y = h, label = lab),
                     hjust = -0.18, vjust = 0.4, size = 1.8, colour = "grey15",
                     fontface = "italic") +
  ggplot2::scale_x_continuous(limits = c(XLO, XHI), labels = lab_signed,
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::scale_y_continuous(limits = c(0, HMAX * 1.06),
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::labs(x = "difference from global Myc scaling  (z)",
                y = "genes (density)") +
  theme_panel(base_size = 6) +
  # NO KEY, and therefore no shape encoding either (author's review, 2026-08-05:
  # the panel is obvious without one). The three marks are identical; WHICH of
  # them was pre-specified and which was found in the scan is a fact about the
  # analysis, not a property of a gene's position, so it belongs in the legend
  # block. An unexplained difference between glyphs on the page would be worse
  # than no difference.
  ggplot2::theme(
    axis.text.y            = ggplot2::element_blank(),
    axis.ticks.y           = ggplot2::element_blank(),
    plot.margin            = ggplot2::margin(1.5, 2.5, 1, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
m       <- function(g, col) MARK[[col]][MARK$gene == g]
scan_of <- function(g, col) cg[[col]][match(g, cg$gene)]
tail05  <- cg[cg$z_resid <= stats::quantile(cg$z_resid, 0.005), ]
near    <- sort(na.omit(tail05$gene))

LEGEND <- panel_legend(
  slot = "Fig. 2H",
  what = paste0(
    "How far each gene's Myc effect at twelve weeks departs from what the global ",
    "dose rescaling predicts, over 8,774 genes: the residual divided by its own ",
    "standard error, signed so that negative means collapsed further than the ",
    "dose explains. Ticks along the foot are every gene in the bottom percentile. ",
    "Three genes are marked -- the two the scan was pre-registered with, and the ",
    "one this panel's sentence is about."),
  detail = c(
    sprintf("n = %s genes (baseMean >= 20 and reported at both ages). The expected twelve-week effect is the six-week effect times %.3f, the global rescaling rate fitted over the 2,648 genes Myc moves (Fig. 1G).",
            format(nrow(cg), big.mark = ","), cmo$defs$global_rate_fitted),
    sprintf("THE CLAIM, AND IT IS SHARP: %s sits at z = %+.2f with %d of %s genes below it (percentile %.2f) and %s at z = %+.2f with %d below it (percentile %.2f) -- and ONLY FOUR GENES LIE BETWEEN THEM.",
            "Foxo3", m("Foxo3", "z"), m("Foxo3", "below"),
            format(nrow(cg), big.mark = ","), m("Foxo3", "pct"),
            "Bbc3", m("Bbc3", "z"), m("Bbc3", "below"), m("Bbc3", "pct")),
    sprintf("THE PRE-REGISTRATION HAD TWO GENES AND THEY SPLIT. Script 44's `pre_specified_genes` are Bbc3 and Bcl2l11, named in advance from the cell experiments. Bbc3 collapses at percentile %.2f; Bcl2l11 sits at percentile %.1f (z %+.2f), the ordinary middle of the distribution, with an interaction p of %.2f. Drawing only the half that worked would be the wrong panel.",
            m("Bbc3", "pct"), m("Bcl2l11", "pct"), m("Bcl2l11", "z"),
            m("Bcl2l11", "int_p")),
    sprintf("The interactions themselves: %s %.4f, %s %.4f, %s %.2f. None survives genome-wide correction (Bbc3's Benjamini-Hochberg value is %.2f), which is the expected outcome for an interaction at n = 6 per cell and is why the licence is pre-specification rather than the p-value.",
            "Bbc3", m("Bbc3", "int_p"), "Foxo3", m("Foxo3", "int_p"),
            "Bcl2l11", m("Bcl2l11", "int_p"), m("Bbc3", "int_padj")),
    sprintf("THE SAME DEPARTURE IS REACHED BY TWO DIFFERENT ROUTES, and the text should keep the distinction. Foxo3: Myc raises it at six weeks (%+.3f) and lowers it at twelve (%+.3f), and the WILD-TYPE gland raises it with age (%+.3f) while the Myc+ gland does not (%+.3f). Bbc3: %+.3f to %+.3f, with the wild-type gland flat (%+.3f) and the Myc+ gland falling (%+.3f).",
            scan_of("Foxo3", "lfc_6W"), scan_of("Foxo3", "lfc_12W"),
            ann$timepoint_neg_log2FC_raw[match("Foxo3", ann$mgi_symbol)],
            ann$timepoint_pos_log2FC_raw[match("Foxo3", ann$mgi_symbol)],
            scan_of("Bbc3", "lfc_6W"), scan_of("Bbc3", "lfc_12W"),
            ann$timepoint_neg_log2FC_raw[match("Bbc3", ann$mgi_symbol)],
            ann$timepoint_pos_log2FC_raw[match("Bbc3", ann$mgi_symbol)]),
    sprintf("For scale at the other end, Myc's own transcript sits at z %+.2f, the %.2fth percentile -- the most retained gene in the transcriptome, which is Fig. S1E's result read on this scan.",
            scan_of("Myc", "z_resid"), scan_of("Myc", "pct_z")),
    sprintf("%d genes lie past the drawn window at z > %.0f and are not shown (%s%s) -- they are genes Myc RETAINS better than the dose predicts, which is the opposite tail from this panel's subject. Every number here is computed on the complete %s.",
            out, XHI, paste(sort(na.omit(cg$gene[cg$z_resid > XHI])), collapse = ", "),
            if (sum(is.na(cg$gene[cg$z_resid > XHI]))) sprintf(", and %d without a current symbol",
                                                               sum(is.na(cg$gene[cg$z_resid > XHI]))) else "",
            format(nrow(cg), big.mark = ","))),
  bounds = c(
    sprintf("FOXO3 WAS NOT PRE-SPECIFIED. It was found in this scan, and %d genes sit below the 0.5th percentile with it -- mostly unrelated (%s). Position alone is therefore not evidence: what makes Foxo3 a lead rather than one of forty-four names is the PRIOR, that FOXO3 is PUMA's canonical p53-independent activator with direct ChIP evidence. The text should introduce it that way round.",
            nrow(tail05), paste(utils::head(near[!near %in% c("Foxo3", "Bbc3")], 6),
                                collapse = ", ")),
    "AND THE ARROW IS NOT THE OBVIOUS ONE. PUMA restrains the mitochondrial pyruvate carrier (Kim, Cancer Cell 2019), so \"PUMA falls, therefore respiration falls\" is backwards; respiration sits upstream (Dey & Moraes), and FOXO3 -> BBC3 closes it into a negative-feedback circuit rather than a linear chain. Nothing on this panel establishes a direction.",
    "A z of -2.5 among 8,774 genes is a RANK STATEMENT, not a test. It says these two are in the far tail of the residual distribution; it does not say the residual is significant, and the interaction p-values above make clear that none of them is after correction.",
    "The residual is measured against a rate fitted on OTHER genes (the 2,648 Myc-responsive ones), so a gene that is itself in that set contributes to its own expectation. Neither Bbc3 nor Foxo3 is: both fall outside the ranking set (their six-week adjusted p-values are 0.22 and 0.21), which is a point in favour of the scan and against reading their six-week effects as established.",
    "Genotype contrasts are clean, but an interaction is a difference of two 6-versus-6 contrasts and is the least powered quantity in the design (median lfcSE 0.333 against 0.233)."),
  source = c(
    "results/collapse_module_ownership.rds (scripts/44_collapse_module_and_ownership.R) -- $collapse_genes, the per-gene departure scan; $defs$global_rate_fitted and $defs$pre_specified_genes",
    "results/interaction_results.rds (scripts/03_deseq_results_qc.R) -- raw interaction p-values",
    "results/combined_df_annotated_raw.rds (scripts/03) -- the symbol-to-Ensembl mapping and the two temporal contrasts quoted in the detail"))

save_panel_p(p, "fig2_departure_from_dose", height = 44)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the three marked genes, with everything the legend quotes
  MARK |> print(row.names = FALSE, digits = 3)

  ## the company Foxo3 keeps: every gene below the 0.5th percentile
  cg[cg$z_resid <= stats::quantile(cg$z_resid, 0.005),
     c("gene", "baseMean", "lfc_6W", "lfc_12W", "z_resid", "pct_z", "myc_induced")] |>
    (\(x) x[order(x$z_resid), ])() |> print(row.names = FALSE, digits = 3)

  ## the four genes between Foxo3 and Bbc3 -- the whole of "closely tracking"
  cg[cg$z_resid > MARK$z[MARK$gene == "Foxo3"] &
     cg$z_resid < MARK$z[MARK$gene == "Bbc3"],
     c("gene", "baseMean", "lfc_6W", "lfc_12W", "z_resid")] |>
    print(row.names = FALSE, digits = 3)

  ## the other tail, for scale
  cg[order(-cg$z_resid), c("gene", "baseMean", "lfc_6W", "lfc_12W", "z_resid",
                           "pct_z")] |> head(10) |> print(row.names = FALSE, digits = 3)

  ## the death genes on this scan -- PUMA is a solo, which is the wider result
  cg[cg$gene %in% c("Bbc3", "Bcl2l11", "Bax", "Bak1", "Bid", "Bcl2l1", "Bcl2",
                    "Mcl1", "Pmaip1", "Bmf", "Bnip3", "Foxo3"),
     c("gene", "lfc_6W", "lfc_12W", "z_resid", "pct_z")] |>
    (\(x) x[order(x$z_resid), ])() |> print(row.names = FALSE, digits = 3)
}
