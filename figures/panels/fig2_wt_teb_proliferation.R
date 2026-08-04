# =============================================================================
# fig2_wt_teb_proliferation.R -- what the normal gland loses between six and
# twelve weeks, and what it keeps
# -----------------------------------------------------------------------------
# SLOT: Fig. 2E.
#
#   "The 6>12W_wt comparison showed that while the TEB signature was lost as
#    expected in the puberty-adult transition, the overall proliferation
#    signalling remained relatively stable (Fig. 2E)."
#
# THE JOB. This is the first panel of section 2 and it establishes the substrate:
# what the wild-type mammary epithelium withdraws from over the pubertal-to-adult
# window is a MORPHOGENETIC programme, not a proliferative one. Everything after
# it depends on that -- the MYC-ER animals never saw Myc, so the change that
# licenses the whole of Figure 2 has to be developmental. It is also the negative
# limb that keeps Fig. 2F honest: had proliferation collapsed here, the
# respiratory withdrawal of 2F would be its shadow rather than a finding.
#
# THE RULER IS THE fGSEA NES on `timepoint_neg` (= 6>12W_wt), which is the same
# instrument Figs. 1D and 1H use, so "enrichment" means one thing across the
# figure set. ONE RANKED LIST, SO THE NORMALISER IS SHARED: the comparison this
# panel makes is WITHIN a single ranking, which is safe in exactly the way Fig.
# 1H's comparison ACROSS two rankings was not (there the flatter twelve-week list
# shrinks the normaliser and lifts every NES). Said again in the legend block.
#
# TWO ROWS, AND NOTHING ELSE (author's review, 2026-08-04). The first version drew
# four rows, a median rule per row and a label on the strongest set, and it read
# as clutter. What is left is the sentence and only the sentence:
#
#   TEB programmes (38)   median NES -1.58, 29 significant, ALL 38 depleted
#   proliferation  (14)   median NES -0.89, NOT ONE significant
#
# THREE THINGS CAME OFF THE PAGE AND STAYED IN THE LEGEND BLOCK, computed and
# asserted below so nothing is lost:
#   - the ductal `_VS_DUCTAL_*_DN` sets, which rise where the TEB sets fall and
#     are the DIRECTION control (author's call: omit);
#   - the 12 `*_TEB_MITO` construction lanes, which do not track the programme
#     they are named for (author: does not add to the interpretation);
#   - the two other rulers -- the gene-level set-average raw log2FC with its
#     expression-matched null, and the purity-adjusted per-sample composite.
# What must NOT happen is folding the first two INTO the TEB row: the `_DN` sets
# have their sign inverted by construction and would cancel part of the effect,
# and the `_MITO` lanes are MitoCarta subsets by build (Fig. 1D's rule). Omitted
# from the page is not the same as merged into the row.
#
# Reads (read-only, no re-run):
#   results/fgsea_percategory.rds              (script 20) -- the drawn ruler
#   results/substrate_specificity_tradeoff.rds (script 43) -- $comparator, $wt_null,
#                                                             $defs (rosters, guard)
#   results/priming_arm_teb.rds                (script 42) -- $teb_signatures
# Output: outputs/figures/panels/fig2_wt_teb_proliferation.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

fg_path <- here::here("results", "fgsea_percategory.rds")
ss_path <- here::here("results", "substrate_specificity_tradeoff.rds")
pa_path <- here::here("results", "priming_arm_teb.rds")
require_fresher_than(fg_path)      # 12:37, newer than gsva_scores.rds at 12:36
require_fresher_than(ss_path)
require_fresher_than(pa_path)

fg <- as.data.frame(readRDS(fg_path)$fgsea)
ss <- readRDS(ss_path)
pa <- readRDS(pa_path)
stopifnot(all(c("ranking", "category", "pathway", "NES", "pval",
                "padj_within_category", "size") %in% names(fg)))

# Scripts 42 and 43 save TIBBLES, and a tibble's `[` returns a tibble rather than
# a scalar -- which silently poisons every sprintf downstream (Fig. 1E's trap).
cmp <- as.data.frame(ss$comparator)
wtn <- as.data.frame(ss$wt_null)
tsg <- as.data.frame(pa$teb_signatures)

# --- guard: a stale results object must fail here, not at review --------------
# Script 43's own positive control, reused verbatim from figures/fig04: the
# wild-type OXPHOS-subunit value reproduces Issue #4's -0.2548. If membership
# resolution has drifted, nothing in this panel is comparable with the committed
# numbers.
ox_obs <- cmp$c_wt_time[cmp$arm == "OXPHOS subunits"]
if (!isTRUE(abs(ox_obs - ss$defs$oxphos_wt_reference) < 0.02))
  stop(sprintf("fig2E: stale results object -- OXPHOS-subunit wild-type %.4f against the reference %.4f. Re-run script 43.",
               ox_obs, ss$defs$oxphos_wt_reference))

# =============================================================================
# the families
# =============================================================================
# `programme_group()` is deliberately NOT the grouping here: it would scatter the
# TEB lanes across five of its rows (mammary development, TF target sets, curated
# mitochondrial ...) and the sentence is about the TEB context, which cuts across
# them. The rosters are therefore named here -- and the proliferation one is the
# ANALYSIS OF RECORD's, asserted identical to script 43's `defs$prolif_sets`, so
# the panel and the arm-level statistic quoted in the legend describe one set.
CONTRAST <- "timepoint_neg"          # = 6>12W_wt, the declared contrast vocabulary
ANCHOR   <- "MG_TEB_VS_DUCTAL_HS_GRAY_UP"

tn <- fg[fg$ranking == CONTRAST, ]
stopifnot(nrow(tn) == 866L)

# Four families are IDENTIFIED; two are DRAWN. Splitting them here rather than at
# the point of drawing is what lets the legend block quote the other two without
# either of them ever touching the TEB row.
fam_of <- function(pw)
  ifelse(pw %in% ss$defs$prolif_sets, "proliferation",
         ifelse(grepl("_MITO$", pw), "TEB x MitoCarta",
                ifelse(grepl("_DN$", pw), "ductal (TEB-down)", "TEB programmes")))

all4 <- tn[grepl("TEB", tn$pathway) | tn$pathway %in% ss$defs$prolif_sets, ]
all4$fam <- fam_of(all4$pathway)
all4$sig <- all4$padj_within_category < 0.05

DRAWN <- c("proliferation", "TEB programmes")     # level order: bottom, top
d <- all4[all4$fam %in% DRAWN, ]
d$fam <- factor(d$fam, levels = DRAWN)

n_all <- table(all4$fam)
stopifnot(
  # the proliferation roster IS script 43's, not a list re-typed here
  setequal(all4$pathway[all4$fam == "proliferation"], ss$defs$prolif_sets),
  n_all[["TEB programmes"]]    == 38L,
  n_all[["ductal (TEB-down)"]] ==  3L,
  n_all[["TEB x MitoCarta"]]   == 12L,
  n_all[["proliferation"]]     == 14L,
  # every TEB programme is depleted -- which is what leaves the top right empty
  # for the key, and is a stronger statement than the median
  max(d$NES[d$fam == "TEB programmes"]) < -1,
  ANCHOR %in% d$pathway,
  !anyNA(d$NES))

# Row labels carry n. Significance is a SHAPE (Fig. 1F's declared encoding) and
# never a fill or a point colour, so the counts stay off the page and go to the
# legend block.
d$row <- factor(sprintf("%s (%d)", as.character(d$fam), n_all[as.character(d$fam)]),
                levels = sprintf("%s (%d)", DRAWN, n_all[DRAWN]))

# =============================================================================
# the panel
# =============================================================================
XR  <- range(d$NES) + c(-1, 1) * diff(range(d$NES)) * 0.06
JIT <- ggplot2::position_jitter(width = 0, height = 0.30, seed = 3)

p <- ggplot2::ggplot(d, ggplot2::aes(NES, row)) +
  # +/-1 is the scale of an UNMOVED set: fGSEA divides the enrichment score by the
  # mean of the same-signed permutation null, so a set that is nowhere in
  # particular lands near |NES| = 1. It is the reference the proliferation row is
  # read against, and it is why no median rule is needed.
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "22",
                      linewidth = 0.25, colour = "grey78") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey45") +
  ggplot2::geom_point(ggplot2::aes(shape = sig), position = JIT,
                      size = 1.3, stroke = 0.32, colour = "grey20",
                      fill = "grey20") +
  ggplot2::scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1),
                              breaks = c(TRUE, FALSE),
                              labels = c("padj < 0.05", "n.s."), name = NULL) +
  ggplot2::scale_x_continuous(limits = XR, labels = lab_signed,
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.50)) +
  ggplot2::labs(x = "NES, 6>12W_wt", y = NULL) +
  ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 1.5))) +
  theme_panel(base_size = 6) +
  # Key inside, TOP RIGHT -- the wedge nothing occupies, because every TEB
  # programme is depleted, so the top row ends well left of zero.
  ggplot2::theme(
    legend.position        = "inside",
    legend.position.inside = c(0.995, 0.99),
    legend.justification   = c(1, 1),
    legend.background      = ggplot2::element_blank(),
    legend.margin          = ggplot2::margin(0, 0, 0, 0),
    legend.key.size        = ggplot2::unit(2.4, "mm"),
    legend.spacing.y       = ggplot2::unit(0.3, "mm"),
    axis.text.y            = ggplot2::element_text(size = 6),
    plot.margin            = ggplot2::margin(1.5, 2.5, 1, 1.5, "mm"))

# =============================================================================
# what is not drawn -- computed and asserted, reported in the legend
# =============================================================================
teb_gene  <- cmp[cmp$arm == "TEB vs ductal (HS)", ]
prol_gene <- cmp[cmp$arm == "PROLIF_* pooled", ]
teb_null  <- wtn[wtn$arm == "TEB vs ductal (HS)", ]
prol_null <- wtn[wtn$arm == "PROLIF_* pooled", ]
teb_samp  <- tsg[tsg$set == ANCHOR, ]
stopifnot(nrow(teb_gene) == 1L, nrow(prol_gene) == 1L, nrow(teb_samp) == 1L,
          # the two script-43 tables must describe the same arms
          abs(teb_null$observed_wt  - teb_gene$c_wt_time)  < 1e-8,
          abs(prol_null$observed_wt - prol_gene$c_wt_time) < 1e-8,
          # gene-level and purity-adjusted per-sample agree on the TEB arm
          abs(teb_samp$wt_time - teb_gene$c_wt_time) < 0.01)

nes_of   <- function(pw) all4$NES[all4$pathway == pw]
fam_line <- function(f) {
  s <- all4[all4$fam == f, ]
  sprintf("%s: %d sets, median NES %+.2f (range %+.2f to %+.2f), %d significant at padj < 0.05.",
          f, nrow(s), stats::median(s$NES), min(s$NES), max(s$NES), sum(s$sig))
}
# a self-contained multiplicity scope over the sets actually drawn
bh_panel <- stats::p.adjust(d$pval, "BH")

LEGEND <- panel_legend(
  slot = "Fig. 2E",
  what = paste0(
    "Gene-set enrichment across the wild-type 6 to 12 week transition. One point ",
    "per gene set, positioned by its fGSEA normalised enrichment score on the ",
    "`6>12W_wt` contrast: the terminal-end-bud programmes above, the curated ",
    "proliferation sets below. Filled points are significant. The dashed guides ",
    "at +/-1 are the scale of a set that has not moved."),
  detail = c(
    sprintf("n = %d gene sets drawn, on one ranked list (the unshrunken Wald statistic of the wild-type temporal contrast, `timepoint_neg`). %s %s",
            nrow(d), fam_line("TEB programmes"), fam_line("proliferation")),
    sprintf("EVERY ONE OF THE %d TEB PROGRAMMES IS DEPLETED (the least of them %+.2f), which is a stronger statement than the median: the loss is the whole family, not its leading members. The strongest is %s at NES %+.2f (p %.1e), which is also the set the two rulers below use.",
            n_all[["TEB programmes"]], max(d$NES[d$fam == "TEB programmes"]),
            ANCHOR, nes_of(ANCHOR), all4$pval[all4$pathway == ANCHOR]),
    sprintf("THREE INDEPENDENT RULERS AGREE ON THE TEB ARM, and two of them are independent of cell composition. Gene-level set-average raw log2FC %+.4f over %d genes, the %.1fth percentile of %d expression-matched random sets (empirical p %.4f); and a purity-adjusted per-sample composite (score ~ timepoint * genotype + epithelial + immune) %+.4f, p %.4f. The two agree to %.3f, so residual stromal or immune contamination does not explain the loss.",
            teb_gene$c_wt_time, teb_gene$n_genes, teb_null$percentile,
            ss$defs$n_set_draws, teb_null$p_emp_lower, teb_samp$wt_time,
            teb_samp$wt_p, abs(teb_samp$wt_time - teb_gene$c_wt_time)),
    sprintf("PROLIFERATION IS SMALL, NOT IMMOBILE, and the sentence should carry the qualifier. Not one of the %d proliferation sets is significant on this ruler, and the family median (%+.2f) sits inside the +/-1 band an unmoved set occupies. But on the gene-level ruler the pooled %d proliferation genes move %+.4f log2, which is the %.1fth percentile of their own matched null (p = %.3f). The honest form of \"relatively stable\" is the ratio: the TEB arm moves %.1f times as far.",
            n_all[["proliferation"]], stats::median(d$NES[d$fam == "proliferation"]),
            prol_gene$n_genes, prol_gene$c_wt_time, prol_null$percentile,
            prol_null$p_emp_lower, abs(teb_gene$c_wt_time / prol_gene$c_wt_time)),
    sprintf("THE DIRECTION IS CONTROLLED, though the control is not drawn (author's call, 2026-08-04). The three `_VS_DUCTAL_*_DN` sets hold the genes lower in the end bud than in the mature duct, so they are the reciprocal of the drawn sets by construction and a ranking that merely drifted downward could not move them the other way. Two of the three lineages show the reciprocal outright -- basal %+.2f against its UP set's %+.2f, hormone-sensing %+.2f against %+.2f -- while the alveolar pair shows neither (%+.2f against %+.2f, neither significant). Only the basal DN set clears padj < 0.05 (%.3f). The gland moves from the end-bud programme toward the ductal one.",
            nes_of("MG_TEB_VS_DUCTAL_BA_GRAY_DN"), nes_of("MG_TEB_VS_DUCTAL_BA_GRAY_UP"),
            nes_of("MG_TEB_VS_DUCTAL_HS_GRAY_DN"), nes_of(ANCHOR),
            nes_of("MG_TEB_VS_DUCTAL_AP_GRAY_DN"), nes_of("MG_TEB_VS_DUCTAL_AP_GRAY_UP"),
            all4$padj_within_category[all4$pathway == "MG_TEB_VS_DUCTAL_BA_GRAY_DN"]),
    sprintf("A further %d TEB-named sets are excluded from the drawn family and from every number above: the `*_TEB_MITO` lanes, whose membership is MitoCarta intersected with a TEB-context regulon (17 to 59 genes each). They are mitochondrial by build, they do not track the programme they are named for (median NES %+.2f, none significant), and they belong to Fig. 2F's question rather than this one. Excluding them is not the same as folding them in, which would have blurred the row.",
            n_all[["TEB x MitoCarta"]],
            stats::median(all4$NES[all4$fam == "TEB x MitoCarta"])),
    "THE GENERIC PROLIFERATION REGULONS AGREE WITH THE CURATED SETS -- TFT_E2F1_DOROTHEA_ABC -1.18 and TFT_E2F1_CHIPATLAS +0.86, neither significant -- but the Gray lineage-context lanes of the same transcription factors fall hard (FOXM1_HS_LE -2.10, MYBL2_HS_LE -2.04, TFDP1_HS_LE -2.08, E2F1_HS_LE -1.72, all padj < 1e-3). Those are lineage-identity lanes, not cell-cycle sets. What the normal gland withdraws from is lineage and morphogenesis; the cell cycle itself does not move.",
    sprintf("Multiplicity, both scopes: drawn significance is script 20's Benjamini-Hochberg WITHIN each library category, the analysis of record. A panel-scope BH over the %d sets drawn here changes nothing -- the same %d TEB programmes survive and the best proliferation set reaches %.3f.",
            nrow(d), sum(bh_panel < 0.05), min(bh_panel[d$fam == "proliferation"]))),
  bounds = c(
    "BATCH = TIMEPOINT. The 6W and 12W cohorts were extracted as two batches, so every value on this panel is DESCRIBED, not claimed. The mitigation is that the panel makes a COMPARISON, not a measurement: both families ride the same two batches, so a batch artefact would have to be TEB-specific to produce this picture.",
    "The comparison is WITHIN one ranked list, so the fGSEA normaliser is shared. That is what makes it safe here and is exactly what was not true of Fig. 1H, where the flatter twelve-week ranking shrinks the normaliser and lifts every NES. No statement about the SIZE of an NES is made across contrasts anywhere in this panel.",
    "fGSEA's gene permutation is anti-conservative for sets of correlated genes, which is most of this library; and the TEB family spans two categories (03_mammary_development, 06_tf_targets) whose BH scopes are separate by script 20's design.",
    "\"The TEB signature was lost as expected\" -- the expectation is external (Gray et al.'s pubertal-to-adult series); nothing on this panel tests it. What the panel shows is that this cohort reproduces it.",
    "The proliferation family is the 14 curated PROLIF_* sets, which is script 43's own roster and the same genes as the pooled arm quoted above. It is not a claim about every proliferative programme in the library."),
  source = c(
    "results/fgsea_percategory.rds (scripts/20_fgsea_percategory.R) -- NES per ranking x category on the unshrunken Wald statistic; the `timepoint_neg` ranking is the contrast drawn",
    "results/substrate_specificity_tradeoff.rds (scripts/43_substrate_specificity_and_tradeoff.R) -- $comparator and $wt_null for the gene-level arms and their expression-matched nulls; $defs$prolif_sets is the drawn proliferation roster",
    "results/priming_arm_teb.rds (scripts/42_priming_arm_and_teb_substrate.R) -- $teb_signatures, the purity-adjusted per-sample composite"))

save_panel_p(p, "fig2_wt_teb_proliferation", height = 30)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## every set of all four families, drawn or not, by rank
  all4[order(all4$fam, all4$NES),
       c("fam", "pathway", "NES", "padj_within_category", "size")] |>
    print(row.names = FALSE, digits = 3)

  ## the two families that came off the page, in one look
  all4[all4$fam %in% c("ductal (TEB-down)", "TEB x MitoCarta"),
       c("fam", "pathway", "NES", "padj_within_category", "size")] |>
    (\(x) x[order(x$fam, -x$NES), ])() |> print(row.names = FALSE, digits = 3)

  ## the TEB programmes by lineage context -- the HS (hormone-sensing) context is
  ## the apoptosis-loaded one (TEB lumen clearance) and it falls furthest
  tb  <- d[d$fam == "TEB programmes", ]
  ctx <- ifelse(grepl("_HS_", tb$pathway), "HS",
         ifelse(grepl("_BA_", tb$pathway), "BA",
         ifelse(grepl("_AP_", tb$pathway), "AP", "other")))
  tapply(tb$NES, ctx, function(v) c(n = length(v), median = stats::median(v))) |> print()

  ## the two rulers the panel does not draw
  cmp[cmp$arm %in% c("TEB vs ductal (HS)", "PROLIF_* pooled"), ] |> print(digits = 4)
  wtn[wtn$arm %in% c("TEB vs ductal (HS)", "PROLIF_* pooled"), ] |> print(digits = 4)

  ## the generic vs lineage-context regulons of the proliferation TFs
  tn[grepl("^TFT_(E2F1|E2F7|FOXM1|MYBL2|CENPA|TFDP1)_", tn$pathway) &
     !grepl("TEB", tn$pathway),
     c("pathway", "NES", "padj_within_category", "size")] |>
    (\(x) x[order(x$NES), ])() |> print(row.names = FALSE, digits = 3)

  ## the mitochondrial arms, which are Fig. 2F's subject, on this same ruler
  tn[tn$category == "01_mitocarta" & grepl("OXPHOS|COMPLEX", tn$pathway),
     c("pathway", "NES", "padj_within_category", "size")] |>
    (\(x) x[order(x$NES), ])() |> print(row.names = FALSE, digits = 3)
}
