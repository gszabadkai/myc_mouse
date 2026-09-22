# =============================================================================
# biogax_regulon_split.R -- a regulon that splits in two directions
# -----------------------------------------------------------------------------
# DISCUSSION PANEL, not a manuscript slot. The `biogax_` prefix keeps it outside
# rebuild_panels.R and panels_to_pdf.R, both of which glob `^fig.*\.R$`: this
# panel takes no slot, does not enter the panel count, and does not appear in
# paper/analysis_record.qmd. It exists to be argued over with collaborators, and
# it is collected by paper/biogenesis_axis_discussion.qmd.
#
# THE DECISIVE TEST. A change in a transcription factor's ACTIVITY acts on its
# REGULON. It cannot act on a functional subset of the regulon. So if the OXPHOS
# members of the ERRa/NRF1/GABP regulon fall while the other members rise, no
# activity change of those factors can be the explanation -- whatever moved the
# respiratory genes moved them AS respiratory genes, not as targets.
#
# WHAT IS DRAWN: for each regulon, the matched-random-set null percentile of its
# OXPHOS-subunit members against that of everything else in the same regulon,
# joined by a segment. The segment IS the result; its length is the dissociation.
#
# WHY PERCENTILE AND NOT THE FOLD CHANGE. Set-mean log2 fold changes are not
# comparable between sets of different size and expression: a -0.2 over 59 genes
# and a -0.2 over 9 are different evidence. The percentile against 2000
# expression-matched random sets of the same size puts every part of every
# regulon on one axis. The content values are in the legend block.
#
# Reads : results/biogenesis_axis_developmental.rds (script 47 PART B)
# Output: outputs/figures/panels/biogax_regulon_split.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

ba_path <- here::here("results", "biogenesis_axis_developmental.rds")
if (!file.exists(ba_path)) stop("run scripts/47_... first")
ba <- readRDS(ba_path)

sp <- as.data.frame(ba$split)
stopifnot(all(c("regulon", "part", "n_genes", "c_wt_time", "wt_null_pct") %in% names(sp)))

# Order: the four regulons that actually contain the respiratory chain first, by
# how much of it they hold, then the three that do not. That ordering is the
# panel's second reading -- the split appears exactly where there is enough
# respiratory content to show one.
# Order is DECLARED, not derived. A derived order (sort by respiratory content)
# is one tie away from rearranging itself on a re-run -- MYC_MITO and
# DEVELOPMENTAL_MITO both hold 9 subunits -- and the reading depends on the four
# testable regulons appearing together at the top. ggplot2 puts the first level at
# the BOTTOM of a discrete y axis, so the vector is written bottom-up.
ORD_TOP_DOWN <- c("CORE_MITO", "ESRRA_MITO", "NRF1_MITO", "GABPA_MITO",
                  "MYC_MITO", "DEVELOPMENTAL_MITO", "E2F1_MITO")
stopifnot(setequal(ORD_TOP_DOWN, unique(sp$regulon)))
ord <- rev(ORD_TOP_DOWN)

n_of <- function(r, prt, fallback = NA_integer_) {
  v <- sp$n_genes[sp$regulon == r & sp$part == prt]
  if (length(v) == 1L) v else fallback
}
lab_of <- function(r) sprintf("%s\n%d OXPHOS / %d other", sub("_MITO$", "", r),
                              n_of(r, "n OXPHOS subunits", 2L), n_of(r, "rest"))

d <- sp[sp$part %in% c("n OXPHOS subunits", "rest"), ]
d$regulon <- factor(d$regulon, levels = ord, labels = vapply(ord, lab_of, character(1)))
d$part <- factor(d$part, levels = c("n OXPHOS subunits", "rest"),
                 labels = c("its OXPHOS subunits", "the rest of the regulon"))

seg <- merge(
  d[d$part == "its OXPHOS subunits",     c("regulon", "wt_null_pct")],
  d[d$part == "the rest of the regulon", c("regulon", "wt_null_pct")],
  by = "regulon", suffixes = c("_ox", "_rest"))

whole <- sp[sp$part == "all", ]
whole$regulon <- factor(whole$regulon, levels = ord,
                        labels = vapply(ord, lab_of, character(1)))

# ASSERTIONS: the reading must not be able to change silently on a re-run.
stopifnot(
  sp$wt_null_pct[sp$regulon == "CORE_MITO" & sp$part == "n OXPHOS subunits"] < 5,
  sp$wt_null_pct[sp$regulon == "CORE_MITO" & sp$part == "rest"]              > 95,
  sp$wt_null_pct[sp$regulon == "CORE_MITO" & sp$part == "all"] > 5,
  sp$wt_null_pct[sp$regulon == "CORE_MITO" & sp$part == "all"] < 95)

# `verdict_cols` is the declared palette for exactly this reading -- withdraws /
# at chance / rises -- so the panel takes it rather than inventing a pair.
part_cols <- c("its OXPHOS subunits"     = unname(verdict_cols[["withdraws"]]),
               "the rest of the regulon" = unname(verdict_cols[["rises"]]))

p <- ggplot2::ggplot(d, ggplot2::aes(x = wt_null_pct, y = regulon)) +
  ggplot2::annotate("rect", xmin = 5, xmax = 95, ymin = -Inf, ymax = Inf,
                    fill = "grey96", colour = NA) +
  ggplot2::geom_vline(xintercept = 50, linewidth = 0.22, colour = "grey65") +
  ggplot2::geom_segment(data = seg,
    ggplot2::aes(x = wt_null_pct_ox, xend = wt_null_pct_rest,
                 y = regulon, yend = regulon),
    inherit.aes = FALSE, linewidth = 0.5, colour = "grey40") +
  ggplot2::geom_point(data = whole,
    ggplot2::aes(x = wt_null_pct, y = regulon),
    inherit.aes = FALSE, shape = 124, size = 2.4,
    colour = unname(verdict_cols[["at chance"]])) +
  ggplot2::geom_point(ggplot2::aes(fill = part), shape = 21, size = 2.3,
                      stroke = 0.25, colour = "white") +
  ggplot2::scale_fill_manual(values = part_cols, name = NULL) +
  ggplot2::scale_x_continuous(
    limits = c(-2, 102), breaks = c(0, 5, 25, 50, 75, 95), expand = c(0, 0),
    name = "percentile of 2000 expression-matched random sets") +
  ggplot2::scale_y_discrete(name = NULL) +
  theme_panel() +
  ggplot2::theme(legend.position = "top",
                 legend.margin = ggplot2::margin(0, 0, -2, 0),
                 panel.grid.major.y = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(lineheight = 0.95))

LEGEND <- panel_legend(
  slot = "Discussion D1",
  what = paste(
    "Within each mitochondrial regulon, the wild-type 6-to-12-week change of its",
    "OXPHOS-subunit members against that of all its other members, expressed as a",
    "percentile of 2000 expression-matched random gene sets. The vertical tick is",
    "the regulon taken whole."),
  detail = c(
    sprintf("CORE_MITO (the ERRa/NRF1/GABP regulon): OXPHOS members %.1f pct, the other %d genes %.1f pct, whole regulon %.1f pct.",
            sp$wt_null_pct[sp$regulon == "CORE_MITO" & sp$part == "n OXPHOS subunits"],
            sp$n_genes[sp$regulon == "CORE_MITO" & sp$part == "rest"],
            sp$wt_null_pct[sp$regulon == "CORE_MITO" & sp$part == "rest"],
            sp$wt_null_pct[sp$regulon == "CORE_MITO" & sp$part == "all"]),
    sprintf("Content values behind those percentiles: %+.3f, %+.3f, %+.3f log2.",
            sp$c_wt_time[sp$regulon == "CORE_MITO" & sp$part == "n OXPHOS subunits"],
            sp$c_wt_time[sp$regulon == "CORE_MITO" & sp$part == "rest"],
            sp$c_wt_time[sp$regulon == "CORE_MITO" & sp$part == "all"]),
    paste("All four regulons that contain the respiratory chain split the same way.",
          "The three that do not split hold 2 to 9 OXPHOS subunits between them, so",
          "the test cannot run: read it as four of four, not four of seven."),
    sprintf("DEVELOPMENTAL_MITO is the one regulon that moves as a UNIT, at the %.1f percentile -- upward. That is the shape an activity change makes.",
            sp$wt_null_pct[sp$regulon == "DEVELOPMENTAL_MITO" & sp$part == "all"])),
  bounds = c(
    "BATCH = TIMEPOINT: the between-age contrast is confounded with cohort. What survives that is precisely this comparison, which is INTERNAL to each regulon -- a shared batch effect cannot move one arm to the floor and the rest to the ceiling.",
    "Regulons are in-silico set memberships (ChIP-Atlas, DoRothEA, MSigDB, Gray CHEA intersected with MitoCarta), not binding measured in these mice.",
    "n = 12 wild-type animals. Ranking plus a negative, not confirmatory inference."),
  source = c("results/biogenesis_axis_developmental.rds (script 47 PART B)",
             "reproduces results/collapse_module_ownership.rds $core_decomp (script 44)"))

save_panel_p(p, "biogax_regulon_split", height = 62)

if (FALSE) {
  print(p)
  sp |> print(n = 25)
  ## the four testable regulons, the reading in one line each
  ba$split_verdict |> print()
}
