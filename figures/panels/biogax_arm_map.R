# =============================================================================
# biogax_arm_map.R -- the gland withdraws from the chain, not from the organelle
# -----------------------------------------------------------------------------
# DISCUSSION PANEL (see biogax_regulon_split.R for why the `biogax_` prefix).
#
# THE ARGUMENT. A withdrawal of mitochondrial BIOGENESIS predicts that the
# respiratory subunits, the factors that assemble them, the mitochondrial
# ribosome and the import channels all move together -- they are one programme.
# What the wild-type gland actually does is move the structural subunits of the
# proton-pumping complexes and nothing else. The assembly factors that build
# those very complexes sit at the centre of their own null.
#
# WHAT IS DRAWN: every MitoCarta arm as a point, positioned by the percentile of
# its wild-type 6-to-12-week change against 2000 expression-matched random sets,
# grouped by what the arm DOES. Class assignment is script 47's `arm_class()`,
# fixed in code and stored in the object, so it can be audited rather than
# believed.
#
# THE mtDNA CLASS IS DRAWN OPEN, AND ON PURPOSE. The thirteen mtDNA-encoded
# subunits sit at the top of the panel, which is the opposite direction to the
# nuclear chain. That axis is three-way confounded in this dataset (mitochondrial
# content x proliferation denominator x prep leak; CLAUDE.md), so it is shown
# because hiding it would be worse, and marked so that no argument can rest on it.
#
# Reads : results/biogenesis_axis_developmental.rds (script 47 PART A)
# Output: outputs/figures/panels/biogax_arm_map.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

ba_path <- here::here("results", "biogenesis_axis_developmental.rds")
if (!file.exists(ba_path)) stop("run scripts/47_... first")
ba <- readRDS(ba_path)

am <- as.data.frame(ba$armmap)
stopifnot(all(c("arm", "class", "n_genes", "c_wt_time", "wt_null_pct") %in% names(am)))

CLASS_LAB <- c(
  chain_structural = "respiratory chain\nstructural subunits",
  chain_assembly   = "chain\nassembly factors",
  organelle_build  = "build and maintain\nthe organelle",
  metabolic        = "mitochondrial\nmetabolism",
  other            = "everything else",
  chain_mtDNA      = "mtDNA-encoded\n(confounded)")
am <- am[am$class %in% names(CLASS_LAB), ]
am$class <- factor(am$class, levels = rev(names(CLASS_LAB)), labels = rev(CLASS_LAB))

# Arms worth naming, with their labels HAND-PLACED (fig2_priming_ratios' idiom):
# several sit in the same class row, so an automatic left/right rule collides
# them. dy moves the label off the point cloud, which the jitter keeps inside
# +/- 0.16 of the row.
LAB <- tibble::tribble(
  ~arm,                                    ~short,                 ~dy,   ~hj,  ~dx,
  "MITOCARTA_OXPHOS_SUBUNITS",             "OXPHOS subunits",     -0.34,  0,     1.5,
  "MITOCARTA_CIV_SUBUNITS",                "complex IV subunits",  0.34,  0,     1.5,
  "MITOCARTA_COMPLEX_II",                  "complex II",           0.34,  0,     2.0,
  "MITOCARTA_ELECTRON_CARRIERS",           "electron carriers",   -0.34,  1,    -2.0,
  "MITOCARTA_OXPHOS_ASSEMBLY_FACTORS",     "OXPHOS assembly",     -0.34,  1,    -2.0,
  "MITOCARTA_MITOCHONDRIAL_RIBOSOME",      "mitoribosome",        -0.34,  0,     2.0,
  "MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA", "central dogma",        0.34,  1,    -2.0,
  "MITOCARTA_FATTY_ACID_OXIDATION",        "fatty-acid oxidation", 0.34,  1,    -2.0)
stopifnot(all(LAB$arm %in% am$arm))
lab <- merge(LAB, am[, c("arm", "class", "wt_null_pct", "c_wt_time", "n_genes")],
             by = "arm")
# merge() reorders; position_nudge() takes vectors POSITIONALLY, so the frame has
# to be put back into the order the nudges were written in or the labels land on
# the wrong points.
lab <- lab[order(match(lab$arm, LAB$arm)), ]
stopifnot(identical(lab$arm, LAB$arm))
# Named points take the declared verdict palette by what they actually did, so a
# riser is never inked as a withdrawal.
lab$verdict <- ifelse(lab$wt_null_pct < 5, "withdraws",
                      ifelse(lab$wt_null_pct > 95, "rises", "at chance"))

# ASSERTIONS -- the two readings the panel exists to make.
cls_med <- tapply(am$wt_null_pct, as.character(am$class), stats::median)
stopifnot(cls_med[[CLASS_LAB[["chain_structural"]]]] < 10,
          cls_med[[CLASS_LAB[["chain_assembly"]]]]   > 40,
          cls_med[[CLASS_LAB[["chain_assembly"]]]]   < 60)

set.seed(7)
p <- ggplot2::ggplot(am, ggplot2::aes(x = wt_null_pct, y = class)) +
  ggplot2::annotate("rect", xmin = 5, xmax = 95, ymin = -Inf, ymax = Inf,
                    fill = "grey96", colour = NA) +
  ggplot2::geom_vline(xintercept = 50, linewidth = 0.22, colour = "grey65") +
  ggplot2::geom_jitter(
    ggplot2::aes(size = n_genes, shape = grepl("confounded", as.character(class))),
    height = 0.16, width = 0, stroke = 0.3,
    fill = "grey30", colour = "grey25", alpha = 0.75) +
  ggplot2::geom_point(data = lab,
                      ggplot2::aes(x = wt_null_pct, y = class, fill = verdict),
                      inherit.aes = FALSE, shape = 21, size = 1.9,
                      colour = "white", stroke = 0.35, show.legend = FALSE) +
  ggplot2::geom_text(data = lab,
                     ggplot2::aes(x = wt_null_pct, y = class, label = short,
                                  hjust = hj),
                     inherit.aes = FALSE, size = 1.8, colour = "grey15",
                     position = ggplot2::position_nudge(x = lab$dx, y = lab$dy)) +
  ggplot2::scale_fill_manual(values = verdict_cols, guide = "none") +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 1), guide = "none") +
  ggplot2::scale_size_area(max_size = 2.6, breaks = c(10, 50, 150),
                           name = "genes in arm") +
  ggplot2::scale_x_continuous(
    limits = c(-16, 116), breaks = c(0, 5, 25, 50, 75, 95), expand = c(0, 0),
    name = "percentile of 2000 expression-matched random sets") +
  ggplot2::scale_y_discrete(name = NULL) +
  theme_panel() +
  ggplot2::theme(legend.position = "bottom",
                 legend.margin = ggplot2::margin(-3, 0, 0, 0),
                 panel.grid.major.y = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(lineheight = 0.95))

pick <- function(a, col) am[[col]][am$arm == a]
LEGEND <- panel_legend(
  slot = "Discussion D2",
  what = paste(
    "Every MitoCarta arm placed by the percentile of its wild-type 6-to-12-week",
    "change against 2000 expression-matched random sets, grouped by function.",
    "Point area is the number of genes in the arm."),
  detail = c(
    sprintf("Class medians: structural subunits %.1f pct, assembly factors %.1f, organelle build %.1f, metabolism %.1f.",
            cls_med[[CLASS_LAB[["chain_structural"]]]], cls_med[[CLASS_LAB[["chain_assembly"]]]],
            cls_med[[CLASS_LAB[["organelle_build"]]]], cls_med[[CLASS_LAB[["metabolic"]]]]),
    sprintf("Complex IV subunits lead the withdrawal (%+.3f log2, %.1f pct); the 87 OXPHOS subunits as a whole are %+.3f (%.1f pct).",
            pick("MITOCARTA_CIV_SUBUNITS", "c_wt_time"), pick("MITOCARTA_CIV_SUBUNITS", "wt_null_pct"),
            pick("MITOCARTA_OXPHOS_SUBUNITS", "c_wt_time"), pick("MITOCARTA_OXPHOS_SUBUNITS", "wt_null_pct")),
    sprintf("The sharp internal control: OXPHOS assembly factors %+.4f log2, %.1f pct -- the centre of their own null, over %d genes.",
            pick("MITOCARTA_OXPHOS_ASSEMBLY_FACTORS", "c_wt_time"),
            pick("MITOCARTA_OXPHOS_ASSEMBLY_FACTORS", "wt_null_pct"),
            pick("MITOCARTA_OXPHOS_ASSEMBLY_FACTORS", "n_genes")),
    sprintf("Two arms of the chain do NOT withdraw: complex II (%+.3f, %.1f pct), the one complex that neither pumps protons nor carries an mtDNA-encoded subunit, and the electron carriers (%+.3f, %.1f pct).",
            pick("MITOCARTA_COMPLEX_II", "c_wt_time"), pick("MITOCARTA_COMPLEX_II", "wt_null_pct"),
            pick("MITOCARTA_ELECTRON_CARRIERS", "c_wt_time"), pick("MITOCARTA_ELECTRON_CARRIERS", "wt_null_pct"))),
  bounds = c(
    "BATCH = TIMEPOINT: the magnitude of the developmental change is confounded with cohort. The ARM-SELECTIVITY is not -- a shared batch effect cannot put one class at the floor and another at the centre.",
    "The mtDNA-encoded class (open symbols) is three-way confounded and is shown for completeness only; no argument in this document rests on it.",
    "Two arms of the build class are low rather than central (chaperones ~9.5 pct, protein import ~8.7 pct). The claim is that the build arms do not FOLLOW the chain, not that nothing else moves.",
    "n = 12 wild-type animals."),
  source = "results/biogenesis_axis_developmental.rds (script 47 PART A)")

save_panel_p(p, "biogax_arm_map", height = 68)

if (FALSE) {
  print(p)
  ba$arm_class_summary |> print()
  am[order(am$wt_null_pct), c("arm", "class", "n_genes", "c_wt_time", "wt_null_pct")] |>
    head(20) |> print()
}
