# =============================================================================
# edfig2b_geneset_library.R -- Extended Data Fig. 2b: the gene-set library, by
# category and by how mitochondrial each set is
# -----------------------------------------------------------------------------
# The text: "... we profiled MECs purified from 6- and 12-week glands ... (ED
# Fig. 2a-c)"; the library is what the pathway axis of ED 2c is scored on. From
# the old Fig. S1B (figures/panels/figS1_geneset_library.R), redrawn here.
#
# Every set in the v1.0 snapshot (data/genesets_from_library/, never rebuilt),
# counted by primary category and filled by its mitochondrial class:
#   MitoCarta                  defined on MitoCarta (or a respiratory metabolism set)
#   Curated (with mitochondria) a transcription-factor target set intersected with
#                              MitoCarta: mitochondrial BY CONSTRUCTION, so it cannot
#                              count as independent mitochondrial coverage
#   Not mitochondrial          everything else
# The classes are script 37's rule (37:445-462), re-derived over the whole library
# and asserted identical to script 37's saved classification on the sets they
# share. Colours are mito_class_cols: purple for the two mitochondrial classes,
# grey for the rest (2026-10-14); blue and orange mean genotype.
#
# FORM: counts per category, three parts -> stacked bars, sorted (data-to-viz,
#   "Order your data"); the key sits inside, bottom right, where the short bars
#   leave room; its middle entry is broken over two lines to fit beside the counts.
#
# Reads:  data/genesets_from_library/provenance_table.csv (v1.0 snapshot)
#         results/pathway_loading.rds (script 37: $mito_classification, assertion)
#         results/gsva_scores.rds     (script 15: $n_sets_scored, legend)
# Output: outputs/natmetab/EDFig2/EDFig2b_geneset_library.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

prov <- utils::read.csv(here::here("data", "genesets_from_library", "provenance_table.csv"),
                        stringsAsFactors = FALSE)
stopifnot(all(c("set_name", "category_primary", "method", "size_mouse", "fgsea_eligible")
              %in% names(prov)))

# script 37's rule (37:445-462), verbatim; `category` is the library's category_primary
mito_class3 <- function(set, category) {
  category  <- ifelse(is.na(category), "", category)
  name_mito <- grepl("_MITO$|_MITO_|MITO_NU|^MITO_|CORE_MITO", set)
  mito_defined <- category == "MitoCarta" | name_mito |
    (category == "Metabolism" & grepl("OXPHOS|KREBS|TCA|ELECTRON|RESPIRAT", set))
  ifelse(category == "MitoCarta", "mitocarta_proper",
    ifelse(name_mito, "construction_MITO",
      ifelse(mito_defined, "mitocarta_proper", "non_mito")))
}
lib <- data.frame(set = prov$set_name, category = prov$category_primary,
                  stringsAsFactors = FALSE)
lib$class3 <- mito_class3(lib$set, lib$category)

ref <- as.data.frame(readRDS(here::here("results", "pathway_loading.rds"))$mito_classification)
chk <- merge(lib[, c("set", "class3")], ref[, c("set", "class3")], by = "set",
             suffixes = c("_mine", "_37"))
stopifnot(nrow(chk) == nrow(ref), identical(chk$class3_mine, chk$class3_37))

cat_lab <- c(TF_targets = "TF regulons", Mammary_development = "Mammary development",
             MitoCarta = "MitoCarta", Metabolism = "Metabolism", Apoptosis = "Apoptosis",
             Biogenesis_apoptosis_intersections = "Biogenesis ×\napoptosis",
             MYC_signatures = "MYC signatures",
             Biogenesis_discrimination = "Biogenesis\ndiscrimination",
             Proliferation = "Proliferation")
stopifnot(setequal(names(cat_lab), unique(lib$category)))

tab <- as.data.frame(table(category = lib$category, class3 = lib$class3), stringsAsFactors = FALSE)
tab <- tab[tab$Freq > 0, ]
tot <- tapply(tab$Freq, tab$category, sum)
stopifnot(sum(tot) == nrow(prov))
lev <- names(sort(tot))                                    # largest at the top
tab$category <- factor(cat_lab[tab$category], levels = cat_lab[lev])
tab$class3   <- factor(tab$class3, levels = names(mito_class_cols))
tots <- data.frame(category = factor(cat_lab[names(tot)], levels = cat_lab[lev]), n = as.vector(tot))

p <- ggplot(tab, aes(Freq, category, fill = class3)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.15,
           position = position_stack(reverse = TRUE)) +
  geom_text(data = tots, aes(n, category, label = n), inherit.aes = FALSE, hjust = -0.25,
            size = pt(NM_TXT), family = NM_FONT) +
  scale_fill_manual(values = mito_class_cols, name = NULL,
                    labels = replace(mito_class_labels, "construction_MITO",
                                     "Curated (with\nmitochondria)")) +   # narrower key
  scale_x_nm(step = 200, top = max(tot) * 1.2, minor = TRUE, labels = lab_num) +  # room for "444"
  labs(x = "Gene sets", y = NULL) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_nm(legend = "inside") +
  theme(legend.position.inside = c(1, 0.03), legend.justification.inside = c(1, 0),
        legend.text = element_text(size = NM_TXT, lineheight = 0.85),
        legend.key.size = unit(2.4, "mm"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(lineheight = 0.8))

save_panel(p, fig = "EDFig2", panel = "b", name = "geneset_library", width = 58, height = 40)

# --- numbers for the legend (printed, never drawn) --------------------------------
n_cls <- table(lib$class3)
cat("\nED Fig. 2b -- for the legend\n")
cat(sprintf("  mammary_geneset_library v1.0: %d sets in %d categories; set size median %d genes (range %d-%d)\n",
            nrow(prov), length(tot), as.integer(stats::median(prov$size_mouse)),
            min(prov$size_mouse), max(prov$size_mouse)))
cat(sprintf("  MitoCarta %d; curated (TF target sets intersected with MitoCarta) %d; not mitochondrial %d\n",
            n_cls[["mitocarta_proper"]], n_cls[["construction_MITO"]], n_cls[["non_mito"]]))
cat(sprintf("  scored per sample by GSVA: %d (results/gsva_scores.rds); fGSEA-eligible: %d\n",
            readRDS(here::here("results", "gsva_scores.rds"))$n_sets_scored, sum(prov$fgsea_eligible)))
cat("  categories overlap by construction, so the counts describe the library, not independent tests\n")

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p)
  tab[order(tab$category), ]
}
