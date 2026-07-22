# scripts/00_setup_packages.R

# List of CRAN packages
cran_packages <- c(
  "here",
  "dplyr", "tibble", "tidyr", "readr", "stringr", "purrr", "magrittr", "forcats",
  "ggplot2", "ggstatsplot", "pheatmap", "RColorBrewer",
  "grid", "gridExtra", "readr",
  "reshape2", "ggrepel",
  "ggbeeswarm",      # quasirandom points for n=6/group panels (figures/)
  "patchwork",       # multi-panel figure composition (figures/)
  "msigdbr",
  "readxl",          # MitoCarta3.0 .xls parsing (08_mitoPPS_analysis.R)
  "splitstackshape"  # cSplit for gene-to-pathway mapping (08_mitoPPS_analysis.R)
)

# List of Bioconductor packages
bioc_packages <- c(
  "DESeq2", "IHW", "ashr", "apeglm", "biomaRt",
  "ComplexHeatmap", "circlize",
  "PoiClaClu", "limma", "edgeR", "variancePartition",
  "fgsea", "gprofiler2",
  "GSVA",            # per-sample gene-set scoring (15_gsva_scoring.R)
  "singscore"        # rank-based cohort-independent set scoring (36_linear_pathway_coupling.R)
)

# Load or install CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load or install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load all packages
library(here)
for (pkg in cran_packages) library(pkg, character.only = TRUE)
for (pkg in bioc_packages) library(pkg, character.only = TRUE)
