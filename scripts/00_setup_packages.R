# scripts/00_setup_packages.R

# List of CRAN packages
cran_packages <- c(
  "here",
  "dplyr", "tibble", "tidyr", "readr", "stringr", "purrr", "magrittr", "forcats",
  "ggplot2", "ggstatsplot", "pheatmap", "RColorBrewer",
  "grid", "gridExtra", "readr",
  "reshape2", "ggrepel",
  "msigdbr"
)

# List of Bioconductor packages
bioc_packages <- c(
  "DESeq2", "IHW", "ashr", "apeglm", "biomaRt",
  "ComplexHeatmap", "circlize",
  "PoiClaClu", "limma", "edgeR", "variancePartition",
  "fgsea", "gprofiler2"               
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
