# scripts/00_setup_packages.R

# List of CRAN packages
cran_packages <- c(
  "dplyr", "tibble", "readr", "stringr", "purrr", "magrittr",
  "ggplot2", "ggstatsplot", "pheatmap", "RColorBrewer",
  "grid", "gridExtra"
)

# List of Bioconductor packages
bioc_packages <- c(
  "DESeq2", "IHW", "ashr", "apeglm", "biomaRt",
  "ComplexHeatmap", "circlize",
  "PoiClaClu", "limma", "edgeR", "variancePartition"
)

# Load or install CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Load or install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
}
