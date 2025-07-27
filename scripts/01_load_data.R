# scripts/01_load_data.R

source("scripts/00_setup_packages.R")

# === Load coldata and reshape ===
coldata <- read.csv("data/coldata.csv", row.names = 1)

coldata$timepoint <- as.factor(coldata$timepoint) %>% relevel('6W')
coldata$myc_status <- as.factor(coldata$myc_status) %>% relevel('neg')
coldata$group <- factor(coldata$group, levels = c('6W_neg', '6W_pos', '12W_neg', '12W_pos'))
coldata <- coldata[order(coldata$group), ]

# === Load count matrix ===
cts <- read.csv("data/FULL.DAT.csv")
cts <- cts[-2]  # remove the second column
cts <- cts %>%
  remove_rownames() %>%
  column_to_rownames(var = "gene_id") %>%
  as.matrix() %>%
  round()
storage.mode(cts) <- "integer"
cts <- cts[, rownames(coldata)]  # match sample order

# === Create DESeqDataSet with interaction design ===
ddsMat_int <- DESeqDataSetFromMatrix(countData = cts,
                                     colData = coldata,
                                     design = ~ timepoint * myc_status)

# Filter genes with at least 10 counts in ≥4 samples
keep <- rowSums(counts(ddsMat_int) >= 10) >= 4
dds_int <- ddsMat_int[keep, ]

# === Load gene sets ===

# Mitocarta
# Convert to named list: pathway → gene vector. Custom selection of gene groups covering all mitocarta3.0
mitocarta_df <- read.csv("data/mitocarta_pathways.csv", header = FALSE, col.names = c("pathway", "symbol_string"))

mitocarta_sets <- setNames(
  lapply(strsplit(mitocarta_df$symbol_string, ","), str_trim),
  mitocarta_df$pathway
)

names(mitocarta_sets) <- paste0("MC_", names(mitocarta_sets))

# MYC .gmx
# Myc signature genes. Using different sets, compiled in Flesher paper: https://doi.org/10.1038/s41388-022-02458-9
# https://github.com/Yenaled/felsher/blob/master/genesets/myc_signature_genesets.gmx

# Read .gmx as column-oriented table
gmx_table <- read.delim("data/myc_signature_genesets.gmx", header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)

# Convert columns into named list of gene sets
gmx_sets <- apply(gmx_table, 2, function(col) {
  col <- na.omit(col[col != ""])  # remove empty entries
  col[-(1:2)]  # remove set name and description
})
names(gmx_sets) <- as.character(gmx_table[1, ])


# Connect to Ensembl mouse mart
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

ortholog_table <- getBM(
  attributes = c(
    "ensembl_gene_id",  # Mouse Ensembl ID
    "external_gene_name",  # Mouse gene symbol
    "hsapiens_homolog_associated_gene_name",  # Human gene symbol
    "hsapiens_homolog_orthology_type"
  ),
  mart = mouse
)

# Clean and filter
ortholog_table <- ortholog_table %>%
  filter(
    external_gene_name != "" & 
      hsapiens_homolog_associated_gene_name != "" &
      hsapiens_homolog_orthology_type %in% c("ortholog_one2one", "ortholog_one2many")
  )

# Create mapping from human → mouse (by gene symbol)
human_to_mouse_map <- deframe(
  ortholog_table[, c("hsapiens_homolog_associated_gene_name", "external_gene_name")]
)

myc_signature_sets <- lapply(gmx_sets, function(genes) {
  # Use the named vector to look up each human gene
  mapped <- human_to_mouse_map[genes]
  unique(na.omit(mapped))
})

# Prefix names for clarity
names(myc_signature_sets) <- paste0("MYC_", names(myc_signature_sets))

# Felsher integrative myc gene set
felsher_human_genes <- read.csv("data/felsher_integrative_signature.csv")$Symbol
felsher_human_genes <- unique(na.omit(felsher_human_genes))

felsher_mouse_genes <- human_to_mouse_map[felsher_human_genes]
felsher_mouse_genes <- unique(na.omit(felsher_mouse_genes))

# Combine into one list
gene_sets_list <- c(mitocarta_sets, myc_signature_sets)

gene_sets_list[["MYC_felsher_integrative_signature"]] <- felsher_mouse_genes


# === Save processed objects ===
saveRDS(cts, "results/count_matrix.rds")
saveRDS(coldata, "results/coldata.rds")
saveRDS(dds_int, "results/dds_int.rds")
saveRDS(gene_sets_list, "results/gene_sets_list.rds")
saveRDS(ortholog_table, "results/ortholog_table.rds")
