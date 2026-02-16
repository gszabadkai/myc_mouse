# scripts/01_load_data.R


source(here::here("scripts", "00_setup_packages.R"))

# === Load coldata and reshape ===
coldata <- read.csv(here("data", "coldata.csv"), row.names = 1)

coldata$timepoint <- as.factor(coldata$timepoint) %>% relevel('6W')
coldata$myc_status <- as.factor(coldata$myc_status) %>% relevel('neg')
coldata$group <- factor(coldata$group, levels = c('6W_neg', '6W_pos', '12W_neg', '12W_pos'))
coldata <- coldata[order(coldata$group), ]

# === Load count matrix ===
cts <- read.csv(here("data", "FULL.DAT.csv"))
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
mitocarta_df <- read.csv(here("data", "mitocarta_pathways.csv"), header = FALSE, col.names = c("pathway", "symbol_string"))

mitocarta_sets <- setNames(
  lapply(strsplit(mitocarta_df$symbol_string, ","), str_trim),
  mitocarta_df$pathway
)

names(mitocarta_sets) <- paste0("MC_", names(mitocarta_sets))

# MYC .gmx
# Myc signature genes. Using different sets, compiled in Flesher paper: https://doi.org/10.1038/s41388-022-02458-9
# https://github.com/Yenaled/felsher/blob/master/genesets/myc_signature_genesets.gmx

# Read .gmx as column-oriented table
gmx_table <- read.delim(here("data", "myc_signature_genesets.gmx"), header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)

# Convert columns into named list of gene sets
gmx_sets <- apply(gmx_table, 2, function(col) {
  col <- na.omit(col[col != ""])  # remove empty entries
  col[-(1:2)]  # remove set name and description
})
names(gmx_sets) <- as.character(gmx_table[1, ])


# === Ortholog mapping (cached) ===
ortholog_cache_path <- here("results", "ortholog_table.rds")

if (file.exists(ortholog_cache_path)) {
  message("Loading cached ortholog table...")
  ortholog_table <- readRDS(ortholog_cache_path)
} else {
  message("Querying biomaRt for ortholog table (this may take a few minutes)...")
  
  # Connect to Ensembl mouse mart
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  ortholog_table <- getBM(
    attributes = c(
      "ensembl_gene_id",
      "external_gene_name",
      "hsapiens_homolog_associated_gene_name",
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
  
  # Save cache for future runs
  saveRDS(ortholog_table, ortholog_cache_path)
}

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
felsher_human_genes <- read.csv(here("data", "felsher_integrative_signature.csv"))$Symbol
felsher_human_genes <- unique(na.omit(felsher_human_genes))

felsher_mouse_genes <- human_to_mouse_map[felsher_human_genes]
felsher_mouse_genes <- unique(na.omit(felsher_mouse_genes))

# Combine into one list
gene_sets_list <- c(mitocarta_sets, myc_signature_sets)

gene_sets_list[["MYC_felsher_integrative_signature"]] <- felsher_mouse_genes

# === Custom apoptosis gene sets (Mitocarta 3.0 -derived) ===

# Lists
apoptosis_pro_symbols <- c(
  "Aifm1","Aifm3","Aifm2","Bad","Bak1","Bax","Bbc3","Bcl2l11","Bid","Bik",
  "Bnip3","Bnip3l","Bok","Casp3","Casp8","Casp9","Cycs","Diablo","Endog",
  "Htra2","Ifi27","Pmaip1","Septin4","Bcl2l13","Sphk2"
)

apoptosis_anti_symbols <- c(
  "Bcl2","Bcl2a1d","Bcl2l1","Bcl2l10","Bcl2l2","Mcl1","Ghitm","Styxl1","Chchd2"
)


# Add to gene_sets_list
gene_sets_list[["MC_Apoptosis_Pro"]]  <- apoptosis_pro_symbols
gene_sets_list[["MC_Apoptosis_Anti"]] <- apoptosis_anti_symbols

# === MSigDB Hallmark pathways ===
# Load Hallmark gene sets from msigdbr package
# These are human gene sets - we convert to mouse orthologs

hallmark_df <- msigdbr(species = "Mus musculus", category = "H")

hallmark_list <- hallmark_df |>
  dplyr::select(gs_name, gene_symbol) |>
  group_by(gs_name) |>
  summarise(genes = list(unique(gene_symbol)), .groups = "drop") |>
  deframe()

# Prefix names for clarity
names(hallmark_list) <- paste0("MSigDB_", names(hallmark_list))

# Combine all gene sets (89 total: 22 MitoCarta + 17 MYC + 50 Hallmark)
gene_sets_list <- c(gene_sets_list, hallmark_list)

message("Gene sets loaded: ", length(gene_sets_list), " total")
message("  MitoCarta: ", sum(grepl("^MC_", names(gene_sets_list))))
message("  MYC signatures: ", sum(grepl("^MYC_", names(gene_sets_list))))
message("  MSigDB Hallmark: ", sum(grepl("^MSigDB_", names(gene_sets_list))))

# === Save processed objects ===
saveRDS(cts, here("results", "count_matrix.rds"))
saveRDS(coldata, here("results", "coldata.rds"))
saveRDS(dds_int, here("results", "dds_int.rds"))
saveRDS(gene_sets_list, here("results", "gene_sets_list.rds"))
