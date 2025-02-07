# for details of choices of some steps see the workflow: https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#introduction


rm(list=ls())

library(DESeq2)
library(dplyr)
library(tibble)
library(ggstatsplot)
library(ggplot2)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(biomaRt)
library(clipr)
library(variancePartition)
library(edgeR)
library(stringr)
library(IHW)
library(gprofiler2)
library(purrr)
library(tidyr)
library(UpSetR)
library(openxlsx)
library(RCy3)



###### Read in and set up coldata ----

coldata <- read.delim("data/FULL.DAT.COL.DATA.txt", sep = "\t", header = TRUE)

#  reshape data, add groupings to detect source of variability - see below
coldata <- coldata %>% remove_rownames %>% column_to_rownames(var="count.file") %>% 
  mutate(myc_status = case_when(
    group == "hyperplastic" ~ "pos",
    group == "WT" ~ "neg"
  )) %>% 
  mutate(timepoint = case_when(
    age == '6 wk' ~ '6W',
    age == '12 wk' ~'12W'
  )) %>% 
  mutate(group = paste(timepoint, myc_status, sep = '_')) %>% 
  mutate(prefix = sub("_.*", "", sample)) %>% 
  mutate(suffix = sub(".*_", "", sample)) %>% 
  mutate(numeric_part = str_extract(suffix, "\\d+")) %>% 
  mutate(alphabetical_part = str_extract(suffix, "[a-zA-Z]+")) %>% 
  mutate(across(c(5, 7:10), as.factor))
coldata <- coldata[-c(3,4)]
coldata$timepoint <- as.factor(coldata$timepoint) %>% relevel('6W')
coldata$group <- factor(coldata$group, levels = c('6W_neg', '6W_pos', '12W_neg', '12W_pos'))
coldata <- coldata[order(coldata$group), ]


###### Read in and set up count file ----

cts <- read.csv("data/FULL.DAT.csv")
cts <- cts[-2]
cts <- cts %>% remove_rownames %>% column_to_rownames(var="gene_id")
cts <- as.matrix(round(cts))
storage.mode(cts) <- "integer"
cts <- cts[, rownames(coldata)]


# Check the matches here before loading the data for further analysis
all(rownames(coldata) %in% colnames(cts)) # This should return TRUE

# look at the total counts
round(colSums(cts) / 1e6, 1 )

# make graph of total counts in different conditions
total_cts_graph_data <- coldata
total_cts_graph_data$tcts <- round(colSums(cts) / 1e6, 1 )

# total count stats
total_cts_graph_data <- total_cts_graph_data %>% 
  mutate(
  Exp_groups = paste(myc_status, timepoint, sep = "_")
  ) %>% 
  mutate(Exp_groups = factor(Exp_groups, levels = c("neg_6W", "pos_6W", "neg_12W", "pos_12W")))

ggbetweenstats(
  data = total_cts_graph_data,
  x = Exp_groups,         
  y = tcts,             
  type = "parametric",  
  pairwise.comparisons = TRUE,  
  pairwise.display = "significant",
  title = "total counts",
  xlab = "myc status ~ timepoint",
  ylab = "total counts (M)"
)

<<<<<<< HEAD
##### create a DESeq object w interaction -----
=======
##### create a DESeq object -----
>>>>>>> e51bf34e25f3ac64d2e1b39c33044662d47df567

ddsMat <- DESeqDataSetFromMatrix(countData = cts,
                                 colData = coldata,
                                 design = ~ timepoint * myc_status) #this expands to ~ timepoint + myc_status + timepoint:myc_status

## filter by counts: 55K rows --> 19K

nrow(ddsMat)

smallestGroupSize <- 4
keep <- rowSums(counts(ddsMat) >= 10) >= smallestGroupSize
dds <- ddsMat[keep,]

nrow(dds)

## normalise and transform to look at variance

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

##### dds data QC -----

total_cts_graph_data$size_factors <- sizeFactors(dds)

ggbetweenstats(
  data = total_cts_graph_data,
  x = Exp_groups,         
  y = size_factors,             
  type = "parametric",  
  pairwise.comparisons = TRUE,  
  pairwise.display = "significant",
  title = "size factors",
  xlab = "myc status ~ timepoint",
  ylab = "size factors"
)

## variance stabilizing transformations 

vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

df_var <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df_var)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df_var$transformation <- factor(df_var$transformation, levels=lvls)

ggplot(df_var, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

## Sample distances
# by vst
sampleDists_vsd <- dist(t(assay(vsd)))
sampleDists_vsd

sampleDistMatrix_vsd <- as.matrix( sampleDists_vsd )
rownames(sampleDistMatrix_vsd) <- paste(vsd$group, vsd$sample, sep = " / " )
colnames(sampleDistMatrix_vsd) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_vsd,
         clustering_distance_rows = sampleDists_vsd,
         clustering_distance_cols = sampleDists_vsd,
         col = colors)

# by rlog
sampleDists_rld <- dist(t(assay(rld)))
sampleDists_rld

sampleDistMatrix_rld <- as.matrix( sampleDists_rld )
rownames(sampleDistMatrix_rld) <- paste(rld$group, rld$sample, sep = " / " )
colnames(sampleDistMatrix_rld) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_rld,
         clustering_distance_rows = sampleDists_rld,
         clustering_distance_cols = sampleDists_rld,
         col = colors)

# by Poisson Distance

poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$group, dds$sample, sep=" / " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# overall there is no clustering by sample groups

## PCA

plotPCA(vsd, intgroup = c("group", "sample"))

# no clear grouping by condition!! also, no batch effect is visible, all seems rather random
## which genes contribute to PCA1 and PCA2 (repeat with changing the numbers), 

pca <- prcomp(t(assay(vsd)))  # Perform PCA on the same data
pc1_loadings <- pca$rotation[, 1]  # Extract PC2 loadings
top_genes <- head(order(abs(pc1_loadings), decreasing = TRUE), n = 50)  # Top 50 contributors

pc1_top_genes <- rownames(pca$rotation)[top_genes]

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
annotation <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = pc1_top_genes,
  mart = ensembl
)

annotation$external_gene_name

write_clip(annotation$external_gene_name)

# examining the resulted gene sets with Gprofiler gives mostly immune/lymphocyte related terms see: https://biit.cs.ut.ee/gprofiler/gost?organism=mmusculus&query=Il7r%0ACd247%0APax5%0AIkzf3%0ACyfip2%0AItk%0AGpr132%0ASatb1%0AMs4a1%0AMs4a6b%0ACd28%0APtprc%0ASell%0ACr2%0ACd2%0ALef1%0ASh2d2a%0ARhoh%0ASt8sia1%0ACd27%0ACd19%0AIl21r%0ACd3e%0ADnah8%0AGalnt6%0ATraf3ip3%0ABank1%0ACcr7%0AGimap3%0ACd53%0AUbash3a%0AGrap2%0ACxcr5%0ASelplg%0AThemis%0AArhgap15%0AP2ry10%0ASkap1%0AGvin3%0AFam169b%0ATrbc2%0ATrac%0AC230085N15Rik%0A5830444F18Rik%0AE430014B02Rik%0AGm37248%0AGm36931%0AIghd%0AA130071D04Rik%0AGm49553&ordered=false&all_results=false&no_iea=false&combined=false&measure_underrepresentation=false&domain_scope=annotated&significance_threshold_method=g_SCS&user_threshold=0.05&numeric_namespace=ENTREZGENE_ACC&sources=GO:MF,GO:CC,GO:BP,KEGG,TF,REAC,MIRNA,HPA,CORUM,HP,WP&background=&highlight=true&no_evidences=false
# this might be a random contamination???

# further exploration with variancePartition
# Prepare count data and metadata
countData <- counts(dds)
metadata <- colData(dds)

# Normalize and transform data
dge <- DGEList(counts = countData)
dge <- calcNormFactors(dge)
v <- voom(dge, design = NULL)

# define formula
formula <- ~ (1|group/myc_status) 

(1|myc_status) + (1|timepoint) 

+ (1|numeric_part) + (1|alphabetical_part) + (1|myc_status) + (1|timepoint) (1|prefix) (1|suffix)

# Fit variance model
varPart <- fitExtractVarPartModel(v$E, formula, metadata)

# Plot variance partitioning results
plotVarPart(varPart)

# only samll part of the variation is explained by any of the parameters - moving on for DGE

###### DGE results / contrasts / lfcShrink /IHW -----

dds <- DESeq(dds)

resultsNames(dds)

#### 1. main effect of timepoints (myc_neg) -----

timepoint_effect_12W_vs_6W_at_myc_neg <- results(dds, name = "timepoint_12W_vs_6W")
timepoint_effect_12W_vs_6W_at_myc_neg.ape <- lfcShrink(dds, coef="timepoint_12W_vs_6W", type="apeglm")
timepoint_effect_12W_vs_6W_at_myc_neg.IHW <- results(dds, name = "timepoint_12W_vs_6W", filterFun=ihw)
# final result with lfcShrink/apeglm and IHW:
timepoint_effect_12W_vs_6W_at_myc_neg.IHW.ape <- lfcShrink(dds, coef="timepoint_12W_vs_6W", type="apeglm", res = timepoint_effect_12W_vs_6W_at_myc_neg.IHW)

# summaries and MA plot
summary(timepoint_effect_12W_vs_6W_at_myc_neg)
summary(timepoint_effect_12W_vs_6W_at_myc_neg.ape) 
DESeq2::plotMA(timepoint_effect_12W_vs_6W_at_myc_neg.ape, ylim=c(-2,2))
DESeq2::plotMA(timepoint_effect_12W_vs_6W_at_myc_neg.IHW.ape, ylim=c(-2,2))

# number of 10% FDR results is less with IHW --> noisy data???
# sum(timepoint_effect_12W_vs_6W_at_myc_neg.ape$padj < 0.1, na.rm = TRUE)
# [1] 1896
# > sum(timepoint_effect_12W_vs_6W_at_myc_neg.IHW.ape$padj < 0.1, na.rm = TRUE)
# [1] 1867

# checking ihw results - using basemean as covariate doesn't add anything --> use timepoint_effect_12W_vs_6W_at_myc_neg.ape in further analysis

res <- results(dds, name = "timepoint_12W_vs_6W", pAdjustMethod="none")
ihw_res <- ihw(
  pvalues    = res$pvalue,
  covariates = res$baseMean,
  alpha      = 0.1
)
res$IHW_padj <- adj_pvalues(ihw_res)

# > sum(res$pvalue < 0.1, na.rm = TRUE)
# [1] 5028
# > sum(res$IHW_padj < 0.1, na.rm = TRUE)
# [1] 1867

plot(ihw_res)
plot(ihw_res, what = "decisionboundary") 
gg <- ggplot(as.data.frame(ihw_res), aes(x = pvalue, y = adj_pvalue, col = group)) + 
  geom_point(size = 0.25) + scale_colour_hue(l = 70, c = 150, drop = FALSE)
gg %+% subset(as.data.frame(ihw_res), adj_pvalue <= 0.1)
  
#### 2. main effect of timepoints (myc_pos) ----

timepoint_effect_12W_vs_6W_at_myc_pos <- results(dds, contrast = list(c("timepoint_12W_vs_6W", "timepoint12W.myc_statuspos")))
    # timepoint_effect_12W_vs_6W_at_myc_pos.normal <- lfcShrink(dds, contrast = list(c("timepoint_12W_vs_6W", "timepoint12W.myc_statuspos")), type="normal") #ofc, lfcShrink work with interactions
summary(timepoint_effect_12W_vs_6W_at_myc_pos)
DESeq2::plotMA(timepoint_effect_12W_vs_6W_at_myc_pos, ylim=c(-2,2))

resMyc <- results(dds, contrast = list(c("timepoint_12W_vs_6W", "timepoint12W.myc_statuspos")), pAdjustMethod="none")
ihw_resMyc <- ihw(
  pvalues    = resMyc$pvalue,
  covariates = resMyc$baseMean,
  alpha      = 0.1
)
resMyc$IHW_padj <- adj_pvalues(ihw_resMyc)

# slight effect of IHW... no worth pursuing
# > sum(resMyc$pvalue < 0.1, na.rm = TRUE)
# [1] 5631
# > # [1] 5028
#   > sum(resMyc$IHW_padj < 0.1, na.rm = TRUE)
# [1] 2623
# > sum(timepoint_effect_12W_vs_6W_at_myc_pos$padj < 0.1, na.rm = TRUE)
# [1] 2636
# > 

plot(ihw_resMyc)
plot(ihw_resMyc, what = "decisionboundary")


##### 3. overall myc effect (6W) ---- 
# (actually it is the same as the 6W effect, due to the use of the interaction term see ?results examples 2 and 3)

myc_effect_overall <- results(dds, name = "myc_status_pos_vs_neg")
myc_effect_overall.ape <- lfcShrink(dds, coef = "myc_status_pos_vs_neg", type = "apeglm")
summary(myc_effect_overall)
summary(myc_effect_overall.ape)

# IHW
resMycO <- results(dds, name = "myc_status_pos_vs_neg", pAdjustMethod="none")
ihw_resMycO <- ihw(
  pvalues    = resMycO$pvalue,
  covariates = resMycO$baseMean,
  alpha      = 0.1
)
resMycO$IHW_padj <- adj_pvalues(ihw_resMycO)

# > sum(resMycO$pvalue < 0.1, na.rm = TRUE)
# [1] 5825
# > sum(resMycO$IHW_padj < 0.1, na.rm = TRUE)
# [1] 2777
# > sum(myc_effect_overall$padj < 0.1, na.rm = TRUE)
# [1] 2704

plot(ihw_resMycO)
plot(ihw_resMycO, what = "decisionboundary") 


# IHW is useful in this case!! integrating with lfcShrink:
myc_effect_overall.ape.IHW <- myc_effect_overall.ape
myc_effect_overall.ape.IHW$padj <- resMycO$IHW_padj

DESeq2::plotMA(myc_effect_overall.ape.IHW, ylim=c(-2,2))

# # myc effect at 6W 
# 
# myc_effect_at_6W <- results(dds, contrast = c("myc_status", "pos", "neg"))
# summary(myc_effect_at_6W)

###### 4. myc effect at 12W ----

myc_effect_at_12W <- results(dds, contrast = list(c("myc_status_pos_vs_neg", "timepoint12W.myc_statuspos")))
summary(myc_effect_at_12W)

resMyc12 <- results(dds, contrast = list(c("myc_status_pos_vs_neg", "timepoint12W.myc_statuspos")), pAdjustMethod="none")
ihw_resMyc12 <- ihw(
  pvalues    = resMyc12$pvalue,
  covariates = resMyc12$baseMean,
  alpha      = 0.1
)
resMyc12$IHW_padj <- adj_pvalues(ihw_resMyc12)

# > sum(resMyc12$pvalue < 0.1, na.rm = TRUE)
# [1] 3048
# > sum(resMyc12$IHW_padj < 0.1, na.rm = TRUE)
# [1] 239
# > sum(myc_effect_at_12W$padj < 0.1, na.rm = TRUE)
# [1] 218

plot(ihw_resMyc12)
plot(ihw_resMyc12, what = "decisionboundary") 

# IHW is useful in this case!! merging
myc_effect_at_12W.IHW <- myc_effect_at_12W
myc_effect_at_12W.IHW$padj <- resMyc12$IHW_padj

DESeq2::plotMA(myc_effect_at_12W.IHW, ylim=c(-2,2))

##### 5. diff myc effect between 6W and 12W ----

myc_effect_12W_vs_6W <- results(dds, name = "timepoint12W.myc_statuspos")
summary(myc_effect_12W_vs_6W)

###### SUM DGE results from interaction model -----

timepoint_effect_12W_vs_6W_at_myc_neg.ape
timepoint_effect_12W_vs_6W_at_myc_pos
myc_effect_overall.ape.IHW
myc_effect_at_12W.IHW
<<<<<<< HEAD
=======




>>>>>>> e51bf34e25f3ac64d2e1b39c33044662d47df567


##### compare LFCs between 6W and 12W -----


# Merge results from the two time points
merged_res_12W_vs_6W <- merge(as.data.frame(myc_effect_overall.ape.IHW), 
                              as.data.frame(myc_effect_at_12W.IHW), 
                              by = "row.names", 
                              suffixes = c("_6W", "_12W"))

# Rename row names column
colnames(merged_res_12W_vs_6W)[1] <- "Ensembl_ID"

# Create scatter plot for ALL genes
ggplot(merged_res_12W_vs_6W, aes(x = log2FoldChange_6W, y = log2FoldChange_12W)) +
  geom_point(alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "log2FC (Myc Effect at 6W)", y = "log2FC (Myc Effect at 12W)", 
       title = "Comparison of Myc Effects at 6W vs 12W (All Data)") +
  theme_minimal()

# Filter for significant genes (padj < 0.1 at both timepoints)

merged_res_12W_vs_6W <- merged_res_12W_vs_6W %>%
  mutate(Significance = case_when(
    padj_6W < 0.1 & padj_12W < 0.1 ~ "Significant at both",
    padj_6W < 0.1 & padj_12W >= 0.1 ~ "Significant at 6W only",
    padj_6W >= 0.1 & padj_12W < 0.1 ~ "Significant at 12W only",
    TRUE ~ "Not significant"
  ))

# Choose ColorBrewer palette
palette_colors <- brewer.pal(n = 4, name = "Dark2")  # "Dark2" is a good categorical choice

# Map significance categories to ColorBrewer colors
color_palette <- c(
  "Significant at both" = palette_colors[1],  # Dark green
  "Significant at 6W only" = palette_colors[2], # Dark orange
  "Significant at 12W only" = palette_colors[3], # Dark purple
  "Not significant" = "gray70"  # Use gray for non-significant
)

# Define alpha values for different categories
alpha_values <- c(
  "Significant at both" = 0.5,
  "Significant at 6W only" = 0.2,
  "Significant at 12W only" = 0.9,
  "Not significant" = 0.02
)

# Create scatter plot with different colors and alpha values
pdf('output/Myc_effect_12W_vs_6W.pdf', width = 14, height = 10)

ggplot(merged_res_12W_vs_6W, 
            aes(x = log2FoldChange_6W, 
                y = log2FoldChange_12W, 
                color = Significance, 
                alpha = Significance)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  scale_color_manual(values = color_palette) +  # Apply ColorBrewer colors
  scale_alpha_manual(values = alpha_values) +   # Apply different alpha values
  xlim(-5, 5) +   # Limit x-axis between -5 and +5
  ylim(-5, 5) +   # Limit y-axis between -5 and +5
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray50") +  # Dashed lines at y = ±1
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +  # Dashed lines at x = ±1
  labs(x = "log2FC (Myc Effect at 6W)", 
       y = "log2FC (Myc Effect at 12W)", 
       title = "Myc Effects 6W vs 12W") +
  theme_minimal() +
  theme(legend.title = element_blank())

dev.off()

# extract gene sets
# Connect to Ensembl biomart for Mouse (Mus musculus)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get gene names for Ensembl IDs
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"), 
  filters = "ensembl_gene_id", 
  values = merged_res_12W_vs_6W$Ensembl_ID, 
  mart = mart
)

# Merge gene names into dataset
merged_res_12W_vs_6W <- merged_res_12W_vs_6W %>%
  left_join(gene_annotations, by = c("Ensembl_ID" = "ensembl_gene_id"))

# Replace missing gene names with Ensembl IDs
merged_res_12W_vs_6W <- merged_res_12W_vs_6W %>%
  mutate(mgi_symbol = ifelse(mgi_symbol == "", Ensembl_ID, mgi_symbol))

# Extract and sort genes for each category, divide positive and negative
genes_both <- merged_res_12W_vs_6W %>%
  filter(Significance == "Significant at both") %>%
  arrange(desc(log2FoldChange_6W)) %>%
  pull(mgi_symbol)

genes_both_pos <- merged_res_12W_vs_6W %>%
  filter(Significance == "Significant at both") %>%
  filter(log2FoldChange_6W > 0) %>% 
  arrange(desc(log2FoldChange_6W)) %>%
  pull(mgi_symbol)

genes_both_neg <- merged_res_12W_vs_6W %>%
  filter(Significance == "Significant at both") %>%
  filter(log2FoldChange_6W < 0) %>% 
  arrange(log2FoldChange_6W) %>%
  pull(mgi_symbol)

genes_6W_only <- merged_res_12W_vs_6W %>%
  filter(Significance == "Significant at 6W only") %>%
  arrange(desc(log2FoldChange_6W)) %>%
  pull(mgi_symbol)

genes_6W_only_pos <- merged_res_12W_vs_6W %>%
  filter(Significance == "Significant at 6W only") %>%
  filter(log2FoldChange_6W > 0) %>% 
  arrange(desc(log2FoldChange_6W)) %>%
  pull(mgi_symbol)

genes_6W_only_neg <- merged_res_12W_vs_6W %>%
  filter(Significance == "Significant at 6W only") %>%
  filter(log2FoldChange_6W < 0) %>% 
  arrange(log2FoldChange_6W) %>%
  pull(mgi_symbol)

genes_6W_all <- merged_res_12W_vs_6W %>%
  filter(Significance %in% c("Significant at both", "Significant at 6W only")) %>%
  arrange(desc(log2FoldChange_6W)) %>%
  pull(mgi_symbol)

genes_6W_all_pos <- merged_res_12W_vs_6W %>%
  filter(Significance %in% c("Significant at both", "Significant at 6W only")) %>%
  filter(log2FoldChange_6W > 0) %>% 
  arrange(desc(log2FoldChange_6W)) %>%
  pull(mgi_symbol)

genes_6W_all_neg <- merged_res_12W_vs_6W %>%
  filter(Significance %in% c("Significant at both", "Significant at 6W only")) %>%
  filter(log2FoldChange_6W < 0) %>% 
  arrange(log2FoldChange_6W) %>%
  pull(mgi_symbol)

genes_12W_only <- merged_res_12W_vs_6W %>%
  filter(Significance == "Significant at 12W only") %>%
  arrange(desc(log2FoldChange_12W)) %>%
  pull(mgi_symbol)

genes_12W_only_pos <- merged_res_12W_vs_6W %>%
  filter(Significance == "Significant at 12W only") %>%
  filter(log2FoldChange_12W > 0) %>% 
  arrange(desc(log2FoldChange_12W)) %>%
  pull(mgi_symbol)

genes_12W_only_neg <- merged_res_12W_vs_6W %>%
  filter(Significance == "Significant at 12W only") %>%
  filter(log2FoldChange_12W < 0) %>% 
  arrange(log2FoldChange_12W) %>%
  pull(mgi_symbol)




##### gprofiler analysis ----

#  for web ui
write_clip(genes_6W_only)
write_clip(genes_12W_only)
write_clip(genes_both)
write_clip(genes_6W_all)

# Background genes for gprofiler = all genes tested in DESeq2
background_genes <- merged_res_12W_vs_6W$mgi_symbol

# Save as a text file for g:Profiler web ui input
write.table(background_genes, "output/background_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Define a function that runs g:Profiler for a given gene set.
run_gprofiler <- function(query_genes, gene_set_name, background) {
  # Run g:Profiler analysis with the following parameters:
  # - ordered_query = TRUE : the query gene list is ordered.
  # - organism = "mmusculus"
  # - custom background provided
  # - sources: GO:BP, GO:CC, KEGG, REAC, TRANSFAC
  res <- gost(query = query_genes,
              organism = "mmusculus",
              ordered_query = TRUE,
              custom_bg = background,
              evcodes = TRUE,  #also lists the intersection genes in a column
              sources = c("GO:BP", "GO:CC", "KEGG", "REAC", "MIRNA", "TF"))
  
  # If there are results, add a column indicating the gene set name.
  if (!is.null(res$result)) {
    res_df <- res$result %>%
      mutate(gene_set = gene_set_name)
  } else {
    # If no significant results were returned, create an empty data frame with the expected structure.
    res_df <- data.frame(gene_set = gene_set_name)
  }
  
  return(res_df)
}

gene_sets <- list(
  '12W_only'     = genes_12W_only,
  'both_n'     = genes_both_neg,
  '6W_only_n'  = genes_6W_only_neg,
  '6W_all_n'   = genes_6W_all_neg,
  'both_p'     = genes_both_pos,
  '6W_only_p'  = genes_6W_only_pos,
  '6W_all_p'   = genes_6W_all_pos)

# genes_12W_only_neg, genes_12W_only_pos doesn't give any significant result, going with genes_12W_only


# Run g:Profiler for each gene set and store the results in a list
gprofiler_results_list <- lapply(names(gene_sets), function(set_name) {
  message("Running g:Profiler for gene set: ", set_name)
  run_gprofiler(query_genes = gene_sets[[set_name]],
                gene_set_name = set_name,
                background = background_genes)
})

# Combine all the results into one data frame
all_pathways <- bind_rows(gprofiler_results_list)

# Reorder columns so that 'gene_set' is the first column
all_pathways <- all_pathways %>% select(gene_set, everything())

# Save the final merged results to a CSV file for further analysis/visualization.
all_pathways_fixed <- all_pathways %>% 
  mutate(across(where(is.list), ~ sapply(., function(x) paste(unlist(x), collapse = ","))))

write.csv(all_pathways_fixed, "output/gprofiler/merged_gprofiler_results.csv", row.names = FALSE)

# Adjust p-values to -log10 scale for visualization
all_pathways$neg_log10_pval <- -log10(all_pathways$p_value)

# Create a summary of pathway occurrence
pathway_counts <- all_pathways %>%
  group_by(term_id, term_name, source) %>%
  summarise(occurrences = n(), .groups = "drop")

# Create a binary presence/absence matrix for UpSet plot
pathway_presence <- all_pathways %>%
  select(term_id, gene_set) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = gene_set, values_from = present, values_fill = 0)

###### UpSet plot for gprofiler results ----

# Convert pathway presence matrix into a format suitable for UpSet
upset_matrix <- pathway_presence %>%
  column_to_rownames("term_id") %>%
  as.matrix()

upset_matrix <- as.data.frame(upset_matrix)

colorsUpSet <- brewer.pal(n = 7, name = "Dark2")

# Create UpSet plot
pdf("output/gprofiler/gprofiler_upset_1.pdf")

upset(
  upset_matrix, 
  sets = colnames(upset_matrix),
  sets.bar.color = colorsUpSet,
  keep.order = TRUE,
  mainbar.y.label = "Number of Shared Pathways",
  sets.x.label = "Number of Pathways in Each Set"
)

dev.off()

# extract pathway names for each category

pathway_presence_merged <- left_join(pathway_presence, 
                                     all_pathways %>% select(term_id, term_name),
                                     by = "term_id")

# Specify the gene set columns
gene_set_columns <- c("both_n", "both_p", 
                      "6W_only_n", "6W_only_p", 
                      "6W_all_n", "6W_all_p", 
                      "12W_only")

# Create an "intersection" identifier for each pathway,
# indicating the gene set columns (from gene_set_columns) where the pathway is present (value == 1).
pathway_presence_merged <- pathway_presence_merged %>%
  rowwise() %>%
  mutate(intersection = {
    # Get the gene set names for which the pathway is present (== 1)
    present_sets <- gene_set_columns[which(c_across(all_of(gene_set_columns)) == 1)]
    if(length(present_sets) == 0) "none" else paste(sort(present_sets), collapse = ";")
  }) %>%
  ungroup()

# Now split the pathway names (term_name) by the unique intersection categories.
# This gives you a list where each element is the vector of pathway names for a given intersection category.
intersection_list <- split(pathway_presence_merged$term_name, pathway_presence_merged$intersection)
intersection_list <- lapply(intersection_list, unique)

# write each path list in a separate sheet of and excel file

wb <- createWorkbook()
for (int_name in names(intersection_list)) {
  # Create a valid sheet name (max 31 characters; truncate if necessary)
  sheet_name <- substr(int_name, 1, 31)
  
  # Create a data frame with the pathway names.
  df <- data.frame(Pathway = intersection_list[[int_name]], stringsAsFactors = FALSE)
  
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, x = df)
}

# Save the Excel workbook.
saveWorkbook(wb, file = "output/gprofiler/upset_pathways.xlsx", overwrite = TRUE)


##### Cytoscape -----















##### fgsea -----

# Create ranked gene list for GSEA
ranks <- merged_res_12W_vs_6W %>%
  arrange(desc(log2FoldChange_6W)) %>%  # Rank by LFC at 6W (or use 12W if needed)
  select(mgi_symbol, log2FoldChange_6W) %>%
  distinct()  # Ensure unique gene names

# Convert to named numeric vector
gene_ranks <- setNames(ranks$log2FoldChange_6W, ranks$mgi_symbol)

# Load gene sets (e.g., MSigDB, KEGG, Reactome)
pathways <- gmtPathways("mouse_msigdb.gmt")  # Replace with the correct gene set file

# Run fgsea
fgsea_results <- fgsea(pathways = pathways, 
                       stats = gene_ranks, 
                       minSize = 15, 
                       maxSize = 500, 
                       nperm = 10000)  # Increase permutations for stability

# Filter significant pathways
fgsea_results <- fgsea_results %>%
  arrange(padj)  # Sort by adjusted p-value

# Save results
write.table(fgsea_results, "GSEA_results_6W.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot top enriched pathways
plotEnrichment(pathways[["KEGG_OXIDATIVE_PHOSPHORYLATION"]], gene_ranks)  # Example pathway







