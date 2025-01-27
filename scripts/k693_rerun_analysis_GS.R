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



###### Read in and set up coldata

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


###### Read in and set up count file

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

##### create a DESeq object

ddsMat <- DESeqDataSetFromMatrix(countData = cts,
                                 colData = coldata,
                                 design = ~ timepoint * myc_status) #this expands to ~ timepoint + myc_status + timepoint:myc_status

##### look at the ddsMat data

## filter by counts: 55K rows --> 19K

nrow(ddsMat)

smallestGroupSize <- 4
keep <- rowSums(counts(ddsMat) >= 10) >= smallestGroupSize
dds <- ddsMat[keep,]

nrow(dds)

## normalise and transform to look at variance

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

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

# MYCF64_4g / 6W_pos is an outlier (also MYCF64_3f / 6W_neg) ???
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

###### DGE results / contrasts

dds <- DESeq(dds)

resultsNames(dds)

# main effect of timepoints (myc_neg)

timepoint_effect_12W_vs_6W_at_myc_neg <- results(dds, name = "timepoint_12W_vs_6W")

# main effect of timepoints (myc_pos)

timepoint_effect_12W_vs_6W_at_myc_pos <- results(dds, contrast = list(c("timepoint_12W_vs_6W", "timepoint12W.myc_statuspos")))

# overall myc effect

myc_effect_overall <- results(dds, name = "myc_status_pos_vs_neg")

# myc effect at 6W

myc_effect_at_6W <- results(dds, contrast = c("myc_status", "pos", "neg"))

# myc effect at 12W

myc_effect_at_12W <- results(dds, contrast = list(c("myc_status_pos_vs_neg", "timepoint12W.myc_statuspos")))

# change in myc effect between 6W and 12W

myc_effect_12W_vs_6W <- results(dds, name = "timepoint12W.myc_statuspos")





