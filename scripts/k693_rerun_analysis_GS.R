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

###### Read in and set up coldata

coldata <- read.delim("data/FULL.DAT.COL.DATA.txt", sep = "\t", header = TRUE)
coldata <- coldata %>% remove_rownames %>% column_to_rownames(var="count.file") %>% 
  mutate(myc_status = case_when(
    group == "hyperplastic" ~ "pos",
    group == "WT" ~ "neg"
  )) %>% 
  mutate(timepoint = case_when(
    age == '6 wk' ~ '6W',
    age == '12 wk' ~'12W'
  )) %>% 
  mutate(group = paste(timepoint, myc_status, sep = '_'))
coldata <- coldata[-c(3,4)]
coldata$timepoint <- as.factor(coldata$timepoint) %>% relevel('6W')
coldata$myc_status <- as.factor(coldata$myc_status)

###### Read in and set up count file

cts <- read.csv("data/FULL.DAT.csv")
cts <- cts[-2]
cts <- cts %>% remove_rownames %>% column_to_rownames(var="gene_id")
cts <- as.matrix(round(cts))
storage.mode(cts) <- "integer"

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

vsd <- vst(dds, blind = TRUE)
head(assay(vsd), 3)

rld <- rlog(dds, blind = TRUE)
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

## PCA

plotPCA(vsd, intgroup = c("group", "sample"))



