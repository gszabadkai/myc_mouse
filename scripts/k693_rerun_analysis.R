rm(list = ls())

setwd("/Users/michaelsworklaptop/Documents/k693_rerun")
# Get the working directory 
getwd()
#BiocManager::install("apeglm")
#BiocManager::install("AnnotationHub")
#packageVersion("msigdbr")
library(stringr)
library(DESeq2)
require(apeglm)
library(dplyr)
library(tibble)
library(ggplot2)
library(enrichplot)
library(pheatmap)
library(ggvenn)
library(RColorBrewer)
library(ggrepel)
library(msigdbr)
library(AnnotationHub)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(cowplot)
library(fgsea)
library(data.table)
library(EnhancedVolcano)
#install.packages("msigdbr")

#load msigd database
all_gene_sets <- msigdbr(species = "Mus musculus")
h_gene_sets <- msigdbr (species = 'mouse' , category = 'H')
msigdbr_collections()

msigdb_mouse_hallmarks <- msigdbr(species = "Mus musculus")

################################################################################
###### Investigating the genome annotation here. no need to run this part ######
###### everytime you run the code as, this is to check the annotation ##########
################################################################################
#columns(org.Mm.eg.db)  # Return the column names here
#keytypes(org.Mm.eg.db) #Filter the database by key, this returns a ENSEMBL KEYS
#keys(org.Mm.eg.db, keytype="ENSEMBL")[1:10]
#select(org.Mm.eg.db, 
#keys="ENSMUSG00000000381",
#keytype = "ENSEMBL",
#columns=c("SYMBOL","GENENAME"))

################################################################################
###### Read the count matrix data in and the Samples information  ##############
################################################################################
################################################################################
#cts <- floor(as.matrix(read.delim("FULL.DAT.txt", sep = "\t", header = TRUE)))
coldata <- read.delim("FULL.DAT.COL.DATA.txt", sep = "\t", header = TRUE)
cts = read.csv("FULL.DAT.csv")

###Remove transcript ID###
cts <- cts[-2]
cts = cts %>% remove_rownames %>% column_to_rownames(var="gene_id")
coldata = coldata %>% remove_rownames %>% column_to_rownames(var="sample")


# Check whether the rows in the sample file are equal to the columns in the count matrix
all(rownames(coldata) %in% colnames(cts)) # This should return true, otherwise needs changing

# Check the matches here before loading the data for further analysis
all(rownames(coldata) %in% colnames(cts)) # This should return TRUE

# Set the conditions as factors in the coldata
coldata$cancer <- as.factor(coldata$age)

# Load the dataset and create a DESeq object for group
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age)
head(dds)

# Pre-filtering => Removes all genes whose row sum is less than 10. 
# This can be adjusted for more inclusivity
# Quality control analysis 
# Normalisation
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# Explore the dataset here
normlzd_dds <- counts(dds, normalized=T)
head(normlzd_dds)
write.csv(normlzd_dds, "Normalised_DESseq_counts.csv")

# Plot the column sums according to size factor
#plot(sizeFactors(dds), colSums(counts(dds)))
#abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))


#######################Visualize data for quality control##########
# Using hierarchical clustering by protocal
#plot(hclust(dist(t(normlzd_dds))), labels=colData(dds)$group)

# Cluster by stage
#plot(hclust(dist(t(normlzd_dds))), labels=colData(dds)$stage)

# Assess the poison noise for low reads
#plot(log(normlzd_dds[,1])+1, log(normlzd_dds[,2])+1, cex =.1)

# Varaiance Stabilizing transformation
#vsd <- vst(dds, blind = T)

# extract the vst matris from the object
#vsd_mat <- assay(vsd)

# compute pairwise correlation values
#vsd_cor <- cor(vsd_mat)

#vsd_cor

# Visualise this using a heatmap
#pheatmap(vsd_cor)

# Perform PCA based on groups and stages here

#plotPCA(vsd, intgroup = "group")
#plotPCA(vsd, intgroup = "cancer")
#plotCounts(dds, gene=which.min(res$padj), intgroup="stage")

## Mean-variance relationship of the first 3 samples
# Calculating mean for each gene
#mean_readCounts <- apply(cts[,1:3], 1, mean)

# Calculating variance for each gene
#var_readCounts <- apply(cts[,1:3], 1, var)

# Visualise the data above here
df <- data.frame(mean_readCounts, var_readCounts)
ggplot(df) +
  geom_point(aes(x=mean_readCounts, y= var_readCounts)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Mean counts per gene") +
  ylab("Variance per gene") +
  labs(title = "DESeq2 model - Dispersion")


##########################################################################
################# Filter some genes out, counts below 10 ##################
###########################################################################
keep <- rowSums(counts(dds)) >= 10    
dds <- dds[keep,]
head(dds)
write.csv(normlzd_dds, "Normalised_DESseq_counts.csv")

###########################################################################
# Differential gene expression on the whole dataframe
dds <- DESeq(dds)

#for multiple compairson eg groups

#POS_v_Neg <- results(dds, contrast=c("group","hyperplastic","WT"))
#twelve_Neg_v_six_Neg <- results(dds, contrast=c("cancer","12 Wk MYC neg","6 Wk MYC neg"))

#### Perform the comprison of all 12 week v 6 week regardless of Myc status
twelve_all_v_six_all <- results(dds, contrast=c("age","12 wk","6 wk"))
head(twelve_all_v_six_all )

res <- twelve_all_v_six_all 
res<-na.omit(res)

head(res)
type(res)
res = as.data.frame(res)
write.csv(res, "twelve_all_v_six_all_deseq2.csv")


##### Plot a volcano plot of DGs #########
#res = read.csv("POS_all_v_NEG_all_deseq2.csv")
EnhancedVolcano(res,
                lab = res,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                title = "All Myc +ive v All Myc -ive",
                )


# resSort is the res file (DESeq results) sorted by log2foldchange
resSort<- res[order(res$log2FoldChange, decreasing = TRUE) ,]
head(resSort)

#Processing for Over Rrepresentation Analysis
#Make a genelist, set the parameters of baseMean>50 and padj < 0.05
gene_list_1 <- resSort[resSort$baseMean>50 ,]
gene_list_2 <- gene_list_1[gene_list_1$padj<0.05,]

#Over represented
gene_list_sig_1 <- gene_list_2[gene_list_2$log2FoldChange >0.5,]
#Underrepresented
gene_list_sig_2<- gene_list_2[gene_list_2$log2FoldChange < 0.5,]
gene_list_sig_2


# 1. make a gene list for ORA
gene_list_sig_for_ORA_pos = rownames(gene_list_sig_1)
head(gene_list_1)

write.csv(gene_list_sig_for_ORA_pos, "POS_all_v_NEG_all_Sig_UP_genelist.csv")

gene_list_sig_for_ORA_neg = rownames(gene_list_sig_2)
gene_list_sig_for_ORA_neg
write.csv(gene_list_sig_for_ORA_neg, "POS_all_v_NEG_all_Sig_down_genelist.csv")

#Perform ORA large v 6wk neg
pos_enriched <- enrichGO(gene     = gene_list_sig_for_ORA_pos,
                         OrgDb         = "org.Mm.eg.db",
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable = TRUE)
as.data.frame(pos_enriched)
write.csv(pos_enriched, "POS_all_v_NEG_all_Sig_UP_ORA.csv")
pos_enriched_plot <- plot(barplot(pos_enriched, showCategory = 15))


neg_enriched <- enrichGO(gene     = gene_list_sig_for_ORA_neg,
                         OrgDb         = "org.Mm.eg.db",
                         keyType       = 'ENSEMBL',
                         ont           = "all",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable = TRUE)
yy = as.data.frame(neg_enriched)
fit = plot(barplot(neg_enriched, showCategory = 15))

write.csv(neg_enriched, "POS_all_v_NEG_all_Sig_down_ORA_genelist.csv")


######## Annotate our data here #########
# Res is the original DEseq output
# Check the row names here and change them to one easily accessible here
res <- as.data.frame(resSort) %>% rownames_to_column("ENSEMBL")
res


## Actually complete the mapping with symbol and Gene name and EntrezID
res_anno <- AnnotationDbi::select(org.Mm.eg.db,keys=res$ENSEMBL,
                                  columns=c("ENSEMBL","SYMBOL","GENENAME","ENTREZID"),
                                  keytype="ENSEMBL") %>% 
  filter(!duplicated(ENSEMBL))

# Have a look at the annotation
head(res_anno)

# Bind the annotation to the original dataset here
res.annotated <- left_join(res_anno, res, by= "ENSEMBL")
res.annotated<-na.omit(res.annotated)
head(res.annotated)


#now we can order the data frame

res.annotated.sort <- res.annotated[order(-res.annotated$log2FoldChange) ,]
res.annotated.sort
write.csv(res.annotated.sort, "POS_all_v_NEG_all_annotated_genelist.csv")


#make a df with log2fold change and entrezid 
res_ENTREZ = res.annotated.sort$log2FoldChange
names(res_ENTREZ) = res.annotated.sort$ENTREZID
head(res_ENTREZ)

#make a df with log2fold change and EnsembleID
res_ENSEMBL = res.annotated.sort$log2FoldChange
names(res_ENSEMBL) = res.annotated.sort$ENSEMBL

#make a df with log2fold change and symbol
res_Symbol = res.annotated.sort$log2FoldChange
names(res_Symbol) = res.annotated.sort$SYMBOL

##### GSEA uses the full DG list
#perform GSEA using Gene Ontology terms with the gse package
## super slow so use the fgsea command below###
#require(DOSE)
#gse_12_WK_NEG_v_neg <- gseGO(geneList=res_ENSEMBL, 
                         #ont ="BP", 
                         #keyType = "ENSEMBL", 
                         #minGSSize = 10, 
                         #maxGSSize = 900, 
                         #pvalueCutoff = 0.05, 
                         #verbose = TRUE, 
                         #eps = 1e-100,
                         #OrgDb = "org.Mm.eg.db", 
                         #pAdjustMethod = "BH")


#gseaplot(gse_12_WK_NEG_v_neg, geneSetID = 24)
#dotplot(gse_12_WK_NEG_v_neg, showCategory=10, split=".sign") + facet_grid(.~.sign)
#as.data.frame(gse_Large_v_neg)

#GSEA using KEGG
#gse_12_WK_NEG_v_6wk_neg <- gseKEGG(geneList     = res_ENTREZ,
#                           organism     = 'mmu',
#                           keyType = "kegg", 
#                           minGSSize    = 15,
#                           maxGSSize  = 800, 
#                           pvalueCutoff = 0.05,
#                           verbose      = TRUE,
#                          pAdjustMethod = "BH",
#                         eps = 1e-100)

#as.data.frame(large_v_6wk_neg)
#dotplot(gse_12_WK_NEG_v_6wk_neg, showCategory=20, split=".sign") + facet_grid(.~.sign)


#GSEA using msigdb and fgsea package
all_gene_sets = msigdbr(species = "Mus musculus")
head(all_gene_sets)

#msigdbr_species()
#You can retrieve data for a specific collection, such as the hallmark gene sets.

h_gene_sets = msigdbr(species = "mouse", category = "H")
#head(h_gene_sets)

cgp_gene_sets = msigdbr(species = "mouse", category = "C2")
#(cgp_gene_sets)

c5_gene_sets = msigdbr(species = "mouse", category = "C5")
#(c5_gene_sets)

c6_gene_sets = msigdbr(species = "mouse", category = "C6")
#(c6_gene_sets)

#There is a helper function to show the available collections.

#msigdbr_collections()
#The msigdbr() function output is a data frame and can be manipulated using more standard methods.

#all_gene_sets %>%
# dplyr::filter(gs_cat == "H") %>%
#head()

#this command we are manipulating the dataframe into the different gene sets as lists
# as EntrezID, use code head(msigdbr_h_list) to see the output

msigdbr_h_list = split(x = h_gene_sets$entrez_gene, f = h_gene_sets$gs_name)

msigdbr_cpg_list = split(x = cgp_gene_sets$entrez_gene, f = cgp_gene_sets$gs_name)

msigdbr_c5_list = split(x = c5_gene_sets$entrez_gene, f = c5_gene_sets$gs_name)

msigdbr_c6_list = split(x = c6_gene_sets$entrez_gene, f = c6_gene_sets$gs_name)


res_h_fgsea <- fgseaMultilevel(pathways = msigdbr_h_list,
                               stats = res_ENTREZ,
                               minSize = 15,
                               maxSize = 900,
)

#Results have all the Hallmarks, filter on p value
res_h_fgsea
res_h_fgsea_filt<- res_h_fgsea[res_h_fgsea$padj<0.05,]

fwrite(res_h_fgsea, "POS_all_v_NEG_all_MSigDB_Hallmarks_full.csv")
fwrite(res_h_fgsea_filt, "POS_all_v_NEG_all_MSigDB_Hallmarks_significant.csv") 

plotEnrichment(msigdbr_h_list[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]],
               res_ENTREZ) + labs(title="12 Wk v 6 Wk - HALLMARK_OXIDATIVE_PHOSPHORYLATION")

res_cpg_fgsea <- fgseaMultilevel(pathways = msigdbr_cpg_list,
                                 stats = res_ENTREZ,
                                 minSize = 15,
                                 maxSize = 800,
)

res_cpg_fgsea_filt<- res_cpg_fgsea[res_cpg_fgsea$padj<0.05,]

fwrite(res_cpg_fgsea, file = "POS_all_v_NEG_all_MSigDB_CPG_full.csv") 
fwrite(res_cpg_fgsea_filt, file = "POS_all_v_NEG_all_MSigDB_CPG_significant.csv") 


res_c6_fgsea <- fgseaMultilevel(pathways = msigdbr_c6_list,
                                stats = res_ENTREZ,
                                minSize = 10,
                                maxSize = 800,
)

res_c6_fgsea_filt<- res_c6_fgsea[res_c6_fgsea$padj<0.05,]

fwrite(res_c6_fgsea, file = "POS_all_v_NEG_all_MSigDB_C6_full.csv") 
fwrite(res_c6_fgsea_filt, file = "POS_all_v_NEG_all_MSigDB_C6_significant.csv") 

res_c5_fgsea <- fgseaMultilevel(pathways = msigdbr_c5_list,
                                stats = res_ENTREZ,
                                minSize = 10,
                                maxSize = 800,
)

res_c5_fgsea_filt<- res_c5_fgsea[res_c5_fgsea$padj<0.05,]
fwrite(res_c5_fgsea, file = "POS_all_v_NEG_all_MSigDB_C5_full.csv") 
fwrite(res_c5_fgsea_filt, file = "POS_all_v_NEG_all_MSigDB_C5_significant.csv") 

mito_pathways_ENTREZ <- gmtPathways("Mouse.MitoCarta3.0.GMT.EntrezIDs.gmt")
mito_pathways_ENTREZ

res_mitoCarta_fgsea <- fgseaMultilevel(pathways = mito_pathways_ENTREZ,
                                       stats = res_ENTREZ,
                                       minSize = 5,
                                       maxSize = 800,
                                       eps = 0,
                                       
                                       
)

res_mitoCarta_fgsea_filt<- res_mitoCarta_fgsea[res_mitoCarta_fgsea$padj<0.05,]

fwrite(res_mitoCarta_fgsea, file = "POS_v_NEG_MSigDB_MitoCarta_full.csv") 
fwrite(res_mitoCarta_fgsea_filt, file = "POS_v_NEG_MSigDB_MitoCarta_significant.csv") 

####################### Importing results so dont re run above

res_h_fgsea_filt = read.csv("POS_all_v_NEG_all_MSigDB_Hallmarks_significant.csv")
res_CPG_fgsea_filt = read.csv("POS_all_v_NEG_all_MSigDB_CPG_significant.csv")
res_CPG_fgsea_filt_reactome = read.csv("POS_all_v_NEG_all_MSigDB_CPG_significant_REACTOME.csv")
res_CPG_fgsea_filt_kegg = read.csv("POS_all_v_NEG_all_MSigDB_CPG_significant_KEGG.csv")
res_C5_fgsea_filt = read.csv("POS_all_v_NEG_all_MSigDB_C5_significant.csv")
res_C5_fgsea_filt_bp = read.csv("POS_all_v_NEG_all_MSigDB_C5_significant_BP.csv")
res_mitoCarta_fgsea_filt = read.csv("POS_v_NEG_MSigDB_MitoCarta_significant.csv")
res_C6_fgsea_filt = read.csv("POS_all_v_NEG_all_MSigDB_C6_significant.csv")

############### To plot halllmarks ##################

ggplot(res_h_fgsea_filt, aes(x = NES, y = reorder(pathway, NES), color = padj, size = size)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Normalized Enrichment Score", y = "Pathway",
       color = "Adjusted p-value", size = "Gene Count") +
  ggtitle("12 Wk Myc Pos V Myc Neg (Isolated epithelial cells) \nMSigDB Hallmark pathways") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank())
ggsave("Pos_V_Neg_Hallmarks_Dotplot.tiff", width = 8, height = 6.5, dpi = 800)

############### To plot MitoCarta ##################
top20_enriched <- head(res_mitoCarta_fgsea_filt[order(-res_mitoCarta_fgsea_filt$NES),], 30)

ggplot(top20_enriched, aes(x = NES, y = reorder(pathway, NES), color = padj, size = size)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Normalized Enrichment Score", y = "Pathway",
       color = "Adjusted p-value", size = "Gene Count") +
  ggtitle("12 Wk Myc Pos V Myc Neg (Isolated epithelial cells) \nTop 20 MitoCarta pathways") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank())
ggsave("Pos_V_Neg_mitocarta_Dotplot.tiff", width = 8, height = 6.5, dpi = 800)

# Filter the data frame to include only the top 10 enriched and top 10 downregulated gene sets
top10_enriched <- head(res_CPG_fgsea_filt_reactome[order(-res_CPG_fgsea_filt_reactome$NES),], 10)
top10_downregulated <- head(res_CPG_fgsea_filt_reactome[order(res_CPG_fgsea_filt_reactome$NES),], 10)
res_top10 <- rbind(top10_enriched, top10_downregulated)

# To plot reactome bar plot
ggplot(res_top10, aes(x = NES, y = reorder(pathway, NES), color = padj, size = size)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Normalized Enrichment Score", y = "Pathway",
       color = "Adjusted p-value", size = "Gene Count") +
  ggtitle("Top 10 enriched and Top 10 \n downregulated - Reactome") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank())
ggsave("Cancer_V_WT_Reactome_dot_plot.tiff", width = 9.5, height = 6, dpi = 800)


ggplot(res_top10, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = as.factor(sign(NES)))) +
  scale_fill_manual(values = c("red", "blue"), labels = c("Underrepresented", "Overrepresented")) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       fill = "NES Sign") + 
  ggtitle("Top 10 enriched and top 10 downregulated - Reactome  ") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank())
ggsave("Cancer_V_WT_Reactome bar pot.tiff", width = 4, height = 6, dpi = 800)

####Plot kegg pathways

ggplot(res_CPG_fgsea_filt_kegg, aes(x = NES, y = reorder(pathway, NES), color = padj, size = size)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Normalized Enrichment Score", y = "Pathway",
       color = "Adjusted p-value", size = "Gene Count") +
  ggtitle("KEGG Terms") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank())
ggsave("Cancer_V_WT_KEGG_dot_plot.tiff", width = 7, height = 6, dpi = 800)

ggplot(res_CPG_fgsea_filt_kegg, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = as.factor(sign(NES)))) +
  scale_fill_manual(values = c("red", "blue"), labels = c("Underrepresented", "Overrepresented")) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       fill = "NES Sign") + 
  ggtitle("KEGG pathways") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank())
ggsave("Cancer_V_WT_KEGG_bar_pot.tiff", width = 12, height = 6, dpi = 800)

#####GOBP plott

top10_enriched <- head(res_C5_fgsea_filt_bp[order(-res_C5_fgsea_filt_bp$NES),], 10)
top10_downregulated <- head(res_C5_fgsea_filt_bp[order(res_C5_fgsea_filt_bp$NES),], 10)
res_BP_top10 <- rbind(top10_enriched, top10_downregulated)


ggplot(res_BP_top10, aes(x = NES, y = reorder(pathway, NES), color = padj, size = size)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Normalized Enrichment Score", y = "Pathway",
       color = "Adjusted p-value", size = "Gene Count") +
  ggtitle("Top 10 enriched and Top 10 \n downregulated - Gene Ontology") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank())
ggsave("Cancer_V_WT_GOBP_top10_dotplot.tiff", width = 10, height = 6, dpi = 800)

