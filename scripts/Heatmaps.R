rm(list=ls())

library(gplots)
library(DESeq2)
require(apeglm)
library(dplyr)
library(tibble)
library(ggplot2)
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

setwd("/Users/michaelsworklaptop/Documents/k693_rerun")

cts = read.csv("FULL.DAT_copy.csv")
colnames(cts) <- sub("X", "", colnames(cts))

###Remove transcript ID###
res <- cts[-2]

res=column_to_rownames(res,var='gene_id')

res = rownames_to_column(res,"ENSEMBL")

res_anno <- AnnotationDbi::select(org.Mm.eg.db,keys=res$ENSEMBL,
                                  columns=c("ENSEMBL","SYMBOL","GENENAME","ENTREZID"),
                                  keytype="ENSEMBL") %>% 
  filter(!duplicated(ENSEMBL))

# Have a look at the annotation
res_anno

# Bind the annotation to the original dataset here
res.annotated <- left_join(res_anno, res, by="ENSEMBL")
res.annotated<-na.omit(res.annotated)
head(res.annotated)


df_2 = res.annotated
head(df_2)

columns_to_remove <- c(1, 3, 4)

df_2 <- res.annotated[, -columns_to_remove]
head(df_2)

write.csv(df_2,"df_2.csv")

# Remove duplicated rows based on the SYMBOL column
df_2 <- df_2 %>% distinct(SYMBOL, .keep_all = TRUE)
# Set column 1 as row names
rownames(df_2) <- df_2[, 1]
head(df_2)
# Remove column 1 from the dataframe
df_2 <- df_2[, -1]

# Remove the specified columns and save as gene_expression
gene_expression <- df_2

#####################################################################
#####################################################################
#########################Glycolysis##################################
#####################################################################
#####################################################################
# Define the list of genes
gene_list <- c("Hk1", "Hk2", "Hk3", "Gpi1", "Pfkm", "Pfkl", "Pfkp", "Aldoa", "Aldob", "Aldoc", "Tpi1", "Gapdh", "Pgk1", "Pgk2", "Pgam1", "Pgam2", "Pgam5", "Eno1", "Eno1b", "Eno2", "Eno3", "Eno4", "Pklr", "Ldha", "Ldhb", "Ldhc", "Ldhd")

# Extract the rows for the genes in the list
selected_genes <- gene_expression[rownames(gene_expression) %in% gene_list,]

# Transform the data to matrix, as required by heatmap.2
selected_genes_matrix <- as.matrix(selected_genes)

# Normalize the data by gene
#selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) x / max(x)))
#t(selected_genes_matrix_norm)
##z-scpre normlaising
selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) scale(x)))
selected_genes_matrix_norm <- na.omit(selected_genes_matrix_norm)
colnames(selected_genes_matrix_norm) <- colnames(selected_genes_matrix)

# Calculate the total expression (e.g., sum) for each sample
total_expression <- colSums(selected_genes_matrix_norm, na.rm = TRUE)

# Order the samples based on total expression (descending order)
ordered_samples <- names(sort(total_expression, decreasing = TRUE))

# Reorder the matrix columns based on expression levels
selected_genes_matrix_norm <- selected_genes_matrix_norm[, ordered_samples]

dev.off()  # Close the current device, if any
dev.new()  # Open a new graphics device

heatmap.2(selected_genes_matrix_norm, 
          main = "Glycolysis", 
          trace = "none", 
          margins = c(5, 10), 
          cexRow = 0.5,
          cexCol = 0.5, 
          col = colorRampPalette(c("blue", "white", "red"))(256),
          key.xlab = "Expression",
          key.title = NA,
          Colv = FALSE,
          dendrogram = "row" # Don't cluster columns
)



#####################################################################
#####################################################################
######################OXPHOS subunits###############################
#####################################################################
#####################################################################

gene_list <- c("Atp5a1", "Atp5b", "Atp5c1", "Atp5d", "Atp5e", "Atp5g1", "Atp5g2", "Atp5g3", "Atp5h", "Atp5j", "Atp5j2", "Atp5k", "Atp5l", "Atp5md", "Atp5mpl", "Atp5o", "Atp5pb", "Atpif1", "Cox4i1", "Cox4i2", "Cox5a", "Cox5b", "Cox6a1", "Cox6a2", "Cox6b1", "Cox6b2", "Cox6c", "Cox7a1", "Cox7a2", "Cox7a2l", "Cox7b", "Cox7c", "Cox8a", "Cox8b", "Cox8c", "Cyc1", "Cycs", "Cyct", "Dmac2l", "Hccs", "mt-Atp6", "mt-Atp8", "mt-Co1", "mt-Co2", "mt-Co3", "mt-Cytb", "mt-Nd1", "mt-Nd2", "mt-Nd3", "mt-Nd4", "mt-Nd4l", "mt-Nd5", "mt-Nd6", "Ndufa1", "Ndufa10", "Ndufa11", "Ndufa12", "Ndufa13", "Ndufa2", "Ndufa3", "Ndufa4", "Ndufa5", "Ndufa6", "Ndufa7", "Ndufa8", "Ndufa9", "Ndufab1", "Ndufb10", "Ndufb11", "Ndufb2", "Ndufb3", "Ndufb4", "Ndufb5", "Ndufb6", "Ndufb7", "Ndufb8", "Ndufb9", "Ndufc1", "Ndufc2", "Ndufs1", "Ndufs2", "Ndufs3", "Ndufs4", "Ndufs5", "Ndufs6", "Ndufs7", "Ndufs8", "Ndufv1", "Ndufv2", "Ndufv3", "Sdha", "Sdhb", "Sdhc", "Sdhd", "Uqcr10", "Uqcr11", "Uqcrb", "Uqcrc1", "Uqcrc2", "Uqcrfs1", "Uqcrh", "Uqcrq")

# Extract the rows for the genes in the list
selected_genes <- gene_expression[rownames(gene_expression) %in% gene_list,]

# Transform the data to matrix, as required by heatmap.2
selected_genes_matrix <- as.matrix(selected_genes)

# Normalize the data by gene
#selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) x / max(x)))
#t(selected_genes_matrix_norm)
selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) scale(x)))
colnames(selected_genes_matrix_norm) <- colnames(selected_genes_matrix)

na_rows <- which(rowSums(is.na(selected_genes_matrix_norm)) > 0)
selected_genes_matrix_norm <- selected_genes_matrix_norm[-na_rows, ]

# Calculate the total expression (e.g., sum) for each sample
total_expression <- colSums(selected_genes_matrix_norm, na.rm = TRUE)

# Order the samples based on total expression (descending order)
ordered_samples <- names(sort(total_expression, decreasing = TRUE))

# Reorder the matrix columns based on expression levels
selected_genes_matrix_norm <- selected_genes_matrix_norm[, ordered_samples]

dev.off()  # Close the current device, if any
dev.new()  # Open a new graphics device

heatmap.2(selected_genes_matrix_norm, 
          main = "OXPHOS Subunits", 
          trace = "none", 
          margins = c(5, 10), 
          cexRow = 0.5,
          cexCol = 0.5, 
          col = colorRampPalette(c("blue", "white", "red"))(256),
          key.xlab = "Expression",
          key.title = NA,
          Colv = FALSE,
          dendrogram = "row" # Don't cluster columns
)




# Add a custom key title above the key
grid.text("Expression", x = unit(0.95, "npc"), y = unit(0.1, "npc"), 
          gp = gpar(fontsize = 12, fontface = "bold"))


#####################################################################
#####################################################################
########################Complex I####################################
#####################################################################
#####################################################################

# Define the list of genes
gene_list <- c("Acad9", "Aifm1", "Dmac1", "Dmac2", "Ecsit", "Foxred1", "Lyrm2", "mt-Nd1", "mt-Nd2", "mt-Nd3", 
               "mt-Nd4", "mt-Nd4l", "mt-Nd5", "mt-Nd6", "Ndufa1", "Ndufa10", "Ndufa11", "Ndufa12", "Ndufa13", "Ndufa2", 
               "Ndufa3", "Ndufa5", "Ndufa6", "Ndufa7", "Ndufa8", "Ndufa9", "Ndufab1", "Ndufaf1", "Ndufaf2", "Ndufaf3", 
               "Ndufaf4", "Ndufaf5", "Ndufaf6", "Ndufaf7", "Ndufaf8", "Ndufb10", "Ndufb11", "Ndufb2", "Ndufb3", "Ndufb4", 
               "Ndufb5", "Ndufb6", "Ndufb7", "Ndufb8", "Ndufb9", "Ndufc1", "Ndufc2", "Ndufs1", "Ndufs2", "Ndufs3", "Ndufs4", 
               "Ndufs5", "Ndufs6", "Ndufs7", "Ndufs8", "Ndufv1", "Ndufv2", "Ndufv3", "Nubpl", "Timmdc1", "Tmem126a", 
               "Tmem126b", "Tmem186", "Tmem70")


# Extract the rows for the genes in the list
selected_genes <- gene_expression[rownames(gene_expression) %in% gene_list,]

# Transform the data to matrix, as required by heatmap.2
selected_genes_matrix <- as.matrix(selected_genes)

# Normalize the data by gene
#selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) x / max(x)))
#t(selected_genes_matrix_norm)
selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) scale(x)))
colnames(selected_genes_matrix_norm) <- colnames(selected_genes_matrix)

# Calculate the total expression (e.g., sum) for each sample
total_expression <- colSums(selected_genes_matrix_norm, na.rm = TRUE)

# Order the samples based on total expression (descending order)
ordered_samples <- names(sort(total_expression, decreasing = TRUE))

# Reorder the matrix columns based on expression levels
selected_genes_matrix_norm <- selected_genes_matrix_norm[, ordered_samples]

dev.off()  # Close the current device, if any
dev.new()  # Open a new graphics device
#tiff("heatmap_complex_I.tiff", width = 1200, height = 1200, units = "px", res = 300)

heatmap_plot = heatmap.2(selected_genes_matrix_norm, 
                         main = "Complex I", 
                         trace = "none", 
                         margins = c(5, 7), 
                         cexRow = 0.5,
                         cexCol = 0.5, 
                         col = colorRampPalette(c("blue", "white", "red"))(300),
                         key.xlab = "Expression",
                         key.title = NA,
                         Colv  =T
) # Don't cluster columns
dev.off()


#####################################################################
#####################################################################
############################TCA cycle################################
#####################################################################
#####################################################################

gene_list <- c("Aco1", "Aco2", "Cs", "Fh1", "Gls", "Gls2", "Glul", "Idh1", "Idh2", "Idh3a", "Idh3b", "Idh3g", "Mdh1", "Mdh2", "Me1", "Me2", "Me3", "Mpc1", "Mpc2", "Ogdh", "Ogdhl", "Pdha1", "Pdhb", "Pdhx", "Pdk1", "Pdk2", "Pdk3", "Pdk4", "Pdp1", "Pdp2", "Pdpk1", "Pdpn", "Pdpr", "Sdha", "Sdhaf1", "Sdhaf2", "Sdhaf3", "Sdhaf4", "Sdhb", "Sdhc", "Sdhd")

# Extract the rows for the genes in the list
selected_genes <- gene_expression[rownames(gene_expression) %in% gene_list,]

# Transform the data to matrix, as required by heatmap.2
selected_genes_matrix <- as.matrix(selected_genes)

# Normalize the data by gene
#selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) x / max(x)))
#t(selected_genes_matrix_norm)
selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) scale(x)))
colnames(selected_genes_matrix_norm) <- colnames(selected_genes_matrix)

dev.off()  # Close the current device, if any
dev.new()  # Open a new graphics device

heatmap_plot = heatmap.2(selected_genes_matrix_norm, 
                         main = "TCA cycle", 
                         trace = "none", 
                         margins = c(5, 10), 
                         cexRow = 0.5,
                         cexCol = 0.5, 
                         col = colorRampPalette(c("blue", "white", "red"))(256),
                         key.xlab = "Expression",
                         key.title = NA,
                         Colv = T
                         
) # Don't cluster columns

#####################################################################
#####################################################################
########################Fatty acid oxidation#########################
#####################################################################
#####################################################################

gene_list <- c("Acaa1a", "Acaa2", "Acacb", "Acad10", "Acad11", "Acad12", "Acadl", "Acadm", "Acads", "Acadsb", "Acadvl", "Acat1", "Acot11", "Acsf2", "Acsl1", "Acsl6", "Acsm1", "Acsm2", "Acsm3", "Acsm4", "Acsm5", "Acss1", "Amacr", "Cpt1a", "Cpt1b", "Cpt1c", "Cpt2", "Crat", "Crot", "Decr1", "Echs1", "Eci1", "Eci2", "Etfa", "Etfb", "Etfdh", "Hadh", "Hadha", "Hadhb", "Hsd17b10", "Mcee", "Mmut", "Pcca", "Pccb", "Slc25a20")

# Extract the rows for the genes in the list
selected_genes <- gene_expression[rownames(gene_expression) %in% gene_list,]

# Transform the data to matrix, as required by heatmap.2
selected_genes_matrix <- as.matrix(selected_genes)

# Normalize the data by gene
#selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) x / max(x)))
#t(selected_genes_matrix_norm)
selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) scale(x)))
colnames(selected_genes_matrix_norm) <- colnames(selected_genes_matrix)

heatmap.2(selected_genes_matrix_norm, 
          main = "Fatty Acid Oxidation", 
          trace = "none", 
          margins = c(5, 10), 
          cexRow = 0.5,
          cexCol = 0.5, 
          col = colorRampPalette(c("blue", "white", "red"))(256),
          key.xlab = "Expression",
          key.title = NA,
          Colv = F
          
) # Don't cluster columns



#####################################################################
#####################################################################
########################Myc targets##################################
#####################################################################
#####################################################################

gene_list <- c("Abca1", "Adam10", "Adam17", "Akt1", "Ap3m1", "Atp11a", "Atp8a1", "Baz1b", "Bex3", "Calcrl", "Canx", "Cdca8", "Ceacam1", "Ceacam2", "Ddx3x", "Dr1", "Fgfr2", "Galnt7", "Ggta1", "Hbp1", "Hif1a", "Ifngr2", "Il13ra1", "Il6st", "Itga8", "Kif11", "Klf6", "Mpp5", "Msn", "Msrb3", "Mtmr1", "Myc", "Nedd4", "Nfatc3", "Nfia", "Pcdh18", "Pcyt1a", "Pitpnb", "Pkn2", "Prcp", "Reep3", "Runx1", "Scamp1", "Sec61a1", "Sema3c", "Sema3e", "Senp2", "Sgk3", "Sgpp1", "Shmt1", "Slc12a2", "Slc12a7", "Slc35a5", "Slc39a8", "Slc6a15", "Slc6a2", "Smpdl3b", "Sqle", "Suz12", "Tbl1x", "Tcf12", "Tfb2m", "Tfrc", "Timp2", "Tiparp", "Tmod3", "Txnip", "Ube3a", "Ulk1", "Utrn", "Zfp36l1", "Zmynd11")

# Extract the rows for the genes in the list
selected_genes <- gene_expression[rownames(gene_expression) %in% gene_list,]

# Transform the data to matrix, as required by heatmap.2
selected_genes_matrix <- as.matrix(selected_genes)

# Normalize the data by gene
#selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) x / max(x)))
#t(selected_genes_matrix_norm)
selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) scale(x)))
colnames(selected_genes_matrix_norm) <- colnames(selected_genes_matrix)

# Calculate the total expression (e.g., sum) for each sample
total_expression <- colSums(selected_genes_matrix_norm, na.rm = TRUE)

# Order the samples based on total expression (descending order)
ordered_samples <- names(sort(total_expression, decreasing = TRUE))

# Reorder the matrix columns based on expression levels
selected_genes_matrix_norm <- selected_genes_matrix_norm[, ordered_samples]

heatmap.2(selected_genes_matrix_norm, 
          main = "Myc Targets", 
          trace = "none", 
          margins = c(5, 10), 
          cexRow = 0.5,
          cexCol = 0.5, 
          col = colorRampPalette(c("blue", "white", "red"))(256),
          key.xlab = "Expression",
          key.title = NA,
          Colv = FALSE,
          
) # Don't cluster columns

#####################################################################
#####################################################################
##########################Apoptosis
#####################################################################
#####################################################################

gene_list <- c("Aifm1", "Aifm2", "Aifm3", "Bad", "Bak1", "Bax", "Bbc3", "Bcl2", "Bcl2a1d", "Bcl2l1", "Bcl2l10", "Bcl2l11", "Bcl2l13", "Bcl2l2", "Bid", "Bik", "Bnip3", "Bnip3l", "Bok", "Casp3", "Casp8", "Casp9", "Chchd2", "Cycs", "Diablo", "Endog", "Ghitm", "Htra2", "Ifi27", "Mcl1", "Pmaip1", "Septin4", "Sphk2", "Styxl1")
# Extract the rows for the genes in the list
selected_genes <- gene_expression[rownames(gene_expression) %in% gene_list,]

# Transform the data to matrix, as required by heatmap.2
selected_genes_matrix <- as.matrix(selected_genes)

# Normalize the data by gene
#selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) x / max(x)))
#t(selected_genes_matrix_norm)
selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) scale(x)))
colnames(selected_genes_matrix_norm) <- colnames(selected_genes_matrix)

# Calculate the total expression (e.g., sum) for each sample
total_expression <- colSums(selected_genes_matrix_norm, na.rm = TRUE)

# Order the samples based on total expression (descending order)
ordered_samples <- names(sort(total_expression, decreasing = TRUE))

# Reorder the matrix columns based on expression levels
selected_genes_matrix_norm <- selected_genes_matrix_norm[, ordered_samples]

dev.off()  # Close the current device, if any
dev.new()  # Open a new graphics device

heatmap_plot = heatmap.2(selected_genes_matrix_norm, 
                         main = "Apoptosis", 
                         trace = "none", 
                         margins = c(5, 10), 
                         cexRow = 0.5,
                         cexCol = 0.5, 
                         col = colorRampPalette(c("blue", "white", "red"))(256),
                         key.xlab = "Expression",
                         key.title = NA,
                         Colv = FALSE,
                         
) # Don't cluster columns

#####################################################################
#####################################################################
####Mitochondrial translation
#####################################################################
#####################################################################

gene_list <- c("2810006K23Rik", "Aars2", "Aurkaip1", "Cars2", "Chchd1", "Coa3", "Cox14", "Dap3", "Dars2", "Ddx28", "Dhx30", "Ears2", "Eral1", "Exd2", "Fars2", "Fastkd2", "Gadd45gip1", "Gars", "Gatb", "Gatc", "Gfm1", "Gfm2", "Grsf1", "Gtpbp10", "Guf1", "Hars2", "Hemk1", "Iars2", "Kars", "Lars2", "Lrpprc", "Malsu1", "Mars2", "Metap1d", "Mettl17", "Mief1", "Mpv17l2", "Mrm2", "Mrm3", "Mrpl1", "Mrpl10", "Mrpl11", "Mrpl12", "Mrpl13", "Mrpl14", "Mrpl15", "Mrpl16", "Mrpl17", "Mrpl18", "Mrpl19", "Mrpl2", "Mrpl20", "Mrpl21", "Mrpl22", "Mrpl23", "Mrpl24", "Mrpl27", "Mrpl28", "Mrpl3", "Mrpl30", "Mrpl32", "Mrpl33", "Mrpl34", "Mrpl35", "Mrpl36", "Mrpl37", "Mrpl38", "Mrpl39", "Mrpl4", "Mrpl40", "Mrpl41", "Mrpl42", "Mrpl43", "Mrpl44", "Mrpl45", "Mrpl46", "Mrpl47", "Mrpl48", "Mrpl49", "Mrpl50", "Mrpl51", "Mrpl52", "Mrpl53", "Mrpl54", "Mrpl55", "Mrpl57", "Mrpl58", "Mrpl9", "Mrps10", "Mrps11", "Mrps12", "Mrps14", "Mrps15", "Mrps16", "Mrps17", "Mrps18a", "Mrps18b", "Mrps18c", "Mrps2", "Mrps21", "Mrps22", "Mrps23", "Mrps24", "Mrps25", "Mrps26", "Mrps27", "Mrps28", "Mrps30", "Mrps31", "Mrps33", "Mrps34", "Mrps35", "Mrps36", "Mrps5", "Mrps6", "Mrps7", "Mrps9", "Mrrf", "Mterf1a", "Mterf3", "Mterf4", "Mtfmt", "Mtg1", "Mtg2", "Mtif2", "Mtif3", "Mtres1", "Mtrf1", "Mtrf1l", "Nars2", "Ngrn", "Noa1", "Nsun4", "Oxa1l", "Pars2", "Pdf", "Ppa2", "Ptcd3", "Pusl1", "Qrsl1", "Rars2", "Rbfa", "Rmnd")
# Extract the rows for the genes in the list

selected_genes <- gene_expression[rownames(gene_expression) %in% gene_list,]

# Transform the data to matrix, as required by heatmap.2
selected_genes_matrix <- as.matrix(selected_genes)

# Normalize the data by gene
#selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) x / max(x)))
#t(selected_genes_matrix_norm)
selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) scale(x)))
colnames(selected_genes_matrix_norm) <- colnames(selected_genes_matrix)

# Calculate the total expression (e.g., sum) for each sample
total_expression <- colSums(selected_genes_matrix_norm, na.rm = TRUE)

# Order the samples based on total expression (descending order)
ordered_samples <- names(sort(total_expression, decreasing = TRUE))

# Reorder the matrix columns based on expression levels
selected_genes_matrix_norm <- selected_genes_matrix_norm[, ordered_samples]

heatmap_plot = heatmap.2(selected_genes_matrix_norm, 
                         main = "Mitochondrial Translation", 
                         trace = "none", 
                         margins = c(5, 10), 
                         cexRow = 0.5,
                         cexCol = 0.5, 
                         col = colorRampPalette(c("blue", "white", "red"))(256),
                         key.xlab = "Expression",
                         key.title = NA,
                         Colv = FALSE,
                         
) # Don't cluster columns

#####################################################################
#####################################################################
######EAA transporter
#####################################################################
#####################################################################

gene_list <- c("Slc6a5", "Slc7a5", "Slc6a20b", "Slc7a11", "Slc7a8", "Slc38a2", "Slc38a4", "Slc38a9", "Slc38a10", "Slc38a3", "Slc43a1", "Slc6a20a", "Slc38a1", "Slc43a2", "Slc38a11", "Slc38a7", "Slc38a8", "Slc38a6", "Slc7a1", "Slc6a13", "Slc6a9", "Slc6a6", "Slc7a2")
selected_genes <- gene_expression[rownames(gene_expression) %in% gene_list,]

# Transform the data to matrix, as required by heatmap.2
selected_genes_matrix <- as.matrix(selected_genes)

# Normalize the data by gene
#selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) x / max(x)))
#t(selected_genes_matrix_norm)
selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) scale(x)))
colnames(selected_genes_matrix_norm) <- colnames(selected_genes_matrix)

# Calculate the total expression (e.g., sum) for each sample
total_expression <- colSums(selected_genes_matrix_norm, na.rm = TRUE)

# Order the samples based on total expression (descending order)
ordered_samples <- names(sort(total_expression, decreasing = TRUE))

# Reorder the matrix columns based on expression levels
selected_genes_matrix_norm <- selected_genes_matrix_norm[, ordered_samples]

# Remove rows with missing or invalid values
selected_genes_matrix_clean <- selected_genes_matrix_norm[complete.cases(selected_genes_matrix_norm), ]

heatmap_plot = heatmap.2(selected_genes_matrix_clean, 
                         main = "EAA transporters", 
                         trace = "none", 
                         margins = c(5, 10), 
                         cexRow = 0.5,
                         cexCol = 0.5, 
                         col = colorRampPalette(c("blue", "white", "red"))(256),
                         key.xlab = "Expression",
                         key.title = NA,
                         Colv = FALSE,
                         
) 
# Don't cluster columns

#####################################################################
#####################################################################
######All AA transporter
#####################################################################
#####################################################################

gene_list <- c("Slc3a1", "Slc3a2", "Slc1a3", "Slc1a2", "Slc1a1", "Slc1a6", "Slc1a7", "Slc6a9", "Slc6a5", "Slc6a7", "Slc6a14", "Slc6a15", "Slc6a17", "Slc6a18", "Slc6a19", "Slc6a20", "Slc38a1", "Slc38a2", "Slc38a3", "Slc38a4", "Slc38a5", "Slc36a1", "Slc36a2", "Slc36a3", "Slc36a4", "Slc1a4", "Slc1a5", "Slc7a5", "Slc7a8", "Slc7a7", "Slc7a6", "Slc7a9", "Slc7a10", "Slc7a11", "Slc7a13", "Slc7a1", "Slc7a2", "Slc7a3", "Slc7a4", "Slc16a10", "Slc43a1", "Slc43a2")
selected_genes <- gene_expression[rownames(gene_expression) %in% gene_list,]

# Transform the data to matrix, as required by heatmap.2
selected_genes_matrix <- as.matrix(selected_genes)

# Normalize the data by gene
#selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) x / max(x)))
#t(selected_genes_matrix_norm)
selected_genes_matrix_norm <- t(apply(selected_genes_matrix, 1, function(x) scale(x)))
colnames(selected_genes_matrix_norm) <- colnames(selected_genes_matrix)

# Calculate the total expression (e.g., sum) for each sample
total_expression <- colSums(selected_genes_matrix_norm, na.rm = TRUE)

# Order the samples based on total expression (descending order)
ordered_samples <- names(sort(total_expression, decreasing = TRUE))

# Reorder the matrix columns based on expression levels
selected_genes_matrix_norm <- selected_genes_matrix_norm[, ordered_samples]

# Remove rows with missing or invalid values
selected_genes_matrix_clean <- selected_genes_matrix_norm[complete.cases(selected_genes_matrix_norm), ]

heatmap_plot = heatmap.2(selected_genes_matrix_clean, 
                         main = "All AA transporters", 
                         trace = "none", 
                         margins = c(5, 10), 
                         cexRow = 0.5,
                         cexCol = 0.5, 
                         col = colorRampPalette(c("blue", "white", "red"))(256),
                         key.xlab = "Expression",
                         key.title = NA,
                         Colv = FALSE,
                         
) 
# Don't cluster columns

