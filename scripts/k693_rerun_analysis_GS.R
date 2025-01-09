rm(list=ls())

library(DESeq2)
library(dplyr)
library(tibble)


###### Read in and set rownames for count matrix data and sample information

coldata <- read.delim("data/FULL.DAT.COL.DATA.txt", sep = "\t", header = TRUE)
cts = read.csv("data/FULL.DAT.csv")
cts <- cts[-2]
cts = cts %>% remove_rownames %>% column_to_rownames(var="gene_id")
coldata = coldata %>% remove_rownames %>% column_to_rownames(var="sample")

# Check the matches here before loading the data for further analysis
all(rownames(coldata) %in% colnames(cts)) # This should return TRUE

# Set the conditions as factors in the coldata
coldata$cancer <- as.factor(coldata$age)
