
# Libraries
library(tidyverse)
library(biomaRt)
library(limma)
library(ArrayTools)

# PRE-RANKED GSEA
# Load data
setwd("/home/mario/Projects/eosinophils")
dataset <- read.csv("selected_genes.txt", sep="\t")

# Change Emsembl by Entrez IDs
mart <- useEnsembl(biomart='ensembl', dataset='mmusculus_gene_ensembl')
translation <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id'),
                    filters = 'ensembl_gene_id',
                    values = dataset$ensembl_gene_id,
                    mart = mart)
translation <- translation[!duplicated(translation$ensembl_gene_id),]
dataset <- merge(x = dataset, y = translation, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")

# Create ranked file
rnk <- dataset[, c(5, 4)]
rnk <- rnk[order(-rnk$logFC), ]
# Drop NA and empty values
rnk <- rnk %>% drop_na()
rnk <- rnk[rnk[, 1] != "", ]
# Write data
write.table(rnk, file="selected_gene.rnk", quote = F, sep = "\t", row.names = F,
            col.names = F)

# # GSEA
# # Load data
# setwd('/home/mario/Projects/eosinophils/GSE62999_RAW')
# cel_files <- list.files(path = getwd(), pattern = '*.CEL.gz', full.names = TRUE)
# setwd('/home/mario/Projects/eosinophils')
# parsed_cels <- oligo::read.celfiles(cel_files, verbose = TRUE) 
# parsed_cels_rma <- oligo::rma(parsed_cels, normalize = TRUE, background = TRUE)  
# 
# # Create the .gct and the .cls
# filename <- "probes"
# output.gct(parsed_cels_rma, filename)
# output.cls(parsed_cels_rma, filename)
