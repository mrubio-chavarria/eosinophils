
# Libraries
library(tidyverse)
library(biomaRt)
library(limma)
library(ArrayTools)
library(BWStest)
library(onewaytests)

# PRE-RANKED GSEA
# Set directory
setwd("/home/mario/Projects/eosinophils")
# Code taken from the other script
expression_data <-  read.csv('GSM1537865_RMA_processed_log2_Ensembl_09112020.csv', row.names=1)
colnames(expression_data) <- gsub('(GSM[0-9]+_MA[0-9]+).*', '\\1', colnames(expression_data))
metadata <- read.table('E-GEOD-62999.sdrf.txt', sep = '\t', header=T)
metadata$Simple_names <- gsub('(GSM[0-9]+_MA[0-9]+).*', '\\1', metadata$Array.Data.File)
metadata <- metadata[match(colnames(expression_data),metadata$Simple_names),]
metadata <- metadata[,c('Simple_names', 'FactorValue..genotype.', 'FactorValue..treatment.')]
colnames(metadata) <- c('Sample_name', 'DUSP5_Genotype', 'IL33_Treatment')
metadata$DUSP5_Genotype<- sapply(strsplit(metadata$DUSP5_Genotype, ' '), '[', 1)
metadata$IL33_Treatment[metadata$IL33_Treatment == 'IL-33'] <- 'Treated'
metadata$DUSP5_Genotype <- factor(metadata$DUSP5_Genotype, levels=c('WT', 'KO'))
metadata$IL33_Treatment <- factor(metadata$IL33_Treatment, levels=c('Untreated', 'Treated'))
model_matrix <- model.matrix(~1+DUSP5_Genotype+IL33_Treatment+DUSP5_Genotype:IL33_Treatment, data=metadata)
limma_fit <- lmFit(expression_data, model_matrix)
limma_fit <- eBayes(limma_fit, trend=TRUE, robust=TRUE)
treatment_results <-topTable(limma_fit, number = Inf, coef = 3, sort.by = 'p')
treatment_results$ensembl_gene_id <- row.names(treatment_results)
# mart <- useEnsembl(biomart='ensembl', dataset='mmusculus_gene_ensembl')
# translation <- getBM(attributes = c('ensembl_gene_id', 'mgi_symbol', 'description'),
#                     filters = 'ensembl_gene_id',  
#                     values = row.names(expression_data),  
#                     mart = mart)
# translation <- translation[!duplicated(translation$ensembl_gene_id),]
treatment_results <- merge(translation, treatment_results, by='ensembl_gene_id')
treatment_results <- treatment_results[order(treatment_results$adj.P.Val),]
# The section of my function
fit <- limma_fit
coef <- 3
results <- topTable(fit, number = Inf, coef = coef, sort.by = 'p')
results$ensembl_gene_id <- row.names(results)
results <- merge(translation, results, by='ensembl_gene_id')
results <- results[order(results$P.Value),]
relevant.b <- results[order(results$B),][1:20, ]  %>%
  dplyr::select(mgi_symbol, logFC, B) %>% mutate(mlog10B = -log(abs(B), 10))
results <- results %>% dplyr::select(ensembl_gene_id, mgi_symbol, logFC, P.Value, B) %>% 
  mutate(adj.P.Val = p.adjust(P.Value, method = 'BH')) %>%
  mutate(mlog10P = -log(P.Value, 10),
         mlog10AdjP = -log(adj.P.Val, 10),
         mlog10B = -log(abs(B), 10))
relevants <- results
relevants <- relevants[c("ensembl_gene_id", "P.Value", "adj.P.Val", "logFC")]
final_data <- merge(x = relevants, y = translation, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
final_data$description <- NULL
# Create ranked file
rnk <- final_data[, c(5, 4)]
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
