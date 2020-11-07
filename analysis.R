
library(tidyverse)
library(dplyr)

# PREPARE
setwd('/home/aetius/Projects/eosinophils/')
raw_data <- read.csv(file = 'GSE62999_RAW/GSM1537865_RMA_processed_log2_Ensembl_05112020.csv')

# TIDY
gene_ids <- raw_data[, 1]
colnames(data) <- substring(colnames(data), 23, 37)

# Wild types
wt_mouses <- list()
for (i in seq(1, 5)) {
  untreated <- data.frame(intensity = raw_data[, 1 + i], 
                          gene_id = gene_ids,
                          mouse_id = i,
                          treated = F,
                          genotype = 'WT')
  treated <- data.frame(intensity = raw_data[, 11 + i],
                        gene_id = gene_ids,
                        mouse_id = i,
                        treated = T,
                        genotype = 'WT')
  mouse <- rbind(untreated, treated)
  wt_mouses[[i]] <- mouse
}
wt_mouses <- do.call(rbind, wt_mouses)

# Knockout
ko_mouses <- list()
for (i in seq(1, 5)) {
  untreated <- data.frame(intensity = raw_data[, 6 + i],
                                gene_id = gene_ids,
                                mouse_id = i,
                                treated = F,
                                genotype = 'KO')
  treated <- data.frame(intensity = raw_data[, 16 + i],
                              gene_id = gene_ids,
                              mouse_id = i,
                              treated = T,
                              genotype = 'KO')
  mouse <- rbind(untreated, treated)
  ko_mouses[[i]] <- mouse
}
ko_mouses <- do.call(rbind, ko_mouses)

eo_data <- rbind(wt_mouses, ko_mouses)

