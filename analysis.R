
library(Biobase)
library(tidyverse)
library(dplyr)
library(limma)


# PREPARE
setwd('/media/aetius/STORE N GO/Projects/eosinophils/')
filename <- 'GSE62999_RAW/GSM1537865_RMA_processed_log2_Ensembl_05112020.csv'
raw_data <- read.csv(file = filename)

# TIDY
gene_ids <- raw_data[, 1]
n_genes <- length(gene_ids)
colnames(raw_data) <- substring(colnames(raw_data), 23, 37)

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

# ANALYSES

# Fix t-test parameters
genes_p_values <- data.frame(gene_id = gene_ids,
                       analysis_1 = seq(1:n_genes),
                       analysis_2 = seq(1:n_genes),
                       analysis_3 = seq(1:n_genes),
                       analysis_4 = seq(1:n_genes))
################################################################################
# Analysis 1
# Legend
# All: untreated
# Differences: wild type vs knockout
# Perform the tests to detect significant genes
for (i in seq(1, n_genes)) {
  gene_data <- filter(eo_data,
                      gene_id == gene_ids[i],
                      !treated)
  group1 <- filter(gene_data, genotype == 'WT')$intensity
  group2 <- filter(gene_data, genotype == 'KO')$intensity
  test <- t.test(group1, group2, alt='two.sided', var.eq = F, paired = F)
  genes_p_values[i, 2] <- test$p.value
}
################################################################################

################################################################################
# Analysis 2
# Legend
# All: treated
# Differences: wild type vs knockout
# Perform the tests to detect significant genes
for (i in seq(1, n_genes)) {
  gene_data <- filter(eo_data,
                      gene_id == gene_ids[i],
                      treated)
  group1 <- filter(gene_data, genotype == 'WT')$intensity
  group2 <- filter(gene_data, genotype == 'KO')$intensity
  test <- t.test(group1, group2, alt='two.sided', var.eq = F, paired = F)
  genes_p_values[i, 3] <- test$p.value
}
################################################################################

################################################################################
# Analysis 3
# Legend
# All: wild type
# Differences: treated vs untreated
# Perform the tests to detect significant genes
for (i in seq(1, n_genes)) {
  gene_data <- filter(eo_data,
                      gene_id == gene_ids[i],
                      genotype == 'WT')
  group1 <- filter(gene_data, treated)$intensity
  group2 <- filter(gene_data, !treated)$intensity
  test <- t.test(group1, group2, alt='two.sided', var.eq = F, paired = T)
  genes_p_values[i, 4] <- test$p.value
}
################################################################################

################################################################################
# Analysis 4
# Legend
# All: knockout
# Differences: treated vs untreated
# Perform the tests to detect significant genes
for (i in seq(1, n_genes)) {
  gene_data <- filter(eo_data,
                      gene_id == gene_ids[i],
                      genotype == 'KO')
  group1 <- filter(gene_data, treated)$intensity
  group2 <- filter(gene_data, !treated)$intensity
  test <- t.test(group1, group2, alt='two.sided', var.eq = F, paired = T)
  genes_p_values[i, 5] <- test$p.value
}
################################################################################

# Create the table with the means of the intensities by gene and group
names <- c('gene_id', 'genotype', 'treated', 'avg_intensity')
avg_intensities <- data.frame(character(), character(), logical(), double())
colnames(avg_intensities) <- names

# Obtain the average intensities for every gene
for (i in seq(1, n_genes)) {
  gene <- gene_ids[i]
  related_rows <- eo_data %>%
    filter(gene_id == gene)
  intensities <- related_rows %>% 
    group_by(genotype, treated) %>%
    summarise(avg_intensity = mean(intensity)) %>%
    ungroup()
  # It is 4 because there is only 4 groups
  for (j in seq(1, 4)) {
    new_row <- c(gene)
    new_row[2:4] <- intensities[j, 1:3]
    new_row <- data.frame(new_row)
    colnames(new_row) <- names
    avg_intensities <- rbind(avg_intensities, new_row)
  }
}

# Create the heatmap with the most important genes
gene_expression_matrix <- raw_data
colnames(gene_expression_matrix)[2:21] <- paste('Sample', seq(1:20))
colnames(gene_expression_matrix)[1] <- 'gene_id'
# Genotype heatmap
# Take the n genes with the highest difference
n <- 40
selected <- genes_p_values[, c(1, 2, 3)] %>%
  transmute(gene_id = gene_id, total = analysis_1 + analysis_2) 
most_diff_genes <- selected[order(selected$total, decreasing = F), ][1:n, 1]
genes_matrix <- gene_expression_matrix %>% filter(gene_id %in% most_diff_genes)
genes_matrix <-as.matrix(genes_matrix[, 2:21])
rownames(genes_matrix) <- paste('Gene', seq(1:n))
coolmap(genes_matrix, cexRow = 1)
# Treatment heatmap
# Take the n genes with the highest difference
n <- 40
selected <- genes_p_values[, c(1, 4, 5)] %>%
  transmute(gene_id = gene_id, total = analysis_3 + analysis_4) 
most_diff_genes <- selected[order(selected$total, decreasing = F), ][1:n, 1]
genes_matrix <- gene_expression_matrix %>% filter(gene_id %in% most_diff_genes)
genes_matrix <-as.matrix(genes_matrix[, 2:21])
rownames(genes_matrix) <- paste('Gene', seq(1:n))
coolmap(genes_matrix, cexRow = 1)
