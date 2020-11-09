
library(Biobase)
library(tidyverse)
library(dplyr)
library(limma)


# PREPARE
setwd('/media/aetius/STORE N GO/Projects/eosinophils/')
filename <- 'GSE62999_RAW/GSM1537865_RMA_processed_log2_Ensembl_09112020.csv'
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

# Obtain the average intensities for every gene and conditions
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

# Calculate the differences between situations, e.g. treated-KO vs utreated-WT 
# and so on
confidence <- 0.95
p_value <- 1 - confidence
n_shown <- 20 # genes to be shown

################################################################################
# Difference 1: 
# Compared: genotypes
# Common: untreated
# Analysis: 1
analysis <- 'analysis_1'
# Calculate the comparisons and relative values
interest_genes <- genes_p_values[, c('gene_id', analysis)] %>%
  filter(analysis_1 <= p_value)
non_interest_genes <- genes_p_values[, c('gene_id', analysis)] %>%
  filter(analysis_1 > p_value)
interest_genes <- interest_genes$gene_id
condition_intensities <- avg_intensities %>% filter(gene_id %in% interest_genes,
                                                    !treated)
condition_intensities$treated <- NULL
group_1 <- condition_intensities %>% filter(genotype == 'WT')
group_2 <- condition_intensities %>% filter(genotype == 'KO')
compared_intensities <- data.frame(gene = unique(condition_intensities$gene_id),
                                   comparison = group_1$avg_intensity / group_2$avg_intensity)
compared_intensities$log_comparison <- log(compared_intensities$comparison)
compared_intensities$squared_log_comparison <- compared_intensities$log_comparison ** 2 
compared_intensities <- compared_intensities[order(compared_intensities$squared_log_comparison,
                                                   decreasing = T), ]
compared_intensities <- head(compared_intensities, n_shown + 1)
compared_intensities$relative <- compared_intensities$comparison / max(compared_intensities$comparison)
compared_intensities$relative <- 100 * compared_intensities$relative
compared_intensities <- compared_intensities[order(compared_intensities$relative, decreasing = T), ]
compared_intensities$expression <- sign(compared_intensities$log_comparison) * compared_intensities$relative
diff_1_genes <- compared_intensities$gene
# Plot the differences
ggplot(data = compared_intensities[2:21, ],
       mapping = aes(x = reorder(gene, -relative), y = expression, fill = relative)) + 
  geom_bar(stat = 'identity') + coord_cartesian(ylim = c(-75, 75)) +
  ylab('Relative expression') + xlab('Genes') + labs(fill = "Fraction (Dusp5)")
################################################################################

################################################################################
# Difference 2: 
# Compared: genotypes
# Common: treated
# Analysis: 2
analysis <- 'analysis_2'
# Calculate the comparisons and relative values
interest_genes <- genes_p_values[, c('gene_id', analysis)] %>%
  filter(analysis_2 <= p_value)
non_interest_genes <- genes_p_values[, c('gene_id', analysis)] %>%
  filter(analysis_2 > p_value)
interest_genes <- interest_genes$gene_id
condition_intensities <- avg_intensities %>% filter(gene_id %in% interest_genes,
                                                    treated)
condition_intensities$treated <- NULL
group_1 <- condition_intensities %>% filter(genotype == 'WT')
group_2 <- condition_intensities %>% filter(genotype == 'KO')
compared_intensities <- data.frame(gene = unique(condition_intensities$gene_id),
                                   comparison = group_1$avg_intensity / group_2$avg_intensity)
compared_intensities$log_comparison <- log(compared_intensities$comparison)
compared_intensities$squared_log_comparison <- compared_intensities$log_comparison ** 2 
compared_intensities <- compared_intensities[order(compared_intensities$squared_log_comparison,
                                                   decreasing = T), ]
compared_intensities <- head(compared_intensities, n_shown + 1)
compared_intensities$relative <- compared_intensities$comparison / max(compared_intensities$comparison)
compared_intensities$relative <- 100 * compared_intensities$relative
compared_intensities <- compared_intensities[order(compared_intensities$relative, decreasing = T), ]
compared_intensities$expression <- sign(compared_intensities$log_comparison) * compared_intensities$relative
diff_2_genes <- compared_intensities$gene
# Plot the differences
ggplot(data = compared_intensities[2:21, ],
       mapping = aes(x = reorder(gene, -relative), y = expression, fill = relative)) + 
  geom_bar(stat = 'identity') + coord_cartesian(ylim = c(-75, 75)) +
  ylab('Relative expression') + xlab('Genes') + labs(fill = "Fraction (Dusp5)")
################################################################################

################################################################################
# Difference 3: 
# Compared: treatments
# Common: wild type
# Analysis: 3
analysis <- 'analysis_3'
# Calculate the comparisons and relative values
interest_genes <- genes_p_values[, c('gene_id', analysis)] %>%
  filter(analysis_3 <= p_value)
non_interest_genes <- genes_p_values[, c('gene_id', analysis)] %>%
  filter(analysis_3 > p_value)
interest_genes <- interest_genes$gene_id
condition_intensities <- avg_intensities %>% filter(gene_id %in% interest_genes,
                                                    genotype == 'WT')
condition_intensities$genotype <- NULL
group_1 <- condition_intensities %>% filter(treated)
group_2 <- condition_intensities %>% filter(!treated)
compared_intensities <- data.frame(gene = unique(condition_intensities$gene_id),
                                   comparison = group_1$avg_intensity / group_2$avg_intensity)
compared_intensities$log_comparison <- log(compared_intensities$comparison)
compared_intensities$squared_log_comparison <- compared_intensities$log_comparison ** 2 
compared_intensities <- compared_intensities[order(compared_intensities$squared_log_comparison,
                                                   decreasing = T), ]
compared_intensities <- head(compared_intensities, n_shown + 1)
compared_intensities$relative <- compared_intensities$comparison / max(compared_intensities$comparison)
compared_intensities$relative <- 100 * compared_intensities$relative
compared_intensities <- compared_intensities[order(compared_intensities$relative, decreasing = T), ]
compared_intensities$expression <- sign(compared_intensities$log_comparison) * compared_intensities$relative
diff_3_genes <- compared_intensities$gene
# Plot the differences
ggplot(data = compared_intensities[1:21, ],
       mapping = aes(x = reorder(gene, -relative), y = expression, fill = relative)) + 
  geom_bar(stat = 'identity') + coord_cartesian(ylim = c(0, 105)) +
  ylab('Relative expression') + xlab('Genes') + labs(fill = "Fraction (Ptgs2)")
################################################################################

################################################################################
# Difference 4: 
# Compared: treatments
# Common: wild type
# Analysis: 4
analysis <- 'analysis_4'
# Calculate the comparisons and relative values
interest_genes <- genes_p_values[, c('gene_id', analysis)] %>%
  filter(analysis_4 <= p_value)
non_interest_genes <- genes_p_values[, c('gene_id', analysis)] %>%
  filter(analysis_4 > p_value)
interest_genes <- interest_genes$gene_id
condition_intensities <- avg_intensities %>% filter(gene_id %in% interest_genes,
                                                    genotype == 'KO')
condition_intensities$genotype <- NULL
group_1 <- condition_intensities %>% filter(treated)
group_2 <- condition_intensities %>% filter(!treated)
compared_intensities <- data.frame(gene = unique(condition_intensities$gene_id),
                                   comparison = group_1$avg_intensity / group_2$avg_intensity)
compared_intensities$log_comparison <- log(compared_intensities$comparison)
compared_intensities$squared_log_comparison <- compared_intensities$log_comparison ** 2 
compared_intensities <- compared_intensities[order(compared_intensities$squared_log_comparison,
                                                   decreasing = T), ]
compared_intensities <- head(compared_intensities, n_shown + 1)
compared_intensities$relative <- compared_intensities$comparison / max(compared_intensities$comparison)
compared_intensities$relative <- 100 * compared_intensities$relative
compared_intensities <- compared_intensities[order(compared_intensities$relative, decreasing = T), ]
compared_intensities$expression <- sign(compared_intensities$log_comparison) * compared_intensities$relative
diff_4_genes <- compared_intensities$gene
# Plot the differences
ggplot(data = compared_intensities[1:21, ],
       mapping = aes(x = reorder(gene, -relative), y = expression, fill = relative)) + 
  geom_bar(stat = 'identity') + coord_cartesian(ylim = c(0, 105)) +
  ylab('Relative expression') + xlab('Genes') + labs(fill = "Fraction (Ptgs2)")
################################################################################

# Prepare data of the genes with the biggest differences
diff_genes <- unique(c(diff_1_genes, diff_2_genes, diff_3_genes, diff_4_genes))
gene_expression_matrix <- raw_data
colnames(gene_expression_matrix)[2:21] <- paste('Sample', seq(1:20))
colnames(gene_expression_matrix)[1] <- 'gene_id'
# Build the heatmap
genes_matrix <- gene_expression_matrix %>% filter(gene_id %in% diff_genes)
genes_matrix <-as.matrix(genes_matrix[, 2:21])
rownames(genes_matrix) <- diff_genes
coolmap(genes_matrix, cexRow = 0.6)
