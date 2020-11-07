
library(tidyverse)
library(dplyr)

# PREPARE
setwd('/media/aetius/STORE N GO/Projects/eosinophils/')
raw_data <- read.csv(file = 'GSE62999_RAW/GSM1537865_RMA_processed_log2_Ensembl_05112020.csv')

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
confidence <- 0.99
p_value <- 1 - confidence

################################################################################
# Analysis 1
# Legend
# All: untreated
# Differences: wild type vs knockout
sig_genes <- c()
no_sig_genes <- c()
# Perform the tests to detect significant genes
for (i in seq(1, n_genes)) {
  gene_data <- filter(eo_data,
                      gene_id == gene_ids[i],
                      !treated)
  group1 <- filter(gene_data, genotype == 'WT')$intensity
  group2 <- filter(gene_data, genotype == 'KO')$intensity
  test <- t.test(group1, group2, alt='two.sided', conf=confidence, var.eq = F, paired = F)
  if (test$p.value <= p_value) {
    sig_genes[i] <- gene_ids[i]
  } else {
    no_sig_genes[i] <- gene_ids[i]
  }
}
unt_wtko_sig_genes <- sig_genes[!is.na(sig_genes)]
unt_wtko_no_sig_genes <- no_sig_genes[!is.na(no_sig_genes)]
################################################################################

################################################################################
# Analysis 2
# Legend
# All: treated
# Differences: wild type vs knockout
sig_genes <- c()
no_sig_genes <- c()
# Perform the tests to detect significant genes
for (i in seq(1, n_genes)) {
  gene_data <- filter(eo_data,
                      gene_id == gene_ids[i],
                      treated)
  group1 <- filter(gene_data, genotype == 'WT')$intensity
  group2 <- filter(gene_data, genotype == 'KO')$intensity
  test <- t.test(group1, group2, alt='two.sided', conf=confidence, var.eq = F, paired = F)
  if (test$p.value <= p_value) {
    sig_genes[i] <- gene_ids[i]
  } else {
    no_sig_genes[i] <- gene_ids[i]
  }
}
t_wtko_sig_genes <- sig_genes[!is.na(sig_genes)]
t_wtko_no_sig_genes <- no_sig_genes[!is.na(no_sig_genes)]
################################################################################

################################################################################
# Analysis 3
# Legend
# All: wild type
# Differences: treated vs untreated
sig_genes <- c()
no_sig_genes <- c()
# Perform the tests to detect significant genes
for (i in seq(1, n_genes)) {
  gene_data <- filter(eo_data,
                      gene_id == gene_ids[i],
                      genotype == 'WT')
  group1 <- filter(gene_data, treated)$intensity
  group2 <- filter(gene_data, !treated)$intensity
  test <- t.test(group1, group2, alt='two.sided', conf=confidence, var.eq = F, paired = T)
  if (test$p.value <= p_value) {
    sig_genes[i] <- gene_ids[i]
  } else {
    no_sig_genes[i] <- gene_ids[i]
  }
}
wt_tunt_sig_genes <- sig_genes[!is.na(sig_genes)]
wt_tunt_no_sig_genes <- no_sig_genes[!is.na(no_sig_genes)]
################################################################################

################################################################################
# Analysis 4
# Legend
# All: knockout
# Differences: treated vs untreated
sig_genes <- c()
no_sig_genes <- c()
# Perform the tests to detect significant genes
for (i in seq(1, n_genes)) {
  gene_data <- filter(eo_data,
                      gene_id == gene_ids[i],
                      genotype == 'KO')
  group1 <- filter(gene_data, treated)$intensity
  group2 <- filter(gene_data, !treated)$intensity
  test <- t.test(group1, group2, alt='two.sided', conf=confidence, var.eq = F, paired = T)
  if (test$p.value <= p_value) {
    sig_genes[i] <- gene_ids[i]
  } else {
    no_sig_genes[i] <- gene_ids[i]
  }
}
ko_tunt_sig_genes <- sig_genes[!is.na(sig_genes)]
ko_tunt_no_sig_genes <- no_sig_genes[!is.na(no_sig_genes)]
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

