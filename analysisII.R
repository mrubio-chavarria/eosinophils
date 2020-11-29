
# The second version started on 23/11/2020

# Libraries
library(limma)
library(dplyr)

# Set working directory and load data
setwd('/media/mario/E788-1B77/Projects/eosinophils')
dataset <- read.csv('dataset.csv')

# Design matrix
# Key:
# Wild type --> 0
# Knockout --> 1
# Untreated --> 0
# Treated --> 1
# wild type untreated
genotype <- rep(0, 5)
treatment <- rep(0, 5)
part_data_1 <- data.frame(genotype, treatment)
# wild type treated
genotype <- rep(0, 5)
treatment <- rep(1, 5)
part_data_2 <- data.frame(genotype, treatment)
# knockout untreated
genotype <- rep(1, 5)
treatment <- rep(0, 5)
part_data_3 <- data.frame(genotype, treatment)
# knockout treated
genotype <- rep(1, 5)
treatment <- rep(1, 5)
part_data_4 <- data.frame(genotype, treatment)
design_data <- rbind(part_data_1, part_data_2, part_data_3, part_data_4)
model <- as.formula("~genotype + treatment + genotype:treatment")
design_matrix <- model.matrix(model, data = design_data)
colnames(design_matrix)[1] <- "Intercept"
colnames(design_matrix)[4] <- "genotype_treatment"

labels <- dataset$Ensembl
values <- dataset[, 2:21]

fit1 <- lmFit(values)

contrast.matrix <- makeContrasts(genotype - treatment,
                                 genotype_treatment,
                                 levels = design_matrix)

fit2 <- contrasts.fit(fit1, contrast.matrix)

fit3 <- eBayes(contrast1)

deg <- topTable(fit3, p.value=0.05)




