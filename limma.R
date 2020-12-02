
library(limma)
library(ggplot2)
library(biomaRt)
library(dplyr)
library(gridExtra)

setwd('/media/mario/2860-9FFC/Projects/eosinophils')

# Load our processed microarray data
expression_data = read.csv('GSM1537865_RMA_processed_log2_Ensembl_09112020.csv', row.names=1)
# Simplify names of samples
colnames(expression_data) <- gsub('(GSM[0-9]+_MA[0-9]+).*', '\\1', colnames(expression_data))

# Load metadata (downloaded from an alternative source: https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-62999/)
metadata = read.table('E-GEOD-62999.sdrf.txt', sep = '\t', header=T)
# make short sample names for matching to expression
metadata$Simple_names = gsub('(GSM[0-9]+_MA[0-9]+).*', '\\1', metadata$Array.Data.File)

# Order metadata as array data and simplify
metadata = metadata[match(colnames(expression_data),metadata$Simple_names),]
metadata = metadata[,c('Simple_names', 'FactorValue..genotype.', 'FactorValue..treatment.')]
colnames(metadata) = c('Sample_name', 'DUSP5_Genotype', 'IL33_Treatment')
metadata$DUSP5_Genotype = sapply(strsplit(metadata$DUSP5_Genotype, ' '), '[', 1)
metadata$IL33_Treatment[metadata$IL33_Treatment == 'IL-33'] = 'Treated'

# Change to factors
metadata$DUSP5_Genotype = factor(metadata$DUSP5_Genotype, levels=c('WT', 'KO'))  # this sets up BASELINE as WT (will be intercept)
metadata$IL33_Treatment = factor(metadata$IL33_Treatment, levels=c('Untreated', 'Treated'))  # this sets up BASELINE as UNTREATED (will be intercept)

# Now ready to set up the factors
model_matrix = model.matrix(~1+DUSP5_Genotype+IL33_Treatment+DUSP5_Genotype:IL33_Treatment, data=metadata)

# run limma with eBayes
limma_fit = lmFit(expression_data, model_matrix)
limma_fit = eBayes(limma_fit, trend=TRUE, robust=TRUE)

# extract results without using contrasts (as we're using intercept already)
# check the ordering of those coefficients to use the correct number
colnames(limma_fit$coefficients)  # e.g. coef=2 is the effect of the DUSP5 KO

KO_results = topTable(limma_fit, number = Inf, coef = 2, sort.by = 'p')
KO_results$ensembl_gene_id = row.names(KO_results)

treatment_results = topTable(limma_fit, number = Inf, coef = 3, sort.by = 'p')
treatment_results$ensembl_gene_id = row.names(treatment_results)

interaction_results = topTable(limma_fit, number = Inf, coef = 4, sort.by = 'p')
interaction_results$ensembl_gene_id = row.names(interaction_results)

# Everything below p adjusted value of 0.1 is typically consider significant
# log2FC cutoff is typically in the 0.5-1.0 range (absolute value since it can be up or down)
# E.g. I would consider the two top genes for the interaction term as significant

# Let's annotated the data with nicer gene names using biomaRt
mart = useEnsembl(biomart='ensembl', dataset='mmusculus_gene_ensembl')
translation = getBM(attributes = c('ensembl_gene_id', 'mgi_symbol', 'description'), # this specifies the fields we want to obtain
                    filters = 'ensembl_gene_id',  # this specifies what field we are selecting rows on (passed on in values attribute)
                    values = row.names(expression_data),  # here we actually say what genes we need data for
                    mart = mart)
translation = translation[!duplicated(translation$ensembl_gene_id),]  # sometimes IDs map to more genes, we can skip it now

# Add the descriptions and names to results
KO_results = merge(translation, KO_results, by='ensembl_gene_id')
KO_results = KO_results[order(KO_results$adj.P.Val),]

treatment_results = merge(translation, treatment_results, by='ensembl_gene_id')
treatment_results = treatment_results[order(treatment_results$adj.P.Val),]

interaction_results = merge(translation, interaction_results, by='ensembl_gene_id')
interaction_results = interaction_results[order(interaction_results$adj.P.Val),]

# Plot selected results and double-check if the values make sense
plot_gene = function(ensembl_id){
  
  # Extract and annotate
  expression_gene = expression_data[row.names(expression_data) == ensembl_id,]
  expression_gene = data.frame(t(expression_gene))
  colnames(expression_gene) = 'Log2_gene_expression'
  expression_gene$Sample = row.names(expression_gene)
  expression_gene = merge(metadata, expression_gene, by.x='Sample_name', by.y='Sample')
  expression_gene$Group = paste(expression_gene$DUSP5_Genotype, expression_gene$IL33_Treatment, sep='_')
  expression_gene$Group = factor(expression_gene$Group, levels=c('WT_Untreated', 'WT_Treated', 'KO_Untreated', 'KO_Treated'))
  expression_gene = expression_gene[order(expression_gene$Group),]
  gene_MGI_name = translation$mgi_symbol[translation$ensembl_gene_id == ensembl_id]
  
  # plot
  ggplot(expression_gene) + geom_point(aes(x=Group, y=Log2_gene_expression, color=Group)) +
    ggtitle(paste('Gene expression values for gene: ', gene_MGI_name, sep=''))
  
}

# Plot the top gene from KO - DUSP5 should have very low signal (background mostly)
plot_gene(KO_results$ensembl_gene_id[3])

# Plot the top gene affected by IL-33 treatment
plot_gene(treatment_results$ensembl_gene_id[1])

# Plot the gene most affected by interaction
plot_gene(interaction_results$ensembl_gene_id[1])


# Another classical plot in gene expression is a volcano plot
# It shows the effects of one factor (e.g. IL33-treatment vs control)
# The X-axis is Log2FC and the Y-axis is -log10(p.adjusted)
# Each dot is one gene (you should use the topTable results for this)
# You can highlight with color genes that are significant, e.g. p.adj < 0.1 and abs(log2FC) > 0.5
plot_volcano = function(fit, coef){
  ## Write your code
  results = topTable(fit, number = Inf, coef = coef, sort.by = 'p')
  results$ensembl_gene_id = row.names(results)
  results = merge(translation, results, by='ensembl_gene_id')
  results = results[order(results$P.Value),]
  relevant.b <- results[order(results$B),][1:20, ]  %>%
    select(mgi_symbol, logFC, B) %>% mutate(mlog10B = -log(abs(B), 10))
  results <- results %>% select(mgi_symbol, logFC, P.Value, B) %>% 
    mutate(adj.P.Val = p.adjust(P.Value, method = 'BH')) %>%
    mutate(mlog10P = -log(P.Value, 10),
           mlog10AdjP = -log(adj.P.Val, 10),
           mlog10B = -log(abs(B), 10))
  # Correction for Dusp5. In this case this point does not provide any 
  # information and it deforms the whole chart.
  if (coef == 2) {
    results <- results %>% filter(mgi_symbol != 'Dusp5')
  }
  
  # Take the relevant genes
  # This varies depending on the coefficient because the efect of the treatment
  # is much bigger compared to the other two
  if (coef == 2) {
    # Criterion:  adjusted p.val < 0.1 and abs(log2FC) > 0.8
    relevants <- results[1:20, ] %>% filter(adj.P.Val < 0.1, abs(logFC) > 0.75)
    title <- 'Effec of the Dusp5 knockout'
  } else if (coef == 3) {
    # Criterion:  adjusted p.val < 0.1 and abs(log2FC) > 1
    relevants <- results[1:20, ] %>% filter(adj.P.Val < 0.1, abs(logFC) > 1) 
    title <- 'Effec of the IL-33 treatment'
    } else {
    # Criterion:  p.value < 0.01 and abs(log2FC) > 0.8
    title <- 'Effect of the knockout-treatment interaction'
    relevants <- results[1:20, ] %>% filter(P.Value < 0.01, abs(logFC) > 0.75) 
  }
  picture <- ggplot(data = results) +
    geom_point(size = 0.1, aes(x = logFC, y = mlog10P, colour = mlog10P >= 1)) +
    geom_text(data = relevants, size = 2, aes(x = logFC, y = mlog10P, label=mgi_symbol), hjust=0, vjust=0) +
    ggtitle(title) + xlab('Log2 Fold Change') + ylab('-Log10(P-value)') +
    labs(colour = 'P-value < 0.1')
  return(picture)
}

fit <- limma_fit

# Plot the most altered genes by the Dusp5 KO
coef <- 2
plot_volcano(fit, coef)

# Plot the most altered genes by the IL-33 treatment
coef <- 3
plot_volcano(fit, coef)

# Plot the most altered genes by the interaction between treatment and genotype
coef <- 4
plot_volcano(fit, coef)

# # The following charts but made with the limma library, to compare. 
# coef <- 2
# # Be noticed that here, because it has not been eliminated, the outstanding 
# # point corresponds to the Dusp5 gene. 
# volcanoplot(fit, coef = coef, style = 'p-value', highlight = 0, 
#             names = fit$genes$ID, hl.col='blue', 
#             xlab = 'Log2 Fold Change', ylab = NULL, pch=16, cex=0.35)
# volcanoplot(fit, coef = coef, style = 'B-statistic', highlight = 0, 
#             names = fit$genes$ID, hl.col='blue', 
#             xlab = 'Log2 Fold Change', ylab = NULL, pch=16, cex=0.35)
# 
# coef <- 3
# volcanoplot(fit, coef = coef, style = 'p-value', highlight = 0, 
#             names = fit$genes$ID, hl.col='blue', 
#             xlab = 'Log2 Fold Change', ylab = NULL, pch=16, cex=0.35)
# volcanoplot(fit, coef = coef, style = 'B-statistic', highlight = 0, 
#             names = fit$genes$ID, hl.col='blue', 
#             xlab = 'Log2 Fold Change', ylab = NULL, pch=16, cex=0.35)
# 
# coef <- 4
# volcanoplot(fit, coef = coef, style = 'p-value', highlight = 0, 
#             names = fit$genes$ID, hl.col='blue', 
#             xlab = 'Log2 Fold Change', ylab = NULL, pch=16, cex=0.35)
# volcanoplot(fit, coef = coef, style = 'B-statistic', highlight = 0, 
#             names = fit$genes$ID, hl.col='blue', 
#             xlab = 'Log2 Fold Change', ylab = NULL, pch=16, cex=0.35)
