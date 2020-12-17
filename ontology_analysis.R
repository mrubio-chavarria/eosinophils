 
# Load libraries
library(topGO)
library(DataCombine)
library(biomaRt)
library(tidyverse)


# Functions
compareMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

topDiffGenes <- function (allScore) {
  # allScore is a vector of p-values
  confidence <- 0.90
  limit = 1 - confidence
  return(allScore < limit)
}

# Read the list of genes
setwd('/home/mario/Projects/eosinophils')
filename <- "selected_genes.txt"
genes <- read.table(filename, header = TRUE)
genes$p_value <- NULL

# Obtain GO terms
# bm <- useMart("ensembl")
# bm <- useDataset("mmusculus_gene_ensembl", mart=bm)
gene2GO <- getBM(mart=bm, attributes=c('ensembl_gene_id', 'go_id', 'external_gene_name'))
genes <- merge(x = genes, y = gene2GO, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
genes <- genes[order(genes$adj.P.Val), ]
# Be noticed that some genes are to be discarded because they do not have GO term
genes_go <- genes[!(is.na(genes$go_id) | genes$go_id==""), ]
print(paste("Entries lost: ", nrow(genes) - nrow(genes_go)))
genes.names <- genes_go$ensembl_gene_id
gene2GO <- sapply(genes_go[, 4], c)
names(gene2GO) <- genes.names

# Adapt data to the format 
genes.names <- genes$ensembl_gene_id
genes.ids <- as.vector(genes$adj.P.Val)
names(genes.ids) <- genes.names
gene2GO <- annFUN.gene2GO(whichOnto = "MF", gene2GO = gene2GO)

sampleGOdata <- new("topGOdata",
                    description = "Ensembl GO enrichment", ontology = "BP",
                    allGenes = genes.ids, geneSel = topDiffGenes,
                    annot = annFUN.gene2GO, gene2GO=gene2GO, nodeSize = 10)

# Fisher test
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

# Kolmogorov-Smirnov test
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

# Show most significant genes
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 19)

# Prepare a display of the results
names(allRes)<-str_replace_all(names(allRes), c(" " = "_"))
allRes <- allRes[order(allRes$Rank_in_classicFisher), ]
allRes %>%
ggplot(aes(y = reorder(Term, -Rank_in_classicFisher), x = Annotated, fill = Rank_in_classicFisher)) + 
  labs(fill = "Rank (Fisher)", y = "Term", x = "Annotated") +
  geom_bar(stat = "identity")

