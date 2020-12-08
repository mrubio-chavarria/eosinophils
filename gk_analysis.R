
# Load libraries
library(limma)
library(tidyverse)

# Functions
findMid <- function(dtf) {
  n_results <- dim(dtf)[1]
  rownames(dtf) <- seq(1, n_results, 1)
  ratios <- dtf$ratio
  m1 <- max(ratios)
  m2 <- min(ratios)
  ratio_level <- mean(c(m1, m2))
  r1 <- rownames(dtf[dtf$ratio >= ratio_level, ])[1]
  r2 <- tail(rownames(dtf[dtf$ratio < ratio_level, ]), 1)
  height <- mean(as.integer(c(r1, r2)))
  return(height)
}

goana_analysis <- function(entrezgene_ids, ont, n_results, limit_n) {
  goana_results <- goana(entrezgene_ids, species = "Mm", Ont = ont) 
  goana_results <- goana_results[goana_results$N <= limit_n, ] 
  results <- goana_results %>% transmute(p.value = P.DE, term = Term,
                                         term.id = rownames(goana_results),
                                         adj.p.value = p.adjust(P.DE, method = "BH"),
                                         ratio = DE / N,
                                         de = DE)
  results <- results[order(results$ratio), ]
  results <- results[results$adj.p.value <= 0.1, ]
  main_results <- tail(results, n_results)
  height <- findMid(main_results)
  main_results %>%
    ggplot(aes(y = reorder(term, ratio), x = de, fill = ratio)) + 
    labs(fill = "Annotated term percentage", y = "GO term description",
         x = "Annotated genes of the set") +
    geom_bar(stat = "identity") + 
    geom_hline(yintercept = height , colour = "red", linetype = "dashed") + 
    ggtitle("GOANA: Biological Processes")
}

kegga_analysis <- function(entrezgene_ids, n_results = 75, limit_n = 200) {
  analysis_results <- kegga(entrezgene_ids, species = "Mm")
  analysis_results <- analysis_results[analysis_results$N <= limit_n, ] 
  results <- analysis_results %>% transmute(p.value = P.DE,
                                            term = Pathway,
                                            term.id = rownames(analysis_results),
                                            adj.p.value = p.adjust(P.DE, method = "BH"),
                                            ratio = DE / N,
                                            de = DE)
  results <- results[order(results$ratio), ]
  results <- results[results$adj.p.value <= 0.1, ]
  main_results <- tail(results, n_results)
  height <- findMid(main_results)
  main_results %>%
    ggplot(aes(y = reorder(term, ratio), x = de, fill = ratio)) + 
    labs(fill = "Annotated term percentage", y = "GO term description",
         x = "Annotated genes of the set") +
    geom_bar(stat = "identity") + 
    geom_hline(yintercept = height , colour = "red", linetype = "dashed") + 
    ggtitle("KEGGA: Pathways")
}

# Load data
setwd('/home/mario/Projects/eosinophils')
filename <- "selected_genes.txt"
ann.dataset <- read.table(filename, header = TRUE)

# Process data
significant <- ann.dataset[ann.dataset$adj.P.Val < 0.1, ]$entrezgene_id
significant <- significant[!(is.na(significant) | significant == "")]


# GOANA
# Biological process
entrezgene_ids <- significant # IDs to perform the test
n_results <- 75 # Show the n_results most relevant terms
ont <- "BP" # Ontology to show
limit_n <- 200 # Upper limit in annotated genes to select term
goana_analysis(entrezgene_ids, ont, n_results, limit_n)

# Cellular compartment
entrezgene_ids <- significant
n_results <- 75
ont <- "CC"
limit_n <- 200
goana_analysis(entrezgene_ids, ont, n_results, limit_n)

# Molecular function
entrezgene_ids <- significant
n_results <- 75
ont <- "MF"
limit_n <- 200
goana_analysis(entrezgene_ids, ont, n_results, limit_n)

# KEGGA
entrezgene_ids <- significant
n_results <- 75
limit_n <- 200
kegga_analysis(entrezgene_ids, n_results, limit_n)

