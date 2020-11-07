library(Biobase)
library(oligo)
library(data.table)
library(arrayQualityMetrics)
library(biomaRt)

setwd('/home/aetius/Projects/eosinophils/GSE62999_RAW/')  # set your working directory
cel_files = list.files(path = getwd(), pattern = '*.CEL.gz', full.names = TRUE)
 
# Parse the CEL files and run QC
# T6his function reads in microarray raw data (CEL files) and creates a container
# object that stores them plus metadata on the microarray platform, species etc.
parsed_cels = oligo::read.celfiles(cel_files, verbose = TRUE)  
# This function runs the main normalisation pipeline (RMA), it performs multiple
# steps to normalise the data, you will read more on it in a paper.
parsed_cels_rma = oligo::rma(parsed_cels, normalize = TRUE, background = TRUE)  

# QC on RMA-normalised files, will help you check if there are outlier microarrays
# This will write a report into your disk (in the working directory which you
# can check with getwd(), open index.html with your browser)
# qc_data = arrayQualityMetrics(parsed_cels_rma)  

# The data is still in the format of (probes x samples), we have to aggregate
# the probes (which have format nnnnn_at) into mouse genes
# However we first need to get an ID mapping table that will tell us which probes
# belong to which genes in the Ensembl Gene ID system, we will do this using
# biomaRt, a very helpful package for bioinformatics.
mart = useEnsembl(biomart='ensembl', dataset='mmusculus_gene_ensembl')  # set up connection to the Ensembl Mouse database
mouse_probes = row.names(parsed_cels_rma@assayData$exprs)
id_translation_table = getBM(attributes = c('affy_mouse430_2', 'ensembl_gene_id'),
                             filters = 'affy_mouse430_2',
                             values = mouse_probes,
                             mart=mart)  # might take a while

# Now add the Ensembl gene IDs
expression_data = parsed_cels_rma@assayData$exprs  # the intensity data is log2-normalised by the 'rma' function already
expression_data = as.data.frame(expression_data)
expression_data$Ensembl = id_translation_table$ensembl_gene_id[match(row.names(expression_data), id_translation_table$affy_mouse430_2)]

# Turn into a faster data.table object
expression_data_dt = as.data.table(expression_data)
# Some genes were not mapped from probes to Ensembl IDs, drop them
expression_data_dt = expression_data_dt[!is.na(expression_data_dt$Ensembl),]
 
# Aggregate id-wise and use mean (multiple probes can map to the same gene!)
sample_columns = row.names(parsed_cels_rma@phenoData@data)
expression_data_dt_agg = expression_data_dt[,lapply(.SD, mean), by=Ensembl, .SDcols=sample_columns]  # this is a bit tricky function call if you don't know the package, I can explain
 
# Sanity check, do global intensity distributions look similar?
boxplot(expression_data_dt_agg[,2:21])
 
# The sample names are a bit ugly, let's format them a bit too
colnames(expression_data_dt_agg) = gsub('_Mouse430v2.CEL.gz', '', colnames(expression_data_dt_agg))  # this command changes the ending into 'nothing' in all the column names

# Round the values a bit
expression_data_dt_agg[,2:21] = round(expression_data_dt_agg[,2:21], 3)

# We're done! Write the data file to disk
write.csv(expression_data_dt_agg, './GSM1537865_RMA_processed_log2_Ensembl_05112020.csv', quote=F, row.names = F)  
# I like to keep the date of data export in the file name and sometimes a database ID (GSM...) to be able to quickly
# find sources in few years, but feel free to come up with your own system
