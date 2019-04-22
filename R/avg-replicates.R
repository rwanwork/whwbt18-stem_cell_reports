#!/usr/bin/env Rscript
#####################################################################
##
##  Perform an average across replicates
##
#####################################################################

##############################
##  Read in libraries and modules

#  Load in libraries
library (org.Hs.eg.db)
library (Cairo)  #  For transparency
library (docopt)  #  See:  https://github.com/docopt/docopt.R
library (ggplot2)

#  Load additional R functions
source ("R/modules/general.R")
source ("R/modules/clean-labels.R")


#####################################################################
##  Constants
#####################################################################

#  Value added before taking the log
LOG_CONSTANT <- 1

#  Any value above 1 would be fine
PVALUE_ERROR <- 2 

#  The maximum number of genes accepted by the DAVID API
DAVID_GENES_LIMIT <- 3000

#  Threshold for the p-value
PVALUE_THRESHOLD <- 0.01

#  Since log base 2 is used, this implies a fold change of 2
LOG_FOLDCHG_CUTOFF <- 1


#####################################################################
##  Functions
#####################################################################

ApplyTTest <- function (dataset, sample1_data, sample2_data) {
  for(i in 1:nrow (sample1_data)) {
    data_a = unlist (sample1_data[i,])
    data_b = unlist (sample2_data[i,])
    
    dataset[i, "pvalue_orig"] <- t.test (data_a, data_b)$p.value
  }
  
  #  Adjust the p-values using FDR
  dataset[["pvalue"]] <- p.adjust (dataset[["pvalue_orig"]], method="fdr")
  
  #  Change p-values of NA to 2 (i.e., a large, non-sensical number)
  dataset <- within (dataset, pvalue[is.na (pvalue)] <- PVALUE_ERROR)
  
  #  Calculate the log of the p-value
  dataset[["logpvalue"]] <- -1 * log10 (dataset[["pvalue"]])

  return (dataset)
}


AugmentData = function (dataset, sample1, sample2) {
  x_mean <- paste (SAMPLE1_ARG, "_mean", sep="")
  y_mean <- paste (SAMPLE2_ARG, "_mean", sep="")

  x_log <- paste (SAMPLE1_ARG, "_log", sep="")
  y_log <- paste (SAMPLE2_ARG, "_log", sep="")

  #  Calculate log of the FPKMs
  dataset[x_log] <- log2 (dataset[x_mean] + LOG_CONSTANT)
  dataset[y_log] <- log2 (dataset[y_mean] + LOG_CONSTANT)
  
  #  Calculate values for the MA plot
  dataset["M"] <- dataset[x_log] - dataset[y_log]
  dataset["A"] <- 0.5 * (dataset[x_log] + dataset[y_log])

  return (dataset)
}


#####################################################################
##  Process arguments using docopt
#####################################################################

"Usage:  avg-replicates.R --input INPUT --allgenes ALLGENES --siggenes SIGGENES --sample1 SAMPLE1 --sample2 SAMPLE2

Options:
 --input INPUT  The input file
 --allgenes ALLGENES  The output file of all genes
 --siggenes SIGGENES  The output file of just the significant genes
 --sample1 SAMPLE1  Sample #1
 --sample2 SAMPLE2  Sample #2
.
" -> options

# Retrieve the command-line arguments
opts <- docopt (options)

# Update values with the program arguments
INPUT_ARG <- opts$input
ALLGENES_ARG <- opts$allgenes
SIGGENES_ARG <- opts$siggenes
SAMPLE1_ARG <- opts$sample1
SAMPLE2_ARG <- opts$sample2


#####################################################################
##  Variables for testing
#####################################################################

# INPUT_ARG <- "Results/combo3/ballgown/genes.tsv"
# OUTPUT_ARG <- "./genes.tsv"
# SAMPLE1_ARG <- "MPC"
# SAMPLE2_ARG <- "ASC"


#####################################################################
##  Read in the data
#####################################################################

in_fn <- INPUT_ARG
genes_all <- read.table (in_fn, sep="\t", header=TRUE)

colnames (genes_all)[1] <- "Gene"

head (genes_all)

dim (genes_all)

#  Work on the data frame
genes_dataframe <- data.frame (genes_all)


#####################################################################
##  Perform a subset of the data so that only the columns we're interested in are included
#####################################################################

#  Take the union of the two samples
sample1_arg_search <- paste (SAMPLE1_ARG, "_\\d", sep="")
sample2_arg_search <- paste (SAMPLE2_ARG, "_\\d", sep="")
x <- grepl (sample1_arg_search, names (genes_dataframe))
y <- grepl (sample2_arg_search, names (genes_dataframe))
z <- x | y

#  Also include the gene names
z[1] <- TRUE

#  Perform a subset
genes_subset <- genes_dataframe[z]


#####################################################################
##  Extract the two samples into separate tables
#####################################################################

x <- genes_subset[grepl (sample1_arg_search, names (genes_subset))]
y <- genes_subset[grepl (sample2_arg_search, names (genes_subset))]

x_len <- length (x)
if (x_len == 0) {
  stop (paste ("Number of samples for ", SAMPLE1_ARG, " is 0!\n", sep=""))
}
y_len <- length (y)
if (y_len == 0) {
  stop (paste ("Number of samples for ", SAMPLE2_ARG, " is 0!\n", sep=""))
}

print (paste ("Number of samples for", SAMPLE1_ARG, ":", x_len))
print (colnames (x))

print (paste ("Number of samples for", SAMPLE2_ARG, ":", y_len))
print (colnames (y))


#####################################################################
##  Calculate the mean across the two samples' replicates
#####################################################################

x_mean <- paste (SAMPLE1_ARG, "_mean", sep="")
genes_subset[[x_mean]] <- rowMeans (x, na.rm=TRUE)

y_mean <- paste (SAMPLE2_ARG, "_mean", sep="")
genes_subset[[y_mean]] <- rowMeans (y, na.rm=TRUE)


#####################################################################
##  Calculate the p-value on a row-by-row basis
#####################################################################

#  Add p-values to the data frame by using a t-test
genes_subset <- ApplyTTest (genes_subset, x, y)

#  Add the M and A values to the data frame
genes_subset <- AugmentData (genes_subset, SAMPLE1_ARG, SAMPLE2_ARG)


#####################################################################
##  Identify significant genes based on a combination of p-value and log fold change
#####################################################################

#  Create new columns based on significance and fold change cutoffs
#    We have to create a new factor value "No2" since duplicates are not allowed
genes_subset[["pvalue_sig"]] <- cut (genes_subset[["pvalue"]], breaks = c(-Inf, 0, 1 * PVALUE_THRESHOLD, Inf), labels = c("No", "Yes", "No2"))
genes_subset[["foldchg_flag"]] <- cut (genes_subset[["M"]], breaks = c(-Inf, -1 * LOG_FOLDCHG_CUTOFF, 1 * LOG_FOLDCHG_CUTOFF, Inf), labels = c("Down", "No", "Up"))

#  Identify genes that have significant p-values *and* have an up or down fold change
pvalue_bool = (genes_subset[["pvalue_sig"]] == "Yes")
foldchg_bool_up = (genes_subset[["foldchg_flag"]] == "Up")
foldchg_bool_down = (genes_subset[["foldchg_flag"]] == "Down")
genes_subset[["up_reg_flag"]] = (pvalue_bool & foldchg_bool_up)
genes_subset[["down_reg_flag"]] = (pvalue_bool & foldchg_bool_down)

#  Combine the above two tests.  Create a new column with a default of 
#    "No" and then set "Up" or "Down" afterwards
genes_subset[["combined"]] = "No"
genes_subset <- within (genes_subset, combined[up_reg_flag == TRUE] <- "Up")
genes_subset <- within (genes_subset, combined[down_reg_flag == TRUE] <- "Down")


#####################################################################
##  Output set of significant genes
#####################################################################

#  Make a copy of the above table, taking only genes that are statistically significant
#    (But they can be either up or down regulated.)
genes_subset_2 <- genes_subset[genes_subset[["combined"]] != "No",]

#  Select only the significant genes (whether they are up or down regulated)
significant_genes <- genes_subset_2[,c("Gene", "M")]
print (paste ("II\tDimensions of table of signficant genes:  ", toString (dim (significant_genes)), sep=""))

#  Check if we have too many up or down-regulated genes
significant_up_genes <- significant_genes[significant_genes[["M"]] >= 0,]
significant_down_genes <- significant_genes[significant_genes[["M"]] < 0,]

print (paste ("II\tDimensions of table of up-regulated signficant genes:  ", toString (dim (significant_up_genes)), sep=""))
print (paste ("II\tDimensions of table of down-regulated signficant genes:  ", toString (dim (significant_down_genes)), sep=""))

if (dim (significant_up_genes)[1] > DAVID_GENES_LIMIT) {
  stop (paste ("Number of up-regulated genes (", dim (significant_up_genes)[1], ") exceeded ", DAVID_GENES_LIMIT, " for ", SAMPLE1_ARG, " and ", SAMPLE2_ARG, ".", sep=""))
}

if (dim (significant_down_genes)[1] > DAVID_GENES_LIMIT) {
  stop (paste ("Number of down-regulated genes (", dim (significant_down_genes)[1], ") exceeded ", DAVID_GENES_LIMIT, " for ", SAMPLE1_ARG, " and ", SAMPLE2_ARG, ".", sep=""))
}

#  Give column names to the newly generated table
names (significant_genes) <- c ("SYMBOL", "M")
  
#  Use "as.character ()" to prevent the warning "'keys' must be a character vector"
gene_mapping <- select (org.Hs.eg.db, as.character (significant_genes[["SYMBOL"]]), "ENTREZID", "SYMBOL")

merged_table <- merge (significant_genes, gene_mapping, by="SYMBOL")
  
merged_table_na <- merged_table[is.na (merged_table[["ENTREZID"]]),]
print (paste ("Set of genes which did not have a matching Entrez ID:", sep=""))
print (merged_table_na)


#####################################################################
##  Export the entire gene table (genes_subset) and the significant 
##    genes to separate files
#####################################################################

out_table <- as.matrix (genes_subset)
out_fn <- ALLGENES_ARG
writeTable (out_table, out_fn)

out_fn <- SIGGENES_ARG
writeTable (merged_table, out_fn)


#####################################################################
##  For debugging purposes, show the head of the table
#####################################################################

head (genes_subset)

