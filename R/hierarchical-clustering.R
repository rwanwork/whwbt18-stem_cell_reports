#!/usr/bin/env Rscript
#####################################################################
##
##  Run hierarchical clustering
##
#####################################################################

##############################
##  Read in libraries and modules

#  Load in libraries
library (Cairo)  #  For transparency
library (genefilter)

#  If we load plyr before dplyr (as indicated by the warning), then summary.R script doesn't work.  So, we should load dplyr first (despite the warning)
library (dplyr)
library (data.table)
library (devtools)
library (amap)  #  Dist function
library (ggplot2)
library (ggfortify)
library (gplots)  #  heatmap.2 function
library (docopt)  #  See:  https://github.com/docopt/docopt.R

#  Load additional R functions
source ("R/modules/general.R")
source ("R/modules/multiplot.R")
source ("R/modules/clustering.R")
source ("R/modules/clean-labels.R")
source ("R/modules/about-data.R")


#####################################################################
##  Constants
#####################################################################

HCLUST_HEIGHT <- 9
HCLUST_WIDTH <- 7

HEATMAP_HEIGHT <- 14
HEATMAP_WIDTH <- 14


#####################################################################
##  Functions
#####################################################################


#####################################################################
##  Process arguments using docopt
#####################################################################

"Usage:  hc.R --input INPUT --graphs GRAPHS [--log]

Options:
 --input INPUT  The input directory.
 --graphs GRAPHS  The graphs directory.
 --log  Apply log base 10 to the FPKM values (i.e., do not apply log by default!)
.
" -> options

# Retrieve the command-line arguments
opts <- docopt (options)

LOG_ARG <- NULL

# Update values with the program arguments
INPUT_ARG <- opts$input
GRAPHS_ARG <- opts$graphs
if (!is.null (opts$log)) {
  LOG_ARG <- opts$log
}


#####################################################################
##  Variables for testing
#####################################################################

# INPUT_ARG <- "Results/combo3/ballgown/genes.tsv"
# GRAPHS_ARG <- "./"


#####################################################################
##  Read in the data
#####################################################################

in_fn <- INPUT_ARG
genes_fpkm <- read.table (in_fn, sep="\t", header=TRUE)

#  Change the first column's name and then remove it
colnames (genes_fpkm)[1] <- "Gene"
genes_fpkm[["Gene"]] <- NULL

#  Convert the data frame to a numerical matrix
genes_fpkm <- (data.matrix (genes_fpkm))

#  We do not export the matrix with the longer labels since it is more convenient to refer to them with shorter labels
genes_fpkm <- shortToLongLabels (genes_fpkm)

head (genes_fpkm)

dim (genes_fpkm)


#####################################################################
##  Perform log transformation
#####################################################################

if (!is.null (LOG_ARG)) {
  genes_fpkm <- genes_fpkm + 1
  genes_fpkm <- log10 (genes_fpkm)
}

head (genes_fpkm)


#####################################################################
##  Plot hierarchical clusterings and heatmaps on all samples
#####################################################################

#  Pearson correlation
plotHC (GRAPHS_ARG, genes_fpkm, "pearson", "Pearson correlation")

#  Spearman rank correlation
plotHC (GRAPHS_ARG, genes_fpkm, "spearman", "Spearman rank correlation")

#  Euclidean distance
plotHC (GRAPHS_ARG, genes_fpkm, "euclidean", "Euclidean distance")


