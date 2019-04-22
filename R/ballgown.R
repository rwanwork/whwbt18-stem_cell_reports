#!/usr/bin/env Rscript
#####################################################################
##
##  Run ballgown
##
#####################################################################

##############################
##  Read in libraries and modules

#  Load in libraries
library (Cairo)  #  For transparency
library (ballgown)
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

FPKM_THRESHOLD <- 5
FOLD_CHANGE_CUTOFF <- 5

GENE_HEATMAP_HEIGHT <- 3.5  # 2.5
GENE_HEATMAP_WIDTH <- 1.6

SCATTERPLOT_FIGURE_WIDTH <- 15
SCATTERPLOT_FIGURE_HEIGHT <- 15


#####################################################################
##  Functions
#####################################################################


#####################################################################
##  Process arguments using docopt
#####################################################################

"Usage: ballgown.R --sample SAMPLE --results RESULTS --graphs GRAPHS [--log]

Options:
 --sample SAMPLE  The sample to use as input.
 --results RESULTS  The results directory.
 --graphs GRAPHS  The graphs directory.
 --log  Apply log base 10 to the FPKM values (i.e., do not apply log by default!)
.
" -> options

# Retrieve the command-line arguments
opts <- docopt (options)

LOG_ARG <- NULL

# Update values with the program arguments
SAMPLE_ARG <- opts$sample
RESULTS_ARG <- opts$results
GRAPHS_ARG <- opts$graphs
if (!is.null (opts$log)) {
  LOG_ARG <- opts$log
}


#####################################################################
##  Variables for testing
#####################################################################

# SAMPLE_ARG <- "combo1"
# RESULTS_ARG <- "./"
# GRAPHS_ARG <- "./"


#####################################################################
##  Set up variables based on parameters
#####################################################################

writeLines (paste ("Sample combination:  ", SAMPLE_ARG, sep=""))
writeLines (paste ("Perform log?:  ", LOG_ARG, sep=""))

INPUT_DIR <- paste ("Data/", SAMPLE_ARG, "/", sep="")
PHENODATA_FN <- paste ("Data/combos/", SAMPLE_ARG, ".csv", sep="")
if (!file.exists (PHENODATA_FN)) {
  stop (paste ("The file ", PHENODATA_FN, "does not exist!", sep=""))
}


#####################################################################
##  Read in the data
#####################################################################

pheno_data <- read.csv (PHENODATA_FN)

order_correct <- all (pheno_data[["ids"]] == list.files (INPUT_DIR))
writeLines (paste ("Order of files (must be TRUE):  ", order_correct, sep=""))
if (isFALSE (order_correct)) {
  print (list.files (INPUT_DIR))
  stop (paste ("Order of files incorrect!  Check the above with ", PHENODATA_FN, "!", sep=""))
}
writeLines ("\n\n\n")

ballgown_obj = ballgown (dataDir = INPUT_DIR, samplePattern = "*", pData = pheno_data)
ballgown_obj


#####################################################################
##  Perform a subset of the data, and clean up the experiment names
#####################################################################

genes_all <- gexpr (ballgown_obj)
genes_all <- cleanUpBasic (genes_all)
dim (genes_all)
out_fn <- paste (RESULTS_ARG, "/", "genes-all.tsv", sep="")
writeTable (genes_all, out_fn)

#  Select rows where the row mean is greater than 1 (i.e., remove genes where expression levels are low throughout all data sets)
genes_mean <- genes_all[rowMeans (genes_all[,]) > 1,]
dim (genes_mean)
out_fn <- paste (RESULTS_ARG, "/", "genes-mean.tsv", sep="")
writeTable (genes_mean, out_fn)

#  Select rows where at least one cell value is greater than the threshold
genes_fpkm <- genes_mean[apply (genes_mean[, 1:ncol(genes_mean)], 1, function( x ) any( x > FPKM_THRESHOLD ) ), ]
dim (genes_fpkm)
out_fn <- paste (RESULTS_ARG, "/", "genes-fpkm.tsv", sep="")
writeTable (genes_fpkm, out_fn)


#####################################################################
##  Create a heatmap using a subset of genes
#####################################################################

head (genes_all)

#  Intersect gene list
intersect_genes <- intersect (rownames (genes_all), HEATMAP_GENES)

#  Get a subset of genes_all
x <- genes_all[intersect_genes,]

#  Order the subset by the order of the HEATMAP_GENES
x <- x[match (HEATMAP_GENES, rownames (x)),]

#  Remove records with NA
x <- x[complete.cases (x), ]

#  Reorder the columns
# x <- RearrangeSamples (SAMPLE_ARG, x)

#  Take the log (base 10) of the values to give a better range for the heat map
x <- log10 (x + 1)

#  Set the locations of each part of the heatmap
#    https://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
mylmat = rbind (c (0, 4), c (2, 1), c (0, 3))
mylhei = c (1, 4, 0.1)
mylwid = c (0.1, 4)

out_fn <- paste (GRAPHS_ARG, "/", "genes-heatmap.eps", sep="")
postscript (file=out_fn, onefile=FALSE, width=GENE_HEATMAP_WIDTH, height=GENE_HEATMAP_HEIGHT, paper="special", horizontal=FALSE)
# heatmap.2 (x, srtCol=45, key=FALSE, dendrogram="none", col=bluered, trace="none", Colv=FALSE, Rowv=FALSE, lwid=c(0.1,4), lhei=c(0.1,4), cexRow=1.25, cexCol=1.25, rowsep=c(16), colsep=c(6), sepcolor="black")
heatmap.2 (x, 
  srtCol=90, 
  dendrogram="none", 
  col=greenred, 
  trace="none", 
  Colv=FALSE, 
  Rowv=FALSE, 
  lmat=mylmat,
  lwid=mylwid,
  lhei=mylhei,
  key=TRUE,
  key.title=NA,
  key.xlab=NA,
  keysize=0.25,
  key.par = list (cex=0.5, mar=c(2, 0.25, 0.5, 2.5)),
  density.info="none",
#   rowsep=c(16), 
#   colsep=c(6), 
#   sepcolor="yellow", 
  margins=c(3, 2.5), # bottom and right
  offsetRow=-0.25, 
  offsetCol=-0.25,
  cexRow=0.75,
  cexCol=0.75
)
dev.off ()


