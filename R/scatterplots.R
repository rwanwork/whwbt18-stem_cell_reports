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
library (dplyr)
library (data.table)
library (devtools)
library (amap)  #  Dist function
library (ggplot2)
library (ggfortify)
library (docopt)  #  See:  https://github.com/docopt/docopt.R

#  Load additional R functions
source ("../modules/general.R")
source ("../modules/multiplot.R")
source ("../modules/clustering.R")
source ("../modules/pca.R")


#####################################################################
##  Constants
#####################################################################

FPKM_THRESHOLD <- 5
FOLD_CHANGE_CUTOFF <- 5

HCLUST_HEIGHT <- 9
HCLUST_WIDTH <- 7
HEATMAP_HEIGHT <- 10
HEATMAP_WIDTH <- 10

SCATTERPLOT_FIGURE_WIDTH <- 15
SCATTERPLOT_FIGURE_HEIGHT <- 15


#####################################################################
##  Functions
#####################################################################

cleanUp = function (gene_matrix) {
  colnames (gene_matrix) <- gsub ("FPKM.", "", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("_s", "", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("asc", "ASC", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("sb", "SB", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("qsc", "QSC", colnames (gene_matrix))

  return (gene_matrix)
}



#####################################################################
##  Process arguments using docopt
#####################################################################

"Usage: scatterplots.R --sample SAMPLE [--log]

Options:
 --sample SAMPLE  The sample to use as input.
 --log  Apply log base 2 to the FPKM values (i.e., do not apply log by default!)
.
" -> options

# Retrieve the command-line arguments
opts <- docopt (options)

# Update values with the program arguments
SAMPLE_ARG <- opts$sample
LOG_ARG <- opts$log


#####################################################################
##  Variables for testing
#####################################################################

# SAMPLE_ARG <- "combo1"


#####################################################################
##  Set up variables based on parameters
#####################################################################

if (SAMPLE_ARG == "combo1") {
  SLURM_ID <- "246858"
  INPUT_PATTERN <- "MPC*"
} else if (SAMPLE_ARG == "combo2") {
  SLURM_ID <- "246859"
  INPUT_PATTERN <- "MPC*"
} else if (SAMPLE_ARG == "combo3") {
  SLURM_ID <- "246860"
  INPUT_PATTERN <- "ASC*|MPC*|QSC*"
} else if (SAMPLE_ARG == "combo4") {
  SLURM_ID <- "246861"
  INPUT_PATTERN <- "ASC*|MPC*|QSC*"
} else if (SAMPLE_ARG == "combo5") {
  SLURM_ID <- "246862"
  INPUT_PATTERN <- "ASC*|MPC*|QSC*"
} else if (SAMPLE_ARG == "combo6") {
  SLURM_ID <- "246863"
  INPUT_PATTERN <- "ASC*|MPC*|QSC*"
} else if (SAMPLE_ARG == "combo7") {
  SLURM_ID <- "246864"
  INPUT_PATTERN <- "*"
} else {
  writeLines ("That sample name is invalid!\n")
}

writeLines (paste ("Sample combination:  ", SAMPLE_ARG, sep=""))
writeLines (paste ("Perform log?:  ", LOG_ARG, sep=""))
writeLines (paste ("Slurm ID:  ", SLURM_ID, sep=""))
writeLines (paste ("Input pattern:  ", INPUT_PATTERN, sep=""))

INPUT_DIR <- paste ("data/", SLURM_ID, "/ballgown/", sep="")
PHENODATA_FN <- paste (SAMPLE_ARG, ".csv", sep="")


#####################################################################
##  Read in the data
#####################################################################

pheno_data <- read.csv (PHENODATA_FN)

order_correct <- all (pheno_data$ids == list.files (INPUT_DIR))
writeLines (paste ("Order of files (must be TRUE):  ", order_correct, sep=""))
if (isFALSE (order_correct)) {
  print (list.files (INPUT_DIR))
  stop (paste ("Order of files incorrect!  Check the above with ", PHENODATA_FN, "!", sep=""))
}
writeLines ("\n\n\n")

ballgown_obj = ballgown (dataDir = INPUT_DIR, samplePattern = INPUT_PATTERN, pData = pheno_data)
ballgown_obj


#####################################################################
##  Perform a subset of the data, and clean up the experiment names
#####################################################################

genes_all <- gexpr (ballgown_obj)
dim (genes_all)
writeTable (genes_all, "genes-all.tsv")

genes_all <- cleanUp (genes_all)

genes_mean <- genes_all[rowMeans (genes_all[,]) > 1,]
dim (genes_mean)
writeTable (genes_mean, "genes-mean.tsv")

#  Select rows where at least one cell value is greater than the threshold
genes_fpkm <- genes_mean[apply (genes_mean[, 1:ncol(genes_mean)], 1, function( x ) any( x > FPKM_THRESHOLD ) ), ]
dim (genes_fpkm)
writeTable (genes_fpkm, "genes-fpkm.tsv")


#####################################################################
##  Perform log transformation
#####################################################################

if (LOG_ARG) {
  genes_fpkm <- genes_fpkm + 1
  genes_fpkm <- log2 (genes_fpkm)
}

head (genes_fpkm)


#####################################################################
##  Create scatter plots
#####################################################################

x <- data.frame (genes_fpkm)

for (i in 1:ncol (x)) {
  for (j in 1:ncol (x)) {
    if (i > j) {
      sample1_name <- names (x[1,])[i]
      sample2_name <- names (x[1,])[j]
    
      sample1 <- x[,i]
      sample2 <- x[,j]
    
      pearson_value <- cor (sample1, sample2, method="pearson")
      spearman_value <- cor (sample1, sample2, method="spearman")
      
      fn_eps <- paste (sample1_name, "-", sample2_name, ".eps", sep="")
      
      my_title <- paste (sample1_name, " vs ", sample2_name, sep="")
      my_subtitle <- paste ("Pearson = ", pearson_value, "; Spearman = ", spearman_value, sep="")

      plot_obj <- qplot (sample1, sample2)
      
      plot_obj <- plot_obj + geom_smooth (method=lm)   # Add linear regression line (by default includes 95% confidence region)
      plot_obj <- plot_obj + ylab (sample1_name) + xlab (sample2_name) 
      plot_obj <- plot_obj + geom_abline (slope=1, intercept=0, colour="red")
      plot_obj <- plot_obj + labs (title=my_title, subtitle=my_subtitle)
      
      ggsave (fn_eps, device="eps")
      
      writeLines (paste (names (x[1,])[i], "\t", names (x[1,])[j], "\t", pearson_value, "\t", spearman_value, "\n"), sep="")
    }
  }
}


