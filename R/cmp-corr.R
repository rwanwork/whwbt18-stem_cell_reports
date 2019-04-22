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
library (genefilter)
library (data.table)
library (ggplot2)
library (docopt)  #  See:  https://github.com/docopt/docopt.R

#  Load additional R functions
# source ("../modules/general.R")
# source ("../modules/multiplot.R")
# source ("../modules/clustering.R")
# source ("../modules/pca.R")


#####################################################################
##  Constants
#####################################################################


#####################################################################
##  Functions
#####################################################################


#####################################################################
##  Process arguments using docopt
#####################################################################

"Usage: cmp-corr.R [--log]

Options:
 --log  Apply log transformation.
.
" -> options

# Retrieve the command-line arguments
opts <- docopt (options)

# Update values with the program arguments
LOG_ARG <- opts$log


#####################################################################
##  Variables for testing
#####################################################################



#####################################################################
##  Set up variables based on parameters
#####################################################################

if (LOG_ARG) {
  LOG_PATH <- "with-log"
  OUTPUT_FN <- "with-log.eps"
  SUBTITLE <- "(with log)"
} else if (!LOG_ARG) {
  LOG_PATH <- "without-log"
  OUTPUT_FN <- "without-log.eps"
  SUBTITLE <- "(without log)"
}

INPUT_FN <- paste ("scatterplots/", LOG_PATH, "/correlations.tsv", sep="")

writeLines (paste ("Log transformation:  ", LOG_ARG, sep=""))
writeLines (paste ("Output filename:  ", OUTPUT_FN, sep=""))
writeLines (paste ("Input filename:  ", INPUT_FN, sep=""))


#####################################################################
##  Read in the data
#####################################################################

x <- read.table (INPUT_FN, sep="\t", header=FALSE)
names (x) <- c("Sample1", "Sample2", "Pearson", "Spearman")

summary (x)


#####################################################################
##  Create scatter plots
#####################################################################

my_title <- paste ("Pearson correlation vs Spearman rank correlation")
my_subtitle <- SUBTITLE

plot_obj <- ggplot (x, aes (y=Pearson, x=Spearman)) + geom_point (shape=1)      # Use hollow circles
plot_obj <- plot_obj + geom_smooth (method=lm)   # Add linear regression line (by default includes 95% confidence region)
plot_obj <- plot_obj + ylab ("Pearson") + xlab ("Spearman") 
#       plot_obj <- plot_obj + ggtitle ()
plot_obj <- plot_obj + geom_abline (slope=1, intercept=0, colour="red")
plot_obj <- plot_obj + labs (title=my_title, subtitle=my_subtitle)

#  Spearman and Pearson correlations are all between [0, 1]; none are between [-1, 0]...not sure why
plot_obj <- plot_obj + scale_x_continuous (limits = c(0, 1))
plot_obj <- plot_obj + scale_y_continuous (limits = c(0, 1))
      
ggsave (OUTPUT_FN, device=cairo_ps)
      
