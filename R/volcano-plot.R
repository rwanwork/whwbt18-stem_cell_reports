#!/usr/bin/env Rscript
#####################################################################
##
##  Create volcano plot after averaging across replicates
##
#####################################################################

##############################
##  Read in libraries and modules

#  Load in libraries
library (Cairo)  #  For transparency
library (docopt)  #  See:  https://github.com/docopt/docopt.R
library (ggplot2)

#  Load additional R functions
source ("R/modules/general.R")
source ("R/modules/multiplot.R")


#####################################################################
##  Constants
#####################################################################

FIGURE_WIDTH <- 20
FIGURE_HEIGHT <- 10
FIGURE_DPI <- 600

#  Data levels ("No" should be first)
DOT_LEVELS <- c("No", "Up", "Down")

#  Colours for each of the dots
DOT_COLOURS <- c("Up" = "red", "Down" = "green", "No" = "gray")


#####################################################################
##  Functions
#####################################################################


#####################################################################
##  Process arguments using docopt
#####################################################################

"Usage:  volcano-plot.R --input INPUT --output OUTPUT --sample1 SAMPLE1 --sample2 SAMPLE2

Options:
 --input INPUT  The input genes file
 --output OUTPUT  The output graph file
 --sample1 SAMPLE1  Sample #1
 --sample2 SAMPLE2  Sample #2

.
" -> options

# Retrieve the command-line arguments
opts <- docopt (options)

# Update values with the program arguments
INPUT_ARG <- opts$input
OUTPUT_ARG <- opts$output
SAMPLE1_ARG <- opts$sample1
SAMPLE2_ARG <- opts$sample2


#####################################################################
##  Variables for testing
#####################################################################

# INPUT_ARG <- "./genes.tsv"
# GENES_ARG<- "./out-genes.tsv"
# GRAPHS_ARG <- "./maplot.eps"
# SAMPLE1_ARG <- "MPC"
# SAMPLE2_ARG <- "ASC"


#####################################################################
##  Read in the data
#####################################################################

in_fn <- INPUT_ARG
genes_table <- read.table (in_fn, sep="\t", header=TRUE)

dim (genes_table)


#####################################################################
##  Create the Volcano plot
#####################################################################

#  Change the significance level into a factor so that we can sort it; this ensures the coloured dots are on top of the gray ones.
genes_table[["combined"]] <- factor (genes_table[["combined"]], levels = DOT_LEVELS)
genes_table <- genes_table[order (genes_table[["combined"]]),]

ggplot_obj <- ggplot (genes_table, aes (y=logpvalue, x=M, color=combined))
ggplot_obj <- ggplot_obj + geom_hline (yintercept=0, colour="gray")
ggplot_obj <- ggplot_obj + geom_vline (xintercept=0, colour="gray")
ggplot_obj <- ggplot_obj + geom_point (shape=20, size=1)
ggplot_obj <- ggplot_obj + ylab (expression ("-" * log[10] ~ " p-value"))
ggplot_obj <- ggplot_obj + xlab (bquote (log[2] ~ "of" ~ .(SAMPLE1_ARG) * "/" * .(SAMPLE2_ARG)))
ggplot_obj <- ggplot_obj + scale_color_manual (values = DOT_COLOURS)
# ggplot_obj <- ggplot_obj + scale_x_continuous (limits = c (-5000, 5000))
ggplot_obj <- ggplot_obj + labs (colour="combined")
ggplot_obj <- ggplot_obj + theme (legend.position="none")  #  Remove legend
  
out_fn <- OUTPUT_ARG
ggsave (out_fn, plot=ggplot_obj, device=cairo_ps, width = FIGURE_WIDTH, height = FIGURE_HEIGHT, dpi = FIGURE_DPI, units = "cm")  


#####################################################################
##  For debugging purposes, show the output of the table
#####################################################################

head (genes_table)


