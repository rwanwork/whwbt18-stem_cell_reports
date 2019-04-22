#!/usr/bin/env Rscript
#####################################################################
##
##  Plot funcchart graphs using existing data
##
#####################################################################

##############################
##  Read in libraries and modules

library (ggplot2)
library (Cairo)  #  For transparency
library (docopt)  #  See:  https://github.com/docopt/docopt.R

#  Load additional R functions
source ("R/modules/multiplot.R")


#####################################################################
##  Constants
#####################################################################

FIGURE_WIDTH <- 14
FIGURE_HEIGHT <- 3
FIGURE_DPI <- 600

NUM_CATEGORIES <- 15

MAX_X_AXIS <- 20


#####################################################################
##  Functions
#####################################################################

FirstUpper <- function (x) {
  substr (x, 1, 1) <- toupper (substr (x, 1, 1))
  return (x)
}

ReadData = function (funcchart_arg) {
  func_chart <- data.frame ()

  #  Read in the data
  in_fn <- paste (funcchart_arg, sep="")
  if (!file.exists (in_fn)) {
    return (func_chart)
  }
  
  message (paste ("Processing:  ", in_fn, "...", sep=""))
  
  #  Disable quoting of characters since some pathway names have a quote (i.e., Huntington's)
  func_chart <- read.table (in_fn, sep="\t", header=TRUE, quote = "")

  #  Remove the GO or KEGG identifier to save some space
  func_chart[["Term"]] <- gsub ("^GO:\\d+~", "", func_chart[["Term"]])
  func_chart[["Term"]] <- gsub ("^mmu\\d+:", "", func_chart[["Term"]])
  func_chart[["Term"]] <- gsub ("^hsa\\d+:", "", func_chart[["Term"]])

  #  Make the first letter in upper case
  func_chart[["Term"]] <- FirstUpper (func_chart[["Term"]])

  #  Trim any trailing whitespace
  func_chart[["Term"]] <- trimws (func_chart[["Term"]])

  #  Minor tweaks to the capitalization
  func_chart[["Term"]] <- gsub ("^MRNA", "mRNA", func_chart[["Term"]])
  
  #  Calculate the log p-value
  func_chart["logPValue"] <- -1 * log10 (func_chart["PValue"])
  
  return (func_chart)
}


RoundUpEven = function (num) {
  result <- 0
  
  tmp <- floor (num)
  
  if (tmp %% 2 == 0) {
    result <- tmp + 2
  } else {
    result <- tmp + 1
  }
  
  return (result)
}


#  Print out the functional chart as a bar plot
SaveFigure = function (outprefix_arg, up_dataset, down_dataset, db, figure_fn) {
  if (!is.null (up_dataset)) {
    ##  Up-regulated genes
  
    #  Round up to the nearest even integer
    max_x_axis_up <- RoundUpEven (max (up_dataset[["logPValue"]]))
  
    up_ggplot_obj <- ggplot (up_dataset, aes (x=reorder (Term, logPValue), y=-logPValue))
    up_ggplot_obj <- up_ggplot_obj + geom_bar (position = position_dodge (), stat="identity", fill = "red")

    #  Flip the x and y axes
    up_ggplot_obj <- up_ggplot_obj + coord_flip () 

    #  Remove the legend
    up_ggplot_obj <- up_ggplot_obj + theme (legend.position="none") 

    #  Add axes labels to graph
    up_ggplot_obj <- up_ggplot_obj + ylab (expression ("-" * log[10] ~ " p-value"))  
    up_ggplot_obj <- up_ggplot_obj + xlab ("")
    up_ggplot_obj <- up_ggplot_obj + scale_x_discrete (position = "bottom")
    up_ggplot_obj <- up_ggplot_obj + scale_y_continuous (limits=c(-max_x_axis_up, 0),
      breaks = seq (0, -max_x_axis_up, by = -2),  # y axis values (before coord_flip) 
      labels = seq (0,  max_x_axis_up, by =  2))  # show non-negative values
  }
                     
  ##  Down-regulated genes
  
  #  Round up to the nearest even integer
  max_x_axis_down <- RoundUpEven (max (down_dataset[["logPValue"]]))
  
  down_ggplot_obj <- ggplot (down_dataset, aes (x=reorder (Term, logPValue), y=logPValue))
  down_ggplot_obj <- down_ggplot_obj + geom_bar (position = position_dodge(), stat="identity", fill = "green")

  #  Flip the x and y axes
  down_ggplot_obj <- down_ggplot_obj + coord_flip () 

  #  Remove the legend
  down_ggplot_obj <- down_ggplot_obj + theme (legend.position="none") 

  #  Add axes labels to graph
  down_ggplot_obj <- down_ggplot_obj + ylab (expression ("-" * log[10] ~ " p-value"))  
  
  down_ggplot_obj <- down_ggplot_obj + xlab ("")
  down_ggplot_obj <- down_ggplot_obj + scale_x_discrete (position = "top")
  down_ggplot_obj <- down_ggplot_obj + scale_y_continuous (limits=c(0, max_x_axis_down), breaks = seq (0, max_x_axis_down, by = 2))

  #  If we're plotting only one panel, then halve the width
  if (is.null (up_dataset)) {
    FIGURE_WIDTH = FIGURE_WIDTH / 2
  }
  
  out_fn <- figure_fn
  postscript (file=out_fn, onefile=FALSE, width=FIGURE_WIDTH, height=FIGURE_HEIGHT, paper="special", horizontal=FALSE)
  if (is.null (up_dataset)) {
    multiplot (down_ggplot_obj, cols=1)
  } else {
    multiplot (up_ggplot_obj, down_ggplot_obj, cols=2)
  }
  dev.off ()
}


#####################################################################
##  Process arguments using docopt
#####################################################################

"Usage:  funcchart.R --upregulate UPCHART --downregulate DOWNCHART --kegg KEGG --go GO

Options:
 --upregulate UPCHART  The input file of up-regulated genes, as generated by DAVID
 --downregulate DOWNCHART  The input file of down-regulated genes, as generated by DAVID
 --kegg KEGG  Filename for the KEGG file
 --go GO  Filename for the Gene Ontology file

.
" -> options

# Retrieve the command-line arguments
opts <- docopt (options)

# Update values with the program arguments
UPCHART_ARG <- opts$upregulate
DOWNCHART_ARG <- opts$downregulate
KEGG_ARG <- opts$kegg
GO_ARG <- opts$go


#####################################################################
##  Testing variables
#####################################################################

# UPCHART_ARG <- "Results/combo3/QSC-MPC_MT/david/david-up.tsv"
# DOWNCHART_ARG <- "Results/combo3/QSC-MPC_MT/david/david-down.tsv"
# KEGG_ARG <- "kegg.eps"
# GO_ARG <- "go.eps"


######################################################################
##  Process a sample at a time
######################################################################

##  Note:  upchart is on the left; downchart is on the right
##         If they are combined, there is only downchart; upchart is NULL
downchart <- ReadData (DOWNCHART_ARG)
if (UPCHART_ARG != DOWNCHART_ARG) {
  upchart <- ReadData (UPCHART_ARG)
} else {
  upchart <- NULL
}

##  Re-order the data frames in decreasing order of p-value
downchart <- downchart[order (downchart[["PValue"]]),] 
if (!is.null (upchart)) {
  upchart <- upchart[order (upchart[["PValue"]]),]   
}

if (is.null (upchart)) {
  if ((is.data.frame (downchart)) && (nrow (downchart) != 0)) {
    downchart_subset <- head (downchart[downchart["Category"] == "KEGG_PATHWAY",], NUM_CATEGORIES)
    SaveFigure (outprefix_arg, NULL, downchart_subset, "kegg", KEGG_ARG)

    downchart_subset <- head (downchart[downchart["Category"] != "KEGG_PATHWAY",], NUM_CATEGORIES)
    SaveFigure (outprefix_arg, NULL, downchart_subset, "go", GO_ARG)
  }  
} else {
  if ((is.data.frame (upchart)) && (nrow (upchart) != 0)) {
    if ((is.data.frame (downchart)) && (nrow (downchart) != 0)) {
      upchart_subset <- head (upchart[upchart["Category"] == "KEGG_PATHWAY",], NUM_CATEGORIES)
      downchart_subset <- head (downchart[downchart["Category"] == "KEGG_PATHWAY",], NUM_CATEGORIES)
      SaveFigure (outprefix_arg, upchart_subset, downchart_subset, "kegg", KEGG_ARG)

      upchart_subset <- head (upchart[upchart["Category"] != "KEGG_PATHWAY",], NUM_CATEGORIES)
      downchart_subset <- head (downchart[downchart["Category"] != "KEGG_PATHWAY",], NUM_CATEGORIES)
      SaveFigure (outprefix_arg, upchart_subset, downchart_subset, "go", GO_ARG)
    }
  }  
}


######################################################################
##  Indicate the date of analysis
######################################################################

format (Sys.time(), "Date of analysis:  %a %b %d %X %Y")


