#####################################################################
##
##  General functions for clustering
##
#####################################################################


##
##  Calculate hierarchical clustering
##
plotHC = function (outdir, genes_local, method, method_full, suffix="") {
  subtext <- method_full
  
  dist_matrix <- Dist (t (genes_local), method=method) 
  fn <- paste (outdir, "/", "heatmap-", method, suffix, ".eps", sep="")
#   par (mar=c(7,4,4,2)+0.1) 
  postscript (file=fn, onefile=FALSE, width=HEATMAP_WIDTH, height=HEATMAP_HEIGHT, paper="special", horizontal=FALSE)
  heatmap (as.matrix (dist_matrix), main=subtext, margins=c(12,8))
  
# heatmap.2(sel, col=redgreen(75), scale="row", ColSideColors=col, key=TRUE, symkey=FALSE, density.info="none",cexRow=1,cexCol=1,margins=c(12,8),trace="none",srtCol=45)
  
  dev.off ()

  hc <- hclust (dist_matrix, "ave")
  fn <- paste (outdir, "/", "hc-", method, suffix, ".eps", sep="")
  postscript (file=fn, onefile=FALSE, width=HCLUST_WIDTH, height=HCLUST_HEIGHT, paper="special", horizontal=FALSE)
  plot (hc, sub=subtext, xlab="")
  dev.off ()
}


##
##  Clean up sample names for PCA
##
cleanCategories = function (categories) {
  categories <- gsub ("_1A", "", categories)
  categories <- gsub ("_1B", "", categories)
  categories <- gsub ("_2A", "", categories)
  categories <- gsub ("_2B", "", categories)
  categories <- gsub ("_1", "", categories)
  categories <- gsub ("_2", "", categories)
  categories <- gsub ("_3", "", categories)
  categories <- gsub ("_A", "", categories)
  categories <- gsub ("_B", "", categories)
  categories <- gsub ("_C", "", categories)
  
  return (categories)
}


filterPublicData = function (categories) {
  categories <- gsub ("^Chon$", "Others", categories)
  categories <- gsub ("^CMC$", "Others", categories)
  categories <- gsub ("^ESC$", "Others", categories)
  categories <- gsub ("^HDF$", "Others", categories)
  categories <- gsub ("^HEP$", "Others", categories)
  categories <- gsub ("^HSC$", "Others", categories)
  categories <- gsub ("^KiEC$", "Others", categories)
  categories <- gsub ("^MB$", "Others", categories)
  categories <- gsub ("^MCSC$", "Others", categories)
  categories <- gsub ("^MT$", "Others", categories)
  categories <- gsub ("^NSC$", "Others", categories)
  categories <- gsub ("^SC$", "Others", categories)
  categories <- gsub ("^SMC$", "Others", categories)
  
  return (categories)
}


##
##  Calculate PCA
##
createPCA = function (genes_local) {
  pca_obj <- prcomp (t (genes_local), center=TRUE, scale=TRUE)
  print (summary (pca_obj))

  return (pca_obj)
} 


showPCAVariance = function (outdir, pca_obj, suffix="") {
  # http://www.sthda.com/english/wiki/principal-component-analysis-in-r-prcomp-vs-princomp-r-software-and-data-mining
  eig <- (pca_obj$sdev)^2
  # Variances in percentage
  variance <- eig*100 / sum (eig)
  # Cumulative variances
  cumvar <- cumsum (variance)
  eigen.pca <- data.frame (eig = eig, variance = variance, cumvariance = cumvar)
  
  ##  Plot the variances
  fn <- paste (outdir, "/", "variances", suffix, ".eps", sep="")
  postscript (file=fn, onefile=FALSE, width=10, height=7, paper="special", horizontal=FALSE)
  # plot (pca, type="l", main=subtext)
  barplot (eigen.pca[, 2], names.arg=1:nrow(eigen.pca), main = "Variances", xlab = "Principal Components", ylab = "Percentage of variances", col ="black")
  #lines (x = 1:nrow(eigen.pca), eigen.pca[, 2], type="b", pch=19, col = "red")
  dev.off ()  
}


showPCAPlot = function (outdir, genes_local, pca_obj, suffix="") {
  ##  Plot the PCA
  fn <- paste (outdir, "/", "pca", suffix, ".eps", sep="")
  postscript (file=fn, onefile=FALSE, width=10, height=7, paper="special", horizontal=FALSE)
  
  # http://stackoverflow.com/questions/16653581/color-pca-depending-on-predefined-groups
  categories_public <- colnames (genes_local)
  categories_public <- cleanCategories (categories_public)
  categories_public <- filterPublicData (categories_public)
  
  ggplot_obj <- ggplot (pca_obj, aes (x=pca_obj$x[,1], y=pca_obj$x[,2], colour=categories_public)) 
  ggplot_obj <- ggplot_obj + geom_point (shape=20, size=3) 
  ggplot_obj <- ggplot_obj + xlab ("Principal Component #1") + ylab ("Principal Component #2")
  ggplot_obj <- ggplot_obj + labs (color="Legend") 
  print (ggplot_obj)

  dev.off ()

  return
}

