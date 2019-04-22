#####################################################################
##
##  General functions
##
#####################################################################

isFALSE = function (x) {
  return (!isTRUE (x))
}

writeTable = function (dataset, filename) {
  if (class (dataset) == "matrix") {
    #  Need to write out the blank cell at the top left
    write.table (dataset, filename, sep='\t', row.names = TRUE, col.names = NA, quote = FALSE)
  } else if (class (dataset) == "data.frame") {
    write.table (dataset, filename, sep='\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
  } else if (class (dataset) == "character") {
    #  List of strings (so actually the separator character is meaningless); drop the column names
    write.table (dataset, filename, sep='\t', row.names = FALSE, col.names = FALSE, quote=FALSE)
  } else {
    stop ("Unknown file format for function call to writeTable ().")
  }
}


