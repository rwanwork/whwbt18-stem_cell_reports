#####################################################################
##
##  General functions for cleaning up labels
##
#####################################################################

cleanUpBasic = function (gene_matrix) {
  colnames (gene_matrix) <- gsub ("FPKM.", "", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("_s", "", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("asc", "ASC", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("sb", "SB", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("qsc", "QSC", colnames (gene_matrix))

  return (gene_matrix)
}


shortToLongLabels = function (gene_matrix) {  
  colnames (gene_matrix) <- gsub ("^HDF_", "Dermal Fibroblasts ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^Chon_", "Chrondocytes ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^MCSC_", "Mesenchymal Stem Cells ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^SMC", "Smooth Muscle Cells", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^CMC", "Cardiomyocytes", colnames (gene_matrix))

  colnames (gene_matrix) <- gsub ("^SC_", "Myosatellite Cells ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^MT_", "Primary Myotubes ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^MB_", "Primary Myoblasts ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^ESC_", "Embryonic Stem Cells ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^KiEC_", "Primary Kidney Epithelial Cells ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^HEP_", "Primary Hepatocytes ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^HSC", "Hematopoietic Stem Cells", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^NSC_", "Neural Stem Cells ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^MPC_MT_", "MPCs 4 day differentiated control ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^MPC_", "MPCs control ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^ASC_SB_", "Activated Muscle Stem Cells P38 treated ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^ASC_", "Activated Muscle Stem Cells ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^QSC_", "Quiescence Muscle Stem Cells ", colnames (gene_matrix))
  
  return (gene_matrix)
}


shortToLongLabelsAggr = function (gene_matrix) {  
  colnames (gene_matrix) <- gsub ("^HDF_", "Dermal Fibroblasts ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^Chon_", "Chrondocytes ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^MCSC_", "Mesenchymal Stem Cells ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^SMC", "Smooth Muscle Cells", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^CMC", "Cardiomyocytes", colnames (gene_matrix))

  colnames (gene_matrix) <- gsub ("^SC_", "Myosatellite Cells ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^MT_", "Primary Myotubes ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^MB_", "Primary Myoblasts ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^ESC_", "Embryonic Stem Cells ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^KiEC_", "Primary Kidney Epithelial Cells ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^HEP_", "Primary Hepatocytes ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^HSC", "Hematopoietic Stem Cells", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^NSC_", "Neural Stem Cells ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^MPC_MT_", "MPCs 4 day differentiated control ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^MPC_", "MPCs control ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^ASC_SB_", "Activated Muscle Stem Cells P38 treated ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^ASC_", "Activated Muscle Stem Cells ", colnames (gene_matrix))
  colnames (gene_matrix) <- gsub ("^QSC_", "Quiescence Muscle Stem Cells ", colnames (gene_matrix))
  
  return (gene_matrix)
}
