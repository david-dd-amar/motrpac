#' A wrapper for preprocessCore's quantile normalization.
#' Comment: we do not use this by default
run_quantile_normalization<-function(x){
  x = as.matrix(x)
  mode(x) = "numeric"
  newx = preprocessCore::normalize.quantiles.robust(x)
  rownames(newx) = rownames(x)
  colnames(newx) = colnames(x)
  return(newx)
}

#' Takes an FPKM matrix, removes lowly expressed genes and log transform
#' the remaining matrix. This can be used on the FPKM output matrix
#' from the RNA-seq pipeline.
#' @return A matrix of log transformed FPKMs
process_fpkm1 <-function(fpkm_matrix, intensity_threshold=0,intensity_pct=0.2){
  lowly_expressed_genes = rowSums(
    fpkm_matrix==intensity_threshold)/ncol(fpkm_matrix) > intensity_pct
  fpkm_matrix = fpkm_matrix[!lowly_expressed_genes,]
  fpkm_matrix = log(fpkm_matrix+1,base = 2)
  return(fpkm_matrix)
}

try(library(NOISeq))
run_tmm_normalization<-function(x,...){
  return(NOISeq::tmm(x,...))
}

