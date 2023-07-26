#' Generation of absolute correlation matrixes
#'
#' @name Absolute_correlation_matrix
#' @aliases Absolute_correlation_matrix
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom stats na.omit
#' 
#' @param input_se A SummarizedExperiment object  
#' @param list_gs A list of gene set collection
#'
#' @return A list of matrixes
#' @export
#'
#' @examples
#' 
#' ## One Gene Set
#' 
#' library(airway)
#' data(airway)
#' 
#' Signature <- list(GS1=names(airway@rowRanges)[1:6])
#' 
#' Mtx_1 <- Absolute_correlation_matrix(airway,Signature)
#' 
#' ## Two or more Gene Sets
#' 
#' data(airway)
#' Signatures <- list(GS1=names(airway@rowRanges)[1:6],GS2=names(airway@rowRanges)[7:12])
#' 
#' Mtx_2 <- Absolute_correlation_matrix(airway,Signatures)
#' 
Absolute_correlation_matrix <- function(input_se = input_se,
                                        list_gs = list_gs){
  
  if (length(assay(input_se)) == 0L) {
    stop("Summarized experiment object is required")
  }
  
  assay_obj <- assay(input_se)
  
  if (length(list_gs) == 0L) {
    stop("Signature list of gene set collection is not provided or it is empty")
  }
  
  abs_mtx_list <- list()
  
  j <- seq(1,length(list_gs),1)
  
  abs_mtx_list <- lapply(list_gs, function(j){
    selected_genes <- c(unlist(list_gs[j]))
    
    subset <- na.omit(assay_obj[selected_genes,])
    subset.t <- t(subset)
    corr_table <- cor(subset.t, method = "spearman")
    corr_table[is.na(corr_table)] <- 0
    corr_table <- abs(corr_table)
    
    abs_mtx_list[[length(abs_mtx_list) + 1]] <- corr_table
  })
  
  return(abs_mtx_list)
  
}
