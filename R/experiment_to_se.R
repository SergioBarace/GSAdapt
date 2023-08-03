#' @name experiment_to_se
#' @aliases experiment_to_se
#' 
#' @title Creation of a SummarizedExperiment object from a gene expression data set and its metadata table
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors SimpleList
#'
#' @param gematrix A matrix
#' @param metadata Metadata of gematrix in data frame format
#' 
#' @usage experiment_to_se(gematrix, metadata)
#' 
#' @details The gematrix is a gene expression data frame where gene Symbols are in the first column. Sample names of gematrix must be in the first row. 
#' The metadata is a data frame where sample names are in the first column.
#'
#' @return summarized experiment object
#' @export
#' 
#' @examples 
#' 
#' library(SummarizedExperiment)
#' 
#' gematrix <- as.data.frame(matrix(rnorm(120),nrow=20))
#' colnames(gematrix) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")
#' gematrix <- cbind(seq(1,20,1), gematrix)
#' metadata <- data.frame(Sample=c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6"),
#' Treatment=c("Ctrl","Ctrl","Ctrl","Treated","Treated","Treated"))
#' 
#' se <- experiment_to_se(gematrix, metadata)

experiment_to_se <- function(gematrix = gematrix, 
                            metadata = metadata) {
  if (length(gematrix) == 0L) {
    stop("a Gene Expression matrix is required")
  }
  rownames(gematrix) <- gematrix[,1]
  assays <- as.matrix(gematrix[,-1])
  rownames(metadata) <- metadata[,1]
  colData <- metadata
  rowData <- data.frame(gene_id = rownames(assays))
  rownames(rowData) <- rowData[,1]
  
  if (identical(colnames(assays), rownames(colData))) {
    se <- SummarizedExperiment(
      assays = SimpleList(counts=assays),
      rowData = rowData,
      colData = colData
    )
    return(se)
  }
}
