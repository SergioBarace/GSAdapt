#' Signature_input_to_GSAdapt
#' 
#' @name Signature_input_to_GSAdapt
#' @aliases Signature_input_to_GSAdapt
#' @title Conversion of gene set collections
#' 
#' @usage Signature_input_to_GSAdapt(Signature)
#' @param Signature Dataframe with genes and signature columns
#' 
#' @description Convert a dataframe in a list using signature name as factor and genes as features
#'
#' @details Function Signature_input_to_GSAdapt generates a list taking a two-column dataframe. It is supposed that features are in first column and signature name in second one. Features can be expressed as Ensembl ID, gene symbol or other annotated gene expression.
#' @return List with signature information
#' @export
#' @examples 
#' 
#' ## One Gene Set
#' Signature <- data.frame(GeneID = c("A","B"),Sig = rep("GS1",2))
#' Signature <- Signature_input_to_GSAdapt(Signature)
#' 
#' ## Two and more Gene Sets
#' Signature <- data.frame(Gene_ID = c("A","B","C","D"),Sig = c(rep("GS1",2),rep("GS2",2)))
#' Signature <- Signature_input_to_GSAdapt(Signature)

Signature_input_to_GSAdapt <- function(Signature = Signature) {
  if (missing(Signature)) {
    stop("A dataframe with genes and annotation columns is required")
  } else {
    Signature[, 2] <- as.factor(Signature[, 2])
    colnames(Signature) <- c("GeneID", "Gene_set")
    
    Signature <- as.list(Signature)
    

    return(Signature)
  }
}
