#' Gene Set Adaptation of one gene set
#' 
#' @name gsadapt_of_one_gs
#' @aliases gsadapt_of_one_gs
#' 
#' @importFrom igraph graph_from_adjacency_matrix degree eigen_centrality simplify
#' @importFrom stats cor
#'
#' @param abs_mtx Absolute correlation matrix
#' @param sig_name The name of the gene set
#' @param CCIM Correlation cut-off
#' @param NCEC Centrality cut-off
#' 
#' @details This function performs a simply gene set adaptation using graph models. CCIM and NCEC correspond to adaptation filters able to remove weak associated genes and poor connected genes into resultant graph g1. 
#'
#' @return Dataframe with adapted genes
#' @export
#'
#' @examples
#' 
#' library("igraph")
#' 
#' Mtx <- matrix(rnorm(25),ncol=5)
#' 
#' Abs_mtx <- abs(cor(Mtx))
#' dimnames(Abs_mtx) <- list(c("Gene1","Gene2","Gene3","Gene4","Gene5"),
#' c("Gene1","Gene2","Gene3","Gene4","Gene5"))
#' 
#' Adapt_df <- gsadapt(abs_mtx=Abs_mtx,sig_name="GS1",CCIM=0.1,NCEC=0.2)
#' 

gsadapt <- function(abs_mtx = abs_mtx,
                    sig_name = sig_name,
                    CCIM = CCIM,
                    NCEC = NCEC){
  
  abs_mtx[abs_mtx < CCIM] <- 0
  
  genes_input <- rownames(abs_mtx)
  
  g1 <- graph_from_adjacency_matrix(abs_mtx, weighted = NULL, mode = "directed", diag = F)
  
  gene_names_filtered <- names(which(degree(g1) > 0))
  
  net <- simplify(g1)
  
  ec <- eigen_centrality(net, directed = F, weights = F)
  
  genes_ec <- ec$vector[ec$vector > NCEC]

  Tpm <- data.frame(Gene_ID = names(genes_ec), Signature = rep(sig_name, length(names(genes_ec))))
  
  return(Tpm)
}
