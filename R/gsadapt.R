#' gsadapt
#'
#' @param se A summarizedExperiment
#' @param gs_df Dataframe with information of gene_symbol and signature name
#' @param mggcc Median gene-gene correlation Cut-off Numeric value between 0 to 1 
#' @param ncec Node Centrality Eigenvalue Cut-off. Numeric value between 0 to 1
#' @param plot_graph Whether to plot graphs in directory
#'
#' @title Gene Set Adaptation
#' 
#' @description Adaptation of gene set collections to transcriptomic data
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom dplyr filter
#' @importFrom tidyr %>%
#' @importFrom stats na.omit cor
#' @importFrom igraph graph_from_adjacency_matrix layout_with_fr degree simplify eigen_centrality delete.vertices
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics par
#' @import airway
#' 
#' @aliases gsadapt 
#'
#' @return dataframe with adapted gene sets
#' @export
#' 
#' @examples 
#' 
#' ### Two or more signatures (permissive filter)
#' library(SummarizedExperiment)
#' library(airway)
#' 
#' data(airway, package="airway")
#' 
#' Signature <- data.frame(Gene_ID = head(rownames(assay(airway)), n=20), 
#' Gene_Set = c(rep("GS1",10),rep("GS2",10))) 
#' 
#' DataMeta <- gsadapt(se = airway, gs_df=Signature, mggcc=0.2, ncec=0.1, plot_graph=FALSE)
#' 
#' ### Two or more signatures (restrictive filter)
#' library(SummarizedExperiment)
#' library(airway)
#' 
#' data(airway, package="airway")
#' 
#' Signature <- data.frame(Gene_ID = head(rownames(assay(airway)), n=20), 
#' Gene_Set = c(rep("GS1",10),rep("GS2",10))) 
#' 
#' DataMeta <- gsadapt(se = airway, gs_df=Signature, mggcc=0.6, ncec=0.3, plot_graph=FALSE)


gsadapt <- function(se, gs_df, mggcc, ncec, plot_graph = c(FALSE, TRUE)) {
    if (length(assay(se)) == 0L) {
      stop("A Summarized experiment object is required")
    }
  if (length(gs_df) == 0L) {
    stop("A Gene Set data frame object is required")
  }
  
    se_mtx <- as.data.frame(assay(se))
    gs_df[,2] <- as.factor(gs_df[,2])
    gs_list <- as.list(levels(gs_df[, 2]))

    gene_symbols <- list()
    gs_name <- list()

    adapted_gs <- lapply(gs_list, function(j) {
      selected_genes <- gs_df %>% dplyr::filter(gs_df[,2] == j)

      subset <- na.omit(se_mtx[selected_genes[, 1], ])
      subset.t <- t(subset)
      corr_table <- cor(subset.t, method = "spearman")
      corr_table[is.na(corr_table)] <- 0
      corr_table <- abs(corr_table)

      corr_table[corr_table < mggcc] <- 0
      original_gs <- rownames(corr_table)

      if (sum(corr_table) > 10) {
        graph1 <- graph_from_adjacency_matrix(corr_table, weighted = T, mode = "undirected", diag = F)


        gene_names_filtered <- names(which(degree(graph1) > 0))
        net <- simplify(graph1)
        ec <- eigen_centrality(net, directed = F, weights = NULL)
        genes_ec <- ec$vector[ec$vector > ncec]
        gene_symbols[[length(gene_symbols) + 1]] <- names(genes_ec)
        gs_name[[length(gs_name) + 1]] <- paste0(j, "_", mggcc, "_", ncec)
        
        if (plot_graph == T) {
          par(mar=c(0,0,0,0)+.1)
          LO <- layout_with_fr(graph1)
          pdf(paste0(j, "_", mggcc, "_graph1", "_", Sys.Date(), "_", ".pdf"), paper = "a4")
          plot(graph1, layout = LO, vertex.label.cex=1.0, vertex.size=8, vertex.color="salmon", vertex.shape="sphere")
          dev.off()
          
          if (length(gene_names_filtered) > 0) {
            if (length(gene_names_filtered) < length(original_gs)) {
              if (length(gene_names_filtered) >= 7) {
                Isolated <- which(degree(graph1) == 0)
                graph2 <- delete.vertices(graph1, Isolated)
                
                if (plot_graph == T) {
                  LO2 <- LO[-Isolated, ]
                  if (length(LO2) > 0) {
                    pdf(paste0(j, "_", mggcc, "_graph2", "_", Sys.Date(), "_", ".pdf"), paper = "a4")
                    plot(graph2, layout = LO2, vertex.label.cex=1.0, vertex.size=8, vertex.color="salmon", vertex.shape="sphere")
                    dev.off()
                  }
                }
                
                x <- exists("graph2")
                if (x == FALSE) {
                  stop("Graph2 is not computed or number of filtered genes is not enough")
                } else {
                  
                  Ext_genes <- names(ec$vector[ec$vector < ncec])
                  degree_Ext_genes <- degree(graph1)[c(Ext_genes,names(Isolated))]
                  LO3 <- LO[-degree_Ext_genes, ]
                  if (length(LO3) > 0) {
                    pdf(paste0(j, "_", mggcc, "_graph3", "_", Sys.Date(), "_", ".pdf"), paper = "a4")
                    plot(graph1, layout = LO3, vertex.label.cex=1.0, vertex.size=8, vertex.color="salmon", vertex.shape="sphere")
                    dev.off()
                  }
                }
                
                
              } else {stop("There were not enough genes that passed the mggc cut-off. The gene set cannot be adapted using the provided gene expression matrix ")}
            } else {stop("All genes in the input gene set were above the mggc cut-off. The adapted gene set is the same as the original gene set")}
          } else {stop("None of the genes had a mggc above the cut-off. This signature cannot be adapted using the provided gene expression matrix")}
          
        }
        
        names(gene_symbols) <- gs_name
        splited.list <- split(gene_symbols, names(gene_symbols))
        
        genes <- c()
        name <- c()
        
        temporal_list <- lapply(splited.list, function(splited.elem) {
          list <- unlist(splited.elem)
          genes <- c(genes, list)
          name <- c(name, rep(names(splited.elem), length(list)))
          temporal_list <- list(GeneID = genes, Signature = name)
          temporal_list <- data.frame(temporal_list)
        })
      

        
      } else {"Correlation matrix size of this gene set does not achieve the minimum sum of correlations."}
    })
    
    
    return(adapted_gs)
  }

