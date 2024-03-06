#' test_rank_ags
#'
#' @param input_df A dataframe with two columns. This table is possible to obtain from gsadapt function
#' @param additional_sig Dataframe containing gene set collection out of GSEA collection
#' @param gs_column Name of the column of additional_sig file containing gene set names out of GSEA collection
#' @param gene_symbol Name of the column of additional_sig file. Gene name.
#' @param species Parameter of msigdbr
#' @param category Parameter of msigdbr
#' @param subcategory Parameter of msigdbr
#' 
#' @title Test of rank position
#' 
#' @description Evaluation of rank position of each adapted gene set using hypergeometric formula
#' 
#' @importFrom msigdbr msigdbr
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom tidyr %>%
#' @importFrom clusterProfiler enricher
#' 
#' @aliases rank_test_ags
#'
#' @return list of dataframes with information about rank, gene symbol and gene size based on hypergeometric function
#' @export
#'
#' @examples
#' 
#' ## One signature from GSEA 
#' library(msigdbr)
#' library(clusterProfiler)
#' 
#' input_df <- msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL)
#' input_df <- input_df %>% dplyr::select(gene_symbol, gs_name) 
#' input_df <- input_df %>% dplyr::filter(gs_name=="HALLMARK_ANGIOGENESIS")
#' 
#' # The first ten genes from angiogenesis gene set was selected as an example
#' input_df <- as.data.frame(input_df[1:10,])
#' 
#' Summary_table <- test_rank_ags(input_df,species="Homo sapiens",category="H",subcategory=NULL,additional_sig=NULL)
#'
#'## Multiple signatures from GSEA
#'
#'input_df <- msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL) 
#'input_df <- input_df %>% dplyr::select(gene_symbol, gs_name) 
#'random_genes <- sample(input_df$gene_symbol,2000)
#'input_df <- input_df %>% dplyr::filter(gene_symbol %in% random_genes)
#'
#'Summary_table <- test_rank_ags(input_df,species="Homo sapiens",category="H",subcategory=NULL,additional_sig=NULL)
#' 

test_rank_ags <- function(input_df,
                          additional_sig=c(NULL,additional_sig),
                          gs_column=c(NULL,gs_column),
                          gene_symbol=c(NULL,gene_symbol),
                          species=c(NULL,"Homo sapiens","Mus musculus"),
                          category=c(NULL,"H","C1","C2","C3","C4","C5","C6","C7","C8"),
                          subcategory=c(NULL,"CGP","CP","MIR","TFT","3CA","CGN","CM","GO","HPO","IMMUNESIGDB","VAX")){
  if (length(input_df) == 0L) {
    stop("input_df argument is empty")
  }
  
  if (ncol(input_df) >2) {
    stop("Dataframe with two columns in required")
  }
  
  colnames(input_df) <- c("GeneID","Signature")
  
  if(length(subcategory) == 0L){
    Msigdb_cat_df <- msigdbr(species = species, category = category)
  } else {
    Msigdb_cat_df <- msigdbr(species = species, category = category, subcategory = subcategory)
  }
   
  m_t2g <- Msigdb_cat_df[,c(3,4)]
  
  gene_set_list <- as.list(levels(as.factor(input_df$Signature)))
  
  if(is.null(additional_sig)){
    message("Additional information is not supplied")
  } else {
    additional_sig <- additional_sig %>% dplyr::select(gs_column,gene_symbol)
    
    m_t2g <- rbind(Msigdb_cat_df,additional_sig)
  }
  
  rank_GS <- vector()
  ags_name <- vector()
  
  columns = c("ags_name","rank_GS") 
  empty_df = data.frame(matrix(nrow = 1, ncol = length(columns))) 
  colnames(empty_df) = columns
  
  summary_ranks_ags <- lapply(gene_set_list, function(i){
    
    summary_test <- enricher(input_df$GeneID, TERM2GENE=m_t2g)
    rank_gsea <- summary_test@result
    rank_gsea_position <- which(rank_gsea$ID == i)
    
    ifelse(length(rank_gsea_position) == 0,
           rank_gsea_position <- "Out of list",
           rank_gsea_position <- rank_gsea_position)
    
    rank_GS <- c(rank_GS,rank_gsea_position)
    ags_name <- c(ags_name,i)
    
    summary_ranks_ags_temp <- data.frame(ags_name, rank_GS)
    summary_ranks_ags <- rbind(empty_df,summary_ranks_ags_temp)
    summary_ranks_ags <- na.omit(summary_ranks_ags)
    
  })
  
  empty_list_genes <- vector()
  gene_set_size <- vector()
  
  summary_ranks_genes_ags <- lapply(gene_set_list, function(i){
    
    Signature <- input_df$Signature
    genes <- input_df %>% dplyr::filter(Signature == i)
    vector_char <- paste(as.character(genes$GeneID),sep = "",collapse = ",")
    empty_list_genes <- c(empty_list_genes,vector_char)
    gene_set_size <- c(gene_set_size, length(genes$GeneID))
    
    summary_ranks_genes_ags <- data.frame(Gene_ID = empty_list_genes,gs_size=gene_set_size)
  })
  
  final_list <- list(summary_ranks_ags,summary_ranks_genes_ags)
  return(final_list)
}
