#' @name Matrix_to_se
#' @aliases Matrix_to_se
#' 
#' @title Conversion of dataframe into SummarizedExperiment object
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors SimpleList
#'
#' @param inp_df A matrix
#' @param chromosome_field The name of the column related to chromosome
#' @param ranges IRanges from each feature
#' @param strand The name of the column related to strand
#' @param id gene_id or feature name
#' @param condition Metadata of inp_df in dataframe format
#' 
#' @usage Dataframe_to_se(inp_df,chromosome_field,ranges,strand,id,condition)
#' 
#' @details Conversion using a expression matrix containing information about gene ranges and experimental condition of each samples of expression matrix. For creating genomic ranges chromosome location, strand, start and end gene locations are required to construct gene ranges. ColData is compulsory and should give enough information about sample condition. Names of samples must be identical between matrix colnames and colData sample names.
#'
#' @return summarized experiment object
#' @export
#' 
#' @examples 
#' 
#' library(GenomicRanges)
#' library(SummarizedExperiment)
#' library(S4Vectors)
#' 
#' matrix <- matrix(rnorm(120),nrow=20)
#' colnames(matrix) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")
#' rownames(matrix) <- seq(1,20,1)
#' 
#' chromosome_field <- rep("chr1",20)
#' ranges <- IRanges(floor(runif(20,1e5, 1e6)),width=100)
#' strand <- sample(c("+", "-","*"), 20, TRUE)
#' id <- rownames(matrix)
#' 
#' condition <- data.frame(Sample=c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6"),
#' Treatment=c("Ctrl","Ctrl","Ctrl","Treated","Treated","Treated"))
#' 
#' se <- Dataframe_to_se(matrix, chromosome_field, ranges, strand, id, condition)

Dataframe_to_se <- function(inp_df = inp_df, 
                            chromosome_field = chromosome_field,
                            ranges = ranges,
                            strand = strand,
                            id = id,
                            condition = condition) {
  if (length(inp_df) == 0L) {
    stop("Dataframe is required")
  }

  assays <- inp_df

  if (identical(rownames(inp_df),id)) {
    rowRanges <- GRanges(seqnames = paste0(chromosome_field,"_",id),
                         ranges = ranges,
                         strand = strand)
  }

  colData <- condition

  if (identical(ncol(assays), nrow(condition))) {
    se <- SummarizedExperiment(
      assays = SimpleList(counts=assays),
      rowRanges = rowRanges,
      colData = colData
    )

    return(se)
   
  }
}
