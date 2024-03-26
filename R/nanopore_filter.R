#' Nanopore_filter
#' @author Auden Block
#' 
#' @param reads A dataframe containing the column names: sequence_length_template, mean_qscore_template, and (if IDs=TRUE) read_ids.
#' @param min_quality (default = 1) Minimum acceptable quality score
#' @param min_length (default = 200) Minimum acceptable read length
#' @param IDs (default = false) By default, dataframe will be returned with the reads that pass filtering. If TRUE, a vector of read_IDs that passed quality scoring will be returned
#' @return A dataframe of reads that pass quality and length filtering, or a character vector of readIDs that pass quality filtering.
#' 
#' @examples 
#' library2_duplex_table <- dorado_table(simplex.df, duplex.df)
#' all.libraries.duplex.table <- as.data.frame(lapply(simplex.list, dorado_table))
#' @export
nanopore_filter <- function(reads, min_quality = 1, min_length = 200, IDs = FALSE) {
  stopifnot("The reads variable must be a dataframe" = is.dataframe(reads))  
  if (isFALSE(IDs)){
    #Remove the rows from the dataframe
    return(reads[which(reads$sequence_length_template >= min_length & reads$mean_qscore_template >= min_quality),])
  } else {
    stopifnot("The read_id column must be present when IDs is true" = )
    return(reads$read_id[which(reads$sequence_length_template >= min_length & reads$mean_qscore_template >= min_quality)])
  }
  
}
