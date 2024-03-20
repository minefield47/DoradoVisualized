#' Subset __X coverage using Nanopore Reads
#' @author Auden Block
#' 
#' @param reads A dataframe or list of dataframes containing the column names: read_id and sequence_length_template.
#' @param genome_size_mb Estimated genome size in megabases.

#' @return A vector of subset read_ids/indexes whose length is equal to the desired coverage of a given genome.
#' 
#' @examples 
#' 2x.subset <- nanopore_subset(reads, 2, 1500, TRUE, 10)

#' @export
coverage <- function(reads, genome_size_mb) {
  if (lapply(reads,is.data.frame) |> unlist() |> all()) {
    
    return((lapply(reads, function(i) {colSums(i['sequence_length_template'])}) |> unlist() |> sum())/(1600*1000000)
           )
  } else {
    if (is.data.frame(reads)){
      return(sum(reads$sequence_length_template))
    } else {
      return(sum(reads))
    }
  }
}