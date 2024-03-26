#' Subset __X coverage using Nanopore Reads
#' @author Auden Block
#' 
#' @param reads A vector or list of dataframes (for users with multiple libraries wishing to calculate the total coverage, the sequence lengths must be labelled "sequence_length_template"). 
#' @param genome_size Estimated genome size.
#' @param units is the genome in bases, kilobases, megabases, or gigabases (default = "gb")

#' @return A vector of subset read_ids/indexes whose length is equal to the desired coverage of a given genome.
#' 
#' @examples 
#' 2x.subset <- nanopore_subset(reads, 2, 1500, TRUE, 10)

#' @export
coverage <- function(reads, genome_size, units = "gb") {
  
  stopifnot("Size parameter must be b, kb, mb, gb" = units %in% c("b", "kb", "mb","gb"))
  size <- switch(
    units, 
    "b" = 1,
    "kb" = 1e3,
    "mb"= 1e6,
    "gb"= 1e9,
  )
  
  if (lapply(reads,is.data.frame) |> unlist() |> all()) {
    
    return((lapply(reads, function(i) {colSums(i['sequence_length_template'])}) |> unlist() |> sum())/(genome_size*size)
           )
  } else {
    stopifnot("Vector input not given." = is.vector(reads))
      return(sum(reads))
  }
}