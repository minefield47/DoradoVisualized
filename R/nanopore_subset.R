#' Subset __X coverage using Nanopore Reads
#' @author Auden Block
#' 
#' @param reads A dataframe or list of dataframes containing the column names: read_id and sequence_length_template.
#' @param coverage The desired coverage.
#' @param genome_size_mb Estimated genome size in megabases.
#' @param ID (Default = TRUE) Return a character vector of readIDs or the index of reads selected, if false. 
#' @return A vector of subset read_ids/indexes whose length is equal to the desired coverage of a given genome.
#' 
#' @examples 
#' 2x.subset <- nanopore_subset(reads.df, 2, 1500, TRUE)

#' @export
nanopore_subset <- function(reads, coverage, genome_size_mb, ID = TRUE) {
  total_bases <- coverage * genome_size_mb*1000000
  
  stopifnot("object reads is not in dataframe format" = lapply(reads,is.data.frame))
  stopifnot("Sum read length is not enough for the given coverage/genome size. Please select a different coverage or genome size." =  total_bases < (do.call(rbind, lapply(reads, function(i)colSums(i['sequence_length_template'])))))
  
  if (is.list(reads)) {
    #If a list is given, need to randomly subset from any element of list. 
    
    #Get some summary stats
    list_length <- length(reads)
    
    #Create a blank list to store index values.
    list_row_index <- vector(mode = "list", length =list_length)
    
    #A random list
    list_index <- sample(1:list_length, 1)
    #Within the list_index, find a random row.
    value <- sample(1:nrow(reads[[list_index]]), 1)
    
    #Store the row value
    list_row_index[[list_index]] <- append(list_row_index[[list_index]], value)
    
    #Store the read length
    read.length <- reads[[list_index]]$sequence_length_template[value]
    
    #Remove the read from being sampled twice.
    reads[[list_index]] <- reads[[list_index]][-value,]
    
    
    #Repeat until coverage has been reached.
    while (read.length < total_bases){
      list_index <- sample(1:list_length, 1)
      value <- sample(1:nrow(reads[[list_index]]), 1)
      list_row_index[[list_index]] <- append(list_row_index[[list_index]], value)
      read.length <- read.length + reads[[list_index]]$sequence_length_template[value]
      reads[[list_index]] <- reads[[list_index]][-value,]
      
    }
    #Print the values 
    if (isTRUE(ID)){
      return(mapply(function(reads,list_row_index) {
        reads$read_id[list_row_index]
      }, reads, list_row_index))
    } else {
      return(do.call("rbind", mapply(function(reads,list_row_index) {
        reads[list_row_index,]
      }, reads, list_row_index, SIMPLIFY = FALSE) ))
    }
    
  } else{
    value <- sample(1:nrow(reads), 1)
    c_row_index <- value
    read.length <- reads$sequence_length_template[value]
    reads <- reads[-value,]
    while (read.length <= total_bases){
      value <- sample(1:nrow(reads), 1)
      c_row_index <- append(c_row_index, value)
      read.length <- read.length + reads$sequence_length_template[value]
      reads <- reads[-value,]
    }
    
    if (isTRUE(ID)){
      return(reads$read_id[c_row_index])
    } else {
      return(reads$read_id[c_row_index])
    }
    
  }
  
  
}
