#' Subset __X coverage using Nanopore Reads
#' @author Auden Block
#' 
#' Export the readIDs for a given dataframe of reads. 
#' This script is typically ran after nanopore_subset to save the randomly generated subset of readIDs for utilization in Shell/Bash. 
#' 
#' @param reads A dataframe containing the read_ids you wish to export. 
#' @param file_name (Default = paste0(as.character(substitute(reads)),"_readIDs.txt")) Name for the file to be created. 
#' @return A txt file of readIDs. 
#' 
#' @examples 
#' 2x.subset <- nanopore_subset(reads, 2, 1500, TRUE)

#' @export
write_readIDs <- function(reads, file_name = "NULL") {
  
  if (is.list(reads)) {
    dir.create(as.character(substitute(reads)))
    
    if (file_name == "NULL") {
      #No filenames provided, lets create our own. 
      if (is.null(names(reads))) { 
        #if the elements of reads are unset, use the run_id. 
        file_name <- paste0("./", 
                            as.character(substitute(reads)), 
                            "/",
                            unlist(sapply(seq_along(reads), 
                                          function(i) paste0(as.character(reads[[i]]$run_id[1]),"_readIDs.txt"))))
        } else {
             #Else, use the element names.
             file_name <- paste0("./", 
                                 as.character(substitute(reads)), 
                                 "/", 
                                 unlist(sapply(names(reads), 
                                               function(i) paste0(i,"_readIDs.txt"))))
       }
    }
    #Error catching
    stopifnot("The length of file_name must be the same length as reads" = identical(length(reads), length(file_name)))
    invisible(mapply(function(reads, file_name) {write.table(reads, file = file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)}, reads = reads, file = file_name))
    
  } else {
    stopifnot("Dataframe or list input must be provided" = is.data.frame(reads))
    if (is.null(file_name)){
    file_name <- paste0(as.character(reads$run_id[1]),"_readIDs.txt")
    }
    write.table(reads$read_id, file = file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
}