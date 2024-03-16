#' Find the duplex reads in which a simplex read was used as a template for one duplex read, and a complement for a different duplex read.
#' @author Auden Block
#' 
#' @param duplex.df A dataframe of duplex reads
#' @return A subset of the input duplex dataframe in which a read is used as both a template and a complement. Useful for finding false positive duplex reads
#' 
#' @examples 
#' duplex.duplicates <- duplex.duplicates(duplex.df);
#' summary.list$duplex.duplicates <- duplex.duplicates(summary.list$duplex);
#' @export
duplex_duplicates <- function(duplex) {
  #duplex = duplex dataframe. 
  stopifnot("The duplex parameter must be a dataframe." = is.data.frame(duplex))
  
    if ("complement_id" %in% colnames(duplex)) {
    #Since duplex_parents has already been run...no need to utilize processing time trying to subset the dataframe.   
      
      #Find the intersect (those that appear in both) between template id and complement id
      duplicates <- intersect(duplex$complement_id, duplex$template_id)
      #return the matches in either the complement id or template id. 
      return(unique( #Because the rbind returns duplicates that link them together, use unique to not oversaturate duplicates
        rbind(
          duplex[which(duplex$template_id %in% duplicates),], 
          duplex[which(duplex$complement_id %in% duplicates),]
              )
            )
          )
    } else {
      #If duplex.parents has not been ran, first need to split: 
      
      #Find the complements
      str.split <- as.data.frame(do.call(rbind, strsplit(duplex$read_id, ";")))
      #Find the row indices of reads are found in both complements and templates. 
      
      duplicates <- intersect(str.split$V1, str.split$V2)
      
     return(unique(
        rbind(
          duplex[which(str.split$V1 %in% duplicates),], 
          duplex[which(str.split$V2 %in% duplicates),]
            )
          ))
      
    }
}



