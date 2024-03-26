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
      
      return(
          duplex[c(which(duplex$template_id %in% duplex$complement_id),which(duplex$complement_id %in% duplex$template_id)),])
    } else {
      stopifnot("Duplex dataframe missing the column read_id" = ("read_id" %in% colnames(duplex)))
      str_split <- as.data.frame(do.call(rbind, strsplit(duplex$read_id, ";")))
      return(
        duplex[c(which(str_split$V1 %in% str_split$V2),which(str_split$V2 %in% str_split$V1)),]
        )
      
    }
}



