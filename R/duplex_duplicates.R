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
duplex.duplicates <- function(duplex) {
  #duplex = duplex dataframe. 
  if (!is.data.frame(duplex)) { #If either the duplex or simplex parameter given is not a dataframe, return an error. 
    stop("The duplex parameter must be a dataframe.")
    
  } else {
  
    if ("complement_id" %in% colnames(duplex)) {
      #If duplex.parents has already been done on the dataframe, the template and complement ID columns already exist for utilization. 
      return(duplex[
        c(which(duplex$template_id %in% duplex$complement_id & duplex$complement_id%in%duplex$template_id), #Find the complements reads also in the template read column
          which(duplex$template_id %in% duplex$complement_id[which(duplex$template_id %in% duplex$complement_id & duplex$complement_id %in% duplex$template_id)])),]) #Find the templates that are also in the complements (This allows both to be returned instead of just the first match)
      
    } else {
      #If duplex.parents has not been ran, first need to split: 
      
      #Find the complements
      complement <- sapply(str_split(duplex$read_id, ";"), `[`,1)
      #Find the templates
      template <- sapply(str_split(duplex$read_id, ";"), `[`,2)
      #Find the row indices of reads are found in both complements and templates. 
      return(duplex[
        c(which(template %in% complement & complement %in% template),
          which(template %in% duplex$complement[which(template %in% complement & complement %in% template)])),])
    }
    
    
    
  }
}


