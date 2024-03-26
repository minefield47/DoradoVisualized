#' Simplex_only
#' @author Auden Block
#' @description
#' Subset the simplex reads dataframe to contain only reads from which no duplex read was generated.
#' 
#' @param duplex A vector or dataframe (containing the column "read_id" or c("complement_id", "template_id")) of duplex reads.
#' @param simplex A vector or dataframe of simplex reads including from those from which the duplex reads spawned.
#' @return A character vector or subsetted dataframe of simplex reads from which NO duplex reads were generated.
#' 
#' @examples
#' simplex_only(summary.list$duplex.df, summary.list$simplex.df)
#' @export simplex_only
simplex_only <- function(duplex, simplex){
  stopifnot("The duplex parameter must be a vector" = is.data.frame(duplex))
  if (is.data.frame(simplex)){
  
  if ("complement_id" %in% colnames(duplex)) {
    #If duplex.parents has already been done on the dataframe, the template and complement ID columns already exist for utilization. 
    
    
    return(simplex[! simplex$read_id %in% c(duplex$template_id, duplex$complement_id)])
    
    
  
  } else {
    #Find the complements/templates when complement_id/template_id not in the duplex dataframe.
    
    return(simplex[! simplex$read_id %in% as.data.frame(do.call(rbind, strsplit(duplex$read_id, ";")))])
    
    }  
  } else {
    stopifnot("Simplex must be a character vector" = is.character(simplex))
    if ("complement_id" %in% colnames(duplex)) {
      #If duplex.parents has already been done on the dataframe, the template and complement ID columns already exist for utilization. 
      
      
      return(simplex[!(simplex %in% c(duplex$template_id, duplex$complement_id))])
      
      
      
    } else {
      #Find the complements/templates when complement_id/template_id not in the duplex dataframe.
      
      return(simplex[! simplex %in% as.character(do.call(rbind, strsplit(duplex, ";")))])
      
    }  
  }
}
