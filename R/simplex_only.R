#' Find the simplex reads from which no duplex read was generated.
#' @author Auden Block
#' 
#' @param duplex.df A dataframe of duplex reads
#' @param simplex.df A dataframe of simplex reads of which the duplex reads spawned
#' @return A dataframe of the simplex reads from which NO duplex reads were generated.
#' 
#' @examples 
#' simplex.only <- simplex.only(duplex.df, simplex.df);
#' summary.list$simplex.only <- simplex.only(summary.list$duplex.df, summary.list$simplex.df);
#' @export
simplex_only <- function(duplex, simplex){
  stopifnot("The duplex parameter must be a dataframe." = is.data.frame(duplex))
  stopifnot("The simplex parameter must be a dataframe." = is.data.frame(simplex))
  
  if ("complement_id" %in% colnames(duplex)) {
    #If duplex.parents has already been done on the dataframe, the template and complement ID columns already exist for utilization. 
    
    
    return(subset(simplex, !(simplex$read_id %in% c(duplex$template_id, duplex$complement_id))))
    
    
  
  } else {
    #Find the complements/templates
    values <- unlist(do.call(rbind, strsplit(duplex$read_id, ";")))

    
    return(subset(simplex, !(simplex$read_id %in% values)))
    
    }
}
