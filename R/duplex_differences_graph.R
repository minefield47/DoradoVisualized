#' Dataframe containing a summary table about a Dorado run.
#' @author Auden Block
#' 
#' @param reads A read dataframe containing the columns "sequence_length_template" and "mean_qscore_template". 
#' @param outliers (Default = FALSE) If true, the longest 1% of reads will not be included as this can create skewing within the data. 
#' @param histogram (Default = FALSE) If true, ggExtra will be used to create a histogram of the x- & y-axis.
#' @param colors If ggplot is active in the namespace, the colors used to differentiate between pass and fail. If using default R, the colors will unfortunately be random.
#' @param theme (Default = theme_classic()) When using ggplot2, what theme should the plot be formatted in.
#' @return A subset of the input duplex dataframe in which a read is used as both a template and a complement. Useful for finding false positive duplex reads
#' @import ggplot2
#' @import ggExtra
#' @examples 
#' library2_duplex_table <- dorado_table(reads.df, duplex.df)
#' all.libraries.duplex.table <- as.data.frame(lapply(reads.list, dorado_table))
#' @export
duplex_differences_graph <- function(duplex, simplex = NULL, difference = "ALL", theme = ggplot2::theme_classic()) {
  stopifnot("The differences parameter must be ALL, TEMPLATE, COMPLEMENT" = difference %in% c("ALL", "TEMPLATE", "COMPLEMENT"))
  stopifnot("The duplex parameter is not a dataframe" = is.data.frame(duplex))
  if (c("duplex_sequence_length", "template_length", "complement_length") %in% colnames(duplex)){
    stopifnot("The Duplex dataframe does not have information about the simplex parents. Please run either duplex_parents() or include the simplex dataframe in duplex_differences_graph()" = !is.null(simplex))
    
    duplex <- duplex_parents(duplex, simplex)
  }
  
  if (difference == "ALL"){
    duplex[,"template_differences"] <- -(duplex$template_length - duplex$duplex_sequence_length)
    duplex[,"complement_differences"] <- (duplex$complement_length - duplex$duplex_sequence_length)
  } else if (difference == "TEMPLATE") {
    duplex[,"template_differences"] <- -(duplex$template_length - duplex$duplex_sequence_length)
  } else {
    duplex[,"complement_differences"] <- (duplex$complement_length - duplex$duplex_sequence_length)
  }
  
  
  
  
  
}
