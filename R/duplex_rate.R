#' Find either the total duplex percentage over a run.
#' @author Auden Block
#' 
#' Dorado will automatically report the duplex rate over the entire run time as a nucleotide rate at the end of basecaling. However, this can be difficult to utilize within R or recalculate if a standard output file is lost. We also included the ability to evaluate the number of duplex reads generated over the experiment.
#' 
#' @param duplex.df A dataframe of duplex reads
#' @param simplex.df A dataframe of simplex reads
#' @param type (default = "nucleotide") The type of rate to be reported: nucleotide or read.
#' @return A dataframe of the run time evaluated in seconds, minutes, and hours, and the duplex rate as a decimal value. 
#' 
#' @examples 
#' duplex_rate.df <- duplex_rate(duplex.df, simplex.df)
#' duplex_rate.df <- duplex_rate(duplex.df, simplex.df, type = "read")
#' @export
duplex_rate <- function(duplex, simplex, type = "nucleotide") {
  #If either the duplex or simplex parameter given is not a dataframe, return an error. 
  stopifnot("The duplex parameter must be a dataframe." = is.data.frame(duplex))
  stopifnot("The simplex parameter must be a dataframe." = is.data.frame(simplex))

  #Type setting.
  stopifnot("Type must be either read or nucleotide" = type %in% c("read", "nucleotide"))
  
  
  if (type == "nucleotide") {
    #This is for the nucleotide version
  
      return(
        data.frame(
          nucleotide.rate = round((sum(duplex$sequence_length_template) * 2) / sum(simplex$sequence_length_template),5),
          n.duplex = nrow(duplex),
          n.simplex = nrow(simplex)
        )
      )
      #Duplex * 2 because each duplex read represents two simplex reads. 

  } else {
      return(
        data.frame(
          read.rate = round((nrow(duplex$sequence_length_template) * 2) / nrow(simplex$sequence_length_template),5),
          n.duplex = nrow(duplex),
          n.simplex = nrow(simplex)
        )
      )
      #Duplex * 2 because each duplex read represents two simplex reads. 

  }
}

