#' Find either the total duplex READ percentage or duplex percentage over specific time intervals.
#' @author Auden Block
#' 
#' Dorado will automatically report the duplex rate over the entire run time as a nucleotide rate into the standard output (std.out). However, this can be difficult to utilize within R or if the std.out files are lost. In addition to a nucleotide rate duplex.rate can output the duplex read rate to evaluate the number of duplex reads generated.
#' 
#' @param duplex.df A dataframe of duplex reads
#' @param simplex.df A dataframe of simplex reads
#' @param interval (Default = 0, the entire experimental run) Interval divisions in MINUTES to report the rate as a decimal value. For instance, 60 will return the rate of duplex reads generated every 60 minutes of run time. 
#' @param type (default = "nucleotide") The type of rate to be reported: nucleotide or read.
#' @return A dataframe of (,1:3) the interval evaluated reported in seconds, minutes, and hours, and (,4) the duplex rate as a decimal during that time interval. 
#' 
#' @examples 
#' simplex.parents <- simplex.parents(duplex.df, simplex.df);
#' summary.list$simplex.parents <- simplex.parents(summary.list$duplex.df, summary.list$simplex.df);
#' @export

duplex.rate <- function(duplex, simplex, interval = 0, type = "nucleotide") {

  #If either the duplex or simplex parameter given is not a dataframe, return an error. 
  stopifnot("The duplex parameter must be a dataframe." = is.data.frame(duplex))
  stopifnot("The simplex parameter must be a dataframe." = is.data.frame(simplex))
  
  #Interval being greater than zero
  stopifnot("Interval must be a number greater than zero" = interval >= 0 )
  
  #Type setting.
  stopifnot("Type must be either read or nucleotide" = type %in% c("read", "nucleotide"))
  
  
  if (type == "read") {
    if (interval == 0) {
      #Since the interval is zero, no need to create divisions and try to analyze each, just use nrows
      
      return(
        data.frame(
          interval.seconds = paste(interval, "-", (which.max(simplex$template_start))),
          interval.minutes = paste(interval, "-", round((which.max(simplex$template_start)/60),3)),
          interval.hours = paste(interval, "-", round((which.max(simplex$template_start)/60/60),3)),
          read.rate = (nrow(duplex) * 2) / nrow(simplex)
          )
        )
      #Duplex * 2 because each duplex read represents two simplex reads. 
    } else {
      i <- seq(from=which.min(simplex$template_start), to=which.max(simplex$template_start), by = (interval*60)) #For convenience, users input the interval in minutes, dorado is in seconds, so convert.
    
      
      
      return()      
      
      
    }
    
  } else {
    if (interval == 0) {
      #Since the interval is zero, no need to create divisions and try to analyze each, just use nrows and then find the maximum start time (idicating the end of the run)
      return(
        data.frame(
          interval.seconds = paste(interval, "-", which.max(simplex$template_start)),
          interval.minutes = paste(interval, "-", round((which.max(simplex$template_start)/60),3)),
          interval.hours = paste(interval, "-", round((which.max(simplex$template_start)/60/60),3)),
          nucleotide.rate = (sum(duplex$sequence_length_template) * 2) / sum(simplex$sequence_length_template)
        )
      )
      #Duplex * 2 because each duplex read represents two simplex reads. 
    } else {
      
    }
  }
}