#' Find either the rolling duplex percentage over specific time intervals.
#' @author Auden Block
#' 
#' Calculate the rolling average duplex nucleotide or read rate over an experiment at given intervals.
#' 
#' @param duplex.df A dataframe of duplex reads
#' @param simplex.df A dataframe of simplex reads
#' @param interval (Default = 60) Interval in MINUTES to generate a rolling average. For instance, 60 will return the rate of duplex reads generated every 60 minutes of run time. 
#' @param type (Default = "nucleotide") The type of rate to be reported: nucleotide or read.
#' @return A dataframe of (,1:3) the interval evaluated reported in seconds, minutes, and hours, and (,4) the duplex rate as a decimal during that time interval. 
#' 
#' @examples 
#' rolling_duplex_rate.df <- rolling_duplex_rate(duplex.df, simplex.df, interval = 120)
#' rolling_duplex_rate.df <- rolling_duplex_rate(duplex.df, simplex.df, type = "read")
#' @export
rolling_duplex_rate <- function(duplex, simplex, interval = 60, type = "nucleotide") {
  #If either the duplex or simplex parameter given is not a dataframe, return an error. 
  stopifnot("The duplex parameter must be a dataframe." = is.data.frame(duplex))
  stopifnot("The simplex parameter must be a dataframe." = is.data.frame(simplex))
  
  #Interval being greater than zero
  stopifnot("Interval must be a number greater than zero" = interval >= 0 )
  
  #Type setting.
  stopifnot("Type must be either read or nucleotide" = type %in% c("read", "nucleotide"))
  
  if (type == "nucleotide") {
      seq_nums <- seq(
        from=0+(interval*60), 
        to=simplex$start_time[which.max(simplex$start_time)], 
        by = ((interval*60)*2)) 
      #For convenience, users input the interval in minutes, dorado is in seconds, so convert.
      
      return(as.data.frame(t(sapply(seq_nums, window_duplex_rate, interval=interval, duplex=duplex, simplex=simplex))))

  } else {
      
    seq_nums <- seq(from=which.min(simplex$template_start), to=which.max(simplex$template_start), by = (interval*60)*2) #For convenience, users input the interval in minutes, dorado is in seconds, so convert.
    
    
    return(as.data.frame(t(sapply(seq_nums, window_duplex_rate, interval = interval, duplex = duplex, simplex = simplex, type = "read"))))

      }
}

