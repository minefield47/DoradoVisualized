#' Find the duplex rate at a specific timepoint.
#' @author Auden Block
#' 
#' #' Calculate the average duplex nucleotide or read rate at a specific timepoint and window size.
#' 
#' @param duplex.df A dataframe of duplex reads
#' @param simplex.df A dataframe of simplex reads
#' @param interval (Default = 60) The window before and after the timepoint to analyze in minutes.
#' @param time_point The timepoint at which to analyze in seconds
#' @param type (Default = "nucleotide") The type of rate to be reported: nucleotide or read.
#' @return A dataframe of (,1:3) the interval evaluated reported in seconds, minutes, and hours, and (,4) the duplex rate as a decimal during that time interval. 
#' 
#' @examples 
#' window_duplex_rate(duplex.df, simplex.df, 3600, 60, type = "read")
#' @export

window_duplex_rate <- function(duplex, simplex, time_point, interval = 60, type = "nucleotide") {
  #If either the duplex or simplex parameter given is not a dataframe, return an error. 
  stopifnot("The duplex parameter must be a dataframe." = is.data.frame(duplex))
  stopifnot("The simplex parameter must be a dataframe." = is.data.frame(simplex))
  
  #Interval being greater than zero
  stopifnot("Interval must be a number greater than zero" = interval >= 0 )
  
  #Type setting.
  stopifnot("Type must be either read or nucleotide" = type %in% c("read", "nucleotide"))

  
  window.min <- time_point - (interval*60)
  window.max <- time_point + (interval*60)
  
  stopifnot("Timepoint minus interval must be a number greater than zero" = window.min >= 0 )
  
  simplex_subset <- subset(simplex, simplex$start_time>window.min & simplex$start_time<window.max)
  duplex_subset <- subset(duplex, duplex$start_time>window.min & duplex$start_time<window.max)
  
  
  
  if (type == "nucleotide") {
    return(
      data.frame(
        point_seconds = time_point,
        interval_seconds = paste(window.min, "-", window.max),
        point_minutes = time_point / 60,
        interval_minutes = paste(round((window.min/60),3), "-", round((window.max/60),3)),
        point_hours = time_point / 60 / 60,
        interval_hours = paste(round((window.min/60/60),3), "-", round((window.max/60/60),3)),
        nucleotide_rate_decimal = round((sum(duplex_subset$sequence_length_template) * 2) / sum(simplex_subset$sequence_length_template),5),
        nucleotide_rate_percentage = round((sum(duplex_subset$sequence_length_template) * 2) / sum(simplex_subset$sequence_length_template) * 100,5),
        n_duplex = nrow(duplex_subset),
        n_simplex = nrow(simplex_subset)
      )
    )
  } else {
    return(
      data.frame(
        interval.seconds = paste(window.min, "-", window.max),
        interval.minutes = paste(round((window.min/60),3), "-", round((window.max/60),3)),
        interval.hours = paste(round((window.min/60/60),3), "-", round((window.max/60/60),3)),
        read.rate = round((nrow(duplex_subset) * 2) / nrow(simplex_subset),5), #As this is read rates, use a count of reads. 
        n.duplex = nrow(duplex_subset),
        n.simplex = nrow(simplex_subset)
      )
    )
  }
}

