#' Dataframe containing a summary table about a Dorado run.
#' @author Auden Block
#' 
#' @param intervals The intervals to graph (such as rates$points_hours from `rolling_duplex_rate()`)
#' @param rates The calculated rate percentage for each interval time point (such as rates$nucleotide_rate_percentage from `rolling_duplex_rate()`)
#' @param units (Default = c("Hours", "%", "nucleotides per hour")) The units used for the interval and rate.
#' @param theme (Default = theme_classic()) When using ggplot2, what theme should the plot be formatted in.
#' @return A graph
#' @import ggplot2
#' @examples 
#' rolling_rate_graph(rates$points_hours, rates$nucleotide_rate_percentage, c("Hours", "% nucleotides per hour"))
#' @export
rolling_rate_graph <- function(intervals, rates, units= c("Hours", "Duplex/Simplex Nucleotides per hour"), theme = ggplot2::theme_classic()) {
  
  rate_dataframe <- data.frame("interval" = unlist(intervals), "rates" = unlist(rates))
  return(ggplot2::ggplot(rate_dataframe, aes(x=.data$interval, y=.data$rates)) + 
    ggplot2::geom_line() +
    ggplot2::labs(x=paste0("Experiment Run (",units[1], ")"), y= paste0("Percent Duplex Rate (",units[2],")")) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) + 
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    theme)
}