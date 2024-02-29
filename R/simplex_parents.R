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