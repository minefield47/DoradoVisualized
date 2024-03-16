#' Find the simplex reads from which no duplex read was generated.
#' @author Auden Block
#' 
#' @param duplicated.df A dataframe of reads with which have been duplicated
#' @return A dataframe of the simplex reads from which NO duplex reads were generated.
#' 
#' @examples 
#' simplex.only <- simplex.only(duplex.df, simplex.df);
#' summary.list$simplex.only <- simplex.only(summary.list$duplex.df, summary.list$simplex.df);
#' @export