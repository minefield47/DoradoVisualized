#' Subset __X coverage using Nanopore Reads
#' @author Auden Block
#' 
#' @param reads A dataframe containing the column names: read_id and sequence_length_template.
#' @param coverage The desired coverage.
#' @return A subset of the input duplex dataframe in which a read is used as both a template and a complement. Useful for finding false positive duplex reads
#' 
#' @examples 
#' library2_duplex_table <- dorado_table(simplex.df, duplex.df)
#' all.libraries.duplex.table <- as.data.frame(lapply(simplex.list, dorado_table))
#' @export