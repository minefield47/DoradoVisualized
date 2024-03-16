#' Dataframe containing a summary table about a Dorado run.
#' @author Auden Block
#' 
#' @param simplex A simplex read dataframe.
#' @param duplex A duplex read dataframe from the same summary as the simplex dataframe. 
#' @return A subset of the input duplex dataframe in which a read is used as both a template and a complement. Useful for finding false positive duplex reads
#' 
#' @examples 
#' library2_duplex_table <- dorado_table(simplex.df, duplex.df)
#' all.libraries.duplex.table <- as.data.frame(lapply(simplex.list, dorado_table))
#' @export
qscore_vs_readlength_graph <- function(simplex, histogram = TRUE, colors = c("#e76f51", "#2a9d8f"), qscore= 10, theme = theme_classic()) {
  if (isTRUE(histogram)) {
    stopifnot("ggplot2 and ggExtra must be loaded to use this function" = (!requireNamespace(c("ggplot2", "ggExtra"), quietly = TRUE)))
    
    scatter.histo <- ggplot2::ggplot(simplex, aes(x=.data$sequence_length_template, y=.data$mean_qscore_template, color=(.data$mean_qscore_template >= qscore))) + 
                     ggplot2::geom_point() +
                     ggplot2::labs(y="Mean Q-Score", x= "Read Length") +
                     ggplot2::scale_color_manual(name=paste0("Q-Score: ", qscore), values=colors, labels=c("Below Q-Score", "Above Q-Score")) 
      
    return(ggExtra::ggMarginal(scatter.histo, theme_classic() + theme(legend.position = "left"), type="histogram", binwidth=1)) #Bin Length has to be 1 for the Quality Score index

    
  } else {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::ggplot(simplex)
  }
  }
}

ggplot(bp.g1.simplex, aes(sequence_length_template, mean_qscore_template)) + geom_point() + theme
  