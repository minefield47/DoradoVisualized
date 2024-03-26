#' Dataframe containing a summary table about a Dorado run.
#' @author Auden Block
#' 
#' @param reads A read dataframe containing the columns "sequence_length_template" and "mean_qscore_template". 
#' @param qscore (Default = 10) What value is an acceptable q-score.
#' @param outliers (Default = FALSE) If true, the longest 1% of reads will not be included as this can create skewing within the data. 
#' @param histogram (Default = FALSE) If true, ggExtra will be used to create a histogram of the x- & y-axis.
#' @param colors If ggplot is active in the namespace, the colors used to differentiate between pass and fail. If using default R, the colors will unfortunately be random.
#' @param theme (Default = theme_classic()) When using ggplot2, what theme should the plot be formatted in.
#' @return A subset of the input duplex dataframe in which a read is used as both a template and a complement. Useful for finding false positive duplex reads
#' @import ggplot2
#' @import ggExtra
#' @examples 
#' qscore_vs_length_plot <- qscore_vs_readlength_graph(reads.df)
#' @export
qscore_vs_readlength_graph <- function(reads, qscore= 10, outliers = FALSE, histogram = FALSE, colors = c("#e76f51", "#2a9d8f"),  theme = ggplot2::theme_classic()) {
  stopifnot("Reads must be a dataframe"= is.data.frame(reads))
  stopifnot("The reads dataframe is missing the column names 'sequence_length_template' or 'mean_qscore_template'" = c("sequence_length_template", "mean_qscore_template") %in% colnames(reads))

  options(scipen = 999)

  if (isTRUE(outliers)) {
    #With outliers = true, we remove the longest 1% of reads from the dataframe to prevent visual skewing. 
    reads <- reads[reads$sequence_length_template > quantile(reads$sequence_length_template, prob = 1-(99/100)),]
  }
  

  if (isTRUE(histogram)) {
    stopifnot("ggExtra must be loaded to use the histogram function" = (requireNamespace("ggExtra", quietly = TRUE)))
    
    scatter.histo <- ggplot2::ggplot(reads, aes(x=.data$sequence_length_template, y=.data$mean_qscore_template, color=(.data$mean_qscore_template >= qscore))) + 
                     ggplot2::geom_point() +
                     ggplot2::labs(y="Mean Q-Score", x= "Read Length") +
                     ggplot2::scale_color_manual(name=paste0("Q-Score: ", qscore), values=colors, labels=c("Below Q-Score", "Above Q-Score"))
      
    return(ggExtra::ggMarginal(scatter.histo, theme + theme(legend.position = "left"), type="histogram", binwidth=1)) #Bin Length has to be 1 for the Quality Score index

    
  } else {
  if (requireNamespace("ggplot", quietly = TRUE)){
    return(ggplot2::ggplot(reads, aes(x=.data$sequence_length_template, y=.data$mean_qscore_template, color=(.data$mean_qscore_template >= qscore))) + 
             ggplot2::geom_point() +
             ggplot2::labs(y="Mean Q-Score", x= "Read Length") +
             ggplot2::scale_color_manual(name=paste0("Q-Score: ", qscore), values=colors, labels=c("Below Q-Score", "Above Q-Score")))
  } else{
    #Since ggplot2 isn't loaded. We will create a simple graph using default R
    
    #
    pass_fail <- ifelse(reads$mean_qscore_template >= qscore, TRUE, FALSE)
    reads[,"passfail"] <- pass_fail
    
    plot <- plot(reads$sequence_length_template, reads$mean_qscore_template, col=as.factor(reads$passfail), pch=19)
    plot <- legend("bottomright", legend = paste("Group", 1:2), col = 1:2, pch = 19, bty = "n")
    return(plot)
    }
  } 
}
