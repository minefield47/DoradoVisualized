#' Export the fastq file into R
#' @author Auden Block
#' 
#' This function is for those wishing to use R to manipulate and evaluate fastq files either through DuplexVisualized or other bioinformatics programs.
#' 
#' @param fastq.df The f
#' @return A fastq file of the fastq dataframe. 
#' 
#' @examples 
#' write_fastq(fastq.df,"<path_to_file>")
#' @export
write_fastq <- function(fastq.df, path){
  stopifnot("A fastq dataframe must be given" = (!is.null(fastq.df)))
  stopifnot("The fastq file must be a dataframe" = (is.data.frame(fastq.df)))
  stopifnot("A file path must be given"=(!is.null(path)))
  #Find the expected fastq length#Find the expected fastq length#Find the expected fastq length
  seq_num <- seq(1, (nrow(fastq.df)*4), by= 4)
  
  #Create a blank dataframe of the predicted length
  df <- data.frame(matrix(NA, ncol=1, nrow = seq_num))
  
  df[seq_num,] <- fastq.df[,"read.id.header"]
  df[seq_num+1,] <- fastq.df[,"sequence"]
  df[seq_num+2,] <- paste("+")
  df[seq_num+3,] <- fastq.df[,"quality"]
  
  #write the dataframe into a fastq.df file.
  write.table(df, path, row.names = FALSE, col.names = FALSE, quote = FALSE)
}