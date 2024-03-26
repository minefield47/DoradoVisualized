#' Export the fastq file into R
#' @author Auden Block
#' 
#' This function is for those wishing to use R to manipulate and evaluate fastq files either through DuplexVisualized or other bioinformatics programs.
#' 
#' @param headers The Headers (with or without the \@ ) to write
#' @param sequences The bases from which to write
#' @param qualities The qualities 
#' @return A fastq file of the fastq dataframe. 
#' 
#' @examples 
#' write_fastq(fastq.df,"<path_to_file>")
#' @export
write_fastq <- function(headers, sequences, qualities, path = "./") {
  
  #Confirm the length of provided vectors are each identical. 
  stopifnot("The header, sequence, or header vector are of different lengths" = (sapply(list(length(sequences), length(qualities)), FUN = identical, length(headers)) |> all()))
    
  #Find the expected fastq length#Find the expected fastq length#Find the expected fastq length
  seq_num <- seq(1, (nrow(fastq.df)*4), by= 4)

  
  #Create a blank dataframe of the predicted length
  df <- data.frame(matrix(NA, ncol=1, nrow = seq_num))
  
  df[seq_num,] <- headers
  df[seq_num+1,] <- sequences
  df[seq_num+2,] <- paste("+")
  df[seq_num+3,] <- qualities
  
  #write the dataframe into a fastq.df file.
  write.table(df, path, row.names = FALSE, col.names = FALSE, quote = FALSE)
}
