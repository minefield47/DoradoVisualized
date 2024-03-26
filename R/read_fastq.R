#' Import a fastq file into R
#' @author Auden Block
#' 
#' This function is for those wishing to use R to manipulate and evaluate fastq files either through DuplexVisualized or other bioinformatics programs.
#' 
#' @param fastq The path to the fastq file to be imported.
#' @return A dataframe containing the raw header line (ReadID, duplex value, etc), sequence, base quality scores, and a column containing the parsed ReadIDs.
#' 
#' @examples 
#' fastq <- read_fastq("<path_to_file>")
#' @export
read_fastq <- function(fastq, TAGS = TRUE){
  #Parsing the fastq file in with R's read.csv but instead of commas it is new lines. 
  fastq <- read.csv(fastq, sep = "\n", header = FALSE)
  
  #create a sequence number since we know all fastq's are composed of 4 lines. 
  seq_num <- seq(1, nrow(fastq), by= 4)
  
  #mutate the csv 
  fastq.final <- data.frame(
    read.id.header = fastq[seq_num,],
    sequence = fastq[(seq_num+1),],
    quality = fastq[(seq_num+3),]
  )
  fastq.final["read_id"] <- gsub("@", "", sapply(strsplit(as.character(fastq.final$read.id.header), "\\s"), "[[", 1))
  if (isTRUE(TAGS)){
    
  }
  

return(fastq.final)
}