#' Import the summary tsv into R and subset by type
#' @author Auden Block
#' 
#' @description
#' Dorado summary creates a summary.tsv which can be read into R using read.delim(<file>, header=T,na.strings=c("","NA")). However, this can only read a single dataframe at a time and for duplex reads, is not subset into smaller dataframes. 
#' For a single file: The function imports the data using read.delim as a dataframe.
#' For a directory of reads: The function imports each .tsv file as a dataframe index of a list, then binds the list into a single dataframe. 
#' If Duplex reads are found: The imported dataframe is converted into a list of (1) duplex/simplex reads, (2) duplex reads, (3) simplex reads (including simplex reads that spawned duplex reads). 
#' If no duplex reads are found, the script returns a dataframe of all imported reads.
#' 
#' This is the SIMPLE version of importation. Designed for users who want to quickly subset their data into the three categories outlined above without additional subsetting of reads, creating either (1) additional overlapping of reads between dataframes, or (2) unnecessary divisions that hinder the user.
#' @param file_or_directory The path to a single duplex summary file or directory of duplex summary files. 
#' @return Depending upon input, either a list containing three dataframes: (1) Both duplex/simplex reads, (2) only duplex reads, (3) only simplex reads, OR if no duplex data is found, a single dataframe of simplex reads.
#' 
#' @examples 
#' summary.list <- import.dorado.tsv("/Users/user/Dorado/summary.tsv");
#' summary.list <- import.dorado.tsv("/Users/user/Dorado/summary_directory");
#' @export
import_dorado_tsv <- function(file_or_directory)  {
    if (dir.exists(file_or_directory)) { #Check if argument is a directory and that a DoradoQCduplex_<file_name> does not exist. If true, lump the directory and execute. 
      all_reads <- do.call("rbind", lapply(list.files(path=file_or_directory, full.names = TRUE),read.delim, header = T, na.strings =c("","NA")))
    } else {
      all_reads <-read.delim(file_or_directory, header = T, na.strings =c("","NA"))
    }
  
    if (any(is.na(all_reads$filename))) {
     #If filenames has NA values, this indicated duplex data, so we want to automatically return a list of the duplex and simplex data.   
       return(list(all_reads=all_reads,
            duplex=subset(all_reads, is.na(all_reads$filename)), #Subset of only duplex values since duplex do not have a filename value.
            all_simplex=subset(all_reads,!is.na(all_reads$filename))
                  )
            )
      }  else {
          return(all_reads)
        }
}  
