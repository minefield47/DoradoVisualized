#' Subset __X coverage using Nanopore Reads
#' @author Auden Block
#' 
#' @param reads A dataframe or list of dataframes containing the column names: read_id and sequence_length_template.
#' @param coverage The desired coverage.
#' @param genome_size_mb Estimated genome size in megabases.
#' @param random This parameter only applies when a list of dataframes is given as input. When true, the script pulls randomly, which could bias against libraries with a low number of reads. If false, the script will try to pull uniformly from each element of the list (for instance a library) until no more reads can be pulled from the element. This option is advantageous where all elements have a sum sequence length greater than the desired coverage as small dataframes will always be present in a given coverage. 
#' @param ID (Default = TRUE) Return a character vector of readIDs or a dataframe of reads selected. 
#' @param tolerance (Default = 0.05) How leinient will you allow the algorithm to overshoot the given coverage. The lower the tolerance, the longer this function will take to run.
#' @param n_tries_limit (Default = 10) The number of tries the script will take to reach a coverage within the specified tolerance level.
#' @return A vector of subset read_ids/indexes whose length is equal to the desired coverage of a given genome.
#' 
#' @examples 
#' 2x.subset <- nanopore_subset(reads, 2, 1500, TRUE, 10)

#' @export
nanopore_subset <- function(reads, coverage, genome_size, units = "gb", random = TRUE, ID = TRUE, tolerance = 0.05, n_tries_limit = 10) {
  
  size <- switch(
    units, 
    "b" = 1,
    "kb" = 1e3,
    "mb"= 1e6,
    "gb"= 1e9,
  )
  
  
  total_bases <- coverage * genome_size*size #1MB = 1000000 (1mil) bases
  stopifnot("object reads is not in dataframe/list format" = is.list(reads))
  
  n_tries <- 0
  
  if (isTRUE(random)) {
    
    if (lapply(reads,is.data.frame) |> unlist() |> all()){
      #Each element of the reads list will be condensed to a dataframe.  
      names <- names(reads)
      
      reads<- do.call("rbind", mapply(cbind, reads, "list_ID"=names, SIMPLIFY=F))
      #We need to convert the list back into a dataframe
      convert_back <- TRUE
    } else {
      #We don't need to convert the list back into a dataframe as it should already be in dataframe format.
      convert_back <- FALSE
    }
    
    stopifnot("Input is not a data frame or list of dataframes" = is.data.frame(reads))
    stopifnot("Input read length is not enough to reach the requested coverage/genome_size." = (sum(as.numeric(reads$sequence_length_template)) > total_bases))
  
    #Find the unweighted average of the reads within each list, then find the average of those averages. 
    average_length <- mean(reads$sequence_length_template) |> unlist() |> mean() |> floor()
  
  repeat{
    #Get a sequence of 1:nreads#Get a sequence of 1:nrow of each element in the list. 
    element_index <- seq(1, nrow(reads), 1)
    
    read_length <- 0
    
    repeat{
      #How many bases do we need to 
      bases_needed <- total_bases - read_length
      
      #The number of reads we could select from each element of the list and only reach 80% (being conservative) of the remaining bases necessary. s
      n_to_select <- ((bases_needed*.8)/average_length) |> floor()
      
      if (n_to_select < 1){
        #Because we take the floor, anything less than 1 will return an n_to_select of zero...causing the loop to lock up.
        break
      }
      
      read_row_index <- sample(element_index, n_to_select)
      
      #Store the currently calculated read_length (current + new)
      read_length <- read_length + sum(reads$sequence_length_template[read_row_index])
      
      element_index <- element_index[! element_index %in% read_row_index]
    }
    
    #At this stage we should only need to iteratively sample one random read from a random element until we reach coverage. 
    while (read_length < total_bases) {
      #Within the selected element, select a random read row from the list of possible elements. 
      read_row_index <- sample(element_index, 1)
      
      #Store the currently calculated read_length (current + new)
      read_length <- read_length + reads$sequence_length_template[read_row_index]
      
      element_index <- element_index[! element_index %in% read_row_index]
    } 
    
    #Tolerance check: If we are above our tolerance: calculate again, else break the repeat loop 
    if (read_length < (total_bases + total_bases*tolerance)) {
      break
    }  
    
    
    if (n_tries == n_tries_limit){
      stop("Set Tolerance might be too stringent (Suggest decreasing tolerance) or n_tries_limit too low (Suggest increasing n_tries_limit if the desired tolerance is vital)")
    }
    n_tries <- n_tries + 1 
  }
  reads <- reads[-element_index,]
  
  
  if (isTRUE(convert_back)) { 
    #The input was a list of dataframes. We now want to give back a list of dataframes. 
    reads <- split(reads, reads$list_ID)
    #remove the list_ID parameter. 
    reads <- lapply(reads, function(x) {x[!(names(x) %in% c("list_ID"))]})
    if (isTRUE(ID)){
      return(lapply(reads, function(x) {x[(names(x) %in% c("read_id"))]}))
    } else {
      return(reads)
    }
  } else {
    #We do not need to convert the reads dataframe back into a list. It already is one. 
    if (isTRUE(ID)){
      return(reads$read_id[-element_index])
    } else {
      return(reads[-element_index])
    }
  } 
 } else {
    #We want to pull uniformly across a list of elements. 
    stopifnot("Elements of reads are not in dataframe format" = lapply(reads,is.data.frame) |> unlist() |> all())
    stopifnot("Sum read length is not enough for the given coverage/genome size." =  total_bases < (lapply(reads, function(i) {colSums(i['sequence_length_template'])}) |> unlist() |> sum()))
    #If a list is given, need to randomly subset from any element of list. 
    #If a list is given, need to randomly subset from any element of list. 
    #This is desired as it allows the user to retain which list element the reads generated came from 
    # without having to mutate the dataframe or search the selected readIDs from the input list. 
    
    
    repeat{
      #Everytime we repeat (i.e. we overshot higher than our tolerance allows), we need to redefine these values:
      element_index_list <- mapply(function(reads, nrows) {seq(1, nrows, 1)}, reads, nrows = lapply(reads, nrow))
      
      #Redefine our sequence and total length. 
      list_length_seq <- seq(1:length(reads))
      list_length <- length(reads)
      
      #Reset read_length to zero. 
      read_length <- 0
      
      rm_list_index <- FALSE
      #This determines the actual reads.
      repeat {
        #Find the unweighted average of the reads within each list, then find the average of those averages. 
        average_length <- lapply(seq_along(reads[list_length_seq]), function(x) {mean(reads[[x]]$sequence_length_template) + sd(reads[[x]]$sequence_length_template)}) |> unlist() |> mean() |> floor()
        #How many bases do we need to 
        
        
        #The number of reads we could select from each element of the list and only reach 80% (being conservative) of the remaining bases necessary. s
        n_to_select <- (((total_bases - read_length)*.8)/average_length/list_length) |> floor()
        
        if (n_to_select <= 1){
          #If we pull from each element, we would overshoot our conservative estimate, so break the repeat function. 
          break
        }
        
        
        if (lapply(seq_along(element_index_list[list_length_seq]), function(i, n_to_select) {isFALSE(n_to_select < length(element_index_list[[i]]))}, n_to_select=n_to_select) |> unlist() |> any()) {
          #We don't have enough reads in one of the elements to pull from.
          
          #Let's find the lowest nrow value, that will be our new n_to_select
          #Find the lowest_index
          list_index <- list_length_seq[lapply(reads[list_length_seq], nrow) |> which.min()]
          #Since which.min finds the lowest of the given seq, it doesn't return the actual ELEMENT POSITION within the list that is being called, just of the sequence which is being called. 
          #So for this, lets just use the NAME of the list instead of the vector value. 
          n_to_select <- nrow(reads[[list_index]][element_index_list[[list_index]],])
          rm_list_index <- TRUE
        } 
        
        
        for (i in list_length_seq) {
          #Within the selected element, select a random read row from the list of possible elements. 
          read_row_index <- sample(element_index_list[[i]], n_to_select)
          #Store the currently calculated read_length (current + new)
          read_length <- read_length + sum(reads[[i]]$sequence_length_template[read_row_index])
          
          element_index_list[[i]] <- element_index_list[[i]][! element_index_list[[i]] %in%  read_row_index]
          
        }
        
        if (isTRUE(rm_list_index)){
          #If the length of the element is zero, we have fully exhausted that element, so remove it from the our sequence of elements to pull from.  
          list_length_seq <- list_length_seq[! list_length_seq %in% list_index]
          list_length <- list_length - 1
          rm_list_index <- FALSE
          
          #However, we can still keep pulling from other samples to boost up towards are goal, so keep calling batchs. 
          # This prevents a file with only a couple of reads from significantly inhibiting read selection. Since you have an equal probibilty 
        }
        
      }
      
      
      
      #Since pulling one random read from each element would theoretically overshoot...sample one random read from a random element until we reach coverage. 
      while (read_length < total_bases) {
        #Select a random element of the list. 
        list_index <- sample(list_length_seq, 1)
        #Within the selected element, select a random read row from the list of possible elements. 
        read_row_index <- sample(element_index_list[[list_index]], 1)
        
        #Store the currently calculated read_length (current + new)
        read_length <- read_length + reads[[list_index]]$sequence_length_template[read_row_index]
        
        element_index_list[[list_index]] <- element_index_list[[list_index]][! element_index_list[[list_index]] %in%  read_row_index]
        
        if (length(element_index_list[[list_index]]) == 0){
          
          #If the length of the element is zero, we have fully exhausted that element, so remove it from the dataset. 
          list_length_seq <- list_length_seq[-list_index]
        }
        
      } 
      
      #Tolerance check: If we are above our tolerance: calculate again, else break the repeat loop 
      if (read_length < (total_bases + total_bases*tolerance)) {
        break
      }  
      
      if (n_tries == n_tries_limit){
        stop("Set Tolerance might be too stringent (Suggest decreasing tolerance) or n_tries_limit too low (Suggest increasing n_tries_limit if the desired tolerance is vital)")
      }
      n_tries <- n_tries + 1 
      print(paste("Trial:",n_tries))
    }
    
    #Return the values.
    if (isTRUE(ID)){
      return( mapply(function(reads,element_index_list, all_rows) {
        reads$read_id[! all_rows %in% element_index_list]
      }, reads, element_index_list, all_rows = mapply(function(reads, nrows) {seq(1, nrows, 1)}, reads, nrows = lapply(reads, nrow))) )
    } else {
      return(mapply(function(reads, element_index_list, all_rows) {
        reads[! all_rows %in% element_index_list,]
      }, reads, element_index_list, mapply(function(reads, nrows) {seq(1, nrows, 1)}, reads, nrows = lapply(reads, nrow)), SIMPLIFY = FALSE))
    }
  }
}

