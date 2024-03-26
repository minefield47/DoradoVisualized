#' Find information about the simplex reads from which a duplex read was created.
#' @author Auden Block
#' 
#' @param duplex A dataframe of duplex reads
#' @param simplex A dataframe of simplex reads from which the duplex reads were generated
#' @param mutate (Default = True) Mutate an existing dataframe to where the duplex parents are in the same row as the duplex read information, or return a dataframe only containing the simplex information. 
#' @return Either the duplex dataframe with each duplex read now having the length and quality of its simplex parents (mutate = true) or a new dataframe of only the duplex parents. 
#' 
#' @examples 
#' duplex.with.parents.df <- duplex_parents(duplex.df, simplex.df);
#' summary.list$duplex <- duplex_parents(summary.list$duplex, summary.list$simplex);
#' duplex.parents.only.df <- duplex_parents(duplex.df, simplex.df, mutate = FALSE)
#' @export
duplex_parents <- function(duplex, simplex, mutate = TRUE) {
  #Value checkers: 
    #Dataframes
    stopifnot("The duplex parameter must be a dataframe." = is.data.frame(duplex))
    stopifnot("The simplex parameter must be a dataframe." = is.data.frame(simplex))
    #Mutate
    stopifnot("Mutate must be a boolean (T/F) value" = is.logical(mutate))
    
  
  if (mutate == TRUE){
    #Create two new columns with names "template_id" and "complement_id" by splitting the duplex read name into the respective simplex parents. 
    duplex[c("template_id", "complement_id")] <- do.call(rbind, strsplit(duplex$read_id, ";"))
    
    #Differentiation of Duplex Length and Q-Score.
    #For comparisons later, need the column names to be different. 
    
    #Change the name of the duplex sequence length column to duplex_length.
    names(duplex)[names(duplex) == 'sequence_length_template'] <- 'duplex_sequence_length'
    #Change the name of Q-Score. 
    names(duplex)[names(duplex) == 'mean_qscore_template'] <- 'duplex_mean_qscore'
    
    
    #Pull the matching reads from the simplex dataframe that match the complement column,
    #remove all but the read id, sequence length, and mean q score, 
    #and merge into the duplex dataframe. 
    duplex <- merge(x=duplex, simplex[simplex$read_id %in% duplex$template_id,c("read_id", "sequence_length_template", "mean_qscore_template")], by.x = "template_id", by.y= "read_id")
    
    #Piped Version here for read-ability: 
    #
    #duplex <- simplex[simplex$read_id %in% duplex$template_id,] %>% 
    #   .[c("read_id", "sequence_length_template", "mean_qscore_template")] %>% 
    #   merge(duplex, ., by.x = "template_id", by.y= "read_id")
    
    
    #Change the name of the sequence_length_template column to be the template length
    names(duplex)[names(duplex) == 'sequence_length_template'] <- 'template_length'
    #Change the name of Q-Score. 
    names(duplex)[names(duplex) == 'mean_qscore_template'] <- 'template_mean_qscore'
    
    
    #Pull the matching reads from the simplex dataframe that match the complement column,
    #remove all but the read id, sequence length, and mean q score, 
    #and merge into the duplex dataframe. 
    duplex <- merge(x=duplex, simplex[simplex$read_id %in% duplex$complement_id,c("read_id", "sequence_length_template", "mean_qscore_template")], by.x = "complement_id", by.y= "read_id")
    
    
    # Piped Version here for readability:
    #   duplex <- simplex[simplex$read_id %in% duplex$complement_id,] %>% 
    #   .[c("read_id", "sequence_length_template", "mean_qscore_template")] %>%
    #   merge(duplex, ., by.x = "complement_id", by.y= "read_id")
    
    
    #Change the name of the sequence_length_template column to be the complement length
    names(duplex)[names(duplex) == 'sequence_length_template'] <- 'complement_length'
    #Change the name of Q-Score. 
    names(duplex)[names(duplex) == 'mean_qscore_template'] <- 'complement_mean_qscore'
    
    return(duplex)
  } else {
  #if mutate ain't true...it false. We just want to 
    
    return(
      subset(simplex, (simplex$read_id %in% as.character(do.call(rbind, strsplit(duplex$read_id, ";")))))
      )
  } 
}
