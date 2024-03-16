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
dorado_table <- function(simplex, duplex = NA) {
  stopifnot("The simplex parameter must be a dataframe" = (is.data.frame(simplex)))
  stopifnot("The duplex parameter must be a dataframe or NA" = (is.data.frame(duplex) | (is.na(duplex))))
  if (is.data.frame(duplex)){
    data.frame(
      run_id = simplex$run_id[1],
      run_length_seconds = simplex$start_time[which.max(simplex$start_time)],
      n_active_channels = length(unique(simplex$channel)),
      n_active_pores = nrow(unique(simplex[,c("channel", "mux")])),
      perc_active_pores = nrow(unique(simplex[,c("channel", "mux")])) / 2048,
      min_simplex_length = simplex$sequence_length_template[which.min(simplex$sequence_length_template)], #Minimum Simplex
      max_simplex_length = simplex$sequence_length_template[which.max(simplex$sequence_length_template)], #Maximum Simplex
      N_50 = n_function(simplex$sequence_length_template),
      L_50 = l_function(simplex$sequence_length_template),
      N_75 = n_function(simplex$sequence_length_template, 75),
      L_75 = l_function(simplex$sequence_length_template, 75),
      N_90 = n_function(simplex$sequence_length_template, 90),
      L_90 = l_function(simplex$sequence_length_template, 90),
      n_simplex_reads = nrow(simplex),
      n_simplex_bases = sum(simplex$sequence_length_template),
      avg_q_score = mean(simplex$mean_qscore_template),
      max_q_score = simplex$mean_qscore_template[which.max(simplex$mean_qscore_template)],
      n_less_than_q10 = length(simplex$mean_qscore_template[which(simplex$mean_qscore_template<10)]),
      perc_less_than_equal_q10 = (length(simplex$mean_qscore_template[which(simplex$mean_qscore_template<10)])/nrow(simplex) * 100), #which returns a character vector, so calculate the length to get the number of rows which can be divided by the total number of rows in the dataframe_ 
      n_greater_than_equal_q10 = length(simplex$mean_qscore_template[which(simplex$mean_qscore_template>=10)]),
      perc_greater_than_equal_q10 = (length(simplex$mean_qscore_template[which(simplex$mean_qscore_template>=10)])/nrow(simplex) * 100),
      n_greater_than_equal_q20 = length(simplex$mean_qscore_template[which(simplex$mean_qscore_template>=20)]),
      perc_greater_than_equal_q20 = (length(simplex$mean_qscore_template[which(simplex$mean_qscore_template>=20)])/nrow(simplex) * 100),
      n_greater_than_equal_q30 = length(simplex$mean_qscore_template[which(simplex$mean_qscore_template>=30)]),
      perc_greater_than_equal_q30 = (length(simplex$mean_qscore_template[which(simplex$mean_qscore_template>=30)])/nrow(simplex) * 100),
      n_greater_than_equal_q40 = length(simplex$mean_qscore_template[which(simplex$mean_qscore_template>=40)]),
      perc_greater_than_equal_q40 = (length(simplex$mean_qscore_template[which(simplex$mean_qscore_template>=40)])/nrow(simplex) * 100),
      n_duplex_reads = nrow(duplex),
      n_duplex_bases = sum(duplex$sequence_length_template),
      duplex_rate_nucleotide = duplex_rate(duplex, simplex),
      duplex_average_qscore = mean(duplex$mean_qscore_template),
      duplex_max_qscore = duplex$mean_qscore_template[which.max(duplex$mean_qscore_template)],
      duplex_n_channels = length(unique(duplex$channel)),
      duplex_n_active_pores = nrow(unique(duplex[,c("channel", "mux")]))
    )
 
  } else {
     #No Duplex data provided. 
    data.frame(
      run_id = simplex$run_id[1],
      run_length_seconds = simplex$start_time[which.max(simplex$start_time)],
      n_active_channels = length(unique(simplex$channel)),
      n_active_pores = nrow(unique(simplex[,c("channel", "mux")])),
      perc_active_pores = nrow(unique(simplex[,c("channel", "mux")])) / 2048,
      min_simplex_length = simplex$sequence_length_template[which.min(simplex$sequence_length_template)], #Minimum Simplex
      max_simplex_length = simplex$sequence_length_template[which.max(simplex$sequence_length_template)], #Maximum Simplex
      N_50 = n_function(simplex$sequence_length_template),
      L_50 = l_function(simplex$sequence_length_template),
      N_75 = n_function(simplex$sequence_length_template, 75),
      L_75 = l_function(simplex$sequence_length_template, 75),
      N_90 = n_function(simplex$sequence_length_template, 90),
      L_90 = l_function(simplex$sequence_length_template, 90),
      n_reads = nrow(simplex),
      n_bases = sum(simplex$sequence_length_template),
      avg_q_score = mean(simplex$mean_qscore_template),
      max_q_score = simplex$mean_qscore_template[which.max(simplex$mean_qscore_template)],
      n_less_than_q10 = length(simplex$mean_qscore_template[which(simplex$mean_qscore_template<10)]),
      perc_less_than_equal_q10 = (length(simplex$mean_qscore_template[which(simplex$mean_qscore_template<10)])/nrow(simplex) * 100), #which returns a character vector, so calculate the length to get the number of rows which can be divided by the total number of rows in the dataframe_ 
      n_greater_than_equal_q10 = length(simplex$mean_qscore_template[which(simplex$mean_qscore_template>=10)]),
      perc_greater_than_equal_q10 = (length(simplex$mean_qscore_template[which(simplex$mean_qscore_template>=10)])/nrow(simplex) * 100),
      n_greater_than_equal_q20 = length(simplex$mean_qscore_template[which(simplex$mean_qscore_template>=20)]),
      perc_greater_than_equal_q20 = (length(simplex$mean_qscore_template[which(simplex$mean_qscore_template>=20)])/nrow(simplex) * 100),
      n_greater_than_equal_q30 = length(simplex$mean_qscore_template[which(simplex$mean_qscore_template>=30)]),
      perc_greater_than_equal_q30 = (length(simplex$mean_qscore_template[which(simplex$mean_qscore_template>=30)])/nrow(simplex) * 100),
      n_greater_than_equal_q40 = length(simplex$mean_qscore_template[which(simplex$mean_qscore_template>=40)]),
      perc_greater_than_equal_q40 = (length(simplex$mean_qscore_template[which(simplex$mean_qscore_template>=40)])/nrow(simplex) * 100)
    )
  
  }
}