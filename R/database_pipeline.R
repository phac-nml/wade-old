#' Generic Profiler : SampleNo_contigs.fasta vs. arg-annot/resfinder2/CARD or VFDB
#'
#' @param org_id Organism to query: GAS or GONO
#' @param samples.df Data frame of user selected samples and the sample paths
#' @param is_vfdb Boolean, determine which database is being used. Specified from main call
#'
#' @details
#' GAS AMR loci
#' dfrF, dfrG, DHFR, ermA, ermB, ermT, folA, folP, gyrA, mefAE, msrD, parC, rpsJ, tetM, tetO
#' note that ermA == ermTR, ermC == ermT, mefA/E includes mefA and mefE; DHFR is the same as folA
#'
#' ARG-ANNOT is a downloaded as a single multifasta
#'
#' RESFINDER is separated into antimicrobial classes
#'
#' VFDB is a downloaded from the web site as a single multifasta
#'
#' @return A table frame containing the results of the query
#' @export

database_pipeline <- function(org_id, samples.df, is_vfdb){
  # TODO ####
  # - Will the sbatch command be used for blastn?

  if(is_vfdb){ # Change DB based on is_vfdb
    databases <- "VFDB"
  } else {
    databases <- c("arg-annot", "resfinder2", "CARD")
  }

  # Initializing ####
  inc_amount <- 1/(length(databases)*nrow(samples.df)) # increment amount, for progress bar
  sample_files <- file.path(samples.df[,"parent_dir"], samples.df[,"subdir_id"], samples.df[,"filename"]) # List of all files and their path
  
  db_samples.df <- databases %>% map_dfr(~ data.frame("db" = .x, 
                                                      "sample" = sample_files, 
                                                      stringsAsFactors = FALSE))
  
  output.df <- map2_df(db_samples.df[,"db"], db_samples.df[,"sample"], ~ execute_blastout(.x, .y, inc_amount)) %>%
    select("SampleNo", "DataBase", "GeneID", "MatchID") # Just rearranging the columns here

  writeLines(paste("Writing output to", here("data", "DB_PIPELINE_OUT.tsv")))
  write.csv(x = output.df,
            file = here("data", "DB_PIPELINE_OUT.tsv"),
            row.names = FALSE)
  
  writeLines("DONE: database_pipeline() finished....")
  output.df
}

execute_blastout <- function(curr_db, sample, inc_amount){
  sample_num <- basename(sample) # Vector of sample names .fasta
  
  #incProgress(amount = inc_amount,
  #            message = paste("Blasting ", sample_num, " against ", curr_db, sep = ""))
  
  writeLines(paste("Executing blastn for:", sample_num, "on db", curr_db))
  headers <- c("SampleNo", "DataBase", "GeneID", "MatchID")
  
  if(curr_db == "VFDB"){
    blast_evalue <- "10e-50"
  } else {
    blast_evalue <- "10e-100"
  }
  
  # DATABASE ACCESS ####
  db_dir <- here("data", "databases", curr_db, paste(curr_db, ".fasta", sep = "")) # data/databases/curr_db/curr_db.fasta
  output_location <- here("data", "databases", curr_db, paste(curr_db, "_blast_out.tsv", sep = "")) # data/curr_db/curr_db.fasta
  
  #-------------
  # using "-outfmt 6" in the blast command provides us with an output table to read from. This there's no need to
  #   parse through a file, line by painful line.
  blast_command <- paste("blastn -db ", db_dir, " -query ",
                         sample, " -outfmt 6 -out ", output_location, " -evalue ", blast_evalue, sep = "")
  sys_command <- blast_command
  
  # If batch command is needed ####  
  # sbatch_command <- "sbatch -p NMLResearch -c 1 --mem=1G -J %u-database_pipeline-%J --wrap=" # Use a better job name
  # sys_command <- paste(sbatch_command, "'", blast_command, "'", sep = "")
  
  try(system(sys_command)) # Blast command call
  
  info <- file.info(output_location)
  
  if(!is.na(info) && info$size > 0){
    # Read the table to get both the GeneID and MatchID
    table <- read.table(file = output_location, stringsAsFactors = FALSE)
    curr_blast_table.df <- select(table, "V1", "V2")
    names(curr_blast_table.df) <- c("GeneID", "MatchID")
    
    # Add the SampleNo and DataBase
    curr_blast_table.df <- mutate(curr_blast_table.df,
                                  "SampleNo" = sample_num,
                                  "DataBase" = curr_db)
    
    curr_blast_table.df[, headers] # Order by headers
  } else {
    curr_blast_table.df <- data.frame()
  }
  curr_blast_table.df
}
