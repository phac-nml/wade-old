#' database_pipeline.R
#' Generic Profiler : SampleNo_contigs.fasta vs. arg-annot/resfinder2/CARD or VFDB
#'
#' @param org_id Organism to query: GAS or GONO
#' @param sample_num Sample number (or list of sample numbers) associated with contig.fasta file
#' @param is_vfdb Boolean, determine which database is being used. Specified from main call
#'
#' @details
#' GAS AMR loci
#' dfrF, dfrG, DHFR, ermA, ermB, ermT, folA, folP, gyrA, mefAE, msrD, parC, rpsJ, tetM, tetO
#' note that ermA == ermTR, ermC == ermT, mefA/E includes mefA and mefE; DHFR is the same as folA
#'
#' ARG-ANNOT is a downloaded as a single multifasta. When a new download is made, on the Windows
#' command prompt, naviate to :
#' "L:\ GC_PAHO\ Whole_Genome_Sequencing\ WGS_Typing\ Resistance_Genes\ arg_annot"
#' Then make BLAST index by running from command prompt: "makeblastdb -in arg-annot.fasta -dbtype nucl"
#'
#' RESFINDER is separated into antimicrobial classes.  Use Windows command prompt to navigate to
#' "L:\ GC_PAHO\ Whole_Genome_Sequencing\ WGS_Typing\ Resistance_Genes\ resfinder2" folder, then run:
#' "copy /b *.fsa resfinder2.fasta" to concatenate all into one multifasta.
#' Then make BLAST index by running from command prompt: "makeblastdb -in resfinder2.fasta -dbtype nucl"
#'
#' VFDB is a downloaded from the web site as a single multifasta. When a new download is made, on the Windows
#' command prompt, naviate to :
#' "L:\ GC_PAHO\ Whole_Genome_Sequencing\ WGS_Typing\ VFDB \ "
#' Then make BLAST index by running from command prompt: "makeblastdb -in VFDB.fasta -dbtype nucl"
#'
#' @return A table frame containing the results of the query
#' @export

database_pipeline <- function(org_id, sample_num, is_vfdb){
  # TODO ####
  # - Error Checking
  #     - Does the sample entered exist?
  #       - We don't want it to crash because a bad sample
  #       - What should be output then?
  #
  # - Will the sbatch command be used for blastn?
  #
  # - When there is only one sample IncProgress is useless, should go as it executes the blast commands

  if(is_vfdb){ # Change DB based on is_vfdb
    databases <- "VFDB"
  } else {
    databases <- c("arg-annot", "resfinder2", "CARD")
  }

  # Initialize Return DF ####
  results.df <- data.frame()

  # Directories ####
  contigs_dir <- here("data", "databases", org_id, "assemblies")

  sample_numbers <- get_samples(org_id, sample_num) %>% select("SampleNo") # List of the samples
  inc_amount <- 1/(length(databases)*nrow(sample_numbers)) # Amount to increment, based on the number of databases and samples

  sample_files <- sample_numbers %>%
    map_df(~ paste(contigs_dir, "/", .x, ".fasta", sep = "")) %>%
    filter(file.exists(SampleNo))
  names(sample_files) <- "path"

  output.df <- unlist(sample_files) %>% map_dfr(function(x){
    curr_sample = x
    out <- databases %>% map_dfr(~ execute_blastout(curr_db = .x,
                                          sample = curr_sample))
  })

  output.df <- select(output.df, "SampleNo", "DataBase", "GeneID", "MatchID")
  print(output.df)

  writeLines(paste("Writing output to", here("data", "DB_PIPELINE_OUT.tsv")))
  write.csv(x = output.df,
            file = here("data", "DB_PIPELINE_OUT.tsv"),
            row.names = FALSE)

  writeLines("DONE: database_pipeline() finished....")
  output.df
}

#-------------
# using "-outfmt 6" in the blast command provides us with an output table to read from. This there's no need to
#   parse through a file, line by painful line.
execute_blastout <- function(curr_db, sample){#, inc_amount){
  sample_num <- basename(sample) # Vector of sample names .fasta

  # incProgress(amount = inc_amount,
  #             message = paste(sample_num, " against ", curr_db, sep = ""))

  writeLines(paste("Executing blastn for:", sample_num))
  headers <- c("SampleNo", "DataBase", "GeneID", "MatchID")

  if(curr_db %in% "VFDB"){
    blast_evalue <- "10e-50"
  } else {
    blast_evalue <- "10e-100"
  }

  db_dir <- here("data", "databases", curr_db, paste(curr_db, ".fasta", sep = "")) # data/databases/curr_db/curr_db.fasta
  output_location <- here("data", "databases", curr_db, paste(curr_db, "_blast_out.tsv", sep = "")) # data/curr_db/curr_db.fasta

  #sbatch_command <- "sbatch -p NMLResearch -c 1 --mem=1G -J %u-database_pipeline-%J --wrap=" # Use a better job name
  blast_command <- paste("blastn -db ", db_dir, " -query ",
                         sample, " -outfmt 6 -out ", output_location, " -evalue ", blast_evalue, sep = "")

  #sys_command <- paste(sbatch_command, "'", blast_command, "'", sep = "")
  sys_command <- blast_command
  try(system(sys_command)) # Blast command call

  # grab the file to make sure it's not empty
  info <- file.info(output_location)

  if(info$size > 0){
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
