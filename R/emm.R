#' emm typing pipeline from WGS assemblies
#'
#' Takes Organism, Sample Number, Locus, and a Variable at queries a contig.fasta file
#' @param org_id Organism to query: GAS, PNEUMO or GONO
#' @param samples.df Data frame of user selected samples and the sample paths
#' @param locus Sample number associated with contig.fasta file
#' @importFrom dplyr filter bind_rows
#' @importFrom purrr is_empty map_df pmap
#' @import here
#' @import utils
#' @import magrittr
#' @return A table frame containing the results of the query
#' @export

emm <- function(org_id, samples.df, locus){
  
  # Directories ####
  
  lookup_dir <- system.file(paste("extdata/databases", org_id, "EMM", "allele_lkup_dna", sep = "/"), package = "wade") # extdata/databases/curr_db/curr_db.fasta
  all_loci <- system.file(paste("extdata/databases", org_id, "EMM", "temp", "loci.csv", sep = "/"), package = "wade") # Location of the loci
  blast_out_file <- paste(out_location, "emm_blast_out.csv", sep = "") # A temporary location to hold the blast output (Is deleted later)

  # ---------------- Variables ----------------
  error <- "Sample_Err"
  blast_names <- c("SampleNo", "Allele", "Ident", "Align",
                   "Mismatches", "Gaps", "SampleStart", "SampleEnd", 
                   "AlleleStart", "AlleleEnd", "eValue", "bit")

  #  ---------------- Initialize DF's ----------------
  blastout.df <- data.frame()
  
  output.df <- data.frame(matrix(ncol=5, nrow=0), stringsAsFactors = FALSE)
  names(output.df) <- c("Sample", "Type", "Subtype_Rep", "Subtype", "bp_id")

  # same as both AMR, and 23SrRNA files
  loci.df <- read.csv(all_loci) # Changed from read_csv
  if(locus != "list") { loci.df <- filter(loci.df, Locus_id == locus) }
  num_loci <- (dim(loci.df))[1]
  
  #  ---------------- Set up Sample list table ----------------
  num_samples <- (dim(samples.df))[1]
  
  all_loci <- loci.df %>% pmap(~ .x) # Get the loci as a list
  all_samples <- samples.df[,'filename']
  
  query_files <- file.path(samples.df[,"parent_dir"], samples.df[,"subdir_id"], samples.df[,"filename"])
  
  # progress_frac <- 1/(length(all_samples)*length(all_loci)) # Determine the incProgress increment amount #progressrelated
  
  loci_dna_lookup <- all_loci %>% pmap(~ paste(lookup_dir, "/", .x, ".fasta", sep = ""))
  
  counter <- 1
  blastout.df <- query_files %>% map_df(function(x){ # Use the fastas on each loci and call emm_blastout
    if(file.exists(x)){
      new_blast <- data.frame()
      # Need to check if the loci.fasta exists
      # incProgress(amount = progress_frac,
      #            message = paste("Executing emm blast on ", basename(x), sep="")) #progressrelated
      
      # Perform and Read blast
      
      blast_command <- paste("blastn -query ", x, " -db ", loci_dna_lookup, " -out ", blast_out_file, " -num_alignments 10 -evalue 10e-50 -outfmt 6")
      info <- file.info(blast(blast_command, blast_out_file))
      
      if(is.na(info$size) || info$size == 0){ # There's nothing to be read!
        writeLines("Blast file not found!")
      } else { # There is file information so read it
        new_blast <- read.csv(blast_out_file,
                              header = FALSE,
                              sep = "\t",
                              stringsAsFactors = FALSE)
        
        names(new_blast) <- blast_names
        new_blast <- cbind(new_blast, "filename"=all_samples[counter])
        counter <<- counter + 1
      }
      file.remove(blast_out_file) # Remove the file after using it
      new_blast
    }
  })
  
  # Check if there is atleast SOME data in blastout.df
  #   If there is none then nothing happened! No data worked!
  if(purrr::is_empty(blastout.df)){
    return(data.frame("Error"="No data succeeded", stringsAsFactors = FALSE))
  }

  # Go through each row, do something and add it to output.df
  for(i in 1:nrow(blastout.df)){
    sub_type <- tolower(blastout.df[i, "Allele"])
    type <- sub_type
    bp <- (180 - blastout.df[i, "Mismatches"])
    type_rep <- get_emm_type_rep(bp, sub_type)
    
    temp.df <- data.frame("Sample" = blastout.df[i, "filename"],
                          "Type" = type,
                          "Subtype_Rep" = type_rep,
                          "Subtype" = sub_type,
                          "bp_id" = bp,
                          stringsAsFactors = FALSE)
    
    output.df <- dplyr::bind_rows(output.df, temp.df)
  }
  
  # A row will have an error when there is a sample that does not exist
  output.df <- filter(output.df, Type != error) # Filter out the rows that contain errors

  if(num_samples == 1){
    write_emm_output(write_blast = TRUE, blastout.df, output.df, org_id)
    return(blastout.df)
  } else {
    write_emm_output(write_blast = FALSE, blastout.df, output.df, org_id)
    return(output.df)
  }
} # end function call

# ------------------------------------
# write_emm_output()
#
# create .csv's based on the given parameters
write_emm_output <- function(write_blast, blast.df, sample.df, org_id, out=out_location){
  datetime <- format(Sys.time(), "%Y-%m-%d")
  
  # Multiple output files
  writeLines(paste("Writing output to", out))
  
  emm_blast_file <- paste(out, paste(Sys.Date(), org_id, "emmBLAST", "WADE.csv", sep = "_"), sep = "")
  emm_file <- paste(out, paste(Sys.Date(), org_id, "emm", "WADE.csv", sep = "_"), sep = "")
  emm_labware_file <- paste(out, paste(Sys.Date(), org_id, "emmLW", "WADE.csv", sep = "_"), sep = "")
  
  if(write_blast){ # We only write the blast when there is one file
    write.csv(blast.df, emm_blast_file, row.names = FALSE)
  } else { # Otherwise just write E V E R Y T H I N G
    write.csv(sample.df, emm_file, row.names = FALSE)
  }
  
  write.csv(sample.df, emm_labware_file, quote = FALSE, row.names = FALSE) # Always write a LabwareOutput
  
  writeLines("DONE: EMM_pipeline()")
}

# ------------------------------------
# get_emm_type_rep()
#
# Determine the emm_type_rep based on the bp_ids and the subtype
# Can bp_ids be above 180?
get_emm_type_rep <- function(bp_ids, subtype){
  if(bp_ids >= 180){
    type_rep <- subtype
  } else if(180 > bp_ids & bp_ids >= 92){
    type_rep <- "NT"
  } else {
    type_rep <- "New Type"
  }
  type_rep
}
