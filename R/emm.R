#' emm typing pipeline from WGS assemblies
#'
#' Takes Organism, Sample Number, Locus, and a Variable at queries a contig.fasta file
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param samples.df Data frame of user selected samples and the sample paths
#' @param locus Sample number associated with contig.fasta file

#' @return A table frame containing the results of the query
#' @export

emm <- function(org_id, samples.df, locus){
  
  # Directories ####
  contigs_dir <- here("data", "databases", org_id, "assemblies")
  lookup_dir <- here("data", "databases", org_id, "EMM", "allele_lkup_dna")
  all_loci <- here("data", "databases", org_id, "EMM", "temp", "loci.csv")
  dest_file <- here("data", "output", "temp", "queryfile.fasta")
  blast_out_file <- here("data", "output", "temp", "blastout.csv")

  # ---------------- Variables ----------------
  Variable <- NA
  error <- "Sample_Err"
  blast_names <- c("SampleNo", "Allele", "Ident", "Align", "Mismatches", "Gaps", "SampleStart", "SampleEnd", "AlleleStart", "AlleleEnd", "eValue", "bit")

  #  ---------------- Initialize DF's ----------------
  sample_output.df <- data.frame()
  blastout.df <- data.frame()

  # same as both AMR, and 23SrRNA files
  loci.df <- read.csv(all_loci) # Changed from read_csv
  if(locus != "list") { loci.df <- filter(loci.df, Locus_id == locus) }
  num_loci <- (dim(loci.df))[1]

  #  ---------------- Set up Sample list table ----------------
  #samples.df <- get_samples(org_id, sample_num)
  num_samples <- (dim(samples.df))[1]

  all_loci <- loci.df %>% pmap(~ .x) # Get the loci as a list
  all_samples <- samples.df[,'filename']
  
  query_files <- file.path(samples.df[,"parent_dir"], samples.df[,"subdir_id"], samples.df[,"filename"])

  progress_frac <- 1/(length(all_samples)*length(all_loci))

  loci_dna_lookup <- all_loci %>% pmap(~ paste(lookup_dir, "/", .x, ".fasta", sep = ""))

  
  blast_info <- query_files %>% map(function(x){ # Use the fastas on each loci and call emm_blastout
    if(file.exists(x)){
      # Need to check if the loci.fasta exists
      incProgress(amount = progress_frac,
                  message = NULL)

      emm_blastout(x, loci_dna_lookup, blast_out_file) # How can we do this if there are multiple loci?
      info = file.info(blast_out_file)

      if(info$size == 0){
        writeLines("Blast file not found!")
      } else {
        blastout.df <<- read.csv(blast_out_file,
                                 header = FALSE,
                                 sep = "\t",
                                 stringsAsFactors = FALSE)
        names(blastout.df) <- blast_names
        blastout.df
      }
    }
  })

  output.df <- map2_dfr(blast_info, all_samples, function(x, y){
    if(!is.null(x)){
      sub_type <- tolower(x$Allele) # A list of alleles for y sample
      type <- sub_type
      bp <- (180 - x$Mismatches) # A list of mismatches for y sample
      type_rep <- unlist(map2(bp, sub_type, ~ get_emm_type_rep(.x, .y))) #

      data.frame(Sample = y,
                 Type = type,
                 Subtype_Rep = type_rep,
                 Subtype = sub_type,
                 bp_id = bp,
                 stringsAsFactors = FALSE)
    }
  })

  output.df <- filter(output.df, Type != error) # Filter out the rows that contain errors
  # A row will have an error when there is a sample that does not exist

  if(num_samples == 1){
    write_emm_output(write_blast = TRUE, blastout.df, sample_output.df, org_id)
    return(blastout.df)
  } else {
    write_emm_output(write_blast = FALSE, blastout.df, sample_output.df, org_id)
    return(output.df)
  }

} # end function call

# ------------------------------------
# write_emm_output()
#
# create .csv's based on the given parameters
write_emm_output <- function(write_blast, blast.df, sample.df, org_id){
  if(write_blast){
    write.csv(blast.df, here("data", "output", "output_profile_emm.csv"), row.names = FALSE)
  } else {
    write.csv(sample.df, here("data", "output", "output_profile_emm.csv"), quote = FALSE, row.names = FALSE)
  }

  write.csv(sample.df, here("data", "output", "LabWareUpload_GAS_emm.csv"), quote = FALSE, row.names = FALSE)

  writeLines("DONE: EMM_pipeline()")
}

# ------------------------------------
# emm_blastout()
#
# Perform a blastout command
emm_blastout <- function(dest, database, output){
  blast_command <- paste("blastn -query ", dest, " -db ", database, " -out ", output, " -num_alignments 10 -evalue 10e-50 -outfmt 6")
  try(system(blast_command))
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


