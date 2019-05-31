# Makeblastdb.R indexes Blast fasta databases

#' Index BLAST fasta datases of lookup data
#'
#' Takes Organism, Test, and Locus and create BLAST databases
#' @param org_id Organism to query: GAS, PNEUMO or GONO
#' @param TestID Sample number associated with contig.fasta file
#' @param locus_id Locus_name to build blast database
#' @return A table frame containing the results of the query
#' @export

Index_pipeline <- function(org_id, test_id, locus_id){
  lookup_dir <- here("data", "databases", org_id, test_id)

  if(locus_id == "list") {
    loci.df <- read.csv(file = paste(lookup_dir, "/temp/loci.csv", sep = ""),
                        header = TRUE,
                        stringsAsFactors = FALSE,
                        sep = ",")
  } else {
    loci.df <- data.frame(locus_id)
  }

  # Iteration ####
  for(x in 1:(dim(loci.df))[1])
  {
    locus <- as.character(loci.df[x,1])
    locus_dna_lookup <- paste(lookup_dir, "/allele_lkup_dna/", locus, ".fasta", sep = "")

    if(file.exists(locus_dna_lookup)) {
      blast_command <- paste("makeblastdb -in ", locus_dna_lookup, " -dbtype nucl", sep = "")
      try(system(blast_command))
    }
  }

  # Output ####
  writeLines("DONE:  BLAST databases indexed!")
  done_signal.df <- data.frame(Output = "Blast databases indexed.", stringsAsFactors = FALSE)

  # Return ####
  done_signal.df
}
