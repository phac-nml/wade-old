#' rRNA.R
#'
#' @param org_id Organism to query: GAS, PNEUMO or GONO
#' @param sample_num Sample number or "list" or "folder" of sample numbers associated with VCF file(s)
#' @return A table frame containing the results of the query
#' @details
#' 23S rRNA pipeline for WGS assemblies to determine number of mutated alleles
#'
#' E.coli: A2059G and C2611T
#'
#' GONO:
#' Run SNP core pipeline with NCCP11945_23S4.fasta file: A2045G or C2597T

#' PNEUMO:
#' Run SNP core pipeline with 23S_R6.fasta file: A2061G or C2613T
#'
#' Takes Organism, Sample Number and queries a contig.fasta
#' #Parses 23s rRNA mutations from VCF files
#'
#' ON GALAXY:
#'
#' Alternative_allele_proporition = 0.1
#' min_coverage = 15
#' min_mean_mapping = 30
#' run_name = 23S
#' in FreeBayes step of phylogeny set ploidy = 4
#'
#' Export VCF from Galaxy:
#' Under "Collection Tools" , use the tool "Export to Warehouse" on "Filter vcf on collection 210".
#' It will bundle them together and makes a history item.  They are in /Warehouse/Temporary/galaxy_exports/RANDOM_NAME.
#' Click the "eye" to view RANDOM_NAME. You can then copy them from there using the regular file explorer into:
#'             " folder
#'
#' @export

rna_23s <- function(org_id, sample_num){

  # ------------------ Organism Variables ------------------
  rna_dir <- here("data/databases", org_id, "23S_rRNA")
  vcf_folder <- paste(rna_dir, "/VCF", sep = "")

  switch(org_id,
         GONO = {
           chrom <- "23S4_NCCP11945"

           rRNA23S_pos1.df <- data.frame(CHROM = chrom,
                                         POS = 2045,
                                         stringsAsFactors = FALSE)
           rRNA23S_pos2.df <- data.frame(CHROM = chrom,
                                         POS = 2597,
                                         stringsAsFactors = FALSE)
         },
         PNEUMO = {
           chrom <- "23S_rRNA_R6_sprr02"

           rRNA23S_pos1.df <- data.frame(CHROM = chrom,
                                         POS = 2061)
           rRNA23S_pos2.df <- data.frame(CHROM = chrom,
                                         POS = 2613)
         }
  )

  # ------------------ Setting up Data Frames ------------------
  if(sample_num == "list") { # Get samples can only
    samples.df <- read.csv(paste(rna_dir, "/list.csv", sep = ""))
  } else {
    samples.df <- data.frame(SampleNo = sample_num, Variable = NA)
  }

  num_samples <- (dim(samples.df))[1] # Number of samples from samples.df

  writeLines("\n23s rRNA counts from VCF") # Output ####

  # Get a list of each sample.vcf
  sample_list <- map(1:num_samples, ~ paste(samples.df[.x, "SampleNo"], ".vcf", sep = "")) # list of all the samples files
  sample_paths <- map(sample_list, ~ paste(vcf_folder,"/", .x, sep = "")) # list of all the samples file paths

  vcf_files.df <- map_df(sample_paths, ~ read_vcf(.x)) # Give path, get a data frame

  # ------------------ Query Data Frames ------------------
  values <- map(1:num_samples, function(x){ # Go through the rows of A SINGLE df
    frame <- vcf_files.df[x,]

    if(!is_empty(frame)){
      a20_values <- get_allele_fraction(curr_vcf_frame = frame, # Get the allele fractions for a20
                                        rna_pos = rRNA23S_pos1.df)

      c26_values <- get_allele_fraction(curr_vcf_frame = frame, # and now for c26
                                        rna_pos = rRNA23S_pos2.df)

      sample <- frame[length(frame)] # Grab the sample
      data.frame(SampleNo = sample, A2059G = a20_values, C2611T = c26_values)
    } else {
      data.frame() # Return an empty data frame into it
    }
  })

  # ------------------ Output ------------------
  output_profile.df <- data.frame() # Initialize Return Var ####
  output_profile.df <- bind_rows(output_profile.df, values)

  #write_csv(output_profile.df, here("data", "output", "output_profile_23S.csv"))
  writeLines("rRNA23S_pipeline() Completed...")

  output_profile.df
  #filter(output_profile.df, !(A2059G == 0 && C2611T == 0)) # Remove all 0, 0 columns
}

# ------------------ get_allele_fraction() ------------------
# Determines number of mutated alleles
#
# Called from rRNA23S_pipeline()
#
# PARAM(S):
#           curr_vcf_frame - a data frame from a .vcf file
#           rna_pos - a data frame of the CHROM and POS
#
# Return: an integer value
#--------------------------
get_allele_fraction <- function(curr_vcf_frame, rna_pos){
  val <- 0L

  # are the chromosomes the same?
  if(rna_pos[1] == curr_vcf_frame[1]){  #2597 & 2045 for GONO; 2061 for pneumo

    if(rna_pos[2] == curr_vcf_frame[2]){ # are the pos's the same?

      unknown_col <- curr_vcf_frame[10]
      vcf_parts <- str_split(unknown_col, ":")
      vcf_parts <- unlist(vcf_parts)

      depth <- as.integer(vcf_parts[2])   #total number of reads (depth of reads) (DP)
      alt_obs <- as.integer(vcf_parts[5])  #alternate observations from reference (AO)

      allele_fraction <- as.integer((alt_obs/depth)*100) # Get the fraction alt_obs/depth

      # Determine a value based on that fraction
      if(allele_fraction < 15) {
        val <- 0L
      } else if ((15 <= allele_fraction) && (allele_fraction < 35)) {
        val <- 1L
      } else if ((35 <= allele_fraction) && (allele_fraction < 65)) {
        val <- 2L
      } else if ((65 <= allele_fraction) && (allele_fraction < 85)) {
        val <- 3L
      } else if (85 <= allele_fraction) {
        val <- 4L
      }
    }
  }

  val
}
