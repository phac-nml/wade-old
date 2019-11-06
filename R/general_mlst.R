#' Create a profile for MLST, NGMAST, or NGSTAR using blastn
#'
#' @param org_id Organism to query: GAS, PNEUMO or GONO
#' @param samples.df Data frame of user selected samples and the sample paths
#' @param locus The locus to query, or enter "list" to use a list of alleles
#' @param seq_type The sequence type to use. Can be MLST, NGMAST, or NSTAR
#'
#' @details
#'
#' @importFrom dplyr filter arrange bind_cols bind_rows select
#' @import here
#' @import utils
#' @return A profile containing the sample number, and the loci tested on
#' @export

general_mlst_pipeline <- function(org_id, samples.df, locus, seq_type){
  if (seq_type %notin% c("MLST", "NGMAST", "NSTAR")) {
    stop('General MLST seq_type must be one of MLST, NGMAST, or NSTAR', call.=FALSE)
  }

  # ------------- Directories ----------------
  db_dir <- system.file("extdata/databases", package = "wade") # extdata/databases/curr_db/curr_db.fasta
  lookup_dir <- paste(db_dir, org_id, seq_type, "allele_lkup_dna", sep = "/")  # Lookups
  temp_dir <- paste(db_dir, org_id, seq_type, "temp", sep="/") # Temporary
  profiles_dir <- paste(temp_dir, "/", "profiles.csv", sep = "")
  loci.list <- paste(temp_dir, "/", "loci.csv", sep = "") # Loci list
  blast_out_file <- paste(out_location, org_id, "_blast_out.tsv", sep = "")
  
  
  if(!file.exists(blast_out_file)){
    file.create(blast_out_file) # If it doesn't exist yet just create it.
  }
  
  # ------------- Initialize Loci ----------------
  writeLines("Initialize Loci....")
  loci.df <- read.csv(file = loci.list, # read loci.csv
                      header = TRUE,
                      sep = ",",
                      stringsAsFactors = FALSE)
  
  if(locus != "list") { loci.df <- filter(loci.df, Locus_id == locus) }
  num_loci <- (dim(loci.df))[1]
  
  # ------------- Initialize Sample ---------------
  writeLines("Initialize Sample....")
  num_samples <- (dim(samples.df))[1]
  
  # ------------- Initialize Profiles ------------------
  writeLines("Initialize Profiles....")
  profiles.df <- read.csv(file = profiles_dir,
                          header = TRUE,
                          sep = ",",
                          stringsAsFactors = FALSE)
  
  # ------------- Initialize Other Variables -------------
  
  sample_out.df <- data.frame()
  sample_mut_out.df <- data.frame()
  
  # ------------ Loop Through Samples -----------------
  for(i in 1:num_samples){
    allele <- ""
    allele_num <- ""
    
    curr_sample_num <- sub("([^.]+)\\.[[:alnum:]]+$", "\\1", samples.df[i, "filename"])
    curr_file <- file.path(samples.df[i, "parent_dir"], samples.df[i, "subdir_id"], samples.df[i, "filename"])
    
    dest_file <- here("data", "queryfile.fasta") # Still need a finalized location
    if(!file.exists(dest_file)){
      file.create(dest_file) # If it doesn't exist yet just create it.
    }
    
    if(!file.copy(curr_file, dest_file, overwrite = TRUE)){
      allele <- "Sample Number Error" # Do we need this?
      allele_num <- "Sample_Err"
    }
    
    # ------------------ Loop Through Loci ----------------------
    for(k in 1:num_loci){
      curr_locus <- as.character(loci.df[k, "Locus_id"])
      curr_locus_length <- as.integer(loci.df[k, "size"])
      
      # Initialize these here to resolve scope issues
      mutations <- ""
      blast_out.df <- ""
      
      if(allele_num != "Sample_Err"){ # This is redundant!! Put it before so nothing can execute!
        
        #incProgress(amount = 1/(num_samples*num_loci),
        #            message = paste(seq_type, " Pipeline: ", curr_sample_num, " blastn on ", curr_locus, sep = "")) #progressrelated
        
        locus_dna_lookup <- paste(lookup_dir, "/", curr_locus, ".fasta", sep = "")
        
        if(!(file.exists(locus_dna_lookup))){ # The dna lookup doesn't exist
          writeLines(paste("Locus error for: ", locus, ".... file not found!"))
          sample_error.df <- data.frame(Locus_ID = locus,
                                        Output = "File Not Found",
                                        stringsAsFactors = FALSE)
          return(sample_error.df)
        }
        
        blast_command <- paste("blastn -db ", locus_dna_lookup, "-query ", dest_file, " -out ", blast_out_file, " -num_alignments 10 -evalue 10e-175 -outfmt 6")
        blast_info <- file.info(blast(blast_command, blast_out_file))
        
        if(blast_info$size == 0){
          allele <- "No gene present"
          allele_num <- "0"
          
          # ---> difference in NGSTAR here <---
          mutations <- "?"
          blast_out.df <- data.frame(curr_sample_num, curr_locus)
        } else {
          blast_out.df <- read.csv(file = blast_out_file,
                                   header = FALSE,
                                   sep = "\t",
                                   stringsAsFactors = FALSE)
          names(blast_out.df) <- c("SampleNo",
                                   "Allele",
                                   "Ident",
                                   "Align",
                                   "Mismatches",
                                   "Gaps",
                                   "SampleStart",
                                   "SampleEnd",
                                   "AlleleStart",
                                   "AlleleEnd",
                                   "eValue",
                                   "bit")
          
          if(seq_type == "NGSTAR"){
            # ---> difference in NGSTAR here <---
            blast_out_100.df <- filter(blast_out.df, Ident == 100 & Mismatches == 0 & Gaps == 0)
          } else {
            blast_out_100.df <- filter(blast_out.df,
                                       Ident == 100 & Align == curr_locus_length)
          }
          
          if(seq_type == "NGMAST"){
            # ---> addition in NGMAST here <---
            if(curr_locus == "porB"){
              str_start <- 5
            } else { str_start <- 6 }
            blast_out_100.df$Allele2 <- as.numeric(substr(blast_out_100.df$Allele, str_start, 10)) # add Allele2 as a column
            blast_out_100.df <- arrange(blast_out_100.df, Allele2)
          }
          
          blast_size <- nrow(blast_out_100.df) # ISSUE ####
          # If neither NGSTAR nor NGMAST are called blast_out_100.df is not created
          # Look into it
          
          if(blast_size > 0){
            allele <- blast_out_100.df[1, 2]
            allele_parts <- unlist(strsplit(allele, "_"))
            allele_num <- allele_parts[2]
            
            mutations <- allele_parts[3] # For NGSTAR, not used otherwise ####
          } else {
            allele <- "Not Found"
            allele_num <- "?"
            
            mutations <- "?" # For NGSTAR, not used otherwise
          }
        }
      } # end if != "Sample Err"
      
      #------------- If it's the first locus -------------
      if(k == 1){
        curr_output.df <- data.frame(curr_sample_num,
                                     allele_num,
                                     stringsAsFactors = FALSE)
        names(curr_output.df) <- c("SampleNo", curr_locus)
        
        # ---> addition in NGSTAR here <---
        curr_mut_output.df <- data.frame(curr_sample_num,
                                         mutations,
                                         stringsAsFactors = FALSE)
        names(curr_mut_output.df) <- c("SampleNo", curr_locus)
      } else {
        locus_output.df <- data.frame(allele_num, stringsAsFactors = FALSE)
        names(locus_output.df) <- c(curr_locus)
        curr_output.df <- bind_cols(curr_output.df, locus_output.df) # Bind a new column to add a new locus!
        
        # ---> addition in NGSTAR here <---
        locus_mut_output.df <- data.frame(mutations, stringsAsFactors = FALSE)
        names(locus_mut_output.df) <- c(curr_locus)
        
        curr_mut_output.df <- bind_cols(curr_mut_output.df, locus_mut_output.df)
        #curr_mut_output.df <- data.frame(curr_mut_output.df, locus_mut_output.df)
      }
    } # End Locus Loop ####
    
    #------ REMOVE FILES -------
    # We remove these files after our looping
    file.remove(blast_out_file) # Delete the blast file
    file.remove(dest_file) # Delete the copied .fasta
    
    if(locus == "list"){
      # ----------- lookup profiles
      profiles_copy.df <- data.frame(profiles.df)
      
      for(m in 1:num_loci){
        profiles_copy.df <- filter(profiles_copy.df, profiles_copy.df[,m] == curr_output.df[1, m+1]) # what are we doing here?
      }
      
      copy_dim <- dim(profiles_copy.df)
      if(is.null(copy_dim) || copy_dim[1] == 0 || copy_dim[2] == 0){
        ST <- NA # Set sequence typing to NA because?
        mlst_type.df <- data.frame(ST, stringsAsFactors = FALSE)
        curr_output.df <- bind_cols(curr_output.df, mlst_type.df)
      } else {
        curr_output.df <- bind_cols(curr_output.df,
                                    select(profiles_copy.df, ST))
        mlst_type.df <- select(profiles_copy.df, ST) # Get the correct ST
      }
      
      #curr_output.df <- bind_rows(curr_output.df, mlst_type.df) # Add the mlst_type to the column (or row because tidy data)
      # CHANGED ####
    } else {
      mlst_type.df <- data.frame(curr_output.df[2])
    }
    
    sample_out.df <- bind_rows(sample_out.df, curr_output.df)
    sample_mut_out.df <- bind_rows(sample_mut_out.df, curr_mut_output.df) # NGSTAR ####
    
    if(locus == "list"){  # if this is done on ALL
      if(seq_type == "NGMAST"){ # NGMAST ####
        writeLines(paste(curr_sample_num, "\t por-", curr_output.df[1,2], "\t tbpB-", curr_output.df[1,3], "\t ST-",  mlst_type.df[1,1], sep = ""))
      } else {
        writeLines(paste(curr_sample_num, "\t", "ST-", mlst_type.df[1,1], sep = ""))
      }
    } else {
      writeLines(paste(curr_sample_num, "\t", locus, "-", mlst_type.df[1,1], sep = ""))
    }
  } # End Sample Loop
  
  writeLines("---- End Sample Loop ----")
  
  if(seq_type == "NGMAST"){
    
    # ---> addition in NGMAST here <----
    sample_out.df$porB <- paste("porB-", sample_out.df$porB, sep = "")
    sample_out.df$tbpB <- paste("tbpB-", sample_out.df$tbpB, sep = "")
    sample_out.df$ST <- paste("porB-", sample_out.df$porB, sep = "")
  }
  
  if(num_loci == 1){ # If we only used one loci - Not in first Galaxy iteration
    out_file <- paste(out_location, paste(Sys.Date(), org_id, "MLST-blast", "WADE.csv", sep = "_"), sep = "")
    write.csv(blast_out.df, out_file, row.names = FALSE)
    writeLines("DONE: general_mlst_pipeline()")
    
    return(blast_out.df)
    
  } else {
    
    # ---------------- File Output ----------------
    
    if(seq_type == "MLST"){
      # --- MLST out
      # These write the exact same things with different names...
      
      write.csv(sample_out.df,
                paste(out_location, paste(Sys.Date(), org_id, "MLST", "WADE.csv", sep = "_"), sep = ""),
                quote = FALSE,  row.names = FALSE)
      
      write.csv(sample_out.df,
                paste(out_location, paste(Sys.Date(), org_id, "MLST-profile", "WADE.csv", sep = "_"), sep = ""),
                quote = FALSE,
                row.names = FALSE)
      writeLines("--- profile_mlst.csv ---")
      
    } else {
      
      write.csv(sample_out.df,
                paste(out_location, paste(Sys.Date(), org_id, "LabwareUpload", "WADE.csv", sep = "_"), sep = ""),
                quote = FALSE,
                row.names = FALSE)
      
      write.csv(sample_out.df,
                paste(out_location, paste(Sys.Date(), org_id, "MLST-profile", "WADE.csv", sep = "_"), sep = ""),
                quote = FALSE,
                row.names = FALSE)
    }

    if(seq_type == "NGSTAR"){
      
      write.csv(sample_mut_out.df,
                paste(out_location, paste(Sys.Date(), org_id, "NGSTAR-mut", "WADE.csv", sep = "_"), sep = ""),
                quote = FALSE,
                row.names = FALSE)
    } else if(seq_type == "NGMAST"){
      # --- NGMAST out
      write.csv(sample_out.df,
                paste(out_location, paste(Sys.Date(), org_id, "NGMAST-profile", "WADE.csv", sep = "_"), sep = ""),
                row.names = FALSE)
    }

    writeLines("DONE: general_mlst_pipeline()")
    return(sample_out.df) # Files saved, why print?
  }

} # general_mlst_pipeline()
