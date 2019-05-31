#' Molecular typing pipeline for WGS assemblies
#'
#' Takes Organism, Sample Number, Locus, to query a contig.fasta file
#' @param org_id Organism to query: GAS, PNEUMO or GONO
#' @param test_id AMR, TOXINS, VIRULENCE, NGSTAR, use MASTER for 16S id
#' @param sample_num Sample number associated with contig.fasta file
#' @param locus_id The locus to query, or enter "list" to use a list of alleles
#' @details GAS AMR loci:
#'
#' dfrF, dfrG, DHFR, ermA, ermB, ermT, folA, folP, gyrA, mefAE, msrD, parC, rpsJ, tetM, tetO
#' note that ermA == ermTR, ermC == ermT, mefA/E includes mefA and mefE; DHFR is the same as folA
#' ermB + = ERY-R/CLI-R; ermT = ERY-R/CLI-?; mefA + = ERY-R/CLI-S; ermA/TR = ERY-R/CLI-Ind.
#'
#' PNEUMO AMR loci:
#' cat ermB, ermTR, folA, folP, gyrA, parC, mefA, mefE, msrD, tetM, tetO, pbp1a, pbp2b, pbp2x
#' note that ermA == ermTR
#'
#' GAS VIRULENCE FACTORS:
#' covR => V128A
#' covS => I332V; E226G (1228 deletion)
#' rocA => V47A; R259Q; V333A; D396N; I404T; P405S
#' rpoB/rgg => V169I
#' nga => promoter mutations: A8G (K3R) or T13G (F5V) or C17T (A6V) (protein): F5V+A6V=pnga3
#' nga => D349G is the D330G mutation
#' cpA, fctA, fctB, lepA and srtC1 all pilin proteins - use only fctA (major pilin protein)
#' fbp54, grab, ideS_mac, mf_spd, speB, scpA, ska always positive?
#' hasA,B,C are the same operon, use only hasA
#' @return A table frame containing the results of the query
#' @export

master_blastr <- function(org_id, test_id, sample_num, locus_id){
  # Variables ####
  Variable <- NA
  Blast_evalue <- "10e-50"            #sets sensitivity of Blast gene match 10e-50 to 10e-150; use 10e-5 for primers

  #--------------------------------------------------------------------------------------------------------

  SampList <- here("data", "assemblies", "list.csv")

  writeLines(paste("SampleNo \t", test_id, sep = ""))

  if (org_id == "GONO" && (locus_id %in% c("penA", "mtrR", "porB", "ponA", "gyrA", "parC", "rRNA23S"))) {
    test_id <- "NGSTAR"
  }

  # Hardcoded Dir(s) ####

  contigs_dir <- here("data", "databases", org_id, "assemblies") # Change depending on Organism!!
  lookup_dir <- here("data", "databases", org_id, test_id)

  # switch(org_id,
  #        GAS={
  #          contigs_dir <- "W:\\Projects\\Project_GC_WalterD\\MiSeq\\Streptococcus\\GAS\\contigs\\"
  #
  #          switch(test_id,
  #                 AMR={lookup_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\GAS\\Wamr_R\\"},
  #                 TOXINS={lookup_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\GAS\\Toxins_R\\"},
  #                 VIRULENCE={lookup_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\GAS\\Virulence_R\\"},
  #                 MASTER={lookup_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\GAS\\Master_BlastR\\"}
  #          )
  #
  #        },
  #
  #        PNEUMO={
  #          contigs_dir <- "W:\\Projects\\Project_GC_WalterD\\MiSeq\\Streptococcus\\Pneumo\\contigs\\"
  #
  #          switch(test_id,
  #                 AMR={lookup_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\PNEUMO\\Wamr_R\\"},
  #                 TOXINS={lookup_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\PNEUMO\\Toxins_R\\"},
  #                 VIRULENCE={lookup_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\PNEUMO\\Virulence_R\\"},
  #                 MASTER={lookup_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\PNEUMO\\Master_Blaster_R\\"}
  #          )
  #        },
  #
  #        GONO={
  #          contigs_dir <- "W:\\Projects\\Project_GC_WalterD\\MiSeq\\Gonorrhoea\\contigs\\"
  #
  #          switch(test_id,
  #                 AMR={lookup_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\GONO\\Wamr_R\\"},  #but use NGSTAR databases for penA, etc
  #                 AMR_LW={lookup_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\GONO\\Wamr_R\\"},  #but use NGSTAR databases for penA, etc
  #
  #                 NGMAST={lookup_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\GONO\\NGMAST_R\\"},
  #                 MASTER={lookup_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\GONO\\Master_Blaster_R\\"},
  #                 NGSTAR={lookup_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\GONO\\NGSTAR_R\\"}
  #
  #          )
  #        }
  #
  # )

  #------------------------------------------------------------------------------------------------------####

  # Hardcoded Dir(s) ####
  unlink("C:\\WGS_Typing\\Output\\output_dna.fasta") #this deletes the file!
  unlink("C:\\WGS_Typing\\Output\\output_dna_notfound.fasta")
  unlink("C:\\WGS_Typing\\Output\\output_aa.fasta")
  # unlink("C:\\Temp\\Typing\\Output\\output_aa_notfound.fasta") #this deletes the file!

  samples.df <- get_samples(org_id, sample_num)
  print(samples.df)
  # if(sample_num == "list") {
  #   samples.df <- read.csv(SampList, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  # } else
  # {
  #   samples.df <- data.frame(sample_num, Variable)
  # }

  num_samples <- (dim(samples.df))[1]

  if(locus_id == "list") {
    #writeLines(paste(lookup_dir, "/temp/loci.csv", sep = ""))
    loci.df <- read.csv(file = paste(lookup_dir, "/temp/loci.csv", sep = ""),
                             header = TRUE,
                             sep = ",",
                             stringsAsFactors = FALSE)
  } else {
    loci.df <- data.frame(locus_id)
  }

  # Variables ####
  num_loci <- (dim(loci.df))[1]

  # Progress ####
  progress_ratio <- num_samples * num_loci

  if (test_id == "MASTER") {
    for(q in 1L:num_loci) {

      locus <- as.character(loci.df[q,1])
      locus_dna_lookup <- paste(lookup_dir, "/allele_lkup_dna/", locus, ".fasta", sep = "")

      if(file.exists(locus_dna_lookup)) {
        blast_format_command <- paste("makeblastdb -in ", locus_dna_lookup, " -dbtype nucl", sep = "")
        try(system(blast_format_command))
      }
    }
  }

  # Iteration ####
  for (m in 1L:num_samples){ # Loop through the samples

    # Variables ####
    curr_sample_num <- as.character(samples.df[m, "SampleNo"])
    curr_sample_var <-as.character(samples.df[m, "Variable"])
    curr_sample.df <- filter(samples.df, SampleNo == curr_sample_num)

    sample_profile <- ""

    for(p in 1L:num_loci) { # Loop through the loci

      locus <- as.character(loci.df[p,1])
      locus_dna_lookup <- paste(lookup_dir, "/allele_lkup_dna/", locus, ".fasta", sep = "")
      locus_file <- paste(lookup_dir, "/wildgenes/", locus, ".fasta", sep = "") # Hardcoded Dir ####

      if(!file.exists(locus_file)) { # Check if the Locus file exists!
        sample_error.df <- data.frame(Loucs_ID = locus,
                                      Output = "File Not Found",
                                      stringsAsFactors = FALSE)
        return(sample_error.df)
      }

      if(file.exists(locus_dna_lookup)) {
        dna_present = TRUE
      } else {
        dna_present = FALSE
      }

      AlleleLine <- ""
      AlleleInfo <- NA
      AlleleInfo[1] <- NA
      AlleleInfo[2] <- NA
      AlleleInfo[3] <- NA
      AlleleInfo[4] <- NA
      IdentLine <- ""


      #  .....................................................................BLAST and parse blastout.txt

      query_file <- paste(contigs_dir, "/", curr_sample_num, ".fasta", sep = "")
      db_dir <- here("data", "output", "temp", "queryfile.fasta") # "C:\\WGS_Typing\\temp\\queryfile.fasta" # Hardcoded Dir ####

      if (!file.copy(query_file, db_dir, overwrite = T)) { # Can we copy the file? Should we copy the file?

        AlleleInfo[1] <- "Sample_Err"
        AlleleInfo[2] <- NA
        AlleleInfo[3] <- NA
        AlleleInfo[4] <- NA
        IdLine <- ""
        IDpercent2 <- ""

      } else { #make blast database of contig file, then blast locus against it
        format_command <- paste("makeblastdb -in ", db_dir, " -dbtype nucl", sep = "")
        try(system(command = format_command,
                   intern = TRUE))

        blast_db_file <- paste(db_dir, "/", curr_sample_num, ".fasta", sep = "")
        output_location <- paste(lookup_dir, "/temp/blastout.txt", sep = "")

        blast_command <- paste("blastn -query ", locus_file, " -db ", db_dir, " -out ", output_location, " -evalue ", Blast_evalue, sep = "")
        try(system(blast_command, intern = TRUE))

        con <- file(output_location, open="r")
        linn <- readLines(con)
        close(con)

        #check if gene was found in BLAST
        blast_result <- NA

        for (i in 1:length(linn)) {
          if (str_detect(linn[i], "No hits found")) {
            blast_result <- "NEG"
            AlleleInfo[1] <- "NEG"
            AlleleInfo[2] <- ""
            AlleleInfo[3] <- ""
            AlleleInfo[4] <- ""
            break()
          } else {
            blast_result <- "POS"
            AlleleInfo[1] <- "POS"
            AlleleInfo[2] <- ""
            AlleleInfo[3] <- ""
            AlleleInfo[4] <- ""
          }
        }

        # Variables ####
        DNASeqLine_str <- ""
        WTDNASeqLine_str <- ""
        TimesThrough <- 0
        WTlengthLine <- ""
        WTlength <- 0
        IdLine <- ""
        IDlength <- ""
        IDpercent <- ""
        IdNum <- ""

        if (blast_result == "POS") {
          # Iteration ####
          for (i in 1:length(linn)) {
            if (str_detect(linn[i], "Score =")) { #set counter to parse only first align block of blastout
              TimesThrough <- TimesThrough + 1
            }

            if ((str_detect(linn[i], "Length=" )) & (WTlength == 0)) { #get length of WT gene
              WTlengthLine <-  unlist(linn[i])
              WTlength <- as.numeric(substr(WTlengthLine, 8, 50))
            }

            if (TimesThrough == 1) {
              if (str_detect(linn[i], "Identities")) {
                IdLine <-  unlist(linn[i])
                #IdLine <- " Identities = 2431/2460 (99%), Gaps = 0/2460 (0%)"
                IdLine <- substr(IdLine, 15, 50)
                IdLineParts <- strsplit(IdLine, "/")
                IDLineParts2 <- unlist(IdLineParts)
                IDlength <- as.numeric(IDLineParts2[1])
                IDcoverage <- ((IDlength / WTlength) * 100)
                IDpercent <- sprintf("%3.1f%%", IDcoverage)
              }

              if (str_detect(linn[i], "Query ")) {
                QueryLine <-  unlist(strsplit(linn[i], " "))
                QueryLine <- QueryLine[QueryLine != ""]
                WTDNASeqLine_str <- paste(WTDNASeqLine_str, QueryLine[3], sep = "")
              }

              if (str_detect(linn[i], "Sbjct ")) {
                SbjctLine <-  unlist(strsplit(linn[i], " "))
                SbjctLine <- SbjctLine[SbjctLine != ""]
                DNASeqLine_str <- paste(DNASeqLine_str, SbjctLine[3], sep = "")
              }
            }
          }

          # Variables ####
          WTDNASeqLine <- DNAString(WTDNASeqLine_str)
          WTDNASeqLine_NoDash_str <- str_replace_all(WTDNASeqLine_str, "-", "")
          WTDNASeqLine_NoDash <- DNAString(WTDNASeqLine_NoDash_str)

          if (locus_id != "list") {
            cat("\n\n>", locus , "(Wildtype)\n", WTDNASeqLine_NoDash_str, "\n", sep ="")
          }

          # Variables ####
          DNASeqLine <- DNAString(DNASeqLine_str)
          DNASeqLine_NoDash_str1 <- str_replace_all(DNASeqLine_str, "-", "")
          DNASeqLine_NoDash_str <- str_replace_all(DNASeqLine_NoDash_str1, "N", "")
          DNASeqLine_NoDash <- DNAString(DNASeqLine_NoDash_str)

          if (locus_id != "list") {
            cat("\n\n>", curr_sample_num, "_", locus , "\n", DNASeqLine_NoDash_str, "\n", sep ="")
          }

          #-------------------------------------------------------------------------------make Protein sequence
          WTAASeqLine <- translate(WTDNASeqLine_NoDash)
          WTAASeqLine_str <- toString(WTAASeqLine)

          if (locus_id != "list") {
            cat("\n\n>", locus , "(Wildtype)\n", WTAASeqLine_str, "\n", sep ="")
          }

          AASeqLine <- translate(DNASeqLine_NoDash)
          AASeqLine_str <- toString(AASeqLine)
          # AASeqLine_str <- str_replace_all(AASeqLine_str, "[*]", "")
          if (locus_id != "list") {
            cat("\n\n>", curr_sample_num, "_", locus , "\n", AASeqLine_str, "\n", sep ="")
          }

          # Variables ####
          #-----------------------------------make sure AASeqLine, WTAASeqLine are listed in correct order!!!!!!!!!!!!
          globalAlign_AA <- pairwiseAlignment(AASeqLine, WTAASeqLine, substitutionMatrix = "BLOSUM50", gapOpening = -2,
                                              gapExtension = -8, scoreOnly = FALSE)

          WTAASeqLine_aln <- subject(globalAlign_AA)
          WTAASeqLine_aln_str <- toString(WTAASeqLine_aln)
          #cat("\n\n>", locus, "_alignment\n", WTAASeqLine_aln_str, sep = "")

          AASeqLine_aln <- pattern(globalAlign_AA)
          AASeqLine_aln_str <- toString(AASeqLine_aln)
          #cat("\n\n>", sample_num, "_", locus , "_alignment\n", AASeqLine_aln_str, sep = "")

          DNASeqLine_aln <- ""

          for (j in 1:str_length(WTDNASeqLine_str)) {
            if (str_sub(WTDNASeqLine_str, j, j) == str_sub(DNASeqLine_str, j, j)) {
              DNASeqLine_aln <- paste(DNASeqLine_aln, ".", sep = "")
            } else {
              DNASeqLine_aln <- paste(DNASeqLine_aln, str_sub(DNASeqLine_str, j, j), sep = "")
            }
          }
          #if (sample_num != "list")
          #{
          #  cat("\n\nDNA Alignment:\n", DNASeqLine_aln, sep = "")
          #}
          if (locus_id != "list") {
            cat("\n\nDNA Alignment:\n", DNASeqLine_aln, "\n", sep = "")
          }

          #------------------------------------------------------------------------------mutations
          motifs <- ""

          # MANY REPEATED IF's ##################################
          if (locus == "pbp1a") {

            motif_1<-substr(AASeqLine_aln_str, 370, 373)
            if (motif_1 == "STMK") {
              motif_1<-"WT"
            }

            motif_2<-substr(AASeqLine_aln_str, 466, 468)
            if (motif_2 == "SSN") {
              motif_2<-"WT"
            }

            motif_3<-substr(AASeqLine_aln_str, 557, 559)
            if (motif_3 == "KTG") {
              motif_3<-"WT"
            }
            motifs <- paste(motif_1, motif_2,motif_3, sep = "/")
          }
          else if (locus == "pbp2b") {

            motif_1<-substr(AASeqLine_aln_str, 386, 389)
            if (motif_1 == "SVVK") {
              motif_1<-"WT"
            }

            motif_2<-substr(AASeqLine_aln_str, 443, 445)
            if (motif_2 == "SSN") {
              motif_2<-"WT"
            }

            motif_3<-substr(AASeqLine_aln_str, 615, 617)
            if (motif_3 == "KTG") {
              motif_3<-"WT"
            }
            motifs <- paste(motif_1, motif_2,motif_3, sep = "/")
          }
          else if (locus == "pbp2x") {

            motif_1<-substr(AASeqLine_aln_str, 337, 340)
            if (motif_1 == "STMK") {
              motif_1<-"WT"
            }

            motif_2<-substr(AASeqLine_aln_str, 395, 397)
            if (motif_2 == "SSN") {
              motif_2<-"WT"
            }

            motif_3<-substr(AASeqLine_aln_str, 546, 549)
            if (motif_3 == "LKSG") {
              motif_3<-"WT"
            }
            motifs <- paste(motif_1, motif_2,motif_3, sep = "/")
          }

          AASeqLine_aln_disp <- ""
          mutations <- ""

          for (k in 1:str_length(WTAASeqLine_aln_str)) {

            if (str_sub(WTAASeqLine_aln_str, k, k) == str_sub(AASeqLine_aln_str, k, k)) {
              AASeqLine_aln_disp <- paste(AASeqLine_aln_disp, ".", sep = "")
            } else {
              AASeqLine_aln_disp <- paste(AASeqLine_aln_disp, str_sub(AASeqLine_str, k, k), sep = "")
              mutations <- paste(mutations, str_sub(WTAASeqLine_aln_str, k, k), k, str_sub(AASeqLine_str, k, k), " ", sep = "")
            }
          }

          # Output ####
          if (locus_id != "list") {
            cat("\n\nProtein Alignment:\n", AASeqLine_aln_disp, "\n", sep = "")
            cat("\nMutations: ", mutations, "\n\n", sep = "")
            cat("\nMotifs: ", motifs, "\n\n", sep = "")
          }
          #----------------------------------------------------------------------------------Lookup Alleles DNA
          # write DNASeqLine_NoDash_str to a file named = querygene.fasta,
          # BLAST vs. lookup table,
          # parse out the allele numbers


          Seq_File <- paste(lookup_dir, "/temp/querygene.fasta", sep = "") # "C:\\WGS_Typing\\temp\\querygene.fasta" # Hardcoded Dir ####
          sink(Seq_File, split=FALSE, append = FALSE)
          cat(">", curr_sample_num, "_", locus , "\n", DNASeqLine_NoDash_str, sep ="")
          sink()

          ExactMatchFound <- FALSE

          if(dna_present) {

            #BLAST lookup table
            # Hardcoded Dir ####
            output_location <- paste(lookup_dir, "/temp/blastout2.txt", sep = "")

            BlastCommand2 <- paste("blastn -query ", Seq_File, " -db ", locus_dna_lookup, " -out ", output_location, " -evalue 10e-85 -num_alignments 1", sep = "")
            try(system(BlastCommand2, intern = TRUE))

            con <- file(output_location, open="r")
            linn <- readLines(con)
            close(con)

            #parse blastout2
            for (i in 1:length(linn)) {

              if (str_detect(linn[i], "Identities")) {
                IdentLine <-  unlist(linn[i])
                IdentLine <- substr(IdentLine, 15, 50)

                if (str_detect(IdentLine, "100%")) {
                  ExactMatchFound <- TRUE
                }
              }

              if (str_detect(linn[i], ">")) {
                AlleleLine <- unlist(linn[i])
                AlleleLine <- substr(AlleleLine, 2, 50)
                AlleleParts <- strsplit(AlleleLine, "_")
                AlleleParts2 <- unlist(AlleleParts)
              }
            }

            if(ExactMatchFound) {   #if only not found in fasta make files here.
              #AlleleParts2[1] holds POS/NEG; [2] allele number; [3] mutations or WT; [4] extra info
              AlleleInfo[2] <- AlleleParts2[2]  #allele number
              AlleleInfo[3] <- AlleleParts2[3]  #mutations
              AlleleInfo[4] <- AlleleParts2[4]  #comments
            } else {
              AlleleInfo[2] <- "NF"
              AlleleInfo[3] <- "???"
              AlleleInfo[4] <- ""
            }
          }#close bracket for DNA lookup file exists check


          #-------------------------------------------write a fasta file of all sequences output_dna.fasta and output_aa.fasta
          sink(file = here("data", "output_dna.fasta"),
               split = FALSE,
               append = TRUE)
          cat(">", locus, "_", curr_sample_num, "_", locus, AlleleInfo[2], "_", AlleleInfo[3], "_", curr_sample_var, "\n", DNASeqLine_NoDash_str, "\n", sep ="")
          sink()

          if (AlleleInfo[2] == "NF") {
            sink(file = here("data", "output_dna_notfound.fasta"),
                 split = FALSE,
                 append = TRUE)
            cat(">", locus, "_", curr_sample_num, "_", locus, AlleleInfo[3], "_", curr_sample_var, "\n", DNASeqLine_NoDash_str, "\n", sep ="")
            sink()
          }

          sink(file = here("data", "output_aa.fasta"),
               split = FALSE,
               append = TRUE)
          cat(">", locus, "_", curr_sample_num, "_", curr_sample_var, "\n", AASeqLine_str, "\n", sep ="")
          sink()

          IDpercent2 <- paste(IDlength, "/", WTlength, " (", IDpercent, ")", sep = "") #made a blast ID line based on full WT gene

          if ((locus == "porB") & (IDcoverage < 90)) {
            AlleleInfo[4] <- "porB1a"
          } else if ((locus == "porB") & (IDcoverage >= 90)) {
            AlleleInfo[4] <- "porB1b"
          }

        } else { #close bracket for BLAST positive
          # cat(curr_sample_num, "\t", locus, "\tBLAST Negative\n", sep = "")
          IdLine <-""
          IDpercent2 <- ""
        }
      } #close bracket for contig file found.

      col1_name <- paste(locus, "_result", sep = "")
      col2_name <- paste(locus, "_allele", sep = "")
      col3_name <- paste(locus, "_mutations", sep = "")
      col4_name <- paste(locus, "_comments", sep = "")
      col5_name <- paste(locus, "_BlastID", sep = "")
      col6_name <- paste(locus, "_PctWT", sep = "")

      headers <- list(col1_name, col2_name, col3_name, col4_name, col5_name, col6_name)

      if(p == 1) {
        OutputLocus.df <- data.frame(AlleleInfo[1], AlleleInfo[2], AlleleInfo[3], AlleleInfo[4], IdLine, IDpercent2, stringsAsFactors = FALSE)
        names(OutputLocus.df) <- headers
      } else {
        OutputLocus2.df <- data.frame(AlleleInfo[1], AlleleInfo[2], AlleleInfo[3], AlleleInfo[4], IdLine, IDpercent2, stringsAsFactors = FALSE)
        names(OutputLocus2.df) <- headers
        OutputLocus.df <- cbind(OutputLocus.df, OutputLocus2.df)
      }



      #} #close bracket for contig file found.


      ProfileEntry <- ""

      if (AlleleInfo[1] == "POS") {
        ProfileEntry<-locus

        if (dna_present) {

          if (AlleleInfo[3] == "" | is.na(AlleleInfo[3]) | AlleleInfo[3] == "???") {
            ProfileEntry <- paste(locus, AlleleInfo[3], sep = " ")
          } else if (AlleleInfo[3]=="WT" | AlleleInfo[3]=="WT/WT" | AlleleInfo[3]=="WT/WT/WT") {
            ProfileEntry <- ""
          } else {
            ProfileEntry <- paste(locus, AlleleInfo[3], sep = " ")
          }
        }

        if (sample_profile == "") {
          sample_profile <- ProfileEntry
        } else {
          if (ProfileEntry!="") {
            ProfileEntry <- paste("-", ProfileEntry, sep = "")
          }
          sample_profile <- paste(sample_profile, ProfileEntry, sep = "")
        }
      }

      cat(curr_sample_num, locus, AlleleInfo[1], AlleleInfo[2], AlleleInfo[3], AlleleInfo[4], IdLine, IDpercent2, "\n", sep = "\t")
      incProgress(amount = 1/progress_ratio,
                message = paste(curr_sample_num, "on", locus)) # Progress ####

    } #end of locus list loop

    sample_profile.df <- data.frame(SampleProfile = sample_profile)
    cat("Molecular Profile:  ", sample_profile, "\n\n", sep = "")
    OutputLocus1.df <- bind_cols(curr_sample.df, OutputLocus.df, sample_profile.df)

    if(m == 1) {  #if first sample make one row profile table, otherwise add new row to table
      OutputProfile.df <- data.frame(OutputLocus1.df)
      headers2 <- names(OutputLocus1.df)
      names(OutputProfile.df) <- headers2
    } else {
      OutputProfile.df <- rbind(OutputProfile.df, OutputLocus1.df)
    }

    #incProgress(amount = 1/num_samples) # PROGRESS ####

  } #close brack for sample list loop

  outfile <- paste(here("data", "input", "labware"), "/output_profile_", org_id, "_", test_id, ".csv", sep = "")

  writeLines(paste("Writing to: ", outfile))
  write.csv(OutputProfile.df, outfile, row.names = F)

  #----------------------------------------------------------------Make LABWARE UPLOAD FILES
  if ((org_id=="GAS") & (test_id=="AMR") & (num_loci > 1)) {
    OutputProfile.df <- labware_gas_amr()
  }

  if (org_id == "GAS" & (test_id=="TOXINS" & num_loci > 1)) {
    OutputProfile.df <- labware_gas_toxins()
  }
  #----------------------------------------------------------------Make LABWARE UPLOAD FILES
  # Output ####
  writeLines("DONE: MasterBlastR finished....")

  # Return ####
  return(OutputProfile.df)

}
