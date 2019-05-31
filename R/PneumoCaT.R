#' PneumoCaT_Batch (Waffles version) :

#' Log in to Waffles
#' Download files from IRIDA using ngsArchiveLinker using folder structures
#' Rename files to remove _001.fastq using: prename 's/_001.fastq/.fastq/' */*.fastq
#' Run command line to send all samples in the "fastqs" folder :
#' $ for sample in `ls fastqs` ; do  sbatch -c 1 -p high --mem 4G --output $sample-%j.out --wrap="PneumoCaT.py -i fastqs/$sample -o Output/$sample -s /share/apps/samtools-0.1.19/samtools" ; done
#' $ watch sq
#' CTRL-C to exit watch.

#' @param sample_num Sample number associated with contig.fasta file

#' @return A table frame containing the results of the query
#' @export

PneumoCaT_pipeline <- function(sample_num){
  # TODO ####
  # - See if there's a way to abstract the innards of the if and else statements (they're similar)
  # - Remove the for loops, if mapping we have to use <<-

  # Directories ####
  # Hardcoded Dir(s) ####
  samplist <- here("data", "input", "assemblies", "list.csv")
  pneumocat_out <- here("data", "databases", "PNEUMO", "Output")

  #--------------------------------------------------------------------------------------Setup sample list table ####
  if(sample_num == "list") {
    samples.df <- read.csv(samplist, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  } else {
    samples.df <- data.frame(sample_num)
  }

  num_samples <- (dim(samples.df))[1]

  # Initialize Return ####
  output_profile.df <- data.frame()

  # Iteration ####
  for (m in 1L:num_samples) {

    # Variables ####
    curr_sample_num <- as.character(samples.df[m, "SampleNo"])

    snp_dir <- paste(pneumocat_out, "/", curr_sample_num, "/SNP_based_serotyping/", sep = "")
    map_dir <- paste(pneumocat_out, "/", curr_sample_num, sep = "")

    file_list <- list.files(snp_dir, "*.results.xml")
    print(file_list)

    typing_method <- ""   # "SNP" or "MAP" based
    Serotype1 <- ""
    QC_coverage1 <- ""
    QC_mean_depth <- ""
    QC_min_depth <- ""
    qc_mean <- ""

    QC_min_depth <- ""
    Serotype2 <- ""
    QC_coverage2 <- ""

    if(length(file_list) > 0) {
      typing_method <- "SNP-based"

      SNP_based_file <- paste(snp_dir, "/", file_list, sep = "")
      con <- file(SNP_based_file, open="r")
      linn <- readLines(con)
      close(con)

      for (i in 1:length(linn)) {

        if (str_detect(linn[i], "Serotype_Distinction\" value=")) {
          curr_ident <-  unlist(linn[i])
          ident_line.list <- unlist(strsplit(curr_ident, "\""))
          Serotype1 <- ident_line.list[4]
        }

        if (str_detect(linn[i], "Total_Hits") && (QC_coverage1 == "")) {
          curr_ident <-  unlist(linn[i])
          ident_line.list <- unlist(strsplit(curr_ident, "\""))
          QC_coverage1 <- ident_line.list[4]
        }

      } # End SNP Based For()
    } else {
      typing_method <- "Mapping-based"

      file_list2 <- list.files(map_dir, "*.results.xml")

      map_file <- paste(map_dir, "/", file_list2, sep = "")
      con <- file(map_file, open = "r")
      linn <- readLines(con)
      close(con)

      # Iteration ####
      for (i in 1:length(linn)) {

        if (str_detect(linn[i], "Serotype")) {
          curr_ident <-  unlist(linn[i])
          ident_line.list <- unlist(strsplit(curr_ident, "\""))
          Serotype1 <- ident_line.list[4]
        }

        if (str_detect(linn[i], "QC_coverage\" value")) {
          curr_ident <-  unlist(linn[i])
          ident_line.list <- unlist(strsplit(curr_ident, "\""))
          QC_coverage1 <- ident_line.list[4]
        }

        if (str_detect(linn[i], "QC_mean_depth"))
        {
          curr_ident <-  unlist(linn[i])
          ident_line.list <- unlist(strsplit(curr_ident, "\""))
          QC_mean_depth <- ident_line.list[4]
        }

        if (str_detect(linn[i], "QC_minimum_depth"))
        {
          curr_ident <-  unlist(linn[i])
          ident_line.list <- unlist(strsplit(curr_ident, "\""))
          QC_min_depth <- ident_line.list[4]
        }

        if (str_detect(linn[i], "QC_meanQ"))
        {
          curr_ident <-  unlist(linn[i])
          ident_line.list <- unlist(strsplit(curr_ident, "\""))
          qc_mean <- ident_line.list[4]
        }

        if (str_detect(linn[i], "\"second_value"))
        {
          curr_ident <-  unlist(linn[i])
          ident_line.list <- unlist(strsplit(curr_ident, "\""))
          Serotype2 <- ident_line.list[4]
        }

        if (str_detect(linn[i], "QC_coverage_second"))
        {
          curr_ident <-  unlist(linn[i])
          ident_line.list <- unlist(strsplit(curr_ident, "\""))
          QC_coverage2 <- ident_line.list[4]
        }
      } # End Mapping based For()
    }

    cov_summ_file <- paste(map_dir, "/coverage_summary.txt", sep = "")
    coverage_summary.df <- read.csv(cov_summ_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    seq_type_1 <- coverage_summary.df[1,1]
    seq_type_2 <- coverage_summary.df[2,1]
    seq_type_3 <- coverage_summary.df[3,1]

    coverage_1 <- coverage_summary.df[1,2]
    coverage_2 <- coverage_summary.df[2,2]
    coverage_3 <- coverage_summary.df[3,2]

    # Output ####
    writeLines(paste(curr_sample_num, "\t",
                     seq_type_1, " (", coverage_1, ")", "\t",
                     seq_type_2, "(", coverage_2, ")", "\t",
                     seq_type_3, "(", coverage_3, ")", sep =""))

    SampleResults.df <- data.frame(curr_sample_num,
                                   Serotype1,
                                   typing_method,
                                   QC_coverage1,
                                   QC_mean_depth,
                                   QC_min_depth,
                                   qc_mean,
                                   Serotype2,
                                   QC_coverage2,
                                   seq_type_1,
                                   coverage_1,
                                   seq_type_2,
                                   coverage_2,
                                   seq_type_3,
                                   coverage_3)

    output_profile.df <- bind_rows(output_profile.df, SampleResults.df)
  }

  serotype_upload.df <- data.frame(SampleNo = output_profile.df$curr_sample_num,
                                   Serotype = output_profile.df$Serotype1)

  # Output ####
  write_to_files(output_profile.df, serotype_upload.df)

  writeLines("DONE: PneumoCaT_pipeline() finished...")

  # Return ####
  output_profile.df
}

# Separate the writing of files to this function
write_to_files <- function(output_profile.df, serotype_upload.df){
  # write.csv(output_profile.df,
  #           here("data", "output", "pneumo_output_profile.csv"),
  #           row.names = FALSE,
  #           quote = FALSE)
  write.csv(output_profile.df,
            here("data", "output", "output_profile_PNEUMO_SEROTYPE.csv"),
            row.names = FALSE,
            quote = FALSE)
  write.csv(serotype_upload.df,
            here("data", "output", "LabWareUpload_PNEUMO_SERO.csv"),
            row.names = FALSE,
            quote = FALSE)

  # write.csv(output_profile.df, "C:\\WGS_Typing\\Output\\output_profile.csv", row.names = F, quote = FALSE) # Don't these two do the EXACT SAME THING?
  # write.csv(output_profile.df, "C:\\WGS_Typing\\Output\\output_profile_PNEUMO_SEROTYPE.csv", row.names = F, quote = FALSE)
  # write.csv(serotype_upload.df, "C:\\WGS_Typing\\Output\\LabWareUpload_PNEUMO_SERO.csv", row.names = F, quote = FALSE)
}

