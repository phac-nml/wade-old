#' Browse to raw downloaded .csv files to create the upload format table for Labware
#' Doesn't matter in what order files are selected
#' Will combine columns into Labware upload ready .csv
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @return A table frame containing the results of the query
#' @export
#'
#'

metrics <- function(org_id) {

  # TODO ####
  # - Do I need to allow it to take more than two files?
  #

  switch(org_id,
         GONO = { curr_dir <- here("data", "Gonorrhoea", "wgs_metrics") },
         GAS = { curr_dir <- here("data", "Streptococcus", "GAS", "wgs_metrics") },
         PNEUMO = { curr_dir <- here("data", "Streptococcus", "PNEUMO", "wgs_metrics")}
  )

  # File Managing ####
  all_files <- fileInput(default = "",
                            caption = "Select Files",
                            multi = TRUE) # list.files in the import directory in Docker!
  # Look into shinyFileChoose in the shinyFiles package!!!!! It (supposedly) works with Docker

  # Initialize DF's ####
  one.df <- read.csv(all_files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  two.df <- read.csv(all_files[2], header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  colnames(one.df)[1] <- "StrainID"
  colnames(two.df)[1] <- "StrainID"

  # Join the data frames together
  combined_files.df <- full_join(one.df, two.df, by = "StrainID") # Combine by StrainID
  combined_files.df  <- combined_files.df[order(combined_files.df$StrainID),] # Order by StrainID

  combined_files.df$StrainID <- sub("\\__.*", "", combined_files.df$StrainID)
  combined_files.df$N50.contig.length <- sub("\\,", "", combined_files.df$N50.contig.length)

  # Initialize Return Var ####
  metrics.df <-  data.frame(cbind(combined_files.df$StrainID,
                                  combined_files.df$X..of.Reads,
                                  combined_files.df$Estimated.Coverage,
                                  combined_files.df$Mean.contig.length,
                                  combined_files.df$N50.contig.length,
                                  combined_files.df$Number.of.contigs))
  names(metrics.df) <- list("Sample", "Num_Reads", "Coverage", "MeanContigLength", "N50ContigLength", "Comments")

  # Output ####
  write.csv(x = metrics.df,
            file = here("data", "output", "LabWareUpload_METRICS.csv"),
            quote = FALSE,
            row.names = FALSE)

  writeLines("DONE .. Output \ LabWareUpload_METRICS.csv has been created...")

  # Return ####
  metrics.df
}
