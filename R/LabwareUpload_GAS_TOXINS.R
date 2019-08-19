#' run MasterBlastR first, then this script outputs csv to LabWareUpload_GAS_TOXINS.csv in Output folder
#'
#' Then run this analysis to combine data the full amr profile to upload to LabWare.
#'
#'
#' @return A table frame containing the results of the query
#' @export
#'
#'
#'

labware_gas_toxins <- function() {

  # TODO ####
  # - Make sure labware.df is actually being created the correct way
  # - Make sure the in.csv found in Output_to_delete/ is the correct format
  #

  output.df <- read.csv(paste(here("data", "input", "labware", "output_profile_GAS_TOXINS.csv")), # What data should we be reading here?
                        header = TRUE,
                        sep = ",",
                        stringsAsFactors = FALSE)

  output.df$SampleProfile[output.df$SampleProfile == ""] <- "Error" # set it to an error?? Or is it experimental error

  # Initialize DF's ####
  # Initialize Return Var ####
  labware.df <- data.frame(
    SampleNo = output.df$SampleNo.filename,
    select(output.df, contains("result")),
    SampleProfile = output.df$SampleProfile
  )

  # Output ####
  write.csv(x = labware.df,
            file = here("data", "output", "LabWareUpload_GAS_TOXINS.csv"),
            quote = FALSE,
            row.names = FALSE)

  writeLines("DONE: labware_gas_toxins()")

  # Return ####
  labware.df
}

