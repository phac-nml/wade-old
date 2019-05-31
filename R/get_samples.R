#' determine_sample
#'
#' This function determines whether the user selected "list" or a specific sample number
#' and creates SampleList.df accordingly.
#' Used in: 23SrRNA.R, amr_databases.R, emm.R, master.R, pneumo_cat.R, vfdb.R
#'
#' @param sample_num the sample the user input
#' @return A data frame containing a read .csv file or containing a single sample

# TODO
#   - Change data location based on the function accessing this function
#   - Sometimes Variable needs to be omitted (PneumoCaT.R)

get_samples <- function(org_id, sample_num) {
  Variable = NA

  if(sample_num == "list") {# If it's multiple
    path <- here("data", "databases", org_id, "assemblies", "list.csv")

    if(file.exists(path)){
      sample_list.df <- read.csv(file = path,
                                 header = TRUE,
                                 sep = ",",
                                 stringsAsFactors = FALSE) # File not found -- Raising errors!
    } else {
      stop(paste("get_samples.R: ", path, " does not exist", sep = ""))
    }
  } else {
    sample_list.df <- data.frame(SampleNo = sample_num, Variable)
  }

  sample_list.df
}
