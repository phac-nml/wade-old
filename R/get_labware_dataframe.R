#' get_labware_dataframe.R
#'
#' Creates and returns an empty data frame with headers following Labware guidelines
#'
#' @return A data frame
#' @export
#'

get_labware_dataframe <- function(){
  genes <- c("ermA",
             "ermB",
             "ermT",
             "mefAE",
             "gyrA",
             "parC",
             "tetM",
             "tetO",
             "cat",
             "dfrF",
             "dfrG",
             "msrD",
             "folA",
             "folP",
             "tetT",
             "catQ")

  genes_info <- map(genes, function(i){ # Attach the info to the genes
    info <- c("result",
              "allele",
              "mutations",
              "comments",
              "BlastID",
              "PctWT")

    curr_gene_info <- map(info, ~ paste(i, "_", .x, sep = "")) # Combine the gene and the info into a list
  })

  all_info <- unlist(genes_info) # get one big list
  labware.df <- data.frame(matrix(ncol = 99)) # create an empty df
  names(labware.df) <- c("SampleNo",
                         "Variable",
                         all_info,
                         "SampleProfile") # Rename the df
}
