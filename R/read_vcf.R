#' Reads a .vcf file ignoring any headers, creates a data frame from the remaining lines
#'
#' 1. Grabs the header from the correct line
#' 2. Creates a data frame from all found mutations
#' 3. Change the column names of the data frame to the headers in 1.
#'
#' @param path The path to file corresponding .vcf file
#' @return A Data Frame of the data from a .vcf
#' @importFrom stringr str_remove str_split
#' @import utils
#'

read_vcf <- function(path){
  # comment.char in read.table can only take one char. So we have to get the header first.

  if(file.exists(path)){

    # Read Header####
    file_lines <- readLines(path) # read the lines
    header_str <- file_lines[grep("#CHROM", file_lines)] # Get the header line
    header <- str_split(header_str, '\t') # convert the header to a list

    # Read Table ####
    in_file.df <- read.table(file = path,
                             comment.char = "#",
                             stringsAsFactors = FALSE,
                             colClasses = "character")
    names(in_file.df) <- header[[1]] # name the columns
    in_file.df$SampleNo <- str_remove(basename(path), ".vcf") # Assign a new column tracking the sample number it's from
    
    # Return ####
    in_file.df

  } 
}
