#' Reads a .vcf file ignoring any headers, creates a data frame from the remaining lines
#'
#' 1. Grabs the header from the correct line
#' 2. Creates a data frame from all found mutations
#' 3. Change the column names of the data frame to the headers in 1.
#'
#' @param path The path to file corresponding .vcf file
#' @return A Data Frame of the data from a .vcf
#'

read_vcf <- function(path){
  # comment.char in read.table can only take one char. So we have to get the header first.

  if(file.exists(path)){

    curr_file <- file(path, open = "r") # open the file
    file_lines <- readLines(curr_file) # read the lines
    header_str <- file_lines[grep("#CHROM", file_lines)] # Get the header line

    header <- str_split(header_str, '\t') # convert the header to a list
    close(curr_file)

    # Read Table ####
    curr_file <- file(path, open = "r") # open the file
    in_file.df <- read.table(file = curr_file,
                             comment.char = "#",
                             stringsAsFactors = FALSE,
                             colClasses = "character")
    close(curr_file)

    names(in_file.df) <- as.list(header[[1]]) # name the columns
    in_file.df$SampleNo <- str_remove(basename(path), ".vcf") # Assign a new column tracking the sample number it's from
    typeof(in_file.df$SampleNo)
    # Return ####
    in_file.df

  } else {
    data.frame() # Return an empty data frame
  }
}
