library(shiny)
library(jsonlite)
library(shinydashboard)
library(shinyWidgets)
library(tidyverse)
library(here)
library(GalaxyConnector)

if(!require(wade)){
  options(repos = c(CRAN = "http://cran.rstudio.com"))
  if (!require(remotes)) { install.packages("remotes") }
  remotes::install_github("phac-nml/wade")
  library(pavian)
}

# if (!require(Rsamtools)) {
#   source("https://bioconductor.org/biocLite.R")
#   biocLite("Rsamtools")
#   library(Rsamtools)
# }

if(!require(GalaxyConnector)){
  options(repos = c(CRAN = "http://cran.rstudio.com")) 
  if(!require(remotes)){ install.packages("remotes") }
  remotes::install_github("justinband/r-galaxy-connector", ref="pavian") # This needs to change to phac-nml/galaxy-connector
  library(GalaxyConnector)
}

if (!dir.exists(rappdirs::user_config_dir("wade", expand = FALSE))) {
  dir.create(rappdirs::user_config_dir("wade", expand = FALSE),
             recursive = TRUE)
}

# Shiny app call
shiny::shinyApp(wade::ui, wade::server, enableBookmarking="server")