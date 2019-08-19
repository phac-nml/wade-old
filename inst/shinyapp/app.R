library(devtools)
library(dplyr)
library(GalaxyConnector)
library(jsonlite)
library(here)
library(purrr)
library(RCurl)
library(shiny)
library(shinydashboard)
library(shinyFiles)
library(shinyWidgets)
library(stringr)

if(!require(wade)){
  options(repos = c(CRAN = "http://cran.rstudio.com"))
  if (!require(remotes)) { install.packages("remotes") }
  remotes::install_github("phac-nml/wade")
  library(pavian)
}

if(!require(GalaxyConnector)){
  options(repos = c(CRAN = "http://cran.rstudio.com")) 
  if(!require(remotes)){ install.packages("remotes") }
  remotes::install_github("phac-nml/r-galaxy-connector", ref="master") # This needs to change to phac-nml/galaxy-connector
  library(GalaxyConnector)
}

if (!dir.exists(rappdirs::user_config_dir("wade", expand = FALSE))) {
  dir.create(rappdirs::user_config_dir("wade", expand = FALSE),
             recursive = TRUE)
}

options(shiny.maxRequestSize = 500*1024^2)

# Shiny app call
shiny::shinyApp(wade::ui, wade::server, enableBookmarking="server")