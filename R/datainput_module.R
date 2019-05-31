selectDataPanel <- function(ns) {
  api <- Sys.getenv("GX_API_KEY")
  url <- Sys.getenv("GX_GALAXY_URL")
  history_id <- Sys.getenv("GX_HISTORY_ID")
  GalaxyConnector::gx_init(API_KEY = api, GALAXY_URL = url, HISTORY_ID = history_id) # Initialize our pkg env
  
  if(api == "" || url == "" || history_id == ""){
    # We should exception handle here
    # Print a message to the user, but for now let's just leave it blank
    user_data <- "Unable to connect to Galaxy"
  } else {
    user_data <- GalaxyConnector::gx_list_history_datasets()['name']
  }

  tabPanel("Select data from history",
           selectizeInput(inputId = ns("select_dataset"),
                          label = "Select by dataset name",
                          choices = user_data,
                          multiple = FALSE,
                          selected = NULL
           ),
           shiny::actionButton(inputId = ns("btn_confirm_selection"),
                               label = "Confirm")
  )
}

dataInputModuleUI <- function(id){
  ns <- shiny::NS(id)
  selectDataPanel(ns)
}
