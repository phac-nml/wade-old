#' Manta Server Function
#'
#' @param input Input object
#' @param output Output object
#' @param session Session object
#'
#' @export
server <- function(input, output, session){
  file_out <- NULL
  final_out.df <- data.frame()
  
  names <- c("parent_dir", "subdir_id", "filename")
  fullpath <- ""
  downloaded_samples.df <- data.frame(matrix(ncol = 3, nrow = 0), stringsAsFactors = FALSE) # Initialize to an empty df
  names(downloaded_samples.df) <- names # Give the df some col names
  
  # ---- Grab the Galaxy info from the Environment ---- #
  api <- Sys.getenv("GX_API_KEY")
  url <- Sys.getenv("GX_GALAXY_URL")
  history_id <- Sys.getenv("GX_HISTORY_ID")
  
  GalaxyConnector::gx_init(API_KEY = api, GALAXY_URL = url, HISTORY_ID = history_id) # Initialize our pkg env
  
  if(api == "" || url == "" || history_id == ""){
    # We should exception handle here
    # Print a message to the user, but for now let's just leave it blank
    user_data <- "Unable to connect to Galaxy"
    all_data <- "Unable to connect to Galaxy"
  } else {
    user_data <- GalaxyConnector::gx_list_history_datasets()['name']
    all_data <- dplyr::filter(GalaxyConnector::gx_list_history_datasets(), deleted == FALSE) # Filter out any deleted dataset
  }
  
  output$selectize <- renderUI({
    selectizeInput(inputId = "select_dataset",
                   label = "Select by dataset name",
                   choices = user_data,
                   multiple = FALSE,
                   selected = NULL
    )
  })
  
  output$confirm_selectize <- renderUI({
    shiny::actionButton(inputId = "btn_confirm_selection",
                        label = "Confirm")
  })
  
  observeEvent(input$btn_confirm_selection, {
    if(input$select_dataset != ''){ # Make sure that we have something selected first!
      
      data_hid.df <- dplyr::filter(all_data, name == input$select_dataset) # Grab the data that's been selected
      
      if(dim(data_hid.df)[1] != 0){
        
        # Take the first piece of data
        #   This could possibly be an issue since with collections there are multiple versions of the same data (by name).
        #   If one of the pieces of data we need are changed then things can possibly go wrong
        
        # Using [1,'hid'] because this will happen for each single data piece.
        # With a collection it only needs the parent directory so grabbing the first hid will also work
        datapath <- GalaxyConnector::gx_get(data_hid.df[1, 'hid']) 
        
        if(!is.null(datapath)){
          fullpath <- base::dirname(datapath) # Get everything up until the name of the data
          parent_dir <- base::dirname(fullpath) # Append dir_num so we can use list.files!
          dir_num <- base::basename(fullpath) # Get only the directory the data exists to use for tidy data
          
          combine.df <- data.frame("parent_dir" = parent_dir, "subdir_id" = dir_num, "filename" = list.files(fullpath), stringsAsFactors = FALSE) %>% 
            dplyr::filter(!(filename %in% unlist(downloaded_samples.df["filename"]))) # Filter out selected data (works with collections here)
          
          if(dim(combine.df)[1] != 0){
            downloaded_samples.df <<- rbind(downloaded_samples.df, combine.df) # We need to change it globally. Could also use a reactive?
            
            # Use our downloaded data frame
            samples <- downloaded_samples.df %>% pull(filename)
            
            sendSweetAlert(session = session,
                           title = "Confirmed",
                           text = "Data was added",
                           type = "success",
                           showCloseButton = FALSE)
            
            output$sample_selection <- renderUI({
              prettyCheckboxGroup(inputId = "samples_check",
                                  label = "",
                                  inline = TRUE,
                                  status = "info",
                                  choices = samples,
                                  selected = samples)
              
            })
          } else {
            sendSweetAlert(session = session,
                           title = "Data Already Added",
                           text = "The data already exists",
                           type = "info")
          }
        } else {
          # Show user 
          sendSweetAlert(session = session,
                         title = "Failed",
                         text = "Data failed to load. Please make sure you have an internet connection. Otherwise try relaunching WADE.",
                         type = "error")
        }
      }
    }
    
  })
  
  observeEvent(input$select_all, {
    samples <- downloaded_samples.df %>% pull(filename)
    updatePrettyCheckboxGroup(session = session,
                              inputId = "samples_check",
                              selected = samples)
  })
  
  observeEvent(input$deselect_all, {
    samples <- "NA"
    updatePrettyCheckboxGroup(session = session,
                              inputId = "samples_check",
                              selected = samples)
  })
  
  # Home tab ####
  ## Change the locus output depending on user input
  getLocus <- reactive({
    locus <- "list"
    if(input$locus != "list"){
      locus <- input$user_locus
    }
    locus
  })
  
  downloadedDim <- reactive({
    dimensions <- dim(downloaded_samples.df)
    dimensions
  })
  
  getSamples <- reactive({ # getSamples() will return a data frame
    samples <- "" # start as empty
    
    if(!is.null(input$samples_check) && !purrr::is_empty(input$samples_check)){
      if(downloadedDim()[1] > 0 && downloadedDim()[2] > 0){ # Is there actually data that we can use?
        samples <- dplyr::filter(downloaded_samples.df, filename %in% input$samples_check)
      }
    }
    samples
  })
  
  output$selected_Org <- renderText({
    paste("Organism: ", input$org_tab_box)
  })
  
  output$selected_test <- renderText({
    test_name <- paste(input$org_tab_box, "_test", sep = "")
    paste("Test:", input[[test_name]])
  })
  
  output$entered_locus <- renderText({
    paste("Locus: ", getLocus())
  })
  
  # Output Tab ####
  output$download_data <- downloadHandler(
    filename = function() { paste("test-", Sys.Date(), ".", input$download_type, sep = "") },
    content = function(file) {
      sep <- switch(input$download_type, "csv" = ",", "tsv" = "\t")
      write.table(x = final_out.df,
                  file,
                  row.names = FALSE,
                  sep = sep,
                  qmethod = "escape")
    }
  )
  
  #------------------
  # CONFIRMATION CHECK
  observeEvent(input$confirm_execution, {
    if(isTRUE(input$confirm_execution)){
      execute()
    }
  })
  
  #---------------------
  # createDownloadButton
  #
  # Creates a download button
  createDownloadButton <- function(){
    output$download <- renderUI({ # Download Button ####
      box(id = "download_table",
          radioButtons(inputId = "download_type", 
                       label = "File type:",
                       choices = c("csv", "tsv")),
          downloadBttn(outputId = "download_data",
                       label = "Download",
                       style = "material-flat",
                       size = "lg")
      )
    })
  }
  
  #------------------
  # createOutputTable
  #
  # Renders a data a given data table
  createOutputTable <- function(output.df){
    renderDataTable(output.df,
                    options = list(scrollX = TRUE,
                                   pageLength = 10))
  }
  
  #---------------
  # createFilters
  #
  # render a UI object containing the filters
  # of a given data frame
  createFilters <- function(output.df){
    output$filter <- renderUI({ # Filter Boxes
      output_cols <- colnames(output.df)
      
      if(output_cols %>% when(str_detect(., "_result")) %>% any()){
        filter_header <- output_cols[grepl("_result", output_cols)]
        filter_header <- gsub("_result", "", filter_header) %>%
          prepend(output_cols[1]) %>%
          append(output_cols[length(output_cols)])
      } else {
        filter_header <- output_cols
      }
      
      box(id = "out_filter",
          prettyCheckboxGroup(inputId = "filter",
                              label = "Filter",
                              status = "info",
                              inline = TRUE,
                              choices = filter_header # Make the filters the column names
          )
      )
    })
  }
  
  #-------------
  # execute
  #
  # The "main" call.
  # switches the session to the output tab to display output.
  #
  # execute() is called from the confirm_execution observeEvent. Only called when
  # confirm is selected.
  execute <- function(){
    updateTabItems(session = session, # After "execute!" is pressed, change tabs to view the output
                   inputId = "side_tabs",
                   selected = "output")
    
    locus <- getLocus()
    samples <- getSamples() # We need to make sure that sample isn't empty.
    #print(samples)
    org <- input$org_tab_box
    test <- input[[paste(input$org_tab_box, "_test", sep = "")]] # Does this cause an issue?
    
    output.df <- execute_analysis(org, samples, locus, test) # EXECUTE ANALYSIS
    
    
    output$profile_table <- createOutputTable(output.df) # Update the output tab
    
    final_out.df <<- output.df # Update our global data frame
    
    # Output ####
    filename <- paste("output_profile_", test, ".csv", sep="")
    file_out <<- here("data", "output", filename)
    write.csv(x = output.df, # Write file
              file = file_out,
              row.names = FALSE)
    
    GalaxyConnector::gx_put(file_out, filename) # Writes the file to Galaxy
    
    createDownloadButton()
    createFilters(output.df)
  }
  
  observeEvent(input$execute, {
    # If there is no data selected
    if(is.null(getSamples()[1,1]) || getSamples()[1,1] == ""){
      shinyWidgets::sendSweetAlert(session = session,
                                   title = "No Data Selected",
                                   text = "Please select some data and try again",
                                   type = "warning",
                                   btn_labels = "Ok")
    } else {
      shinyWidgets::confirmSweetAlert(session = session,
                                      inputId = "confirm_execution",
                                      title = "Confirm Execution?",
                                      text = "Once confirmed it cannot be stopped!",
                                      type = "question",
                                      btn_labels = c("Cancel", "Confirm"))
    }
  }) # End Main Call ####
  
  # Output Filter ####
  observeEvent(input$filter, {
    vals <- input$filter # This is a list
    headers <- colnames(x = final_out.df) # Headers of the data table
    filters <- unlist(vals) %>% map(~ grep(.x, headers)) %>% unlist() # Get the columns that contain the vals from the filter
    
    if(!is.null(vals)){ # Are there filters selected??
      output$profile_table <- renderDataTable(select(final_out.df, filters), # Do a select on the table
                                              options = list(scrollX = TRUE, # This render can probably be made into a small module
                                                             pageLength = 10))
    } else { # No filters, just do regular
      output$profile_table <- renderDataTable(final_out.df,
                                              options = list(scrollX = TRUE,
                                                             pageLength = 10))
    }
  }, ignoreNULL = FALSE) # Needed to detect any deselection to NULL
  
  # Called from Main server call
  # Calls the correct analysis
  #
  # org: string
  # samples: dataframe of Filename and the path
  # locus: string
  # test: string
  execute_analysis <- function(org, samples, locus, test){
    withProgress(message = "Work in progress...",
                 value = 0, {
                   switch(test,
                          AMR_DB = { database_pipeline(org, samples, FALSE) },
                          AMR_LW = { labware_gono_amr() },
                          EMM = { emm(org, samples, locus) },
                          MASTER = { master_blastr(org, test, samples, locus) },
                          MLST = { general_mlst_pipeline(org, samples, locus, test) },
                          NGSTAR = { general_mlst_pipeline(org, samples, locus, test) },
                          NGMAST = { general_mlst_pipeline(org, samples, locus, test) },
                          rRNA23S = { rna_23s(org, samples) },
                          SERO = { PneumoCaT_pipeline(samples) },
                          VFDB = { database_pipeline(org, samples, TRUE) },
                          { master_blastr(org, test, samples, locus) }
                   )
                 })
  }

  # Settings Tab ####
  observeEvent(input$make_blast_db, { Index_pipeline(org, test, locus) })
}
