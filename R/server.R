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
  
  names <- c("subdir_id", "filename")
  downloaded_fullpath <- ""
  downloaded_samples.df <- data.frame(subdir_id, filename)
  
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

      all_data <- dplyr::filter(GalaxyConnector::gx_list_history_datasets(), deleted == FALSE) # Filter out any deleted dataset
      data_hid.df <- dplyr::filter(all_data, name == input$select_dataset) # Grab the data that's been selected

      # Take the first piece of data
      #   This could possibly be an issue since with collections there are multiple versions of the same data (by name).
      #   If one of the pieces of data we need are changed then things can possibly go wrong
      
      # Let's check if the data is already downloaded.
      datapath <- GalaxyConnector::gx_get(data_hid.df[1, 'hid'])
      
      if(!is.null(datapath)){
        if(downloaded_fullpath == ""){ # Grab the path earlier so we can use tidy data
          downloaded_fullpath <- base::basename(base::dirname(datapath))
        }
        
        downloaded_samples.df <- addDirFiles(base::dirname(datapath), downloaded_samples.df) # DATAFRAME OF ALL THE CONCATENATED DATA
        
        downloaded_data <- list.files(base::dirname(datapath)) # Take this and make it look really nice.
        output$sample_selection <- renderUI({
          prettyCheckboxGroup(inputId = "samples_check",
                              label = "",
                              inline = TRUE,
                              status = "info",
                              choices = downloaded_data,
                              selected = downloaded_data)
        })
      } else {
        # RAISE ERROR ####
      }
    }
    
  })
  
  addDirFiles <- function(dir, df){
    subdir_id <- basename(dir) # Will be a number representing the data/collection number from Galaxy
    all_files <- list.files(dir) # All the files in the subdir
    
    dir_files.df <- data.frame(subdir_id, all_files)
    names(dir_files.df) <- names
    df <- rbind(df, dir_files.df)
    df # return df!
  }
  

  # Home tab ####
  ## Change the locus output depending on user input
  getLocus <- reactive({
    locus <- "list"
    if(input$locus != "list"){
      locus <- input$user_locus
    }
    locus
  })
  
  getSamples <- reactive({
    samples <- "" # start as empty
    if(!is.null(input$samples_check) && !purrr::is_empty(input$samples_check)){
      samples <- input$samples_check
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
    filename = function() { paste("test-", Sys.Date(), ".csv", sep = "") },
    content = function(file) {
      write.table(x = read_table(file_out),
                  file)
    }
  )

  # Main Call ####
  observeEvent(input$execute, {
    updateTabItems(session = session, # After "execute!" is pressed, change tabs to view the output
                   inputId = "side_tabs",
                   selected = "output")

    locus <- getLocus()
    sample <- getSamples() # We need to make sure that sample isn't empty.
    org <- input$org_tab_box
    test <- input[[paste(input$org_tab_box, "_test", sep = "")]] # Does this cause an issue?

    output.df <- execute_analysis(org, sample, locus, test)
    output$profile_table <- renderDataTable(output.df,
                                            options = list(scrollX = TRUE,
                                                           pageLength = 10)
    )

    final_out.df <<- output.df

    # Output ####
    file_out <<- here("data", "output", paste("output_profile_", test, ".csv", sep = "")) # Do we need this for downloading???
    write.csv(x = output.df,
              file = file_out,
              row.names = F)



    # Change Output UI ####
    output$download <- renderUI({ # Download Button ####
      box(id = "download_table",
          radioButtons("download_type", "File type:",
                       choices = c("csv", "tsv")),
          downloadBttn(outputId = "download_data",
                       label = "Download",
                       style = "material-flat",
                       size = "lg")
      )
    })

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
  # sample: dataframe of Filename and the path
  # locus: string
  # test: string
  execute_analysis <- function(org, sample, locus, test){
    withProgress(message = "Work in progress...",
                 value = 0, {

                   switch(test,
                          AMR_DB = { database_pipeline(org, sample, FALSE) },
                          AMR_LW = { labware_gono_amr() },
                          EMM = { emm(org, sample, locus) },
                          MASTER = { master_blastr(org, test, sample, locus) },
                          MLST = { general_mlst_pipeline(org, sample, locus, test) },
                          NGSTAR = { general_mlst_pipeline(org, sample, locus, test) },
                          NGMAST = { general_mlst_pipeline(org, sample, locus, test) },
                          rRNA23S = { rna_23s(org, sample) },
                          SERO = { PneumoCaT_pipeline(sample) },
                          VFDB = { database_pipeline(org, sample, TRUE) },
                          { master_blastr(org, test, sample, locus) }
                   )
                 })
  }

  # Settings Tab ####
  observeEvent(input$make_blast_db, { Index_pipeline(org, test, locus) })
}
