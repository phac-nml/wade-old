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

  # Home tab ####

  ## Change the locus output depending on user input
  getLocus <- reactive({
    locus <- "list"
    if(input$locus != "list"){
      locus <- input$user_locus
    }
    locus
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

  # Sample Selection Checkboxes
  output$sample_selection <- renderUI({
    files <- list.files(path = here("data", "databases", input$org_tab_box, "assemblies"),
                        pattern = "*.fasta") # Grab assemblies of selected org
    chosen <- "NULL"
    if(input$sample == "list"){ # Is the sample list selected?
      chosen <- files # Update the chosen checkboxes
    }

    prettyCheckboxGroup(inputId = "samples_check",
                        label = "Select a Sample...",
                        inline = TRUE,
                        status = "info",
                        choices = files,
                        selected  = chosen)
  })

  observeEvent(input$samples_check, {
    samples <- input$samples_check
    if(is.null(samples)){
      writeLines("NULL")
    } else {
      print(input$samples_check)
    }
  })

  # Output Tab ####
  output$download_data <- downloadHandler(
    filename = function() { paste("test-", Sys.Date(), ".csv", sep = "") },
    content = function(file) {
      print(file_out)
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
    sample <- "list" #user_sample()
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
