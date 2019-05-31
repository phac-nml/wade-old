#' Function building dashboard UI, used in Shiny app
#'
#' @export
ui <- dashboardPage(skin = "blue",
                    dashboardHeader(title = "WADE"),
                    dashboardSidebar(
                      tags$head(
                        tags$style(
                          "#execute { width: 90%; margin-left: auto; margin-right: auto }
                          .shiny-notification { width = 100% } ",
                          HTML("dashboard{ margin-bottom: 100%; }
                               .shiny-notification {
                               width: 80%;
                               position:fixed;
                               top: calc(90%);;
                               left: calc(15%);;
                               font-size: 250%;
                               text-align: center;
                               }
                               "))
                      ),

                      sidebarMenu(id = "side_tabs",
                                  # Home Tab ####
                                  menuItem(text = "Home",
                                           tabName = "home",
                                           icon = icon("dashboard")
                                  ),
                                  # Output Tab ####
                                  menuItem(text = "Output",
                                           tabName = "output",
                                           icon = icon("file")
                                  ),
                                  # Settings Tab ####
                                  menuItem(text = "Settings",
                                           tabName = "settings",
                                           icon = icon("cog")
                                  )
                      ),
                      # execute Button ####
                      actionBttn(inputId = "execute",
                                 label = "execute!",
                                 color = "success",
                                 style = "material-flat",
                                 size = "lg",
                                 block = "FALSE")
                    ),

                    dashboardBody(
                      tags$head(tags$style("#selected_output { font-size: 175%; }
                                           #org_tab_box { font-size: 150%; }
                                           #sample_tab_box { font-size: 150%; }
                                           #update-db { width: 49%; margin-left: auto; margin-right: auto }
                                           #make_blast_db { width: 49%; margin-left: auto; margin-right: auto }")),
                      tabItems(
                        # Home tab content ####
                        tabItem(tabName = "home",
                                fluidRow(
                                  # Organism Box ####
                                  tabBox(id = "org_tab_box",
                                         tabPanel("GAS", # GAS Tab panel
                                                  prettyRadioButtons("GAS_test", h3("Choose an analysis:"),
                                                                     choices = sort(c("AMR Profile" = "AMR",
                                                                                      "Toxin Profile" = "TOXINS",
                                                                                      "MLST Type" = "MLST",
                                                                                      "Virulence Factors" = "VIRULENCE",
                                                                                      "EMM Typing" = "EMM",
                                                                                      "MasterBlastR*" = "MASTER",
                                                                                      "ARG-ANNOT/Resfinder/CARD" = "AMR_DB",
                                                                                      "VFDB (Virulence Factor Database)" = "VFDB"
                                                                                      #"LabWare Metrics" = "LW_METRICS"

                                                                     ), decreasing = FALSE),
                                                                     selected = "AMR",
                                                                     animation = "smooth"
                                                  ),
                                                  value = "GAS"
                                         ),
                                         tabPanel("GONO", # GONO Tab Panel
                                                  prettyRadioButtons("GONO_test", h3("Choose an analysis:"),
                                                                     choices = sort(c("AMR Profile" = "AMR",
                                                                                      "NG-STAR Type" = "NGSTAR",
                                                                                      "23S rRNA Alleles" = "rRNA23S",
                                                                                      "LabWare AMR profile*" = "AMR_LW",
                                                                                      "MLST Type" = "MLST",
                                                                                      "NG-MAST Type" = "NGMAST",
                                                                                      "MasterBlastR*" = "MASTER"
                                                                                      #"LabWare Metrics" = "LW_METRICS"

                                                                     ), decreasing = FALSE),
                                                                     selected = "AMR",
                                                                     animation = "smooth"
                                                  ),
                                                  value = "GONO"
                                         ),
                                         tabPanel("PNEUMO", # PNEUMO Tab Panel
                                                  prettyRadioButtons("PNEUMO_test",
                                                                     h3("Choose an analysis:"),
                                                                     choices = sort(c("AMR Profile" = "AMR",
                                                                                      "23S rRNA Alleles" = "rRNA23S",
                                                                                      "Serotyping*" = "SERO",
                                                                                      "MLST Type" = "MLST",
                                                                                      "Virulence Factors" = "VIRULENCE",
                                                                                      "MasterBlastR*" = "MASTER",
                                                                                      "ARG-ANNOT/Resfinder/CARD" = "AMR_DB",
                                                                                      "VFDB (Virulence Factor Database)" = "VFDB"
                                                                                      #"LabWare Metrics" = "LW_METRICS"
                                                                     ), decreasing = FALSE),
                                                                     selected = "SERO",
                                                                     animation = "smooth"
                                                  ),
                                                  value = "PNEUMO"
                                         )
                                  ),
                                  # Loci and Sample Box ####
                                  tabBox(id = "sample_tab_box",
                                         tabPanel("Select Data",
                                                  dataInputModuleUI("data_input")),
                                         tabPanel("Loci",
                                                  # Select Locus ####
                                                  prettyRadioButtons(inputId = "locus",
                                                                     label = h3("Select a locus to query"),
                                                                     choices = list("Default loci list" = "list",
                                                                                    "Input locus query" = "input_loci"),
                                                                     selected = "list",
                                                                     animation = "smooth"

                                                  ),
                                                  conditionalPanel( # selecting to input a loci will show the text input box
                                                    condition = "input.locus == 'input_loci'",
                                                    textInput(inputId = "user_locus", label = "Enter a locus to query", value = "")
                                                    # Should this maybe be textAreaInput? That way they can expand it if they wish?
                                                  )
                                         )
                                  )
                                ), # End fluidRow

                                fluidRow(
                                  box(id = "selected_output",
                                      textOutput("selected_Org"),
                                      textOutput("selected_test"),
                                      textOutput("entered_locus"),
                                      textOutput("entered_sample"),
                                      textOutput("galaxy_access")
                                  ),
                                  box(id = "samples",
                                      collapsible = TRUE,
                                      collapsed = TRUE,
                                      title = h2("Samples"),
                                      uiOutput("sample_selection")
                                  )
                                )
                        ), # End Home tab

                        # Output tab content
                        tabItem(tabName = "output",
                                dataTableOutput(outputId = "profile_table"),
                                uiOutput("output_title"),
                                uiOutput("download"),
                                uiOutput("filter")
                        ),

                        # Settings tab content
                        tabItem(tabName = "settings",
                                h2("Settings"),
                                h3("Database Versions:"),
                                actionBttn(inputId = "update-db",
                                           label = "Update Databases",
                                           icon = shiny::icon("exchange-alt"),
                                           color = "primary",
                                           style = "material-flat",
                                           block = "FALSE"
                                ),
                                actionBttn(inputId = "make_blast_db",
                                           label = "Make BlastDB",
                                           icon = shiny::icon("database"),
                                           color = "primary",
                                           style = "material-flat",
                                           block = "FALSE"
                                )
                        )
                      )
                    )
)
