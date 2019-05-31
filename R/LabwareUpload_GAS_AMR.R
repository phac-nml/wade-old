#' Run AMR first
#' Then run this analysis to combine data the full amr profile to upload to LabWare.
#'
#' @return A table frame containing the results of the query
#' @export
#'

labware_gas_amr <- function() {

  test_file <- here("data", "input", "labware", "output_profile_GAS_AMR.csv")

  # Initialize DF's ####
  Output.df <- read.csv(file = test_file,
                        header = TRUE,
                        sep = ",",
                        stringsAsFactors = FALSE)

  # Variables ####
  sample_num <- (dim(Output.df))[1]
  delim_dash <- " - " # delimiter
  m <- 1 # More descriptive name

  # Iteration ####
  for (m in 1L:sample_num)
  {
    # Variables ####
    molec_profile <- ""
    pos = "POS"
    lw_CurrSampleNo <- as.character(Output.df[m, "SampleNo"])

    # Similar Vars ####
    lw_ermA <- as.character(Output.df[m, "ermA_result"]) # What is ermA??
    lw_ermB <- as.character(Output.df[m, "ermB_result"])
    lw_ermT <- as.character(Output.df[m, "ermT_result"])
    lw_mefAE <- as.character(Output.df[m, "mefAE_result"])
    lw_tetM <- as.character(Output.df[m, "tetM_result"])
    lw_tetO <- as.character(Output.df[m, "tetO_result"])
    lw_cat <- as.character(Output.df[m, "cat_result"])
    lw_dfrF <- as.character(Output.df[m, "dfrF_result"])
    lw_dfrG <- as.character(Output.df[m, "dfrG_result"])

    # Consecutive IF's ####
    if (lw_ermA == pos) { molec_profile <- get_labware_profile(molec_profile == "", molec_profile, "ermA", delim_dash) }

    if (lw_ermB == pos) { molec_profile <- get_labware_profile(molec_profile == "", molec_profile, "ermB", delim_dash) }

    if (lw_ermT == pos) { molec_profile <- get_labware_profile(molec_profile == "", molec_profile, "ermT", delim_dash) }

    if (lw_mefAE == pos) { molec_profile <- get_labware_profile(molec_profile == "", molec_profile, "mefAE", delim_dash) }

    if (lw_tetM == pos) { molec_profile <- get_labware_profile(molec_profile == "", molec_profile, "tetM", delim_dash) }

    if (lw_tetO == pos) { molec_profile <- get_labware_profile(molec_profile == "", molec_profile, "tetO", delim_dash) }

    if (lw_cat == pos){ molec_profile <- get_labware_profile(molec_profile == "", molec_profile, "cat", delim_dash) }

    if (lw_dfrF == pos) { molec_profile <- get_labware_profile(molec_profile == "", molec_profile, "dfrF", delim_dash) }

    if (lw_dfrG == pos) { molec_profile <- get_labware_profile(molec_profile == "", molec_profile, "dfrG", delim_dash) }

    #--------------------------------

    lw_gyrA <-   as.character(Output.df[m, "gyrA_mutations"])


    if (lw_gyrA == "???" | lw_gyrA == "x" | lw_gyrA == "" | is.na(lw_gyrA)){
      lw_gyrA <- "Err"
    }
    if (lw_gyrA != "WT")
    {
      if (molec_profile == "")
      {molec_profile <- paste("gyrA ", lw_gyrA, sep = "")}else
      {molec_profile <- paste(molec_profile, delim_dash, "gyrA ", lw_gyrA, sep = "")}
    }

    #---------------------------------

    lw_parC <-   as.character(Output.df[m, "parC_mutations"])

    if (lw_parC == "???" | lw_parC == "x" | lw_parC == "" | is.na(lw_parC)) {
      lw_parC <- "Err/Err/Err"
    }
    if (lw_parC != "WT")  # Everything is here is similar to everything else....
    {
      lw_parC_parts <- unlist(strsplit(lw_parC, "/"))
      lw_parC_D78 <-  lw_parC_parts[1]
      lw_parC_S79 <- lw_parC_parts[2]
      lw_parC_D83 <- lw_parC_parts[3]
      lw_parC_prof <- NA

      if (lw_parC_D78 != "WT")
      {
        if (is.na(lw_parC_prof))
        {lw_parC_prof <- paste("parC ", lw_parC_D78, sep = "")} else
        {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_D78, sep = "")}
      }
      if (lw_parC_S79 != "WT")
      {
        if (is.na(lw_parC_prof))
        {lw_parC_prof <- paste("parC ", lw_parC_S79, sep = "")} else
        {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_S79, sep = "")}
      }
      if (lw_parC_D83 != "WT")
      {
        if (is.na(lw_parC_prof))
        {lw_parC_prof <- paste("parC ", lw_parC_D83, sep = "")} else
        {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_D83, sep = "")}
      }

      if (!is.na(lw_parC_prof))
      {
        if (molec_profile == "")
        {molec_profile <- lw_parC_prof} else
        {molec_profile <- paste(molec_profile, delim_dash, lw_parC_prof, sep = "")}
      }

    }

    lw_folA <-   as.character(Output.df[m, "folA_mutations"])
    if (lw_folA == "???" | lw_folA == "x" | lw_folA == "" | is.na(lw_folA)) {lw_folA <- "Err"}
    if (lw_folA != "WT")
    {
      if (molec_profile == "")
      {molec_profile <- paste("folA ", lw_folA, sep = "")}else
      {molec_profile <- paste(molec_profile, delim_dash, "folA ", lw_folA, sep = "")}
    }

    lw_folP <-   as.character(Output.df[m, "folP_mutations"])
    if (lw_folP == "???" | lw_folP == "x" | lw_folP == "" | is.na(lw_folP)) {lw_folP <- "Err"}
    if (lw_folP != "WT")
    {
      if (molec_profile == "")
      {molec_profile <- paste("folP ", lw_folP, sep = "")}else
      {molec_profile <- paste(molec_profile, delim_dash, "folP ", lw_folP, sep = "")}
    }

    if (molec_profile == "") {molec_profile <- "Wild Type"}


    #----------------------------------------------------------------------  INTERPRETATIONS
    amr_profile <- "Susceptible"

    if (str_detect(molec_profile, paste(c("ermA", "ermB", "ermT", "mefAE"),collapse = '|')))
    {
      ery <- "Resistant"
    }else{ery <- "Susceptible"}

    if (str_detect(molec_profile, "cat"))
    {
      chl <- "Resistant"
    }else{chl <- "Susceptible"}

    if (str_detect(molec_profile, "ermB"))
    {
      cli <- "Resistant"
    }else if (str_detect(molec_profile, paste(c("ermA", "ermT"),collapse = '|')))
    {
      cli <- "Inducible"
    } else {cli <-"Susceptible"}

    # Variables ####
    cip <- "Undetermined"
    lev <- "Undetermined"

    if ( (str_detect(molec_profile, "gyrA Err")) | (str_detect(molec_profile, "parC Err")) )
    {cip <- "Error"
    amr_profile <- "Error"
    }else
    {
      if (str_detect(molec_profile, paste(c("gyrA", "parC"),collapse = '|')))
      {
        cip <- "Resistant"
      }else{cip <- "Susceptible"}

      if (str_detect(molec_profile, "gyrA S81F"))
      {
        lev <- "Resistant"
      }else{lev <- "Susceptible"}
    }


    if (str_detect(molec_profile, paste(c("tetM", "tetO"),collapse = '|')))
    {
      tet <- "Resistant"
    }else{tet <- "Susceptible"}

    if ( (str_detect(molec_profile, "folA Err")) | (str_detect(molec_profile, "folP Err")) )
    {sxt <- "Error"
    amr_profile <- "Error"
    } else { # else if
      if (str_detect(molec_profile, paste(c("dfrG", "dfrF", "folA", "folP"),collapse = '|'))) # this can be an else if
      {

        if (str_detect(molec_profile, paste(c("dfrG", "dfrF"),collapse = '|')))
        {
        sxt <- "Resistant"
        } else if (str_detect(molec_profile, paste(c("folA", "folP"),collapse = '|')))
        {
          sxt <- "Intermediate"
          if ((str_detect(molec_profile, "folA")) & (str_detect(molec_profile, "folP")))
          {
          sxt <- "Resistant"
          }
        }
      }else{sxt <- "Susceptible"} # the else, after the if and else if

    }

    #--------------------------------------------------------------  MAKE AMR PROFILE
    if (amr_profile != "Error") {

      amr_profile <- "Susceptible"
      delim_slash <- "/"

      ### Consecutive IF's ###
      if (ery == "Resistant") { amr_profile <- get_labware_profile(amr_profile == "Susceptible", amr_profile, "ERY-R", delim_slash) }

      if (cli == "Resistant") { amr_profile <- get_labware_profile(amr_profile == "Susceptible", amr_profile, "CLI-R", delim_slash) }

      if (cli == "Inducible") { amr_profile <- get_labware_profile(amr_profile == "Susceptible", amr_profile, "CLI-Ind", delim_slash) }

      if (chl == "Resistant") { amr_profile <- get_labware_profile(amr_profile == "Susceptible", amr_profile, "CHL-R", delim_slash) }

      if (cip == "Resistant") { amr_profile <- get_labware_profile(amr_profile == "Susceptible", amr_profile, "CIP-R", delim_slash) }

      if (lev == "Resistant") { amr_profile <- get_labware_profile(amr_profile == "Susceptible", amr_profile, "LEV-R", delim_slash) }

      if (tet == "Resistant") { amr_profile <- get_labware_profile(amr_profile == "Susceptible", amr_profile, "TET-R", delim_slash) }

      if (sxt == "Resistant") { amr_profile <- get_labware_profile(amr_profile == "Susceptible", amr_profile, "SXT-R", delim_slash) }

      if (sxt == "Intermediate") {
        amr_profile <- get_labware_profile(amr_profile == "Susceptible", amr_profile, "SXT-I", delim_slash)
        # if (amr_profile == "Susceptible")
        # {amr_profile <- "SXT-I"} else
        # {amr_profile <- paste(amr_profile, sepr2, "SXT-I", sep = "")}
      }

    }

    #--------------------------------------------------------------

    labware_sample.df <- data.frame(
                                    lw_CurrSampleNo,
                                    "ND",
                                    "ND",
                                    lw_ermA,
                                    lw_ermB,
                                    lw_ermT,
                                    lw_mefAE,
                                    lw_gyrA,
                                    lw_parC,
                                    lw_tetM,
                                    lw_tetO,
                                    lw_cat,
                                    lw_dfrF,
                                    lw_dfrG,
                                    lw_folA,
                                    lw_folP,
                                    molec_profile,
                                    ery,
                                    chl,
                                    cip,
                                    cli,
                                    tet,
                                    sxt,
                                    amr_profile
                                    )
    headers <- list("SampleNo", "A2059G", "C2611T", "ermA", "ermB", "ermT", "mefAE", "gyrA", "parC","tetM", "tetO",
                    "cat", "dfrF", "dfrG", "folA", "folP", "Molec Profile",  "ERY", "CHL", "CIP", "CLI", "TET", "SXT", "AMR Profile")

    names(labware_sample.df) <- headers

    if(m==1) {  #if first sample make one row profile table, otherwise add new row to table
      labware.df <- data.frame(labware_sample.df)
      names(labware.df) <- headers
    } else {
      labware.df <- bind_rows(labware.df, labware_sample.df) # How does labware.df work here if it DNE?
    }

  } # End for loop

  # Output ####
  write.csv(labware.df, here("data", "output", "LabWareUpload_GAS_AMR.csv") , quote = FALSE,  row.names = FALSE)

  # Return ####
  # return(labware.df)
  labware.df
}
