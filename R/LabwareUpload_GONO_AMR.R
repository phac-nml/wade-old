#' Run AMR first, then run the 23S allele counts, and then the NGSTAR-MLST analyses
#' Then run this analysis to combine data from AMR, 23S rRNA and NG-STAR
#' to prepare full amr profile to upload to LabWare.
#'
#' @return A table frame containing the results of the query
#' @export
#'
#'

labware_gono_amr <- function(amrDF, ngstarDF, rnaDF) {
  # TODO ####
  # - Line 104, The comparator is NA, we should handle this case for everything
  
  AMR_Output.df <- amrDF
  NGSTAR_Output.df <- ngstarDF
  rRNA23S_Output.df <- rnaDF
  
    
  # Initialize DF's ####
  # AMR_Output.df <- read.csv(here("data", "input", "labware", "output_profile_GONO_AMR.csv"),
  #                           header = TRUE,
  #                           sep = ",",
  #                           stringsAsFactors = FALSE)
  # NGSTAR_Output.df <- read.csv(here("data", "output", "GONO", "output_profile_mut.csv"),
  #                              header = TRUE,
  #                              sep = ",",
  #                              stringsAsFactors = FALSE)
  # rRNA23S_Output.df <- read.csv(here("data", "output", "output_profile_23S.csv"),
  #                               header = TRUE,
  #                               sep = ",",
  #                               stringsAsFactors = FALSE)
  rRNA23S_Output.df <- rRNA23S_Output.df %>% mutate(SampleNo = as.character(SampleNo)) # Convert the SampleNo in rRNA23S to a character

  Combined_Output_first.df <- full_join(AMR_Output.df, NGSTAR_Output.df, by = "SampleNo")
  Combined_Output.df <- full_join(Combined_Output_first.df, rRNA23S_Output.df, by = "SampleNo")

  NumSamples <- (dim(Combined_Output.df))[1]

  # Variables ####
  delim_dash <- " - "
  delim_slash <- "/"
  m <- 1
  ErrorFound <- FALSE

  # Iteration ####
  for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Build Molecular Profile
  {
    molec_profile <- NA

    lw_CurrSampleNo <- as.character(Combined_Output.df[m, "SampleNo"])

    # -------- ermB --------
    lw_ermB <- as.character(Combined_Output.df[m, "ermB_result"])
    if (lw_ermB == "POS")
    {
      molec_profile <- "ermB"
      V_ermB <- 1L
    } else (v_ermB <- 0L)

    # -------- ermC --------
    lw_ermC <- as.character(Combined_Output.df[m, "ermC_result"])
    if (lw_ermC == "POS") {
      v_ermC <- 1L
      molec_profile <- get_labware_profile(compare = is.na(molec_profile),
                                           profile = molec_profile,
                                           gene = "ermC",
                                           delimiter = delim_dash)
    } else { v_ermC <- 0L }

    # -------- rpsJ --------
    lw_rpsJ <- as.character(Combined_Output.df[m, "rpsJ_mutations"])
    if (is.na(lw_rpsJ)) {
      lw_rpsJ <- "Err"
    }

    if (lw_rpsJ != "WT")
    {
      if (is.na(molec_profile)) {
        molec_profile <- paste("rpsJ ", lw_rpsJ, sep = "")
      } else {
        molec_profile <- paste(molec_profile, delim_dash, "rpsJ ", lw_rpsJ, sep = "")
      }
    }

    # -------- tetM --------
    lw_tetM <- as.character(Combined_Output.df[m, "tetM_result"])
    if (lw_tetM == "POS") {
      molec_profile <- get_labware_profile(is.na(molec_profile), molec_profile, "tetM", delim_dash)
    }

    # -------- bla--------
    lw_bla <- as.character(Combined_Output.df[m, "bla_result"])
    if (lw_bla == "POS") {
      v_bla <- 1L
      molec_profile <- get_labware_profile(is.na(molec_profile), molec_profile, "bla", delim_dash)
    } else {v_bla <- 0L}

    #--------------------------------------------------------------------------------------gyrA
    # Variables ####
    lw_16S <-   as.character(Combined_Output.df[m, "X16S_mutations"])

    if (lw_16S == "?" | lw_16S == "x") {
      lw_16S <- "Err/Err"
    }

    lw_16S_parts <- unlist(strsplit(lw_16S, "/"))
    lw_16S_C1192T <-  lw_16S_parts[1]
    lw_16S_T1458C <- lw_16S_parts[2]
    lw_16S_prof <- NA

    if (lw_16S_C1192T != "WT")
    {
      if (is.na(lw_16S_prof))
      {lw_16S_prof <- paste("16S rRNA ", lw_16S_C1192T, sep = "")} else
      {lw_16S_prof <- paste(lw_16S_prof, "/", lw_16S_C1192T, sep = "")}
    }
    if (lw_16S_T1458C != "WT")
    {
      if (is.na(lw_16S_prof))
      {lw_16S_prof <- paste("23S rRNA ", lw_16S_T1458C, sep = "")} else
      {lw_16S_prof <- paste(lw_16S_prof, "/", lw_16S_T1458C, sep = "")}
    }

    if (!is.na(lw_16S_prof)) {
      molec_profile <- get_labware_profile(is.na(molec_profile), molec_profile, lw_16S_prof, delim_dash)
    }

    #---------------------------------------------------------------------------  NG-STAR penA
    # Variables ####
    lw_penA <-   as.character(Combined_Output.df[m, "penA"])
    lw_penA_prof <- NA
    if (lw_penA == "?" | lw_penA == "x") {
      lw_penA <- "Err/Err/Err/Err/Err"
    }

    lw_penA_parts <- unlist(strsplit(lw_penA, "/"))
    lw_penA_A311V <-  lw_penA_parts[1]
    lw_penA_A501 <- lw_penA_parts[2]
    lw_penA_N513Y <- lw_penA_parts[3]
    lw_penA_A517G <- lw_penA_parts[4]
    lw_penA_G543S <- lw_penA_parts[5]

    # -------- penA_A311V --------
    if (lw_penA_A311V != "WT") { lw_penA_prof <- get_labware_profile(compare = is.na(lw_penA_prof),
                                                                     profile = lw_penA_prof,
                                                                     gene = lw_penA_A311V,
                                                                     delimiter = delim_slash)
    }

    # -------- penA_A501 --------
    if (lw_penA_A501 != "WT") {
      lw_penA_prof <- get_labware_profile(is.na(lw_penA_prof), lw_penA_prof, lw_penA_A501, delim_slash)
    }

    # -------- penA_N513Y --------
    if (lw_penA_N513Y != "WT") {
      lw_penA_prof <- get_labware_profile(is.na(lw_penA_prof), lw_penA_prof, lw_penA_N513Y, delim_slash)
    }

    # -------- penA_A517G --------
    if (lw_penA_A517G != "WT") {
      lw_penA_prof <- get_labware_profile(is.na(lw_penA_prof), lw_penA_prof, lw_penA_A517G, delim_slash)
    }

    # -------- penA_G543S --------
    if (lw_penA_G543S != "WT") {
      lw_penA_prof <- get_labware_profile(is.na(lw_penA_prof), lw_penA_prof, lw_penA_G543S, delim_slash)
    }

    # -------- penA_prof --------
    if (is.na(lw_penA_prof))
    {
      lw_penA_prof <- "WT"
      if (is.na(molec_profile))
      {molec_profile <- lw_penA_prof} else
      {molec_profile <- paste(molec_profile, delim_dash, lw_penA_prof, sep = "")}
    }

    ### TOO MUCH! CONFUSING! ###
    ifelse(str_detect(lw_penA_prof, "A311V"), (v_penA_A311V <- 1L), (v_penA_A311V <- 0L))
    ifelse(str_detect(lw_penA_prof, "A501P"), (v_penA_A501P <- 1L), (v_penA_A501P <- 0L))
    ifelse(str_detect(lw_penA_prof, "A501V"), (v_penA_A501V <- 1L), (v_penA_A501V <- 0L))
    ifelse(str_detect(lw_penA_prof, "A501T"), (v_penA_A501T <- 1L), (v_penA_A501T <- 0L))
    ifelse(str_detect(lw_penA_prof, "N513Y"), (v_penA_N513Y <- 1L), (v_penA_N513Y <- 0L))
    ifelse(str_detect(lw_penA_prof, "A517G"), (v_penA_A517G <- 1L), (v_penA_A517G <- 0L))
    ifelse(str_detect(lw_penA_prof, "G543S"), (v_penA_G543S <- 1L), (v_penA_G543S <- 0L))

    #---------------------------------------------------------------------------  mtrR
    # Variables ####
    lw_mtrR <-   as.character(Combined_Output.df[m, "mtrR"])
    lw_mtrR_prof <- NA

    if (lw_mtrR == "?" | lw_mtrR == "x") {lw_mtrR <- "Err/Err/Err"}

    lw_mtrR_parts <- unlist(strsplit(lw_mtrR, "/"))
    lw_mtrR_p <-  lw_mtrR_parts[1]
    lw_mtrR_A39 <- lw_mtrR_parts[2]
    lw_mtrR_G45 <- lw_mtrR_parts[3]

    # -------- mtrR_p --------
    if (lw_mtrR_p != "WT")
    {
      if (is.na(lw_mtrR_prof))
      {lw_mtrR_prof <- paste("mtrR ", lw_mtrR_p, sep = "")} else
      {lw_mtrR_prof <- paste(lw_mtrR_prof, "mtrR ", lw_mtrR_p, sep = "")}

    }

    # -------- mtrR_A39 --------
    if (lw_mtrR_A39 == "A39T")
    {
      if (is.na(lw_mtrR_prof))
      {lw_mtrR_prof <- "mtrR A39T"} else
      {lw_mtrR_prof <- paste(lw_mtrR_prof, "/A39T", sep = "")}
    }

    # -------- mtrR_G45 --------
    if (lw_mtrR_G45 == "G45D")
    {
      if (is.na(lw_mtrR_prof))
      {lw_mtrR_prof <- "mtrR G45D"} else
      {lw_mtrR_prof <- paste(lw_mtrR_prof, "/G45D", sep = "")}
    }

    # -------- mtrR_prof --------
    if (!is.na(lw_mtrR_prof))
    {
      if (is.na(molec_profile))
      {molec_profile <- lw_mtrR_prof} else
      {molec_profile <- paste(molec_profile, delim_dash, lw_mtrR_prof, sep = "")}
    }

    ### TOO MUCH! CONFUSING! ###
    if (is.na(lw_mtrR_prof)) {lw_mtrR_prof <- ""}
    if (str_detect(lw_mtrR_prof, "-35Adel")){v_mtrR35Adel <- 1L}else {v_mtrR35Adel <- 0L}
    if (str_detect(lw_mtrR_prof, "MEN")){v_mtrRMEN <- 1L}else {v_mtrRMEN <- 0L}
    if (str_detect(lw_mtrR_prof, "Disrupted")){v_mtrR35Disrupted <- 1L}else {v_mtrRDisrupted <- 0L}
    if (str_detect(lw_mtrR_prof, "G45D")){v_mtrR_G45D <- 1L}else {v_mtrR_G45D <- 0L}



    #--------------------------------------------------------------------------------------porB
    lw_porB <-   as.character(Combined_Output.df[m, "porB"])
    if (lw_porB == "?" | lw_porB == "x")
      {
        lw_porB <- "Err/Err"
        lw_porB_struct <- "Err"
        lw_porB_G120 <- "Err"
        lw_porB_A121 <- "Err"
      } else {
        lw_porB_parts <- unlist(strsplit(lw_porB, "/"))

        if (lw_porB_parts[1] == "porB1a")
        {
          lw_porB_struct <- lw_porB_parts[1]
          lw_porB_G120 <- NA
          lw_porB_A121 <- NA
        } else {
          lw_porB_struct <- "porB1b"
          lw_porB_G120 <- lw_porB_parts[1]
          lw_porB_A121 <- lw_porB_parts[2]
        }
      }

    lw_porB_prof <- NA

    if (lw_porB_struct == "porB1a")
    {
      lw_porB_prof <- "porB1a"
    } else {

      if (lw_porB_G120 != "WT")
      {
        if (is.na(lw_porB_prof))
        {lw_porB_prof <- paste("porB ", lw_porB_G120, sep = "")} else
        {lw_porB_prof <- paste(lw_porB_prof, "/", lw_porB_G120, sep = "")}

        v_porB_G120 <- 1L
      } else {v_porB_G120 <- 0L}

      if (lw_porB_A121 != "WT")
      {
        v_porB_A121 <- 1L
        if (is.na(lw_porB_prof))
        {lw_porB_prof <- paste("porB ", lw_porB_A121, sep = "")} else
        {lw_porB_prof <- paste(lw_porB_prof, "/", lw_porB_A121, sep = "")}
      }else {v_porB_A121 <- 0L}
    }

    if (!is.na(lw_porB_prof)) {
      if (is.na(molec_profile)) {
        molec_profile <- lw_porB_prof
      } else {
        molec_profile <- paste(molec_profile, delim_dash, lw_porB_prof, sep = "")
      }
    }

    #--------------------------------------------------------------------------------------ponA
    lw_ponA <-   as.character(Combined_Output.df[m, "ponA"])

    if (lw_ponA == "?" | lw_ponA == 'x') {lw_ponA <- "Err"}

    if (lw_ponA != "WT") {
      if (is.na(molec_profile)) {
        molec_profile <- paste("ponA ", lw_ponA, sep = "")
      } else {
        molec_profile <- paste(molec_profile, delim_dash, "ponA ", lw_ponA, sep = "")
      }

      v_ponA_L421P <- 1L
    } else {
      v_ponA_L421P <- 0L
    }

    #--------------------------------------------------------------------------------------gyrA
    lw_gyrA <-   as.character(Combined_Output.df[m, "gyrA"])
    if (lw_gyrA == "?" | lw_gyrA == "x") {
      lw_gyrA <- "Err/Err"
    }

    lw_gyrA_parts <- unlist(strsplit(lw_gyrA, "/"))
    lw_gyrA_S91 <-  lw_gyrA_parts[1]
    lw_gyrA_D95 <- lw_gyrA_parts[2]
    lw_gyrA_prof <- NA

    if (lw_gyrA_S91 != "WT") {

      v_gyrA_S91 <- 1L

      if (is.na(lw_gyrA_prof)) {
        lw_gyrA_prof <- paste("gyrA ", lw_gyrA_S91, sep = "")
      } else {
        lw_gyrA_prof <- paste(lw_gyrA_prof, "/", lw_gyrA_S91, sep = "")
      }
    } else {v_gyrA_S91 <- 0L}

    if (lw_gyrA_D95 != "WT") {

      v_gyrA_D95 <- 1L

      if (is.na(lw_gyrA_prof)) {
        lw_gyrA_prof <- paste("gyrA ", lw_gyrA_D95, sep = "")
      } else {
        lw_gyrA_prof <- paste(lw_gyrA_prof, "/", lw_gyrA_D95, sep = "")
      }
    } else {
      v_gyrA_D95 <- 0L
    }

    if (!is.na(lw_gyrA_prof))
    {
      if (is.na(molec_profile))
      {molec_profile <- lw_gyrA_prof} else
      {molec_profile <- paste(molec_profile, delim_dash, lw_gyrA_prof, sep = "")}
    }


    #--------------------------------------------------------------------------------------parC
    # Variables ####
    lw_parC <-   as.character(Combined_Output.df[m, "parC"])
    if (lw_parC == "?" | lw_parC == "x") {
      lw_parC <- "Err/Err/Err"
    }

    lw_parC_parts <- unlist(strsplit(lw_parC, "/"))
    lw_parC_D86 <-  lw_parC_parts[1]
    lw_parC_S87 <- lw_parC_parts[2]
    lw_parC_S88 <- lw_parC_parts[3]
    lw_parC_prof <- NA

    if (lw_parC_D86 != "WT") # parC_D86 ####
    {
      v_parC_D86 <- 1L

      if (is.na(lw_parC_prof)) {
        lw_parC_prof <- paste("parC ", lw_parC_D86)
      } else {
        lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_D86, sep = "")
      }
    } else {
      v_parC_D86 <- 0L
    }


    if (lw_parC_S87 != "WT") # parC_S87 ####
    {
      v_parC_S87 <- 1L

      if (is.na(lw_parC_prof)) {
        lw_parC_prof <- paste("parC ", lw_parC_S87)
      } else {
        lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_S87, sep = "")
      }
    } else {
      v_parC_S87 <- 0L
    }

    if (lw_parC_S88 != "WT") # parC_S88 ####
    {
      V_parC_S88 <- 1L
      if (is.na(lw_parC_prof))
      {lw_parC_prof <- paste("parC ", lw_parC_S88)} else
      {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_S88, sep = "")}
    } else {v_parC_S88 <- 0L}

    if (!is.na(lw_parC_prof))
    {
      if (is.na(molec_profile))
      {molec_profile <- lw_parC_prof} else
      {molec_profile <- paste(molec_profile, delim_dash, lw_parC_prof, sep = "")}
    }


    # Variables ####
    #lw_rRNA23S <-   as.character(Combined_Output.df[m, "rRNA23S"]) #use allele counts instead
    #lw_rRNA23S_parts <- unlist(strsplit(lw_rRNA23S, "/"))
    lw_23S_A2059G <- as.integer(Combined_Output.df[m, "A2059G"])
    lw_23S_C2611T <- as.integer(Combined_Output.df[m, "C2611T"])
    lw_23S_prof <- NA

    if (is.na(lw_23S_A2059G)) {lw_23S_A2059G <- "Err"}
    if (is.na(lw_23S_C2611T)) {lw_23S_C2611T <- "Err"}

    if (lw_23S_A2059G != 0L)
    {
      if (is.na(lw_23S_prof))
      {lw_23S_prof <- paste("23S rRNA A2059G ", lw_23S_A2059G, "/4", sep = "")} else
      {lw_23S_prof <- paste(ls_23S_prof, "/A20159G ", lw_23S_A2059G, "/4",  sep = "")}
    }

    if (lw_23S_C2611T != 0L)
    {
      if (is.na(lw_23S_prof))
      {lw_23S_prof <- paste("23S rRNA C2611T ", lw_23S_C2611T, "/4", sep = "")} else
      {lw_23S_prof <- paste(lw_23S_prof, "/C2611T ", lw_23S_C2611T, "/4",  sep = "")}
    }

    if (!is.na(lw_23S_prof))
    {
      if (is.na(molec_profile))
      {molec_profile <- lw_23S_prof} else
      {molec_profile <- paste(molec_profile, delim_dash, lw_23S_prof, sep = "")}
    }

    #----------------------------------------------------------------------------------------------------
    #                             Assign Resistance to AZI, CEPH, CIP and TET, later assign MIC equation
    #----------------------------------------------------------------------------------------------------
    if (is.na(molec_profile))
    {
      molec_profile <- ""
    }

    # Variables ####
    amr_profile <- "Susceptible"
    sepr2 <- "/"

    #---------------------------------------------------------------------------------------------- AZI
    azi <- "Susceptible"

    if (str_detect(molec_profile, paste(c("mtrR Err", "A2059 Err", "C2611T Err"),collapse = '|')))
    {
      azi_interp <- "Error"
    } else {
      azi_MIC_inc <- round(1.27+
                            (2.48*lw_23S_A2059G)+
                            (1.10*lw_23S_C2611T)+
                            (0.82*v_mtrR35Adel)+
                            (2.86*v_mtrRMEN)+
                            (2.73*v_mtrRDisrupted)+
                            (2.41*v_ermB)+
                            (3.57*v_ermC)
                          )

      if (azi_MIC_inc > 12) {azi_MIC_inc <- 12L}

      # Hardcoded Dir ####
      azi_mic.df <- read.csv(here("data", "databases", "NGSTAR", "temp", "inc_mic_azi.csv"), # MAKE ME WORK WITH system.file
                             header = TRUE,
                             sep = ",",
                             stringsAsFactors = FALSE)
      azi <- paste(azi_mic.df$MIC[azi_mic.df$Inc == azi_MIC_inc], " ug/ml", sep = "")
      azi_interp <- azi_mic.df$Interp[azi_mic.df$Inc == azi_MIC_inc]

      if (azi_interp != "S" )
      {
        if (amr_profile == "Susceptible")
        {amr_profile <- paste("AZI-", azi_interp, sep = ""  )} else
        {amr_profile <- paste(amr_profile, sepr2, "AZI-", azi_interp, sep = "")}
      }
    }

    #---------------------------------------------------------------------------------------------- penicillin
    ### DF ###
    pen_MIC_inc <- round(5.98+
                          (0.69*v_mtrRMEN)+
                          (1.12*v_porB_G120)+
                          (1.71*v_ponA_L421P)+
                          (0.19*v_porB_A121)
                        )
    if (pen_MIC_inc > 16) {cfm_MIC_inc <- 16L}

    # Hardcoded Dir ####
    pen_mic.df <- read.csv(here("data", "databases", "NGSTAR", "temp", "inc_mic_penicillin.csv"),
                           header = TRUE,
                           sep = ",",
                           stringsAsFactors = FALSE)
    pen <- paste(pen_mic.df$MIC[pen_mic.df$Inc == pen_MIC_inc], " ug/ml", sep = "")
    pen_interp <- pen_mic.df$Interp[pen_mic.df$Inc == pen_MIC_inc]

    if (pen_interp != "S" )
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- paste("PEN-", pen_interp, sep = ""  )} else
      {amr_profile <- paste(amr_profile, sepr2, "PEN-", pen_interp, sep = "")}
    }

    #---------------------------------------------------------------------------------------------- cephalosporins
    if (str_detect(molec_profile, "penA Err"))
    {
      cfx_interp <- "Error"
      cfm_interp <- "Error"
    } else {

    cfx_MIC_inc <- round(0.98+
                          (0.69*v_mtrRMEN)+
                          (0.93*v_porB_G120)+
                          (0.24*v_ponA_L421P)+
                          (3.53*v_penA_A311V)+
                          (5.03*v_penA_A501P)+
                          (1.37*v_penA_A501T)+
                          (1.81*v_penA_A501V)+
                          (2.83*v_penA_N513Y)+
                          (0.65*v_penA_A517G)+
                          (0.74*v_penA_G543S)+
                          (0.36*v_mtrR_G45D)
                        )

    if (cfx_MIC_inc > 10) {cfm_MIC_inc <- 10L}

    # Hardcoded Dir ####
    cfx_mic.df <- read.csv(here("data", "databases", "NGSTAR", "temp", "inc_mic_ceftriaxone.csv"),
                           header = TRUE,
                           sep = ",",
                           stringsAsFactors = FALSE)
    cfx <- paste(cfx_mic.df$MIC[cfx_mic.df$Inc == cfx_MIC_inc], " ug/ml", sep = "")
    cfx_interp <- cfx_mic.df$Interp[cfx_mic.df$Inc == cfx_MIC_inc]

    cfm_MIC_inc <- round(0.98+
                          (0.68*v_mtrRMEN)+
                          (0.93*v_porB_G120)+
                          (0.24*v_ponA_L421P)+
                          (4.02*v_penA_A311V)+
                          (5.03*v_penA_A501P)+
                          (1.37*v_penA_A501T)+
                          (1.81*v_penA_A501V)+
                          (2.83*v_penA_N513Y)+
                          (0.65*v_penA_A517G)+
                          (0.74*v_penA_G543S)+
                          (0.36*v_mtrR_G45D)
                          )
    if (cfm_MIC_inc > 10) {cfm_MIC_inc <- 10L}

    # Hardcoded Dir ####
    cfm_mic.df <- read.csv(here("data", "databases", "NGSTAR", "temp", "inc_mic_cefixime.csv"),
                           header = TRUE,
                           sep = ",",
                           stringsAsFactors = FALSE)
    cfm <- paste(cfm_mic.df$MIC[cfm_mic.df$Inc == cfm_MIC_inc], " ug/ml", sep = "")
    cfm_interp <- cfm_mic.df$Interp[cfm_mic.df$Inc == cfm_MIC_inc]
    }

    if (cfm_interp != "S" )
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- paste("CE-", cfm_interp, sep = ""  )} else
      {amr_profile <- paste(amr_profile, sepr2, "CE-", cfm_interp, sep = "")}
    }

    if (cfx_interp != "S" )
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- paste("CX-", cfx_interp, sep = ""  )} else
      {amr_profile <- paste(amr_profile, sepr2, "CX-", cfx_interp, sep = "")}
    }

    #---------------------------------------------------------------------------------------------- fluorquilones
    cip_MIC_inc <- round(1.08+
                          (4.03*v_gyrA_S91)+
                          (4.68*v_gyrA_D95)+
                          (0.33*v_parC_D86)+
                          (0.81*v_parC_S87)+
                          (1.45*v_parC_S88)
                        )
    if (cip_MIC_inc > 13) {cfm_MIC_inc <- 13L}

    # Hardcoded Dir ####
    cip_mic.df <- read.csv(here("data", "databases", "NGSTAR", "temp", "inc_mic_ciprofloxacin.csv"),
                           header = TRUE,
                           sep = ",",
                           stringsAsFactors = FALSE)
    cip <- paste(cip_mic.df$MIC[cip_mic.df$Inc == cip_MIC_inc], " ug/ml", sep = "")
    cip_interp <- cip_mic.df$Interp[cip_mic.df$Inc == cip_MIC_inc]

    if (cip_interp != "S" )
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- paste("CIP-", cip_interp, sep = ""  )} else
      {amr_profile <- paste(amr_profile, sepr2, "CIP-", cip_interp, sep = "")}
    }

    #---------------------------------------------------------------------------------------------- tetracycline
    tet <- "<= 1 ug/ml"
    if (str_detect(molec_profile, "rpsJ Err"))
      {tet <- "Error"
      amr_profile <- "Error"
      }else
    {

      if (lw_rpsJ == "V57M")
      {
        tet <- "4 ug/ml"
      }
      if (lw_tetM == "POS")
      {
        tet <- ">= 16 ug/ml"
      }

      if (tet != "<= 1 ug/ml")
      {
        if (amr_profile == "Susceptible")
        {amr_profile <- "TET-R"} else
        {amr_profile <- paste(amr_profile, sepr2, "TET-R", sep = "")}
      }

    }

    #-----------------------------------------------------------------------------------  spectinomycin
    # Variables ####
    spe <- "<= 8 ug/ml"

    if (str_detect(molec_profile, "16S rRNA Err"))
    {
      spe <- "Error"
      amr_profile <- "Error"
    } else {

      if (lw_16S_C1192T == "C1192T")
      {
        spe <- ">=128 ug/ml"
      }
      if (lw_16S_T1458C == "T1485C")
      {
        spe <- "16 ug/ml"
      }

      if (spe == ">=128 ug/ml" | spe == "64 ug/ml")
      {
        if (amr_profile == "Susceptible")
        {amr_profile <- "SPE-R"} else
        {amr_profile <- paste(amr_profile, sepr2, "SPE-R", sep = "")}
      }

    }

    #--------------------------------------------------------------------------------------
    # Used to create the return var
    sample_data.df <- data.frame(lw_CurrSampleNo, lw_ermB, lw_ermC, lw_rpsJ, lw_tetM, lw_bla, lw_penA_prof,
                                lw_mtrR_p, lw_mtrR_A39, lw_mtrR_G45,
                                lw_porB_struct, lw_porB_G120, lw_porB_A121, lw_ponA,
                                lw_gyrA_S91, lw_gyrA_D95, lw_parC_D86, lw_parC_S87, lw_parC_S88,
                                lw_23S_A2059G, lw_23S_C2611T, lw_16S_C1192T, lw_16S_T1458C, amr_profile,
                                azi, cfx, cfm, cip, tet, pen, spe, stringsAsFactors = FALSE)

    # Initialize Return Var ####
    if(m == 1) {
      lw_Output.df <- data.frame(sample_data.df, stringsAsFactors = FALSE)
    } else {
      lw_Output.df <- bind_rows(lw_Output.df, sample_data.df)
    }

  } # End For Loop

  # Output ####
  write.csv(x = lw_Output.df,
            #here("data", "input", "labware", "LabWareUpload_GONO_AMR.csv"),
            paste(out_location, paste(Sys.Date(), org_id, "LabwareGONOAMR", "WADE.csv", sep = "_"), sep = ""),
            quote = FALSE,
            row.names = FALSE)

  # Return ####
  return(lw_Output.df)
} # End Function

# I want it to just change the objects to change sub
change_enz <- function(enz, enz_value, enz_profile, enz_string){
  if(enz != "WT") {
    enz_value <- 1

    if(is.na(enz_profile)){
      enz_profile <- paste(enz_string, " ", enz)
    } else {
      enz_profile <- paste(enz_prof, "/", enz, sep = "")
    }
  } else {
    enz_value <- 0
  }
}
