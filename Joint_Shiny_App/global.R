####------------------ Covid Predictor: Global Shiny Fine ------------------####

####-- 1. Required Libraries -----------------------------------------------####
  
  # List of all required packages for the application to run
  pkgs <- c('shiny', 'mlbench', 'shinythemes', 'dplyr', 'shinyWidgets',
            'ggplot2', 'shinydashboard', 'data.table', 'optimx', 'Matrix',
            'shinycssloaders', 'DT', 'deSolve', 'plotly', 'readr', 'tidyverse')
  
  # Handles any packages that are not installed, and installs them locally
  pkgs_to_install <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if (length(pkgs_to_install)) {
    install.packages(pkgs_to_install)
  } 
  
  # Load all packages - continues if all loaded safely, errors if any failures
  pkg_status <- sapply(pkgs, function(x) {
    suppressPackageStartupMessages(
      invisible(require(x, character.only = TRUE))
    )}
  )
  
  if (sum(pkg_status) != length(pkgs)) {
    stop("Error loading some of the packages required - please check!")
  }

####-- 2. Required Functions -----------------------------------------------####
  
  # Locate all files in the code directory
  func_files <- list.files(path = 'code', full.names = TRUE)
  
  # Load all function files in the code directory
  func_status <- sapply(func_files, source)
  
####-- 3. Required Data ----------------------------------------------------####
  
  # Load Data Objects as a List
  data_sources <- load_data(data_path = 'data/', county = 'Dublin')
  
  # Create Data Frames from the Data Object List
  list2env(data_sources, envir = .GlobalEnv)
  
  # Remove the List containing the Data Objects
  rm(data_sources)
  
####-- 4. Extra App Data  --------------------------------------------------####
  
  # Create Default Age Groups
  def_age_groups <- c(paste0(seq(0, 70, 5), " - ", seq(4, 74, 5)), "75+")
  
  # Change the population data age labels
  dub_population$agegroup <- def_age_groups
  irl_population$agegroup <- def_age_groups
  
  # Vector of compartment names
  comp_vec <- c("S", "Ev", "Ip", "IA", "Ii", "It", "Iti", "Iq", "R")
  
  # Compartment extract function
  comp_extract <- function(dat, comp) {
    dat[grepl(paste0(comp, "_"), names(dat))]
  }
  
  # Summary plotting function
  summary_plot <- function(dat, xval, comp_vec, group, I_type){
    
    # Extract Compartments
    comps <- Map(comp_extract, comp = comp_vec,
                 MoreArgs = list(dat = dat))
    names(comps)[3:8] <- c("I_Pr", "I_As", "I_Im", "I_Aw", "I_Is", "I_No")
    comps$I_Sy <- comps$I_Im + comps$I_Aw + comps$I_Is + comps$I_No
    comps$I_Al <- comps$I_Pr + comps$I_As + comps$I_Sy 
    comps$I <- comps[[paste0("I_", substr(I_type, start = 1, stop = 2))]]
    
    # Sum across desired age groups
    comp_groups <- lapply(comps[c("S", "Ev", "I", "R")],
                          "[", i = group)
    comp_draw <- lapply(comp_groups, rowSums)
    
    # Draw the Plot
    plot_ly(x = ~xval, y = ~comp_draw$S, name = 'Susceptible', type = 'scatter', mode = 'lines',
            line = list(color = 'rgb(69, 95, 245)')) %>%
      add_trace(y = ~comp_draw$Ev, name = 'Exposed', mode = 'lines', line = list(color = 'rgb(214, 122, 17)')) %>% 
      add_trace(y = ~comp_draw$I, name = 'Infected', mode = 'lines', line = list(color = 'rgb(186, 24, 19)')) %>%
      add_trace(y = ~comp_draw$R, name = 'Removed', mode = 'lines', line = list(color = 'rgb(23, 191, 26)')) %>% 
      layout(xaxis = list(title = "Time"), 
             yaxis = list(title = "Compartment Size"), 
             title = list(text = "Overall Output"))
    
  }
  
  comp_plot <- function(dat, xval, comp, group, inp){
    
    comp_dat <- dat[grepl(comp, names(dat))][group]
    
    ### Draw the Plot
    p <- plot_ly(x = ~xval, y = NA, type = 'scatter', mode = 'lines')
    
    for ( i in seq_along(inp$age_sel) ) {
      p <- p %>% add_trace(y = comp_dat[[i]], name = inp$age_sel[i], mode = 'lines')
    }
    
    p 
    
  }
  
  td <- as.Date(Sys.time())
  
  # Lockdown Info
  l_mean <- apply(boot_lockdown_scalars, 2, mean)
  l_q_025 <- apply(boot_lockdown_scalars, 2, quantile, 0.025)
  l_q_975 <- apply(boot_lockdown_scalars, 2, quantile, 0.975)
  
  interventions_info$l_mean <- l_mean
  interventions_info$l_q_025 <- l_q_025
  interventions_info$l_q_975 <- l_q_975
  
  linfo_mean <- interventions_info[c(1, 2, 4)]
  linfo_lower <- interventions_info[c(1, 2, 5)]
  linfo_upper <- interventions_info[c(1, 2, 6)]
  
  