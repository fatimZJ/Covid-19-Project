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
  
####-- 4. Extra App Data ----------------------------------------------------####
  
  # Create Default Age Groups
  def_age_groups <- c(paste0(seq(0, 70, 5), " - ", seq(4, 74, 5)), "75+")
  
  