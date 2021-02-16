####------------------ Covid Predictor: Global Shiny Fine ------------------####

####-- 1. Required Libraries -----------------------------------------------####
  
  # List of all required packages for the application to run
  pkgs <- c('shiny', 'mlbench', 'shinythemes', 'dplyr', 'shinyWidgets',
            'ggplot2', 'shinydashboard', 'data.table', 'optimx', 'Matrix',
            'shinycssloaders', 'DT', 'deSolve', 'plotly', 'readr', 'tidyverse',
            'parallel', 'doParallel')
  
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
  #func_status <- sapply(func_files, source)
  sapply(func_files, source)
  
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
  
  # Vector of lockdown measures
  lockdown_measures <- c('No Intervention', paste0('Lockdown Level ', 1:5))
  
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
  
  # Starting Values
  dub_xstart <- c(dub_population$popage - 15/16,
                  rep(1, 16),
                  rep(1/16, 16),
                  rep(0, 16 * 6))
  irl_xstart <- c(irl_population$popage - 15/16,
                  rep(1, 16),
                  rep(1/16, 16),
                  rep(0, 16 * 6))
  names(dub_xstart) <- names(irl_xstart) <- c(paste0('S_',1:16),
                                              paste0('Ev_',1:16),
                                              paste0('Ip_',1:16),
                                              paste0('IA_',1:16),
                                              paste0('Ii_',1:16),
                                              paste0('It_',1:16),
                                              paste0('Iti_',1:16),
                                              paste0('Iq_',1:16),
                                              paste0('R_',1:16))
  
  def_pars <- c(3.7, 5.8, 11.6, 0.55, 0.2, 0.8, 0.1, 0.05, 7)
  names(def_pars) <- c('L', 'Cv', 'Dv', 'h', 'f', 'tv', 'q', 'k', 'TT')
  
  dub_def_beta <- getbeta(3.7, pars = def_pars, p_age = dub_population$propage,
            CONTACTMATRIX = contacts)
  
  N_known <- 13
  N_full <- 68
  
  # Create Forecast Components
  comp_sel <- function(x, y) {
    comps <- Map(comp_extract, comp = comp_vec, 
                 MoreArgs = list(dat = x))
    names(comps) <- c("Su", "Ex", "Pr", "As", "I_Im", 
                      "I_Aw", "I_Is", "I_No", "Re")
    comps$Sy <- comps$I_Im + comps$I_Aw + comps$I_Is + comps$I_No
    comps$Al <- comps$Pr + comps$As + comps$Sy
    comps_sel <- comps[[substr(y, start = 1, stop = 2)]]
    N <- nrow(comps_sel)
    rowSums(comps_sel[(N-N_full+1):N, ])
  }
  
  intervention_adjust <- function(dat, x) {
    
    old_end <- dat$end[nrow(dat)]
    missing_interventions <- as.numeric(difftime(x, old_end, units = "days"))
    if (missing_interventions <= 1)
      return(dat)
    rbind(dat, data.frame(start = old_end + 1, 
                          end = x, 
                          policy = "No Intervention"))
    
  }
  
  ### Alter datasets to include lockdown measures
  L1 <- mean(c(optim_res$optim_res[optim_res$policy == 'No Intervention'],
               optim_res$optim_res[optim_res$policy == 'Lockdown Level 2']))
  L4 <- mean(c(optim_res$optim_res[optim_res$policy == 'Lockdown Level 3'],
               optim_res$optim_res[optim_res$policy == 'Lockdown Level 5']))
  optim_res <- rbind(optim_res, data.frame(policy = paste0('Lockdown Level ', c(1, 4)),
                                           optim_res = c(L1, L4)))
  
  boot_lockdown_scalars['Lockdown Level 1'] <- (boot_lockdown_scalars$`No Intervention` + boot_lockdown_scalars$`Lockdown Level 2`)/2
  boot_lockdown_scalars['Lockdown Level 4'] <- (boot_lockdown_scalars$`Lockdown Level 3` + boot_lockdown_scalars$`Lockdown Level 5`)/2
  
  boot_scales_t <- cbind( policy = names(boot_lockdown_scalars), 
                          as.data.frame(t(boot_lockdown_scalars)) )
  rm(boot_lockdown_scalars)
    
  ### Get quantiles
  df_quant <- function(x, p) {
    sapply(x, quantile, probs = p)
  }
  
  ### Parallelisation
  n_cores <- detectCores()-1
  
  ### For debugging
  #boot_scales_t <- boot_scales_t[1:11]
  #boot_lockdown_scalars <- boot_lockdown_scalars[1:10, ]
  
  