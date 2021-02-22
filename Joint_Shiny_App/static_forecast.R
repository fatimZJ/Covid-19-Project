##### Create and Save Projections for 8 weeks after 1st December 2020

### Load in libraries and data
setwd("Joint_Shiny_App/")
source("global.R")

### Solve SEIR ODEs
# Set dates  
start_date <- as.Date("2020-12-01")
date_start_seq <- seq.Date(start_date, start_date + 55L, 14L)
tmax <- as.numeric( difftime(start_date + 55, as.Date('2020-02-29'), units = "days") )

# Set lockdown levels
all_levs <- vector("list", 5)
all_levs[[1]] <- rep(optim_res$policy[1], 4)
all_levs[[2]] <- rep(optim_res$policy[8], 4)
all_levs[[3]] <- rep(optim_res$policy[9], 4)
all_levs[[4]] <- c( rep(optim_res$policy[8], 2), rep(optim_res$policy[9], 2) )
all_levs[[5]] <- c( rep(optim_res$policy[7], 2), rep(optim_res$policy[9], 2) )

# Create lockdown forecast data
lockdown_forecast <- optim_dat <- all_linfo <- vector("list", 5)
for (i in 1:5) {
  
  lockdown_forecast[[i]] <- rbind(interventions_info,
                                  data.frame(start = date_start_seq, 
                                             end = date_start_seq + 13L,
                                             policy = all_levs[[i]]))
    optim_dat[[i]] <- merge(lockdown_forecast[[i]], optim_res, by = "policy")
    all_linfo[[i]] <- merge(optim_dat[[i]], boot_scales_t, by = "policy")
  
}

rm(lockdown_forecast, optim_dat)

# Parallel initiation
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(
           cl, varlist=c("SEIR_model_simulation", "def_pars", 
                         "contacts", "dub_xstart", "dub_population",
                         "dub_def_beta", "lsoda", "SEIR_model_D", 
                         "tmax", "all_linfo"),
           envir=environment())

# Run SEIR Models
UL <- MID <- LL <- vector("list", 5)
for (i in 1:5) {
  linfo <- all_linfo[[i]]
  All_Runs <- foreach(i = 5:ncol(linfo), .combine = rbind) %dopar% {
    SEIR_model_simulation(pars = def_pars,
                          contacts_ireland = contacts,
                          dateStart = as.Date('2020-02-29'),
                          startval = dub_xstart,
                          POP = dub_population,
                          beta = dub_def_beta,
                          tmax = tmax, 
                          lockdown_information = linfo[, c(2:3, i)])$solution
  }
  
  All_Runs <- split(All_Runs, All_Runs$time)
  UL[[i]] <- foreach(x = All_Runs, .combine = rbind,
                     .final = as.data.frame) %dopar% {
    sapply(x, quantile, probs = 0.975, names = FALSE)
  }
  LL[[i]] <- foreach(x = All_Runs, .combine = rbind,
                .final = as.data.frame) %dopar% {
    sapply(x, quantile, probs = 0.025, names = FALSE)
  }
  rm(All_Runs)
  
  MID[[i]] <- SEIR_model_simulation(pars = def_pars,
                                    contacts_ireland = contacts,
                                    dateStart = as.Date('2020-02-29'),
                                    startval = dub_xstart,
                                    POP = dub_population,
                                    beta = dub_def_beta,
                                    tmax = tmax, 
                                    lockdown_information = linfo[, 2:4])$solution
    
}

stopCluster(cl)
forecast_fits <- list(MID = MID, UL = UL, LL = LL)
rm(MID, UL, LL)

# Write results to file
save(forecast_fits, file = "data/forecast_fits.Rdata")

### Repeat the above but reducing the contacts of the 75+ population to 0
set_to_zero <- function(x, y = 16) {
  x[y, ] <- x[, y] <- 0
  x
}
contacts_isolated75 <- lapply(contacts, set_to_zero)

dub_def_beta_isolated75 <- getbeta(3.7, pars = def_pars, p_age = dub_population$propage,
            CONTACTMATRIX = contacts_isolated75)

cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(
           cl, varlist=c("SEIR_model_simulation", "def_pars", 
                         "contacts_isolated75", "dub_xstart", "dub_population",
                         "dub_def_beta_isolated75", "lsoda", "SEIR_model_D", 
                         "tmax", "all_linfo"),
           envir=environment())

# Run SEIR Models
UL <- MID <- LL <- vector("list", 5)
for (i in 1:5) {
  linfo <- all_linfo[[i]]
  All_Runs <- foreach(i = 5:ncol(linfo), .combine = rbind) %dopar% {
    SEIR_model_simulation(pars = def_pars,
                          contacts_ireland = contacts_isolated75,
                          dateStart = as.Date('2020-02-29'),
                          startval = dub_xstart,
                          POP = dub_population,
                          beta = dub_def_beta,
                          tmax = tmax, 
                          lockdown_information = linfo[, c(2:3, i)],
                          isolated = TRUE,
                          isolated_contacts = contacts_isolated75,
                          isolated_beta = dub_def_beta_isolated75)$solution
  }
  
  All_Runs <- split(All_Runs, All_Runs$time)
  UL[[i]] <- foreach(x = All_Runs, .combine = rbind,
                     .final = as.data.frame) %dopar% {
    sapply(x, quantile, probs = 0.975, names = FALSE)
  }
  LL[[i]] <- foreach(x = All_Runs, .combine = rbind,
                .final = as.data.frame) %dopar% {
    sapply(x, quantile, probs = 0.025, names = FALSE)
  }
  rm(All_Runs)
  
  MID[[i]] <- SEIR_model_simulation(pars = def_pars,
                                    contacts_ireland = contacts,
                                    dateStart = as.Date('2020-02-29'),
                                    startval = dub_xstart,
                                    POP = dub_population,
                                    beta = dub_def_beta,
                                    tmax = tmax, 
                                    lockdown_information = linfo[, 2:4],
                                    isolated = TRUE,
                                    isolated_contacts = contacts_isolated75,
                                    isolated_beta = dub_def_beta_isolated75)$solution
    
}

stopCluster(cl)
forecast_fits_isolated75 <- list(MID = MID, UL = UL, LL = LL)
rm(MID, UL, LL)

# Write results to file
save(forecast_fits_isolated75, file = "data/forecast_fits_isolated75.Rdata")


