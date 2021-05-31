##### Create and Save Projections for 8 weeks after 1st December 2020

### Load in libraries and data
setwd("Joint_Shiny_App/")
source("global.R")
library(xtable)

### Solve SEIR ODEs
# Set dates  
start_date <- as.Date("2021-02-01")
date_start_seq <- c(start_date, start_date + 30)
date_start_end <- c(start_date + 29, start_date + 55)
tmax <- as.numeric( difftime(start_date + 55, as.Date('2020-02-29'), units = "days") )
start_date - date_start_end[2] - 1

# Set lockdown levels
all_levs <- vector("list", 7)
str <- 'Avg. Lockdown Level '
all_levs[[1]] <- paste0(str, rep(0, 2))
all_levs[[2]] <- paste0(str, rep(3, 2))
all_levs[[3]] <- paste0(str, rep(5, 2))
all_levs[[4]] <- paste0(str, c(3, 5))
all_levs[[5]] <- paste0(str, c(2, 5))
all_levs[[6]] <- paste0(str, c(1, 5))
all_levs[[7]] <- paste0(str, c(0, 5))
len <- length(all_levs)

# Find estimated cost
diff_date <- as.numeric( difftime(date_start_end, date_start_seq) )
est_costs$Level[1:6] <- paste0(str, 0:5) 
est_costs$Estimated_Costs[which(est_costs$Level %in% all_levs[[2]])]
cost_extract <- function(x) {
  y <- which(est_costs$Level %in% x)
  if (length(y) < 2) { y <- rep(y, 2) }
  est_costs$Estimated_Costs[y]
}

all_costs <- lapply(all_levs, cost_extract) 
costs <- sapply(all_costs, "%*%", diff_date)
nms <- c( 'No intervention', 'Level 3', 'Level 5',
          'Level 3 > Level 5', 'Level 2 > Level 5',
          'Level 1 > Level 5', 'No Intervention > Level 5')
costs_df <- data.frame(Policy = nms, Cost = costs)

pdf('projected_costs.pdf', height = 10)
par(mar = c(8, 4, 4, 2) + 0.1)
barplot(height = costs_df$Cost[-7], las = 2, ylim = c(-200, 400),
        col = c('orange', rep('dodgerblue', 5)), 
        names.arg = costs_df$Policy[-7], 
        main = 'Projected Costs (in billions of euros)')
abline(h = seq(-200, 400, 100), lty = 3, col = 'grey')
dev.off()

# Create lockdown forecast data
lockdown_forecast <- optim_dat <- all_linfo <- vector("list", 7)
for (i in 1:len) {
  
  lockdown_forecast[[i]] <- rbind(interventions_info,
                                  data.frame(start = date_start_seq, 
                                             end = date_start_end,
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
UL <- MID <- LL <- vector("list", len)
for (j in 1:len) {
  linfo <- all_linfo[[j]]
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
  UL[[j]] <- foreach(x = All_Runs, .combine = rbind,
                     .final = as.data.frame) %dopar% {
    sapply(x, quantile, probs = 0.975, names = FALSE)
  }
  LL[[j]] <- foreach(x = All_Runs, .combine = rbind,
                .final = as.data.frame) %dopar% {
    sapply(x, quantile, probs = 0.025, names = FALSE)
  }
  rm(All_Runs)
  
  MID[[j]] <- SEIR_model_simulation(pars = def_pars,
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

### Repeat the above but reducing the contacts of the 70+ population to 0
set_to_zero <- function(x, y = 15:16) {
  x[y, ] <- x[, y] <- 0
  x
}
contacts_isolated70 <- lapply(contacts, set_to_zero)

dub_def_beta_isolated70 <- getbeta(3.7, pars = def_pars, p_age = dub_population$propage,
            CONTACTMATRIX = contacts_isolated70)

cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(
           cl, varlist=c("SEIR_model_simulation", "def_pars", 
                         "contacts_isolated70", "dub_xstart", "dub_population",
                         "dub_def_beta_isolated70", "lsoda", "SEIR_model_D", 
                         "tmax", "all_linfo"),
           envir=environment())

# Run SEIR Models
UL <- MID <- LL <- vector("list", len)
for (j in 1:len) {
  linfo <- all_linfo[[j]]
  All_Runs <- foreach(i = 5:ncol(linfo), .combine = rbind) %dopar% {
    SEIR_model_simulation(pars = def_pars,
                          contacts_ireland = contacts,
                          dateStart = as.Date('2020-02-29'),
                          startval = dub_xstart,
                          POP = dub_population,
                          beta = dub_def_beta,
                          tmax = tmax, 
                          lockdown_information = linfo[, c(2:3, i)],
                          isolated = TRUE,
                          isolated_contacts = contacts_isolated70,
                          isolated_beta = dub_def_beta_isolated70)$solution
  }
  
  All_Runs <- split(All_Runs, All_Runs$time)
  UL[[j]] <- foreach(x = All_Runs, .combine = rbind,
                     .final = as.data.frame) %dopar% {
    sapply(x, quantile, probs = 0.975, names = FALSE)
  }
  LL[[j]] <- foreach(x = All_Runs, .combine = rbind,
                .final = as.data.frame) %dopar% {
    sapply(x, quantile, probs = 0.025, names = FALSE)
  }
  rm(All_Runs)
  
  MID[[j]] <- SEIR_model_simulation(pars = def_pars,
                                    contacts_ireland = contacts,
                                    dateStart = as.Date('2020-02-29'),
                                    startval = dub_xstart,
                                    POP = dub_population,
                                    beta = dub_def_beta,
                                    tmax = tmax, 
                                    lockdown_information = linfo[, 2:4],
                                    isolated = TRUE,
                                    isolated_contacts = contacts_isolated70,
                                    isolated_beta = dub_def_beta_isolated70)$solution
    
}

stopCluster(cl)
forecast_fits_isolated70 <- list(MID = MID, UL = UL, LL = LL)
rm(MID, UL, LL)

# Write results to file
save(forecast_fits_isolated70, file = "data/forecast_fits_isolated70.Rdata")


