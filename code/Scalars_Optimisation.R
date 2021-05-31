## Load packages
library(tidyverse)
library(deSolve)
library(Matrix)
library(optimx)
library(parallel)

## Load the data
source('code/function_load_data.r')
#Loading the required data
data_list <- load_data()
list2env(data_list, envir = .GlobalEnv)
rm(data_list)

## Load the get beta function
source("code/function_getbeta.R")

## Load the model function
source("code/function_SEIR_Model.R")

## load the simulation function 
source("code/function_SEIR_Simulation.R")

## Define objective function
#scalars <- scalars_init
optim_fun <- function(scalars, cumulative_obs_cases, population_df, contacts_list, interventions_df) {
  
  scalars_est <- exp(scalars)  
  
  names(scalars_est) <- unique(interventions_df$policy)
  
  Sim <- try(simulation_SEIR_model(R0t = 3.4,
                                   POP = population_df,
                                   contacts_ireland = contacts_list,
                                   interventions = interventions_df, 
                                   dt = 1,
                                   dateStart = interventions_df$start[1],
                                   dateEnd = cumulative_obs_cases$date[dim(cumulative_obs_cases)[1]],
                                   scalars = scalars_est), silent = TRUE)
  
  if (class(Sim) != "try-error") {
    
    Est_Cc <- rowSums(Sim$sol_out[grepl('Cc_',names(Sim$sol_out))])
    
    ss <- try(sum((cumulative_obs_cases$cases - (Est_Cc))^2), silent = TRUE)
    
    if ( !is.na(ss) & class(ss) != "try-error"){
      RSS <- ss
    } else {
      RSS <- Inf}
  } else {
    RSS <- Inf
  }
  
  RSS
  
}


Rand_start <- function(i) {
  
  inters <- length(unique(interventions_info$policy))
  
  scalars_init <- c(runif(1, 1.5, 3), runif((inters - 1), 0, 1))
  
  startt <- Sys.time()
  ests_optimx <-  optimx(log(scalars_init), 
                         optim_fun,
                         cumulative_obs_cases = cumulative_cases, 
                         itnmax = 5000,
                         interventions_df = interventions_info, 
                         population_df = population,
                         contacts_list = contacts, 
                         method = c("Nelder-Mead"),
                         control = list(maximize = FALSE, kkt = FALSE))
  endt <- Sys.time()
  
  ests <- c(as.numeric(exp(ests_optimx[1:9])), ests_optimx$value, 
            ests_optimx$convcode, ests_optimx$fevals,
            as.numeric(endt - startt))
  names(ests) <- c(paste0("Scalar", 1:inters), "TSS", "Convergence", "NO. Calls", "time")
  return(ests)
  
}

nRandStarts <- 2



ncores <- detectCores()
cl <- makeCluster(ncores - 4) 

clusterEvalQ(cl, {
  library(deSolve)
  library(optimx)
  library(tidyverse)
  library(Matrix)
})

clusterExport(cl, c("optim_fun", "simulation_SEIR_model", "SEIR_model", "getbeta", 
                    "contacts", "cumulative_cases", "interventions_info", "population"))

startt <- Sys.time()

estsList <- parLapply(cl, 1:nRandStarts, Rand_start) 

estsMat <- as.data.frame(do.call(rbind, parestsMat))
names(estsMat) <- names(estsList[[1]])

Endt <- Sys.time()

stopCluster(cl)

write.csv(file =  "MultipleInitialResults_PARALLEL_test.csv", estsMat)

Endt - startt

estsMat[which.min((estsMat$TSS)),]

#optim_fun(log(as.numeric(mutlistarts[which.min((mutlistarts$TSS)),1:inters])),
#          interventions_df = interventions_info, 
#          population_df = population, contacts_list = contacts, 
#          cumulative_obs_cases = cumulative_cases)


