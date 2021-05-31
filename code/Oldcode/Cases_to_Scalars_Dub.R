## Load packages
library(deSolve)
#library(DEoptim)
library(optimx)
library(numDeriv)
library(parallel)
## Load the data
#source('code/1_loadData.r')

source('code/01-load_data.r')
data_list <- load_data()

#data <- load_data(data_path = 'data/', county = 'Dublin')

# Create Data Frames from the Data Object List
list2env(data_list, envir = .GlobalEnv)

# Remove the List containing the Data Objects
rm(data_list)

## Load the get beta function
source("code/getbeta.R")

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
      tss <- ss
    } else {
      tss <- Inf}
  } else {
    tss <- Inf
  }
  
  tss
  
}

## Nlm

scalars_init <- c(2.175648, 0.7734942, 0.2349296, 0.05953305, 0.1763019, 0.3924366, 0.4126519,
                 0.4111283, 0.1557194)

startt <- Sys.time()
ests_nlm <- nlm(optim_fun,log(scalars_init),
                interventions_df = interventions_info, population_df = population,
                contacts_list = contacts, #stepmax = 0.5,
                iterlim = 1000, cumulative_obs_cases = cumulative_cases)
endt <- Sys.time()

exp(ests_nlm$estimate)
ests_nlm$min
ests_nlm$code

nlm_est <- c(2.27619341, 0.94894146, 0.21845530, 0.01503771, 0.31684429, 0.40807037, 0.49130460,
             0.34598300, 0.20890798)
#31219054

nlm_est <- c(2.4242360, 0.9130557, 0.2019896, 0.0153258, 0.3168725, 0.4115858,
             0.5386053, 0.3075835, 0.2579594)
#53007074
nlm_est <- c(2.34845285, 0.90252470, 0.22250976, 0.01596992, 0.30632324, 0.37813130, 0.58107912,
             0.30148189, 0.24569175)
#38799530

optim_fun(log(nlm_est), interventions_df = interventions_info, 
        population_df = population, contacts_list = contacts, 
        cumulative_obs_cases = cumulative_cases)

## Optimx- NM

startt <- Sys.time()
ests_optimx <-  optimx(log(scalars_init), 
                       optim_fun,
                       #lower = 0, upper = 3,
                       cumulative_obs_cases = cumulative_cases, itnmax = 5000,
                       interventions_df = interventions_info, population_df = population,
                       contacts_list = contacts, 
                       method = c("Nelder-Mead"),
                       control = list(maximize = FALSE, kkt = FALSE))
endt <- Sys.time()


endt - startt
c(as.numeric(exp(ests_optimx[1:9])), ests_optimx$value, ests_optimx$convcode, ests_optimx$fevals,
  endt - startt)

c(2.033784e+00, 1.093872e+00, 2.143144e-01, 3.064505e-02, 3.295509e-01, 3.654255e-01,
  5.485281e-01, 3.371615e-01, 2.035889e-01)#, 3.064904e+07, 1.000000e+01, 2.362000e+03, 1.95 hours)
## DEoptim

lower <- rep(0, length(scalars_init))
upper <- rep(4, length(scalars_init))

startt <- Sys.time()
ests_DEoptim <- DEoptim(optim_fun, cumulative_obs_cases = cumulative_cases,
                        lower, upper, DEoptim.control(itermax = 200, reltol = 1e-8,
                                                      steptol = 50))

endt <- Sys.time()

################################################################################

## Multiple random starts 

## Nlm
 
n <- 20 

inters <- length(unique(interventions_info$policy))
ests <- matrix( NA, nrow = n, ncol = inters + 4 )
for (i in 1:n) {
  
  scalars_init <- c(runif(1, 1.5, 3), runif((inters - 1), 0, 1))
  
  startt <- Sys.time()
  
  ests_nlm <- nlm(optim_fun,log(scalars_init),
                  #stepmax = 0.5,
                  interventions_df = interventions_info, population_df = population,
                  contacts_list = contacts, 
                  iterlim = 1000, cumulative_obs_cases = cumulative_cases)
  endt <- Sys.time()
  ests[i,] <- c(exp(ests_nlm$estimate), ests_nlm$minimum, ests_nlm$code, ests_nlm$iterations,
                as.numeric(endt - startt))
  
  write.csv(file =  "MultipleInitialResults_nlm.csv", ests)
  
}


## Optimx- NM
#set.seed(97532468)
n <- 20 

inters <- length(unique(interventions_info$policy))
ests <- matrix( NA, nrow = n, ncol = inters + 4 )

colnames(ests) <- c(paste0("Scalar", 1:inters), "TSS", "Convergence", "NO. Calls", "time")

for (i in 1:n) {
  
  scalars_init <- c(runif(1, 1.5, 3), runif((inters - 1), 0, 1))
  
  startt <- Sys.time()
  ests_optimx <-  optimx(log(scalars_init), 
                         optim_fun,
                         cumulative_obs_cases = cumulative_cases, itnmax = 5000,
                         interventions_df = interventions_info, population_df = population,
                         contacts_list = contacts, 
                         method = c("Nelder-Mead"),
                         control = list(maximize = FALSE, kkt = FALSE))
  endt <- Sys.time()
  
  ests[i,] <- c(as.numeric(exp(ests_optimx[1:9])), ests_optimx$value, 
                ests_optimx$convcode, ests_optimx$fevals,
                as.numeric(endt - startt) )  
  
  write.csv(file =  "MultipleInitialResults_NM5.csv", ests)
  
}

mutlistarts <- read.csv(file =  "MultipleInitialResults_NM.csv", row.names = NULL)[-1]

colnames(mutlistarts) <- c(paste0("Scalar", 1:inters), "TSS", "Convergence", "NO. Calls", "time")

mutlistarts[which.min((mutlistarts$TSS)),]

optim_fun(log(as.numeric(mutlistarts[which.min((mutlistarts$TSS)),1:inters])),
interventions_df = interventions_info, 
population_df = population, contacts_list = contacts, 
cumulative_obs_cases = cumulative_cases)


HESS <- hessian(optim_fun, log(as.numeric(mutlistarts[which.min((mutlistarts$TSS)),1:inters])),
        interventions_df = interventions_info, 
        population_df = population, contacts_list = contacts, 
        cumulative_obs_cases = cumulative_cases)

solve(HESS)
sqrt(diag(solve(HESS)))
#summaryRprof(filename="Profile.out")

# negative hessian may indicate poor identification of the model parameters 
#or could be just a problem of bad scaling of variables.
# assuming the derivatives defining the entries off the main diagonal are continuous at that point,
#  the special case of a smooth function at a local minimum, the Hessian will be
# a suitable covariance matrix (or inverse covariance matrix) for a Gaussian distribution.


## Running optimisation in parallel
library(parallel)


Rand_start <- function(i) {
  
  data_list <- load_data()#data_path = data_path, county = county)

  # Create Data Frames from the Data Object List
  list2env(data_list, envir = .GlobalEnv)
  
  # Remove the List containing the Data Objects
  rm(data_list)
  
  inters <- length(unique(interventions_info$policy))
  
  scalars_init <- c(runif(1, 1.5, 3), runif((inters - 1), 0, 1))
  
  startt <- Sys.time()
  ests_optimx <-  optimx(log(scalars_init), 
                         optim_fun,
                         cumulative_obs_cases = cumulative_cases, 
                         itnmax = 5000,
                         interventions_df = interventions_info, population_df = population,
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

Rand_start()#data_path = 'data/', county = 'Dublin')

nRandStarts <- 20

## Parallel 
Unlist <- function(list){
  mat <- matrix(NA, nrow = length(list), ncol = length(list[[1]]))
  for (i in 1:length(list)) {
    mat[i, ] <- unlist(list[[i]])
  }
  mat
}

ncores <- detectCores()
cl <- makeCluster(ncores - 2) 

clusterEvalQ(cl, {
  library(deSolve)
  library(optimx)
  library(tidyverse)
  library(Matrix)
})

clusterExport(cl, c("optim_fun", "simulation_SEIR_model", "SEIR_model", "getbeta", "load_data"))

startt <- Sys.time()
    parestsMat <- parLapply(cl, 1:nRandStarts, Rand_start) 
    
    estsMat <- as.data.frame(Unlist(parestsMat))
    
    colnames(estsMat) <- c(paste0("Scalar", 1:9), "TSS", "Convergence", "NO. Calls", "time")
    write.csv(file =  "MultipleInitialResults_PARALLEL.csv", estsMat)
    
Endt <- Sys.time()
    
stopCluster(cl)

Endt - startt

estsMat[which.min((estsMat$TSS)),]

optim_fun(log(as.numeric(mutlistarts[which.min((mutlistarts$TSS)),1:inters])),
          interventions_df = interventions_info, 
          population_df = population, contacts_list = contacts, 
          cumulative_obs_cases = cumulative_cases)


