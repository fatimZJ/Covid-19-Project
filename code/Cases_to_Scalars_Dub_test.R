## Load packages
library(deSolve)
library(DEoptim)
library(optimx)

## Load packages
library(deSolve)

## Load the data
source('code/1_loadData.r')

## Load the get beta function
source("code/getbeta.R")

## Load the model function
source("code/function_SEIR_Model.R")

## load the simulation function 
source("code/function_SEIR_Simulation.R")

## Define objective function

nlm_fun <- function(scalars, Actual_Cc, population_df, contacts_list, interventions_df) {
  
  scalars_est <- exp(scalars)  
  
  names(scalars_est) <- unique(interventions_df$policy)
  
  Sim <- try(simulation_SEIR_model(R0t = 3.4,
                                POP = population_df,
                                contacts_ireland = contacts_list,
                                interventions = interventions_df, 
                                dt = 1,
                                dateStart = interventions_df$start[1],
                                dateEnd = cumulative_cases$date[dim(cumulative_cases)[1]],# Time step (days)
                                scalars = scalars_est), silent = TRUE)
  
  if (class(Sim) != "try-error") {
  
  Est_Cc <- rowSums(Sim$sol_out[grepl('Cc_',names(Sim$sol_out))])
  
  tss <- sum((Actual_Cc - (Est_Cc))^2) } else {
    
    tss <- Inf
  }
  
  tss
  
}

## Nlm

scalars_init <- c(2.45199033, 0.920820963, 0.19089041, 0.01511209, 0.32126002, 
                  0.42580101, 0.52733763,0.30954179, 0.26381850)

startt <- Sys.time()
ests_nlm <- nlm(nlm_fun,log(scalars_init),
                interventions_df = interventions_info, population_df = population,
                contacts_list = contacts, stepmax = 0.5,
                iterlim = 1000, Actual_Cc = cumulative_cases$cases)
endt <- Sys.time()


scalars <- log(scalars_init)

exp(ests_nlm$estimate)
ests_nlm$min

## Optimx- NM

startt <- Sys.time()
ests_optimx <-  optimx(log(scalars_init), 
                       nlm_fun,
                       #lower = 0, upper = 3,
                       Actual_Cc = cumulative_cases$cases, itnmax = 5000,
                       interventions_df = interventions_info, population_df = population,
                       contacts_list = contacts, 
                       method = c("Nelder-Mead"),
                       control = list(maximize = FALSE, kkt = FALSE))
endt <- Sys.time()


endt - startt
c(as.numeric(exp(ests_optimx[1:9])), ests_optimx$value, ests_optimx$convcode, ests_optimx$fevals)

## DEoptim

lower <- rep(0, length(scalars_init))
upper <- rep(4, length(scalars_init))

startt <- Sys.time()
ests_DEoptim <- DEoptim(nlm_fun, Actual_Cc = cumulative_cases$cases,
                        lower, upper, DEoptim.control(itermax = 200, reltol = 1e-8,
                                                      steptol = 50))

endt <- Sys.time()

################################################################################

## Multiple random starts 

## Nlm
n <- 20 

inters <- length(unique(interventions_info$policy))
ests <- matrix( NA, nrow = n, ncol = inters + 3 )
for (i in 1:n) {
  
  scalars_init <- runif(inters, 0, 1)
  
  ests_nlm <- nlm(nlm_fun,log(scalars_init),
                  stepmax = 0.5,interventions_df = interventions_info, population_df = population,
                  contacts_list = contacts, 
                  iterlim = 1000, Actual_Cc = cumulative_cases$cases)
  
  ests[i,] <- c(exp(ests_nlm$estimate), ests_nlm$minimum, ests_nlm$code, ests_nlm$iterations)
  
  write.csv(file =  "MultipleInitialResults_nlm.csv", ests)
  
}



## Optimx- NM

n <- 20 

inters <- length(unique(interventions_info$policy))
ests <- matrix( NA, nrow = n, ncol = inters + 4 )

for (i in 1:n) {
  
  scalars_init <- runif(inters, 0, 1)
  
  startt <- Sys.time()
  ests_optimx <-  optimx(log(scalars_init), 
                         nlm_fun,
                         Actual_Cc = cumulative_cases$cases, itnmax = 5000,
                         interventions_df = interventions_info, population_df = population,
                         contacts_list = contacts, 
                         method = c("Nelder-Mead"),
                         control = list(maximize = FALSE, kkt = FALSE))
  endt <- Sys.time()
  
  ests[i,] <- length(c(as.numeric(exp(ests_optimx[1:9])), ests_optimx$value, 
                       ests_optimx$convcode, ests_optimx$fevals,
                       as.numeric(endt - startt) )  )
  
  write.csv(file =  "MultipleInitialResults_NM.csv", ests)
  
}


#summaryRprof(filename="Profile.out")
