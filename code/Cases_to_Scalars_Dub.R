## Load packages
library(deSolve)
#library(DEoptim)
library(optimx)

## Load the data
source('code/1_loadData.r')

## Load the get beta function
source("code/getbeta.R")

## Load the model function
source("code/function_SEIR_Model.R")

## load the simulation function 
source("code/function_SEIR_Simulation.R")

## Define objective function
#scalars <- scalars_init
nlm_fun <- function(scalars, Actual_Cc, population_df, contacts_list, interventions_df) {
  
  scalars_est <- exp(scalars)  
  
  names(scalars_est) <- unique(interventions_df$policy)
  
  Sim <- try(simulation_SEIR_model(R0t = 3.4,
                                POP = population_df,
                                contacts_ireland = contacts_list,
                                interventions = interventions_df, 
                                dt = 1,
                                dateStart = interventions_df$start[1],
                                dateEnd = cumulative_cases$date[dim(cumulative_cases)[1]],
                                scalars = scalars_est), silent = TRUE)
  
  if (class(Sim) != "try-error") {
  
  Est_Cc <- rowSums(Sim$sol_out[grepl('Cc_',names(Sim$sol_out))])
  
  ss <- try(sum((Actual_Cc - (Est_Cc))^2), silent = TRUE)
  
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

scalars_init <- c(2.45199033, 0.920820963, 0.19089041, 0.01511209, 0.32126002, 
                  0.42580101, 0.52733763,0.30954179, 0.26381850)

startt <- Sys.time()
ests_nlm <- nlm(nlm_fun,log(scalars_init),
                interventions_df = interventions_info, population_df = population,
                contacts_list = contacts, #stepmax = 0.5,
                iterlim = 1000, Actual_Cc = cumulative_cases$cases)
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

nlm_fun(log(nlm_est), interventions_df = interventions_info, 
        population_df = population, contacts_list = contacts, 
        Actual_Cc = cumulative_cases$cases)

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
c(as.numeric(exp(ests_optimx[1:9])), ests_optimx$value, ests_optimx$convcode, ests_optimx$fevals,
  endt - startt)

c(2.033784e+00, 1.093872e+00, 2.143144e-01, 3.064505e-02, 3.295509e-01, 3.654255e-01,
   5.485281e-01, 3.371615e-01, 2.035889e-01)#, 3.064904e+07, 1.000000e+01, 2.362000e+03, 1.95 hours)
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
ests <- matrix( NA, nrow = n, ncol = inters + 4 )
for (i in 1:n) {
  
  scalars_init <- c(runif(1, 1.5, 3), runif((inters - 1), 0, 1))
  
  startt <- Sys.time()
  
  ests_nlm <- nlm(nlm_fun,log(scalars_init),
                  #stepmax = 0.5,
                  interventions_df = interventions_info, population_df = population,
                  contacts_list = contacts, 
                  iterlim = 1000, Actual_Cc = cumulative_cases$cases)
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

for (i in 1:n) {
  
  scalars_init <- c(runif(1, 1.5, 3), runif((inters - 1), 0, 1))
  
  startt <- Sys.time()
  ests_optimx <-  optimx(log(scalars_init), 
                         nlm_fun,
                         Actual_Cc = cumulative_cases$cases, itnmax = 5000,
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


#summaryRprof(filename="Profile.out")
