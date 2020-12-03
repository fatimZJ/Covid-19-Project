
## Load the data
source('code/1_loadData.r')

## Load the get beta function
source("code/getbeta.R")

## Load packages
library(deSolve)
library(tidyverse)

num_inf <- 0.947286/16
num_exp <- 14.5344/16
groups <- 16
xstart <- c(Irlpop$popage - num_inf - num_exp,
            rep(num_exp, groups),
            rep(num_inf, groups),
            rep(0, groups),
            rep(0, groups),
            rep(0, groups),
            rep(0, groups),
            rep(0, groups),
            rep(0, groups))

names(xstart) <- c(paste0('S_',1:groups),
                   paste0('Ev_',1:groups),
                   paste0('Ip_',1:groups),
                   paste0('IA_',1:groups),
                   paste0('Ii_',1:groups),
                   paste0('It_',1:groups),
                   paste0('Iti_',1:groups),
                   paste0('Iq_',1:groups),
                   paste0('R_',1:groups))

## Define the model for lsoda
SEIR_model <- function (t_ind, x, parms) {
  
  # Initialise the time-dependent variables
  S <- x[grepl('S_',names(x))]
  Ev <- x[grepl('Ev_',names(x))]
  Ip <- x[grepl('Ip_',names(x))]
  IA <- x[grepl('IA_',names(x))]
  Ii <- x[grepl('Ii_',names(x))]
  It <- x[grepl('It_',names(x))]
  Iti <- x[grepl('Iti_',names(x))]
  Iq <- x[grepl('Iq_',names(x))]
  R <- x[grepl('R_',names(x))]
  
  # Define model parameter values
  L <- parms[["L"]]
  Cv <- parms[["Cv"]]
  Dv <- parms[["Dv"]]
  h <- parms[["h"]]
  f <- parms[["f"]]
  tv <- parms[["tv"]]
  q <- parms[["q"]]
  TT <- parms[["TT"]]
  
  beta <- parms[["beta"]]
  N_age <- parms[["N_age"]]
  linfo <- parms[["linfo"]]
  intervention_scales <- parms[["intervention_scales"]]
  
  scale_ind <- (t_ind >= linfo[[1]]) & (t_ind <= linfo[[2]])
  C <- parms[["C1"]] * ifelse( !any(scale_ind), 1L, intervention_scales[scale_ind] )
  
  # calculate the number of infections and recoveries between time t_ind and t_ind + dt
  dSdt <- -S*beta*(C%*%(Ip + h*IA + It + Iq) + parms[["H"]]%*%(Ii + Iti))/N_age
  dEvdt <- -Ev/L - dSdt
  dIpdt <- -Ip/(Cv - L) + (1 - f)*Ev/L
  dIAdt <- -IA/Dv + f*Ev/L
  dIidt <- -Ii/(Dv - Cv + L) + q*Ip/(Cv - L)
  dItdt <- -It/TT + tv*Ip/(Cv - L)
  dItidt <- -Iti/(Dv - Cv + L - TT) + It/TT
  dIqdt <- -Iq/(Dv - Cv + L) + (1 - q - tv)*Ip/(Cv - L)
  dRdt <-IA/Dv + Ii/(Dv - Cv + L) + Iq/(Dv - Cv + L) + Iti/(Dv - Cv + L - TT)
  
  list(c(dSdt, dEvdt, dIpdt, dIAdt, dIidt, dItdt, dItidt, dIqdt, dRdt))
  
}

simulation_SEIR_model <- function(pars = c(4.9, 5.9, 7.0, 0.25, 0.5, 0.75, 0.13, 3.6),
                                  dateStart = as.Date('2020-02-28'),
                                  lockdown_information = NULL,
                                  POP = Irlpop,
                                  contacts_ireland = contacts,
                                  beta = 0.1816126,
                                  startval = xstart,
                                  dt = 1,  
                                  tmax = 225 )   
{
  
  ## Load population information
  p_age <- POP$propage
  N_age <- POP$popage
  
  ## Initialising the comparments
  groups <- dim(contacts_ireland[[1]])[2]
  
  ## Setting time scale                     
  numSteps <- tmax/dt;
  times <- seq(from = 0, to = tmax, by = dt)
  dateEnd <- dateStart + (tmax - 1)
  
  ## Defining model parameters
  names(pars) <- c("L","Cv","Dv","h","f","tv","q","TT")
  
  ## Defining time points at which interventions come in
  ds <- as.numeric(dateStart)
  int_begin <- difftime(lockdown_information[[1]], dateStart, units = "days")
  int_end <- difftime(lockdown_information[[2]], dateStart, units = "days")
  linfo <- data.frame(c1 = int_begin, c2 = int_end)
  
  ## defining all parameters required for solving model equations
  parms <- list(L = pars["L"],Cv = pars["Cv"],Dv =  pars["Dv"],h = pars["h"],
                f = pars["f"],tv = pars["tv"], q = pars["q"], TT = pars["TT"], 
                beta = beta, N_age = N_age, C1 = contacts_ireland[[5]],
                H = contacts_ireland[[1]], linfo = linfo, 
                intervention_scales = lockdown_information[[3]])
  
  ## Solving the equations and returning result
  sol <- lsoda(startval, times, SEIR_model, parms)
  list( sol_out = as.data.frame(sol), N_age = N_age, beta = beta,
        dateStart = dateStart, dateEnd = dateEnd,
        lockdown_information = lockdown_information )
  
}

