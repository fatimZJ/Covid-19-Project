
## Load the data
source('code/1_loadData.r')

## Load the get beta function
source("code/getbeta.R")

## Load packages
library(deSolve)
library(tidyverse)


simulation_SEIR_model <- function(R0t = 3.4, pars = c(4.9, 5.9, 7.0, 0.25, 0.5, 0.75, 0.13, 3.6),
                                  dateStart = as.Date('2020-02-28'), #start date for epidemic in Ireland
                                  lockdown_information = NULL,
                                  POP = Irlpop,
                                  numWeekStagger = c(3,6,9,12,15),
                                  contacts_ireland = contacts,
                                  beta = 0.1816126,
                                  dt = 1,  # Time step (days)
                                  tmax = 225 )  # Time horizon (days))   
{
  
  ## Load population information
  p_age <- POP$propage
  N_age <- POP$popage
  
  ## Initialising the comparments
  groups <- dim(contacts_ireland[[1]])[2]
  
  num_inf <- 0.947286/groups
  num_exp <- 14.5344/groups
  
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
  
  ## create a data frame for initial values (input must be a named vector,
  ##                                         order the same as the order of the equations,
  ##                                         required by solver)
  
  xstart <- c(N_age - num_inf - num_exp,
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
  
  ## Solving the equations and returning result
  sol <- lsoda(xstart, times, SEIR_model, parms)
  list( sol_out = as.data.frame(sol), N_age = N_age, R0t = R0t, beta = beta,
        dateStart = dateStart, dateEnd = dateEnd,
        lockdown_information = lockdown_information )
  
}

Test <- FALSE

if (Test) {
  Base <- simulation_SEIR_model(R0t = 3.4,
                                dateStartSchoolClosure = as.Date('2020-03-12') , #schools closed before imtense lockdown
                                dateStartIntenseIntervention = as.Date('2020-03-27') , #Intense intervention: starts at Wuhan Lockdown
                                dateEndIntenseIntervention = as.Date('2020-05-18'), #date we begin relaxing intense intervention
                                dateStart = as.Date('2020-02-28'), #start date for epidemic in Ireland
                                POP = Irlpop,
                                numWeekStagger = c(3,6,9,12,15),
                                contacts_ireland = contacts,
                                intervention_scales = c(1.10233162, 0.49031872, 0.10325630, 0.03780625, 0.11139420, 0.20539652, 0.22703778,
                                            0.17883774, 0.22092837), # l2 norm for 2019 pop
                                dt = 0.1,  
                                tmax = 225 ) 
  
  # doNothing <- simulation_SEIR_model(R0t = 3.65,
  #                                    dateStartSchoolClosure = as.Date('2020-03-10') , #schools closed before imtense lockdown
  #                                    dateStartIntenseIntervention = as.Date('2020-03-10') , #Intense intervention: starts at Wuhan Lockdown
  #                                    dateEndIntenseIntervention = as.Date('2020-05-10'), #date we begin relaxing intense intervention
  #                                    dateStart = as.Date('2020-02-28'), #start date for epidemic in Ireland
  #                                    POP = Irlpop,
  #                                    numWeekStagger = c(0,0,0,0,0),
  #                                    contacts_ireland = contacts,
  #                                    intervention_scales = c(1,1,1,1,1,1,1,1,1),
  #                                    dt = 1,  
  #                                    tmax = 225 ) 
  
  ## Plotting the solution
  
  S <- Base$sol_out[grepl('S_',names(Base$sol_out))]
  Ev <- Base$sol_out[grepl('Ev_',names(Base$sol_out))]
  Ip <- Base$sol_out[grepl('Ip_',names(Base$sol_out))]
  IA <- Base$sol_out[grepl('IA_',names(Base$sol_out))]
  Ii <- Base$sol_out[grepl('Ii_',names(Base$sol_out))]
  It <- Base$sol_out[grepl('It_',names(Base$sol_out))]
  Iti <- Base$sol_out[grepl('Iti_',names(Base$sol_out))]
  Iq <- Base$sol_out[grepl('Iq_',names(Base$sol_out))]
  R <- Base$sol_out[grepl('R_',names(Base$sol_out))]
  
  ## Adapted from https://gist.github.com/jonocarroll/b17ce021b0637a31f584ed08a1fbe733
  read.tcsv = function(file, header=TRUE, sep=",", ...) {
    
    n = max(count.fields(file, sep=sep), na.rm=TRUE)
    x = readLines(file)[-1]
    
    .splitvar = function(x, sep, n) {
      var = unlist(strsplit(x, split=sep))
      length(var) = n
      return(var)
    }
    
    x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
    x = apply(x, 1, paste, collapse=sep)
    ## empty strings are converted to NA
    out = read.csv(text=x, sep=sep, header=header, na.strings = "", ...)
    return(out)
  }
  
  jg_dat <- read.tcsv("data/dat_seir_code.csv")
  jg_dat$Date <- as.Date(jg_dat$Date, format = "%a %d %b %Y") # Reformat date column
  
  
  
  plot(jg_dat$Infected[1:225], type = "l", lwd = 2, xlab ="Time(days)",
       ylab = "Daily no. of infections",
       panel.first = rect(c(Base$tStartSchoolClosure, Base$tStartIntenseIntervention,Base$tEndIntenseIntervention,
                            Base$tRelaxIntervention1, Base$tRelaxIntervention2,Base$tRelaxIntervention3,
                            Base$tRelaxIntervention4), -1e6,
                          c(Base$tStartIntenseIntervention,Base$tEndIntenseIntervention, Base$tRelaxIntervention1,
                            Base$tRelaxIntervention2,Base$tRelaxIntervention3,
                            Base$tRelaxIntervention4,Base$tRelaxIntervention5 ), 1e6,
                          col=c('gray63','gray48', 'gray70', 'gray75','gray80','gray85', 'gray90'), border=NA))
  
  lines(Base$sol_out$time[-1], (rowSums(cbind(Ip, IA, Ii,It,Iti,Iq)))[-1],lwd=2,col='tomato')
  legend(120, 12000,legend = c("JG's model", "Our Model"),
         col = c("black", "red"), bty = 'n',lty = c(1,1),lwd = c(2,2), cex = 1)
}
