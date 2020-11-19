
## Load the data
source('code/1_loadData.r')

## Load the get beta function
source("code/getbeta.R")

## Load packages
library(deSolve)
library(tidyverse)


simulation_SEIR_model <- function(R0t = 3.4,
                                  #R0tpostoutbreak = 1.5,
                                  #pWorkOpen = c(0.1,0.15,0.25,0.40,0.55,0.7), # pWorkOpen: proportion of the work force that is working (will be time-varying)
                                  dateStartSchoolClosure = as.Date('2020-03-12') , # Schools closed before Intense lockdown
                                  dateStartIntenseIntervention = as.Date('2020-03-27') , #Intense intervention: lockdown
                                  dateEndIntenseIntervention = as.Date('2020-05-18'), #date we begin relaxing intense intervention
                                  dateStart = as.Date('2020-02-28'), #start date for epidemic in Ireland
                                  POP = Irlpop,
                                  numWeekStagger = c(3,6,9,12,15),
                                  contacts_ireland = contacts,
                                  scalars = c(1.10233162, 0.49031872, 0.10325630, 0.03780625, 0.11139420, 0.20539652, 0.22703778,
                                              0.17883774, 0.22092837),
                                  beta = 0.1816126,
                                  dt = 1,  # Time step (days)
                                  tmax = 225 )  #nrow(jg_dat)         # Time horizon (days))   
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
  pars <- c(4.9, 5.9, 7.0, 0.25, 0.05, 0.05, 0.5, 0.75, 0.13, 3.6)
  names(pars) <- c("L","Cv","Dv","h","i","j","f","tv","q","TT")
  
  ## Estimating Beta
  Beta <-  beta #getbeta(R0t = R0t, pars = pars, p_age = Irlpop$propage, CONTACTMATRIX = contacts_ireland)
  
  ## Defining time points at which interventions come in
  tStartSchoolClosure <- (as.vector(dateStartSchoolClosure - dateStart) + 1) #Time point to add the school closure effect
  tStartIntenseIntervention <- (as.vector(dateStartIntenseIntervention - dateStart) + 1)#Time point to add the intense lockdown effect
  tEndIntenseIntervention <- (as.vector(dateEndIntenseIntervention - dateStart) + 1)     
  tRelaxIntervention1 <- tEndIntenseIntervention + (numWeekStagger[1]*7)                               
  tRelaxIntervention2 <- tEndIntenseIntervention + (numWeekStagger[2]*7)                               
  tRelaxIntervention3 <- tEndIntenseIntervention + (numWeekStagger[3]*7)
  tRelaxIntervention4 <- tEndIntenseIntervention + (numWeekStagger[4]*7)                              
  tRelaxIntervention5 <- tEndIntenseIntervention + (numWeekStagger[5]*7)                              
  
  ## defining all parameters required for solving model equations
  parms <- list(L = pars["L"],Cv = pars["Cv"],Dv =  pars["Dv"],h = pars["h"],
                i = pars["i"],j = pars["j"],f = pars["f"],tv = pars["tv"],
                q = pars["q"],TT = pars["TT"], beta = Beta, N_age = N_age, contacts_ireland = contacts_ireland,
                tStartSchoolClosure = tStartSchoolClosure,
                tStartIntenseIntervention = tStartIntenseIntervention,
                tEndIntenseIntervention = tEndIntenseIntervention,
                tRelaxIntervention1 = tRelaxIntervention1,
                tRelaxIntervention2 = tRelaxIntervention2,
                tRelaxIntervention3 = tRelaxIntervention3,
                tRelaxIntervention4 = tRelaxIntervention4,
                tRelaxIntervention5 = tRelaxIntervention5,
                scalars = scalars)
  
  ## create a data frame for initial values (input must be a named vector,
  ##                                         order the same as the order of the equations,
  ##                                         required by solver)
  
  xstart <- c(c(N_age - num_inf - num_exp),
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
  SEIR_model <- function (t, x, parms) {
    x <- as_vector(x)
    #browser()
    # Initialise the time-dependent variables
    S <- (x[grepl('S_',names(x))])
    Ev <- (x[grepl('Ev_',names(x))])
    Ip <- (x[grepl('Ip_',names(x))])
    IA <- (x[grepl('IA_',names(x))])
    Ii <- (x[grepl('Ii_',names(x))])
    It <- (x[grepl('It_',names(x))])
    Iti <- (x[grepl('Iti_',names(x))])
    Iq <- (x[grepl('Iq_',names(x))])
    R <- (x[grepl('R_',names(x))])
    
    # Define model parameter values
    L <- parms[["L"]]
    Cv <- (parms[["Cv"]])
    Dv <- (parms[["Dv"]])
    h <- (parms[["h"]])
    i <- (parms[["i"]])
    j <- (parms[["j"]])
    f <- (parms[["f"]])
    tv <- (parms[["tv"]])
    q <- (parms[["q"]])
    TT <- (parms[["TT"]])
    
    beta <- parms[["beta"]]
    N_age <- (parms[["N_age"]])
    
    contacts_ireland <- parms[["contacts_ireland"]]
    
    tStartSchoolClosure <- parms[["tStartSchoolClosure"]]
    tStartIntenseIntervention <- parms[["tStartIntenseIntervention"]]
    tEndIntenseIntervention <- parms[["tEndIntenseIntervention"]]
    tRelaxIntervention1 <- parms[["tRelaxIntervention1"]]
    tRelaxIntervention2 <- parms [["tRelaxIntervention2"]]
    tRelaxIntervention3 <- parms[["tRelaxIntervention3"]]
    tRelaxIntervention4 <- parms[["tRelaxIntervention4"]]
    tRelaxIntervention5 <- parms [["tRelaxIntervention5"]]
    scalars <- parms[["scalars"]]
    
    if(t < tStartSchoolClosure)
    {
      INTERVENTION <- scalars[1]
      CONSTRAINT  <- diag(INTERVENTION,16,16)
    }
    # Schools closed before lockdown period
    if(t >= tStartSchoolClosure & t < tStartIntenseIntervention)
    {
      INTERVENTION <- scalars[2]
      CONSTRAINT  <- diag(INTERVENTION,16,16)
    }
    #Intense intervention- lock down
    if(t >= tStartIntenseIntervention & t < tEndIntenseIntervention)
    {
      INTERVENTION <- scalars[3]
      CONSTRAINT  <- diag(INTERVENTION,16,16)
    }
    #Relaxing interventions - Phase 1
    if(t >= tEndIntenseIntervention & t < tRelaxIntervention1)
    {
      INTERVENTION <- scalars[4]
      CONSTRAINT  <- diag(INTERVENTION,16,16)
    }
    #Relaxing interventions - Phase 2
    if(t >= tRelaxIntervention1 & t < tRelaxIntervention2)
    {
      INTERVENTION <- scalars[5]
      CONSTRAINT  <- diag(INTERVENTION,16,16)
    }
    #Relaxing interventions - Phase 3
    if(t >= tRelaxIntervention2 & t < tRelaxIntervention3)
    {
      INTERVENTION <- scalars[6]
      CONSTRAINT  <- diag(INTERVENTION,16,16)
    }
    #Relaxing interventions - Phase 4
    if(t >= tRelaxIntervention3 & t < tRelaxIntervention4)
    {
      INTERVENTION <- scalars[7]
      CONSTRAINT  <- diag(INTERVENTION,16,16)
    }
    #Relaxing interventions - Phase 5
    if(t >= tRelaxIntervention4 & t < tRelaxIntervention5)
    {
      INTERVENTION <- scalars[8]
      CONSTRAINT  <- diag(INTERVENTION,16,16)
    }
    #Post lockdown, no intervention
    if(t >= tRelaxIntervention5)
    {
      INTERVENTION <- scalars[9]
      CONSTRAINT  <- diag(INTERVENTION,16,16)
    }
    
    #print(INTERVENTION)
    #print(t)
    
    # total contacts matrix (work + school + household + other)
    C <- CONSTRAINT%*%contacts_ireland[[5]]
    
    # calculate the number of infections and recoveries between time t and t + dt
    dSdt <- -(S*(beta*(as.matrix(C)%*%as.matrix((Ip + (h*IA) + (i*Ii) + (It) +
                                                   (j*Iti) + (Iq))/N_age))))
    
    dEvdt <- -(Ev/L) +(S*(beta*(as.matrix(C)%*%as.matrix((Ip + (h*IA) + (i*Ii) +
                                                            (It) + (j*Iti) + (Iq))/N_age))))
    
    dIpdt <- - (Ip/(Cv - L)) + (((1 - f)*Ev)/L)
    
    dIAdt <- - (IA/Dv) + ((f*Ev)/L)
    
    dIidt <- - (Ii/(Dv - Cv + L)) + ((q*Ip)/(Cv - L))
    
    dItdt <- - (It/TT) + ((tv*Ip)/(Cv - L))
    
    dItidt <- - (Iti/(Dv - Cv + L - TT)) + (It/TT)
    
    dIqdt <- - (Iq/(Dv - Cv + L)) + (((1 - q - tv)*Ip)/(Cv - L))
    
    dRdt <- (IA/Dv) + (Ii/(Dv - Cv + L)) + (Iq/(Dv - Cv + L)) + (Iti/(Dv - Cv + L - TT))
    
    
    dxdt <- list(c(dSdt, dEvdt, dIpdt, dIAdt, dIidt, dItdt, dItidt, dIqdt, dRdt))
    
    return(dxdt)
  }
  
  ## Solving the equations
  sol <- lsoda(xstart, times, SEIR_model, parms)
  
  sol_out <- as.data.frame(sol)
  
  output = list( sol_out = sol_out, N_age = N_age, R0t = R0t, beta = beta,
                 dateStart = dateStart, dateEnd = dateEnd,
                 dateStartIntenseIntervention = dateStartIntenseIntervention, 
                 dateEndIntenseIntervention = dateEndIntenseIntervention, 
                 dateStartSchoolClosure = dateStartSchoolClosure,
                 tStartSchoolClosure =  tStartSchoolClosure,
                 tStartIntenseIntervention = tStartIntenseIntervention,
                 tEndIntenseIntervention = tEndIntenseIntervention,
                 tRelaxIntervention1 = tRelaxIntervention1,
                 tRelaxIntervention2 = tRelaxIntervention2,
                 tRelaxIntervention3 = tRelaxIntervention3,
                 tRelaxIntervention4 = tRelaxIntervention4, 
                 tRelaxIntervention5 = tRelaxIntervention5)
  
  return(output)
  
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
                              scalars = c(1.10233162, 0.49031872, 0.10325630, 0.03780625, 0.11139420, 0.20539652, 0.22703778,
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
#                                    scalars = c(1,1,1,1,1,1,1,1,1),
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
