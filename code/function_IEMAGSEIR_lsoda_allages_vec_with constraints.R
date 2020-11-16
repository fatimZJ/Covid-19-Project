
## Load the data
source('code/1_loadData.r')
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



## Specifying interventions and other model parameters
loadInterventions <- function(p_workopen)
{
  list(
    # constraints under a DO-NOTHING scenario 
    base =list(home = diag(1,16,16),
               work = diag(1,16,16),
               school = diag(1,16,16),
               others = diag(1,16,16)),
    # constraints under school closure + some social distancing for school-age going children but 100% workplace
    schcloseonly = list(home = diag(c(rep(1,4),rep(1,12))),
                        work = diag(1,16,16),
                        school = diag(0,16,16),
                        others = diag(c(rep(0.5,4),rep(1,12)))), 
    
    # constraints under work place distancing + schoolclosure 
    schcloseworkplacedist = list(home = diag(1,16,16),
                                 work = diag(p_workopen,16,16),
                                 school = diag(0,16,16),
                                 others = diag(c(rep(0.1,4),rep(0.1,12)))),
    
    # Post Outbeak, people still cautious 
    postoutbreak = list(home = diag(1,16,16),
                        work = diag(0.8,16,16),
                        school = diag(1.0,16,16),
                        others = diag(c(rep(0.7,4),rep(0.5,12)))))
  
}


## Load the get beta function 
source("code/getbeta.R")

## Load packages 
library(deSolve)
library(tidyverse)
x <- xstart
t = 1
SEIR_model <- function (t, x, params) {
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
  
  constraintsIntervention = loadInterventions(p_workopen = 1)
  
  if(t < tStartSchoolClosure)  
  {
    CONSTRAINT = constraintsIntervention$base
  }
  # I1:  When school winter break but before lockdown period, use 'schcloseonly'
  if(t >= tStartSchoolClosure & t < tStartIntenseIntervention) 
  {
    INTERVENTION = "schcloseonly"   
    CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
  }  
  # I2:  Intense intervention- lock down
  if(t >= tStartIntenseIntervention & t < tEndIntenseIntervention) 
  { 
    INTERVENTION = "schcloseworkplacedist"   
    CONSTRAINT = loadInterventions(p_workopen = 0.06)[[INTERVENTION]]
  }  
  if(t >= tEndIntenseIntervention & t < tRelaxIntervention1) 
  { 
    INTERVENTION = "schcloseworkplacedist"   
    CONSTRAINT = loadInterventions(p_workopen = 0.15)[[INTERVENTION]]
    #CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
  }  
  
  if(t >= tRelaxIntervention1 & t < tRelaxIntervention2) 
  {
    INTERVENTION = "schcloseworkplacedist"   
    CONSTRAINT = loadInterventions(p_workopen = 0.25)[[INTERVENTION]] 
  }  
  
  if(t >= tRelaxIntervention2 & t < tRelaxIntervention3) 
  {
    INTERVENTION = "schcloseworkplacedist"   
    CONSTRAINT = loadInterventions(p_workopen = 0.40)[[INTERVENTION]]
  }  
  
  if(t >= tRelaxIntervention3 & t < tRelaxIntervention4) 
  {
    INTERVENTION = "schcloseworkplacedist"   
    CONSTRAINT = loadInterventions(p_workopen = 0.55)[[INTERVENTION]]
  }  
  
  if(t >= tRelaxIntervention4 & t < tRelaxIntervention5) 
  {
    INTERVENTION = "schcloseworkplacedist"   
    CONSTRAINT = loadInterventions(p_workopen = 0.70)[[INTERVENTION]]
  }  
  
  if(t >= tRelaxIntervention5)  
  {
    CONSTRAINT = constraintsIntervention$postoutbreak
  }
  
  # total contacts matrix (work + school + household + other)
  C = CONSTRAINT[[1]]%*%contacts_ireland[[1]]+
    CONSTRAINT[[2]]%*%contacts_ireland[[2]]+
    CONSTRAINT[[3]]%*%contacts_ireland[[3]]+
    CONSTRAINT[[4]]%*%contacts_ireland[[4]]
  
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

## Load population information
POP <- dubpop

## Fix population
Nv <- sum(POP$popage)
#POP$popage[1] <- sum(POP$popage[1], 4.9e6 - Nv)
#sum(POP$popage)
p_age <- POP$propage
N_age <- POP$popage  

## setting contacts matrix to 0
cc <- rep(list(matrix(0, nrow = 16, ncol = 16)), 5)
cc[[4]] <- diag(16)
names(cc) <- c("home", "work", "school", "others", "all")

contacts_ireland <-  contacts #cc

Avg_contacts <- contacts_ireland[[1]] + contacts_ireland[[2]] + 
  contacts_ireland[[3]] + contacts_ireland[[4]]

## Initialising the comparments 
groups <- dim(contacts_ireland[[1]])[2]

num_inf <- 0.947286/groups    
num_exp <- 14.5344/groups

#"Mean relative difference: 0.0236159"

## Setting time scale
dt <- 1;                           # Time step (days)
tmax <- 225#nrow(jg_dat)#               # Time horizon (days) 
numSteps <- tmax/dt; 
times <- seq(from = 0, to = tmax, by = dt)
length(times)

## Estimating Beta
pars <- c(4.9, 5.9, 7.0, 0.25, 0.05, 0.05, 0.5, 0.75, 0.13, 3.6)
names(pars) <- c("L","Cv","Dv","h","i","j","f","tv","q","TT")

Beta <- getbeta(R0t = 3.65, pars = pars, p_age = dubpop$propage, CONTACTMATRIX = contacts_ireland)

#plot(jg_dat$Beta[1:tmax], type = "l", lwd = 2)
#lines(beta[1:tmax], lwd = 2, col="red")

## Defining Model parameters (As list)
#pWorkOpen = c(0.1,0.25,0.5,0.9)

dateStartSchoolClosure = as.Date('2020-03-12') # School closure
dateStartIntenseIntervention = as.Date('2020-03-27')  #Intense intervention
dateEndIntenseIntervention = as.Date('2020-05-18') #date we begin relaxing intense intervention 
dateStart = as.Date('2020-02-28') #start date for epidemic in Ireland

numWeekStagger = c(3,6,9,12,15)

tStartSchoolClosure = as.vector(dateStartSchoolClosure - dateStart)+1 #Time point to add the school closure effect
tStartIntenseIntervention = as.vector(dateStartIntenseIntervention - dateStart)+1 # #Time point to add the intense lockdown effect
tEndIntenseIntervention = as.vector(dateEndIntenseIntervention - dateStart)+1     # for pw = 0.1
tRelaxIntervention1 = tEndIntenseIntervention + numWeekStagger[1]*7                               # for pw = 0.25
tRelaxIntervention2 = tEndIntenseIntervention + numWeekStagger[2]*7                               # for pw = 0.5
tRelaxIntervention3 = tEndIntenseIntervention + numWeekStagger[3]*7
tRelaxIntervention4 = tEndIntenseIntervention + numWeekStagger[4]*7                               # for pw = 0.25
tRelaxIntervention5 = tEndIntenseIntervention + numWeekStagger[5]*7                               # for pw = 0.5
 

parms <- list(L = 4.9,Cv = 5.9,Dv = 7.0,h = 0.25,i = 0.05,j = 0.05,f = 0.5,tv = 0.75,
              q = 0.13,TT = 3.6, beta = Beta, N_age = N_age, contacts_ireland = contacts_ireland,
              tStartSchoolClosure = tStartSchoolClosure, 
              tStartIntenseIntervention = tStartIntenseIntervention,
              tEndIntenseIntervention = tEndIntenseIntervention,
              tRelaxIntervention1 = tRelaxIntervention1,
              tRelaxIntervention2 = tRelaxIntervention2,
              tRelaxIntervention3 = tRelaxIntervention3,
              tRelaxIntervention4 = tRelaxIntervention4,                              
              tRelaxIntervention5 = tRelaxIntervention5                     
) 


## create a data frame for initial values (input must be a named vector, 
##                                         order the same as the order of the eqations,
##                                         required by solver)

xstart <- c(c(N_age - num_inf - num_exp),
            rep(num_exp, groups),        
            rep(num_inf, groups), 
            rep(0, groups),
            rep(0, groups),
            rep(0, groups),
            rep(0, groups),
            rep(0, groups),
            rep(0, groups))#rep(num_rec, groups))

names(xstart) <- c(paste0('S_',1:groups),
                   paste0('Ev_',1:groups), 
                   paste0('Ip_',1:groups),
                   paste0('IA_',1:groups),
                   paste0('Ii_',1:groups), 
                   paste0('It_',1:groups),
                   paste0('Iti_',1:groups),
                   paste0('Iq_',1:groups), 
                   paste0('R_',1:groups))

length(xstart)

## Solving the equations 
sol <- lsoda(xstart, times, SEIR_model, parms)#, rtol = 1e-14, atol = 1e-14)

class(sol)
diagnostics(sol)

## Plotting the solution
sol_out <- as.data.frame(sol)
str(sol_out)

S <- sol_out[grepl('S_',names(sol_out))]
Ev <- sol_out[grepl('Ev_',names(sol_out))] 
Ip <- sol_out[grepl('Ip_',names(sol_out))] 
IA <- sol_out[grepl('IA_',names(sol_out))] 
Ii <- sol_out[grepl('Ii_',names(sol_out))]  
It <- sol_out[grepl('It_',names(sol_out))]  
Iti <- sol_out[grepl('Iti_',names(sol_out))] 
Iq <- sol_out[grepl('Iq_',names(sol_out))] 
R <- sol_out[grepl('R_',names(sol_out))] 

rowSums(sol_out[-1])

plot((rowSums(cbind(Ip, IA, Ii,It,Iti,Iq)))[1:170], lwd = 2, type = "l",
     xlab ="Time(days)", ylab = "Daily no. of infections",
     panel.first = rect(c(tStartSchoolClosure, tStartIntenseIntervention,tEndIntenseIntervention,
                          tRelaxIntervention1, tRelaxIntervention2,tRelaxIntervention3,
                          tRelaxIntervention4), -1e6, 
                        c(tStartIntenseIntervention,tEndIntenseIntervention, tRelaxIntervention1,
                          tRelaxIntervention2,tRelaxIntervention3,
                          tRelaxIntervention4,tRelaxIntervention5 ), 1e6, 
                        col=c('gray63','gray48', 'gray70', 'gray75','gray80','gray85', 'gray90'), border=NA))
#abline(v = tStartSchoolClosure, col = "red", lwd = 2)
#abline(v = tStartIntenseIntervention, col = "blue", lwd = 2)
#abline(v = tEndIntenseIntervention, col = "green", lwd = 2)


plot((rowSums(cbind(Ip, IA, Ii,It,Iti,Iq))), lwd = 2, type = "l",
     xlab ="Time(days)", ylab = "Daily no. of infections", col = "red",
     panel.first = rect(c(tStartSchoolClosure, tStartIntenseIntervention,tEndIntenseIntervention,
                          tRelaxIntervention1, tRelaxIntervention2,tRelaxIntervention3,
                          tRelaxIntervention4), -1e6, 
                        c(tStartIntenseIntervention,tEndIntenseIntervention, tRelaxIntervention1,
                          tRelaxIntervention2,tRelaxIntervention3,
                          tRelaxIntervention4,tRelaxIntervention5 ), 1e6, 
                        col=c('gray63','gray48', 'gray70', 'gray75','gray80','gray85', 'gray90'), border=NA))
lines(jg_dat$Infected[1:tmax], lwd = 2)
#abline(v = tStartSchoolClosure, col = "red", lwd = 2)
#abline(v = tStartIntenseIntervention, col = "blue", lwd = 2)
#abline(v = tEndIntenseIntervention, col = "green", lwd = 2)


#plot((rowSums(cbind(Ip, Ii,It,Iti,Iq)))[1:170], lwd = 2, type = "l",
#     xlab ="Time(days)", ylab = "Daily no. of infections")
#abline(v = tStartSchoolClosure, col = "red", lwd = 2)
#abline(v = tStartIntenseIntervention, col = "blue", lwd = 2)
#abline(v = tEndIntenseIntervention, col = "green", lwd = 2)


# plot((rowSums(cbind(Ip, IA, Ii,It,Iti,Iq))), lwd = 2, type = "l",
#      xlab ="Time(days)", ylab = "Daily no. of infections")
# lines(jg_dat$Infected[1:tmax], lwd = 2, col = "red")
# abline(v = tStartSchoolClosure, col = "red", lwd = 2)
# abline(v = tStartIntenseIntervention, col = "blue", lwd = 2)
# abline(v = tEndIntenseIntervention, col = "green", lwd = 2)


# plot(jg_dat$Infected[1:tmax], lwd = 2, type = "l",
#      xlab ="Time(days)", ylab = "Daily no. of infections")
# abline(v = tStartSchoolClosure, col = "red", lwd = 2)
# abline(v = tStartIntenseIntervention, col = "blue", lwd = 2)
# abline(v = tEndIntenseIntervention, col = "green", lwd = 2)
# 
