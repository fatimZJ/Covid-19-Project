## Possible sources of errors
# 1. For one, each method has discretization error, which causes a difference between them.
# 2. It's a very long timespan, and so you need very low tolerances to not have a substantial buildup of error.
# the initial condition was slightly inaccurate as well, which builds up over time to cause an error in the first digit. 

## Load the data
source('code/1_loadData.r')
par(mfrow = c(1,1))

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

## Load the get beta function 
source("code/getbeta.R")
## Load packages 
library(deSolve)
library(tidyverse)
#x <- xstart
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
  beta <- (parms[["beta"]][t])
  print(beta)
  N_age <- (parms[["N_age"]])
  
  # total contacts matrix (work + school + household + other)
  C <- parms[["C"]]
  
  # calculate the number of infections and recoveries between time t and t+dt
  
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
POP$popage[1] <- sum(POP$popage[1], 4.9e6 - Nv)
#sum(POP$popage)
p_age <- POP$propage
N_age <- POP$popage  

## setting contacts matrix to 0
cc <- rep(list(matrix(0, nrow = 16, ncol = 16)), 5)
cc[[4]] <- diag(16)
names(cc) <- c("home", "work", "school", "others", "all")

contacts_ireland <- cc # contacts

groups <- dim(contacts_ireland[[1]])[2]

pInfected <- 0.00002
num_inf <- jg_dat$Infected[1]/groups #pInfected*sum(N_age)/groups
num_exp <- jg_dat$Exposed[1]/groups
num_rec <- jg_dat$Recovered[1]/groups

Avg_contacts <- contacts_ireland[[1]] + contacts_ireland[[2]] + 
                contacts_ireland[[3]] + contacts_ireland[[4]]

## Setting time scale
dt <- 1;                           # Time step (days)
tmax <- nrow(jg_dat)#-370 #365;      # Time horizon (days) 
numSteps <- tmax/dt; 

times <- seq(from = 1, to = tmax, by = dt)
length(times)

## interpolating R
R0_int <- approx(jg_dat$Rt, xout = times)$y
length(R0_int)
#is.na(R0_int)

## Estimating Beta
pars <- c(4.9, 5.9, 7.0, 0.25, 0.05, 0.05, 0.5, 0.75, 0.13, 3.6)
names(pars) <- c("L","Cv","Dv","h","i","j","f","tv","q","TT")

beta <- getbeta(R0t = jg_dat$Rt, pars = pars, p_age = dubpop$propage, CONTACTMATRIX = contacts_ireland)
#beta <- getbeta(R0t = R0_int, pars = pars, p_age = dubpop$propage, CONTACTMATRIX = contacts_ireland)

plot(jg_dat$Beta[1:tmax], type = "l", lwd = 2)
lines(beta[1:tmax], lwd = 2, col="red")

## Defining Model parameters (As list)
parms <- list(L = 4.9,Cv = 5.9,Dv = 7.0,h = 0.25,i = 0.05,j = 0.05,f = 0.5,tv=0.75,
              q = 0.13,TT = 3.6, beta = beta, N_age = N_age, C = Avg_contacts)
#names(parms)
#class(parms)

## create a data frame for initial values (input must be a named vector, 
##                                         order the same as the order of the eqations,
##                                         required by solver)

xstart <- c(c(N_age - num_inf),# - num_rec)# - num_exp),
            rep(0, groups), #rep(num_exp, groups),# 
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

length(xstart)

## Solving the equations 
sol <- lsoda(xstart, times, SEIR_model, parms, rtol = 1e-14, atol = 1e-14)

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

plot(x = sol_out$time, jg_dat$Infected[1:tmax], lwd = 2, type = "l")
lines((rowSums(cbind(Ip, IA, Ii,It,Iti,Iq))), lwd = 2, col = "red")
#lines(xx, I, col = "blue", lwd = 2)

plot(x = sol_out$time, jg_dat$Susceptible[1:tmax], type = "l", lwd = 2)
lines(rowSums(S), col = "red", lwd = 2)

plot(x = sol_out$time, jg_dat$Recovered[1:tmax], type = "l", lwd = 2)
lines(rowSums(R), col = "red", lwd = 2)

plot(x = sol_out$time, jg_dat$Exposed[1:tmax], type = "l", lwd = 2)
lines(rowSums(Ev), col = "red", lwd = 2)

plot(abs(jg_dat$Infected[1:tmax] - (rowSums(cbind(Ip, IA, Ii,It,Iti,Iq)))), lwd = 2, type = "l")

all.equal(jg_dat$Infected , (rowSums(cbind(Ip, IA, Ii,It,Iti,Iq))))

#plot(jg_dat$Infected[1:tmax], lwd = 2, type = "l")
#xx <- seq_along(Ip[,1])*dt
#lines(xx,(rowSums(cbind(Ip, IA, Ii,It,Iti,Iq))), lwd = 2, col = "red")

##################################################################################

#type: "lsoda", "lsode", "lsodes","lsodar","vode", "daspk", "euler", "rk4", "ode23", 
#"ode45", "radau", "bdf", "bdf_d", "adams", "impAdams" or "impAdams_d" ,"iteration"

sol <- ode(xstart, times, SEIR_model, parms ,"lsode")

class(sol)

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

plot(x = sol_out$time, jg_dat$Beta, type = "l", lwd = 2)
lines(beta, lwd = 2, col="red")

plot(x = sol_out$time, jg_dat$Infected, lwd = 2, type = "l")
lines((rowSums(cbind(Ip, IA, Ii,It,Iti,Iq))), lwd = 2, col = "red")
#lines(xx, I, col = "blue", lwd = 2)

plot(x = sol_out$time, jg_dat$Susceptible, type = "l", lwd = 2)
lines(rowSums(S), col = "red", lwd = 2)

plot(x = sol_out$time, jg_dat$Recovered, type = "l", lwd = 2)
lines(rowSums(R), col = "red", lwd = 2)

plot(abs(jg_dat$Infected - (rowSums(cbind(Ip, IA, Ii,It,Iti,Iq)))), lwd = 2, type = "l")
plot(abs(jg_dat$Susceptible - (rowSums(S))), lwd = 2, type = "l")
lines(abs(jg_dat$Recovered - (rowSums(R))), lwd = 2, type = "l")
lines(abs(jg_dat$Infected - (rowSums(cbind(Ip, IA, Ii,It,Iti,Iq)))), lwd = 2, type = "l")



all.equal(jg_dat$Infected , (rowSums(cbind(Ip, IA, Ii,It,Iti,Iq))))

