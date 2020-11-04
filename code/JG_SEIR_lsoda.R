#Load the data
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

#install.packages("deSolve")
library(deSolve)
library(tidyverse)

## Define function to return derivatives 
SEIR_model <- function (t, x, parms) {
  x <- as_vector(x)
  # browser()
  # Initialise the time-dependent variables
  S <- x['S']
  Ev <- x['Ev']
  Ip <- x['Ip'] 
  IA <- x['IA']
  Ii <- x['Ii']
  It <- x['It']
  Iti <- x['Iti']
  Iq <- x['Iq']
  R <- x['R']
  
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
  #print(beta)
  Nv <- (parms[["Nv"]])
  
  # calculate the number of infections and recoveries between time t and t+dt
  
  dSdt <- -(S*beta*((Ip + (h*IA) + (i*Ii) + (It) + (j*Iti) + (Iq))/Nv))
  
  dEvdt <- -(Ev/L) +(S*beta*((Ip + (h*IA) + (i*Ii) + (It) + (j*Iti) +
                                                  (Iq))/Nv))
  
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

## vectorising beta
dt <- 1;                      # Time step (days)
tmax <- nrow(jg_dat);                    # Time horizon (days) 
numSteps <- tmax/dt; 

times <- seq(from = 1, to = tmax, by = dt)

Nv <- (4.9*(10^6))
Infected <- 0.00002
num_inf <- 1.1222#1.1#0.947286#pInfected*Nv (1.55322# 20)
num_exp <- 22#14.5344

xstart <- c((Nv - num_inf - num_exp), num_exp, num_inf, 0, 0, 0, 0, 0, 0)
names(xstart) <- c("S", "Ev", "Ip", "IA", "Ii", "It", "Iti", "Iq","R")

length(xstart)

Avg_contacts <-  as.matrix(contacts_ireland[[1]] + contacts_ireland[[2]] + contacts_ireland[[3]] +
                             contacts_ireland[[4]])[1:groups,]

beta <- getbeta(R0t = jg_dat$Rt, pars = pars, p_age = dubpop$propage, CONTACTMATRIX = contacts_ireland)

parms <- list(L = 4.9,Cv = 5.9,Dv = 7.0,h = 0.25,i = 0.05,j = 0.05,f = 0.5,tv=0.75,
              q = 0.13,TT = 3.6, beta = jg_dat$Beta, Nv = Nv)
names(parms)
class(parms)

sol <- lsoda(xstart, times, SEIR_model, parms )

class(sol)
#diagnostics(sol)

sol_out <- as.data.frame(sol)

str(sol_out)

plot(x = sol_out$time, jg_dat$Infected, type = "l", lwd = 2)
lines(rowSums(sol_out[,4:9]), col = "red", lwd = 2)
all.equal(jg_dat$Infected, rowSums(sol_out[,4:9]))
#type: "lsoda", "lsode", "lsodes","lsodar","vode", "daspk", "euler", "rk4", "ode23", 
#"ode45", "radau", "bdf", "bdf_d", "adams", "impAdams" or "impAdams_d" ,"iteration"

sol <- ode(xstart, times, SEIR_model, parms ,"impAdams_d")

sol_out <- as.data.frame(sol)

str(sol_out)

plot(x = sol_out$time, jg_dat$Infected, type = "l", lwd = 2)
lines(rowSums(sol_out[,4:9]), col = "red", lwd = 2)

plot(abs(jg_dat$Infected - rowSums(sol_out[,4:9])), type = "l", lwd = 2)

all.equal(jg_dat$Infected, rowSums(sol_out[,4:9])) # 0.01853623