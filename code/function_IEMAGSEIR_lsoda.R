#install.packages("deSolve")
library(deSolve)
library(tidyverse)

t <- times
params <- pars
x <- xstart

SEIR_model <- function (t, params, x) {
# Initialise the time-dependent variables
  
Ev <- as.numeric(x[grepl('Ev_',names(x))])
Ip <- as.numeric(x[grepl('Ip_',names(x))])
IA <- as.numeric(x[grepl('IA_',names(x))])
Ii <- as.numeric(x[grepl('Ii_',names(x))])
It <- as.numeric(x[grepl('It_',names(x))])
Iti <- as.numeric(x[grepl('Iti_',names(x))])
Iq <- as.numeric(x[grepl('Iq_',names(x))])
R <- as.numeric(x[grepl('R_',names(x))])
S <- as.numeric(x[grepl('S_',names(x))])

L <- as.numeric(params[[1]]["L"])
Cv <- as.numeric(params[[1]]["Cv"])
Dv <- as.numeric(params[[1]]["Dv"])
h <- as.numeric(params[[1]]["h"])
i <- as.numeric(params[[1]]["i"])
j <- as.numeric(params[[1]]["j"])
f <- as.numeric(params[[1]]["f"])
tv <- as.numeric(params[[1]]["tv"])
q <- as.numeric(params[[1]]["q"])
TT <- as.numeric(params[[1]]["TT"])
beta <- as.numeric(params[[1]]["beta"])

N_age <- as.numeric(params[[2]])
## total contacts matrix
C <- as.matrix(params[[3]])

# calculate the number of infections and recoveries between time t and t+dt
# Derivatives

dSdt <- -(S*(beta*(C%*%((as.matrix(Ip) + (h*as.matrix(IA)) + 
                                      (i*as.matrix(Ii)) + (as.matrix(It)) + 
                                      (j*as.matrix(Iti)) +
                                      (as.matrix(Iq)))/N_age))))
#dim(C)

#dim(((as.matrix(as.numeric(Ip)) + (h*as.matrix(as.numeric(IA))) + 
#        (i*as.matrix(as.numeric(Ii))) + (as.matrix(as.numeric(It))) + 
#        (as.numeric(j)*as.matrix(as.numeric(Iti))) +
#        (as.matrix(as.numeric(Iq))))/N_age))

dEvdt <- -(Ev/L) +(S*(beta*(C%*%((as.matrix(Ip) + (h*as.matrix(IA)) + 
                                    (i*as.matrix(Ii)) + (as.matrix(It)) + 
                                    (j*as.matrix(Iti)) +
                                    (as.matrix(Iq)))/N_age))))

dIpdt <- (((1 - f)*Ev)/L) - (Ip/(Cv - L))

dIAdt <- ((f*Ev)/L) - (IA/Dv)

dIidt <- ((q*Ip)/(Cv - L)) - (Ii/(Dv - Cv + L))

dItdt <- ((tv*Ip)/(Cv - L)) - (It/TT)

dItidt <- (It/TT) - (Iti/(Dv - Cv + L - TT))

dIqdt <- (((1 - q - tv)*Ip)/(Cv - L)) - (Iq/(Dv - Cv + L))

dRdt <- (IA/Dv) + (Ii/(Dv - Cv + L)) + 
  (Iq/(Dv - Cv + L)) + (Iti/(Dv - Cv + L - TT))


dxdt <- c(dSdt, dEvdt, dIpdt, dIAdt, dIidt, dItdt, dItidt, dIqdt, dRdt)

return(list(dxdt))

}

dt <- 1;                                                # Time step (days)
tmax <- 365;                                            # Time horizon (days) 366 days in 2020 cause of leap year
numSteps <- tmax/dt; 

times <- seq(from = 0, to = tmax, by = dt)
dt <- 0.01 # step size
times <- seq(from = 0, to = 30, by = dt) # times we want

# Load population information

p_age <- dubpop$propage
N_age <- dubpop$popage   

contacts_ireland <- contacts

groups <- dim(contacts_ireland[[1]])[2]
pInfected <- 0.0002
num_inf <- pInfected*sum(N_age)/groups

# create a data frame for initial values
xstart_df <- tibble(
  names = c(paste0('Ev_',1:groups), 
            paste0('Ip_',1:groups),
            paste0('IA_',1:groups),
            paste0('Ii_',1:groups), 
            paste0('It_',1:groups),
            paste0('Iti_',1:groups),
            paste0('Iq_',1:groups), 
            paste0('R_',1:groups),
            paste0('S_',1:groups)), 
  value = c(rep(0, groups),
            rep(num_inf, groups ),
            rep(0, groups),
            rep(0, groups),
            rep(0, groups),
            rep(0, groups),
            rep(0, groups),
            rep(0, groups),
            c(N_age - num_inf)) %>% 
    ceiling()
)
# turn it into a named vactor (required input by desolve)
xstart <- c(xstart_df$value)
names(xstart) <- xstart_df$names

#lambda <- as.numeric(beta)*(as.matrix(C)%*%(as.matrix(Ip) + (h*as.matrix(IA)) + 
#                                              (i*as.matrix(Ii)) + (as.matrix(It)) + 
#                                              (j*as.matrix(Iti)) + (as.matrix(Iq))));

pars <- c(4.9, 5.9, 7.0, 0.25, 0.05, 0.05, 0.5, 0.75, 0.13, 3.6, 0.05)
names(pars) <- c("L","Cv","Dv","h","i","j","f","tv","q","TT", "beta")

pars <- list(pars)

pars[[2]] <- N_age
pars[[3]] <-  contacts_ireland[[1]] + contacts_ireland[[2]] + contacts_ireland[[3]] +
  contacts_ireland[[4]]

params = pars

sol <- lsoda(xstart, times, SEIR_model, pars )

#sol_out <- lsoda(xstart, times, closed.sir.model_n_pops, params) 
