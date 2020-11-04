#Load the data
source('code/1_loadData.r')

#install.packages("deSolve")
library(deSolve)
library(tidyverse)

t <- times
params <- pars
x <- xstart

SEIR_model <- function (t, x, parms) {
  x <- as_vector(x)
  # browser()
  #Initialise the time-dependent variables
  S <- x['S']#(x[grepl('S_',names(x))])
  Ev <- x['Ev']#(x[grepl('Ev_',names(x))])
  Ip <- x['Ip'] #(x[grepl('Ip_',names(x))])
  IA <- x['IA']#(x[grepl('IA_',names(x))])
  Ii <- x['Ii']#(x[grepl('Ii_',names(x))])
  It <- x['It']#(x[grepl('It_',names(x))])
  Iti <- x['Iti']#(x[grepl('Iti_',names(x))])
  Iq <- x['Iq']#(x[grepl('Iq_',names(x))])
  R <- x['R']#(x[grepl('R_',names(x))])
  
  
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
  beta <- (parms[["beta"]])
  
  N_age <- (parms[["N_age"]])
  ## total contacts matrix
  C <- parms[["C"]]
  
  # calculate the number of infections and recoveries between time t and t+dt
  # Derivatives
  
  dSdt <- -(S*(beta*sum(C%*%as.matrix((Ip + (h*IA) + (i*Ii) + (It) + (j*Iti) + (Iq))/N_age))))
  
  dEvdt <- -(Ev/L) +(S*(beta*sum(C%*%as.matrix((Ip + (h*IA) + (i*Ii) + (It) + (j*Iti) +
                                        (Iq))/N_age))))
  
  dIpdt <- - (Ip/(Cv - L)) + (((1 - f)*Ev)/L) 
  
  dIAdt <- - (IA/Dv) + ((f*Ev)/L) 
  
  dIidt <- - (Ii/(Dv - Cv + L)) + ((q*Ip)/(Cv - L)) 
  
  dItdt <- - (It/TT) + ((tv*Ip)/(Cv - L)) 
  
  dItidt <- - (Iti/(Dv - Cv + L - TT)) + (It/TT) 
  
  dIqdt <- - (Iq/(Dv - Cv + L)) + (((1 - q - tv)*Ip)/(Cv - L)) 
  
  dRdt <- (IA/Dv) + (Ii/(Dv - Cv + L)) + (Iq/(Dv - Cv + L)) + (Iti/(Dv - Cv + L - TT))
  
  
  dxdt <- list(c(dSdt, dEvdt, dIpdt, dIAdt, dIidt, dItdt, dItidt, dIqdt, dRdt))
  #browser()
  return(dxdt)
  
}

dt <- 0.1;                                                # Time step (days)
tmax <- 365;                                            # Time horizon (days) 366 days in 2020 cause of leap year
numSteps <- tmax/dt; 

times <- seq(from = 1, to = tmax, by = dt)

# Load population information
groups <- 1#dim(contacts_ireland[[1]])[2]
p_age <- dubpop$propage[1:groups]
N_age <- dubpop$popage[1:groups]   
contacts_ireland <- contacts#[1:groups,]

pInfected <- 0.0002
num_inf <- pInfected*sum(N_age)/groups


# create a data frame for initial values
# xstart_df <- tibble(
#   names = c(paste0('S_',1:groups),
#             paste0('Ev_',1:groups), 
#             paste0('Ip_',1:groups),
#             paste0('IA_',1:groups),
#             paste0('Ii_',1:groups), 
#             paste0('It_',1:groups),
#             paste0('Iti_',1:groups),
#             paste0('Iq_',1:groups), 
#             paste0('R_',1:groups)), 
#   value = c(c(N_age - num_inf),
#             rep(0, groups),
#             rep(num_inf, groups),
#             rep(0, groups),
#             rep(0, groups),
#             rep(0, groups),
#             rep(0, groups),
#             rep(0, groups),
#             rep(0, groups)) %>% 
#     ceiling()
# )
# # turn it into a named vactor (required input by desolve)
# dim(xstart_df)
# xstart <- c(xstart_df$value)
xstart <- c((N_age - num_inf), 0, num_inf, 0, 0, 0, 0, 0, 0)
names(xstart) <- c("S", "Ev", "Ip", "IA", "Ii", "It", "Iti", "Iq","R")

length(xstart)
#lambda <- as.numeric(beta)*(as.matrix(C)%*%(as.matrix(Ip) + (h*as.matrix(IA)) + 
#                                              (i*as.matrix(Ii)) + (as.matrix(It)) + 
#                                              (j*as.matrix(Iti)) + (as.matrix(Iq))));

#pars <- c(4.9, 5.9, 7.0, 0.25, 0.05, 0.05, 0.5, 0.75, 0.13, 3.6, 0.05)
#names(pars) <- c("L","Cv","Dv","h","i","j","f","tv","q","TT", "beta")


#pars[[2]] <- N_age
Avg_contacts <-  as.matrix(contacts_ireland[[1]] + contacts_ireland[[2]] + contacts_ireland[[3]] +
  contacts_ireland[[4]])[1:groups,]

#pars[[3]] <- Avg_contacts
#parms = pars

parms <- list(L = 4.9,Cv = 5.9,Dv = 7.0,h = 0.25,i = 0.05,j = 0.05,f = 0.5,tv=0.75,
               q = 0.13,TT = 3.6, beta = 0.05, N_age = 331515, C = Avg_contacts)
names(parms)
class(parms)
#class(parms[["L"]])
sol <- lsoda(xstart, times, SEIR_model, parms )
sol
diagnostics(sol)
#ode(func = SEIR_model, y = xstart, parms = parms, times =  times)
#sol_out <- lsoda(xstart, times, closed.sir.model_n_pops, params) 
class(sol)


#sol_put <- sol %>% as_tibble() %>%
#  mutate_all(as.numeric)

sol_out <- as.data.frame(sol)

str(sol_out)

plot(x = sol_out$time, sol_out$S, type = "l", ylim = c(0, max(sol_out$S) + 2))

lines(x = sol_out$time, (rowSums(sol_out[,4:9])), col = "red")
lines(x = sol_out$time, sol_out$R, col = "blue")

plot(rowSums(sol_out[,4:9]), type = "l")
par(mfrow = c(1,1))
#options(error=recover)
#sol %>% as.tibble() %>%
#  gather(variable,value,-time) %>%
#  ggplot(aes(x=time,y=value,color=variable))+
#  geom_line(size=2)+
#  labs(x='time (yr)',y='number of individuals')
