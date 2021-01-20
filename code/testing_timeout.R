## Load packages
library(deSolve)
library(plotly)
## Load the data
source('code/1_loadData.r')

## Load the get beta function
source("code/getbeta.R")

## Load the model function
source("code/function_SEIR_Model.R")

## load the simulation function 
#source("code/function_SEIR_Simulation.R")

R0t = 3.4
POP = population
contacts_ireland = contacts
interventions = interventions_info 
dt = 1
dateStart = interventions_info$start[1]
dateEnd = cumulative_cases$date[dim(cumulative_cases)[1]]# Time step (days)
#interventions_info$end[dim(interventions_info)[1]]

#scalars =  c(2.033784, 1.093872, 0.2143144, 0.03064505, 0.3295509, 0.3654255,
#                            0.5485281, 0.3371615, 0.2035889)

scalars <- c(2.27619341, 0.94894146, 0.21845530, 0.01503771, 0.31684429, 0.40807037, 0.49130460,
             0.34598300, 0.20890798)

#scalars <-  c(2.45199033, 0.920820963, 0.19089041, 0.01511209, 0.32126002, 
#                             0.42580101, 0.52733763,0.30954179, 0.26381850)

#scalars =  c(1, 0.3, 40, 50, 100, 3.654255e-01,
#            1, 0, 0)

names(scalars) <- unique(interventions_info$policy)

## Load population information
p_age <- POP$propage
N_age <- POP$popage

## Initialising the comparments
groups <- dim(contacts_ireland[[1]])[2]

num_inf <- 0.947286/groups
num_exp <- 14.5344/groups

## Setting time scale                     
tmax <- as.numeric(difftime(dateEnd, dateStart))
numSteps <- tmax/dt;
times <- seq(from = 0, to = tmax, by = dt)

## Defining model parameters
pars <- c(4.9, 5.9, 7.0, 0.25, 0.05, 0.05, 0.5, 0.75, 0.13, 3.6)
names(pars) <- c("L","Cv","Dv","h","i","j","f","tv","q","TT")

## Estimating Beta
Beta <- getbeta(R0t = R0t, pars = pars, p_age = POP$propage, CONTACTMATRIX = contacts_ireland)

## Rescaling contact matrices to ensure reprocity of contacts
Csym <- lapply(contacts_ireland, function(x, p_age) (x + t(x)*((p_age)%*%t(1/p_age)))/2, p_age) 

## Defining time points at which interventions come in

tstart_intervention <- as.numeric(difftime(interventions$start, dateStart, units = "days"))
names(tstart_intervention) <- interventions$policy

tend_intervention <- as.numeric(difftime(interventions$end, dateStart, units = "days"))
names(tend_intervention) <-  interventions$policy


## defining all parameters required for solving model equations
params <- list(L = pars["L"],Cv = pars["Cv"],Dv =  pars["Dv"],h = pars["h"],
               i = pars["i"],j = pars["j"],f = pars["f"],tv = pars["tv"],
               q = pars["q"],TT = pars["TT"], beta = Beta, N_age = N_age, contacts_ireland = Csym,#contacts_ireland,
               tstart_intervention = tstart_intervention,
               tend_intervention = tend_intervention,
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
                   paste0('R_',1:groups),
                   paste0('Cc_',1:groups))

## Solving the SEIR model equations
system.time(lsoda(xstart, times, SEIR_model, params, timeout = 2))

system.time(lsoda(xstart, times, SEIR_model, params))

startt <- Sys.time()

sol <- lsoda(xstart, times, SEIR_model, params, timeout = 2)

endt <- Sys.time()

endt - startt


sol <- lsoda(xstart, times, SEIR_model, params)#, timeout = 2)

attr(sol, "timeout")

sol <- as.data.frame(sol)

sol$time
Cc <- sol[grepl('Cc_',names(sol))]

length(rowSums(Cc))

plot(cumulative_cases$cases, type = "l", lwd = 2, xlab ="Time(days)",
     ylab = "Daily no. of cumulative cases")
lines(sol$time, (rowSums(Cc)),lwd=2,col='tomato')

sum((cumulative_cases$cases - (rowSums(Cc)))^2)

