## Load packages
library(deSolve)

## Load the data
source('code/1_loadData.r')

## Load the get beta function
source("code/getbeta.R")

## Define the SEIR model for lsoda  
SEIR_model <- function (t, x, parms) {
  #x <- as_vector(x)
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
  Cc <- (x[grepl('Cc_',names(x))])
  
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
  
  tstart_intervention <- parms[["tstart_intervention"]]
  
  tend_intervention <- parms[["tend_intervention"]]
  
  scalars <- parms[["scalars"]]
  
 
  
  intervention_ind <- (t >= tstart_intervention) & (t < (tend_intervention + 1))  ## rounding t or not rounding t?
  
  INTERVENTION <- noquote(names(intervention_ind[intervention_ind == 1]))
  
  CONSTRAINT <- ifelse( !any(intervention_ind), scalars["No Intervention"], scalars[INTERVENTION])
  
  print(INTERVENTION)
  print(CONSTRAINT)
 # print(t)
  
  # applying the intervention to the contact matrices
  C <-  diag(CONSTRAINT, length(N_age))%*%contacts_ireland[[5]]
  
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
  
  dCcdt <- (It/TT)
  
  dxdt <- list(c(dSdt, dEvdt, dIpdt, dIAdt, dIidt, dItdt, dItidt, dIqdt, dRdt,dCcdt))
  
  return(dxdt)
}


simulation_SEIR_model <- function(R0t = 3.4,
                                  POP = population,
                                  contacts_ireland = contacts,
                                  interventions = interventions_info, 
                                  dt = 1, 
                                  dateStart = interventions_info$start[1],
                                  dateEnd = interventions_info$end[dim(interventions_info)[1]],# Time step (days)
                                  scalars)    
{
  
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
  parms <- list(L = pars["L"],Cv = pars["Cv"],Dv =  pars["Dv"],h = pars["h"],
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
  sol <- lsoda(xstart, times, SEIR_model, parms)
  
  ## Defining the output
  output <- list( sol_out = as.data.frame(sol) , N_age = N_age, R0t = R0t, beta = beta,
                  dateStart = dateStart, dateEnd = dateEnd,
                  tstart_intervention = tstart_intervention, 
                  tend_intervention = tend_intervention, scalars = scalars)
  
  return(output)
}

Test <- TRUE #
#scalars_test <- c(1.362216848, 1.50116460, 0.23050230, 0.0105608, 0.0122922232,
#                  0.597667229, 0.30907, 0.251967661, 1)

#scalars_test <- c(1.839589e+00, 0.510116460, 0.09050230, 0.008395608, 0.122922232,
#                  0.207667229, 0.308520907, 0.101967661, 0.10010276)

#scalars_test <- c(0.003476602, 1.929626439, 0.101822647, 0.001159886, 0.003668096,
#                     0.338307651, 0.214795326, 0.176496238, 0.097322440)
                   
#scalars_test <- c(3.36625700, 0.35436691, 0.30468075, 0.01002319, 0.25632503, 
#                  0.52871643, 0.04051897, 0.17398844, 0.06356881)

#scalars_test <- c(3.234072e+00, 2.348389e-01, 2.136190e-01, 3.272638e-02, 2.050354e-01, 7.639543e-01,
#4.721503e-01, 3.246413e-01, 5.988597e-06)

#scalars_test <- c(0.003476602, 1.929626439, 0.101822647, 0.001159886, 0.003668096,
#                 0.338307651, 0.214795326, 0.176496238, 0.097322440)

## NM Random initial values
#scalars_test <- c(0.7811568, 0.9227932, 0.1053569, 0.02633053, 0.04398501, 0.2400942,
#0.2515661, 0.1412688, 0.1473584)

## NM
#scalars_test <- c(1.391473, 0.5718681, 0.1083594, 0.01474348, 0.08742271, 0.2026824, 0.2653538,
# 0.1644852, 0.1034355)

#scalars_test <- c(1.332476753,	0.595615258,	0.09840976,	0.162825863,	0.172910861,	0.074046275,	
#  0.144918579,	0.321139668,	0.217734612)


#scalars_test <- c(1.839589e+00, 0.510116460, 0.09050230, 0.008395608, 0.122922232,
#                  0.207667229, 0.308520907, 0.101967661, 0.10010276)

#scalars_test <- c(6.858719e-02, 1.839589e+00, 9.594844e-02, 4.283682e-04, 2.721637e-06,
#                   3.610250e-01, 2.099779e-01, 1.736524e-01, 9.436696e-02)

scalars_test <- c(1.362216848, 0.510116460, 0.106050230, 0.008395608, 0.122922232,
                  0.197667229, 0.248520907, 0.171967661, 0.091010276)

names(scalars_test) <- unique(interventions_info$policy)


#scalars <- scalars_test

if (Test) {
  
  Base <- simulation_SEIR_model(scalars = scalars_test) 
  
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
  Cc <- Base$sol_out[grepl('Cc_',names(Base$sol_out))]
  
  plot(cumulative_cases$cases, type = "l", lwd = 2, xlab ="Time(days)",
       ylab = "Daily no. of cumulative cases",
       panel.first = rect(c(Base$tstart_intervention[[2]], Base$tstart_intervention[[3]],
                            Base$tstart_intervention[[4]], Base$tstart_intervention[[5]], 
                            Base$tstart_intervention[[6]],Base$tstart_intervention[[7]],
                            Base$tstart_intervention[[8]], Base$tstart_intervention[[9]]), -1e6,
                          c(Base$tend_intervention[[2]], Base$tend_intervention[[3]],
                            Base$tend_intervention[[4]], Base$tend_intervention[[5]], 
                            Base$tend_intervention[[6]],Base$tend_intervention[[7]],
                            Base$tend_intervention[[8]], Base$tend_intervention[[9]]), 1e6,
                          col=c('gray88','gray48', 'gray70', 'gray75','gray80','gray85', 'gray70','gray54'), border=NA))
  
  lines(Base$sol_out$time[-1], (rowSums(Cc))[-1],lwd=2,col='tomato')
  legend(40, 24000,legend = c("Actual data", "Our Model"),
         col = c("black", "red"), bty = 'n',lty = c(1,1),lwd = c(2,2), cex = 1)
  
}

#plot(Base$sol_out$time, (rowSums(Cc)),lwd=2,col='tomato', type = "l")
