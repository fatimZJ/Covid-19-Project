## Load packages
library(deSolve)

## Load the data
source('code/1_loadData.r')

## Load the get beta function
source("code/getbeta.R")

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
  tmax <- (as.vector(dateEnd - dateStart) + 1)
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
  
  tstart_intervention <- list()
  tend_intervention <- list()
  
  for ( i in 1:length(interventions$policy)) {
    
    tstart_intervention[[i]] <- as.vector(interventions$start[i] - dateStart) + 1
    names(tstart_intervention)[[i]] <- interventions$policy[i]
    tend_intervention[[i]] <- as.vector(interventions$end[i] - dateStart) + 1
    names(tend_intervention)[[i]] <- interventions$policy[i]
  }
  
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
  
  #x = xstart
  #parms = parms
  #t = 1
  
  ## Define the SEIR model for lsoda  
  SEIR_model <- function (t, x, parms) {
    #x <- as_vector(x)
    
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
    
    tstart_intervention <- as.list(parms[["tstart_intervention"]])
    
    tend_intervention <- as.list(parms[["tend_intervention"]])

    scalars <- parms[["scalars"]]
    
    # defining the intervention
    # for( i in 1:length(tstart_intervention)){
    #   if( round(t,4) > (tstart_intervention[[i]] - 1) & round(t,4) <= tend_intervention[[i]])
    #   {
    #     INTERVENTION <- noquote(names(tstart_intervention)[i])
    #     CONSTRAINT  <- as.numeric(scalars[INTERVENTION])
    #   }
    # }

    ## rounding t
    for( i in 1:length(tstart_intervention)){
      if (round(t) == 0){
        INTERVENTION <- noquote(names(tstart_intervention)[1])
        CONSTRAINT  <- as.numeric(scalars[INTERVENTION])
      }
      else if( (round(t) > ((tstart_intervention[[i]]) - 1)) & (round(t) <= (tend_intervention[[i]])))
      {
        INTERVENTION <- noquote(names(tstart_intervention)[i])
        CONSTRAINT  <- as.numeric(scalars[INTERVENTION])
      }
    }
    
    print(INTERVENTION)
   print(round(t))
    #print(round(t,4))
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
  
  ## Solving the SEIR model equations
  sol <- lsoda(xstart, times, SEIR_model, parms)
  
  ## Defining the output
  output <- list( sol_out = as.data.frame(sol), N_age = N_age, R0t = R0t, beta = beta,
                 dateStart = dateStart, dateEnd = dateEnd,
                 tstart_intervention = tstart_intervention, 
                 tend_intervention = tend_intervention, scalars = scalars)
  
  return(output)
}

Test <- TRUE#FALSE
#scalars_test <- c(1.362216848, 0.510116460, 0.106050230, 0.008395608, 0.122922232,
#                  0.197667229, 0.248520907, 0.171967661, 0.091010276)

scalars_test <- c(1.839589e+00, 0.510116460, 0.09050230, 0.008395608, 0.122922232,
                  0.207667229, 0.308520907, 0.101967661, 0.10010276)

scalars_test <- c(6.858719e-02, 1.839589e+00, 9.594844e-02, 4.283682e-04, 2.721637e-06,
                   3.610250e-01, 2.099779e-01, 1.736524e-01, 9.436696e-02)

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
