## SEIR Model simulation
## Inputs 
# R0t: the basic reproduction number
# pars: names vector of model parameters ("L","Cv","Dv","h","i","j","f","tv","q","TT")
# POP: population age group structure 
# contacts_ireland: list of contact matrices
# interventions: start, end date and policy or no policy
# dt: time step  
# dateStart: start date of simulation, set to first start date on interventions
# dateEnd: end date of simulation
# scalars: vector of scalars, length should be same as the number of interventions


simulation_SEIR_model <- function(R0t, pars, POP, contacts_ireland, interventions, 
                                  dt, dateStart, dateEnd, scalars)    
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
  
  ## Defining model parameters (as per document)
  if (missing(pars)){
  pars <- c(3.6, 5.8, 13.0, 0.55, 0.05, 0.05, 0.21, 0.8, 0.1, 2)
  names(pars) <- c("L","Cv","Dv","h","i","j","f","tv","q","TT")
  }
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
  sol <- lsoda(xstart, times, SEIR_model, params, timeout = 5)
  
  ## Defining the output
  output <- list( sol_out = as.data.frame(sol) , N_age = N_age, R0t = R0t, Beta = Beta,
                  dateStart = dateStart, dateEnd = dateEnd,
                  tstart_intervention = tstart_intervention, 
                  tend_intervention = tend_intervention, scalars = scalars)
  
  return(output)
}