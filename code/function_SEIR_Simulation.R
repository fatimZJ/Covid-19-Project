
simulation_SEIR_model <- function(R0t, POP, contacts_ireland, interventions, 
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
  sol <- lsoda(xstart, times, SEIR_model, params)
  
  ## Defining the output
  output <- list( sol_out = as.data.frame(sol) , N_age = N_age, R0t = R0t, beta = beta,
                  dateStart = dateStart, dateEnd = dateEnd,
                  tstart_intervention = tstart_intervention, 
                  tend_intervention = tend_intervention, scalars = scalars)
  
  return(output)
}