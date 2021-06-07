## Define the model for lsoda
SEIR_model_D <- function (t_ind, x, parms) {
  
  #browser()
  # Initialise the time-dependent variables
  S <- x[1:16]
  Ev <- x[17:32]
  Ip <- x[33:48]
  IA <- x[49:64]
  Ii <- x[65:80]
  It <- x[81:96]
  Iti <- x[97:112]
  Iq <- x[113:128]
  R <- x[129:144]
  Cc <- x[145:160]
  
  # Define model parameter values
  L <- parms[["L"]]
  Cv <- parms[["Cv"]]
  #denom_1 <- parms[["denom_1"]]
  #denom_2 <- parms[["denom_2"]]
  Dv <- parms[["Dv"]]
  f <- parms[["f"]]
  k <- parms[["k"]]
  tv <- parms[["tv"]]
  q <- parms[["q"]]
  h <- parms[["h"]]
  TT <- parms[["TT"]]
  beta <- parms[["beta"]]
  linfo <- parms[["linfo"]]
  N_age <- parms[["N_age"]]
  
  # applying the intervention to the contact matrices
  C <- parms[["C1"]] * parms[["intervention_scales"]][(t_ind >= linfo[[1]]) & (t_ind < (linfo[[2]] + 1))]
  
  # calculate the number of infections and recoveries between time t_ind and t_ind + dt
  #dSdt <- -S*parms[["beta"]]*C%*%(Ip + parms[["h"]]*IA + It + Iq + k*Ii + k*Iti)/parms[["N_age"]]
  #dEvdt <- -Ev/L - dSdt
  #dIpdt <- -Ip/denom_1 + (1 - f)*Ev/L
  #dIAdt <- -IA/Dv + f*Ev/L
  #dIidt <- -Ii/denom_2 + q*Ip/denom_1
  #dItdt <- -It/TT + tv*Ip/denom_1
  #dItidt <- -Iti/(denom_2 - TT) + It/TT
  #dIqdt <- -Iq/denom_2 + (1 - q - tv)*Ip/denom_1
  #dRdt <-IA/Dv + Ii/denom_2 + Iq/denom_2 + Iti/(denom_2 - TT)
  
  dSdt <- -S*(beta*C %*% ((Ip + (h*IA) + (k*Ii) + (It) + (k*Iti) + (Iq))/N_age))
  dEvdt <- -(Ev/L) - dSdt
  dIpdt <- -(Ip/(Cv - L)) + (((1 - f)*Ev)/L)
  dIAdt <- -(IA/Dv) + ((f*Ev)/L)
  dIidt <- -(Ii/(Dv - Cv + L)) + ((q*Ip)/(Cv - L))
  dItdt <- -(It/TT) + ((tv*Ip)/(Cv - L))
  dItidt <- -(Iti/(Dv - Cv + L - TT)) + (It/TT)
  dIqdt <- -(Iq/(Dv - Cv + L)) + (((1 - q - tv)*Ip)/(Cv - L))
  dRdt <- (IA/Dv) + (Ii/(Dv - Cv + L)) + (Iq/(Dv - Cv + L)) + (Iti/(Dv - Cv + L - TT))
  dCcdt <- It/TT
  
  list(c(dSdt, dEvdt, dIpdt, dIAdt, dIidt, dItdt, dItidt, dIqdt, dRdt, dCcdt))
  
}

## Define the model for lsoda for when we want to isolate over 75s
SEIR_model_I <- function (t_ind, x, parms) {
  
  # Initialise the time-dependent variables
  S <- x[1:16]
  Ev <- x[17:32]
  Ip <- x[33:48]
  IA <- x[49:64]
  Ii <- x[65:80]
  It <- x[81:96]
  Iti <- x[97:112]
  Iq <- x[113:128]
  R <- x[129:144]
  Cc <- x[145:160]
  
  # Define model parameter values
  L <- parms[["L"]]
  Cv <- parms[["Cv"]]
  #denom_1 <- parms[["denom_1"]]
  #denom_2 <- parms[["denom_2"]]
  Dv <- parms[["Dv"]]
  f <- parms[["f"]]
  k <- parms[["k"]]
  tv <- parms[["tv"]]
  q <- parms[["q"]]
  h <- parms[["h"]]
  TT <- parms[["TT"]]
  beta <- parms[["beta"]]
  linfo <- parms[["linfo"]]
  N_age <- parms[["N_age"]]
  
  if (t_ind >= 277) {
    C1 <- parms[["C2"]]
    beta <- parms[["beta2"]]
  }
  else {
    C1 <- parms[["C1"]]
    beta <- parms[["beta"]]
  }
  
  C <- C1 * parms[["intervention_scales"]][(t_ind >= linfo[[1]]) & (t_ind < (linfo[[2]] + 1))]
  
  # calculate the number of infections and recoveries between time t_ind and t_ind + dt
  #dSdt <- -S*beta*C%*%(Ip + parms[["h"]]*IA + It + Iq + k*Ii + k*Iti)/parms[["N_age"]]
  #dEvdt <- -Ev/L - dSdt
  #dIpdt <- -Ip/denom_1 + (1 - f)*Ev/L
  #dIAdt <- -IA/Dv + f*Ev/L
  #dIidt <- -Ii/denom_2 + q*Ip/denom_1
  #dItdt <- -It/TT + tv*Ip/denom_1
  #dItidt <- -Iti/(denom_2 - TT) + It/TT
  #dIqdt <- -Iq/denom_2 + (1 - q - tv)*Ip/denom_1
  #dRdt <-IA/Dv + Ii/denom_2 + Iq/denom_2 + Iti/(denom_2 - TT)
  
  dSdt <- -S*(beta*C %*% ((Ip + (h*IA) + (k*Ii) + (It) + (k*Iti) + (Iq))/N_age))
  dEvdt <- -(Ev/L) - dSdt
  dIpdt <- -(Ip/(Cv - L)) + (((1 - f)*Ev)/L)
  dIAdt <- -(IA/Dv) + ((f*Ev)/L)
  dIidt <- -(Ii/(Dv - Cv + L)) + ((q*Ip)/(Cv - L))
  dItdt <- -(It/TT) + ((tv*Ip)/(Cv - L))
  dItidt <- -(Iti/(Dv - Cv + L - TT)) + (It/TT)
  dIqdt <- -(Iq/(Dv - Cv + L)) + (((1 - q - tv)*Ip)/(Cv - L))
  dRdt <- (IA/Dv) + (Ii/(Dv - Cv + L)) + (Iq/(Dv - Cv + L)) + (Iti/(Dv - Cv + L - TT))
  dCcdt <- It/TT
  
  list(c(dSdt, dEvdt, dIpdt, dIAdt, dIidt, dItdt, dItidt, dIqdt, dRdt, dCcdt))
  
}

SEIR_model_simulation <- function(pars,
                                  dateStart = as.Date('2020-02-29'),
                                  lockdown_information = NULL,
                                  POP = population,
                                  contacts_ireland = contacts,
                                  beta = 0.1816126,
                                  startval,
                                  dt = 1,  
                                  tmax = 225,
                                  isolated = FALSE,
                                  isolated_contacts = contacts,
                                  isolated_beta = 0.1816126)   
{
  
  ## Load population information
  N_age <- POP$popage
  p_age <- POP$propage
  
  ## Initialising the compartments
  groups <- dim(contacts_ireland[[1]])[2]
  
  ## Setting time scale                     
  numSteps <- tmax/dt;
  times <- seq(from = 0, to = tmax, by = dt)
  dateEnd <- dateStart + (tmax - 1)
  
  ## Defining time points at which interventions come in
  linfo <- data.frame(c1 = difftime(lockdown_information[[1]], dateStart, units = "days"),
                      c2 = difftime(lockdown_information[[2]], dateStart, units = "days"))
  
  Csym <- lapply(contacts_ireland, function(x, p_age) (x + t(x)*((p_age) %*% t(1/p_age)))/2, p_age) 
  
  ## defining all parameters required for solving model equations
  parms <- list(L = pars[1], Cv = pars[2], Dv = pars[3],
                h = pars[4], f = pars[5], tv = pars[6], 
                q = pars[7], k = pars[8], TT = pars[9], beta = beta, 
                N_age = N_age, C1 = Csym[[5]], linfo = linfo, 
                intervention_scales = lockdown_information[[3]],
                C2 = isolated_contacts[[5]], beta2 = isolated_beta)
  
  ## Solving the equations and returning result
  if (isolated) { sol <- lsoda(startval, times, SEIR_model_I, parms) }
  else { sol <- lsoda(startval, times, SEIR_model_D, parms) }
  
  list( solution = as.data.frame(sol), N_age = N_age, beta = beta,
        dateStart = dateStart, dateEnd = dateEnd,
        lockdown_information = lockdown_information )
  
}

