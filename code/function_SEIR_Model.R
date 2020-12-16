
## Define the SEIR model for lsoda  

# x must be a named vector of initial values

SEIR_model <- function (t, x, params) {
  
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
  L <- params[["L"]]
  Cv <- (params[["Cv"]])
  Dv <- (params[["Dv"]])
  h <- (params[["h"]])
  i <- (params[["i"]])
  j <- (params[["j"]])
  f <- (params[["f"]])
  tv <- (params[["tv"]])
  q <- (params[["q"]])
  TT <- (params[["TT"]])
  
  beta <- params[["beta"]]
  N_age <- (params[["N_age"]])
  
  contacts_ireland <- params[["contacts_ireland"]]
  
  tstart_intervention <- params[["tstart_intervention"]]
  
  tend_intervention <- params[["tend_intervention"]]
  
  scalars <- params[["scalars"]]
  
  
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
