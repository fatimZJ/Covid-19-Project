### get beta

### Using R0 and parameters
getbeta1 = function(R0t, pars)
{
  
  # Extract Parameters
  h <- pars["h"]
  i <- pars["i"]
  j <- pars["j"]
  Nv <- pars["Nv"]
  L <- pars["L"]
  Cv <- pars["Cv"]
  Dv <- pars["Dv"]
  f <- pars["f"]
  q <- pars["q"]
  tv <- pars["tv"]
  TT <- pars["TT"]
  Nv <- pars["Nv"]  
  g <- 1 - f
  
  # Compartment Contributions
  IaC <- prod(f, Dv, h)
  IpC <- prod(g, Cv - L)
  IiC <- prod(g, q, Dv - Cv + L, i)
  It1C <- prod(g, tv, TT)
  It2C <- prod(g, tv, Dv - Cv + L - TT, j)
  InC <- prod(g, 1 - q - tv, Dv - Cv + L)
  allC <- sum(IaC, IpC, IiC, It1C, It2C, InC)
  
  R0t / allC
  
}

### Using compartments and parameters
getbeta2 = function(compartments, pars)
{
  
  # Extract Compartments
  S <- compartments[["S"]]
  len <- length(S) - 1
  Ip <- compartments[["Ip"]][len]
  Ia <- compartments[["Ia"]][len]
  Ii <- compartments[["Ii"]][len]
  It1 <- compartments[["It1"]][len]
  It2 <- compartments[["It2"]][len]
  In <- compartments[["In"]][len]
  dSdt <- diff(S)
  S <- S[len]
  
  # Extract Parameters
  h <- pars["h"]
  i <- pars["i"]
  j <- pars["j"]
  Nv <- pars["Nv"]
  
  (Nv / S) * dSdt * (-1)/sum(Ip, h*Ia, i*Ii, It1, j*It2, In)
  
}