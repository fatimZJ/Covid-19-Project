## functions summarising the simulations 

# total number per day (across all agegroups)
summariseSimulations = function(VAR,CI,SIMS)
{
  temp = lapply(SIMS,FUN = function(x) rowSums(x[[VAR]])) # total per day per simulation run list of length #sim runs
  temp1 = do.call(cbind.data.frame, temp)  # cbind total per day per simulation- dataframe of size #days x #sim runs
  var_p_median = apply(temp1,1,function(x) quantile(x,0.5)) # median of the simulation runs for each day
   var_p_lci = apply(temp1,1,function(x) quantile(x,(1-CI/100)/2)) # confidence bounds
  var_p_uci = apply(temp1,1,function(x) quantile(x,1-(1-CI/100)/2))
  SUMMARY = list(median = var_p_median,lci = var_p_lci,uci=var_p_uci)
  rm(temp,temp1,var_p_median,var_p_lci,var_p_uci)
  
  
  results = list(summary=SUMMARY, Sim1 = SIMS[[1]])
  
  return(results)
  
}
# get peak incidence per day and the corresponding t
summarisePeakTimePeakSize= function(VAR = 'incidence',SIMS)
{
  time = SIMS[[1]]$time 
  temp = lapply(SIMS,FUN = function(x) rowSums(x[[VAR]])) #total incidence per day per simulation run
  temp1 = do.call(cbind.data.frame, temp) # cbind simulation runs of total per day- dataframe of size #days x #sim runs
  peaksize = as.numeric(apply(temp1,2,max)) # identify the peak (max total incidence in a day) in each simulation run
  peaktime = time[as.numeric(apply(temp1,2,function(x) which.max(x)))] # find t which corresponds to the peak
  results = list(time = time, peaktime = peaktime,peaksize = peaksize)
  return(results)
}

# "Median" simulation! 
summariseSimulations_mid = function(CI,SIMS)
{
  temp = lapply(SIMS,FUN = function(x) rowSums(x[['S']])) #total S per day per simulation run list of length #sim runs
  temp1 = do.call(cbind.data.frame, temp) # cbind simulation runs of total S per day- dataframe of size #days x #sim runs
  i = which.min(abs(as.numeric(temp1[nrow(temp1),])-as.numeric(quantile(temp1[nrow(temp1),],0.5)))) # median of total S per day in last sim
  j = which.min(abs(as.numeric(temp1[nrow(temp1),])-as.numeric(quantile(temp1[nrow(temp1),],0.25)))) # Confidence bounds
  k = which.min(abs(as.numeric(temp1[nrow(temp1),])-as.numeric(quantile(temp1[nrow(temp1),],0.75))))
  
  S = data.frame(med = temp[[i]],lci = temp[[k]],uci = temp[[j]],time = SIMS[[1]]$time) # get "median" simulation + ci total S
  S_age = list(med = SIMS[[i]]$S,lci = SIMS[[k]]$S,uci = SIMS[[j]]$S) # median simulation results- S per age group
  inc = list(med = SIMS[[i]]$incidence,lci = SIMS[[k]]$incidence,uci = SIMS[[j]]$incidence) # median simulation results- Incidence per age group
  time = SIMS[[1]]$time # time
  N_age = SIMS[[1]]$N_age # total in each age group
  results = list(S=S,inc = inc, time = time,N_age=N_age,S_age = S_age)
  return(results)
}

# median value per day over all the simulation runs for each age group
summariseSimulationsAGE = function(VAR,CI,SIMS)
{
  var_p_median = var_p_lci = var_p_uci = array(NA,c(length(SIMS[[1]]$time),16))
  for(age in 1:16)
    age = 1
  {
    temp = lapply(SIMS,FUN = function(x) (x[[VAR]][,age])) # select incidence or S per day for each age group from every simulation
    temp1 = do.call(cbind.data.frame, temp) # cbind the simulation runs per age group
    var_p_median[,age] = apply(temp1,1,function(x) quantile(x,0.5)) # get the median for everyday
    var_p_lci[,age] = apply(temp1,1,function(x) quantile(x,(1-CI/100)/2)) # confidence intervals
    var_p_uci[,age] = apply(temp1,1,function(x) quantile(x,1-(1-CI/100)/2))
    rm(temp,temp1)
  }
  
  SUMMARY = list(median = var_p_median,lci = var_p_lci,uci=var_p_uci)
  rm(var_p_median,var_p_lci,var_p_uci)
  
  
  results = list(summary=SUMMARY, Sim1 = SIMS[[1]])
  return(results)
  
}