## functions to simulate the SEIR and SEIcIscR outbreak

## load data: Population age structure 

source("code/getbeta.R")
source('code/1_loadData.r')

#### Loading in James Gleason's data

## Adapted from https://gist.github.com/jonocarroll/b17ce021b0637a31f584ed08a1fbe733
read.tcsv = function(file, header=TRUE, sep=",", ...) {
  
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)[-1]
  
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep) 
  ## empty strings are converted to NA
  out = read.csv(text=x, sep=sep, header=header, na.strings = "", ...)
  return(out)
  
}

jg_dat <- read.tcsv("data/dat_seir_code.csv")
jg_dat$Date <- as.Date(jg_dat$Date, format = "%a %d %b %Y") # Reformat date column

####

loadPopInfo = function(POP)
{
  pop = list()
  pop$N = sum(POP$popage)    # Population size 
  pop$p_age = POP$propage    # Population age structure of China
  return(pop)
}

loadInterventions = function(p_workopen)
{
  list(
    # constraints under a DO-NOTHING scenario 
    base =list(home = diag(1,16,16),
               work = diag(1,16,16),
               school = diag(1,16,16),
               others = diag(1,16,16)),
    # Dublin's lockdown--assume from XX march to May
    wuhanlockdown = list(home = diag(1,16,16),
                         work = diag(0.1,16,16),
                         school = diag(0,16,16),
                         others = diag(c(rep(0.1,4),rep(0.1,12)))),
    # constraints under school closure + some social distancing for school-age going children but 100% workplace
    schcloseonly = list(home = diag(c(rep(1,4),rep(1,12))),
                        work = diag(1,16,16),
                        school = diag(0,16,16),
                        others = diag(c(rep(0.5,4),rep(1,12)))), 
    # constraints under work place distancing only (MAYBE UNREALISTIC, should close schools too)
    workplacedistonly = list(home = diag(1,16,16),
                             work = diag(0.5,16,16),
                             school = diag(1,16,16),
                             others = diag(0.1,16,16)) ,
    
    # constraints under work place distancing + schoolclosure 
    schcloseworkplacedist = list(home = diag(1,16,16),
                                 work = diag(p_workopen,16,16),
                                 school = diag(0,16,16),
                                 others = diag(c(rep(0.1,4),rep(0.1,12)))),
    
    # Post Outbeak, people still cautious 
    postoutbreak = list(home = diag(1,16,16),
                        work = diag(1,16,16),
                        school = diag(1.0,16,16),
                        others = diag(c(rep(1.0,4),rep(1.0,12)))))
  
}



# Children less infectious and as susceptible 
simulateOutbreakSEIcIscR = function(R0t = 3.5,
                                    R0tpostoutbreak = 1.5,
                                    pWorkOpen = c(0.1,0.25,0.5,0.9), # pWorkOpen: proportion of the work force that is working (will be time-varying)
                                    dateStartSchoolClosure = as.Date('2020-03-12') , # cause winter term break 
                                    dateStartIntenseIntervention = as.Date('2020-03-27') , #Intense intervention: lockdown
                                    dateEndIntenseIntervention = as.Date('2020-05-18'), #date we begin relaxing intense intervention 
                                    dateStart = as.Date('2020-02-28'),
                                    POP = dubpop,
                                    numWeekStagger = c(4,8,12),
                                    pInfected = 0.00002,
                                    contacts_ireland = contacts,
                                    tmax = 365*2,
                                    dt = 1,
                                    dat)
{
  # debug dateStartIntenseIntervention = as.Date('2020-01-23')  
  # debug dateEndIntenseIntervention = as.Date('2020-03-01')
  # debug R0est = rep(2,3660) 
  # debug rho = rep(0.8,3660) 
  # debug pWorkOpen =  c(0.1,0.25,0.5,1)
  
  # Load population information
  pop = list()
  pop$N = sum(POP$popage)
  pop$p_age = POP$propage
  N_age = POP$popage                                # Population age structure (in numbers)
  
  # Specify epi info
  
  pars <- c(4.9, 5.9, 7.0, 0.25, 0.05, 0.05, 0.5, 0.75, 0.13, 3.6)
  names(pars) <- c("L","Cv","Dv","h","i","j","f","tv","q","TT")
  
  L <- pars["L"]
  Cv <- pars["Cv"]
  Dv <- pars["Dv"]
  h <- pars["h"]
  i <- pars["i"]
  j <- pars["j"]
  f <- pars["f"]
  tv <- pars["tv"]
  q <- pars["q"]
  TT <- pars["TT"]
  
  
  numSteps = tmax/dt;  	                                 # Total number of simulation time steps
  dateEnd = dateStart+(tmax-1)
  
  # Declare the state variables and related variables:
  # The values of these variables change over time
  # S = E = Isc = Ic = R = array(0,c(numSteps,length(pop$p_age)))
  
  S = Ev = Ip = IA = Ii = It = Iti = Iq = R = Cr <- array(0,c(numSteps,length(pop$p_age)))
  
  lambda = incidence = cumulativeIncidence <- array(0,c(numSteps,length(pop$p_age)))
  
  time = array(0,numSteps)
  
  # Initialise the time-dependent variables, i.e. setting the values of the variables at time 0
  #Ev[1,] = pInfected*5*sum(N_age)/16
  Ev[1,] = 22/16#dat$Exposed[1]/16
  #Ip[1,] = pInfected*sum(N_age)/16
  Ip[1,] = 1.1222/16 #dat$Infected[1]/16
  IA[1,] = 0
  Ii[1,] = 0
  It[1,] = 0
  Iti[1,] = 0
  Iq[1,] = 0
  #R[1,] = 0
  R[1,] = dat$Recovered[1]/16
  Cr[1,] = pInfected*sum(N_age)/16 # Initial number of infected people # Assign 100 infected person in each age group (TODO RELAX?)
  S[1,] = N_age - Ev[1,] - Ip[1,] -  IA[1,] - Ii[1,] - It[1,] - Iti[1,] - Iq[1,] - R[1,]
  
  incidence[1,] = 0;
  
  time[1] = 0;

  ## INTERVENTIONS 
  # School closed 2020-03-12, 
  #lockdown (intense intervention) started 2020-03-27, end of intense intervention: user-specified 
  # note that intense intervention is time-varying control by pWorkOpen: proportion of the work force that is working
  # debug pWorkOpen = c(0.1,0.25,0.5,1)
  tStartSchoolClosure = as.vector(dateStartSchoolClosure - dateStart)+1   #Time to add the school closure effect
  tStartIntenseIntervention = as.vector(dateStartIntenseIntervention - dateStart)+1 # or pw = 0.1
  tEndIntenseIntervention = as.vector(dateEndIntenseIntervention - dateStart)+1     # for pw = 0.1
  tRelaxIntervention1 = tEndIntenseIntervention + numWeekStagger[1]*7                               # for pw = 0.25
  tRelaxIntervention2 = tEndIntenseIntervention + numWeekStagger[2]*7                               # for pw = 0.5
  tRelaxIntervention3 = tEndIntenseIntervention + numWeekStagger[3]*7                               # for pw = 1
  # tStartEndClosure = as.vector(dateEndSchoolClosure - dateStart)+1
  pwork = array(1,numSteps)
  pwork[1:tRelaxIntervention3] =c(rep(1,(tStartIntenseIntervention-0)), # dont know there is outbreak 
                                  rep(pWorkOpen[1],(tEndIntenseIntervention-tStartIntenseIntervention)),
                                  rep(pWorkOpen[2],(tRelaxIntervention1-tEndIntenseIntervention)),
                                  rep(pWorkOpen[3],(tRelaxIntervention2-tRelaxIntervention1)),
                                  rep(pWorkOpen[4],(tRelaxIntervention3-tRelaxIntervention2)))
  
  R0tpostoutbreak = R0t #overwrites the default reduction in R0 post-outbreak 
  # bring r0t to levels pre lockdown so that the only thing that changes is the contacts

  ################################################################################################################### 

  constraintsIntervention = loadInterventions(p_workopen = pwork[1]) # Put here because we need it for the beta function to work, not sure if it actually makes sense? - Daniel
  beta <- getbeta(R0t = R0t, constraints = constraintsIntervention$base, pars = pars, p_age = POP$popage, CONTACTMATRIX = contacts_ireland)
  
  if (length(beta) < numSteps) {
    warning("Not enough R0t values for full time granularity. Elements of R0t will be vectorised.")
    beta <- rep(beta, length.out = numSteps)
  }
  
  if (length(beta) > numSteps) {
    warning("Too many R0t values for time granularity. Vector will be truncated.")
    beta <- beta[1:numSteps]
  }
  
  if(pWorkOpen[2]<1) beta_postfirstwave = getbeta(R0t = R0t,constraints = constraintsIntervention$base,pars = pars,p_age = POP$popage)
  if(pWorkOpen[2]>=1) beta_postfirstwave = beta#getbeta(R0t = R0t[2],constraints = constraintsIntervention$base,gamma = gamma,p_age = pop$p_age)
  
  ###################################################################################################################
  
  ## solve for S, E, I and R at different time points
  
  for (stepIndex in 1: (numSteps-1))
  { 
    
    # load plausible intervetions 
    constraintsIntervention = loadInterventions(p_workopen = pwork[stepIndex])
    
    ## Age- and location-specific contact rates for the given interventions 
    
    # I0: before school winter break intervention period, use base-case
    if(time[stepIndex] < tStartSchoolClosure)  
    {
      CONSTRAINT = constraintsIntervention$base
    }
    # I1:  When school winter break but before lockdown period, use 'schcloseonly'
    if(time[stepIndex] >= tStartSchoolClosure & time[stepIndex] < tStartIntenseIntervention) 
    {
      INTERVENTION = "schcloseonly"   
      CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
    }  
    # I2:  Intense intervention- lock down
    if(time[stepIndex] >= tStartIntenseIntervention & time[stepIndex] < tRelaxIntervention3) 
    {
      INTERVENTION = "schcloseworkplacedist"   
      CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
    }  
    # I3: post outbreak 
    if(time[stepIndex] >= tRelaxIntervention3)  
    {
      CONSTRAINT = constraintsIntervention$postoutbreak
    }
    # 
    
    ## total contacts matrix
    C = CONSTRAINT[[1]]%*%contacts_ireland[[1]]+
      CONSTRAINT[[2]]%*%contacts_ireland[[2]]+
      CONSTRAINT[[3]]%*%contacts_ireland[[3]]+
      CONSTRAINT[[4]]%*%contacts_ireland[[4]]
    
    # calculate the force of infection
    
    # beta = getbeta(R0t = R0t[stepIndex],constraints = constraintsIntervention$base,gamma = gamma,p_age = pop$p_age)
    if(time[stepIndex] < tEndIntenseIntervention+0) lambda[stepIndex,] = as.numeric(beta[stepIndex])*(C%*%(as.matrix(Ip[stepIndex,]) + 
                                                                                                                        (h*as.matrix(IA[stepIndex,])) + 
                                                                                                                        (i*as.matrix(Ii[stepIndex,])) + 
                                                                                                                        (as.matrix(It[stepIndex,])) + 
                                                                                                                        (j*as.matrix(Iti[stepIndex,])) + 
                                                                                                                        (as.matrix(Iq[stepIndex,]))))/N_age
    
    if(time[stepIndex] >= tEndIntenseIntervention+0)lambda[stepIndex,] = as.numeric(beta_postfirstwave[stepIndex])*(C%*%(as.matrix(Ip[stepIndex,]) + 
                                                                                                                                      (h*as.matrix(IA[stepIndex,])) + 
                                                                                                                                      (i*as.matrix(Ii[stepIndex,])) + 
                                                                                                                                      (as.matrix(It[stepIndex,])) + 
                                                                                                                                      (j*as.matrix(Iti[stepIndex,])) + 
                                                                                                                                      (as.matrix(Iq[stepIndex,]))))/N_age
    # calculate the number of infections and recoveries between time t and t+dt
    # Derivatives
    dSdt <- -lambda[stepIndex,]*(S[stepIndex,])#-(beta*S[stepIndex,]*(Ip[stepIndex,] + h*IA[stepIndex,] + i*Ii[stepIndex,] + It[stepIndex,] + j*Iti[stepIndex,] + Iq[stepIndex,]))/N_age
    
    dEvdt <- -(Ev[stepIndex,]/L) + (lambda[stepIndex,]*(S[stepIndex,]))#(beta*S[stepIndex,]*(Ip[stepIndex,] + h*IA[stepIndex,] + i*Ii[stepIndex,] + It[stepIndex,] + j*Iti[stepIndex,] + Iq[stepIndex,]))/N_age
    
    dIpdt <- (((1 - f)*Ev[stepIndex,])/L) - (Ip[stepIndex,]/(Cv - L))
    
    dIAdt <- ((f*Ev[stepIndex,])/L) - (IA[stepIndex,]/Dv)
    
    dIidt <- ((q*Ip[stepIndex,])/(Cv - L)) - (Ii[stepIndex,]/(Dv - Cv + L))
    
    dItdt <- ((tv*Ip[stepIndex,])/(Cv - L)) - (It[stepIndex,]/TT)
    
    dItidt <- (It[stepIndex,]/TT) - (Iti[stepIndex,]/(Dv - Cv + L - TT))
    
    dIqdt <- (((1 - q - tv)*Ip[stepIndex,])/(Cv - L)) - (Iq[stepIndex,]/(Dv - Cv + L))
    
    dRdt <- (IA[stepIndex,]/Dv) + (Ii[stepIndex,]/(Dv - Cv + L)) + (Iq[stepIndex,]/(Dv - Cv + L)) + (Iti[stepIndex,]/(Dv - Cv + L - TT))
    
    dCrdt <- It[stepIndex,]/TT
    
    S[stepIndex + 1,] = S[stepIndex,] + (dSdt*dt)
    Ev[stepIndex + 1,] = Ev[stepIndex,] + (dEvdt*dt)
    Ip[stepIndex + 1,] = Ip[stepIndex,] + (dIpdt*dt)
    IA[stepIndex + 1,] = IA[stepIndex,] + (dIAdt*dt)
    Ii[stepIndex + 1,] = Ii[stepIndex,] + (dIidt*dt)
    It[stepIndex + 1,] = It[stepIndex,] + (dItdt*dt)
    Iti[stepIndex + 1,] = Iti[stepIndex,] + (dItidt*dt)
    Iq[stepIndex + 1,] = Iq[stepIndex,] + (dIqdt*dt)
    R[stepIndex + 1,] = R[stepIndex,] + (dRdt*dt)
    Cr[stepIndex + 1,] = Cr[stepIndex,] + (dCrdt*dt)
    
    incidence[stepIndex+1,] = Ip[stepIndex+1,];     # Only the clinical cases are included in the incidence per day
    time[stepIndex+1] = time[stepIndex] + dt;
  }
  
  output = list( S = S, Ev = Ev, Ip = Ip, IA = IA, Ii = Ii, It = It, Iti = Iti, Iq = Iq, R = R,
                 Cr = Cr, time = time, lambda=lambda, incidence = incidence, N_age= N_age, R0t = R0t,#rho = rho,
                 dateStart = dateStart, dateEnd = dateEnd, dateStartIntenseIntervention = dateStartIntenseIntervention, 
                 dateEndIntenseIntervention = dateEndIntenseIntervention, dateStartSchoolClosure = dateStartSchoolClosure,
                 beta = beta, lambda = lambda)
  
  return(output)
}

CHECKMODEL  = TRUE

if(CHECKMODEL)
{
  # Quick checks: Simulate an outbreak for sanity checks
  set.seed(666)
  
  
  # test for an R0 value of 2.2
  #R0est <- 2.3
  R0est <- jg_dat$Rt
  # R0est = sample(x = r0posterior,size = 100)
  
  ### Set other function parameters
  dateStart <- jg_dat$Date[1]
  
  POP <- dubpop
  ### Fix population
  Nv <- sum(POP$popage)
  POP$popage[1] <- sum(POP$popage[1], 4.9e6 - Nv)
  sum(POP$popage)
  
  tmax <- nrow(jg_dat)
  dt <- 0.1
  
  R0_int <- approx(R0est, xout = seq(1, tmax, dt))$y
  R0_int <- c(R0_int, rep(R0_int[length(R0_int)], (1/dt)-1))
  
  nsim = 1
  epi_doNothing = vector('list',nsim)
  epi_base = vector('list',nsim)
  
  cc <- rep(list(matrix(0, nrow = 16, ncol = 16)), 5)
  #cc[[4]] <- matrix(1, nrow = 16, ncol = 16)
  cc[[4]] <- diag(16)
  names(cc) <- c("home", "work", "school", "others", "all")
  
  for(sim in 1:nsim)
  {
    
    epi_doNothing[[sim]] = simulateOutbreakSEIcIscR(R0t =R0_int,dateStartSchoolClosure = as.Date('2020-03-10'),
                                                    dateStartIntenseIntervention =as.Date('2020-03-10'),
                                                    dateEndIntenseIntervention = as.Date('2020-03-10'),
                                                    pWorkOpen = c(1,1,1,1),numWeekStagger = c(0,0,0,0),
                                                    dateStart = dateStart, POP = POP, tmax = tmax,
                                                    contacts_ireland = cc, pInfected = 0.000001, dt = dt,
                                                    dat = jg_dat)
    
    #epi_base[[sim]] = simulateOutbreakSEIcIscR(R0t = R0est ,dateEndIntenseIntervention = as.Date('2020-05-18'))
    
    
  }
  par(mfrow=c(2,1))
  
  # incidence over time
  agegp = 10
  plot(epi_doNothing[[1]]$time, epi_doNothing[[1]]$incidence[,agegp], type='l', lwd=2,
       main=paste0("Incidence for age [",(agegp-1)*5,',',agegp*5,')'),
       xlab="Time(days)", ylab="Daily no. of infections");
  #lines(x=epi_base[[1]]$time,y=epi_base[[1]]$incidence[,agegp],lwd=2,col='tomato')
  
  
  
  # cumulative incidence over time
  plot(epi_doNothing[[1]]$time, (epi_doNothing[[1]]$N_age[agegp]-epi_doNothing[[1]]$S[,agegp])/epi_doNothing[[1]]$N_age[agegp], lwd=2,type='l', 
       main=paste0("Cum incidence for age [",(agegp-1)*5,',',agegp*5,')'),
       xlab="Time(days)", ylab="Cum incidence",ylim = c(0,1));
  #lines(epi_base[[1]]$time, (epi_base[[1]]$N_age[agegp]-epi_base[[1]]$S[,agegp])/epi_base[[1]]$N_age[agegp],lwd=2,col='tomato')
  #legend(0.25, 0.98, legend=c("Do Nothing", "Base"),
  #       col=c("black", "tomato"), bty='n',lty=c(1,1),lwd=c(2,2), cex=0.7)
  
  
  # agegp = 4
  # plot(epi_doNothing[[1]]$time, epi_doNothing[[1]]$S[,agegp], type='l', lwd=2,
  #      main=paste0("Incidence for age [",(agegp-1)*5,',',agegp*5,')'),
  #      xlab="Time(days)", ylab="Daily no. of infections");
  # lines(x=epi_base[[1]]$time,y=epi_base[[1]]$S[,agegp],lwd=2,col='tomato')
  # lines(x=epi_base[[1]]$time,y=epi_base[[1]]$incidence[,agegp],lwd=2,col='grey')
  # 
}

### Validation
par(mfrow = c(1, 1))

# Beta
plot(jg_dat$Beta, type = "l")
lines(xx, epi_doNothing[[1]]$beta, col = "red")

# Infected
Nv <- sum(POP$popage)
I <- Nv - rowSums(epi_doNothing[[1]]$S)- rowSums(epi_doNothing[[1]]$Ev)- rowSums(epi_doNothing[[1]]$R)
xx <- seq_along(I)*dt
plot(jg_dat$Infected, type = "l")
lines(xx, I, col = "red")

# Susceptiple
plot(jg_dat$Susceptible, type = "l")
lines(xx, rowSums(epi_doNothing[[1]]$S), col = "red")

# Recovered
plot(jg_dat$Recovered, type = "l")
lines(xx, rowSums(epi_doNothing[[1]]$R), col = "red")

#FOI
plot(xx, rowSums(epi_doNothing[[1]]$lambda), col = "red", type = "l")
lines(jg_dat$ForceInfection, type = "l")

#all.equal(jg_dat$Infected , I)

